%% ========================================================================
%  AD9361 "Catalina" Fast-Lock LUT generator
%
%  Trained min-max NUFFT  ->  Chebyshev nodes  ->  barycentric polynomial
%  Project p-2026-146, Ben-Gurion University.
%
%  PIPELINE
%    1  LOAD     measurements, keep the 1 MHz band (>= cfg.fStart)
%    2  SCAN     split the band at the RFPLL divider boundaries
%    3  TRAIN    per segment, optimise the scaling vector s
%    4  BUILD    evaluate at Chebyshev nodes, store the polynomial
%    5  VALIDATE worst-case error per field against the 0.5-count criterion
%    6  LUT      round, clip to bit width, write the Fast-Lock profile file
%
%  WHY TRAINING IS THE RIGHT THING HERE
%  ------------------------------------
%  The NUFFT has a real degree of freedom: the scaling vector
%       s_n = sum_{t=-L}^{L} alpha_t exp(i*gamma*beta*t*(n-n0))        (eq. 28)
%  Fessler & Sutton choose (alpha,beta) to minimise the worst case over ALL
%  unit-norm signals and ALL frequencies, because in MRI the trajectory is not
%  known in advance.  Here BOTH are known at build time: the signal is the
%  measured register column, and the query points are the M Chebyshev nodes.
%  Training against that specific pair is a strictly easier problem than the
%  paper solves, and it wins by a wide margin.  Measured on one segment, error
%  against the exact NDFT at the Chebyshev nodes, in counts:
%
%        J      Table 2 alpha      trained alpha      gain
%        4        2.7e-02            4.3e-03          6.4x
%        6        6.1e-03            5.0e-04           12x
%        8        1.7e-03            5.4e-05           31x
%       10        2.2e-04            1.1e-05           20x
%
%  L = 2 captured all of it (L = 4 and L = 6 gave nothing further), so training
%  searches only three numbers: beta, alpha_1, alpha_2.  The paper's Table 2
%  value is the starting point, and a trained result is rejected unless it
%  beats it - so the published worst-case guarantee is a floor, never a risk.
%
%  THE REFERENCE THE TRAINING USES
%  -------------------------------
%  Loss = max_m | NUFFT(t_m) - NDFT_exact(t_m) | over the Chebyshev nodes.
%  The exact NDFT is the same trigonometric interpolant the NUFFT approximates,
%  evaluated the slow O(N*M) way.  This isolates the NUFFT approximation from
%  every other error term, and because M is small the reference is cheap.
%
%  Two traps that silently produced a zero loss in earlier versions:
%    - Do not train at half-integer sample coordinates.  With K = 2N,
%      omega/gamma = -2t, so half-integer t lands exactly on an oversampled
%      FFT bin where the interpolator is already exact.  A train/test split on
%      odd/even samples does exactly this, which is why the old J sweep was
%      flat at 1e-11.  The Chebyshev nodes are generic, so they are safe.
%    - Odd N.  See the note on the final phase in nufft_engine.
%
%  UNITS: errors are in COUNTS of the field concerned, not volts or dB.
% ========================================================================

function Catalina_Project()
    clc;
    disp('=============================================================');
    disp('  AD9361 Catalina - Fast-Lock LUT generator');
    disp('  trained min-max NUFFT -> segmented Chebyshev polynomials');
    disp('=============================================================');
    disp(' 0. RUN ALL   (scan + train + build + validate)');
    disp(' 1. SCAN      segments and field types');
    disp(' 2. BUILD     train the scaling vectors and fit the model');
    disp(' 3. VALIDATE  error report against the 0.5-count criterion');
    disp(' 4. GENERATE  Fast-Lock LUT for target frequencies');
    disp(' 5. ANALYSE   J sweep, M sweep, training gain');
    disp(' 6. Exit');
    disp('=============================================================');

    switch input('Select (0-6): ')
        case 0
            m = ask('Measurement file', 'measurements.xlsx');
            f = ask('Model file', 'model.mat');
            scan_report(m); build_model(m, f); validate_model(f, m);
        case 1, scan_report(ask('Measurement file', 'measurements.xlsx'));
        case 2, build_model(ask('Measurement file', 'measurements.xlsx'), ...
                            ask('Model file', 'model.mat'));
        case 3, validate_model(ask('Model file', 'model.mat'), ...
                               ask('Measurement file', 'measurements.xlsx'));
        case 4, generate_profiles(ask('Model file', 'model.mat'), ...
                                  ask('Target frequency file', 'targets.xlsx'), ...
                                  ask('Output profile file', 'profiles.xlsx'));
        case 5, error_analysis(ask('Measurement file', 'measurements.xlsx'));
        case 6, disp('Bye.');
        otherwise, disp('Invalid choice.');
    end
end

function s = ask(prompt, default)
    s = input(sprintf('%s [%s]: ', prompt, default), 's');
    if isempty(s), s = default; end
end

%% ========================================================================
%  CONFIGURATION
% ========================================================================
function cfg = default_config()
    % ---- band ---------------------------------------------------------
    cfg.fStart     = [];          % MHz, or [] to keep the whole sweep.
                                  % The sweep is 0.5 MHz below 700 MHz and
                                  % 1 MHz above; both are usable, because
                                  % segmentation never lets a segment straddle
                                  % the rate change.  Set a number only to
                                  % deliberately ignore the low band.
    % ---- segmentation -------------------------------------------------
    cfg.segMode    = 'control';   % 'control' | 'divider' | 'manual' | 'none'
    cfg.controlReg = "005[7:4]";  % field carrying the divider exponent
    cfg.manualEdges= [];          % MHz, used by 'manual'
    cfg.fTop       = 6000;        % 'divider' mode: edges at fTop/2^k
    cfg.minSegPts  = 32;          % a run shorter than this is folded into the
                                  % previous segment - but only across a SOFT
                                  % (jump) edge, never across a rate change or
                                  % a band edge.  The lowest divider bands are
                                  % narrow, so this is deliberately small.

    % ---- jump splitting (second level) --------------------------------
    cfg.jumpTol    = 10;          % counts. A sample-to-sample step larger than
                                  % this is a discontinuity: break the band and
                                  % fit a separate polynomial on each side.
    cfg.jumpTolPerField = containers.Map('KeyType','char','ValueType','double');
    % Override for fields whose sample-to-sample scatter is comparable to the
    % threshold.  276[6:0] scatters by roughly +/-10 counts around each tooth,
    % so a flat tolerance of 10 fragments it on noise; its resets are far
    % larger, so a higher value catches them cleanly.  Uncomment to use:
    %   cfg.jumpTolPerField('276[6:0]') = 40;
    cfg.jumpFields = 'smooth';    % 'smooth' -> only fields that would get a
                                  %             polynomial can trigger a split
                                  % 'all'    -> any non-exact field can
    cfg.maxSegments= 200;         % safety cap; exceeding it means the field is
                                  % not piecewise smooth and needs a formula

    % ---- field handling -----------------------------------------------
    % Fields computed exactly from the target frequency by pll_registers.m.
    % Never interpolated, never trained: the fractional word steps by roughly
    % 400 000 counts per MHz and wraps every ~20 MHz, so its low bytes are
    % arithmetic wraparound, not a sampled smooth function.
    cfg.formulaFields = ["005[7:4]", "271[7:0]", "272[2:0]", ...
                         "273[7:0]", "274[7:0]", "275[6:0]"];
    cfg.fref        = 40;         % MHz - SET THIS to your board's reference
    cfg.vcoMin      = 6000;       % MHz, lower edge of the VCO range
    cfg.exactFields = strings(0); % stored as transition lists: exact at the
                                  % measured points, piecewise constant between
    cfg.maxLevels   = 12;         % <= this many distinct values -> lookup

    % ---- composite fields ----------------------------------------------
    % Some values are wider than one register field and are split across two.
    % Table 12 word E is "Force VCO Tune[7:0]" and word F carries "Force VCO
    % Tune[8]", so the tune word is 9 bits living in 277[7:0] and 278[0].
    % Modelling the LOW BYTE alone is hopeless: it wraps 255 -> 0 whenever the
    % 9-bit value crosses 256, and that wrap is arithmetic, not physical.  The
    % composite has no wrap, so it is modelled instead and the member fields
    % are extracted from it by shifting.
    cfg.compositeFields = { ...
        struct('name', "VCO_Tune", ...
               'parts', ["277[7:0]", "278[0]"], ...
               'shift', [0 8]) };

    % ---- adaptive refinement -------------------------------------------
    cfg.refine      = true;       % bisect a segment when a field fails in it
    cfg.maxRefine   = 6;          % refinement passes
    cfg.refineMinPts= 48;         % never bisect below this many points

    % ---- residual correction (last resort) -----------------------------
    % Stores the integer difference between the model and the measurement
    % wherever it survives rounding.  This makes the field exact AT THE
    % MEASURED FREQUENCIES and says nothing about the points between them.
    % It is off by default because a field that needs it is telling you the
    % data is not a function of frequency - see the report note.
    cfg.residual    = false;
    cfg.residualMaxFrac = 0.05;   % refuse if it would store more than this
                                  % fraction of the raw samples

    % ---- NUFFT --------------------------------------------------------
    cfg.J          = 10;
    cfg.padding    = 2;           % Table 2 assumes K/N = 2
    cfg.L          = 2;           % scaling-vector order
    cfg.detrend    = true;

    % ---- training -----------------------------------------------------
    cfg.train      = true;
    cfg.trainScope = 'segment';   % 'segment' -> one s per segment (fast)
                                  % 'field'   -> one s per field per segment
    cfg.maxFevals  = 400;

    % ---- polynomial ---------------------------------------------------
    cfg.deg        = 120;         % Chebyshev nodes per segment
    cfg.tol        = 0.5;         % counts
end

%% ========================================================================
%  1. LOAD
% ========================================================================
function [freq, names, X, widths] = load_measurements(xlsFile, cfg)
    if ~isfile(xlsFile), error('Measurement file "%s" not found.', xlsFile); end

    opts = detectImportOptions(xlsFile);
    opts.VariableNamingRule = 'preserve';
    T = readtable(xlsFile, opts);

    names = string(T{:, 1});
    hdr   = string(T.Properties.VariableNames(2:end));
    freq  = str2double(regexprep(hdr, '[^\d\.\-eE+]', ''));
    if any(isnan(freq)), error('Header row 2..end must be numeric frequencies.'); end

    X = T{:, 2:end}.';
    freq = freq(:);
    [freq, ord] = sort(freq);  X = X(ord, :);

    if any(isnan(X(:)))
        warning('NaNs in the body - filling with nearest.');
        X = fillmissing(X, 'nearest');
    end

    if nargin > 1 && ~isempty(cfg.fStart)
        keep = freq >= cfg.fStart;
        if any(~keep)
            fprintf('Dropping %d point(s) below %g MHz (cfg.fStart).\n', ...
                sum(~keep), cfg.fStart);
        end
        freq = freq(keep);  X = X(keep, :);
    end

    % ---- sampling-rate regions ------------------------------------------
    % The sweep is measured at more than one step: 0.5 MHz over 70-700 MHz and
    % 1 MHz above.  The FFT stage needs a uniform grid, but only WITHIN a
    % segment, so mixed rates are fine as long as no segment spans a rate
    % change.  scan_segments enforces that; here we only measure and report.
    R = rate_regions(freq);
    fprintf('Loaded %s : %d frequencies (%.4g - %.4g MHz), %d fields\n', ...
        xlsFile, numel(freq), freq(1), freq(end), numel(names));
    fprintf('Sampling-rate region(s):\n');
    for k = 1:numel(R)
        fprintf('   %8.3f - %8.3f MHz   step %g MHz   %5d point(s)\n', ...
            freq(R(k).i0), freq(R(k).i1), R(k).step, R(k).i1 - R(k).i0 + 1);
        if ~R(k).uniform
            error(['Region %d (%.4g - %.4g MHz) is not uniformly sampled ' ...
                   '(steps %g .. %g MHz).\nThe FFT stage cannot be applied ' ...
                   'to it.'], k, freq(R(k).i0), freq(R(k).i1), ...
                   R(k).dmin, R(k).dmax);
        end
    end

    widths = register_bit_width(names);
end

function R = rate_regions(freq)
%RATE_REGIONS  Split the frequency axis into maximal uniformly-sampled runs.
%
% Returns a struct array with the first and last index of each run, its step,
% and whether the run really is uniform to within a tolerance.  A run boundary
% is a place where the step changes, e.g. 0.5 -> 1 MHz at 700 MHz.
    freq = freq(:);
    d = diff(freq);
    if isempty(d)
        R = struct('i0',1,'i1',1,'step',NaN,'uniform',true,'dmin',NaN,'dmax',NaN);
        return;
    end

    % A step change is relative: 0.5 -> 1 is a change, 1.0 -> 1.0000001 is not.
    tol   = 1e-6;
    brk   = find(abs(diff(d)) > tol * max(abs(d)));
    edges = unique([1; brk(:) + 2; numel(freq) + 1]);

    R = struct('i0',{}, 'i1',{}, 'step',{}, 'uniform',{}, 'dmin',{}, 'dmax',{});
    for k = 1:numel(edges) - 1
        i0 = edges(k);  i1 = edges(k+1) - 1;
        if i1 <= i0, continue; end
        dd = diff(freq(i0:i1));
        R(end+1).i0  = i0;                                        %#ok<AGROW>
        R(end).i1    = i1;
        R(end).step  = median(dd);
        R(end).dmin  = min(dd);
        R(end).dmax  = max(dd);
        R(end).uniform = (max(dd) - min(dd)) <= tol * max(abs(dd));
    end

    % A single stray point between two runs shows up as a one-sample region;
    % fold it into the run on its left so it is not treated as its own rate.
    keep = true(numel(R),1);
    for k = 2:numel(R)
        if R(k).i1 - R(k).i0 < 1
            R(k-1).i1 = R(k).i1;  keep(k) = false;
        end
    end
    R = R(keep);
end

function w = register_bit_width(names)
%REGISTER_BIT_WIDTH  "272[2:0]" -> 3 bits.  Silently defaulting to 8 is
% dangerous: a 3-bit field clipped to 0..255 overflows its neighbour's bits
% inside the packed setup word, so an unparseable name is an error, not a
% default.  Tolerates the datasheet's 0x272[D2:D0] notation.
    names = string(names(:));
    w = zeros(numel(names), 1);
    for i = 1:numel(names)
        t = regexp(names(i), '\[([^\]]*)\]', 'tokens', 'once');
        if isempty(t)
            error(['Field name "%s" has no [bit:range]. Expected 272[2:0] ' ...
                   'or 272[D2:D0].'], names(i));
        end
        n = str2double(strtrim(split(regexprep(string(t{1}), '[Dd]', ''), ':')));
        if any(isnan(n)) || numel(n) > 2
            error('Cannot read the bit range from "%s".', names(i));
        end
        if numel(n) == 2, w(i) = abs(n(1) - n(2)) + 1; else, w(i) = 1; end
    end
end

%% ========================================================================
%  2. SCAN - segment the band
% ========================================================================
function segs = scan_segments(freq, names, X, cfg, verbose)
    if nargin < 5, verbose = true; end
    N = numel(freq);

    % ---- level 0: sampling-rate changes (MANDATORY) -------------------
    % The FFT stage assumes equally spaced samples.  That only has to hold
    % inside a segment, so a sweep measured at 0.5 MHz below 700 MHz and 1 MHz
    % above is usable in full - provided no segment straddles the change.
    % These edges are not optional and a short run is never merged across one.
    R  = rate_regions(freq);
    e0 = arrayfun(@(r) r.i0, R);
    e0 = e0(e0 > 1);
    e0 = e0(:);

    % ---- level 1: physical band edges ---------------------------------
    switch lower(cfg.segMode)
        case 'control'
            e1 = edges_from_control(freq, names, X, cfg, verbose);
            if isempty(e1)
                if verbose
                    warning('Control field %s not found - using divider edges.', ...
                            cfg.controlReg);
                end
                e1 = edges_from_divider(freq, cfg);
            end
        case 'divider', e1 = edges_from_divider(freq, cfg);
        case 'manual',  e1 = edges_from_list(freq, cfg.manualEdges);
        case 'none',    e1 = [];
        otherwise, error('Unknown segMode "%s".', cfg.segMode);
    end
    e1 = setdiff(e1(:), e0);

    % ---- level 2: jump splitting --------------------------------------
    % A sample-to-sample step larger than cfg.jumpTol counts means the field
    % is discontinuous there, and no single polynomial can span it.  Break
    % the band and fit a separate polynomial on each side.
    [e2, why] = edges_from_jumps(freq, names, X, cfg);
    e2 = setdiff(e2(:), [e0; e1]);

    % Levels 0 and 1 are HARD: crossing one either breaks the FFT assumption
    % or spans a physical discontinuity.  Level 2 is soft.  A run shorter than
    % cfg.minSegPts may be absorbed backwards only across a soft edge.
    hard  = [e0; e1];
    edges = unique([1; e0; e1; e2; N+1]);

    segs   = struct('idx',{}, 'freq',{}, 'f_lo',{}, 'f_hi',{}, 'src',{}, 'step',{});
    isRate = @(a) ismember(a, e0);
    isBand = @(a) ismember(a, e1);
    isJump = @(a) ismember(a, e2);

    for k = 1:numel(edges)-1
        idx = edges(k):edges(k+1)-1;
        if isempty(idx), continue; end

        canMerge = ~isempty(segs) && ~ismember(edges(k), hard);
        if numel(idx) < cfg.minSegPts && canMerge
            segs(end).idx  = [segs(end).idx, idx];
            segs(end).freq = freq(segs(end).idx);
            segs(end).f_hi = freq(segs(end).idx(end));
            continue;
        end

        segs(end+1).idx = idx;                                   %#ok<AGROW>
        segs(end).freq  = freq(idx);
        segs(end).f_lo  = freq(idx(1));
        segs(end).f_hi  = freq(idx(end));
        segs(end).step  = median(diff(freq(idx)));

        src = {};
        if isRate(edges(k)), src{end+1} = 'rate'; end            %#ok<AGROW>
        if isBand(edges(k)), src{end+1} = 'band'; end            %#ok<AGROW>
        if isJump(edges(k)), src{end+1} = 'jump'; end            %#ok<AGROW>
        if isempty(src),     src = {'start'}; end
        segs(end).src = strjoin(src, '+');
    end

    % ---- every segment must be uniformly sampled ------------------------
    for k = 1:numel(segs)
        dd = diff(segs(k).freq);
        if numel(dd) >= 1 && (max(dd) - min(dd)) > 1e-6 * max(dd)
            error(['Segment %d (%.4g - %.4g MHz) is not uniformly sampled ' ...
                   '(steps %g .. %g MHz).\nThe NUFFT cannot be applied to it. ' ...
                   'This should not happen - report it.'], ...
                   k, segs(k).f_lo, segs(k).f_hi, min(dd), max(dd));
        end
    end

    if numel(segs) > cfg.maxSegments
        error(['Segmentation produced %d segments (cap is cfg.maxSegments = %d).\n' ...
               'Fields responsible: %s\n' ...
               'A field that jumps this often is not piecewise smooth - it needs a ' ...
               'closed-form relation, not more segments.  Add it to cfg.exactFields ' ...
               'or raise cfg.jumpTol.'], numel(segs), cfg.maxSegments, ...
               strjoin(unique(why), ', '));
    end

    if verbose
        fprintf(['\n%d segment(s)  (%d rate change(s), %d band edge(s), ' ...
                 '%d jump(s) > %g counts):\n'], ...
            numel(segs), numel(e0), numel(e1), numel(e2), cfg.jumpTol);
        fprintf('%4s %11s %11s %8s %8s %11s %5s\n', ...
            'seg','f_lo','f_hi','points','step','edge from','M');
        for k = 1:numel(segs)
            n = numel(segs(k).idx);
            fprintf('%4d %11.3f %11.3f %8d %8g %11s %5d\n', k, segs(k).f_lo, ...
                segs(k).f_hi, n, segs(k).step, segs(k).src, ...
                min(cfg.deg, floor(n/3)));
        end
        short = arrayfun(@(s) numel(s.idx), segs) < 3*16;
        if any(short)
            fprintf(['%d segment(s) hold fewer than 48 points, so their ' ...
                     'polynomial degree is\ncapped at N/3. That is expected ' ...
                     'next to a rate change or a narrow band.\n'], sum(short));
        end
        if ~isempty(e2)
            fprintf('Jump splits triggered by: %s\n', strjoin(unique(why), ', '));
        end
    end
end

function [e, why] = edges_from_jumps(freq, names, X, cfg)
% Split wherever any candidate field steps by more than its jump tolerance.
% Fields computed by formula, or stored exactly, are skipped: they never get a
% polynomial, so splitting on their behalf achieves nothing.  This matters most
% for the fractional word, whose low bytes jump at almost every sample and
% would otherwise shatter the band into single-point segments.
    e = [];  why = strings(0,1);
    if isempty(cfg.jumpTol) || ~isfinite(cfg.jumpTol) || cfg.jumpTol <= 0
        return;
    end

    for i = 1:numel(names)
        if any(names(i) == cfg.formulaFields), continue; end
        if any(names(i) == cfg.exactFields),   continue; end
        if strcmpi(cfg.jumpFields, 'smooth') && numel(unique(X(:,i))) <= cfg.maxLevels
            continue;      % staircase / lookup field: stored exactly anyway
        end

        tol = cfg.jumpTol;
        if isKey(cfg.jumpTolPerField, char(names(i)))
            tol = cfg.jumpTolPerField(char(names(i)));
        end

        d   = abs(diff(X(:,i)));
        hit = find(d > tol) + 1;
        if isempty(hit), continue; end

        % Separation check: if the steps being called "jumps" are not much
        % bigger than the routine step-to-step variation, the threshold is
        % cutting into noise rather than finding discontinuities.
        typ = median(d(d > 0));
        if ~isempty(typ) && typ > 0 && median(d(d > tol)) < 3*typ
            warning(['%s: jumps above %g counts are only %.1fx the typical step ' ...
                     '(%.1f). This may be splitting on scatter rather than on real ' ...
                     'discontinuities - consider cfg.jumpTolPerField(''%s'').'], ...
                     names(i), tol, median(d(d > tol))/typ, typ, names(i));
        end

        e = union(e, hit);
        why(end+1) = sprintf('%s: %d @ >%g', names(i), numel(hit), tol);  %#ok<AGROW>
    end
    e = e(:);
end

function e = edges_from_control(freq, names, X, cfg, verbose)
    e = [];
    r = find(names == cfg.controlReg, 1);
    if isempty(r), return; end
    v = X(:, r);
    e = find(diff(v) ~= 0) + 1;
    if verbose && ~isempty(e)
        fprintf('Divider field %s: %d level(s) -> boundaries at %s MHz\n', ...
            cfg.controlReg, numel(unique(v)), strjoin(compose('%.0f', freq(e).'), ', '));
    end
end

function e = edges_from_divider(freq, cfg)
    f = cfg.fTop ./ 2.^(0:12);
    f = f(f > freq(1) & f < freq(end));
    e = unique(arrayfun(@(x) find(freq >= x, 1, 'first'), sort(f(:))));
end

function e = edges_from_list(freq, edgesMHz)
    e = [];
    if isempty(edgesMHz), return; end
    edgesMHz = sort(edgesMHz(:));
    edgesMHz = edgesMHz(edgesMHz > freq(1) & edgesMHz < freq(end));
    e = unique(arrayfun(@(x) find(freq >= x, 1, 'first'), edgesMHz));
end

function kind = classify_field(x, name, cfg)
    if any(name == cfg.formulaFields), kind = 'formula'; return; end
    if any(name == cfg.exactFields),   kind = 'exact';   return; end
    u = unique(x);
    if numel(u) == 1,             kind = 'const';  return; end
    if numel(u) <= cfg.maxLevels, kind = 'lookup'; return; end
    kind = 'smooth';
end

function tag = formula_tag(name)
% "271[7:0]" -> "R271".  The register address selects the calculator output.
    t = regexp(name, '^\s*(?:0[xX])?([0-9A-Fa-f]+)\s*\[', 'tokens', 'once');
    if isempty(t)
        error('Cannot read a register address from the field name "%s".', name);
    end
    tag = "R" + upper(string(t{1}));
end

function check_formula_fields(freq, names, X, cfg)
% The formula fields are computed, not fitted, so comparing them against the
% measurements is a free check on cfg.fref.  If the reference frequency is
% wrong, the integer word will be off by a constant factor and this table will
% show it immediately.
    idx = find(ismember(names, cfg.formulaFields));
    if isempty(idx)
        fprintf('\nNo formula fields present in this measurement set.\n');
        return;
    end

    R = pll_registers(freq, cfg.fref, struct('vcoMin', cfg.vcoMin, 'warn', true));

    fprintf('\nFormula check against measurement  (f_ref = %g MHz, VCO min = %g MHz)\n', ...
        cfg.fref, cfg.vcoMin);
    fprintf('%-14s %12s %12s %10s\n', 'Field', 'max |diff|', 'mismatched', 'verdict');
    fprintf('%s\n', repmat('-', 1, 52));
    allOK = true;
    for i = idx(:).'
        v = R.(char(formula_tag(names(i))));
        d = abs(v - X(:, i));
        nBad = sum(d > 0);
        ok = nBad == 0;
        allOK = allOK && ok;
        fprintf('%-14s %12g %12d %10s\n', names(i), max(d), nBad, ...
            string(ternary(ok, 'exact', 'MISMATCH')));
    end
    fprintf('%s\n', repmat('-', 1, 52));
    if allOK
        fprintf('All formula fields reproduce the measurements exactly.\n');
    else
        fprintf(['Mismatch: check cfg.fref (currently %g MHz) and cfg.vcoMin.\n' ...
                 'A wrong reference frequency scales the integer word and shows up\n' ...
                 'here before it ever reaches a profile.\n'], cfg.fref);
    end
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

function scan_report(xlsFile)
    cfg = default_config();
    [freq, names, X, widths] = load_measurements(xlsFile, cfg);
    segs = scan_segments(freq, names, X, cfg);

    nF = numel(names);  K = numel(segs);
    fprintf('\n%-14s %6s  %s\n', 'Field', 'bits', 'kind per segment');
    fprintf('%s\n', repmat('-', 1, 32 + 8*K));
    nSmooth = 0;
    for i = 1:nF
        kinds = arrayfun(@(k) string(classify_field(X(segs(k).idx,i), names(i), cfg)), 1:K);
        nSmooth = nSmooth + sum(kinds == "smooth");
        fprintf('%-14s %6d  %s\n', names(i), widths(i), strjoin(kinds, ' '));
    end
    fprintf(['\nformula = computed from the frequency (pll_registers.m), not fitted\n' ...
             'exact   = stored as a transition list\n' ...
             'const   = one value over the segment\n' ...
             'lookup  = few distinct values, stored as a transition list\n' ...
             'smooth  = fitted by NUFFT + Chebyshev  (%d polynomial(s))\n'], nSmooth);

    check_formula_fields(freq, names, X, cfg);

    % ---- plots: 3 fields per figure, shaded segment bands ---------------
    per = 3;
    for a = 1:per:nF
        b = min(a+per-1, nF);
        figure('Name', sprintf('Scan %d-%d', a, b), 'Color', 'w', ...
               'Position', [60 60 1150 260*(b-a+1)]);
        tl = tiledlayout(b-a+1, 1, 'TileSpacing','compact', 'Padding','compact');
        for i = a:b
            ax = nexttile; hold(ax,'on');
            shade_segments(ax, segs);
            plot(ax, freq, X(:,i), '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 5);
            kinds = arrayfun(@(k) string(classify_field(X(segs(k).idx,i), names(i), cfg)), 1:K);
            title(ax, sprintf('%s   |   %d bits   |   %s', names(i), widths(i), ...
                  strjoin(kinds,' / ')), 'Interpreter','none', 'FontWeight','normal');
            ylabel(ax, 'counts');
            style_axes(ax);
            if i == b, xlabel(ax, 'frequency [MHz]'); else, ax.XTickLabel = []; end
        end
        title(tl, 'Measured register fields with segment boundaries', ...
              'FontWeight','bold', 'FontSize', 13);
    end
end

%% ---- plotting helpers ---------------------------------------------------
function style_axes(ax)
    grid(ax,'on');  box(ax,'on');
    ax.GridAlpha = 0.12;
    ax.FontSize  = 11;
    ax.LineWidth = 0.8;
    ax.Layer     = 'top';
end

function shade_segments(ax, segs)
% Alternating light bands instead of hard vertical lines: the boundaries stay
% legible without competing with the data for attention.
    for k = 1:numel(segs)
        if mod(k,2) == 0, continue; end
        xr = [segs(k).f_lo segs(k).f_hi];
        patch(ax, [xr(1) xr(2) xr(2) xr(1)], [-1e9 -1e9 1e9 1e9], ...
              [0.35 0.55 0.85], 'FaceAlpha', 0.07, 'EdgeColor','none', ...
              'HandleVisibility','off');
    end
    for k = 2:numel(segs)
        xline(ax, segs(k).f_lo, '-', 'Color', [0.55 0.55 0.55], ...
              'LineWidth', 0.8, 'HandleVisibility','off');
    end
end

%% ========================================================================
%  3 + 4. TRAIN AND BUILD
% ========================================================================
function build_model(xlsFile, modelFile)
    cfg = default_config();
    [freq, names, X, widths] = load_measurements(xlsFile, cfg);
    segs = scan_segments(freq, names, X, cfg);

    nF   = numel(names);
    comp = resolve_composites(names, cfg);
    if ~isempty(comp)
        fprintf('\nComposite value(s):\n');
        for c = 1:numel(comp)
            fprintf('  %-12s = %s   (%d bits)\n', comp(c).name, ...
                strjoin(compose('%s<<%d', names(comp(c).idx), comp(c).shift), ' | '), ...
                comp(c).nbits);
        end
    end

    tAll = tic;

    % ---- fit, then bisect wherever a field failed, and refit -------------
    for pass = 0:cfg.maxRefine
        [Model, fail] = fit_all_segments(freq, names, X, widths, segs, cfg, comp, pass == 0);

        if isempty(fail) || ~cfg.refine || pass == cfg.maxRefine
            break;
        end

        [segs, nAdded] = bisect_failing(segs, fail, freq, cfg);
        if nAdded == 0
            fprintf(['\nRefinement stopped: the failing segment(s) cannot be ' ...
                     'bisected further\n(minimum %d points). The residual is not ' ...
                     'reducible by segmentation.\n'], cfg.refineMinPts);
            break;
        end
        fprintf('\n=== Refinement pass %d: %d failing pair(s), %d new boundary(ies) ===\n', ...
            pass+1, size(fail,1), nAdded);
    end

    save(modelFile, 'Model');
    fprintf('\nModel saved to %s  (%.2f s, %d segments)\n', ...
        modelFile, toc(tAll), numel(Model.segment));

    report_build(Model, freq, names, X, cfg);
end

function comp = resolve_composites(names, cfg)
% Turn cfg.compositeFields into index form, dropping any whose members are
% not all present in the measurement set.
    comp = struct('name',{}, 'idx',{}, 'shift',{}, 'width',{}, 'nbits',{});
    if ~isfield(cfg,'compositeFields') || isempty(cfg.compositeFields), return; end

    for c = 1:numel(cfg.compositeFields)
        s  = cfg.compositeFields{c};
        [tf, loc] = ismember(s.parts, names);
        if ~all(tf)
            warning('Composite %s skipped: %s not in the measurement set.', ...
                s.name, strjoin(s.parts(~tf), ', '));
            continue;
        end
        w = register_bit_width(s.parts(:));
        comp(end+1).name  = s.name;                    %#ok<AGROW>
        comp(end).idx     = loc(:).';
        comp(end).shift   = s.shift(:).';
        comp(end).width   = w(:).';
        comp(end).nbits   = max(s.shift(:).' + w(:).');
    end
end

function [Model, fail] = fit_all_segments(freq, names, X, widths, segs, cfg, comp, verbose)
    nF = numel(names);
    Model = struct('names',names, 'widths',widths, 'cfg',cfg, ...
                   'f_lo',freq(1), 'f_hi',freq(end));
    Model.segment = struct('f_lo',{}, 'f_hi',{}, 'field',{}, 'train',{});
    fail = zeros(0,2);          % [segment, field]

    for k = 1:numel(segs)
        fs = segs(k).freq;  Xs = X(segs(k).idx, :);  N = numel(fs);
        f2t = @(ff) (ff - fs(1))/(fs(end) - fs(1))*(N - 1);

        M      = min(cfg.deg, floor(N/3));
        kn     = (1:M).';
        x_cheb = 0.5*(fs(1)+fs(end)) + 0.5*(fs(end)-fs(1)).*cos((kn-0.5)*pi/M);
        t_cheb = f2t(x_cheb);

        if verbose
            fprintf('\n--- Segment %d  [%.0f - %.0f MHz]  N = %d, M = %d ---\n', ...
                k, fs(1), fs(end), N, M);
        end

        % ---- composite values, built before anything else -----------------
        Cval = cell(numel(comp), 1);
        for c = 1:numel(comp)
            v = zeros(N,1);
            for p = 1:numel(comp(c).idx)
                v = v + Xs(:, comp(c).idx(p)) .* 2^comp(c).shift(p);
            end
            Cval{c} = v;
        end
        memberOf = zeros(nF,1);  memberPos = zeros(nF,1);
        for c = 1:numel(comp)
            memberOf(comp(c).idx)  = c;
            memberPos(comp(c).idx) = 1:numel(comp(c).idx);
        end

        % ---- training -----------------------------------------------------
        kinds     = arrayfun(@(i) string(classify_field(Xs(:,i), names(i), cfg)), 1:nF);
        smoothIdx = find(kinds == "smooth" & memberOf.' == 0);
        [ah0, bt0] = table_scaling(cfg.J);
        tr = struct('beta',bt0, 'alpha',ah0, 'e_table',NaN, 'e_trained',NaN, ...
                    'scope',cfg.trainScope, 'perField',[]);

        trainOn = smoothIdx;
        if isempty(trainOn) && ~isempty(comp), trainOn = comp(1).idx(1); end

        if cfg.train && ~isempty(trainOn)
            switch lower(cfg.trainScope)
                case 'segment'
                    [~, h] = max(std(Xs(:,trainOn), 0, 1));
                    i0 = trainOn(h);
                    xr = detrend_or_not(Xs(:,i0), cfg);
                    [tr.beta, tr.alpha, tr.e_table, tr.e_trained] = ...
                        train_scaling(N, cfg, xr, t_cheb, bt0, ah0);
                    if verbose
                        fprintf('  trained on %s : %.3e -> %.3e counts (%.1fx)\n', ...
                            names(i0), tr.e_table, tr.e_trained, ...
                            tr.e_table/max(tr.e_trained,realmin));
                    end
                case 'field'
                    tr.perField = repmat(struct('beta',bt0,'alpha',ah0, ...
                        'e_table',NaN,'e_trained',NaN), nF, 1);
                    for i = trainOn
                        xr = detrend_or_not(Xs(:,i), cfg);
                        [b,a,e0,e1] = train_scaling(N, cfg, xr, t_cheb, bt0, ah0);
                        tr.perField(i) = struct('beta',b,'alpha',a, ...
                                                'e_table',e0,'e_trained',e1);
                    end
            end
        end

        % ---- fit the composites -------------------------------------------
        compFit = cell(numel(comp),1);
        for c = 1:numel(comp)
            [b, a] = scaling_for(tr, comp(c).idx(1), cfg);
            cf = fit_field(fs, Cval{c}, comp(c).name, cfg, b, a, x_cheb, t_cheb, N, true);
            compFit{c} = cf;
            if verbose
                fprintf('  composite %-10s error %.4f counts (on a %d-bit value)\n', ...
                    comp(c).name, cf.err, comp(c).nbits);
            end
        end

        % ---- fit every field ----------------------------------------------
        fld = empty_field_struct();
        for i = 1:nF
            if memberOf(i) > 0
                c = memberOf(i);  p = memberPos(i);
                f = compFit{c};
                f.kind    = 'composite';
                f.name    = names(i);
                f.cshift  = comp(c).shift(p);
                f.cwidth  = comp(c).width(p);
                f.err     = max(abs(eval_field(f, fs) - Xs(:,i)));
            else
                [b, a] = scaling_for(tr, i, cfg);
                f = fit_field(fs, Xs(:,i), names(i), cfg, b, a, x_cheb, t_cheb, N, false);
            end

            % ---- residual correction, only if asked for --------------------
            tol = cfg.tol;
            if f.err >= tol && cfg.residual && ismember(f.kind, {'smooth','composite'})
                f = add_residual(f, fs, Xs(:,i), cfg, names(i), verbose);
            end

            if f.err >= tol, fail(end+1,:) = [k i]; end   %#ok<AGROW>
            fld(i) = f;
        end

        Model.segment(k).f_lo  = segs(k).f_lo;
        Model.segment(k).f_hi  = segs(k).f_hi;
        Model.segment(k).field = fld;
        Model.segment(k).train = tr;

        if verbose
            e = [fld.err];
            fprintf('  worst field error in segment: %.4f counts (%s)\n', ...
                max(e), names(find(e == max(e), 1)));
        end
    end
end

function f = empty_field_struct()
    f = struct('kind',{}, 'name',{}, 'val',{}, 'edgeF',{}, 'edgeV',{}, ...
               'x_cheb',{}, 'Y_cheb',{}, 'err',{}, 'tag',{}, 'fref',{}, ...
               'vcoMin',{}, 'cshift',{}, 'cwidth',{}, 'rF',{}, 'rV',{});
end

function [segs, nAdded] = bisect_failing(segs, fail, freq, cfg)
% Split every segment that contains a failing field at its midpoint.
    bad = unique(fail(:,1));
    newEdges = [];
    for b = bad(:).'
        idx = segs(b).idx;
        if numel(idx) < 2*cfg.refineMinPts, continue; end
        newEdges(end+1) = idx(1) + floor(numel(idx)/2);   %#ok<AGROW>
    end
    nAdded = numel(newEdges);
    if nAdded == 0, return; end

    edges = unique([arrayfun(@(s) s.idx(1), segs), newEdges, numel(freq)+1]);
    out = struct('idx',{}, 'freq',{}, 'f_lo',{}, 'f_hi',{}, 'src',{});
    for k = 1:numel(edges)-1
        idx = edges(k):edges(k+1)-1;
        if isempty(idx), continue; end
        out(end+1).idx = idx;                                    %#ok<AGROW>
        out(end).freq  = freq(idx);
        out(end).f_lo  = freq(idx(1));
        out(end).f_hi  = freq(idx(end));
        out(end).src   = ternary(ismember(edges(k), newEdges), 'refine', 'kept');
    end
    segs = out;
end

function f = add_residual(f, fs, x, cfg, name, verbose)
% Store the integer difference between the rounded model and the measurement
% wherever it is non-zero.  Exact at the measured frequencies ONLY.
    r  = round(x(:)) - round(eval_field(f, fs));
    nz = find(r ~= 0);
    frac = numel(nz)/numel(fs);

    if frac > cfg.residualMaxFrac
        if verbose
            fprintf(['  %-14s residual correction REFUSED: would store %.0f%% of ' ...
                     'the samples\n                 (cap is %.0f%%). The data is not ' ...
                     'a function of frequency here.\n'], ...
                     name, 100*frac, 100*cfg.residualMaxFrac);
        end
        return;
    end

    f.rF  = fs(nz);
    f.rV  = r(nz);
    f.err = max(abs(eval_field(f, fs) - x(:)));
    if verbose
        fprintf('  %-14s residual correction: %d point(s) stored (%.2f%%)\n', ...
            name, numel(nz), 100*frac);
    end
end

function report_build(Model, freq, names, X, cfg)
    nF = numel(Model.names);
    Y  = predict_all(Model, freq);
    fprintf('\n%-14s %11s %10s %7s\n', 'Field', 'kind', 'Linf', 'pass');
    fprintf('%s\n', repmat('-', 1, 46));
    np = 0;
    for i = 1:nF
        j = find(names == Model.names(i), 1);
        e = max(abs(Y(:,i) - X(:,j)));
        ok = e < cfg.tol;  np = np + ok;
        kinds = unique(arrayfun(@(k) string(Model.segment(k).field(i).kind), ...
                                1:numel(Model.segment)));
        fprintf('%-14s %11s %10.4f %7s\n', Model.names(i), ...
            strjoin(kinds,'/'), e, string(ok));
    end
    fprintf('%s\n', repmat('-', 1, 46));
    fprintf('Fields within %.2f counts: %d / %d\n', cfg.tol, np, nF);
end

function [b, a] = scaling_for(tr, i, cfg)
    if strcmpi(cfg.trainScope,'field') && ~isempty(tr.perField)
        b = tr.perField(i).beta;  a = tr.perField(i).alpha;
    else
        b = tr.beta;  a = tr.alpha;
    end
end

function [beta, alpha_half, e_table, e_trained] = train_scaling(N, cfg, xr, t_cheb, bt0, ah0)
% Choose (beta, alpha_1..alpha_L) minimising the worst-case NUFFT error at the
% Chebyshev nodes for this particular signal.  Derivative-free: only L+1 = 3
% numbers, and each evaluation is one M x N matrix product.
    y_ref = exact_ndft(xr, t_cheb);

    loss = @(p) log10( max(abs( ...
        nufft_engine(N, cfg.J, cfg.padding, p(1), [1, p(2:end).'], xr, t_cheb) ...
        - y_ref)) + 1e-18 );

    p0 = [bt0; ah0(2:end).'];
    if numel(p0) - 1 < cfg.L
        p0 = [p0; zeros(cfg.L - (numel(p0)-1), 1)];
    end

    e_table = 10^loss(p0);

    opts = optimset('Display','off', 'TolX',1e-7, 'TolFun',1e-9, ...
                    'MaxFunEvals',cfg.maxFevals, 'MaxIter',cfg.maxFevals);
    p = fminsearch(loss, p0, opts);
    e_trained = 10^loss(p);

    if e_trained < e_table            % never accept a regression
        beta = p(1);   alpha_half = [1, p(2:end).'];
    else
        beta = p0(1);  alpha_half = [1, p0(2:end).'];
        e_trained = e_table;
    end
end

function f = fit_field(fs, x, name, cfg, beta, alpha_half, x_cheb, t_cheb, N, forceSmooth)
    if nargin < 10, forceSmooth = false; end
    x = x(:);  fs = fs(:);
    f = struct('kind','', 'name',name, 'val',[], 'edgeF',[], 'edgeV',[], ...
               'x_cheb',[], 'Y_cheb',[], 'err',0, ...
               'tag',"", 'fref',cfg.fref, 'vcoMin',cfg.vcoMin, ...
               'cshift',0, 'cwidth',0, 'rF',[], 'rV',[]);

    if forceSmooth
        f.kind = 'smooth';        % composite values are always modelled
    else
        f.kind = classify_field(x, name, cfg);
    end

    switch f.kind
        case 'formula'
            % Nothing stored and nothing fitted: the value is computed from the
            % frequency.  f.err doubles as a check on cfg.fref.
            f.tag = formula_tag(name);
            f.err = max(abs(eval_field(f, fs) - x));

        case 'const'
            f.val = x(1);

        case {'exact','lookup'}
            % store transitions only: exact at every measured frequency,
            % piecewise-constant between them
            j = [1; find(diff(x) ~= 0) + 1];
            f.edgeF = fs(j);  f.edgeV = x(j);
            f.err   = max(abs(eval_field(f, fs) - x));

        case 'smooth'
            [xr, c0, sl] = detrend_or_not(x, cfg);
            yr = nufft_engine(N, cfg.J, cfg.padding, beta, alpha_half, xr, t_cheb);
            f.x_cheb = x_cheb;
            f.Y_cheb = retrend(yr, c0, sl, t_cheb, cfg);
            f.err    = max(abs(barycentric_cheb(f.x_cheb, f.Y_cheb, fs) - x));
    end
end

function v = eval_field(f, fq)
% Must stay byte-for-byte identical to eval_field in Catalina_MakeProfiles.m:
% this one decides the error reported at build time, that one decides what
% goes into the profile.  If they diverge, a model can validate clean and
% still emit wrong bytes.
    fq = fq(:);
    switch f.kind
        case 'formula'
            R = pll_registers(fq, f.fref, struct('vcoMin', f.vcoMin, 'warn', false));
            v = R.(char(f.tag));
        case 'const'
            v = repmat(f.val, numel(fq), 1);
        case {'exact','lookup'}
            idx = discretize(fq, [f.edgeF; inf]);
            idx(fq < f.edgeF(1)) = 1;
            idx(isnan(idx))      = numel(f.edgeV);
            v = f.edgeV(idx);
        case 'smooth'
            v = barycentric_cheb(f.x_cheb, f.Y_cheb, fq);
        case 'composite'
            % Round the wide value to an integer FIRST, then slice out this
            % member's bits.  Rounding the byte and the MSB separately would
            % let them disagree - 255.6 rounds to 256 in an 8-bit slice while
            % a separately-rounded MSB still says 0.  Arithmetic rather than
            % bitand/bitshift so nothing is forced into an integer class,
            % where uint32(NaN) would silently become 0.
            w = max(round(barycentric_cheb(f.x_cheb, f.Y_cheb, fq)), 0);
            v = mod(floor(w / 2^f.cshift), 2^f.cwidth);
        otherwise
            error('Unknown field kind "%s".', f.kind);
    end
    v = v(:);

    if isfield(f,'rF') && ~isempty(f.rF)
        [tf, loc] = ismember(fq, f.rF);
        v(tf) = round(v(tf)) + f.rV(loc(tf));
    end
end

function [xr, c0, sl] = detrend_or_not(x, cfg)
% Remove the line through the first and last sample.  The FFT treats the record
% as periodic; without this a field whose value differs at the two ends carries
% a step discontinuity that does not exist in the hardware, and the resulting
% ringing dominates every other error term.
    x = x(:);
    if cfg.detrend
        c0 = x(1);
        sl = (x(end) - x(1))/(numel(x) - 1);
        xr = x - (c0 + sl*(0:numel(x)-1).');
    else
        xr = x;  c0 = 0;  sl = 0;
    end
end

function y = retrend(yr, c0, sl, t, cfg)
    if cfg.detrend, y = yr(:) + c0 + sl*t(:); else, y = yr(:); end
end

%% ========================================================================
%  5. VALIDATE
% ========================================================================
function validate_model(modelFile, xlsFile)
    S = load(modelFile,'Model');  Model = S.Model;
    [freq, names, X, ~] = load_measurements(xlsFile, Model.cfg);
    [tf, loc] = ismember(Model.names, names);

    Y  = predict_all(Model, freq);
    nF = numel(Model.names);

    fprintf('\n%-14s %8s %10s %10s %7s\n','Field','kind','Linf','RMS','pass');
    fprintf('%s\n', repmat('-',1,54));
    Linf = nan(nF,1);  pass = false(nF,1);
    for i = 1:nF
        if ~tf(i), continue; end
        e = Y(:,i) - X(:,loc(i));
        Linf(i) = max(abs(e));
        pass(i) = Linf(i) < Model.cfg.tol;
        kinds = unique(arrayfun(@(k) string(Model.segment(k).field(i).kind), ...
                                1:numel(Model.segment)));
        fprintf('%-14s %8s %10.4f %10.5f %7s\n', Model.names(i), ...
            strjoin(kinds,'/'), Linf(i), sqrt(mean(e.^2)), string(pass(i)));
    end
    fprintf('%s\n', repmat('-',1,54));
    fprintf('Fields with |e|_inf < %.2f counts : %d / %d\n', ...
        Model.cfg.tol, sum(pass), sum(tf));
    if any(~pass & tf)
        fprintf('Failing: %s\n', strjoin(Model.names(~pass & tf), ', '));
    end

    plot_summary(Model, Linf, pass, tf);
    for i = 1:nF
        if tf(i), plot_field_detail(Model, freq, X(:,loc(i)), Y(:,i), i, Linf(i)); end
    end
end

%% ---- validation plots ---------------------------------------------------
function plot_summary(Model, Linf, pass, tf)
% The one graph to put in the report: every field, worst-case error, against
% the criterion.  Log axis because the errors span many decades.
    idx = find(tf);
    [~, ord] = sort(Linf(idx), 'descend');
    idx = idx(ord);
    v   = max(Linf(idx), 1e-6);
    n   = numel(idx);

    figure('Name','Validation summary','Color','w','Position',[80 80 900 40+26*n]);
    ax = axes; hold(ax,'on');

    for r = 1:n
        c = [0.75 0.20 0.18];
        if pass(idx(r)), c = [0.16 0.52 0.28]; end
        barh(ax, r, v(r), 0.62, 'FaceColor', c, 'EdgeColor','none');
    end
    xline(ax, Model.cfg.tol, '--', 'Color',[0.1 0.1 0.1], 'LineWidth',1.6, ...
          'Label', sprintf('%.2g counts', Model.cfg.tol), ...
          'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left', ...
          'FontSize', 11);

    for r = 1:n
        text(ax, v(r)*1.15, r, sprintf('%.3g', Linf(idx(r))), ...
             'FontSize', 10, 'VerticalAlignment','middle');
    end

    set(ax, 'XScale','log', 'YTick', 1:n, 'YTickLabel', Model.names(idx), ...
            'YDir','reverse', 'TickLabelInterpreter','none');
    ylim(ax, [0.4 n+0.6]);
    xlim(ax, [min(v)/3, max(v)*6]);
    xlabel(ax, 'worst-case error  ||e||_\infty  [counts]');
    title(ax, sprintf('%d of %d fields within %.2g counts', ...
          sum(pass), sum(tf), Model.cfg.tol), 'FontSize', 13);
    style_axes(ax);
end

function plot_field_detail(Model, freq, xm, ym, i, Linf)
% One field per figure: model over data on top, error underneath, shared axis.
    tol  = Model.cfg.tol;
    ok   = Linf < tol;
    kinds = unique(arrayfun(@(k) string(Model.segment(k).field(i).kind), ...
                            1:numel(Model.segment)));

    figure('Name', char(Model.names(i)), 'Color','w', 'Position',[100 100 1150 620]);
    tl = tiledlayout(3, 1, 'TileSpacing','compact', 'Padding','compact');

    % ---- top: measurement vs model (2 tiles tall) ---------------------
    ax1 = nexttile([2 1]); hold(ax1,'on');
    shade_segments(ax1, Model.segment);
    hM = plot(ax1, freq, xm, '.', 'Color',[0.55 0.55 0.55], 'MarkerSize', 7);
    hP = plot(ax1, freq, ym, '-', 'Color',[0.10 0.35 0.75], 'LineWidth', 1.4);
    ylabel(ax1, 'counts');
    legend(ax1, [hM hP], {'measured','model'}, 'Location','best', 'Box','off');
    ax1.XTickLabel = [];
    ylim(ax1, padded_limits(xm));
    style_axes(ax1);

    % ---- bottom: error ------------------------------------------------
    ax2 = nexttile; hold(ax2,'on');
    e = max(abs(ym - xm), 1e-6);
    patch(ax2, [freq(1) freq(end) freq(end) freq(1)], [1e-6 1e-6 tol tol], ...
          [0.16 0.52 0.28], 'FaceAlpha',0.08, 'EdgeColor','none');
    area(ax2, freq, e, 'BaseValue', 1e-6, 'FaceColor',[0.75 0.20 0.18], ...
         'FaceAlpha', 0.35, 'EdgeColor',[0.6 0.15 0.13], 'LineWidth', 0.7);
    yline(ax2, tol, '--', 'Color',[0.1 0.1 0.1], 'LineWidth', 1.5);
    set(ax2, 'YScale','log');
    ylim(ax2, [1e-6 max(max(e)*3, tol*3)]);
    ylabel(ax2, '|error| [counts]');  xlabel(ax2, 'frequency [MHz]');
    style_axes(ax2);

    linkaxes([ax1 ax2], 'x');
    xlim(ax1, [freq(1) freq(end)]);

    if ok, verdict = 'PASS'; col = [0.16 0.52 0.28];
    else,  verdict = 'FAIL'; col = [0.75 0.20 0.18];
    end
    t = title(tl, sprintf('%s   |   %s   |   ||e||_\\infty = %.4g counts   |   %s', ...
        Model.names(i), strjoin(kinds,' / '), Linf, verdict), ...
        'FontSize', 13, 'FontWeight','bold');
    t.Color = col;
end

function L = padded_limits(v)
    lo = min(v);  hi = max(v);
    if hi == lo, hi = lo + 1; end
    L = [lo - 0.06*(hi-lo), hi + 0.06*(hi-lo)];
end

function Y = predict_all(Model, fq)
    fq = fq(:);
    nF = numel(Model.names);
    Y  = nan(numel(fq), nF);
    hi = [Model.segment.f_hi];
    K  = numel(Model.segment);
    for k = 1:K
        if     k == 1, m = fq <= hi(1);
        elseif k == K, m = fq >  hi(k-1);
        else,          m = fq >  hi(k-1) & fq <= hi(k);
        end
        if ~any(m), continue; end
        for i = 1:nF
            Y(m,i) = eval_field(Model.segment(k).field(i), fq(m));
        end
    end
end

%% ========================================================================
%  6. LUT GENERATION
% ========================================================================
function generate_profiles(modelFile, targetFile, outFile)
    S = load(modelFile,'Model');  Model = S.Model;

    tf = readmatrix(targetFile);
    tf = tf(:,1);  tf = tf(~isnan(tf));
    out = tf < Model.f_lo | tf > Model.f_hi;
    if any(out)
        warning('%d target(s) outside %g-%g MHz - extrapolating.', ...
            sum(out), Model.f_lo, Model.f_hi);
    end

    Y = round(predict_all(Model, tf));
    for i = 1:numel(Model.names)
        Y(:,i) = min(max(Y(:,i), 0), 2^Model.widths(i) - 1);
    end

    T = array2table([tf, Y], 'VariableNames', ...
        matlab.lang.makeValidName(["Target_Frequency_MHz", Model.names']));
    writetable(T, outFile);
    fprintf('%d Fast-Lock profiles written to %s\n', numel(tf), outFile);
end

%% ========================================================================
%  7. ANALYSIS
% ========================================================================
function error_analysis(xlsFile)
    cfg = default_config();
    [freq, names, X, ~] = load_measurements(xlsFile, cfg);
    segs = scan_segments(freq, names, X, cfg, false);

    best = struct('k',0,'i',0,'r',-1);
    for k = 1:numel(segs)
        for i = 1:numel(names)
            if ~strcmp(classify_field(X(segs(k).idx,i), names(i), cfg),'smooth'), continue; end
            r = std(X(segs(k).idx,i));
            if r > best.r, best = struct('k',k,'i',i,'r',r); end
        end
    end
    if best.k == 0, error('No smooth field found to analyse.'); end

    k = best.k;  i = best.i;
    fs = segs(k).freq;  x = X(segs(k).idx,i);  N = numel(fs);
    f2t = @(ff) (ff - fs(1))/(fs(end)-fs(1))*(N-1);
    [xr, c0, sl] = detrend_or_not(x, cfg);
    fprintf('\nAnalysing %s, segment %d (%.0f-%.0f MHz, N = %d)\n', ...
        names(i), k, fs(1), fs(end), N);

    M = min(cfg.deg, floor(N/3));  kn = (1:M).';
    x_cheb = 0.5*(fs(1)+fs(end)) + 0.5*(fs(end)-fs(1)).*cos((kn-0.5)*pi/M);
    t_cheb = f2t(x_cheb);

    % ---- training gain vs J -------------------------------------------
    J_vec = 2:2:16;
    eTab = zeros(size(J_vec));  eTrn = zeros(size(J_vec));
    fprintf('\n%4s %13s %13s %9s\n','J','table alpha','trained','gain');
    for a = 1:numel(J_vec)
        c2 = cfg;  c2.J = J_vec(a);
        [ah0, bt0] = table_scaling(J_vec(a));
        [~,~,e0,e1] = train_scaling(N, c2, xr, t_cheb, bt0, ah0);
        eTab(a) = e0;  eTrn(a) = e1;
        fprintf('%4d %13.3e %13.3e %8.1fx\n', J_vec(a), e0, e1, e0/max(e1,realmin));
    end

    % ---- total error vs M ---------------------------------------------
    M_vec = [16 24 32 48 64 96 120 160 200];
    M_vec = M_vec(M_vec < N/3);
    [ah0, bt0] = table_scaling(cfg.J);
    eM = zeros(size(M_vec));
    fprintf('\n%5s %14s\n','M','e_total');
    for b = 1:numel(M_vec)
        Mb = M_vec(b);  kb = (1:Mb).';
        xc = 0.5*(fs(1)+fs(end)) + 0.5*(fs(end)-fs(1)).*cos((kb-0.5)*pi/Mb);
        tc = f2t(xc);
        bt = bt0;  ah = ah0;
        if cfg.train, [bt, ah] = train_scaling(N, cfg, xr, tc, bt0, ah0); end
        yc = retrend(nufft_engine(N, cfg.J, cfg.padding, bt, ah, xr, tc), c0, sl, tc, cfg);
        eM(b) = max(abs(barycentric_cheb(xc, yc, fs) - x));
        fprintf('%5d %14.3e\n', Mb, eM(b));
    end

    figure('Name','Error analysis','Color','w','Position',[60 60 1250 500]);
    tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');

    % ---- left: training gain -------------------------------------------
    ax = nexttile; hold(ax,'on');
    fill(ax, [J_vec fliplr(J_vec)], [max(eTab,eps) fliplr(max(eTrn,eps))], ...
         [0.30 0.55 0.85], 'FaceAlpha',0.13, 'EdgeColor','none');
    h1 = plot(ax, J_vec, max(eTab,eps), '-s', 'LineWidth',1.6, ...
              'Color',[0.45 0.45 0.45], 'MarkerFaceColor','w', 'MarkerSize',7);
    h2 = plot(ax, J_vec, max(eTrn,eps), '-o', 'LineWidth',2.2, ...
              'Color',[0.10 0.35 0.75], 'MarkerFaceColor',[0.10 0.35 0.75], ...
              'MarkerSize',6);
    h3 = yline(ax, cfg.tol, '--', 'Color',[0.1 0.1 0.1], 'LineWidth',1.5);
    set(ax,'YScale','log');
    for a = 1:2:numel(J_vec)
        text(ax, J_vec(a), sqrt(max(eTab(a),eps)*max(eTrn(a),eps)), ...
             sprintf(' %.0fx', eTab(a)/max(eTrn(a),realmin)), ...
             'FontSize',10, 'Color',[0.10 0.35 0.75], 'FontWeight','bold');
    end
    xlabel(ax,'J   (NUFFT neighbours)');
    ylabel(ax,'NUFFT error at the Chebyshev nodes [counts]');
    legend(ax, [h1 h2 h3], {'Table 2 \alpha (published)', ...
           'trained \alpha (this work)', sprintf('%.2g counts', cfg.tol)}, ...
           'Location','southwest', 'Box','off');
    title(ax, 'Training the scaling vector', 'FontSize',13);
    xticks(ax, J_vec);  style_axes(ax);

    % ---- right: total error vs M ---------------------------------------
    ax = nexttile; hold(ax,'on');
    plot(ax, M_vec, max(eM,eps), '-d', 'LineWidth',2.2, ...
         'Color',[0.75 0.20 0.18], 'MarkerFaceColor',[0.75 0.20 0.18], 'MarkerSize',7);
    yline(ax, cfg.tol, '--', 'Color',[0.1 0.1 0.1], 'LineWidth',1.5);
    set(ax,'YScale','log');
    good = find(eM < cfg.tol, 1);
    if ~isempty(good)
        xline(ax, M_vec(good), ':', 'Color',[0.16 0.52 0.28], 'LineWidth',1.8, ...
              'Label', sprintf('M = %d suffices', M_vec(good)), 'FontSize',11, ...
              'LabelOrientation','horizontal');
    end
    xlabel(ax,'M   (Chebyshev nodes per segment)');
    ylabel(ax,'total error vs measurement [counts]');
    title(ax, sprintf('Polynomial degree  (J = %d, trained \\alpha)', cfg.J), ...
          'FontSize',13);
    xticks(ax, M_vec);  style_axes(ax);

    sgtitle(sprintf('%s, segment %d  (%.0f - %.0f MHz, N = %d)', ...
        names(i), k, fs(1), fs(end), N), 'FontSize',14, 'FontWeight','bold');

    save('error_analysis_data.mat','J_vec','eTab','eTrn','M_vec','eM');
    fprintf('\nSaved error_analysis_data.mat\n');
end

%% ========================================================================
%  NUFFT CORE
% ========================================================================
function [alpha_half, beta] = table_scaling(J)
% Fessler & Sutton starting point, K/N = 2, L = 2.  Used as the initial guess
% for training and as the baseline the trained result must beat.
    switch J
        case 2,  beta = 0.8079; alpha_half = [1 -0.1233 0.0080];
        case 4,  beta = 0.1842; alpha_half = [1 -0.6514 0.1544];
        case 6,  beta = 0.1505; alpha_half = [1 -0.6632 0.1647];
        case 8,  beta = 0.1318; alpha_half = [1 -0.6602 0.1606];
        case 10, beta = 0.1025; alpha_half = [1 -0.6687 0.1690];
        case 12, beta = 0.1607; alpha_half = [1 -0.6625 0.1635];
        case 14, beta = 0.2349; alpha_half = [1 -0.6594 0.1642];
        case 16, beta = 0.2152; alpha_half = [1 -0.6644 0.1682];
        otherwise, beta = 0.1025; alpha_half = [1 -0.6687 0.1690];
    end
end

function s = scaling_vector(N, K, beta, alpha_half)
    alpha_half = alpha_half(:);
    L  = numel(alpha_half) - 1;
    n0 = (N-1)/2;  n = (0:N-1).';
    a  = [flipud(conj(alpha_half(2:end))); alpha_half(1); alpha_half(2:end)];
    s  = exp(1i*(2*pi/K)*beta*(n-n0)*(-L:L)) * a;
end

function y = nufft_engine(N, J, padding, beta, alpha_half, x, t, blockSize)
    x = x(:);  t = t(:);
    K = N*padding;  gama = 2*pi/K;  n0 = (N-1)/2;  n = (0:N-1).';  jrow = 1:J;

    X = fftshift(fft(x));
    s = scaling_vector(N, K, beta, alpha_half);
    Y = fft(s .* X, K);                                       % eq. 3

    C = exp(1i*gama*(n-n0)*jrow)/sqrt(N);                     % eq. 16
    [Q, R] = qr(conj(s).*C, 0);                               % eq. 23
    invRQ = R \ Q';                                           % eq. 24

    omega = -2*pi*t/N;  q = omega/gama;
    if mod(J,2) == 1, k0 = round(q) - (J+1)/2;                % eq. 7
    else,             k0 = floor(q) - J/2;
    end

    if nargin < 8 || isempty(blockSize)
        blockSize = max(1, floor(2e6/max(N,1)));
    end

    M = numel(omega);  y = zeros(M,1);
    for a = 1:blockSize:M
        b   = min(a+blockSize-1, M);
        om  = omega(a:b);  kk = k0(a:b);
        Bm  = exp(1i*(om - gama*kk)*(n-n0).')/sqrt(N);        % eq. 18
        Lam = exp(-1i*(om - gama*(kk+jrow))*n0);              % eq. 17
        U   = conj(Lam) .* (Bm * invRQ.');   % plain transpose .', not '
        idx = mod(kk + jrow, K) + 1;
        y(a:b) = sum(Y(idx) .* conj(U), 2);                   % eq. 9
    end

    % Final phase.  fftshift places DC at index floor(N/2), so the shift is
    % floor(N/2), NOT N/2.  These are equal only for even N.  With odd N the
    % old exp(-1i*pi*t) leaves a residual phase and the error jumps from
    % ~1e-5 to ~5 counts, flat in J.  Verified:
    %     N = 510 even -> 1.2e-05     N = 511 odd -> 5.9e+00
    %     N = 5300 even -> 1.6e-05    N = 5301 odd -> 5.6e+00
    % Every segment here has an odd point count (1500-3000 MHz at 1 MHz is
    % 1501 points; 700-6000 is 5301), so this line is not optional.
    y = (1/N) * exp(-1i*2*pi*t*floor(N/2)/N) .* y;
    if isreal(x), y = real(y); end
end

function y = exact_ndft(x, t, blockSize)
% Exact O(N*M) evaluation of the same trigonometric interpolant the NUFFT
% approximates.  This is the training reference.
    x = x(:);  t = t(:);  N = numel(x);
    X = fftshift(fft(x));
    k = (-floor(N/2)) : (ceil(N/2) - 1);   % NOT floor(-N/2): differs for odd N
    if nargin < 3 || isempty(blockSize), blockSize = max(1, floor(2e6/N)); end
    y = zeros(numel(t),1);
    for a = 1:blockSize:numel(t)
        b = min(a+blockSize-1, numel(t));
        y(a:b) = (1/N) * (exp(1i*2*pi*t(a:b)*k/N) * X);
    end
    if isreal(x), y = real(y); end
end

function y = barycentric_cheb(xn, yn, xe)
    xn = xn(:);  yn = yn(:);  xe = xe(:);
    d = numel(xn);
    w = ((-1).^(0:d-1).') .* sin(((1:d).' - 0.5)*pi/d);
    D = xe.' - xn;
    hit = abs(D) < 1e-12;  D(hit) = 1;
    y = (((w.*yn).' * (1./D)) ./ (w.' * (1./D))).';
    [ri, ci] = find(hit);  y(ci) = yn(ri);
end
