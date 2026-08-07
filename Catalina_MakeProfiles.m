%% ========================================================================
%  AD9361 Tx Fast Lock profile generator
%
%  Model + target frequency list  ->  16 setup words per frequency, in hex.
%
%  Usage
%  -----
%     Catalina_MakeProfiles
%     Catalina_MakeProfiles('model.mat','targets.xlsx','profiles.xlsx')
%
%  Chain
%  -----
%     predict 21 register fields at each target frequency
%       -> ROUND to integer
%       -> clip to each field's bit width
%       -> reassemble the source register bytes
%       -> pack into the 16 setup words per UG-570 Table 12
%       -> UNPACK again and check the round trip
%       -> write hex
%
%  Output columns W0..WF are written as text with a 0x prefix.  Without the
%  prefix Excel reads "00", "12", "99" as numbers and "0A", "B5" as text, so
%  a single column comes back mixed and hex2dec on it fails.  Set
%  cfg.hexPrefix = '' if your loader wants bare hex, and read the file with
%  readtable(..., 'TextType','string') so the values stay text.
%
%  Register addresses (verified against UG-570 and the ADI no-OS driver):
%      Rx  0x25A SETUP  0x25B INIT_DELAY  0x25C ADDR  0x25D DATA
%          0x25E READ   0x25F CTRL
%      Tx  the same block offset by 0x40:  0x29A .. 0x29F
%  Table 12's header says the Tx setup words go to 0x25D; that is the Rx data
%  register and appears carried over from Table 11.  Table 12's own source
%  registers (0x271-0x291) are Tx, and UG-570's body text places the Tx fast
%  lock block at 0x29C-0x29F.
% ========================================================================

function Catalina_MakeProfiles(modelFile, freqFile, outFile)
    if nargin >= 1 && (isstring(modelFile) || ischar(modelFile)) && ...
       strcmpi(string(modelFile), "selftest")
        selftest();
        return;
    end
    if nargin < 1 || isempty(modelFile), modelFile = 'model.mat';     end
    if nargin < 2 || isempty(freqFile),  freqFile  = 'targets.xlsx';  end
    if nargin < 3 || isempty(outFile),   outFile   = 'profiles.xlsx'; end

    cfg = profile_config();

    % ---- 1. model and targets -------------------------------------------
    if ~isfile(modelFile), error('Model file "%s" not found.', modelFile); end
    S = load(modelFile, 'Model');  Model = S.Model;
    check_model_compatible(Model, modelFile);

    fq = read_target_frequencies(freqFile);
    fprintf('Model  : %s   (%d fields, %.4g - %.4g MHz)\n', modelFile, ...
        numel(Model.names), Model.f_lo, Model.f_hi);
    fprintf('Targets: %s   (%d frequencies, %.4g - %.4g MHz)\n', freqFile, ...
        numel(fq), min(fq), max(fq));

    out = fq < Model.f_lo | fq > Model.f_hi;
    if any(out)
        warning(['%d target(s) lie outside the modelled band and will be ' ...
                 'EXTRAPOLATED. Their profiles should not be trusted.'], sum(out));
    end

    if ~isempty(cfg.frefOverride)
        Model = override_fref(Model, cfg.frefOverride);
    end
    report_formula_fields(Model);

    % ---- 2. predict, round, clip ----------------------------------------
    Yraw = predict_all(Model, fq);

    bad = any(isnan(Yraw), 2);
    if any(bad)
        error(['%d target frequency(ies) produced NaN, e.g. %.6g MHz. ' ...
               'No segment covers them.'], sum(bad), fq(find(bad,1)));
    end

    Y = round(Yraw);        % round-half-away-from-zero; values are non-negative
                            % integers, so this is the plain rounding rule the
                            % 0.5-count criterion is built around

    % ---- out-of-range handling ------------------------------------------
    % A register field physically cannot hold a value outside 0 .. 2^w-1, so
    % something has to happen.  But clipping silently is how a modelling
    % failure gets written into a profile and never noticed, so the overshoot
    % is measured and reported, and cfg.clipMode decides whether it is fatal.
    over = struct('name',{}, 'n',{}, 'worst',{}, 'at',{});
    for i = 1:numel(Model.names)
        hi  = 2^Model.widths(i) - 1;
        d   = max(max(0 - Y(:,i), Y(:,i) - hi), 0);   % distance outside range
        idx = find(d > 0);
        if ~isempty(idx)
            [wv, wj] = max(d);
            over(end+1) = struct('name', Model.names(i), 'n', numel(idx), ...
                                 'worst', wv, 'at', fq(wj));       %#ok<AGROW>
        end
    end

    if ~isempty(over)
        fprintf('\nValues outside their field range:\n');
        fprintf('  %-14s %8s %12s %14s\n', 'Field','count','worst','at [MHz]');
        for k = 1:numel(over)
            fprintf('  %-14s %8d %12.3f %14.6g\n', over(k).name, over(k).n, ...
                over(k).worst, over(k).at);
        end
        worstOver = max([over.worst]);
        switch lower(cfg.clipMode)
            case 'error'
                error(['%d field(s) predict values their register cannot hold ' ...
                       '(worst overshoot %.3f counts).\nThat is a modelling ' ...
                       'failure, not a formatting one - clipping would write a ' ...
                       'wrong\nvalue into the profile and hide it. Fix the model, ' ...
                       'or set cfg.clipMode\nto ''saturate'' if you have decided ' ...
                       'the overshoot is acceptable.'], numel(over), worstOver);
            case 'saturate'
                fprintf(['Saturating to the field range (cfg.clipMode = ' ...
                         '''saturate'').\n']);
                if worstOver >= cfg.tol
                    warning(['Worst overshoot %.3f counts is at or above the ' ...
                             '%.2f-count criterion,\nso at least one profile is ' ...
                             'wrong by more than a rounding step.'], ...
                             worstOver, cfg.tol);
                end
                for i = 1:numel(Model.names)
                    Y(:,i) = min(max(Y(:,i), 0), 2^Model.widths(i) - 1);
                end
            otherwise
                error('cfg.clipMode must be ''error'' or ''saturate''.');
        end
    end

    margin = max(abs(Yraw - Y), [], 1);
    worst  = max(margin);
    fprintf('Largest rounding distance: %.3f counts (%s)%s\n', worst, ...
        Model.names(find(margin == worst, 1)), ...
        ternary(worst < 0.5, '   [within the 0.5 criterion]', ...
                             '   [OUTSIDE 0.5 - that profile may be wrong]'));

    % ---- 3. parse the field names ONCE ----------------------------------
    P = parse_fields(Model.names);
    if cfg.showParseTable, print_parse_table(P); end
    W = fastlock_map();
    check_coverage(P, W);
    W = bind_map(W, P);          % resolve every part to (address, hi, lo) now

    % ---- 4. pack --------------------------------------------------------
    nF = numel(fq);
    Bytes = zeros(nF, 16);
    for r = 1:nF
        reg        = assemble_registers(Y(r,:), P);
        Bytes(r,:) = pack_setup_words(reg, W, cfg);
    end

    if any(Bytes(:) < 0 | Bytes(:) > 255 | Bytes(:) ~= fix(Bytes(:)))
        error('Packing produced a value outside 0..255. Check fastlock_map().');
    end

    % ---- 4b. diagnostics ------------------------------------------------
    % An all-zero result can come from either stage, and the two have entirely
    % different causes, so say which one produced it.
    diagnose_stage('predicted field values', Y,     Model.names);
    diagnose_stage('packed setup words',     Bytes, "W" + string(dec2hex(0:15,1)).');

    if all(Y(:) == 0)
        error(['Every predicted field value is zero. The packing stage is not ' ...
               'at fault -\nthe model returned zeros. Load model.mat and check ' ...
               'Model.segment(1).field(1),\nor rebuild with Catalina_Project ' ...
               'option 2.']);
    end
    if all(Bytes(:) == 0)
        error(['Field values are non-zero but every packed byte is zero. That ' ...
               'is a packing\nfault, not a model fault. Run ' ...
               'Catalina_MakeProfiles(''selftest'') to check the\npacker against ' ...
               'known values.']);
    end

    % ---- 5. round trip: unpack the words and compare --------------------
    verify_roundtrip(Bytes, Y, Model.names, P, W, cfg);

    % ---- 6. write -------------------------------------------------------
    write_profile_table(outFile, fq, Bytes, cfg);
    if cfg.writeFieldDump
        write_field_dump(replace_ext(outFile, '_fields.xlsx'), fq, Y, Model.names);
    end
    print_preview(fq, Bytes, W, cfg);
end

%% ========================================================================
%  CONFIGURATION
% ========================================================================
function cfg = profile_config()
    cfg.frefOverride = [];      % MHz; [] = use the value stored in the model

    % ---- Init N values (words 6, 8, 9, B) -------------------------------
    cfg.initMode   = 'copy';    % 'copy' | 'zero' | 'manual'
    cfg.initValues = struct('CP', 0, 'R3', 0, 'C3', 0, 'R1', 0);

    % ---- out-of-range handling ------------------------------------------
    cfg.clipMode = 'error';     % 'error'    stop if a value cannot fit its field
                                % 'saturate' clamp to the range and report how
                                %            far out it was
    cfg.tol      = 0.5;         % counts; the rounding criterion

    % ---- output ---------------------------------------------------------
    cfg.hexPrefix      = '0x';  % '' for bare hex (see the note in the header)
    cfg.numProfiles    = 8;     % Program Address[7:4] holds 0..7
    cfg.writeFieldDump = true;
    cfg.showParseTable = true;  % print how each field name was interpreted
    cfg.previewRows    = 5;
    cfg.verifyWrite    = true;  % read the spreadsheet back and compare
    cfg.verifyStrict   = true;  % error rather than warn on a verification failure
end

%% ========================================================================
%  UG-570 TABLE 12
% ========================================================================
function W = fastlock_map()
% Each setup word is one or more parts.  A part is either a measured field,
% named exactly as in the measurement file, or an Init N value tagged with a
% leading '#'.  dst is the [hi lo] destination bit range inside the word.
%
% The datasheet's "shift left/right by n" are relative to where the field sits
% in its SOURCE register, which is why word A shifts one part up and the other
% down: C1 is at 0x27E[3:0] and has to reach [7:4], C2 is at 0x27E[7:4] and has
% to reach [3:0].  Destination ranges are written out here so there is nothing
% to infer.
    W = struct('idx',{},'desc',{},'src',{},'dst',{});

    W = add(W, 0,  'Synth Integer Word[7:0]',        {"271[7:0]",[7 0]});
    W = add(W, 1,  'Synth Integer Word[10:8]',       {"272[2:0]",[2 0]});
    W = add(W, 2,  'Synth Fractional Word[7:0]',     {"273[7:0]",[7 0]});
    W = add(W, 3,  'Synth Fractional Word[15:8]',    {"274[7:0]",[7 0]});
    W = add(W, 4,  'Synth Fractional Word[22:16]',   {"275[6:0]",[6 0]});
    W = add(W, 5,  'VCO Bias Ref | VCO Varactor',    {"282[2:0]",[6 4]}, ...
                                                     {"279[3:0]",[3 0]});
    W = add(W, 6,  'VCO Bias Tcf | CP Current Init', {"282[4:3]",[7 6]}, ...
                                                     {"#CP",     [5 0]});
    W = add(W, 7,  'Charge Pump Current',            {"27B[5:0]",[5 0]});
    W = add(W, 8,  'LF R3 | LF R3 Init',             {"280[3:0]",[7 4]}, ...
                                                     {"#R3",     [3 0]});
    W = add(W, 9,  'LF C3 | LF C3 Init',             {"27F[3:0]",[7 4]}, ...
                                                     {"#C3",     [3 0]});
    W = add(W, 10, 'LF C1 | LF C2',                  {"27E[3:0]",[7 4]}, ...
                                                     {"27E[7:4]",[3 0]});
    W = add(W, 11, 'LF R1 | LF R1 Init',             {"27F[7:4]",[7 4]}, ...
                                                     {"#R1",     [3 0]});
    W = add(W, 12, 'VCO Varactor Ref Tcf | Rx Div',  {"290[6:4]",[6 4]}, ...
                                                     {"005[7:4]",[3 0]});
    % Word D.  UG-570 says "VCO Cal Offset[3:0] shift left by 1", which would
    % place it at bits [4:1] and overlap VCO Varactor Reference[3:0] at [3:0].
    % That cannot be right, so [7:4] (shift left by 4) is used here - the only
    % non-overlapping reading, and consistent with words 8, 9 and B which all
    % pair a high nibble with a low nibble.  If your part disagrees, this is
    % the line to change.
    W = add(W, 13, 'VCO Cal Offset | Varactor Ref',  {"278[6:3]",[7 4]}, ...
                                                     {"291[3:0]",[3 0]});
    W = add(W, 14, 'Force VCO Tune[7:0]',            {"277[7:0]",[7 0]});
    W = add(W, 15, 'Force ALC | Force VCO Tune[8]',  {"276[6:0]",[7 1]}, ...
                                                     {"278[0]",  [0 0]});

    % ---- destination bits must not overlap inside a word -----------------
    for k = 1:numel(W)
        used = 0;
        for j = 1:numel(W(k).src)
            n    = W(k).dst{j}(1) - W(k).dst{j}(2) + 1;
            span = bitshift(2^n - 1, W(k).dst{j}(2));
            if bitand(used, span)
                error('fastlock_map: word %X has overlapping destination bits.', k-1);
            end
            used = bitor(used, span);
        end
    end
end

function W = add(W, idx, desc, varargin)
    n = numel(W) + 1;
    W(n).idx  = idx;
    W(n).desc = desc;
    W(n).src  = cellfun(@(c) c{1}, varargin, 'UniformOutput', false);
    W(n).dst  = cellfun(@(c) c{2}, varargin, 'UniformOutput', false);
end

function W = bind_map(W, P)
% Resolve every source string to (address, hi, lo) once, so the packing loop
% contains no string handling at all.  Also removes ~32 regex calls per
% target frequency.
    for k = 1:numel(W)
        n = numel(W(k).src);
        W(k).addr = strings(1,n);
        W(k).hi   = zeros(1,n);
        W(k).lo   = zeros(1,n);
        W(k).init = strings(1,n);
        for j = 1:n
            s = string(W(k).src{j});
            if startsWith(s, "#")
                W(k).init(j) = extractAfter(s, 1);
            else
                [a, hi, lo] = split_field_name(s);
                W(k).addr(j) = a;  W(k).hi(j) = hi;  W(k).lo(j) = lo;
                if ~any(P.addr == a)
                    error('Word %X needs register 0x%s, which the model lacks.', ...
                        W(k).idx, a);
                end
            end
        end
    end
end

%% ========================================================================
%  FIELD -> REGISTER -> SETUP WORD
% ========================================================================
function [a, hi, lo] = split_field_name(name)
%SPLIT_FIELD_NAME  "282[4:3]" -> "282", 4, 3.   "278[0]" -> "278", 0, 0.
%
% The bracket contents are captured whole and split afterwards, rather than
% matched with an optional capture group.  An optional trailing group is the
% subtle part: when it does not participate, MATLAB may return a shorter token
% list, and a `numel(t) >= 3` guard then silently falls through to lo = hi.
% Every multi-bit field would come out one bit wide, mod() would discard the
% rest, and the packed words would collapse toward zero - producing a
% spreadsheet of 0x00 with no error anywhere.
%
% Also tolerates the datasheet's own notation, 0x271[D7:D0], and a reversed
% range such as [0:7].
    name = strtrim(string(name));

    t = regexp(name, '^\s*(?:0[xX])?([0-9A-Fa-f]+)\s*\[([^\]]*)\]\s*$', ...
               'tokens', 'once');
    if isempty(t) || numel(t) < 2
        error(['Cannot parse the field name "%s".\n' ...
               'Expected forms: 278[6:3], 278[0], 0x278[D6:D3].'], name);
    end

    a = upper(string(t{1}));

    bits = regexprep(string(t{2}), '[Dd]', '');   % "D7:D0" -> "7:0"
    bits = strtrim(bits);
    parts = strtrim(split(bits, ':'));

    n = str2double(parts);
    if any(isnan(n)) || isempty(n) || numel(n) > 2
        error(['Cannot read the bit range from "%s" (bracket contents "%s").'], ...
              name, t{2});
    end

    if numel(n) == 2
        hi = max(n);  lo = min(n);
    else
        hi = n(1);    lo = n(1);
    end

    if hi < 0 || lo < 0 || hi > 31
        error('Field "%s" has an implausible bit range [%d:%d].', name, hi, lo);
    end
end

function P = parse_fields(names)
    n = numel(names);
    P.name = string(names(:));
    P.addr = strings(n,1);  P.hi = zeros(n,1);  P.lo = zeros(n,1);
    for i = 1:n
        [P.addr(i), P.hi(i), P.lo(i)] = split_field_name(P.name(i));
    end
    P.uaddr = unique(P.addr, 'stable');
    P.width = P.hi - P.lo + 1;

    % A field name containing ':' must parse to more than one bit.  If it does
    % not, the parser has failed and every value in that field is about to be
    % truncated to a single bit.
    suspect = contains(P.name, ':') & P.width == 1;
    if any(suspect)
        error(['%d field name(s) contain '':'' but parsed to a single bit, ' ...
               'e.g. "%s" -> [%d:%d].\nThe bit range was not read correctly ' ...
               'and every value would be truncated.'], ...
               sum(suspect), P.name(find(suspect,1)), ...
               P.hi(find(suspect,1)), P.lo(find(suspect,1)));
    end
end

function print_parse_table(P)
    fprintf('\nField name parsing:\n');
    fprintf('  %-14s %8s %6s %6s %7s\n', 'Name','Register','hi','lo','bits');
    for i = 1:numel(P.name)
        fprintf('  %-14s %8s %6d %6d %7d\n', P.name(i), "0x"+P.addr(i), ...
            P.hi(i), P.lo(i), P.width(i));
    end
    fprintf('  %d field(s) across %d register(s), %d total bit(s)\n', ...
        numel(P.name), numel(P.uaddr), sum(P.width));
end

function reg = assemble_registers(row, P)
% Rebuild each source register byte by OR-ing its fields back into their
% original bit positions.  Returns a containers.Map of char address -> double.
%
% Everything is kept in DOUBLE and combined with plus rather than bitor:
% the fields of one register occupy disjoint bit ranges by construction, so
% addition and bitwise-or agree, and no integer-class value ever reaches the
% map.  (A containers.Map declared ValueType 'double' rejects a uint32, which
% is what the previous version tried to store.)
    reg = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(P.name)
        w = P.hi(i) - P.lo(i) + 1;
        v = mod(floor(row(i)), 2^w) * 2^P.lo(i);
        a = char(P.addr(i));
        if isKey(reg, a), reg(a) = reg(a) + v;
        else,             reg(a) = v;
        end
    end
end

function v = field_value_bound(reg, addr, hi, lo)
    a = char(addr);
    if ~isKey(reg, a)
        error('Register 0x%s is required by Table 12 but is not in the model.', a);
    end
    v = mod(floor(reg(a) / 2^lo), 2^(hi - lo + 1));
end

function words = pack_setup_words(reg, W, cfg)
    words = zeros(1,16);
    INIT  = init_values(reg, cfg);
    for k = 1:numel(W)
        b = 0;
        for j = 1:numel(W(k).src)
            n = W(k).dst{j}(1) - W(k).dst{j}(2) + 1;
            if strlength(W(k).init(j)) > 0
                v = INIT.(char(W(k).init(j)));
            else
                v = field_value_bound(reg, W(k).addr(j), W(k).hi(j), W(k).lo(j));
            end
            b = b + mod(v, 2^n) * 2^W(k).dst{j}(2);
        end
        words(W(k).idx + 1) = b;
    end
end

function INIT = init_values(reg, cfg)
% The four "Set per Init N calculation" values.  They differ from steady state
% only when Wide BW mode widens the loop during the transition; with Wide BW
% off, the ADI driver sets Init = steady state, which is cfg.initMode='copy'.
    switch lower(cfg.initMode)
        case 'zero'
            INIT = struct('CP',0,'R3',0,'C3',0,'R1',0);
        case 'manual'
            INIT = cfg.initValues;
        case 'copy'
            INIT.CP = field_value_bound(reg, "27B", 5, 0);
            INIT.R3 = field_value_bound(reg, "280", 3, 0);
            INIT.C3 = field_value_bound(reg, "27F", 3, 0);
            INIT.R1 = field_value_bound(reg, "27F", 7, 4);
        otherwise
            error('Unknown initMode "%s".', cfg.initMode);
    end
end

%% ========================================================================
%  ROUND-TRIP VERIFICATION
% ========================================================================
function verify_roundtrip(Bytes, Y, names, P, W, cfg)
%VERIFY_ROUNDTRIP  Read-only check on the packing step.
%
% Unpacks the 16 generated bytes back into field values and compares them
% against the values that went in.  It modifies nothing, rounds nothing and
% discards nothing - Y is already rounded and range-checked by the time this
% runs.  It is a mirror held up to pack_setup_words.
%
% Two things make it fail, and both mean the hex would be wrong:
%
%   1. The Table 12 map is wrong - a bad shift or bit range.  Every field in
%      the affected word comes back different.
%   2. A field value does not fit the destination slot in the setup word, so
%      packing truncated it with mod().  Only that field comes back different,
%      and the value read back is the input modulo the slot width.
%
% Case 2 is the one worth caring about: it is the packing stage quietly
% mangling a value, which is exactly what a profile generator must never do
% silently.  cfg.verifyStrict = false downgrades the stop to a warning.
    nF   = numel(names);
    seen = false(nF,1);
    bad  = strings(0,1);

    for k = 1:numel(W)
        for j = 1:numel(W(k).src)
            if strlength(W(k).init(j)) > 0, continue; end
            i = find(P.addr == W(k).addr(j) & P.hi == W(k).hi(j) & ...
                     P.lo == W(k).lo(j), 1);
            if isempty(i), continue; end
            seen(i) = true;

            n   = W(k).dst{j}(1) - W(k).dst{j}(2) + 1;
            got = mod(floor(Bytes(:, W(k).idx+1) / 2^W(k).dst{j}(2)), 2^n);
            if any(got ~= Y(:,i))
                r = find(got ~= Y(:,i), 1);
                why = ternary(Y(r,i) >= 2^n, ...
                    sprintf('value needs %d bits, word %X gives it %d', ...
                            ceil(log2(Y(r,i)+1)), W(k).idx, n), ...
                    sprintf('word %X bits [%d:%d] disagree', W(k).idx, ...
                            W(k).dst{j}(1), W(k).dst{j}(2)));
                bad(end+1) = sprintf('%s (row %d: sent %d, read back %d - %s)', ...
                    names(i), r, Y(r,i), got(r), why);              %#ok<AGROW>
            end
        end
    end

    fprintf('\nRound trip: %d of %d field(s) carried by Table 12', sum(seen), nF);
    if any(~seen)
        fprintf(', unused: %s', strjoin(names(~seen), ', '));
    end
    fprintf('\n');

    if isempty(bad)
        fprintf('All carried fields unpack bit-exact.\n');
    else
        msg = sprintf('Round-trip mismatch:\n  %s', strjoin(bad, sprintf('\n  ')));
        if cfg.verifyStrict, error('%s', msg); else, warning('%s', msg); end
    end
end

function check_coverage(P, W)
    needed = strings(0,1);
    for k = 1:numel(W)
        for j = 1:numel(W(k).src)
            s = string(W(k).src{j});
            if startsWith(s, "#"), continue; end
            needed(end+1) = extractBefore(upper(s), "[");   %#ok<AGROW>
        end
    end
    needed = unique(needed);
    have   = unique(P.addr);
    miss   = setdiff(needed, have);

    fprintf('\nTable 12 needs %d source register(s); model provides %d.\n', ...
        numel(needed), numel(have));
    if ~isempty(miss)
        error(['Missing register(s) required by Table 12: %s\n' ...
               'The profile cannot be assembled without them.'], ...
               strjoin("0x" + miss, ', '));
    end
    fprintf('All required registers are present.\n');
    extra = setdiff(have, needed);
    if ~isempty(extra)
        fprintf('Model also holds %s (not used by Tx Fast Lock).\n', ...
            strjoin("0x" + extra, ', '));
    end
end

%% ========================================================================
%  OUTPUT
% ========================================================================
function write_profile_table(outFile, fq, Bytes, cfg)
    nF   = numel(fq);
    pre  = string(cfg.hexPrefix);
    slot = mod((0:nF-1).', cfg.numProfiles);

    T = table(fq, slot, 'VariableNames', {'Frequency_MHz','Profile'});
    for k = 1:16
        % dec2hex returns an n-by-2 CHAR MATRIX; string() turns each row into
        % one element, giving an n-by-1 string array.
        hx = string(dec2hex(Bytes(:,k), 2));
        hx = hx(:);
        if numel(hx) ~= nF
            error('Hex conversion produced %d value(s) for %d row(s).', ...
                numel(hx), nF);
        end
        T.(sprintf('W%X', k-1)) = pre + hx;
    end

    blob = strings(nF,1);
    for r = 1:nF
        blob(r) = string(sprintf('%02X', Bytes(r,:)));
    end
    T.Profile_Hex32 = blob;

    safe_writetable(T, outFile);

    % ---- read the file back and confirm it holds what we sent ------------
    if cfg.verifyWrite
        verify_written_file(outFile, Bytes, cfg);
    end

    fprintf('\nProfiles written to %s  (%d rows, 16 setup words each, hex)\n', ...
        outFile, nF);
    if isempty(cfg.hexPrefix)
        fprintf(['Note: bare hex. Excel will read values like 00 or 12 as ' ...
                 'numbers.\n      Set cfg.hexPrefix = ''0x'' if that is a ' ...
                 'problem.\n']);
    end
    if nF > cfg.numProfiles
        fprintf(['Note: the part holds %d profiles at a time. The Profile ' ...
                 'column cycles 0..%d;\n      the host must load the right ' ...
                 'row before each retune.\n'], cfg.numProfiles, cfg.numProfiles-1);
    end
end

function safe_writetable(T, outFile)
% writetable to an EXISTING .xlsx overwrites cell by cell and leaves whatever
% was in the sheet beyond the new extent.  If the previous run had more rows,
% or different columns, the leftovers stay and the file reads back as a mix of
% old and new.  Deleting first is the only way to get a clean sheet.
%
% It also fails if the file is open in Excel, which produces a confusing error
% at the very end of a long run, so that case is named explicitly.
    if isfile(outFile)
        [ok, msg] = delete_file(outFile);
        if ~ok
            error(['Cannot replace "%s": %s\n' ...
                   'If the file is open in Excel, close it and run again. ' ...
                   'Writing over an\nexisting sheet leaves stale rows behind, ' ...
                   'so this is not skipped.'], outFile, msg);
        end
    end

    try
        writetable(T, outFile);
    catch ME
        error(['Could not write "%s": %s\n' ...
               'The usual cause is the file being open in Excel.'], ...
               outFile, ME.message);
    end
end

function [ok, msg] = delete_file(f)
    ok = true;  msg = '';
    try
        delete(f);
    catch ME
        ok = false;  msg = ME.message;
    end
    if isfile(f)
        ok = false;
        if isempty(msg), msg = 'the file is locked by another program'; end
    end
end

function verify_written_file(outFile, Bytes, cfg)
% Read the spreadsheet back and compare against what was sent.  This is the
% only check that covers the file layer itself: string columns that Excel
% coerced, a stale sheet, a truncated write.
    try
        T = readtable(outFile, 'VariableNamingRule','preserve', 'TextType','string');
    catch ME
        warning('Wrote "%s" but could not read it back: %s', outFile, ME.message);
        return;
    end

    vn = string(T.Properties.VariableNames);
    got = zeros(height(T), 16);
    for k = 1:16
        c = find(vn == "W" + string(dec2hex(k-1,1)), 1);
        if isempty(c)
            warning('Column W%X is missing from the written file.', k-1);
            return;
        end
        col = T{:, c};
        if isnumeric(col)
            got(:,k) = double(col);
        else
            s = erase(erase(strtrim(string(col)), "0x"), "0X");
            got(:,k) = arrayfun(@(x) hex2dec(char(x)), s);
        end
    end

    if height(T) ~= size(Bytes,1)
        warning(['"%s" has %d row(s) but %d were sent. The sheet may hold ' ...
                 'stale rows from an earlier run.'], outFile, height(T), size(Bytes,1));
        return;
    end

    d = got - Bytes;
    if any(d(:) ~= 0)
        [r, c] = find(d ~= 0, 1);
        msg = sprintf(['Readback mismatch in "%s": row %d, W%X sent 0x%02X, ' ...
                       'file holds 0x%02X.'], outFile, r, c-1, Bytes(r,c), got(r,c));
        if cfg.verifyStrict, error('%s', msg); else, warning('%s', msg); end
    else
        fprintf('Readback check: the file matches the generated bytes.\n');
    end
end

function diagnose_stage(label, M, names)
    nz  = nnz(M);
    tot = numel(M);
    fprintf('%-24s  %d/%d non-zero, range %g .. %g\n', ...
        label, nz, tot, min(M(:)), max(M(:)));
    if nz < tot
        allZeroCols = find(all(M == 0, 1));
        if ~isempty(allZeroCols) && numel(allZeroCols) < numel(names)
            fprintf('%-24s  all-zero: %s\n', '', ...
                strjoin(names(allZeroCols), ', '));
        end
    end
end

function write_field_dump(dumpFile, fq, Y, names)
    T = array2table([fq, double(Y)], 'VariableNames', ...
        matlab.lang.makeValidName(["Frequency_MHz", names']));
    safe_writetable(T, dumpFile);
    fprintf('Rounded field values written to %s\n', dumpFile);
end

%% ========================================================================
%  SELF-TEST
% ========================================================================
function selftest()
%SELFTEST  Run the packer on known values, with no model and no data files.
%
% Separates "the code is broken" from "my model is broken". If this passes,
% the packing, the Table 12 map and the spreadsheet writer all work, and an
% all-zero profiles.xlsx is coming from the model.
    fprintf('Catalina_MakeProfiles self-test\n');
    fprintf('%s\n', repmat('=', 1, 46));

    names = ["271[7:0]";"272[2:0]";"273[7:0]";"274[7:0]";"275[6:0]"; ...
             "282[2:0]";"282[4:3]";"279[3:0]";"27B[5:0]";"280[3:0]"; ...
             "27F[3:0]";"27F[7:4]";"27E[3:0]";"27E[7:4]";"290[6:4]"; ...
             "005[7:4]";"278[6:3]";"291[3:0]";"277[7:0]";"276[6:0]";"278[0]"];
    vals  = [0x2C 0x01 0x66 0x66 0x06 0x05 0x02 0x09 0x1E 0x0C ...
             0x07 0x0B 0x0D 0x04 0x03 0x01 0x0A 0x0E 0x93 0x5A 0x01];

    cfg = profile_config();
    P   = parse_fields(names);
    W   = fastlock_map();
    check_coverage(P, W);
    W   = bind_map(W, P);

    reg   = assemble_registers(vals, P);
    Bytes = pack_setup_words(reg, W, cfg);

    fprintf('\nAssembled registers:\n');
    for a = sort(string(reg.keys))
        fprintf('  0x%-4s = 0x%02X\n', a, reg(char(a)));
    end

    fprintf('\nPacked words:\n');
    for k = 1:16
        fprintf('  W%X = 0x%02X   %s\n', k-1, Bytes(k), W(k).desc);
    end

    % ---- hand-computed expectations ------------------------------------
    exp = struct( ...
        'w', {1, 6, 11, 16}, ...
        'v', {0x2C, ...                                   % 271 straight through
              bitor(bitshift(0x05,4), 0x09), ...          % BiasRef<<4 | Varactor
              bitor(bitshift(0x0D,4), 0x04), ...          % C1<<4 | C2  (nibble swap)
              bitor(bitshift(0x5A,1), 0x01)}, ...         % ALC<<1 | Tune[8]
        'why', {'271 passes through', 'W5 nibble pair', ...
                'WA nibble swap of 0x27E', 'WF ALC shifted left by 1'});

    fprintf('\nHand checks:\n');
    ok = true;
    for e = exp
        got = Bytes(e.w);
        pass = got == e.v;
        ok = ok && pass;
        fprintf('  W%X  expect 0x%02X  got 0x%02X  %-6s  %s\n', ...
            e.w-1, e.v, got, string(ternary(pass,'PASS','FAIL')), e.why);
    end

    verify_roundtrip(Bytes, vals, names, P, W, cfg);

    % ---- spreadsheet layer ----------------------------------------------
    tmp = [tempname '.xlsx'];
    write_profile_table(tmp, 2400, Bytes, cfg);
    delete(tmp);

    fprintf('\n%s\n', repmat('=', 1, 46));
    if ok
        fprintf(['Self-test PASSED. Packing, the Table 12 map and the writer ' ...
                 'all work.\nIf your profiles.xlsx is zeros, the zeros come ' ...
                 'from the model:\ncheck profiles_fields.xlsx, or rebuild with ' ...
                 'Catalina_Project option 2.\n']);
    else
        fprintf('Self-test FAILED. The packing stage itself is wrong.\n');
    end
end

function print_preview(fq, Bytes, W, cfg)
    n = min(cfg.previewRows, numel(fq));
    fprintf('\nPreview (first %d target(s)):\n\n%12s ', n, 'f [MHz]');
    for k = 1:16, fprintf(' W%X', k-1); end
    fprintf('\n%s\n', repmat('-', 1, 12 + 16*3 + 1));
    for r = 1:n
        fprintf('%12.4f ', fq(r));
        fprintf(' %02X', Bytes(r,:));
        fprintf('\n');
    end
    fprintf('\nWord map:\n');
    for k = 1:numel(W)
        fprintf('  W%X  %s\n', W(k).idx, W(k).desc);
    end
end

%% ========================================================================
%  MODEL EVALUATION
% ========================================================================
function Model = override_fref(Model, fref)
    for k = 1:numel(Model.segment)
        for i = 1:numel(Model.segment(k).field)
            if strcmp(Model.segment(k).field(i).kind, 'formula')
                Model.segment(k).field(i).fref = fref;
            end
        end
    end
    fprintf('Reference frequency overridden to %g MHz.\n', fref);
end

function report_formula_fields(Model)
    kinds = arrayfun(@(i) string(Model.segment(1).field(i).kind), ...
                     1:numel(Model.names));
    isF = kinds == "formula";
    if ~any(isF)
        fprintf(['\nNo formula fields in this model: the synthesizer words ' ...
                 'will be interpolated,\nwhich is unlikely to be correct. See ' ...
                 'cfg.formulaFields in Catalina_Project.m\n']);
        return;
    end
    fref = Model.segment(1).field(find(isF,1)).fref;
    fprintf('\nComputed from the equation (f_ref = %g MHz): %s\n', ...
        fref, strjoin(Model.names(isF), ', '));
    fprintf('Interpolated from the model: %s\n', strjoin(Model.names(~isF), ', '));
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

function v = eval_field(f, fq)
% Must stay identical to eval_field in Catalina_Project.m.
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
            % round the wide value first, then slice out this member's bits
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

function y = barycentric_cheb(xn, yn, xe)
    xn = xn(:);  yn = yn(:);  xe = xe(:);
    d = numel(xn);
    w = ((-1).^(0:d-1).') .* sin(((1:d).' - 0.5)*pi/d);
    D = xe.' - xn;
    hit = abs(D) < 1e-12;  D(hit) = 1;
    y = (((w.*yn).' * (1./D)) ./ (w.' * (1./D))).';
    [ri, ci] = find(hit);  y(ci) = yn(ri);
end

%% ========================================================================
%  HELPERS
% ========================================================================
function check_model_compatible(Model, modelFile)
%CHECK_MODEL_COMPATIBLE  Refuse a model this file cannot evaluate.
%
% The model format grew: formula fields (tag/fref/vcoMin), composite fields
% (cshift/cwidth) and residual corrections (rF/rV) were all added after the
% first version.  A model saved before those exist loads without complaint and
% then fails deep inside eval_field, or worse, silently evaluates a field by
% the wrong rule.  Check up front and say what to do.
    need = {'kind','name','val','edgeF','edgeV','x_cheb','Y_cheb','err', ...
            'tag','fref','vcoMin','cshift','cwidth','rF','rV'};

    if ~isfield(Model,'segment') || isempty(Model.segment) || ...
       ~isfield(Model.segment,'field')
        error(['"%s" does not look like a model built by Catalina_Project. ' ...
               'Rebuild it with option 2.'], modelFile);
    end

    have    = fieldnames(Model.segment(1).field);
    missing = setdiff(need, have);
    if ~isempty(missing)
        error(['"%s" was built by an older version of Catalina_Project and is ' ...
               'missing:\n    %s\nRebuild it (Catalina_Project, option 2). ' ...
               'The stored bit widths may also be\nwrong, so reusing the old ' ...
               'file is not safe even if these were patched in.'], ...
               modelFile, strjoin(missing, ', '));
    end

    % every kind present must be one eval_field knows
    known = {'formula','const','lookup','exact','smooth','composite'};
    kinds = {};
    for k = 1:numel(Model.segment)
        kinds = [kinds, {Model.segment(k).field.kind}];            %#ok<AGROW>
    end
    bad = setdiff(unique(kinds), known);
    if ~isempty(bad)
        error('"%s" contains unknown field kind(s): %s', ...
              modelFile, strjoin(bad, ', '));
    end

    % a formula field is useless without a reference frequency
    for k = 1:numel(Model.segment)
        for i = 1:numel(Model.segment(k).field)
            f = Model.segment(k).field(i);
            if strcmp(f.kind,'formula') && (isempty(f.fref) || ~isfinite(f.fref))
                error(['Formula field %s in segment %d has no reference ' ...
                       'frequency.\nSet cfg.fref in Catalina_Project and ' ...
                       'rebuild.'], f.name, k);
            end
        end
    end

    nk = numel(unique(kinds));
    fprintf('Model format OK: %d segment(s), %d field kind(s) in use (%s).\n', ...
        numel(Model.segment), nk, strjoin(unique(kinds), ', '));
end

function fq = read_target_frequencies(freqFile)
    if ~isfile(freqFile), error('Target file "%s" not found.', freqFile); end
    raw = readmatrix(freqFile);
    if isempty(raw), error('No numeric data found in "%s".', freqFile); end
    fq = raw(:,1);
    fq = fq(~isnan(fq));
    if isempty(fq)
        error('First column of "%s" contains no frequencies.', freqFile);
    end
    fq = sort(fq);
end

function s = replace_ext(f, newExt)
    [p, n, ~] = fileparts(f);
    s = fullfile(p, [n newExt]);
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
