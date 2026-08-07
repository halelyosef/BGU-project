%% ========================================================================
%  AD9361 Fast Lock profile inspector
%
%  Reads a profile table (or a raw hex blob, or a generated SPI script) and
%  says what is actually in it: decodes the 16 setup words back into named
%  register fields, reconstructs the LO frequency each profile will tune to,
%  and flags anything suspicious.
%
%  Usage
%  -----
%     Catalina_InspectProfiles                                % profiles.xlsx
%     Catalina_InspectProfiles('profiles.xlsx')
%     Catalina_InspectProfiles('profiles.xlsx', 40)           % f_ref = 40 MHz
%     Catalina_InspectProfiles('fastlock_load.txt')           % an SPI script
%     Catalina_InspectProfiles('F000CACC0C501810BBFF3CEE71E95DAE')
%
%  WHY THIS IS WORTH RUNNING
%  -------------------------
%  The frequency check is the point.  Words 0-4 carry the synthesizer integer
%  and fractional words and word C carries the VCO divider, so the frequency a
%  profile will actually produce can be computed from the bytes alone:
%
%      f_RF = (INT + FRAC/8388593) * f_ref / 2^(d+2)
%
%  That is completely independent of the model, the interpolation, the
%  Chebyshev stage and the packing code.  If it agrees with the frequency the
%  profile is labelled with, then the whole chain got the synthesizer right.
%  If it disagrees, the profile will tune somewhere other than where you think,
%  and no amount of checking the analog fields would have told you.
%
%  Everything else here is a consistency check on the remaining fields: fields
%  stuck at a constant, fields pinned at 0 or at their maximum, values that do
%  not fit their bit width, and jumps between adjacent profiles.
%
%  A NOTE ON THE UNPACKING MAP
%  ---------------------------
%  unpack_map() below is written independently from fastlock_map() in
%  Catalina_MakeProfiles.m, from the same Table 12.  That is deliberate: two
%  independent transcriptions that agree are evidence the transcription is
%  right, whereas a shared one would simply repeat any mistake.  If you edit
%  one you must edit the other, and the frequency check will usually catch it
%  if you forget.
% ========================================================================

function Catalina_InspectProfiles(src, fref)
    if nargin < 1 || isempty(src),  src  = 'profiles.xlsx'; end
    if nargin < 2 || isempty(fref), fref = 40;              end

    cfg = inspect_config();
    cfg.fref = fref;

    [fq, Bytes, srcKind] = read_any(src);
    nP = size(Bytes,1);

    fprintf('%s\n', repmat('=', 1, 74));
    fprintf('  AD9361 Fast Lock profile inspection\n');
    fprintf('  source: %s   (%s)   %d profile(s)   f_ref = %g MHz\n', ...
        src, srcKind, nP, cfg.fref);
    fprintf('%s\n', repmat('=', 1, 74));

    U = unpack_map();
    F = decode_all(Bytes, U);            % struct of [nP x 1] field vectors

    report_frequency(fq, F, cfg);
    report_fields(F, U, cfg);
    report_anomalies(fq, F, U, cfg);
    report_raw(fq, Bytes, cfg);

    if cfg.plot && nP >= 3
        plot_fields(fq, F, U);
    end
end

%% ========================================================================
%  CONFIGURATION
% ========================================================================
function cfg = inspect_config()
    cfg.fref      = 40;        % MHz
    cfg.modulus   = 8388593;   % AD9361 fractional-N modulus
    cfg.freqTolHz = 100;       % decoded vs labelled disagreement worth flagging
    cfg.jumpFrac  = 0.5;       % jump > this fraction of a field's range is odd
    cfg.rawRows   = 8;         % rows of raw hex to print
    cfg.plot      = true;
end

%% ========================================================================
%  TABLE 12, INVERTED
% ========================================================================
function U = unpack_map()
%UNPACK_MAP  Where each field lives inside the 16 setup words.
%
% Independent transcription of UG-570 Table 12.  Each entry says: field NAME
% occupies bits [hi lo] of setup word WORD, and is WIDTH bits wide.
%
% Word D uses [7:4] for the VCO cal offset.  UG-570 says "shift left by 1",
% which would put it at [4:1] and overlap VCO Varactor Reference at [3:0].
% [7:4] is the only non-overlapping reading and matches words 8, 9 and B,
% which all pair a high nibble with a low nibble.
    U = struct('name',{}, 'word',{}, 'hi',{}, 'lo',{}, 'meaning',{});

    U = u(U, "271[7:0]",  0,  7, 0, 'Synth integer word [7:0]');
    U = u(U, "272[2:0]",  1,  2, 0, 'Synth integer word [10:8]');
    U = u(U, "273[7:0]",  2,  7, 0, 'Synth fractional word [7:0]');
    U = u(U, "274[7:0]",  3,  7, 0, 'Synth fractional word [15:8]');
    U = u(U, "275[6:0]",  4,  6, 0, 'Synth fractional word [22:16]');
    U = u(U, "282[2:0]",  5,  6, 4, 'VCO bias ref');
    U = u(U, "279[3:0]",  5,  3, 0, 'VCO varactor');
    U = u(U, "282[4:3]",  6,  7, 6, 'VCO bias Tcf');
    U = u(U, "#CP",       6,  5, 0, 'Charge pump current (Init)');
    U = u(U, "27B[5:0]",  7,  5, 0, 'Charge pump current');
    U = u(U, "280[3:0]",  8,  7, 4, 'Loop filter R3');
    U = u(U, "#R3",       8,  3, 0, 'Loop filter R3 (Init)');
    U = u(U, "27F[3:0]",  9,  7, 4, 'Loop filter C3');
    U = u(U, "#C3",       9,  3, 0, 'Loop filter C3 (Init)');
    U = u(U, "27E[3:0]", 10,  7, 4, 'Loop filter C1');
    U = u(U, "27E[7:4]", 10,  3, 0, 'Loop filter C2');
    U = u(U, "27F[7:4]", 11,  7, 4, 'Loop filter R1');
    U = u(U, "#R1",      11,  3, 0, 'Loop filter R1 (Init)');
    U = u(U, "290[6:4]", 12,  6, 4, 'VCO varactor ref Tcf');
    U = u(U, "005[7:4]", 12,  3, 0, 'VCO divider');
    U = u(U, "278[6:3]", 13,  7, 4, 'VCO cal offset');
    U = u(U, "291[3:0]", 13,  3, 0, 'VCO varactor reference');
    U = u(U, "277[7:0]", 14,  7, 0, 'Force VCO tune [7:0]');
    U = u(U, "276[6:0]", 15,  7, 1, 'Force ALC word');
    U = u(U, "278[0]",   15,  0, 0, 'Force VCO tune [8]');
end

function U = u(U, name, word, hi, lo, meaning)
    n = numel(U) + 1;
    U(n).name    = name;
    U(n).word    = word;
    U(n).hi      = hi;
    U(n).lo      = lo;
    U(n).meaning = meaning;
end

function F = decode_all(Bytes, U)
% Pull every field out of the setup words.  Returned as a struct keyed by a
% sanitised field name, each holding an [nP x 1] column.
    F = struct();
    for k = 1:numel(U)
        n = U(k).hi - U(k).lo + 1;
        v = mod(floor(Bytes(:, U(k).word + 1) / 2^U(k).lo), 2^n);
        F.(key(U(k).name)) = v;
    end
end

function s = key(name)
    s = char(matlab.lang.makeValidName(name));
end

%% ========================================================================
%  THE FREQUENCY CHECK
% ========================================================================
function report_frequency(fq, F, cfg)
    INT  = F.(key("271[7:0]")) + 256*F.(key("272[2:0]"));
    FRAC = F.(key("273[7:0]")) + 256*F.(key("274[7:0]")) + 65536*F.(key("275[6:0]"));
    d    = F.(key("005[7:4]"));

    N     = INT + FRAC / cfg.modulus;
    f_vco = N * cfg.fref / 2;
    f_rf  = f_vco ./ 2.^(d + 1);

    fprintf('\nFrequency decoded from the setup words\n');
    fprintf('  f_RF = (INT + FRAC/%d) * %g / 2^(d+2)\n\n', cfg.modulus, cfg.fref);
    fprintf('  %3s %13s %13s %12s %9s %5s %6s %10s\n', ...
        '#','labelled MHz','decoded MHz','error Hz','ppm','div','INT','FRAC');
    fprintf('  %s\n', repmat('-', 1, 78));

    haveLabel = ~isempty(fq) && all(~isnan(fq));
    bad = false(numel(f_rf),1);
    for r = 1:numel(f_rf)
        if haveLabel
            eHz  = (f_rf(r) - fq(r)) * 1e6;
            ppm  = 1e6 * (f_rf(r) - fq(r)) / max(fq(r), eps);
            bad(r) = abs(eHz) > cfg.freqTolHz;
            fprintf('  %3d %13.6f %13.6f %12.1f %9.3f %5d %6d %10d%s\n', ...
                r, fq(r), f_rf(r), eHz, ppm, d(r), INT(r), FRAC(r), ...
                ternary(bad(r), '  <-- MISMATCH', ''));
        else
            fprintf('  %3d %13s %13.6f %12s %9s %5d %6d %10d\n', ...
                r, '-', f_rf(r), '-', '-', d(r), INT(r), FRAC(r));
        end
    end
    fprintf('  %s\n', repmat('-', 1, 78));

    if ~haveLabel
        fprintf(['  No labelled frequencies in the source, so only the decoded ' ...
                 'values are shown.\n']);
        return;
    end

    if any(bad)
        fprintf(['  %d profile(s) decode to a different frequency than they are ' ...
                 'labelled with.\n'], sum(bad));
        fprintf(['  This is independent of the model and of the packing code, so ' ...
                 'one of three\n  things is wrong: f_ref is not %g MHz, the ' ...
                 'divider convention differs, or\n  the synthesizer words were ' ...
                 'not written correctly.\n'], cfg.fref);
        checkFref = median(fq ./ max(f_rf, eps));
        if abs(checkFref - 1) > 1e-6 && abs(checkFref - round(checkFref*8)/8) < 1e-4
            fprintf(['  The decoded frequencies are a consistent factor %.4f off. ' ...
                     'Try f_ref = %g MHz.\n'], checkFref, cfg.fref*checkFref);
        end
    else
        fprintf(['  All %d profile(s) decode to their labelled frequency within ' ...
                 '%g Hz.\n'], numel(f_rf), cfg.freqTolHz);
        fprintf(['  The synthesizer path is correct end to end: model, packing ' ...
                 'and hex all agree.\n']);
    end
end

%% ========================================================================
%  FIELD SUMMARY
% ========================================================================
function report_fields(F, U, cfg) %#ok<INUSD>
    fprintf('\nField values across all profiles\n');
    fprintf('  %-12s %-28s %4s %6s %6s %7s %8s\n', ...
        'Field','Meaning','bits','min','max','unique','span');
    fprintf('  %s\n', repmat('-', 1, 80));

    for k = 1:numel(U)
        v  = F.(key(U(k).name));
        w  = U(k).hi - U(k).lo + 1;
        nu = numel(unique(v));
        sp = max(v) - min(v);
        fprintf('  %-12s %-28s %4d %6d %6d %7d %8d\n', ...
            U(k).name, U(k).meaning, w, min(v), max(v), nu, sp);
    end
    fprintf('  %s\n', repmat('-', 1, 80));
end

%% ========================================================================
%  ANOMALIES
% ========================================================================
function report_anomalies(fq, F, U, cfg)
    fprintf('\nChecks\n');
    issues = strings(0,1);
    nP = numel(F.(key(U(1).name)));

    for k = 1:numel(U)
        v    = F.(key(U(k).name));
        w    = U(k).hi - U(k).lo + 1;
        vmax = 2^w - 1;
        nm   = U(k).name;

        % ---- pinned at a rail --------------------------------------------
        if all(v == 0) && nP > 1
            issues(end+1) = sprintf(['%-12s is zero in every profile. Either the ' ...
                'field is genuinely unused, or\n               it was never ' ...
                'written.'], nm);                                   %#ok<AGROW>
        elseif all(v == vmax) && nP > 1
            issues(end+1) = sprintf(['%-12s is at its maximum (%d) in every ' ...
                'profile, which usually means a\n               value overflowed ' ...
                'and saturated.'], nm, vmax);                       %#ok<AGROW>
        end

        % ---- does not fit its own width ----------------------------------
        if any(v > vmax)
            issues(end+1) = sprintf(['%-12s holds %d which does not fit %d ' ...
                'bits.'], nm, max(v), w);                           %#ok<AGROW>
        end

        % ---- discontinuity between adjacent profiles ---------------------
        % Only meaningful if the profiles are ordered by frequency and the
        % field is analog rather than a discrete control.
        if nP >= 3 && ~isempty(fq) && issorted(fq)
            rng_ = max(v) - min(v);
            if rng_ > 4
                dv = abs(diff(v));
                big = find(dv > cfg.jumpFrac * rng_);
                if ~isempty(big)
                    issues(end+1) = sprintf(['%-12s jumps by %d between %.4g and ' ...
                        '%.4g MHz (range is %d).\n               Expected at a ' ...
                        'band edge, otherwise suspect.'], nm, max(dv(big)), ...
                        fq(big(1)), fq(big(1)+1), rng_);            %#ok<AGROW>
                end
            end
        end
    end

    % ---- Init vs steady state ------------------------------------------
    pairs = {"#CP","27B[5:0]"; "#R3","280[3:0]"; "#C3","27F[3:0]"; "#R1","27F[7:4]"};
    same  = true(size(pairs,1),1);
    for p = 1:size(pairs,1)
        same(p) = isequal(F.(key(pairs{p,1})), F.(key(pairs{p,2})));
    end
    if all(same)
        fprintf(['  Init values equal steady-state values in every profile: ' ...
                 'Wide BW is off.\n']);
    elseif any(same)
        issues(end+1) = sprintf(['Init values match steady state for some fields ' ...
            'but not others (%s differ).\n               That mixture is unusual ' ...
            '- Wide BW is normally all or nothing.'], ...
            strjoin(string(pairs(~same,1)), ', '));                 %#ok<AGROW>
    else
        fprintf('  Init values differ from steady state: Wide BW mode is in use.\n');
    end

    % ---- duplicate profiles ---------------------------------------------
    B = cell2mat(arrayfun(@(k) F.(key(U(k).name)), 1:numel(U), 'UniformOutput', false));
    [~, ia, ic] = unique(B, 'rows', 'stable');
    if numel(ia) < nP
        dup = find(accumarray(ic,1) > 1);
        issues(end+1) = sprintf(['%d profile(s) are byte-identical to another. ' ...
            'Two target frequencies\n               that produce the same profile ' ...
            'is normal only if they are very close.'], nP - numel(ia)); %#ok<AGROW>
        if ~isempty(dup) && ~isempty(fq)
            g = find(ic == dup(1));
            issues(end) = issues(end) + sprintf('\n               First pair: %.6g and %.6g MHz.', ...
                fq(g(1)), fq(g(2)));
        end
    end

    if isempty(issues)
        fprintf('  No anomalies found.\n');
    else
        fprintf('\n');
        for i = 1:numel(issues)
            fprintf('  [!] %s\n', issues(i));
        end
    end
end

%% ========================================================================
%  RAW BYTES
% ========================================================================
function report_raw(fq, Bytes, cfg)
    n = min(cfg.rawRows, size(Bytes,1));
    fprintf('\nRaw setup words (first %d)\n', n);
    fprintf('  %12s ', 'f [MHz]');
    for k = 1:16, fprintf(' W%X', k-1); end
    fprintf('\n  %s\n', repmat('-', 1, 12 + 16*3 + 1));
    for r = 1:n
        if isempty(fq) || isnan(fq(r))
            fprintf('  %12s ', '-');
        else
            fprintf('  %12.4f ', fq(r));
        end
        fprintf(' %02X', Bytes(r,:));
        fprintf('\n');
    end
end

%% ========================================================================
%  PLOT
% ========================================================================
function plot_fields(fq, F, U)
    keep = true(numel(U),1);
    for k = 1:numel(U)
        keep(k) = numel(unique(F.(key(U(k).name)))) > 1;
    end
    idx = find(keep);
    if isempty(idx), return; end

    per = 6;
    for a = 1:per:numel(idx)
        b = min(a+per-1, numel(idx));
        figure('Name','Profile field values','Color','w', ...
               'Position',[70 70 1150 150*(b-a+1)]);
        tiledlayout(b-a+1, 1, 'TileSpacing','compact', 'Padding','compact');
        for i = a:b
            k = idx(i);
            ax = nexttile;
            v = F.(key(U(k).name));
            if isempty(fq) || any(isnan(fq))
                plot(ax, v, '.-', 'LineWidth', 1.2, 'MarkerSize', 12, ...
                     'Color', [0.10 0.35 0.75]);
                xlabel(ax, 'profile index');
            else
                plot(ax, fq, v, '.-', 'LineWidth', 1.2, 'MarkerSize', 12, ...
                     'Color', [0.10 0.35 0.75]);
                if i == b, xlabel(ax, 'frequency [MHz]'); end
            end
            title(ax, sprintf('%s  -  %s', U(k).name, U(k).meaning), ...
                  'Interpreter','none', 'FontWeight','normal');
            ylabel(ax, 'counts');
            grid(ax,'on');  box(ax,'on');
            ax.GridAlpha = 0.12;  ax.FontSize = 10;
        end
    end
end

%% ========================================================================
%  INPUT - spreadsheet, hex blob, or a generated SPI script
% ========================================================================
function [fq, Bytes, kind] = read_any(src)
    src = string(src);

    % ---- a bare hex blob passed directly ---------------------------------
    if ~isfile(src)
        s = char(upper(regexprep(src, '[^0-9A-Fa-f]', '')));
        if numel(s) == 32
            fq    = NaN;
            Bytes = hex2dec(reshape(s, 2, 16).').';
            kind  = 'hex blob';
            return;
        end
        error(['"%s" is not a file, and not a 32-character hex blob either.'], src);
    end

    [~,~,ext] = fileparts(src);
    switch lower(ext)
        case {'.xlsx','.xls','.csv'}
            [fq, Bytes] = read_profile_table(src);
            kind = 'profile table';
        case {'.txt','.log'}
            [fq, Bytes] = read_spi_script(src);
            kind = 'SPI script';
        otherwise
            error('Unrecognised file type "%s".', ext);
    end
end

function [fq, Bytes] = read_profile_table(profileFile)
    T  = readtable(profileFile, 'VariableNamingRule','preserve', 'TextType','string');
    vn = string(T.Properties.VariableNames);

    fcol = find(contains(lower(vn), 'freq'), 1);
    if isempty(fcol)
        fq = nan(height(T),1);
    else
        fq = T{:, fcol};
        if ~isnumeric(fq), fq = str2double(string(fq)); end
    end

    wcol  = arrayfun(@(k) find(vn == "W" + string(dec2hex(k,1)), 1), 0:15, ...
                     'UniformOutput', false);
    if all(~cellfun(@isempty, wcol))
        Bytes = zeros(height(T), 16);
        for k = 1:16
            Bytes(:,k) = hexcol_to_num(T{:, wcol{k}});
        end
    else
        bcol = find(contains(lower(vn),'hex32') | contains(lower(vn),'blob'), 1);
        if isempty(bcol)
            error('"%s" has neither W0..WF nor a Profile_Hex32 column.', profileFile);
        end
        blob  = string(T{:, bcol});
        Bytes = zeros(numel(blob), 16);
        for r = 1:numel(blob)
            s = char(strtrim(blob(r)));
            if numel(s) ~= 32
                error('Row %d: expected 32 hex characters, found %d.', r, numel(s));
            end
            Bytes(r,:) = hex2dec(reshape(s, 2, 16).');
        end
    end

    keep  = ~any(isnan(Bytes), 2);
    fq    = fq(keep);
    Bytes = Bytes(keep,:);
end

function [fq, Bytes] = read_spi_script(txtFile)
% Replay a generated SPIwrite script and recover the words it would program.
% Useful for checking a script by hand, or for inspecting one that came from
% somewhere other than this toolchain.
    txt = string(splitlines(fileread(txtFile)));

    tok = regexp(txt, '^\s*\w+\s+(?:0[xX])?([0-9A-Fa-f]+)\s*[,\s]\s*(?:0[xX])?([0-9A-Fa-f]+)', ...
                 'tokens', 'once');
    have = ~cellfun(@isempty, tok);
    if ~any(have)
        error('No "<cmd> reg,val" lines found in "%s".', txtFile);
    end

    regs = zeros(sum(have),1);  vals = zeros(sum(have),1);
    t = tok(have);
    for i = 1:numel(t)
        regs(i) = hex2dec(t{i}{1});
        vals(i) = hex2dec(t{i}{2});
    end

    % Tx block 0x29C/0x29D, Rx block 0x25C/0x25D
    for base = [hex2dec('29C'), hex2dec('25C')]
        aIdx = regs == base;
        dIdx = regs == base + 1;
        if any(aIdx) && any(dIdx), break; end
    end
    if ~any(aIdx) || ~any(dIdx)
        error('No PROGRAM_ADDR/PROGRAM_DATA pairs found in "%s".', txtFile);
    end

    % walk the file, remembering the last address written before each data byte
    slots = containers.Map('KeyType','double','ValueType','any');
    curAddr = NaN;
    for i = 1:numel(regs)
        if regs(i) == base
            curAddr = vals(i);
        elseif regs(i) == base + 1 && ~isnan(curAddr)
            p = floor(curAddr / 16);
            w = mod(curAddr, 16);
            if ~isKey(slots, p), slots(p) = nan(1,16); end
            tmp = slots(p);  tmp(w+1) = vals(i);  slots(p) = tmp;
        end
    end

    ks = sort(cell2mat(slots.keys));
    Bytes = zeros(numel(ks), 16);
    for i = 1:numel(ks)
        Bytes(i,:) = slots(ks(i));
    end
    if any(isnan(Bytes(:)))
        warning('%d setup word(s) were never written in "%s".', ...
            sum(isnan(Bytes(:))), txtFile);
        Bytes(isnan(Bytes)) = 0;
    end

    % frequencies, if the script carries them as comments
    fq = nan(numel(ks),1);
    fc = regexp(txt, 'profile\s+(\d+)\s*:\s*([\d.]+)\s*MHz', 'tokens', 'once');
    fc = fc(~cellfun(@isempty, fc));
    for i = 1:numel(fc)
        p = str2double(fc{i}{1});
        j = find(ks == p, 1);
        if ~isempty(j), fq(j) = str2double(fc{i}{2}); end
    end
end

function v = hexcol_to_num(col)
    if isnumeric(col)
        v = double(col(:));  return;
    end
    s = erase(erase(strtrim(string(col)), "0x"), "0X");
    v = zeros(numel(s),1);
    for i = 1:numel(s)
        if strlength(s(i)) == 0, v(i) = NaN;
        else,                    v(i) = hex2dec(char(s(i)));
        end
    end
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
