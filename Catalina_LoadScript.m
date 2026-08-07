%% ========================================================================
%  AD9361 Fast Lock loader - SPIWrite / WAIT script generator
%
%  Reads profiles.xlsx (from Catalina_MakeProfiles) and emits a flat command
%  script in the form:
%
%      SPIWrite 0x29D,0x2C
%      WAIT 10                      ; duration in milliseconds
%
%  Usage
%  -----
%     Catalina_LoadScript                                  % Tx, all rows
%     Catalina_LoadScript('profiles.xlsx','load.txt')
%     Catalina_LoadScript('profiles.xlsx','load.txt', 1:8)
%
%  Set cfg.path to 'tx', 'rx' or 'both'.
%
%  BOTH-PATH MODE
%  --------------
%  Rx and Tx have independent fast lock blocks and independent sets of 8
%  profile slots, so a pair of frequencies can be resident at once - one on
%  each path.  cfg.pairMode decides how the rows of the profile table are
%  paired:
%
%     'alternate'  (default)  rows 1,2 -> Tx slot 0, Rx slot 0
%                             rows 3,4 -> Tx slot 1, Rx slot 1  ...
%                             i.e. the first of each pair is Tx, the second Rx
%     'split'                 first half of the rows to Tx, second half to Rx
%     'files'                 Rx frequencies come from cfg.rxProfileFile,
%                             paired row by row with the Tx file
%
%  A CAVEAT ABOUT Rx PROFILES BUILT FROM Tx MEASUREMENTS
%  -----------------------------------------------------
%  The setup word FORMAT is the same for both paths: Table 11 (Rx) and
%  Table 12 (Tx) assign identical meanings to Program Address[3:0]; only the
%  source registers differ (0x231-0x251 for Rx, 0x271-0x291 for Tx).  So a
%  word packed from Tx measurements is structurally valid in the Rx block.
%
%  It is not necessarily CORRECT there.  Words 0-4 and word C - the integer
%  word, the fractional word and the VCO divider - come from the frequency by
%  formula and are identical for either synthesizer.  Everything else (VCO
%  varactor, bias, cal offset, ALC, charge pump, loop filter) is the result of
%  calibrating one specific synthesizer, and the Rx and Tx PLLs are separate
%  circuits that will not calibrate to the same values.  Reusing Tx values on
%  Rx will usually still lock, but the loop dynamics and the VCO operating
%  point will be off, and it may fail at band edges or over temperature.
%
%  For a correct Rx profile, measure the Rx registers and build a second model
%  with cfg.exactFields / cfg.formulaFields pointed at the Rx addresses, then
%  use cfg.pairMode = 'files'.  This script warns when it is asked to send
%  Tx-derived words to the Rx block.
%
%  REGISTERS  (UG-570 body text + ADI no-OS driver)
%  ------------------------------------------------
%      Rx   0x25A SETUP   0x25B INIT_DELAY   0x25C ADDR
%           0x25D DATA    0x25E READ         0x25F CTRL
%      Tx   the same block offset by 0x40:   0x29A .. 0x29F
%
%  THE FULL PROCEDURE THIS SCRIPT EMITS
%  ------------------------------------
%  STEP 1  Enable fast lock mode.  FAST_LOCK_SETUP with D0 set, D1 and D2
%          clear.  Stop the program clocks (CTRL = 0).
%  STEP 2  Mask the synthesizer-ready handshake in ENSM_CONFIG_2 (0x015).
%          THIS IS THE STEP THAT MAKES A RECALL ACTUALLY TAKE EFFECT.  The
%          ENSM normally waits for the synthesizer to report "ready" before
%          it will change state and apply a new frequency.  A fast lock recall
%          does not run a normal lock sequence, so that handshake never
%          completes: the ENSM never advances, the profile sits in memory
%          unused, and the part stays on the old frequency.  ADI's own
%          ad9361_fastlock_prepare() ends with exactly this write.
%  STEP 3  For each profile, for each of the 16 setup words:
%              SPIWrite ADDR,(profile<<4)|word
%              SPIWrite DATA,value
%              SPIWrite CTRL,0x03          ; WRITE | CLOCK_ENABLE
%          then once after the sixteenth word:
%              SPIWrite CTRL,0x00          ; stop clocks
%  STEP 4  Recall, bracketed by a forced ALERT: the synthesizer reprograms
%          itself on a state transition, so selecting a profile while the part
%          sits in FDD/TX/RX may not retune it until the next transition.
%              force ALERT -> write FAST_LOCK_SETUP -> release ALERT
%
%  No VCO calibration is needed anywhere in this sequence.  Words E and F of
%  every profile carry Force VCO Tune and Force ALC, and per ADI, if those are
%  written manually then Profile Init (D2) need not be set and no calibration
%  has to run.  That is the entire point of generating profiles offline.
%
%  FAST_LOCK_SETUP  (0x25A Rx / 0x29A Tx)
%  --------------------------------------
%      [D7:D5]  profile number 0..7            (cfg.profileShift)
%      [D2]     Fast Lock Profile Init  - set ONLY while creating a profile
%               with an on-chip VCO calibration; CLEAR to use saved profiles.
%               Per ADI: if the VCO Tune and ALC words (Program Address 0xE
%               and 0xF) are written manually, this bit need not be set and no
%               VCO calibration has to run.  This toolchain always writes them.
%      [D1]     Fast Lock Profile Pin Select - CTRL_IN0..CTRL_IN2 select the
%               profile instead of SPI.  Leave CLEAR for SPI control.
%      [D0]     Fast Lock Mode Enable   - set while creating profiles and in
%               normal operation while using them.
%
%  The profile number does NOT begin at D1.  An earlier version built the
%  recall word as (profile << 1) | 1.  For profile 1 that sets D1 and hands
%  selection to the CTRL_IN pins; for profile 2 it sets D2 and puts the part
%  into profile-creation mode.  Either way the profile you asked for is not
%  the one in use, and the symptom is a profile that appears not to load.
%
%  STILL UNVERIFIED
%  ----------------
%  cfg.ctrlStrobe (0x03), the stop-clocks placement, and the ENSM_CONFIG_2 bit
%  positions (cfg.rxSynthReadyBit, cfg.txSynthReadyBit, cfg.forceAlertBit) all
%  come from driver headers rather than UG-570.  Register 0x015 itself and the
%  existence of the ready-mask step are confirmed.  If a profile still does not
%  take, in order: check cfg.ensmConfig2Base holds whatever else your board
%  needs in 0x015, then the ready-mask bit positions, then try
%  cfg.stopClocksPerWord = true.
% ========================================================================

function Catalina_LoadScript(profileFile, outFile, rows)
    if nargin < 1 || isempty(profileFile), profileFile = 'profiles.xlsx';     end
    if nargin < 2 || isempty(outFile),     outFile     = 'fastlock_load.txt'; end

    cfg = load_config();

    [fq, Words] = read_profile_table(profileFile);
    if nargin < 3 || isempty(rows), rows = 1:size(Words,1); end
    rows = rows(:).';
    rows = rows(rows >= 1 & rows <= size(Words,1));
    if isempty(rows), error('No valid rows selected.'); end
    fq = fq(rows);  Words = Words(rows,:);

    fqRx = [];  WordsRx = [];
    if strcmpi(cfg.path,'both') && strcmpi(cfg.pairMode,'files')
        if isempty(cfg.rxProfileFile)
            error('pairMode ''files'' needs cfg.rxProfileFile.');
        end
        [fqRx, WordsRx] = read_profile_table(cfg.rxProfileFile);
        fprintf('Rx profiles from %s (%d row(s))\n', cfg.rxProfileFile, numel(fqRx));
    end

    A = build_assignments(fq, Words, fqRx, WordsRx, cfg);

    nTx = sum(strcmp({A.path}, 'tx'));
    nRx = sum(strcmp({A.path}, 'rx'));
    fprintf('Read %s : %d Tx profile(s), %d Rx profile(s)\n', ...
        profileFile, nTx, nRx);
    for a = 1:numel(A)
        fprintf('   %-2s slot %d  <-  %.6g MHz\n', upper(A(a).path), ...
            A(a).slot, A(a).freq);
    end

    if nRx > 0 && ~strcmpi(cfg.pairMode,'files')
        warning(['%d Rx profile(s) will be programmed with words packed from ' ...
                 'Tx measurements.\nWords 0-4 and C are frequency formulas and ' ...
                 'are fine; the VCO cal,\nvaractor, ALC, charge pump and loop ' ...
                 'filter values belong to the Tx\nsynthesizer. See the header ' ...
                 'note.'], nRx);
    end

    check_alc_words(A, cfg);
    write_script(outFile, A, cfg);
    if cfg.writeRecallFile
        write_recall(replace_ext(outFile, '_recall.txt'), A, cfg);
    end

    per = 16*3 + 1;
    fprintf('\n%d SPIWrite command(s) written to %s\n', numel(A)*per, outFile);
    print_preview(A, cfg);
end

%% ========================================================================
%  CONFIGURATION
% ========================================================================
function cfg = load_config()
    cfg.path      = 'tx';        % 'tx' | 'rx' | 'both'
    cfg.pairMode  = 'alternate'; % 'alternate' | 'split' | 'files'
    cfg.rxProfileFile = '';      % used only by pairMode 'files'
    cfg.numSlots  = 8;           % per path; Program Address[7:4] holds 0..7

    % ---- command syntax --------------------------------------------------
    cfg.cmd       = 'SPIWrite';
    cfg.sep       = ',';
    cfg.space     = ' ';
    cfg.radix     = 'hex';       % 'hex' | 'dec'
    cfg.prefix    = '0x';        % '' for bare hex
    cfg.regDigits = 3;
    cfg.valDigits = 2;
    cfg.comment   = ';';         % '' to suppress comments

    % ---- wait ------------------------------------------------------------
    % WAIT takes a duration in milliseconds: "WAIT 10" waits 10 ms.  Each
    % place that needs one has its own value, because they are not the same
    % kind of pause: a settling wait after a synthesizer transition has to be
    % long enough for the PLL, while a pause between register writes only has
    % to satisfy the SPI controller.  A value of 0 emits nothing.
    cfg.waitCmd      = 'WAIT';
    cfg.waitWordMs   = 0;   % between setup words - normally unnecessary
    cfg.waitProfileMs= 1;   % after a whole profile has been written
    cfg.waitSetupMs  = 1;   % after enabling fast lock / masking synth ready
    cfg.waitAlertMs  = 1;   % after forcing ALERT, before selecting a profile
    cfg.waitRecallMs = 2;   % after releasing ALERT - PLL settling

    % ---- programming sequence -------------------------------------------
    cfg.ctrlStrobe        = hex2dec('03');   % WRITE | CLOCK_ENABLE
    cfg.ctrlIdle          = hex2dec('00');   % stop clocks
    cfg.stopClocksPerWord = false;
    cfg.addrBeforeData    = true;

    % ---- FAST_LOCK_SETUP register (0x25A Rx / 0x29A Tx) ------------------
    % Bit map, from UG-570 and the ADI register description:
    %   D0        Fast Lock Mode Enable   - set when creating profiles AND in
    %                                       normal operation when using them
    %   D1        Fast Lock Profile Pin Select - CTRL_IN0..2 select the profile
    %                                       instead of SPI.  Leave CLEAR for
    %                                       SPI control.
    %   D2        Fast Lock Profile Init  - set only while CREATING a profile
    %                                       with an on-chip VCO calibration.
    %                                       Must be CLEAR to use saved profiles.
    %   profile   the 3-bit profile number, at cfg.profileShift
    %
    % The profile number does NOT start at D1.  An earlier version of this
    % file used (profile << 1), which for profile 1 sets D1 and hands profile
    % selection to the CTRL_IN pins, and for profile 2 sets D2 and puts the
    % part into profile-creation mode.  Either way the profile you asked for
    % is not the one that gets used - the symptom is a profile that appears
    % not to load at all.
    cfg.profileShift = 5;        % profile occupies [D7:D5]
    cfg.modeEnable   = hex2dec('01');   % D0
    cfg.pinSelect    = false;    % true sets D1 and hands control to CTRL_IN0..2
    cfg.profileInit  = false;    % D2 - see the note in write_script below

    cfg.emitSetupHeader = true;
    cfg.initDelayNs     = 0;     % >0 writes INIT_DELAY = ns/250 (wide BW only)
    cfg.emitRecallBlock = true;
    cfg.writeRecallFile = true;

    % ---- ENSM (this is what makes a recalled profile actually take) -------
    % ADI's own ad9361_fastlock_prepare() ends with:
    %     ad9361_spi_writef(REG_ENSM_CONFIG_2, ready_mask, 1);
    %     ad9361_trx_vco_cal_control(phy, tx, false);
    % REG_ENSM_CONFIG_2 is 0x015 (confirmed in ad9361.h).
    %
    % Why it matters: the ENSM normally waits for the synthesizer to report
    % "ready" before it will move state and apply a new frequency.  A fast
    % lock recall does not run a normal lock sequence, so that handshake never
    % completes, the ENSM never advances, and the profile sits in memory
    % unused - the part keeps transmitting on the old frequency.  Setting the
    % synth-ready MASK tells the ENSM to stop waiting.
    %
    % Bit positions below are from driver headers, not from UG-570. Verify.
    cfg.emitEnsmReady   = true;
    cfg.ensmConfig2Reg  = hex2dec('015');
    cfg.ensmConfig2Base = hex2dec('00');   % other bits of 0x015 on your board.
                                           % You have no SPI read, so the whole
                                           % byte is written; put whatever else
                                           % belongs in 0x015 here or it is lost.
    cfg.rxSynthReadyBit = hex2dec('08');   % D3
    cfg.txSynthReadyBit = hex2dec('04');   % D2

    % ---- ENSM state around a recall --------------------------------------
    % The synthesizer reprograms itself on a state transition.  Recalling a
    % profile while the part sits in FDD/TX/RX may not retune it until the
    % next transition, so the recall is bracketed: force ALERT, select the
    % profile, release ALERT.
    cfg.emitAlertBracket = true;
    cfg.forceAlertBit    = hex2dec('02');  % D1 of 0x015, FORCE_ALERT_STATE

    % ---- ALC word workaround ---------------------------------------------
    % ADI's ad9361_fastlock_recall() carries a workaround labelled "Lock
    % problem with same ALC word": if the profile being recalled has the same
    % word F as the one currently selected, the PLL does not re-lock.  The fix
    % is to rewrite word F as (value - 1) first, so the synthesizer always sees
    % it change.  Off by default because it costs four extra writes per recall;
    % turn it on if two profiles share an ALC word, which the generator warns
    % about.
    cfg.emitAlcNudge = false;

    % ---- SPIRead --------------------------------------------------------
    % With a read available, the profile memory can be verified instead of
    % assumed.  ADI's ad9361_fastlock_readval() reads a stored word with:
    %       write ADDR  = (profile<<4)|word
    %       write CTRL  = CLOCK_ENABLE | READ
    %       read  PROGRAM_READ                  (0x25E Rx / 0x29E Tx)
    %       write CTRL  = 0
    % The expected value is printed beside each read as a comment, so the log
    % can be compared line by line.
    cfg.readCmd      = 'SPIRead';
    cfg.readFmt      = '%s%s%s';   % cmd, space, register
    cfg.ctrlRead     = hex2dec('06');  % CLOCK_ENABLE | READ - VERIFY
    cfg.emitReadback = true;   % read all 16 words back after writing a profile
    cfg.emitStateReads = true; % read the ENSM state and SETUP around a recall
    cfg.stateReg     = hex2dec('017');  % REG_STATE, confirmed in ad9361.h

    % With a read available, 0x015 no longer has to be written blind.  The
    % script emits a read of it before the first write so you can put the
    % actual value into cfg.ensmConfig2Base once, instead of guessing which
    % other bits your board needs there.
    cfg.emitEnsmRead = true;
end

function R = register_map(path)
    R.SETUP      = hex2dec('25A');
    R.INIT_DELAY = hex2dec('25B');
    R.ADDR       = hex2dec('25C');
    R.DATA       = hex2dec('25D');
    R.READ       = hex2dec('25E');
    R.CTRL       = hex2dec('25F');

    switch lower(path)
        case 'tx'
            for f = string(fieldnames(R)).'
                R.(char(f)) = R.(char(f)) + hex2dec('040');
            end
        case 'rx'
            % base map is already Rx
        otherwise
            error('path must be ''tx'' or ''rx'' here.');
    end
end

%% ========================================================================
%  ASSIGNMENT OF ROWS TO PATHS AND SLOTS
% ========================================================================
function A = build_assignments(fq, Words, fqRx, WordsRx, cfg)
    A = struct('path',{}, 'R',{}, 'slot',{}, 'freq',{}, 'words',{});

    switch lower(cfg.path)
        case {'tx','rx'}
            n = numel(fq);
            if n > cfg.numSlots
                warning(['%d rows but only %d slots on the %s path. ' ...
                         'Keeping the first %d.'], n, cfg.numSlots, ...
                         upper(cfg.path), cfg.numSlots);
                n = cfg.numSlots;
            end
            for i = 1:n
                A = push(A, cfg.path, i-1, fq(i), Words(i,:));
            end

        case 'both'
            switch lower(cfg.pairMode)
                case 'alternate'
                    % rows 1,2 -> Tx slot 0, Rx slot 0;  rows 3,4 -> slot 1 ...
                    n = numel(fq);
                    if mod(n,2) == 1
                        warning(['Odd number of rows (%d): the last frequency ' ...
                                 '(%.6g MHz) has no partner and goes to Tx ' ...
                                 'only.'], n, fq(end));
                    end
                    np = ceil(n/2);
                    if np > cfg.numSlots
                        warning('%d pairs but only %d slots. Keeping the first %d.', ...
                            np, cfg.numSlots, cfg.numSlots);
                        np = cfg.numSlots;
                    end
                    for p = 1:np
                        iTx = 2*p - 1;
                        iRx = 2*p;
                        A = push(A, 'tx', p-1, fq(iTx), Words(iTx,:));
                        if iRx <= n
                            A = push(A, 'rx', p-1, fq(iRx), Words(iRx,:));
                        end
                    end

                case 'split'
                    n  = numel(fq);
                    h  = floor(n/2);
                    nT = min(h, cfg.numSlots);
                    nR = min(n - h, cfg.numSlots);
                    for i = 1:nT
                        A = push(A, 'tx', i-1, fq(i), Words(i,:));
                    end
                    for i = 1:nR
                        j = h + i;
                        A = push(A, 'rx', i-1, fq(j), Words(j,:));
                    end

                case 'files'
                    nT = min(numel(fq),   cfg.numSlots);
                    nR = min(numel(fqRx), cfg.numSlots);
                    for i = 1:nT
                        A = push(A, 'tx', i-1, fq(i), Words(i,:));
                    end
                    for i = 1:nR
                        A = push(A, 'rx', i-1, fqRx(i), WordsRx(i,:));
                    end

                otherwise
                    error('Unknown pairMode "%s".', cfg.pairMode);
            end

        otherwise
            error('cfg.path must be ''tx'', ''rx'' or ''both''.');
    end

    if isempty(A), error('No profiles to program.'); end
end

function A = push(A, path, slot, freq, words)
    n = numel(A) + 1;
    A(n).path  = lower(path);
    A(n).R     = register_map(path);
    A(n).slot  = slot;
    A(n).freq  = freq;
    A(n).words = words;
end

%% ========================================================================
%  SCRIPT GENERATION
% ========================================================================
function v = setup_word(slot, cfg, init)
%SETUP_WORD  Build the value for FAST_LOCK_SETUP (0x25A Rx / 0x29A Tx).
%
%   [D7:D5] profile number      (cfg.profileShift)
%   [D2]    Profile Init        - only while creating a profile on-chip
%   [D1]    Profile Pin Select  - CTRL_IN pins instead of SPI
%   [D0]    Fast Lock Mode Enable
    if nargin < 3, init = false; end
    v = bitshift(bitand(slot, 7), cfg.profileShift);
    v = bitor(v, cfg.modeEnable);
    if cfg.pinSelect, v = bitor(v, 2); end
    if init,          v = bitor(v, 4); end
end

function write_script(outFile, A, cfg)
    fid = fopen(outFile, 'w');
    if fid < 0, error('Cannot open %s for writing.', outFile); end

    banner(fid, cfg, A);

    % =====================================================================
    % STEP 1 - enable fast lock mode, and only fast lock mode
    % =====================================================================
    if cfg.emitSetupHeader
        cmt(fid, cfg, '=== STEP 1: enable fast lock mode ===');
        cmt(fid, cfg, 'D0 = mode enable, D1 = pin select (clear for SPI control),');
        cmt(fid, cfg, 'D2 = profile init (clear - we write VCO tune and ALC ourselves).');
        for p = unique({A.path}, 'stable')
            R = register_map(p{1});
            if cfg.initDelayNs > 0
                spiw(fid, cfg, R.INIT_DELAY, min(255, round(cfg.initDelayNs/250)), ...
                     sprintf('%s init delay %g ns / 250', upper(p{1}), cfg.initDelayNs));
            end
            spiw(fid, cfg, R.SETUP, setup_word(0, cfg, false), ...
                 sprintf('%s fast lock mode enable, profile 0', upper(p{1})));
            spiw(fid, cfg, R.CTRL, cfg.ctrlIdle, ...
                 sprintf('%s program clocks stopped', upper(p{1})));
        end
        waitc(fid, cfg, cfg.waitSetupMs, 'let fast lock mode settle');
        fprintf(fid, '\n');
    end

    % =====================================================================
    % STEP 2 - let the ENSM proceed without a synthesizer-ready handshake
    % =====================================================================
    if cfg.emitEnsmReady
        cmt(fid, cfg, '=== STEP 2: mask the synthesizer-ready handshake ===');
        cmt(fid, cfg, 'A fast lock recall does not run a normal lock sequence, so the');
        cmt(fid, cfg, 'synth never reports "ready". Without this the ENSM waits forever,');
        cmt(fid, cfg, 'never applies the profile, and the part stays on the old frequency.');
        if cfg.emitEnsmRead
            spir(fid, cfg, cfg.ensmConfig2Reg, ...
                 'read 0x015 FIRST, then put that value in cfg.ensmConfig2Base');
            cmt(fid, cfg, 'The write below replaces the whole byte. Until the read value is');
            cmt(fid, cfg, 'in cfg.ensmConfig2Base, any other bits 0x015 held are lost.');
        else
            cmt(fid, cfg, 'The whole byte is written - put any other 0x015 bits your board');
            cmt(fid, cfg, 'needs into cfg.ensmConfig2Base, or they are lost.');
        end
        spiw(fid, cfg, cfg.ensmConfig2Reg, ensm_ready_value(A, cfg), ...
             'ENSM config 2: synth ready mask set');
        if cfg.emitStateReads
            spir(fid, cfg, cfg.ensmConfig2Reg, sprintf('expect %s', ...
                 fmt(cfg, ensm_ready_value(A, cfg), cfg.valDigits)));
        end
        waitc(fid, cfg, cfg.waitSetupMs, 'synth ready mask applied');
        fprintf(fid, '\n');
    end

    % =====================================================================
    % STEP 3 - write the 16 setup words of each profile
    % =====================================================================
    cmt(fid, cfg, '=== STEP 3: write the profile setup words ===');
    fprintf(fid, '\n');

    for a = 1:numel(A)
        R = A(a).R;
        cmt(fid, cfg, sprintf('---- %s profile %d : %.6g MHz ----', ...
            upper(A(a).path), A(a).slot, A(a).freq));

        for k = 1:16
            addr = bitor(bitshift(A(a).slot, 4), k-1);

            if cfg.addrBeforeData
                spiw(fid, cfg, R.ADDR, addr,          sprintf('word %X address', k-1));
                spiw(fid, cfg, R.DATA, A(a).words(k), sprintf('word %X data', k-1));
            else
                spiw(fid, cfg, R.DATA, A(a).words(k), sprintf('word %X data', k-1));
                spiw(fid, cfg, R.ADDR, addr,          sprintf('word %X address', k-1));
            end

            spiw(fid, cfg, R.CTRL, cfg.ctrlStrobe, 'write strobe');

            if cfg.stopClocksPerWord
                spiw(fid, cfg, R.CTRL, cfg.ctrlIdle, 'stop clocks');
            end
            waitc(fid, cfg, cfg.waitWordMs);
        end

        if ~cfg.stopClocksPerWord
            spiw(fid, cfg, R.CTRL, cfg.ctrlIdle, 'stop clocks');
        end
        waitc(fid, cfg, cfg.waitProfileMs, 'profile written');
        fprintf(fid, '\n');

        if cfg.emitReadback
            emit_readback(fid, cfg, A, a);
        end
    end

    % =====================================================================
    % STEP 4 - recall
    % =====================================================================
    if cfg.emitRecallBlock
        cmt(fid, cfg, '=== STEP 4: recall - keep the block you need ===');
        cmt(fid, cfg, sprintf('profile number sits at [D%d:D%d]; D0 stays set.', ...
            cfg.profileShift+2, cfg.profileShift));
        if cfg.emitAlertBracket
            cmt(fid, cfg, 'The synthesizer reloads on a state transition, so the');
            cmt(fid, cfg, 'selection is bracketed by a forced ALERT.');
        end
        fprintf(fid, '\n');

        for s = unique([A.slot])
            sel = find([A.slot] == s);
            lbl = strjoin(arrayfun(@(a) sprintf('%s %.6g MHz', ...
                    upper(A(a).path), A(a).freq), sel, ...
                    'UniformOutput', false), '  +  ');
            cmt(fid, cfg, sprintf('-- slot %d : %s --', s, lbl));
            emit_recall_block(fid, cfg, A, sel);
            fprintf(fid, '\n');
        end
    end

    fclose(fid);
end

function check_alc_words(A, cfg)
%CHECK_ALC_WORDS  Warn when two profiles share a Force ALC word.
%
% ADI's own ad9361_fastlock_recall() carries a workaround labelled
% "Lock problem with same ALC word": if the profile being recalled has the
% same word F as the profile currently selected, the PLL does not re-lock and
% the part keeps running on the old frequency.  The driver's fix is to write
% the new profile's word F as (value - 1) before selecting it, which forces a
% change the synthesizer notices.
%
% Word F is (ALC << 1) | VCO_tune[8], so two nearby frequencies with the same
% ALC result and the same tune MSB collide.  That is not rare: ALC is a coarse
% amplitude setting that changes slowly across a band.
    for p = unique({A.path})
        sel = find(strcmp({A.path}, p{1}));
        if numel(sel) < 2, continue; end

        wF = arrayfun(@(a) A(a).words(16), sel);
        [u, ~, ic] = unique(wF);
        for k = 1:numel(u)
            g = sel(ic == k);
            if numel(g) < 2, continue; end
            warning(['%s profiles %s share Force ALC word 0x%02X ' ...
                     '(%s MHz).\nRecalling one straight after another may not ' ...
                     're-lock: the synthesizer sees no\nchange in word F. Set ' ...
                     'cfg.emitAlcNudge = true to emit ADI''s workaround, or\n' ...
                     'avoid recalling these two consecutively.'], ...
                     upper(p{1}), mat2str([A(g).slot]), u(k), ...
                     strjoin(compose('%.6g', [A(g).freq]), ', '));
        end
    end
end

function emit_alc_nudge(fid, cfg, A, a)
% ADI's workaround: rewrite word F as (value - 1) before selecting the
% profile, so the synthesizer always sees the ALC word change.
    R = A(a).R;
    v = max(A(a).words(16) - 1, 0);
    cmt(fid, cfg, 'ALC nudge: force word F to change so the PLL re-locks');
    addr = bitor(bitshift(A(a).slot, 4), 15);
    spiw(fid, cfg, R.ADDR, addr, 'word F address');
    spiw(fid, cfg, R.DATA, v,    sprintf('word F data (0x%02X - 1)', A(a).words(16)));
    spiw(fid, cfg, R.CTRL, cfg.ctrlStrobe, 'write strobe');
    spiw(fid, cfg, R.CTRL, cfg.ctrlIdle,   'stop clocks');
end

function spir(fid, cfg, reg, note)
%SPIR  Emit a read command.
    s = sprintf(cfg.readFmt, cfg.readCmd, cfg.space, fmt(cfg, reg, cfg.regDigits));
    if ~isempty(cfg.comment) && nargin >= 4 && ~isempty(note)
        fprintf(fid, '%-30s %s %s\n', s, cfg.comment, note);
    else
        fprintf(fid, '%s\n', s);
    end
end

function emit_readback(fid, cfg, A, a)
%EMIT_READBACK  Read all 16 words of one profile out of the part again.
%
% This is the check that separates "the profile never reached the memory" from
% "the profile is in the memory but the part is not applying it".  Those two
% have completely different fixes, and without a read there is no way to tell
% them apart.  The expected value is printed beside each read, so the log can
% be compared line by line.
    R = A(a).R;
    cmt(fid, cfg, sprintf('---- verify %s profile %d (%.6g MHz) ----', ...
        upper(A(a).path), A(a).slot, A(a).freq));

    for k = 1:16
        addr = bitor(bitshift(A(a).slot, 4), k-1);
        spiw(fid, cfg, R.ADDR, addr,        sprintf('word %X address', k-1));
        spiw(fid, cfg, R.CTRL, cfg.ctrlRead, 'clock enable | read');
        spir(fid, cfg, R.READ, sprintf('word %X: expect %s', k-1, ...
             fmt(cfg, A(a).words(k), cfg.valDigits)));
        spiw(fid, cfg, R.CTRL, cfg.ctrlIdle, 'stop clocks');
    end
    fprintf(fid, '\n');
end

function emit_recall_block(fid, cfg, A, sel)
% Force ALERT, select the profile on every path in this slot, release ALERT.
    if cfg.emitStateReads
        spir(fid, cfg, cfg.stateReg, 'ENSM state before');
    end

    if cfg.emitAlertBracket
        spiw(fid, cfg, cfg.ensmConfig2Reg, ...
             bitor(ensm_ready_value(A, cfg), cfg.forceAlertBit), 'force ALERT');
        waitc(fid, cfg, cfg.waitAlertMs, 'reach ALERT before selecting');
        if cfg.emitStateReads
            spir(fid, cfg, cfg.stateReg, 'ENSM state - expect ALERT');
        end
    end

    for a = sel
        if cfg.emitAlcNudge, emit_alc_nudge(fid, cfg, A, a); end
        spiw(fid, cfg, A(a).R.SETUP, setup_word(A(a).slot, cfg, false), ...
             sprintf('%s select profile %d', upper(A(a).path), A(a).slot));
        if cfg.emitStateReads
            spir(fid, cfg, A(a).R.SETUP, sprintf('%s SETUP: expect %s', ...
                 upper(A(a).path), ...
                 fmt(cfg, setup_word(A(a).slot, cfg, false), cfg.valDigits)));
        end
    end
    waitc(fid, cfg, cfg.waitAlertMs, 'profile selected');

    if cfg.emitAlertBracket
        spiw(fid, cfg, cfg.ensmConfig2Reg, ensm_ready_value(A, cfg), ...
             'release ALERT - the profile takes effect here');
        waitc(fid, cfg, cfg.waitRecallMs, 'PLL settling');
        if cfg.emitStateReads
            spir(fid, cfg, cfg.stateReg, 'ENSM state after - expect TX/RX/FDD');
        end
    end
end

function v = ensm_ready_value(A, cfg)
% ENSM config 2 with the synth-ready mask set for whichever paths are in use.
    v = cfg.ensmConfig2Base;
    p = unique({A.path});
    if any(strcmp(p, 'rx')), v = bitor(v, cfg.rxSynthReadyBit); end
    if any(strcmp(p, 'tx')), v = bitor(v, cfg.txSynthReadyBit); end
end

function write_recall(outFile, A, cfg)
% One block per slot.  In both-path mode the Tx and Rx halves of a pair share
% a slot number, so the pair is emitted together and can be recalled as a unit.
    fid = fopen(outFile, 'w');
    if fid < 0, warning('Cannot open %s', outFile); return; end

    banner(fid, cfg, A);
    cmt(fid, cfg, 'Recall snippets. Each block selects one profile (or pair).');
    fprintf(fid, '\n');

    for s = unique([A.slot])
        sel = find([A.slot] == s);
        lbl = strjoin(arrayfun(@(a) sprintf('%s %.6g MHz', ...
                upper(A(a).path), A(a).freq), sel, 'UniformOutput', false), '  +  ');
        cmt(fid, cfg, sprintf('---- slot %d : %s ----', s, lbl));
        emit_recall_block(fid, cfg, A, sel);
        fprintf(fid, '\n');
    end

    fclose(fid);
    fprintf('Recall snippets written to %s\n', outFile);
end

%% ========================================================================
%  COMMAND FORMATTING
% ========================================================================
function spiw(fid, cfg, reg, val, note)
    s = sprintf('%s%s%s%s%s', cfg.cmd, cfg.space, ...
        fmt(cfg, reg, cfg.regDigits), cfg.sep, fmt(cfg, val, cfg.valDigits));
    if ~isempty(cfg.comment) && nargin >= 5 && ~isempty(note)
        fprintf(fid, '%-30s %s %s\n', s, cfg.comment, note);
    else
        fprintf(fid, '%s\n', s);
    end
end

function waitc(fid, cfg, ms, note)
%WAITC  Emit "WAIT <ms>".  A duration of 0 emits nothing.
    if nargin < 3 || isempty(ms) || ms <= 0, return; end
    s = sprintf('%s%s%g', cfg.waitCmd, cfg.space, ms);
    if ~isempty(cfg.comment) && nargin >= 4 && ~isempty(note)
        fprintf(fid, '%-30s %s %s\n', s, cfg.comment, note);
    else
        fprintf(fid, '%s\n', s);
    end
end

function s = fmt(cfg, v, digits)
    if strcmpi(cfg.radix, 'dec')
        s = sprintf('%d', v);
    else
        s = sprintf('%s%s', cfg.prefix, dec2hex(v, digits));
    end
end

function cmt(fid, cfg, txt)
    if ~isempty(cfg.comment)
        fprintf(fid, '%s %s\n', cfg.comment, txt);
    end
end

function banner(fid, cfg, A)
    cmt(fid, cfg, sprintf('AD9361 Fast Lock - generated %s', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS')));                     %#ok<TNOW1,DATST>
    for p = unique({A.path}, 'stable')
        R = register_map(p{1});
        cmt(fid, cfg, sprintf('%s block: SETUP 0x%s  ADDR 0x%s  DATA 0x%s  CTRL 0x%s', ...
            upper(p{1}), dec2hex(R.SETUP,3), dec2hex(R.ADDR,3), ...
            dec2hex(R.DATA,3), dec2hex(R.CTRL,3)));
    end
    cmt(fid, cfg, sprintf('%d profile(s), %.6g - %.6g MHz', ...
        numel(A), min([A.freq]), max([A.freq])));
    fprintf(fid, '\n');
end

%% ========================================================================
%  INPUT
% ========================================================================
function [fq, Words] = read_profile_table(profileFile)
% Accepts the output of Catalina_MakeProfiles: a frequency column plus either
% W0..WF (hex text with or without a 0x prefix, or decimal) or a single
% 32-character Profile_Hex32 column.
    if ~isfile(profileFile)
        error('Profile file "%s" not found.', profileFile);
    end

    T  = readtable(profileFile, 'VariableNamingRule','preserve', ...
                   'TextType','string');
    vn = string(T.Properties.VariableNames);

    fcol = find(contains(lower(vn), 'freq'), 1);
    if isempty(fcol), error('No frequency column in "%s".', profileFile); end
    fq = T{:, fcol};
    if ~isnumeric(fq), fq = str2double(string(fq)); end

    wcol  = arrayfun(@(k) find(vn == "W" + string(dec2hex(k,1)), 1), 0:15, ...
                     'UniformOutput', false);
    haveW = all(~cellfun(@isempty, wcol));

    if haveW
        Words = zeros(height(T), 16);
        for k = 1:16
            Words(:,k) = hexcol_to_num(T{:, wcol{k}});
        end
    else
        bcol = find(contains(lower(vn),'hex32') | contains(lower(vn),'blob'), 1);
        if isempty(bcol)
            error('"%s" has neither W0..WF nor a Profile_Hex32 column.', profileFile);
        end
        blob  = string(T{:, bcol});
        Words = zeros(numel(blob), 16);
        for r = 1:numel(blob)
            s = char(strtrim(blob(r)));
            if numel(s) ~= 32
                error('Row %d: expected 32 hex characters, found %d.', r, numel(s));
            end
            Words(r,:) = hex2dec(reshape(s, 2, 16).');
        end
    end

    keep  = ~isnan(fq) & ~any(isnan(Words), 2);
    fq    = fq(keep);
    Words = Words(keep, :);

    if isempty(fq), error('No usable rows in "%s".', profileFile); end
    if any(Words(:) < 0 | Words(:) > 255)
        error('A setup word is outside the 0..255 byte range.');
    end
end

function v = hexcol_to_num(col)
% A column may arrive as text ("0x2C" or "2C") or, if Excel decided a column
% of all-numeric-looking hex was a number, as a double.  Handle both.
    if isnumeric(col)
        v = double(col(:));
        return;
    end
    s = strtrim(string(col));
    s = erase(erase(s, "0x"), "0X");
    v = zeros(numel(s), 1);
    for i = 1:numel(s)
        if strlength(s(i)) == 0
            v(i) = NaN;
        else
            v(i) = hex2dec(char(s(i)));
        end
    end
end

%% ========================================================================
%  PREVIEW
% ========================================================================
function print_preview(A, cfg)
    R = A(1).R;
    fprintf('\nFirst commands (%s slot %d, %.6g MHz):\n\n', ...
        upper(A(1).path), A(1).slot, A(1).freq);
    addr = bitor(bitshift(A(1).slot,4), 0);
    fprintf('  %s%s%s%s%s\n', cfg.cmd, cfg.space, fmt(cfg,R.ADDR,3), cfg.sep, ...
        fmt(cfg,addr,2));
    fprintf('  %s%s%s%s%s\n', cfg.cmd, cfg.space, fmt(cfg,R.DATA,3), cfg.sep, ...
        fmt(cfg,A(1).words(1),2));
    fprintf('  %s%s%s%s%s\n', cfg.cmd, cfg.space, fmt(cfg,R.CTRL,3), cfg.sep, ...
        fmt(cfg,cfg.ctrlStrobe,2));
    fprintf('  ... 15 more words, then stop clocks\n');

    fprintf('\nBytes: ');  fprintf('%02X ', A(1).words);  fprintf('\n');

    fprintf(['\nBefore hardware use, confirm cfg.ctrlStrobe (0x%02X) and the ' ...
             'stop-clocks\nplacement against your driver. Neither is stated ' ...
             'in UG-570.\n'], cfg.ctrlStrobe);
end

function s = replace_ext(f, newExt)
    [p, n, ~] = fileparts(f);
    s = fullfile(p, [n newExt]);
end
