function R = pll_registers(freq, fref, opts)
%PLL_REGISTERS  AD9361 RFPLL synthesizer registers from the target frequency.
%
%   R = pll_registers(freq, fref)
%   R = pll_registers(freq, fref, opts)
%
%   freq  target RF frequency in MHz (scalar or vector)
%   fref  reference frequency in MHz
%
%   Returns a struct of column vectors:
%       R.R005        VCO divider code        -> 0x005[7:4]
%       R.R271        integer word [7:0]      -> 0x271[7:0]
%       R.R272        integer word [10:8]     -> 0x272[2:0]
%       R.R273        fractional word [7:0]   -> 0x273[7:0]
%       R.R274        fractional word [15:8]  -> 0x274[7:0]
%       R.R275        fractional word [22:16] -> 0x275[6:0]
%       R.divider     VCO divider code (same as R005)
%       R.f_vco       resulting VCO frequency in MHz
%       R.integer     synthesizer integer word
%       R.fractional  synthesizer fractional word
%
%   These six fields are computed exactly from the target frequency.  They are
%   NOT interpolated: the fractional word steps by roughly 400 000 counts per
%   1 MHz and wraps every ~20 MHz, so its low bytes are arithmetic wraparound
%   rather than a sampled smooth function, and no interpolation scheme of any
%   kind can reproduce them.
%
%   Faithful vectorised version of the reference implementation
%   pll_cac / GET_Synthesizer_REGISTER_VAL, with three additions:
%     - vectorised over freq;
%     - the fractional word is wrapped into the integer word if rounding
%       pushes it up to the modulus (otherwise 0x275 can overflow 7 bits);
%     - range checks on the integer and fractional words.

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'vcoMin'),  opts.vcoMin  = 6000;    end  % MHz
    if ~isfield(opts, 'modulus'), opts.modulus = 8388593; end  % AD9361 Fract-N
    if ~isfield(opts, 'warn'),    opts.warn    = true;    end

    freq = double(freq(:));
    if any(freq <= 0), error('pll_registers: frequencies must be positive.'); end
    if isempty(fref) || fref <= 0
        error('pll_registers: fref must be a positive reference frequency in MHz.');
    end

    % ---- VCO divider: double until the VCO clears vcoMin -----------------
    % Reference loop:  freq = freq*2;  while freq <= 6000, d = d+1; freq = freq*2; end
    f_vco = freq * 2;
    d     = zeros(size(freq));
    below = f_vco <= opts.vcoMin;
    guard = 0;
    while any(below)
        d(below)     = d(below) + 1;
        f_vco(below) = f_vco(below) * 2;
        below        = f_vco <= opts.vcoMin;
        guard = guard + 1;
        if guard > 64, error('pll_registers: VCO divider search did not converge.'); end
    end

    % ---- integer and fractional words -----------------------------------
    N   = f_vco / fref * 2;
    INT = floor(N);
    FRC = round(opts.modulus * (N - INT));

    % rounding can land exactly on the modulus; carry it into the integer word
    carry = FRC >= opts.modulus;
    INT(carry) = INT(carry) + 1;
    FRC(carry) = 0;

    % ---- range checks -----------------------------------------------------
    if opts.warn
        if any(INT > 2047 | INT < 0)
            warning(['pll_registers: integer word out of the 11-bit range ' ...
                     '(min %d, max %d). Check fref = %g MHz.'], ...
                     min(INT), max(INT), fref);
        end
        if any(FRC > 8388607 | FRC < 0)
            warning('pll_registers: fractional word out of the 23-bit range.');
        end
        if any(d > 15)
            warning('pll_registers: VCO divider code exceeds the 4-bit field.');
        end
    end

    % ---- byte split (Table 12 word order) ---------------------------------
    R272 = floor(INT / 256);
    R271 = INT - R272 * 256;

    R275 = floor(FRC / 65536);
    R274 = floor(FRC / 256) - R275 * 256;
    R273 = FRC - R274 * 256 - R275 * 65536;

    R = struct( ...
        'R005', d, 'R271', R271, 'R272', R272, ...
        'R273', R273, 'R274', R274, 'R275', R275, ...
        'divider', d, 'f_vco', f_vco, 'integer', INT, 'fractional', FRC);
end
