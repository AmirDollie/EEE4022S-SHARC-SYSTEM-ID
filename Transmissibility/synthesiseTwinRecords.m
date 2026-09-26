function varargout = synthesiseTwinRecords(twin, varargin)
%SYNTHESISETWINRECORDS Stochastic multi-sensor twin records from the EMM,
%with a cached node FRF so Monte Carlo never re-solves the Forward Model.
%
%   TWIN = SYNTHESISETWINRECORDS(SPEC)                     % prepare only
%   TWIN = SYNTHESISETWINRECORDS(SPEC, 'CacheDir', D, 'Verbose', tf)
%   [Y, T, INFO, TWIN] = SYNTHESISETWINRECORDS(TWIN_OR_SPEC, N, DT, Name, Value)
%
%   PURPOSE AND SCOPE
%   Generates the DATA side of the output-only twin: Y(t) from one common
%   incident-wave input through the EMM FRF, plus independent sensor noise.
%   It never defines the inversion feature grid. Two grids are kept apart:
%     node grid  : where the EMM is solved (cached, 'nodeSpacing' apart);
%     fine grid  : the FFT bins of the record, omega_k = k 2 pi / (N DT),
%                  H interpolated from the nodes; used ONLY for synthesis;
%     Welch grid : set by extractTransmissibility (S2) and evaluated on the
%                  model side by transmissibilityFeatureVector (S1).
%   The fine grid changes with the record length, so the expensive part
%   that is cached is the node FRF, not H on a particular fine grid.
%
%   SPEC (struct; a prepared TWIN has the same fields plus the results)
%     p            [beta gamma R], full parameter vector (required for EMM)
%     sensors      nS x 2 [r theta], non-dimensional r (required)
%     band         [wLo wHi] rad/s synthesis band. Default: the validated
%                  EMM band, alpha in ValidatedAlphaRange mapped through
%                  omega = sqrt(alpha g / depth) (about [2.978 8.501]).
%     nodeSpacing  target EMM node spacing (rad/s), default 0.1 (57 EMM
%                  solves over the default band, about 1 min in MATLAB).
%                  Spline error scales as spacing^4: 0.25 already gives
%                  about 3e-5 relative error in T on the EMM (tester group 6)
%     frfOptions   cell of name-value pairs for computeSensorFRF
%                  (e.g. {'Truncation', [60 12 12]}); default {}
%     checkPoints  number of node midpoints at which the interpolated FRF
%                  is compared with a direct solve, default 3 (first, middle
%                  and last gap); 0 disables the check
%     checkReference  sensor used as reference for the T check (default nS)
%     frf          OPTIONAL function handle @(omega) -> nS x nW ACCELERATION
%                  FRF. Replaces the EMM (analytical tests). Not cacheable.
%
%   NODE FRF AND INTERPOLATION
%   Nodes: linspace(band(1), band(2), ceil(diff(band)/nodeSpacing) + 1), so
%   both band edges are nodes and the fine grid never extrapolates. The
%   displacement-type FRF H_a / (-omega^2) (smoother than H_a) is
%   interpolated with a cubic spline in its real and imaginary parts, then
%   multiplied back by -omega^2. twin.interpCheck reports the relative
%   error of H and of T = H_j / H_ref at the check midpoints (worst case,
%   halfway between nodes). twin.Hacc(omega) evaluates the interpolated
%   acceleration FRF anywhere inside the band: this is the exact truth the
%   synthesised records follow (S4 bias should be measured against it,
%   and against the S1 EMM features to include the interpolation error).
%
%   DISK CACHE ('CacheDir')
%   Each cache file stores key = {source, p, sensors, band, nodeSpacing,
%   frfOptions, checkPoints, checkReference} together with omegaNodes, the
%   node FRFs and interpCheck. checkPoints / checkReference do not change
%   the FRF but do change the stored interpCheck, so they are in the key
%   (a different diagnostic setting writes a new file rather than
%   returning a stale check). On a request
%   every file's key is loaded and compared (p and sensors to relative
%   1e-12, the rest exactly); only a matching file is reused, otherwise
%   the nodes are solved and a new file is written. The sampling period is
%   not part of the key: the node FRF does not depend on it (Nyquist is
%   checked at synthesis). twin.cacheStatus is 'computed', 'loaded' or
%   'memory' (a prepared TWIN passed back in).
%
%   RECORD MODEL (per realisation, noise-free part via stochasticSensorSynthesis)
%     y(t) = Re sum_k w_k sqrt(S_u(omega_k) dOmega) H(omega_k) xi_k e^{i omega_k t}
%   with xi_k complex Gaussian and SHARED by all sensors (single input).
%   H is the ACCELERATION FRF per metre of incident amplitude, so S_u is
%   the incident wave elevation spectrum (m^2 per rad/s, one-sided) and Y
%   is in m/s^2. The record is exactly periodic with period N*DT
%   (frequency-domain synthesis); Welch segments shorter than the record
%   see a continuous spectrum. The amplitude taper w is 1 inside the band
%   and falls to 0 at the band edges with raised-cosine ramps of width
%   'TaperWidth' INSIDE the band, so the unvalidated EMM range is never
%   excited. Welch bins inside the ramps receive less input: expect
%   larger variance there when noise is present.
%   Sensor noise: independent white Gaussian per channel, one-sided PSD
%   sigma^2 DT / pi per rad/s (the extractTransmissibility convention, so
%   info.noisePSD(ref) can be passed as its 'NoisePSD'). Noise is drawn from
%   its own stream ('NoiseSeed'), so the wave realisation for a given Seed
%   does not change when the noise level changes (common random numbers).
%
%   RECORD OPTIONS
%     'Seed'            [] (global stream) or nonnegative integer (local stream)
%     'NoiseSeed'       default Seed + 1e6 ([] if Seed is [])
%     'NumRealisations' default 1
%     'Spectrum'        'jonswap' (default), 'white', or a handle @(omega)
%                       returning S_u (m^2 s/rad) at the given frequencies
%     'Hs'              JONSWAP significant wave height (m), default 0.05.
%                       Normalised over the FULL spectrum (4 sqrt(m0) = Hs),
%                       before band restriction; info.HsInBand is what is
%                       actually synthesised.
%     'PeakFrequency'   JONSWAP omega_p (rad/s), default 5.0
%     'PeakEnhancement' JONSWAP gamma_J, default 3.3 (sigma 0.07 / 0.09)
%     'WhiteLevel'      flat S_u for 'white' (m^2 s/rad), default 1e-4
%     'TaperWidth'      ramp width at each band edge (rad/s), default 0.2;
%                       0 gives a hard edge
%     'NoiseStd'        absolute sensor noise std (m/s^2), scalar or 1 x nS,
%                       default 0
%     'NoiseFraction'   noise std as a fraction of each channel's expected
%                       signal std, scalar or 1 x nS, default 0. Use either
%                       NoiseStd or NoiseFraction, not both.
%     'CacheDir', 'Verbose'   as in prepare mode
%
%   OUTPUTS
%     Y     nS x N x nR records (m/s^2), including noise
%     T     1 x N sample times
%     INFO  struct: omega (fine grid), H (fine-grid acceleration FRF),
%           dOmega, N, dt, recordLength, inputPSD (untapered S_u on the fine
%           grid), taper, spectrum (handle for S_u), spectrumType, Hs,
%           HsInBand, expectedPSD / expectedCrossSpectrum /
%           expectedCovariance (signal only, from stochasticSensorSynthesis),
%           signalStd (1 x nS expected), noiseStd (1 x nS), noisePSD (1 x nS),
%           seed, noiseSeed, numRealisations, cacheStatus, interpCheck
%     TWIN  prepared twin (pass it back in to skip preparation)
%
%   ERROR IDENTIFIERS (synthesiseTwinRecords:<id>)
%     badSpec, badOption, badN, badDt, aboveNyquist, tooFewBins, badSpectrum,
%     badNoise, badFRF, cacheNotAllowed
%   Errors from computeSensorFRF / stochasticSensorSynthesis propagate.
%
%   Lives in Transmissibility/. Needs SSI/ (stochasticSensorSynthesis) and
%   Inverse Mapping/ (computeSensorFRF). Tester:
%   Transmissibility/Unit Tests Transmissibility/synthesiseTwinRecordsTester.m

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'SSI'));
addpath(fullfile(thisDir, '..', 'Inverse Mapping'));

%% ---- Mode ------------------------------------------------------------------
prepareOnly = isempty(varargin) || ischar(varargin{1}) || isa(varargin{1}, 'string');
if prepareOnly
    opts = varargin;
else
    if numel(varargin) < 2
        error('synthesiseTwinRecords:badOption', 'Record mode needs N and DT.');
    end
    N = varargin{1}; dt = varargin{2}; opts = varargin(3:end);
end

opt = struct('CacheDir', '', 'Verbose', false, 'Seed', [], 'NoiseSeed', NaN, 'NumRealisations', 1, ...
    'Spectrum', 'jonswap', 'Hs', 0.05, 'PeakFrequency', 5.0, 'PeakEnhancement', 3.3, ...
    'WhiteLevel', 1e-4, 'TaperWidth', 0.2, 'NoiseStd', 0, 'NoiseFraction', 0);
prepNames = {'CacheDir', 'Verbose'};
if mod(numel(opts), 2) ~= 0
    error('synthesiseTwinRecords:badOption', 'Options must be name-value pairs.');
end
names = fieldnames(opt);
for a = 1:2:numel(opts)
    key = opts{a};
    if isa(key, 'string'), key = char(key); end
    hit = ischar(key) & strcmpi(key, names);
    if ~any(hit)
        error('synthesiseTwinRecords:badOption', 'Unknown option.');
    end
    if prepareOnly && ~any(strcmpi(names{hit}, prepNames))
        error('synthesiseTwinRecords:badOption', ...
            'Option ''%s'' only applies in record mode (pass N and DT).', names{hit});
    end
    opt.(names{hit}) = opts{a + 1};
end
cacheDir = opt.CacheDir;
if isa(cacheDir, 'string'), cacheDir = char(cacheDir); end
if ~ischar(cacheDir)
    error('synthesiseTwinRecords:badOption', 'CacheDir must be a folder name.');
end
if ~((islogical(opt.Verbose) || isnumeric(opt.Verbose)) && isscalar(opt.Verbose))
    error('synthesiseTwinRecords:badOption', 'Verbose must be true or false.');
end

%% ---- Prepare (or reuse) the node FRF ---------------------------------------------
twin = prepareTwin(twin, cacheDir, logical(opt.Verbose));
if prepareOnly
    varargout = {twin};
    return
end

%% ---- Record inputs ----------------------------------------------------------------
if ~(isnumeric(N) && isscalar(N) && isfinite(N) && N == round(N) && N >= 16)
    error('synthesiseTwinRecords:badN', 'N must be an integer >= 16.');
end
if ~(isnumeric(dt) && isscalar(dt) && isreal(dt) && isfinite(dt) && dt > 0)
    error('synthesiseTwinRecords:badDt', 'DT must be a positive finite scalar.');
end
if twin.band(2) >= pi / dt
    error('synthesiseTwinRecords:aboveNyquist', 'Band top %.4g rad/s is not below Nyquist %.4g rad/s.', ...
        twin.band(2), pi / dt);
end
nS = size(twin.sensors, 1);
seed = opt.Seed;
if ~isempty(seed) && ~(isnumeric(seed) && isscalar(seed) && isfinite(seed) && seed >= 0 && seed == round(seed))
    error('synthesiseTwinRecords:badOption', 'Seed must be [] or a nonnegative integer.');
end
nseed = opt.NoiseSeed;
if isnumeric(nseed) && isscalar(nseed) && isnan(nseed)
    if isempty(seed), nseed = []; else, nseed = seed + 1e6; end
end
if ~isempty(nseed) && ~(isnumeric(nseed) && isscalar(nseed) && isfinite(nseed) && nseed >= 0 && nseed == round(nseed))
    error('synthesiseTwinRecords:badOption', 'NoiseSeed must be [] or a nonnegative integer.');
end
nR = opt.NumRealisations;
if ~(isnumeric(nR) && isscalar(nR) && isfinite(nR) && nR >= 1 && nR == round(nR))
    error('synthesiseTwinRecords:badOption', 'NumRealisations must be a positive integer.');
end
tw = opt.TaperWidth;
if ~(isnumeric(tw) && isscalar(tw) && isfinite(tw) && tw >= 0 && 2 * tw < diff(twin.band))
    error('synthesiseTwinRecords:badOption', 'TaperWidth must be >= 0 and less than half the band width.');
end
nStd = checkNoise(opt.NoiseStd, nS, 'NoiseStd');
nFrac = checkNoise(opt.NoiseFraction, nS, 'NoiseFraction');
if any(nStd > 0) && any(nFrac > 0)
    error('synthesiseTwinRecords:badNoise', 'Give NoiseStd or NoiseFraction, not both.');
end

%% ---- Fine grid, FRF, spectrum, taper ---------------------------------------------
dW = 2 * pi / (N * dt);
k = ceil(twin.band(1) / dW - 1e-9):floor(twin.band(2) / dW + 1e-9);
k = k(k * dW >= twin.band(1) & k * dW <= twin.band(2));
if numel(k) < 2
    error('synthesiseTwinRecords:tooFewBins', ...
        'Only %d FFT bin(s) of spacing %.4g rad/s fall in the band; lengthen the record.', numel(k), dW);
end
w = k * dW;
H = twin.Hacc(w);
[Sfun, spec] = buildSpectrum(opt);
Su = Sfun(w);
if ~(isreal(Su) && isequal(size(Su), size(w)) && all(isfinite(Su)) && all(Su >= 0))
    error('synthesiseTwinRecords:badSpectrum', 'Spectrum must return real, finite, nonnegative values, one per frequency.');
end
taper = edgeTaper(w, twin.band, tw);

%% ---- Signal (shared input) and noise (own stream) -------------------------------
[Y, t, sinfo] = stochasticSensorSynthesis(w, H, dt, 'Seed', seed, 'NumRealisations', nR, ...
    'InputPSD', Su, 'Taper', taper);
signalStd = sqrt(diag(sinfo.expectedCovariance)).';
if any(nFrac > 0)
    nStd = nFrac .* signalStd;
end
if any(nStd > 0)
    Y = Y + nStd(:) .* drawStandardNormal([nS, N, nR], nseed);
end

%% ---- Info -------------------------------------------------------------------------------
info = struct();
info.omega = w;
info.binIndex = k;
info.H = H;
info.dOmega = dW;
info.N = N;
info.dt = dt;
info.recordLength = N * dt;
info.band = twin.band;
info.inputPSD = Su;
info.taper = taper;
info.spectrum = Sfun;
info.spectrumType = spec.type;
info.spectrumParameters = spec;
info.Hs = spec.Hs;
info.HsInBand = 4 * sqrt(sum(Su .* taper.^2) * dW);
info.expectedPSD = sinfo.expectedPSD;
info.expectedCrossSpectrum = sinfo.expectedCrossSpectrum;
info.expectedCovariance = sinfo.expectedCovariance;
info.signalStd = signalStd;
info.noiseStd = nStd;
info.noisePSD = nStd.^2 * dt / pi;
info.seed = seed;
info.noiseSeed = nseed;
info.numRealisations = nR;
info.cacheStatus = twin.cacheStatus;
info.interpCheck = twin.interpCheck;
if strcmp(twin.cacheStatus, 'computed') || strcmp(twin.cacheStatus, 'loaded')
    twin.cacheStatus = 'memory';                    % next call with this twin reuses it
end
varargout = {Y, t, info, twin};
end

%% ================================================================================================
function twin = prepareTwin(twin, cacheDir, verbose)
if ~isstruct(twin) || ~isscalar(twin)
    error('synthesiseTwinRecords:badSpec', 'First argument must be a SPEC or TWIN struct.');
end
if isfield(twin, 'prepared') && isequal(twin.prepared, true)
    twin.cacheStatus = 'memory';
    return
end
% ---- spec fields and defaults ----
useHandle = isfield(twin, 'frf') && ~isempty(twin.frf);
if ~isfield(twin, 'sensors') || ~(isnumeric(twin.sensors) && size(twin.sensors, 2) == 2 && size(twin.sensors, 1) >= 2)
    error('synthesiseTwinRecords:badSpec', 'SPEC.sensors must be nS x 2 with nS >= 2.');
end
nS = size(twin.sensors, 1);
if ~useHandle
    if ~isfield(twin, 'p') || ~(isnumeric(twin.p) && numel(twin.p) == 3 && all(isfinite(twin.p)) && all(twin.p > 0))
        error('synthesiseTwinRecords:badSpec', 'SPEC.p must be three positive finite values [beta gamma R].');
    end
    twin.p = twin.p(:).';
elseif ~isa(twin.frf, 'function_handle')
    error('synthesiseTwinRecords:badSpec', 'SPEC.frf must be a function handle.');
elseif ~isfield(twin, 'p')
    twin.p = [];
end
if ~isfield(twin, 'frfOptions') || isempty(twin.frfOptions), twin.frfOptions = {}; end
if ~iscell(twin.frfOptions)
    error('synthesiseTwinRecords:badSpec', 'SPEC.frfOptions must be a cell of name-value pairs.');
end
if ~isfield(twin, 'band') || isempty(twin.band)
    twin.band = validatedBand(twin.frfOptions);
end
b = twin.band;
if ~(isnumeric(b) && numel(b) == 2 && all(isfinite(b)) && b(1) > 0 && b(2) > b(1))
    error('synthesiseTwinRecords:badSpec', 'SPEC.band must be [wLo wHi] with 0 < wLo < wHi.');
end
twin.band = b(:).';
if ~isfield(twin, 'nodeSpacing') || isempty(twin.nodeSpacing), twin.nodeSpacing = 0.1; end
if ~(isnumeric(twin.nodeSpacing) && isscalar(twin.nodeSpacing) && isfinite(twin.nodeSpacing) && twin.nodeSpacing > 0)
    error('synthesiseTwinRecords:badSpec', 'SPEC.nodeSpacing must be a positive scalar.');
end
if ~isfield(twin, 'checkPoints') || isempty(twin.checkPoints), twin.checkPoints = 3; end
cp = twin.checkPoints;
if ~(isnumeric(cp) && isscalar(cp) && cp >= 0 && cp == round(cp))
    error('synthesiseTwinRecords:badSpec', 'SPEC.checkPoints must be a nonnegative integer.');
end
if ~isfield(twin, 'checkReference') || isempty(twin.checkReference), twin.checkReference = nS; end
cr = twin.checkReference;
if ~(isnumeric(cr) && isscalar(cr) && cr == round(cr) && cr >= 1 && cr <= nS)
    error('synthesiseTwinRecords:badSpec', 'SPEC.checkReference must be a sensor index in 1..%d.', nS);
end
if useHandle && ~isempty(cacheDir)
    error('synthesiseTwinRecords:cacheNotAllowed', 'A function-handle FRF cannot be cached to disk.');
end

nNodes = max(2, ceil(diff(twin.band) / twin.nodeSpacing - 1e-9) + 1);
wN = linspace(twin.band(1), twin.band(2), nNodes);
if useHandle
    source = 'handle';
else
    source = 'emm';
end
key = struct('source', source, 'p', twin.p, 'sensors', twin.sensors, 'band', twin.band, ...
    'nodeSpacing', twin.nodeSpacing, 'frfOptions', {twin.frfOptions}, ...
    'checkPoints', twin.checkPoints, 'checkReference', twin.checkReference);

% ---- cache lookup ----
cacheFile = '';
loaded = false;
if ~isempty(cacheDir)
    if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
    d = dir(fullfile(cacheDir, 'twinFRF_*.mat'));
    for i = 1:numel(d)
        f = fullfile(cacheDir, d(i).name);
        try
            S = load(f, 'key');
        catch
            continue
        end
        if isfield(S, 'key') && keysMatch(S.key, key)
            S = load(f);
            wN = S.omegaNodes; HaN = S.HaNodes; chk = S.interpCheck; solveSec = S.solveSeconds;
            cacheFile = f; loaded = true;
            if verbose, fprintf('synthesiseTwinRecords: loaded node FRF from %s\n', f); end
            break
        end
    end
end

% ---- solve ----
if ~loaded
    t0 = tic;
    if verbose, fprintf('synthesiseTwinRecords: solving %d nodes in [%.4g %.4g] rad/s\n', nNodes, wN(1), wN(end)); end
    HaN = evalFRF(twin, source, wN, nS);
    solveSec = toc(t0);
end
HdN = HaN ./ (-wN.^2);
Hacc = @(w) interpAcc(wN, HdN, w);

if ~loaded
    chk = interpolationCheck(twin, source, wN, Hacc, nS, cp, cr);
    if ~isempty(cacheDir)
        stamp = datestr(now, 'yyyymmdd_HHMMSS');
        n = 1;
        cacheFile = fullfile(cacheDir, sprintf('twinFRF_%s_%02d.mat', stamp, n));
        while exist(cacheFile, 'file')
            n = n + 1;
            cacheFile = fullfile(cacheDir, sprintf('twinFRF_%s_%02d.mat', stamp, n));
        end
        omegaNodes = wN; HaNodes = HaN; interpCheck = chk; solveSeconds = solveSec; %#ok<NASGU>
        save(cacheFile, 'key', 'omegaNodes', 'HaNodes', 'interpCheck', 'solveSeconds', '-v7');
        if verbose, fprintf('synthesiseTwinRecords: saved node FRF to %s\n', cacheFile); end
    end
end

twin.prepared = true;
twin.source = source;
twin.key = key;
twin.omegaNodes = wN;
twin.HaNodes = HaN;
twin.HdNodes = HdN;
twin.Hacc = Hacc;
twin.interpCheck = chk;
twin.solveSeconds = solveSec;
twin.cacheFile = cacheFile;
if loaded, twin.cacheStatus = 'loaded'; else, twin.cacheStatus = 'computed'; end
end

%% ================================================================================================
function Ha = evalFRF(twin, source, w, nS)
if strcmp(source, 'handle')
    Ha = twin.frf(w);
    if ~(isnumeric(Ha) && isequal(size(Ha), [nS, numel(w)]) && all(isfinite(Ha(:))))
        error('synthesiseTwinRecords:badFRF', 'SPEC.frf must return a finite %d x nW matrix.', nS);
    end
else
    Ha = computeSensorFRF(w, twin.sensors, twin.p, 'Output', 'acceleration', twin.frfOptions{:});
end
end

function H = interpAcc(wN, HdN, w)
w = w(:).';
if any(w < wN(1) - 1e-12 * wN(1)) || any(w > wN(end) + 1e-12 * wN(end))
    error('synthesiseTwinRecords:badOption', 'Interpolated FRF requested outside the node band [%.6g %.6g].', wN(1), wN(end));
end
w = min(max(w, wN(1)), wN(end));
if numel(wN) >= 4
    method = 'spline';
else
    method = 'linear';
end
Hr = interp1(wN(:), real(HdN).', w(:), method);   % nW x nS
Hi = interp1(wN(:), imag(HdN).', w(:), method);
H = -(w.^2) .* complex(Hr, Hi).';
end

function chk = interpolationCheck(twin, source, wN, Hacc, nS, cp, cr)
chk = struct('omega', [], 'relErrH', [], 'relErrT', [], 'maxRelErrH', NaN, 'maxRelErrT', NaN, 'reference', cr);
nGap = numel(wN) - 1;
if cp == 0 || nGap < 1, return; end
g = unique(round(linspace(1, nGap, min(cp, nGap))));
wc = 0.5 * (wN(g) + wN(g + 1));
Hex = evalFRF(twin, source, wc, nS);
Hin = Hacc(wc);
others = setdiff(1:nS, cr);
relH = max(abs(Hin - Hex), [], 1) ./ max(abs(Hex), [], 1);
Tex = Hex(others, :) ./ Hex(cr, :);
Tin = Hin(others, :) ./ Hin(cr, :);
relT = max(abs(Tin - Tex), [], 1) ./ max(abs(Tex), [], 1);
chk.omega = wc; chk.relErrH = relH; chk.relErrT = relT;
chk.maxRelErrH = max(relH); chk.maxRelErrT = max(relT);
end

function tf = keysMatch(a, b)
tf = false;
try
    if ~strcmp(a.source, b.source), return; end
    if ~isequal(size(a.p), size(b.p)) || any(abs(a.p - b.p) > 1e-12 * abs(b.p)), return; end
    if ~isequal(size(a.sensors), size(b.sensors)) || any(abs(a.sensors(:) - b.sensors(:)) > 1e-12 * max(1, abs(b.sensors(:)))), return; end
    if ~isequal(a.band, b.band) || ~isequal(a.nodeSpacing, b.nodeSpacing), return; end
    if ~isequal(a.frfOptions, b.frfOptions), return; end
    if ~isequal(a.checkPoints, b.checkPoints) || ~isequal(a.checkReference, b.checkReference), return; end
    tf = true;
catch
    tf = false;
end
end

function b = validatedBand(frfOptions)
depth = 1.88; g = 9.81; vr = [1.7 13.85];
for a = 1:2:numel(frfOptions) - 1
    switch lower(char(frfOptions{a}))
        case 'waterdepth', depth = frfOptions{a + 1};
        case 'gravity', g = frfOptions{a + 1};
        case 'validatedalpharange', vr = frfOptions{a + 1};
    end
end
b = sqrt(vr(:).' * g / depth) .* [1 + 1e-9, 1 - 1e-9];   % just inside, so no validated-range warning
end

%% ================================================================================================
function [Sfun, spec] = buildSpectrum(opt)
s = opt.Spectrum;
if isa(s, 'string'), s = char(s); end
spec = struct('type', '', 'Hs', NaN, 'peakFrequency', NaN, 'peakEnhancement', NaN, 'whiteLevel', NaN, 'alpha', NaN);
if isa(s, 'function_handle')
    Sfun = s; spec.type = 'custom';
elseif ischar(s) && strcmpi(s, 'white')
    L = opt.WhiteLevel;
    if ~(isnumeric(L) && isscalar(L) && isfinite(L) && L >= 0)
        error('synthesiseTwinRecords:badSpectrum', 'WhiteLevel must be a nonnegative scalar.');
    end
    Sfun = @(w) L * ones(size(w)); spec.type = 'white'; spec.whiteLevel = L;
elseif ischar(s) && strcmpi(s, 'jonswap')
    Hs = opt.Hs; wp = opt.PeakFrequency; gJ = opt.PeakEnhancement;
    pos = @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0;
    if ~(pos(Hs) && pos(wp) && pos(gJ) && gJ >= 1)
        error('synthesiseTwinRecords:badSpectrum', 'Hs and PeakFrequency must be positive; PeakEnhancement >= 1.');
    end
    shape = @(w) jonswapShape(w, wp, gJ);
    wq = wp * logspace(log10(0.2), log10(20), 20000);   % full-spectrum m0 (shape is negligible outside)
    m0 = trapz(wq, shape(wq));
    alpha = (Hs / 4)^2 / m0;
    Sfun = @(w) alpha * jonswapShape(w, wp, gJ);
    spec.type = 'jonswap'; spec.Hs = Hs; spec.peakFrequency = wp; spec.peakEnhancement = gJ; spec.alpha = alpha;
else
    error('synthesiseTwinRecords:badSpectrum', 'Spectrum must be ''jonswap'', ''white'' or a function handle.');
end
end

function S = jonswapShape(w, wp, gJ)
sig = 0.07 * ones(size(w)); sig(w > wp) = 0.09;
S = w.^-5 .* exp(-1.25 * (wp ./ w).^4) .* gJ.^exp(-(w - wp).^2 ./ (2 * sig.^2 * wp^2));
end

function tpr = edgeTaper(w, band, width)
% Amplitude taper: 1 inside [band(1)+width, band(2)-width], raised-cosine
% ramps to 0 at the band edges.
tpr = ones(size(w));
if width <= 0, return; end
lo = w < band(1) + width;
tpr(lo) = 0.5 * (1 - cos(pi * (w(lo) - band(1)) / width));
hi = w > band(2) - width;
tpr(hi) = 0.5 * (1 - cos(pi * (band(2) - w(hi)) / width));
tpr = min(max(tpr, 0), 1);
end

function v = checkNoise(v, nS, name)
if ~(isnumeric(v) && isreal(v) && all(isfinite(v(:))) && all(v(:) >= 0) && (isscalar(v) || numel(v) == nS))
    error('synthesiseTwinRecords:badNoise', '%s must be a nonnegative scalar or a 1 x %d vector.', name, nS);
end
v = v(:).' .* ones(1, nS);
end

function g = drawStandardNormal(sz, seed)
if isempty(seed)
    g = randn(sz);
    return
end
if exist('OCTAVE_VERSION', 'builtin') > 0
    saved = randn('state'); %#ok<RAND>
    randn('state', seed);   %#ok<RAND>
    g = randn(sz);
    randn('state', saved);  %#ok<RAND>
else
    s = RandStream('mt19937ar', 'Seed', seed);
    g = randn(s, sz);
end
end