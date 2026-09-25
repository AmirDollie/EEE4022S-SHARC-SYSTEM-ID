function [y, t, info] = stochasticSensorSynthesis(omega, H, dt, varargin)
%STOCHASTICSENSORSYNTHESIS Multi-channel stochastic time series from a
%sampled frequency response, by frequency-domain synthesis (Y = H*U).
%
%   [Y, T, INFO] = STOCHASTICSENSORSYNTHESIS(OMEGA, H, DT)
%   [Y, T, INFO] = STOCHASTICSENSORSYNTHESIS(..., Name, Value)
%
%   OMEGA - 1 x nW (or nW x 1) angular frequencies (rad/s), uniformly
%           spaced with spacing dOmega, lying EXACTLY on the FFT grid of
%           the record (see GRID RULES). Positive; DC is never excited.
%   H     - nS x nW complex frequency response: H(j,k) is sensor j's
%           response at OMEGA(k) per unit input (e.g. acceleration FRF
%           H_a = -omega.^2 .* H_eta). nS = 1 may also be given as nW x 1.
%   DT    - sampling period (s).
%
%   Y     - nS x N x nR real time series (nR = 'NumRealisations').
%           Y(:,:,r) is exactly the (channels x samples) shape
%           buildHankelMatrix.m expects.
%   T     - 1 x N sample times, T = (0:N-1)*DT.
%   INFO  - struct: grid, taper, drawn input, and the EXPECTED second-order
%           statistics (see OUTPUT STATISTICS).
%
%   MODEL
%   One complex input coefficient per frequency, SHARED by all sensors:
%       A_k = w_k * sqrt(S_u(omega_k) * dOmega) * H(:,k) * xi_k,
%       y(t) = Re{ sum_k A_k exp(i*omega_k*t) }.
%   Sharing xi_k across sensors preserves the spatial coherence produced
%   by the forward model (single incident wave per frequency). Drawing
%   independent noise per sensor would destroy it, and with it the mode
%   shape information SSI is supposed to recover.
%     'gaussian'    : xi_k = a + i*b, a, b ~ N(0,1) independent, so
%                     E|xi_k|^2 = 2 (circular complex Gaussian; |xi_k| is
%                     Rayleigh). This is a stationary Gaussian process
%                     sample: the output-only / white-input assumption SSI
%                     is built on.
%     'randomPhase' : xi_k = sqrt(2) * exp(i*phi_k), phi_k ~ U[0, 2pi).
%                     Same expected power, deterministic amplitudes. This
%                     is the older forced-tone style, kept only as a
%                     comparison case for level4A.
%   Synthesis: A_k/2 is placed in FFT bin k and conj(A_k)/2 in bin N-k
%   (explicit Hermitian symmetry), then y = N * ifft(.). DC and all bins
%   outside OMEGA are zero, so every realisation has zero mean exactly.
%
%   GRID RULES (checked; errors if violated)
%     dOmega = OMEGA(2) - OMEGA(1), uniform to relative tolerance 'Tol'.
%     N = 2*pi / (dOmega * DT) must be an integer: the record length is
%         exactly one period, T_rec = N*DT = 2*pi/dOmega.
%     OMEGA / dOmega must be integers >= 1 (bins aligned with the FFT grid).
%     max(OMEGA) must lie strictly below the Nyquist frequency pi/DT.
%   Build a compatible grid like this:
%       N = 4096; dt = 0.1; dOmega = 2*pi/(N*dt);
%       omega = (ceil(4/dOmega):floor(8.4/dOmega)) * dOmega;
%
%   CIRCULARITY (important for interpretation)
%   The record is exactly periodic with period T_rec = 2*pi/dOmega: the
%   excitation is a finite set of tones spaced dOmega apart. Refining the
%   grid (smaller dOmega) lengthens the record. Whether SSI then sees
%   damped poles or undamped tones depends on how many bins fall across
%   each half-power bandwidth, N_BW = 2*zeta*omega_n / dOmega. That
%   transition is what level4A_broadbandKnownOscillator.m measures; this
%   function makes no claim about it.
%
%   NAME-VALUE OPTIONS
%   'NumRealisations'   positive integer (default 1).
%   'Seed'              [] (default: draw from the current global stream)
%                       or a nonnegative integer. A seed draws from a LOCAL
%                       stream, so the caller's global RNG state is left
%                       untouched.
%   'InputPSD'          one-sided input PSD S_u per rad/s: scalar (flat,
%                       default 1) or 1 x nW nonnegative vector (e.g. a
%                       wave spectrum later).
%   'Taper'             'none' (default), 'cosine', or a 1 x nW AMPLITUDE
%                       taper w in [0, 1]. The PSD is multiplied by w.^2.
%   'TaperFraction'     for 'cosine': fraction of the nW bins ramped at
%                       EACH band edge, in [0, 0.5] (default 0.1). Ramp is a
%                       raised cosine in amplitude. Purpose: avoid a hard
%                       band edge, which can itself be fitted as a pole.
%   'AmplitudeMode'     'gaussian' (default) or 'randomPhase'.
%   'InputCoefficients' nW x nR complex xi to use instead of drawing (for
%                       testing / common random numbers across parameter
%                       perturbations). Overrides Seed and AmplitudeMode.
%   'Tol'               relative tolerance for the grid checks (1e-9).
%
%   OUTPUT STATISTICS (all one-sided, per rad/s)
%     info.expectedPSD          nS x nW:  G_j(omega_k) = S_u w^2 |H_jk|^2
%     info.expectedCrossSpectrum nS x nS x nW:  C_k = S_u w^2 H_k H_k^H
%                               (rank 1 at every frequency: single input)
%     info.expectedCovariance   nS x nS:  E[y y^T] = dOmega * sum_k Re(C_k)
%   Because each record is exactly one period, fft(Y(j,:,r)) returns
%   (N/2) * A_k in bin k+1 with no leakage; the tester uses this.
%
%   ERROR IDENTIFIERS (stochasticSensorSynthesis:<id>)
%     badOmega, nonUniformGrid, gridNotFFTCompatible, offGridFrequency,
%     aboveNyquist, badDt, badH, sizeMismatch, badOption, badTaper,
%     badPSD, badCoefficients, synthesisNotReal
%
%   Lives in SSI/. Tester: SSI/Unit Tests SSI/stochasticSensorSynthesisTester.m

%% ---- Required inputs ---------------------------------------------------
if ~(isnumeric(omega) && isreal(omega) && isvector(omega) && numel(omega) >= 2 ...
        && all(isfinite(omega)) && all(omega > 0) && all(diff(omega(:)) > 0))
    error('stochasticSensorSynthesis:badOmega', ...
        'OMEGA must be a real, finite, positive, strictly increasing vector with at least 2 entries.');
end
omega = omega(:).';
nW = numel(omega);

if ~(isnumeric(dt) && isreal(dt) && isscalar(dt) && isfinite(dt) && dt > 0)
    error('stochasticSensorSynthesis:badDt', 'DT must be a positive finite real scalar.');
end

if ~(isnumeric(H) && ismatrix(H) && ~isempty(H) && all(isfinite(H(:))))
    error('stochasticSensorSynthesis:badH', 'H must be a finite numeric (complex) matrix.');
end
if iscolumn(H) && numel(H) == nW && nW > 1
    H = H.';                               % single sensor given as nW x 1
end
if size(H, 2) ~= nW
    error('stochasticSensorSynthesis:sizeMismatch', ...
        'H must be nS x nW with nW = numel(OMEGA) = %d; got %d x %d.', nW, size(H, 1), size(H, 2));
end
nS = size(H, 1);

%% ---- Options -------------------------------------------------------------
opt = struct('NumRealisations', 1, 'Seed', [], 'InputPSD', 1, 'Taper', 'none', ...
    'TaperFraction', 0.1, 'AmplitudeMode', 'gaussian', 'InputCoefficients', [], 'Tol', 1e-9);
if mod(numel(varargin), 2) ~= 0
    error('stochasticSensorSynthesis:badOption', 'Options must be name-value pairs.');
end
names = fieldnames(opt);
for a = 1:2:numel(varargin)
    key = varargin{a};
    if isa(key, 'string'), key = char(key); end
    if ~ischar(key)
        error('stochasticSensorSynthesis:badOption', 'Option names must be text.');
    end
    hit = strcmpi(key, names);
    if ~any(hit)
        error('stochasticSensorSynthesis:badOption', 'Unknown option ''%s''.', key);
    end
    opt.(names{hit}) = varargin{a + 1};
end

nR = opt.NumRealisations;
if ~(isnumeric(nR) && isscalar(nR) && isfinite(nR) && nR >= 1 && nR == round(nR))
    error('stochasticSensorSynthesis:badOption', 'NumRealisations must be a positive integer.');
end
tol = opt.Tol;
if ~(isnumeric(tol) && isscalar(tol) && isfinite(tol) && tol > 0 && tol < 1e-3)
    error('stochasticSensorSynthesis:badOption', 'Tol must be a positive scalar below 1e-3.');
end
seed = opt.Seed;
if ~isempty(seed) && ~(isnumeric(seed) && isscalar(seed) && isfinite(seed) && seed >= 0 && seed == round(seed))
    error('stochasticSensorSynthesis:badOption', 'Seed must be [] or a nonnegative integer.');
end
ampMode = opt.AmplitudeMode;
if isa(ampMode, 'string'), ampMode = char(ampMode); end
if ~(ischar(ampMode) && any(strcmpi(ampMode, {'gaussian', 'randomPhase'})))
    error('stochasticSensorSynthesis:badOption', 'AmplitudeMode must be ''gaussian'' or ''randomPhase''.');
end
ampMode = lower(ampMode);

%% ---- Grid checks -----------------------------------------------------------
dOmega = omega(2) - omega(1);
if max(abs(diff(omega) - dOmega)) > tol * dOmega * max(1, nW)
    error('stochasticSensorSynthesis:nonUniformGrid', 'OMEGA must be uniformly spaced.');
end
dOmega = (omega(end) - omega(1)) / (nW - 1);      % least round-off estimate

Nreal = 2 * pi / (dOmega * dt);
N = round(Nreal);
if abs(Nreal - N) > tol * max(1, Nreal) || N < 2
    error('stochasticSensorSynthesis:gridNotFFTCompatible', ...
        ['2*pi/(dOmega*DT) = %.10g is not an integer, so OMEGA is not on an FFT grid of this DT. ' ...
         'For this DT the nearest valid dOmega is %.10g (N = %d); for this dOmega the nearest valid DT is %.10g.'], ...
        Nreal, 2 * pi / (max(N, 2) * dt), max(N, 2), 2 * pi / (max(N, 2) * dOmega));
end
dOmega = 2 * pi / (N * dt);                       % exact grid spacing

kReal = omega / dOmega;
k = round(kReal);
if any(abs(kReal - k) > tol * max(1, kReal)) || any(k < 1)
    error('stochasticSensorSynthesis:offGridFrequency', ...
        'OMEGA must be integer multiples (>= 1) of dOmega = %.10g; the grid is offset from the FFT bins.', dOmega);
end
if any(2 * k >= N)
    error('stochasticSensorSynthesis:aboveNyquist', ...
        'max(OMEGA) = %.6g rad/s is not strictly below the Nyquist frequency pi/DT = %.6g rad/s.', ...
        omega(end), pi / dt);
end

%% ---- Input PSD and taper -----------------------------------------------------
Su = opt.InputPSD;
if ~(isnumeric(Su) && isreal(Su) && all(isfinite(Su(:))) && all(Su(:) >= 0) && (isscalar(Su) || numel(Su) == nW))
    error('stochasticSensorSynthesis:badPSD', 'InputPSD must be a nonnegative scalar or a vector with numel(OMEGA) entries.');
end
Su = Su(:).' .* ones(1, nW);

w = buildTaper(opt.Taper, opt.TaperFraction, nW);

%% ---- Input coefficients xi (nW x nR) --------------------------------------------
if ~isempty(opt.InputCoefficients)
    xi = opt.InputCoefficients;
    if ~(isnumeric(xi) && isequal(size(xi), [nW, nR]) && all(isfinite(xi(:))))
        error('stochasticSensorSynthesis:badCoefficients', ...
            'InputCoefficients must be a finite nW x NumRealisations array (%d x %d).', nW, nR);
    end
    ampMode = 'supplied';
else
    g = drawStandardNormal([nW, 2 * nR], seed);
    switch ampMode
        case 'gaussian'
            xi = g(:, 1:nR) + 1i * g(:, nR + 1:end);
        case 'randomphase'
            % Map the first normal block to a uniform phase via the normal CDF:
            % keeps a single draw path (and hence a single seed contract).
            u = 0.5 * erfc(-g(:, 1:nR) / sqrt(2));
            xi = sqrt(2) * exp(1i * 2 * pi * u);
    end
end

%% ---- Synthesis ---------------------------------------------------------------
scale = w .* sqrt(Su * dOmega);                   % 1 x nW
y = zeros(nS, N, nR);
for r = 1:nR
    A = H .* (scale .* xi(:, r).');               % nS x nW complex amplitudes
    Z = zeros(N, nS);
    Z(k + 1, :) = A.' / 2;
    Z(N - k + 1, :) = conj(A.') / 2;              % explicit Hermitian symmetry
    yr = N * ifft(Z, [], 1);                       % N x nS
    magRe = max(abs(real(yr(:))));
    magIm = max(abs(imag(yr(:))));
    if magIm > 1e-10 * max(magRe, realmin)
        error('stochasticSensorSynthesis:synthesisNotReal', ...
            'Inverse FFT is not real (max|imag| = %.3g vs max|real| = %.3g).', magIm, magRe);
    end
    y(:, :, r) = real(yr).';
end
t = (0:N - 1) * dt;

%% ---- Info ----------------------------------------------------------------------------
if nargout > 2
    C = zeros(nS, nS, nW);
    for q = 1:nW
        C(:, :, q) = Su(q) * w(q)^2 * (H(:, q) * H(:, q)');
    end
    info = struct();
    info.dt = dt;
    info.N = N;
    info.recordLength = N * dt;
    info.dOmega = dOmega;
    info.omega = k * dOmega;
    info.binIndex = k;                              % FFT bin = binIndex + 1
    info.omegaNyquist = pi / dt;
    info.numRealisations = nR;
    info.amplitudeMode = ampMode;
    info.seed = seed;
    info.inputPSD = Su;
    info.taper = w;                                 % amplitude taper
    info.psdTaper = w.^2;
    info.inputCoefficients = xi;                    % nW x nR
    info.expectedPSD = (Su .* w.^2) .* abs(H).^2;   % nS x nW
    info.expectedCrossSpectrum = C;                 % nS x nS x nW
    info.expectedCovariance = dOmega * real(sum(C, 3));
end
end

%% ================================================================================
function w = buildTaper(spec, frac, nW)
if isa(spec, 'string'), spec = char(spec); end
if ischar(spec)
    switch lower(spec)
        case 'none'
            w = ones(1, nW);
        case 'cosine'
            if ~(isnumeric(frac) && isscalar(frac) && isfinite(frac) && frac >= 0 && frac <= 0.5)
                error('stochasticSensorSynthesis:badTaper', 'TaperFraction must be in [0, 0.5].');
            end
            nRamp = round(frac * nW);
            w = ones(1, nW);
            if nRamp > 0
                ramp = 0.5 * (1 - cos(pi * (1:nRamp) / (nRamp + 1)));
                w(1:nRamp) = ramp;
                w(end - nRamp + 1:end) = fliplr(ramp);
            end
        otherwise
            error('stochasticSensorSynthesis:badTaper', 'Taper must be ''none'', ''cosine'' or a numeric vector.');
    end
elseif isnumeric(spec) && isreal(spec) && numel(spec) == nW && all(isfinite(spec(:))) ...
        && all(spec(:) >= 0) && all(spec(:) <= 1)
    w = spec(:).';
else
    error('stochasticSensorSynthesis:badTaper', ...
        'A numeric Taper must have numel(OMEGA) entries in [0, 1].');
end
end

%% ================================================================================
function g = drawStandardNormal(sz, seed)
% Standard normal draws. With a seed, use a local stream so the caller's
% global RNG state is not disturbed.
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