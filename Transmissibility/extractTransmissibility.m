function est = extractTransmissibility(Y, dt, varargin)
%EXTRACTTRANSMISSIBILITY Output-only transmissibility estimate from
%multi-channel time series (Welch cross-spectra, H1-type ratio).
%
%   EST = EXTRACTTRANSMISSIBILITY(Y, DT)
%   EST = EXTRACTTRANSMISSIBILITY(Y, DT, Name, Value)
%
%   Y  - nS x n real time series (channels x samples), nS >= 2, e.g. the
%        four IMU acceleration records. Same layout as buildHankelMatrix.
%   DT - sampling period (s).
%
%   ESTIMATOR
%   For reference channel r and every other channel j,
%       T_hat_j(omega) = S_jr(omega) / S_rr(omega),
%       S_jr = E[Y_j Y_r^*]   (orientation matters: with y_j = H_j u this
%                              gives H_j H_r^* S_u, so the ratio is H_j/H_r;
%                              the conjugate orientation would return
%                              conj(T) or its reciprocal).
%   Spectra by Welch's method: segments of 'SegmentLength' samples,
%   fractional 'Overlap', periodic Hann window, per-segment mean removed.
%   The FFT kernel exp(-i omega t) makes the coefficient at +omega
%   proportional to the model's H(omega) under the exp(+i omega t)
%   convention of the Forward Model and stochasticSensorSynthesis, so T_hat
%   estimates exactly the model's H_j/H_r.
%
%   Frequencies: omega_k = 2 pi k / (SegmentLength * DT), k = 1 .. L/2 - 1
%   (DC and Nyquist excluded), optionally restricted to 'Band'. This
%   vector (EST.omega) is the grid on which the model side
%   (transmissibilityFeatureVector) must be evaluated.
%
%   PSD NORMALISATION: one-sided, per rad/s, so the integral over
%   [0, pi/DT] equals the variance. White noise of std sigma sampled at DT
%   has PSD sigma^2 DT / pi. 'NoisePSD' must use the same convention.
%
%   REFERENCE-NOISE CORRECTION
%   Independent sensor noise does not bias S_jr, but it biases the
%   reference auto-spectrum: S_rr_hat = S_rr + S_nn. With 'NoisePSD' (the
%   reference sensor's noise PSD) the denominator becomes S_rr_hat - S_nn.
%   Bins where the corrected denominator is below 'MinDenominatorFraction'
%   of the uncorrected one (or non-positive) are flagged invalid and their
%   corrected T is NaN; they are never silently used. EST.Tuncorrected is
%   always returned.
%
%   OUTPUT STRUCT EST
%     omega          1 x nF   frequencies (rad/s)
%     binIndex       1 x nF   FFT bin indices k
%     reference      r
%     others         1 x m    channel indices of the rows of T (all but r, ascending)
%     T              m x nF   transmissibility (corrected if NoisePSD given)
%     Tuncorrected   m x nF
%     valid          1 x nF   denominator check passed
%     Sjr            m x nF   cross-spectra E[Y_j Y_r^*]
%     Srr            1 x nF   reference auto-spectrum (uncorrected)
%     Sjj            m x nF   auto-spectra of the other channels
%     coherence      m x nF   |S_jr|^2 / (S_jj S_rr)   (uncorrected spectra)
%     nSegments, nEffective   raw and effective number of averages
%                             (Welch overlap correction for the window used)
%     varTdiag       m x nF   ANALYTICAL DIAGNOSTIC ONLY: E|T_hat - T|^2 ~
%                             |T|^2 (1 - coh) / (coh nEffective) (Bendat and
%                             Piersol, single-input H1). It ignores the
%                             correlation between the m ratios (they share
%                             the reference channel); the real covariance
%                             comes from Monte Carlo.
%     settings                options used
%
%   NAME-VALUE OPTIONS
%     'Reference'               channel index r (default 1)
%     'SegmentLength'           even integer L (default 512)
%     'Overlap'                 fraction in [0, 1) (default 0.5)
%     'Band'                    [wLo wHi] rad/s, or [] for all bins (default [])
%     'NoisePSD'                [] (default), a scalar (white), a vector with
%                               one value per returned frequency, or a
%                               function handle @(omega)
%     'MinDenominatorFraction'  default 0.1
%
%   ERROR IDENTIFIERS (extractTransmissibility:<id>)
%     badY, badDt, badOption, tooShort, emptyBand, badNoisePSD
%
%   Lives in Transmissibility/. Tester:
%   Transmissibility/Unit Tests Transmissibility/extractTransmissibilityTester.m

%% ---- Inputs ---------------------------------------------------------------
if ~(isnumeric(Y) && isreal(Y) && ismatrix(Y) && size(Y, 1) >= 2 && all(isfinite(Y(:))))
    error('extractTransmissibility:badY', 'Y must be a finite real nS x n matrix with nS >= 2 channels.');
end
if ~(isnumeric(dt) && isscalar(dt) && isreal(dt) && isfinite(dt) && dt > 0)
    error('extractTransmissibility:badDt', 'DT must be a positive finite scalar.');
end
[nS, n] = size(Y);

opt = struct('Reference', 1, 'SegmentLength', 512, 'Overlap', 0.5, 'Band', [], ...
    'NoisePSD', [], 'MinDenominatorFraction', 0.1);
if mod(numel(varargin), 2) ~= 0
    error('extractTransmissibility:badOption', 'Options must be name-value pairs.');
end
names = fieldnames(opt);
for a = 1:2:numel(varargin)
    key = varargin{a};
    if isa(key, 'string'), key = char(key); end
    hit = ischar(key) & strcmpi(key, names);
    if ~any(hit)
        error('extractTransmissibility:badOption', 'Unknown option.');
    end
    opt.(names{hit}) = varargin{a + 1};
end
r = opt.Reference;
if ~(isnumeric(r) && isscalar(r) && r == round(r) && r >= 1 && r <= nS)
    error('extractTransmissibility:badOption', 'Reference must be a channel index in 1..%d.', nS);
end
L = opt.SegmentLength;
if ~(isnumeric(L) && isscalar(L) && L == round(L) && L >= 8 && mod(L, 2) == 0)
    error('extractTransmissibility:badOption', 'SegmentLength must be an even integer >= 8.');
end
ov = opt.Overlap;
if ~(isnumeric(ov) && isscalar(ov) && isfinite(ov) && ov >= 0 && ov < 1)
    error('extractTransmissibility:badOption', 'Overlap must be in [0, 1).');
end
mdf = opt.MinDenominatorFraction;
if ~(isnumeric(mdf) && isscalar(mdf) && isfinite(mdf) && mdf >= 0 && mdf < 1)
    error('extractTransmissibility:badOption', 'MinDenominatorFraction must be in [0, 1).');
end
if n < L
    error('extractTransmissibility:tooShort', 'Record length %d is shorter than SegmentLength %d.', n, L);
end

%% ---- Frequency grid -----------------------------------------------------------
k = 1:L / 2 - 1;
omega = 2 * pi * k / (L * dt);
if ~isempty(opt.Band)
    b = opt.Band;
    if ~(isnumeric(b) && numel(b) == 2 && all(isfinite(b)) && b(2) > b(1))
        error('extractTransmissibility:badOption', 'Band must be [wLo wHi] with wHi > wLo.');
    end
    keep = omega >= b(1) & omega <= b(2);
    if ~any(keep)
        error('extractTransmissibility:emptyBand', 'No Welch bin lies inside Band [%g %g] rad/s.', b(1), b(2));
    end
    k = k(keep); omega = omega(keep);
end
nF = numel(k);

%% ---- Welch cross-spectra ------------------------------------------------------------
w = 0.5 * (1 - cos(2 * pi * (0:L - 1).' / L));      % periodic Hann
U = sum(w.^2);
hop = max(1, round(L * (1 - ov)));
starts = 1:hop:(n - L + 1);
K = numel(starts);
scale = 2 * dt / (U * 2 * pi);                        % one-sided, per rad/s
others = setdiff(1:nS, r);
m = numel(others);
Sjr = complex(zeros(m, nF)); Sjj = zeros(m, nF); Srr = zeros(1, nF);
for s = 1:K
    seg = Y(:, starts(s):starts(s) + L - 1);
    seg = seg - mean(seg, 2);
    X = fft(seg .* w.', [], 2);                       % nS x L
    X = X(:, k + 1);
    Xr = X(r, :);
    Sjr = Sjr + X(others, :) .* conj(Xr);             % E[Y_j Y_r^*]
    Sjj = Sjj + abs(X(others, :)).^2;
    Srr = Srr + abs(Xr).^2;
end
Sjr = scale * Sjr / K; Sjj = scale * Sjj / K; Srr = scale * Srr / K;

%% ---- Ratios, correction, validity ------------------------------------------------------
Tunc = Sjr ./ Srr;
coh = abs(Sjr).^2 ./ (Sjj .* Srr);
valid = true(1, nF);
T = Tunc;
if ~isempty(opt.NoisePSD)
    Snn = noiseAtBins(opt.NoisePSD, omega);
    den = Srr - Snn;
    valid = den > mdf * Srr;
    T = Sjr ./ den;
    T(:, ~valid) = NaN;
end

% effective number of averages for overlapped windows (Welch 1967)
nEff = K;
if K > 1
    rho2 = 0;
    for l = 1:K - 1
        sh = l * hop;
        if sh >= L, break; end
        c = sum(w(1:L - sh) .* w(1 + sh:L)) / U;
        rho2 = rho2 + (1 - l / K) * c^2;
    end
    nEff = K / (1 + 2 * rho2);
end
cohSafe = min(max(coh, eps), 1);
varTdiag = abs(Tunc).^2 .* (1 - cohSafe) ./ (cohSafe * nEff);

%% ---- Output -------------------------------------------------------------------------------
est = struct();
est.omega = omega;
est.binIndex = k;
est.reference = r;
est.others = others;
est.T = T;
est.Tuncorrected = Tunc;
est.valid = valid;
est.Sjr = Sjr;
est.Srr = Srr;
est.Sjj = Sjj;
est.coherence = coh;
est.nSegments = K;
est.nEffective = nEff;
est.varTdiag = varTdiag;
est.settings = struct('dt', dt, 'SegmentLength', L, 'Overlap', ov, 'hop', hop, 'Band', opt.Band, ...
    'NoisePSD', opt.NoisePSD, 'MinDenominatorFraction', mdf, 'window', 'periodic Hann');
end

%% ================================================================================================
function Snn = noiseAtBins(spec, omega)
if isa(spec, 'function_handle')
    Snn = spec(omega);
elseif isnumeric(spec) && isscalar(spec)
    Snn = spec * ones(size(omega));
elseif isnumeric(spec) && numel(spec) == numel(omega)
    Snn = reshape(spec, size(omega));
else
    error('extractTransmissibility:badNoisePSD', ...
        'NoisePSD must be a scalar, a vector with one value per returned frequency (%d), or a function handle.', numel(omega));
end
if ~(isreal(Snn) && all(isfinite(Snn)) && all(Snn >= 0) && isequal(size(Snn), size(omega)))
    error('extractTransmissibility:badNoisePSD', 'NoisePSD must be real, finite and nonnegative at every returned frequency.');
end
end