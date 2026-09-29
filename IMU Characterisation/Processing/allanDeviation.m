function out = allanDeviation(X, fs, varargin)
%ALLANDEVIATION  Overlapping Allan deviation of stationary records, with noise-term fits.
%   out = allanDeviation(X, fs, Name, Value, ...)
%
%   X    N x C, one column per channel, uniformly sampled, raw units (no filtering, no
%        integration by the caller: the rate/acceleration samples are passed as they are).
%   fs   sampling rate (Hz); for the shell logs the MEASURED received rate. [] = from 'Time'.
%
%   ESTIMATOR (overlapping, IEEE 952 / 1293 form)
%     theta_k = dt * cumsum(x), tau = m / fs,
%     sigma^2(tau) = sum_k (theta_{k+2m} - 2 theta_{k+m} + theta_k)^2 / (2 tau^2 (N - 2m + 1)),
%     m on a log grid from 1 to floor((N-1)/2).
%
%   NOISE TERMS, in THIS PROJECT'S conventions (the same as synthesiseImuNoise):
%     white        one-sided density N (units/sqrt(Hz)):  sigma_A = N / sqrt(2 tau)
%                  If N_IEEE is defined by sigma_A = N_IEEE / sqrt(tau) (IEEE 952, two-sided
%                  PSD), then N_IEEE = N / sqrt(2) = sigma_A(1 s); returned as whiteIEEE.
%                  Datasheet noise densities (ug/sqrt(Hz), mdps/sqrt(Hz)) are one-sided
%                  amplitude densities: compare them with 'white', not 'whiteIEEE'.
%     flicker      one-sided PSD B^2 / (2 pi f):  plateau sigma_A = sqrt(ln2 / pi) B = 0.470 B
%                  (bias instability quoted with the IEEE factor 0.664 is also returned)
%     random walk  increments K sqrt(dt):  sigma_A = K sqrt(tau / 3)
%   Fits use only the tau ranges where the local log-log slope is within SlopeTol of the
%   term's slope (-1/2, 0, +1/2). A term with no such range is returned as NaN: absence of
%   a region is information, not a fit failure to hide.
%
%   CONFIDENCE (approximate): equivalent degrees of freedom for white noise with overlapping
%   estimates, edf = (3 (N-1) / (2 m) - 2 (N-2) / N) * 4 m^2 / (4 m^2 + 5) (Howe, Allan and
%   Barnes), chi-square interval by Wilson-Hilferty; for other noise types it is indicative.
%
%   Name-value options
%     'PointsPerDecade'  tau grid density (default 10)
%     'SlopeTol'         slope tolerance for the term fits (default 0.15)
%     'MinEdf'           taus with fewer equivalent dof are plotted but not used by the fits
%                        or the plateau (default 4: the last few taus are very uncertain)
%     'Time'             N x 1 sample times: uniformity check (as computeNoisePSD) and fs if []
%     'MaxJitter' (0.05), 'MaxGapPeriods' (1.5)
%     'ChannelNames', 'Units' (default 'm/s^2'), 'Plot' (default false)
%
%   Output out: tau (nT x 1, s), m, adev (nT x C), avar, nTerms, edf, adevLo, adevHi (95%),
%     localSlope (nT-1 x C, midpoints tauMid), terms(c) with fields white, whiteIEEE, flicker
%     (B = plateau / 0.470), biasInstability (plateau / 0.664, the IEEE form),
%     randomWalk (K), and the tau range each was fitted over; plateau = mean sigma_A over the
%     flat run (NaN if none), minimum and tauAtMin (over trusted taus only); fs, N, units,
%     channelNames, timing.
%
% Lives in IMU Characterisation/Processing/.

[N, C] = size(X);
opt = struct('PointsPerDecade', 10, 'SlopeTol', 0.15, 'MinEdf', 4, 'Time', [], 'MaxJitter', 0.05, 'MaxGapPeriods', 1.5, ...
    'ChannelNames', {{}}, 'Units', 'm/s^2', 'Plot', false);
if mod(numel(varargin), 2) ~= 0, error('allanDeviation:args', 'Name-value pairs expected.'); end
for i = 1:2:numel(varargin)
    if ~isfield(opt, varargin{i}), error('allanDeviation:args', 'Unknown option ''%s''.', varargin{i}); end
    opt.(varargin{i}) = varargin{i + 1};
end
X = double(X);
if any(~isfinite(X(:))), error('allanDeviation:nonfinite', 'X contains NaN or Inf.'); end
if N < 9, error('allanDeviation:short', 'Need at least 9 samples.'); end

timing = struct('checked', false);
if ~isempty(opt.Time)
    t = double(opt.Time(:));
    if numel(t) ~= N, error('allanDeviation:time', 'Time must have one entry per row of X.'); end
    dt = diff(t); md = median(dt);
    timing = struct('checked', true, 'fsFromTime', (N - 1) / (t(end) - t(1)), 'medianDt', md, ...
        'jitter', std(dt) / md, 'maxGapPeriods', max(dt) / md);
    if isempty(fs), fs = timing.fsFromTime; end
    if timing.jitter > opt.MaxJitter || timing.maxGapPeriods > opt.MaxGapPeriods
        error('allanDeviation:irregular', ['Sampling is not uniform: std(dt)/median(dt) = %.3f (max %.3f), ' ...
            'max dt = %.2f periods (max %.2f).'], timing.jitter, opt.MaxJitter, timing.maxGapPeriods, opt.MaxGapPeriods);
    end
end
if isempty(fs), error('allanDeviation:fs', 'fs is required when Time is not given.'); end
tau0 = 1 / fs;

% ---- tau grid and overlapping estimator ------------------------------------------------------------
mMax = floor((N - 1) / 2);
m = unique(round(logspace(0, log10(mMax), max(2, ceil(opt.PointsPerDecade * log10(mMax)) + 1))));
m = m(:);
theta = [zeros(1, C); cumsum(X, 1)] * tau0;                          % (N+1) x C
nT = numel(m);
avar = zeros(nT, C); nTerms = N - 2 * m + 1;
for i = 1:nT
    k = m(i);
    d = theta(1 + 2 * k:end, :) - 2 * theta(1 + k:end - k, :) + theta(1:end - 2 * k, :);
    tau = k * tau0;
    avar(i, :) = sum(d.^2, 1) / (2 * tau^2 * (N - 2 * k + 1));
end
tau = m * tau0;
adev = sqrt(avar);

% ---- approximate confidence --------------------------------------------------------------------------
edf = (3 * (N - 1) ./ (2 * m) - 2 * (N - 2) / N) .* (4 * m.^2 ./ (4 * m.^2 + 5));
edf = max(edf, 1);
z = 1.959964;
chiq = @(nu, q) nu .* (1 - 2 ./ (9 * nu) + q * sqrt(2 ./ (9 * nu))).^3;
lo = sqrt(edf ./ max(chiq(edf, z), eps)); hi = sqrt(edf ./ max(chiq(edf, -z), eps));
adevLo = adev .* lo; adevHi = adev .* hi;

% ---- local slopes and noise-term fits -----------------------------------------------------------------
lt = log10(tau); la = log10(adev);
slope = diff(la, 1, 1) ./ diff(lt);
tauMid = 10.^((lt(1:end - 1) + lt(2:end)) / 2);
terms = struct([]);
use = edf >= opt.MinEdf;                                              % taus trusted by the fits
nu = find(use, 1, 'last'); if isempty(nu), nu = 2; end
for c = 1:C
    s = slope(1:nu - 1, c); tt = tau(1:nu); aa = adev(1:nu, c);
    [wN, wR] = fitTerm(tt, aa, s, -0.5, opt.SlopeTol, @(x, y) y .* sqrt(2 * x));
    [kK, kR] = fitTerm(tt, aa, s, 0.5, opt.SlopeTol, @(x, y) y .* sqrt(3 ./ x));
    [pl, fR] = fitTerm(tt, aa, s, 0, opt.SlopeTol, @(x, y) y);      % plateau: mean over the flat run
    [amin, iMin] = min(aa);
    if isnan(pl), fl = NaN; bi = NaN; else, fl = pl / sqrt(log(2) / pi); bi = pl / 0.664; end
    terms = [terms, struct('white', wN, 'whiteIEEE', wN / sqrt(2), 'whiteRange', wR, ...
        'flicker', fl, 'biasInstability', bi, 'plateau', pl, 'flatRange', fR, 'minimum', amin, ...
        'tauAtMin', tt(iMin), 'randomWalk', kK, 'randomWalkRange', kR)]; %#ok<AGROW>
end

names = opt.ChannelNames;
if isempty(names), names = arrayfun(@(c) sprintf('ch%d', c), 1:C, 'UniformOutput', false); end
out = struct('tau', tau, 'm', m, 'adev', adev, 'avar', avar, 'nTerms', nTerms, 'edf', edf, ...
    'adevLo', adevLo, 'adevHi', adevHi, 'localSlope', slope, 'tauMid', tauMid, 'terms', terms, ...
    'fs', fs, 'N', N, 'units', opt.Units, 'channelNames', {names}, 'timing', timing);
if opt.Plot, plotADEV(out); end
end

function [val, rangeTau] = fitTerm(tau, a, s, target, tol, conv)
% Mean of conv(tau, a) over the longest run of consecutive grid intervals whose local slope
% is within tol of target (both end points of each interval included).
ok = abs(s - target) <= tol;
best = [0 0]; i = 1; n = numel(ok);
while i <= n
    if ok(i)
        j = i; while j < n && ok(j + 1), j = j + 1; end
        if j - i + 1 > best(2) - best(1) + (best(2) > 0), best = [i j]; end
        i = j + 1;
    else
        i = i + 1;
    end
end
if best(2) == 0 || best(2) - best(1) + 1 < 2
    val = NaN; rangeTau = [NaN NaN]; return
end
idx = best(1):best(2) + 1;
val = exp(mean(log(conv(tau(idx), a(idx)))));
rangeTau = tau(idx([1 end])).';
end

function plotADEV(out)
figure('Color', 'w', 'Position', [60 60 760 520]); hold on;
h = zeros(1, size(out.adev, 2));
for c = 1:size(out.adev, 2)
    h(c) = plot(out.tau, out.adev(:, c), '-', 'LineWidth', 1.5);
    col = get(h(c), 'Color');
    plot(out.tau, out.adevLo(:, c), ':', 'Color', col); plot(out.tau, out.adevHi(:, c), ':', 'Color', col);
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on;
xlabel('\tau (s)'); ylabel(['Allan deviation (', out.units, ')']);
title(sprintf('Overlapping Allan deviation, N = %d, fs = %.4g Hz (dotted: approx. 95%%)', out.N, out.fs));
legend(h, out.channelNames, 'Interpreter', 'none', 'Location', 'southwest');
end