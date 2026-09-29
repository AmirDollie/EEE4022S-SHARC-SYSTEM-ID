function out = computeNoisePSD(X, fs, varargin)
%COMPUTENOISEPSD  One-sided Welch PSD / ASD, cross-spectra and coherence of stationary records.
%   out = computeNoisePSD(X, fs, Name, Value, ...)
%
%   X    N x C, one column per channel, uniformly sampled (e.g. the same axis of three
%        co-located IMUs, aligned tick by tick with alignImuSensors). Raw units, no filtering.
%   fs   sampling rate (Hz). For the shell logs use the MEASURED received rate (~9.75 Hz at a
%        commanded 10 Hz), not the command; [] = estimate from 'Time'.
%
%   CONVENTION (fixed, matches the transmissibility estimator)
%     one-sided, per Hz:        S(f) [units^2/Hz], f = 0 .. fs/2, integral over f = variance
%     amplitude density:        A(f) = sqrt(S(f)) [units/sqrt(Hz)]
%     per rad/s (inversion):    S_w(w) = S(f) / (2 pi), w = 2 pi f; white noise of std sigma
%                               sampled at dt gives S_w = sigma^2 dt / pi (the twin's S_nn)
%     Welch: segments of L samples, Hann window (or rectangular), overlap fraction, mean
%     removed per segment (or linear trend), P = 2 |FFT(w x)|^2 / (fs sum w^2), DC and Nyquist
%     not doubled.
%
%   Name-value options
%     'SegmentLength'    L (default: largest power of 2 <= N/4, i.e. >= 7 segments at 50%)
%     'Overlap'          fraction in [0, 1) (default 0.5)
%     'Window'           'hann' (default) or 'rect'
%     'Detrend'          'mean' (default), 'linear' or 'none', per segment
%     'Pairs'            P x 2 channel pairs for cross-spectra / coherence (default: all pairs
%                        when C <= 6, none otherwise); [] = none
%     'ChannelNames'     cellstr (default ch1..chC);  'Units' string (default 'm/s^2')
%     'Bands'            struct of [fLo fHi] (Hz) for band summaries (default basin
%                        [0.47 1.35] and miz [0.08 0.17])
%     'Time'             N x 1 sample times (s): uniformity check, and fs estimate if fs = []
%     'MaxJitter'        max std(dt)/median(dt) before refusing (default 0.05)
%     'MaxGapPeriods'    max dt / median(dt) before refusing (default 1.5): Welch needs a
%                        uniform grid; a broken record (e.g. the 50 Hz shell log) is refused
%     'Plot'             true = figure: ASD with the bands shaded, and coherence (default false)
%
%   Output out: f (nF x 1, Hz), omega (rad/s), Pxx (nF x C), ASD, PxxPerRadS, units,
%     unitsPSD, fs, N, segmentLength, overlap, window, nSegments, df, enbw (Hz),
%     dofEffective (equivalent chi-square degrees of freedom per PSD value), ci95 ([lo hi]
%     multipliers on Pxx), pairs, Pxy (nF x P complex), coherence (nF x P),
%     coherenceNoiseFloor (EXPECTED gamma^2 for independent channels, 1 / K_eff = 2 / dof),
%     coherence95Threshold (single-bin 95% significance level for independent channels,
%     1 - 0.05^(1/(K_eff - 1)); a bin above it is significantly coherent, one below it is
%     not evidence of independence, only of no detectable common content),
%     channelNames, bands (per band: range, nBins, meanPSD, rms, meanCoherence, note), timing.
%
% Lives in IMU Characterisation/Processing/.

[N, C] = size(X);
if N < C, warning('computeNoisePSD:shape', 'X has fewer rows than columns: channels must be columns.'); end
opt = struct('SegmentLength', [], 'Overlap', 0.5, 'Window', 'hann', 'Detrend', 'mean', 'Pairs', 'default', ...
    'ChannelNames', {{}}, 'Units', 'm/s^2', 'Bands', struct('basin', [0.47 1.35], 'miz', [0.08 0.17]), ...
    'Time', [], 'MaxJitter', 0.05, 'MaxGapPeriods', 1.5, 'Plot', false);
if mod(numel(varargin), 2) ~= 0, error('computeNoisePSD:args', 'Name-value pairs expected.'); end
for i = 1:2:numel(varargin)
    if ~isfield(opt, varargin{i}), error('computeNoisePSD:args', 'Unknown option ''%s''.', varargin{i}); end
    opt.(varargin{i}) = varargin{i + 1};
end
X = double(X);
if any(~isfinite(X(:))), error('computeNoisePSD:nonfinite', 'X contains NaN or Inf (fill or cut gaps first).'); end

% ---- timing: uniform grid required --------------------------------------------------------
timing = struct('checked', false);
if ~isempty(opt.Time)
    t = double(opt.Time(:));
    if numel(t) ~= N, error('computeNoisePSD:time', 'Time must have one entry per row of X.'); end
    dt = diff(t); md = median(dt);
    timing = struct('checked', true, 'fsFromTime', (N - 1) / (t(end) - t(1)), 'medianDt', md, ...
        'jitter', std(dt) / md, 'maxGapPeriods', max(dt) / md);
    if isempty(fs), fs = timing.fsFromTime; end
    if timing.jitter > opt.MaxJitter || timing.maxGapPeriods > opt.MaxGapPeriods
        error('computeNoisePSD:irregular', ['Sampling is not uniform enough for Welch: std(dt)/median(dt) = %.3f ' ...
            '(max %.3f), max dt = %.2f median periods (max %.2f). Use a loss-free record, or cut at the gaps.'], ...
            timing.jitter, opt.MaxJitter, timing.maxGapPeriods, opt.MaxGapPeriods);
    end
    if abs(timing.fsFromTime / fs - 1) > 0.01
        warning('computeNoisePSD:fs', 'fs = %.4g Hz differs from the rate implied by Time (%.4g Hz) by more than 1%%.', ...
            fs, timing.fsFromTime);
    end
end
if isempty(fs), error('computeNoisePSD:fs', 'fs is required when Time is not given.'); end

% ---- segments and window ------------------------------------------------------------------------
L = opt.SegmentLength;
if isempty(L), L = 2^floor(log2(max(8, N / 4))); end
L = min(L, N);
hop = max(1, round(L * (1 - opt.Overlap)));
starts = 1:hop:(N - L + 1);
K = numel(starts);
switch lower(opt.Window)
    case 'hann', w = 0.5 - 0.5 * cos(2 * pi * (0:L - 1).' / L);        % periodic Hann
    case 'rect', w = ones(L, 1);
    otherwise, error('computeNoisePSD:window', 'Window must be ''hann'' or ''rect''.');
end
U = sum(w.^2);
nF = floor(L / 2) + 1;
f = (0:nF - 1).' * fs / L;
scale = 2 * ones(nF, 1) / (fs * U); scale(1) = 1 / (fs * U);
if mod(L, 2) == 0, scale(end) = 1 / (fs * U); end

if ischar(opt.Pairs) && strcmp(opt.Pairs, 'default')
    if C >= 2 && C <= 6, pairs = nchoosek(1:C, 2); else, pairs = zeros(0, 2); end
else
    pairs = opt.Pairs; if isempty(pairs), pairs = zeros(0, 2); end
end
P = size(pairs, 1);

% ---- Welch accumulation -----------------------------------------------------------------------------
Pxx = zeros(nF, C); Pxy = complex(zeros(nF, P));
n = (0:L - 1).';
for k = 1:K
    seg = X(starts(k):starts(k) + L - 1, :);
    switch lower(opt.Detrend)
        case 'linear'
            A = [ones(L, 1), n - mean(n)];
            seg = seg - A * (A \ seg);
        case 'mean'
            seg = seg - mean(seg, 1);
        case 'none'
        otherwise
            error('computeNoisePSD:detrend', 'Detrend must be ''mean'', ''linear'' or ''none''.');
    end
    Z = fft(seg .* w);
    Z = Z(1:nF, :);
    Pxx = Pxx + abs(Z).^2;
    for p = 1:P
        Pxy(:, p) = Pxy(:, p) + Z(:, pairs(p, 1)) .* conj(Z(:, pairs(p, 2)));
    end
end
Pxx = Pxx .* scale / K;
Pxy = Pxy .* scale / K;
coh = zeros(nF, P);
for p = 1:P
    coh(:, p) = abs(Pxy(:, p)).^2 ./ (Pxx(:, pairs(p, 1)) .* Pxx(:, pairs(p, 2)));
end

% ---- statistics: equivalent degrees of freedom (Welch 1967 overlap correction) ------------------------
rho = 0;
for j = 1:K - 1
    lag = j * hop;
    if lag >= L, break; end
    c = (w(1:L - lag).' * w(1 + lag:L)) / U;
    rho = rho + 2 * (1 - j / K) * c^2;
end
dof = 2 * K / (1 + rho);
z = 1.959964;
chiq = @(q) dof * (1 - 2 / (9 * dof) + q * sqrt(2 / (9 * dof)))^3;   % Wilson-Hilferty chi-square quantile
ci95 = [dof / chiq(z), dof / chiq(-z)];                              % S_true in [lo hi] x S_hat
enbw = fs * U / sum(w)^2;
Keff = dof / 2;                                                       % effective independent averages
coh95 = 1 - 0.05^(1 / max(Keff - 1, eps));                           % P(gamma^2 > coh95 | independent) = 5%

% ---- band summaries ------------------------------------------------------------------------------------
bands = struct();
bn = fieldnames(opt.Bands);
df = fs / L;
for b = 1:numel(bn)
    r = opt.Bands.(bn{b});
    in = f >= r(1) & f <= r(2) & f > 0;
    s = struct('range', r, 'nBins', nnz(in), 'meanPSD', NaN(1, C), 'rms', NaN(1, C), ...
        'meanCoherence', NaN(1, P), 'note', '');
    if r(2) > fs / 2
        s.note = 'band extends above Nyquist';
    end
    if nnz(in) == 0
        s.note = sprintf('no bins: resolution df = %.3g Hz is too coarse (need longer segments)', df);
    else
        s.meanPSD = mean(Pxx(in, :), 1);
        s.rms = sqrt(sum(Pxx(in, :), 1) * df);
        if P > 0, s.meanCoherence = mean(coh(in, :), 1); end
        if nnz(in) < 3, s.note = sprintf('only %d bins: coarse; use longer segments', nnz(in)); end
    end
    bands.(bn{b}) = s;
end

names = opt.ChannelNames;
if isempty(names), names = arrayfun(@(c) sprintf('ch%d', c), 1:C, 'UniformOutput', false); end
out = struct('f', f, 'omega', 2 * pi * f, 'Pxx', Pxx, 'ASD', sqrt(Pxx), 'PxxPerRadS', Pxx / (2 * pi), ...
    'units', opt.Units, 'unitsPSD', ['(', opt.Units, ')^2/Hz'], 'fs', fs, 'N', N, 'segmentLength', L, ...
    'overlap', opt.Overlap, 'window', lower(opt.Window), 'detrend', lower(opt.Detrend), 'nSegments', K, ...
    'df', df, 'enbw', enbw, 'dofEffective', dof, 'ci95', ci95, 'pairs', pairs, 'Pxy', Pxy, ...
    'coherence', coh, 'coherenceNoiseFloor', min(1, 2 / dof), 'coherence95Threshold', coh95, ...
    'channelNames', {names}, 'bands', bands, ...
    'timing', timing);

if opt.Plot, plotPSD(out, opt); end
end

function plotPSD(out, opt)
fig = figure('Color', 'w', 'Position', [60 60 1200 440]); %#ok<NASGU>
ax1 = subplot(1, 2, 1); hold(ax1, 'on');
k = out.f > 0;
a = out.ASD(k, :); yl = [min(a(:)) / 2, max(a(:)) * 2];
shadeBands(ax1, opt.Bands, yl);
h = loglog(ax1, out.f(k), out.ASD(k, :), 'LineWidth', 1.3);
set(ax1, 'XScale', 'log', 'YScale', 'log'); ylim(ax1, yl); xlim(ax1, [out.f(2), out.f(end)]); grid(ax1, 'on');
xlabel(ax1, 'f (Hz)'); ylabel(ax1, ['ASD (', out.units, '/\surdHz)']);
title(ax1, sprintf('Welch ASD: L = %d, %d segments, df = %.3g Hz', out.segmentLength, out.nSegments, out.df));
legend(ax1, h, out.channelNames, 'Interpreter', 'none', 'Location', 'southwest');
if ~isempty(out.pairs)
    ax2 = subplot(1, 2, 2); hold(ax2, 'on');
    shadeBands(ax2, opt.Bands, [0 1]);
    lab = cell(1, size(out.pairs, 1));
    for p = 1:size(out.pairs, 1)
        plot(ax2, out.f(k), out.coherence(k, p), 'LineWidth', 1.2);
        lab{p} = [out.channelNames{out.pairs(p, 1)}, ' / ', out.channelNames{out.pairs(p, 2)}];
    end
    plot(ax2, out.f([2 end]), out.coherenceNoiseFloor * [1 1], 'k:');
    plot(ax2, out.f([2 end]), out.coherence95Threshold * [1 1], 'k--');
    set(ax2, 'XScale', 'log'); ylim(ax2, [0 1]); xlim(ax2, [out.f(2), out.f(end)]); grid(ax2, 'on');
    xlabel(ax2, 'f (Hz)'); ylabel(ax2, '\gamma^2');
    title(ax2, 'Coherence (dotted: expected if independent; dashed: 95% significance)');
    legend(ax2, [lab, {'expected (indep.)', '95% threshold'}], 'Interpreter', 'none', 'Location', 'northeast');
end
end

function shadeBands(ax, bands, yl)
bn = fieldnames(bands);
cols = [0.85 0.92 1.0; 1.0 0.92 0.85];
for b = 1:numel(bn)
    r = bands.(bn{b});
    patch(ax, [r(1) r(2) r(2) r(1)], [yl(1) yl(1) yl(2) yl(2)], cols(1 + mod(b - 1, 2), :), ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(ax, r(1), yl(2), [' ', bn{b}], 'VerticalAlignment', 'top', 'FontSize', 8);
end
end