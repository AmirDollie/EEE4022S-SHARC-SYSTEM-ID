%% extractTransmissibilityTester.m
% Unit tests for extractTransmissibility.m (Transmissibility/). No EMM:
% signals come from an analytical modal acceleration model passed through
% stochasticSensorSynthesis.m (SSI/), or from exact delays. About a minute.
%
% Groups
%   1  Orientation: a pure delay y_j(t) = y_r(t - tau) gives T = exp(-i omega tau)
%      (the conjugate orientation would give the opposite phase slope)
%   2  Known modal system, white common input, no noise: T_hat = H_j/H_r in
%      magnitude AND phase; coherence ~ 1
%   3  PSD normalisation: white noise std sigma -> S_rr = sigma^2 dt / pi
%   4  Coloured (JONSWAP-shaped) vs white common input, with sensor noise:
%      same mean, larger variance in low-energy bins
%   5  Noise on the reference only: T_hat biased by S_rr/(S_rr+S_nn) as
%      predicted; the NoisePSD correction removes it
%   6  Noise on a non-reference channel only: no bias
%   7  Invalid-denominator bins flagged and set to NaN; uncorrected kept
%   8  Analytical variance diagnostic vs Monte Carlo (order of magnitude)
%   9  nEffective (Welch overlap correction); frequency grid and Band
%  10  Reference choice and row ordering
%  11  Error handling
%
% Lives in Transmissibility/Unit Tests Transmissibility/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));                     % extractTransmissibility
addpath(fullfile(thisDir, '..', '..', 'SSI'));        % stochasticSensorSynthesis
clearvars -except thisDir
close all, clc

nPass = 0;
if exist('OCTAVE_VERSION', 'builtin') > 0, randn('state', 7); rand('state', 7); else, rng(7); end %#ok<RAND>

%% common test system: 4 channels, 2 well-damped modes, acceleration output
dt = 0.1; N = 8192; L = 512;
dW = 2 * pi / (N * dt);
kS = ceil(3 / dW):floor(9 / dW); wS = kS * dW;                  % synthesis band
sys.wn = [5.0 7.0]; sys.zeta = [0.15 0.20];
sys.Phi = [1.0 0.4; 0.6 -0.8; -0.3 1.0; 0.9 0.5];              % 4 sensors x 2 modes
Hfun = @(w) modalAcc(sys, w);
band = [3.6 8.4];                                                % estimator band, inside the taper
ref = 4;

%% 1  orientation via exact delay
fprintf('=== 1: orientation (pure delay) ===\n');
d = 3; tau = d * dt;
x = randn(1, N + d);
Y = [x(1 + d:end); x(1:end - d)];                                 % row 2 = row 1 delayed by tau
e = extractTransmissibility(Y, dt, 'Reference', 1, 'SegmentLength', L);
Tex = exp(-1i * e.omega * tau);
phErr = max(abs(angle(e.T ./ Tex)));
phWrong = max(abs(angle(conj(e.T) ./ Tex)));
assert(phErr < 0.05, sprintf('1: phase error %.3f rad', phErr));
assert(phWrong > 1, '1: conjugate orientation should fail this test');
assert(max(abs(abs(e.T) - 1)) < 0.05, '1: delay magnitude should be 1');
fprintf('PASS: max phase error %.1e rad (conjugate would be %.2f rad)\n', phErr, phWrong); nPass = nPass + 1;

%% 2  known modal system, white input, noise-free
fprintf('\n=== 2: known modal system, noise-free ===\n');
Y = stochasticSensorSynthesis(wS, Hfun(wS), dt, 'Seed', 11, 'Taper', 'cosine');
e = extractTransmissibility(Y, dt, 'Reference', ref, 'SegmentLength', L, 'Band', band);
Htrue = Hfun(e.omega);
Ttrue = Htrue(e.others, :) ./ Htrue(ref, :);
relErr = abs(e.T - Ttrue) ./ abs(Ttrue);
% Without noise the only error is window leakage, largest where a channel
% is near a response zero (channel 3 near 3.8 rad/s here): coherence dips
% slightly below 1 there even though no noise is present.
assert(median(relErr(:)) < 0.01 && max(relErr(:)) < 0.06, ...
    sprintf('2: T error median %.2e max %.2e', median(relErr(:)), max(relErr(:))));
assert(min(e.coherence(:)) > 0.95 && median(e.coherence(:)) > 0.998, '2: coherence should be ~1 without noise');
assert(all(e.valid), '2: all bins valid without correction');
fprintf('PASS: relative error median %.1e, max %.1e (window leakage); min coherence %.4f\n', ...
    median(relErr(:)), max(relErr(:)), min(e.coherence(:))); nPass = nPass + 1;

%% 3  PSD normalisation
fprintf('\n=== 3: PSD normalisation ===\n');
sig = 2.0;
Y = sig * randn(2, 200000);
e = extractTransmissibility(Y, dt, 'SegmentLength', L);
ratio = mean(e.Srr) / (sig^2 * dt / pi);
assert(abs(ratio - 1) < 0.02, sprintf('3: mean S_rr / (sigma^2 dt / pi) = %.3f', ratio));
fprintf('PASS: mean S_rr / (sigma^2 dt/pi) = %.4f\n', ratio); nPass = nPass + 1;

%% common Monte Carlo helper settings
nR = 60;
Hs = Hfun(wS);
jon = jonswapShape(wS, 5.5, 3.3); jon = jon / max(jon);         % coloured input PSD (peak 1)
sigSig = std(reshape(stochasticSensorSynthesis(wS, Hs, dt, 'Seed', 1), 1, []));
noiseStd = 0.1 * sigSig;                                         % sensor noise, 10% of white-case signal std

%% 4  coloured vs white input
fprintf('\n=== 4: white vs JONSWAP-shaped common input (sensor noise on all channels) ===\n');
[TW, omegaMC] = mcT(nR, wS, Hs, dt, L, band, ref, 1, noiseStd * [1 1 1 1], [], 100);
TJ = mcT(nR, wS, Hs, dt, L, band, ref, jon, noiseStd * [1 1 1 1] * sqrt(mean(jon)), [], 200);
Ht = Hfun(omegaMC); Tt = Ht(setdiff(1:4, ref), :) ./ Ht(ref, :);
biasW = abs(mean(TW, 3) - Tt) ./ abs(Tt);
biasJ = abs(mean(TJ, 3) - Tt) ./ abs(Tt);
varW = mean(abs(TW - mean(TW, 3)).^2, 3) ./ abs(Tt).^2;
varJ = mean(abs(TJ - mean(TJ, 3)).^2, 3) ./ abs(Tt).^2;
jonAt = interp1(wS, jon, omegaMC);
low = jonAt < 0.05; high = jonAt > 0.5;
fprintf('  median relative bias: white %.3f, JONSWAP %.3f (noise biases both through S_rr)\n', ...
    median(biasW(:)), median(biasJ(:)));
fprintf('  median relative variance, low-energy bins: white %.2e, JONSWAP %.2e\n', ...
    median(reshape(varW(:, low), 1, [])), median(reshape(varJ(:, low), 1, [])));
fprintf('  median relative variance, high-energy bins: white %.2e, JONSWAP %.2e\n', ...
    median(reshape(varW(:, high), 1, [])), median(reshape(varJ(:, high), 1, [])));
assert(any(low) && any(high), '4: band does not contain both low- and high-energy bins');
assert(median(reshape(varJ(:, low), 1, [])) > 5 * median(reshape(varW(:, low), 1, [])), ...
    '4: coloured input should raise variance in low-energy bins');
fprintf('PASS\n'); nPass = nPass + 1;

%% 5  reference-noise bias and correction
fprintf('\n=== 5: noise on the reference only ===\n');
sn = 2.0 * sigSig;                                % white over [0, pi/dt]: in-band SNR stays moderate
Snn = sn^2 * dt / pi;
[Tu, ~, Tc, SrrMean] = mcT(nR, wS, Hs, dt, L, band, ref, 1, [0 0 0 sn], Snn, 300);
Srr0 = SrrMean - Snn;                                            % noise-free S_rr estimate
factorPred = Srr0 ./ (Srr0 + Snn);
factorObs = mean(abs(mean(Tu, 3) ./ Tt), 1);
okBins = all(all(isfinite(Tc), 3), 1);                          % bins valid in every realisation
biasC = abs(mean(Tc(:, okBins, :), 3) - Tt(:, okBins)) ./ abs(Tt(:, okBins));
fprintf('  uncorrected |mean T_hat| / |T|: observed median %.3f, predicted median %.3f\n', ...
    median(factorObs), median(factorPred));
fprintf('  corrected: %d of %d bins valid in every realisation; median relative bias %.3f, max %.3f\n', ...
    sum(okBins), numel(okBins), median(biasC(:)), max(biasC(:)));
assert(sum(okBins) > 0.8 * numel(okBins), '5: too many bins flagged');
assert(abs(median(factorObs) - median(factorPred)) < 0.03, '5: bias factor does not match S_rr/(S_rr+S_nn)');
assert(median(factorPred) < 0.9, '5: test not sensitive (noise too small)');
biasU = abs(mean(Tu(:, okBins, :), 3) - Tt(:, okBins)) ./ abs(Tt(:, okBins));
fprintf('  median relative bias: uncorrected %.3f -> corrected %.3f (residual = finite-sample ratio bias\n', ...
    median(biasU(:)), median(biasC(:)));
fprintf('  at low reference SNR; it shrinks with record length: a G2 measurement, not a bug)\n');
assert(median(biasC(:)) < 0.2 * median(biasU(:)) && median(biasC(:)) < 0.05, '5: correction did not remove most of the bias');
fprintf('PASS\n'); nPass = nPass + 1;

%% 6  noise on a non-reference channel only
fprintf('\n=== 6: noise on a non-reference channel only ===\n');
T6 = mcT(nR, wS, Hs, dt, L, band, ref, 1, [0.5 * sigSig 0 0 0], [], 400);
bias6 = abs(mean(T6(1, :, :), 3) - Tt(1, :)) ./ abs(Tt(1, :));
fprintf('  channel 1 (noisy) median relative bias %.3f\n', median(bias6));
assert(median(bias6) < 0.02, '6: non-reference noise should not bias T');
fprintf('PASS\n'); nPass = nPass + 1;

%% 7  invalid-denominator flag
fprintf('\n=== 7: invalid denominator ===\n');
Y = stochasticSensorSynthesis(wS, Hs, dt, 'Seed', 5);
e0 = extractTransmissibility(Y, dt, 'Reference', ref, 'SegmentLength', L, 'Band', band);
big = 0.95 * e0.Srr; big(1:2:end) = 0;                           % every other bin ~ fully "noise"
e = extractTransmissibility(Y, dt, 'Reference', ref, 'SegmentLength', L, 'Band', band, 'NoisePSD', big);
assert(isequal(e.valid, big == 0), '7: wrong bins flagged');
assert(all(all(isnan(e.T(:, ~e.valid)))) && all(all(isfinite(e.T(:, e.valid)))), '7: NaN pattern wrong');
assert(isequal(e.Tuncorrected, e0.T), '7: uncorrected estimate must be unchanged');
e = extractTransmissibility(Y, dt, 'Reference', ref, 'SegmentLength', L, 'Band', band, ...
    'NoisePSD', @(w) 0 * w);
assert(all(e.valid) && max(abs(e.T(:) - e0.T(:))) < 1e-12 * max(abs(e0.T(:))), '7: zero noise handle should change nothing');
fprintf('PASS\n'); nPass = nPass + 1;

%% 8  analytical variance diagnostic vs Monte Carlo
fprintf('\n=== 8: analytical variance diagnostic vs Monte Carlo ===\n');
[T8, ~, ~, ~, vDiag] = mcT(nR, wS, Hs, dt, L, band, ref, 1, noiseStd * [1 1 1 1], [], 500);
vMC = mean(abs(T8 - mean(T8, 3)).^2, 3);
rat = vMC ./ vDiag;
fprintf('  MC variance / diagnostic: median %.2f (IQR %.2f-%.2f)\n', median(rat(:)), ...
    quantileSimple(rat(:), 0.25), quantileSimple(rat(:), 0.75));
assert(median(rat(:)) > 0.5 && median(rat(:)) < 2, '8: diagnostic variance off by more than 2x');
fprintf('PASS (diagnostic only; S4 estimates the real covariance)\n'); nPass = nPass + 1;

%% 9  nEffective, grid, Band
fprintf('\n=== 9: nEffective, frequency grid, Band ===\n');
Y = randn(2, 20 * L);
e50 = extractTransmissibility(Y, dt, 'SegmentLength', L, 'Overlap', 0.5);
e0v = extractTransmissibility(Y, dt, 'SegmentLength', L, 'Overlap', 0);
K = e50.nSegments;
nPred = K / (1 + 2 * (1 - 1 / K) * 0.1667^2);
assert(abs(e50.nEffective - nPred) / nPred < 0.01, '9: Hann 50% nEffective wrong');
assert(e0v.nEffective == e0v.nSegments, '9: no overlap should give nEffective = nSegments');
assert(max(abs(e50.omega - 2 * pi * (1:L / 2 - 1) / (L * dt))) < 1e-12, '9: frequency grid wrong');
eb = extractTransmissibility(Y, dt, 'SegmentLength', L, 'Band', [3 8.5]);
assert(all(eb.omega >= 3 & eb.omega <= 8.5) && numel(eb.omega) == sum(e50.omega >= 3 & e50.omega <= 8.5), '9: Band subset wrong');
assert(isequal(eb.T, e50.T(:, e50.omega >= 3 & e50.omega <= 8.5)), '9: Band must only select bins');
fprintf('PASS: K = %d, nEffective = %.2f; %d bins in [3, 8.5] rad/s at L = %d, dt = %.1f\n', ...
    K, e50.nEffective, numel(eb.omega), L, dt); nPass = nPass + 1;

%% 10  reference choice and row ordering
fprintf('\n=== 10: reference and row ordering ===\n');
Y = stochasticSensorSynthesis(wS, Hs, dt, 'Seed', 21);
for rr = 1:4
    e = extractTransmissibility(Y, dt, 'Reference', rr, 'SegmentLength', L, 'Band', band);
    assert(isequal(e.others, setdiff(1:4, rr)) && e.reference == rr, '10: others/reference wrong');
    H10 = Hfun(e.omega);
    err = abs(e.T - H10(e.others, :) ./ H10(rr, :)) ./ abs(H10(e.others, :) ./ H10(rr, :));
    assert(median(err(:)) < 0.01, sprintf('10: reference %d wrong', rr));
end
fprintf('PASS\n'); nPass = nPass + 1;

%% 11  errors
fprintf('\n=== 11: error handling ===\n');
Y = randn(3, 2000);
cases = {
    'badY',        {randn(1, 2000), dt}
    'badY',        {[Y(:, 1:end - 1), NaN(3, 1)], dt}
    'badY',        {Y + 1i, dt}
    'badDt',       {Y, 0}
    'badOption',   {Y, dt, 'Reference', 4}
    'badOption',   {Y, dt, 'SegmentLength', 511}
    'badOption',   {Y, dt, 'Overlap', 1}
    'badOption',   {Y, dt, 'Band', [5 3]}
    'badOption',   {Y, dt, 'MinDenominatorFraction', 1}
    'badOption',   {Y, dt, 'Nope', 1}
    'tooShort',    {Y(:, 1:100), dt}
    'emptyBand',   {Y, dt, 'Band', [40 50]}
    'badNoisePSD', {Y, dt, 'NoisePSD', [1 2 3]}
    'badNoisePSD', {Y, dt, 'NoisePSD', -1}
    };
for c = 1:size(cases, 1)
    id = ['extractTransmissibility:', cases{c, 1}];
    try
        extractTransmissibility(cases{c, 2}{:});
        error('tester:noError', '11.%d: expected %s but no error was thrown', c, id);
    catch err
        assert(strcmp(err.identifier, id), sprintf('11.%d: expected %s, got %s (%s)', c, id, err.identifier, err.message));
    end
end
fprintf('PASS: %d error cases\n', size(cases, 1)); nPass = nPass + 1;

fprintf('\nAll %d groups passed.\n', nPass);

%% ================================================================================================
function H = modalAcc(sys, w)
% Acceleration FRF of a modal system under exp(+i w t): sum_m Phi(:,m) (-w^2)/(wn^2 - w^2 + 2i zeta wn w).
H = zeros(size(sys.Phi, 1), numel(w));
for m = 1:numel(sys.wn)
    H = H + sys.Phi(:, m) * (-w.^2 ./ (sys.wn(m)^2 - w.^2 + 2i * sys.zeta(m) * sys.wn(m) * w));
end
end

function S = jonswapShape(w, wp, gJ)
sig = 0.07 * ones(size(w)); sig(w > wp) = 0.09;
S = w.^-5 .* exp(-1.25 * (wp ./ w).^4) .* gJ.^exp(-(w - wp).^2 ./ (2 * sig.^2 * wp^2));
end

function [Tu, omega, Tc, SrrMean, vDiag] = mcT(nR, wS, Hs, dt, L, band, ref, Su, noiseStd, Snn, seed0)
% Monte Carlo of T_hat: common input with PSD Su, independent white sensor
% noise with per-channel std noiseStd, optional NoisePSD correction Snn.
Tu = []; Tc = []; SrrAcc = 0; vAcc = 0;
for r = 1:nR
    Y = stochasticSensorSynthesis(wS, Hs, dt, 'Seed', seed0 * 1000 + r, 'InputPSD', Su, 'Taper', 'cosine');
    Y = Y + noiseStd(:) .* randn(size(Y));
    e = extractTransmissibility(Y, dt, 'Reference', ref, 'SegmentLength', L, 'Band', band);
    Tu(:, :, r) = e.T; %#ok<AGROW>
    SrrAcc = SrrAcc + e.Srr;
    vAcc = vAcc + e.varTdiag;
    if ~isempty(Snn)
        ec = extractTransmissibility(Y, dt, 'Reference', ref, 'SegmentLength', L, 'Band', band, 'NoisePSD', Snn);
        Tc(:, :, r) = ec.T; %#ok<AGROW>
    end
end
omega = e.omega;
SrrMean = SrrAcc / nR;
vDiag = vAcc / nR;
end

function q = quantileSimple(x, p)
x = sort(x(~isnan(x)));
q = x(max(1, round(p * numel(x))));
end