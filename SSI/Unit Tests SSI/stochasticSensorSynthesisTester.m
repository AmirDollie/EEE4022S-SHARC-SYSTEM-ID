%% stochasticSensorSynthesisTester.m
% Unit tests for stochasticSensorSynthesis.m (SSI/). No EMM solves: every
% H here is a small analytical or random complex matrix, so the whole
% tester runs in seconds.
%
% Groups
%   1  Shapes, time vector, info fields, zero mean
%   2  Exact regression against the direct tone sum y = Re{sum A_k e^{iwt}}
%   3  Record is exactly one period (circularity)
%   4  Seeds: reproducible, distinct, global RNG untouched
%   5  FFT recovers A_k exactly (no leakage) and Hermitian synthesis is real
%   6  Statistics: cross-spectrum S_u w^2 H H^H, common input (coherence 1),
%      covariance, E|xi|^2 = 2 in both amplitude modes
%   7  Contrast: independent per-sensor inputs would NOT give coherence 1
%   8  Taper and InputPSD are applied as documented (amplitude vs PSD)
%   9  Analytical damped oscillator FRF: PSD shape recovered
%  10  Error handling
%
% Lives in SSI/Unit Tests SSI/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clearvars -except thisDir
close all, clc

nPass = 0;
maxAll = @(X) max(X(:));   % portable max over all elements

%% Common grid used by most groups
N  = 1024;
dt = 0.1;
dW = 2*pi/(N*dt);
k  = 30:60;                               % bins -> omega ~ 1.84..3.68 rad/s
omega = k*dW;
nW = numel(omega);
nS = 3;
Hrand = (1 + 0.5*randn(nS, nW)) .* exp(1i*2*pi*rand(nS, nW));

%% 1  Shapes, time vector, info fields, zero mean
fprintf('=== 1: shapes, time vector, info, zero mean ===\n');
[y, t, info] = stochasticSensorSynthesis(omega, Hrand, dt, 'NumRealisations', 4, 'Seed', 1);
assert(isequal(size(y), [nS, N, 4]), '1: wrong Y size');
assert(isreal(y), '1: Y must be real');
assert(isequal(size(t), [1, N]) && max(abs(t - (0:N-1)*dt)) < 1e-12, '1: wrong T');
assert(info.N == N && abs(info.dOmega - dW) < 1e-14 && isequal(info.binIndex, k), '1: grid info wrong');
assert(abs(info.recordLength - 2*pi/dW) < 1e-9, '1: record length should be 2*pi/dOmega');
assert(abs(info.omegaNyquist - pi/dt) < 1e-12, '1: Nyquist wrong');
assert(isequal(size(info.inputCoefficients), [nW, 4]), '1: xi size wrong');
assert(isequal(size(info.expectedCrossSpectrum), [nS, nS, nW]), '1: cross-spectrum size wrong');
assert(maxAll(abs(mean(y, 2))) < 1e-10 * max(abs(y(:))), '1: DC must be exactly zero');
% single sensor given as nW x 1
y1 = stochasticSensorSynthesis(omega, Hrand(1, :).', dt, 'Seed', 1);
assert(isequal(size(y1), [1, N]), '1: single-sensor column H not accepted');
fprintf('PASS\n'); nPass = nPass + 1;

%% 2  Exact regression against the direct tone sum
fprintf('\n=== 2: regression against direct sum ===\n');
xi = (randn(nW, 2) + 1i*randn(nW, 2));
Su = 0.7;
[y, t] = stochasticSensorSynthesis(omega, Hrand, dt, 'NumRealisations', 2, ...
    'InputCoefficients', xi, 'InputPSD', Su);
worst = 0;
for r = 1:2
    A = Hrand .* (sqrt(Su*dW) * xi(:, r).');           % nS x nW
    yDirect = real(A * exp(1i * omega(:) * t));         % nS x N
    worst = max(worst, maxAll(abs(y(:, :, r) - yDirect)) / max(abs(yDirect(:))));
end
assert(worst < 1e-10, sprintf('2: FFT synthesis differs from direct sum (rel %.2e)', worst));
fprintf('PASS: max relative difference %.2e\n', worst); nPass = nPass + 1;

%% 3  Record is exactly one period
fprintf('\n=== 3: circularity (period = 2*pi/dOmega) ===\n');
A = Hrand .* (sqrt(Su*dW) * xi(:, 1).');
yShift = real(A * exp(1i * omega(:) * (t + info.recordLength)));
d = maxAll(abs(y(:, :, 1) - yShift)) / max(abs(yShift(:)));
assert(d < 1e-9, '3: record is not periodic with period 2*pi/dOmega');
fprintf('PASS: y(t) = y(t + T_rec) to %.2e (record is ONE period of a tone set)\n', d); nPass = nPass + 1;

%% 4  Seeds
fprintf('\n=== 4: seeds ===\n');
ya = stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 42, 'NumRealisations', 3);
yb = stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 42, 'NumRealisations', 3);
yc = stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 43, 'NumRealisations', 3);
assert(isequal(ya, yb), '4: same seed must reproduce exactly');
assert(~isequal(ya, yc), '4: different seeds must differ');
assert(~isequal(ya(:, :, 1), ya(:, :, 2)), '4: realisations must differ');
% global RNG untouched by a seeded call
if exist('OCTAVE_VERSION', 'builtin') > 0
    s0 = randn('state'); %#ok<RAND>
    stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 7);
    assert(isequal(s0, randn('state')), '4: seeded call changed the global RNG state'); %#ok<RAND>
else
    s0 = rng;
    stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 7);
    s1 = rng;
    assert(isequal(s0.State, s1.State), '4: seeded call changed the global RNG state');
end
fprintf('PASS\n'); nPass = nPass + 1;

%% 5  FFT recovers A_k exactly
fprintf('\n=== 5: no leakage, exact bin recovery ===\n');
[y, ~, info] = stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 3, 'Taper', 'cosine');
Y = fft(y, [], 2);                                       % nS x N
Arec = Y(:, k + 1) * 2 / N;
Atrue = Hrand .* (info.taper .* sqrt(info.inputPSD * dW) .* info.inputCoefficients(:, 1).');
errA = maxAll(abs(Arec - Atrue)) / max(abs(Atrue(:)));
outBand = setdiff(1:floor(N/2), k + 1);
leak = maxAll(abs(Y(:, outBand))) / max(abs(Y(:)));
assert(errA < 1e-10, '5: in-band bins not recovered exactly');
assert(leak < 1e-10, '5: energy found outside the excited bins');
fprintf('PASS: bin error %.1e, out-of-band %.1e\n', errA, leak); nPass = nPass + 1;

%% 6  Statistics
fprintf('\n=== 6: cross-spectrum, coherence, covariance, E|xi|^2 ===\n');
nR = 4000;
SuVec = linspace(0.5, 2, nW);
[y, ~, info] = stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 11, ...
    'NumRealisations', nR, 'InputPSD', SuVec);
Chat = zeros(nS, nS, nW);
for r = 1:nR
    Yr = fft(y(:, :, r), [], 2);
    Ar = Yr(:, k + 1) * 2 / N;
    for q = 1:nW
        Chat(:, :, q) = Chat(:, :, q) + Ar(:, q) * Ar(:, q)';
    end
end
Chat = Chat / (nR * 2 * dW);                            % E[A A^H] = 2 dW C
C = info.expectedCrossSpectrum;
relC = norm(Chat(:) - C(:)) / norm(C(:));
% common input => each averaged matrix is exactly rank 1 => coherence 1
coh = zeros(nS, nS, nW);
for q = 1:nW
    Dq = sqrt(real(diag(Chat(:, :, q))));
    coh(:, :, q) = abs(Chat(:, :, q)).^2 ./ (Dq * Dq').^2;
end
cohErr = max(abs(coh(:) - 1));
% covariance over time and realisations
Yall = reshape(y, nS, []);
Rhat = (Yall * Yall.') / size(Yall, 2);
relR = norm(Rhat - info.expectedCovariance, 'fro') / norm(info.expectedCovariance, 'fro');
% per-bin sampling error ~ 1/sqrt(nR) = 0.016, so 0.05 is ~3 sigma
assert(relC < 0.05, sprintf('6: cross-spectrum off by %.3f', relC));
assert(cohErr < 1e-10, '6: coherence is not 1 (input is not shared across sensors)');
assert(relR < 0.03, sprintf('6: covariance off by %.3f', relR));
% E|xi|^2 = 2 in both modes; randomPhase has deterministic |xi|
[~, ~, iG] = stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 5, 'NumRealisations', 2000);
[~, ~, iP] = stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 5, 'NumRealisations', 2000, ...
    'AmplitudeMode', 'randomPhase');
mG = mean(abs(iG.inputCoefficients(:)).^2);
assert(abs(mG - 2) < 0.03, '6: gaussian E|xi|^2 should be 2');
assert(max(abs(abs(iP.inputCoefficients(:)) - sqrt(2))) < 1e-12, '6: randomPhase |xi| must be sqrt(2)');
ph = angle(iP.inputCoefficients(:));
assert(abs(mean(exp(1i*ph))) < 0.02, '6: randomPhase phases not uniform');
fprintf('PASS: cross-spectrum rel err %.3f, max|coherence-1| %.1e, covariance rel err %.3f, E|xi|^2 %.3f\n', ...
    relC, cohErr, relR, mG); nPass = nPass + 1;

%% 7  Contrast: independent per-sensor inputs lose coherence
fprintf('\n=== 7: contrast (independent inputs per sensor) ===\n');
nR7 = 400;
Cind = zeros(nS, nS);
q7 = 10;
for r = 1:nR7
    a = Hrand(:, q7) .* (randn(nS, 1) + 1i*randn(nS, 1));   % WRONG construction on purpose
    Cind = Cind + a * a';
end
Dq = sqrt(real(diag(Cind)));
cohInd = abs(Cind(1, 2))^2 / (Dq(1) * Dq(2))^2;
assert(cohInd < 0.2, '7: contrast case unexpectedly coherent');
fprintf('PASS: independent-input coherence %.3f (vs 1 for the shared input)\n', cohInd); nPass = nPass + 1;

%% 8  Taper and InputPSD
fprintf('\n=== 8: taper (amplitude) and InputPSD ===\n');
[~, ~, i0] = stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 2, 'Taper', 'none');
assert(all(i0.taper == 1), '8: none should be all ones');
[~, ~, iC] = stochasticSensorSynthesis(omega, Hrand, dt, 'Seed', 2, 'Taper', 'cosine', 'TaperFraction', 0.2);
nRamp = round(0.2 * nW);
assert(all(iC.taper(nRamp+1:end-nRamp) == 1), '8: taper interior must be 1');
assert(iC.taper(1) < 0.2 && iC.taper(end) < 0.2 && all(diff(iC.taper(1:nRamp+1)) > 0), '8: ramp shape wrong');
assert(max(abs(iC.taper - fliplr(iC.taper))) < 1e-15, '8: taper not symmetric');
assert(max(abs(iC.psdTaper - iC.taper.^2)) == 0, '8: PSD taper must be the amplitude taper squared');
assert(maxAll(abs(iC.expectedPSD - iC.taper.^2 .* abs(Hrand).^2)) < 1e-14, '8: expected PSD wrong');
wUser = linspace(0.1, 1, nW);
[y8, ~, iU] = stochasticSensorSynthesis(omega, Hrand, dt, 'Taper', wUser, 'InputPSD', 3, ...
    'AmplitudeMode', 'randomPhase', 'Seed', 9);
Y8 = fft(y8, [], 2);
ampRec = abs(Y8(:, k + 1)) * 2 / N;
ampExp = wUser .* sqrt(3 * dW * 2) .* abs(Hrand);        % |xi| = sqrt(2)
assert(maxAll(abs(ampRec - ampExp)) / max(ampExp(:)) < 1e-10, '8: amplitude not w*sqrt(Su*dW)*|H|*|xi|');
assert(isequal(iU.taper, wUser), '8: user taper not honoured');
fprintf('PASS\n'); nPass = nPass + 1;

%% 9  Damped oscillator: acceleration PSD shape
fprintf('\n=== 9: analytical damped oscillator PSD ===\n');
N9 = 8192; dt9 = 0.1; dW9 = 2*pi/(N9*dt9);
k9 = ceil(3/dW9):floor(9/dW9);
w9 = k9 * dW9;
wn = 6; zeta = 0.05;
Ha = -w9.^2 ./ (wn^2 - w9.^2 + 2i*zeta*wn*w9);          % acceleration FRF
[y9, ~, i9] = stochasticSensorSynthesis(w9, Ha, dt9, 'Seed', 21, 'NumRealisations', 200);
P = zeros(1, numel(k9));
for r = 1:200
    Y9 = fft(y9(1, :, r));
    P = P + abs(Y9(k9 + 1) * 2 / N9).^2;
end
P = P / (200 * 2 * dW9);                                 % one-sided PSD estimate
[~, iPk] = max(movmean(P, 9));
[~, iTh] = max(i9.expectedPSD);
nBW = 2*zeta*wn / dW9;
relP = norm(movmean(P, 21) - movmean(i9.expectedPSD, 21)) / norm(i9.expectedPSD);
assert(abs(w9(iPk) - w9(iTh)) < 0.1, '9: PSD peak misplaced');
assert(relP < 0.05, sprintf('9: PSD shape off by %.3f', relP));
fprintf('PASS: peak %.3f vs %.3f rad/s, smoothed PSD rel err %.3f, N_BW = %.0f bins\n', ...
    w9(iPk), w9(iTh), relP, nBW); nPass = nPass + 1;

%% 10  Error handling
fprintf('\n=== 10: error handling ===\n');
f = @(varargin) stochasticSensorSynthesis(varargin{:});
cases = {
    'badOmega',             {[-1 1], Hrand(:, 1:2), dt}
    'badOmega',             {omega(1), Hrand(:, 1), dt}
    'badOmega',             {fliplr(omega), Hrand, dt}
    'badDt',                {omega, Hrand, 0}
    'badH',                 {omega, [Hrand(:, 1:end-1), NaN(nS, 1)], dt}
    'sizeMismatch',         {omega, Hrand(:, 1:end-1), dt}
    'nonUniformGrid',       {[omega(1:end-1), omega(end) + 0.3*dW], Hrand, dt}
    'gridNotFFTCompatible', {omega, Hrand, dt * 1.01}
    'offGridFrequency',     {omega + 0.5*dW, Hrand, dt}
    'aboveNyquist',         {(N/2 - 5:N/2) * dW, ones(1, 6), dt}
    'badOption',            {omega, Hrand, dt, 'Nope', 1}
    'badOption',            {omega, Hrand, dt, 'Seed'}
    'badOption',            {omega, Hrand, dt, 'NumRealisations', 0}
    'badOption',            {omega, Hrand, dt, 'Seed', -1}
    'badOption',            {omega, Hrand, dt, 'AmplitudeMode', 'uniform'}
    'badTaper',             {omega, Hrand, dt, 'Taper', 'hann'}
    'badTaper',             {omega, Hrand, dt, 'Taper', ones(1, nW - 1)}
    'badTaper',             {omega, Hrand, dt, 'Taper', 'cosine', 'TaperFraction', 0.6}
    'badPSD',               {omega, Hrand, dt, 'InputPSD', -1}
    'badPSD',               {omega, Hrand, dt, 'InputPSD', ones(1, 3)}
    'badCoefficients',      {omega, Hrand, dt, 'InputCoefficients', ones(nW, 2)}
    };
for c = 1:size(cases, 1)
    id = ['stochasticSensorSynthesis:', cases{c, 1}];
    try
        f(cases{c, 2}{:});
        error('tester:noError', '10.%d: expected %s but no error was thrown', c, id);
    catch err
        assert(strcmp(err.identifier, id), sprintf('10.%d: expected %s, got %s (%s)', ...
            c, id, err.identifier, err.message));
    end
end
fprintf('PASS: %d error cases\n', size(cases, 1)); nPass = nPass + 1;

fprintf('\nAll %d groups passed.\n', nPass);