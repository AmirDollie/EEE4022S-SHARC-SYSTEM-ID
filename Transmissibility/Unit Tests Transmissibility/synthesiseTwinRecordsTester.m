%% synthesiseTwinRecordsTester.m
% Unit tests for synthesiseTwinRecords.m (Transmissibility/). Groups 1-5 use
% an analytical modal FRF (no EMM); group 6 uses the EMM on a small band
% (7 solves, about 10 s in MATLAB) to test the disk cache.
%
% Groups
%   1  Node FRF from a handle: nodes span the band, exact at nodes,
%      interpolation check small, prepared twin reused from memory
%   2  Fine grid and exact single-input structure: record bins sit on the
%      FFT grid inside the band, and fft(Y_j)/fft(Y_r) = H_j/H_r bin by bin
%   3  Input spectrum: JONSWAP Hs normalisation and peak, white level,
%      custom handle, edge taper shape, HsInBand
%   4  Statistics: ensemble variance vs expected covariance; Welch T
%      (extractTransmissibility) vs the interpolated truth twin.Hacc
%   5  Noise: absolute and relative levels, noise PSD convention,
%      common random numbers (signal unchanged by noise), seeds, global RNG
%   6  EMM node FRF and disk cache: computed / loaded / key mismatch,
%      agreement with transmissibilityFeatureVector (S1) at the nodes
%   7  Error handling
%
% Lives in Transmissibility/Unit Tests Transmissibility/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));                                   % synthesiseTwinRecords, extractTransmissibility
addpath(fullfile(thisDir, '..', '..', 'SSI'));
addpath(fullfile(thisDir, '..', '..', 'Inverse Mapping'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
clearvars -except thisDir
close all, clc

nPass = 0;
tStart = tic;
maxAll = @(X) max(abs(X(:)));

% analytical 4-sensor, 2-mode acceleration FRF (as in extractTransmissibilityTester)
sys.wn = [5.0 7.0]; sys.zeta = [0.15 0.20];
sys.Phi = [1.0 0.4; 0.6 -0.8; -0.3 1.0; 0.9 0.5];
Hfun = @(w) modalAcc(sys, w);
sensorsToy = [0.1 0; 0.1 pi; 0.2 pi; 0.3 0];                         % recorded only
spec = struct('frf', Hfun, 'sensors', sensorsToy, 'band', [3 9], 'nodeSpacing', 0.05);
dt = 0.1; ref = 4;

%% 1  node FRF from a handle
fprintf('=== 1: node FRF from a handle ===\n');
twin = synthesiseTwinRecords(spec);
wN = twin.omegaNodes;
assert(abs(wN(1) - 3) < 1e-12 && abs(wN(end) - 9) < 1e-12, '1: nodes must include both band edges');
assert(max(diff(wN)) <= 0.05 + 1e-12, '1: node spacing above target');
assert(maxAll(twin.Hacc(wN) - Hfun(wN)) < 1e-12 * maxAll(Hfun(wN)), '1: interpolant not exact at nodes');
wm = 0.5 * (wN(1:end - 1) + wN(2:end));
eH = maxAll(twin.Hacc(wm) - Hfun(wm)) / maxAll(Hfun(wm));
assert(eH < 1e-4, sprintf('1: midpoint interpolation error %.2e', eH));
assert(twin.interpCheck.maxRelErrH < 1e-4 && twin.interpCheck.maxRelErrT < 1e-4, '1: interpCheck too large');
assert(strcmp(twin.cacheStatus, 'computed') && strcmp(twin.source, 'handle'), '1: status / source');
[~, ~, i1, twin2] = synthesiseTwinRecords(twin, 1024, dt, 'Seed', 1);
assert(strcmp(i1.cacheStatus, 'memory') && strcmp(twin2.cacheStatus, 'memory'), '1: prepared twin not reused');
tDef = synthesiseTwinRecords(struct('frf', Hfun, 'sensors', sensorsToy));   % default band = validated EMM band
bExp = sqrt([1.7 13.85] * 9.81 / 1.88);
assert(max(abs(tDef.band - bExp)) < 1e-6 && all(1.88 * tDef.band.^2 / 9.81 >= 1.7) ...
    && all(1.88 * tDef.band.^2 / 9.81 <= 13.85) && abs(tDef.nodeSpacing - 0.1) < eps, '1: default band / spacing');
fprintf('PASS: %d nodes, midpoint error %.1e (H), check %.1e (T)\n', numel(wN), eH, twin.interpCheck.maxRelErrT);
nPass = nPass + 1;

%% 2  fine grid and exact single-input structure
fprintf('\n=== 2: fine grid and exact bin ratios ===\n');
N = 4096;
[Y, t, info] = synthesiseTwinRecords(twin, N, dt, 'Seed', 2, 'Spectrum', 'white');
dW = 2 * pi / (N * dt);
assert(isequal(size(Y), [4, N]) && abs(t(2) - dt) < 1e-15 && abs(info.dOmega - dW) < 1e-15, '2: sizes / grid');
assert(all(info.omega >= 3 & info.omega <= 9) && maxAll(info.omega - info.binIndex * dW) < 1e-12, '2: off-grid or out-of-band bins');
assert(info.binIndex(1) == ceil(3 / dW) && info.binIndex(end) == floor(9 / dW), '2: band bins missing');
assert(maxAll(info.H - twin.Hacc(info.omega)) == 0, '2: fine-grid H is not twin.Hacc');
X = fft(Y, [], 2);
live = info.taper > 1e-3;
kk = info.binIndex(live) + 1;
Tbin = X(1:3, kk) ./ X(ref, kk);
Ttrue = info.H(1:3, live) ./ info.H(ref, live);
e2 = maxAll(Tbin - Ttrue) / maxAll(Ttrue);
assert(e2 < 1e-9, sprintf('2: bin ratio differs from H_j/H_r (%.2e)', e2));
outBins = setdiff(1:N / 2 - 1, info.binIndex);
assert(maxAll(X(:, outBins + 1)) < 1e-9 * maxAll(X), '2: energy outside the band');
fprintf('PASS: %d bins, per-bin ratio error %.1e\n', numel(info.omega), e2); nPass = nPass + 1;

%% 3  input spectrum
fprintf('\n=== 3: spectrum and taper ===\n');
[~, ~, iJ] = synthesiseTwinRecords(twin, 8192, dt, 'Seed', 3, 'Hs', 0.08, 'PeakFrequency', 4.5, 'TaperWidth', 0.3);
wq = linspace(0.5, 100, 400000);
HsNum = 4 * sqrt(trapz(wq, iJ.spectrum(wq)));
assert(abs(HsNum / 0.08 - 1) < 2e-3, sprintf('3: JONSWAP Hs %.5f vs 0.08', HsNum));
[~, ip] = max(iJ.spectrum(wq));
assert(abs(wq(ip) - 4.5) < 0.01, '3: JONSWAP peak not at PeakFrequency');
assert(iJ.HsInBand < 0.08 && iJ.HsInBand > 0.5 * 0.08, '3: HsInBand implausible');
tp = iJ.taper; w3 = iJ.omega;
assert(all(tp >= 0 & tp <= 1) && all(tp(w3 >= 3.3 & w3 <= 8.7) == 1), '3: taper not 1 in the interior');
assert(tp(1) < 1e-2 && tp(end) < 1e-2, '3: taper not ~0 at band edges');
lo = w3 < 3.3; hi = w3 > 8.7;
assert(all(diff(tp(lo)) > 0) && all(diff(tp(hi)) < 0), '3: ramps not monotone');
[~, ~, iW] = synthesiseTwinRecords(twin, 4096, dt, 'Spectrum', 'white', 'WhiteLevel', 2e-3, 'Seed', 3);
assert(all(abs(iW.inputPSD - 2e-3) < 1e-15), '3: white level');
[~, ~, iC] = synthesiseTwinRecords(twin, 4096, dt, 'Spectrum', @(w) 1e-3 ./ w, 'Seed', 3);
assert(maxAll(iC.inputPSD - 1e-3 ./ iC.omega) < 1e-15 && strcmp(iC.spectrumType, 'custom'), '3: custom handle');
fprintf('PASS: Hs %.5f (target 0.08), HsInBand %.4f\n', HsNum, iJ.HsInBand); nPass = nPass + 1;

%% 4  statistics
fprintf('\n=== 4: ensemble statistics and Welch transmissibility ===\n');
nR = 200;
[Ye, ~, ie] = synthesiseTwinRecords(twin, 2048, dt, 'Seed', 4, 'NumRealisations', nR, 'Spectrum', 'white', 'WhiteLevel', 1e-3);
v = squeeze(mean(Ye.^2, 2));                                         % nS x nR sample variances
rel = mean(v, 2) ./ diag(ie.expectedCovariance);
assert(all(abs(rel - 1) < 0.08), sprintf('4: ensemble variance ratio %s', mat2str(rel.', 3)));
assert(maxAll(ie.signalStd.^2 - diag(ie.expectedCovariance).') < 1e-15, '4: signalStd');
[YL, ~, iL] = synthesiseTwinRecords(twin, 65536, dt, 'Seed', 5, 'Spectrum', 'white');
est = extractTransmissibility(YL, dt, 'Reference', ref, 'SegmentLength', 512, 'Band', [3.6 8.4]);
Tt = twin.Hacc(est.omega); Tt = Tt(1:3, :) ./ Tt(ref, :);
eW = abs(est.T - Tt) ./ abs(Tt);
assert(median(eW(:)) < 0.01 && max(eW(:)) < 0.06, sprintf('4: Welch T error median %.3f max %.3f', median(eW(:)), max(eW(:))));
assert(min(est.coherence(:)) > 0.95, '4: coherence should be ~1 noise-free');
fprintf('PASS: variance ratio %s; Welch T error median %.2e, max %.2e\n', mat2str(rel.', 3), median(eW(:)), max(eW(:)));
nPass = nPass + 1;

%% 5  noise, seeds, common random numbers
fprintf('\n=== 5: noise and seeds ===\n');
N5 = 32768; sig = 0.02;
[Y0, ~, i0] = synthesiseTwinRecords(twin, N5, dt, 'Seed', 6);
[Yn, ~, in] = synthesiseTwinRecords(twin, N5, dt, 'Seed', 6, 'NoiseStd', sig);
D = Yn - Y0;                                                         % pure noise if CRN holds
sd = std(D, 0, 2).';
assert(all(abs(sd / sig - 1) < 0.03), sprintf('5: noise std %s', mat2str(sd, 3)));
c = corrcoef(D(1, :), Y0(1, :));
assert(abs(c(1, 2)) < 0.03, '5: noise correlated with signal');
cn = corrcoef(D(1, :), D(2, :));
assert(abs(cn(1, 2)) < 0.03, '5: noise correlated across channels');
assert(all(abs(in.noisePSD - sig^2 * dt / pi) < 1e-18) && isequal(in.noiseStd, sig * ones(1, 4)), '5: noisePSD / noiseStd');
eN = extractTransmissibility(D, dt, 'Reference', ref, 'SegmentLength', 512);
assert(abs(mean(eN.Srr) / in.noisePSD(ref) - 1) < 0.05, '5: noise PSD convention differs from extractTransmissibility');
[Yf, ~, iF] = synthesiseTwinRecords(twin, N5, dt, 'Seed', 6, 'NoiseFraction', [0.1 0.2 0.3 0.4]);
assert(maxAll(iF.noiseStd - [0.1 0.2 0.3 0.4] .* i0.signalStd) < 1e-15, '5: NoiseFraction scaling');
sdf = std(Yf - Y0, 0, 2).' ./ i0.signalStd;
assert(all(abs(sdf ./ [0.1 0.2 0.3 0.4] - 1) < 0.03), '5: relative noise level');
% seeds and global RNG
[Ya] = synthesiseTwinRecords(twin, 1024, dt, 'Seed', 9, 'NoiseStd', sig);
[Yb] = synthesiseTwinRecords(twin, 1024, dt, 'Seed', 9, 'NoiseStd', sig);
[Yc] = synthesiseTwinRecords(twin, 1024, dt, 'Seed', 10, 'NoiseStd', sig);
[Yd] = synthesiseTwinRecords(twin, 1024, dt, 'Seed', 9, 'NoiseStd', sig, 'NoiseSeed', 123);
assert(isequal(Ya, Yb) && ~isequal(Ya, Yc) && ~isequal(Ya, Yd), '5: seed contract');
if exist('OCTAVE_VERSION', 'builtin') > 0
    randn('state', 42); g1 = randn(1, 3); randn('state', 42); %#ok<RAND>
    synthesiseTwinRecords(twin, 1024, dt, 'Seed', 9, 'NoiseStd', sig); g2 = randn(1, 3);
else
    rng(42); g1 = randn(1, 3); rng(42);
    synthesiseTwinRecords(twin, 1024, dt, 'Seed', 9, 'NoiseStd', sig); g2 = randn(1, 3);
end
assert(isequal(g1, g2), '5: seeded call disturbed the global RNG');
fprintf('PASS: noise std %s (target %.3f), relative levels %s\n', mat2str(sd, 3), sig, mat2str(sdf, 3));
nPass = nPass + 1;

%% 6  EMM node FRF and disk cache
fprintf('\n=== 6: EMM node FRF and disk cache ===\n');
p0 = [4.6985e-5, 1.4548e-3, 0.3830]; R0 = p0(3);
sensors = [0.3*R0, 0; 0.3*R0, pi; 0.5*R0, pi; 0.9*R0, 0];
cdir = fullfile(tempdir, sprintf('twinCacheTest_%d', floor(1e6 * rem(now, 1))));
if exist(cdir, 'dir'), rmdir(cdir, 's'); end
specE = struct('p', p0, 'sensors', sensors, 'band', [4.5 5.5], 'nodeSpacing', 0.25, 'checkPoints', 2);
t6 = tic; tA = synthesiseTwinRecords(specE, 'CacheDir', cdir); sA = toc(t6);
assert(strcmp(tA.cacheStatus, 'computed') && exist(tA.cacheFile, 'file') == 2, '6: first call should compute and save');
assert(numel(tA.omegaNodes) == 5, '6: expected 5 nodes');
t6 = tic; tB = synthesiseTwinRecords(specE, 'CacheDir', cdir); sB = toc(t6);
assert(strcmp(tB.cacheStatus, 'loaded') && strcmp(tB.cacheFile, tA.cacheFile), '6: second call should load');
assert(isequal(tB.HaNodes, tA.HaNodes) && isequal(tB.omegaNodes, tA.omegaNodes), '6: loaded FRF differs');
assert(sB < 0.5 * sA, sprintf('6: loading (%.2f s) not faster than solving (%.2f s)', sB, sA));
[fS1, oS1] = transmissibilityFeatureVector(log(p0([1 3])), tA.omegaNodes, sensors); %#ok<ASGLU>
Tn = tA.HaNodes(1:3, :) ./ tA.HaNodes(4, :);
e6 = maxAll(Tn - oS1.T) / maxAll(oS1.T);
assert(e6 < 1e-12, sprintf('6: node T differs from S1 (%.2e)', e6));
fprintf('  interpolation check (spacing 0.25): H %.2e, T %.2e\n', tA.interpCheck.maxRelErrH, tA.interpCheck.maxRelErrT);
assert(isfinite(tA.interpCheck.maxRelErrT) && tA.interpCheck.maxRelErrT < 1e-2, '6: interpolation check');
specM = specE; specM.p = [p0(1:2), 0.99 * R0];                       % different parameter point
tM = synthesiseTwinRecords(specM, 'CacheDir', cdir);
assert(strcmp(tM.cacheStatus, 'computed') && ~strcmp(tM.cacheFile, tA.cacheFile), '6: key mismatch reused a file');
specO = specE; specO.frfOptions = {'Truncation', [40 10 10]};        % different truncation
tO = synthesiseTwinRecords(specO, 'CacheDir', cdir);
assert(strcmp(tO.cacheStatus, 'computed'), '6: truncation not part of the key');
specR = specE; specR.checkReference = 2;                            % different diagnostic setting
tR = synthesiseTwinRecords(specR, 'CacheDir', cdir);
assert(strcmp(tR.cacheStatus, 'computed') && tR.interpCheck.reference == 2, '6: stale interpCheck reused');
assert(isequal(tR.HaNodes, tA.HaNodes), '6: FRF should not depend on the check settings');
assert(numel(dir(fullfile(cdir, 'twinFRF_*.mat'))) == 4, '6: expected four cache files');
[Y6, ~, i6] = synthesiseTwinRecords(tB, 4096, dt, 'Seed', 1);
assert(all(isfinite(Y6(:))) && strcmp(i6.cacheStatus, 'memory'), '6: synthesis from loaded twin');
rmdir(cdir, 's');
fprintf('PASS: solve %.1f s, load %.2f s, node T vs S1 %.1e\n', sA, sB, e6); nPass = nPass + 1;

%% 7  errors
fprintf('\n=== 7: error handling ===\n');
cases = {
    'badSpec',          {struct('sensors', sensorsToy)}
    'badSpec',          {struct('frf', Hfun, 'sensors', sensorsToy(1, :))}
    'badSpec',          {struct('frf', Hfun, 'sensors', sensorsToy, 'band', [5 4])}
    'badSpec',          {struct('frf', 3, 'sensors', sensorsToy)}
    'badFRF',           {struct('frf', @(w) ones(3, numel(w)), 'sensors', sensorsToy)}
    'cacheNotAllowed',  {spec, 'CacheDir', tempdir}
    'badOption',        {spec, 'Seed', 1}
    'badOption',        {twin, 1024, dt, 'Nope', 1}
    'badOption',        {twin, 1024}
    'badN',             {twin, 10.5, dt}
    'badDt',            {twin, 1024, -1}
    'aboveNyquist',     {twin, 1024, 0.5}
    'tooFewBins',       {synthesiseTwinRecords(struct('frf', Hfun, 'sensors', sensorsToy, 'band', [4 4.5])), 16, dt}
    'badOption',        {twin, 1024, dt, 'TaperWidth', 4}
    'badOption',        {twin, 1024, dt, 'Seed', -1}
    'badNoise',         {twin, 1024, dt, 'NoiseStd', [1 2]}
    'badNoise',         {twin, 1024, dt, 'NoiseStd', 0.1, 'NoiseFraction', 0.1}
    'badSpectrum',      {twin, 1024, dt, 'Spectrum', 'pm'}
    'badSpectrum',      {twin, 1024, dt, 'Hs', -1}
    'badSpectrum',      {twin, 1024, dt, 'Spectrum', @(w) -ones(size(w))}
    };
for c = 1:size(cases, 1)
    id = ['synthesiseTwinRecords:', cases{c, 1}];
    try
        if size(cases{c, 2}, 2) >= 3 && isnumeric(cases{c, 2}{2})
            [~, ~, ~, ~] = synthesiseTwinRecords(cases{c, 2}{:});
        else
            synthesiseTwinRecords(cases{c, 2}{:});
        end
        error('tester:noError', '7.%d: expected %s but no error was thrown', c, id);
    catch err
        assert(strcmp(err.identifier, id), sprintf('7.%d: expected %s, got %s (%s)', c, id, err.identifier, err.message));
    end
end
fprintf('PASS: %d cases\n', size(cases, 1)); nPass = nPass + 1;

fprintf('\nAll %d groups passed in %.1f s.\n', nPass, toc(tStart));

%% ================================================================================================
function H = modalAcc(sys, w)
% Acceleration FRF of a modal system under exp(+i w t): sum_m Phi(:,m) (-w^2)/(wn^2 - w^2 + 2i zeta wn w).
w = w(:).';
H = zeros(size(sys.Phi, 1), numel(w));
for m = 1:numel(sys.wn)
    H = H + sys.Phi(:, m) * (-w.^2 ./ (sys.wn(m)^2 - w.^2 + 2i * sys.zeta(m) * sys.wn(m) * w));
end
end