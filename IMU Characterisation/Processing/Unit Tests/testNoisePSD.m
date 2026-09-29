%% testNoisePSD.m
% Known-answer tests for computeNoisePSD.m (and alignImuSensors.m on a synthetic record).
%
% Lives in IMU Characterisation/Processing/Unit-Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
nPass = 0; nFail = 0;
fs = 10;

% ---- 1  white noise of known one-sided density ------------------------------------------------------
Nw = 2e-3;
rec = synthesiseImuNoise(4 * 3600, fs, 'Accel', [Nw 0 0], 'Gyro', [1e-4 0 0], 'Seed', 21);
x = rec.acc(:, 1);
o = computeNoisePSD(x, fs, 'SegmentLength', 1024);
in = o.f > 0.05 * fs & o.f < 0.45 * fs;
lev = mean(o.Pxx(in)) / Nw^2;
[nPass, nFail] = check(abs(lev - 1) < 0.02, sprintf('white: mean PSD / N^2 = %.4f', lev), nPass, nFail);
sig2 = var(x);
[nPass, nFail] = check(abs(mean(o.PxxPerRadS(in)) / (sig2 / fs / pi) - 1) < 0.02, ...
    'white: per rad/s level = sigma^2 dt / pi (the twin''s S_nn convention)', nPass, nFail);
% scatter of individual bins consistent with the equivalent dof
r = o.Pxx(in) / Nw^2;
[nPass, nFail] = check(abs(std(r) / sqrt(2 / o.dofEffective) - 1) < 0.15, ...
    sprintf('white: bin scatter %.4f vs predicted sqrt(2/dof) %.4f', std(r), sqrt(2 / o.dofEffective)), nPass, nFail);
inCI = mean(Nw^2 >= o.ci95(1) * o.Pxx(in) & Nw^2 <= o.ci95(2) * o.Pxx(in));
[nPass, nFail] = check(inCI > 0.92 && inCI < 0.98, sprintf('white: 95%% interval covers truth in %.1f%% of bins', 100 * inCI), nPass, nFail);

% ---- 2  Parseval ---------------------------------------------------------------------------------------
[nPass, nFail] = check(abs(sum(o.Pxx) * o.df / sig2 - 1) < 0.02, sprintf('Parseval: sum(P) df / var = %.4f', ...
    sum(o.Pxx) * o.df / sig2), nPass, nFail);
o2 = computeNoisePSD(x, fs, 'SegmentLength', 1024, 'Window', 'rect', 'Overlap', 0);
[nPass, nFail] = check(abs(sum(o2.Pxx) * o2.df / sig2 - 1) < 0.02, 'Parseval with a rectangular window', nPass, nFail);

% ---- 3  frequency axis, DC / Nyquist, a bin-centred sinusoid ----------------------------------------------------
L = 256; t = (0:L * 64 - 1).' / fs; f0 = 37 * fs / L; Amp = 0.05;
s = Amp * sin(2 * pi * f0 * t);
o = computeNoisePSD(s, fs, 'SegmentLength', L);
[~, ipk] = max(o.Pxx);
[nPass, nFail] = check(numel(o.f) == L / 2 + 1 && o.f(1) == 0 && abs(o.f(end) - fs / 2) < 1e-12 && ...
    abs(o.df - fs / L) < 1e-12, 'axis: 0 .. fs/2 in L/2+1 bins of fs/L', nPass, nFail);
[nPass, nFail] = check(abs(o.f(ipk) - f0) < 1e-12, sprintf('sinusoid peak at %.4f Hz', o.f(ipk)), nPass, nFail);
[nPass, nFail] = check(abs(sum(o.Pxx) * o.df / (Amp^2 / 2) - 1) < 0.01, 'sinusoid power A^2/2', nPass, nFail);
o = computeNoisePSD(Amp * cos(pi * fs * t), fs, 'SegmentLength', L, 'Window', 'rect', 'Detrend', 'none');
[nPass, nFail] = check(abs(o.Pxx(end) * o.df / Amp^2 - 1) < 1e-9, 'Nyquist bin not doubled', nPass, nFail);

% ---- 4  coherence ------------------------------------------------------------------------------------------------
rng(22); M = 2^15;
common = randn(M, 1);
Y = [common + 0.1 * randn(M, 1), common + 0.1 * randn(M, 1), randn(M, 1)];
o = computeNoisePSD(Y, fs, 'SegmentLength', 512);
k = o.f > 0;
[nPass, nFail] = check(mean(o.coherence(k, 1)) > 0.97, sprintf('coherence: shared signal %.3f', mean(o.coherence(k, 1))), nPass, nFail);
indep = mean(o.coherence(k, 2:3), 1);
[nPass, nFail] = check(all(indep < 3 * o.coherenceNoiseFloor), sprintf('coherence: independent %.4f %.4f (floor %.4f)', ...
    indep, o.coherenceNoiseFloor), nPass, nFail);
[nPass, nFail] = check(isequal(o.pairs, [1 2; 1 3; 2 3]), 'default pairs for 3 channels', nPass, nFail);
% significance threshold: independent channels exceed it in about 5% of bins (short record, few averages)
rng(26); Z = randn(720, 2);
exceed = [];
for rep = 1:200
    Z = randn(720, 2);
    oz = computeNoisePSD(Z, fs, 'SegmentLength', 128);
    kk = oz.f > 0 & oz.f < fs / 2;
    exceed(end + 1) = mean(oz.coherence(kk) > oz.coherence95Threshold); %#ok<SAGROW>
end
[nPass, nFail] = check(mean(exceed) > 0.03 && mean(exceed) < 0.08, sprintf(['coherence 95%% threshold %.3f (K_eff %.1f): ' ...
    'independent noise exceeds it in %.1f%% of bins'], oz.coherence95Threshold, oz.dofEffective / 2, 100 * mean(exceed)), nPass, nFail);

% ---- 5  coloured noise slopes (also validates the synthesiser) ---------------------------------------------
rec = synthesiseImuNoise(4 * 3600, fs, 'Accel', [0 1e-3 0], 'Gyro', [1e-4 0 0], 'Seed', 23);
o = computeNoisePSD(rec.acc(:, 1), fs, 'SegmentLength', 8192);
in = o.f > 0.01 & o.f < 1;
pf = polyfit(log10(o.f(in)), log10(o.Pxx(in)), 1);
[nPass, nFail] = check(abs(pf(1) + 1) < 0.1, sprintf('flicker slope %.3f (expect -1)', pf(1)), nPass, nFail);
lvl = mean(o.Pxx(in) .* o.f(in)) / (1e-3^2 / (2 * pi));
[nPass, nFail] = check(abs(lvl - 1) < 0.1, sprintf('flicker level S f / (B^2/2pi) = %.3f', lvl), nPass, nFail);
rec = synthesiseImuNoise(4 * 3600, fs, 'Accel', [0 0 1e-4], 'Gyro', [1e-4 0 0], 'Seed', 24);
o = computeNoisePSD(rec.acc(:, 1), fs, 'SegmentLength', 8192, 'Detrend', 'linear');
in = o.f > 0.01 & o.f < 0.5;
pf = polyfit(log10(o.f(in)), log10(o.Pxx(in)), 1);
[nPass, nFail] = check(abs(pf(1) + 2) < 0.15, sprintf('random-walk slope %.3f (expect -2)', pf(1)), nPass, nFail);

% ---- 6  band summaries ---------------------------------------------------------------------------------------------
rec = synthesiseImuNoise(2 * 3600, fs, 'Accel', [Nw 0 0], 'Gyro', [1e-4 0 0], 'Seed', 25);
o = computeNoisePSD(rec.acc(:, 1), fs, 'SegmentLength', 4096);
bw = sum(o.f >= 0.47 & o.f <= 1.35) * o.df;
[nPass, nFail] = check(abs(o.bands.basin.rms / (Nw * sqrt(bw)) - 1) < 0.05, sprintf('basin band rms %.3g (expect %.3g)', ...
    o.bands.basin.rms, Nw * sqrt(bw)), nPass, nFail);
o = computeNoisePSD(rec.acc(1:600, 1), fs, 'SegmentLength', 64);          % df = 0.156 Hz: at most 1 MIZ bin
[nPass, nFail] = check(o.bands.miz.nBins < 3 && ~isempty(o.bands.miz.note), ['coarse resolution flagged: ', o.bands.miz.note], ...
    nPass, nFail);

% ---- 7  timing guard ---------------------------------------------------------------------------------------------
tt = (0:5999).' / fs + 1e-3 * randn(6000, 1);
ok = true;
try, computeNoisePSD(randn(6000, 1), [], 'Time', tt); catch, ok = false; end
[nPass, nFail] = check(ok, 'mild jitter (1 ms at 10 Hz) accepted, fs taken from Time', nPass, nFail);
tt2 = tt; tt2(3000:end) = tt2(3000:end) + 0.5;
refused = false;
try, computeNoisePSD(randn(6000, 1), fs, 'Time', tt2); catch err, refused = strcmp(err.identifier, 'computeNoisePSD:irregular'); end
[nPass, nFail] = check(refused, 'a 0.5 s gap is refused', nPass, nFail);

% ---- 8  alignImuSensors on a synthetic 3-sensor record with one lost line ------------------------------------------
t0 = (0:99).' * 0.1;
mk = @(nm, dtk, drop) struct('name', nm, 't', setdiff1(t0 + dtk, drop), 'acc', zeros(100 - numel(drop), 3), ...
    'gyr', zeros(100 - numel(drop), 3), 'nRows', 100 - numel(drop));
recS.sensors = [mk('A', 0, []), mk('B', 0.001, 50), mk('C', 0.002, [])];
A = alignImuSensors(recS);
[nPass, nFail] = check(numel(A.t) == 99 && A.nDropped == 1, 'align: a tick missing in one sensor is dropped from all', nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testNoisePSD:failed', '%d test(s) failed.', nFail); end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end

function t = setdiff1(t, drop)
t(drop) = [];
end