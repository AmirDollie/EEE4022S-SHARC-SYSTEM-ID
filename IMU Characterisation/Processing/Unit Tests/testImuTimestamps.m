%% testImuTimestamps.m
% Unit tests for synthesiseImuNoise.m (timing part) and checkImuTimestamps.m: every
% defect is injected with a known size and must be recovered.
%
% Lives in IMU Characterisation/Processing/Unit-Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
nPass = 0; nFail = 0;
fs = 100;
wrap32us = 2^32 * 1e-6;

% ---- A  counter, drops, duplicates, counter wraps, timestamp wrap, host clock -------------------------
drops = [1000 3; 50000 1; 120000 40];
dups = [5000 9000 150000];
rec = synthesiseImuNoise(1800, fs, 'OdrPpm', 25, 'TsPpm', -12, 'TsJitter', 20e-6, 'TsResolution', 1e-6, ...
    'TsOffset', 4000, 'TsWrap', wrap32us, 'CounterBits', 16, 'HostPpm', 40, 'Drop', drops, 'Duplicate', dups, ...
    'Seed', 11);
tr = rec.truth;
o = checkImuTimestamps(rec.tau, 'NominalRate', fs, 'Counter', rec.counter, 'CounterBits', 16, ...
    'TsWrap', wrap32us, 'Host', rec.host, 'Verbose', true);
[nPass, nFail] = check(o.nDuplicates == 3 && isequal(sort(o.index(o.duplicateRows)), sort(tr.duplicatedIndex(:))), ...
    'A: duplicates found at the right samples', nPass, nFail);
[nPass, nFail] = check(o.nGaps == 3 && o.nMissing == 44 && isequal(o.gaps(:, 2:3), [drops(:, 1) - 1, drops(:, 2)]), ...
    'A: gaps found with the right start and length', nPass, nFail);
[nPass, nFail] = check(o.nWraps == 1 && o.nBackwardJumps == 0, 'A: one timestamp wrap unwrapped', nPass, nFail);
[nPass, nFail] = check(abs(o.rateErrorPpm - tr.rateErrorPpm) < 0.01, ...
    sprintf('A: rate error %.4f ppm (truth %.4f)', o.rateErrorPpm, tr.rateErrorPpm), nPass, nFail);
jTrue = sqrt((20e-6)^2 + (1e-6)^2 / 12);
[nPass, nFail] = check(abs(o.jitterStd / jTrue - 1) < 0.03, sprintf('A: jitter std %.2f us (truth %.2f)', ...
    1e6 * o.jitterStd, 1e6 * jTrue), nPass, nFail);
[nPass, nFail] = check(abs(o.host.relativePpm - tr.hostVsDevicePpm) < 0.5, sprintf('A: host vs device %.2f ppm (truth %.2f)', ...
    o.host.relativePpm, tr.hostVsDevicePpm), nPass, nFail);
[nPass, nFail] = check(abs(o.driftPpm) < 0.05, sprintf('A: no drift injected, estimate %.4f ppm', o.driftPpm), nPass, nFail);

% ---- B  same record without the counter -----------------------------------------------------------
o = checkImuTimestamps(rec.tau, 'NominalRate', fs, 'TsWrap', wrap32us, 'Verbose', false);
[nPass, nFail] = check(o.nDuplicates == 3 && o.nGaps == 3 && o.nMissing == 44 && abs(o.rateErrorPpm - tr.rateErrorPpm) < 0.01, ...
    'B: no counter: same duplicates, gaps and rate from timestamps alone', nPass, nFail);

% ---- C  a gap longer than the counter period (counter alone would alias) ------------------------------
rec = synthesiseImuNoise(1200, fs, 'TsJitter', 10e-6, 'CounterBits', 16, 'Drop', [20000 70000], 'Seed', 12);
o = checkImuTimestamps(rec.tau, 'NominalRate', fs, 'Counter', rec.counter, 'CounterBits', 16, 'Verbose', false);
[nPass, nFail] = check(o.nGaps == 1 && o.nMissing == 70000, sprintf('C: 70000-sample gap past a 16-bit counter wrap (found %d)', ...
    o.nMissing), nPass, nFail);

% ---- D  clock drift: device clock rate changes by 3 ppm over the record ---------------------------------
rec = synthesiseImuNoise(3600, fs, 'TsPpm', 5, 'TsDriftPpm', 3, 'TsJitter', 20e-6, 'Seed', 13);
o = checkImuTimestamps(rec.tau, 'NominalRate', fs, 'Counter', rec.counter, 'Verbose', false);
[nPass, nFail] = check(abs(o.driftPpm - (-3)) < 0.05, sprintf('D: drift %.3f ppm (truth -3: faster device clock, lower apparent rate)', ...
    o.driftPpm), nPass, nFail);
[nPass, nFail] = check(o.windowRatePpm(1) > o.windowRatePpm(end), 'D: windowed rate falls across the record', nPass, nFail);

% ---- D2 drift with a very asymmetric missing-data pattern (first 40% mostly lost) ---------------------------
drop = [(1000:200:140000).' , 150 * ones(numel(1000:200:140000), 1)];
rec = synthesiseImuNoise(3600, fs, 'TsPpm', 5, 'TsDriftPpm', 3, 'TsJitter', 20e-6, 'Drop', drop, 'Seed', 18);
o = checkImuTimestamps(rec.tau, 'NominalRate', fs, 'Counter', rec.counter, 'Verbose', false);
[nPass, nFail] = check(abs(o.driftPpm - (-3)) < 0.05 && o.nMissing == sum(drop(:, 2)), ...
    sprintf('D2: asymmetric gaps (%d missing): drift %.3f ppm', o.nMissing, o.driftPpm), nPass, nFail);

% ---- E  resolution and a coarse millisecond timestamp -----------------------------------------------------
rec = synthesiseImuNoise(600, fs, 'TsJitter', 0.3e-3, 'TsResolution', 1e-3, 'Seed', 14);
o = checkImuTimestamps(rec.tau, 'NominalRate', fs, 'Counter', rec.counter, 'Verbose', false);
[nPass, nFail] = check(abs(o.apparentResolution - 1e-3) < 1e-9, sprintf('E: apparent resolution %.3g ms', 1e3 * o.apparentResolution), nPass, nFail);
% no jitter, exact nominal rate, 1 us timestamps: every increment is exactly 10 ms, so the
% quantum cannot be seen and must be reported as undetermined (NaN), not guessed
rec = synthesiseImuNoise(60, fs, 'TsResolution', 1e-6, 'Seed', 16);
o = checkImuTimestamps(rec.tau, 'NominalRate', fs, 'Counter', rec.counter, 'Verbose', false);
[nPass, nFail] = check(isnan(o.apparentResolution), 'E: ideal timestamps: resolution undetermined (NaN)', nPass, nFail);
% no jitter but a 37 ppm rate offset: increments step between adjacent 1 us levels
rec = synthesiseImuNoise(600, fs, 'TsPpm', 37, 'TsResolution', 1e-6, 'Seed', 17);
o = checkImuTimestamps(rec.tau, 'NominalRate', fs, 'Counter', rec.counter, 'Verbose', false);
[nPass, nFail] = check(abs(o.apparentResolution - 1e-6) < 1e-12, sprintf('E: no jitter, rate offset: resolution %.3g us', ...
    1e6 * o.apparentResolution), nPass, nFail);

% ---- F  white-noise level of the synthesiser (one-sided density convention) --------------------------------
rec = synthesiseImuNoise(600, fs, 'Accel', [2e-3 0 0], 'Gyro', [1e-4 0 0], 'Seed', 15);
sdA = std(rec.acc(:, 1)); sdG = std(rec.gyr(:, 2));
[nPass, nFail] = check(abs(sdA / (2e-3 * sqrt(fs / 2)) - 1) < 0.02 && abs(sdG / (1e-4 * sqrt(fs / 2)) - 1) < 0.02, ...
    'F: white std = N sqrt(fs/2) (one-sided density)', nPass, nFail);
[nPass, nFail] = check(abs(mean(rec.acc(:, 3)) - 9.81) < 1e-3, 'F: gravity on z', nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testImuTimestamps:failed', '%d test(s) failed.', nFail); end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end