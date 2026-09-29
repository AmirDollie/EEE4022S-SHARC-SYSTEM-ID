%% testAllanDeviation.m
% Known-answer tests for allanDeviation.m, using synthesiseImuNoise.m records.
% Settles the flicker constant under this project's convention (sqrt(ln2/pi) B = 0.470 B).
%
% Lives in IMU Characterisation/Processing/Unit-Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
nPass = 0; nFail = 0;
fs = 10; T = 6 * 3600;

% ---- 1  estimator equals the brute-force definition -------------------------------------------------------
rng(31); x = randn(500, 1); dt = 1 / fs;
o = allanDeviation(x, fs);
err = 0;
for i = 1:numel(o.m)
    k = o.m(i); acc = 0; nn = 0;
    for j = 1:(500 - 2 * k + 1)
        y1 = mean(x(j:j + k - 1)); y2 = mean(x(j + k:j + 2 * k - 1));
        acc = acc + (y2 - y1)^2; nn = nn + 1;
    end
    err = max(err, abs(sqrt(acc / (2 * nn)) / o.adev(i) - 1));
end
[nPass, nFail] = check(err < 1e-10, sprintf('matches brute-force overlapping definition (max rel err %.1e)', err), nPass, nFail);

% ---- 2  white noise: slope -1/2, level N / sqrt(2 tau) -----------------------------------------------------
Nw = 2e-3;
rec = synthesiseImuNoise(T, fs, 'Accel', [Nw 0 0], 'Gyro', [1e-4 0 0], 'Seed', 32);
o = allanDeviation(rec.acc(:, 1), fs);
k = o.tau >= 0.5 & o.tau <= 100;
pf = polyfit(log10(o.tau(k)), log10(o.adev(k)), 1);
[nPass, nFail] = check(abs(pf(1) + 0.5) < 0.03, sprintf('white: slope %.3f (expect -0.5)', pf(1)), nPass, nFail);
[nPass, nFail] = check(abs(o.terms.white / Nw - 1) < 0.03, sprintf('white: N = %.4g (truth %.4g, one-sided)', o.terms.white, Nw), ...
    nPass, nFail);
[nPass, nFail] = check(abs(o.terms.whiteIEEE / (Nw / sqrt(2)) - 1) < 0.03, ...
    sprintf('white: IEEE coefficient = one-sided N / sqrt(2) = sigma_A(1 s): %.4g (expect %.4g)', o.terms.whiteIEEE, Nw / sqrt(2)), nPass, nFail);
[nPass, nFail] = check(abs(o.terms.whiteIEEE / interp1(log(o.tau), o.adev, 0) - 1) < 0.03, ...
    'white: whiteIEEE equals the Allan deviation read at tau = 1 s', nPass, nFail);
[nPass, nFail] = check(isnan(o.terms.randomWalk), 'white: no random-walk region reported', nPass, nFail);
% confidence intervals cover the theory curve at most taus
theory = Nw ./ sqrt(2 * o.tau);
cover = mean(theory >= o.adevLo & theory <= o.adevHi);
[nPass, nFail] = check(cover > 0.85, sprintf('white: approx. 95%% intervals cover theory at %.0f%% of taus', 100 * cover), nPass, nFail);

% ---- 3  random walk: slope +1/2, level K sqrt(tau/3) ------------------------------------------------------------
K = 1e-4;
rec = synthesiseImuNoise(T, fs, 'Accel', [0 0 K], 'Gyro', [1e-4 0 0], 'Seed', 33);
o = allanDeviation(rec.acc(:, 1), fs);
k = o.tau >= 1 & o.tau <= 600;
pf = polyfit(log10(o.tau(k)), log10(o.adev(k)), 1);
[nPass, nFail] = check(abs(pf(1) - 0.5) < 0.05, sprintf('random walk: slope %.3f (expect +0.5)', pf(1)), nPass, nFail);
[nPass, nFail] = check(abs(o.terms.randomWalk / K - 1) < 0.1, sprintf('random walk: K = %.4g (truth %.4g)', o.terms.randomWalk, K), ...
    nPass, nFail);

% ---- 4  flicker: flat plateau at sqrt(ln2/pi) B = 0.470 B (this project's convention) ------------------------------
B = 1e-3; lev = zeros(1, 4);
for s = 1:4
    rec = synthesiseImuNoise(T, fs, 'Accel', [0 B 0], 'Gyro', [1e-4 0 0], 'Seed', 33 + s);
    o = allanDeviation(rec.acc(:, 1), fs);
    k = o.tau >= 1 & o.tau <= 300;
    lev(s) = mean(o.adev(k)) / B;
end
[nPass, nFail] = check(abs(mean(lev) / sqrt(log(2) / pi) - 1) < 0.07, sprintf(['flicker: plateau / B = %.3f (4 seeds; ' ...
    'predicted sqrt(ln2/pi) = %.3f, NOT the IEEE 0.664)'], mean(lev), sqrt(log(2) / pi)), nPass, nFail);
pf = polyfit(log10(o.tau(k)), log10(o.adev(k)), 1);
[nPass, nFail] = check(abs(pf(1)) < 0.08, sprintf('flicker: slope %.3f (expect 0)', pf(1)), nPass, nFail);
[nPass, nFail] = check(abs(o.terms.flicker / B - 1) < 0.15, sprintf('flicker: B from the plateau = %.4g (truth %.4g)', ...
    o.terms.flicker, B), nPass, nFail);

% ---- 5  all three terms together, over 5 records (white dominates < ~40 s, flicker plateau
%          ~100-400 s, random walk beyond ~1000 s). A single record is not a fair test of the
%          random-walk term: it only dominates at the last trusted taus, where few independent
%          averages exist, so whether a clean +1/2 run appears varies from draw to draw. The
%          test asks for it in at least 3 of 5 records (about 90% of draws find it).
Nw = 2e-3; B = 5e-4; K = 1e-5; R = NaN(5, 3);
for sd = 1:5
    rec = synthesiseImuNoise(24 * 3600, fs, 'Accel', [Nw B K], 'Gyro', [1e-4 0 0], 'Seed', 40 + sd);
    t = allanDeviation(rec.acc(:, 1), fs).terms;
    R(sd, :) = [t.white / Nw, t.flicker / B, t.randomWalk / K];
end
fprintf('    combined (ratio to truth, 5 records): N %s  B %s  K %s\n', mat2str(R(:, 1).', 3), mat2str(R(:, 2).', 3), ...
    mat2str(R(:, 3).', 3));
[nPass, nFail] = check(all(abs(R(:, 1) - 1) < 0.1), 'combined: white density within 10% in every record', nPass, nFail);
[nPass, nFail] = check(abs(median(R(:, 2)) - 1) < 0.3, ...
    'combined: flicker B median within 30% (white and random walk lift the plateau)', nPass, nFail);
kf = R(isfinite(R(:, 3)), 3);
[nPass, nFail] = check(numel(kf) >= 3 && abs(median(kf) - 1) < 0.5, sprintf(['combined: random walk found in %d of 5 ' ...
    'records, median within 50%%'], numel(kf)), nPass, nFail);

% ---- 6  grid, multi-channel, units scaling and timing guard ----------------------------------------------------------
rng(41); Y = randn(20000, 3) .* [1 2 3];
o = allanDeviation(Y, fs);
[nPass, nFail] = check(o.tau(1) == 1 / fs && o.m(end) == floor((20000 - 1) / 2) && all(diff(o.m) > 0), ...
    'grid: tau from 1/fs to (N-1)/2 samples, strictly increasing', nPass, nFail);
kk = o.edf >= 20;
[nPass, nFail] = check(max(abs(o.adev(kk, 2) ./ o.adev(kk, 1) - 2)) < 0.6 && size(o.adev, 2) == 3, ...
    'multi-channel: columns handled independently (scale 2 channel ~ 2x)', nPass, nFail);
o2 = allanDeviation(Y(:, 1) * 9.81, fs);
[nPass, nFail] = check(max(abs(o2.adev ./ o.adev(:, 1) - 9.81)) < 1e-9, 'linear in the input (unit scaling)', nPass, nFail);
tt = (0:5999).' / fs; tt(3000:end) = tt(3000:end) + 0.5;
refused = false;
try, allanDeviation(randn(6000, 1), fs, 'Time', tt); catch err, refused = strcmp(err.identifier, 'allanDeviation:irregular'); end
[nPass, nFail] = check(refused, 'a 0.5 s gap is refused', nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testAllanDeviation:failed', '%d test(s) failed.', nFail); end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end