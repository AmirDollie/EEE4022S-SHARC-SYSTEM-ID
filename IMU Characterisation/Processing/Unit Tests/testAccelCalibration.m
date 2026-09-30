%% testAccelCalibration.m
% Known-answer tests for calibrateAccelSixPosition.m on synthetic static holds.
%
% Lives in IMU Characterisation/Processing/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
nPass = 0; nFail = 0;
g = 9.796;
rng(7);
sTrue = [1.0021 0.9964 1.0043];            % calibration scale (sensor reads 1/s times true)
bTrue = [0.30 -0.12 0.78];                 % offsets, m/s^2 (up to 80 mg, as observed)
noise = 2e-4;                              % m/s^2 on each hold mean (about 20 ug)
raw = @(G, M, b) (M \ G.').' + b;          % true gravity vectors -> raw readings

% the seven holds actually logged: up, down, sides at 0/90/180 deg, two in-between rolls
roll = [0 90 180 -34 135] * pi / 180;
G7 = [0 0 g; 0 0 -g; g * [sin(roll(:)) -cos(roll(:)) 0.05 * ones(5, 1)]];
G7(3:end, :) = g * G7(3:end, :) ./ sqrt(sum(G7(3:end, :).^2, 2));

% 1  diag model on the real 7-hold geometry
A = raw(G7, diag(sTrue), bTrue) + noise * randn(7, 3);
c = calibrateAccelSixPosition(A, g, 'TwoPointAxis', [1 2 3]);
eS = max(abs(c.scale ./ sTrue - 1)); eB = max(abs(c.offset - bTrue));
[nPass, nFail] = check(c.converged && eS < 2e-4 && eB < 2e-3, sprintf(['diag fit, 7 real-geometry holds: scale err %.1e, ' ...
    'offset err %.1e m/s^2, dof %d'], eS, eB, c.dof), nPass, nFail);

% 2  the sphere fit needs no orientation knowledge: random 2 deg tilts of every hold change nothing
Gt = G7;
for h = 1:7
    ax = randn(3, 1); ax = ax / norm(ax); th = 2 * pi / 180;
    K = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];
    Gt(h, :) = (expm(th * K) * G7(h, :).').';
end
At = raw(Gt, diag(sTrue), bTrue) + noise * randn(7, 3);
ct = calibrateAccelSixPosition(At, g, 'TwoPointAxis', [1 2 3]);
eFit = abs(ct.scale(3) / sTrue(3) - 1); eTwo = abs(ct.twoPoint.gain * sTrue(3) - 1);
[nPass, nFail] = check(eFit < 2e-4, sprintf('tilted holds: z scale from the fit still within %.1e (two-point: %.1e)', eFit, eTwo), ...
    nPass, nFail);

% 3  two-point estimate on aligned holds
e2 = abs(c.twoPoint.gain * sTrue(3) - 1);
[nPass, nFail] = check(e2 < 1e-4 && abs(c.twoPoint.offset - bTrue(3)) < 1e-3, sprintf('two-point z on aligned holds: gain err %.1e', e2), ...
    nPass, nFail);

% 4  leave-one-out needs >= 8 holds: unavailable with the 7 logged, complete once 3 tilted holds are added
[nPass, nFail] = check(~c.loo.available && c.loo.nDetermined == 0, ...
    '7 holds: leave-one-out correctly unavailable (each refit would have no residual)', nPass, nFail);
tilt = [45 45 45] * pi / 180; az = [0 120 240] * pi / 180;
G10 = [G7; g * [sin(tilt(:)) .* cos(az(:)), sin(tilt(:)) .* sin(az(:)), cos(tilt(:))]];
A10 = raw(G10, diag(sTrue), bTrue) + noise * randn(10, 3);
c10 = calibrateAccelSixPosition(A10, g);
[nPass, nFail] = check(c10.loo.available && c10.loo.nDetermined == 10 && max(c10.looSpreadScale ./ c10.scale) < 5e-4, ...
    sprintf('10 holds (+3 tilted): all 10 refits determined, max half-range %.1e, dof %d', ...
    max(c10.looSpreadScale ./ c10.scale), c10.dof), nPass, nFail);

% 4b hold-count guards at the boundary: diag refuses 6, full refuses 9
threw6 = false; threw9 = false;
try, calibrateAccelSixPosition(A(1:6, :), g); catch, threw6 = true; end
V9 = randn(9, 3); A9 = raw(g * V9 ./ sqrt(sum(V9.^2, 2)), diag(sTrue), bTrue);
try, calibrateAccelSixPosition(A9, g, 'Model', 'full'); catch, threw9 = true; end
[nPass, nFail] = check(threw6 && threw9, 'guards: diag refuses 6 holds, full refuses 9 (need more holds than parameters)', ...
    nPass, nFail);

% 5  full model: 7 holds refused, 14 well-spread holds recover non-orthogonality
Mtrue = [sTrue(1) 0.004 -0.003; 0 sTrue(2) 0.002; 0 0 sTrue(3)];
threw = false;
try, calibrateAccelSixPosition(A, g, 'Model', 'full'); catch, threw = true; end
V = randn(14, 3); G14 = g * V ./ sqrt(sum(V.^2, 2));
A14 = raw(G14, Mtrue, bTrue) + noise * randn(14, 3);
cf = calibrateAccelSixPosition(A14, g, 'Model', 'full');
eM = max(abs(cf.M(:) - Mtrue(:)));
[nPass, nFail] = check(threw && cf.converged && eM < 1e-3, sprintf('full model: refused with 7 holds; 14 holds recover M within %.1e', eM), ...
    nPass, nFail);

% 6  a diag fit on data WITH non-orthogonality leaves a visible residual (the model limit is detectable)
A12 = raw(G14(1:12, :), Mtrue, bTrue);
cd12 = calibrateAccelSixPosition(A12, g);
[nPass, nFail] = check(cd12.rmsResidual > 1e-3, sprintf('diag fit on skewed axes: residual %.1e m/s^2 flags the missing terms', ...
    cd12.rmsResidual), nPass, nFail);

% 7  stillness: a slow 0.3 deg roll at constant |a| is WARNED (drift) but kept, because it does not
%    corrupt the mean; a 15 deg swing within the hold is EXCLUDED by the averaging-loss test
n = 1600; t = (0:n - 1).' / 10;
aStill = [zeros(n, 1), zeros(n, 1), g * ones(n, 1)] + 0.02 * randn(n, 3);
th = 0.3 * pi / 180 * t / t(end);
aRoll = g * [zeros(n, 1), sin(th), cos(th)] + 0.02 * randn(n, 3);
th2 = 15 * pi / 180 * t / t(end);
aSwing = g * [zeros(n, 1), sin(th2), cos(th2)] + 0.02 * randn(n, 3);
[s1, m1] = isStaticHold(aStill, 5e-4 * randn(n, 3));
[s2, m2] = isStaticHold(aRoll, []);
[s3, m3] = isStaticHold(aSwing, []);
ok = s1 && isempty(m1.warnings) && s2 && any(strcmp(m2.warnings, 'drift')) && m2.magStd <= 0.05 && ...
    ~s3 && any(strcmp(m3.fails, 'averaging'));
[nPass, nFail] = check(ok, sprintf(['stillness: still hold clean; 0.3 deg roll warned (drift %.3f, |a| std %.3f) but kept; ' ...
    '15 deg swing excluded (averaging loss %.1e)'], m2.drift, m2.magStd, m3.avgLoss), nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testAccelCalibration:failed', '%d test(s) failed.', nFail); end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end