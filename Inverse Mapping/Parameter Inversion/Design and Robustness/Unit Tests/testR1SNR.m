%% testR1SNR.m
% Pre-flight checks for study R1. Must pass (all) before runR1SNR is run. Uses the cached p0 field (no EMM
% solve); the boundary refinement is r1Boundary, the function the runner uses.
%
%   1  41-point log grid from 1e-8 to 1e-4, uniform in log10
%   2  every marked noise case exists in cfg.noise and lies inside the swept interval; R1 sea = E1 sea
%   3  the R1 path at LSM6DSV16X reproduces the frozen evaluateScenario default bitwise
%   4  only S_nn changes: f0, J and the Welch feature bias b_W are identical at both ends of the axis; F is not
%   5  SNR equivalence: (S_nn x 4, Hs x 2) gives the same sigma, d_W and F as (S_nn, Hs)
%   6  no hidden change: every resolved scenario field except Snn equals the E1 default
%   7  r1Boundary on a known function returns the root to tolerance, inside its bracket
%   8  r1Boundary on the real sigma_lnbeta curve: the refined point lies in its bracket and the metric there
%      matches the threshold to within the bracket's own variation
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
P = setupStudyPaths();
cfg = defineStudyScenarios();
nPass = 0; nFail = 0;
Snn = cfg.R1.Snn;
ev = @(s, sea) evaluateScenario(struct('Snn', s, 'sea', sea, 'cacheDir', P.cache), cfg);

%% ---- 1-2  configuration ---------------------------------------------------------------------------------------
dx = diff(log10(Snn));
ok = numel(Snn) == 41 && abs(Snn(1) - 1e-8) < 1e-20 && abs(Snn(end) - 1e-4) < 1e-16 && max(abs(dx - 0.1)) < 1e-12;
[nPass, nFail] = check(ok, sprintf('1  %d S_nn from %.3g to %.3g, step %.3f decade', numel(Snn), Snn(1), Snn(end), dx(1)), ...
    nPass, nFail);

inside = cellfun(@(n) isfield(cfg.noise, n) && cfg.noise.(n) >= Snn(1) && cfg.noise.(n) <= Snn(end), cfg.R1.marked);
ok = all(inside) && isequal(cfg.R1.sea, cfg.sea);
[nPass, nFail] = check(ok, sprintf('2  %d marked noise cases, all inside the sweep (%.3g to %.3g); R1 sea = E1 sea', ...
    numel(cfg.R1.marked), min(cellfun(@(n) cfg.noise.(n), cfg.R1.marked)), max(cellfun(@(n) cfg.noise.(n), cfg.R1.marked))), ...
    nPass, nFail);

%% ---- 3-6  the E1 chain ----------------------------------------------------------------------------------------------
fprintf('  ... evaluations on the cached p0 field\n');
sL = cfg.noise.lsm6dsv16x;
oL = ev(sL, cfg.R1.sea);
e1 = evaluateScenario(struct('Snn', sL, 'cacheDir', P.cache), cfg);
ok = isequal(oL.sigma, e1.sigma) && isequal(oL.F, e1.F) && oL.dW == e1.dW && strcmp(oL.class, e1.class) && ...
    isequal(oL.validMask, e1.validMask);
[nPass, nFail] = check(ok, sprintf('3  R1 path at LSM6DSV16X = evaluateScenario default: sigma [%.4f%% %.4f%%], d_W %.4f, %s (bitwise)', ...
    100 * oL.sigma, oL.dW, oL.class), nPass, nFail);

oA = ev(Snn(1), cfg.R1.sea); oB = ev(Snn(end), cfg.R1.sea);
ok = isequal(oA.f0, oB.f0) && isequal(oA.Jall, oB.Jall) && isequal(oA.bW, oB.bW) && ~isequal(oA.F, oB.F);
[nPass, nFail] = check(ok, sprintf(['4  S_nn %.0e vs %.0e: f0, J and b_W identical; F differs (sigma_lnbeta %.4f%% vs %.3f%%, ' ...
    'd_W %.3f vs %.3f, invalid bins %.0f%% vs %.0f%%)'], Snn(1), Snn(end), 100 * oA.sigma(1), 100 * oB.sigma(1), oA.dW, oB.dW, ...
    100 * oA.invalidFraction, 100 * oB.invalidFraction), nPass, nFail);

c = 4;
sea2 = cfg.R1.sea; sea2.Hs = sea2.Hs * sqrt(c);
oC = ev(sL * c, sea2);
rs = max(abs(oC.sigma(:) ./ oL.sigma(:) - 1)); rd = abs(oC.dW / oL.dW - 1); rF = norm(oC.F - oL.F) / norm(oL.F);
ok = rs < 1e-9 && rd < 1e-9 && rF < 1e-9 && isequal(oC.validMask, oL.validMask);
[nPass, nFail] = check(ok, sprintf('5  SNR equivalence: (S_nn x %g, Hs x %g) = (S_nn, Hs): rel. diff sigma %.1e, d_W %.1e, F %.1e', ...
    c, sqrt(c), rs, rd, rF), nPass, nFail);

fields = {'p', 'sensors', 'layout', 'ref', 'estIdx', 'nuisanceIdx', 'sea', 'N', 'dt', 'L', 'overlap', 'binStep', 'bandEst', ...
    'bandUse', 'nodeSpacing', 'frfOptions', 'step'};
bad = fields(~cellfun(@(f) isequal(oA.scn.(f), e1.scn.(f)), fields));
ok = isempty(bad) && oA.scn.Snn == Snn(1);
msg = 'only Snn differs from the E1 default (p0, Level 2C, s26, sea, N, dt, L, bands, node grid, truncation)';
if ~ok, msg = [msg ': DIFFERS in ' strjoin(bad, ', ')]; end
[nPass, nFail] = check(ok, ['6  ' msg], nPass, nFail);

%% ---- 7-8  boundary refinement ---------------------------------------------------------------------------------------
[x7, i7] = r1Boundary(@(x) 10.^(2 * x), 10^(2 * -5.2345), [-6 -4], 1e-3);
ok = abs(x7 + 5.2345) <= 1e-3 && x7 >= -6 && x7 <= -4;
[nPass, nFail] = check(ok, sprintf('7  r1Boundary on a known curve: root -5.2345, found %.4f (%d evaluations)', x7, i7.nEval), ...
    nPass, nFail);

% real curve: threshold halfway (geometric) between sigma_lnbeta at 1e-6 and 1e-5, so a crossing is guaranteed
xa = -6; xb = -5;
sigB = @(x) firstSigma(ev(10^x, cfg.R1.sea));     % MATLAB cannot index a call result directly
sa = sigB(xa); sb = sigB(xb);
t = sqrt(sa * sb);
[x8, i8] = r1Boundary(sigB, t, [xa xb], 1e-2);
ok = x8 > xa && x8 < xb && (i8.fLo - t) * (i8.fHi - t) <= 0 && (i8.xHi - i8.xLo) <= 1e-2;
[nPass, nFail] = check(ok, sprintf(['8  r1Boundary on sigma_lnbeta(S_nn): threshold %.3f%% crossed at S_nn = %.4g, inside [1e-6, 1e-5]; ' ...
    'final bracket %.4f decade with values %.4f%% / %.4f%% (%d evaluations)'], 100 * t, 10^x8, i8.xHi - i8.xLo, 100 * i8.fLo, ...
    100 * i8.fHi, i8.nEval), nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testR1SNR:failed', '%d check(s) failed: do not run R1.', nFail); end

%% ================================================================================================
function [nPass, nFail] = check(ok, msg, nPass, nFail)
if ok
    fprintf('PASS  %s\n', msg); nPass = nPass + 1;
else
    fprintf('FAIL  %s\n', msg); nFail = nFail + 1;
end
end

function v = firstSigma(out)
v = out.sigma(1);
end