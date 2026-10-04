%% testV1ModelMismatch.m
% Pre-flight checks for study V1. All must pass before runV1ModelMismatch is run. Uses cached fields where they
% exist; checks 4 to 6 and 10 each need one EMM field (29 nodes, about 30 s); check 11 uses an analytic field.
%
%   1  cfg.V1: registered values unchanged (inverse [50 10 10], depth +-2/5%, radial +-1/2 cm), clarification logged
%   2  stage 0: the newest truncation ladder selected cfg.V1.truncationTruth by the registered rule (consecutive
%      d < 0.1 at p0); its P check is below 0.1 or a logged P amendment is in place; truncationTruthSelected set
%   3  the p0 nominal scenario reproduces E1/R1 (sigma, d_W, class at LSM6DSV16X)
%   4  zero perturbation: the nominal model through v1Tools gives f_W = engine f0 + b_W, b_sys = 0, d_sys = 0
%   5  nominal truncation as "truth" ([50 10 10] passed explicitly): zero discrepancy
%   6  depth eps = 0 ('WaterDepth' 1.88 passed explicitly): zero discrepancy
%   7  depth scaling: beta ~ H^-4, gamma ~ H^-1, R ~ H^-1, sensors' physical radius and angle fixed; alpha range
%      at +5% (top retained bin about 13.87, flagged) and at 0 (inside, not flagged)
%   8  radial +1 cm: every physical radius moves by exactly 0.01 m, angles unchanged
%   9  the projector equals the engine's: dtheta_sys, d_sys and d_tot of an arbitrary feature error equal
%      evaluateScenario with scn.sysBias
%  10  one real case (radial +1 cm): b_sys is Welch minus Welch; the inverse (nominal out) is unchanged by the
%      evaluation; the V1 class equals evaluateScenario's class with that b_sys
%  11  nonlinear check (analytic field): on the nominal Welch features the G3 solver from theta_true converges
%      to theta_true + dth_W within 0.05 sigma
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
P = setupStudyPaths();
cfg = defineStudyScenarios();
T = v1Tools(); V = cfg.V1;
nPass = 0; nFail = 0;
rel = @(a, b) max(abs(a(:) - b(:))) / max(abs(b(:)));

%% ---- 1  configuration ----------------------------------------------------------------------------------------------
ok = isequal(V.truncationInverse, [50 10 10]) && isequal(V.depthErrorFrac, [-0.05 -0.02 0.02 0.05]) && ...
    isequal(V.sensorRadialError, [-0.02 -0.01 0.01 0.02]) && ...
    any(arrayfun(@(c) ~isempty(strfind(c.field, 'V1.truncationTruthSelected')), cfg.meta.changes)) && ...
    any(arrayfun(@(c) ~isempty(strfind(c.field, 'V1.truncationLadderM')), cfg.meta.changes));
[nPass, nFail] = check(ok, sprintf(['1  cfg.V1: inverse [%d %d %d], depth %s, radial %s m; ladder and study ' ...
    'clarifications logged'], V.truncationInverse, mat2str(V.depthErrorFrac), mat2str(V.sensorRadialError)), nPass, nFail);

%% ---- 2  stage 0 selection -------------------------------------------------------------------------------------------
Lst = dir(fullfile(P.results, 'V1', 'diagnostics', 'truncationLadder_*.mat'));
ok = false; msg = '2  stage 0: no truncation ladder found (run diagnoseV1Truncation first)';
if ~isempty(Lst)
    [~, k] = max([Lst.datenum]); S = load(fullfile(Lst(k).folder, Lst(k).name));
    if isfield(S, 'ladder')
        Ld = S.ladder; dSel = Ld.levels(Ld.chosen).d;
        pOk = Ld.pCheck.d < V.truncationConvergedD || V.truncationTruth(2) >= Ld.pCheck.to(2);
        ok = V.truncationTruthSelected && Ld.converged && isequal(Ld.truth, V.truncationTruth) && ...
            dSel < V.truncationConvergedD && pOk;
        msg = sprintf(['2  stage 0 (%s): truth [%d %d %d] = cfg %s, last step d %.4f (< %.2g %s), P check d %.4f %s; ' ...
            'large-kR last step d %s; selected flag %d'], Lst(k).name, Ld.truth, tf(isequal(Ld.truth, V.truncationTruth)), ...
            dSel, V.truncationConvergedD, tf(Ld.converged), Ld.pCheck.d, tf(pOk), mat2str(Ld.dLastStep(2:end), 3), ...
            V.truncationTruthSelected);
    else
        msg = sprintf('2  stage 0: %s holds a partial ladder (the run did not finish)', Lst(k).name);
    end
end
[nPass, nFail] = check(ok, msg, nPass, nFail);

%% ---- 3  nominal p0 reproduces E1 / R1 --------------------------------------------------------------------------------
fprintf('  ... nominal p0 (cached)\n');
out0 = evaluateScenario(struct('cacheDir', P.cache), cfg);
cs = resolveM1Scenarios(cfg, P); src = cs(4).source;
ok = max(abs(out0.sigma(:).' ./ src.sigma - 1)) < 1e-12 && abs(out0.dW - src.dW) < 1e-12 && strcmp(out0.class, src.class);
[nPass, nFail] = check(ok, sprintf('3  nominal p0: sigma [%.4f%% %.4f%%], d_W %.4f, %s = R1 source (%s)', 100 * out0.sigma, ...
    out0.dW, out0.class, src.class), nPass, nFail);

%% ---- 4-6  zero discrepancies -------------------------------------------------------------------------------------------
b0 = T.base(out0);
fprintf('  ... nominal Welch features through v1Tools (one EMM field, about 30 s the first time)\n');
fW0 = T.welchFeatures(T.truthModel(b0, 'none', [], cfg), out0, P.cache);
r0 = T.evaluate(b0, out0, fW0, 'none', [], cfg, P.cache);
e4 = rel(fW0, out0.f0(:) + out0.bW(:));
ok = e4 < 1e-10 && all(r0.bSys == 0) && r0.dSys == 0 && r0.betaPct == 0 && strcmp(r0.class, out0.class) && ...
    abs(r0.dTot - out0.dW) < 1e-12 && b0.H == cfg.physics.Hbasin;
[nPass, nFail] = check(ok, sprintf(['4  zero perturbation: f_W(v1Tools) = f0 + b_W to %.1e; b_sys = 0, d_sys = %g, ' ...
    'd_tot = d_W, class %s, H0 = %.2f m'], e4, r0.dSys, r0.class, b0.H), nPass, nFail);

fprintf('  ... [50 10 10] passed explicitly as the truth truncation (one EMM field)\n');
r5 = T.evaluate(b0, out0, fW0, 'truncation', V.truncationInverse, cfg, P.cache);
ok = max(abs(r5.bSys)) <= 1e-12 * max(abs(fW0)) && r5.dSys < 1e-9;
[nPass, nFail] = check(ok, sprintf('5  nominal truncation as truth: max |b_sys| %.1e (relative), d_sys %.1e', ...
    max(abs(r5.bSys)) / max(abs(fW0)), r5.dSys), nPass, nFail);

fprintf('  ... depth eps = 0 with WaterDepth passed explicitly (one EMM field)\n');
m6 = T.truthModel(b0, 'depth', 0, cfg);
r6 = T.evaluate(b0, out0, fW0, 'depth', 0, cfg, P.cache);
ok = isequal(m6.p, b0.p) && isequal(m6.sensors, b0.sensors) && max(abs(r6.bSys)) <= 1e-12 * max(abs(fW0)) && ...
    r6.dSys < 1e-9 && ~r6.alpha.extrapolated && ~r6.alpha.retainedOutside;
[nPass, nFail] = check(ok, sprintf(['6  depth eps = 0: p and sensors unchanged, max |b_sys| %.1e, d_sys %.1e; alpha ' ...
    'support [%.3f %.3f] inside'], max(abs(r6.bSys)) / max(abs(fW0)), r6.dSys, r6.alpha.minSupport, ...
    r6.alpha.maxSupport), nPass, nFail);

%% ---- 7  depth scaling and alpha -----------------------------------------------------------------------------------------
e = 0.05; m7 = T.truthModel(b0, 'depth', e, cfg); s = 1 + e;
law = abs(m7.p ./ b0.p - [s^-4, s^-1, s^-1]);
phys = abs(m7.sensors(:, 1) * m7.H - b0.sensors(:, 1) * b0.H);
a7 = T.alphaRange(m7.H, out0.omega(out0.validMask), m7.band, cfg);
ok = max(law) < 1e-14 && max(phys) < 1e-14 && isequal(m7.sensors(:, 2), b0.sensors(:, 2)) && ...
    abs(m7.H - 1.88 * s) < 1e-14 && isequal(m7.band, b0.band) && abs(a7.maxRetained - 13.866) < 0.005 && ...
    a7.extrapolated && any(strcmp(m7.frfOptions, 'WaterDepth'));
[nPass, nFail] = check(ok, sprintf(['7  depth +5%%: p ratios [%.6f %.6f %.6f] (laws H^-4, H^-1, H^-1), physical radii ' ...
    'fixed (%.1e), band unchanged; alpha retained max %.3f, support [%.3f %.3f], flagged'], m7.p ./ b0.p, max(phys), ...
    a7.maxRetained, a7.minSupport, a7.maxSupport), nPass, nFail);

%% ---- 8  radial placement -----------------------------------------------------------------------------------------------
m8 = T.truthModel(b0, 'radial', 0.01, cfg);
dr = (m8.sensors(:, 1) - b0.sensors(:, 1)) * b0.H;
ok = max(abs(dr - 0.01)) < 1e-14 && isequal(m8.sensors(:, 2), b0.sensors(:, 2)) && isequal(m8.p, b0.p) && ...
    all(b0.sensors(:, 1) > 0);
[nPass, nFail] = check(ok, sprintf('8  radial +1 cm: physical radii change by %s m, angles unchanged, floe unchanged', ...
    mat2str(dr.', 6)), nPass, nFail);

%% ---- 9  projector = engine --------------------------------------------------------------------------------------------
bArb = out0.Jall * [0.004; -0.0007] + 0.1 * out0.bW;
q9 = T.project(out0, bArb); qt9 = T.project(out0, out0.bW + bArb);
oS = evaluateScenario(struct('cacheDir', P.cache, 'sysBias', bArb), cfg);
ok = max(abs(q9.dtheta - oS.dthetaSys)) < 1e-14 && abs(q9.d - oS.dSys) < 1e-12 && abs(qt9.d - oS.dTot) < 1e-12 && ...
    strcmp(T.classify(out0, qt9.d, cfg), oS.class);
[nPass, nFail] = check(ok, sprintf(['9  projector = engine: dtheta_sys [%+.5f %+.5f], d_sys %.4f, d_tot %.4f, class %s ' ...
    '(engine %.4f, %.4f, %s)'], q9.dtheta, q9.d, qt9.d, T.classify(out0, qt9.d, cfg), oS.dSys, oS.dTot, oS.class), ...
    nPass, nFail);

%% ---- 10  one real case -----------------------------------------------------------------------------------------------
fprintf('  ... radial +1 cm truth (one EMM field)\n');
F0 = out0.F; s0 = out0.sigma;
r10 = T.evaluate(b0, out0, fW0, 'radial', 0.01, cfg, P.cache);
fW10 = T.welchFeatures(T.truthModel(b0, 'radial', 0.01, cfg), out0, P.cache);
o10 = evaluateScenario(struct('cacheDir', P.cache, 'sysBias', r10.bSys), cfg);
ok = isequal(r10.bSys, fW10 - fW0) && isequal(out0.F, F0) && isequal(out0.sigma, s0) && ...
    strcmp(r10.class, o10.class) && abs(r10.dSys - o10.dSys) < 1e-12 && abs(r10.dTot - o10.dTot) < 1e-12;
[nPass, nFail] = check(ok, sprintf(['10 radial +1 cm: b_sys = f_W(truth) - f_W(nominal); d_sys %.3f, dbeta/beta %+.3f%%, ' ...
    'd_tot %.3f, %s = engine with that b_sys; nominal inverse unchanged'], r10.dSys, r10.betaPct, r10.dTot, r10.class), ...
    nPass, nFail);

%% ---- 11  nonlinear check (analytic field) ---------------------------------------------------------------------------
oA = evaluateScenario(struct('Snn', cfg.noise.lsm6dsv16x, 'frf', analyticFRF(cfg.p0.vec)), cfg);
thLin = log(oA.scn.p([1 3])).' + oA.dthetaW(:);
nl = T.nonlinear(oA, oA.f0(:) + oA.bW(:), thLin, cfg);
ok = nl.converged && max(abs(nl.gapSigma)) < 0.05;
[nPass, nFail] = check(ok, sprintf(['11 nonlinear (analytic field): %s in %d it; NL - (theta_true + dth_W) = [%+.4f %+.4f] ' ...
    'sigma (d_W %.3f)'], tf(nl.converged), nl.it, nl.gapSigma, oA.dW), nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testV1ModelMismatch:failed', '%d check(s) failed: do not run V1.', nFail); end

%% ================================================================================================
function [nPass, nFail] = check(ok, msg, nPass, nFail)
if ok
    fprintf('PASS  %s\n', msg); nPass = nPass + 1;
else
    fprintf('FAIL  %s\n', msg); nFail = nFail + 1;
end
end

function s = tf(b)
if b, s = 'ok'; else, s = 'NO'; end
end

function frf = analyticFRF(p0)
% smooth, non-degenerate test field (as testM1): acceleration FRF at the Level 2C points depending on beta and R
frf = @(w, pp) fieldAt(w, pp, p0);
end

function H = fieldAt(w, pp, p0)
persistent sens
if isempty(sens)
    cfg = defineStudyScenarios(); sens = cfg.layout.level2C.frac .* [cfg.p0.R 1];
end
w = w(:).'; b = (pp(1) / p0(1))^0.35; x = sens(:, 1) / pp(3);
H = -(ones(size(sens, 1), 1) * w.^2) .* (1 + 0.6 * (x.^2 .* cos(sens(:, 2))) * ones(1, numel(w)) .* ...
    exp(1i * (0.9 * b * x) * w) + 0.25 * b * (x * w));
end