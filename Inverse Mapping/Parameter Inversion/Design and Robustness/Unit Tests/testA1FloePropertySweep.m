%% testA1FloePropertySweep.m
% Pre-flight checks for study A1. Must pass (all) before runA1FloePropertySweep is launched.
% The numerical core (a1Tools) is the one the runner uses, so these checks test the batch itself.
%
%   1-7   scenario construction, conversions, bands, nu and sensor placement (no EMM, seconds)
%   8     the A1 production path at p0 reproduces evaluateScenario (the frozen E1 chain) exactly
%   9     node-grid check machinery: 0.2 vs 0.1 rad/s at p0 gives a tiny d (validated earlier)
%   10    truncation gate at p0: [50 10 10] vs [70 15 15] metrics converge (sigma, d_W, class); feature bias reported
%   11    checkpoint round trip, and a changed configuration makes a checkpoint stale
%   12    a singular F is Not identifiable with the model checks NOT ASSESSED (analytic field, no EMM)
%
% COST. Checks 8-10 run the full A1 evaluation of basin_p0 once: the production field is cached from
% earlier studies; the 0.1 rad/s and [70 15 15] fields are solved the first time (several minutes) and
% cached, so the batch's basin_p0 then loads them. Nothing is written to Results/A1.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
P = setupStudyPaths();
cfg = defineStudyScenarios();
A = a1Tools();
nPass = 0; nFail = 0;
p0 = cfg.p0.vec;

%% ---- 1  counts ------------------------------------------------------------------------------------------
sc = A.buildScenarios(cfg);
g = {sc.group};
cnt = [nnz(strcmp(g, 'ring')), nnz(strcmp(g, 'corner')), nnz(~ismember(g, {'ring', 'corner'}))];
ok = isequal(cnt, [8 4 12]) && numel(unique({sc.name})) == 24 && isequal([sc.index], 1:24);
[nPass, nFail] = check(ok, sprintf('1  %d ring + %d corner + %d physical = %d scenarios, unique names', cnt, numel(sc)), ...
    nPass, nFail);

%% ---- 2  ring and corner ratios ----------------------------------------------------------------------------
ring = sc(strcmp(g, 'ring')); corner = sc(strcmp(g, 'corner'));
rb = [ring.betaRatio]; rr = [ring.RRatio];
near = @(x, y) isequal(size(x), size(y)) && max(abs(x(:) - y(:))) < 1e-12;
okRing = near(rb(1:4), cfg.A1.ringBetaRatio) && near(rr(1:4), ones(1, 4)) && near(rb(5:8), ones(1, 4)) && ...
    near(rr(5:8), cfg.A1.ringRRatio);
okCorner = near([[corner.betaRatio].', [corner.RRatio].'], cfg.A1.cornerRatios);
okGamma = all(abs([ring.gammaRatio, corner.gammaRatio] - 1) < 1e-12) && all([ring.nu, corner.nu] == 0.3);
[nPass, nFail] = check(okRing && okCorner && okGamma, ...
    '2  ring and corner (beta/beta0, R/R0) exactly as cfg.A1; gamma = gamma0, nu = 0.3', nPass, nFail);

%% ---- 3  basin_p0 round trip ---------------------------------------------------------------------------------
b0 = sc(strcmp({sc.name}, 'basin_p0'));
ok = max(abs(b0.p ./ p0 - 1)) < 1e-12 && isempty(b0.frfOptions) && b0.scale == 1;
[nPass, nFail] = check(ok, sprintf('3  basin_p0 -> p = %s (p0 to %.1e), nu = 0.3 passes nothing to the EMM', ...
    mat2str(b0.p, 5), max(abs(b0.p ./ p0 - 1))), nPass, nFail);

%% ---- 4  polypropylene p0 analogue ------------------------------------------------------------------------------
q = sc(strcmp({sc.name}, 'pp3mm_R072'));
ok = abs(q.betaRatio - 0.698) < 2e-3 && abs(q.gammaRatio - 1) < 0.01 && abs(q.RRatio - 0.72 / cfg.p0.Rphys) < 1e-12;
[nPass, nFail] = check(ok, sprintf('4  pp3mm_R072: beta = %.3f beta0, gamma = %.3f gamma0, R = %.4f R0 (expected 0.70, ~1, 1)', ...
    q.betaRatio, q.gammaRatio, q.RRatio), nPass, nFail);

%% ---- 5  finite-depth bands ----------------------------------------------------------------------------------------
fin = sc(~[sc.deepWater]);
okB = all(arrayfun(@(s) isequal(s.band, cfg.band.basin) && s.Hphys == cfg.physics.Hbasin && s.scale == 1, fin));
aOk = true;
for s = fin
    w = welchBins(cfg.welch.L, cfg.acq.dt, cfg.band.estimate, cfg.welch.binStep, s.band);
    G = A.dimensionlessGroups(s, w, cfg);
    aOk = aOk && G.alphaValid && numel(w) == 84;
end
[nPass, nFail] = check(okB && aOk, sprintf(['5  %d finite-depth cases: depth 1.88 m, band exactly %s rad/s, 84 bins, ' ...
    'alpha inside [%.2f %.2f]'], numel(fin), mat2str(cfg.band.basin), cfg.emm.alphaValidated), nPass, nFail);

%% ---- 6  MIZ: effective depth, Froude-scaled band ----------------------------------------------------------------
miz = sc([sc.deepWater]);
ok = numel(miz) == 4; txt = {};
for s = miz
    T = cfg.A1.bandRule.mizPeriods;
    scr = screenFloeRegime(struct('R', s.Rphys, 'h', s.h, 'E', s.E, 'nu', s.nu, 'rhoIce', s.rhoI), T, ...
        'Depth', Inf, 'RhoWater', s.rhoW);
    kMin = min(scr.k(2), scr.kappa(2));                       % at the longest period, 12.5 s
    okDepth = abs(kMin * s.Hphys - pi) < 1e-9;
    okP = max(abs(s.p ./ [scr.emm.beta scr.emm.gamma scr.emm.Rnd] - 1)) < 1e-12;
    % the band, mapped back to full-scale periods, must lie inside 6-12.5 s and inside the basin band
    okBand = s.band(1) >= cfg.band.basin(1) - 1e-12 && s.band(2) <= cfg.band.basin(2) + 1e-12 && ...
        s.bandPeriods(1) >= min(T) - 1e-9 && s.bandPeriods(2) <= max(T) + 1e-9 && ...
        abs(s.scale - s.Hphys / cfg.physics.Hbasin) < 1e-12;
    w = welchBins(cfg.welch.L, cfg.acq.dt, cfg.band.estimate, cfg.welch.binStep, s.band);
    G = A.dimensionlessGroups(s, w, cfg);
    % Froude invariance: alpha of the retained bins at the basin equals alpha at H_eff of the full-scale omega
    wFull = w / sqrt(s.scale);
    okAlpha = max(abs(s.Hphys * wFull.^2 / cfg.physics.g - cfg.physics.Hbasin * w.^2 / cfg.physics.g)) < 1e-9 && G.alphaValid;
    ok = ok && okDepth && okP && okBand && okAlpha;
    txt{end + 1} = sprintf('%s H_eff %.0f m, %.1f-%.1f s, %d bins', s.name, s.Hphys, s.bandPeriods, numel(w)); %#ok<SAGROW>
end
[nPass, nFail] = check(ok, ['6  MIZ: min(k,kappa) H_eff = pi at 12.5 s, p from screenFloeRegime, band inside 6-12.5 s ' ...
    'and the basin band, alpha Froude-invariant: ' strjoin(txt, '; ')], nPass, nFail);

%% ---- 7  nu and sensor placement (analytic field, no EMM) -----------------------------------------------------------
pp = sc(strcmp({sc.group}, 'tank'));
okNu = all(arrayfun(@(s) isequal(s.frfOptions, {'Nu', 0.4}), pp)) && ...
    all(arrayfun(@(s) isempty(s.frfOptions), sc(abs([sc.nu] - 0.3) < 1e-12)));
s = sc(strcmp({sc.name}, 'ring_R_x2'));
scn = A.buildScn(s, cfg, struct());
okPos = max(max(abs(scn.sensors - cfg.layout.level2C.frac .* [s.p(3) 1]))) == 0;
% analytic field H_j = 1 + r_j / R: at FIXED points dH/dlnR = -r_j / R; if the points moved with R
% (r = xi R) the field would not depend on R at all and the derivative would be 0
frf = @(w, pp) (1 + scn.sensors(:, 1) / pp(3)) * ones(1, numel(w));
fld = emmField(scn.p, scn.sensors, 'ThetaIdx', [1 3], 'FRF', frf);
w = linspace(3.5, 8, 5);
dH = fld.dH(w);
okFix = max(max(abs(dH(:, :, 2) - (-scn.sensors(:, 1) / s.p(3)) * ones(1, 5)))) < 1e-3 && max(max(abs(dH(:, :, 1)))) < 1e-12;
[nPass, nFail] = check(okNu && okPos && okFix, sprintf(['7  nu = 0.4 passed to the EMM for the %d polypropylene/model-floe ' ...
    'cases only; sensors = Level 2C fractions x true R (ring_R_x2); dH/dlnR at fixed points = -r/R (max error %.1e)'], ...
    numel(pp), max(max(abs(dH(:, :, 2) - (-scn.sensors(:, 1) / s.p(3)) * ones(1, 5))))), nPass, nFail);

%% ---- 8-10  the A1 path at p0, with both model checks -------------------------------------------------------------------
fprintf('  ... basin_p0 through the full A1 path (0.1 rad/s and [70 15 15] fields are solved the first time: minutes)\n');
Snn = cfg.noise.(cfg.A1.noiseCase);
tic; r = A.evaluate(b0, cfg, struct('cacheDir', P.cache, 'verbose', true)); tA1 = toc;
e1 = evaluateScenario(struct('Snn', Snn, 'cacheDir', P.cache), cfg);
rs = max(abs(r.sigma(:) ./ e1.sigma(:) - 1)); rd = abs(r.dW / e1.dW - 1); rF = norm(r.F - e1.F) / norm(e1.F);
ok = rs < 1e-10 && rd < 1e-10 && rF < 1e-10 && strcmp(r.inverseClass, e1.class) && r.nodeSpacingUsed == cfg.emm.nodeSpacing;
[nPass, nFail] = check(ok, sprintf(['8  A1 path at p0 = evaluateScenario: sigma [%.4f%% %.4f%%], d_W %.4f, %s ' ...
    '(rel. diff sigma %.1e, d_W %.1e, F %.1e; %.0f s)'], 100 * r.sigma, r.dW, r.inverseClass, rs, rd, rF, tA1), nPass, nFail);

ok = r.nodeGrid.done && isfinite(r.nodeGrid.d) && r.nodeGrid.d < 0.01 && r.nodeGrid.relT < 1e-3 && ~r.nodeGrid.fallback;
[nPass, nFail] = check(ok, sprintf('9  node grid 0.2 vs 0.1 rad/s at p0: d = %.2e (< 0.01), max |dT|/|T| = %.2e, no fallback', ...
    r.nodeGrid.d, r.nodeGrid.relT), nPass, nFail);

ok = r.truncation.done && r.truncationValid && r.truncation.relSigma < 1e-3 && r.truncation.deltaDW < 1e-3 && ...
    r.truncation.sameClass && r.modelValid && strcmp(r.finalClass, r.inverseClass) && isfinite(r.truncation.d);
[nPass, nFail] = check(ok, sprintf(['10 truncation [50 10 10] vs [70 15 15] at p0, metric gate: sigma %.1e rel, d_W %.1e abs, ' ...
    'class %s (gate %.0f%%, %.2g); feature bias d = %.3f reported (diagnostic: 0.57); model valid = %d, final class %s'], ...
    r.truncation.relSigma, r.truncation.deltaDW, r.truncation.classHi, 100 * cfg.A1.truncationGate.maxRelSigma, ...
    cfg.A1.truncationGate.maxDeltaDW, r.truncation.d, r.modelValid, r.finalClass), nPass, nFail);

%% ---- 11  checkpoints ------------------------------------------------------------------------------------------------------
tmp = fullfile(tempdir, sprintf('a1test_%s', datestr(now, 'HHMMSSFFF')));
key = A.checkpointKey(b0, cfg);
file = A.checkpointPath(tmp, b0);
A.saveCheckpoint(file, r, key);
[r2, ok1] = A.loadCheckpoint(file, key);
cfg2 = cfg; cfg2.emm.truncationCheck = [80 15 15];
[~, ok2] = A.loadCheckpoint(file, A.checkpointKey(b0, cfg2));
[~, ok3] = A.loadCheckpoint(fullfile(tmp, 'missing.mat'), key);
ok = ok1 && isequaln(r2.sigma, r.sigma) && ~ok2 && ~ok3;
[nPass, nFail] = check(ok, '11 checkpoint saves and reloads; a changed configuration makes it stale; a missing file is not loaded', ...
    nPass, nFail);
delete(file); rmdir(tmp);

%% ---- 12  singular F: checks not assessable, never a silent pass (analytic field, no EMM) --------------------
% a field that does not depend on beta at all: the ln beta column of J is zero, so F is singular
s = sc(strcmp({sc.name}, 'ring_R_x2'));
pos = A.placeSensors(s, cfg); rr = pos(:, 1); tt = pos(:, 2);
frfSing = @(w, pp) -(ones(numel(rr), 1) * w.^2) .* (1 + 0.5 * (rr / pp(3)).^2 * ones(1, numel(w))) .* ...
    exp(-0.3i * (rr .* cos(tt)) * w / pp(3)^0.2);
rS = A.evaluate(s, cfg, struct('frf', frfSing));
ok = rS.singular && ~rS.modelAssessed && ~rS.modelValid && ~rS.nodeGrid.done && ~rS.truncation.done && ...
    ~rS.nodeGridValid && ~rS.truncationValid && strcmp(rS.finalClass, 'Not identifiable') && ...
    ~isempty(strfind(rS.modelReason, 'not assessable'));
[nPass, nFail] = check(ok, sprintf(['12 singular F (field independent of beta): final class %s, model checks not run and ' ...
    'not recorded as passed (assessed %d, valid %d; "%s")'], rS.finalClass, rS.modelAssessed, rS.modelValid, ...
    rS.modelReason), nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testA1FloePropertySweep:failed', '%d check(s) failed: do not launch the A1 batch.', nFail); end

%% ================================================================================================
function [nPass, nFail] = check(ok, msg, nPass, nFail)
if ok
    fprintf('PASS  %s\n', msg); nPass = nPass + 1;
else
    fprintf('FAIL  %s\n', msg); nFail = nFail + 1;
end
end