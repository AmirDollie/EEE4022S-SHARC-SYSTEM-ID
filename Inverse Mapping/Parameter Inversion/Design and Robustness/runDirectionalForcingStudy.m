%% runDirectionalForcingStudy.m
% V2: directional forcing. The data come from a sea with a second, independent incident direction; the inverse keeps
% the frozen single-direction model. How large is the parameter bias, and does coherence warn first?
%
% WHY. The output-only inverse rests on the single-input identity S_jr / S_rr = H_j / H_r (rank-one spectral matrix).
% With two independent components y = H_1 u_1 + H_2 u_2 the matrix is rank two and the estimated transmissibility is
% an energy-weighted mixture, which the inverse reads as a change of beta and R. V2 measures that bias and asks
% whether the coherence drop, which is observable, reveals the failed assumption before the bias matters.
%
% METHOD (programme rev. 7 and cfg.V2; clarification logged 2026-10-05 before the run)
%   Forcing  primary component at 0 deg (the EMM incidence, the Level 2C axis), secondary at +phi, both with the
%            reference JONSWAP shape; total energy fixed: S_u1 = (1 - f2) S_u, S_u2 = f2 S_u. A component from
%            heading phi is the nominal field at the sensors rotated by -phi (circular floe): ONE EMM field at the
%            rotated Level 2C points serves every case. S = sum_c f_c H_c H_c^H S_u (Welch-smoothed, engine kernel)
%            + diag(S_nn): directionalSpectralMatrix, shared with D3.
%   V2a bias b_dir = f_W(directional) - f_W(single), projected with the NOMINAL J, Sigma, F on the nominal accepted
%            bins; d_sys, dbeta/beta, dR/R; d_tot from the vector sum with dth_W; class. phi in cfg.V2.angleDeg,
%            f2 in cfg.V2.secondEnergyFraction (4 x 4).
%   V2b coherence: expected MSC per pair, single-direction (noisy) baseline vs the directional sea; Fisher z; the drop
%            of the bin-mean D_j over the accepted bins, in units of its null scatter, CALIBRATED by Monte Carlo of the
%            real estimator on single-direction records (bins not assumed independent; amendment 2026-10-05);
%            detectable when max_j D_j >= z(1 - 0.05/m) + z(0.80) (80%-power detectability threshold).
%   Local regime: d_sys / f2 at f2 = cfg.V2.linearityF2 per phi (constant when b_dir ~ f2).
%   Crossings per phi (first, in f2, bisection in log10 f2 to 0.005 decade, not beyond 0.5): d_sys = 0.5, 1,
%            |dbeta/beta| = 1%, and detection; boundary refinements inside the bracket between the exact f2 = 0 baseline
%            and the first registered point. Status exact | below_floor | not_crossed. EARLY WARNING at phi if
%            f_detect < f(d_sys = 0.5); unresolved if both lie below the floor.
%   Nonlinear diagnostic: at the d_sys = 1 crossing of the most sensitive angle (else the registered case with d_sys
%            closest to 1), G3 solver from theta_true on its noiseless directional Welch features, against
%            theta_true + dth_W + dth_dir.
%   Diagnostics (not criteria): sigma the estimator would have under the directional sea (its larger scatter), and
%            how many bins the estimator's acceptance rule would change.
%
% PRECONDITIONS. testDirectionalForcingStudy passes (12/12); no uncommitted .m files.
% COST. One EMM field at 5 x 4 rotated points (29 nodes, about 40 s); the coherence null Monte Carlo (400 records,
% about 1 min); linear algebra (seconds); the nonlinear diagnostic about 5 to 15 min (one EMM field per evaluation).
% OUTPUT. Results/V2/ (log, .mat, summary, V2_grid.csv, V2_crossings.csv, V2_nonlinear.csv). Figures:
% plotDirectionalForcingStudy (after the run).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

ALLOW_DIRTY = false;    % true: run with uncommitted .m files (provenance is then NOT the logged commit)

P = setupStudyPaths();
cfg = defineStudyScenarios();
T = v2Tools(); T1 = v1Tools(); Tm = m1Tools();
V = cfg.V2;
run = startStudy('V2', 'Directional forcing', P.results);

try
%% ---- 0  provenance ---------------------------------------------------------------------------------------------------
[st0, gitOut] = system(sprintf('git -C "%s" status --porcelain', P.root));
if st0 ~= 0
    fprintf('  WARNING: git status failed (%s); the code state cannot be checked.\n', strtrim(gitOut));
else
    [clean, dirty] = Tm.cleanTree(gitOut);
    if ~clean
        fprintf('  uncommitted .m files:\n'); fprintf('    %s\n', dirty{:});
        if ~ALLOW_DIRTY
            error('V2:dirtyTree', '%d uncommitted .m file(s): commit first (or set ALLOW_DIRTY = true).', numel(dirty));
        end
        fprintf('  *** ALLOW_DIRTY: this run is NOT reproducible from commit %s ***\n', run.commit);
    else
        fprintf('  working tree clean (.m files): the code that runs is commit %s\n', run.commit);
    end
end
phis = V.angleDeg * pi / 180; f2s = V.secondEnergyFraction;
fprintf('  operating point: %s\n  phi = %s deg, f2 = %s (total energy fixed)\n', V.operatingPoint, ...
    mat2str(V.angleDeg), mat2str(f2s));

%% ---- 1  nominal p0 and the rotated-point field --------------------------------------------------------------------------
out0 = evaluateScenario(struct('cacheDir', P.cache), cfg);
fprintf('\nStep 1: nominal p0: sigma [%.4f%% %.4f%%], d_W %.4f, %s; %d accepted bins, N_eff %.1f\n', 100 * out0.sigma, ...
    out0.dW, out0.class, sum(out0.validMask), out0.nEff);
ctx = T.context(out0, [V.primaryHeadingDeg, V.angleDeg] * pi / 180, cfg, P.cache);
zc = [max(abs(ctx.fW0 - (out0.f0(:) + out0.bW(:)))) / max(abs(ctx.fW0)), ...
    max(abs(ctx.coh0(:) - out0.coherence(:))), max(abs(ctx.blocks0(:) - out0.blocksAll(:))) / max(abs(out0.blocksAll(:)))];
fprintf(['  field at %d rotated points (%.0f s): the single direction through the directional path reproduces the ' ...
    'engine: f_W %.1e, coherence %.1e, Sigma %.1e\n'], size(ctx.pts.sensors, 1), ctx.fieldSeconds, zc);
if any(zc > 1e-10) || ~isequal(ctx.valid0, out0.validMask)
    error('V2:baseline', 'The single-direction baseline does not reproduce the engine (%.2e %.2e %.2e).', zc);
end
m = size(ctx.coh0, 1);
cohB = mean(ctx.coh0(:, ctx.valid), 2);
fprintf(['  single-direction coherence (mean over accepted bins) per pair: %s; 80%%-power detectability threshold ' ...
    'D >= %.3f (m = %d; the 5%% rejection threshold is %.3f)\n'], mat2str(cohB.', 5), T.normInv(1 - V.coherence.alpha / m) + ...
    T.normInv(V.coherence.power), m, T.normInv(1 - V.coherence.alpha / m));
fprintf('\nStep 1b: coherence null Monte Carlo (%d single-direction records, the real estimator)\n', V.coherence.nullRecords);
nc = T.nullCoherence(out0, cfg, V.coherence.nullRecords, V.coherence.nullSeedBase);
ctx.sdNull = nc.sd;
fprintf(['  per pair: sd of the bin-mean Fisher z %s (analytic, independent bins %.4g); variance ratio %s; effective ' ...
    'independent bins %s of %d; mean offset %s sd (%.0f s)\n'], mat2str(nc.sd.', 4), nc.sdAnalytic, ...
    mat2str(nc.varianceRatio.', 3), mat2str(round(nc.nBinsEff.'), 4), nc.nBins, mat2str((nc.mean ./ nc.sd).', 3), nc.seconds);

%% ---- 2  the registered grid -------------------------------------------------------------------------------------------
fprintf('\nStep 2: %d x %d grid (phi x f2); class under the nominal sigma ruler\n', numel(phis), numel(f2s));
G = struct('row', {});
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
evalAt = @(phi, f2) cachedCase(cache, T, ctx, out0, phi, f2, cfg);
for i = 1:numel(phis)
    for f2 = f2s
        r = evalAt(phis(i), f2);
        G(end + 1).row = r; %#ok<SAGROW>
        fprintf(['  phi %2.0f  f2 %.2f  d_sys %8.3f  dbeta/beta %+8.3f%%  dR/R %+8.4f%%  d_tot %8.3f  %-16s  D_max %7.2f ' ...
            '%-3s  sigma x [%.3f %.3f]%s\n'], r.phiDeg, f2, r.dSys, r.betaPct, r.RPct, r.dTot, r.class, r.Dmax, ...
            ternary(r.detect, 'DET', '-'), r.sigmaRatio, ternary(r.nValidChanged > 0, sprintf('  (%d bins change acceptance)', ...
            r.nValidChanged), ''));
    end
end

%% ---- 2b  local regime ------------------------------------------------------------------------------------------------------
fprintf('\nStep 2b: local regime, d_sys / f2 at f2 = %s (constant if b_dir ~ f2)\n', mat2str(V.linearityF2));
LIN = struct('phiDeg', {}, 'f2', {}, 'dPerF2', {}, 'betaPerF2', {}, 'spread', {});
for i = 1:numel(phis)
    Li = T.linearity(ctx, out0, phis(i), V.linearityF2, cfg);
    LIN(end + 1) = Li; %#ok<SAGROW>
    r05 = evalAt(phis(i), min(f2s));
    fprintf('  phi %2.0f: d_sys/f2 %s (spread %.2f%%); at f2 = %g: %.4g (%.2f x the smallest-f2 value)\n', ...
        V.angleDeg(i), mat2str(Li.dPerF2, 5), 100 * Li.spread, min(f2s), r05.dSys / min(f2s), ...
        r05.dSys / min(f2s) / Li.dPerF2(1));
end

%% ---- 3  crossings in f2 per phi -------------------------------------------------------------------------------------------
fprintf('\nStep 3: first crossing in f2 per phi (bisection in log10 f2 to %.3g decade)\n', V.toleranceTolDecades);
crit = struct('name', {'d_sys', 'd_sys', '|dbeta/beta|', 'detection'}, 'field', {'dSys', 'dSys', 'absBeta', 'detect'}, ...
    'threshold', {V.thresholds.dSys(1), V.thresholds.dSys(2), V.thresholds.betaPct, NaN});
X = struct('phiDeg', {}, 'criterion', {}, 'threshold', {}, 'status', {}, 'crossed', {}, 'f2', {}, 'nEval', {}, ...
    'note', {});
floorOf = min(f2s) * V.toleranceFloorFraction;
for i = 1:numel(phis)
    xs = [0, sort(f2s)];
    for c = 1:numel(crit)
        thr = crit(c).threshold;
        if strcmp(crit(c).field, 'detect'), thr = 0; end            % metric D_max - threshold crosses 0
        met = @(f2) metricOf(evalAt(phis(i), f2), crit(c).field);
        vals = arrayfun(met, xs);
        k = find(vals(1:end - 1) < thr & vals(2:end) >= thr, 1);
        t = struct('phiDeg', V.angleDeg(i), 'criterion', crit(c).name, 'threshold', crit(c).threshold, ...
            'status', 'not_crossed', 'crossed', false, 'f2', NaN, 'nEval', 0, 'note', '');
        if isempty(k)
            t.note = sprintf('not crossed within f2 <= %g', xs(end));
        else
            xa = xs(k); xb = xs(k + 1); nExtra = 0;
            if xa == 0
                xa = floorOf; nExtra = 1;
                if met(xa) >= thr
                    t.status = 'below_floor'; t.nEval = nExtra;
                    t.note = sprintf('below the floor: already met at f2 = %.3g', xa);
                end
            end
            if ~strcmp(t.status, 'below_floor')
                [ub, info] = r1Boundary(@(u) met(10^u), thr, log10([xa xb]), V.toleranceTolDecades);
                t.status = 'exact'; t.crossed = true; t.f2 = 10^ub; t.nEval = info.nEval + nExtra;
                t.note = sprintf('f2 = %.4g', t.f2);
            end
        end
        X(end + 1) = t; %#ok<SAGROW>
        fprintf('  phi %2.0f  %-13s %-4s: %s\n', V.angleDeg(i), crit(c).name, numStr(crit(c).threshold), t.note);
    end
end
fprintf('\n  Early warning (coherence detectable before d_sys = %.1f):\n', V.thresholds.dSys(1));
warnTab = struct('phiDeg', num2cell(V.angleDeg), 'fDetect', NaN, 'fBias', NaN, 'early', false, 'verdict', '');
for i = 1:numel(phis)
    xd = X([X.phiDeg] == V.angleDeg(i) & strcmp({X.criterion}, 'detection'));
    xb = X([X.phiDeg] == V.angleDeg(i) & strcmp({X.criterion}, 'd_sys') & [X.threshold] == V.thresholds.dSys(1));
    sd = xd.status; sb = xb.status; fd = xd.f2; fb = xb.f2; early = false;
    if strcmp(sd, 'below_floor') && strcmp(sb, 'below_floor')
        v = sprintf('ordering unresolved (both below f2 = %.3g)', floorOf);
    elseif strcmp(sd, 'below_floor') || (strcmp(sd, 'exact') && (strcmp(sb, 'not_crossed') || (strcmp(sb, 'exact') && fd < fb)))
        v = 'EARLY WARNING (coherence first)'; early = true;
    elseif strcmp(sb, 'below_floor') || strcmp(sb, 'exact')
        v = 'NO WARNING (bias first)';
    else
        v = 'neither within f2 <= 0.5';
    end
    warnTab(i).fDetect = fd; warnTab(i).fBias = fb; warnTab(i).early = early; warnTab(i).verdict = v;
    fprintf('  phi %2.0f: f_detect %s (%s), f(d_sys = 0.5) %s (%s): %s', V.angleDeg(i), numStr(fd), sd, numStr(fb), sb, v);
    if ~isnan(fd) && ~isnan(fb), fprintf(' (ratio f_bias / f_detect %.2f)', fb / fd); end
    fprintf('\n');
end

%% ---- 4  nonlinear diagnostic -------------------------------------------------------------------------------------------------
% at the d_sys = 1 crossing of the most sensitive angle (smallest crossing f2), where the tolerance conclusions are
% drawn; if no angle crosses, the registered case with d_sys closest to 1 (cfg.V2.nonlinear.rule)
rows = [G.row];
x1 = X(strcmp({X.criterion}, 'd_sys') & [X.threshold] == V.thresholds.dSys(2) & [X.crossed]);
if ~isempty(x1)
    [~, j] = min([x1.f2]);
    rN = evalAt(x1(j).phiDeg * pi / 180, x1(j).f2);
    how = sprintf('the d_sys = %g crossing of the most sensitive angle', V.thresholds.dSys(2));
else
    [~, iN] = min(abs(log10([rows.dSys])));
    rN = rows(iN);
    how = 'no crossing: the registered case with d_sys closest to 1';
end
thLin = log(out0.scn.p([1 3])).' + out0.dthetaW(:) + rN.dthetaSys(:);
fprintf('\nStep 4 (diagnostic): nonlinear check at %s: phi %.0f, f2 %.4g (d_sys %.3f)\n', how, rN.phiDeg, rN.f2, rN.dSys);
cfgN = cfg; cfgN.V1.nonlinear.solver = V.nonlinear.solver;
nl = T1.nonlinear(out0, rN.fW, thLin, cfgN);
nl.caseName = sprintf('phi %.0f f2 %.4g', rN.phiDeg, rN.f2); nl.phiDeg = rN.phiDeg; nl.f2 = rN.f2; nl.dSysLin = rN.dSys;
fprintf(['  %s in %d it (%d model evals, %.0f s): dbeta/beta %+.3f%% (linear %+.3f%%), dR/R %+.4f%% (linear %+.4f%%)\n' ...
    '  NL - (theta_true + dth_W + dth_dir) = [%+.3f %+.3f] sigma, d = %.3f (reference %.1f sigma, not a gate)\n'], ...
    ternary(nl.converged, 'converged', 'NOT converged'), nl.it, nl.evals, nl.seconds, nl.betaPct, nl.linBetaPct, ...
    nl.RPct, nl.linRPct, nl.gapSigma, nl.dGap, V.nonlinear.referenceGapSigma);

%% ---- 5  tables, summary ------------------------------------------------------------------------------------------------------
rowsG = arrayfun(@(r) {r.phiDeg, r.f2, r.dSys, r.betaPct, r.RPct, r.dTot, r.class, r.Dmax, r.threshold, r.detect, ...
    max(r.perBinFrac), min(r.meanCoh), r.sigmaRatio(1), r.sigmaRatio(2), r.nValidChanged}, rows, 'UniformOutput', false);
writeStudyCsv(run, 'grid', {'phi_deg', 'f2', 'd_sys', 'dbeta_beta_pct', 'dR_R_pct', 'd_tot', 'class', 'D_max', ...
    'D_threshold', 'detected', 'max_perbin_fraction', 'min_pair_mean_coherence', 'sigma_ratio_lnbeta', ...
    'sigma_ratio_lnR', 'bins_changing_acceptance'}, vertcat(rowsG{:}));
rowsX = arrayfun(@(t) {t.phiDeg, t.criterion, t.threshold, t.status, t.f2, t.nEval, t.note}, X, 'UniformOutput', false);
writeStudyCsv(run, 'crossings', {'phi_deg', 'criterion', 'threshold', 'status', 'f2', 'evaluations', 'note'}, ...
    vertcat(rowsX{:}));
writeStudyCsv(run, 'nonlinear', {'case', 'converged', 'iterations', 'model_evals', 'd_sys_linear', ...
    'dbeta_beta_pct_NL', 'dbeta_beta_pct_linear', 'dR_R_pct_NL', 'dR_R_pct_linear', 'gap_lnbeta_sigma', ...
    'gap_lnR_sigma', 'd_gap'}, {nl.caseName, nl.converged, nl.it, nl.evals, nl.dSysLin, nl.betaPct, nl.linBetaPct, ...
    nl.RPct, nl.linRPct, nl.gapSigma(1), nl.gapSigma(2), nl.dGap});

Ls = {};
Ls{end + 1} = 'V2: directional forcing (second independent incident direction, frozen single-direction inverse)';
Ls{end + 1} = sprintf('run %s, %s, code commit %s', run.stamp, run.platform, run.commit);
Ls{end + 1} = sprintf('operating point: %s; nominal d_W %.3f, %s', V.operatingPoint, out0.dW, out0.class);
Ls{end + 1} = '';
Ls{end + 1} = sprintf(['Coherence null (MC, %d records): sd of the bin-mean Fisher z %s vs analytic %.4g; effective ' ...
    'independent bins %s of %d.'], nc.nRec, mat2str(nc.sd.', 4), nc.sdAnalytic, mat2str(round(nc.nBinsEff.'), 4), nc.nBins);
Ls{end + 1} = 'Local regime d_sys / f2 at small f2 (spread):';
for i = 1:numel(LIN)
    Ls{end + 1} = sprintf('  phi %2.0f: %s (%.2f%%)', LIN(i).phiDeg, mat2str(LIN(i).dPerF2, 5), 100 * LIN(i).spread); %#ok<SAGROW>
end
Ls{end + 1} = 'Grid (phi deg, f2): d_sys, dbeta/beta, dR/R, d_tot, class (nominal sigma), coherence D_max (detected):';
for k = 1:numel(rows)
    r = rows(k);
    Ls{end + 1} = sprintf('  %2.0f  %.2f  d_sys %8.3f  %+8.3f%%  %+8.4f%%  d_tot %8.3f  %-16s  D_max %7.2f %s', r.phiDeg, ...
        r.f2, r.dSys, r.betaPct, r.RPct, r.dTot, r.class, r.Dmax, ternary(r.detect, '(detected)', '')); %#ok<SAGROW>
end
Ls{end + 1} = 'First crossings in f2:';
for k = 1:numel(X)
    Ls{end + 1} = sprintf('  phi %2.0f  %-13s %-4s: %s', X(k).phiDeg, X(k).criterion, numStr(X(k).threshold), ...
        X(k).note); %#ok<SAGROW>
end
Ls{end + 1} = 'Early warning:';
for i = 1:numel(warnTab)
    Ls{end + 1} = sprintf('  phi %2.0f: f_detect %s, f(d_sys = 0.5) %s: %s', warnTab(i).phiDeg, numStr(warnTab(i).fDetect), ...
        numStr(warnTab(i).fBias), warnTab(i).verdict); %#ok<SAGROW>
end
Ls{end + 1} = sprintf(['Nonlinear check (%s): %s, NL - linear total [%+.3f %+.3f] sigma (d %.3f); dbeta/beta %+.3f%% ' ...
    'vs linear %+.3f%%.'], nl.caseName, ternary(nl.converged, 'converged', 'NOT converged'), nl.gapSigma, nl.dGap, ...
    nl.betaPct, nl.linBetaPct);

strip = @(r) rmfield(r, {'fW', 'b'});                          % keep the .mat small
results = struct('grid', {arrayfun(strip, rows, 'UniformOutput', false)}, 'crossings', X, 'warning', warnTab, ...
    'nonlinear', nl, 'nominal', struct('sigma', out0.sigma, 'dW', out0.dW, 'dthetaW', out0.dthetaW, 'class', out0.class, ...
    'omega', out0.omega, 'validMask', out0.validMask, 'coherence', out0.coherence, 'nEff', out0.nEff), ...
    'baselineCoherence', cohB, 'nullCoherence', nc, 'linearity', LIN, 'nPairs', m, 'nDof', ctx.nDof, 'angleDeg', V.angleDeg, 'f2', f2s, 'tools', T.version);
finishStudy(run, results, cfg, Ls);
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================
function r = cachedCase(cache, T, ctx, out0, phi, f2, cfg)
key = sprintf('%.12e_%.12e', phi, f2);
if isKey(cache, key), r = cache(key); return; end
r = T.evaluate(ctx, out0, phi, f2, cfg);
cache(key) = r; %#ok<NASGU>
end

function v = metricOf(r, field)
switch field
    case 'absBeta', v = abs(r.betaPct);
    case 'detect', v = r.Dmax - r.threshold;
    otherwise, v = r.(field);
end
end

function s = numStr(x)
if isnan(x), s = '-'; else, s = sprintf('%.4g', x); end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end