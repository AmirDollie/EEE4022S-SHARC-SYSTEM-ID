%% runV1ModelMismatch.m
% V1: forward-model mismatch. The data come from a perturbed TRUTH model; the inverse keeps the frozen nominal
% model ([50 10 10], nominal depth, nominal sensor positions). How large is the parameter bias?
%
% WHY. Every recovery study so far generated and inverted with the same EMM (the "inverse crime"). M1 showed the
% estimator behaves as predicted when the model is right; V1 quantifies what happens when it is not, and separates
% estimator error (G3, M1) from forward-model discrepancy.
%
% METHOD (programme rev. 6 and cfg.V1; clarifications logged 2026-10-04 before the run)
%   b_sys = f_W(truth) - f_W(nominal), both Welch expectations from the same function (v1Tools.welchFeatures);
%   dth_sys = F^-1 J' Sigma^-1 b_sys and d_sys with the NOMINAL J, Sigma, F; d_tot from dth_W + dth_sys; class.
%   V1a  truncation: [50 10 10] inverse against cfg.V1.truncationTruth (selected by diagnoseV1Truncation), at p0
%        and the large-kR A1 floes cfg.V1.truncationFloes
%   V1b  depth: true H = H0 (1 + eps), eps in cfg.V1.depthErrorFrac, physically consistent (beta ~ H^-4, gamma, R,
%        r_i ~ H^-1, alpha ~ H), physical sea and band fixed; alpha range reported per case
%   V1c  radial sensor error: every sensor (r > 0) moved by dr in cfg.V1.sensorRadialError (m), angles unchanged
%   Tolerances (V1b, V1c): first crossings of d_sys = 0.5, 1 and |dbeta/beta| = 1% per side, bracketed by the
%   registered points and refined by bisection (r1Boundary); not extended beyond the registered range.
%   Nonlinear check (diagnostic): the case with the largest d_sys, G3 solver from theta_true on its noiseless
%   Welch features, against theta_true + dth_W + dth_sys.
%   p0 cases use the E1 default (Level 2C, reference s26, LSM6DSV16X, setting F).
%
% PRECONDITIONS. diagnoseV1Truncation has run and cfg.V1.truncationTruth / truncationTruthSelected are logged;
% testV1ModelMismatch passes; no uncommitted .m files.
% COST. One EMM field (29 nodes, no derivatives) per case: 4 truncation (slow: high truncation, but cached by
% stage 0 for p0 and the floes), 4 depth, 4 radial, up to about 6 x 6 bisection evaluations, one nonlinear
% inversion (about 15 model evaluations). Roughly 30 to 60 min; each evaluation prints its time.
% OUTPUT. Results/V1/ (log, .mat, summary, V1_truncation.csv, V1_mismatch.csv, V1_tolerances.csv,
% V1_nonlinear.csv). Figures: plotV1ModelMismatch (after the run).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

ALLOW_DIRTY = false;    % true: run with uncommitted .m files (provenance is then NOT the logged commit)

P = setupStudyPaths();
cfg = defineStudyScenarios();
T = v1Tools(); A = a1Tools(); Tm = m1Tools();
V = cfg.V1;
if ~V.truncationTruthSelected
    error('V1:truth', ['cfg.V1.truncationTruth is not yet selected: run diagnoseV1Truncation, paste and log the ' ...
        'result, set cfg.V1.truncationTruthSelected = true, commit.']);
end
run = startStudy('V1', 'Forward-model mismatch', P.results);

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
            error('V1:dirtyTree', '%d uncommitted .m file(s): commit first (or set ALLOW_DIRTY = true).', numel(dirty));
        end
        fprintf('  *** ALLOW_DIRTY: this run is NOT reproducible from commit %s ***\n', run.commit);
    else
        fprintf('  working tree clean (.m files): the code that runs is commit %s\n', run.commit);
    end
end
fprintf('  inverse model [%d %d %d], truth truncation [%d %d %d]\n', V.truncationInverse, V.truncationTruth);

%% ---- 1  nominal p0 -----------------------------------------------------------------------------------------------------
out0 = evaluateScenario(struct('cacheDir', P.cache), cfg);
b0 = T.base(out0);
tc = tic; fW0 = T.welchFeatures(T.truthModel(b0, 'none', [], cfg), out0, P.cache);
zeroCheck = max(abs(fW0 - (out0.f0(:) + out0.bW(:)))) / max(abs(fW0));
fprintf(['\nStep 1: nominal p0: sigma [%.4f%% %.4f%%], d_W %.4f, %s; nominal Welch features (%.0f s) reproduce the ' ...
    'engine''s f0 + b_W to %.1e\n'], 100 * out0.sigma, out0.dW, out0.class, toc(tc), zeroCheck);
if zeroCheck > 1e-10, error('V1:nominal', 'The nominal Welch features do not reproduce the engine (%.2e).', zeroCheck); end

%% ---- 2  V1a truncation ----------------------------------------------------------------------------------------------------
fprintf('\nStep 2 (V1a): [%d %d %d] inverse against truth [%d %d %d]\n', V.truncationInverse, V.truncationTruth);
sc = A.buildScenarios(cfg);
floes = [{'basin_p0'}, V.truncationFloes(:).'];
TR = struct('floe', floes, 'row', []);
outs = cell(1, numel(floes)); bases = cell(1, numel(floes)); fWn = cell(1, numel(floes));
for i = 1:numel(floes)
    if i == 1
        outs{i} = out0; bases{i} = b0; fWn{i} = fW0;
    else
        s = sc(strcmp({sc.name}, floes{i}));
        outs{i} = evaluateScenario(A.buildScn(s, cfg, struct('cacheDir', P.cache)), cfg);
        bases{i} = T.base(outs{i});
        fWn{i} = T.welchFeatures(T.truthModel(bases{i}, 'none', [], cfg), outs{i}, P.cache);
    end
    tc = tic;
    r = T.evaluate(bases{i}, outs{i}, fWn{i}, 'truncation', V.truncationTruth, cfg, P.cache);
    TR(i).row = r;
    fprintf(['  %-24s d_sys %.3f  dbeta/beta %+.3f%%  dR/R %+.4f%%  | d_W %.3f, d_tot %.3f: %s (production class ' ...
        '%s)  (%.0f s)\n'], floes{i}, r.dSys, r.betaPct, r.RPct, outs{i}.dW, r.dTot, r.class, outs{i}.class, toc(tc));
end

%% ---- 3-4  V1b depth, V1c radial ---------------------------------------------------------------------------------------------
axesV = struct('kind', {'depth', 'radial'}, 'values', {V.depthErrorFrac, V.sensorRadialError}, ...
    'unit', {'%', 'mm'}, 'scale', {100, 1000});
MM = struct('kind', {}, 'row', {});
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
evalAt = @(kind, x) cachedCase(cache, T, b0, out0, fW0, kind, x, cfg, P.cache);
for a = 1:numel(axesV)
    ax = axesV(a);
    fprintf('\nStep %d (%s): %s = %s %s\n', 2 + a, ternary(a == 1, 'V1b', 'V1c'), ax.kind, ...
        mat2str(ax.scale * ax.values), ax.unit);
    for x = ax.values
        tc = tic; r = evalAt(ax.kind, x);
        MM(end + 1) = struct('kind', ax.kind, 'row', r); %#ok<SAGROW>
        al = '';
        if strcmp(ax.kind, 'depth')
            al = sprintf('  alpha retained [%.3f %.3f], support [%.3f %.3f]%s', r.alpha.minRetained, ...
                r.alpha.maxRetained, r.alpha.minSupport, r.alpha.maxSupport, ternary(r.alpha.extrapolated, ...
                ' EXTRAPOLATED', ''));
        end
        fprintf('  %+6.2f %-2s  d_sys %.3f  dbeta/beta %+.3f%%  dR/R %+.4f%%  d_tot %.3f  %s  (%.0f s)%s\n', ...
            ax.scale * x, ax.unit, r.dSys, r.betaPct, r.RPct, r.dTot, r.class, toc(tc), al);
    end
end

%% ---- 5  tolerances --------------------------------------------------------------------------------------------------------
fprintf('\nStep 5: tolerances (first crossing per side; bisection in log10 |error| to %.3g decade)\n', ...
    V.toleranceTolDecades);
crit = struct('name', {'d_sys', 'd_sys', '|dbeta/beta|'}, 'field', {'dSys', 'dSys', 'absBeta'}, ...
    'threshold', {V.toleranceThresholds.dSys(1), V.toleranceThresholds.dSys(2), V.toleranceThresholds.betaPct});
TOL = struct('kind', {}, 'side', {}, 'criterion', {}, 'threshold', {}, 'crossed', {}, 'x', {}, 'nEval', {}, 'note', {});
for a = 1:numel(axesV)
    ax = axesV(a);
    for side = [-1 1]
        xs = sort(ax.values(sign(ax.values) == side), 'ascend');
        if side < 0, xs = fliplr(xs); end                                 % outward from 0
        xs = [0, xs]; %#ok<AGROW>
        for c = 1:numel(crit)
            met = @(x) metricOf(evalAt(ax.kind, x), crit(c).field);
            vals = arrayfun(met, xs);
            k = find(vals(1:end - 1) < crit(c).threshold & vals(2:end) >= crit(c).threshold, 1);
            t = struct('kind', ax.kind, 'side', side, 'criterion', crit(c).name, 'threshold', crit(c).threshold, ...
                'crossed', ~isempty(k), 'x', NaN, 'nEval', 0, 'note', '');
            if isempty(k)
                t.note = sprintf('not crossed within %s%g %s (max %.3f)', ternary(side < 0, '-', '+'), ...
                    ax.scale * max(abs(xs)), ax.unit, max(vals));
            else
                % bisection in u = log10 |x| (logged pre-run correction): crossings can lie far below the first
                % registered point; from 0, the lower bracket is toleranceFloorFraction of that point
                xa = xs(k); xb0 = xs(k + 1); nExtra = 0;
                if xa == 0
                    xa = xb0 * V.toleranceFloorFraction; nExtra = 1;
                    if met(xa) >= crit(c).threshold
                        t.note = sprintf('already exceeded at %+.3g %s: tolerance below that', ax.scale * xa, ax.unit);
                        t.x = NaN; t.nEval = nExtra;
                    end
                end
                if isempty(t.note)
                    fu = @(u) met(side * 10^u);
                    [ub, info] = r1Boundary(fu, crit(c).threshold, log10(abs([xa xb0])), V.toleranceTolDecades);
                    t.x = side * 10^ub; t.nEval = info.nEval + nExtra;
                    t.note = sprintf('%+.4g %s', ax.scale * t.x, ax.unit);
                end
            end
            TOL(end + 1) = t; %#ok<SAGROW>
            fprintf('  %-7s %s  %-13s = %-4g: %s\n', ax.kind, ternary(side < 0, '-', '+'), crit(c).name, ...
                crit(c).threshold, t.note);
        end
    end
end

%% ---- 6  nonlinear check (diagnostic) -------------------------------------------------------------------------------------
allRows = [TR.row, MM.row];                                      % same fields (v1Tools.evaluate)
allOut = [outs, repmat({out0}, 1, numel(MM))];
allName = [cellfun(@(f) ['truncation ' f], floes, 'UniformOutput', false), ...
    arrayfun(@(m) sprintf('%s %+g', m.kind, m.row.value), MM, 'UniformOutput', false)];
[~, iMax] = max([allRows.dSys]);
rN = allRows(iMax); oN = allOut{iMax};
thLin = log(oN.scn.p([1 3])).' + oN.dthetaW(:) + rN.dthetaSys(:);
fprintf('\nStep 6 (diagnostic): nonlinear check of the largest d_sys: %s (d_sys %.3f, d_tot %.3f)\n', allName{iMax}, ...
    rN.dSys, rN.dTot);
nl = T.nonlinear(oN, rN.fW, thLin, cfg);
fprintf(['  %s in %d it (%d model evals, %.0f s): dbeta/beta %+.3f%% (linear %+.3f%%), dR/R %+.4f%% (linear %+.4f%%)\n' ...
    '  NL - (theta_true + dth_W + dth_sys) = [%+.3f %+.3f] sigma, d = %.3f (reference %.1f sigma, not a gate)\n'], ...
    ternary(nl.converged, 'converged', 'NOT converged'), nl.it, nl.evals, nl.seconds, nl.betaPct, nl.linBetaPct, ...
    nl.RPct, nl.linRPct, nl.gapSigma, nl.dGap, V.nonlinear.referenceGapSigma);
nl.caseName = allName{iMax};

%% ---- 7  tables, summary --------------------------------------------------------------------------------------------------
rowsT = arrayfun(@(t) {t.floe, sprintf('[%d %d %d]', V.truncationInverse), sprintf('[%d %d %d]', V.truncationTruth), ...
    t.row.dSys, t.row.betaPct, t.row.RPct, outs{strcmp(floes, t.floe)}.dW, t.row.dTot, t.row.class, ...
    outs{strcmp(floes, t.floe)}.class}, TR, 'UniformOutput', false);
writeStudyCsv(run, 'truncation', {'floe', 'inverse', 'truth', 'd_sys', 'dbeta_beta_pct', 'dR_R_pct', 'd_W', 'd_tot', ...
    'class_with_mismatch', 'class_matched'}, vertcat(rowsT{:}));
rowsM = arrayfun(@(m) {m.kind, m.row.value, m.row.dSys, m.row.betaPct, m.row.RPct, m.row.dTot, m.row.class, ...
    m.row.alpha.minRetained, m.row.alpha.maxRetained, m.row.alpha.minSupport, m.row.alpha.maxSupport, ...
    m.row.alpha.extrapolated}, MM, 'UniformOutput', false);
writeStudyCsv(run, 'mismatch', {'kind', 'value', 'd_sys', 'dbeta_beta_pct', 'dR_R_pct', 'd_tot', 'class', ...
    'alpha_min_retained', 'alpha_max_retained', 'alpha_min_support', 'alpha_max_support', 'alpha_extrapolated'}, ...
    vertcat(rowsM{:}));
rowsL = arrayfun(@(t) {t.kind, t.side, t.criterion, t.threshold, t.crossed, t.x, t.nEval, t.note}, TOL, ...
    'UniformOutput', false);
writeStudyCsv(run, 'tolerances', {'kind', 'side', 'criterion', 'threshold', 'crossed', 'value', 'evaluations', 'note'}, ...
    vertcat(rowsL{:}));
writeStudyCsv(run, 'nonlinear', {'case', 'converged', 'iterations', 'model_evals', 'dbeta_beta_pct_NL', ...
    'dbeta_beta_pct_linear', 'dR_R_pct_NL', 'dR_R_pct_linear', 'gap_lnbeta_sigma', 'gap_lnR_sigma', 'd_gap'}, ...
    {nl.caseName, nl.converged, nl.it, nl.evals, nl.betaPct, nl.linBetaPct, nl.RPct, nl.linRPct, nl.gapSigma(1), ...
    nl.gapSigma(2), nl.dGap});

Ls = {};
Ls{end + 1} = 'V1: forward-model mismatch (data from a perturbed truth model, frozen nominal inverse)';
Ls{end + 1} = sprintf('run %s, %s, code commit %s', run.stamp, run.platform, run.commit);
Ls{end + 1} = '';
Ls{end + 1} = sprintf('V1a truncation: inverse [%d %d %d], truth [%d %d %d]', V.truncationInverse, V.truncationTruth);
for i = 1:numel(TR)
    r = TR(i).row;
    Ls{end + 1} = sprintf('  %-24s d_sys %.3f, dbeta/beta %+.3f%%, dR/R %+.4f%%; d_tot %.3f: %s (matched: %s)', floes{i}, ...
        r.dSys, r.betaPct, r.RPct, r.dTot, r.class, outs{i}.class); %#ok<SAGROW>
end
Ls{end + 1} = 'V1b depth and V1c radial position (p0, Level 2C):';
for k = 1:numel(MM)
    r = MM(k).row; ax = axesV(strcmp({axesV.kind}, MM(k).kind));
    Ls{end + 1} = sprintf('  %-6s %+6.2f %-2s d_sys %.3f, dbeta/beta %+.3f%%, dR/R %+.4f%%, d_tot %.3f: %s%s', MM(k).kind, ...
        ax.scale * r.value, ax.unit, r.dSys, r.betaPct, r.RPct, r.dTot, r.class, ...
        ternary(strcmp(MM(k).kind, 'depth') && r.alpha.extrapolated, ' (alpha extrapolated)', '')); %#ok<SAGROW>
end
Ls{end + 1} = 'Tolerances (first crossing per side):';
for k = 1:numel(TOL)
    t = TOL(k);
    Ls{end + 1} = sprintf('  %-6s %s %-13s %-4g: %s', t.kind, ternary(t.side < 0, '-', '+'), t.criterion, t.threshold, ...
        t.note); %#ok<SAGROW>
end
Ls{end + 1} = sprintf(['Nonlinear check (%s): %s, NL - linear total [%+.3f %+.3f] sigma (d %.3f); dbeta/beta %+.3f%% ' ...
    'vs linear %+.3f%%.'], nl.caseName, ternary(nl.converged, 'converged', 'NOT converged'), nl.gapSigma, nl.dGap, ...
    nl.betaPct, nl.linBetaPct);

strip = @(r) rmfield(r, {'fW', 'bSys'});                        % keep the .mat small
results = struct('truncation', struct('floes', {floes}, 'rows', {arrayfun(@(t) strip(t.row), TR, 'UniformOutput', false)}, ...
    'truth', V.truncationTruth, 'inverse', V.truncationInverse), ...
    'mismatch', {arrayfun(@(m) strip(m.row), MM, 'UniformOutput', false)}, 'tolerances', TOL, 'nonlinear', nl, ...
    'nominal', struct('sigma', out0.sigma, 'dW', out0.dW, 'dthetaW', out0.dthetaW, 'class', out0.class, ...
    'omega', out0.omega, 'F', out0.F), 'axes', axesV, 'tools', T.version);
finishStudy(run, results, cfg, Ls);
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================
function r = cachedCase(cache, T, b0, out0, fW0, kind, x, cfg, cacheDir)
key = sprintf('%s_%.12e', kind, x);
if isKey(cache, key), r = cache(key); return; end
if x == 0
    r = T.evaluate(b0, out0, fW0, 'none', [], cfg, cacheDir);
else
    r = T.evaluate(b0, out0, fW0, kind, x, cfg, cacheDir);
end
cache(key) = r; %#ok<NASGU>
end

function v = metricOf(r, field)
if strcmp(field, 'absBeta'), v = abs(r.betaPct); else, v = r.(field); end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end