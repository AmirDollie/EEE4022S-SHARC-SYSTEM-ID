%% runD1SensorLayout.m
% STUDY D1 (MUST): how many sensors, where, and which one is the reference?
%
% QUESTION. What does each additional sensor buy in parameter identifiability, measured with the
% Fisher information of the actual transmissibility inverse, and how does the earlier Gramian
% placement compare?
%
% METHOD (pre-registered, cfg.D1). Nominal floe p0, production setting F, the 33-point Level 2C
% candidate grid, noise cases cfg.D1.noiseCases (the first, LSM6DSV16X, is the headline).
%   1. ONE emmField for all 33 points (sensors fixed at their physical positions at the true R).
%   2. Per reference r (33 of them): transmissibilities to the other 32 points, Jacobian, Welch bias
%      and the full 64 x 64 shared-reference covariance per accepted bin (d1LayoutTools.precompute).
%   3. GATE: the Level 2C layout [s2 s6 s14 s26] with reference s26, extracted from step 2, must
%      reproduce evaluateScenario (the frozen E1 chain) before any search runs.
%   4. Exhaustive search over all layouts for n_s in cfg.D1.sensorCounts, EVERY sensor of each layout
%      tried as reference. Criteria: A = tr(F^-1) (primary, minimise), D = log det F and
%      E = lambda_min(F) (secondary, maximise); each criterion picks its own reference.
%   5. Admissibility: radius <= 0.9 R, buoy wholly on the floe, centre spacing >= cfg.D1.minSpacing
%      (0.10 m, primary); 0.05 and 0.15 m and the unconstrained optimum reported alongside.
%   6. Comparisons: Level 2C (reference s26, and reference re-optimised) placed in the n_s = 4 A
%      ranking; the Gramian lambda_min optimum per n_s (recomputed exactly as level2C_finalCheck),
%      its Fisher rank and efficiency eta = A_opt / A_Gramian.
%   7. Sensor count: gain G = 1 - sigma_lnbeta(n+1) / sigma_lnbeta(n) at the A-optima ("worth it" if
%      >= cfg.D1.worthItFraction); recommended count = smallest n_s whose optimum is Identifiable at
%      the headline noise. Distribution of sigma_lnbeta (best reference per layout) per n_s.
%   Singular F is a result (A = Inf, Not identifiable), never regularised.
%
% RUN ORDER. testD1SensorLayout must pass first. Then QUICK = true (n_s = 2, headline noise only, a few
% seconds after the field) to check the ranking and reference bookkeeping; then QUICK = false.
%
% COST. One 33-point EMM field (5 evaluations, cached in Results/twinCache), then about 181 000
% layout-reference evaluations per noise case, small matrices only.
%
% OUTPUT. Results/D1/ (log, .mat with everything needed to redraw, summary, CSV tables; figures via
% plotD1SensorLayout once it exists).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

QUICK = false; % true: n_s = 2 and the headline noise case only (bookkeeping check)
RUN_GRAMIAN = true;             % the Gramian comparison (cheap; needs the SSI and Forward Model code)

P = setupStudyPaths();
cfg = defineStudyScenarios();
T = d1LayoutTools();
run = startStudy('D1', 'Sensor count, layout and reference', P.results);

try
p0 = cfg.p0.vec;
nsList = cfg.D1.sensorCounts(:).';
noiseCases = cfg.D1.noiseCases;
if QUICK, nsList = nsList(1); noiseCases = noiseCases(1); fprintf('QUICK run: n_s = %d, %s only\n', nsList, noiseCases{1}); end
headline = noiseCases{1};
nN = numel(noiseCases); nK = numel(nsList);
spacings = [cfg.D1.minSpacing, cfg.D1.minSpacingSensitivity(:).', 0];   % primary, sensitivity, unconstrained
spacingNames = [{'primary'}, arrayfun(@(s) sprintf('spacing %.2f m', s), cfg.D1.minSpacingSensitivity, ...
    'UniformOutput', false), {'unconstrained'}];

%% ---- 1  geometry and admissibility --------------------------------------------------------------------
geom = T.buildGeometry(cfg, p0);
dd = geom.dist + diag(Inf(size(geom.dist, 1), 1));
fprintf('Step 1: %d candidate points, R = %.3f m, closest pair %.3f m; Level 2C = %s, reference %s\n', ...
    size(geom.xy, 1), geom.Rphys, min(dd(:)), mat2str(geom.level2C), geom.names{geom.level2CRef});
lay = cell(1, nK); okAll = cell(1, nK);
for k = 1:nK
    okAll{k} = false(nchoosek(size(geom.xy, 1), nsList(k)), numel(spacings));
    for s = 1:numel(spacings)
        [lay{k}, okAll{k}(:, s)] = T.enumerateLayouts(geom, nsList(k), cfg, spacings(s));
    end
    fprintf('  n_s = %d: %d layouts; admissible: %s\n', nsList(k), size(lay{k}, 1), ...
        strjoin(arrayfun(@(s) sprintf('%s %d', spacingNames{s}, nnz(okAll{k}(:, s))), 1:numel(spacings), ...
        'UniformOutput', false), ', '));
end

%% ---- 2  the 33-point field ------------------------------------------------------------------------------
fprintf('\nStep 2: EMM field at all %d points (5 evaluations; cached after the first run)\n', size(geom.xy, 1));
full = T.buildFullField(cfg, p0, geom, P.cache, true);
fprintf('  %d retained bins, N_eff = %.3f\n', numel(full.omega), full.nEff);

%% ---- 3-4  per noise case: precompute, gate, search -----------------------------------------------------
search = cell(nN, nK); gate = struct('noise', noiseCases, 'relSigma', [], 'relDW', [], 'relF', []);
for c = 1:nN
    Snn = cfg.noise.(noiseCases{c});
    fprintf('\nStep 3 (%s): precompute 33 references\n', noiseCases{c});
    pre = T.precompute(full, Snn, cfg);
    fprintf('  accepted bins per reference: %d to %d of %d; full covariance blocks positive definite for %d of 33\n', ...
        min([pre.nValid]), max([pre.nValid]), numel(full.omega), nnz([pre.fullPD]));
    % GATE: Level 2C with s26 against the frozen E1 chain
    e1 = evaluateScenario(struct('Snn', Snn, 'cacheDir', P.cache), cfg);
    qd = T.evaluateLayout(pre, geom.level2C, geom.level2CRef, cfg, 'direct');
    qf = T.evaluateLayout(pre, geom.level2C, geom.level2CRef, cfg, 'fast');
    gate(c).relSigma = max(abs([qd.sigma; qf.sigma] ./ [e1.sigma; e1.sigma] - 1));
    gate(c).relDW = max(abs([qd.dW qf.dW] / e1.dW - 1));
    gate(c).relF = max(norm(qd.F - e1.F), norm(qf.F - e1.F)) / norm(e1.F);
    fprintf('  GATE Level 2C + s26 vs evaluateScenario: sigma %.2e, d_W %.2e, F %.2e (relative)\n', ...
        gate(c).relSigma, gate(c).relDW, gate(c).relF);
    if gate(c).relSigma > 1e-6 || gate(c).relDW > 1e-6 || gate(c).relF > 1e-6
        error('runD1:gate', 'Level 2C is not reproduced (%s): stop and check the indexing.', noiseCases{c});
    end
    fprintf('Step 4 (%s): exhaustive search\n', noiseCases{c});
    for k = 1:nK
        search{c, k} = T.searchLayouts(pre, lay{k}, cfg, noiseCases{c});
        search{c, k}.admissible = okAll{k};
    end
    clear pre
end

%% ---- 5  optima, distribution, Level 2C ---------------------------------------------------------------------
nameOf = @(L) strjoin(geom.names(double(L)), ' ');
best = struct();
for c = 1:nN
    for k = 1:nK
        S = search{c, k};
        for s = 1:numel(spacings)
            ok = S.admissible(:, s);
            best(c, k, s).A = pickBest(S, ok, 'A');
        end
        ok = S.admissible(:, 1);
        best(c, k, 1).D = pickBest(S, ok, 'D');
        best(c, k, 1).E = pickBest(S, ok, 'E');
        % distribution over admissible layouts (best reference per layout)
        a = S.A(ok); sb = S.sB(ok); fin = isfinite(a);
        Aopt = min(a);
        dist(c, k) = struct('nAdmissible', nnz(ok), 'nSingular', nnz(~fin), ...
            'pctSigmaB', prctileSafe(sb(fin), [0 5 25 50 75 95 100]), ...
            'within5', nnz(a <= 1.05 * Aopt), 'within10', nnz(a <= 1.10 * Aopt), ...
            'refMatters', median(max(S.Aall(ok, :), [], 2) ./ min(S.Aall(ok, :), [], 2), 'omitnan')); %#ok<SAGROW>
    end
end
% sensor count: gain and recommendation (primary admissibility, A-optima)
countRes = struct('noise', noiseCases);
for c = 1:nN
    sb = arrayfun(@(k) best(c, k, 1).A.sB, 1:nK);
    cls = arrayfun(@(k) best(c, k, 1).A.classCode, 1:nK);
    G = [NaN, 1 - sb(2:end) ./ sb(1:end - 1)];
    rec = find(cls == 1, 1, 'first');
    countRes(c).ns = nsList; countRes(c).sigmaB = sb; countRes(c).gain = G;
    countRes(c).worthIt = G >= cfg.D1.worthItFraction; countRes(c).classCode = cls;
    countRes(c).recommended = NaN; if ~isempty(rec), countRes(c).recommended = nsList(rec); end
end
% Level 2C in the n_s = 4 ranking
L2 = struct('noise', noiseCases);
k4 = find(nsList == numel(geom.level2C));
for c = 1:nN
    if isempty(k4), break; end
    S = search{c, k4};
    row = find(all(double(S.layouts) == geom.level2C, 2));
    ok = S.admissible(:, 1);
    pos26 = find(geom.level2C == geom.level2CRef);
    A26 = S.Aall(row, pos26);
    L2(c).row = row;
    L2(c).A_s26 = A26; L2(c).A_bestRef = S.A(row); L2(c).bestRef = double(S.refA(row));
    L2(c).rank_s26 = 1 + nnz(S.A(ok) < A26);          % where s26 would sit among best-reference layouts
    L2(c).rank_bestRef = 1 + nnz(S.A(ok) < S.A(row));
    L2(c).eff_s26 = best(c, k4, 1).A.A / A26;
    L2(c).eff_bestRef = best(c, k4, 1).A.A / S.A(row);
    L2(c).sB_bestRef = S.sB(row); L2(c).nAdmissible = nnz(ok);
end

%% ---- 6  Gramian comparison --------------------------------------------------------------------------------
gram = [];
if RUN_GRAMIAN
    try
        fprintf('\nStep 6: Gramian lambda_min optima (as level2C_finalCheck, N = 4001)\n');
        gram = T.gramianLayouts(cfg, nsList);
        for k = 1:nK
            for c = 1:nN
                S = search{c, k}; ok = S.admissible(:, 1);
                tied = gram(k).tied; r = zeros(size(tied, 1), 1); eta = r; sb = r;
                for t = 1:size(tied, 1)
                    row = find(all(double(S.layouts) == tied(t, :), 2));
                    r(t) = 1 + nnz(S.A(ok) < S.A(row)); eta(t) = best(c, k, 1).A.A / S.A(row); sb(t) = S.sB(row);
                end
                gram(k).fisher(c) = struct('noise', noiseCases{c}, 'rank', r, 'eta', eta, 'sigmaB', sb, ...
                    'nAdmissible', nnz(ok));
            end
            fprintf('  n_s = %d: Gramian optimum %s (%d tied); Fisher rank %s of %d, eta %s (%s)\n', nsList(k), ...
                nameOf(gram(k).layout), size(gram(k).tied, 1), mat2str(gram(k).fisher(1).rank.'), ...
                gram(k).fisher(1).nAdmissible, mat2str(round(gram(k).fisher(1).eta.', 3)), headline);
        end
    catch errG
        warning('runD1:gramian', 'Gramian comparison failed (%s); the rest of D1 is unaffected.', errG.message);
        gram = [];
    end
end

%% ---- 7  tables ----------------------------------------------------------------------------------------------
rowsT = {};
crit = {'A', 'D', 'E'};
for c = 1:nN
    for k = 1:nK
        for s = 1:numel(spacings)
            b = best(c, k, s).A;
            rowsT(end + 1, :) = {noiseCases{c}, nsList(k), ['A, ' spacingNames{s}], nameOf(b.layout), ...
                geom.names{b.ref}, 100 * b.sB, 100 * b.sR, b.dW, b.kappa, b.A, cfg.class.names{b.classCode}}; %#ok<AGROW>
        end
        for j = 2:3
            b = best(c, k, 1).(crit{j});
            rowsT(end + 1, :) = {noiseCases{c}, nsList(k), [crit{j} ', primary'], nameOf(b.layout), ...
                geom.names{b.ref}, 100 * b.sB, 100 * b.sR, b.dW, b.kappa, b.A, cfg.class.names{b.classCode}}; %#ok<AGROW>
        end
    end
end
writeStudyCsv(run, 'optimum', {'noise_case', 'n_s', 'criterion', 'layout', 'reference', 'sigma_lnbeta_pct', ...
    'sigma_lnR_pct', 'd_W', 'kappa_F', 'A_trace_Finv', 'class'}, rowsT);
rowsT = {};
for c = 1:nN
    for k = 1:nK
        rowsT(end + 1, :) = {noiseCases{c}, nsList(k), 100 * countRes(c).sigmaB(k), 100 * countRes(c).gain(k), ...
            countRes(c).worthIt(k), cfg.class.names{countRes(c).classCode(k)}, countRes(c).recommended}; %#ok<AGROW>
    end
end
writeStudyCsv(run, 'sensor_count', {'noise_case', 'n_s', 'sigma_lnbeta_opt_pct', 'gain_vs_previous_pct', ...
    'worth_it', 'class', 'recommended_n_s'}, rowsT);
rowsT = {};
for c = 1:nN
    for k = 1:nK
        d = dist(c, k);
        rowsT(end + 1, :) = [{noiseCases{c}, nsList(k), d.nAdmissible, d.nSingular}, num2cell(100 * d.pctSigmaB), ...
            {d.within5, d.within10, d.refMatters}]; %#ok<AGROW>
    end
end
writeStudyCsv(run, 'distribution', {'noise_case', 'n_s', 'n_admissible', 'n_singular', 'sigma_lnbeta_min_pct', ...
    'p5_pct', 'p25_pct', 'median_pct', 'p75_pct', 'p95_pct', 'max_pct', 'n_within_5pct_of_Aopt', ...
    'n_within_10pct_of_Aopt', 'median_worst_over_best_reference_A'}, rowsT);
rowsT = {};
for c = 1:nN
    if ~isempty(k4) && isfield(L2, 'row')
        rowsT(end + 1, :) = {noiseCases{c}, 'Level 2C, reference s26', nameOf(geom.level2C), ...
            geom.names{geom.level2CRef}, L2(c).A_s26, L2(c).rank_s26, L2(c).nAdmissible, L2(c).eff_s26}; %#ok<AGROW>
        rowsT(end + 1, :) = {noiseCases{c}, 'Level 2C, reference optimised', nameOf(geom.level2C), ...
            geom.names{L2(c).bestRef}, L2(c).A_bestRef, L2(c).rank_bestRef, L2(c).nAdmissible, L2(c).eff_bestRef}; %#ok<AGROW>
    end
    for k = 1:numel(gram)
        f = gram(k).fisher(c);
        for t = 1:numel(f.rank)
            S = search{c, k}; row = find(all(double(S.layouts) == gram(k).tied(t, :), 2));
            rowsT(end + 1, :) = {noiseCases{c}, sprintf('Gramian lambda_min, n_s = %d', nsList(k)), ...
                nameOf(gram(k).tied(t, :)), geom.names{S.refA(row)}, S.A(row), f.rank(t), f.nAdmissible, f.eta(t)}; %#ok<AGROW>
        end
    end
end
writeStudyCsv(run, 'comparison', {'noise_case', 'design', 'layout', 'reference', 'A_trace_Finv', 'fisher_rank', ...
    'n_admissible', 'efficiency_Aopt_over_A'}, rowsT);

%% ---- 8  summary -----------------------------------------------------------------------------------------------
L = {};
L{end + 1} = sprintf('Setting: p0, setting F (N = %d, L = %d), %d-point grid, R = %.3f m; admissible: r <= %.1f R, buoy %.2f m.', ...
    cfg.acq.Nref, cfg.welch.L, size(geom.xy, 1), geom.Rphys, cfg.D1.maxRadiusFrac, cfg.D1.buoyDiameter);
L{end + 1} = sprintf('Gate (Level 2C + s26 vs evaluateScenario), worst over noise cases: sigma %.1e, d_W %.1e, F %.1e.', ...
    max([gate.relSigma]), max([gate.relDW]), max([gate.relF]));
L{end + 1} = '';
L{end + 1} = sprintf('HEADLINE (%s): A-optimal designs, primary admissibility', headline);
for k = 1:nK
    b = best(1, k, 1).A;
    L{end + 1} = sprintf('  n_s = %d: %-16s ref %-4s sigma [%.3f%% %.4f%%], d_W %.3f, kappa %.0f, %s; gain %s, worth it: %s', ...
        nsList(k), nameOf(b.layout), geom.names{b.ref}, 100 * b.sB, 100 * b.sR, b.dW, b.kappa, ...
        cfg.class.names{b.classCode}, pctStr(countRes(1).gain(k)), yn(countRes(1).worthIt(k), k)); %#ok<AGROW>
end
L{end + 1} = sprintf('  Recommended n_s (smallest Identifiable optimum): %s', num2str(countRes(1).recommended));
L{end + 1} = '';
L{end + 1} = 'Secondary criteria (each with its own best reference):';
for k = 1:nK
    L{end + 1} = sprintf('  n_s = %d: D-opt %s (ref %s, sigma_lnbeta %.3f%%) | E-opt %s (ref %s, sigma_lnbeta %.3f%%)', ...
        nsList(k), nameOf(best(1, k, 1).D.layout), geom.names{best(1, k, 1).D.ref}, 100 * best(1, k, 1).D.sB, ...
        nameOf(best(1, k, 1).E.layout), geom.names{best(1, k, 1).E.ref}, 100 * best(1, k, 1).E.sB); %#ok<AGROW>
end
L{end + 1} = '';
L{end + 1} = 'Admissibility sensitivity (A-optimum, headline):';
for k = 1:nK
    s = arrayfun(@(j) sprintf('%s: %s (%.3f%%)', spacingNames{j}, nameOf(best(1, k, j).A.layout), ...
        100 * best(1, k, j).A.sB), 1:numel(spacings), 'UniformOutput', false);
    L{end + 1} = sprintf('  n_s = %d | %s', nsList(k), strjoin(s, ' | ')); %#ok<AGROW>
end
L{end + 1} = '';
L{end + 1} = 'How forgiving is placement (headline, best reference per layout)?';
for k = 1:nK
    d = dist(1, k);
    L{end + 1} = sprintf('  n_s = %d: %d admissible (%d singular); sigma_lnbeta min %.3f%%, median %.3f%%, p95 %.3f%%; %d within 5%% and %d within 10%% of the optimal A; reference choice changes A by x%.2f (median)', ...
        nsList(k), d.nAdmissible, d.nSingular, 100 * d.pctSigmaB(1), 100 * d.pctSigmaB(4), 100 * d.pctSigmaB(6), ...
        d.within5, d.within10, d.refMatters); %#ok<AGROW>
end
if ~isempty(k4) && isfield(L2, 'row')
    L{end + 1} = '';
    L{end + 1} = sprintf('Level 2C (%s) at %s: with s26 rank %d of %d (efficiency %.3f); best reference %s, rank %d (efficiency %.3f)', ...
        nameOf(geom.level2C), headline, L2(1).rank_s26, L2(1).nAdmissible, L2(1).eff_s26, geom.names{L2(1).bestRef}, ...
        L2(1).rank_bestRef, L2(1).eff_bestRef);
end
if ~isempty(gram)
    L{end + 1} = '';
    L{end + 1} = sprintf('Gramian lambda_min optimum (N = %d) in the Fisher A ranking (%s):', gram(1).N, headline);
    for k = 1:nK
        L{end + 1} = sprintf('  n_s = %d: %s (%d tied): rank %s of %d, eta %s', nsList(k), nameOf(gram(k).layout), ...
            size(gram(k).tied, 1), mat2str(gram(k).fisher(1).rank.'), gram(k).fisher(1).nAdmissible, ...
            mat2str(round(gram(k).fisher(1).eta.', 3))); %#ok<AGROW>
    end
end
if nN > 1
    L{end + 1} = '';
    L{end + 1} = 'Other noise cases (A-optimum, primary):';
    for c = 2:nN
        s = arrayfun(@(k) sprintf('n_s %d: %s ref %s %.3f%% %s', nsList(k), nameOf(best(c, k, 1).A.layout), ...
            geom.names{best(c, k, 1).A.ref}, 100 * best(c, k, 1).A.sB, cfg.class.names{best(c, k, 1).A.classCode}), ...
            1:nK, 'UniformOutput', false);
        L{end + 1} = sprintf('  %-11s %s | recommended %s', noiseCases{c}, strjoin(s, ' | '), num2str(countRes(c).recommended)); %#ok<AGROW>
    end
end

%% ---- results (everything the figures need, no recomputation) -------------------------------------------------
results = struct('quick', QUICK, 'ns', nsList, 'noise', {noiseCases}, 'headline', headline, ...
    'spacings', spacings, 'spacingNames', {spacingNames}, 'omega', full.omega, 'nEff', full.nEff, 'gate', gate);
results.geometry = struct('xy', geom.xy, 'Rphys', geom.Rphys, 'frac', geom.frac, 'names', {geom.names}, ...
    'level2C', geom.level2C, 'level2CRef', geom.level2CRef, 'maxRadiusFrac', cfg.D1.maxRadiusFrac, ...
    'buoyDiameter', cfg.D1.buoyDiameter);
results.search = search;            % per noise case x n_s: every layout, A per reference, best-reference metrics
results.best = best;                % (noise, n_s, spacing).A / (noise, n_s, 1).D / .E: layout, ref, metrics
results.count = countRes; results.distribution = dist; results.level2C = L2; results.gramian = gram;
if exist('plotD1SensorLayout', 'file')
    try
        plotD1SensorLayout(results, cfg, run);
    catch errPlot
        warning('runD1:plot', 'Figures failed (%s); results are still saved. Redraw with plotD1SensorLayout.', errPlot.message);
    end
end
finishStudy(run, results, cfg, L);
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================
function b = pickBest(S, ok, crit)
% best admissible layout under criterion crit, each criterion with ITS OWN best reference; the
% metrics reported are those of that layout with that reference
idx = find(ok);
switch crit
    case 'A', [~, i] = min(S.A(idx)); ref = S.refA(idx(i));
    case 'D', [~, i] = max(S.D(idx)); ref = S.refD(idx(i));
    case 'E', [~, i] = max(S.E(idx)); ref = S.refE(idx(i));
end
r = idx(i);
lay = double(S.layouts(r, :));
j = find(lay == double(ref), 1);
b = struct('row', r, 'layout', lay, 'ref', double(ref), 'A', S.Aall(r, j), 'sB', S.sBall(r, j), ...
    'sR', S.sRall(r, j), 'dW', S.dWall(r, j), 'kappa', S.kappaAll(r, j), 'classCode', double(S.classAll(r, j)), ...
    'D', S.D(r), 'E', S.E(r));
end

function p = prctileSafe(x, q)
% percentiles without the Statistics Toolbox (nearest-rank on the sorted values)
x = sort(x(:)); n = numel(x);
if n == 0, p = NaN(size(q)); return; end
p = x(max(1, min(n, round(q / 100 * (n - 1)) + 1))).';
end

function s = pctStr(g)
if isnan(g), s = '-'; else, s = sprintf('%.1f%%', 100 * g); end
end

function s = yn(tf, k)
if k == 1, s = '-'; elseif tf, s = 'yes'; else, s = 'no'; end
end