%% runDirectionRobustLayoutStudy.m
% D3: direction-robust sensor count and layout. D1 assumed one KNOWN wave heading (every optimum lies on the incidence
% axis). D3 removes that assumption in two separate ways and asks whether two sensors still suffice:
%   D3a  unknown but single heading (model correct; heading psi a nuisance parameter), 24 headings
%   D3b  multi-directional sea (frozen single-direction model wrong; V2's bias projected per layout), 7 seas x 24 headings
%
% METHOD (programme rev. 4; cfg.D3, clarification logged 2026-10-05 before the run)
%   Candidates, admissibility, reference optimisation, noise cases and class rules as D1. Headings 7.5 + 15k deg (no
%   heading on a grid axis: see the log). One emmField at the D1 grid and its copies rotated by delta = 7.5, 22.5, 37.5
%   deg and delta +- 0.5 deg; every heading is a 45 deg index permutation of one of the three delta (d3Tools).
%   D3a  per heading, theta = [ln beta, ln R, psi]: F_eff (Schur), A = tr F_eff^-1, sigma, d_W, kappa, class, best
%        reference per heading; per layout the MEAN and WORST-CASE A over headings (never A of a mean F), also with
%        the heading known (the nuisance penalty). Heading 0 known-heading search reproduces D1 (gate: n_s = 2 optimum
%        s6 + s18, reference s18 at LSM6DSV16X).
%   D3b  headline noise: seas cfg.D3.d3b (30/60/90 deg at 10/25%, cos^2s spread s = 10), each rotated with the primary
%        heading; Sigma, F and weights under the directional sea; d_tot from b_W + b_dir. Reference per condition: the
%        smallest A among references with d_tot <= 0.5 (gate), the smallest d_tot (least-bad diagnostic), and the
%        precision-optimal one (reported). Gate: d_tot <= 0.5 in every sea and heading; counts per n_s and per sea;
%        the gate is never relaxed. Feasible optima: mean A (primary ranking) and worst-case A (count rule).
%   Count rule (as D1, worst-case): the worst-case optimum (min over layouts of max A over headings) per n_s; a sensor
%   is worth adding if its worst-case sigma_lnbeta falls by >= 20%; recommended n_s = smallest whose worst-case optimum
%   is Identifiable at every heading (D3a), and with D3b also bias-feasible.
%   Earlier designs: the D1 optimum per n_s (heading 0, known) and Level 2C under every criterion, the D1 optimum's
%   efficiency against heading, and the heading-nuisance penalty sigma_lnbeta(psi estimated) / sigma_lnbeta(psi known).
%
% PRECONDITIONS. testDirectionRobustLayoutStudy passes; no uncommitted .m files.
% COST. One EMM field at 10 grid copies (330 points, 5 solves, cached); D3a 4 grid searches per noise case (each about
% D1's search), D3b 21 searches at the headline noise. Roughly 20 to 40 min in MATLAB.
% OUTPUT. Results/D3/ (log, .mat, summary, D3_d3a.csv, D3_d3b.csv, D3_counts.csv). Figures:
% plotDirectionRobustLayoutStudy (after the run).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

ALLOW_DIRTY = false;    % true: run with uncommitted .m files (provenance is then NOT the logged commit)

P = setupStudyPaths();
cfg = defineStudyScenarios();
T = d3Tools(); Tm = m1Tools();
V = cfg.D3;
run = startStudy('D3', 'Direction-robust sensor count and layout', P.results);

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
            error('D3:dirtyTree', '%d uncommitted .m file(s): commit first (or set ALLOW_DIRTY = true).', numel(dirty));
        end
        fprintf('  *** ALLOW_DIRTY: this run is NOT reproducible from commit %s ***\n', run.commit);
    else
        fprintf('  working tree clean (.m files): the code that runs is commit %s\n', run.commit);
    end
end

%% ---- 1  geometry, layouts, field --------------------------------------------------------------------------------------
p0 = cfg.p0.vec;
geom = T.geometry(cfg, p0);
nsList = V.sensorCounts; nK = numel(nsList);
Ls = cell(1, nK);
for k = 1:nK
    Ls{k} = T.layoutSet(geom, nsList(k), cfg);
    fprintf('  n_s = %d: %d admissible layouts (closed under 45 deg rotation)\n', nsList(k), Ls{k}.n);
end
hDeg = V.headingsDeg; nH = numel(hDeg);
deltas = unique(round(mod(hDeg, 45) * 1e9) / 1e9);
angles = T.copyAngles(cfg);
fprintf('\nStep 1: one field at %d grid copies (%s deg), %d points\n', numel(angles), mat2str(angles), ...
    numel(angles) * geom.perm(1, end));
F = T.buildField(cfg, p0, geom, angles, P.cache, [], true);
fprintf('  field ready (%.0f s); %d retained bins, N_eff %.1f\n', F.seconds, numel(F.omega), F.nEff);
L2C = sort(geom.level2C); k4 = find(nsList == numel(L2C), 1);
iL2C = []; if ~isempty(k4), iL2C = find(Ls{k4}.key == L2C * (34.^(0:numel(L2C) - 1)).', 1); end
names = @(lay) strjoin(geom.names(lay), '+');

%% ---- 2  D3a per noise case ------------------------------------------------------------------------------------------------
noiseCases = V.noiseCases; nN = numel(noiseCases);
A3 = struct('noise', noiseCases, 'ns', [], 'perNs', []);
for c = 1:nN
    Snn = cfg.noise.(noiseCases{c});
    fprintf('\nStep 2 (D3a, %s): heading 0 known (D1), then delta = %s deg with psi as nuisance\n', noiseCases{c}, ...
        mat2str(deltas));
    pre = T.precomputeA(F, 0, Snn, false);
    S0 = cell(1, nK);
    for k = 1:nK, S0{k} = T.searchA(pre, Ls{k}, cfg, false, [noiseCases{c} ' heading 0']); end
    perD = struct('deltaDeg', num2cell(deltas), 'S', []);
    SA = cell(numel(deltas), nK);
    for d = 1:numel(deltas)
        pre = T.precomputeA(F, deltas(d), Snn, true);
        for k = 1:nK, SA{d, k} = T.searchA(pre, Ls{k}, cfg, true, sprintf('%s delta %.1f', noiseCases{c}, deltas(d))); end
    end
    clear pre
    perNs = struct('ns', num2cell(nsList));
    for k = 1:nK
        L = Ls{k};
        for d = 1:numel(deltas), perD(d).S = SA{d, k}; end
        An = T.headingTable(perD, 'An', L, F, hDeg); Ak = T.headingTable(perD, 'Ak', L, F, hDeg);
        sBn = T.headingTable(perD, 'sBn', L, F, hDeg); sBk = T.headingTable(perD, 'sBk', L, F, hDeg);
        sRn = T.headingTable(perD, 'sRn', L, F, hDeg);
        clsN = T.headingTable(perD, 'clsN', L, F, hDeg); dWn = T.headingTable(perD, 'dWn', L, F, hDeg);
        kapN = T.headingTable(perD, 'kapN', L, F, hDeg);
        agg = struct('AbarN', mean(An, 2), 'AworstN', max(An, [], 2), 'sBworstN', max(sBn, [], 2), ...
            'sBmeanN', mean(sBn, 2), 'clsWorstN', max(clsN, [], 2), 'dWworstN', max(dWn, [], 2), ...
            'kapWorstN', max(kapN, [], 2), 'AbarK', mean(Ak, 2), 'AworstK', max(Ak, [], 2), 'sBworstK', max(sBk, [], 2));
        [~, iMean] = min(agg.AbarN); [~, iWorst] = min(agg.AworstN);
        [~, iD1] = min(S0{k}.Ak);                       % the D1 optimum: heading 0, known, A-best
        sel = struct('name', {'robustMean', 'robustWorst', 'D1'}, 'row', {iMean, iWorst, iD1});
        if k == k4 && ~isempty(iL2C), sel(end + 1) = struct('name', 'Level2C', 'row', iL2C); end %#ok<AGROW>
        for s = 1:numel(sel)
            i = sel(s).row;
            sel(s).layout = L.rows(i, :); %#ok<SAGROW>
            sel(s).AbarN = agg.AbarN(i); sel(s).AworstN = agg.AworstN(i); sel(s).sBworstN = agg.sBworstN(i);
            sel(s).sBmeanN = agg.sBmeanN(i); sel(s).sRworstN = max(sRn(i, :)); sel(s).clsWorstN = agg.clsWorstN(i);
            sel(s).dWworstN = agg.dWworstN(i); sel(s).kapWorstN = agg.kapWorstN(i);
            sel(s).AbarK = agg.AbarK(i); sel(s).AworstK = agg.AworstK(i); sel(s).sBworstK = agg.sBworstK(i);
            sel(s).effMean = agg.AbarN(iMean) / agg.AbarN(i); sel(s).effWorst = agg.AworstN(iWorst) / agg.AworstN(i);
            sel(s).perHeading = struct('An', An(i, :), 'Ak', Ak(i, :), 'sBn', sBn(i, :), 'sBk', sBk(i, :), ...
                'cls', clsN(i, :), 'effKnown', min(Ak, [], 1) ./ Ak(i, :));
            sel(s).heading0 = struct('A', S0{k}.Ak(i), 'sB', S0{k}.sBk(i), 'sR', S0{k}.sRk(i), 'dW', S0{k}.dWk(i), ...
                'cls', S0{k}.clsK(i), 'ref', S0{k}.refK(i));
        end
        perNs(k).sel = sel;
        perNs(k).nIdentAllHeadings = nnz(agg.clsWorstN == 1);
        perNs(k).worstSigmaDist = [min(agg.sBworstN), median(agg.sBworstN), max(agg.sBworstN)];
        fprintf('  n_s = %d (%d layouts, %d Identifiable at every heading):\n', nsList(k), L.n, perNs(k).nIdentAllHeadings);
        for s = 1:numel(sel)
            q = sel(s);
            fprintf(['    %-11s %-16s mean A %.3e, worst A %.3e | sigma_lnb worst %.3f%% (known %.3f%%), mean ' ...
                '%.3f%%; worst class %s; eff mean %.3f, worst %.3f\n'], q.name, names(q.layout), q.AbarN, q.AworstN, ...
                100 * q.sBworstN, 100 * q.sBworstK, 100 * q.sBmeanN, cfg.class.names{q.clsWorstN}, q.effMean, q.effWorst);
        end
    end
    % gate (headline): the heading-0 known search reproduces D1
    if c == 1
        q = perNs(nsList == numel(V.d1Gate.layout)).sel(3);
        ok = isequal(q.layout, V.d1Gate.layout) && q.heading0.ref == V.d1Gate.ref;
        fprintf('  D1 gate: heading-0 n_s = 2 optimum %s ref %s (D1: %s ref %s) %s\n', names(q.layout), ...
            geom.names{q.heading0.ref}, names(V.d1Gate.layout), geom.names{V.d1Gate.ref}, ternary(ok, 'OK', 'MISMATCH'));
        if ~ok, error('D3:d1Gate', 'The heading-0 known-heading search does not reproduce the D1 optimum.'); end
    end
    % count rule on the worst-case optimum
    sw = arrayfun(@(k) perNs(k).sel(2).sBworstN, 1:nK);
    cw = arrayfun(@(k) perNs(k).sel(2).clsWorstN, 1:nK);
    worth = [NaN, 1 - sw(2:end) ./ sw(1:end - 1) >= V.worthItFraction];
    rec = nsList(find(cw == 1, 1));
    if isempty(rec), rec = NaN; end
    A3(c).ns = nsList; A3(c).perNs = perNs;
    A3(c).count = struct('sBworst', sw, 'clsWorst', cw, 'gain', [NaN, 1 - sw(2:end) ./ sw(1:end - 1)], 'worth', worth, ...
        'recommended', rec);
    fprintf('  count rule (worst-case optimum): sigma_lnb worst %s %%, gain %s, recommended n_s (D3a) = %s\n', ...
        mat2str(100 * sw, 4), mat2str(A3(c).count.gain, 3), num2str(rec));
end

%% ---- 3  D3b at the headline noise -------------------------------------------------------------------------------------
Snn = cfg.noise.(V.d3bNoise);
seas = V.d3b; nS = numel(seas);
seaName = arrayfun(@(s) seaLabel(s), seas, 'UniformOutput', false);
fprintf('\nStep 3 (D3b, %s, ruler %s): %d seas x %d headings, gate d_tot <= %.2g\n', V.d3bNoise, V.d3bRuler, nS, nH, ...
    V.biasGate);
B3 = struct('ns', num2cell(nsList));
acc = cell(1, nK);
for k = 1:nK
    n = Ls{k}.n;
    acc{k} = struct('dTotWorstSea', NaN(n, nS), 'dTotAbestWorstSea', NaN(n, nS), 'dSysWorstSea', NaN(n, nS), ...
        'Asum', zeros(n, 1), 'Aworst', -Inf(n, 1), 'sBworst', -Inf(n, 1), 'clsWorst', zeros(n, 1), ...
        'Afsum', zeros(n, 1), 'Afworst', -Inf(n, 1), 'sBfworst', -Inf(n, 1), 'clsfWorst', zeros(n, 1), 'nCond', 0);
end
for s = 1:nS
    perD = struct('deltaDeg', num2cell(deltas), 'S', []);
    SB = cell(numel(deltas), nK);
    for d = 1:numel(deltas)
        pre = T.precomputeB(F, deltas(d), seas(s), Snn, V.d3bRuler);
        for k = 1:nK, SB{d, k} = T.searchB(pre, Ls{k}, cfg, sprintf('%s delta %.1f', seaName{s}, deltas(d))); end
    end
    clear pre
    for k = 1:nK
        for d = 1:numel(deltas), perD(d).S = SB{d, k}; end
        tb = @(f) T.headingTable(perD, f, Ls{k}, F, hDeg);
        a = acc{k};
        a.dTotWorstSea(:, s) = max(tb('dTotMin'), [], 2);        % least-biased reference per condition: the gate
        a.dTotAbestWorstSea(:, s) = max(tb('dTot'), [], 2);      % precision-optimal reference (reported)
        a.dSysWorstSea(:, s) = max(tb('dSys'), [], 2);
        A = tb('A'); a.Asum = a.Asum + sum(A, 2); a.Aworst = max(a.Aworst, max(A, [], 2));
        a.sBworst = max(a.sBworst, max(tb('sB'), [], 2)); a.clsWorst = max(a.clsWorst, max(tb('cls'), [], 2));
        Af = tb('Af'); a.Afsum = a.Afsum + sum(Af, 2); a.Afworst = max(a.Afworst, max(Af, [], 2));
        a.sBfworst = max(a.sBfworst, max(tb('sBf'), [], 2)); a.clsfWorst = max(a.clsfWorst, max(tb('clsf'), [], 2));
        a.nCond = a.nCond + nH;
        acc{k} = a;
        fprintf(['    %-22s n_s = %d: min over layouts of worst d_tot %.3g (best reference per condition; %.3g with ' ...
            'the precision-optimal one); %d layouts pass this sea\n'], seaName{s}, nsList(k), min(a.dTotWorstSea(:, s)), ...
            min(a.dTotAbestWorstSea(:, s)), nnz(a.dTotWorstSea(:, s) <= V.biasGate));
    end
end
for k = 1:nK
    a = acc{k}; L = Ls{k};
    worst = max(a.dTotWorstSea, [], 2);
    feas = worst <= V.biasGate;                                   % every sea and heading has a reference that passes
    [~, iLB] = min(worst);
    B3(k).nFeasible = nnz(feas); B3(k).nPassPerSea = sum(a.dTotWorstSea <= V.biasGate, 1);
    B3(k).leastBad = struct('layout', L.rows(iLB, :), 'dTotWorst', worst(iLB), 'dTotPerSea', a.dTotWorstSea(iLB, :), ...
        'dTotAbestPerSea', a.dTotAbestWorstSea(iLB, :), 'Abar', a.Asum(iLB) / a.nCond, 'sBworst', a.sBworst(iLB));
    B3(k).minWorstPerSea = min(a.dTotWorstSea, [], 1);
    B3(k).minWorstAbestPerSea = min(a.dTotAbestWorstSea, [], 1);
    B3(k).bestFeasibleMean = []; B3(k).bestFeasibleWorst = [];
    if any(feas)
        idx = find(feas);
        mk = @(i) struct('layout', L.rows(i, :), 'Abar', a.Afsum(i) / a.nCond, 'Aworst', a.Afworst(i), ...
            'sBworst', a.sBfworst(i), 'clsWorst', a.clsfWorst(i), 'dTotWorst', worst(i));
        [~, j] = min(a.Afsum(idx)); B3(k).bestFeasibleMean = mk(idx(j));     % primary D3b ranking
        [~, j] = min(a.Afworst(idx)); B3(k).bestFeasibleWorst = mk(idx(j));  % count rule
    end
    rowsSel = [A3(1).perNs(k).sel.row];
    namesSel = {A3(1).perNs(k).sel.name};
    for j = 1:numel(rowsSel)
        i = rowsSel(j);
        B3(k).sel(j) = struct('name', namesSel{j}, 'layout', L.rows(i, :), 'dTotWorst', worst(i), ...
            'dTotPerSea', a.dTotWorstSea(i, :), 'dTotAbestPerSea', a.dTotAbestWorstSea(i, :), ...
            'dSysPerSea', a.dSysWorstSea(i, :), 'Abar', a.Asum(i) / a.nCond, 'sBworst', a.sBworst(i), 'clsWorst', a.clsWorst(i));
    end
    fprintf(['  n_s = %d: %d of %d layouts bias-feasible in every sea and heading; per sea %s; least bad %s ' ...
        '(worst d_tot %.3g)\n'], nsList(k), B3(k).nFeasible, L.n, mat2str(B3(k).nPassPerSea), names(B3(k).leastBad.layout), ...
        B3(k).leastBad.dTotWorst);
    if ~isempty(B3(k).bestFeasibleWorst)
        q = B3(k).bestFeasibleWorst;
        fprintf('    feasible worst-case optimum %s: sigma_lnb worst %.3f%%, class %s\n', names(q.layout), ...
            100 * q.sBworst, cfg.class.names{q.clsWorst});
    end
end
% count rule with D3b: the bias-feasible WORST-CASE optimum (20% rule on its worst-case sigma_lnbeta)
swB = NaN(1, nK); clB = NaN(1, nK);
for k = 1:nK
    if ~isempty(B3(k).bestFeasibleWorst), swB(k) = B3(k).bestFeasibleWorst.sBworst; clB(k) = B3(k).bestFeasibleWorst.clsWorst; end
end
countB = struct('sBworst', swB, 'clsWorst', clB, 'gain', [NaN, 1 - swB(2:end) ./ swB(1:end - 1)]);
countB.worth = countB.gain >= V.worthItFraction;
recB = NaN;
for k = 1:nK
    if clB(k) == 1 && A3(1).count.clsWorst(k) == 1, recB = nsList(k); break; end
end
fprintf('\n  Recommended n_s: D3a (unknown heading) %s; with D3b (directional seas) %s\n', num2str(A3(1).count.recommended), ...
    ternary(isnan(recB), 'none up to 4 (no bias-feasible worst-case optimum is Identifiable)', num2str(recB)));

%% ---- 4  tables, summary --------------------------------------------------------------------------------------------------
rowsA = {};
for c = 1:nN
    for k = 1:nK
        for s = 1:numel(A3(c).perNs(k).sel)
            q = A3(c).perNs(k).sel(s);
            rowsA(end + 1, :) = {noiseCases{c}, nsList(k), q.name, names(q.layout), q.AbarN, q.AworstN, q.sBmeanN, ...
                q.sBworstN, q.sBworstK, q.sRworstN, q.dWworstN, cfg.class.names{q.clsWorstN}, q.effMean, q.effWorst, ...
                q.heading0.sB}; %#ok<SAGROW>
        end
    end
end
writeStudyCsv(run, 'd3a', {'noise', 'n_s', 'design', 'layout', 'A_mean', 'A_worst', 'sigma_lnbeta_mean', ...
    'sigma_lnbeta_worst', 'sigma_lnbeta_worst_known_heading', 'sigma_lnR_worst', 'd_W_worst', 'class_worst', ...
    'efficiency_mean', 'efficiency_worst', 'sigma_lnbeta_heading0_known'}, rowsA);
rowsB = {};
for k = 1:nK
    for j = 1:numel(B3(k).sel)
        q = B3(k).sel(j);
        rowsB(end + 1, :) = [{nsList(k), q.name, names(q.layout), q.dTotWorst, max(q.dTotAbestPerSea), q.sBworst, ...
            cfg.class.names{q.clsWorst}}, num2cell(q.dTotPerSea)]; %#ok<SAGROW>
    end
    q = B3(k).leastBad;
    rowsB(end + 1, :) = [{nsList(k), 'leastBad', names(q.layout), q.dTotWorst, max(q.dTotAbestPerSea), q.sBworst, '-'}, ...
        num2cell(q.dTotPerSea)]; %#ok<SAGROW>
    for f = {'bestFeasibleMean', 'bestFeasibleWorst'}
        q = B3(k).(f{1});
        if ~isempty(q)
            rowsB(end + 1, :) = [{nsList(k), f{1}, names(q.layout), q.dTotWorst, NaN, q.sBworst, ...
                cfg.class.names{q.clsWorst}}, num2cell(NaN(1, nS))]; %#ok<SAGROW>
        end
    end
end
writeStudyCsv(run, 'd3b', [{'n_s', 'design', 'layout', 'd_tot_worst_best_reference', 'd_tot_worst_precision_reference', ...
    'sigma_lnbeta_worst', 'class_worst'}, ...
    cellfun(@(x) ['d_tot_' x], matlab_safe(seaName), 'UniformOutput', false)], rowsB);
rowsC = {};
for c = 1:nN
    for k = 1:nK
        rowsC(end + 1, :) = {noiseCases{c}, nsList(k), A3(c).count.sBworst(k), A3(c).count.gain(k), A3(c).count.worth(k), ...
            cfg.class.names{A3(c).count.clsWorst(k)}, ternary(c == 1, B3(k).nFeasible, NaN), ...
            ternary(c == 1, countB.sBworst(k), NaN)}; %#ok<SAGROW>
    end
end
writeStudyCsv(run, 'counts', {'noise', 'n_s', 'worst_case_optimum_sigma_lnbeta_worst', 'gain', 'worth_it', ...
    'class_worst', 'n_bias_feasible_D3b', 'D3b_feasible_worst_case_optimum_sigma_lnbeta_worst'}, rowsC);

Ls2 = {};
Ls2{end + 1} = 'D3: direction-robust sensor count and layout';
Ls2{end + 1} = sprintf('run %s, %s, code commit %s', run.stamp, run.platform, run.commit);
Ls2{end + 1} = sprintf('headings %s deg; D3b seas: %s (%s, ruler %s)', mat2str(hDeg), strjoin(seaName, '; '), V.d3bNoise, ...
    V.d3bRuler);
Ls2{end + 1} = '';
for c = 1:nN
    Ls2{end + 1} = sprintf('D3a %s: recommended n_s %s; worst-case optimum sigma_lnbeta worst %s %%', noiseCases{c}, ...
        num2str(A3(c).count.recommended), mat2str(100 * A3(c).count.sBworst, 4)); %#ok<SAGROW>
    for k = 1:nK
        for s = 1:numel(A3(c).perNs(k).sel)
            q = A3(c).perNs(k).sel(s);
            Ls2{end + 1} = sprintf(['  n_s %d %-11s %-16s sigma_lnb worst %.3f%% (known %.3f%%), mean %.3f%%, class %s, ' ...
                'eff mean %.3f worst %.3f'], nsList(k), q.name, names(q.layout), 100 * q.sBworstN, 100 * q.sBworstK, ...
                100 * q.sBmeanN, cfg.class.names{q.clsWorstN}, q.effMean, q.effWorst); %#ok<SAGROW>
        end
    end
end
Ls2{end + 1} = sprintf('D3b (%s): bias-feasible layouts per n_s %s; recommended with D3b: %s', V.d3bNoise, ...
    mat2str([B3.nFeasible]), ternary(isnan(recB), 'none up to 4', num2str(recB)));
for k = 1:nK
    Ls2{end + 1} = sprintf('  n_s %d: per-sea passes %s; min worst d_tot per sea %s; least bad %s (%.3g)', nsList(k), ...
        mat2str(B3(k).nPassPerSea), mat2str(B3(k).minWorstPerSea, 3), names(B3(k).leastBad.layout), ...
        B3(k).leastBad.dTotWorst); %#ok<SAGROW>
end
results = struct('d3a', {A3}, 'd3b', {B3}, 'countD3b', countB, 'recommendedD3b', recB, 'headingsDeg', hDeg, 'deltas', deltas, ...
    'seaNames', {seaName}, 'nsList', nsList, 'geom', struct('frac', geom.frac, 'names', {geom.names}, ...
    'level2C', geom.level2C, 'Rphys', geom.Rphys, 'xy', geom.xy), 'nLayouts', cellfun(@(L) L.n, Ls), ...
    'fieldSeconds', F.seconds, 'tools', T.version);
finishStudy(run, results, cfg, Ls2);
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================
function s = seaLabel(sea)
if strcmp(sea.kind, 'two')
    s = sprintf('phi %g f2 %.2f', sea.phiDeg, sea.f2);
else
    s = sprintf('spread s %g', sea.s);
end
end

function c = matlab_safe(c)
c = cellfun(@(x) regexprep(x, '[^A-Za-z0-9]+', '_'), c, 'UniformOutput', false);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end