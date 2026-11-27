%% reclassifyStudies.m  --  F10 step 5: every saved study point under the thesis AND the accuracy class rule
%
% WHY. The thesis rule classifies on Fisher distance (bias / sigma). At p0 sigma_lnbeta is 0.05%, so a bias that
% is physically negligible (9 um of sensor position, 0.008 mm of draught) already fails it. The accuracy rule
% (Engine/identifiabilityClass, cfg.classAccuracy) classifies on e_beta = sqrt(sigma^2 + bias^2) and keeps
% d_tot <= 0.5 as a calibration flag. This script applies BOTH rules to the raw records the studies saved
% (rawOf, 2026-11-26 onwards), with no EMM solve, and reports where they differ.
%
% WHAT. For each study it loads the NEWEST Results/<ID>/<ID>_*.mat and
%   1  recomputes the THESIS class of every saved point and checks it against the saved class (regression:
%      a mismatch means the record is incomplete, and is reported, never patched);
%   2  computes the ACCURACY class, e_beta, the 1/5/10/20% ladder level and the calibration flag;
%   3  re-derives the study's headline quantity under the accuracy rule:
%        R1  noise window                    R3, R4  tolerance on the noise floor / gamma for e_beta <= 1, 5, 20%
%        V1  tolerance on depth / radial position (from every evaluated point, incl. the bisection steps)
%        V2  second-direction energy fraction f2 per angle
%        D3  D3a layouts accurate at every heading; D3b layouts accurate in every sea and heading
%        A1  class of each floe (V1a truncation rows use the A1 ruler of the same floe)
%   Crossings between evaluated points are interpolated in log|x| vs log e_beta and printed with their bracket.
%
% NOT INCLUDED. D1 (layout ranking does not use the class), D2 and R2 (no raw bias vectors saved), M1.
%
% OUTPUT. Printed, and Results/F10/reclassify_<stamp>.log, .mat (all tables) and _points.csv (one row per point).
% Read-only on the study results. Runs in seconds.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.   (F10 step 5, 2026-11-27)

P = setupStudyPaths();
cfg = defineStudyScenarios();
A = cfg.classAccuracy;
outDir = fullfile(P.results, 'F10');
if ~exist(outDir, 'dir'), mkdir(outDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS');
logFile = fullfile(outDir, sprintf('reclassify_%s.log', stamp));
diary(logFile); diary('on');
try
fprintf('\n==== F10: thesis vs accuracy classes from the saved study records (%s) ====\n', stamp);
fprintf(['  accuracy rule: e_beta = sqrt(sigma_lnbeta^2 + bias_lnbeta^2); Identifiable <= %g%%, Marginal <= %g%%, ' ...
    'kappa <= %g;\n  calibrated = d_tot <= %g; ladder %s%%\n'], 100 * A.identifiableMaxE, 100 * A.marginalMaxE, ...
    A.maxCondF, A.calibrationMaxD, mat2str(100 * A.ladder));
PTS = {};                                     % one row per point: study, point, saved, thesis, accuracy, e, ladder, cal, dTot
OUT = struct();
lv = [A.ladder(1), A.identifiableMaxE, A.marginalMaxE];   % thresholds reported for tolerances: 1%, 5%, 20%

%% ---- A1 (first: V1a borrows its rulers) -----------------------------------------------------------------------
[S, f] = newest(P, 'A1');
a1 = struct('name', {}, 'raw', {});
if ~isempty(S)
    sectionHead('A1', f);
    R = S.results.res;
    fprintf('  %-26s %-17s %-17s %8s %7s %5s %8s  %s\n', 'scenario', 'thesis', 'accuracy', 'e_beta', 'ladder', 'cal', ...
        'd_tot', 'model check');
    for i = 1:numel(R)
        r = R{i}; if isempty(r) || ~isfield(r, 'sigma'), continue; end
        raw = struct('sigma', r.sigma(:), 'singular', r.singular, 'condF', r.condF, 'sMax', r.sMax, 'dTot', r.dTot, ...
            'dW', r.dW, 'dthetaW', r.dthetaW(:));
        [PTS, c] = addPoint(PTS, 'A1', r.name, r.inverseClass, raw, cfg);
        a1(end + 1) = struct('name', r.name, 'raw', raw); %#ok<SAGROW>
        note = ''; if ~r.modelValid, note = 'FAILED (class robust: inverse class shown)'; end
        fprintf('  %-26s %-17s %-17s %7.3f%% %7s %5s %8.3f  %s\n', r.name, c.thesis, c.acc, 100 * c.e, ladStr(c.ladder), ...
            yn(c.cal), r.dTot, note);
    end
    OUT.A1 = summarise(PTS, 'A1', cfg);
end

%% ---- R1 -----------------------------------------------------------------------------------------------------------
[S, f] = newest(P, 'R1');
if ~isempty(S) && isfield(S.results.grid, 'raw')
    sectionHead('R1', f);
    G = S.results.grid; Snn = S.results.Snn;
    e = NaN(1, numel(G.raw)); code = NaN(1, numel(G.raw));
    for k = 1:numel(G.raw)
        [PTS, c] = addPoint(PTS, 'R1', sprintf('S_nn %.4g', Snn(k)), G.raw{k}.class, G.raw{k}, cfg);
        e(k) = c.e; code(k) = c.accCode;
    end
    fprintf('  accuracy-Identifiable for %d of %d noise levels (thesis: %d); e_beta from %.3f%% to %.3f%%\n', ...
        nnz(code == 1), numel(code), nnz(strcmp(cellfun(@(r) r.class, G.raw, 'UniformOutput', false), cfg.class.names{1})), ...
        100 * min(e), 100 * max(e));
    OUT.R1 = summarise(PTS, 'R1', cfg);
    OUT.R1.eBeta = e; OUT.R1.Snn = Snn;
end

%% ---- R3 and R4: tolerances, exact on the fine grid ----------------------------------------------------------------
for id = {'R3', 'R4'}
    id = id{1}; %#ok<FXSET>
    [S, f] = newest(P, id);
    if isempty(S) || ~isfield(S.results.res, 'raw'), continue; end
    sectionHead(id, f);
    res = S.results.res; epsGrid = S.results.epsGrid;
    if strcmp(id, 'R3')
        epsFine = linspace(cfg.R3.toleranceSearch(1), cfg.R3.toleranceSearch(2), cfg.R3.toleranceN);
        what = 'noise-floor error eps';
    else
        epsFine = linspace(cfg.R4.toleranceSearch(1), cfg.R4.toleranceSearch(2), cfg.R4.toleranceN);
        what = 'gamma (draught) error eps';
    end
    fprintf('  tolerance on the %s (smaller |eps| of the two signs; "-" = not reached in the search range)\n', what);
    fprintf('  %-11s %10s %10s %10s | %10s %10s %10s\n', 'noise', 'e<=1%', 'e<=5%', 'e<=20%', '|db|<=1%', 'd<=0.5', 'd<=1');
    tolTab = struct('noise', {}, 'tol', {}, 'thesis', {});
    for n = 1:numel(res)
        r = res(n); raw = r.raw;
        for k = 1:numel(epsGrid)
            rk = raw; rk.dTot = r.dTot(k); rk.dthetaTot = r.dthetaTot(:, k);
            PTS = addPoint(PTS, id, sprintf('%s eps %+.0f%%', r.noise, 100 * epsGrid(k)), r.class{k}, rk, cfg);
        end
        if strcmp(id, 'R3')
            if ~isfield(r, 'dthetaSysFine'), continue; end
            bFine = raw.dthetaW(1) + r.dthetaSysFine(1, :);
        else
            bFine = raw.dthetaW(1) + r.slope(1) * log(1 + epsFine);
        end
        eFine = sqrt(raw.sigma(1)^2 + bFine.^2);
        tol = arrayfun(@(t) firstCrossingAbs(epsFine, eFine, t), lv);
        tolTab(end + 1) = struct('noise', r.noise, 'tol', tol, 'thesis', [r.tolBeta1 r.tolD05 r.tolD1]); %#ok<SAGROW>
        fprintf('  %-11s %10s %10s %10s | %10s %10s %10s\n', r.noise, pct(tol(1)), pct(tol(2)), pct(tol(3)), ...
            pct(r.tolBeta1), pct(r.tolD05), pct(r.tolD1));
    end
    OUT.(id) = summarise(PTS, id, cfg);
    OUT.(id).tolerances = tolTab;
end

%% ---- V1: depth and radial position, from every evaluated point ----------------------------------------------------
[S, f] = newest(P, 'V1');
if ~isempty(S) && isfield(S.results, 'evaluated')
    sectionHead('V1', f);
    nom = S.results.nominal.raw;
    ev = S.results.evaluated;
    kind = cell(1, numel(ev.keys)); x = NaN(1, numel(ev.keys)); e = x;
    for j = 1:numel(ev.keys)
        rec = ev.records{j};
        u = find(ev.keys{j} == '_', 1, 'last');
        kind{j} = ev.keys{j}(1:u - 1); x(j) = str2double(ev.keys{j}(u + 1:end));
        rk = nom; rk.dTot = rec.dTot; rk.dthetaTot = rec.dthetaTot(:);
        [PTS, c] = addPoint(PTS, 'V1', sprintf('%s %+.6g', kind{j}, x(j)), rec.class, rk, cfg);
        e(j) = c.e;
    end
    % V1a: each floe's truncation bias on that floe's own ruler (from A1)
    if ~isempty(a1)
        fprintf('  V1a truncation (ruler = the A1 record of the same floe):\n');
        rows = S.results.truncation.rows; fl = S.results.truncation.floes;
        for i = 1:numel(rows)
            m = find(strcmp({a1.name}, fl{i}), 1);
            if isempty(m), fprintf('    %-24s no A1 record: skipped\n', fl{i}); continue; end
            rk = a1(m).raw; rk.dTot = rows{i}.dTot; rk.dthetaTot = rows{i}.dthetaTot(:);
            [PTS, c] = addPoint(PTS, 'V1a', fl{i}, rows{i}.class, rk, cfg);
            fprintf('    %-24s thesis %-17s accuracy %-17s e_beta %.3f%%, calibrated %s\n', fl{i}, c.thesis, c.acc, ...
                100 * c.e, yn(c.cal));
        end
    end
    sc = struct('depth', struct('unit', '%', 'scale', 100), 'radial', struct('unit', 'mm', 'scale', 1000));
    fprintf('  tolerance for e_beta <= 1 / 5 / 20%% (first crossing per side, log-interpolated between evaluated points)\n');
    tolV1 = struct('kind', {}, 'side', {}, 'tol', {}, 'bracket', {});
    for kd = {'depth', 'radial'}
        kd = kd{1}; %#ok<FXSET>
        for side = [-1 1]
            sel = strcmp(kind, kd) & (sign(x) == side | x == 0);
            [tol, br] = arrayfun(@(t) crossOutward(abs(x(sel)), e(sel), t), lv, 'UniformOutput', false);
            tol = cell2mat(tol);
            tolV1(end + 1) = struct('kind', kd, 'side', side, 'tol', side * tol, 'bracket', {br}); %#ok<SAGROW>
            fprintf('    %-6s %s  ', kd, ternary(side < 0, '-', '+'));
            for t = 1:numel(lv)
                fprintf('e<=%2.0f%%: %s   ', 100 * lv(t), tolStr(side * tol(t) * sc.(kd).scale, sc.(kd).unit, ...
                    sc.(kd).scale * br{t}));
            end
            fprintf('\n');
        end
    end
    T = S.results.tolerances;
    fprintf('  thesis tolerances (saved): ');
    for t = 1:numel(T)
        if T(t).crossed, fprintf('%s %s %s %g: %s; ', T(t).kind, ternary(T(t).side < 0, '-', '+'), T(t).criterion, ...
                T(t).threshold, T(t).note); end
    end
    fprintf('\n');
    OUT.V1 = summarise(PTS, 'V1', cfg);
    OUT.V1.tolerances = tolV1;
end

%% ---- V2: directional forcing ---------------------------------------------------------------------------------------
[S, f] = newest(P, 'V2');
if ~isempty(S) && isfield(S.results, 'evaluated')
    sectionHead('V2', f);
    nom = S.results.nominal.raw;
    ev = S.results.evaluated;
    phi = NaN(1, numel(ev.keys)); f2 = phi; e = phi; dB = phi;
    for j = 1:numel(ev.keys)
        rec = ev.records{j};
        v = sscanf(strrep(ev.keys{j}, '_', ' '), '%f %f');
        phi(j) = round(v(1) * 180 / pi * 1e6) / 1e6; f2(j) = v(2);
        rk = nom; rk.dTot = rec.dTot; rk.dthetaTot = rec.dthetaTot(:);
        [PTS, c] = addPoint(PTS, 'V2', sprintf('phi %g f2 %.4g', phi(j), f2(j)), rec.class, rk, cfg);
        e(j) = c.e; dB(j) = rec.betaPct;
    end
    fprintf('  largest second-direction energy fraction f2 with e_beta <= 1 / 5 / 20%% (thesis: d_sys <= 0.5, |dbeta| <= 1%%)\n');
    X = S.results.crossings;
    tolV2 = struct('phiDeg', {}, 'tol', {}, 'bracket', {});
    for ph = unique(phi(isfinite(phi)))
        sel = phi == ph & f2 > 0;
        [tol, br] = arrayfun(@(t) crossOutward(f2(sel), e(sel), t), lv, 'UniformOutput', false);
        tol = cell2mat(tol);
        tolV2(end + 1) = struct('phiDeg', ph, 'tol', tol, 'bracket', {br}); %#ok<SAGROW>
        th = X([X.phiDeg] == ph);
        fprintf('    phi %3g:  e<=1%%: %-10s e<=5%%: %-10s e<=20%%: %-10s | thesis d_sys 0.5: %-10s |db| 1%%: %s\n', ph, ...
            numStr(tol(1)), numStr(tol(2)), numStr(tol(3)), crossNote(th, 'd_sys', 0.5), crossNote(th, '|dbeta/beta|', 1));
    end
    OUT.V2 = summarise(PTS, 'V2', cfg);
    OUT.V2.tolerances = tolV2;
end

%% ---- D3: layouts (worst over headings / seas), vectorised on the saved per-layout records ---------------------------
[S, f] = newest(P, 'D3');
if ~isempty(S) && isfield(S.results.d3b, 'perLayout')
    sectionHead('D3', f);
    ns = S.results.nsList;
    fprintf('  D3a (unknown heading): layouts Identifiable at every heading, thesis vs accuracy (best worst-heading e_beta)\n');
    for c = 1:numel(S.results.d3a)
        a3 = S.results.d3a(c);
        for k = 1:numel(ns)
            L = a3.perNs(k).perLayout;
            nTh = nnz(L.clsWorstN == 1);
            reg = nTh == a3.perNs(k).nIdentAllHeadings;
            cAcc = accCodes(L.eBworstN, L.kapWorstN, A);
            [eMin, iMin] = min(L.eBworstN);
            fprintf('    %-11s n_s = %d: thesis %5d, accuracy %5d of %5d (calibrated as well: %5d); best e_beta %.3f%% %s\n', ...
                a3.noise, ns(k), nTh, nnz(cAcc == 1), numel(cAcc), nnz(cAcc == 1 & L.dWworstN <= A.calibrationMaxD), ...
                100 * eMin, ternary(reg, '', '[THESIS COUNT MISMATCH]'));
            OUT.D3.d3a(c, k) = struct('noise', a3.noise, 'ns', ns(k), 'thesis', nTh, 'accuracy', nnz(cAcc == 1), ...
                'accuracyCalibrated', nnz(cAcc == 1 & L.dWworstN <= A.calibrationMaxD), 'eBest', eMin, 'iBest', iMin, ...
                'regressionOK', reg);
        end
    end
    fprintf(['  D3b (directional seas, every sea and heading): thesis gate d_tot <= 0.5 (best reference per condition) vs ' ...
        'accuracy\n  (e_beta at the A-best reference; kappa is not stored per layout for D3b, D3a kappa worst <= %.3g)\n'], ...
        max(arrayfun(@(c) max(arrayfun(@(k) max(S.results.d3a(c).perNs(k).perLayout.kapWorstN), 1:numel(ns))), ...
        1:numel(S.results.d3a))));
    for k = 1:numel(ns)
        B = S.results.d3b(k).perLayout;
        worst = max(B.dTotWorstSea, [], 2);
        nTh = nnz(worst <= cfg.D3.biasGate);
        reg = nTh == S.results.d3b(k).nFeasible;
        cAcc = accCodes(B.eBworst, zeros(size(B.eBworst)), A);
        [eMin, iMin] = min(B.eBworst);
        fprintf('    n_s = %d: thesis feasible %5d, accuracy Identifiable %5d, Marginal %5d of %5d; best e_beta %.3f%% %s\n', ...
            ns(k), nTh, nnz(cAcc == 1), nnz(cAcc == 2), numel(cAcc), 100 * eMin, ternary(reg, '', '[THESIS COUNT MISMATCH]'));
        OUT.D3.d3b(k) = struct('ns', ns(k), 'thesis', nTh, 'accuracy', nnz(cAcc == 1), 'accuracyMarginal', nnz(cAcc == 2), ...
            'eBest', eMin, 'iBest', iMin, 'regressionOK', reg);
    end
end

%% ---- summary ----------------------------------------------------------------------------------------------------------
fprintf('\nSUMMARY  (points per study; "changed" = accuracy class differs from the thesis class)\n');
fprintf('  %-5s %6s %18s %18s %8s %10s %12s\n', 'study', 'points', 'thesis I/M/N', 'accuracy I/M/N', 'changed', ...
    'calibrated', 'thesis regr.');
ids = unique(PTS(:, 1), 'stable');
for i = 1:numel(ids)
    q = PTS(strcmp(PTS(:, 1), ids{i}), :);
    cnt = @(col) arrayfun(@(k) nnz(strcmp(q(:, col), cfg.class.names{k})), 1:3);
    t = cnt(4); a = cnt(5);
    fprintf('  %-5s %6d %18s %18s %8d %10d %12s\n', ids{i}, size(q, 1), sprintf('%d/%d/%d', t), sprintf('%d/%d/%d', a), ...
        nnz(~strcmp(q(:, 4), q(:, 5))), nnz([q{:, 8}]), ...
        ternary(all(strcmp(q(:, 3), q(:, 4)) | cellfun(@isempty, q(:, 3))), 'OK', ...
        sprintf('%d MISMATCH', nnz(~strcmp(q(:, 3), q(:, 4)) & ~cellfun(@isempty, q(:, 3))))));
end

hdr = {'study', 'point', 'savedClass', 'thesisClass', 'accuracyClass', 'eBeta', 'ladder', 'calibrated', 'dTot'};
Tpts = cell2table(PTS, 'VariableNames', hdr);
writetable(Tpts, fullfile(outDir, sprintf('reclassify_%s_points.csv', stamp)));
save(fullfile(outDir, sprintf('reclassify_%s.mat', stamp)), 'OUT', 'Tpts', 'cfg', '-v7');
fprintf('\n  saved %s (.log, .mat, _points.csv)\n', fullfile(outDir, sprintf('reclassify_%s', stamp)));
diary('off');
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================================
function [S, f] = newest(P, id)
L = dir(fullfile(P.results, id, [id '_20*.mat']));
S = []; f = '';
if isempty(L), fprintf('\n---- %s: no saved run\n', id); return; end
[~, k] = max([L.datenum]); f = L(k).name;
S = load(fullfile(L(k).folder, f), 'results');
end

function sectionHead(id, f)
fprintf('\n---- %s  (%s) ----\n', id, f);
end

function [PTS, c] = addPoint(PTS, study, label, saved, raw, cfg)
[c.thesis] = identifiabilityClass(raw, cfg, 'thesis');
[c.acc, c.accCode, info] = identifiabilityClass(raw, cfg, 'accuracy');
c.e = info.eBeta; c.ladder = info.ladder; c.cal = info.calibrated;
if isempty(saved), saved = ''; end
PTS(end + 1, :) = {study, label, saved, c.thesis, c.acc, c.e, c.ladder, c.cal, info.dTot};
end

function s = summarise(PTS, study, cfg)
q = PTS(strcmp(PTS(:, 1), study), :);
s = struct('n', size(q, 1), 'thesis', arrayfun(@(k) nnz(strcmp(q(:, 4), cfg.class.names{k})), 1:3), ...
    'accuracy', arrayfun(@(k) nnz(strcmp(q(:, 5), cfg.class.names{k})), 1:3), ...
    'changed', nnz(~strcmp(q(:, 4), q(:, 5))), 'regressionMismatch', nnz(~strcmp(q(:, 3), q(:, 4)) & ~cellfun(@isempty, q(:, 3))));
end

function c = accCodes(e, kappa, A)
% the accuracy rule, vectorised (identical to identifiabilityClass for finite inputs; NaN/Inf e -> 3)
c = 3 * ones(size(e));
ok = isfinite(e) & kappa <= A.maxCondF;
c(ok & e <= A.marginalMaxE) = 2;
c(ok & e <= A.identifiableMaxE) = 1;
end

function t = firstCrossingAbs(x, v, thr)
% smallest |x| at which v first reaches thr going outward from x = 0 on either side (linear interpolation)
t = NaN; best = Inf;
[~, i0] = min(abs(x));
for dir = [-1 1]
    idx = i0:dir:ternary(dir < 0, 1, numel(x));
    for j = 2:numel(idx)
        a = idx(j - 1); b = idx(j);
        if v(a) < thr && v(b) >= thr
            xc = x(a) + (thr - v(a)) * (x(b) - x(a)) / (v(b) - v(a));
            if abs(xc) < best, best = abs(xc); t = abs(xc); end
            break
        end
    end
end
end

function [t, br] = crossOutward(ax, v, thr)
% first crossing of thr by v as |x| grows (ax = |x| >= 0), log-log interpolated between evaluated points
[ax, o] = sort(ax(:).'); v = v(o);
t = NaN; br = [NaN NaN];
for j = 2:numel(ax)
    if v(j - 1) < thr && v(j) >= thr
        br = [ax(j - 1) ax(j)];
        if ax(j - 1) > 0 && v(j - 1) > 0
            t = exp(log(ax(j - 1)) + (log(thr) - log(v(j - 1))) * (log(ax(j)) - log(ax(j - 1))) / (log(v(j)) - log(v(j - 1))));
        else
            t = ax(j - 1) + (thr - v(j - 1)) * (ax(j) - ax(j - 1)) / (v(j) - v(j - 1));
        end
        return
    end
end
if ~isempty(v) && all(v >= thr), t = 0; br = [0 ax(1)]; end   % already exceeded at the smallest evaluated point
end

function s = crossNote(th, crit, thr)
s = '-';
k = find(strcmp({th.criterion}, crit) & [th.threshold] == thr, 1);
if ~isempty(k) && th(k).crossed, s = numStr(th(k).f2); end
end

function s = tolStr(v, unit, br)
if ~isfinite(v), s = 'not reached'; return; end
s = sprintf('%+.3g %s [%.3g, %.3g]', v, unit, abs(br(1)), abs(br(2)));
end

function s = pct(v)
if ~isfinite(v), s = '-'; else, s = sprintf('%.2f%%', 100 * v); end
end

function s = numStr(v)
if ~isfinite(v), s = '-'; else, s = sprintf('%.4g', v); end
end

function s = ladStr(v)
if ~isfinite(v), s = '>20%'; else, s = sprintf('%g%%', 100 * v); end
end

function s = yn(t)
if t, s = 'yes'; else, s = 'no'; end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
