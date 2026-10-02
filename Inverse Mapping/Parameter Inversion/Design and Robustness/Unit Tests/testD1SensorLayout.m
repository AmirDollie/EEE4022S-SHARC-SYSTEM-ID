%% testD1SensorLayout.m
% Pre-flight checks for study D1. Must pass (all) before runD1SensorLayout is run.
% The numerical core (d1LayoutTools) is the one the runner uses, so these checks test the search itself.
%
% The decisive check is 5: the Level 2C layout [s2 s6 s14 s26] with reference s26, EXTRACTED from the
% 33-point machinery, must reproduce evaluateScenario (the frozen E1 chain). If it does not, the
% exhaustive search would be ranking the wrong covariance or Jacobian sub-matrices.
%
% COST. The first run solves the 33-point EMM field (5 evaluations, a few minutes in MATLAB; cached in
% Results/twinCache, so later runs and runD1SensorLayout load it). The rest takes seconds.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
P = setupStudyPaths();
cfg = defineStudyScenarios();
T = d1LayoutTools();
nPass = 0; nFail = 0;
p0 = cfg.p0.vec;

%% ---- geometry and admissibility (no EMM) ---------------------------------------------------------------
geom = T.buildGeometry(cfg, p0);
ok = size(geom.xy, 1) == 33 && isequal(geom.level2C, [2 6 14 26]) && geom.level2CRef == 26;
[nPass, nFail] = check(ok, sprintf('1  33 candidate points; Level 2C found on the grid as %s, reference s%d', ...
    mat2str(geom.level2C), geom.level2CRef), nPass, nFail);

n = zeros(1, 3); nPrim = n; n15 = n;
for k = 1:3
    [L, okU] = T.enumerateLayouts(geom, k + 1, cfg, 0);
    [~, okP] = T.enumerateLayouts(geom, k + 1, cfg, cfg.D1.minSpacing);
    [~, ok15] = T.enumerateLayouts(geom, k + 1, cfg, 0.15);
    n(k) = nnz(okU); nPrim(k) = nnz(okP); n15(k) = nnz(ok15);
end
[nPass, nFail] = check(isequal(n, [528 5456 40920]), sprintf('2  unconstrained layout counts %s (528, 5456, 40920)', ...
    mat2str(n)), nPass, nFail);
dd = geom.dist + diag(Inf(33, 1));
ok = isequal(nPrim, n) && abs(min(dd(:)) - 0.144) < 5e-4 && all(n15 < n);
[nPass, nFail] = check(ok, sprintf(['3  primary 0.10 m spacing removes no layout at p0 (closest pair %.3f m); ' ...
    '0.15 m does bind (%s admissible)'], min(dd(:)), mat2str(n15)), nPass, nFail);

%% ---- the 33-point field and the per-reference precompute -----------------------------------------------
fprintf('  ... 33-point EMM field (5 evaluations; "computed" takes minutes the first time, then "loaded")\n');
full = T.buildFullField(cfg, p0, geom, P.cache, true);
Snn = cfg.noise.(cfg.D1.noiseCases{1});
tic; pre = T.precompute(full, Snn, cfg); tPre = toc;
ok = numel(pre) == 33 && all([pre.ref] == 1:33) && all(arrayfun(@(q) numel(q.others) == 32, pre));
[nPass, nFail] = check(ok, sprintf('4  precompute: 33 references x 32 transmissibilities, %d to %d accepted bins (%.1f s)', ...
    min([pre.nValid]), max([pre.nValid]), tPre), nPass, nFail);

%% ---- GATE: Level 2C reproduces the frozen E1 chain ----------------------------------------------------
e1 = evaluateScenario(struct('Snn', Snn, 'cacheDir', P.cache), cfg);
qd = T.evaluateLayout(pre, geom.level2C, 26, cfg, 'direct');
qf = T.evaluateLayout(pre, geom.level2C, 26, cfg, 'fast');
rs = max(abs([qd.sigma; qf.sigma] ./ [e1.sigma; e1.sigma] - 1));
rd = max(abs([qd.dW qf.dW] / e1.dW - 1));
rF = max(norm(qd.F - e1.F), norm(qf.F - e1.F)) / norm(e1.F);
ok = rs < 1e-6 && rd < 1e-6 && rF < 1e-6 && strcmp(qf.class, e1.class);
[nPass, nFail] = check(ok, sprintf(['5  GATE Level 2C + s26 = evaluateScenario: sigma [%.5f %.6f] vs [%.5f %.6f], d_W %.4f vs %.4f ' ...
    '(rel. diff sigma %.1e, d_W %.1e, F %.1e), %s'], qf.sigma, e1.sigma, qf.dW, e1.dW, rs, rd, rF, qf.class), nPass, nFail);

%% ---- bookkeeping ---------------------------------------------------------------------------------------
qp = T.evaluateLayout(pre, [14 2 26 6], 26, cfg, 'fast');
rp = norm(qp.F - qf.F) / norm(qf.F);
[nPass, nFail] = check(rp < 1e-12, sprintf('6  permuting the non-reference sensors leaves F unchanged (rel. diff %.1e)', rp), ...
    nPass, nFail);

C = qf.F \ eye(2);
ok = abs(qf.A - trace(C)) < 1e-12 * qf.A && abs(qf.A - sum(qf.sigma.^2)) < 1e-12 * qf.A && ...
    abs(qf.D - log(det(qf.F))) < 1e-9 * abs(qf.D) && abs(qf.E - min(eig(qf.F))) < 1e-9 * qf.E;
[nPass, nFail] = check(ok, sprintf('7  A = tr(F^-1) = sigma_lnbeta^2 + sigma_lnR^2 = %.3e; D = log det F; E = lambda_min(F)', qf.A), ...
    nPass, nFail);

[L3, ~] = T.enumerateLayouts(geom, 3, cfg, 0);
rows = [1:15, 2000:2010].';
S = T.searchLayouts(pre, L3(rows, :), cfg, 'test');
okRef = all(~isnan(S.sBall(:))) && size(S.Aall, 2) == 3 && all(abs(S.A - min(S.Aall, [], 2)) <= 1e-15 * S.A) && ...
    all(arrayfun(@(i) any(L3(rows(i), :) == S.refA(i)), 1:numel(rows)));
row4 = find(all(nchoosek(1:33, 4) == [2 6 14 26], 2));
S4 = T.searchLayouts(pre, [2 6 14 26], cfg, 'test');
okRef = okRef && abs(S4.Aall(1, 4) - qf.A) < 1e-12 * qf.A && ~isempty(row4);
[nPass, nFail] = check(okRef, sprintf(['8  reference loop: exactly n_s references evaluated per layout, best = argmin A; ' ...
    'Level 2C A by reference [%s] x 1e-5, best s%d'], sprintf(' %.3f', 1e5 * S4.Aall), S4.refA), nPass, nFail);

preS = pre; preS(26).dT(:, :, 1) = 0;              % no beta information at all
qs = T.evaluateLayout(preS, geom.level2C, 26, cfg, 'fast');
qs2 = T.evaluateLayout(preS, geom.level2C, 26, cfg, 'direct');
ok = qs.singular && isinf(qs.A) && qs.D == -Inf && qs.E == 0 && strcmp(qs.class, cfg.class.names{3}) && qs2.singular;
[nPass, nFail] = check(ok, '9  a layout with no beta information: singular F, A = Inf, Not identifiable (fast and direct), not regularised', ...
    nPass, nFail);

rng(3); worst = 0; worstD = 0; tic; nE = 0;
for ns = 2:4
    [Ln, ~] = T.enumerateLayouts(geom, ns, cfg, 0);
    pick = randperm(size(Ln, 1), 20);
    for i = pick
        r = Ln(i, randi(ns));
        a = T.evaluateLayout(pre, Ln(i, :), r, cfg, 'fast');
        b = T.evaluateLayout(pre, Ln(i, :), r, cfg, 'direct');
        nE = nE + 1;
        if a.singular ~= b.singular, worst = Inf; continue; end
        if ~a.singular
            worst = max(worst, norm(a.F - b.F) / norm(b.F));
            worstD = max(worstD, abs(a.dW / b.dW - 1));
        end
    end
end
[nPass, nFail] = check(worst < 1e-9 && worstD < 1e-8, sprintf(['10 fast (paged) and direct (fisherFromBlocks) evaluations agree on %d random ' ...
    'layouts: F %.1e, d_W %.1e'], nE, worst, worstD), nPass, nFail);

% speed estimate for the full search
L4 = nchoosek(1:33, 4);
tic; Ssp = T.searchLayouts(pre, L4(1:2000, :), cfg, 'timing'); %#ok<NASGU>
t4 = toc;
fprintf('\n  timing: 2000 four-sensor layouts x 4 references in %.1f s; full D1 (3 noise cases, n_s = 2,3,4) about %.0f min\n', ...
    t4, 3 * t4 * (528 * 2 + 5456 * 3 + 40920 * 4) / (2000 * 4) / 60);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testD1SensorLayout:failed', '%d check(s) failed: do NOT run the D1 search.', nFail); end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end