%% diagnoseV1Truncation.m
% V1 STAGE 0 (a diagnostic, not the V1 study): which EMM truncation is converged enough to serve as the "truth"
% when V1 asks what the frozen [50 10 10] inverse model gets wrong?
%
% WHY. A1 found that [50 10 10] gives converged METRICS but not a converged MODEL: at p0 the features differ
% from [70 15 15] by d = 0.57 and [70 15 15] from [90 20 20] by d = 0.25, dominated by the number of vertical
% modes M. [70 15 15] is therefore not a valid truth. The rule (programme rev. 5, cfg.V1): extend M until one
% consecutive step changes the features by d < cfg.V1.truncationConvergedD; the higher member of that first
% pair becomes cfg.V1.truncationTruth, logged in defineStudyScenarios BEFORE runV1ModelMismatch.
%
% WHAT IT COMPUTES (Level 2C, reference s26, LSM6DSV16X, setting F, node grid 0.2 rad/s; every difference is
% measured with the SAME ruler: the scenario's own F, J, Sigma at the production truncation [50 10 10],
% d = sqrt(dth' F dth), dth = F^-1 J' Sigma^-1 b, b = exact-feature difference)
%   1. The ladder at p0 (the registered selection case): [90 20 20] (cached from diagnoseTruncationP0), then
%      [M P N] with M = cfg.V1.truncationLadderM(2:end) and [P N] = cfg.V1.truncationLadderPN, stopping at the
%      first consecutive step with d < cfg.V1.truncationConvergedD (or at the end of the list: not converged).
%   2. Diagnostic only (does not change the selection): the same last step at the large-kR A1 floes
%      cfg.V1.truncationFloes, so a truth that is converged at p0 but not there is reported, not hidden.
%   3. Diagnostic only: the P index at the selected M, [M* P N] -> [M* cfg.V1.truncationCheckP N] at p0.
%   4. Preview of V1a: [50 10 10] against the selected truth at every floe (exact features; the V1 study itself
%      uses Welch expectation minus Welch expectation).
%
% COST. One nominal field (29 nodes, no derivatives) per level and floe. Each level is timed and the time of
% the next is projected (the EMM cost grows with M). Every field is cached in Results/twinCache, so a rerun
% after an interruption repeats only the level that was running.
% OUTPUT. Printed, plus Results/V1/diagnostics/truncationLadder_<stamp>.log/.mat, and the line to paste into
% defineStudyScenarios.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();
A = a1Tools();
V = cfg.V1;
dDir = fullfile(P.results, 'V1', 'diagnostics');
if ~exist(dDir, 'dir'), mkdir(dDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS');
diary(fullfile(dDir, sprintf('truncationLadder_%s.log', stamp))); diary('on');
fprintf('\n==== V1 stage 0: truncation ladder (%s) ====\n', stamp);
fprintf('  rule: first consecutive step in M with d < %.2g (ruler: F, J, Sigma at [%d %d %d])\n\n', ...
    V.truncationConvergedD, cfg.emm.truncation);

try
sc = A.buildScenarios(cfg);
opts = struct('cacheDir', P.cache, 'verbose', false);
names = [{'basin_p0'}, V.truncationFloes(:).'];
fl = struct('name', names, 'sc', [], 'scn', [], 'out', []);
for i = 1:numel(fl)
    s = sc(strcmp({sc.name}, fl(i).name));
    if numel(s) ~= 1, error('diagnoseV1Truncation:floe', 'A1 scenario %s not found exactly once.', fl(i).name); end
    fl(i).sc = s; fl(i).scn = A.buildScn(s, cfg, opts);
    fl(i).out = evaluateScenario(fl(i).scn, cfg);                     % the ruler (cached A1 fields)
    fprintf('  %-24s production [%d %d %d]: sigma [%.4f%% %.4f%%], d_W %.4f, %s\n', fl(i).name, cfg.emm.truncation, ...
        100 * fl(i).out.sigma, fl(i).out.dW, fl(i).out.class);
end
feat = @(i, tr) A.features(fl(i).sc, fl(i).scn, fl(i).out.omega, cfg.emm.nodeSpacing, ...
    [fl(i).sc.frfOptions, {'Truncation', tr}], opts);
dist = @(i, b) A.projectBias(fl(i).out, b);
lab = @(t) sprintf('[%d %d %d]', t);

%% ---- 1  ladder at p0 ---------------------------------------------------------------------------------------------
PN = V.truncationLadderPN;
Ms = V.truncationLadderM;
fprintf('\nStep 1: ladder at p0, M = %s with [P N] = [%d %d]\n', mat2str(Ms), PN);
L = struct('tr', {}, 'f', {}, 'seconds', {}, 'd', {}, 'dtheta', {}, 'relT', {});
m = size(fl(1).out.blocksAll, 1) / 2;
relT = @(fa, fb) max(max(abs(unstackTransmissibility(fa, m) - unstackTransmissibility(fb, m)) ./ ...
    abs(unstackTransmissibility(fb, m))));
chosen = [];
for k = 1:numel(Ms)
    tr = [Ms(k) PN];
    tc = tic;
    f = feat(1, tr);
    L(k).tr = tr; L(k).f = f; L(k).seconds = toc(tc);
    if k == 1
        L(k).d = NaN; L(k).dtheta = [NaN; NaN]; L(k).relT = NaN;
        fprintf('  %-14s %7.0f s\n', lab(tr), L(k).seconds);
    else
        q = dist(1, f - L(k - 1).f);
        L(k).d = q.d; L(k).dtheta = q.dtheta; L(k).relT = relT(f(:), L(k - 1).f(:));
        fprintf('  %-14s %7.0f s   step %s -> %s: d = %.4f  (dlnbeta %+.4f%%, dlnR %+.5f%%, max |dT|/|T| %.1e)%s\n', ...
            lab(tr), L(k).seconds, lab(L(k - 1).tr), lab(tr), q.d, 100 * q.dtheta, L(k).relT, ...
            ternaryStr(q.d < V.truncationConvergedD, '  < threshold', ''));
        if q.d < V.truncationConvergedD, chosen = k; break; end
    end
    if k < numel(Ms) && L(k).seconds > 60
        fprintf('      next level (M = %d) projected about %.0f s if the cost scales with M^2\n', Ms(k + 1), ...
            L(k).seconds * (Ms(k + 1) / Ms(k))^2);
    end
    save(fullfile(dDir, sprintf('truncationLadder_%s.mat', stamp)), 'L', '-v7');   % partial result, every level
    if exist('OCTAVE_VERSION', 'builtin'), fflush(stdout); else, drawnow; end
end
converged = ~isempty(chosen);
if converged
    truth = L(chosen).tr; prevTr = L(chosen - 1).tr;
    fprintf('  SELECTED: %s (step from %s, d = %.4f < %.2g)\n', lab(truth), lab(prevTr), L(chosen).d, ...
        V.truncationConvergedD);
else
    chosen = numel(L); truth = L(end).tr; prevTr = L(end - 1).tr;
    fprintf('  NOT CONVERGED within M <= %d: last step d = %.4f. %s is the best available truth, NOT converged;\n', ...
        Ms(end), L(end).d, lab(truth));
    fprintf('  report V1a as a lower bound on the truncation bias, or extend cfg.V1.truncationLadderM (logged).\n');
end
fT = cell(1, numel(fl)); fT{1} = L(chosen).f; fPrev = cell(1, numel(fl)); fPrev{1} = L(chosen - 1).f;

%% ---- 2  the same last step at the large-kR floes (diagnostic) ------------------------------------------------------
fprintf('\nStep 2 (diagnostic): the last step %s -> %s at the large-kR A1 floes\n', lab(prevTr), lab(truth));
dLast = NaN(1, numel(fl)); dLast(1) = L(chosen).d;
for i = 2:numel(fl)
    tc = tic;
    fPrev{i} = feat(i, prevTr); fT{i} = feat(i, truth);
    q = dist(i, fT{i} - fPrev{i}); dLast(i) = q.d;
    fprintf('  %-24s d = %.4f%s  (%.0f s)\n', fl(i).name, q.d, ternaryStr(q.d < V.truncationConvergedD, '', ...
        '  >= threshold: the p0-selected truth is NOT converged here'), toc(tc));
end

%% ---- 3  the P index at the selected M (diagnostic) -------------------------------------------------------------------
trP = [truth(1), V.truncationCheckP, truth(3)];
fprintf('\nStep 3 (diagnostic): P index at p0, %s -> %s\n', lab(truth), lab(trP));
tc = tic;
fP = feat(1, trP);
qP = dist(1, fP - fT{1});
fprintf('  d = %.4f%s  (%.0f s)\n', qP.d, ternaryStr(qP.d < V.truncationConvergedD, '', ...
    '  >= threshold: P is not converged at the selected M; say so in V1'), toc(tc));

%% ---- 4  preview of V1a: production model against the selected truth ------------------------------------------------
fprintf('\nStep 4 (preview of V1a, exact features): [%d %d %d] against %s\n', cfg.emm.truncation, lab(truth));
pre = struct('name', names, 'd', NaN, 'dtheta', [], 'betaPct', NaN, 'RPct', NaN, 'dTot', NaN);
for i = 1:numel(fl)
    b = fT{i} - fl(i).out.f0;
    q = dist(i, b); qt = dist(i, fl(i).out.bW + b);
    pre(i).d = q.d; pre(i).dtheta = q.dtheta; pre(i).dTot = qt.d;
    pre(i).betaPct = 100 * (exp(q.dtheta(1)) - 1); pre(i).RPct = 100 * (exp(q.dtheta(2)) - 1);
    fprintf('  %-24s d_sys %.3f, dbeta/beta %+.3f%%, dR/R %+.4f%%; with the Welch bias d_tot %.3f (d_W alone %.3f)\n', ...
        names{i}, q.d, pre(i).betaPct, pre(i).RPct, qt.d, fl(i).out.dW);
end

%% ---- verdict --------------------------------------------------------------------------------------------------------
fprintf('\nVERDICT\n');
if converged
    fprintf('  Paste into defineStudyScenarios (replacing the TO REVISE line), log it in cfg.meta.changes, commit, then V1:\n');
else
    fprintf('  NOT converged: if you accept the best available truth, paste (and log it as not converged):\n');
end
fprintf('    cfg.V1.truncationTruth = [%d %d %d];\n', truth);
fprintf('  large-kR floes, last step d: %s; P check d = %.4f\n', mat2str(dLast(2:end), 3), qP.d);

ladder = struct('levels', L, 'chosen', chosen, 'converged', converged, 'truth', truth, 'previous', prevTr, ...
    'floes', {names}, 'dLastStep', dLast, 'pCheck', struct('to', trP, 'd', qP.d, 'dtheta', qP.dtheta), ...
    'preview', pre, 'rule', V.truncationConvergedD, 'production', cfg.emm.truncation); %#ok<NASGU>
for k = 1:numel(ladder.levels), ladder.levels(k).f = []; end                          % keep the file small
save(fullfile(dDir, sprintf('truncationLadder_%s.mat', stamp)), 'ladder', '-v7');
fprintf('\n  saved %s\n', fullfile(dDir, sprintf('truncationLadder_%s.mat', stamp)));
diary('off');
catch err
    diary('off');
    rethrow(err);
end

function s = ternaryStr(c, a, b)
if c, s = a; else, s = b; end
end