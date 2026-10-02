%% diagnoseTruncationP0.m
% DIAGNOSTIC (not a study): EMM truncation convergence of the transmissibility FEATURES at p0.
%
% WHY. testA1FloePropertySweep check 10 failed: [50 10 10] vs [70 15 15] at p0 gives d = 0.57 (threshold
% 0.1), with only 0.2% change in T. [50 10 10] is the frozen truncation of every study so far, so
% before A1 is changed or launched we need to know: is [70 15 15] itself converged, which index drives
% the difference, where is it in frequency, and do the A1 METRICS (sigma, d_W) move, or only the bias?
%
% WHAT IT COMPUTES (Level 2C, reference s26, LSM6DSV16X, setting F, node grid 0.2 rad/s)
%   1. Ladder of nominal features: [50 10 10] (production), [70 10 10] (M only), [50 15 10] (P only),
%      [50 10 15] (N only), [70 15 15], [90 20 20]. Every difference is projected with the SAME ruler,
%      F at [50 10 10]: d = sqrt(dth' F dth), dth = F^-1 J' Sigma^-1 b. Key number: d(70 -> 90).
%   2. Where: per accepted bin, the standardised feature error s_k = sqrt(b_k' Sigma_k^-1 b_k) for
%      b = f70 - f50; top bins with omega, the largest |dT|/|T| and which transmissibility, and
%      |H_s26| relative to its in-band maximum (a reference dip would amplify the ratio error).
%   3. Metrics: a full evaluateScenario at [70 15 15] (Jacobian included) against [50 10 10]:
%      sigma, d_W, kappa. If they agree to well under 1%, A1's classes do not depend on the truncation
%      and the failure is purely a MODEL-DISCREPANCY BIAS.
%   4. What a real inversion with the [50 10 10] model would carry if the truth were [70 15 15]:
%      d_tot from the Welch bias plus the truncation bias, added as vectors.
%
% COST. New nominal fields at [70 10 10], [50 15 10], [50 10 15], [90 20 20] (29 nodes each) and the four
% derivative fields at [70 15 15]: roughly 10 to 20 min the first time; everything is cached.
% Output: printed, and Results/A1/diagnostics/truncationP0_<stamp>.log/.mat
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();
A = a1Tools();
dDir = fullfile(P.results, 'A1', 'diagnostics');
if ~exist(dDir, 'dir'), mkdir(dDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS');
diary(fullfile(dDir, sprintf('truncationP0_%s.log', stamp))); diary('on');
fprintf('\n==== Truncation diagnostic at p0 (%s) ====\n\n', stamp);

try
sc = A.buildScenarios(cfg);
b0 = sc(strcmp({sc.name}, 'basin_p0'));
opts = struct('cacheDir', P.cache, 'verbose', true);
scn = A.buildScn(b0, cfg, opts);
out50 = evaluateScenario(scn, cfg);
w = out50.omega; v = out50.validMask; pB = size(out50.blocksAll, 1); m = pB / 2;
fprintf('Production [50 10 10]: sigma = [%.4f%% %.4f%%], d_W = %.4f, kappa = %.3g, %d of %d bins accepted\n\n', ...
    100 * out50.sigma, out50.dW, out50.condF, nnz(v), numel(v));

%% ---- 1  the ladder ----------------------------------------------------------------------------------------
tr = {[50 10 10], [70 10 10], [50 15 10], [50 10 15], [70 15 15], [90 20 20]};
lab = cellfun(@(t) sprintf('[%d %d %d]', t), tr, 'UniformOutput', false);
f = cell(size(tr)); f{1} = out50.f0;
for i = 2:numel(tr)
    fprintf('features at %s\n', lab{i});
    tc = tic;
    f{i} = A.features(b0, scn, w, cfg.emm.nodeSpacing, {'Truncation', tr{i}}, opts);
    fprintf('  %.0f s\n', toc(tc));
end
relT = @(fa, fb) max(max(abs(unstackTransmissibility(fa, m) - unstackTransmissibility(fb, m)) ./ ...
    abs(unstackTransmissibility(fb, m))));
pairs = [1 2; 1 3; 1 4; 1 5; 1 6; 5 6];
pairNote = {'M only', 'P only', 'N only', 'the A1 check', 'production vs best', 'IS [70 15 15] CONVERGED?'};
fprintf('\nLadder (common ruler: F at [50 10 10]; threshold %.2g)\n', cfg.A1.maxTruncationBiasDistance);
fprintf('  %-28s %10s %14s   %s\n', 'difference', 'd', 'max |dT|/|T|', '');
D = struct('from', {}, 'to', {}, 'd', {}, 'relT', {}, 'dtheta', {});
for k = 1:size(pairs, 1)
    a = pairs(k, 1); b = pairs(k, 2);
    q = A.projectBias(out50, f{b} - f{a});
    rt = relT(f{a}(:), f{b}(:));
    D(k) = struct('from', lab{a}, 'to', lab{b}, 'd', q.d, 'relT', rt, 'dtheta', q.dtheta);
    fprintf('  %-28s %10.4f %14.2e   %s (dlnbeta %+.3f%%, dlnR %+.4f%%)\n', [lab{a} ' -> ' lab{b}], q.d, rt, ...
        pairNote{k}, 100 * q.dtheta);
end

%% ---- 2  where in frequency ------------------------------------------------------------------------------------
b = f{5} - f{1};
idx = find(v);
sK = zeros(size(idx)); worstJ = zeros(size(idx)); relK = zeros(size(idx));
T50 = unstackTransmissibility(f{1}, m); T70 = unstackTransmissibility(f{5}, m);
for n = 1:numel(idx)
    k = idx(n); rows = (k - 1) * pB + (1:pB);
    L = chol(out50.blocksAll(:, :, k), 'lower');
    sK(n) = norm(L \ b(rows));
    [relK(n), worstJ(n)] = max(abs(T70(:, k) - T50(:, k)) ./ abs(T70(:, k)));
end
Hall = out50.field.H(w);
Href = abs(Hall(scn.layout(scn.ref), :)); HrefRel = Href / max(Href);
names = {'s2/s26', 's6/s26', 's14/s26'};
[~, ord] = sort(sK, 'descend');
share = cumsum(sK(ord).^2) / sum(sK.^2);
fprintf('\nWhere ([70 15 15] - [50 10 10]): standardised error per accepted bin, s_k = sqrt(b_k'' Sigma_k^-1 b_k)\n');
fprintf('  total sum s_k^2 = %.3f over %d bins; top 5 bins carry %.0f%%, top 10 %.0f%%\n', sum(sK.^2), numel(sK), ...
    100 * share(min(5, end)), 100 * share(min(10, end)));
fprintf('  %8s %8s %10s %10s %10s\n', 'omega', 's_k', 'max dT/T', 'worst T', '|H_s26|/max');
for n = reshape(ord(1:min(10, end)), 1, [])
    k = idx(n);
    fprintf('  %8.3f %8.3f %10.2e %10s %10.3f\n', w(k), sK(n), relK(n), names{worstJ(n)}, HrefRel(k));
end
[~, kMin] = min(HrefRel);
fprintf('  |H_s26| minimum in the band: %.3f of max at omega = %.3f rad/s\n', HrefRel(kMin), w(kMin));
lowBand = w(idx) < 5; fprintf('  share of sum s_k^2 below 5 rad/s: %.0f%%; above 7 rad/s: %.0f%%\n', ...
    100 * sum(sK(lowBand).^2) / sum(sK.^2), 100 * sum(sK(w(idx) > 7).^2) / sum(sK.^2));

%% ---- 3  do the A1 metrics move? -------------------------------------------------------------------------------
fprintf('\nFull evaluateScenario at [70 15 15] (Jacobian included)\n');
scn70 = scn; scn70.frfOptions = {'Truncation', [70 15 15]};
out70 = evaluateScenario(scn70, cfg);
relS = out70.sigma(:).' ./ out50.sigma(:).' - 1;
fprintf('  sigma [%.4f%% %.4f%%] vs [%.4f%% %.4f%%]: relative change [%+.2e %+.2e]\n', 100 * out70.sigma, ...
    100 * out50.sigma, relS);
fprintf('  d_W %.4f vs %.4f (%+.2e); kappa %.4g vs %.4g; class %s vs %s\n', out70.dW, out50.dW, ...
    out70.dW / out50.dW - 1, out70.condF, out50.condF, out70.class, out50.class);

%% ---- 4  what a real inversion with [50 10 10] would carry ------------------------------------------------------------
qTot = A.projectBias(out50, out50.bW + (f{5} - f{1}));
qTot90 = A.projectBias(out50, out50.bW + (f{6} - f{1}));
fprintf('\nIf the truth is [70 15 15] ([90 20 20]) and the inverse model is [50 10 10]:\n');
fprintf('  d_W alone %.3f; truncation alone %.3f (%.3f); Welch + truncation as vectors d_tot = %.3f (%.3f)\n', ...
    out50.dW, D(4).d, D(5).d, qTot.d, qTot90.d);
fprintf('  (class threshold d_tot <= 0.5: p0 would be %s with the [50 10 10] model)\n', ...
    ternaryStr(qTot.d <= 0.5, 'still Identifiable', 'NO LONGER Identifiable'));

%% ---- verdict -------------------------------------------------------------------------------------------------------
fprintf('\nVERDICT\n');
if D(6).d < cfg.A1.maxTruncationBiasDistance
    fprintf('  [70 15 15] -> [90 20 20]: d = %.4f < %.2g: [70 15 15] is converged for these features; [50 10 10] is not.\n', ...
        D(6).d, cfg.A1.maxTruncationBiasDistance);
else
    fprintf('  [70 15 15] -> [90 20 20]: d = %.4f >= %.2g: [70 15 15] is NOT converged either; extend the ladder.\n', ...
        D(6).d, cfg.A1.maxTruncationBiasDistance);
end
fprintf('  A1 metrics at [70 15 15] vs [50 10 10]: max |relative change in sigma| = %.2e, d_W %.2e.\n', ...
    max(abs(relS)), abs(out70.dW / out50.dW - 1));

trunc = struct('ladder', D, 'truncations', {tr}, 'omega', w, 'validMask', v, 'sK', sK, 'binIdx', idx, ...
    'HrefRel', HrefRel, 'sigma50', out50.sigma, 'sigma70', out70.sigma, 'dW50', out50.dW, 'dW70', out70.dW, ...
    'dTot50vs70', qTot.d, 'dTot50vs90', qTot90.d); %#ok<NASGU>
save(fullfile(dDir, sprintf('truncationP0_%s.mat', stamp)), 'trunc', '-v7');
fprintf('\n  saved %s\n', fullfile(dDir, sprintf('truncationP0_%s.mat', stamp)));
diary('off');
catch err
    diary('off');
    rethrow(err);
end

function s = ternaryStr(c, a, b)
if c, s = a; else, s = b; end
end