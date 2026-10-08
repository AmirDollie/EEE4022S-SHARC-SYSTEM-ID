%% validateAlphaMigration.m
% Checks after applying the incident-amplitude fix and migrating twinCache.
%   1  patched computeSensorFRF is on the path
%   2  NODES: migrated cached FRF equals a fresh patched EMM solve at 3 node frequencies (exact)
%   3  MIDPOINTS: migrated, interpolated FRF vs a fresh patched solve at 3 midpoints
%      (interpolation error, same order as the original interpCheck)
%   4  CONSISTENCY with checkAlphaScaling: migrated Hacc vs alpha(w) .* old Hacc at the analysis
%      bins (the two differ only because interpolation and the alpha factor do not commute)
%   5  PRODUCTION PATH: evaluateScenario through the normal cached chain (no FRF handle) at p0,
%      Level 2C, L 2048, lsm6dsv16x, compared with the 'new' row of checkAlphaScaling
% About 6 EMM solves. Prints PASS / CHECK lines; send the output.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();
p0 = cfg.p0.vec;
sens = cfg.layout.level2C.frac .* [p0(3) 1];
g = cfg.physics.g; depth = cfg.physics.Hbasin;
alphaOf = @(w) depth * (w(:).') .^ 2 / g;

% ---- 1 ----
ws = warning('off', 'computeSensorFRF:outsideValidatedAlpha');
[~, gi] = computeSensorFRF(6.0, [0 0], p0, 'Truncation', [14 4 3]);
warning(ws);
assert(isfield(gi, 'HdeflectionRaw'), '1: computeSensorFRF is NOT the patched version.');
fprintf('PASS 1: patched computeSensorFRF on the path\n');

% ---- load the migrated p0 Level 2C twin (same key emmField uses at the nominal point) ----
spec = struct('p', p0, 'sensors', sens, 'nodeSpacing', cfg.emm.nodeSpacing, 'frfOptions', {{}}, ...
    'band', [], 'checkPoints', 3);
tw = synthesiseTwinRecords(spec, 'CacheDir', P.cache);
assert(strcmp(tw.cacheStatus, 'loaded'), ['p0 Level 2C twin was not loaded from the migrated cache ' ...
    '(status %s). Was it in the old cache? If it was just solved, checks 2-4 compare fresh with fresh.'], tw.cacheStatus);
S = load(tw.cacheFile, 'alphaMigration');
assert(isfield(S, 'alphaMigration'), 'loaded cache file has no alphaMigration tag');

% ---- 2 nodes ----
wN = tw.omegaNodes;
iN = round(linspace(2, numel(wN) - 1, 3));
Hfresh = computeSensorFRF(wN(iN), sens, p0);
e2 = max(abs(tw.HaNodes(:, iN) - Hfresh), [], 1) ./ max(abs(Hfresh), [], 1);
fprintf('%s 2: migrated nodes vs fresh patched solve, max rel %.2e\n', passStr(all(e2 < 1e-9)), max(e2));

% ---- 3 midpoints ----
wM = 0.5 * (wN(iN) + wN(iN + 1));
HfreshM = computeSensorFRF(wM, sens, p0);
e3 = max(abs(tw.Hacc(wM) - HfreshM), [], 1) ./ max(abs(HfreshM), [], 1);
fprintf('%s 3: migrated interpolation vs fresh patched solve at midpoints, max rel %.2e (old interpCheck H %.2e)\n', ...
    passStr(all(e3 < 1e-3)), max(e3), tw.interpCheck.maxRelErrH);

% ---- 4 consistency with checkAlphaScaling's handle ----
[omega, ~] = welchBins(cfg.welch.L, cfg.acq.dt, cfg.band.estimate, cfg.welch.binStep, cfg.band.basin);
old = load(tw.cacheFile);                                   % node values are migrated; rebuild old interpolant
HdOld = (old.HaNodes ./ alphaOf(wN)) ./ (-wN .^ 2);
HoldScaled = alphaOf(omega) .* (-(omega .^ 2)) .* interpComplex(wN, HdOld, omega);
e4 = max(abs(tw.Hacc(omega) - HoldScaled), [], 1) ./ max(abs(HoldScaled), [], 1);
fprintf('%s 4: migrated vs alpha .* old interpolant at %d analysis bins, max rel %.2e\n', ...
    passStr(max(e4) < 1e-3), numel(omega), max(e4));

% ---- 5 production path ----
out = evaluateScenario(struct('Snn', cfg.noise.lsm6dsv16x, 'cacheDir', P.cache), cfg);
fprintf('\nProduction path, p0, Level 2C, L %d, lsm6dsv16x:\n', cfg.welch.L);
fprintf('  sigma lnB %.4f %%  sigma lnR %.4f %%  dW %.3f  dTot %.3f  minCoh %.3f  class %s\n', ...
    100 * out.sigma(1), 100 * out.sigma(2), out.dW, out.dTot, out.minCoherence, out.class);
d = dir(fullfile(P.results, 'E1', 'alphaScalingCheck_*.mat'));
if ~isempty(d)
    [~, j] = max([d.datenum]); R = load(fullfile(d(j).folder, d(j).name));
    k = find(strcmp({R.res.noise}, 'lsm6dsv16x'), 1);
    ref = R.res(k).new;
    rs = out.sigma(:) ./ ref.sigma(:);
    fprintf('  checkAlphaScaling NEW: sigma lnB %.4f %%  sigma lnR %.4f %%  dW %.3f  class %s\n', ...
        100 * ref.sigma(1), 100 * ref.sigma(2), ref.dW, ref.class);
    fprintf('%s 5: production / handle sigma ratios %.4f %.4f, dW difference %.4f, same class: %d\n', ...
        passStr(all(abs(rs - 1) < 0.01) && abs(out.dW - ref.dW) < 0.01 && strcmp(out.class, ref.class)), ...
        rs(1), rs(2), out.dW - ref.dW, strcmp(out.class, ref.class));
else
    fprintf('CHECK 5: no alphaScalingCheck_*.mat found in Results/E1 to compare against\n');
end

function s = passStr(ok)
if ok, s = 'PASS'; else, s = 'CHECK'; end
end

function Z = interpComplex(x, Y, xq)
% Same scheme as synthesiseTwinRecords/interpAcc (spline on real and imaginary parts).
xq = min(max(xq(:).', x(1)), x(end));
Z = complex(interp1(x(:), real(Y).', xq(:), 'spline'), interp1(x(:), imag(Y).', xq(:), 'spline')).';
end