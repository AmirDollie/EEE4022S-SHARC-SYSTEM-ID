%% runEMMTwinInversion.m
% First nonlinear inversion: the NOISE-FREE EMM twin test.
%
%   p_true --(EMM, predicted anchors)--> f_obs --(solveInverse)--> p_hat
%
% Deliberately an "inverse crime": the SAME forward map (same EMM
% truncation, same predicted-anchor feature extractor, same prediction
% model) generates f_obs and is inverted. It is an integration test of the
% whole inversion chain on the real physics, answering two questions:
%   1. Can the nonlinear inverse recover all three perturbed parameters?
%   2. Do different starting points converge to the same solution?
% No noise, no feature subsets, gamma not fixed: those come after this works.
%
% The truth is generated ONCE, outside the optimiser; solveInverse only ever
% sees f_obs, never p_true. Recovery errors are computed here, afterwards.
%
% COST: each forward call is ~35-45 s. One lsqnonlin iteration with central
% differences costs ~7 calls (~4-5 min); an exact twin typically needs
% ~5-10 iterations, so ~30-50 min per start, capped by MAX_FEVALS below.
% Results are saved after EVERY start, so stopping early keeps finished
% starts. Use START_INDICES to run a subset (e.g. 1 first, the rest later).
%
% OUTPUT (in Parameter Inversion/Results/)
%   EMMTwin_noiseFree_<timestamp>.mat : twin definition + every solveInverse
%                                       result + recovery summary
%   EMMTwin_noiseFree_<timestamp>.txt : console log of the whole run
%
% Lives in Inverse Mapping/Parameter Inversion/EMM Twin/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));         % solveInverse.m, inverseResidual.m (Parameter Inversion)
addpath(fullfile(thisDir, '..', '..'));   % modalFeatureVector.m, predictFeatureAnchors.m (Inverse Mapping)

clearvars -except thisDir
close all, clc

%% ======================= EXPERIMENT DEFINITION ===========================
% Reference point (reference/test parameter set, NOT the physical floe).
beta0 = 4.6985e-5; gamma0 = 1.4548e-3; R0 = 0.3830;
pRef = [beta0; gamma0; R0];
paramNames = {'beta', 'gamma', 'R'};

% Hidden truth, scaled units x = p./pRef. All three genuinely perturbed.
xTrue = [1.08; 0.93; 1.04];

% Optimiser bounds: conservative, well inside the predicted-anchor box
% verified by runProbeValidRegionEMM (probeValidRegion_pred_20260923_185905:
% all 8 corners of [0.84 0.84 0.92]..[1.16 1.16 1.08] valid).
lb = [0.90; 0.85; 0.95];
ub = [1.10; 1.15; 1.05];
VERIFIED_LB = [0.84; 0.84; 0.92];
VERIFIED_UB = [1.16; 1.16; 1.08];

% Starts (columns): the reference point plus three mixed perturbations.
xStarts = [1.00  0.95  1.05  1.08
           1.00  1.08  0.90  1.02
           1.00  1.03  0.97  0.96];
START_INDICES = 2:4;     % reference-point start first; then 2:4

% Minimum distance from any bound, as a fraction of that coordinate's box
% width, for the truth and for every start.
MIN_CLEARANCE = 0.05;

% Solver budget. The solveInverse default (1000 evaluations) would allow
% ~10 h per start at ~38 s per call; 150 caps a start at ~1.5-2 h.
MAX_FEVALS = 150;
MAX_ITERS = 30;

modelFile = fullfile(thisDir, '..', '..', 'anchorPredictionModel.mat');
EXPECTED_ERROR_IDS = {'modalFeatureVector:invalidFeatureFit'};
%% =========================================================================

%% Checks before anything expensive runs
assert(isfile(modelFile), 'runEMMTwinInversion:missingModel', ...
    'Prediction model not found: %s', modelFile);
model = load(modelFile);
assert(max(abs(model.p0(:) ./ pRef - 1)) < 1e-12, 'runEMMTwinInversion:modelRefMismatch', ...
    'anchorPredictionModel.p0 = %s does not match pRef = %s.', mat2str(model.p0(:).', 6), mat2str(pRef.', 6));
assert(all(lb >= VERIFIED_LB) && all(ub <= VERIFIED_UB), 'runEMMTwinInversion:boundsOutsideVerified', ...
    'Bounds [%s]..[%s] extend outside the verified box [%s]..[%s].', ...
    num2str(lb.'), num2str(ub.'), num2str(VERIFIED_LB.'), num2str(VERIFIED_UB.'));
width = ub - lb;
clearance = @(x) min([x - lb, ub - x], [], 2) ./ width;   % per coordinate, fraction of width
assert(all(clearance(xTrue) >= MIN_CLEARANCE - 1e-12), 'runEMMTwinInversion:truthNearBound', ...
    'xTrue = %s is within %.0f%% of a bound (clearance %s).', mat2str(xTrue.', 4), ...
    100 * MIN_CLEARANCE, mat2str(clearance(xTrue).', 3));
for k = 1:size(xStarts, 2)
    assert(all(clearance(xStarts(:, k)) >= MIN_CLEARANCE - 1e-12), 'runEMMTwinInversion:startNearBound', ...
        'start %d = %s is within %.0f%% of a bound.', k, mat2str(xStarts(:, k).', 4), 100 * MIN_CLEARANCE);
end
assert(all(ismember(START_INDICES, 1:size(xStarts, 2))), 'runEMMTwinInversion:badStartIndices', ...
    'START_INDICES must index columns of xStarts.');

%% Output files and log
resultsDir = fullfile(thisDir, '..', 'Results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['EMMTwin_noiseFree_', stamp]);
diary([baseName, '.txt']);
try   % closes the log even if the run errors (onCleanup does not fire at the end of a script)

fwd = @(p) modalFeatureVector(p(1), p(2), p(3), 'Strict', true, ...
    'AnchorMode', 'predicted', 'AnchorModel', modelFile);

fprintf('runEMMTwinInversion  (%s)   noise-free EMM twin, predicted anchors\n', stamp);
fprintf('model: %s  (built at eps = %g)\n', modelFile, model.epsilon);
fprintf('pRef   = [%.6g %.6g %.6g]\n', pRef);
fprintf('bounds (scaled): lb = %s   ub = %s\n', mat2str(lb.', 4), mat2str(ub.', 4));
fprintf('budget per start: MaxFunctionEvaluations = %d, MaxIterations = %d\n\n', MAX_FEVALS, MAX_ITERS);

%% Hidden truth -> synthetic observation (once)
pTrue = xTrue .* pRef;
t0 = tic;
[fObs, truthDetails] = fwd(pTrue);
fprintf('TRUTH (hidden from the solver)\n');
fprintf('  xTrue = %s\n  pTrue = %s\n', mat2str(xTrue.', 4), mat2str(pTrue.', 6));
fprintf('  fObs  = %s   (%.1fs)\n', mat2str(fObs(:).', 8), toc(t0));
truthMargins = arrayfun(@(d) min(d.omegaStar - min(d.omega3), max(d.omega3) - d.omegaStar) / ...
    (max(d.omega3) - min(d.omega3)), truthDetails);
fprintf('  anchor-bracket margins at the truth: %s\n\n', mat2str(truthMargins(:).', 3));

twin = struct('stamp', stamp, 'pRef', pRef, 'paramNames', {paramNames}, 'xTrue', xTrue, ...
    'pTrue', pTrue, 'fObs', fObs, 'truthDetails', {truthDetails}, 'lb', lb, 'ub', ub, ...
    'verifiedLB', VERIFIED_LB, 'verifiedUB', VERIFIED_UB, 'xStarts', xStarts, ...
    'startIndices', START_INDICES, 'maxFevals', MAX_FEVALS, 'maxIters', MAX_ITERS, ...
    'expectedErrorIDs', {EXPECTED_ERROR_IDS}, 'modelFile', modelFile, 'model', model);
results = cell(1, size(xStarts, 2));
summary = struct('start', {}, 'x0', {}, 'status', {}, 'exitflag', {}, 'xHat', {}, 'pHat', {}, ...
    'relErr', {}, 'maxAbsRelErr', {}, 'resnorm', {}, 'nIterations', {}, 'funcCount', {}, ...
    'runtimeMin', {}, 'leftValidRegion', {});

%% Solve from each start
for k = START_INDICES
    x0 = xStarts(:, k);
    fprintf('==== START %d of %d: x0 = %s ====\n', k, size(xStarts, 2), mat2str(x0.', 4));
    res = solveInverse(fwd, fObs, pRef, x0, lb, ub, ...
        'ExpectedErrorIDs', EXPECTED_ERROR_IDS, ...
        'MaxFunctionEvaluations', MAX_FEVALS, ...
        'MaxIterations', MAX_ITERS, ...
        'Display', 'iter');
    % solverOptions holds a handle to solveInverse's nested OutputFcn, which
    % would drag its whole workspace into the .mat; res.settings already
    % records every option as a plain struct.
    results{k} = rmfield(res, 'solverOptions');

    relErr = res.pHat ./ pTrue - 1;                         % = xHat./xTrue - 1
    s = struct('start', k, 'x0', x0, 'status', res.status, 'exitflag', res.exitflag, ...
        'xHat', res.xHat, 'pHat', res.pHat, 'relErr', relErr, 'maxAbsRelErr', max(abs(relErr)), ...
        'resnorm', res.resnorm, 'nIterations', res.nIterations, 'funcCount', res.funcCount, ...
        'runtimeMin', res.runtimeSec / 60, 'leftValidRegion', strcmp(res.status, 'leftValidRegion'));
    summary(end+1) = s; %#ok<AGROW>

    fprintf('\n  status %s (exitflag %g), %d iterations, %d forward calls, %.1f min\n', ...
        res.status, res.exitflag, res.nIterations, res.funcCount, res.runtimeSec / 60);
    if s.leftValidRegion
        fprintf('  LEFT VALID REGION at x = %s: %s\n', mat2str(res.failedX.', 6), res.error.message);
    else
        fprintf('  xHat   = %s   (xTrue %s)\n', mat2str(res.xHat.', 8), mat2str(xTrue.', 4));
        fprintf('  relErr = [%s]   (beta gamma R)\n', sprintf(' %+.2e', relErr));
        fprintf('  resnorm = %.3e\n', res.resnorm);
        h = res.history;
        fprintf('  iterate history (scaled):\n    %4s  %10s %10s %10s   %10s\n', 'it', 'beta', 'gamma', 'R', 'resnorm');
        for i = 1:numel(h.iteration)
            fprintf('    %4d  %10.6f %10.6f %10.6f   %10.3e\n', h.iteration(i), h.x(i, :), h.resnorm(i));
        end
    end
    fprintf('\n');

    save([baseName, '.mat'], 'twin', 'results', 'summary');   % after every start
end

%% Summary across starts
fprintf('==== SUMMARY ====\n');
fprintf('  xTrue = %s\n', mat2str(xTrue.', 4));
fprintf('  %5s  %-16s %-36s %-36s %10s %5s %6s %7s\n', 'start', 'status', 'xHat (beta gamma R)', ...
    'relErr (beta gamma R)', 'resnorm', 'iters', 'calls', 'min');
for s = summary
    fprintf('  %5d  %-16s %-36s %-36s %10.2e %5d %6d %7.1f\n', s.start, s.status, ...
        sprintf('%.6f ', s.xHat), sprintf('%+.1e ', s.relErr), s.resnorm, s.nIterations, s.funcCount, s.runtimeMin);
end

% Only starts with status 'solved' are aggregated: 'budgetExhausted',
% 'solverStopped' and 'leftValidRegion' end wherever they stopped, which is
% not a solution (they still appear, with their status, in the table above).
solved = strcmp({summary.status}, 'solved');
if nnz(solved) >= 2
    X = [summary(solved).xHat];
    spread = max(X, [], 2) - min(X, [], 2);
    fprintf('\n  spread of xHat across %d solved starts: [%s]   (beta gamma R)\n', ...
        nnz(solved), sprintf(' %.2e', spread));
end

% Local identifiability at the solution (solved starts only), from the
% Jacobian lsqnonlin returned (no extra forward calls). res.jacobian = dr/dx =
% diag(1./fObs)*J*diag(pRef); multiplying columns by xHat gives the scaled
% Jacobian at p_hat, comparable to (but not expected to equal) the
% reference-point J_s (kappa ~ 293, weak vector ~ gamma).
for k = find(~cellfun(@isempty, results))
    res = results{k};
    if ~strcmp(res.status, 'solved'), continue; end
    JsHat = res.jacobian * diag(res.xHat);
    [~, S, V] = svd(JsHat, 'econ');
    sv = diag(S);
    fprintf('  start %d: scaled Jacobian at p_hat: sigma = [%s], kappa = %.0f, weak vector = [%s]\n', ...
        k, sprintf(' %.4g', sv), sv(1) / sv(end), sprintf(' %+.4f', V(:, end)));
end

save([baseName, '.mat'], 'twin', 'results', 'summary');
fprintf('\nSaved: %s.mat\n', baseName);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');