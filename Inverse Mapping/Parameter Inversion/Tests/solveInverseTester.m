%% solveInverseTester.m
% Validation tests for solveInverse.m, using cheap analytic toy forward
% maps (never the EMM), so the file runs in a few seconds. solveInverse.m
% is forward-model-agnostic, so its own logic (option wiring, bounds,
% error classification, logging, result assembly) is fully testable
% without the physics. The EMM only enters at runEMMTwinInversion.m.
%
% Requires the Optimization Toolbox (lsqnonlin, optimoptions).
%
% TEST LIST
%   1. Input validation: bad handle / pRef / x0 / bounds / fObs / options,
%      and an underdetermined FeatureMask, all rejected before any solve.
%   2. Exact toy twin recovery: status 'solved', parameters recovered,
%      pHat = xHat.*pRef, residual ~0, history logged and consistent,
%      returned Jacobian matches the analytic dr/dx at the solution.
%   3. Finite-difference settings actually reach lsqnonlin: the points the
%      solver evaluates around x0 show central differences at step ~1e-4
%      by default (NOT lsqnonlin's ~1.5e-8 default), and forward
%      differences at 1e-3 when requested.
%   4. Error classification:
%      (4a) an expected error mid-solve -> status 'leftValidRegion',
%           recorded with the failing point, not thrown;
%      (4b) the same error WITHOUT ExpectedErrorIDs -> re-thrown;
%      (4c) an UNEXPECTED identifier -> re-thrown even when other IDs are
%           expected (a bug is never recorded as a finding);
%      (4d) an expected error at the very first evaluation (x0 itself
%           invalid) -> still recorded (lsqnonlin may wrap first-evaluation
%           failures; the cause chain is searched);
%      (4e) ExpectedErrorIDs given as a single char string also works.
%   5. FeatureMask is forwarded: a NaN feature masked out still solves,
%      with residual/Jacobian sized to the selected features.
%   6. Budget exhaustion: MaxIterations = 1 -> status 'budgetExhausted'.
%
% Lives in Inverse Mapping/Parameter Inversion/Tests/; solveInverse.m and
% inverseResidual.m are one level up in Inverse Mapping/Parameter Inversion/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));   % solveInverse.m, inverseResidual.m

clearvars -except thisDir
clc

fprintf('================================================================\n');
fprintf(' solveInverseTester.m\n');
fprintf('================================================================\n\n');

% Shared toy problem: 3 parameters -> 5 features (see toyMap below),
% truth within +/-10%% of the reference point, like the EMM twin will be.
pRef  = [2; 3; 5];
pTrue = [2.2; 2.7; 5.5];              % x_true = [1.1; 0.9; 1.1]
xTrue = pTrue ./ pRef;
fObs  = toyMap(pTrue);
x0 = [1; 1; 1];
lb = [0.5; 0.5; 0.5];
ub = [2; 2; 2];
toyFwd = @toyMap;

%% TEST 1: input validation
fprintf('--- TEST 1: input validation -------------------------------\n');

cases = {
    'forwardFcn not a handle',  @() solveInverse('toyMap', fObs, pRef, x0, lb, ub),                   'solveInverse:badForwardFcn'
    'pRef has a zero',          @() solveInverse(toyFwd, fObs, [2; 0; 5], x0, lb, ub),               'solveInverse:badPRef'
    'x0 wrong length',          @() solveInverse(toyFwd, fObs, pRef, [1; 1], lb, ub),                'solveInverse:badX0'
    'x0 outside bounds',        @() solveInverse(toyFwd, fObs, pRef, [1; 1; 3], lb, ub),             'solveInverse:badX0'
    'lb has a zero',            @() solveInverse(toyFwd, fObs, pRef, x0, [0; 0.5; 0.5], ub),         'solveInverse:badBounds'
    'lb >= ub',                 @() solveInverse(toyFwd, fObs, pRef, x0, lb, [2; 0.5; 2]),           'solveInverse:badBounds'
    'ub has Inf',               @() solveInverse(toyFwd, fObs, pRef, x0, lb, [2; Inf; 2]),           'solveInverse:badBounds'
    'fObs has NaN',             @() solveInverse(toyFwd, [1; NaN; 1; 1; 1], pRef, x0, lb, ub),       'solveInverse:badFObs'
    'unknown option',           @() solveInverse(toyFwd, fObs, pRef, x0, lb, ub, 'Strict', true),   'solveInverse:badOption'
    'odd option count',         @() solveInverse(toyFwd, fObs, pRef, x0, lb, ub, 'Display'),        'solveInverse:badOption'
    'bad ExpectedErrorIDs',     @() solveInverse(toyFwd, fObs, pRef, x0, lb, ub, 'ExpectedErrorIDs', {1}), 'solveInverse:badOption'
    'bad FD step',              @() solveInverse(toyFwd, fObs, pRef, x0, lb, ub, 'FiniteDifferenceStepSize', 0), 'solveInverse:badOption'
    'bad FD type',              @() solveInverse(toyFwd, fObs, pRef, x0, lb, ub, 'FiniteDifferenceType', 'backward'), 'solveInverse:badOption'
    'underdetermined mask',     @() solveInverse(toyFwd, fObs, pRef, x0, lb, ub, 'FeatureMask', [1 3]), 'solveInverse:underdetermined'
};
for c = 1:size(cases, 1)
    gotId = '<no error>';
    try
        cases{c, 2}();
    catch err
        gotId = err.identifier;
    end
    assert(strcmp(gotId, cases{c, 3}), sprintf('Test1:case%d', c), ...
        '%s: expected %s, got %s', cases{c, 1}, cases{c, 3}, gotId);
end
fprintf('  (1) all %d invalid-input cases throw their documented identifier: PASS\n\n', size(cases, 1));

%% TEST 2: exact toy twin recovery
fprintf('--- TEST 2: exact toy twin recovery ------------------------\n');

res = solveInverse(toyFwd, fObs, pRef, x0, lb, ub);
relErr = abs(res.pHat - pTrue) ./ pTrue;

assert(strcmp(res.status, 'solved') && res.exitflag > 0, 'Test2:status', ...
    'status = %s, exitflag = %g; expected solved with exitflag > 0', res.status, res.exitflag);
assert(max(relErr) < 1e-8, 'Test2:recovery', 'max relative parameter error %.3g, expected < 1e-8', max(relErr));
assert(max(abs(res.pHat - res.xHat .* pRef)) < 1e-14, 'Test2:unscale', 'pHat ~= xHat.*pRef');
assert(res.resnorm < 1e-16 && abs(res.resnorm - sum(res.residual.^2)) < 1e-20, 'Test2:resnorm', ...
    'resnorm %.3g not ~0 or inconsistent with residual', res.resnorm);
assert(isempty(res.error) && isempty(res.failedX), 'Test2:noError', 'error fields set on a clean solve');
fprintf('  (2a) toy truth recovered, max relative error %.2g, resnorm %.2g: PASS\n', max(relErr), res.resnorm);

h = res.history;
assert(~isempty(h.iteration) && size(h.x, 1) == numel(h.iteration) && size(h.x, 2) == 3, 'Test2:historySize', ...
    'history is empty or inconsistently sized');
assert(all(diff(h.iteration) >= 0), 'Test2:historyOrder', 'history iterations are not in order');
assert(all(all(h.x >= lb.' & h.x <= ub.')), 'Test2:historyBounds', 'a logged iterate lies outside [lb, ub]');
assert(max(abs(h.p(end, :).' - h.x(end, :).' .* pRef)) < 1e-14, 'Test2:historyP', 'history.p ~= history.x.*pRef');
assert(norm(h.x(end, :).' - res.xHat) < 1e-6, 'Test2:historyEnd', 'last logged iterate is not the returned solution');
assert(res.funcCount > 0 && res.runtimeSec > 0 && res.nIterations >= 1, 'Test2:bookkeeping', ...
    'funcCount/runtime/nIterations not populated');
fprintf('  (2b) %d iterations logged in order, within bounds, ending at xHat; %d forward calls: PASS\n', ...
    numel(h.iteration), res.funcCount);

% Analytic dr/dx at the solution: diag(1./fObs) * df/dp * diag(pRef).
p = res.pHat;
dfdp = [1, 1, 0; p(3), 0, p(1); 0, 2 * p(2), 0; 1, 0, 1; 0, 2, 0];
Jexact = diag(1 ./ fObs) * dfdp * diag(pRef);
assert(isequal(size(res.jacobian), [5, 3]), 'Test2:jacSize', 'jacobian is %s, expected 5x3', mat2str(size(res.jacobian)));
jacErr = max(abs(res.jacobian(:) - Jexact(:))) / max(abs(Jexact(:)));
assert(jacErr < 1e-6, 'Test2:jacobian', 'returned Jacobian differs from analytic dr/dx by %.3g (relative)', jacErr);
fprintf('  (2c) returned Jacobian matches analytic dr/dx at the solution (rel. error %.2g): PASS\n\n', jacErr);

%% TEST 3: finite-difference settings reach lsqnonlin
fprintf('--- TEST 3: finite-difference settings ---------------------\n');

% (3a) Default: central differences, step 1e-4.
recFwd('reset');
res = solveInverse(@(p) recFwd('eval', p), fObs, pRef, x0, lb, ub);
X = recFwd('get') ./ pRef.';          % evaluated points, in scaled units
[plusStep, minusStep] = fdStepsAround(X, x0);
assert(all(isfinite(plusStep)) && all(isfinite(minusStep)), 'Test3a:central', ...
    'did not find both +h and -h evaluations around x0 in every coordinate (central differences not in use?)');
assert(all(plusStep > 0.5e-4 & plusStep < 2e-4) && all(minusStep > 0.5e-4 & minusStep < 2e-4), 'Test3a:step', ...
    'FD steps around x0 are %s / %s, expected ~1e-4 (lsqnonlin default would be ~1e-8)', ...
    mat2str(plusStep.', 3), mat2str(minusStep.', 3));
assert(res.settings.FiniteDifferenceStepSize == 1e-4 && strcmp(res.settings.FiniteDifferenceType, 'central'), ...
    'Test3a:settings', 'result.settings does not record the default FD settings');
assert(res.solverOptions.FiniteDifferenceStepSize == 1e-4, 'Test3a:solverOptions', ...
    'optimoptions object does not carry FiniteDifferenceStepSize = 1e-4');
fprintf('  (3a) default: central differences at step ~1e-4 around x0 (+%s / -%s): PASS\n', ...
    mat2str(plusStep.', 2), mat2str(minusStep.', 2));

% (3b) Forward differences at 1e-3 when requested.
recFwd('reset');
solveInverse(@(p) recFwd('eval', p), fObs, pRef, x0, lb, ub, ...
    'FiniteDifferenceType', 'forward', 'FiniteDifferenceStepSize', 1e-3);
X = recFwd('get') ./ pRef.';
[plusStep, minusStep] = fdStepsAround(X, x0);
assert(all(plusStep > 0.5e-3 & plusStep < 2e-3), 'Test3b:step', ...
    'forward FD steps around x0 are %s, expected ~1e-3', mat2str(plusStep.', 3));
assert(all(~isfinite(minusStep)), 'Test3b:forwardOnly', ...
    'found -h evaluations around x0 although forward differences were requested');
fprintf('  (3b) FiniteDifferenceType/StepSize options take effect (forward, ~1e-3): PASS\n\n');

%% TEST 4: error classification
fprintf('--- TEST 4: error classification ---------------------------\n');

% Toy "valid region": x(1) <= 1.3. The truth sits outside it at x(1) = 1.5,
% so the solver must step past the wall to converge.
wall = 1.3;
fObsWall = toyMap([3.0; 2.7; 5.5]);   % x_true = [1.5; 0.9; 1.1]
wallFwd = @(p) wallMap(p, pRef, wall, 'toyModel:leftValidRegion');

% (4a) Expected error mid-solve: recorded, not thrown.
res = solveInverse(wallFwd, fObsWall, pRef, x0, lb, ub, ...
    'ExpectedErrorIDs', {'toyModel:leftValidRegion'});
assert(strcmp(res.status, 'leftValidRegion'), 'Test4a:status', 'status = %s, expected leftValidRegion', res.status);
assert(~isempty(res.failedX) && res.failedX(1) > wall, 'Test4a:failedX', ...
    'failedX does not record the point beyond the wall');
assert(max(abs(res.failedP - res.failedX .* pRef)) < 1e-14, 'Test4a:failedP', 'failedP ~= failedX.*pRef');
assert(any(strcmp([{res.error.identifier}, res.error.causeIdentifiers], 'toyModel:leftValidRegion')), ...
    'Test4a:errorRecord', 'recorded error does not carry the expected identifier');
assert(all(isnan(res.xHat)) && isnan(res.resnorm) && res.funcCount > 0, 'Test4a:fields', ...
    'aborted solve should have NaN xHat/resnorm and a positive funcCount');
fprintf('  (4a) expected error mid-solve recorded as leftValidRegion at x = %s: PASS\n', mat2str(res.failedX.', 4));

% (4b) Same error, not declared expected: re-thrown with its identifier.
gotId = '<no error>';
try
    solveInverse(wallFwd, fObsWall, pRef, x0, lb, ub);
catch err
    gotId = err.identifier;
end
assert(strcmp(gotId, 'toyModel:leftValidRegion'), 'Test4b:rethrow', ...
    'undeclared error was not re-thrown with its identifier (got %s)', gotId);
fprintf('  (4b) the same error is re-thrown when not declared expected: PASS\n');

% (4c) Unexpected identifier: re-thrown even though other IDs are expected.
bugFwd = @(p) wallMap(p, pRef, wall, 'toyModel:genuineBug');
gotId = '<no error>';
try
    solveInverse(bugFwd, fObsWall, pRef, x0, lb, ub, 'ExpectedErrorIDs', {'toyModel:leftValidRegion'});
catch err
    gotId = err.identifier;
end
assert(strcmp(gotId, 'toyModel:genuineBug'), 'Test4c:unexpected', ...
    'an unexpected error identifier was not re-thrown (got %s)', gotId);
fprintf('  (4c) an unexpected identifier is re-thrown, never recorded as a finding: PASS\n');

% (4d) x0 itself outside the valid region: fails on the first evaluation.
x0bad = [1.4; 1; 1];
res = solveInverse(wallFwd, fObsWall, pRef, x0bad, lb, ub, ...
    'ExpectedErrorIDs', {'toyModel:leftValidRegion'});
assert(strcmp(res.status, 'leftValidRegion') && isequal(res.failedX, x0bad), 'Test4d:initial', ...
    'failure at the initial point was not recorded as leftValidRegion at x0');
fprintf('  (4d) failure at the very first evaluation (x0) is recorded: PASS\n');

% (4e) ExpectedErrorIDs as a plain char string.
res = solveInverse(wallFwd, fObsWall, pRef, x0, lb, ub, 'ExpectedErrorIDs', 'toyModel:leftValidRegion');
assert(strcmp(res.status, 'leftValidRegion'), 'Test4e:charID', 'char ExpectedErrorIDs not accepted');
fprintf('  (4e) ExpectedErrorIDs accepted as a single char string: PASS\n\n');

%% TEST 5: FeatureMask forwarded
fprintf('--- TEST 5: FeatureMask -----------------------------------\n');

nanFwd = @(p) setElem(toyMap(p), 2, NaN);   % feature 2 never available
res = solveInverse(nanFwd, fObs, pRef, x0, lb, ub, 'FeatureMask', [1 3 4 5]);
relErr = abs(res.pHat - pTrue) ./ pTrue;
assert(strcmp(res.status, 'solved') && max(relErr) < 1e-8, 'Test5:recovery', ...
    'masked solve did not recover the truth (status %s, max rel. error %.3g)', res.status, max(relErr));
assert(isequal(size(res.residual), [4, 1]) && isequal(size(res.jacobian), [4, 3]), 'Test5:sizes', ...
    'residual/Jacobian not sized to the 4 selected features');
fprintf('  (5) masked-out NaN feature: solves, residual 4x1, Jacobian 4x3: PASS\n\n');

%% TEST 6: budget exhaustion
fprintf('--- TEST 6: budget exhaustion ------------------------------\n');

res = solveInverse(toyFwd, fObs, pRef, x0, lb, ub, 'MaxIterations', 1);
assert(strcmp(res.status, 'budgetExhausted') && res.exitflag == 0, 'Test6:budget', ...
    'MaxIterations = 1 gave status %s / exitflag %g, expected budgetExhausted / 0', res.status, res.exitflag);
fprintf('  (6) MaxIterations = 1 -> budgetExhausted (exitflag 0): PASS\n\n');

fprintf('================================================================\n');
fprintf(' ALL solveInverseTester.m TESTS PASSED\n');
fprintf('================================================================\n');

%% ------------------------------------------------------------------------
function f = toyMap(p)
% 3 parameters -> 5 features, nonlinear and identifiable.
f = [p(1) + p(2); p(1) * p(3); p(2)^2; p(1) + p(3); 2 * p(2)];
end

function f = wallMap(p, pRef, wall, id)
% toyMap, but throws error `id` when x(1) = p(1)/pRef(1) exceeds `wall`
% (a stand-in for modalFeatureVector leaving its valid region under Strict).
if p(1) / pRef(1) > wall
    error(id, 'toy forward model: x(1) = %.4f is beyond the valid region (%.2f)', p(1) / pRef(1), wall);
end
f = toyMap(p);
end

function out = recFwd(cmd, p)
% toyMap that also records every point it is evaluated at.
persistent log
switch cmd
    case 'reset'
        log = zeros(0, 3); out = [];
    case 'eval'
        log(end+1, :) = p(:).'; out = toyMap(p);
    case 'get'
        out = log;
end
end

function [plusStep, minusStep] = fdStepsAround(X, x0)
% For each coordinate j, find evaluated points equal to x0 in every other
% coordinate and offset in coordinate j only; return the smallest positive
% and negative offsets found (NaN if none), i.e. the FD steps used at x0.
n = numel(x0);
plusStep = NaN(n, 1); minusStep = NaN(n, 1);
for j = 1:n
    others = setdiff(1:n, j);
    onAxis = all(abs(X(:, others) - x0(others).') < 1e-12, 2);
    d = X(onAxis, j) - x0(j);
    if any(d > 1e-12),  plusStep(j)  = min(d(d > 1e-12)); end
    if any(d < -1e-12), minusStep(j) = min(-d(d < -1e-12)); end
end
end

function v = setElem(v, k, val)
% Return v with element k replaced by val.
v(k) = val;
end