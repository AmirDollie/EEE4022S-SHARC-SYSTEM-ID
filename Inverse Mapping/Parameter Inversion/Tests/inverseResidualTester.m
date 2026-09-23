%% inverseResidualTester.m
% Validation tests for inverseResidual.m.
%
% Everything here uses cheap analytic toy forward maps, never the EMM, so
% the whole file runs in well under a second. inverseResidual.m is
% forward-model-agnostic by design, so nothing about its OWN logic
% (unscaling, relative residual, masking, validation, error propagation)
% needs the real physics to be tested. The EMM only enters at
% runEMMTwinInversion.m.
%
% TEST LIST
%   1. Definition: forwardFcn receives exactly p = x.*pRef; the residual
%      is (f - fObs)./fObs; zero at the truth; hand-computed values;
%      row/column orientation; invariance to rescaling the features.
%   2. FeatureMask: logical and index masks agree with the matching subset
%      of the full residual; index order does not matter; NaN / zero
%      entries on UNSELECTED features are tolerated, on selected ones not.
%   3. Error propagation: an error thrown by forwardFcn reaches the
%      caller with its original identifier (so solveInverse.m can
%      classify it), not wrapped or swallowed.
%   4. Input validation: every guard throws its documented identifier.
%   5. Fail-cheap ordering: invalid inputs are rejected BEFORE forwardFcn
%      is called (matters when one forward call is an EMM solve).
%   6. Smoke test: the residual, used as a least-squares objective by a
%      basic optimiser (fminsearch, base MATLAB/Octave), drives a toy
%      inverse problem to its known truth. This is NOT the solver test
%      (that is solveInverseTester.m); it only confirms the residual is a
%      well-formed objective.
%
% Lives in Inverse Mapping/Parameter Inversion/Tests/;
% inverseResidual.m is one level up in Inverse Mapping/Parameter Inversion/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));   % inverseResidual.m

clearvars -except thisDir
clc

fprintf('================================================================\n');
fprintf(' inverseResidualTester.m\n');
fprintf('================================================================\n\n');

% Toy forward map used throughout: 3 parameters -> 5 features, nonlinear,
% and identifiable (p2 from f3, p1 from f1, p3 from f4), mimicking the
% shape of the EMM problem ([beta;gamma;R] -> 5 features).
toyFwd = @(p) [p(1) + p(2); p(1) * p(3); p(2)^2; p(1) + p(3); 2 * p(2)];

%% TEST 1: definition
fprintf('--- TEST 1: definition ------------------------------------\n');

% (1a) Identity forward map with fObs = pRef gives r = x - 1 exactly, which
%      can only happen if forwardFcn receives p = x.*pRef (elementwise, in
%      order) and the residual divides by fObs.
pRef = [2; 3; 5];
x    = [1.1; 0.9; 1.5];
r = inverseResidual(x, pRef, pRef, @(p) p);
assert(isequal(size(r), [3, 1]), 'Test1a:size', 'r is %s, expected 3x1', mat2str(size(r)));
assert(max(abs(r - (x - 1))) < 1e-15, 'Test1a:unscale', ...
    'identity map did not give r = x - 1; forwardFcn is not receiving x.*pRef');
fprintf('  (1a) forwardFcn receives p = x.*pRef; r = x - 1 for the identity map: PASS\n');

% (1b) Zero residual at the truth.
pTrue = [2.2; 2.7; 5.5];
fTrue = toyFwd(pTrue);
r = inverseResidual(pTrue ./ pRef, pRef, fTrue, toyFwd);
assert(isequal(size(r), [5, 1]), 'Test1b:size', 'r is %s, expected 5x1', mat2str(size(r)));
assert(norm(r) < 1e-14, 'Test1b:zeroAtTruth', 'residual at the truth is %.3g, expected ~0', norm(r));
fprintf('  (1b) r = 0 at x = pTrue./pRef: PASS\n');

% (1c) Hand-computed value. x = [1;2;1], pRef = [2;3;5] -> p = [2;6;5]
%      -> toyFwd(p) = [8; 10; 36; 7; 12]. With fObs = [10;10;30;7;10]:
%      r = [(8-10)/10; 0; (36-30)/30; 0; (12-10)/10] = [-0.2; 0; 0.2; 0; 0.2].
r = inverseResidual([1; 2; 1], pRef, [10; 10; 30; 7; 10], toyFwd);
assert(max(abs(r - [-0.2; 0; 0.2; 0; 0.2])) < 1e-15, 'Test1c:handValue', ...
    'r = %s, expected [-0.2;0;0.2;0;0.2]', mat2str(r.', 6));
fprintf('  (1c) matches hand-computed residual: PASS\n');

% (1d) Orientation: row x / pRef / fObs, and a forward map returning a
%      ROW vector, give the identical column result.
toyFwdRow = @(p) toyFwd(p).';
rCol = inverseResidual([1; 2; 1], pRef, [10; 10; 30; 7; 10], toyFwd);
rRow = inverseResidual([1, 2, 1], pRef.', [10, 10, 30, 7, 10], toyFwdRow);
assert(isequal(rCol, rRow) && iscolumn(rRow), 'Test1d:orientation', ...
    'row/column inputs gave a different (or non-column) residual');
fprintf('  (1d) row/column orientation of inputs and forward output does not matter: PASS\n');

% (1e) Relative residual: scaling BOTH the forward output and fObs by the
%      same constant (a change of feature units) leaves r unchanged.
c = 37.3;
fObs = [10; 10; 30; 7; 10];
r1 = inverseResidual([1; 2; 1], pRef, fObs, toyFwd);
r2 = inverseResidual([1; 2; 1], pRef, c * fObs, @(p) c * toyFwd(p));
assert(max(abs(r1 - r2)) < 1e-14, 'Test1e:unitInvariance', ...
    'residual changed when features were rescaled; it is not a relative residual');
fprintf('  (1e) residual is invariant to rescaling feature units: PASS\n\n');

%% TEST 2: FeatureMask
fprintf('--- TEST 2: FeatureMask -----------------------------------\n');

xT = [1.05; 0.97; 1.2];
fObs = toyFwd(pTrue);
rFull = inverseResidual(xT, pRef, fObs, toyFwd);

% (2a) Logical mask = matching subset of the full residual.
lmask = logical([1 0 1 0 1]);
rL = inverseResidual(xT, pRef, fObs, toyFwd, 'FeatureMask', lmask);
assert(isequal(rL, rFull(lmask(:))), 'Test2a:logical', 'logical mask did not return the matching subset');
fprintf('  (2a) logical mask returns the matching subset of the full residual: PASS\n');

% (2b) Index masks, in any order, give the same result as the logical mask,
%      always in original feature order.
rI1 = inverseResidual(xT, pRef, fObs, toyFwd, 'FeatureMask', [1 3 5]);
rI2 = inverseResidual(xT, pRef, fObs, toyFwd, 'FeatureMask', [5 1 3]);
rI3 = inverseResidual(xT, pRef, fObs, toyFwd, 'featuremask', [3; 5; 1]);   % case-insensitive name, column indices
assert(isequal(rI1, rL) && isequal(rI2, rL) && isequal(rI3, rL), 'Test2b:index', ...
    'index masks did not match the logical mask / original feature order');
fprintf('  (2b) index masks match the logical mask regardless of index order: PASS\n');

% (2c) NaN on an UNSELECTED feature is tolerated; on a selected one it is
%      rejected with inverseResidual:nonFiniteForward.
fwdNaN2 = @(p) setElem(toyFwd(p), 2, NaN);   % toy map with feature 2 = NaN
r = inverseResidual(xT, pRef, fObs, fwdNaN2, 'FeatureMask', [1 3 4 5]);
assert(all(isfinite(r)) && numel(r) == 4, 'Test2c:nanUnselected', 'NaN on an unselected feature was not tolerated');
threw = false;
try
    inverseResidual(xT, pRef, fObs, fwdNaN2);
catch err
    threw = strcmp(err.identifier, 'inverseResidual:nonFiniteForward');
end
assert(threw, 'Test2c:nanSelected', 'NaN on a selected feature did not throw inverseResidual:nonFiniteForward');
fprintf('  (2c) non-finite forward output: tolerated if masked out, rejected if selected: PASS\n');

% (2d) fObs = 0 on an UNSELECTED feature is tolerated; on a selected one
%      it is rejected (it is a divisor).
fObsZ = fObs; fObsZ(2) = 0;
r = inverseResidual(xT, pRef, fObsZ, toyFwd, 'FeatureMask', [1 3 4 5]);
assert(all(isfinite(r)), 'Test2d:zeroUnselected', 'fObs = 0 on an unselected feature was not tolerated');
threw = false;
try
    inverseResidual(xT, pRef, fObsZ, toyFwd);
catch err
    threw = strcmp(err.identifier, 'inverseResidual:badFObs');
end
assert(threw, 'Test2d:zeroSelected', 'fObs = 0 on a selected feature did not throw inverseResidual:badFObs');
fprintf('  (2d) fObs = 0: tolerated if masked out, rejected if selected: PASS\n\n');

%% TEST 3: error propagation from forwardFcn
fprintf('--- TEST 3: error propagation ------------------------------\n');

threw = false;
try
    inverseResidual(xT, pRef, fObs, @(p) toyThrow('toyModel:leftValidRegion'));
catch err
    threw = strcmp(err.identifier, 'toyModel:leftValidRegion');
end
assert(threw, 'Test3:propagation', ...
    'an error thrown by forwardFcn did not reach the caller with its original identifier');
fprintf('  (3) forwardFcn errors propagate with their original identifier (not wrapped): PASS\n\n');

%% TEST 4: input validation
fprintf('--- TEST 4: input validation -------------------------------\n');

f5 = fObs;
% {description, call, expected identifier}
cases = {
    'x contains NaN',             @() inverseResidual([1; NaN; 1], pRef, f5, toyFwd),            'inverseResidual:badX'
    'x has a zero',               @() inverseResidual([1; 0; 1], pRef, f5, toyFwd),              'inverseResidual:badX'
    'x has a negative',           @() inverseResidual([1; -1; 1], pRef, f5, toyFwd),             'inverseResidual:badX'
    'x complex',                  @() inverseResidual([1; 1+1i; 1], pRef, f5, toyFwd),           'inverseResidual:badX'
    'x is a matrix',              @() inverseResidual(ones(3, 2), pRef, f5, toyFwd),             'inverseResidual:badX'
    'x wrong length',             @() inverseResidual([1; 1], pRef, f5, toyFwd),                 'inverseResidual:sizeMismatch'
    'pRef has a zero',            @() inverseResidual([1; 1; 1], [2; 0; 5], f5, toyFwd),         'inverseResidual:badPRef'
    'pRef has a negative',        @() inverseResidual([1; 1; 1], [2; -3; 5], f5, toyFwd),        'inverseResidual:badPRef'
    'pRef has Inf',               @() inverseResidual([1; 1; 1], [2; Inf; 5], f5, toyFwd),       'inverseResidual:badPRef'
    'fObs has NaN',               @() inverseResidual([1; 1; 1], pRef, [1; NaN; 1; 1; 1], toyFwd), 'inverseResidual:badFObs'
    'fObs not numeric',           @() inverseResidual([1; 1; 1], pRef, 'abcde', toyFwd),         'inverseResidual:badFObs'
    'forwardFcn not a handle',    @() inverseResidual([1; 1; 1], pRef, f5, 'toyFwd'),            'inverseResidual:badForwardFcn'
    'forward returns 4 features', @() inverseResidual([1; 1; 1], pRef, f5, @(p) ones(4, 1)),     'inverseResidual:forwardSizeMismatch'
    'forward returns a matrix',   @() inverseResidual([1; 1; 1], pRef, f5, @(p) ones(5, 2)),     'inverseResidual:forwardSizeMismatch'
    'forward returns text',       @() inverseResidual([1; 1; 1], pRef, f5, @(p) 'abcde'),        'inverseResidual:forwardSizeMismatch'
    'logical mask wrong length',  @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 'FeatureMask', true(1, 4)),    'inverseResidual:badMask'
    'logical mask all false',     @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 'FeatureMask', false(1, 5)),   'inverseResidual:badMask'
    'index mask has 0',           @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 'FeatureMask', [0 1]),         'inverseResidual:badMask'
    'index mask out of range',    @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 'FeatureMask', [1 6]),         'inverseResidual:badMask'
    'index mask non-integer',     @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 'FeatureMask', [1 2.5]),       'inverseResidual:badMask'
    'index mask duplicate',       @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 'FeatureMask', [1 1 3]),       'inverseResidual:badMask'
    'mask is a cell',             @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 'FeatureMask', {1, 2}),        'inverseResidual:badMask'
    'odd number of options',      @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 'FeatureMask'),                'inverseResidual:badOption'
    'unknown option name',        @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 'Strict', true),               'inverseResidual:badOption'
    'non-text option name',       @() inverseResidual([1; 1; 1], pRef, f5, toyFwd, 3, true),                      'inverseResidual:badOption'
};
for c = 1:size(cases, 1)
    threw = false; gotId = '<no error>';
    try
        cases{c, 2}();
    catch err
        gotId = err.identifier;
        threw = strcmp(gotId, cases{c, 3});
    end
    assert(threw, sprintf('Test4:case%d', c), '%s: expected %s, got %s', cases{c, 1}, cases{c, 3}, gotId);
end
fprintf('  (4) all %d invalid-input cases throw their documented identifier: PASS\n\n', size(cases, 1));

%% TEST 5: invalid inputs are rejected before forwardFcn is called
fprintf('--- TEST 5: fail-cheap ordering ----------------------------\n');

tripwire = @(p) toyThrow('tester:forwardWasCalled');
cheapCases = {
    @() inverseResidual([1; NaN; 1], pRef, f5, tripwire),                          'inverseResidual:badX'
    @() inverseResidual([1; 1; 1], [2; 0; 5], f5, tripwire),                       'inverseResidual:badPRef'
    @() inverseResidual([1; 1; 1], pRef, [1; 0; 1; 1; 1], tripwire),              'inverseResidual:badFObs'
    @() inverseResidual([1; 1; 1], pRef, f5, tripwire, 'FeatureMask', [1 9]),      'inverseResidual:badMask'
};
for c = 1:size(cheapCases, 1)
    gotId = '<no error>';
    try
        cheapCases{c, 1}();
    catch err
        gotId = err.identifier;
    end
    assert(strcmp(gotId, cheapCases{c, 2}), sprintf('Test5:case%d', c), ...
        'expected %s before any forward call, got %s', cheapCases{c, 2}, gotId);
end
fprintf('  (5) bad x / pRef / fObs / mask are rejected without calling forwardFcn: PASS\n\n');

%% TEST 6: smoke test as a least-squares objective
fprintf('--- TEST 6: residual as an objective (smoke test) ---------\n');

pTrue = [2.2; 2.7; 5.5];
fObs  = toyFwd(pTrue);
obj = @(x) sum(inverseResidual(x, pRef, fObs, toyFwd).^2);
opts = optimset('TolX', 1e-12, 'TolFun', 1e-20, 'MaxFunEvals', 1e4, 'MaxIter', 1e4, 'Display', 'off');
xHat = fminsearch(obj, [1; 1; 1], opts);
pHat = xHat .* pRef;
relErr = abs(pHat - pTrue) ./ pTrue;
assert(max(relErr) < 1e-6, 'Test6:recovery', ...
    'fminsearch on the residual did not recover the toy truth (max relative error %.3g)', max(relErr));
fprintf('  (6) toy truth recovered from x0 = [1;1;1], max relative error %.2g: PASS\n\n', max(relErr));

fprintf('================================================================\n');
fprintf(' ALL inverseResidualTester.m TESTS PASSED\n');
fprintf('================================================================\n');

%% ------------------------------------------------------------------------
function v = setElem(v, k, val)
% Return v with element k replaced by val (MATLAB has no chained indexing
% on a function-call result, so this cannot be done inline in a handle).
v(k) = val;
end

function f = toyThrow(id) %#ok<STOUT>
% Stand-in for a forward model that fails (e.g. modalFeatureVector under
% 'Strict',true leaving its valid region). Declares an output it never
% sets so that it can be called as f = forwardFcn(p): Octave checks the
% output count before running the body, MATLAB would reach error() first.
error(id, 'toy forward model failure (%s)', id);
end