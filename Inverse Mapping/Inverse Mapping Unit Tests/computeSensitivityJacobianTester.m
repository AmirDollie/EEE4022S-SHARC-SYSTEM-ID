%% computeSensitivityJacobianTester.m
% Validation tests for computeSensitivityJacobian.m and its pure helper
% centralDiffJacobian.m.
%
% TEST LIST
%   1. centralDiffJacobian.m  -- pure finite-difference arithmetic
%                                 (column indexing, the 2*deltaP
%                                 normalisation), no EMM involved.
%   2. computeSensitivityJacobian.m input validation -- bad beta/gamma/R/
%                                 epsilon, and the reserved 'Strict'
%                                 option, all rejected before any EMM
%                                 call is made.
%   3. Forwarding + fail-path integration against the real Forward Model
%      -- confirms options are actually passed through to
%      modalFeatureVector.m, and that an invalid baseline fit correctly
%      aborts the whole Jacobian (propagating
%      modalFeatureVector:invalidFeatureFit) rather than silently
%      producing a contaminated J.
%
% WHAT THIS TESTER DELIBERATELY DOES NOT COVER: a full, successful 7-call
% Jacobian (baseline + 3 parameters x 2 sides, all genuinely clean under
% Strict) against real EMM output. Finding a synthetic (beta,gamma,R)
% where all 5 features have a genuine, non-degenerate, in-bracket local
% extremum requires essentially re-running a scoped-down
% modalRAOReconnaissance.m/modalFeatureRefinement.m by hand -- expensive,
% and not needed to validate this function's OWN logic, which is fully
% covered by: the pure arithmetic (TEST 1), the guards (TEST 2), and
% real-EMM-backed proof that a bad fit is correctly propagated and aborts
% the computation (TEST 3c). The genuine "does the full Jacobian run
% cleanly" check belongs at the one point you have already validated for
% real: your actual (beta0,gamma0,R0) reference point, which
% modalFeatureVectorTester.m/your own manual modalFeatureVector run has
% already confirmed is clean for all 5 features under 'Strict',true.
% Run computeSensitivityJacobian there directly -- see the note at the
% end of this file.
%
% Lives in Inverse Mapping/Inverse Mapping Unit Tests/, alongside
% modalFeatureVectorTester.m; computeSensitivityJacobian.m and
% centralDiffJacobian.m themselves stay one level up in
% Inverse Mapping/, per this repo's convention (see
% modalFeatureVectorTester.m's own header for the fuller rationale).
% Reuses makeSyntheticFeatureRefinement.m from that same tester (also in
% this folder) for TEST 3's synthetic anchor file.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);                                              % makeSyntheticFeatureRefinement.m
addpath(fullfile(thisDir, '..'));                              % computeSensitivityJacobian.m, centralDiffJacobian.m, modalFeatureVector.m, etc.
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));

clearvars -except thisDir
close all, clc

scratchDir = tempname;
mkdir(scratchDir);
cleanupObj = onCleanup(@() rmdir(scratchDir, 's'));

fprintf('================================================================\n');
fprintf(' computeSensitivityJacobianTester.m\n');
fprintf(' scratch folder: %s\n', scratchDir);
fprintf('================================================================\n\n');

%% TEST 1: centralDiffJacobian.m -- pure arithmetic
fprintf('--- TEST 1: centralDiffJacobian.m ---------------------------\n');

% (1a) Known synthetic case, hand-computed: 3 features, 2 parameters.
fPlus  = [1.20, 5.00; 2.40, 5.10; 0.90, 5.30];   % 3x2
fMinus = [1.00, 5.02; 2.00, 5.06; 1.10, 5.20];   % 3x2
deltaP = [0.10, 0.02];                            % 1x2
J = centralDiffJacobian(fPlus, fMinus, deltaP);
expectedJ = (fPlus - fMinus) ./ (2 * deltaP);     % same formula, computed independently here
assert(isequal(size(J), [3, 2]), 'Test1a:size', 'J is %s, expected 3x2', mat2str(size(J)));
assert(max(abs(J(:) - expectedJ(:))) < 1e-14, 'Test1a:values', 'J does not match hand-computed expectation');
% Spot-check one entry by hand: column 1, row 1: (1.20-1.00)/(2*0.10) = 1.0
assert(abs(J(1,1) - 1.0) < 1e-14, 'Test1a:spotCheck', 'J(1,1) = %.10f, expected 1.0', J(1,1));
fprintf('  (1a) known synthetic case matches hand computation exactly: PASS\n');

% (1b) deltaP given as a column vector must work identically (orientation
%      should not matter).
Jcol = centralDiffJacobian(fPlus, fMinus, deltaP(:));
assert(isequal(J, Jcol), 'Test1b:orientation', 'centralDiffJacobian gave a different result for deltaP as a column vs row vector');
fprintf('  (1b) deltaP orientation (row vs column) does not affect the result: PASS\n');

% (1c) Error handling: size mismatch between fPlus/fMinus.
threw = false;
try
    centralDiffJacobian([1 2; 3 4], [1 2 3; 4 5 6], [0.1, 0.1]);
catch err
    threw = strcmp(err.identifier, 'centralDiffJacobian:sizeMismatch');
end
assert(threw, 'Test1c:sizeMismatch', 'did not throw centralDiffJacobian:sizeMismatch for mismatched fPlus/fMinus sizes');

% (1d) Error handling: deltaP with the wrong number of elements.
threw = false;
try
    centralDiffJacobian([1 2; 3 4], [1 2; 3 4], [0.1, 0.1, 0.1]);
catch err
    threw = strcmp(err.identifier, 'centralDiffJacobian:deltaPSizeMismatch');
end
assert(threw, 'Test1d:deltaPSize', 'did not throw centralDiffJacobian:deltaPSizeMismatch for a wrongly-sized deltaP');

% (1e) Error handling: a non-positive or non-finite deltaP entry.
threw = false;
try
    centralDiffJacobian([1 2], [1 2], [0.1, 0]);
catch err
    threw = strcmp(err.identifier, 'centralDiffJacobian:badDeltaP');
end
assert(threw, 'Test1e:zeroDeltaP', 'did not throw centralDiffJacobian:badDeltaP for a zero deltaP entry');

threw = false;
try
    centralDiffJacobian([1 2], [1 2], [0.1, -0.05]);
catch err
    threw = strcmp(err.identifier, 'centralDiffJacobian:badDeltaP');
end
assert(threw, 'Test1e:negDeltaP', 'did not throw centralDiffJacobian:badDeltaP for a negative deltaP entry');
fprintf('  (1c-e) error handling (size mismatches, zero/negative deltaP): PASS\n\n');

%% TEST 2: computeSensitivityJacobian.m -- input validation (no EMM cost;
%  every case here must fail before any modalFeatureVector call is made)
fprintf('--- TEST 2: computeSensitivityJacobian.m input validation ---\n');

badCases = {
    'badBeta',  {-1,    0.1,  1.0, 1e-3}, 'computeSensitivityJacobian:badBeta';
    'badBeta2', {0,     0.1,  1.0, 1e-3}, 'computeSensitivityJacobian:badBeta';
    'badGamma', {1e-2,  -0.1, 1.0, 1e-3}, 'computeSensitivityJacobian:badGamma';
    'badR',     {1e-2,  0.1,  -1,  1e-3}, 'computeSensitivityJacobian:badR';
    'badEps0',  {1e-2,  0.1,  1.0, 0},    'computeSensitivityJacobian:badEpsilon';
    'badEpsNeg',{1e-2,  0.1,  1.0, -1e-3},'computeSensitivityJacobian:badEpsilon';
    'badEpsGE1',{1e-2,  0.1,  1.0, 1.5},  'computeSensitivityJacobian:badEpsilon';
};
for r = 1:size(badCases, 1)
    args = badCases{r, 2};
    expectedId = badCases{r, 3};
    threw = false;
    try
        computeSensitivityJacobian(args{1}, args{2}, args{3}, args{4});
    catch err
        threw = strcmp(err.identifier, expectedId);
    end
    assert(threw, sprintf('Test2:%s', badCases{r,1}), ...
        'case ''%s'' did not throw %s', badCases{r,1}, expectedId);
end
fprintf('  (2a) bad beta/gamma/R/epsilon are all rejected with the correct identifiers: PASS\n');

% 'Strict' is reserved and must be rejected regardless of its value.
threw = false;
try
    computeSensitivityJacobian(1e-2, 0.1, 1.0, 1e-3, 'Strict', true);
catch err
    threw = strcmp(err.identifier, 'computeSensitivityJacobian:strictReserved');
end
assert(threw, 'Test2:strictTrue', 'did not throw computeSensitivityJacobian:strictReserved for ''Strict'',true');

threw = false;
try
    computeSensitivityJacobian(1e-2, 0.1, 1.0, 1e-3, 'Strict', false);
catch err
    threw = strcmp(err.identifier, 'computeSensitivityJacobian:strictReserved');
end
assert(threw, 'Test2:strictFalse', 'did not throw computeSensitivityJacobian:strictReserved for ''Strict'',false (still reserved regardless of value)');
fprintf('  (2b) ''Strict'' is rejected in either case, before any EMM call: PASS\n\n');

%% TEST 3: forwarding + fail-path integration against the real Forward
%  Model (never touches the real modalFeatureRefinementResults.mat).
fprintf('--- TEST 3: forwarding and fail-path integration -------------\n');

[refined, modeList, nu, M, P, N, H, g] = makeSyntheticFeatureRefinement(); %#ok<ASGLU>
smokeFile = fullfile(scratchDir, 'smoke.mat');
save(smokeFile, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g');

betaTest = 1e-2; gammaTest = 0.1; Rtest = 1.0;

% (3a) A bad 'RefinementResultsFile' must propagate modalFeatureVector's
%      own error unchanged (proves the option is actually forwarded, and
%      cheaply -- no EMM solve happens before this file is checked).
threw = false;
try
    computeSensitivityJacobian(betaTest, gammaTest, Rtest, 1e-3, ...
        'RefinementResultsFile', fullfile(scratchDir, 'does_not_exist.mat'));
catch err
    threw = strcmp(err.identifier, 'modalFeatureVector:missingRefinementResults');
end
assert(threw, 'Test3a:missingFile', 'a missing RefinementResultsFile did not propagate modalFeatureVector:missingRefinementResults');
fprintf('  (3a) ''RefinementResultsFile'' is forwarded (missing-file error propagates): PASS\n');

% (3b) A non-integer 'NTheta' must be rejected by modalFeatureVector's
%      own inputParser -- proves NTheta is forwarded too, and cheaply
%      (fails at argument-parsing time, before any EMM solve).
threw = false;
try
    computeSensitivityJacobian(betaTest, gammaTest, Rtest, 1e-3, ...
        'RefinementResultsFile', smokeFile, 'NTheta', 45.5);
catch %#ok<CTCH>
    threw = true; % identifier differs between MATLAB/Octave inputParser versions -- see modalFeatureVectorTester.m TEST 6
end
assert(threw, 'Test3b:badNTheta', 'a non-integer NTheta was not rejected -- ''NTheta'' does not appear to be forwarded to modalFeatureVector');
fprintf('  (3b) ''NTheta'' is forwarded (non-integer value rejected): PASS\n');

% (3c) Fail-path integration against the REAL Forward Model: the default
%      makeSyntheticFeatureRefinement() anchors are arbitrary and
%      uninformed for this cheap test floe -- modalFeatureVectorTester.m
%      (TEST 4/5) already demonstrated that this exact configuration
%      genuinely fails 'Strict',true at these parameters (extrapolated
%      and/or kindMismatch on real EMM output, not a mocked value). That
%      means the very first (baseline) call computeSensitivityJacobian
%      makes here is expected to throw, before any perturbed evaluation
%      -- confirming the whole call chain (argument forwarding, the
%      Strict-always-on baseline call, and error propagation) is wired
%      correctly against real physics, end to end.
ticS = tic;
threw = false;
try
    computeSensitivityJacobian(betaTest, gammaTest, Rtest, 1e-3, 'RefinementResultsFile', smokeFile);
catch err
    threw = strcmp(err.identifier, 'modalFeatureVector:invalidFeatureFit');
end
fprintf('  (real EMM baseline call completed in %.1fs)\n', toc(ticS));
assert(threw, 'Test3c:failPath', ...
    'an invalid baseline fit did not propagate modalFeatureVector:invalidFeatureFit -- a bad fit may be getting silently absorbed');
fprintf('  (3c) an invalid baseline correctly aborts the WHOLE Jacobian (real EMM, ''Strict'',true propagates): PASS\n\n');

fprintf('================================================================\n');
fprintf(' ALL computeSensitivityJacobianTester.m TESTS PASSED\n');
fprintf('================================================================\n\n');
fprintf('NOTE: this tester does not exercise a full, successful 7-call\n');
fprintf('Jacobian -- see the header comment. Run it directly at your own\n');
fprintf('already-validated reference point instead, e.g.:\n');
fprintf('  [J, info] = computeSensitivityJacobian(beta0, gamma0, R0, 1e-3);\n');