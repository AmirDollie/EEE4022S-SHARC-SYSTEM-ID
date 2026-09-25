%% featureInformationTester.m
% Unit tests for featureInformation.m (Inverse Mapping/). Pure linear
% algebra, no EMM; runs in well under a second.
%
% Groups
%   1  Invertible linear transform of the reference, noise propagated
%      consistently -> lambda = 1 exactly, stds identical
%   2  Row subset / rank-reducing linear map (consistent propagation)
%      -> 0 <= lambda <= 1, stds can only grow
%   3  Feature rescaling (JF -> D JF, SF -> D SF D') leaves everything unchanged
%   4  Known analytic case: diagonal information, lambda and V known exactly
%   5  Inconsistent (understated) feature noise -> lambda > 1, flagged + warning
%   6  Sigma forms: scalar, vector and matrix give the same answer
%   7  A direction completely lost -> lambda = 0, stdF = Inf
%   8  Generalised eigen relation F_f V = F_r V diag(lambda) holds
%   9  Error handling
%
% Lives in Inverse Mapping/Inverse Mapping Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clearvars -except thisDir
close all, clc

nPass = 0;
rng0 = 12345;
if exist('OCTAVE_VERSION', 'builtin') > 0, randn('state', rng0); else, rng(rng0); end %#ok<RAND>

%% common random reference: 40 observables, 2 parameters, correlated noise
mR = 40; nP = 2;
Jr = randn(mR, nP);
B = randn(mR); Sr = B * B.' / mR + 0.1 * eye(mR);        % SPD reference covariance

%% 1  invertible transform -> lambda = 1
fprintf('=== 1: invertible linear transform ===\n');
A = randn(mR) + 3 * eye(mR);
o = featureInformation(A * Jr, A * Sr * A.', Jr, Sr);
assert(max(abs(o.lambda - 1)) < 1e-8, sprintf('1: lambda should be 1 (got %s)', mat2str(o.lambda.', 6)));
assert(max(abs(o.stdRatio - 1)) < 1e-8, '1: stds should be identical');
assert(abs(o.detRatio - 1) < 1e-7 && o.consistent, '1: detRatio / consistent wrong');
fprintf('PASS: lambda = %s\n', mat2str(o.lambda.', 10)); nPass = nPass + 1;

%% 2  subsets and rank-reducing maps -> 0 <= lambda <= 1
fprintf('\n=== 2: information cannot increase ===\n');
worst = -Inf;
for t = 1:200
    k = randi([nP, mR - 1]);
    A = randn(k, mR);                                    % generic rank-reducing map
    if mod(t, 2) == 0
        idx = randperm(mR, k); A = eye(mR); A = A(idx, :);   % plain row subset
    end
    o = featureInformation(A * Jr, A * Sr * A.', Jr, Sr);
    assert(all(o.lambda >= -1e-10) && all(o.lambda <= 1 + 1e-8), sprintf('2.%d: lambda outside [0,1]: %s', t, mat2str(o.lambda.', 6)));
    assert(all(o.stdRatio >= 1 - 1e-8), '2: a std decreased');
    assert(o.consistent, '2: flagged inconsistent');
    worst = max(worst, max(o.lambda));
end
fprintf('PASS: 200 random maps, max lambda %.6f\n', worst); nPass = nPass + 1;

%% 3  rescaling invariance
fprintf('\n=== 3: feature rescaling ===\n');
A = randn(10, mR);
Jf = A * Jr; Sf = A * Sr * A.';
D = diag(10.^(4 * rand(10, 1) - 2));
o1 = featureInformation(Jf, Sf, Jr, Sr);
o2 = featureInformation(D * Jf, D * Sf * D, Jr, Sr);
assert(max(abs(o1.lambda - o2.lambda)) < 1e-9 * max(o1.lambda), '3: lambda changed under rescaling');
assert(max(abs(o1.stdF - o2.stdF) ./ o1.stdF) < 1e-9, '3: stds changed under rescaling');
assert(max(abs(abs(o1.V(:)) - abs(o2.V(:)))) < 1e-8, '3: eigenvectors changed under rescaling');
fprintf('PASS\n'); nPass = nPass + 1;

%% 4  analytic case
fprintf('\n=== 4: analytic case ===\n');
% reference: iid unit noise, F_r = diag(4, 9); features keep half of p1's
% information and a ninth of p2's -> lambda = [0.5; 1/9], V = I
JrA = [2 0; 0 3]; JfA = [sqrt(2) 0; 0 1];
o = featureInformation(JfA, 1, JrA, 1, 'ParamNames', {'lnBeta', 'lnR'});
assert(max(abs(o.lambda - [0.5; 1/9])) < 1e-12, '4: lambda wrong');
assert(max(max(abs(abs(o.V) - eye(2)))) < 1e-12, '4: V wrong');
assert(max(max(abs(abs(o.Vscaled) - eye(2)))) < 1e-12, '4: Vscaled wrong');
assert(max(abs(o.stdR - [1/2; 1/3])) < 1e-12 && max(abs(o.stdF - [1/sqrt(2); 1])) < 1e-12, '4: stds wrong');
assert(max(abs(o.stdRatio - 1 ./ sqrt([0.5; 1/9]))) < 1e-12, '4: stdRatio wrong');
assert(isequal(o.paramNames, {'lnBeta', 'lnR'}), '4: names not stored');
fprintf('PASS: lambda = %s\n', mat2str(o.lambda.', 6)); nPass = nPass + 1;

%% 5  inconsistent noise flagged
fprintf('\n=== 5: understated feature noise ===\n');
A = randn(mR) + 3 * eye(mR);
lastwarn('', '');
fprintf('(the warning printed next is EXPECTED)\n');
o = featureInformation(A * Jr, 0.25 * (A * Sr * A.'), Jr, Sr);   % noise understated 4x
[~, id] = lastwarn;
assert(~o.consistent && max(o.lambda) > 3.9, '5: inconsistency not detected');
assert(strcmp(id, 'featureInformation:lambdaAboveOne'), '5: warning not raised');
lastwarn('', '');
o = featureInformation(A * Jr, 0.25 * (A * Sr * A.'), Jr, Sr, 'WarnInconsistent', false);
[~, id] = lastwarn;
assert(isempty(id) && ~o.consistent, '5: WarnInconsistent=false should silence the warning only');
fprintf('PASS: lambda = %s\n', mat2str(o.lambda.', 4)); nPass = nPass + 1;

%% 6  Sigma forms agree
fprintf('\n=== 6: scalar / vector / matrix Sigma ===\n');
v = 0.3 + rand(mR, 1);
oS = featureInformation(Jr, 2.0, Jr, v);
oV = featureInformation(Jr, 2.0 * ones(mR, 1), Jr, diag(v));
oM = featureInformation(Jr, 2.0 * eye(mR), Jr, v);
assert(max(abs(oS.lambda - oV.lambda)) < 1e-12 && max(abs(oS.lambda - oM.lambda)) < 1e-12, '6: Sigma forms disagree');
assert(max(max(abs(oS.Fr - oV.Fr))) < 1e-10 * max(abs(oS.Fr(:))), '6: Fr differs between forms');
fprintf('PASS\n'); nPass = nPass + 1;

%% 7  completely lost direction
fprintf('\n=== 7: lost direction ===\n');
% consistent map that is blind to p2: rows orthogonal to Jr(:, 2)
Nb = null(Jr(:, 2).');
A = Nb(:, 1:5).';
o = featureInformation(A * Jr, A * Sr * A.', Jr, Sr);
assert(min(o.lambda) < 1e-10 && all(isinf(o.stdF)) && o.consistent, '7: lost direction not reported');
fprintf('PASS: lambda = %s\n', mat2str(o.lambda.', 4)); nPass = nPass + 1;

%% 8  generalised eigen relation
fprintf('\n=== 8: F_f V = F_r V diag(lambda) ===\n');
A = randn(15, mR);
o = featureInformation(A * Jr, A * Sr * A.', Jr, Sr);
res = norm(o.Ff * o.V - o.Fr * o.V * diag(o.lambda)) / norm(o.Ff);
assert(res < 1e-10, '8: eigen relation violated');
assert(max(abs(vecnorm(o.V, 2, 1) - 1)) < 1e-12, '8: V columns not unit norm');
fprintf('PASS: relative residual %.1e\n', res); nPass = nPass + 1;

%% 9  error handling
fprintf('\n=== 9: error handling ===\n');
cases = {
    'badJ',                     {[], 1, Jr, 1}
    'badJ',                     {Jr + 1i, 1, Jr, 1}
    'sizeMismatch',             {Jr(:, 1), 1, Jr, 1}
    'sizeMismatch',             {Jr, ones(mR - 1, 1), Jr, 1}
    'badSigma',                 {Jr, NaN, Jr, 1}
    'badSigma',                 {Jr, triu(ones(mR)) + eye(mR), Jr, 1}
    'sigmaNotPositiveDefinite', {Jr, -1, Jr, 1}
    'sigmaNotPositiveDefinite', {Jr, [-1; ones(mR - 1, 1)], Jr, 1}
    'sigmaNotPositiveDefinite', {Jr, -eye(mR), Jr, 1}
    'referenceNotInformative',  {Jr, 1, [Jr(:, 1), 2 * Jr(:, 1)], 1}
    'badOption',                {Jr, 1, Jr, 1, 'Nope', 1}
    'badOption',                {Jr, 1, Jr, 1, 'ParamNames', {'a'}}
    'badOption',                {Jr, 1, Jr, 1, 'Tol'}
    };
for c = 1:size(cases, 1)
    id = ['featureInformation:', cases{c, 1}];
    try
        featureInformation(cases{c, 2}{:});
        error('tester:noError', '9.%d: expected %s but no error was thrown', c, id);
    catch err
        assert(strcmp(err.identifier, id), sprintf('9.%d: expected %s, got %s (%s)', c, id, err.identifier, err.message));
    end
end
fprintf('PASS: %d error cases\n', size(cases, 1)); nPass = nPass + 1;

fprintf('\nAll %d groups passed.\n', nPass);