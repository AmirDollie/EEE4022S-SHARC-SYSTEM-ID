%% greedySensorSelectionTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

omega1 = 2.70; omega2 = 6.00; dt = 0.02;
A = buildKnownOscillatorA(omega1, omega2, dt);
N = 500;

%% Test 1: greedy matches brute-force optimum on a small, checkable case
fprintf('=== Test 1: greedy matches brute-force optimum on a small case ===\n');
phis = [1+0i, 0+0i;
        0+0i, 1+0i;
        0.1+0i, 0.1+0i;
        0.05+0i, 0.05+0i];

candidateRows = zeros(4,4);
for k = 1:4
    candidateRows(k,:) = buildCandidateCrows(phis(k,1), phis(k,2));
end

traceMetricFn = @(Wo) gramianMetrics(Wo).trace;

[selectedIdx, scoreHistory] = greedySensorSelection(candidateRows, A, N, traceMetricFn, 2);

pairs = nchoosek(1:4, 2);
bruteScores = zeros(size(pairs,1),1);
for k = 1:size(pairs,1)
    C_pair = candidateRows(pairs(k,:), :);
    bruteScores(k) = traceMetricFn(finiteObservabilityGramian(A, C_pair, N));
end
[bruteBestScore, bruteBestIdx] = max(bruteScores);
bruteBestPair = sort(pairs(bruteBestIdx,:));

fprintf('Greedy selected: %s (score=%.4f)\n', mat2str(sort(selectedIdx)), scoreHistory(end));
fprintf('Brute-force best: %s (score=%.4f)\n', mat2str(bruteBestPair), bruteBestScore);

assert(isequal(sort(selectedIdx), bruteBestPair), 'Test 1 failed: greedy does not match brute-force optimum.');
assert(abs(scoreHistory(end) - bruteBestScore) < 1e-10, 'Test 1 failed: greedy final score mismatch.');
fprintf('PASS: greedy selection exactly matches brute-force optimum (trace is modular, as expected).\n');

%% Test 2: for a monotone maximised objective (trace), the cumulative
% score should be non-decreasing as sensors are added. This is NOT
% claimed for every objective -- e.g. -kappa can genuinely decrease when
% a new sensor makes the configuration worse-conditioned.
fprintf('\n=== Test 2: score history is non-decreasing for the (monotone) trace objective ===\n');
assert(all(diff(scoreHistory) >= -1e-10), 'Test 2 failed: trace-based score should never decrease as sensors are added.');
fprintf('PASS: score history %s is non-decreasing\n', mat2str(scoreHistory));

%% Test 3: marginal-gain-consistent (not independent) scoring changes
% the result for a non-modular objective (lambdaMin)
fprintf('\n=== Test 3: total-score-after-addition scoring changes the result for lambdaMin ===\n');
phisRedundant = [1+0i, 0+0i;
                 1+0i, 0+0i;
                 0+0i, 0.9+0i];

candidateRowsRedundant = zeros(3,4);
for k = 1:3
    candidateRowsRedundant(k,:) = buildCandidateCrows(phisRedundant(k,1), phisRedundant(k,2));
end

lambdaMinMetricFn = @(Wo) gramianMetrics(Wo).lambdaMin;
[selIdx3b, scoreHist3b] = greedySensorSelection(candidateRowsRedundant, A, N, lambdaMinMetricFn, 2);
fprintf('lambdaMin objective: selected %s, scores %s\n', mat2str(sort(selIdx3b)), mat2str(scoreHist3b));
assert(any(selIdx3b == 3), ...
    'Test 3 failed: expected candidate 3 (mode-2 coverage) to be selected, not two redundant mode-1 candidates.');
fprintf('PASS: lambdaMin-based greedy correctly avoids picking two fully redundant candidates.\n');

%% Test 4: minimising kappa via sign flip works correctly
fprintf('\n=== Test 4: kappa minimisation via negated metric function ===\n');
kappaMinFn = @(Wo) -gramianMetrics(Wo).kappa;
[selIdx4, scoreHist4] = greedySensorSelection(candidateRows, A, N, kappaMinFn, 2);
C_final = candidateRows(selIdx4, :);
finalKappa = gramianMetrics(finiteObservabilityGramian(A, C_final, N)).kappa;
fprintf('Selected for MIN kappa: %s, final kappa=%.4f\n', mat2str(selIdx4), finalKappa);
assert(abs(scoreHist4(end) - (-finalKappa)) < 1e-10, 'Test 4 failed: negated score/kappa mismatch.');
fprintf('PASS: kappa minimisation via negation is internally consistent.\n');

%% Test 5: error handling for too few candidates
fprintf('\n=== Test 5: error when requesting more sensors than candidates ===\n');
try
    greedySensorSelection(candidateRows, A, N, traceMetricFn, 10);
    error('Test 5 failed: should have thrown for numSensors > numCandidates.');
catch ME
    assert(strcmp(ME.identifier, 'greedySensorSelection:tooFewCandidates'), 'Test 5 failed: wrong error identifier.');
    fprintf('PASS: correctly threw greedySensorSelection:tooFewCandidates\n');
end
%% Test 6: every candidate scoring -Inf at a step must not break selection
fprintf('\n=== Test 6: all-candidates-score--Inf case is handled deterministically ===\n');

% Artificial metric used ONLY to test the -Inf tie-handling logic, not a
% real observability measure. All one-sensor configurations in this
% test have trace < 750 (confirmed from Test 1's own data: individual
% candidates score well below this), while at least one two-sensor
% configuration (pair [1 2], trace=1000, per Test 1) exceeds 750 --
% forcing exactly one all--Inf step followed by one finite step.
infUntilInformative = @(Wo) testInfMetric(Wo);

[selIdx6, scoreHist6] = greedySensorSelection(candidateRows, A, N, infUntilInformative, 2);

fprintf('Selected: %s, scores: %s\n', mat2str(selIdx6), mat2str(scoreHist6));

assert(scoreHist6(1) == -Inf, 'Test 6 failed: first step was supposed to force an all--Inf tie.');
assert(~any(selIdx6 <= 0), 'Test 6 failed: selection left an invalid index.');
assert(isfinite(scoreHist6(2)), 'Test 6 failed: second step should escape the forced -Inf regime.');

fprintf('PASS: all candidates tied at -Inf on step 1, a valid candidate was still\n');
fprintf('selected deterministically, and step 2 escaped to a finite score.\n');

function score = testInfMetric(Wo)
    tr = trace(Wo);
    if tr < 750
        score = -Inf;
    else
        score = tr;
    end
end