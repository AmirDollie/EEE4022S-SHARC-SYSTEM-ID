%% exhaustiveSensorSelectionTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

omega1 = 2.70; omega2 = 6.00; dt = 0.02;
A = buildKnownOscillatorA(omega1, omega2, dt);
N = 100;

candidateRows = [buildCandidateCrows(1+0.5i, 0.3-0.2i);
                 buildCandidateCrows(-0.2+0.9i, 0.7-0.1i);
                 buildCandidateCrows(0.4-0.6i, -0.5+0.3i);
                 buildCandidateCrows(0.1+0.1i, 0.9+0.9i);
                 buildCandidateCrows(-0.8+0.2i, 0.2-0.7i)];

%% Test 1: additive shortcut matches the DIRECT (slow, obviously correct)
% computation -- build the full candidate C matrix for a given subset
% and compute Wo directly, versus summing precomputed single-sensor
% contributions. These must agree exactly (up to roundoff), since this
% is just linearity of the defining sum, not an approximation.
fprintf('=== Test 1: additive shortcut matches direct computation ===\n');
traceFn = @(Wo) gramianMetrics(Wo).trace;

testSubsets = {[1,2], [1,3,5], [2,3,4,5], [1,2,3,4,5]};
for t = 1:length(testSubsets)
    subset = testSubsets{t};
    C_direct = candidateRows(subset, :);
    Wo_direct = finiteObservabilityGramian(A, C_direct, N);

    Wo_additive = zeros(4,4);
    for j = subset
        Wo_additive = Wo_additive + finiteObservabilityGramian(A, candidateRows(j,:), N);
    end

    err = max(abs(Wo_direct(:) - Wo_additive(:)));
    fprintf('  subset %s: max|direct-additive| = %.3e\n', mat2str(subset), err);
    assert(err < 1e-10, 'Test 1 failed: additive shortcut does not match direct computation for subset %s.', mat2str(subset));
end
fprintf('PASS: additive shortcut matches direct computation exactly for all tested subsets.\n');

%% Test 2: exhaustiveSensorSelection itself matches brute-force-by-hand
% on a small, fully checkable case
fprintf('\n=== Test 2: exhaustiveSensorSelection matches independent brute force ===\n');
[bestIdx, bestScore, allScores] = exhaustiveSensorSelection(candidateRows, A, N, traceFn, 2);

combos = nchoosek(1:5, 2);
independentScores = zeros(size(combos,1),1);
for c = 1:size(combos,1)
    Wo_c = finiteObservabilityGramian(A, candidateRows(combos(c,:), :), N);
    independentScores(c) = traceFn(Wo_c);
end
[indepBest, indepBestRow] = max(independentScores);
indepBestIdx = sort(combos(indepBestRow,:));

fprintf('exhaustiveSensorSelection: %s (score=%.4f)\n', mat2str(bestIdx), bestScore);
fprintf('Independent brute force:   %s (score=%.4f)\n', mat2str(indepBestIdx), indepBest);
assert(isequal(bestIdx, indepBestIdx), 'Test 2 failed: exhaustiveSensorSelection does not match independent brute force.');
assert(abs(bestScore - indepBest) < 1e-10, 'Test 2 failed: score mismatch.');
fprintf('PASS: exhaustiveSensorSelection matches an independently-written brute-force search.\n');

%% Test 3: for the MODULAR trace objective, exhaustive must agree with
% greedy exactly (greedy is provably optimal for trace -- if these ever
% disagree, something is wrong in one of the two implementations)
fprintf('\n=== Test 3: exhaustive agrees with greedy for the modular trace objective ===\n');
greedyIdx = sort(greedySensorSelection(candidateRows, A, N, traceFn, 2));
assert(isequal(bestIdx, greedyIdx), 'Test 3 failed: exhaustive and greedy disagree for trace, which is provably modular -- greedy should be exactly optimal here.');
fprintf('PASS: exhaustive (%s) matches greedy (%s) for trace, as required by modularity.\n', mat2str(bestIdx), mat2str(greedyIdx));

%% Test 4: for a NON-modular objective (lambdaMin), exhaustive score must
% be AT LEAST as good as greedy's (never worse -- exhaustive checks
% everything, so it can only match or beat a heuristic)
fprintf('\n=== Test 4: exhaustive is never worse than greedy for a non-modular objective ===\n');
lambdaMinFn = @(Wo) gramianMetrics(Wo).lambdaMin;
[exhaustiveIdxLM, exhaustiveScoreLM] = exhaustiveSensorSelection(candidateRows, A, N, lambdaMinFn, 3);
[greedyIdxLM, greedyScoreHistLM] = greedySensorSelection(candidateRows, A, N, lambdaMinFn, 3);
greedyScoreLM = greedyScoreHistLM(end);

fprintf('Exhaustive lambdaMin: %s (score=%.6f)\n', mat2str(exhaustiveIdxLM), exhaustiveScoreLM);
fprintf('Greedy lambdaMin:     %s (score=%.6f)\n', mat2str(sort(greedyIdxLM)), greedyScoreLM);
assert(exhaustiveScoreLM >= greedyScoreLM - 1e-10, 'Test 4 failed: exhaustive should never score worse than greedy.');
fprintf('PASS: exhaustive score is >= greedy score, as required (exhaustive checks every combination).\n');

fprintf('\nAll 4 tests passed.\n');