function [bestIdx, bestScore, allScores] = exhaustiveSensorSelection(candidateRows, A, N, metricFn, numSensors)
%EXHAUSTIVESENSORSELECTION Brute-force optimal sensor selection,
%exploiting the fact that per-sensor Gramian contributions ADD:
%Wo(S) = sum_{j in S} Wo_j, where Wo_j is sensor j's own single-row
%contribution. Precomputing all single-sensor contributions once means
%evaluating each candidate 4-subset only costs summing 4 small matrices,
%not rebuilding the Gramian from scratch each time.
%
%   Use this to check whether GREEDY (for lambdaMin/kappa specifically,
%   which have no submodularity guarantee) is finding the true global
%   optimum or a path-dependent local one -- exactly the question raised
%   after Level 2C's initial greedy-only run.

    numCandidates = size(candidateRows, 1);
    n = size(A,1);

    % Precompute each candidate's own single-sensor Gramian contribution once
    singleWo = zeros(n, n, numCandidates);
    for j = 1:numCandidates
        singleWo(:,:,j) = finiteObservabilityGramian(A, candidateRows(j,:), N);
    end

    combos = nchoosek(1:numCandidates, numSensors);
    allScores = zeros(size(combos,1), 1);

    for c = 1:size(combos,1)
        Wo_combo = zeros(n,n);
        for j = combos(c,:)
            Wo_combo = Wo_combo + singleWo(:,:,j);
        end
        allScores(c) = metricFn(Wo_combo);
    end

    [bestScore, bestRow] = max(allScores);
    bestIdx = sort(combos(bestRow, :));
end