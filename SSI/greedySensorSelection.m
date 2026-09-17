function [selectedIdx, scoreHistory] = greedySensorSelection(candidateRows, A, N, metricFn, numSensors)
%GREEDYSENSORSELECTION Generic greedy forward selection of sensor rows,
%parameterized by any scalar Gramian-based objective function.
%
%   [SELECTEDIDX, SCOREHISTORY] = GREEDYSENSORSELECTION(CANDIDATEROWS, A, N, METRICFN, NUMSENSORS)
%
%   CANDIDATEROWS - M x 4 matrix, one row per candidate sensor location
%       (e.g. from stacking buildCandidateCrows.m outputs), M >= NUMSENSORS
%   A             - the known state matrix (buildKnownOscillatorA.m)
%   N             - finite horizon for the Gramian (finiteObservabilityGramian.m)
%   METRICFN      - function handle taking a Gramian matrix Wo and
%       returning a SCALAR SCORE TO BE MAXIMISED. For an objective that
%       is naturally minimised (e.g. kappa, the condition number),
%       METRICFN must already encode the sign flip -- e.g. pass
%       @(Wo) -gramianMetrics(Wo).kappa, not the raw kappa.
%   NUMSENSORS    - how many sensors to select (e.g. 4)
%
%   Returns SELECTEDIDX, a 1 x NUMSENSORS vector of row indices into
%   CANDIDATEROWS, in the order they were added, and SCOREHISTORY, the
%   best achievable score after each addition (length NUMSENSORS).
%
%   ALGORITHM: at each step, evaluate the objective on every remaining
%   candidate's trial Gramian, COMBINED WITH the sensors already
%   selected (not the candidate alone), and add whichever single
%   candidate produces the best resulting total score. For a fixed
%   current selection, maximising the resulting total score is
%   equivalent to maximising the marginal improvement, so this achieves
%   the same outcome as an explicit marginal-gain computation
%   F(S union {j}) - F(S), without computing that difference directly.
%   Scoring candidates independently of the current selection, then
%   taking the top-N, would be a materially different (and generally
%   worse) procedure than this, and would not match what the reviewed
%   literature (Summers & Lygeros; Bartos & Kerkez; Yamada et al.)
%   actually describes or its guarantees actually cover.
%
%   HANDLING -Inf SCORES (rank-deficient trial Gramians): gramianMetrics
%   deliberately returns logDet=-Inf for a rank-deficient Wo -- entirely
%   expected at early steps (e.g. with only 1 of 4 sensors selected, the
%   state may not yet be fully observable at all). If every remaining
%   candidate scores -Inf at some step, the naive comparison
%   "score > bestScore" with bestScore initialised to -Inf never fires
%   (since -Inf > -Inf is false in MATLAB), leaving no candidate
%   selected. This is handled by always accepting the FIRST candidate
%   evaluated at each step regardless of its score, then only replacing
%   it if a STRICTLY better score is found -- so a step where every
%   candidate ties at -Inf deterministically keeps the first one, and
%   later steps (with more sensors, hence more likely to be rank
%   sufficient) can still escape rank deficiency and discriminate
%   normally.
%
%   GREEDY GUARANTEES ARE OBJECTIVE-SPECIFIC, NOT UNIVERSAL. Per the
%   Level 2C literature review:
%     trace       - MODULAR (Summers & Lygeros): greedy is exactly optimal.
%     rank        - SUBMODULAR (Bartos & Kerkez, citing Cortesi et al.):
%                   greedy provably within 63% of optimal.
%     logDet      - submodular under stated conditions (Yamada et al.):
%                   greedy or SDP relaxation both reasonable.
%     lambdaMin, kappa - NO general submodularity guarantee. Greedy is
%                   used here as a heuristic only for these two, with no
%                   claim of near-optimality.

    numCandidates = size(candidateRows, 1);
    if numSensors > numCandidates
        error('greedySensorSelection:tooFewCandidates', ...
            'Requested %d sensors but only %d candidates available.', numSensors, numCandidates);
    end

    selectedIdx = zeros(1, numSensors);
    scoreHistory = zeros(1, numSensors);
    remaining = 1:numCandidates;

    for step = 1:numSensors
        bestScore = -Inf;
        bestCandidate = -1;
        bestCandidatePos = -1;

        for k = 1:length(remaining)
            candidateIdx = remaining(k);
            trialIdx = [selectedIdx(1:step-1), candidateIdx];
            C_trial = candidateRows(trialIdx, :);

            Wo_trial = finiteObservabilityGramian(A, C_trial, N);
            score = metricFn(Wo_trial);

            % Always accept the first candidate evaluated this step,
            % regardless of score (handles the case where every
            % candidate ties at -Inf, e.g. all trial Gramians are
            % rank-deficient); otherwise only replace on a strict
            % improvement.
            if bestCandidate == -1 || score > bestScore
                bestScore = score;
                bestCandidate = candidateIdx;
                bestCandidatePos = k;
            end
        end

        selectedIdx(step) = bestCandidate;
        scoreHistory(step) = bestScore;
        remaining(bestCandidatePos) = [];
    end
end