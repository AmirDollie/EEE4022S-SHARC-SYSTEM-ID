function persistence = computeBranchPersistence(branches, windowSize, minStableInWindow)
%COMPUTEBRANCHPERSISTENCE Applies a sliding-window persistence criterion
%to reconstructed modal branches (buildModalBranches.m), rather than
%requiring strict unbroken consecutiveness.
%
%   PERSISTENCE = COMPUTEBRANCHPERSISTENCE(BRANCHES, WINDOWSIZE, MINSTABLEINWINDOW)
%
%   BRANCHES           - struct array from buildModalBranches.m
%   WINDOWSIZE         - W, number of consecutive order-transitions to
%                        look at as one sliding window
%   MINSTABLEINWINDOW  - k, minimum number of class-3 transitions
%                        required within any window of size W for that
%                        window to count as "persistent"
%
%   Returns PERSISTENCE, a struct array of the same length as BRANCHES.
%   persistence(b) has fields:
%   isPersistent    - true if AT LEAST ONE window of size WINDOWSIZE
%                     within this branch's class history contains at
%                     least MINSTABLEINWINDOW class-3 entries
%   bestWindowStart - the order (branches(b).order value, not an index)
%                     at which the best (most class-3-dense) window
%                     starts, or NaN if the branch is too short to form
%                     even one window
%   bestWindowCount - the number of class-3 transitions in that best
%                     window
%
%   DEFINITION: isPersistent answers "does at least one qualifying
%   window exist anywhere in this branch's history" -- NOT "is this
%   branch stable for its entire lifetime." A branch that is solidly
%   class-3 for its first several orders and then degrades completely
%   afterward still registers isPersistent=true, since the qualifying
%   window occurred. This is a deliberate, narrower definition suited to
%   CANDIDATE DISCOVERY (does evidence of stability exist at all) rather
%   than a claim about stability throughout. Formally:
%       isPersistent(branch) <=> max over all windows W of N_3(W) >= k
%   where N_3(W) is the count of class-3 transitions in window W.
%
%   WHY A SLIDING WINDOW, NOT STRICT CONSECUTIVENESS: the Pass 1
%   integration test demonstrated directly that even a single non-class-3
%   transition sitting inside an otherwise long run of class-3 links
%   breaks a strict "N consecutive" requirement, even when the underlying
%   signal is genuinely stable -- see the noisy-data run where all 3 true
%   modes had visible, repeated class-3 runs that a strict consecutive-run
%   check nonetheless reported as "not found". A sliding window with a
%   small tolerance for gaps (k out of W, not W out of W) is deliberately
%   more robust to exactly that failure mode, at the cost of being a
%   softer criterion -- this tradeoff is a deliberate design response to
%   that empirical finding, not an arbitrary relaxation.
%
%   NOTE: branches(b).class(1) is always NaN (a branch's first
%   appearance has no previous-order transition to classify -- see
%   buildModalBranches.m). This function correctly ignores that first
%   NaN entry when forming windows, since it is not a transition at all.

    nBranches = length(branches);
    persistence = repmat(struct('isPersistent', false, ...
                                 'bestWindowStart', NaN, ...
                                 'bestWindowCount', 0), nBranches, 1);

    for b = 1:nBranches
        classHistory = branches(b).class(2:end);   % drop the leading NaN (not a real transition)
        orderHistory = branches(b).order(2:end);    % order at which each transition landed

        nTransitions = length(classHistory);
        if nTransitions < windowSize
            continue;   % too short to ever form a full window; stays not persistent
        end

        bestCount = -1;
        bestStart = NaN;
        for startIdx = 1:(nTransitions - windowSize + 1)
            windowClasses = classHistory(startIdx : startIdx+windowSize-1);
            countInWindow = sum(windowClasses == 3);
            if countInWindow > bestCount
                bestCount = countInWindow;
                bestStart = orderHistory(startIdx);
            end
        end

        persistence(b).bestWindowStart = bestStart;
        persistence(b).bestWindowCount = bestCount;
        persistence(b).isPersistent = bestCount >= minStableInWindow;
    end
end