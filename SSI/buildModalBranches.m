function branches = buildModalBranches(results, matches)
%BUILDMODALBRANCHES Reconstructs continuous modal "branches" (a single
%pole tracked across many consecutive model orders) from the per-order
%results and per-transition matches already computed by
%sweepModelOrders.m and classifyPoleStability.m.
%
%   BRANCHES = BUILDMODALBRANCHES(RESULTS, MATCHES)
%
%   RESULTS - 1 x nMax struct array from sweepModelOrders.m
%   MATCHES - 1 x nMax cell array; MATCHES{n} is the output of
%             classifyPoleStability(results(n-1), results(n), ...) for
%             n=2..nMax. MATCHES{1} is unused (no previous order to
%             compare against).
%
%   Returns BRANCHES, a struct array. Each branches(b) represents one
%   continuous chain and has fields:
%   order      - row vector of model orders this branch appears at
%   frequency  - matching row vector of frequencies at those orders
%   damping    - matching row vector of damping ratios
%   class      - matching row vector of the classifyPoleStability class
%                that linked this branch INTO that order (class(1) is
%                NaN, since the branch's first appearance has no
%                previous-order link by definition)
%
%   WHY THIS IS ITS OWN FUNCTION: classifyPoleStability only ever looks
%   at ONE pair of adjacent orders at a time -- it has no notion of a
%   pole's history beyond the immediately preceding order. Persistence
%   (how long has this SAME pole been showing up) can only be assessed
%   once these order-by-order links are threaded together into actual
%   chains, which is exactly what this function does and nothing more --
%   deciding what counts as "persistent enough" is left to
%   computeBranchPersistence.m.

    nMax = length(results);

    branches = struct('order', {}, 'frequency', {}, 'damping', {}, 'class', {});
    activeBranchOfPole = [];   % activeBranchOfPole(p) = branch index currently ending at pole p of the previous order

    for n = 1:nMax
        k = length(results(n).frequency);
        newActiveBranchOfPole = zeros(k,1);

        for c = 1:k
            linkedBranch = 0;
            linkedClass = NaN;

            if n > 1
                prevIdx = matches{n}(c).previousIndex;
                if prevIdx > 0 && prevIdx <= length(activeBranchOfPole) && activeBranchOfPole(prevIdx) > 0
                    linkedBranch = activeBranchOfPole(prevIdx);
                    linkedClass = matches{n}(c).class;
                end
            end

            if linkedBranch > 0
                b = linkedBranch;
                branches(b).order(end+1)     = n;
                branches(b).frequency(end+1) = results(n).frequency(c);
                branches(b).damping(end+1)   = results(n).damping(c);
                branches(b).class(end+1)     = linkedClass;
            else
                branches(end+1) = struct( ...
                    'order', n, ...
                    'frequency', results(n).frequency(c), ...
                    'damping', results(n).damping(c), ...
                    'class', NaN); %#ok<AGROW>
                b = length(branches);
            end

            newActiveBranchOfPole(c) = b;
        end

        activeBranchOfPole = newActiveBranchOfPole;
    end
end