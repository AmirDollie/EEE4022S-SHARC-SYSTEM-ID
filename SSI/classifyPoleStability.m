function match = classifyPoleStability(prevResults, currResults, freqTol, dampTol, macTol)
%CLASSIFYPOLESTABILITY Matches poles at model order n against order n-1
%and assigns a hierarchical stability class (Peeters & De Roeck 1999,
%Section 8; classification hierarchy per project discussion).
%
%   usage MATCH = CLASSIFYPOLESTABILITY(PREVRESULTS, CURRRESULTS, FREQTOL, DAMPTOL, MACTOL)
%
%   PREVRESULTS, CURRRESULTS - single struct entries from sweepModelOrders
%       output, i.e. PREVRESULTS = results(n-1), CURRRESULTS = results(n).
%   FREQTOL - relative frequency tolerance (e.g. 0.01 for 1%)
%   DAMPTOL - relative damping tolerance (e.g. 0.05 for 5%)
%   MACTOL  - MAC threshold (e.g. 0.99)
%
%   Returns MATCH, a struct array of length k = length(currResults.frequency),
%   one entry per pole at the CURRENT order, with fields:
%   previousIndex   - index into prevResults arrays this pole matched to,
%                     or 0 if unmatched
%   frequencyStable - true if the matched pair passed the frequency gate
%   dampingStable   - true if it also passed the damping gate
%   macStable       - true if it also passed the MAC gate
%   class           - 0 (no match), 1 (frequency only), 2 (frequency +
%                     damping), or 3 (frequency + damping + MAC, "fully
%                     stable")
%
%   MATCHING PROCEDURE: frequency is a hard gate (candidates outside
%   freqTol are never considered at all, regardless of damping/MAC).
%   Among candidates passing the frequency gate, damping and MAC are
%   evaluated but do NOT further restrict candidacy for matching itself
%   -- a candidate can still be matched (as class 1 or 2) even if damping
%   or MAC fail, since a partial match is still useful diagnostic
%   information (per the project discussion: "don't throw away partially
%   stable poles"). Assignment itself is a GREEDY one-to-one match: among
%   all (previous, current) pairs passing the frequency gate, repeatedly
%   pick the single best-MAC pair, remove both indices from further
%   consideration, and repeat. This guarantees no previous-order pole is
%   claimed by two different current-order poles (or vice versa).
%
%   Ranking by MAC (rather than by frequency closeness) for the greedy
%   selection is deliberate: MAC is the strictest, most specific
%   similarity measure available, so it is the best tie-breaker when
%   multiple candidates pass the frequency gate.

    kCurr = length(currResults.frequency);
    kPrev = length(prevResults.frequency);

    match = repmat(struct('previousIndex', 0, 'frequencyStable', false, ...
                           'dampingStable', false, 'macStable', false, ...
                           'class', 0), kCurr, 1);

    if kPrev == 0 || kCurr == 0
        return;   % nothing to match against; all entries stay class 0
    end

    % Build candidate table: for every (prev,curr) pair passing the
    % frequency gate, compute damping/MAC and whether each passes.
    freqOK = false(kPrev, kCurr);
    dampOK = false(kPrev, kCurr);
    macVal = zeros(kPrev, kCurr);

    for p = 1:kPrev
        for c = 1:kCurr
            fRel = abs(currResults.frequency(c) - prevResults.frequency(p)) / prevResults.frequency(p);
            if fRel < freqTol
                freqOK(p,c) = true;
                zRel = abs(currResults.damping(c) - prevResults.damping(p)) / abs(prevResults.damping(p));
                dampOK(p,c) = zRel < dampTol;
                macVal(p,c) = computeMAC(prevResults.modeShapes(:,p), currResults.modeShapes(:,c));
            end
        end
    end

    % Greedy one-to-one assignment among frequency-gated candidates,
    % ranked by MAC (highest first).
    availablePrev = true(kPrev, 1);
    availableCurr = true(kCurr, 1);

    candidateList = [];
    for p = 1:kPrev
        for c = 1:kCurr
            if freqOK(p,c)
                candidateList = [candidateList; p, c, macVal(p,c)]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(candidateList)
        [~, order] = sort(candidateList(:,3), 'descend');
        candidateList = candidateList(order, :);

        for row = 1:size(candidateList,1)
            p = candidateList(row,1);
            c = candidateList(row,2);
            if availablePrev(p) && availableCurr(c)
                availablePrev(p) = false;
                availableCurr(c) = false;

                match(c).previousIndex = p;
                match(c).frequencyStable = true;
                match(c).dampingStable = dampOK(p,c);
                match(c).macStable = macVal(p,c) >= macTol;

                if match(c).macStable
                    match(c).class = 3;
                elseif match(c).dampingStable
                    match(c).class = 2;
                else
                    match(c).class = 1;
                end
            end
        end
    end
end