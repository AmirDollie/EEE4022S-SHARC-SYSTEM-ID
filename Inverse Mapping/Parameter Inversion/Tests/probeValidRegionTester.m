%% probeValidRegionTester.m
% Validation tests for probeValidRegion.m, using a cheap analytic toy
% forward map with a KNOWN valid region (never the EMM), so every expected
% number below can be worked out by hand and the file runs in well under
% a second.
%
% TOY VALID REGION (in scaled units x = p./pRef), chosen so that the axis
% probes succeed individually but a combined corner fails, forcing a
% shrink round:
%     0.93 <= x1 <= 1.12,   0.85 <= x2 <= 1.30,   0.97 <= x3 <= 1.04,
%     (x1 - 1) + (x2 - 1) <= 0.12          (coupling constraint)
% Outside it the toy map throws 'toyModel:outsideRegion', standing in for
% modalFeatureVector:invalidFeatureFit under 'Strict',true.
%
% HAND-DERIVED EXPECTATIONS (default perturbations 1, 2, 5, 10, 20%):
%   axis:   x1 verified -5% .. +10%  (first failures -10% / +20%)
%           x2 verified -10% .. +10% (first failures -20% / +20%)
%           x3 verified -2% .. +2%   (first failures -5% / +5%)
%   round 1 box (BoundFraction 0.8): x1 [0.96, 1.08], x2 [0.92, 1.08],
%           x3 [0.984, 1.016]. Corners with x1 = x2 = 1.08 violate the
%           coupling (0.16 > 0.12): 2 of 8 fail.
%   round 2 box (ShrinkFactor 0.5): x1 [0.98, 1.04], x2 [0.96, 1.04],
%           x3 [0.992, 1.008]: all 8 corners valid -> suggested bounds.
%   evaluations: 1 reference + 25 axis (5 skipped after failures)
%                + 16 corners = 42.
%
% TEST LIST
%   1. Input validation, including a malformed forward output.
%   2. Axis probes, candidate box, corner shrink and suggested bounds
%      match the hand-derived values exactly; record bookkeeping.
%   3. StopAtFirstFailure = false probes every delta (no skips), same box.
%   4. Error classification: reference failure is always re-thrown (even
%      with a declared-expected ID); an unexpected ID mid-probe is
%      re-thrown; expected failures are recorded with their identifier.
%   5. relShift = (f - f0)./f0 on valid points.
%   6. CaptureDetails: anchor-bracket margins computed from
%      omega3/omegaStar, near-edge points counted; off by default.
%   7. Degenerate outcomes: no verified neighbourhood for a parameter
%      ('noAxisNeighbourhood'); corners never all valid ('cornersFailed').
%   8. Non-monotone validity is detected (StopAtFirstFailure = false), and
%      complex/non-finite forward output is recorded as invalid with its
%      own identifier instead of aborting the run.
%
% Lives in Inverse Mapping/Parameter Inversion/Tests/; probeValidRegion.m
% is one level up in Inverse Mapping/Parameter Inversion/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));   % probeValidRegion.m

clearvars -except thisDir
clc

fprintf('================================================================\n');
fprintf(' probeValidRegionTester.m\n');
fprintf('================================================================\n\n');

pRef = [2; 3; 5];
ID = 'toyModel:outsideRegion';
box = [0.93, 1.12; 0.85, 1.30; 0.97, 1.04];   % per-parameter [lo, hi]
fwd = @(p) regionMap(p, pRef, box, 0.12, ID);
common = {'ExpectedErrorIDs', {ID}, 'ParamNames', {'b', 'g', 'R'}, 'Verbose', false};
tol = 1e-12;

%% TEST 1: input validation
fprintf('--- TEST 1: input validation -------------------------------\n');

cases = {
    'forwardFcn not a handle',   @() probeValidRegion('regionMap', pRef),                               'probeValidRegion:badForwardFcn'
    'pRef has a zero',           @() probeValidRegion(fwd, [2; 0; 5]),                                  'probeValidRegion:badPRef'
    'perturbation of 0',         @() probeValidRegion(fwd, pRef, 'Perturbations', [0 0.1]),             'probeValidRegion:badOption'
    'perturbation of 1',         @() probeValidRegion(fwd, pRef, 'Perturbations', [0.1 1]),             'probeValidRegion:badOption'
    'empty perturbations',       @() probeValidRegion(fwd, pRef, 'Perturbations', []),                  'probeValidRegion:badOption'
    'unknown option',            @() probeValidRegion(fwd, pRef, 'Strict', true),                       'probeValidRegion:badOption'
    'odd option count',          @() probeValidRegion(fwd, pRef, 'Verbose'),                            'probeValidRegion:badOption'
    'bad ExpectedErrorIDs',      @() probeValidRegion(fwd, pRef, 'ExpectedErrorIDs', {1}),              'probeValidRegion:badOption'
    'ParamNames wrong count',    @() probeValidRegion(fwd, pRef, 'ParamNames', {'a', 'b'}),             'probeValidRegion:badOption'
    'BoundFraction 0',           @() probeValidRegion(fwd, pRef, 'BoundFraction', 0),                   'probeValidRegion:badOption'
    'BoundFraction > 1',         @() probeValidRegion(fwd, pRef, 'BoundFraction', 1.5),                 'probeValidRegion:badOption'
    'ShrinkFactor 1',            @() probeValidRegion(fwd, pRef, 'ShrinkFactor', 1),                    'probeValidRegion:badOption'
    'MaxShrinks non-integer',    @() probeValidRegion(fwd, pRef, 'MaxShrinks', 1.5),                    'probeValidRegion:badOption'
    'forward output is text',    @() probeValidRegion(@(p) 'abc', pRef, 'Verbose', false),              'probeValidRegion:badForwardOutput'
    'forward output changes size', @() probeValidRegion(@(p) ones(4 + (p(1) ~= pRef(1)), 1), pRef, 'Verbose', false), 'probeValidRegion:badForwardOutput'
};
for c = 1:size(cases, 1)
    gotId = '<no error>';
    try
        cases{c, 2}();
    catch err
        gotId = err.identifier;
    end
    assert(strcmp(gotId, cases{c, 3}), sprintf('Test1:case%d', c), ...
        '%s: expected %s, got %s', cases{c, 1}, cases{c, 3}, gotId);
end
fprintf('  (1) all %d invalid-input cases throw their documented identifier: PASS\n\n', size(cases, 1));

%% TEST 2: full probe on the toy region (printed, so the summary format is visible)
fprintf('--- TEST 2: axis probes, corner shrink, suggested bounds --\n');

probe = probeValidRegion(fwd, pRef, 'ExpectedErrorIDs', {ID}, 'ParamNames', {'b', 'g', 'R'});

assert(max(abs([probe.axis.verifiedLow] - [0.05 0.10 0.02])) < tol && ...
       max(abs([probe.axis.verifiedHigh] - [0.10 0.10 0.02])) < tol, 'Test2:verified', ...
    'verified extents %s / %s, expected [0.05 0.1 0.02] / [0.1 0.1 0.02]', ...
    mat2str([probe.axis.verifiedLow]), mat2str([probe.axis.verifiedHigh]));
assert(max(abs([probe.axis.firstFailLow] - [0.10 0.20 0.05])) < tol && ...
       max(abs([probe.axis.firstFailHigh] - [0.20 0.20 0.05])) < tol, 'Test2:firstFail', ...
    'first failures %s / %s, expected [0.1 0.2 0.05] / [0.2 0.2 0.05]', ...
    mat2str([probe.axis.firstFailLow]), mat2str([probe.axis.firstFailHigh]));
fprintf('  (2a) axis verified extents and first failures match hand derivation: PASS\n');

assert(numel(probe.combined) == 2 && probe.combined(1).nValid == 6 && probe.combined(2).nValid == 8, ...
    'Test2:rounds', 'expected 2 corner rounds with 6/8 then 8/8 valid');
r1 = probe.points(strcmp({probe.points.stage}, 'combined') & [probe.points.round] == 1);
bad = r1(strcmp({r1.status}, 'invalid'));
badX = [bad.x];
assert(numel(bad) == 2 && all(abs(badX(1, :) - 1.08) < tol) && all(abs(badX(2, :) - 1.08) < tol), ...
    'Test2:failedCorners', 'the failing round-1 corners should be exactly those with x1 = x2 = 1.08');
fprintf('  (2b) round 1: exactly the two x1 = x2 = +8%% corners fail; round 2: all 8 valid: PASS\n');

assert(strcmp(probe.boxStatus, 'verified'), 'Test2:status', 'boxStatus = %s, expected verified', probe.boxStatus);
assert(max(abs(probe.suggestedLB - [0.98; 0.96; 0.992])) < tol && ...
       max(abs(probe.suggestedUB - [1.04; 1.04; 1.008])) < tol, 'Test2:bounds', ...
    'suggested bounds %s .. %s, expected [0.98 0.96 0.992] .. [1.04 1.04 1.008]', ...
    mat2str(probe.suggestedLB.'), mat2str(probe.suggestedUB.'));
assert(max(abs(probe.suggestedLBp - probe.suggestedLB .* pRef)) < tol && ...
       max(abs(probe.suggestedUBp - probe.suggestedUB .* pRef)) < tol, 'Test2:boundsP', ...
    'physical-unit bounds are not suggested*.pRef');
fprintf('  (2c) suggested bounds [0.98 0.96 0.992] .. [1.04 1.04 1.008], physical units consistent: PASS\n');

st = {probe.points.status};
assert(probe.nEvaluations == 42 && sum(strcmp(st, 'skipped')) == 5 && numel(probe.points) == 47, ...
    'Test2:counts', 'expected 42 evaluations, 5 skipped, 47 records; got %d, %d, %d', ...
    probe.nEvaluations, sum(strcmp(st, 'skipped')), numel(probe.points));
assert(strcmp(probe.points(1).stage, 'reference') && isequal(probe.reference.x, [1; 1; 1]), ...
    'Test2:reference', 'first record is not the reference point');
pts = probe.points;
for k = 1:numel(pts)
    assert(max(abs(pts(k).p - pts(k).x .* pRef)) < tol, 'Test2:pConsistency', 'record %d: p ~= x.*pRef', k);
    switch pts(k).status
        case 'valid'
            assert(numel(pts(k).f) == 5 && isempty(pts(k).errorIdentifier), 'Test2:validRecord', ...
                'record %d: valid but f/error fields inconsistent', k);
        case 'invalid'
            assert(isempty(pts(k).f) && strcmp(pts(k).errorIdentifier, ID) && ~isempty(pts(k).errorMessage), ...
                'Test2:invalidRecord', 'record %d: invalid but error fields not recorded', k);
        case 'skipped'
            assert(isempty(pts(k).f) && isnan(pts(k).runtimeSec), 'Test2:skippedRecord', ...
                'record %d: skipped point appears to have been evaluated', k);
    end
end
fprintf('  (2d) 42 evaluations, 5 skipped; every record internally consistent: PASS\n\n');

%% TEST 3: StopAtFirstFailure = false
fprintf('--- TEST 3: StopAtFirstFailure = false ---------------------\n');

probeAll = probeValidRegion(fwd, pRef, common{:}, 'StopAtFirstFailure', false);
assert(probeAll.nEvaluations == 47 && ~any(strcmp({probeAll.points.status}, 'skipped')), 'Test3:noSkips', ...
    'expected 47 evaluations and no skips, got %d', probeAll.nEvaluations);
assert(isequal(probeAll.suggestedLB, probe.suggestedLB) && isequal(probeAll.suggestedUB, probe.suggestedUB), ...
    'Test3:sameBox', 'probing every delta changed the suggested box for a monotone region');
fprintf('  (3) every delta probed (47 evaluations, no skips), identical suggested box: PASS\n\n');

%% TEST 4: error classification
fprintf('--- TEST 4: error classification ---------------------------\n');

% (4a) Reference outside the region: re-thrown even though the ID is expected.
fwdRefBad = @(p) regionMap(p, pRef, [1.01, 1.2; 0.85, 1.3; 0.97, 1.04], 0.12, ID);
gotId = '<no error>';
try
    probeValidRegion(fwdRefBad, pRef, common{:});
catch err
    gotId = err.identifier;
end
assert(strcmp(gotId, ID), 'Test4a:reference', 'a failing reference was not re-thrown (got %s)', gotId);
fprintf('  (4a) failure at the reference point is re-thrown, not recorded: PASS\n');

% (4b) Unexpected identifier mid-probe: re-thrown.
fwdBug = @(p) regionMap(p, pRef, box, 0.12, 'toyModel:genuineBug');
gotId = '<no error>';
try
    probeValidRegion(fwdBug, pRef, common{:});
catch err
    gotId = err.identifier;
end
assert(strcmp(gotId, 'toyModel:genuineBug'), 'Test4b:unexpected', ...
    'an unexpected identifier was not re-thrown (got %s)', gotId);
fprintf('  (4b) an unexpected error identifier mid-probe is re-thrown: PASS\n');

% (4c) Expected failures recorded (already checked per record in 2d);
%      also accepted as a single char ID.
probeChar = probeValidRegion(fwd, pRef, 'ExpectedErrorIDs', ID, 'Verbose', false);
assert(strcmp(probeChar.boxStatus, 'verified'), 'Test4c:charID', 'char ExpectedErrorIDs not accepted');
fprintf('  (4c) expected failures recorded; ExpectedErrorIDs accepted as a char string: PASS\n\n');

%% TEST 5: relShift
fprintf('--- TEST 5: relative feature shifts ------------------------\n');

f0 = toyMap(pRef);
k = find(strcmp({probe.points.stage}, 'axis') & [probe.points.paramIndex] == 1 & ...
    abs([probe.points.delta] - 0.05) < tol, 1);
expected = (toyMap(pRef .* [1.05; 1; 1]) - f0) ./ f0;
assert(max(abs(probe.points(k).relShift - expected)) < tol, 'Test5:relShift', 'relShift ~= (f - f0)./f0');
assert(isequal(probe.reference.relShift, zeros(5, 1)), 'Test5:refShift', 'reference relShift is not zero');
fprintf('  (5) relShift = (f - f0)./f0 at x1 = +5%%, zero at the reference: PASS\n\n');

%% TEST 6: CaptureDetails and anchor-bracket margins
fprintf('--- TEST 6: CaptureDetails / bracket margins ---------------\n');

% Toy details: every feature has anchors [1 2 3] and vertex 2 + 9*(x1 - 1),
% so margin = min(w* - 1, 3 - w*)/2: 0.5 at x1 = 1, 0.275 at +5%, 0.05 at
% +10% (the only valid point below the 0.1 warning level).
pd = probeValidRegion(fwd, pRef, common{:}, 'CaptureDetails', true);
assert(abs(pd.reference.minBracketMargin - 0.5) < tol, 'Test6:refMargin', ...
    'reference margin %.4f, expected 0.5', pd.reference.minBracketMargin);
k5 = find(strcmp({pd.points.stage}, 'axis') & [pd.points.paramIndex] == 1 & abs([pd.points.delta] - 0.05) < tol, 1);
k10 = find(strcmp({pd.points.stage}, 'axis') & [pd.points.paramIndex] == 1 & abs([pd.points.delta] - 0.10) < tol, 1);
assert(abs(pd.points(k5).minBracketMargin - 0.275) < tol && abs(pd.points(k10).minBracketMargin - 0.05) < tol, ...
    'Test6:margins', 'margins at +5%%/+10%% are %.4f/%.4f, expected 0.275/0.05', ...
    pd.points(k5).minBracketMargin, pd.points(k10).minBracketMargin);
assert(numel(pd.points(k5).bracketMargin) == 5 && numel(pd.points(k5).details) == 5, 'Test6:perFeature', ...
    'details/margins not stored per feature');
assert(pd.nNearEdge == 1, 'Test6:nearEdge', 'nNearEdge = %d, expected 1', pd.nNearEdge);
assert(isnan(probe.reference.minBracketMargin) && isempty(probe.reference.details) && probe.nNearEdge == 0, ...
    'Test6:offByDefault', 'details captured although CaptureDetails was not set');
fprintf('  (6) bracket margins 0.5 / 0.275 / 0.05 as derived; 1 near-edge point; off by default: PASS\n\n');

%% TEST 7: degenerate outcomes
fprintf('--- TEST 7: degenerate outcomes ----------------------------\n');

% (7a) x3 valid only within [0.995, 1.004]: the smallest (1%) probe fails
%      both ways -> no box can be formed, no corner test is run.
fwdTight = @(p) regionMap(p, pRef, [0.93, 1.12; 0.85, 1.3; 0.995, 1.004], 0.12, ID);
pt = probeValidRegion(fwdTight, pRef, common{:});
assert(strcmp(pt.boxStatus, 'noAxisNeighbourhood') && all(isnan(pt.suggestedLB)) && isempty(pt.combined), ...
    'Test7a:noNeighbourhood', 'expected noAxisNeighbourhood, NaN bounds and no corner rounds');
assert(pt.axis(3).verifiedLow == 0 && pt.axis(3).verifiedHigh == 0, 'Test7a:zeroExtent', ...
    'x3 verified extent should be 0 in both directions');
fprintf('  (7a) smallest probe failing both ways -> noAxisNeighbourhood, NaN bounds: PASS\n');

% (7b) No shrinking allowed and round 1 fails -> cornersFailed.
pc = probeValidRegion(fwd, pRef, common{:}, 'MaxShrinks', 0);
assert(strcmp(pc.boxStatus, 'cornersFailed') && all(isnan(pc.suggestedUB)) && numel(pc.combined) == 1, ...
    'Test7b:cornersFailed', 'expected cornersFailed with NaN bounds after a single round');
fprintf('  (7b) corners never all valid -> cornersFailed, NaN bounds: PASS\n\n');

%% TEST 8: non-monotone validity and non-finite output
fprintf('--- TEST 8: non-monotone / non-finite -----------------------\n');

% holeMap: x1 in [1.07, 1.15] throws (a "hole"), x1 < 0.85 returns NaN;
% x2, x3 unrestricted. So along x1: +10% fails but +20% is valid again,
% and -20% runs but produces a NaN feature.
fwdHole = @(p) holeMap(p, pRef, ID);
ph = probeValidRegion(fwdHole, pRef, common{:}, 'StopAtFirstFailure', false);
a1 = ph.axis(1);
assert(abs(a1.verifiedHigh - 0.05) < tol && abs(a1.firstFailHigh - 0.10) < tol && a1.nonMonotoneHigh, ...
    'Test8a:nonMonotone', 'expected x1 verified +5%%, first failure +10%%, non-monotone on the + side');
assert(~a1.nonMonotoneLow && ~any([ph.axis(2:3).nonMonotoneLow, ph.axis(2:3).nonMonotoneHigh]), ...
    'Test8a:falsePositive', 'non-monotone flagged where validity is monotone');
fprintf('  (8a) +10%% fails, +20%% valid again -> nonMonotoneHigh flagged, verified extent stays +5%%: PASS\n');

phStop = probeValidRegion(fwdHole, pRef, common{:});
assert(~phStop.axis(1).nonMonotoneHigh, 'Test8b:stopHides', ...
    'with early stopping the +20% point is skipped, so non-monotonicity cannot be (and must not be) claimed');
fprintf('  (8b) with early stopping (default) the flag stays false, as documented: PASS\n');

k = find(strcmp({ph.points.stage}, 'axis') & [ph.points.paramIndex] == 1 & abs([ph.points.delta] + 0.20) < tol, 1);
assert(strcmp(ph.points(k).status, 'invalid') && ...
    strcmp(ph.points(k).errorIdentifier, 'probeValidRegion:nonFiniteOutput'), 'Test8c:nonFinite', ...
    'NaN forward output at x1 = -20%% was not recorded as invalid with probeValidRegion:nonFiniteOutput');
assert(abs(a1.verifiedLow - 0.10) < tol && abs(a1.firstFailLow - 0.20) < tol, 'Test8c:extent', ...
    'the non-finite point should end the verified extent at -10%%');
fprintf('  (8c) NaN output recorded as invalid (probeValidRegion:nonFiniteOutput), run not aborted: PASS\n\n');

fprintf('================================================================\n');
fprintf(' ALL probeValidRegionTester.m TESTS PASSED\n');
fprintf('================================================================\n');

%% ------------------------------------------------------------------------
function f = toyMap(p)
% 3 parameters -> 5 features (same toy map as the other Parameter
% Inversion testers).
f = [p(1) + p(2); p(1) * p(3); p(2)^2; p(1) + p(3); 2 * p(2)];
end

function f = holeMap(p, pRef, id)
% toyMap with a hole in validity at x1 in [1.07, 1.15], and a NaN
% feature (not an error) for x1 < 0.85.
x1 = p(1) / pRef(1);
if x1 >= 1.07 && x1 <= 1.15
    error(id, 'toy forward model: x1 = %.3f is in the invalid hole', x1);
end
f = toyMap(p);
if x1 < 0.85
    f(3) = NaN;
end
end

function [f, det] = regionMap(p, pRef, box, coupling, id)
% toyMap inside a known valid region, error `id` outside it. Also returns
% modalFeatureVector-style details (omega3, omegaStar) for 5 features.
x = p(:) ./ pRef(:);
inside = all(x >= box(:, 1) - 1e-12 & x <= box(:, 2) + 1e-12) && (x(1) - 1) + (x(2) - 1) <= coupling + 1e-12;
if ~inside
    error(id, 'toy forward model: x = %s is outside the valid region', mat2str(x.', 5));
end
f = toyMap(p);
det = repmat(struct('omega3', [1 2 3], 'omegaStar', 2 + 9 * (x(1) - 1)), 1, 5);
end