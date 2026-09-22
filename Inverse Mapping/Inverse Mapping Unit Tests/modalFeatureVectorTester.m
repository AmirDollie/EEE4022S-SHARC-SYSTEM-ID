%% modalFeatureVectorTester.m
% Validation tests for modalFeatureVector.m and its three helper
% functions: quadFit3.m (exact 3-point quadratic fit), loadFeatureAnchors.m
% (reads the fixed sub-grid interpolation anchors from a
% modalFeatureRefinement.m results file), and fitFeatureFrom3Points.m (the
% pure degenerate/extrapolated/kindMismatch decision logic).
%
% TEST LIST
%   1. quadFit3.m                     -- exact recovery, degenerate case,
%                                         error handling.
%   2. fitFeatureFrom3Points.m        -- direct, EMM-free tests of the
%                                         degenerate/extrapolated/
%                                         kindMismatch flags and the
%                                         omegaStar fallback, using known
%                                         synthetic (omega3,y3,kind)
%                                         triples.
%   3. loadFeatureAnchors.m           -- bracket extraction, constant
%                                         pass-through, and every error
%                                         path (missing file, malformed
%                                         file, edge-of-window feature,
%                                         featureOmega/grid mismatch,
%                                         wrong label order).
%   4. modalFeatureVector.m           -- end-to-end smoke test against the
%                                         real Forward Model at cheap
%                                         sanity parameters, via a full
%                                         5-feature synthetic anchor file.
%   5. 'Strict' option                -- confirms Strict=true throws
%                                         modalFeatureVector:invalidFeatureFit
%                                         whenever Strict=false would have
%                                         only warned, and otherwise
%                                         returns the identical fVec.
%   6. Integer validation             -- M, P, N, NTheta reject
%                                         non-integer values.
%   7. Override-file independence     -- two different
%                                         'RefinementResultsFile' values,
%                                         called back to back, each
%                                         resolve to their own correct
%                                         anchors (no cross-call
%                                         contamination).
%
% SAFETY NOTE: this tester never reads or writes the real
% modalFeatureRefinementResults.mat next to modalFeatureVector.m. Every
% test here builds its own small synthetic results file in a scratch
% folder (via makeSyntheticFeatureRefinement.m) and points
% modalFeatureVector.m/loadFeatureAnchors.m at it via the
% 'RefinementResultsFile' override -- so this is safe to run in the real
% repo even before (or after) modalFeatureRefinement.m has been run for
% real, and it will never clobber real results. Because of that, this
% tester does NOT exercise the default-path behaviour (the case with no
% 'RefinementResultsFile' given, which reads
% modalFeatureRefinementResults.mat next to modalFeatureVector.m) --
% constructing a fake file at that real location, even temporarily, risks
% colliding with or masking a real results file if one is present. That
% path is exercised naturally the first time modalFeatureVector.m is
% actually used against the real refinement output (see the "production
% nominal call" note at the end of this file's companion delivery notes).
%
% NOTE ON modalFeatureVector.m's anchor loading: as of the current
% revision, loadFeatureAnchors.m is always called fresh on every call (no
% persistent cache), so "staleness" in the old sense is no longer
% possible by construction. TEST 7 below is kept as a regression check
% that two different override files, called back to back, still resolve
% independently and correctly -- i.e. that removing the cache did not
% introduce some other form of cross-call state leakage.
%
% NOTE ON TEST 5 ('Strict'): this test exercises the real Forward
% Model/projection pipeline (unlike TEST 2, which tests the
% degenerate/extrapolated/kindMismatch DECISION LOGIC directly with no
% EMM dependency at all). The five synthetic anchor windows used here are
% NOT the real, validated feature locations for the cheap sanity
% parameters used elsewhere in this file -- they are arbitrary,
% uninformed omega windows -- so which, if any, of the three flags
% actually fires is a genuine property of the forward model at those
% frequencies and cannot be hardcoded in advance. TEST 5 is written to
% adapt to whichever case actually occurs on a given run (see its
% comments) rather than assume one outcome.
%
% This file lives in Inverse Mapping/Inverse Mapping Unit Tests/,
% alongside modalRAOReconnaissanceTester.m, per this repo's convention of
% keeping testers in a Unit-Tests-style subfolder (mirroring
% Forward Model/Unit-Tests/) rather than flat next to the source they
% test. modalFeatureVector.m, loadFeatureAnchors.m, quadFit3.m and
% fitFeatureFrom3Points.m themselves stay one level up in
% Inverse Mapping/, so (unlike modalRAOReconnaissanceTester.m, which only
% ever loads a sibling .mat file and never calls modalRAOReconnaissance.m
% as a function) this tester also needs Inverse Mapping/ itself on the
% path, since it calls those functions directly -- see the extra addpath
% below. makeSyntheticFeatureRefinement.m is TEST-ONLY, so it lives here
% instead, alongside this tester, not in Inverse Mapping/.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);                                              % makeSyntheticFeatureRefinement.m
addpath(fullfile(thisDir, '..'));                              % modalFeatureVector.m, loadFeatureAnchors.m, quadFit3.m, fitFeatureFrom3Points.m
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));

clearvars -except thisDir
close all, clc

scratchDir = tempname;
mkdir(scratchDir);
cleanupObj = onCleanup(@() rmdir(scratchDir, 's'));

fprintf('================================================================\n');
fprintf(' modalFeatureVectorTester.m\n');
fprintf(' scratch folder: %s\n', scratchDir);
fprintf('================================================================\n\n');

%% TEST 1: quadFit3.m -- exact recovery, degenerate case, error handling
fprintf('--- TEST 1: quadFit3.m -------------------------------------\n');

% (1a) Known parabola, UNEQUALLY spaced points: y = 2*(x-5)^2 + 3
%      => a=2, b=-20, c=53, vertex at x=5, y=3.
xTest = [4.7, 4.9, 5.4];
yTest = 2*(xTest - 5).^2 + 3;
[aFit, bFit, cFit] = quadFit3(xTest, yTest);
assert(abs(aFit - 2)  < 1e-9, 'Test1a:a',  'quadFit3 recovered a=%.10f, expected 2', aFit);
assert(abs(bFit - -20) < 1e-8, 'Test1a:b',  'quadFit3 recovered b=%.10f, expected -20', bFit);
assert(abs(cFit - 53) < 1e-7, 'Test1a:c',  'quadFit3 recovered c=%.10f, expected 53', cFit);
vertexFit = -bFit / (2*aFit);
assert(abs(vertexFit - 5) < 1e-9, 'Test1a:vertex', 'quadFit3 vertex=%.10f, expected 5', vertexFit);
fprintf('  (1a) exact recovery of a known parabola (unequal spacing): PASS (a=%.6f, b=%.6f, c=%.6f, vertex=%.6f)\n', aFit, bFit, cFit, vertexFit);

% (1b) Flat/degenerate points: a must come back exactly 0.
[aFlat, bFlat, cFlat] = quadFit3([1, 2, 3], [5, 5, 5]);
assert(aFlat == 0, 'Test1b:a', 'quadFit3 on flat points gave a=%.3e, expected exactly 0', aFlat);
assert(bFlat == 0, 'Test1b:b', 'quadFit3 on flat points gave b=%.3e, expected exactly 0', bFlat);
assert(cFlat == 5, 'Test1b:c', 'quadFit3 on flat points gave c=%.10f, expected 5', cFlat);
fprintf('  (1b) flat/degenerate points give a=0 exactly: PASS\n');

% (1c) A parabola whose vertex falls WAY outside the sampled bracket
%      (this is the shape the 'extrapolated' flag is built to catch):
%      y=(x-100)^2 sampled at x=[1,2,3].
[aOut, bOut, ~] = quadFit3([1, 2, 3], (( [1,2,3] - 100).^2));
vertexOut = -bOut / (2*aOut);
assert(abs(vertexOut - 100) < 1e-6, 'Test1c:vertex', 'quadFit3 vertex=%.6f, expected 100', vertexOut);
assert(vertexOut < 1 || vertexOut > 3, 'Test1c:outside', 'expected vertex outside [1,3]');
fprintf('  (1c) far-outside-bracket vertex correctly recovered (%.4f, outside [1,3]): PASS\n', vertexOut);

% (1d) Error handling: repeated x, and wrong-sized input.
threw = false;
try
    quadFit3([2, 2, 3], [1, 2, 3]);
catch err
    threw = strcmp(err.identifier, 'quadFit3:repeatedX');
end
assert(threw, 'Test1d:repeatedX', 'quadFit3 did not throw quadFit3:repeatedX on repeated x');

threw = false;
try
    quadFit3([1, 2], [1, 2]);
catch err
    threw = strcmp(err.identifier, 'quadFit3:badInput');
end
assert(threw, 'Test1d:badInput', 'quadFit3 did not throw quadFit3:badInput on a 2-element input');
fprintf('  (1d) error handling (repeated x, wrong size): PASS\n\n');

%% TEST 2: fitFeatureFrom3Points.m -- direct, EMM-free flag tests
fprintf('--- TEST 2: fitFeatureFrom3Points.m -------------------------\n');

% (2a) Clean 'max' feature: y = -(x-5)^2 + 10, sampled at x=[4.7,5.0,5.6].
%      True peak at x=5, well inside the bracket, clear curvature (a=-1),
%      correctly downward-opening for a 'max' feature.
x2a = [4.7, 5.0, 5.6];
y2a = -(x2a - 5).^2 + 10;
[wStar, deg, extr, kMis, aC, ~, ~] = fitFeatureFrom3Points(x2a, y2a, 'max');
assert(abs(wStar - 5) < 1e-9, 'Test2a:omegaStar', 'expected omegaStar=5, got %.10f', wStar);
assert(~deg && ~extr && ~kMis, 'Test2a:flags', 'expected all flags false for a clean max feature, got degenerate=%d extrapolated=%d kindMismatch=%d', deg, extr, kMis);
assert(aC < 0, 'Test2a:curvature', 'expected a<0 for a max feature, got a=%.4f', aC);
fprintf('  (2a) clean ''max'' feature: omegaStar=%.6f, all flags false: PASS\n', wStar);

% (2b) Clean 'zero' feature: y = (x-5)^2, same x -- a local minimum
%      (upward-opening), which is the correct shape for a near-zero
%      modal-coefficient feature.
y2b = (x2a - 5).^2;
[wStar, deg, extr, kMis, aC, ~, ~] = fitFeatureFrom3Points(x2a, y2b, 'zero');
assert(abs(wStar - 5) < 1e-9, 'Test2b:omegaStar', 'expected omegaStar=5, got %.10f', wStar);
assert(~deg && ~extr && ~kMis, 'Test2b:flags', 'expected all flags false for a clean zero feature, got degenerate=%d extrapolated=%d kindMismatch=%d', deg, extr, kMis);
assert(aC > 0, 'Test2b:curvature', 'expected a>0 for a zero feature, got a=%.4f', aC);
fprintf('  (2b) clean ''zero'' feature: omegaStar=%.6f, all flags false: PASS\n', wStar);

% (2c) Degenerate: exactly flat points. omegaStar must fall back to the
%      anchor with the extreme y-value (here, all equal -> first index).
[wStar, deg, extr, kMis, ~, ~, ~] = fitFeatureFrom3Points([1, 2, 3], [5, 5, 5], 'max');
assert(deg, 'Test2c:degenerate', 'expected degenerate=true for exactly flat points');
assert(~extr && ~kMis, 'Test2c:forcedFalse', 'extrapolated and kindMismatch must both be forced false when degenerate');
assert(wStar == 1, 'Test2c:fallback', 'expected fallback omegaStar=1 (first of three equal max candidates), got %.10f', wStar);
fprintf('  (2c) degenerate (flat) case: correctly detected, fallback omegaStar=%.4f: PASS\n', wStar);

% (2d) Extrapolated: true vertex far outside the anchor bracket.
%      y=(x-100)^2 at x=[1,2,3], kind='zero' (a>0 is the correct shape,
%      so this is a non-degenerate, non-kindMismatch, purely
%      out-of-bracket case).
[wStar, deg, extr, kMis, ~, ~, ~] = fitFeatureFrom3Points([1, 2, 3], ([1,2,3]-100).^2, 'zero');
assert(~deg, 'Test2d:notDegenerate', 'expected degenerate=false (real curvature present)');
assert(extr, 'Test2d:extrapolated', 'expected extrapolated=true for a vertex at 100 with bracket [1,3]');
assert(~kMis, 'Test2d:notMismatch', 'expected kindMismatch=false (a>0 is correct for a zero feature)');
assert(abs(wStar - 100) < 1e-6, 'Test2d:omegaStar', 'expected omegaStar=100, got %.6f', wStar);
fprintf('  (2d) extrapolated (vertex outside bracket) case: correctly detected, omegaStar=%.4f: PASS\n', wStar);

% (2e) Kind mismatch: y=(x-5)^2 (upward-opening, true local minimum) but
%      declared as a 'max' feature.
[wStar, deg, extr, kMis, aC, ~, ~] = fitFeatureFrom3Points(x2a, y2b, 'max');
assert(~deg && ~extr, 'Test2e:otherFlags', 'expected degenerate=false, extrapolated=false for this well-conditioned in-bracket case');
assert(kMis, 'Test2e:kindMismatch', 'expected kindMismatch=true: a=%.4f (>0) declared as a ''max'' feature', aC);
fprintf('  (2e) kind-mismatch (upward parabola declared ''max''): correctly detected: PASS\n');

% (2f) Error handling: invalid kind string, and wrong-sized input.
threw = false;
try
    fitFeatureFrom3Points([1, 2, 3], [1, 2, 3], 'bogus');
catch err
    threw = strcmp(err.identifier, 'fitFeatureFrom3Points:badKind');
end
assert(threw, 'Test2f:badKind', 'did not throw fitFeatureFrom3Points:badKind for an invalid kind string');

threw = false;
try
    fitFeatureFrom3Points([1, 2], [1, 2], 'max');
catch err
    threw = strcmp(err.identifier, 'fitFeatureFrom3Points:badInput');
end
assert(threw, 'Test2f:badInput', 'did not throw fitFeatureFrom3Points:badInput for a 2-element input');
fprintf('  (2f) error handling (bad kind, wrong size): PASS\n\n');

%% TEST 3: loadFeatureAnchors.m -- bracket extraction and error handling
fprintf('--- TEST 3: loadFeatureAnchors.m ----------------------------\n');

% (3a) Well-formed, correctly-ordered 5-feature file (default centres
%      [2.0 2.5 3.0 3.5 4.0], dOmega=0.1 -> anchors are centre +/- 0.1).
[refined, modeList, nu, M, P, N, H, g] = makeSyntheticFeatureRefinement(); %#ok<ASGLU>
goodFile = fullfile(scratchDir, 'good.mat');
save(goodFile, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g');

expectedAnchors = { [1.9, 2.0, 2.1], [2.4, 2.5, 2.6], [2.9, 3.0, 3.1], ...
    [3.4, 3.5, 3.6], [3.9, 4.0, 4.1] };
Anc = loadFeatureAnchors(goodFile);
assert(numel(Anc.features) == 5, 'Test3a:count', 'expected 5 features, got %d', numel(Anc.features));
for k = 1:5
    assert(isequal(Anc.features(k).omega3, expectedAnchors{k}), 'Test3a:bracket', ...
        'feature %d (%s): expected omega3=[%s], got [%s]', k, Anc.features(k).label, ...
        mat2str(expectedAnchors{k}), mat2str(Anc.features(k).omega3));
end
assert(Anc.M == 10 && Anc.P == 6 && Anc.N == 10 && Anc.nu == 0.3 && Anc.H == 1.88 && Anc.g == 9.81, ...
    'Test3a:passthrough', 'forward-model constants did not pass through loadFeatureAnchors correctly');
assert(isequal(Anc.modeList, modeList), 'Test3a:modeList', 'modeList did not pass through correctly');
fprintf('  (3a) correct bracket extraction (all 5 features) + constant pass-through: PASS\n');

% (3b) Missing file.
threw = false;
try
    loadFeatureAnchors(fullfile(scratchDir, 'does_not_exist.mat'));
catch err
    threw = strcmp(err.identifier, 'modalFeatureVector:missingRefinementResults');
end
assert(threw, 'Test3b:missingFile', 'did not throw modalFeatureVector:missingRefinementResults for a missing file');
fprintf('  (3b) missing-file error: PASS\n');

% (3c) Malformed file (missing a required field).
malformedFile = fullfile(scratchDir, 'malformed.mat');
save(malformedFile, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H'); % 'g' omitted
threw = false;
try
    loadFeatureAnchors(malformedFile);
catch err
    threw = strcmp(err.identifier, 'modalFeatureVector:malformedRefinementResults');
end
assert(threw, 'Test3c:malformed', 'did not throw modalFeatureVector:malformedRefinementResults when a field is missing');
fprintf('  (3c) malformed-file error: PASS\n');

% (3d) Feature at its window edge (idx=1): must error, not silently clip.
refinedEdge = refined;
refinedEdge(1).featureOmega = refinedEdge(1).omega(1); % force idx=1
edgeFile = fullfile(scratchDir, 'edge.mat');
refined = refinedEdge; %#ok<NASGU>
save(edgeFile, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g');
threw = false;
try
    loadFeatureAnchors(edgeFile);
catch err
    threw = strcmp(err.identifier, 'modalFeatureVector:featureAtWindowEdge');
end
assert(threw, 'Test3d:atEdge', 'did not throw modalFeatureVector:featureAtWindowEdge for an edge-of-window feature');
fprintf('  (3d) feature-at-window-edge error: PASS\n');

% (3e) featureOmega does not exactly match any grid point.
refinedMismatch = refinedEdge;
refinedMismatch(1).featureOmega = 2.0 + 1e-6; % off the grid
mismatchFile = fullfile(scratchDir, 'mismatch.mat');
refined = refinedMismatch; %#ok<NASGU>
save(mismatchFile, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g');
threw = false;
try
    loadFeatureAnchors(mismatchFile);
catch err
    threw = strcmp(err.identifier, 'modalFeatureVector:featureOmegaMismatch');
end
assert(threw, 'Test3e:mismatch', 'did not throw modalFeatureVector:featureOmegaMismatch for an off-grid featureOmega');
fprintf('  (3e) featureOmega/grid mismatch error: PASS\n');

% (3f) Wrong label order: swap the first two features' labels (and their
%      row/kind, so it is still an internally-consistent struct, just in
%      the wrong sequence) -- must be refused outright.
[refinedSwap, modeListSwap, nuS, MS, PS, NS, HS, gS] = makeSyntheticFeatureRefinement(); %#ok<ASGLU>
tmpLabel = refinedSwap(1).label; tmpRow = refinedSwap(1).row; tmpKind = refinedSwap(1).kind;
refinedSwap(1).label = refinedSwap(2).label; refinedSwap(1).row = refinedSwap(2).row; refinedSwap(1).kind = refinedSwap(2).kind;
refinedSwap(2).label = tmpLabel; refinedSwap(2).row = tmpRow; refinedSwap(2).kind = tmpKind;
swapFile = fullfile(scratchDir, 'swapped_order.mat');
refined = refinedSwap; modeList = modeListSwap; nu = nuS; M = MS; P = PS; N = NS; H = HS; g = gS; %#ok<NASGU>
save(swapFile, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g');
threw = false;
try
    loadFeatureAnchors(swapFile);
catch err
    threw = strcmp(err.identifier, 'modalFeatureVector:unexpectedFeatureOrder');
end
assert(threw, 'Test3f:wrongOrder', 'did not throw modalFeatureVector:unexpectedFeatureOrder when the first two feature labels were swapped');
fprintf('  (3f) wrong feature-label-order error: PASS\n\n');

%% TEST 4: modalFeatureVector.m -- end-to-end smoke test at cheap sanity
%  parameters, via the 'RefinementResultsFile' override (never touches
%  the real modalFeatureRefinementResults.mat), against the real
%  Forward Model.
fprintf('--- TEST 4: modalFeatureVector.m end-to-end smoke test ------\n');

% Cheap sanity forward-model parameters (NOT the real validated nominal
% M=50,P=10,N=10 -- this test is about wiring/shape/consistency
% correctness, not physical accuracy) baked into the synthetic results
% file itself, so the call below exercises the "defaults pulled from the
% frozen anchors file" path with no overrides. Full 5-feature synthetic
% set from makeSyntheticFeatureRefinement.m.
[refined, modeList, nu, M, P, N, H, g] = makeSyntheticFeatureRefinement(); %#ok<ASGLU>
smokeFile = fullfile(scratchDir, 'smoke.mat');
save(smokeFile, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g');

betaTest = 1e-2; gammaTest = 0.1; Rtest = 1.0;
ticS = tic;
[fVecLoose, detailsLoose] = modalFeatureVector(betaTest, gammaTest, Rtest, 'RefinementResultsFile', smokeFile);
fprintf('  modalFeatureVector call completed in %.1fs\n', toc(ticS));

assert(isequal(size(fVecLoose), [5, 1]), 'Test4:size', 'fVec is %s, expected 5x1', mat2str(size(fVecLoose)));
assert(all(isfinite(fVecLoose)), 'Test4:finite', 'fVec contains non-finite values: %s', mat2str(fVecLoose));
assert(numel(detailsLoose) == 5, 'Test4:detailsSize', 'details has %d entries, expected 5', numel(detailsLoose));

expectedAnchors = { [1.9, 2.0, 2.1], [2.4, 2.5, 2.6], [2.9, 3.0, 3.1], ...
    [3.4, 3.5, 3.6], [3.9, 4.0, 4.1] };
anyFlagged = false;
for k = 1:5
    assert(isequal(detailsLoose(k).omega3, expectedAnchors{k}), 'Test4:anchors', ...
        'feature %d (%s): unexpected anchor bracket [%s]', k, detailsLoose(k).label, mat2str(detailsLoose(k).omega3));
    assert(detailsLoose(k).omegaStar == fVecLoose(k), 'Test4:consistency', ...
        'details(%d).omegaStar does not match fVec(%d)', k, k);
    assert(max(abs(detailsLoose(k).y3 - abs(detailsLoose(k).A3).^2)) < 1e-12, 'Test4:y3', ...
        'details(%d).y3 does not equal |A3|^2', k);

    % Both-bounds bracket sanity (reviewer feedback: the original version
    % of this test only checked the lower bound).
    if detailsLoose(k).extrapolated
        assert(detailsLoose(k).omegaStar < detailsLoose(k).omega3(1) - 10*eps || ...
               detailsLoose(k).omegaStar > detailsLoose(k).omega3(3) + 10*eps, ...
            'Test4:extrapSanity', 'feature %d: extrapolated flag is set but the vertex is actually inside its bracket', k);
    else
        assert(detailsLoose(k).omegaStar >= detailsLoose(k).omega3(1) - 10*eps && ...
               detailsLoose(k).omegaStar <= detailsLoose(k).omega3(3) + 10*eps, ...
            'Test4:bracketSanity', 'feature %d: vertex falls outside its anchor bracket [%.4f, %.4f] without the extrapolated flag set', ...
            k, detailsLoose(k).omega3(1), detailsLoose(k).omega3(3));
    end

    if detailsLoose(k).degenerate || detailsLoose(k).extrapolated || detailsLoose(k).kindMismatch
        anyFlagged = true;
        fprintf('  NOTE: feature %s flagged under this arbitrary synthetic anchor set (degenerate=%d, extrapolated=%d, kindMismatch=%d) -- expected/benign here, these anchors are NOT real feature locations.\n', ...
            detailsLoose(k).label, detailsLoose(k).degenerate, detailsLoose(k).extrapolated, detailsLoose(k).kindMismatch);
    end
end
fprintf('  (4a) shapes, internal consistency, and correct fixed anchors (both-bounds bracket check): PASS\n');
fprintf('       fVec = [%.4f; %.4f; %.4f; %.4f; %.4f]\n', fVecLoose(1), fVecLoose(2), fVecLoose(3), fVecLoose(4), fVecLoose(5));

threw = false;
try
    modalFeatureVector(betaTest, gammaTest, Rtest, 'RefinementResultsFile', fullfile(scratchDir, 'nope.mat'));
catch err
    threw = strcmp(err.identifier, 'modalFeatureVector:missingRefinementResults');
end
assert(threw, 'Test4b:missingFile', 'modalFeatureVector did not propagate the missing-results-file error');
fprintf('  (4b) missing RefinementResultsFile propagates correctly: PASS\n\n');

%% TEST 5: 'Strict' option
fprintf('--- TEST 5: ''Strict'' option --------------------------------\n');

% Reuses the loose (Strict=false, default) call from TEST 4 above:
% detailsLoose already tells us whether ANY feature was flagged under
% this synthetic anchor set. Because the anchors here are arbitrary and
% uninformed relative to the real feature locations, which case actually
% occurs is a genuine property of the forward model at these frequencies
% and is not hardcoded -- this test adapts to whichever case is true on
% this run rather than assuming one.
if anyFlagged
    threw = false;
    try
        modalFeatureVector(betaTest, gammaTest, Rtest, 'RefinementResultsFile', smokeFile, 'Strict', true);
    catch err
        threw = strcmp(err.identifier, 'modalFeatureVector:invalidFeatureFit');
    end
    assert(threw, 'Test5:strictThrows', ...
        '''Strict'',true did not throw modalFeatureVector:invalidFeatureFit despite an invalid-fit flag being set under Strict=false');
    fprintf('  (5a) at least one feature was flagged under Strict=false; ''Strict'',true correctly throws modalFeatureVector:invalidFeatureFit: PASS\n\n');
else
    [fVecStrict, ~] = modalFeatureVector(betaTest, gammaTest, Rtest, 'RefinementResultsFile', smokeFile, 'Strict', true);
    assert(isequal(fVecStrict, fVecLoose), 'Test5:strictConsistent', ...
        '''Strict'',true changed the returned fVec relative to Strict=false when no flags were set');
    fprintf('  (5a) no feature was flagged under Strict=false; ''Strict'',true correctly does NOT throw and returns the identical fVec: PASS\n\n');
end

%% TEST 6: integer validation for M, P, N, NTheta
fprintf('--- TEST 6: integer validation (M, P, N, NTheta) -------------\n');

badParams = {'M', 10.4; 'P', 6.5; 'N', 10.2; 'NTheta', 90.5};
for r = 1:size(badParams, 1)
    threw = false;
    try
        modalFeatureVector(betaTest, gammaTest, Rtest, 'RefinementResultsFile', smokeFile, badParams{r,1}, badParams{r,2});
    catch %#ok<CTCH>
        % Identifier is whatever this MATLAB/Octave version's inputParser
        % validation-function failure uses (it differs between them);
        % what matters is that a non-integer value is rejected at all.
        threw = true;
    end
    assert(threw, sprintf('Test6:%s', badParams{r,1}), ...
        '''%s'',%.1f (non-integer) was NOT rejected -- expected a positive-integer validation error', badParams{r,1}, badParams{r,2});
end
fprintf('  (6a) non-integer M, P, N, NTheta are all rejected: PASS\n\n');

%% TEST 7: override-file independence (no cross-call contamination)
fprintf('--- TEST 7: override-file independence -----------------------\n');

[refinedA, modeListA, nuA, MA, PA, NA, HA, gA] = makeSyntheticFeatureRefinement(); % default centres
fileA = fullfile(scratchDir, 'variantA.mat');
refined = refinedA; modeList = modeListA; nu = nuA; M = MA; P = PA; N = NA; H = HA; g = gA; %#ok<NASGU>
save(fileA, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g');

[refinedB, modeListB, nuB, MB, PB, NB, HB, gB] = makeSyntheticFeatureRefinement('Centers', [5.0, 5.5, 6.0, 6.5, 7.0]);
fileB = fullfile(scratchDir, 'variantB.mat');
refined = refinedB; modeList = modeListB; nu = nuB; M = MB; P = PB; N = NB; H = HB; g = gB; %#ok<NASGU>
save(fileB, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g');

[~, detA] = modalFeatureVector(betaTest, gammaTest, Rtest, 'RefinementResultsFile', fileA);
[~, detB] = modalFeatureVector(betaTest, gammaTest, Rtest, 'RefinementResultsFile', fileB);
assert(isequal(detA(1).omega3, [1.9, 2.0, 2.1]), 'Test7:variantA', 'variant A returned unexpected anchors: %s', mat2str(detA(1).omega3));
assert(isequal(detB(1).omega3, [4.9, 5.0, 5.1]), 'Test7:variantB', 'variant B returned unexpected anchors: %s', mat2str(detB(1).omega3));
fprintf('  (7a) two different override files in a row give their own, independent, correct anchors: PASS\n\n');

fprintf('================================================================\n');
fprintf(' ALL modalFeatureVectorTester.m TESTS PASSED\n');
fprintf('================================================================\n');