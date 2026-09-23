%% predictFeatureAnchorsTester.m
% Validation tests for predictFeatureAnchors.m and for the
% 'AnchorMode','predicted' option it adds to modalFeatureVector.m.
%
% TEST LIST
%   1. predictFeatureAnchors.m, pure arithmetic (no EMM, milliseconds):
%      (1a) p = p0 -> anchors, fPred and shift reproduce the reference
%           EXACTLY (bit for bit);
%      (1b) hand-computed log-linear prediction for a single-parameter move;
%      (1c) brackets translate rigidly: spacing unchanged, all three anchors
%           shift by exactly fPred - f0;
%      (1d) first-order agreement with J_s for a 1e-4 relative step;
%      (1e) symmetry in ln p: fPred(p0*e^d)/f0 = f0/fPred(p0*e^-d);
%      (1f) smoothness: central-difference derivatives match the analytic
%           derivative, and second differences along a 401-point sweep
%           match the analytic second derivative (no jumps anywhere);
%      (1g) every validation error, incl. the model/anchor consistency check.
%   2. modalFeatureVector.m option handling (errors raised before any EMM
%      solve): bad AnchorMode, predicted without a model, a model in fixed
%      mode, mismatched labels, missing model file.
%   3. modalFeatureVector.m against the REAL Forward Model, on the cheap
%      synthetic anchors of makeSyntheticFeatureRefinement.m (~8 calls):
%      (3a) predicted mode at p = p0 is identical to fixed mode (model
%           passed as a FILE path);
%      (3b) at p ~= p0 (model passed as a STRUCT): the anchors used are
%           exactly predictFeatureAnchors' output, spacing is unchanged, and
%           the returned feature is the vertex of the quadratic fitted at
%           those anchors, NOT the predicted centre;
%      (3c) 'extrapolated' is judged against the MOVED bracket, and
%           'Strict',true still throws in predicted mode whenever a check
%           fails;
%      (3d) central finite differences through predicted mode show no
%           jumps: D(h) and D(2h) agree (h = 1e-4 relative, in R).
%      The synthetic anchors are arbitrary (not real extrema), so 3a-3d
%      test WIRING and SMOOTHNESS, not physical accuracy.
%   4. OPTIONAL regression at the real reference point, run only if the real
%      modalFeatureRefinementResults.mat sits next to modalFeatureVector.m
%      (~4 calls, a few minutes): fixed mode still reproduces the recorded
%      f0; predicted mode at p0 is identical to it; and the predicted-mode
%      R-derivative is compared with the J_s R-column of the Jacobian
%      report (reported, not asserted: two estimators of one derivative).
%
% Lives in Inverse Mapping/Inverse Mapping Unit Tests/, alongside
% modalFeatureVectorTester.m (whose makeSyntheticFeatureRefinement.m it
% reuses); predictFeatureAnchors.m and modalFeatureVector.m are one level up.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);                                              % makeSyntheticFeatureRefinement.m
addpath(fullfile(thisDir, '..'));                              % predictFeatureAnchors.m, modalFeatureVector.m, ...
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));

clearvars -except thisDir
close all, clc

scratchDir = tempname;
mkdir(scratchDir);
cleanupObj = onCleanup(@() rmdir(scratchDir, 's'));

fprintf('================================================================\n');
fprintf(' predictFeatureAnchorsTester.m\n');
fprintf('================================================================\n\n');

labels = {'A_{2,0}', 'A_{0,1}', 'A_{1,1}', 'A_{2,1}', 'A_{0,0}'};

% Pure-arithmetic fixture: the real anchors, reference point, f0 and J_s
% from the Jacobian report (Sections 6, 7 and 9).
omega3Ref = [5.70 5.80 5.90; 5.50 5.65 5.80; 6.45 6.60 6.75; 4.90 4.95 5.00; 7.80 7.85 7.90];
model.p0 = [4.6985e-5; 1.4548e-3; 0.3830];
model.f0 = [5.836159; 5.676641; 6.647953; 4.959282; 7.846157];
model.Js = [-0.0210  0.0003 -0.4184
            -0.0578  0.0000 -0.2644
            -0.0568  0.0014 -0.2580
             0.0449 -0.0004 -0.6927
             0.0739 -0.0070 -0.7965];
model.labels = labels;
p0 = model.p0; f0 = model.f0; Js = model.Js;

%% TEST 1: predictFeatureAnchors.m -- pure arithmetic
fprintf('--- TEST 1: predictFeatureAnchors.m (no EMM) ----------------\n');

% (1a) Exact reduction at p0.
[w3, fP, sh] = predictFeatureAnchors(omega3Ref, p0, model);
assert(isequal(w3, omega3Ref) && isequal(fP, f0) && all(sh == 0), 'Test1a:identity', ...
    'at p = p0 the prediction does not reproduce the reference exactly');
fprintf('  (1a) p = p0: anchors, fPred and shift reproduce the reference bit for bit: PASS\n');

% (1b) Hand-computed: +10% beta only -> fPred = f0 .* 1.1.^Js(:,1).
p = p0 .* [1.1; 1; 1];
[w3, fP] = predictFeatureAnchors(omega3Ref, p, model);
expected = f0 .* 1.1 .^ Js(:, 1);
assert(max(abs(fP - expected) ./ expected) < 1e-14, 'Test1b:value', 'fPred ~= f0.*1.1.^Js(:,1)');
fprintf('  (1b) +10%% beta: fPred = f0 .* 1.1.^Js(:,beta) exactly: PASS\n');

% (1c) Rigid translation, spacing unchanged.
p = p0 .* [1.07; 0.8; 1.05];
[w3, fP, sh] = predictFeatureAnchors(omega3Ref, p, model);
assert(max(max(abs(diff(w3, 1, 2) - diff(omega3Ref, 1, 2)))) < 1e-12, 'Test1c:spacing', 'anchor spacing changed');
assert(max(max(abs((w3 - omega3Ref) - repmat(sh, 1, 3)))) < 1e-12 && max(abs(sh - (fP - f0))) < 1e-15, ...
    'Test1c:rigid', 'anchors did not all move by exactly shift = fPred - f0');
fprintf('  (1c) brackets translate rigidly by fPred - f0; spacing unchanged: PASS\n');

% (1d) First order: for a relative step of 1e-4 in parameter j,
%      shift ~= f0 .* Js(:,j) * log(1 + 1e-4), error O(step^2).
for j = 1:3
    e = ones(3, 1); e(j) = 1 + 1e-4;
    [~, ~, sh] = predictFeatureAnchors(omega3Ref, p0 .* e, model);
    lin = f0 .* Js(:, j) * log(1 + 1e-4);
    assert(max(abs(sh - lin)) <= 1e-7 * max(abs(f0)) + 1e-15, 'Test1d:firstOrder', ...
        'parameter %d: shift deviates from the first-order J_s prediction', j);
end
fprintf('  (1d) 1e-4 relative steps move each centre by f0.*Js(:,j)*dlnp to first order: PASS\n');

% (1e) Symmetry in ln p.
d = [0.08; -0.15; 0.03];
[~, fPlus] = predictFeatureAnchors(omega3Ref, p0 .* exp(d), model);
[~, fMinus] = predictFeatureAnchors(omega3Ref, p0 .* exp(-d), model);
assert(max(abs((fPlus ./ f0) - (f0 ./ fMinus))) < 1e-14, 'Test1e:symmetry', 'prediction is not symmetric in ln p');
fprintf('  (1e) fPred(p0*e^d)/f0 = f0/fPred(p0*e^-d): PASS\n');

% (1f) Smoothness. Analytic: d fPred / d ln p_j = fPred .* Js(:,j).
pOff = p0 .* [1.05; 0.9; 1.02];
[~, fOff] = predictFeatureAnchors(omega3Ref, pOff, model);
h = 1e-4;
for j = 1:3
    e = zeros(3, 1); e(j) = h;
    [~, fp] = predictFeatureAnchors(omega3Ref, pOff .* exp(e), model);
    [~, fm] = predictFeatureAnchors(omega3Ref, pOff .* exp(-e), model);
    Dfd = (fp - fm) / (2 * h);
    Dan = fOff .* Js(:, j);
    assert(max(abs(Dfd - Dan)) <= 1e-7 * max(abs(Dan)) + 1e-10, 'Test1f:derivative', ...
        'parameter %d: FD derivative of fPred disagrees with the analytic derivative', j);
end
t = linspace(-0.2, 0.2, 401); ht = t(2) - t(1);
S = zeros(5, numel(t));
for i = 1:numel(t)
    [~, ~, S(:, i)] = predictFeatureAnchors(omega3Ref, p0 .* [1; 1; exp(t(i))], model);
end
D2 = diff(S, 2, 2) / ht^2;
D2an = f0 .* Js(:, 3).^2 .* exp(Js(:, 3) * t(2:end-1));
assert(max(abs(D2(:) - D2an(:))) <= 1e-3 * max(abs(D2an(:))), 'Test1f:noJumps', ...
    'second differences along a +/-20%% R sweep do not match the analytic second derivative (jump or kink?)');
fprintf('  (1f) FD derivatives match analytic; 401-point R sweep has no jumps or kinks: PASS\n');

% (1g) Validation.
bad = model; bad.f0(5) = 7.95;                 % outside [7.80, 7.90]
badJs = model; badJs.Js = Js(:, 1:2);
noJs = rmfield(model, 'Js');
negF = model; negF.f0(1) = -1;
cases = {
    'non-ascending anchors',  @() predictFeatureAnchors(omega3Ref(:, [2 1 3]), p0, model), 'predictFeatureAnchors:badAnchors'
    'anchors K x 2',          @() predictFeatureAnchors(omega3Ref(:, 1:2), p0, model),     'predictFeatureAnchors:badAnchors'
    'p has a zero',           @() predictFeatureAnchors(omega3Ref, [p0(1); 0; p0(3)], model), 'predictFeatureAnchors:badP'
    'p has NaN',              @() predictFeatureAnchors(omega3Ref, [NaN; p0(2:3)], model), 'predictFeatureAnchors:badP'
    'model missing Js',       @() predictFeatureAnchors(omega3Ref, p0, noJs),              'predictFeatureAnchors:badModel'
    'model f0 negative',      @() predictFeatureAnchors(omega3Ref, p0, negF),              'predictFeatureAnchors:badModel'
    'Js wrong size',          @() predictFeatureAnchors(omega3Ref, p0, badJs),             'predictFeatureAnchors:sizeMismatch'
    'p wrong length',         @() predictFeatureAnchors(omega3Ref, p0(1:2), model),        'predictFeatureAnchors:sizeMismatch'
    'f0 outside its bracket', @() predictFeatureAnchors(omega3Ref, p0, bad),               'predictFeatureAnchors:modelAnchorMismatch'
    'anchor driven <= 0',     @() predictFeatureAnchors(omega3Ref, p0 .* [1; 1; 1e4], model), 'predictFeatureAnchors:nonPositiveAnchor'
};
for c = 1:size(cases, 1)
    gotId = '<no error>';
    try
        cases{c, 2}();
    catch err
        gotId = err.identifier;
    end
    assert(strcmp(gotId, cases{c, 3}), sprintf('Test1g:case%d', c), ...
        '%s: expected %s, got %s', cases{c, 1}, cases{c, 3}, gotId);
end
fprintf('  (1g) all %d invalid-input cases throw their documented identifier: PASS\n\n', size(cases, 1));

%% Synthetic anchor file + synthetic model for TESTS 2-3
[refined, modeList, nu, M, P, N, H, g] = makeSyntheticFeatureRefinement(); %#ok<ASGLU>
smokeFile = fullfile(scratchDir, 'smoke.mat');
save(smokeFile, 'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g');
omega3Syn = cell2mat(arrayfun(@(r) r.omega(2:4), refined(:), 'UniformOutput', false));  % [c-0.1 c c+0.1]

pT = [1e-2; 0.1; 1.0];                          % cheap test point (as in modalFeatureVectorTester)
synModel.p0 = pT;
synModel.f0 = [2.0; 2.5; 3.0; 3.5; 4.0];         % bracket centres: inside, as required
synModel.Js = Js;                                % realistic magnitudes; exact values irrelevant here
synModel.labels = labels;
synModelFile = fullfile(scratchDir, 'synModel.mat');
tmp = synModel; save(synModelFile, '-struct', 'tmp');
mfv = @(p, varargin) modalFeatureVector(p(1), p(2), p(3), 'RefinementResultsFile', smokeFile, varargin{:});

%% TEST 2: option handling in modalFeatureVector.m (no EMM solve reached)
fprintf('--- TEST 2: modalFeatureVector option handling --------------\n');

wrongLabels = synModel; wrongLabels.labels = labels([2 1 3 4 5]);
cases = {
    'unknown AnchorMode',        @() mfv(pT, 'AnchorMode', 'adaptive'),                               ''
    'predicted without a model', @() mfv(pT, 'AnchorMode', 'predicted'),                              'modalFeatureVector:missingAnchorModel'
    'model file missing',        @() mfv(pT, 'AnchorMode', 'predicted', 'AnchorModel', fullfile(scratchDir, 'nope.mat')), 'modalFeatureVector:missingAnchorModel'
    'model given in fixed mode', @() mfv(pT, 'AnchorModel', synModel),                                'modalFeatureVector:anchorModelWithoutPredictedMode'
    'model labels mismatch',     @() mfv(pT, 'AnchorMode', 'predicted', 'AnchorModel', wrongLabels),  'modalFeatureVector:anchorModelMismatch'
};
for c = 1:size(cases, 1)
    threw = false; gotId = '<no error>';
    try
        cases{c, 2}();
    catch err
        gotId = err.identifier;
        threw = isempty(cases{c, 3}) || strcmp(gotId, cases{c, 3});   % '' = any error (inputParser IDs differ MATLAB/Octave)
    end
    assert(threw, sprintf('Test2:case%d', c), '%s: expected %s, got %s', cases{c, 1}, cases{c, 3}, gotId);
end
fprintf('  (2) all %d option errors raised before any EMM solve: PASS\n\n', size(cases, 1));

%% TEST 3: modalFeatureVector.m, real Forward Model, synthetic anchors
fprintf('--- TEST 3: predicted mode against the real Forward Model --\n');
fprintf('  (%d-ish EMM-backed modalFeatureVector calls; warnings about the synthetic\n', 8);
fprintf('   anchors failing their checks are expected and are part of the test)\n');

% (3a) Predicted at p0 == fixed.
t3 = tic;
[fFix, dFix] = mfv(pT);
[fPre, dPre] = mfv(pT, 'AnchorMode', 'predicted', 'AnchorModel', synModelFile);
assert(isequal(fFix, fPre), 'Test3a:identical', 'predicted mode at p0 differs from fixed mode');
for k = 1:5
    assert(isequal(dFix(k).omega3(:).', dPre(k).omega3(:).') && dPre(k).anchorShift == 0 && ...
           isequal(dFix(k).y3, dPre(k).y3), 'Test3a:details', 'feature %d: anchors/field differ at p0', k);
    assert(strcmp(dFix(k).anchorMode, 'fixed') && strcmp(dPre(k).anchorMode, 'predicted'), 'Test3a:mode', ...
        'details.anchorMode not recorded');
end
fprintf('  (3a) predicted mode at p0 (model from FILE) identical to fixed mode: PASS  [%.0fs]\n', toc(t3));

% (3b) Off-reference: anchors used, spacing, value comes from the fit.
t3 = tic;
p1 = pT .* [1.05; 1; 1.02];
[fP1, dP1] = mfv(p1, 'AnchorMode', 'predicted', 'AnchorModel', synModel);
w3Expected = predictFeatureAnchors(omega3Syn, p1, synModel);
for k = 1:5
    assert(isequal(dP1(k).omega3, w3Expected(k, :)), 'Test3b:anchors', ...
        'feature %d: anchors used are not predictFeatureAnchors'' output', k);
    assert(max(abs(diff(dP1(k).omega3) - diff(omega3Syn(k, :)))) < 1e-12, 'Test3b:spacing', ...
        'feature %d: anchor spacing changed', k);
    assert(max(abs(dP1(k).y3 - abs(dP1(k).A3).^2)) < 1e-12, 'Test3b:y3', 'feature %d: y3 ~= |A3|^2', k);
    wFit = fitFeatureFrom3Points(dP1(k).omega3, dP1(k).y3, dP1(k).kind);
    assert(isequal(fP1(k), wFit, dP1(k).omegaStar), 'Test3b:fromFit', ...
        'feature %d: returned value is not the vertex fitted at the moved anchors', k);
end
assert(any(fP1(:) ~= w3Expected(:, 2)), 'Test3b:notCentre', 'returned features equal the predicted centres');
fprintf('  (3b) moved anchors = predictFeatureAnchors output, spacing kept, value = fitted vertex (not centre): PASS  [%.0fs]\n', toc(t3));

% (3c) Checks judged against the moved bracket; Strict still throws.
for k = 1:5
    outside = dP1(k).omegaStar < dP1(k).omega3(1) - 10*eps || dP1(k).omegaStar > dP1(k).omega3(3) + 10*eps;
    assert(dP1(k).extrapolated == outside, 'Test3c:extrapolated', ...
        'feature %d: extrapolated flag not consistent with the MOVED bracket', k);
end
anyFlag = any([dP1.degenerate] | [dP1.extrapolated] | [dP1.kindMismatch]);
t3 = tic; threw = false;
try
    mfv(p1, 'AnchorMode', 'predicted', 'AnchorModel', synModel, 'Strict', true);
catch err
    threw = strcmp(err.identifier, 'modalFeatureVector:invalidFeatureFit');
end
assert(threw == anyFlag, 'Test3c:strict', ...
    'Strict in predicted mode: threw = %d but a check failed = %d', threw, anyFlag);
fprintf('  (3c) flags judged against the moved bracket; Strict throws iff a check fails (%d): PASS  [%.0fs]\n', anyFlag, toc(t3));

% (3d) No jumps under central differences (relative step 1e-4 in R).
t3 = tic;
h = 1e-4;
fR = @(s) mfv(p1 .* [1; 1; 1 + s], 'AnchorMode', 'predicted', 'AnchorModel', synModel);
Fm2 = fR(-2*h); Fm1 = fR(-h); Fp1 = fR(h); Fp2 = fR(2*h);
Dh = (Fp1 - Fm1) / (2*h);
D2h = (Fp2 - Fm2) / (4*h);
relDiff = abs(Dh - D2h) ./ max(abs(Dh), 1e-12);
% Threshold 1e-2: a jump from anchor placement would make D(h) and D(2h)
% differ by O(1); the ~1e-4..1e-3 seen here is EMM round-off amplified by the
% (deliberately meaningless) synthetic features' flat vertices.
assert(all(relDiff < 1e-2), 'Test3d:smooth', ...
    'D(h) and D(2h) disagree by up to %.3g (relative): a jump or kink in predicted mode?', max(relDiff));
fprintf('  (3d) central FD through predicted mode is smooth: max |D(h)-D(2h)|/|D| = %.2g: PASS  [%.0fs]\n\n', ...
    max(relDiff), toc(t3));

%% TEST 4 (optional): real reference point regression
fprintf('--- TEST 4: real reference point (optional) -----------------\n');
realFile = fullfile(thisDir, '..', 'modalFeatureRefinementResults.mat');
if ~isfile(realFile)
    fprintf('  SKIPPED: %s not found.\n\n', realFile);
else
    t4 = tic;
    [f0Fix, d0Fix] = modalFeatureVector(p0(1), p0(2), p0(3));
    assert(max(abs(f0Fix - f0)) < 5e-7, 'Test4:fixedRegression', ...
        'fixed mode no longer reproduces the recorded f0 = %s (got %s)', mat2str(f0.', 7), mat2str(f0Fix.', 7));
    realModel = model; realModel.f0 = f0Fix;   % exact f0, report J_s
    [f0Pre, d0Pre] = modalFeatureVector(p0(1), p0(2), p0(3), 'AnchorMode', 'predicted', 'AnchorModel', realModel);
    assert(isequal(f0Fix, f0Pre) && all([d0Pre.anchorShift] == 0), 'Test4:predictedAtP0', ...
        'predicted mode at the real p0 differs from fixed mode');
    fprintf('  (4a) fixed mode reproduces the recorded f0; predicted mode at p0 identical: PASS  [%.0fs]\n', toc(t4));

    t4 = tic;
    h = 1e-4;
    fPlus  = modalFeatureVector(p0(1), p0(2), p0(3) * (1 + h), 'AnchorMode', 'predicted', 'AnchorModel', realModel, 'Strict', true);
    fMinus = modalFeatureVector(p0(1), p0(2), p0(3) * (1 - h), 'AnchorMode', 'predicted', 'AnchorModel', realModel, 'Strict', true);
    JsR = (fPlus - fMinus) ./ (2 * h * f0Fix);    % scaled R-column measured through predicted mode
    relDev = abs(JsR - Js(:, 3)) ./ abs(Js(:, 3));
    fprintf('       J_s R-column, report (fixed anchors) vs measured (predicted anchors):\n');
    for k = 1:5
        fprintf('         %-8s %9.4f  %9.4f   (%.2g%%)\n', labels{k}, Js(k, 3), JsR(k), 100 * relDev(k));
    end
    % Reported, not asserted: the two are different estimators of the same
    % derivative. With fixed anchors the fitted vertex's position bias
    % changes as the feature drifts inside the bracket; with moving anchors
    % it does not. A small difference is expected; a large one says the
    % fixed-anchor J_s carries that bias (worth knowing, not a bug here).
    if all(relDev < 0.05)
        fprintf('  (4b) predicted-mode R-derivative agrees with the report J_s (max %.2g%%): PASS  [%.0fs]\n\n', ...
            100 * max(relDev), toc(t4));
    else
        fprintf(['  (4b) NOTE: predicted-mode R-derivative differs from the report J_s by up to %.2g%%.\n' ...
                 '       Not a failure of this code; report this number (see comment in the tester).  [%.0fs]\n\n'], ...
            100 * max(relDev), toc(t4));
    end
end

fprintf('================================================================\n');
fprintf(' ALL predictFeatureAnchorsTester.m TESTS PASSED\n');
fprintf('================================================================\n');