function [fVec, details] = modalFeatureVector(beta, gamma, R, varargin)
%MODALFEATUREVECTOR  Continuous, sub-grid-refined feature-location vector
%   f(beta,gamma,R) for the five candidate modal-response features found
%   by modalRAOReconnaissance.m / modalFeatureRefinement.m.
%
% DESIGN (resolves the open question left at the end of
% modalFeatureRefinement.m's header): feature locations are estimated by
% FIXED-ANCHOR sub-grid parabolic interpolation, not by re-searching a
% grid at every trial (beta,gamma,R). Concretely, for each of the five
% features this function:
%   1. Uses three FIXED omega anchors -- the three points, from the
%      nominal-parameter fine grid already computed by
%      modalFeatureRefinement.m, that bracket that feature's discrete
%      argmax/argmin (i.e. omega(idx-1), omega(idx), omega(idx+1) around
%      the already-extracted refined(w).featureOmega). These anchors are
%      fixed for the whole life of this function (see loadFeatureAnchors.m)
%      and do NOT move when (beta,gamma,R) changes.
%   2. Re-solves the EMM field and projects onto the relevant single mode
%      at exactly those three omegas, at the TRIAL (beta,gamma,R).
%   3. Fits the unique quadratic y(omega) = a*omega^2 + b*omega + c
%      through the three (omega, |A|^2) points and takes the vertex
%      omega* = -b/(2a) as the continuous, sub-grid feature location
%      (fitFeatureFrom3Points.m; quadFit3.m does the underlying fit).
%
% (Since the parameter-inversion stage there is also an optional
% 'AnchorMode','predicted' in which the three anchors are moved smoothly
% with (beta,gamma,R) by a Jacobian prediction; see NAME-VALUE OPTIONS.
% The default remains the fixed anchors described here, unchanged.)
%
% WHY FIXED ANCHORS, NOT RE-SEARCH: a discrete argmax/argmin recomputed
% on a fixed grid at every trial is a step function of (beta,gamma,R) --
% it can sit dead flat (derivative artificially zero) until a parameter
% perturbation crosses a grid-cell boundary, then jump by a full grid
% step. That is exactly the "broad-maximum caveat" already documented in
% the companion writeup. Fixing the three omegas and only ever
% re-evaluating the (smooth) complex field A(omega) at those omegas, then
% reading off a continuously-varying vertex, gives an observable that
% varies smoothly with the parameters -- which is what a finite-difference
% Jacobian needs.
%
% USAGE
%   fVec = modalFeatureVector(beta, gamma, R)
%   [fVec, details] = modalFeatureVector(beta, gamma, R, 'Name', Value, ...)
%
%   fVec    : 5x1 vector of sub-grid-refined omega locations, in the
%             FIXED order [A_{2,0}; A_{0,1}; A_{1,1}; A_{2,1}; A_{0,0}].
%             loadFeatureAnchors.m refuses to load a results file whose
%             refined(:) labels are not in exactly this order, so this
%             ordering cannot silently drift.
%   details : 1x5 struct array with per-feature diagnostics: label, row,
%             kind, omega3 (the anchors actually used at THIS trial: the
%             fixed ones, or the predicted ones in 'predicted' mode), A3
%             (complex field there at THIS trial), y3 (|A3|^2), the fitted
%             [a,b,c], omegaStar, three flags -- degenerate, extrapolated,
%             kindMismatch (see below) -- and anchorMode, omega3Ref (the
%             fixed reference anchors) and anchorShift (omega3 - omega3Ref,
%             0 in fixed mode).
%
% NAME-VALUE OPTIONS (all optional; default to the values frozen in
% modalFeatureRefinementResults.mat, i.e. the validated nominal settings)
%   'nu'                 : override, positive scalar.
%   'M', 'P', 'N'         : override, each a positive INTEGER (truncation
%                          orders are counts, not continuous quantities;
%                          a non-integer value is rejected here rather
%                          than silently truncated or passed through to
%                          the forward model).
%   'H', 'g'              : override, each a positive scalar.
%       Physically only (beta,gamma,R) are the unknowns being recovered;
%       these six overrides exist mainly so a truncation-convergence- or
%       geometry-sensitivity-style check can be run without editing this
%       file.
%   'NTheta'              : angular quadrature points for
%                          projectEMMOntoPlateModes.m (default 90, the
%                          value used throughout this project). Must be
%                          a positive integer.
%   'RefinementResultsFile' : path to the .mat file the fixed anchors are
%                          loaded from (default: modalFeatureRefinementResults.mat
%                          next to this file). Exists so
%                          modalFeatureVectorTester.m can point at a
%                          synthetic file without touching the real one;
%                          not intended for normal use.
%   'AnchorMode'          : 'fixed' (default) or 'predicted'.
%       'fixed'     : the original behaviour, unchanged -- the three anchors
%                     per feature are the frozen ones from the refinement
%                     file. Use this for Jacobian / identifiability work.
%       'predicted' : each feature's three anchors are translated rigidly
%                     (same spacing) to where the local scaled Jacobian
%                     predicts that feature has moved at this trial
%                     (predictFeatureAnchors.m, log-linear prediction).
%                     Everything else is identical: same EMM solves at the
%                     (moved) anchors, same quadratic fit, same three
%                     validity checks -- now judged against the MOVED
%                     bracket, so 'extrapolated' means "the prediction
%                     missed". The returned value is the fitted vertex,
%                     never the predicted centre. At (beta,gamma,R) = the
%                     model's p0 the shift is exactly zero and the result
%                     is identical to fixed mode. Built for the nonlinear
%                     inversion, where the fixed brackets are left after
%                     ~0.8% in R.
%   'AnchorModel'         : required with 'predicted' (and rejected with
%                          'fixed', so it can never be silently ignored).
%                          Either a struct with fields p0 (3x1), f0 (5x1),
%                          Js (5x3) and preferably labels (1x5 cell), or the
%                          path to a .mat file saved from such a struct (see
%                          predictFeatureAnchors.m for how to build one from
%                          computeSensitivityJacobian's output). If labels
%                          are present they must match the anchor file's
%                          feature labels exactly.
%   'Strict'              : false (default) returns a value and warns
%                          (via fprintf) whenever degenerate, extrapolated
%                          or kindMismatch is true for a feature. true
%                          instead throws immediately on any of the three
%                          -- use this inside computeSensitivityJacobian.m
%                          (or any other automated caller) so a single
%                          bad perturbation cannot silently contaminate a
%                          finite-difference derivative; the interactive
%                          default is meant for exploratory command-line
%                          use where a human reads the warnings.
%
% FLAGS IN `details` (see fitFeatureFrom3Points.m for the exact decision
% logic; always surfaced here, not silently resolved, because this
% function is meant to be called deep inside a finite-difference loop
% where nobody is watching the command window unless 'Strict' is set)
%   degenerate    : the three (omega,|A|^2) points are numerically flat
%                   (no usable curvature) -- omegaStar falls back to
%                   whichever of the three anchor omegas has the extreme
%                   y-value, which is NOT a sub-grid estimate.
%   extrapolated  : the fitted vertex omega* falls outside the anchor
%                   bracket -- i.e. this trial's perturbation has
%                   plausibly moved the true feature far enough that the
%                   FIXED anchors no longer bracket it. The value is
%                   still returned (not clipped) when 'Strict' is false,
%                   since silently clipping would hide exactly the
%                   situation a Jacobian step-size study needs to see.
%   kindMismatch  : the fitted parabola opens the "wrong" way for this
%                   feature's kind (a>0 for a 'max' feature, or a<0 for a
%                   'zero' feature) -- a strong warning sign that the
%                   local shape has changed character under this
%                   perturbation, not just shifted.
%
% Also fails fast (modalFeatureVector:nonFiniteField), before ever
% reaching the quadratic fit, if the EMM/projection pipeline returns a
% NaN or Inf at any anchor: a non-finite value should be diagnosed at its
% source, not silently propagated into quadFit3 and disguised as a
% numeric-looking (but meaningless) omegaStar.
%
% DEPENDS ON modalFeatureRefinement.m having already been run in this
% folder (modalFeatureRefinementResults.mat must exist next to this
% file): the fixed anchors are read from it, never re-hardcoded here, so
% this function cannot silently drift out of sync with the refinement
% stage it builds on. Anchors are reloaded from disk on every call (no
% persistent cache): the .mat file is tiny next to the cost of the 15 EMM
% solves this function makes, so caching bought negligible speed at the
% cost of a real footgun -- a rerun of modalFeatureRefinement.m producing
% new anchors would otherwise go unnoticed for the rest of the MATLAB
% session.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', 'Forward Model', 'Animation'));

%% Parse inputs
isPosIntScalar = @(x) isnumeric(x) && isscalar(x) && x > 0 && isfinite(x) && mod(x, 1) == 0;

p = inputParser;
addRequired(p, 'beta',  @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'gamma', @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'R',     @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'nu', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'M',  [], @(x) isempty(x) || isPosIntScalar(x));
addParameter(p, 'P',  [], @(x) isempty(x) || isPosIntScalar(x));
addParameter(p, 'N',  [], @(x) isempty(x) || isPosIntScalar(x));
addParameter(p, 'H',  [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'g',  [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'NTheta', 90, @(x) isPosIntScalar(x));
addParameter(p, 'RefinementResultsFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Strict', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'AnchorMode', 'fixed', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    any(strcmpi(char(x), {'fixed', 'predicted'})));
addParameter(p, 'AnchorModel', [], @(x) isempty(x) || isstruct(x) || ischar(x) || (isstring(x) && isscalar(x)));
parse(p, beta, gamma, R, varargin{:});
opts = p.Results;
anchorMode = lower(char(opts.AnchorMode));

if isempty(opts.RefinementResultsFile)
    refFile = fullfile(thisDir, 'modalFeatureRefinementResults.mat');
else
    refFile = char(opts.RefinementResultsFile);
end

%% Load the fixed anchors (always fresh -- see header note on caching)
A = loadFeatureAnchors(refFile);

nu = opts.nu; if isempty(nu), nu = A.nu; end
M  = opts.M;  if isempty(M),  M  = A.M;  end
P  = opts.P;  if isempty(P),  P  = A.P;  end
N  = opts.N;  if isempty(N),  N  = A.N;  end
H  = opts.H;  if isempty(H),  H  = A.H;  end
g  = opts.g;  if isempty(g),  g  = A.g;  end

%% Anchors used at this trial: fixed (default) or Jacobian-predicted
omega3Ref = cell2mat(arrayfun(@(f) f.omega3(:).', A.features(:), 'UniformOutput', false));  % nFeat x 3, frozen anchors (row per feature, whatever orientation the file stores)
switch anchorMode
    case 'fixed'
        if ~isempty(opts.AnchorModel)
            error('modalFeatureVector:anchorModelWithoutPredictedMode', ...
                ['''AnchorModel'' was given but ''AnchorMode'' is ''fixed'', so it would be ignored. ' ...
                 'Pass ''AnchorMode'',''predicted'' to use it, or drop it.']);
        end
        omega3Used = omega3Ref;
    case 'predicted'
        model = resolveAnchorModel(opts.AnchorModel, {A.features.label});
        omega3Used = predictFeatureAnchors(omega3Ref, [beta; gamma; R], model);
end

%% Evaluate each feature at its fixed anchors, at THIS trial (beta,gamma,R)
nFeat = numel(A.features);
fVec = nan(nFeat, 1);
details = repmat(struct('label', '', 'row', 0, 'kind', '', 'omega3', [], ...
    'A3', [], 'y3', [], 'a', NaN, 'b', NaN, 'c', NaN, 'omegaStar', NaN, ...
    'degenerate', false, 'extrapolated', false, 'kindMismatch', false, ...
    'anchorMode', '', 'omega3Ref', [], 'anchorShift', 0), 1, nFeat);

for k = 1:nFeat
    feat = A.features(k);
    if strcmp(anchorMode, 'fixed')
        omega3 = feat.omega3;          % exactly as before this option existed
    else
        omega3 = omega3Used(k, :);
    end
    A3 = complex(zeros(1, 3));
    for a3 = 1:3
        alpha = H * omega3(a3)^2 / g;
        data = precomputeDeflectionData(alpha, beta, gamma, R, nu, M, P, N);
        wFun = @(r, theta) evaluateDeflection(data, r, theta);
        A3(a3) = projectEMMOntoPlateModes(wFun, A.modeList(feat.row, :), R, nu, 'NTheta', opts.NTheta);
    end

    if any(~isfinite(A3))
        badIdx = find(~isfinite(A3), 1);
        error('modalFeatureVector:nonFiniteField', ...
            ['%s: the EMM/projection pipeline returned a non-finite value (%s) at ' ...
             'omega=%.6g for trial (beta=%.6g, gamma=%.6g, R=%.6g) -- failing at the ' ...
             'source rather than feeding NaN/Inf into the quadratic fit.'], ...
            feat.label, mat2str(A3(badIdx)), omega3(badIdx), beta, gamma, R);
    end

    y3 = abs(A3).^2;
    [omegaStar, degenerate, extrapolated, kindMismatch, a2, b2, c2] = ...
        fitFeatureFrom3Points(omega3, y3, feat.kind);

    if opts.Strict && (degenerate || extrapolated || kindMismatch)
        badFlags = {};
        if degenerate,    badFlags{end+1} = 'degenerate';    end %#ok<AGROW>
        if extrapolated,  badFlags{end+1} = 'extrapolated';  end %#ok<AGROW>
        if kindMismatch,  badFlags{end+1} = 'kindMismatch';  end %#ok<AGROW>
        error('modalFeatureVector:invalidFeatureFit', ...
            ['%s: sub-grid fit failed its validity check(s) [%s] at trial ' ...
             '(beta=%.6g, gamma=%.6g, R=%.6g) with ''Strict'',true -- refusing to ' ...
             'return a value that could silently contaminate a derivative.'], ...
            feat.label, strjoin(badFlags, ', '), beta, gamma, R);
    end

    if degenerate
        fprintf('  WARNING: %s sub-grid fit is numerically flat (no usable curvature) at this trial; falling back to the discrete anchor extremum, omega = %.4f.\n', feat.label, omegaStar);
    end
    if extrapolated
        fprintf('  WARNING: %s sub-grid vertex (omega=%.4f) falls OUTSIDE its %s anchor bracket [%.4f, %.4f]; this trial may have moved the feature past what these anchors can resolve.\n', ...
            feat.label, omegaStar, anchorMode, min(omega3), max(omega3));
    end
    if kindMismatch
        fprintf('  WARNING: %s parabola opens the "wrong" way for a %s feature at this trial (a=%.3e); the local shape may have changed character, not just shifted.\n', ...
            feat.label, feat.kind, a2);
    end

    fVec(k) = omegaStar;
    details(k).label = feat.label;
    details(k).row = feat.row;
    details(k).kind = feat.kind;
    details(k).omega3 = omega3;
    details(k).A3 = A3;
    details(k).y3 = y3;
    details(k).a = a2; details(k).b = b2; details(k).c = c2;
    details(k).omegaStar = omegaStar;
    details(k).degenerate = degenerate;
    details(k).extrapolated = extrapolated;
    details(k).kindMismatch = kindMismatch;
    details(k).anchorMode = anchorMode;
    details(k).omega3Ref = omega3Ref(k, :);
    details(k).anchorShift = omega3(2) - omega3Ref(k, 2);
end

end % modalFeatureVector

%% ------------------------------------------------------------------------
function model = resolveAnchorModel(modelIn, labels)
% Load / validate the prediction model for 'AnchorMode','predicted'. Only
% presence, source and label consistency are checked here; sizes, signs and
% the f0-inside-bracket consistency check are predictFeatureAnchors.m's.
if isempty(modelIn)
    error('modalFeatureVector:missingAnchorModel', ...
        ['''AnchorMode'',''predicted'' needs ''AnchorModel'' (a struct with p0, f0, Js, or a .mat ' ...
         'file saved from one) -- see predictFeatureAnchors.m for how to build it.']);
end
if isstruct(modelIn)
    model = modelIn;
else
    modelFile = char(modelIn);
    if ~isfile(modelFile)
        error('modalFeatureVector:missingAnchorModel', 'AnchorModel file %s not found.', modelFile);
    end
    model = load(modelFile);
end
if isfield(model, 'labels') && ~isequal(model.labels(:).', labels(:).')
    error('modalFeatureVector:anchorModelMismatch', ...
        'AnchorModel labels {%s} do not match the anchor file''s features {%s}.', ...
        strjoin(model.labels, ', '), strjoin(labels, ', '));
end
end