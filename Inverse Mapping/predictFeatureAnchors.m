function [omega3Pred, fPred, shift] = predictFeatureAnchors(omega3Ref, p, model)
%PREDICTFEATUREANCHORS  Jacobian-predicted anchor brackets for
%   modalFeatureVector.m's 'AnchorMode','predicted'.
%
% The fixed anchors of loadFeatureAnchors.m bracket each feature at the
% reference point p0 only; the reference-point probe showed they are left
% after ~0.8% in R. This function moves each feature's three anchors to
% where the local scaled Jacobian says that feature should now be, so the
% same three-point quadratic fit can follow it over a much larger region.
% The prediction only decides WHERE TO LOOK: the returned feature value is
% still the vertex of the quadratic fitted to the EMM field at the moved
% anchors (modalFeatureVector.m), never the predicted centre itself.
%
% USAGE
%   [omega3Pred, fPred, shift] = predictFeatureAnchors(omega3Ref, p, model)
%
% PREDICTION (log-linear, consistent with J_s ~ d ln f / d ln p):
%   fPred      = f0 .* exp( Js * log(p ./ p0) )          (K x 1)
%   shift      = fPred - f0                              (K x 1)
%   omega3Pred = omega3Ref + shift   (row k shifted rigidly by shift(k))
%
%   - Each bracket is translated, never stretched: anchor spacing is
%     exactly that of omega3Ref.
%   - At p = p0 the shift is exactly zero, so omega3Pred == omega3Ref bit
%     for bit and predicted mode reproduces fixed mode exactly.
%   - The map is smooth (C-infinity) in p, so finite-difference derivatives
%     of the extracted feature see no jumps from the anchor placement.
%   - The exponential form keeps fPred > 0 and treats +delta and -delta in
%     ln p symmetrically (fPred(p0*e^d)/f0 = f0/fPred(p0*e^-d)).
%
% INPUTS
%   omega3Ref : K x 3 reference anchors, one row per feature, strictly
%               ascending within each row (loadFeatureAnchors' omega3 rows).
%   p         : n x 1 trial parameters, finite and > 0 (n = 3 for
%               [beta; gamma; R]).
%   model     : struct with
%                 p0 : n x 1 reference parameters (> 0)
%                 f0 : K x 1 reference feature vector (> 0), i.e.
%                      modalFeatureVector(p0) in FIXED mode
%                 Js : K x n scaled Jacobian, Js = diag(1./f0)*J*diag(p0)
%               Other fields (labels, provenance) are ignored here.
%
% CONSISTENCY CHECK: each f0(k) must lie inside its own reference bracket
% [omega3Ref(k,1), omega3Ref(k,3)]. It must, if the model was built at the
% point the anchors were validated at; if not, the model and the anchor
% file do not belong together (e.g. a model from a different refinement
% run), and predicting with it would silently misplace every bracket.
%
% ERRORS: predictFeatureAnchors:badAnchors, :badP, :badModel,
%   :sizeMismatch, :modelAnchorMismatch, :nonPositiveAnchor (a prediction so
%   far out that an anchor would be <= 0 rad/s).
%
% BUILDING A MODEL from computeSensitivityJacobian's output at the
% reference point (fixed mode):
%   [J, info]    = computeSensitivityJacobian(beta0, gamma0, R0, 1e-4);
%   model.p0     = info.p0;
%   model.f0     = info.f0;
%   model.Js     = diag(1 ./ info.f0) * J * diag(info.p0);
%   model.labels = {info.details0.label};
%   model.epsilon = info.epsilon;                 % provenance
%   save('anchorPredictionModel.mat', '-struct', 'model');

%% Validate
if ~(isnumeric(omega3Ref) && isreal(omega3Ref) && ismatrix(omega3Ref) && size(omega3Ref, 2) == 3 && ...
        size(omega3Ref, 1) >= 1 && all(isfinite(omega3Ref(:))) && all(all(diff(omega3Ref, 1, 2) > 0)))
    error('predictFeatureAnchors:badAnchors', ...
        'omega3Ref must be a real, finite K x 3 matrix with strictly ascending rows.');
end
K = size(omega3Ref, 1);

if ~(isnumeric(p) && isvector(p) && isreal(p) && all(isfinite(p)) && all(p > 0))
    error('predictFeatureAnchors:badP', 'p must be a real vector of finite, strictly positive values.');
end
p = p(:);
n = numel(p);

if ~(isstruct(model) && isscalar(model) && all(isfield(model, {'p0', 'f0', 'Js'})))
    error('predictFeatureAnchors:badModel', 'model must be a struct with fields p0, f0 and Js.');
end
p0 = model.p0; f0 = model.f0; Js = model.Js;
if ~(isnumeric(p0) && isvector(p0) && isreal(p0) && all(isfinite(p0)) && all(p0 > 0))
    error('predictFeatureAnchors:badModel', 'model.p0 must be a real vector of finite, strictly positive values.');
end
if ~(isnumeric(f0) && isvector(f0) && isreal(f0) && all(isfinite(f0)) && all(f0 > 0))
    error('predictFeatureAnchors:badModel', 'model.f0 must be a real vector of finite, strictly positive values.');
end
if ~(isnumeric(Js) && ismatrix(Js) && isreal(Js) && all(isfinite(Js(:))))
    error('predictFeatureAnchors:badModel', 'model.Js must be a real, finite matrix.');
end
p0 = p0(:); f0 = f0(:);
if numel(p0) ~= n || numel(f0) ~= K || ~isequal(size(Js), [K, n])
    error('predictFeatureAnchors:sizeMismatch', ...
        ['Sizes do not agree: p has %d elements, model.p0 %d, omega3Ref %d rows, model.f0 %d, ' ...
         'model.Js is %s (expected %d x %d).'], n, numel(p0), K, numel(f0), mat2str(size(Js)), K, n);
end

outside = f0 < omega3Ref(:, 1) | f0 > omega3Ref(:, 3);
if any(outside)
    k = find(outside, 1);
    error('predictFeatureAnchors:modelAnchorMismatch', ...
        ['model.f0(%d) = %.6g lies outside its reference anchor bracket [%.6g, %.6g]: this model was ' ...
         'not built at the point these anchors were validated at, so it cannot be used to move them.'], ...
        k, f0(k), omega3Ref(k, 1), omega3Ref(k, 3));
end

%% Predict
fPred = f0 .* exp(Js * log(p ./ p0));
shift = fPred - f0;
omega3Pred = omega3Ref + shift;       % implicit expansion: rigid shift per row

if any(omega3Pred(:) <= 0)
    k = find(any(omega3Pred <= 0, 2), 1);
    error('predictFeatureAnchors:nonPositiveAnchor', ...
        'Predicted anchors for feature %d reach %.6g rad/s (<= 0): p = %s is far outside the prediction''s range.', ...
        k, min(omega3Pred(k, :)), mat2str(p.', 6));
end
end