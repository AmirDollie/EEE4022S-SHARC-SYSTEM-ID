function [J, info] = computeSensitivityJacobian(beta, gamma, R, epsilon, varargin)
%COMPUTESENSITIVITYJACOBIAN  Central-difference sensitivity Jacobian of
%   modalFeatureVector.m's feature map f(beta,gamma,R) at a reference
%   point (beta,gamma,R), with respect to p = [beta; gamma; R].
%
% USAGE
%   [J, info] = computeSensitivityJacobian(beta, gamma, R, epsilon)
%   [J, info] = computeSensitivityJacobian(beta, gamma, R, epsilon, Name, Value, ...)
%
% J is 5x3 (rows = the 5 features, in modalFeatureVector's fixed order
% [A_{2,0}; A_{0,1}; A_{1,1}; A_{2,1}; A_{0,0}]; columns = [beta,gamma,R]):
%
%   J(:,j) = ( f(p0 + deltaP(j)*e_j) - f(p0 - deltaP(j)*e_j) ) / (2*deltaP(j))
%
% using a RELATIVE step deltaP(j) = epsilon * p0(j), p0 = [beta;gamma;R]
% -- NOT an absolute step, since beta, gamma and R differ by orders of
% magnitude and a single absolute step could not sensibly perturb all
% three.
%
% EVERY modalFeatureVector call this function makes -- the baseline and
% all 6 perturbed evaluations -- always uses 'Strict',true, and this is
% NOT overridable (passing 'Strict' here throws
% computeSensitivityJacobian:strictReserved). If ANY perturbation moves a
% feature outside its fixed anchor bracket, or flips its curvature
% character (degenerate/extrapolated/kindMismatch -- see
% fitFeatureFrom3Points.m), the whole computation errors out immediately
% by propagating modalFeatureVector's own modalFeatureVector:invalidFeatureFit
% error. This is deliberate: a single bad perturbation must not be
% allowed to silently contaminate a derivative that later feeds an SVD.
%
% NAME-VALUE OPTIONS: anything modalFeatureVector.m itself accepts (nu,
% M, P, N, H, g, NTheta, RefinementResultsFile) is forwarded verbatim to
% every one of the 7 modalFeatureVector calls made here, so all 7 read
% the exact same fixed anchors and forward-model constants. Those
% options are NOT re-validated here -- modalFeatureVector.m already does
% that, and duplicating it here would risk the two silently drifting
% apart. 'Strict' is the one exception (see above): it is reserved by
% this function and cannot be passed through.
%
% info FIELDS
%   p0            : 3x1, [beta;gamma;R] as given
%   paramNames    : {'beta','gamma','R'}, in the column order of J/p0
%   epsilon       : the relative step passed in, echoed back
%   deltaP        : 3x1, the actual absolute perturbation used per
%                   parameter (epsilon*p0(j))
%   f0            : 5x1, baseline feature vector at p0 (for reference/
%                   plots only -- NOT used in the central-difference
%                   formula, which only uses the +/- evaluations)
%   fPlus, fMinus : 5x3, column j is f(p0) with parameter j perturbed by
%                   +deltaP(j) / -deltaP(j)
%   details0      : 1x5 struct array, modalFeatureVector's own
%                   diagnostics at the baseline p0
%   detailsPlus, detailsMinus : 1x3 cell arrays of 1x5 struct arrays --
%                   detailsPlus{j}/detailsMinus{j} are modalFeatureVector's
%                   diagnostics for the +deltaP(j)/-deltaP(j) evaluation
%                   of parameter j. Kept so that if a run at a LARGER
%                   epsilon later fails, you can trace exactly which
%                   feature, side and parameter moved, without
%                   re-running anything.
%
% The central-difference ASSEMBLY itself (the boxed formula above) is
% delegated to centralDiffJacobian.m, a pure function with no EMM
% dependency, specifically so that arithmetic can be unit-tested
% directly -- see computeSensitivityJacobianTester.m.
%
% This function does NOT, by itself, tell you whether epsilon is small
% enough for the central difference to have converged, nor whether
% (beta,gamma,R) is a physically realistic floe. The intended workflow
% is to call this at a few epsilon values (e.g. 1e-3, 1e-4, 1e-5) and
% compare the resulting J across them for consistency before trusting
% any single one, then only afterwards look at the (scaled) SVD.

%% Validate the physical parameters and epsilon ourselves -- these are
% NOT forwarded to modalFeatureVector, so nothing else checks them.
if ~(isnumeric(beta) && isscalar(beta) && isfinite(beta) && beta > 0)
    error('computeSensitivityJacobian:badBeta', 'beta must be a finite, positive scalar.');
end
if ~(isnumeric(gamma) && isscalar(gamma) && isfinite(gamma) && gamma > 0)
    error('computeSensitivityJacobian:badGamma', 'gamma must be a finite, positive scalar.');
end
if ~(isnumeric(R) && isscalar(R) && isfinite(R) && R > 0)
    error('computeSensitivityJacobian:badR', 'R must be a finite, positive scalar.');
end
if ~(isnumeric(epsilon) && isscalar(epsilon) && isfinite(epsilon) && epsilon > 0 && epsilon < 1)
    error('computeSensitivityJacobian:badEpsilon', ...
        'epsilon must be a finite scalar strictly between 0 and 1 (it is a RELATIVE step, deltaP = epsilon*p0).');
end

%% 'Strict' is reserved -- this function always sets it itself (see header)
isStrictName = cellfun(@(c) (ischar(c) || isstring(c)) && strcmpi(char(c), 'Strict'), varargin);
if any(isStrictName)
    error('computeSensitivityJacobian:strictReserved', ...
        ['''Strict'' cannot be passed to computeSensitivityJacobian -- every modalFeatureVector call ' ...
         'this function makes always uses ''Strict'',true by design, so that a perturbation failing ' ...
         'its validity check fails the whole Jacobian rather than silently degrading to a warning.']);
end

paramNames = {'beta', 'gamma', 'R'};
p0 = [beta; gamma; R];

%% Baseline (also always Strict -- a baseline that is itself invalid
% must not be used to build a Jacobian either)
[f0, details0] = modalFeatureVector(p0(1), p0(2), p0(3), varargin{:}, 'Strict', true);
nFeat = numel(f0);

%% Perturbed evaluations, one parameter at a time
deltaP = epsilon * p0;
fPlus  = nan(nFeat, 3);
fMinus = nan(nFeat, 3);
detailsPlus  = cell(1, 3);
detailsMinus = cell(1, 3);

for j = 1:3
    pPlus  = p0; pPlus(j)  = pPlus(j)  + deltaP(j);
    pMinus = p0; pMinus(j) = pMinus(j) - deltaP(j);

    if pMinus(j) <= 0
        % Unreachable given epsilon in (0,1) and p0(j) > 0 checked above
        % (pMinus(j) = p0(j)*(1-epsilon) > 0 always) -- kept as a
        % defensive backstop rather than removed, consistent with this
        % project's fail-at-the-source convention (see
        % modalFeatureVector.m's nonFiniteField check).
        error('computeSensitivityJacobian:nonPositivePerturbation', ...
            ['perturbing %s downward by epsilon=%.3g (deltaP=%.6g) would make it non-positive ' ...
             '(p0(%s)=%.6g) -- reduce epsilon.'], paramNames{j}, epsilon, deltaP(j), paramNames{j}, p0(j));
    end

    [fPlus(:,j), detailsPlus{j}]   = modalFeatureVector(pPlus(1),  pPlus(2),  pPlus(3),  varargin{:}, 'Strict', true);
    [fMinus(:,j), detailsMinus{j}] = modalFeatureVector(pMinus(1), pMinus(2), pMinus(3), varargin{:}, 'Strict', true);
end

%% Assemble (pure arithmetic, unit-tested separately -- see header)
J = centralDiffJacobian(fPlus, fMinus, deltaP);

info.p0 = p0;
info.paramNames = paramNames;
info.epsilon = epsilon;
info.deltaP = deltaP;
info.f0 = f0;
info.fPlus = fPlus;
info.fMinus = fMinus;
info.details0 = details0;
info.detailsPlus = detailsPlus;
info.detailsMinus = detailsMinus;

end