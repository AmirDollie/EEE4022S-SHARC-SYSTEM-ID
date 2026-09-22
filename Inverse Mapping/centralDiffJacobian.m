function J = centralDiffJacobian(fPlus, fMinus, deltaP)
%CENTRALDIFFJACOBIAN  Pure central-difference Jacobian assembly.
%   J = centralDiffJacobian(fPlus, fMinus, deltaP)
%
%   fPlus, fMinus : nFeat x nParam. Column j of each is the feature
%                   vector evaluated with the j-th parameter perturbed
%                   by +deltaP(j) / -deltaP(j) respectively, all other
%                   parameters held fixed at the reference point.
%   deltaP        : nParam x 1 (or 1 x nParam) -- the ACTUAL perturbation
%                   used per parameter (already scaled, e.g.
%                   epsilon*p0(j) -- not epsilon itself).
%
%   J             : nFeat x nParam,
%                   J(:,j) = (fPlus(:,j) - fMinus(:,j)) ./ (2*deltaP(j)).
%
% Deliberately has no idea what beta/gamma/R, features, or epsilon are,
% and makes no modalFeatureVector/EMM call itself. Factored out of
% computeSensitivityJacobian.m specifically so the finite-difference
% ARITHMETIC (column indexing, the 2*deltaP normalisation) can be
% unit-tested directly against known synthetic numbers, without paying
% for any EMM solves -- see computeSensitivityJacobianTester.m. This
% mirrors why fitFeatureFrom3Points.m was factored out of
% modalFeatureVector.m for the same reason.

if ~(isnumeric(fPlus) && isnumeric(fMinus))
    error('centralDiffJacobian:badInput', 'fPlus and fMinus must both be numeric.');
end
if ~isequal(size(fPlus), size(fMinus))
    error('centralDiffJacobian:sizeMismatch', ...
        'fPlus is %s but fMinus is %s -- they must be the same size.', ...
        mat2str(size(fPlus)), mat2str(size(fMinus)));
end

nParam = size(fPlus, 2);
deltaP = deltaP(:).';
if numel(deltaP) ~= nParam
    error('centralDiffJacobian:deltaPSizeMismatch', ...
        'deltaP has %d element(s), expected %d (one per column of fPlus/fMinus).', ...
        numel(deltaP), nParam);
end
if any(~isfinite(deltaP)) || any(deltaP <= 0)
    error('centralDiffJacobian:badDeltaP', ...
        'deltaP must be all finite and strictly positive, got [%s].', mat2str(deltaP));
end

J = (fPlus - fMinus) ./ (2 * deltaP);
end