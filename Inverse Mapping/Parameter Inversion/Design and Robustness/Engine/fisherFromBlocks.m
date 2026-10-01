function q = fisherFromBlocks(J, blocks, bias, varargin)
%FISHERFROMBLOCKS  Fisher information, parameter covariance and bias projection with a
%   block-diagonal feature covariance (one 2m x 2m block per retained, decorrelated bin).
%   q = fisherFromBlocks(J, blocks)
%   q = fisherFromBlocks(J, blocks, bias)          bias: feature vector(s), one per column
%
%   F = sum_k J_k' Sigma_k^-1 J_k (Cholesky whitening per block)
%   q.F, covTheta = F^-1, sigma, corr (theta_1, theta_2), cond (kappa(F)), perBin (nTheta x nF:
%   diagonal information per bin), and for each bias column:
%   dtheta = F^-1 J' Sigma^-1 b,  d = sqrt(dtheta' F dtheta).
%
%   NON-IDENTIFIABILITY IS A RESULT, NOT AN ERROR. A block whose Cholesky factorisation fails
%   (not positive definite: a genuinely degenerate bin) is excluded and listed in q.degenerateBins;
%   it is never regularised, which could manufacture information. If F is singular or numerically
%   so (rcond(F) < 'RcondMin', default 1e-12, or no usable bins), q.singular = true and sigma =
%   Inf, cond = Inf, covTheta = Inf, dtheta = NaN, d = NaN.
%   q = fisherFromBlocks(J, blocks, bias, 'RcondMin', 1e-12)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Engine/.
if nargin < 3, bias = []; end
opt = struct('RcondMin', 1e-12);
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i + 1}; end
[p, ~, nF] = size(blocks);
if ndims(blocks) == 2, nF = 1; end
if isempty(blocks), nF = 0; end
nT = size(J, 2);
F = zeros(nT);
g = zeros(nT, max(1, size(bias, 2)));
perBin = zeros(nT, nF);
degenerate = false(1, nF);
for k = 1:nF
    ix = (k - 1) * p + (1:p);
    [Lk, flag] = chol(blocks(:, :, k), 'lower');
    if flag ~= 0, degenerate(k) = true; continue; end
    Jw = Lk \ J(ix, :);
    Fk = Jw.' * Jw;
    F = F + Fk;
    perBin(:, k) = diag(Fk);
    if ~isempty(bias)
        g = g + Jw.' * (Lk \ bias(ix, :));
    end
end
q.F = F; q.perBin = perBin; q.degenerateBins = find(degenerate);
q.nUsed = nF - nnz(degenerate);
q.singular = q.nUsed == 0 || ~all(isfinite(F(:))) || rcond(F) < opt.RcondMin;
nb = max(1, size(bias, 2));
if q.singular
    q.covTheta = Inf(nT); q.sigma = Inf(nT, 1); q.corr = NaN; q.cond = Inf;
    q.dtheta = NaN(nT, nb); q.d = NaN(1, nb);
    if isempty(bias), q.dtheta = []; q.d = []; end
    return
end
C = F \ eye(nT);
q.covTheta = C;
q.sigma = sqrt(diag(C));
q.corr = NaN; if nT >= 2, q.corr = C(1, 2) / sqrt(C(1, 1) * C(2, 2)); end
q.cond = cond(F);
q.dtheta = []; q.d = [];
if ~isempty(bias)
    q.dtheta = F \ g;
    q.d = sqrt(sum(q.dtheta .* (F * q.dtheta), 1));
end
end