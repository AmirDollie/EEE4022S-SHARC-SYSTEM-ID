function out = featureInformation(Jf, Sf, Jr, Sr, varargin)
%FEATUREINFORMATION How much of a reference observable's local parameter
%information a candidate feature set retains.
%
%   OUT = FEATUREINFORMATION(JF, SF, JR, SR)
%   OUT = FEATUREINFORMATION(..., Name, Value)
%
%   JF - mF x nP real Jacobian of the candidate features w.r.t. the
%        parameters theta (e.g. [ln beta, ln R]).
%   SF - covariance of the candidate features: mF x mF SPD matrix, an
%        mF-vector of variances (diagonal), or a scalar (iid).
%   JR, SR - the same for the REFERENCE observables (e.g. the complex
%        sensor FRFs as [Re; Im] rows), mR x nP and mR x mR / vector / scalar.
%
%   Local (Gauss-Newton / Fisher) information:
%       F_f = JF' SF^-1 JF,       F_r = JR' SR^-1 JR.
%   Retention is measured by the generalised eigenproblem
%       F_f v = lambda F_r v,
%   solved symmetrically via F_r = L L' (Cholesky): lambda are the
%   eigenvalues of L^-1 F_f L^-T. lambda_i is the fraction of the reference
%   information kept along direction v_i (1 = nothing lost, 0 = everything
%   lost). Information is matrix-valued, so no single percentage is
%   reported; the eigenpairs plus per-parameter stds are the summary.
%
%   CONSISTENCY: if the features are a (locally linear) function of the
%   reference observables, f ~ A y, and SF is propagated consistently,
%   SF = A SR A', then 0 <= lambda <= 1 (no processing can add
%   information). lambda > 1 + 'Tol' means the two noise models are
%   inconsistent (e.g. feature noise understated); OUT.consistent is then
%   false and the warning featureInformation:lambdaAboveOne is raised.
%   Rescaling features (JF -> D JF, SF -> D SF D') leaves lambda unchanged.
%
%   OUT fields
%     Ff, Fr            nP x nP information matrices
%     lambda            nP x 1, descending
%     V                 nP x nP generalised eigenvectors (columns, unit norm)
%     Vscaled           V expressed in units of the reference stds
%                       (V ./ stdR, columns renormalised): which parameters a
%                       direction involves RELATIVE TO HOW WELL THE REFERENCE
%                       KNOWS THEM. Raw V is dominated by poorly known
%                       parameters simply because a unit change in them
%                       carries little information; read Vscaled instead.
%     stdF, stdR        joint stds sqrt(diag(inv(F))) (all parameters free)
%     fixedStdF, fixedStdR   1 ./ sqrt(diag(F)) (the other parameters known)
%     stdRatio          stdF ./ stdR (>= 1 when consistent)
%     detRatio          det(F_f) / det(F_r) = prod(lambda)
%     consistent        all(lambda <= 1 + Tol)
%     paramNames
%
%   NAME-VALUE OPTIONS
%     'ParamNames'    cellstr, default {'p1', 'p2', ...}
%     'Tol'           tolerance on lambda <= 1, default 1e-8
%     'WarnInconsistent'  true (default) / false
%
%   ERROR IDENTIFIERS (featureInformation:<id>)
%     badJ, sizeMismatch, badSigma, sigmaNotPositiveDefinite,
%     referenceNotInformative, badOption
%
%   Lives in Inverse Mapping/. Tester:
%   Inverse Mapping/Inverse Mapping Unit Tests/featureInformationTester.m

%% ---- Options ---------------------------------------------------------------
opt = struct('ParamNames', {{}}, 'Tol', 1e-8, 'WarnInconsistent', true);
if mod(numel(varargin), 2) ~= 0
    error('featureInformation:badOption', 'Options must be name-value pairs.');
end
names = fieldnames(opt);
for a = 1:2:numel(varargin)
    key = varargin{a};
    if isa(key, 'string'), key = char(key); end
    hit = ischar(key) & strcmpi(key, names);
    if ~any(hit)
        error('featureInformation:badOption', 'Unknown option.');
    end
    opt.(names{hit}) = varargin{a + 1};
end
if ~(isnumeric(opt.Tol) && isscalar(opt.Tol) && isfinite(opt.Tol) && opt.Tol >= 0)
    error('featureInformation:badOption', 'Tol must be a nonnegative scalar.');
end

%% ---- Inputs --------------------------------------------------------------------
checkJ(Jf, 'JF'); checkJ(Jr, 'JR');
nP = size(Jf, 2);
if size(Jr, 2) ~= nP
    error('featureInformation:sizeMismatch', 'JF and JR must have the same number of columns (parameters).');
end
if isempty(opt.ParamNames)
    opt.ParamNames = arrayfun(@(k) sprintf('p%d', k), 1:nP, 'UniformOutput', false);
elseif ~(iscellstr(opt.ParamNames) && numel(opt.ParamNames) == nP) %#ok<ISCLSTR>
    error('featureInformation:badOption', 'ParamNames must be a cellstr with one name per parameter.');
end

Ff = infoMatrix(Jf, Sf, 'SF');
Fr = infoMatrix(Jr, Sr, 'SR');

[L, flag] = chol(Fr, 'lower');
if flag ~= 0 || rcond(Fr) < 1e-14
    error('featureInformation:referenceNotInformative', ...
        'The reference information matrix is singular: the reference does not identify all %d parameters.', nP);
end

%% ---- Generalised eigenproblem ------------------------------------------------------------
M = L \ (Ff / L.');                       % L^-1 Ff L^-T
M = (M + M.') / 2;
[U, D] = eig(M);
[lambda, ord] = sort(diag(D), 'descend');
U = U(:, ord);
V = L.' \ U;                              % F_f V = F_r V diag(lambda)
V = V ./ vecnorm(V, 2, 1);
for k = 1:nP                              % sign convention: largest component positive
    [~, im] = max(abs(V(:, k)));
    if V(im, k) < 0, V(:, k) = -V(:, k); end
end

consistent = all(lambda <= 1 + opt.Tol);
if ~consistent && opt.WarnInconsistent
    warning('featureInformation:lambdaAboveOne', '%s', sprintf( ...
        'max lambda = %.6g > 1: the feature and reference noise models are inconsistent.', max(lambda)));
end

%% ---- Summary -----------------------------------------------------------------------------
out = struct();
out.Ff = Ff;
out.Fr = Fr;
out.lambda = lambda;
out.V = V;
out.stdR = sqrt(diag(inv(Fr)));
Vs = V ./ out.stdR;
Vs = Vs ./ vecnorm(Vs, 2, 1);
for k = 1:nP
    [~, im] = max(abs(Vs(:, k)));
    if Vs(im, k) < 0, Vs(:, k) = -Vs(:, k); end
end
out.Vscaled = Vs;
if rcond(Ff) > 1e-14
    out.stdF = sqrt(diag(inv(Ff)));
else
    out.stdF = Inf(nP, 1);                % a direction is completely lost
end
out.fixedStdF = 1 ./ sqrt(diag(Ff));
out.fixedStdR = 1 ./ sqrt(diag(Fr));
out.stdRatio = out.stdF ./ out.stdR;
out.detRatio = prod(lambda);
out.consistent = consistent;
out.paramNames = opt.ParamNames;
end

%% ================================================================================================
function checkJ(J, name)
if ~(isnumeric(J) && isreal(J) && ismatrix(J) && ~isempty(J) && all(isfinite(J(:))))
    error('featureInformation:badJ', '%s must be a nonempty finite real matrix.', name);
end
end

function F = infoMatrix(J, S, name)
m = size(J, 1);
if ~(isnumeric(S) && isreal(S) && all(isfinite(S(:))))
    error('featureInformation:badSigma', '%s must be real and finite.', name);
end
if isscalar(S)
    if S <= 0
        error('featureInformation:sigmaNotPositiveDefinite', '%s must be positive.', name);
    end
    F = (J.' * J) / S;
elseif isvector(S) && numel(S) == m
    if any(S(:) <= 0)
        error('featureInformation:sigmaNotPositiveDefinite', '%s variances must be positive.', name);
    end
    F = J.' * (J ./ S(:));
elseif isequal(size(S), [m m])
    if max(max(abs(S - S.'))) > 1e-10 * max(abs(S(:)))
        error('featureInformation:badSigma', '%s must be symmetric.', name);
    end
    [C, flag] = chol((S + S.') / 2, 'lower');
    if flag ~= 0
        error('featureInformation:sigmaNotPositiveDefinite', '%s is not positive definite.', name);
    end
    W = C \ J;                            % whitened Jacobian
    F = W.' * W;
else
    error('featureInformation:sizeMismatch', ...
        '%s must be a scalar, a %d-vector or a %d x %d matrix.', name, m, m, m);
end
F = (F + F.') / 2;
end