function result = solveInverse(forwardFcn, fObs, pRef, x0, lb, ub, varargin)
%SOLVEINVERSE  Bounded nonlinear least-squares parameter inversion,
%   min_x || r(x) ||^2  with  r = inverseResidual(x, pRef, fObs, forwardFcn),
%   solved by lsqnonlin in SCALED parameters x = p ./ pRef. Generic: knows
%   nothing about the EMM or SSI, only about the forward map p -> f.
%
% USAGE
%   result = solveInverse(forwardFcn, fObs, pRef, x0, lb, ub)
%   result = solveInverse(..., Name, Value, ...)
%
% The hidden truth is deliberately NOT an input. solveInverse only sees
% what a real inversion would see (the observed features); recovery errors
% against a known truth are computed by the calling twin-test script.
%
% INPUTS
%   forwardFcn : f = forwardFcn(p), p the UNSCALED n-vector. For the EMM
%                twin, e.g.
%                  @(p) modalFeatureVector(p(1), p(2), p(3), 'Strict', true)
%                Must be DETERMINISTIC for the whole solve. Any random
%                forcing phases / injected noise (SSI, JONSWAP) must be
%                drawn ONCE outside forwardFcn and held fixed; redrawing
%                per call makes the objective change between evaluations.
%   fObs       : m-vector of observed features (the inversion target).
%   pRef       : n-vector of positive reference values; x = p./pRef.
%   x0, lb, ub : n-vectors in SCALED units. Requires 0 < lb < ub and
%                lb <= x0 <= ub. lb must be strictly positive:
%                inverseResidual rejects x <= 0, so a zero lower bound
%                would let the solver probe an invalid point and fail on
%                a bounds mistake rather than on physics.
%
% NAME-VALUE OPTIONS (defaults in brackets)
%   'FeatureMask'              ([] = all)  forwarded to inverseResidual.
%   'ExpectedErrorIDs'         ({})  error identifiers that mean "the
%                              solver stepped outside the region where the
%                              forward map is valid", e.g.
%                              {'modalFeatureVector:invalidFeatureFit'}.
%                              Such an error ends the solve but is RECORDED
%                              (status 'leftValidRegion'), not thrown. Any
%                              other error is re-thrown unchanged, so a bug
%                              can never be mistaken for a finding. The
%                              error's cause chain is searched too, since
%                              lsqnonlin may wrap a failure in its first
%                              objective evaluation.
%   'FiniteDifferenceStepSize' (1e-4)  see STEP SIZE below.
%   'FiniteDifferenceType'     ('central')  'central' or 'forward'.
%   'FunctionTolerance'        (1e-10)
%   'StepTolerance'            (1e-10)
%   'OptimalityTolerance'      (1e-10)
%   'MaxIterations'            (100)
%   'MaxFunctionEvaluations'   (1000)  counts forward-model calls; with
%                              central differences each Jacobian costs 2n
%                              calls, so budget roughly (2n+1) per iteration.
%   'Display'                  ('off')  passed to lsqnonlin ('iter' to watch).
%
%   The tolerances are much tighter than lsqnonlin's defaults (1e-6) on
%   purpose: in an exact twin test the solver should run until it genuinely
%   stalls, so that any remaining error (particularly along the weak gamma
%   direction) reflects the problem, not an early stopping rule. If the
%   forward map's own numerical noise floor is reached first, lsqnonlin
%   simply stops on a step/function tolerance (exitflag 2 or 3).
%
% STEP SIZE
%   lsqnonlin's DEFAULT finite-difference step is ~sqrt(eps) = 1.5e-8
%   (forward), far below where the EMM feature map's derivatives were
%   validated (Section 8 of the Jacobian report: 1e-3 to 1e-5, with signs
%   of a numerical floor for gamma already at 1e-5). Left at the default,
%   the solver would steer on derivatives dominated by solver noise. MATLAB
%   uses delta_j = v * max(|x_j|, TypicalX_j) with TypicalX = 1 by default,
%   so with v = 1e-4 and x ~ 1 the step is ~1e-4 in scaled units, i.e. a
%   ~1e-4 relative step in p: the validated setting.
%
% OUTPUT: result struct
%   status       : 'solved'           exitflag > 0 (a tolerance was met; this
%                                     does NOT by itself prove the truth was
%                                     reached, check resnorm/recovery)
%                  'budgetExhausted'  exitflag 0 (MaxIterations or
%                                     MaxFunctionEvaluations hit)
%                  'solverStopped'    exitflag < 0
%                  'leftValidRegion'  an ExpectedErrorIDs error occurred
%   exitflag, message        : from lsqnonlin (NaN / '' if leftValidRegion)
%   xHat, pHat               : n-by-1 solution, scaled / unscaled (NaN if
%                              leftValidRegion)
%   resnorm, residual        : ||r||^2 and r at xHat (NaN if leftValidRegion)
%   jacobian                 : dr/dx at xHat from lsqnonlin (full matrix),
%                              i.e. diag(1/fObs) * (df/dp) * diag(pRef) on
%                              the selected rows. At x = [1;1;1] with
%                              fObs = f(pRef) this is the scaled Jacobian
%                              J_s of the identifiability analysis.
%   nIterations              : lsqnonlin iterations (NaN if leftValidRegion)
%   funcCount                : forward-model calls made (counted here, so
%                              also available when the solve was aborted)
%   runtimeSec               : wall-clock time of the solve
%   history                  : per-iteration log from the OutputFcn:
%                              .iteration (k-by-1), .x (k-by-n), .p (k-by-n),
%                              .resnorm (k-by-1), .funccount (k-by-1)
%   failedX, failedP         : the point whose evaluation raised the
%                              expected error ([] otherwise)
%   error                    : [] or struct(identifier, message,
%                              causeIdentifiers) for the recorded error;
%                              causeIdentifiers lists the identifiers of
%                              the error itself and everything in its
%                              cause chain, top-level first
%   x0, lb, ub, pRef, fObs, featureMask, settings : inputs as used
%                              (settings = the option values above, as a
%                              plain struct, so results files are
%                              self-describing and portable)
%   solverOptions            : the optimoptions object actually passed

%% ---------------- Parse options --------------------------------------
s = struct( ...
    'FeatureMask', [], ...
    'ExpectedErrorIDs', {{}}, ...
    'FiniteDifferenceStepSize', 1e-4, ...
    'FiniteDifferenceType', 'central', ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-10, ...
    'MaxIterations', 100, ...
    'MaxFunctionEvaluations', 1000, ...
    'Display', 'off');
names = fieldnames(s);

if mod(numel(varargin), 2) ~= 0
    error('solveInverse:badOption', 'Name-value options must come in pairs.');
end
for k = 1:2:numel(varargin)
    nm = varargin{k};
    if ~(ischar(nm) || (isstring(nm) && isscalar(nm)))
        error('solveInverse:badOption', 'Option names must be text.');
    end
    hit = strcmpi(char(nm), names);
    if ~any(hit)
        error('solveInverse:badOption', 'Unknown option ''%s''.', char(nm));
    end
    s.(names{hit}) = varargin{k+1};
end

%% ---------------- Validate ------------------------------------------
if ~isa(forwardFcn, 'function_handle')
    error('solveInverse:badForwardFcn', 'forwardFcn must be a function handle, f = forwardFcn(p).');
end
if ~(isnumeric(pRef) && isvector(pRef) && isreal(pRef) && all(isfinite(pRef)) && all(pRef > 0))
    error('solveInverse:badPRef', 'pRef must be a real vector of finite, strictly positive values.');
end
n = numel(pRef);
pRef = pRef(:);

vecOK = @(v) isnumeric(v) && isvector(v) && isreal(v) && all(isfinite(v)) && numel(v) == n;
if ~vecOK(x0), error('solveInverse:badX0', 'x0 must be a real, finite vector with %d elements.', n); end
if ~vecOK(lb), error('solveInverse:badBounds', 'lb must be a real, finite vector with %d elements.', n); end
if ~vecOK(ub), error('solveInverse:badBounds', 'ub must be a real, finite vector with %d elements.', n); end
x0 = x0(:); lb = lb(:); ub = ub(:);
if any(lb <= 0)
    error('solveInverse:badBounds', ...
        'lb must be strictly positive (scaled parameters must stay > 0 for inverseResidual).');
end
if any(lb >= ub)
    error('solveInverse:badBounds', 'Every lb(j) must be strictly less than ub(j).');
end
if any(x0 < lb | x0 > ub)
    error('solveInverse:badX0', 'x0 must lie within [lb, ub].');
end

if ~(isnumeric(fObs) && isvector(fObs) && isreal(fObs) && all(isfinite(fObs)))
    error('solveInverse:badFObs', 'fObs must be a real vector of finite values.');
end
fObs = fObs(:);
m = numel(fObs);

% Number of features that will actually enter the residual. (Full mask
% validation is inverseResidual's job; this is only the count.)
mask = s.FeatureMask;
if isempty(mask) && ~islogical(mask)
    nSel = m;
elseif islogical(mask)
    nSel = nnz(mask);
elseif isnumeric(mask)
    nSel = numel(unique(mask(:)));
else
    error('solveInverse:badOption', 'FeatureMask must be a logical vector or a vector of indices.');
end
if nSel < n
    error('solveInverse:underdetermined', ...
        ['Only %d feature(s) selected for %d parameters. lsqnonlin''s trust-region-reflective ' ...
         'algorithm needs at least as many residuals as unknowns, and the parameters could not ' ...
         'all be identified anyway.'], nSel, n);
end

ids = s.ExpectedErrorIDs;
if ischar(ids) || (isstring(ids) && isscalar(ids))
    ids = {char(ids)};
end
if ~(iscell(ids) && all(cellfun(@(c) ischar(c) && ~isempty(c), ids)))
    error('solveInverse:badOption', 'ExpectedErrorIDs must be a cell array of non-empty error identifiers.');
end
s.ExpectedErrorIDs = ids(:).';

if ~(isnumeric(s.FiniteDifferenceStepSize) && isscalar(s.FiniteDifferenceStepSize) && ...
        isfinite(s.FiniteDifferenceStepSize) && s.FiniteDifferenceStepSize > 0)
    error('solveInverse:badOption', 'FiniteDifferenceStepSize must be a finite positive scalar.');
end
if ~(ischar(s.FiniteDifferenceType) && any(strcmpi(s.FiniteDifferenceType, {'central', 'forward'})))
    error('solveInverse:badOption', 'FiniteDifferenceType must be ''central'' or ''forward''.');
end
s.FiniteDifferenceType = lower(s.FiniteDifferenceType);

%% ---------------- Shared state (nested functions below) --------------
evalCount = 0;
lastX = [];
histIter = zeros(0, 1);
histX = zeros(0, n);
histResnorm = zeros(0, 1);
histFunccount = zeros(0, 1);

opts = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'FiniteDifferenceType', s.FiniteDifferenceType, ...
    'FiniteDifferenceStepSize', s.FiniteDifferenceStepSize, ...
    'FunctionTolerance', s.FunctionTolerance, ...
    'StepTolerance', s.StepTolerance, ...
    'OptimalityTolerance', s.OptimalityTolerance, ...
    'MaxIterations', s.MaxIterations, ...
    'MaxFunctionEvaluations', s.MaxFunctionEvaluations, ...
    'Display', s.Display, ...
    'OutputFcn', @recordIteration);

%% ---------------- Solve -----------------------------------------------
result = struct();
result.status = '';
result.exitflag = NaN;
result.message = '';
result.xHat = NaN(n, 1);
result.pHat = NaN(n, 1);
result.resnorm = NaN;
result.residual = NaN(nSel, 1);
result.jacobian = NaN(nSel, n);
result.nIterations = NaN;
result.funcCount = NaN;
result.runtimeSec = NaN;
result.history = [];
result.failedX = [];
result.failedP = [];
result.error = [];

t0 = tic;
try
    [xHat, resnorm, residual, exitflag, output, ~, jac] = ...
        lsqnonlin(@objective, x0, lb, ub, opts);

    result.exitflag = exitflag;
    if isfield(output, 'message'), result.message = output.message; end
    result.xHat = xHat(:);
    result.pHat = xHat(:) .* pRef;
    result.resnorm = resnorm;
    result.residual = residual(:);
    result.jacobian = full(jac);
    if isfield(output, 'iterations'), result.nIterations = output.iterations; end

    if exitflag > 0
        result.status = 'solved';
    elseif exitflag == 0
        result.status = 'budgetExhausted';
    else
        result.status = 'solverStopped';
    end
catch err
    [isExpected, causeIDs] = matchErrorChain(err, s.ExpectedErrorIDs);
    if ~isExpected
        rethrow(err);
    end
    result.status = 'leftValidRegion';
    result.failedX = lastX;
    result.failedP = lastX .* pRef;
    result.error = struct('identifier', err.identifier, 'message', err.message, ...
        'causeIdentifiers', {causeIDs});
end
result.runtimeSec = toc(t0);
result.funcCount = evalCount;

result.history = struct( ...
    'iteration', histIter, ...
    'x', histX, ...
    'p', histX .* pRef.', ...
    'resnorm', histResnorm, ...
    'funccount', histFunccount);

result.x0 = x0;
result.lb = lb;
result.ub = ub;
result.pRef = pRef;
result.fObs = fObs;
result.featureMask = s.FeatureMask;
result.settings = s;
result.solverOptions = opts;

%% ---------------- Nested functions -----------------------------------
    function r = objective(x)
        % Counts forward calls and remembers the point being evaluated, so
        % that if forwardFcn throws, we know exactly where it failed.
        x = x(:);
        lastX = x;
        evalCount = evalCount + 1;
        r = inverseResidual(x, pRef, fObs, forwardFcn, 'FeatureMask', s.FeatureMask);
    end

    function stop = recordIteration(x, optimValues, state)
        stop = false;
        if strcmp(state, 'iter')
            histIter(end+1, 1) = getOr(optimValues, 'iteration', NaN);
            histX(end+1, :) = x(:).';
            histResnorm(end+1, 1) = getOr(optimValues, 'resnorm', NaN);
            histFunccount(end+1, 1) = getOr(optimValues, 'funccount', NaN);
        end
    end
end

%% ------------------------------------------------------------------------
function v = getOr(st, field, default)
if isstruct(st) && isfield(st, field)
    v = st.(field);
else
    v = default;
end
end

function [tf, allIDs] = matchErrorChain(err, ids)
% True if err, or anything in its cause chain, has an identifier in ids.
% Works with MATLAB MException objects (which have .cause) and Octave's
% struct-valued caught errors (which do not).
allIDs = {};
queue = {err};
while ~isempty(queue)
    e = queue{1};
    queue(1) = [];
    allIDs{end+1} = e.identifier; %#ok<AGROW>
    try
        c = e.cause;
    catch
        c = {};
    end
    if iscell(c)
        queue = [queue, c(:).']; %#ok<AGROW>
    end
end
tf = any(ismember(allIDs, ids));
end