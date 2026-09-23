function r = inverseResidual(x, pRef, fObs, forwardFcn, varargin)
%INVERSERESIDUAL  Relative feature residual for the parameter inversion,
%   evaluated at SCALED parameters x. Generic: knows nothing about the EMM
%   or SSI, only about a supplied forward map p -> f.
%
% USAGE
%   r = inverseResidual(x, pRef, fObs, forwardFcn)
%   r = inverseResidual(x, pRef, fObs, forwardFcn, 'FeatureMask', mask)
%
% DEFINITION
%   p = x .* pRef                          (unscale, elementwise)
%   f = forwardFcn(p)                      (predicted feature vector)
%   r = (f(sel) - fObs(sel)) ./ fObs(sel)  (relative error, selected features)
%
% so that r is O(feature relative error) and dimensionless, and x is O(1)
% in every component (x = [1;1;1] is the reference point pRef). This is
% the function handed to lsqnonlin by solveInverse.m, as
%   @(x) inverseResidual(x, pRef, fObs, forwardFcn, ...)
%
% INPUTS
%   x          : n-vector of scaled parameters (e.g. [beta/beta0;
%                gamma/gamma0; R/R0]). Finite and strictly positive: the
%                relative scaling only makes sense for positive physical
%                parameters, and a non-positive trial value is always a
%                bounds mistake, never something to pass to the EMM.
%   pRef       : n-vector of reference parameter values used for scaling.
%                Finite and strictly positive. Same length as x.
%   fObs       : m-vector of observed (target) features, same length and
%                order as forwardFcn's output. Must be finite everywhere
%                and non-zero on the selected features (it is a divisor).
%   forwardFcn : function handle, f = forwardFcn(p), taking the UNSCALED
%                n-vector p (column) and returning an m-vector f. For the
%                EMM twin, e.g.
%                  @(p) modalFeatureVector(p(1), p(2), p(3), 'Strict', true)
%                forwardFcn MUST be deterministic for the whole solve: if
%                a later forward map involves random forcing phases or
%                injected noise (SSI/JONSWAP), that realisation must be
%                fixed once, outside forwardFcn, not redrawn per call.
%                Otherwise the objective itself changes between
%                evaluations and the optimiser is chasing a moving target.
%
% NAME-VALUE OPTIONS
%   'FeatureMask' : which features enter the residual. Either a logical
%                   m-vector, or a vector of distinct integer indices in
%                   1..m. Default: all m features. Index masks are
%                   converted to a logical mask, so r is ALWAYS returned
%                   in the original feature order, regardless of the
%                   order the indices were given in. At least one feature
%                   must be selected. (Whether enough features are
%                   selected to identify all n parameters is NOT checked
%                   here; that is a property of the inverse problem, not
%                   of the residual.)
%
% OUTPUT
%   r : column vector, one entry per selected feature.
%
% ERRORS
%   Errors thrown BY forwardFcn itself (e.g. modalFeatureVector's
%   'modalFeatureVector:invalidFeatureFit' under 'Strict',true) are NOT
%   caught or wrapped here: they propagate with their original identifier,
%   so solveInverse.m can decide which identifiers count as "left the
%   valid region" and which are genuine bugs.
%
%   This function's own checks throw:
%     inverseResidual:badX, :badPRef, :sizeMismatch, :badFObs,
%     :badForwardFcn, :badMask, :badOption,
%     :forwardSizeMismatch  (forwardFcn returned the wrong number of features)
%     :nonFiniteForward     (forwardFcn returned NaN/Inf on a SELECTED feature)
%
%   Non-finite forward output on UNSELECTED features is deliberately
%   allowed: a later SSI forward map may legitimately return NaN for a
%   mode it failed to recover, and masking that feature out should be
%   enough to proceed.

%% Options
mask = [];
if mod(numel(varargin), 2) ~= 0
    error('inverseResidual:badOption', 'Name-value options must come in pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~(ischar(name) || (isstring(name) && isscalar(name)))
        error('inverseResidual:badOption', 'Option names must be text.');
    end
    switch lower(char(name))
        case 'featuremask'
            mask = varargin{k+1};
        otherwise
            error('inverseResidual:badOption', 'Unknown option ''%s''.', char(name));
    end
end

%% Validate inputs (all before calling forwardFcn: fail cheaply)
if ~(isnumeric(pRef) && isvector(pRef) && isreal(pRef) && all(isfinite(pRef)) && all(pRef > 0))
    error('inverseResidual:badPRef', 'pRef must be a real vector of finite, strictly positive values.');
end
if ~(isnumeric(x) && isvector(x) && isreal(x) && all(isfinite(x)) && all(x > 0))
    error('inverseResidual:badX', 'x must be a real vector of finite, strictly positive scaled parameters.');
end
if numel(x) ~= numel(pRef)
    error('inverseResidual:sizeMismatch', 'x has %d elements but pRef has %d.', numel(x), numel(pRef));
end
if ~(isnumeric(fObs) && isvector(fObs) && isreal(fObs) && all(isfinite(fObs)))
    error('inverseResidual:badFObs', 'fObs must be a real vector of finite values.');
end
if ~isa(forwardFcn, 'function_handle')
    error('inverseResidual:badForwardFcn', 'forwardFcn must be a function handle, f = forwardFcn(p).');
end

m = numel(fObs);
sel = resolveMask(mask, m);

fObs = fObs(:);
if any(fObs(sel) == 0)
    error('inverseResidual:badFObs', ...
        'fObs is zero on a selected feature (index %d); the relative residual divides by fObs.', ...
        find(sel & fObs == 0, 1));
end

%% Forward evaluation (errors from forwardFcn propagate unchanged)
p = x(:) .* pRef(:);
f = forwardFcn(p);

if ~(isnumeric(f) && isvector(f) && numel(f) == m)
    error('inverseResidual:forwardSizeMismatch', ...
        'forwardFcn returned a %s %s, expected a numeric vector with %d elements (matching fObs).', ...
        mat2str(size(f)), class(f), m);
end
f = f(:);
if ~all(isfinite(f(sel)))
    error('inverseResidual:nonFiniteForward', ...
        'forwardFcn returned a non-finite value on selected feature %d at p = %s.', ...
        find(sel & ~isfinite(f), 1), mat2str(p.', 8));
end

%% Residual
r = (f(sel) - fObs(sel)) ./ fObs(sel);

end

%% ------------------------------------------------------------------------
function sel = resolveMask(mask, m)
% Convert a logical or index mask into a logical m-by-1 selection vector.
if isempty(mask) && ~islogical(mask)
    sel = true(m, 1);   % default: all features ([] means "not given")
    return
end
if islogical(mask)
    if ~(isvector(mask) && numel(mask) == m)
        error('inverseResidual:badMask', ...
            'A logical FeatureMask must have exactly %d elements (one per feature), got %d.', m, numel(mask));
    end
    sel = mask(:);
elseif isnumeric(mask) && isvector(mask) && isreal(mask)
    if ~all(isfinite(mask) & mask == round(mask) & mask >= 1 & mask <= m)
        error('inverseResidual:badMask', ...
            'An index FeatureMask must contain integers in 1..%d.', m);
    end
    if numel(unique(mask)) ~= numel(mask)
        error('inverseResidual:badMask', 'An index FeatureMask must not contain duplicate indices.');
    end
    sel = false(m, 1);
    sel(mask) = true;
else
    error('inverseResidual:badMask', 'FeatureMask must be a logical vector or a vector of indices.');
end
if ~any(sel)
    error('inverseResidual:badMask', 'FeatureMask selects no features.');
end
end