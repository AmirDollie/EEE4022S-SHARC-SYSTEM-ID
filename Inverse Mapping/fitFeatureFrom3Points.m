function [omegaStar, degenerate, extrapolated, kindMismatch, a, b, c] = ...
    fitFeatureFrom3Points(omega3, y3, kind)
%FITFEATUREFROM3POINTS  Pure sub-grid feature-location decision logic.
%   Given three (omega, y) points bracketing a feature (y = |A_{n,j}|^2
%   in modalFeatureVector.m, but this function has no idea about that --
%   it only ever sees numbers), fits the exact quadratic through them
%   (via quadFit3.m) and decides on a sub-grid feature location, plus the
%   three diagnostic flags modalFeatureVector.m surfaces per feature.
%
%   Factored out of modalFeatureVector.m specifically so the
%   degenerate / extrapolated / kindMismatch decision logic can be unit
%   tested directly against known synthetic (omega3,y3) data, without
%   paying for an EMM solve -- see modalFeatureVectorTester.m.
%
%   [omegaStar, degenerate, extrapolated, kindMismatch, a, b, c] = ...
%       fitFeatureFrom3Points(omega3, y3, kind)
%
%   omega3, y3 : each 1x3 (or 3x1), the three anchor points and their
%                y-values (three distinct omega3 values required).
%   kind       : 'max' or 'zero' -- which extremum this feature is.
%
%   omegaStar     : the sub-grid feature location. Equal to the fitted
%                   parabola's vertex -b/(2a) unless degenerate is true,
%                   in which case it falls back to whichever of the
%                   three anchor omegas has the extreme (max or min, per
%                   `kind`) y-value -- NOT a sub-grid estimate.
%   degenerate    : true if the three points are numerically too flat to
%                   fit a usable parabola (|a| negligible relative to
%                   the curvature scale needed to explain the y-span).
%   extrapolated  : true if the fitted vertex falls strictly outside the
%                   anchor bracket [min(omega3), max(omega3)]. Always
%                   false when degenerate is true (the fallback value is
%                   one of the anchors by construction).
%   kindMismatch  : true if the fitted parabola opens the "wrong" way
%                   for the stated kind (a>0 for 'max', a<0 for 'zero').
%                   Always false when degenerate is true.
%   a, b, c       : the fitted quadratic's coefficients, from quadFit3.m.

if ~(ischar(kind) || isstring(kind)) || ~ismember(char(kind), {'max', 'zero'})
    if ischar(kind) || isstring(kind)
        badKindDesc = char(kind);
    else
        badKindDesc = sprintf('<%s>', class(kind)); % avoids MATLAB's string() ctor, unsupported in Octave
    end
    error('fitFeatureFrom3Points:badKind', ...
        'kind must be ''max'' or ''zero'', got ''%s''.', badKindDesc);
end
omega3 = omega3(:).'; y3 = y3(:).';
if numel(omega3) ~= 3 || numel(y3) ~= 3
    error('fitFeatureFrom3Points:badInput', 'omega3 and y3 must each have exactly 3 elements.');
end

[a, b, c] = quadFit3(omega3, y3);

omegaLo = min(omega3); omegaHi = max(omega3);
yspan = max(y3) - min(y3);
xspan = omegaHi - omegaLo;
degenerate = (yspan == 0) || (abs(a) * xspan^2 < 1e-10 * max(yspan, eps));

if degenerate
    if strcmp(kind, 'max')
        [~, mi] = max(y3);
    else
        [~, mi] = min(y3);
    end
    omegaStar = omega3(mi);
    extrapolated = false;
    kindMismatch = false;
else
    omegaStar = -b / (2 * a);
    extrapolated = (omegaStar < omegaLo) || (omegaStar > omegaHi);
    kindMismatch = (strcmp(kind, 'max') && a > 0) || (strcmp(kind, 'zero') && a < 0);
end

end