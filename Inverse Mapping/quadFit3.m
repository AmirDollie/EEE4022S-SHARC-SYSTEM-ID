function [a, b, c] = quadFit3(x, y)
%QUADFIT3  Exact quadratic y = a*x^2 + b*x + c through three points.
%   [a,b,c] = quadFit3(x,y), with x,y each 1x3 (or 3x1), returns the
%   unique quadratic passing exactly through the three points, via Newton
%   divided differences. Valid for unequally spaced x. Three collinear
%   (or otherwise exactly flat) points give a = 0 exactly -- callers that
%   need the vertex -b/(2a) must check for a == 0 (or, for near-flat real
%   data, a small-a case) themselves rather than dividing by it; this
%   function does not guess a fallback.
%
% The vertex of the returned parabola, x* = -b/(2*a) (only meaningful
% when a ~= 0), is the sub-grid interpolated extremum location used by
% modalFeatureVector.m.

x = x(:).'; y = y(:).';
if numel(x) ~= 3 || numel(y) ~= 3
    error('quadFit3:badInput', 'x and y must each have exactly 3 elements.');
end
if numel(unique(x)) < 3
    error('quadFit3:repeatedX', 'x must have three distinct values.');
end

f12  = (y(2) - y(1)) / (x(2) - x(1));
f23  = (y(3) - y(2)) / (x(3) - x(2));
f123 = (f23 - f12) / (x(3) - x(1));

a = f123;
b = f12 - f123 * (x(1) + x(2));
c = y(1) - f12 * x(1) + f123 * x(1) * x(2);
end