function [x, info] = r1Boundary(metric, threshold, bracket, tolX)
%R1BOUNDARY  Locate where a metric crosses a threshold, by bisection in x (R1: x = log10 S_nn).
%   [x, info] = r1Boundary(metric, threshold, [xLo xHi], tolX)
%
%   metric     handle x -> scalar value of the criterion (e.g. s_max, d_W or kappa at S_nn = 10^x)
%   threshold  the class threshold the metric crosses
%   bracket    [xLo xHi] with the crossing inside: sign(metric - threshold) differs at the two ends
%   tolX       stop when the bracket is narrower than tolX (decades for R1)
%
%   x is the midpoint of the final bracket. info: xLo, xHi (final bracket), fLo, fHi (metric at its ends),
%   nEval, and the end values of the input bracket. Bisection rather than a secant/fzero step: the metrics
%   can jump slightly where a bin is accepted or rejected by the estimator's denominator rule, and
%   bisection only needs the sign change.
%
%   Shared by runR1SNR and testR1SNR (that is the only reason it is a separate file).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

a = bracket(1); b = bracket(2);
fa = metric(a); fb = metric(b); n = 2;
sa = sign(fa - threshold); sb = sign(fb - threshold);
if sa == sb || sa == 0 || sb == 0
    error('r1Boundary:bracket', 'No sign change of metric - %.4g on [%.4g, %.4g] (%.4g, %.4g).', threshold, a, b, fa, fb);
end
info = struct('inputBracket', [a b], 'inputValues', [fa fb]);
while (b - a) > tolX
    m = 0.5 * (a + b);
    fm = metric(m); n = n + 1;
    sm = sign(fm - threshold);
    if sm == 0, a = m; b = m; fa = fm; fb = fm; break; end
    if sm == sa
        a = m; fa = fm;
    else
        b = m; fb = fm;
    end
end
x = 0.5 * (a + b);
info.xLo = a; info.xHi = b; info.fLo = fa; info.fHi = fb; info.nEval = n;
end