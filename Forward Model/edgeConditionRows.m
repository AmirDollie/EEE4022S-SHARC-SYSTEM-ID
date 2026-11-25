function [row1, row2] = edgeConditionRows(kappaVal, R, n, nu)
%EDGECONDITIONROWS The two bracketed expressions from (A.3a)-(A.3b).
%Extracted from edgeMap.m (previously a private local function there)
%since the end-to-end residual validation needs to call this too.
%edgeMap.m should be updated to call this external version instead of
%its own private copy -> same reasoning as dispersionFunction.m and
%scaledBesselDerivative.m: one formula, one place.
% F1 (2026-11-25): I_{n-1}/I_n from the exponentially scaled besseli (the scaling cancels in the
% ratio), so it no longer returns Inf/Inf = NaN for |Re(kappa R)| > ~700. Mathematically identical.
    x = kappaVal * R;
    IhatP = besseli(n-1, x, 1)/besseli(n, x, 1) - n/x;
    row1 = x^2 - (1-nu)*(x*IhatP - n^2);
    row2 = x^3*IhatP - n^2*(1-nu)*(x*IhatP - 1);
end
