function [row1, row2] = edgeConditionRows(kappaVal, R, n, nu)
%EDGECONDITIONROWS The two bracketed expressions from (A.3a)-(A.3b).
%Extracted from edgeMap.m (previously a private local function there)
%since the end-to-end residual validation needs to call this too.
%edgeMap.m should be updated to call this external version instead of
%its own private copy -> same reasoning as dispersionFunction.m and
%scaledBesselDerivative.m: one formula, one place.
    x = kappaVal * R;
    IhatP = besseli(n-1, x)/besseli(n, x) - n/x;
    row1 = x^2 - (1-nu)*(x*IhatP - n^2);
    row2 = x^3*IhatP - n^2*(1-nu)*(x*IhatP - 1);
end