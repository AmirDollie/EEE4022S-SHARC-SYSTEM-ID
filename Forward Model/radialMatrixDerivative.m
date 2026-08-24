function Enprime = radialMatrixDerivative(n, R, kappa, V, Ln)
%RADIALMATRIXDERIVATIVE d/dr[E_n(r)] at r=R, the derivative of the
%corrected radial matrix from radialMatrix.m.
%
% usage: ENPRIME = RADIALMATRIXDERIVATIVE(N, R, KAPPA, V, LN)
% Note it has the exact same form of radialMatrix.m, due to the linearity
% of the derivative:

%Since V and LN don't depend on r, d/dr[E_n(r)] = d/dr[Escript_n(r)]
% + V * d/dr[EscriptTilde_n(r)] * LN -- same structure as
%radialMatrix.m, just differentiating the diagonal entries.

    kappaTower = kappa(3:end);
    M1 = length(kappaTower);
    kappaDamped = [kappa(1), kappa(2)];

    EscriptPrimeDiag = zeros(M1,1);
    for m = 1:M1
        EscriptPrimeDiag(m) = scaledBesselDerivative('I', n, kappaTower(m), R, R);
    end
    EscriptPrime = diag(EscriptPrimeDiag);

    EscriptTildePrimeDiag = zeros(2,1);
    for j = 1:2
        EscriptTildePrimeDiag(j) = scaledBesselDerivative('I', n, kappaDamped(j), R, R);
    end
    EscriptTildePrime = diag(EscriptTildePrimeDiag);
    
    %apply linearity of d/dr
    Enprime = EscriptPrime + V*EscriptTildePrime*Ln;
end