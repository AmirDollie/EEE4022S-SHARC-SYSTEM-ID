% Utilizes the edge conditions to reduce the unknown damped amplitudes into
% a matrix multiplied by the undamped amplitudes.
function Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa)
%For more detail, see Montiel's appendix A1
%Need to find Ln such that:
% Ae_n = Ln * A_n, Ae_n = (A_{n,-2}, A_{n,-1})^T
% which reduces our original (M+3) unknowns to (M+1)
%kappa needs to be the output of the dispersion roots for under the
%disc
%Ln is obviously a 2x(M+1)
%TRANSCRIPTION ERROR AGAIN: as mentioned before in edgeDerivative.m
% they evaluate at phi'_m(z = 0) not phi'_m(z = -gamma). It makes sense
% to evaluate at the actual edge though lol...
%ACTUAL TRANSCRIPTION NOTE WORTH BRINGING TO ATTENTION:
%Semiiiii ambiguous as to whether I'(kappa_m*R) is d/dx(I(x)) eval at
%x = kappa_m*R OR
% I'(kappa_m*R) os d/dr(I(kappa_m*r)) with kappa fixed and solve thru
% chain rule...
%I am almost 100% sure it is option 1.
%Note that the kappa vector is already ordered nicely such that no
%further reshuffling (unlike bennettsCoefficients.m)
%
% NOTE: edgeConditionRows (the A.3a/A.3b bracket expressions) used to
% live here as a private local function. It's now its own file,
% edgeConditionRows.m, since scatteringMatricesTester.m needs to call it
% too for the end-to-end residual check -- same "one formula, one place"
% reasoning as dispersionFunction.m and scaledBesselDerivative.m.
    kappaTower = kappa(3:end);
    M1 = length(kappaTower);
%build curly L (Lscript)
%each column is a tower (m>0) mode
    Lscript = zeros(2, M1);
for m = 1:M1
        km = kappaTower(m);
        dphi = edgeDerivative(km, gamma);
        [row1, row2] = edgeConditionRows(km, R, n, nu);
        Lscript(1,m) = dphi * row1;
        Lscript(2,m) = dphi * row2;
end
%build curly L tilde (LscriptTilde)
%cols = kappa_{-2} and kappa_{-1}
    LscriptTilde = zeros(2,2);
    kappaDamped = [kappa(1), kappa(2)];   % [kappa_{-2}, kappa_{-1}]
for col = 1:2
        kCol = kappaDamped(col);
        dphi = edgeDerivative(kCol, gamma);
        [row1, row2] = edgeConditionRows(kCol, R, n, nu);
        LscriptTilde(1,col) = dphi * row1;
        LscriptTilde(2,col) = dphi * row2;
end
%Invert: Ln = -(LscriptTilde)\(Lscript) equ A.4
    Ln = -(LscriptTilde) \ (Lscript);
end