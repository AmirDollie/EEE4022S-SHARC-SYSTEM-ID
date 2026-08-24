%% Helper function to compute d/dr (scaled Bessel)
function val = scaledBesselDerivative(besselType, n, K, r, R)

%SCALEDBESSELDERIVATIVE d/dr[Xhat_n(Kr)] evaluated at r, where
%Xhat_n(Kr) = X_n(Kr)/X_n(KR) is the scaled Bessel function (eq. 6), for

%Note there have been 2 types:
%X = I (besselType='I') or X = K (besselType='K').

% By the chain rule: d/dr[X_n(Kr)/X_n(KR)] = K * X_n'(Kr) / X_n(KR).
% NAMING CLASH (inherited from the physics literature, not introduced
% here):
%  K is the wavenumber (kappa_m or k_m); K_n is the Bessel
% function. 

    x  = K*r;
    xR = K*R;

    switch besselType
        case 'I'
            Xn_x  = besseli(n, x);
            Xn_xR = besseli(n, xR);
            logDeriv = besseli(n-1, x)/Xn_x - n/x;
        case 'K'
            Xn_x  = besselk(n, x);
            Xn_xR = besselk(n, xR);
            logDeriv = -besselk(n-1, x)/Xn_x - n/x;
        otherwise
            error('scaledBesselDerivative:badType', "besselType must be 'I' or 'K'.");
    end

    val = K * Xn_x * logDeriv / Xn_xR;


end