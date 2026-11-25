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

% F1 (2026-11-25): every Bessel value is now MATLAB's exponentially scaled form, so nothing
% overflows for |Re(K R)| > ~700.
%   besseli(n,x,1) = I_n(x) exp(-|Re x|)   ->  I_n(x)/I_n(xR) = Is(x)/Is(xR) exp(|Re x| - |Re xR|)
%   besselk(n,x,1) = K_n(x) exp(x)         ->  K_n(x)/K_n(xR) = Ks(x)/Ks(xR) exp(xR - x)
% The log-derivatives X_{n-1}/X_n use the same scaling in numerator and denominator, which cancels.
% Mathematically identical to the previous version.

    x  = K*r;
    xR = K*R;

    switch besselType
        case 'I'
            ratio = besseli(n, x, 1) / besseli(n, xR, 1) * exp(abs(real(x)) - abs(real(xR)));
            logDeriv = besseli(n-1, x, 1)/besseli(n, x, 1) - n/x;
        case 'K'
            ratio = besselk(n, x, 1) / besselk(n, xR, 1) * exp(xR - x);
            logDeriv = -besselk(n-1, x, 1)/besselk(n, x, 1) - n/x;
        otherwise
            error('scaledBesselDerivative:badType', "besselType must be 'I' or 'K'.");
    end

    val = K * ratio * logDeriv;


end
