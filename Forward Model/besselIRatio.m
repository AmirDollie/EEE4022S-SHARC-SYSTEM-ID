function val = besselIRatio(n, z, zR)
%BESSELIRATIO  I_n(z) / I_n(zR) without overflow.   (F1, 2026-11-25)
%   Uses MATLAB's exponentially scaled Bessel function, besseli(n, x, 1) = I_n(x) .* exp(-abs(real(x))),
%   so the ratio is formed from O(1) numbers and the exponential factor is applied once:
%       I_n(z) / I_n(zR) = [Is_n(z) / Is_n(zR)] .* exp(abs(real(z)) - abs(real(zR))).
%   For r <= R the factor is exp(Re(kappa) (r - R)) <= 1, which underflows gracefully to 0.
%   The unscaled form overflows to Inf/Inf = NaN once |Re(kappa R)| passes about 700.
%   Elementwise in z and zR (implicit expansion).
    val = besseli(n, z, 1) ./ besseli(n, zR, 1) .* exp(abs(real(z)) - abs(real(zR)));
end
