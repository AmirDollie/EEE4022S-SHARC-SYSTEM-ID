function val = gegenbauerBasisFunction(p, gammaDraught, y)
%GEGENBAUERBASISFUNCTION Evaluates the weighted Gegenbauer basis function
%C_p^(delta)(y) (eq. 8), delta=1/6 fixed.
%
%   VAL = GEGENBAUERBASISFUNCTION(P, GAMMADRAUGHT, Y)
%
%   Validation-only utility: used to numerically cross-check
%   gegenbauerCrossMatrix.m via direct quadrature, not called anywhere
%   in the main pipeline itself. Requires the Symbolic Math Toolbox
%   (gegenbauerC).
%
%   Uses builtin('gamma', ...) throughout for the same reason documented
%   in gegenbauerCrossMatrix.m: "gamma" is used elsewhere in this
%   codebase as the draught parameter, so a bare gamma(x) call inside
%   any function taking a gamma-named argument would not call the
%   built-in Gamma function.


    delta = 1/6;
    GammaDelta = builtin('gamma', delta);

    coeff = (2/(1-gammaDraught)) * (factorial(p)*(p+delta)*GammaDelta^2) ...
            / (pi * 2^(1-2*delta) * builtin('gamma', p+2*delta));

    val = coeff .* (1-y.^2).^(delta-0.5) .* gegenbauerC(p, delta, y);
end