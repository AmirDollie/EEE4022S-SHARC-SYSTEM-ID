%% Here are the D and D' dispersion functions
function [D, Dp] = dispersionFunction(xi, alpha, beta, gamma)

    z = xi .* (1 - gamma);
    A = 1 - alpha*gamma + beta .* xi.^4;

    %as per notes (Montiel Chapter 2 (just above equ 2.16)
    D = A .* xi .* sin(z) + alpha .* cos(z);

    if nargout > 1
        % analytical form of its derivative (Wolfram alpha verified)
        Dp = (1-gamma) .* A .* xi .* cos(z) + (1 - alpha + 5*beta.*xi.^4) .* sin(z);
    end

end