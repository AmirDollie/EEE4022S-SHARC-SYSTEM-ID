function [omega, zeta, Phi, lambda_c, isOscillatory] = modalParameters(A, C, dt)
    n = size(A,1);
    [Psi, Lambda] = eig(A);
    lambda = diag(Lambda);

    lambda_c = log(lambda) / dt;

    omega = abs(lambda_c);
    zeta = -real(lambda_c) ./ omega;

    % A discrete eigenvalue is non-oscillatory if it's real -- INCLUDING
    % negative real, which log() maps to a spurious +i*pi/dt via the
    % complex branch cut, not a genuine oscillation. Checking imag(lambda)
    % (the discrete eigenvalue) rather than imag(lambda_c) (post-log)
    % avoids being fooled by that artifact.
    isOscillatory = abs(imag(lambda)) >= 1e-10 * max(abs(lambda));

    if any(~isOscillatory)
        warning('modalParameters:realEigenvalue', ...
            '%d of %d identified poles are real (non-oscillatory) -- likely overdamped or spurious, not a vibration mode in the usual sense.', ...
            sum(~isOscillatory), n);
    end

    Phi = C * Psi;
end