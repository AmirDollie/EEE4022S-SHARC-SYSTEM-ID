function [omega, zeta, Phi, lambda_c] = modalParameters(A, C, dt)
%MODALPARAMETERS Converts an identified discrete-time state matrix into
%physical modal parameters: natural frequency, damping ratio, and mode
%shapes (Peeters & De Roeck 1999, Section 7.1).
%
%  usgae: [OMEGA, ZETA, PHI, LAMBDA_C] = MODALPARAMETERS(A, C, DT)
%
%   A  - n x n identified discrete state matrix (from extractStateSpace.m)
%   C  - l x n identified output matrix (from extractStateSpace.m)
%   DT - sampling period (seconds)
%
%   Returns, for each of the n discrete eigenvalues of A:
%   OMEGA    - n x 1 natural (undamped) angular frequency, rad/s.
%              Real-valued (physically must be non-negative).
%   ZETA     - n x 1 damping ratio (dimensionless).
%   PHI      - l x n mode shapes at the sensor locations (complex, columns
%              matching the order of omega/zeta/lambda_c).
%   LAMBDA_C - n x 1 continuous-time poles themselves (complex), returned
%              in case the raw pole location is wanted directly rather
%              than the (omega,zeta) decomposition.
%
%   Method: 
% eigendecompose A = Psi*Lambda*Psi^-1. 
% Each discrete eigenvalue
%   lambda_q converts to a continuous-time pole via lambda_c = ln(lambda)/dt.
%   For a genuine oscillatory mode, these occur in complex-conjugate pairs
%   lambda_c, lambda_c* = -zeta*omega +/- i*omega*sqrt(1-zeta^2), giving
%   omega = |lambda_c| and zeta = -Re(lambda_c)/|lambda_c|. Mode shapes at
%   the sensor locations are Phi = C*Psi.
%
%   NOTE ON REAL EIGENVALUES: a real (non-complex) discrete eigenvalue
%   corresponds to an overdamped/non-oscillatory pole (zeta>=1, or a pole
%   on the real axis with no companion). This function still computes
%   omega=|lambda_c|, zeta=-real(lambda_c)/abs(lambda_c) for these -> the
%   formula is still algebraically valid, but such poles do not
%   represent an oscillatory structural mode in the usual sense, and are
%   flagged via a warning rather than silently returned as if they were
%   ordinary vibration modes.

    n = size(A,1);
    [Psi, Lambda] = eig(A);
    lambda = diag(Lambda);

    lambda_c = log(lambda)/dt;

    omega = abs(lambda_c);
    zeta = -real(lambda_c) ./ omega;

    realMask = abs(imag(lambda)) < 1e-10 * max(abs(lambda));
    if any(realMask)
        warning('modalParameters:realEigenvalue', ...
            '%d of %d identified poles are real (non-oscillatory) -- likely overdamped or spurious, not a vibration mode in the usual sense.', ...
            sum(realMask), n);
    end

    Phi = C * Psi;
end