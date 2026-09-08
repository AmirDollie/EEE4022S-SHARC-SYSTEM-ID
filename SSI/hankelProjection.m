function [P, Q1, R] = hankelProjection(Yp_ref, Yf)

% Orthogonal projection of future outputs onto the row space of past
% reference outputs. See Peeters and De Roeck equations 32, ..., 35

% The goal of this is simple:
%   "Of everything in YF, how much is explainable purely from YP_REF?
%   By projecting YF onto YP_REF's row space, it produces the best possible
%   pewdiction of YF - using only linear combination of elements of YP_REF!

%   Note that whatever cannot be predicted is chosen to depict NOISE, or in
%   SSI terms, the Kalman filter innovations...

% Uses outputs from buildHankelMatrix.m:

%   YP_REF - (r*i) x j "past" block from buildHankelMatrix.m
%   YF     - (l*i) x j "future" block from buildHankelMatrix.m

%   Returns:
%   P  - (l*i) x j projection matrix, P = Yf/Yp_ref (the row space of Yf
%        projected onto the row space of Yp_ref). This is what gets fed
%        into the SVD step next (extractStateSpace.m).
%   Q1 - j x (r*i), the orthonormal-column factor used to build P.
%        Returned mainly for testing/inspection; not needed downstream.
%   R  - the full (r*i+l*i) x (r*i+l*i) lower-triangular factor from the
%        QR step, also returned for inspection.

%   Computed via QR rather than the direct projection formula
%   (P = Yf*Yp_ref'*pinv(Yp_ref*Yp_ref')*Yp_ref) because that direct
%   formula squares up Yp_ref*Yp_ref', which is numerically worse and,
%   for real data, far more expensive than a single QR factorization of
%   H itself. Both routes are mathematically identical -> confirmed in
%   the accompanying tester by comparing this against the direct formula
%   on a small case, not assumed.

%   As usual, a lot of yap for a simple program:
    r_i = size(Yp_ref, 1);
    j = size(Yp_ref, 2);

    H = [Yp_ref; Yf];   % (r_i + l_i) x j

    % Economy QR of H' (j x (r_i+l_i)): Qfull has orthonormal columns,
    % Rfull is upper triangular. Since H = (Qfull*Rfull)' = Rfull'*Qfull',
    % transposing gives H = R*Q' with R lower triangular.
    [Qfull, Rfull] = qr(H', 0);
    R = Rfull';
    Q1 = Qfull(:, 1:r_i);

    R21 = R(r_i+1:end, 1:r_i);

    P = R21 * Q1';
end