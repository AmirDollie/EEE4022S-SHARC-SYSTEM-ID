function Wo = finiteObservabilityGramian(A, C, N)
%FINITEOBSERVABILITYGRAMIAN Computes the finite-horizon discrete-time
%observability Gramian, valid for !ANY! A (no stability requirement),
%unlike the infinite-horizon Lyapunov-equation form used by several
%papers in the Level 2C literature review.
%
%   WO = FINITEOBSERVABILITYGRAMIAN(A, C, N)
%
%   A - n x n state matrix
%   C - l x n output matrix (l sensor rows, e.g. from stacking
%       buildCandidateCrows.m outputs)
%   N - number of time steps to sum over
%
%   Returns WO, the n x n matrix
%       Wo = sum_{k=0}^{N-1} (A^k)' * C' * C * A^k
%
%   DELIBERATELY MINIMAL: this function only computes the sum for
%   whatever N it is given. It does NOT decide what N is appropriate!

%   That is the orchestration script's job (level2C_sensorPlacementStudy.m),
%   since for an undamped A (eigenvalues on the unit circle, as
%   buildKnownOscillatorA.m produces), Wo genuinely does NOT converge to
%   a fixed matrix as N grows: each additional cycle keeps adding
%   observable energy, with nothing to bound the sum. What should be
%   checked instead is whether CANDIDATE RANKINGS under a given metric
%   stabilise as N grows, not whether Wo itself does, see
%   finiteObservabilityGramianTester.m for a direct demonstration of why.

    n = size(A,1);
    Wo = zeros(n,n);
    Ak = eye(n);
    for k = 0:(N-1)
        Wo = Wo + Ak' * (C'*C) * Ak;
        Ak = A * Ak;
    end
end