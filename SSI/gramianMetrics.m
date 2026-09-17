function metrics = gramianMetrics(Wo)
%GRAMIANMETRICS Computes the four scalar observability measures from
%Section 2.3 of the Level 2C literature review, given a Gramian matrix.
%
%   METRICS = GRAMIANMETRICS(WO)
%
%   WO - an n x n symmetric positive semi-definite Gramian matrix (e.g.
%        from finiteObservabilityGramian.m)
%
%   Returns METRICS, a struct with fields:
%   logDet    - log(det(Wo)), D-optimality, computed via sum(log(eigenvalues))
%               rather than log(det(Wo)) directly, to avoid numerical
%               underflow/overflow when det(Wo) itself is extremely small
%               or large -- particularly relevant to the potentially
%               weak-observability regime under this project's
%               four-sensor constraint (per Qi et al.'s finding that
%               determinant-based measures become unreliable there).
%   trace     - trace(Wo), A-optimality-related
%   lambdaMin - smallest eigenvalue, E-optimality
%   kappa     - condition number lambda_max/lambda_min. A SMALL kappa is
%               good (well-conditioned); this is the opposite direction
%               from the other three metrics, which are all maximised
%               for a good configuration. Callers doing a greedy search
%               on kappa should therefore MINIMISE it, or equivalently
%               maximise -kappa or 1/kappa.
%   eigenvalues - all eigenvalues of Wo, sorted ascending, returned for
%               inspection.
%
%   NUMERICAL RANK HANDLING: only eigenvalues that are numerically
%   negligible relative to the Gramian's own scale (|eigenvalue| <=
%   size(Wo)*eps(lambdaMax)) are treated as exactly zero, giving
%   logDet=-Inf and kappa=Inf -- correctly reflecting a direction that is
%   observable in theory but not resolvable in floating-point practice.
%   A SUBSTANTIALLY negative eigenvalue is different in kind: it means
%   Wo is not actually a valid Gramian (an upstream bug), and this is
%   surfaced as an error, not silently clipped to zero and hidden. An
%   earlier version of this function clipped ALL negative eigenvalues
%   unconditionally, which would have silently disguised a genuinely
%   broken (indefinite) input as an ordinary rank-deficient one -- caught
%   during review before this function was used anywhere else.

    n = size(Wo,1);

    Wo_sym = (Wo + Wo')/2;
    eigenvalues = sort(eig(Wo_sym), 'ascend');

    lambdaMax = max(abs(eigenvalues));
    eigTol = n * eps(max(lambdaMax, 1));

    if eigenvalues(1) < -eigTol
        error('gramianMetrics:notPSD', ...
            'Wo is not positive semi-definite: minimum eigenvalue = %.3e (tolerance = %.3e). This indicates an upstream error in how Wo was constructed, not floating-point roundoff.', ...
            eigenvalues(1), eigTol);
    end

    eigenvalues(abs(eigenvalues) <= eigTol) = 0;

    metrics.eigenvalues = eigenvalues;
    metrics.trace = sum(eigenvalues);
    metrics.lambdaMin = eigenvalues(1);

    if eigenvalues(1) <= 0
        metrics.logDet = -Inf;
    else
        metrics.logDet = sum(log(eigenvalues));
    end

    if eigenvalues(1) <= 0
        metrics.kappa = Inf;
    else
        metrics.kappa = eigenvalues(n) / eigenvalues(1);
    end
end