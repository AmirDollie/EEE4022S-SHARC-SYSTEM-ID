function data = precomputeDeflectionDataSplit(alpha, beta, gamma, R, nu, M, P, N, Mo)
%PRECOMPUTEDEFLECTIONDATASPLIT  As precomputeDeflectionData, with the open-water tower truncated
%separately at Mo modes (typically Mo >> M).   (F1, 2026-11-25)
%
%   DATA = PRECOMPUTEDEFLECTIONDATASPLIT(ALPHA, BETA, GAMMA, R, NU, M, P, N, MO)
%   returns the same struct as precomputeDeflectionData (alpha, R, N, kappa, dphiFull, AfullMat), for
%   evaluateDeflection, plus data.truncation = [M P N Mo].
%
%   WHY. diag_F1_floor.m and diag_F1_Mo_ladder.m showed the truncation error has two independent parts:
%   the plate tower (M) converges as M^-p with p ~ 1.5-1.75, and the open-water tower (Mo) converges as a
%   staircase with steps every ~1/gamma modes (envelope ~ Mo^-1.4). With Mo = M, the open-water error
%   (~1e-4 at M = 200-400 for p0) masks the plate convergence.
%
%   HOW. The two towers couple only through the P+1 Gegenbauer coefficients U. The open-water modes are
%   orthogonal (D0 diagonal), so MIn and M0 are elementwise and the open-water tower reduces to the
%   (P+1)x(P+1) sum G0 = sum_m D0g(m,:).' D0g(m,:) / (d0_m E0'_m). Per angular order n:
%       Mmat = (D Enprime) \ Dg,   Aeq = Dg.' En_R Mmat - G0,
%       U    = Aeq \ [D0g(1,:).' (1 + MIn_00) A^(In)_0],   A_n = Mmat U,   A~_n = L_n A_n.
%   This is the production algebra (matchingMatrices + scatteringMatrices) with the diagonal structure
%   used explicitly; it agrees with precomputeDeflectionData at Mo = M to ~1e-15
%   (diag_F1_Mo_ladder.m Part 0; Unit-Tests/precomputeDeflectionDataSplitTester.m). Cost is linear in Mo.
%   Open-water quantities use exponentially scaled Bessel functions; d0 uses the root identity
%   k tan k = -alpha, so nothing loses precision at k ~ 1e5.
%
%   The plate side (M) uses the production functions unchanged.
    addpath('..');

    kappa = dispersionRoots(alpha, beta, gamma, M);
    k = openWaterRoots(alpha, Mo);
    V = bennettsCoefficients(alpha, beta, gamma, kappa);
    k0 = k(1);

    % plate tower (production)
    D  = innerProductMatrix(kappa(3:end), gamma);
    Dg = gegenbauerCrossMatrix(kappa(3:end), gamma, P, 1-gamma);
    numModes = length(kappa);
    dphiFull = zeros(numModes, 1);
    for idx = 1:numModes
        dphiFull(idx) = edgeDerivative(kappa(idx), gamma);
    end

    % open-water tower (reduced)
    D0g = gegenbauerColumns(k, gamma, P);               % (Mo+1) x (P+1), as gegenbauerCrossMatrix(k, gamma, P, 1)
    t = -alpha ./ k;                                    % tan(k_m) from the dispersion relation
    d0 = 0.5 * (t ./ k + 1 + t.^2);                     % diag of innerProductMatrix(k, 0)
    x = k * R; x0 = x(1);

    AfullMat = zeros(numModes, 2*N+1);
    for nIdx = 1:(2*N+1)
        n = nIdx - N - 1;

        Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa);
        Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);
        Mmat = (D * Enprime) \ Dg;
        En_R = radialMatrix(n, R, R, kappa, V, Ln);

        E0p = k .* (-besselk(n-1, x, 1) ./ besselk(n, x, 1) - n ./ x);   % d/dr Khat_n(k r) at r = R
        G0 = D0g.' * (D0g ./ (d0 .* E0p));
        EIp0 = k0 * (besseli(n-1, x0, 1) / besseli(n, x0, 1) - n / x0); % d/dr Ihat_n(k0 r) at r = R
        MIn00 = -EIp0 / E0p(1);
        AIn0 = besseli(n, k0*R) / (1i*alpha);                            % incidentAmplitude, m = 0

        U = (Dg.' * En_R * Mmat - G0) \ (D0g(1, :).' * (1 + MIn00) * AIn0);
        An = Mmat * U;
        Aen = Ln * An;
        AfullMat(:, nIdx) = [Aen(1); Aen(2); An];
    end

    data.alpha = alpha;
    data.R = R;
    data.N = N;
    data.kappa = kappa;
    data.dphiFull = dphiFull;
    data.AfullMat = AfullMat;
    data.truncation = [M P N Mo];
end

function Dg = gegenbauerColumns(r, gam, P)
% Vectorised gegenbauerCrossMatrix(r, gam, P, 1): same closed form, delta = 1/6.
    delta = 1/6; G = builtin('gamma', delta);
    r = r(:); x = r * (1 - gam);
    Dg = zeros(numel(r), P + 1);
    for j = 0:P
        Dg(:, j+1) = G * (2 ./ x).^delta * (2*j + delta) * (-1)^j .* besselj(2*j + delta, x) ./ cos(r);
    end
end
