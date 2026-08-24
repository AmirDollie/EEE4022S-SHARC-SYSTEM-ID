function data = precomputeDeflectionData(alpha, beta, gamma, R, nu, M, P, N)
%PRECOMPUTEDEFLECTIONDATA Runs the full pipeline (Steps 0-6) once, for
%every angular order n=-N..N, and packages the result for fast repeated
%evaluation at many (r,theta) points via evaluateDeflection.m.
    addpath('..');

    kappa = dispersionRoots(alpha, beta, gamma, M);
    k = dispersionRoots(alpha, 0, 0, M);
    V = bennettsCoefficients(alpha, beta, gamma, kappa);
    k0 = k(1);

    D0 = innerProductMatrix(k, 0);
    D = innerProductMatrix(kappa(3:end), gamma);
    D0g = gegenbauerCrossMatrix(k, gamma, P, 1);
    Dg = gegenbauerCrossMatrix(kappa(3:end), gamma, P, 1-gamma);

    numModes = length(kappa);
    dphiFull = zeros(numModes, 1);
    for idx = 1:numModes
        dphiFull(idx) = edgeDerivative(kappa(idx), gamma);
    end

    AfullMat = zeros(numModes, 2*N+1);
    for nIdx = 1:(2*N+1)
        n = nIdx - N - 1;

        Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa);
        [EIprime, E0prime] = openWaterRadialDerivatives(n, R, k);
        Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);
        [MIn, M0, Mmat] = matchingMatrices(EIprime, E0prime, Enprime, D0, D, D0g, Dg);
        En_R = radialMatrix(n, R, R, kappa, V, Ln);
        [~, Sn, ~] = scatteringMatrices(D0g, Dg, En_R, MIn, M0, Mmat);

        AIn = incidentAmplitude(alpha, n, k0, R, M);
        An = Sn*AIn;
        Aen = Ln*An;
        AfullMat(:, nIdx) = [Aen(1); Aen(2); An];
    end

    data.alpha = alpha;
    data.R = R;
    data.N = N;
    data.kappa = kappa;
    data.dphiFull = dphiFull;
    data.AfullMat = AfullMat;
end