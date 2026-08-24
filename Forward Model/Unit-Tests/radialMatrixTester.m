%% radialMatrixTester.m

%Here I am trying to test the identity I derived:
% E_n(R) = I + V*L_n

% It does NOT validate the untested Bennetts coefficients though
% As at the boundary they could somehow still satisfy this condition.
addpath('..');
clear all, close all, clc

cases = {
    10, 1e-2, 0.1,  20;
    10, 1e-4, 1e-3, 20;
    5,  1e-3, 0.05, 30;
};

for c = 1:size(cases,1)
    alpha = cases{c,1}; beta = cases{c,2}; gamma = cases{c,3}; M = cases{c,4};
    R = 1; nu = 0.3;
    kappa = dispersionRoots(alpha, beta, gamma, M);
    V = bennettsCoefficients(alpha, beta, gamma, kappa);

    fprintf('\n=== Case %d: alpha=%.4g, beta=%.4g, gamma=%.4g, M=%d ===\n', c, alpha, beta, gamma, M);

    for n = [0, 1, 5]
        Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa);
        En_at_R = radialMatrix(n, R, R, kappa, V, Ln);

        expected = eye(M+1) + V*Ln;
        maxDiff = max(abs(En_at_R(:) - expected(:)));

        fprintf('  n=%3d: max|En(R) - (I+V*Ln)| = %.3e\n', n, maxDiff);
    end
end