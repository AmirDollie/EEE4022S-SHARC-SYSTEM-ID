%% Simply checks dimensions are correct
addpath('..');
clear all, close all, clc

alpha = 10; beta = 1e-2; gamma = 0.1; M = 6; P = 4; 
%Clearly from above we want (M+1) x (M+1) = 7x7
%and (M+1) x (P+1) = 7x5
R = 1; nu = 0.3; n = 3;

kappa = dispersionRoots(alpha, beta, gamma, M);
k = dispersionRoots(alpha, 0, 0, M);
V = bennettsCoefficients(alpha, beta, gamma, kappa);
Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa);

[EIprime, E0prime] = openWaterRadialDerivatives(n, R, k);
Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);

D0 = innerProductMatrix(k, 0);
D  = innerProductMatrix(kappa(3:end), gamma);
D0g = gegenbauerCrossMatrix(k, gamma, P, 1);
Dg  = gegenbauerCrossMatrix(kappa(3:end), gamma, P, 1-gamma);

[MIn, M0, Mmat] = matchingMatrices(EIprime, E0prime, Enprime, D0, D, D0g, Dg);


fprintf('MIn:  size=%dx%d (expect %dx%d), finite=%d\n', ...
    size(MIn,1), size(MIn,2), M+1, M+1, all(isfinite(MIn(:))));
fprintf('M0:   size=%dx%d (expect %dx%d), finite=%d\n', ...
    size(M0,1), size(M0,2), M+1, P+1, all(isfinite(M0(:))));
fprintf('Mmat: size=%dx%d (expect %dx%d), finite=%d\n', ...
    size(Mmat,1), size(Mmat,2), M+1, P+1, all(isfinite(Mmat(:))));