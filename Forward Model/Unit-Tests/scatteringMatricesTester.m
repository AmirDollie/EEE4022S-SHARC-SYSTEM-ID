%% This is a HUGE TEST! First closed loop test
% Ensures that so far, every file together solves the problem I have 
% instructed it to

% BUT doesn't say that the problem was set up correctly... validated l8r
addpath('..');
clear all, close all, clc

alpha = 10; beta = 1e-2; gamma = 0.1; M = 6; P = 4;
R = 1; nu = 0.3; n = 3;

kappa = dispersionRoots(alpha, beta, gamma, M);
k = dispersionRoots(alpha, 0, 0, M);
V = bennettsCoefficients(alpha, beta, gamma, kappa);
Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa);

[EIprime, E0prime] = openWaterRadialDerivatives(n, R, k);
Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);

D0 = innerProductMatrix(k, 0);
D = innerProductMatrix(kappa(3:end), gamma);
D0g = gegenbauerCrossMatrix(k, gamma, P, 1);
Dg = gegenbauerCrossMatrix(kappa(3:end), gamma, P, 1-gamma);

[MIn, M0, Mmat] = matchingMatrices(EIprime, E0prime, Enprime, D0, D, D0g, Dg);
En_R = radialMatrix(n, R, R, kappa, V, Ln);
[MU, Sn, Sn0] = scatteringMatrices(D0g, Dg, En_R, MIn, M0, Mmat);

%ensure matrix dimensions are consistent
fprintf('MU: %dx%d, Sn: %dx%d, Sn0: %dx%d, all finite: %d\n', ...
    size(MU,1), size(MU,2), size(Sn,1), size(Sn,2), size(Sn0,1), size(Sn0,2), ...
    all(isfinite(MU(:))) && all(isfinite(Sn(:))) && all(isfinite(Sn0(:))));

%% THE highest-value check: recover (An, Aen), plug back into (A.3a)-(A.3b)
%
k0 = k(1);
AIn = incidentAmplitude(alpha, n, k0, R, M);

An = Sn*AIn;      % disc-covered tower amplitudes
Aen = Ln*An;      % damped amplitudes [A_{-2}; A_{-1}]

AfullVec = [Aen(1); Aen(2); An];   % matches kappa ordering exactly

residA = 0; residB = 0;
for idx = 1:length(kappa)
    kv = kappa(idx);
    dphi = edgeDerivative(kv, gamma);
    [row1, row2] = edgeConditionRows(kv, R, n, nu);
    residA = residA + dphi*row1*AfullVec(idx);
    residB = residB + dphi*row2*AfullVec(idx);
end

%Need the residual to be low.
%This would confirm that the calculated response still satisfies the edge
%conditions! No where is that enforced, thus it serves as a pretty good
%check.

fprintf('\nResidual of (A.3a): %.3e\n', abs(residA));
fprintf('Residual of (A.3b): %.3e\n', abs(residB));