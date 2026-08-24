%% Dg check
addpath('..');
clear all, close all, clc

alpha = 10; beta = 1e-2; gamma = 0.1; M = 6; P = 4;
kappa = dispersionRoots(alpha, beta, gamma, M);
kappaTower = kappa(3:end);
k = dispersionRoots(alpha, 0, 0, M);

Dg = gegenbauerCrossMatrix(kappaTower, gamma, P, 1-gamma);

fprintf('Quadrature cross-check for Dg (disc-covered):\n');
for pair = [1 0; 3 1; 5 2]'
    i = pair(1); j = pair(2);
    ki = kappaTower(i);
    integrand = @(z) verticalEigenfunction(z, ki, gamma) .* ...
                     gegenbauerBasisFunction(2*j, gamma, (z+1)/(1-gamma));
    numeric = integral(integrand, -1, -gamma, 'ArrayValued', true);
    closedForm = Dg(i, j+1);
    fprintf('  Dg(%d,%d): closed-form=%.6f%+.6fi, quadrature=%.6f%+.6fi, diff=%.3e\n', ...
        i, j, real(closedForm), imag(closedForm), real(numeric), imag(numeric), abs(closedForm-numeric));
end

% D0g check:
D0g = gegenbauerCrossMatrix(k, gamma, P, 1);

fprintf('\nQuadrature cross-check for D0g (open water):\n');
for pair = [1 0; 3 1; 5 2]'
    i = pair(1); j = pair(2);
    ki = k(i);
    integrand = @(z) verticalEigenfunction(z, ki, 0) .* ...
                     gegenbauerBasisFunction(2*j, gamma, (z+1)/(1-gamma));
    numeric = integral(integrand, -1, -gamma, 'ArrayValued', true);
    closedForm = D0g(i, j+1);
    fprintf('  D0g(%d,%d): closed-form=%.6f%+.6fi, quadrature=%.6f%+.6fi, diff=%.3e\n', ...
        i, j, real(closedForm), imag(closedForm), real(numeric), imag(numeric), abs(closedForm-numeric));
end