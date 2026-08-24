%% D and D_0 tester
%Tests in 2 ways:

% 1: numerical quadrature cross-check against the derived CF solns
% 2: orthogonality of D_0 ie should be exactly diagonal!
addpath('..');
clear all, close all, clc

alpha = 10; beta = 1e-2; gamma = 0.1; M = 6;
kappa = dispersionRoots(alpha, beta, gamma, M);
kappaTower = kappa(3:end);
k = dispersionRoots(alpha, 0, 0, M);   % open-water tower, for D0

D0 = innerProductMatrix(k, 0);
D  = innerProductMatrix(kappaTower, gamma);

% Check 1: D0 should be exactly diagonal (orthogonality of V0)
offDiagD0 = D0 - diag(diag(D0));
%max off diagonal element should be ~0!
fprintf('max|off-diagonal of D0| = %.3e (should be ~machine precision)\n', max(abs(offDiagD0(:))));

% Check 2: D should NOT be diagonal (V is not orthogonal)
offDiagD = D - diag(diag(D));
% Max off diagonal element should NOT be ~0!!!
fprintf('max|off-diagonal of D|  = %.3e (should be clearly nonzero)\n', max(abs(offDiagD(:))));

% Check 3: numerical quadrature cross-check against closed form, a few entries
% literally just checks my closed form equation vs Matlab's quadrature
% formulation
fprintf('\nNumerical quadrature cross-check:\n');
for pair = [1 1; 1 3; 2 4]'
    i = pair(1); j = pair(2);
    ki = kappaTower(i); kj = kappaTower(j);
    integrand = @(z) verticalEigenfunction(z, ki, gamma) .* verticalEigenfunction(z, kj, gamma);
    numeric = integral(integrand, -1, -gamma, 'ArrayValued', true);
    fprintf('  D(%d,%d): closed-form=%.6f%+.6fi, quadrature=%.6f%+.6fi, diff=%.3e\n', ...
        i, j, real(D(i,j)), imag(D(i,j)), real(numeric), imag(numeric), abs(D(i,j)-numeric));
end