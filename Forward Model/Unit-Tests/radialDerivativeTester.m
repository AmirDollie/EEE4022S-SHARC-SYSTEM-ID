%% Compares the radial derivative closed forms to finite differences

% Checks scaledBesselDerivative, openWaterRadialDerivatives, and
% radialMatrixDerivative all against finite-difference derivatives of
% their "value" counterparts -> same style of check as edgeDerivative.m

addpath('..');
clear all, close all, clc

alpha = 10; beta = 1e-2; gamma = 0.1; M = 6;
R = 1; nu = 0.3;
n = 3;   %non-trivial angular order
h = 1e-6;

kappa = dispersionRoots(alpha, beta, gamma, M);
k = dispersionRoots(alpha, 0, 0, M);
V = bennettsCoefficients(alpha, beta, gamma, kappa);
Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa);

%% Test 1: scaledBesselDerivative vs finite difference, both I and K types
fprintf('=== Test 1: scaledBesselDerivative ===\n');
Kval = kappa(4);   % some real tower kappa, arbitrary choice

for besselType = ["I", "K"]
    if besselType == "I"
        scaledFn = @(r) besseli(n, Kval*r)/besseli(n, Kval*R);
    else
        scaledFn = @(r) besselk(n, Kval*r)/besselk(n, Kval*R);
    end

    fd = (scaledFn(R+h) - scaledFn(R-h)) / (2*h);
    closedForm = scaledBesselDerivative(char(besselType), n, Kval, R, R);

    fprintf('  type=%s: finite-diff=%.6f%+.6fi, closed-form=%.6f%+.6fi, diff=%.3e\n', ...
        besselType, real(fd), imag(fd), real(closedForm), imag(closedForm), abs(fd-closedForm));
end

%% Test 2: openWaterRadialDerivatives vs finite difference
fprintf('\n=== Test 2: openWaterRadialDerivatives ===\n');
[EIprime, E0prime] = openWaterRadialDerivatives(n, R, k);

for m = [1, 3, 5]
    km = k(m);

    EIscaled = @(r) besseli(n, km*r)/besseli(n, km*R);
    E0scaled = @(r) besselk(n, km*r)/besselk(n, km*R);

    fdI = (EIscaled(R+h) - EIscaled(R-h)) / (2*h);
    fdK = (E0scaled(R+h) - E0scaled(R-h)) / (2*h);

    fprintf('  m=%d: EIprime diag diff=%.3e, E0prime diag diff=%.3e\n', ...
        m, abs(fdI - EIprime(m,m)), abs(fdK - E0prime(m,m)));
end

%% Test 3: radialMatrixDerivative vs finite difference of radialMatrix
fprintf('\n=== Test 3: radialMatrixDerivative ===\n');
Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);

En_plus  = radialMatrix(n, R+h, R, kappa, V, Ln);
En_minus = radialMatrix(n, R-h, R, kappa, V, Ln);
Enprime_fd = (En_plus - En_minus) / (2*h);

maxDiff = max(abs(Enprime(:) - Enprime_fd(:)));
fprintf('  max|Enprime - finite-diff(radialMatrix)| = %.3e\n', maxDiff);