%% deflectionTester.m
addpath('..');
clear all, close all, clc

alpha = 10; beta = 1e-2; gamma = 0.1; M = 10; P = 6; N = 10;
R = 1; nu = 0.3;

% Check 1: eta(0, theta) must be theta-independent
etaVals = zeros(4,1);
thetas = [0, pi/4, pi/2, pi];
for i = 1:4
    etaVals(i) = deflection(0, thetas(i), alpha, beta, gamma, R, nu, M, P, N);
end
fprintf('eta(0,theta) at theta = 0, pi/4, pi/2, pi:\n');
disp(etaVals);
fprintf('max spread = %.3e (should be ~0)\n', max(abs(etaVals - etaVals(1))));

% Check 2: convergence in N at a fixed interior point
r_test = 0.5; theta_test = pi/3;
fprintf('\nConvergence in N at r=%.2f, theta=%.2f:\n', r_test, theta_test);
for Ntest = [2, 5, 10, 15]
    etaN = deflection(r_test, theta_test, alpha, beta, gamma, R, nu, M, P, Ntest);
    fprintf('  N=%2d: eta = %.6f%+.6fi\n', Ntest, real(etaN), imag(etaN));
end