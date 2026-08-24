%% deflectionRefactorTester.m
addpath('..');
addpath('../Animation');

clear all, close all, clc

alpha = 10; beta = 1e-2; gamma = 0.1; M = 10; P = 6; N = 10; R = 1; nu = 0.3;

data = precomputeDeflectionData(alpha, beta, gamma, R, nu, M, P, N);

testPoints = [0.3, 0.5; 0.7, 1.2; 0.9, 3.0];
for i = 1:size(testPoints,1)
    r = testPoints(i,1); theta = testPoints(i,2);
    etaSlow = deflection(r, theta, alpha, beta, gamma, R, nu, M, P, N);
    etaFast = evaluateDeflection(data, r, theta);
    fprintf('r=%.2f, theta=%.2f: slow=%.6f%+.6fi, fast=%.6f%+.6fi, diff=%.3e\n', ...
        r, theta, real(etaSlow), imag(etaSlow), real(etaFast), imag(etaFast), abs(etaSlow-etaFast));
end