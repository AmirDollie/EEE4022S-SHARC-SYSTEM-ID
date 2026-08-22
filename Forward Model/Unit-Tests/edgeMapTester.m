addpath('..');
clear all, close all, clc

% Test 1: basic sanity across a range of n and (alpha,beta,gamma,M)
% Should output 2 x (M+1) regardless of the parameters
% Moreover, everything should be FINITE
%Doesn't check physics, only ensuring it doesn't output GARBAGE
cases = {
    10, 1e-2, 0.1,  20;
    10, 1e-4, 1e-3, 20;
    5,  1e-3, 0.05, 30;
};

for c = 1:size(cases,1)
    alpha = cases{c,1}; beta = cases{c,2}; gamma = cases{c,3}; M = cases{c,4};
    kappa = dispersionRoots(alpha, beta, gamma, M);
    R = 1;
    nu = 0.3;

    for n = [0, 1, -1, 5, -5]
        Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa);
        sizeOK = isequal(size(Ln), [2, M+1]);
        anyBad = any(~isfinite(Ln(:)));
        fprintf('case %d, n=%3d: size=%dx%d (expect 2x%d), finite=%d\n', ...
            c, n, size(Ln,1), size(Ln,2), M+1, ~anyBad);
    end
end

% Test 2: confirm n=0 (besseli(-1,x)) doesn't error and gives a sane value
% testing the overflow
n = 0;
kappa = dispersionRoots(10, 1e-2, 0.1, 20);
Ln = edgeMap(n, 10, 1e-2, 0.1, 1, 0.3, kappa);
fprintf('\nn=0 case: any NaN/Inf = %d\n', any(~isfinite(Ln(:))));