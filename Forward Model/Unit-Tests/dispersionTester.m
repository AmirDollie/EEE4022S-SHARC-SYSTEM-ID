clear all, close all, clc
% Open water case:
%Note these parameters were chosen to match the ones present in Figure 2.3
%in chapter 2 of Montiel's thesis, to check accuracy.
alpha = 10;
M = 6;
k = dispersionRoots(alpha, 0, 0, M)
% k is 7x1: k(1) = k0 (imaginary), k(2:7) = k1..k6 (real, ascending)
% Checked outputs, can confirm they're correct! Read off graph.

%The above is validated directly from Montiels paper, exact same values

%% Disc Covered (gamma, beta > 0)

%Set of test cases 
cases = {
    10, 1e-4, 1e-3;
    10, 1e-2, 0.1;
    10, 1.0,  0.2;
    5,  1e-3, 0.05;
    20, 1e-3, 0.15;
    2,  0.5,  0.3;
};

for i = 1:size(cases,1)
    alpha = cases{i,1}; beta = cases{i,2}; gamma = cases{i,3};
    kappa = dispersionRoots(alpha, beta, gamma, 5);
    fprintf('case %d: kappa0=%.4fi  damped=%.4f%+.4fi  first evanescent=%.4f\n', ...
        i, imag(kappa(3)), real(kappa(1)), imag(kappa(1)), kappa(4));
end

%He never supplied gamma and beta (only alpha and M)
%So although i can plug these into the dispersion equ and show they are
%roots
%the question remains... are they the CORRECT roots.

%Perhaps contour integration would reveal.
