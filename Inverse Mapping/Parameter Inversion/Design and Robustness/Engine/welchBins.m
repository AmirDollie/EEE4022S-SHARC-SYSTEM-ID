function [omega, k] = welchBins(L, dt, bandEst, binStep, bandUse)
%WELCHBINS  Retained Welch bins, built exactly as the U1/L2048/G3 scripts build them:
%   bins of the L-point grid inside the estimator band, every binStep-th from the first,
%   then trimmed to the use band.
%   [omega, k] = welchBins(L, dt, bandEst, binStep, bandUse)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Engine/.
k = 1:L / 2 - 1;
w = 2 * pi * k / (L * dt);
keep = w >= bandEst(1) & w <= bandEst(2);
k = k(keep); w = w(keep);
k = k(1:binStep:end); w = w(1:binStep:end);
use = w >= bandUse(1) & w <= bandUse(2);
k = k(use); omega = w(use);
end