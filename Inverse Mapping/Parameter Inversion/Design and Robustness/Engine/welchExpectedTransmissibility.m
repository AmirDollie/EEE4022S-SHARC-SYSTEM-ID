function T = welchExpectedTransmissibility(omegaBins, info, L, ref, others, SnnDen)
%WELCHEXPECTEDTRANSMISSIBILITY  Large-sample (no Monte Carlo noise) mean of the corrected Welch
%   transmissibility of a twin record: the spectral-window average of Saturday's estimator note,
%   identical to welchExpectedZ in runTransmissibilityCovariance.m.
%   T = welchExpectedTransmissibility(omegaBins, info, L, ref, others, SnnDen)
%
%   info    the INFO struct of synthesiseTwinRecords (fine grid omega, H, inputPSD, taper, N, dt)
%   SnnDen  noise PSD left in the denominator after correction (0 when the correction is exact)
%   T       numel(others) x nBins complex
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Engine/.
wq = info.omega;
G = info.inputPSD .* info.taper.^2;
Hq = info.H;
U = 3 * L / 8;
K = abs(hannDTFT((omegaBins(:) - wq) * info.dt, L)).^2 / (info.N * U);
Srr = K * (G .* abs(Hq(ref, :)).^2).' + SnnDen;
Sjr = K * (G .* Hq(others, :) .* conj(Hq(ref, :))).';
T = (Sjr ./ Srr).';
end

function W = hannDTFT(nu, L)
W = 0.5 * dirichletSum(nu, L) - 0.25 * dirichletSum(nu - 2 * pi / L, L) - 0.25 * dirichletSum(nu + 2 * pi / L, L);
end

function D = dirichletSum(x, L)
den = 1 - exp(-1i * x);
D = (1 - exp(-1i * x * L)) ./ den;
small = abs(den) < 1e-9;
D(small) = L;
end