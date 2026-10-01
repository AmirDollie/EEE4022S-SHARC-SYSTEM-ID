function S = welchSmoothedSpectrum(omegaBins, info, L, rows)
%WELCHSMOOTHEDSPECTRUM  Expected Welch (periodic Hann, length L) spectral matrix of the signal
%   part of a twin record at the given bins: the kernel-weighted average of G H H^H over the fine
%   grid, with the same kernel as welchExpectedTransmissibility. At coarse resolution this
%   matrix is no longer rank one (T varies across the window), which lowers coherence and adds
%   variance; the point-frequency H H^H Su misses that.
%   S = welchSmoothedSpectrum(omegaBins, info, L, rows)
%
%   info   synthesiseTwinRecords INFO (omega, H, inputPSD, taper, N, dt)
%   rows   sensor rows of info.H to include (layout order)
%   S      nR x nR x nBins, one-sided per rad/s (sum of the kernel over the band is ~1)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Engine/.
wq = info.omega;
G = info.inputPSD .* info.taper.^2;
Hq = info.H(rows, :);
U = 3 * L / 8;
K = abs(hannDTFT((omegaBins(:) - wq) * info.dt, L)).^2 / (info.N * U);   % nB x nQ
nR = numel(rows); nB = numel(omegaBins);
S = complex(zeros(nR, nR, nB));
for a = 1:nR
    for b = a:nR
        v = K * (G .* Hq(a, :) .* conj(Hq(b, :))).';
        S(a, b, :) = v; S(b, a, :) = conj(v);
    end
end
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