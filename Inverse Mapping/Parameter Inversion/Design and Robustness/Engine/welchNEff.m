function [nEff, K] = welchNEff(N, L, overlap)
%WELCHNEFF  Segments and effective number of averages of a periodic-Hann Welch estimate,
%   exactly as extractTransmissibility computes them (Welch 1967).
%   [nEff, K] = welchNEff(N, L, overlap)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Engine/.
w = 0.5 * (1 - cos(2 * pi * (0:L - 1).' / L));
U = sum(w.^2);
hop = max(1, round(L * (1 - overlap)));
K = numel(1:hop:(N - L + 1));
nEff = K;
if K > 1
    rho2 = 0;
    for l = 1:K - 1
        sh = l * hop;
        if sh >= L, break; end
        c = sum(w(1:L - sh) .* w(1 + sh:L)) / U;
        rho2 = rho2 + (1 - l / K) * c^2;
    end
    nEff = K / (1 + 2 * rho2);
end
end