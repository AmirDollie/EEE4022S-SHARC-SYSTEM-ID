function specData = precomputeSpectralData(Hs, w0, gammaJONSWAP, wMin, wMax, numBins, H, beta, gamma, R, nu, M, P, N)
% discretizes JONSWAP and runs unmodified forward model once per freak bin
% Hs, w0, GAMMAJONSWAP, WMIN, WMAX, NUMBINS: passed straight to
% discretizeSpectrum.m (JONSWAP shape and discretization).
%   H = water depth (m): physical, used only to convert each bin's
%   angular frequency into that bin's non-dim alpha.
%   BETA, GAMMA, R, NU, M, P, N - the SAME floe parameters and truncation
%   for every bin (fixed, per the earlier truncation check).
%Since model is an LTI system, the superposition of each frequency bin will
% output the actual deflection
% Returns a struct SPECTDATA holding:
% omega, a, epsilon, alpha  - Nx1 vectors, one entry per bin
% H, R                      - carried through for evaluateSpectralDeflection
% binData                  - Nx1 cell array; binData{i} is the
%                                precomputeDeflectionData output for bin i
%Requires all forward model functions AND the precompute in Animation!
%
% addpath uses mfilename('fullpath') to anchor to THIS file's actual
% location on disk, not to whatever folder you happened to run from --
% plain addpath('..') breaks if this file is ever called from a
% different working directory (e.g. from a nested Unit-Tests folder).
    thisDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(thisDir, '..'));
    addpath(fullfile(thisDir, '..', 'Animation'));
%call the discretize spectrum function to get frequ bins, mags and
%phases
    [omega, a, epsilon] = discretizeSpectrum(Hs, w0, gammaJONSWAP, wMin, wMax, numBins);
%convert w to non-dim (alpha)
    g = 9.81;
    alphaVec = H * omega.^2 / g;
%assign elements to struct
    specData.omega = omega;
    specData.a = a;
    specData.epsilon = epsilon;
    specData.alpha = alphaVec;
    specData.H = H;
    specData.R = R;
    binData = cell(numBins, 1);
    ticStart = tic;
%precompute parameters for each alpha (omega)
for i = 1:numBins
        binData{i} = precomputeDeflectionData(alphaVec(i), beta, gamma, R, nu, M, P, N);
if mod(i, 10) == 0 || i == numBins
            fprintf('  bin %3d/%d done (alpha=%.4f), elapsed=%.1fs\n', ...
                i, numBins, alphaVec(i), toc(ticStart));
end
end
    specData.binData = binData;
    fprintf('precomputeSpectralData: %d bins, total %.1fs (%.2fs/bin avg)\n', ...
        numBins, toc(ticStart), toc(ticStart)/numBins);
end