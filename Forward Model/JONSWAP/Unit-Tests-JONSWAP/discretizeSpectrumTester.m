%% Checks that the discrete set of comps carry the total var the cont spectrum says it should
addpath('..');
clear all, close all, clc

Hs = 10; w0 = 0.8; gammaJONSWAP = 3.3;
wMin = 0.1; wMax = 5;

fprintf('Checking variance match as numBins increases:\n');
for numBins = [20, 100, 500, 2000]
    [omega, a, epsilon] = discretizeSpectrum(Hs, w0, gammaJONSWAP, wMin, wMax, numBins);

    variance_discrete = sum(a.^2) / 2;
    Hs_recovered = 4 * sqrt(variance_discrete);

    fprintf('  numBins=%5d: variance=%.4f, Hs_recovered=%.4f (target=%.1f)\n', ...
        numBins, variance_discrete, Hs_recovered, Hs);
end

%Independent check: continuous integral of S(omega) over the SAME range,
%using a much finer grid, should match the total (Hs/4)^2 as wMax -> large
%enough to capture essentially all the spectral energy
wFine = linspace(wMin, wMax, 200000)';
Sfine = waveSpectrum(7, [Hs, w0, gammaJONSWAP], wFine, 0);
m0_continuous = trapz(wFine, Sfine);
fprintf('\nContinuous integral over [%.2f, %.2f]: Hs_check = %.4f\n', ...
    wMin, wMax, 4*sqrt(m0_continuous));

%Check phases are genuinely random and cover [0, 2*pi)
[~, ~, epsilonCheck] = discretizeSpectrum(Hs, w0, gammaJONSWAP, wMin, wMax, 1000);
fprintf('\nPhase range: min=%.4f, max=%.4f (should span close to [0, %.4f])\n', ...
    min(epsilonCheck), max(epsilonCheck), 2*pi);