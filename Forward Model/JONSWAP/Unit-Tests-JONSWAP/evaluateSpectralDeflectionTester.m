%% evaluateSpectralDeflectionTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..'));
addpath(fullfile(thisDir, '..', '..', 'Animation'));
clear all, close all, clc

H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;
rTest = 0.7*R; thetaTest = pi/4;

%% Test 1: near-single-spike spectrum should collapse onto deflection.m directly
% A very narrow, single-bin "spectrum" isolates one frequency almost
% exactly -> if the spectral machinery is wired correctly, this should
% match calling deflection.m directly at that one alpha, scaled by that
% bin's own amplitude and phase.
fprintf('=== Test 1: single-bin collapse check ===\n');
Hs = 0.1; w0 = 5.5; gammaJONSWAP = 3.3;
wMin = 5.49; wMax = 5.51; numBins = 1;   % one very narrow bin at essentially w0

specData1 = precomputeSpectralData(Hs, w0, gammaJONSWAP, wMin, wMax, numBins, ...
    H, beta, gamma, R, nu, M, P, N);

alpha1 = specData1.alpha(1);
a1 = specData1.a(1);
eps1 = specData1.epsilon(1);
omega1 = specData1.omega(1);

etaDirect = deflection(rTest, thetaTest, alpha1, beta, gamma, R, nu, M, P, N);

tVec = linspace(0, 10, 50);
zetaSpectral = evaluateSpectralDeflection(specData1, rTest, thetaTest, tVec);
zetaExpected = a1 * real(etaDirect * exp(1i*(omega1*tVec + eps1)));

fprintf('max|zetaSpectral - zetaExpected| = %.3e (should be ~0)\n', ...
    max(abs(zetaSpectral - zetaExpected)));

%% Test 2: incident-wave-only statistical check (no floe physics at all)
% Isolates whether the spectral synthesis itself (amplitudes, random
% phases, summation) is correct, decoupled from any floe response ->
% synthesize the raw sea surface directly from (omega,a,epsilon) and
% check its variance/Hs over a long synthetic time series.
%
% Repeated across several trials with FRESH random phases each time,
% since Hs_synthesized = 4*std(...) is a Monte-Carlo estimate from a
% finite, randomly-phased record -> a single run scattering a bit above
% or below the target is expected statistical noise, not necessarily a
% bug. What matters is whether it scatters AROUND the target across
% trials, or sits CONSISTENTLY below/above it.
fprintf('\n=== Test 2: raw incident wave statistics (5 trials, fresh phases each) ===\n');
Hs2 = 0.1; w0_2 = 5.5; gammaJONSWAP2 = 3.3;
wMin2 = 3.0; wMax2 = 8.5; numBins2 = 200;

tLong = linspace(0, 500, 20000);   % long synthetic record for stable statistics

HsResults = zeros(5,1);
for trial = 1:5
    [omega2, a2, eps2] = discretizeSpectrum(Hs2, w0_2, gammaJONSWAP2, wMin2, wMax2, numBins2);
    seaSurface = zeros(size(tLong));
    for i = 1:numBins2
        seaSurface = seaSurface + a2(i)*cos(omega2(i)*tLong + eps2(i));
    end
    HsResults(trial) = 4*std(seaSurface);
    fprintf('trial %d: Hs_synthesized = %.4f\n', trial, HsResults(trial));
end

fprintf('\nTarget Hs = %.4f\n', Hs2);
fprintf('Mean across trials = %.4f, std across trials = %.4f\n', ...
    mean(HsResults), std(HsResults));