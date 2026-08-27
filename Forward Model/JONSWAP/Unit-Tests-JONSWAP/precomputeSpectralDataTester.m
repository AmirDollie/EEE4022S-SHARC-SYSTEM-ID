%% precomputeSpectralDataTester:
% The actual main point of this is to gauge the time it will take to run a
% large scale (high bin count) simulation.

% Output appears to be ~1s/bin. Thus 100s/100bins etc
% addpath anchored to this file's actual location, not the working
% directory -- this file lives in JONSWAP/Unit-Tests-JONSWAP/, two
% levels below Forward Model/, so it needs two '..' to reach it.
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));                    % JONSWAP/
addpath(fullfile(thisDir, '..', '..'));              % Forward Model/
addpath(fullfile(thisDir, '..', '..', 'Animation'));  % Forward Model/Animation/
clear all, close all, clc
% Settled test configuration (small numBins for a quick first run)


Hs = 0.1; w0 = 5.5; gammaJONSWAP = 3.3;
wMin = 3.0; wMax = 8.5; numBins = 5;
H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;
specData = precomputeSpectralData(Hs, w0, gammaJONSWAP, wMin, wMax, numBins, ...
    H, beta, gamma, R, nu, M, P, N);
fprintf('\n--- Shape/finiteness checks ---\n');
fprintf('numBins requested = %d, binData entries = %d\n', numBins, length(specData.binData));
fprintf('alpha range: [%.4f, %.4f]\n', min(specData.alpha), max(specData.alpha));
fprintf('all alpha finite: %d\n', all(isfinite(specData.alpha)));
fprintf('all amplitudes finite and non-negative: %d\n', ...
    all(isfinite(specData.a)) && all(specData.a >= 0));
fprintf('\n--- Regression check: bin 3 vs direct call ---\n');
testBin = 3;
alphaTest = specData.alpha(testBin);
directData = precomputeDeflectionData(alphaTest, beta, gamma, R, nu, M, P, N);
rTest = 0.7*R; thetaTest = pi/4;
etaFromSpec = evaluateDeflection(specData.binData{testBin}, rTest, thetaTest);
etaDirect = evaluateDeflection(directData, rTest, thetaTest);
fprintf('From precomputeSpectralData: %.6f%+.6fi\n', real(etaFromSpec), imag(etaFromSpec));
fprintf('From direct call:            %.6f%+.6fi\n', real(etaDirect), imag(etaDirect));
fprintf('diff = %.3e (should be ~0, exact same computation)\n', abs(etaFromSpec - etaDirect));