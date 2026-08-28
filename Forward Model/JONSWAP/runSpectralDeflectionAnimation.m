%% 3D plot, along with optional saving
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', 'Animation'));
clear all, close all, clc
%% Settled, validated wave-basin-scale test configuration
Hs = 0.1; w0 = 5.5; gammaJONSWAP = 3.3;
wMin = 3.0; wMax = 8.5; numBins = 30;
H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;
fprintf('Precomputing spectral data (%d bins)...\n', numBins);
specData = precomputeSpectralData(Hs, w0, gammaJONSWAP, wMin, wMax, numBins, ...
    H, beta, gamma, R, nu, M, P, N);
%% Animate
nr = 50; ntheta = 50; nFrames = 100;   % bumped up from 60 -- smoother for a saved video on repeat playback
T0 = 2*pi/w0;        % peak period, for choosing a sensible animation length
tMax = 5*T0;         % show ~5 peak-wave cycles


%alter the arguments here to choose between saving or not
%see function description

%to do no saving: 
% plotSpectralDeflectionSurface(specData, nr, ntheta, nFrames, tMax);
plotSpectralDeflectionSurface(specData, nr, ntheta, nFrames, tMax, [], 'jonswapDeflection');