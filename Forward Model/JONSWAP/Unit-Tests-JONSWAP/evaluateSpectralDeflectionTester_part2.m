%% evaluateSpectralDeflectionTester_part2.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..'));
addpath(fullfile(thisDir, '..', '..', 'Animation'));
clear all, close all, clc

H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;

%% Test 3: genuine multi-bin regression (not just a single bin)
fprintf('=== Test 3: multi-bin sum regression ===\n');
Hs = 0.1; w0 = 5.5; gammaJONSWAP = 3.3;
wMin = 4.5; wMax = 6.5; numBins = 3;

specData3 = precomputeSpectralData(Hs, w0, gammaJONSWAP, wMin, wMax, numBins, ...
    H, beta, gamma, R, nu, M, P, N);

rTest = 0.7*R; thetaTest = pi/4;
tVec = linspace(0, 8, 40);

% Build the expected sum by hand, bin by bin
zetaExpected = zeros(size(tVec));
for i = 1:numBins
    eta_i = evaluateDeflection(specData3.binData{i}, rTest, thetaTest);
    zetaExpected = zetaExpected + specData3.a(i) * ...
        real(eta_i * exp(1i*(specData3.omega(i)*tVec + specData3.epsilon(i))));
end

zetaActual = evaluateSpectralDeflection(specData3, rTest, thetaTest, tVec);
fprintf('max|actual - hand-summed expected| = %.3e (should be ~0)\n', ...
    max(abs(zetaActual - zetaExpected)));

%% Test 4: theta-independence at r=0, carried over to the time-domain spectral case
fprintf('\n=== Test 4: eta(0,theta,t) theta-independence, multi-bin ===\n');
thetas = [0, pi/3, pi, 3*pi/2];
zAtCentre = zeros(length(thetas), length(tVec));
for k = 1:length(thetas)
    zAtCentre(k,:) = evaluateSpectralDeflection(specData3, 0, thetas(k), tVec);
end
maxSpread = max(max(abs(zAtCentre - zAtCentre(1,:))));
fprintf('max spread across theta at r=0 (all t): %.3e (should be ~0)\n', maxSpread);

%% Test 5: confirm the r>R warning actually fires
fprintf('\n=== Test 5: r>R warning check ===\n');
lastwarn('');   % clear any prior warning
evaluateSpectralDeflection(specData3, 0.7, thetaTest, 0);   % r=0.7 > R=0.383
[warnMsg, warnId] = lastwarn;
fprintf('Warning fired: %d (id: %s)\n', ~isempty(warnMsg), warnId);