%% Tester for synthetic sensor data to be fed into Hankel stuff
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'JONSWAP'));
clear all, close all, clc

%% Build a small precomputed spectral dataset to test against
Hs = 0.1; w0 = 5.5; gammaJONSWAP = 3.3;
wMin = 5.49; wMax = 5.51; numBins = 1;   % near-single-frequency, matches earlier twin test

H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;

specData = precomputeSpectralData(Hs, w0, gammaJONSWAP, wMin, wMax, numBins, ...
    H, beta, gamma, R, nu, M, P, N);

tVec = 0:0.05:20;

%% Test 1: single-sensor regression against calling evaluateSpectralDeflection directly
fprintf('=== Test 1: single-sensor regression ===\n');
sensorLocations1 = [0.7*R, pi/4];
yDirect = evaluateSpectralDeflection(specData, 0.7*R, pi/4, tVec);
ySynth = syntheticSensorData(specData, sensorLocations1, tVec);

assert(isequal(size(ySynth), [1, length(tVec)]), 'Test 1 failed: wrong output shape for 1 sensor.');
maxDiff = max(abs(ySynth - yDirect));
assert(maxDiff == 0, 'Test 1 failed: single-sensor output does not match direct call exactly.');
fprintf('PASS: shape correct, max|synth - direct| = %.3e\n', maxDiff);

%% Test 2: multiple sensors -- shape and per-row regression against direct calls
fprintf('\n=== Test 2: multi-sensor shape and per-row regression ===\n');
sensorLocations2 = [0, 0;
                    0.5*R, 0;
                    0.9*R, 0;
                    0.9*R, pi];

yMulti = syntheticSensorData(specData, sensorLocations2, tVec);
assert(isequal(size(yMulti), [4, length(tVec)]), 'Test 2 failed: wrong output shape for 4 sensors.');

allRowsMatch = true;
for i = 1:4
    r = sensorLocations2(i,1); theta = sensorLocations2(i,2);
    yRowDirect = evaluateSpectralDeflection(specData, r, theta, tVec);
    rowDiff = max(abs(yMulti(i,:) - yRowDirect));
    if rowDiff ~= 0
        allRowsMatch = false;
        fprintf('  row %d MISMATCH: diff = %.3e\n', i, rowDiff);
    end
end
assert(allRowsMatch, 'Test 2 failed: at least one row did not match its direct-call equivalent.');
fprintf('PASS: shape correct (4 sensors x %d samples), all 4 rows match direct calls exactly\n', length(tVec));

%% Test 3: rows are genuinely different from each other (sanity check that
% different sensor locations actually produce different signals, not the
% same value copied 4 times due to an indexing bug)
fprintf('\n=== Test 3: distinct sensor locations give distinct signals ===\n');
pairwiseDiffs = zeros(4,4);
for i = 1:4
    for j = 1:4
        pairwiseDiffs(i,j) = max(abs(yMulti(i,:) - yMulti(j,:)));
    end
end
offDiag = pairwiseDiffs(~eye(4));
assert(all(offDiag > 1e-10), 'Test 3 failed: at least two different sensor locations produced identical signals.');
fprintf('PASS: all 4 sensor signals are genuinely distinct (min pairwise diff = %.3e)\n', min(offDiag));

%% Test 4: output feeds cleanly into buildHankelMatrix.m
fprintf('\n=== Test 4: output is directly usable by buildHankelMatrix.m ===\n');
i_hankel = 30;
[Yp_ref, Yf, H_full] = buildHankelMatrix(yMulti, i_hankel);
fprintf('PASS: buildHankelMatrix ran without error on syntheticSensorData output. Yp_ref: %dx%d, Yf: %dx%d\n', ...
    size(Yp_ref,1), size(Yp_ref,2), size(Yf,1), size(Yf,2));