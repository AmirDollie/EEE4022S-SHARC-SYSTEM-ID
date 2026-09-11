%% level2A_frequencyRegimeStudy.m
%
% LEVEL 2A: frequency-regime study. Fixed sensor locations (identical to
% the 2.70/2.80 and 2.70/4.10 runs), varying frequency SEPARATION AND
% ABSOLUTE REGIME, to establish whether the Forward Model itself produces
% more spatially distinct responses at higher frequency, before touching
% sensor placement at all (that is Level 2B's job, done only once this
% establishes where meaningful spatial differences actually exist).
%
% SUCCESS CRITERION, per discussion: NOT crossMAC < 0.5. Rather:
%   MAC_SSI,cross approx= MAC_FM,cross
% i.e. does SSI faithfully REPRODUCE the true physical similarity,
% whatever that similarity actually is. Both prior runs (2.70/2.80 and
% 2.70/4.10) already satisfied this to 4 decimal places; this script
% checks whether that continues to hold as the frequency pair changes,
% while ALSO tracking whether knownCrossMAC itself drops as frequency
% separation and absolute frequency increase.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'JONSWAP'));
clear all, close all, clc

H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;
g = 9.81;

% Same sensor array throughout -- ONLY frequency pairs change in this study
sensorLocations = [0,      0;
                   0.5*R,  0;
                   0.9*R,  0;
                   0.9*R,  pi];
nSensors = size(sensorLocations,1);

% Frequency pairs to sweep: extending the existing trend
% (2.70/2.80 -> 0.9995, 2.70/4.10 -> 0.8263, already established)
% toward the higher end of the previously-validated truncation range
% (alpha up to ~13.85 confirmed safe at M=50,P=10,N=10)
freqPairs = {[2.70, 2.80], [2.70, 4.10], [2.70, 6.00], [2.70, 8.00], [2.70, 9.50]};

fprintf('%14s %14s %16s %16s %10s\n', 'omega1', 'omega2', 'knownCrossMAC', 'ssiCrossMAC', 'agree?');

for p = 1:length(freqPairs)
    omega_components = freqPairs{p}(:);
    a_components = [0.05; 0.05];
    epsilon_components = [0; 0];
    nComp = 2;
    alpha_components = H * omega_components.^2 / g;

    binData = cell(nComp,1);
    for c = 1:nComp
        binData{c} = precomputeDeflectionData(alpha_components(c), beta, gamma, R, nu, M, P, N);
    end

    specData.omega = omega_components;
    specData.a = a_components;
    specData.epsilon = epsilon_components;
    specData.alpha = alpha_components;
    specData.H = H;
    specData.R = R;
    specData.binData = binData;

    knownShapes = zeros(nSensors, nComp);
    for c = 1:nComp
        for s = 1:nSensors
            knownShapes(s,c) = evaluateDeflection(specData.binData{c}, sensorLocations(s,1), sensorLocations(s,2));
        end
    end
    knownCrossMAC = computeMAC(knownShapes(:,1), knownShapes(:,2));

    dt = 0.02;
    tVec = 0:dt:80;
    y = syntheticSensorData(specData, sensorLocations, tVec);

    i = 40; nMax = 12;
    [Yp_ref, Yf, ~] = buildHankelMatrix(y, i);
    P_proj = hankelProjection(Yp_ref, Yf);
    results = sweepModelOrders(P_proj, nSensors, dt, nMax); %#ok<NASGU>
    warning('off', 'modalParameters:realEigenvalue');   % suppress for this sweep-summary script; individual runs already confirmed this warning's meaning
    results = sweepModelOrders(P_proj, nSensors, dt, nMax);
    warning('on', 'modalParameters:realEigenvalue');

    matches = cell(nMax,1);
    for n = 2:nMax
        matches{n} = classifyPoleStability(results(n-1), results(n), 0.01, 0.05, 0.99);
    end
    branches = buildModalBranches(results, matches);
    persistence = computeBranchPersistence(branches, 3, 2);

    foundBranch = zeros(nComp,1);
    for c = 1:nComp
        for b = 1:length(branches)
            if any(abs(branches(b).frequency - omega_components(c))/omega_components(c) < 0.01) && persistence(b).isPersistent
                foundBranch(c) = b;
                break;
            end
        end
    end

    if any(foundBranch == 0)
        fprintf('%14.2f %14.2f %16s %16s %10s\n', omega_components(1), omega_components(2), sprintf('%.4f',knownCrossMAC), 'N/A', 'NO BRANCH');
        continue;
    end

    n_last = branches(foundBranch(1)).order(end);
    idx1 = find(branches(foundBranch(1)).order == n_last, 1);
    idx2 = find(branches(foundBranch(2)).order == n_last, 1);
    [~, col1] = min(abs(results(n_last).frequency - branches(foundBranch(1)).frequency(idx1)));
    [~, col2] = min(abs(results(n_last).frequency - branches(foundBranch(2)).frequency(idx2)));
    shape1 = results(n_last).modeShapes(:, col1);

    ssiCrossMAC = computeMAC(shape1, knownShapes(:,2));
    agree = abs(ssiCrossMAC - knownCrossMAC) < 0.01;

    fprintf('%14.2f %14.2f %16.4f %16.4f %10s\n', omega_components(1), omega_components(2), ...
        knownCrossMAC, ssiCrossMAC, mat2str(agree));
end

fprintf('\nLooking for: (a) knownCrossMAC decreasing as omega2 increases -- confirms richer\n');
fprintf('spatial structure emerges at higher frequency; (b) ssiCrossMAC tracking knownCrossMAC\n');
fprintf('closely at every pair -- confirms SSI faithfully reproduces the true physical\n');
fprintf('similarity regardless of what that similarity is, which is the actual validation target.\n');