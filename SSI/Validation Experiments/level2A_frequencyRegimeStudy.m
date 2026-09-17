%% level2A_frequencyRegimeStudy.m
%
% LEVEL 2A: frequency-regime study. Fixed sensor locations, omega1 held
% fixed at 2.70 rad/s, omega2 varied across several pairs. For each pair:
% compute the TRUE (Forward-Model-only) cross-MAC between the two known
% mode shapes, then run the full SSI chain and compute SSI's OWN
% cross-MAC between its two recovered shapes, and compare the two.
%
% BUG FIX (found by external code review): a previous version of this
% script computed ssiCrossMAC as computeMAC(shape1, knownShapes(:,2)) --
% comparing SSI's recovered shape 1 against the TRUE shape 2, never
% against SSI's own recovered shape 2. Since matching MAC (a recovered
% shape vs its OWN true shape) was already near-perfect everywhere,
% ssiCrossMAC was strongly biased toward knownCrossMAC whenever shape 1
% was recovered accurately, regardless of how well SSI actually recovered
% shape 2. This version extracts BOTH shape1 and shape2 from SSI's own
% output and compares them to EACH OTHER, which is the correct,
% non-circular test.
%
% SECOND FIX (also found by review): the original common-order logic
% assumed both branches end at the same model order, which silently
% produced a false "NO COMMON ORDER" whenever they didn't. Now uses the
% actual intersection of both branches' orders and picks the latest one
% genuinely shared by both.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'JONSWAP'));
clear all, close all, clc

H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;
g = 9.81;

sensorLocations = [0, 0; 0.5*R, 0; 0.9*R, 0; 0.9*R, pi];
nSensors = size(sensorLocations,1);

dt = 0.02; tVec = 0:dt:80;
i = 40; nMax = 12;
freqTol = 0.01; dampTol = 0.05; macTol = 0.99;

omega1 = 2.70;
omega2List = [2.80, 4.10, 6.00, 8.00, 9.50];

fprintf('%10s %10s %16s %16s %10s\n', 'omega1', 'omega2', 'knownCrossMAC', 'ssiCrossMAC', 'agree?');

for p = 1:length(omega2List)
    omega2 = omega2List(p);
    omega_components = [omega1; omega2];
    a_components = [0.05; 0.05];
    epsilon_components = [0; 0];
    alpha_components = H * omega_components.^2 / g;

    binData = cell(2,1);
    for c = 1:2
        binData{c} = precomputeDeflectionData(alpha_components(c), beta, gamma, R, nu, M, P, N);
    end
    specData.omega = omega_components;
    specData.a = a_components;
    specData.epsilon = epsilon_components;
    specData.alpha = alpha_components;
    specData.H = H;
    specData.R = R;
    specData.binData = binData;

    knownShapes = zeros(nSensors, 2);
    for c = 1:2
        for s = 1:nSensors
            knownShapes(s,c) = evaluateDeflection(specData.binData{c}, sensorLocations(s,1), sensorLocations(s,2));
        end
    end
    knownCrossMAC = computeMAC(knownShapes(:,1), knownShapes(:,2));

    y = syntheticSensorData(specData, sensorLocations, tVec);
    [Yp_ref, Yf, ~] = buildHankelMatrix(y, i);
    P_proj = hankelProjection(Yp_ref, Yf);

    warning('off', 'modalParameters:realEigenvalue');
    results = sweepModelOrders(P_proj, nSensors, dt, nMax);
    warning('on', 'modalParameters:realEigenvalue');

    matches = cell(nMax,1);
    for n = 2:nMax
        matches{n} = classifyPoleStability(results(n-1), results(n), freqTol, dampTol, macTol);
    end
    branches = buildModalBranches(results, matches);
    persistence = computeBranchPersistence(branches, 3, 2);

    foundBranch = zeros(2,1);
    for c = 1:2
        for b = 1:length(branches)
            if any(abs(branches(b).frequency - omega_components(c))/omega_components(c) < freqTol) && persistence(b).isPersistent
                foundBranch(c) = b;
                break;
            end
        end
    end

    if any(foundBranch == 0)
        fprintf('%10.2f %10.2f %16.4f %16s %10s\n', omega1, omega2, knownCrossMAC, 'N/A', 'NO BRANCH');
        continue;
    end

    b1 = foundBranch(1);
    b2 = foundBranch(2);

    % Defensive check: the two target frequencies must correspond to
    % distinct modal branches.
    if b1 == b2
        fprintf('%10.2f %10.2f %16.4f %16s %10s\n', omega1, omega2, knownCrossMAC, 'N/A', 'SAME BRANCH');
        continue;
    end

    % Use the highest model order at which BOTH branches are present,
    % not just wherever branch 1 happens to end.
    commonOrders = intersect(branches(b1).order, branches(b2).order);

    if isempty(commonOrders)
        fprintf('%10.2f %10.2f %16.4f %16s %10s\n', omega1, omega2, knownCrossMAC, 'N/A', 'NO COMMON ORDER');
        continue;
    end

    n_last = max(commonOrders);

    idx1 = find(branches(b1).order == n_last, 1);
    idx2 = find(branches(b2).order == n_last, 1);

    [~, col1] = min(abs(results(n_last).frequency - branches(b1).frequency(idx1)));
    [~, col2] = min(abs(results(n_last).frequency - branches(b2).frequency(idx2)));

    shape1 = results(n_last).modeShapes(:, col1);
    shape2 = results(n_last).modeShapes(:, col2);

    ssiCrossMAC = computeMAC(shape1, shape2);

    agree = abs(ssiCrossMAC - knownCrossMAC) < 0.01;

    fprintf('%10.2f %10.2f %16.4f %16.4f %10s\n', omega1, omega2, knownCrossMAC, ssiCrossMAC, mat2str(agree));
end

fprintf('\nCompare knownCrossMAC (true physics) against ssiCrossMAC (SSI''s own recovered\n');
fprintf('shapes, compared to EACH OTHER, not to the true shapes) -- this is the corrected,\n');
fprintf('non-circular version of this check.\n');