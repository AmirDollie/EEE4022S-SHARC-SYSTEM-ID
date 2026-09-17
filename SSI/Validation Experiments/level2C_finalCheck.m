%% level2C_finalCheck.m
%
% Final targeted checks before freezing Level 2C, per review.
% N counts MEASUREMENT INSTANTS in the finite-horizon Gramian sum
% Wo^(N) = sum_{k=0}^{N-1} (A^k)'C'CA^k, not state-transition intervals --
% tVec=0:dt:80 has 4001 measurement instants (4000 transitions between
% them), so N_record=4001 correctly matches the actual SSI record.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'JONSWAP'));
clear all, close all, clc

H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N_trunc = 10;
g = 9.81;
dt = 0.02;

omega1 = 2.70; omega2 = 6.00;
alpha1 = H*omega1^2/g; alpha2 = H*omega2^2/g;

A = buildKnownOscillatorA(omega1, omega2, dt);
binData1 = precomputeDeflectionData(alpha1, beta, gamma, R, nu, M, P, N_trunc);
binData2 = precomputeDeflectionData(alpha2, beta, gamma, R, nu, M, P, N_trunc);

rGrid = [0.3, 0.5, 0.7, 0.9] * R;
thetaGrid = linspace(0, 2*pi, 9); thetaGrid(end) = [];
candidateLocations = [0, 0];
for rr = rGrid
    for th = thetaGrid
        candidateLocations = [candidateLocations; rr, th]; %#ok<AGROW>
    end
end
numCandidates = size(candidateLocations,1);

phi1_raw = zeros(numCandidates,1);
phi2_raw = zeros(numCandidates,1);
for c = 1:numCandidates
    phi1_raw(c) = evaluateDeflection(binData1, candidateLocations(c,1), candidateLocations(c,2));
    phi2_raw(c) = evaluateDeflection(binData2, candidateLocations(c,1), candidateLocations(c,2));
end
norm1 = norm(phi1_raw); norm2 = norm(phi2_raw);
phi1_tilde = phi1_raw/norm1; phi2_tilde = phi2_raw/norm2;

candidateRows = zeros(numCandidates, 4);
for c = 1:numCandidates
    candidateRows(c,:) = buildCandidateCrows(phi1_tilde(c), phi2_tilde(c));
end
knownShapesAll = [phi1_tilde, phi2_tilde];

kappaFn = @(Wo) -gramianMetrics(Wo).kappa;
lambdaMinFn = @(Wo) gramianMetrics(Wo).lambdaMin;

%% ===================================================================
%  CHECK 1: cross-MAC (and quick SSI check) for the TRUE kappa optimum
%  at N=582, with a defensive common-order guard restored
%  ===================================================================
fprintf('=== Check 1: true kappa optimum at N=582 -- cross-MAC and SSI check ===\n');
N_main = 582;

[trueKappaIdx, trueKappaScore] = exhaustiveSensorSelection(candidateRows, A, N_main, kappaFn, 4);
fprintf('True kappa-optimal configuration: %s (kappa=%.6f)\n', mat2str(trueKappaIdx), -trueKappaScore);

trueKappaCrossMAC = computeMAC(knownShapesAll(trueKappaIdx,1), knownShapesAll(trueKappaIdx,2));
fprintf('Cross-MAC of the TRUE kappa optimum: %.4f (previously reported, on greedy''s near-miss: 0.3806)\n', trueKappaCrossMAC);

sensorLocations = candidateLocations(trueKappaIdx, :);
a0 = 0.05;
a_components = a0 ./ [norm1; norm2];
specData.omega = [omega1; omega2];
specData.a = a_components;
specData.epsilon = [0; 0];
specData.alpha = [alpha1; alpha2];
specData.H = H;
specData.R = R;
specData.binData = {binData1; binData2};

tVec = 0:dt:80;
y = syntheticSensorData(specData, sensorLocations, tVec);
[Yp_ref, Yf, ~] = buildHankelMatrix(y, 40);
P_proj = hankelProjection(Yp_ref, Yf);
warning('off', 'modalParameters:realEigenvalue');
results = sweepModelOrders(P_proj, 4, dt, 12);
warning('on', 'modalParameters:realEigenvalue');
matches = cell(12,1);
for n = 2:12
    matches{n} = classifyPoleStability(results(n-1), results(n), 0.01, 0.05, 0.99);
end
branches = buildModalBranches(results, matches);
persistence = computeBranchPersistence(branches, 3, 2);

foundBranch = zeros(2,1);
for c = 1:2
    for b = 1:length(branches)
        if any(abs(branches(b).frequency - specData.omega(c))/specData.omega(c) < 0.01) && persistence(b).isPersistent
            foundBranch(c) = b;
            break;
        end
    end
end

if all(foundBranch > 0)
    b1 = foundBranch(1);
    b2 = foundBranch(2);

    if b1 == b2
        fprintf('WARNING: both target frequencies were assigned to the same branch.\n');
    else
        commonOrders = intersect(branches(b1).order, branches(b2).order);
        if isempty(commonOrders)
            fprintf('WARNING: persistent branches found, but with no common model order.\n');
        else
            n_last = max(commonOrders);
            idx1 = find(branches(b1).order == n_last, 1);
            idx2 = find(branches(b2).order == n_last, 1);
            [~, col1] = min(abs(results(n_last).frequency - branches(b1).frequency(idx1)));
            [~, col2] = min(abs(results(n_last).frequency - branches(b2).frequency(idx2)));
            shape1 = results(n_last).modeShapes(:, col1);
            shape2 = results(n_last).modeShapes(:, col2);
            ssiCrossMAC_trueKappa = computeMAC(shape1, shape2);
            fprintf('SSI cross-MAC (true kappa optimum): %.4f (agrees with FM cross-MAC: %s)\n', ...
                ssiCrossMAC_trueKappa, mat2str(abs(ssiCrossMAC_trueKappa-trueKappaCrossMAC)<0.01));
        end
    end
else
    fprintf('WARNING: SSI did not find persistent branches for both frequencies on this configuration.\n');
end

%% ===================================================================
%  CHECK 2: exhaustive lambdaMin/kappa at N=4001 (matches the actual
%  80s SSI validation record: 4001 measurement instants), with a proper
%  three-way comparison against both previously-tested nearby horizons
%  for BOTH objectives
%  ===================================================================
fprintf('\n=== Check 2: exhaustive lambdaMin/kappa at N=4001 (matches 80s SSI record) ===\n');
N_record = length(tVec);   % 4001 measurement instants
fprintf('Using N=%d (measurement instants in the 80s SSI validation record).\n', N_record);

[lambdaMinIdx_record, lambdaMinScore_record] = exhaustiveSensorSelection(candidateRows, A, N_record, lambdaMinFn, 4);
[kappaIdx_record, kappaScore_record] = exhaustiveSensorSelection(candidateRows, A, N_record, kappaFn, 4);

fprintf('lambdaMin optimum at N=%d: %s (score=%.6f)\n', N_record, mat2str(sort(lambdaMinIdx_record)), lambdaMinScore_record);
fprintf('kappa optimum at N=%d:     %s (kappa=%.6f)\n', N_record, mat2str(sort(kappaIdx_record)), -kappaScore_record);

fprintf('\n--- lambdaMin: three-way comparison ---\n');
lambda5 = [26 27 29 30];
lambda10 = [26 27 30 33];
lambdaRecord = sort(lambdaMinIdx_record);
if isequal(lambdaRecord, lambda10)
    fprintf('  -> matches the N=1164 (10-cycle) configuration.\n');
    fprintf('     This is CONSISTENT WITH the longer-horizon configuration persisting\n');
    fprintf('     to the actual SSI record length (endpoints checked, not continuity).\n');
elseif isequal(lambdaRecord, lambda5)
    fprintf('  -> matches the N=582 (5-cycle) configuration instead: the optimum has\n');
    fprintf('     reverted relative to the 10-cycle result, confirming placement\n');
    fprintf('     remains horizon-dependent even beyond 10 cycles.\n');
else
    fprintf('  -> differs from BOTH the 5-cycle and 10-cycle configurations.\n');
    fprintf('     lambdaMin placement is clearly horizon-specific over this range.\n');
end

fprintf('\n--- kappa: three-way comparison ---\n');
kappaLong = sort([2 6 14 26]);    % N=582, 1164
kappaShort = sort([2 14 22 26]);  % N=116, 233
kappaRecord = sort(kappaIdx_record);
if isequal(kappaRecord, kappaLong)
    fprintf('  -> matches the N=582/1164 (5- and 10-cycle) configuration %s.\n', mat2str(kappaLong));
    fprintf('     Consistent with that longer-horizon regime persisting to the actual record.\n');
elseif isequal(kappaRecord, kappaShort)
    fprintf('  -> reverts to the N=116/233 (1- and 2-cycle) configuration %s.\n', mat2str(kappaShort));
else
    fprintf('  -> differs from both previously observed regimes.\n');
end

recordCrossMAC_lambdaMin = computeMAC(knownShapesAll(lambdaMinIdx_record,1), knownShapesAll(lambdaMinIdx_record,2));
recordCrossMAC_kappa = computeMAC(knownShapesAll(kappaIdx_record,1), knownShapesAll(kappaIdx_record,2));
fprintf('\nCross-MAC at the actual record length: lambdaMin=%.4f, kappa=%.4f\n', recordCrossMAC_lambdaMin, recordCrossMAC_kappa);

%% ===================================================================
%  CHECK 3: are lambdaMin's tied optima at N=116, 233, 582 actually the
%  SAME four configurations -- using a strict numerical-tie tolerance
%  ===================================================================
fprintf('\n=== Check 3: are the 4 tied lambdaMin optima IDENTICAL across N=116,233,582? ===\n');
Ns_short = [116, 233, 582];
combos = nchoosek(1:numCandidates, 4);

optimalSets = cell(length(Ns_short),1);
for k = 1:length(Ns_short)
    [~, bestScore, allScores] = exhaustiveSensorSelection(candidateRows, A, Ns_short(k), lambdaMinFn, 4);
    tol = 1e-10 * max(1, abs(bestScore));   % strict: testing numerical ties, not near-optimality
    tiedRows = find(abs(allScores - bestScore) <= tol);
    optimalSets{k} = sort(combos(tiedRows,:), 2);
    fprintf('  N=%3d: %d tied optima:\n', Ns_short(k), length(tiedRows));
    for r = 1:length(tiedRows)
        fprintf('    %s\n', mat2str(optimalSets{k}(r,:)));
    end
end

setsMatch = isequal(sortrows(optimalSets{1}), sortrows(optimalSets{2})) && ...
            isequal(sortrows(optimalSets{2}), sortrows(optimalSets{3}));
if setsMatch
    fprintf('\n-> CONFIRMED: the exact same 4 configurations are tied-optimal at N=116, 233, AND 582.\n');
    fprintf('   "Stable from 1-5 cycles" is verified directly, not merely inferred from matching counts.\n');
else
    fprintf('\n-> The tied-optimal SETS differ across N=116,233,582 despite matching counts/representatives --\n');
    fprintf('   "stable from 1-5 cycles" is NOT actually supported; treat this range as unverified.\n');
end

fprintf('\n=== All checks complete. Ready to write up Level 2C using these verified figures. ===\n');