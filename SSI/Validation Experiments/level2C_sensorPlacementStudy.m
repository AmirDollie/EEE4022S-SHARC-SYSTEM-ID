%% level2C_sensorPlacementStudy.m
%
% LEVEL 2C: given a maximum of four sensors, which of the four candidate
% Gramian-based objectives (trace, logDet, lambdaMin, kappa) actually
% predicts good SSI-based frequency separation and spatial recovery?
%
% Frequencies fixed at 2.70/6.00 rad/s: widely separated and already
% robustly recovered in Levels 2A/2B, so this experiment isolates
% GEOMETRY, not frequency proximity (Level 3A) or amplitude imbalance
% (Level 3B).
%
% TERMINOLOGY: phi1, phi2 below are the TRUE SPATIAL RESPONSE FIELDS at
% omega1, omega2 -- forced steady-state responses, not established
% natural hydroelastic modes (see Level 1's interpretive note).

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

%% ===================================================================
%  PART 1: known oscillator A, and the two true spatial response fields
%  ===================================================================
fprintf('=== Part 1: known oscillator and true spatial response fields ===\n');
A = buildKnownOscillatorA(omega1, omega2, dt);

binData1 = precomputeDeflectionData(alpha1, beta, gamma, R, nu, M, P, N_trunc);
binData2 = precomputeDeflectionData(alpha2, beta, gamma, R, nu, M, P, N_trunc);
fprintf('A built (4x4, eigenvalues on unit circle by construction).\n');
fprintf('True spatial response fields ready at alpha1=%.4f, alpha2=%.4f.\n', alpha1, alpha2);

%% ===================================================================
%  PART 2: candidate spatial grid (INCLUDING the centre, so the Level 2B
%  baseline configuration can be reproduced within this grid), and the
%  two response fields evaluated on it
%  ===================================================================
fprintf('\n=== Part 2: candidate grid and true field evaluation ===\n');
rGrid = [0.3, 0.5, 0.7, 0.9] * R;
thetaGrid = linspace(0, 2*pi, 9); thetaGrid(end) = [];

candidateLocations = [0, 0];   % centre, added explicitly
for rr = rGrid
    for th = thetaGrid
        candidateLocations = [candidateLocations; rr, th]; %#ok<AGROW>
    end
end
numCandidates = size(candidateLocations,1);
fprintf('Candidate grid: %d locations (centre + %d radii x %d angles)\n', numCandidates, length(rGrid), length(thetaGrid));

phi1_raw = zeros(numCandidates,1);
phi2_raw = zeros(numCandidates,1);
for c = 1:numCandidates
    phi1_raw(c) = evaluateDeflection(binData1, candidateLocations(c,1), candidateLocations(c,2));
    phi2_raw(c) = evaluateDeflection(binData2, candidateLocations(c,1), candidateLocations(c,2));
end
fprintf('True spatial response fields evaluated at every candidate location.\n');

%% ===================================================================
%  PART 3: explicit normalisation (isolate geometric distinguishability
%  from amplitude imbalance, which is Level 3B's separate question).
%  norm1, norm2 are kept for Part 6's amplitude compensation.
%  ===================================================================
fprintf('\n=== Part 3: explicit spatial field normalisation ===\n');
norm1 = norm(phi1_raw);
norm2 = norm(phi2_raw);
phi1_tilde = phi1_raw / norm1;
phi2_tilde = phi2_raw / norm2;
fprintf('Both fields normalised to unit norm across the candidate grid.\n');
fprintf('(||phi1_raw||=%.4f, ||phi2_raw||=%.4f before normalisation)\n', norm1, norm2);

candidateRows = zeros(numCandidates, 4);
for c = 1:numCandidates
    candidateRows(c,:) = buildCandidateCrows(phi1_tilde(c), phi2_tilde(c));
end

%% ===================================================================
%  PART 4: N-sensitivity / ranking-stability check for ALL FOUR
%  objectives, not just trace -- lambdaMin and kappa have no
%  submodularity guarantee, so their stability across N cannot be
%  assumed from trace's behaviour and must be checked independently.
%  ===================================================================
fprintf('\n=== Part 4: N-sensitivity check, ALL FOUR objectives ===\n');
T_slow = 2*pi/omega1;
Ns_to_check = round([1, 2, 5, 10] * T_slow / dt);
fprintf('Checking N corresponding to 1, 2, 5, 10 cycles of the slower frequency: N = %s\n', mat2str(Ns_to_check));

objectiveFns = struct( ...
    'name', {'trace', 'logDet', 'lambdaMin', 'kappa (minimised)'}, ...
    'fn', {@(Wo) gramianMetrics(Wo).trace, ...
           @(Wo) gramianMetrics(Wo).logDet, ...
           @(Wo) gramianMetrics(Wo).lambdaMin, ...
           @(Wo) -gramianMetrics(Wo).kappa} ...
);

stableAcrossN = true(1, length(objectiveFns));
for o = 1:length(objectiveFns)
    fprintf('\n  Objective: %s\n', objectiveFns(o).name);
    selectionsThisObjective = cell(length(Ns_to_check),1);
    for k = 1:length(Ns_to_check)
        selectionsThisObjective{k} = sort(greedySensorSelection(candidateRows, A, Ns_to_check(k), objectiveFns(o).fn, 4));
        fprintf('    N=%5d (%.0f cycles): selected = %s\n', Ns_to_check(k), Ns_to_check(k)*dt/T_slow, mat2str(selectionsThisObjective{k}));
    end
    for k = 2:length(selectionsThisObjective)
        if ~isequal(selectionsThisObjective{k}, selectionsThisObjective{1})
            stableAcrossN(o) = false;
        end
    end
    if stableAcrossN(o)
        fprintf('    STABLE across all tested N.\n');
    else
        fprintf('    WARNING: selection CHANGES with N for this objective.\n');
    end
end

if all(stableAcrossN)
    fprintf('\nAll four objectives'' selections are stable across N -- proceeding with N=%d.\n', Ns_to_check(3));
else
    fprintf('\nWARNING: at least one objective''s selection is N-dependent (see above).\n');
    fprintf('Proceeding with N=%d regardless, but results for any unstable objective should be\n', Ns_to_check(3));
    fprintf('interpreted with that instability explicitly in mind, not treated as settled.\n');
end
N_horizon = Ns_to_check(3);

%% ===================================================================
%  PART 5: greedy selection under all four objectives at the chosen,
%  now-checked horizon
%  ===================================================================
fprintf('\n=== Part 5: greedy selection under all four objectives (N=%d) ===\n', N_horizon);

finalists = struct('name', {}, 'selectedIdx', {}, 'scoreHistory', {}, 'locations', {});
for o = 1:length(objectiveFns)
    [selIdx, scoreHist] = greedySensorSelection(candidateRows, A, N_horizon, objectiveFns(o).fn, 4);
    finalists(o).name = objectiveFns(o).name;
    finalists(o).selectedIdx = sort(selIdx);
    finalists(o).scoreHistory = scoreHist;
    finalists(o).locations = candidateLocations(finalists(o).selectedIdx, :);
    fprintf('%-20s: sensors %s\n', objectiveFns(o).name, mat2str(finalists(o).selectedIdx));
    for s = 1:4
        fprintf('    (r=%.3f, theta=%.3f)\n', finalists(o).locations(s,1), finalists(o).locations(s,2));
    end
end

%% ===================================================================
%  PART 6: expensive validation -- run the REAL SSI pipeline on each
%  DISTINCT finalist configuration. Incident amplitudes are compensated
%  (a_components = a0./[norm1;norm2]) so the generated signal matches
%  the NORMALISED fields the Gramian search actually optimised over --
%  otherwise this step would silently reintroduce the amplitude
%  imbalance Part 3 deliberately removed. Failed configurations are
%  KEPT in the results table (as NaN with a status flag), not dropped,
%  since a theoretically-good configuration failing outright is
%  itself an important result.
%  ===================================================================
fprintf('\n=== Part 6: real SSI validation on finalist configurations ===\n');

a0 = 0.05;
a_components = a0 ./ [norm1; norm2];   % amplitude-compensated: matches the normalised search
fprintf('Amplitude-compensated a_components = %s (matches normalised Gramian search assumption)\n', mat2str(a_components));

omega_components = [omega1; omega2];
alpha_components = [alpha1; alpha2];
epsilon_components = [0; 0];
specData.omega = omega_components;
specData.a = a_components;
specData.epsilon = epsilon_components;
specData.alpha = alpha_components;
specData.H = H;
specData.R = R;
specData.binData = {binData1; binData2};

tVec = 0:dt:80;
i_hankel = 40; nMax = 12;
freqTol = 0.01; dampTol = 0.05; macTol = 0.99;

uniqueConfigs = {};
uniqueNames = {};
for o = 1:length(finalists)
    isNew = true;
    for u = 1:length(uniqueConfigs)
        if isequal(uniqueConfigs{u}, finalists(o).selectedIdx)
            uniqueNames{u} = [uniqueNames{u}, ', ', finalists(o).name]; %#ok<AGROW>
            isNew = false;
        end
    end
    if isNew
        uniqueConfigs{end+1} = finalists(o).selectedIdx; %#ok<AGROW>
        uniqueNames{end+1} = finalists(o).name; %#ok<AGROW>
    end
end
fprintf('%d distinct configuration(s) to validate:\n', length(uniqueConfigs));

resultsTable = struct('objectives', {}, 'success', {}, 'status', {}, ...
    'knownCrossMAC', {}, 'ssiCrossMAC', {}, 'matchMAC1', {}, 'matchMAC2', {}, ...
    'freqErr1', {}, 'freqErr2', {}, 'persistent1', {}, 'persistent2', {});

for u = 1:length(uniqueConfigs)
    sensorLocations = candidateLocations(uniqueConfigs{u}, :);
    nSensors = size(sensorLocations,1);
    fprintf('\n--- Configuration for [%s]: sensors %s ---\n', uniqueNames{u}, mat2str(uniqueConfigs{u}));

    knownShapes = zeros(nSensors, 2);
    for s = 1:nSensors
        knownShapes(s,1) = evaluateDeflection(binData1, sensorLocations(s,1), sensorLocations(s,2));
        knownShapes(s,2) = evaluateDeflection(binData2, sensorLocations(s,1), sensorLocations(s,2));
    end
    knownCrossMAC = computeMAC(knownShapes(:,1), knownShapes(:,2));

    y = syntheticSensorData(specData, sensorLocations, tVec);
    [Yp_ref, Yf, ~] = buildHankelMatrix(y, i_hankel);
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
    isPersistent = false(2,1);
    for c = 1:2
        for b = 1:length(branches)
            if any(abs(branches(b).frequency - omega_components(c))/omega_components(c) < freqTol) && persistence(b).isPersistent
                foundBranch(c) = b;
                isPersistent(c) = true;
                break;
            end
        end
    end

    if any(foundBranch == 0)
        fprintf('  FAILED: no persistent branch found for at least one frequency.\n');
        resultsTable(end+1) = struct('objectives', uniqueNames{u}, 'success', false, ...
            'status', 'NO PERSISTENT BRANCH', 'knownCrossMAC', knownCrossMAC, 'ssiCrossMAC', NaN, ...
            'matchMAC1', NaN, 'matchMAC2', NaN, 'freqErr1', NaN, 'freqErr2', NaN, ...
            'persistent1', isPersistent(1), 'persistent2', isPersistent(2)); %#ok<SAGROW>
        continue;
    end

    b1 = foundBranch(1); b2 = foundBranch(2);
    commonOrders = intersect(branches(b1).order, branches(b2).order);
    if isempty(commonOrders)
        fprintf('  FAILED: no common order between the two branches.\n');
        resultsTable(end+1) = struct('objectives', uniqueNames{u}, 'success', false, ...
            'status', 'NO COMMON ORDER', 'knownCrossMAC', knownCrossMAC, 'ssiCrossMAC', NaN, ...
            'matchMAC1', NaN, 'matchMAC2', NaN, 'freqErr1', NaN, 'freqErr2', NaN, ...
            'persistent1', true, 'persistent2', true); %#ok<SAGROW>
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
    matchMAC1 = computeMAC(shape1, knownShapes(:,1));
    matchMAC2 = computeMAC(shape2, knownShapes(:,2));
    freqErr1 = abs(branches(b1).frequency(idx1) - omega1)/omega1;
    freqErr2 = abs(branches(b2).frequency(idx2) - omega2)/omega2;

    fprintf('  SUCCESS. knownCrossMAC=%.4f, ssiCrossMAC=%.4f, matchMAC=[%.4f,%.4f], freqErr=[%.2e,%.2e]\n', ...
        knownCrossMAC, ssiCrossMAC, matchMAC1, matchMAC2, freqErr1, freqErr2);

    resultsTable(end+1) = struct('objectives', uniqueNames{u}, 'success', true, 'status', 'OK', ...
        'knownCrossMAC', knownCrossMAC, 'ssiCrossMAC', ssiCrossMAC, ...
        'matchMAC1', matchMAC1, 'matchMAC2', matchMAC2, ...
        'freqErr1', freqErr1, 'freqErr2', freqErr2, ...
        'persistent1', true, 'persistent2', true); %#ok<SAGROW>
end

%% ===================================================================
%  PART 7: final summary table, including failures
%  ===================================================================
fprintf('\n=== Part 7: summary (including any failures) ===\n');
fprintf('%-25s %-18s %12s %12s %10s %10s %10s %10s\n', ...
    'Objective(s)', 'Status', 'knownXMAC', 'ssiXMAC', 'match1', 'match2', 'fErr1', 'fErr2');
for r = 1:length(resultsTable)
    fprintf('%-25s %-18s %12.4f %12.4f %10.4f %10.4f %10.2e %10.2e\n', ...
        resultsTable(r).objectives, resultsTable(r).status, resultsTable(r).knownCrossMAC, ...
        resultsTable(r).ssiCrossMAC, resultsTable(r).matchMAC1, resultsTable(r).matchMAC2, ...
        resultsTable(r).freqErr1, resultsTable(r).freqErr2);
end
fprintf('\nLower knownCrossMAC/ssiCrossMAC = better geometric distinguishability. A NaN row\n');
fprintf('with a non-OK status is a real, reportable failure, not a missing data point.\n');