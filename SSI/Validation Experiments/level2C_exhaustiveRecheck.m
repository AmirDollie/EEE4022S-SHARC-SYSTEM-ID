%% level2C_exhaustiveRecheck.m
%
% Re-checks Level 2C's initial greedy-only result exhaustively, per
% external review, COMPARING OBJECTIVE SCORES throughout, not sensor
% index sets -- on a grid with real angular symmetry, two different
% configurations can be tied or near-tied, in which case comparing
% indices would wrongly report "different optimum" for what is actually
% a tie. Greedy results are recomputed fresh here, not hard-coded, so
% this script stays correct if anything upstream changes.
%
% Part B uses the FULL set of tied-optimal configurations at each
% horizon (via exhaustiveSensorSelection's allScores output) and
% intersects across all tested horizons -- an earlier draft's
% diagonal-vs-row-max check was logically unable to detect genuine
% horizon dependence, since each horizon's own exhaustive-search result
% is, by construction, already that row's maximum.

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

%% Rebuild exactly the same candidate grid, A, and normalised fields as
% the original Level 2C run
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
knownShapesAll = [phi1_tilde, phi2_tilde];   % for cheap FM cross-MAC lookups if needed

objectiveFns = struct( ...
    'name', {'trace', 'logDet', 'lambdaMin', 'kappa (minimised)'}, ...
    'fn', {@(Wo) gramianMetrics(Wo).trace, ...
           @(Wo) gramianMetrics(Wo).logDet, ...
           @(Wo) gramianMetrics(Wo).lambdaMin, ...
           @(Wo) -gramianMetrics(Wo).kappa} ...
);

%% ===================================================================
%  PART A: exhaustive vs freshly-recomputed greedy at N=582, COMPARING
%  OBJECTIVE SCORES (not sensor index sets)
%  ===================================================================
fprintf('=== Part A: exhaustive vs greedy at N=582, comparing OBJECTIVE SCORES ===\n');
N_main = 582;

fprintf('%-25s %-16s %-16s %14s %14s %14s\n', 'Objective', 'Greedy set', 'Exhaustive set', 'GreedyScore', 'BestScore', 'Gap');
for o = 1:length(objectiveFns)
    objLabel = objectiveFns(o).name;
    if strcmp(objLabel, 'kappa (minimised)')
        objLabel = 'kappa (score=-kappa)';   % explicit: larger/less-negative score = SMALLER true kappa
    end

    [greedyIdx, ~] = greedySensorSelection(candidateRows, A, N_main, objectiveFns(o).fn, 4);
    greedyIdx = sort(greedyIdx);
    Wo_greedy = finiteObservabilityGramian(A, candidateRows(greedyIdx,:), N_main);
    greedyScore = objectiveFns(o).fn(Wo_greedy);

    [exhIdx, exhScore] = exhaustiveSensorSelection(candidateRows, A, N_main, objectiveFns(o).fn, 4);

    gap = exhScore - greedyScore;
    tol = 1e-10 * max(1, abs(exhScore));
    tiedOrEqual = abs(gap) <= tol;

    gapStr = iif(tiedOrEqual, 'TIED', sprintf('%.3e', gap));   % scientific notation: a small
                                                                  % but real gap should not print as 0.0000
    fprintf('%-25s %-16s %-16s %14.4f %14.4f %14s\n', objLabel, mat2str(greedyIdx), mat2str(exhIdx), ...
        greedyScore, exhScore, gapStr);
end
fprintf('\n"TIED" means greedy achieved the SAME SCORE as the global grid optimum (over this\n');
fprintf('33-point candidate grid, not all continuous locations); a numeric Gap means greedy\n');
fprintf('genuinely left objective value on the table. For "kappa (score=-kappa)", a larger\n');
fprintf('printed score means a SMALLER (better) true condition number.\n');

%% ===================================================================
%  PART B: intersection of GLOBALLY OPTIMAL configuration sets across
%  horizons, using the FULL set of tied optima per horizon (allScores),
%  not just one representative index -- directly answers whether a
%  configuration is globally optimal at every tested horizon
%  ===================================================================
fprintf('\n=== Part B: intersection of GLOBALLY OPTIMAL configuration sets across horizons ===\n');
T_slow = 2*pi/omega1;
Ns_to_check = round([1, 2, 5, 10] * T_slow / dt);
nH = length(Ns_to_check);

combos = nchoosek(1:numCandidates, 4);

for objIdx = [3, 4]
    fprintf('\nObjective: %s\n', objectiveFns(objIdx).name);
    optimalMask = false(size(combos,1), nH);

    for j = 1:nH
        [exhIdx, bestScore, allScores] = exhaustiveSensorSelection(candidateRows, A, Ns_to_check(j), objectiveFns(objIdx).fn, 4);
        tol = 1e-6 * max(1, abs(bestScore));
        optimalMask(:,j) = abs(allScores - bestScore) <= tol;
        fprintf('  N=%5d: representative optimum = %s (score=%.6f), tied optima at this horizon = %d\n', ...
            Ns_to_check(j), mat2str(exhIdx), bestScore, sum(optimalMask(:,j)));
    end

    commonOptimal = all(optimalMask, 2);
    if any(commonOptimal)
        commonRows = find(commonOptimal);
        fprintf('  -> %d configuration(s) are globally optimal at ALL %d tested horizons.\n', length(commonRows), nH);
        fprintf('     Example common optimum: %s\n', mat2str(combos(commonRows(1),:)));
        fprintf('     Apparent index changes across horizons are therefore a TIE/SYMMETRY ARTIFACT,\n');
        fprintf('     not genuine horizon-dependence of the true (grid) optimum.\n');
    else
        fprintf('  -> NO configuration is globally optimal at all %d tested horizons.\n', nH);
        fprintf('     This is direct evidence the %s objective genuinely has a horizon-dependent\n', objectiveFns(objIdx).name);
        fprintf('     global optimum over this 33-point candidate grid -- not a greedy artifact,\n');
        fprintf('     and not merely a tie/symmetry effect.\n');
    end
end

fprintf('\n=== Done. Part A''s Gap column and Part B''s common-optimum intersection are the\n');
fprintf('two results that settle whether greedy lost value, and whether the earlier\n');
fprintf('lambdaMin/kappa horizon-dependence is real or an artifact. ===\n');

function out = iif(cond, a, b)
    if cond, out = a; else, out = b; end
end