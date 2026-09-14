%% level3B_amplitudeImbalanceStudy.m
%
% LEVEL 3B: amplitude imbalance study, isolating amplitude as a SINGLE
% variable. Uses 2.70/6.00 rad/s i.e. the pair already established (Level
% 2A/2B, Level 3A) as easy on both the spatial axis (MAC_FM,cross=0.1314,
% far from the near-degenerate 0.9995 seen for 2.70/2.80) and the
% temporal axis (well above the ~0.003 rad/s resolution limit found in
% 3A). 
% 
% This is deliberate: starting from a pair that is NOT already
% difficult on any other axis means a disappearing component here can be
% attributed to amplitude specifically. A combined-difficulty version
% (2.70/2.80, amplitude imbalance ON TOP OF an already near-degenerate
% spatial pair) is a genuinely different experiment, reserved for Level
% 3C, not folded in here.
%
% Ratio swept MULTIPLICATIVELY (amplitude imbalance is inherently a ratio
% effect), coarse-then-refine staging matching 3A (refine LOGARITHMICALLY
% if a transition is bracketed, since the sweep itself is logarithmic).
%
% Tracks frequency error for BOTH branches at every ratio (not just
% found/not-found), since the possible progression:

%   amplitude decreases -> frequency estimate biased -> spatial shape
%   degrades -> branch loses persistence

% has three observable stages before the final "found" flag flips; only
% checking "found" would only ever reveal the last of these.
%
% INTERPRETIVE NOTE: matchMAC/crossMAC are DIAGNOSTICS, only computable
% while both branches are still found (once a branch disappears, its
% shape cannot be evaluated at all). Any observed MAC degradation before
% branch loss is evidence that spatial identification degrades before
% frequency detection fails at this configuration -- it does not
% establish MAC degradation as the CAUSE of the eventual branch loss,
% only as an earlier-observable symptom.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'JONSWAP'));

H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;
g = 9.81;

% The usual sensor locations in a line
sensorLocations = [0, 0; 0.5*R, 0; 0.9*R, 0; 0.9*R, pi];
nSensors = size(sensorLocations,1);

dt = 0.02; tVec = 0:dt:80;
i = 40; nMax = 12;
freqTol = 0.01; dampTol = 0.05; macTol = 0.99;

omega_components = [2.70; 6.00];
alpha_components = H * omega_components.^2 / g;
binData = cell(2,1);
for c = 1:2
    binData{c} = precomputeDeflectionData(alpha_components(c), beta, gamma, R, nu, M, P, N);
end

knownShapes = zeros(nSensors, 2);
for c = 1:2
    for s = 1:nSensors
        knownShapes(s,c) = evaluateDeflection(binData{c}, sensorLocations(s,1), sensorLocations(s,2));
    end
end
knownCrossMAC = computeMAC(knownShapes(:,1), knownShapes(:,2));
fprintf('Baseline (this pair, equal amplitude): known FM cross-MAC = %.4f\n\n', knownCrossMAC);

function row = runOneRatio(a1, a2, omega_components, alpha_components, binData, knownShapes, ...
        H, R, sensorLocations, nSensors, dt, tVec, i, nMax, freqTol, dampTol, macTol)
    a_components = [a1; a2];
    epsilon_components = [0; 0];
    specData.omega = omega_components; specData.a = a_components;
    specData.epsilon = epsilon_components; specData.alpha = alpha_components;
    specData.H = H; specData.R = R; specData.binData = binData;

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

    row.found1 = false; row.found2 = false;
    row.err1 = NaN; row.err2 = NaN;
    row.matchMAC1 = NaN; row.matchMAC2 = NaN; row.crossMAC = NaN;
    b1 = 0; b2 = 0;
    for b = 1:length(branches)
        if any(abs(branches(b).frequency - omega_components(1))/omega_components(1) < freqTol) && persistence(b).isPersistent
            row.found1 = true; b1 = b;
        end
        if any(abs(branches(b).frequency - omega_components(2))/omega_components(2) < freqTol) && persistence(b).isPersistent
            row.found2 = true; b2 = b;
        end
    end
    if b1 > 0
        row.err1 = abs(branches(b1).frequency(end) - omega_components(1)) / omega_components(1);
    end
    if b2 > 0
        row.err2 = abs(branches(b2).frequency(end) - omega_components(2)) / omega_components(2);
    end
    if b1 > 0 && b2 > 0
        n_last = branches(b1).order(end);
        idx1 = find(branches(b1).order == n_last, 1);
        idx2 = find(branches(b2).order == n_last, 1);
        if ~isempty(idx1) && ~isempty(idx2)
            [~, col1] = min(abs(results(n_last).frequency - branches(b1).frequency(idx1)));
            [~, col2] = min(abs(results(n_last).frequency - branches(b2).frequency(idx2)));
            s1 = results(n_last).modeShapes(:, col1);
            s2 = results(n_last).modeShapes(:, col2);
            row.matchMAC1 = computeMAC(s1, knownShapes(:,1));
            row.matchMAC2 = computeMAC(s2, knownShapes(:,2));
            row.crossMAC = computeMAC(s1, s2);
        end
    end
end

%% Forward direction: omega1 (2.70) dominant, omega2 (6.00, weaker) shrinking
ratiosForward = [1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, ...
                  1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12];
a1 = 0.05;

fprintf('=== Forward: a1 fixed (omega1=2.70), a2/a1 shrinking (omega2=6.00 weakening) ===\n');
fprintf('%10s %10s %10s %8s %8s %11s %11s %10s %10s %10s\n', ...
    'a2/a1', 'a1', 'a2', 'found1', 'found2', 'err1', 'err2', 'matchMAC1', 'matchMAC2', 'crossMAC');
for r = 1:length(ratiosForward)
    a2 = a1 * ratiosForward(r);
    row = runOneRatio(a1, a2, omega_components, alpha_components, binData, knownShapes, ...
        H, R, sensorLocations, nSensors, dt, tVec, i, nMax, freqTol, dampTol, macTol);
    fprintf('%10.4f %10.4f %10.6f %8s %8s %11.3e %11.3e %10.4f %10.4f %10.4f\n', ...
        ratiosForward(r), a1, a2, mat2str(row.found1), mat2str(row.found2), ...
        row.err1, row.err2, row.matchMAC1, row.matchMAC2, row.crossMAC);
end

%% Reciprocal control: omega2 (6.00) dominant, omega1 (2.70, weaker) shrinking
fprintf('\n=== Reciprocal: a2 fixed (omega2=6.00), a1/a2 shrinking (omega1=2.70 weakening) ===\n');
fprintf('%10s %10s %10s %8s %8s %11s %11s %10s %10s %10s\n', ...
    'a1/a2', 'a1', 'a2', 'found1', 'found2', 'err1', 'err2', 'matchMAC1', 'matchMAC2', 'crossMAC');
a2_fixed = 0.05;
for r = 1:length(ratiosForward)
    a1_weak = a2_fixed * ratiosForward(r);
    row = runOneRatio(a1_weak, a2_fixed, omega_components, alpha_components, binData, knownShapes, ...
        H, R, sensorLocations, nSensors, dt, tVec, i, nMax, freqTol, dampTol, macTol);
    fprintf('%10.4f %10.6f %10.4f %8s %8s %11.3e %11.3e %10.4f %10.4f %10.4f\n', ...
        ratiosForward(r), a1_weak, a2_fixed, mat2str(row.found1), mat2str(row.found2), ...
        row.err1, row.err2, row.matchMAC1, row.matchMAC2, row.crossMAC);
end

fprintf('\nLooking for: the ratio where found1/found2 flips false, AND whether err (of the\n');
fprintf('weakening component) grows noticeably in the ratios just before that flip -- a\n');
fprintf('visible error-growth stage before disappearance would mirror the clean, orderly\n');
fprintf('degradation pattern already seen in 3A. If a transition is bracketed, refine\n');
fprintf('LOGARITHMICALLY (e.g. logspace) between the last-found and first-not-found ratios,\n');
fprintf('not linearly, matching the logarithmic spacing of the sweep itself.\n');