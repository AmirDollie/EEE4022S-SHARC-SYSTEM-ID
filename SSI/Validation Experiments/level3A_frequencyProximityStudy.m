
% LEVEL 3A: frequency resolution study.
%  As omega2 approaches omega1, at
% what separation does SSI stop producing two DISTINCT persistent
% branches? 
%
% NOTE ON THE REFERENCE SCALE: 2*pi/T (T = record duration) is a
% Fourier-bin / finite-record reference scale used to BRACKET this
% experiment, not a predicted SSI resolution limit. Parametric/subspace
% methods like SSI can, under favourable conditions, resolve frequencies
% more finely than this naive bound -> the actual empirical question is
% whether the observed transition occurs near, above, or below this
% scale, not an assumption that it must occur there.
%
% CONFOUND, acknowledged rather than eliminated: omega2 changing also
% changes alpha2, and therefore the Forward Model's spatial response at
% omega2 (consistent with the real physics, per Level 2A). This
% experiment therefore tests frequency proximity under physically
% consistent, simultaneously-changing spatial responses, NOT a
% mathematically pure measurement of temporal resolution alone in
% isolation from spatial content. This is arguably the more relevant test
% for the actual application, but should not be read as isolating
% temporal resolution as a single independent variable.
%
% Uses freqTolTight = 0.001 (vs. the standard 0.01) specifically so the
% branch-matching gate's own width does not become wider than the
% separations being tested, which would confuse "can extraction resolve
% two frequencies" with "does the matching gate wrongly fuse them".

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'JONSWAP'));
clear all, close all, clc

% Set up forward model inputs (the same floe as usualllll)
H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;
g = 9.81;

% The same basic sensor locations as well (reach row is ordered (r, theta))
sensorLocations = [0, 0; 0.5*R, 0; 0.9*R, 0; 0.9*R, pi];
nSensors = size(sensorLocations,1);

% Chosen an 80s long record for bow
dt = 0.02;
tVec = 0:dt:80;
T_record = tVec(end);
% If i used a Fourier method, this would be the max resolution. SSI MAY
% perform better though, we will see:
referenceScale = 2*pi/T_record;
fprintf('Record duration T=%.1fs -> Fourier-bin reference scale 2*pi/T = %.4f rad/s\n', T_record, referenceScale);
fprintf('(bracketing reference only -> NOT assumed to be where SSI actually transitions)\n\n');

%base frequency, will build other ones around it
omega1 = 2.70;
freqTolTight = 0.001;
dampTol = 0.05; macTol = 0.99;
i = 40; nMax = 12;

% Helper functions:
function row = runOneDelta(delta, omega1, H, beta, gamma, R, nu, M, P, N, g, ...
        sensorLocations, nSensors, dt, tVec, i, nMax, freqTolTight, dampTol, macTol)
    omega2 = omega1 + delta;
    omega_components = [omega1; omega2];
    a_components = [0.05; 0.05];
    epsilon_components = [0; 0];
    alpha_components = H * omega_components.^2 / g;

    binData = cell(2,1);
    for c = 1:2
        binData{c} = precomputeDeflectionData(alpha_components(c), beta, gamma, R, nu, M, P, N);
    end
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
        matches{n} = classifyPoleStability(results(n-1), results(n), freqTolTight, dampTol, macTol);
    end
    branches = buildModalBranches(results, matches);
    persistence = computeBranchPersistence(branches, 3, 2);

    found1 = false; found2 = false; b1 = 0; b2 = 0;
    for b = 1:length(branches)
        if any(abs(branches(b).frequency - omega1)/omega1 < freqTolTight) && persistence(b).isPersistent
            found1 = true; b1 = b;
        end
        if any(abs(branches(b).frequency - omega2)/omega2 < freqTolTight) && persistence(b).isPersistent
            found2 = true; b2 = b;
        end
    end
    distinct = found1 && found2 && (b1 ~= b2);

    row.delta = delta; row.omega2 = omega2;
    row.found1 = found1; row.found2 = found2; row.distinct = distinct;
    row.err1 = NaN; row.err2 = NaN; row.orderRange1 = [NaN NaN]; row.orderRange2 = [NaN NaN];
    if b1 > 0
        row.err1 = abs(branches(b1).frequency(end) - omega1)/omega1;
        row.orderRange1 = [branches(b1).order(1), branches(b1).order(end)];
    end
    if b2 > 0
        row.err2 = abs(branches(b2).frequency(end) - omega2)/omega2;
        row.orderRange2 = [branches(b2).order(1), branches(b2).order(end)];
    end
end

function printRows(deltas, headerLabel)
    fprintf('\n%s\n', headerLabel);
    fprintf('%8s %10s %8s %8s %10s %14s %14s %12s %12s\n', ...
        'delta', 'omega2', 'found1', 'found2', 'distinct', 'err1', 'err2', 'orders1', 'orders2');
end

%% Stage 1: coarse sweep
%deltasCoarse = [0.20, 0.15, 0.10, 0.075, 0.05, 0.03, 0.02, 0.01];
% The above passed with 100% accuracy, making it smaller now:
deltasCoarseExtended = [0.20, 0.15, 0.10, 0.075, 0.05, 0.03, 0.02, 0.01, ...
                        0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001];
fprintf('=== Stage 1: coarse sweep ===\n');
fprintf('%8s %10s %8s %8s %10s %14s %14s %14s %14s\n', ...
    'delta', 'omega2', 'found1', 'found2', 'distinct', 'err1', 'err2', 'orders1', 'orders2');
lastDistinct = NaN; firstNonDistinct = NaN;
for d = 1:length(deltasCoarseExtended)
    row = runOneDelta(deltasCoarseExtended(d), omega1, H, beta, gamma, R, nu, M, P, N, g, ...
        sensorLocations, nSensors, dt, tVec, i, nMax, freqTolTight, dampTol, macTol);
    fprintf('%8.4f %10.4f %8s %8s %10s %14.3e %14.3e %14s %14s\n', ...
        row.delta, row.omega2, mat2str(row.found1), mat2str(row.found2), mat2str(row.distinct), ...
        row.err1, row.err2, mat2str(row.orderRange1), mat2str(row.orderRange2));
    if row.distinct, lastDistinct = row.delta; end
    if ~row.distinct && isnan(firstNonDistinct), firstNonDistinct = row.delta; end
end

fprintf('\nCoarse transition bracket: last distinct delta ~ %.4f, first non-distinct delta ~ %.4f\n', ...
    lastDistinct, firstNonDistinct);

%% Stage 2: refine around the observed transition, if one was found
if ~isnan(lastDistinct) && ~isnan(firstNonDistinct) && firstNonDistinct < lastDistinct
    deltasFine = linspace(firstNonDistinct, lastDistinct, 6);
    fprintf('\n=== Stage 2: refined sweep between %.4f and %.4f ===\n', firstNonDistinct, lastDistinct);
    fprintf('%8s %10s %8s %8s %10s %14s %14s %14s %14s\n', ...
        'delta', 'omega2', 'found1', 'found2', 'distinct', 'err1', 'err2', 'orders1', 'orders2');
    for d = 1:length(deltasFine)
        row = runOneDelta(deltasFine(d), omega1, H, beta, gamma, R, nu, M, P, N, g, ...
            sensorLocations, nSensors, dt, tVec, i, nMax, freqTolTight, dampTol, macTol);
        fprintf('%8.4f %10.4f %8s %8s %10s %14.3e %14.3e %14s %14s\n', ...
            row.delta, row.omega2, mat2str(row.found1), mat2str(row.found2), mat2str(row.distinct), ...
            row.err1, row.err2, mat2str(row.orderRange1), mat2str(row.orderRange2));
    end
else
    fprintf('\nNo clean transition bracketed in the coarse sweep -- consider widening deltasCoarse.\n');
end

fprintf('\nCompare the empirical transition point against the reference scale %.4f rad/s printed above.\n', referenceScale);