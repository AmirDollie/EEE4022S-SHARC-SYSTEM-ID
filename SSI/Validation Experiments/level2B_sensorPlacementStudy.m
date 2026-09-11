%% level2B_sensorPlacementStudy.m
%
% LEVEL 2B: given a frequency pair the Forward Model has ALREADY shown to
% be genuinely spatially distinguishable (2.70/6.00 rad/s,
% knownCrossMAC=0.1314 from Level 2A), does sensor PLACEMENT preserve or
% enhance that distinguishability? Fixed frequencies throughout; only
% sensor locations vary.
%
% This deliberately follows Level 2A rather than preceding it, per the
% established sequencing: searching for a good sensor placement is only
% well-posed once it's confirmed the underlying physics CAN be
% distinguished somewhere -> otherwise a placement search could spend
% effort hunting for something that doesn't exist at any placement.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'JONSWAP'));
clear all, close all, clc

H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P = 10; N = 10;
g = 9.81;

omega_components = [2.70; 6.00];
a_components = [0.05; 0.05];
epsilon_components = [0; 0];
alpha_components = H * omega_components.^2 / g;

binData = cell(2,1);
for c = 1:2
    binData{c} = precomputeDeflectionData( ...
        alpha_components(c), beta, gamma, R, nu, M, P, N);
end

specData.omega = omega_components;
specData.a = a_components;
specData.epsilon = epsilon_components;
specData.alpha = alpha_components;
specData.H = H;
specData.R = R;
specData.binData = binData;

%% Candidate sensor configurations, all with 4 sensors (Robyn's count),
% varying ONLY placement. Configuration 1 is the ORIGINAL layout used in
% every prior test, kept as the baseline for comparison.

configs = {
    [0, 0; ...
     0.5*R, 0; ...
     0.9*R, 0; ...
     0.9*R, pi], ...                              % 1: original baseline

    [0.9*R, 0; ...
     0.9*R, pi/2; ...
     0.9*R, pi; ...
     0.9*R, 3*pi/2], ...                          % 2: all near-edge, spread in theta

    [0.3*R, 0; ...
     0.6*R, 0; ...
     0.9*R, 0; ...
     0.95*R, pi/4], ...                           % 3: radial spread, one theta, one off-axis

    [0.9*R, 0; ...
     0.9*R, pi/4; ...
     0.9*R, pi/2; ...
     0.9*R, 3*pi/4]                               % 4: near-edge, narrower theta spread
};

configNames = {
    'Original (baseline)', ...
    'Edge, 4 quadrants', ...
    'Radial spread + 1 off-axis', ...
    'Edge, narrow theta spread'
};

fprintf('%-30s %16s %16s %10s\n', ...
    'Configuration', 'FM cross-MAC', 'SSI cross-MAC', 'agree?');

dt = 0.02;
tVec = 0:dt:80;
i = 40;
nMax = 12;

for cfg = 1:length(configs)

    sensorLocations = configs{cfg};
    nSensors = size(sensorLocations,1);

    %% Forward Model spatial response at each sensor

    knownShapes = zeros(nSensors, 2);

    for c = 1:2
        for s = 1:nSensors
            knownShapes(s,c) = evaluateDeflection( ...
                specData.binData{c}, ...
                sensorLocations(s,1), ...
                sensorLocations(s,2));
        end
    end

    % Ground-truth spatial distinguishability for this sensor layout
    knownCrossMAC = computeMAC( ...
        knownShapes(:,1), ...
        knownShapes(:,2));

    %% Generate synthetic multi-sensor response

    y = syntheticSensorData( ...
        specData, sensorLocations, tVec);

    %% SSI

    [Yp_ref, Yf, ~] = buildHankelMatrix(y, i);
    P_proj = hankelProjection(Yp_ref, Yf);

    warning('off', 'modalParameters:realEigenvalue');

    results = sweepModelOrders( ...
        P_proj, nSensors, dt, nMax);

    warning('on', 'modalParameters:realEigenvalue');

    %% Pole stability and branch construction

    matches = cell(nMax,1);

    for n = 2:nMax
        matches{n} = classifyPoleStability( ...
            results(n-1), ...
            results(n), ...
            0.01, ...     % frequency tolerance
            0.05, ...     % damping tolerance
            0.99);        % MAC threshold
    end

    branches = buildModalBranches( ...
        results, matches);

    persistence = computeBranchPersistence( ...
        branches, ...
        3, ...            % minimum number of orders
        2);              % maximum order gap

    %% Find persistent branches corresponding to known frequencies

    foundBranch = zeros(2,1);

    for c = 1:2

        for b = 1:length(branches)

            frequencyMatch = any( ...
                abs(branches(b).frequency - omega_components(c)) ...
                / omega_components(c) < 0.01);

            if frequencyMatch && persistence(b).isPersistent
                foundBranch(c) = b;
                break;
            end

        end

    end

    %% Check that both frequencies were successfully identified

    if any(foundBranch == 0)

        fprintf('%-30s %16.4f %16s %10s\n', ...
            configNames{cfg}, ...
            knownCrossMAC, ...
            'N/A', ...
            'NO BRANCH');

        continue;
    end

    %% Find a common model order for both persistent branches
    %
    % This avoids assuming that the two branches necessarily terminate
    % at exactly the same model order.

    commonOrders = intersect( ...
        branches(foundBranch(1)).order, ...
        branches(foundBranch(2)).order);

    if isempty(commonOrders)

        fprintf('%-30s %16.4f %16s %10s\n', ...
            configNames{cfg}, ...
            knownCrossMAC, ...
            'N/A', ...
            'NO COMMON ORDER');

        continue;
    end

    % Evaluate both modes at the highest model order where BOTH
    % persistent branches are present.

    n_eval = commonOrders(end);

    %% Locate the branch entries at the common model order

    idx1 = find( ...
        branches(foundBranch(1)).order == n_eval, ...
        1);

    idx2 = find( ...
        branches(foundBranch(2)).order == n_eval, ...
        1);

    %% Match each branch frequency to the corresponding SSI mode

    [~, col1] = min(abs( ...
        results(n_eval).frequency ...
        - branches(foundBranch(1)).frequency(idx1)));

    [~, col2] = min(abs( ...
        results(n_eval).frequency ...
        - branches(foundBranch(2)).frequency(idx2)));

    %% Extract the two SSI spatial response vectors

    shapeSSI1 = results(n_eval).modeShapes(:, col1);
    shapeSSI2 = results(n_eval).modeShapes(:, col2);

    %% Frequency-specific SSI-vs-Forward-Model matching MAC
    %
    % These test whether SSI recovers the correct spatial response
    % associated with each known frequency.

    matchingMAC1 = computeMAC( ...
        shapeSSI1, ...
        knownShapes(:,1));

    matchingMAC2 = computeMAC( ...
        shapeSSI2, ...
        knownShapes(:,2));

    %% SSI cross-MAC
    %
    % This directly compares the two SSI-recovered spatial signatures.

    ssiCrossMAC = computeMAC( ...
        shapeSSI1, ...
        shapeSSI2);

    %% Compare SSI spatial distinguishability with Forward Model truth

    agree = abs(ssiCrossMAC - knownCrossMAC) < 0.01;

    %% Main summary output

    fprintf('%-30s %16.4f %16.4f %10s\n', ...
        configNames{cfg}, ...
        knownCrossMAC, ...
        ssiCrossMAC, ...
        mat2str(agree));

    %% Additional validation information

    fprintf('    Model order = %d\n', n_eval);

    fprintf(['    Matching MAC: omega=%.2f -> %.4f, ' ...
             'omega=%.2f -> %.4f\n'], ...
        omega_components(1), ...
        matchingMAC1, ...
        omega_components(2), ...
        matchingMAC2);

end

%% Interpretation

fprintf('\nLooking for:\n');

fprintf(['(a) meaningful variation in knownCrossMAC across sensor ' ...
         'configurations --\n']);

fprintf(['    demonstrates that sensor placement changes the ' ...
         'distinguishability of\n']);

fprintf(['    the two frequency-dependent spatial response signatures;\n']);

fprintf(['(b) matching MACs close to 1 -- confirms SSI recovers the ' ...
         'correct spatial\n']);

fprintf(['    response at each identified frequency;\n']);

fprintf(['(c) ssiCrossMAC tracking knownCrossMAC -- confirms SSI ' ...
         'faithfully preserves\n']);

fprintf(['    the physical spatial distinguishability for each sensor ' ...
         'configuration.\n']);