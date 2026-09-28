%% runSensorFRFSensitivity.m
% File 5b of the SSI bridge: how much do beta, gamma and R change the
% sensor acceleration FRFs themselves, before any feature is extracted?
% No feature derived from these FRFs (poles, damping, shape ratios, ...)
% can carry more parameter information than the FRFs do, so this is an
% upper bound on what any SSI-compatible sensor feature can achieve, and
% it shows WHERE in frequency and in which spatial component the
% information lives. File 6 waits on this result.
%
% METHOD
%   Central differences in log-parameters (the three live on very
%   different scales):
%       dH/dln p_i ~ [H(p_i e^eps) - H(p_i e^-eps)] / (2 eps),
%   at every frequency and sensor, for each eps in EPS_LIST (the smaller
%   one is primary; the other checks the step).
%
% THREE REPRESENTATIONS of the same 4-sensor complex response
%   raw     : H itself.
%   affine  : P_A H, the projection onto span{1, x_j} (all sensors lie on
%             the x-axis), i.e. the best fit c0 + c1 x_j. c0 is heave-LIKE
%             and c1 pitch-LIKE: they are least-squares affine coefficients,
%             not motions measured at the floe centre.
%   nonAff  : (I - P_A) H, the non-affine remainder (flexure plus following
%             the incident-wave curvature).
%   P_A and I - P_A are orthogonal, so ||P_A H||^2 + ||(I-P_A) H||^2 = ||H||^2
%   and all three are in the same (sensor-space) units.
%
% DIMENSIONLESS STACKED JACOBIAN (per representation)
%   For each frequency, rows [Re; Im] of (P dH/dln p) / ||H(omega)||,
%   stacked over frequency and multiplied by sqrt(dOmega), so that
%   F = J'J ~ integral of J(omega)' J(omega) d omega: FREQUENCY-INTEGRATED
%   information, approximately independent of how densely the band is
%   sampled. Its SVD gives a STRUCTURAL identifiability diagnostic:
%   sigma_1..3, kappa, the weak right singular vector, and
%     joint std     = sqrt(diag(inv(F))): std of ln p_i when all three are
%                     estimated together,
%     fixed-other std = 1 ./ sqrt(diag(F)): std of ln p_i if the other two
%                     were known,
%   both per unit relative FRF noise density (per sqrt(rad/s)); their ratio
%   measures parameter coupling. These are NOT the Table 8 metric
%   (different observable, different scaling); compare qualitatively
%   (e.g. "is the weak direction still almost pure gamma?"), not numerically.
%
% CUMULATIVE INFORMATION
%   F(w_max) ~ integral from omega_min to w_max of J' J d omega; the joint
%   std sqrt(diag(inv(F))) versus w_max answers "from which frequency does
%   beta become observable?" and "does extending the band keep adding
%   information, in particular an independent gamma direction?". The value
%   at w_max = omega_max equals the full-band result of section 4.
%
% OPTIMISTIC BOUND: the complex FRF assumes access to amplitude AND phase
% relative to the incident wave. Output-only SSI has less than that, so a
% parameter that is weak here will not become identifiable after feature
% extraction; one that is strong here may still be lost.
%
% HEAVE-DIP CHECK (the kR ~ 3.83 hypothesis)
%   For a rigid disc in deep water the plane-wave average over the disc is
%   2 J1(kR)/(kR), zero at kR = 3.8317, i.e. omega = sqrt(g * 3.8317 / R_phys)
%   (~7.22 rad/s at p0; the reconnaissance found a heave-like dip at 7.275).
%   If the dip is this geometric zero, d omega_dip / d ln R ~ -omega/2 and
%   the dip barely moves with beta and gamma. A local sweep measures the
%   dip location at p0 and at +-DIP_EPS in each ln p.
%
% COST: about 1 s per EMM solve in MATLAB. Default: 45 frequencies x
%   (1 + 2 eps x 3 params x 2 signs) = 585 solves, plus the dip study
%   (DIP_N x 7 = 147 solves). Roughly 12 min.
%
% OUTPUT (in Inverse Mapping/Parameter Inversion/Results/)
%   sensorFRFSensitivity_<stamp>.mat / .txt / _sensitivity.png / _cumulative.png
%
% Lives in Inverse Mapping/Parameter Inversion/SSI Twin/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));      % computeSensorFRF.m (Inverse Mapping)
clearvars -except thisDir
close all, clc

%% ---- Configuration -----------------------------------------------------------
P0 = [4.6985e-5, 1.4548e-3, 0.3830];         % [beta gamma R], reference (not the physical floe)
PARAM_NAMES = {'beta', 'gamma', 'R'};
R0 = P0(3);
SENSOR_IDX = [2 6 14 26];                    % level2C_finalCheck candidate indices (label only)
SENSORS = [0.3*R0, 0;  0.3*R0, pi;  0.5*R0, pi;  0.9*R0, 0];     % [r theta], FIXED in absolute
                                             % non-dimensional coordinates (they do not move when R
                                             % is perturbed; the outermost sensor stays inside the floe)
DEPTH = 1.88; GRAVITY = 9.81;
OMEGA = linspace(3.0, 8.5, 45);              % trusted band, dOmega = 0.125 rad/s
EPS_LIST = [1e-2 1e-3];                      % ln-parameter steps; the LAST one is primary
TRUNCATION = [50 10 10];

DIP_BAND = [6.9 7.6];                        % local sweep around the heave-like dip
DIP_N = 21;                                  % points per sweep (dOmega = 0.035)
DIP_EPS = 0.05;                              % ln-parameter step for the dip study (large: the dip
                                             % must move by more than the fitting resolution)
USE_PARALLEL = false;

%% ---- Output ------------------------------------------------------------------------
resultsDir = fullfile(thisDir, '..', 'Results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['sensorFRFSensitivity_', stamp]);
diary([baseName, '.txt']);
try
tStart = tic;
frfArgs = {'Output', 'acceleration', 'Truncation', TRUNCATION, 'WaterDepth', DEPTH, ...
    'Gravity', GRAVITY, 'UseParallel', USE_PARALLEL};
nW = numel(OMEGA); nS = size(SENSORS, 1); nP = 3; nE = numel(EPS_LIST);
if SENSORS(end, 1) > R0 * exp(-max([EPS_LIST, DIP_EPS]))
    error('runSensorFRFSensitivity:sensorOutside', ...
        'The outermost sensor would lie outside the floe when R is reduced by the largest step.');
end

fprintf('runSensorFRFSensitivity  (%s)\n', stamp);
fprintf('p0 = [beta gamma R] = [%.6g %.6g %.6g]\n', P0);
fprintf('sensors (candidate indices %s): %s\n', mat2str(SENSOR_IDX), mat2str(SENSORS, 4));
fprintf('%d frequencies in [%.2f %.2f] rad/s; ln-parameter steps %s (primary %.0e)\n\n', ...
    nW, OMEGA(1), OMEGA(end), mat2str(EPS_LIST), EPS_LIST(end));

%% ---- 1  Sweeps ------------------------------------------------------------------------
fprintf('=== 1: sweeps ===\n');
[H0, info0] = computeSensorFRF(OMEGA, SENSORS, P0, frfArgs{:});
fprintf('p0: %d solves, %.1f s\n', nW, sum(info0.runtimeSec));
dH = complex(zeros(nS, nW, nP, nE));                   % dH/dln p, per eps
for e = 1:nE
    for i = 1:nP
        pp = P0; pp(i) = P0(i) * exp(EPS_LIST(e));
        pm = P0; pm(i) = P0(i) * exp(-EPS_LIST(e));
        Hp = computeSensorFRF(OMEGA, SENSORS, pp, frfArgs{:});
        Hm = computeSensorFRF(OMEGA, SENSORS, pm, frfArgs{:});
        dH(:, :, i, e) = (Hp - Hm) / (2 * EPS_LIST(e));
        fprintf('  eps %.0e, %-5s done  [%.1f min elapsed]\n', EPS_LIST(e), PARAM_NAMES{i}, toc(tStart) / 60);
    end
end
dHp = dH(:, :, :, end);                                % primary
fprintf('\n');

%% ---- 2  Step check -------------------------------------------------------------------
fprintf('=== 2: finite-difference step check (relative change of dH/dln p between steps) ===\n');
stepCheck = NaN(1, nP);
if nE > 1
    for i = 1:nP
        a = dH(:, :, i, end); b = dH(:, :, i, end - 1);
        stepCheck(i) = norm(a(:) - b(:)) / norm(a(:));
        fprintf('  %-5s  eps %.0e vs %.0e: %.2e\n', PARAM_NAMES{i}, EPS_LIST(end), EPS_LIST(end - 1), stepCheck(i));
    end
else
    fprintf('  only one step in EPS_LIST: not checked\n');
end
fprintf('\n');

%% ---- 3  Representations and frequency-resolved sensitivity --------------------------
xS = SENSORS(:, 1) .* cos(SENSORS(:, 2));
Xa = [ones(nS, 1), xS];
PA = Xa * pinv(Xa);                                    % affine projector (4 x 4)
reps = struct('name', {'raw', 'affine', 'nonAff'}, 'P', {eye(nS), PA, eye(nS) - PA});
nR = numel(reps);
normH = vecnorm(H0, 2, 1);
S = zeros(nR, nP, nW);                                 % ||P dH/dln p|| / ||H||
for rr = 1:nR
    for i = 1:nP
        S(rr, i, :) = vecnorm(reps(rr).P * dHp(:, :, i), 2, 1) ./ normH;
    end
end
fracOfH = zeros(nR, nW);
for rr = 1:nR, fracOfH(rr, :) = vecnorm(reps(rr).P * H0, 2, 1) ./ normH; end

% affine coefficients themselves: relative sensitivity of c0 (heave-like) and c1 (pitch-like)
C0 = pinv(Xa) * H0;
Sc = zeros(2, nP, nW);
for i = 1:nP
    dC = pinv(Xa) * dHp(:, :, i);
    Sc(:, i, :) = abs(dC) ./ abs(C0);
end

fprintf('=== 3: frequency-resolved sensitivity S = ||P dH/dln p|| / ||H|| ===\n');
kk = unique(round(linspace(1, nW, 7)));
fprintf('  %-7s %-6s %s\n', 'rep', 'param', sprintf('%9.2f', OMEGA(kk)));
for rr = 1:nR
    for i = 1:nP
        fprintf('  %-7s %-6s %s\n', reps(rr).name, PARAM_NAMES{i}, sprintf('%9.2e', squeeze(S(rr, i, kk))));
    end
end
fprintf('  (share of ||H||: raw 1, affine %s, nonAff %s at the same frequencies)\n', ...
    mat2str(fracOfH(2, kk), 2), mat2str(fracOfH(3, kk), 2));
fprintf('  band-averaged S (rms over frequency):\n');
for rr = 1:nR
    fprintf('    %-7s beta %.2e   gamma %.2e   R %.2e\n', reps(rr).name, sqrt(mean(S(rr, :, :).^2, 3)));
end
fprintf('\n');

%% ---- 4  Stacked Jacobian: SVD and conditional std --------------------------------------
dOmega = OMEGA(2) - OMEGA(1);
fprintf('=== 4: dimensionless stacked Jacobian (frequency-integrated, F ~ int J''J d omega) ===\n');
fprintf('(joint std = sqrt(diag(inv(F))): all parameters estimated together; fixed-other std =\n');
fprintf(' 1/sqrt(diag(F)): the other two known. Per unit relative FRF noise density. Structural\n');
fprintf(' diagnostic only; not the Table 8 metric.)\n');
svdRes = struct();
for rr = 1:nR
    J = stackJ(reps(rr).P, dHp, normH) * sqrt(dOmega);
    [~, Sg, V] = svd(J, 'econ');
    sg = diag(Sg);
    F = J.' * J;
    jointStd = sqrt(diag(pinv(F))).';                  % all parameters estimated jointly
    fixedStd = 1 ./ sqrt(diag(F)).';                    % other parameters known
    [~, iw] = max(abs(V(:, 3)));
    svdRes(rr).name = reps(rr).name;
    svdRes(rr).sigma = sg.';
    svdRes(rr).kappa = sg(1) / sg(3);
    svdRes(rr).V = V;
    svdRes(rr).jointStd = jointStd;
    svdRes(rr).fixedStd = fixedStd;
    svdRes(rr).coupling = jointStd ./ fixedStd;
    fprintf('  %-7s sigma = %s  kappa = %.3g\n', reps(rr).name, mat2str(sg.', 3), svdRes(rr).kappa);
    fprintf('          weak direction (beta gamma R) = %s  (dominated by %s)\n', mat2str(V(:, 3).', 3), PARAM_NAMES{iw});
    fprintf('          joint std       beta %.3g  gamma %.3g  R %.3g\n', jointStd);
    fprintf('          fixed-other std beta %.3g  gamma %.3g  R %.3g\n', fixedStd);
    fprintf('          coupling ratio  %s\n', mat2str(svdRes(rr).coupling, 3));
end
fprintf('\n');

%% ---- 5  Cumulative information vs upper frequency ------------------------------------------
fprintf('=== 5: cumulative (frequency-integrated) joint std vs upper frequency ===\n');
cum = struct('omegaMax', OMEGA, 'jointStd', NaN(nR, nP, nW), 'sigmaMin', NaN(nR, nW));
for rr = 1:nR
    Jall = stackJ(reps(rr).P, dHp, normH);              % 2 nS nW x 3, rows grouped by frequency
    for k = 1:nW
        Jk = sqrt(dOmega) * Jall(1:2 * nS * k, :);      % integrated up to omega_k
        F = Jk.' * Jk;
        sv = svd(Jk);
        cum.sigmaMin(rr, k) = sv(end);
        if sv(end) > 1e-12 * sv(1)
            cum.jointStd(rr, :, k) = sqrt(diag(inv(F))).';
        end
    end
end
for rr = 1:nR
    fprintf('  %s\n  %-9s %s\n', reps(rr).name, 'omega_max', sprintf('%11s', 'beta', 'gamma', 'R'));
    for k = kk
        fprintf('  %9.2f %s\n', OMEGA(k), sprintf('%11.3g', squeeze(cum.jointStd(rr, :, k))));
    end
end
fprintf('(NaN = fewer independent directions than parameters so far)\n\n');

%% ---- 6  Heave-dip check -----------------------------------------------------------------
fprintf('=== 6: heave-like dip vs the kR = 3.8317 geometric zero ===\n');
wDip = linspace(DIP_BAND(1), DIP_BAND(2), DIP_N);
dip = struct('p0', NaN, 'plus', NaN(1, nP), 'minus', NaN(1, nP));
dip.p0 = locateDip(wDip, computeSensorFRF(wDip, SENSORS, P0, frfArgs{:}), Xa);
for i = 1:nP
    pp = P0; pp(i) = P0(i) * exp(DIP_EPS);
    pm = P0; pm(i) = P0(i) * exp(-DIP_EPS);
    dip.plus(i) = locateDip(wDip, computeSensorFRF(wDip, SENSORS, pp, frfArgs{:}), Xa);
    dip.minus(i) = locateDip(wDip, computeSensorFRF(wDip, SENSORS, pm, frfArgs{:}), Xa);
end
dip.dOmega_dlnp = (dip.plus - dip.minus) / (2 * DIP_EPS);
Rphys = R0 * DEPTH;
kdisc = 3.8317 / Rphys;
dip.bessel = sqrt(GRAVITY * kdisc * tanh(kdisc * DEPTH));
dip.besselSlopeR = -dip.bessel / 2;                    % deep water: omega ~ R^(-1/2)
fprintf('  dip at p0: %.4f rad/s   rigid-disc deep-water prediction: %.4f rad/s\n', dip.p0, dip.bessel);
fprintf('  d omega_dip / d ln p:  beta %+.4f   gamma %+.4f   R %+.4f   (prediction for R: %+.4f)\n', ...
    dip.dOmega_dlnp, dip.besselSlopeR);
fprintf('  (dip located by a parabola through the three smallest-|c0| points of a %d-point sweep)\n\n', DIP_N);

%% ---- 7  Save ----------------------------------------------------------------------------
sens = struct('p0', P0, 'paramNames', {PARAM_NAMES}, 'sensors', SENSORS, 'sensorIdx', SENSOR_IDX, ...
    'omega', OMEGA, 'H0', H0, 'dH', dH, 'epsList', EPS_LIST, 'truncation', TRUNCATION, ...
    'xSensors', xS, 'affineProjector', PA, 'S', S, 'repNames', {{reps.name}}, 'shareOfH', fracOfH, ...
    'Scoeff', Sc, 'stepCheck', stepCheck, 'svd', svdRes, 'cumulative', cum, 'dip', dip, ...
    'depth', DEPTH, 'gravity', GRAVITY);
save([baseName, '.mat'], 'sens');
fprintf('Saved: %s.mat\n', baseName);

%% ---- 8  Figures ------------------------------------------------------------------------------
colours = {[42 120 214]/255, [235 104 52]/255, [27 175 122]/255};   % validated slots 1-3
ink = [0.35 0.35 0.33]; gridInk = [0.85 0.85 0.83];
markers = {'o', 's', '^'};
mEvery = max(1, round(nW / 10));
repTitles = {'(a) raw sensors', '(b) affine part (heave-like + pitch-like)', '(c) non-affine remainder'};

fig1 = figure('Color', 'w', 'Position', [60 60 1250 420]);
yAll = S(:); yAll = yAll(yAll > 0);
yl = [10^floor(log10(min(yAll))), 10^ceil(log10(max(yAll)))];
for rr = 1:nR
    ax = subplot(1, nR, rr); hold(ax, 'on');
    h = gobjects(1, nP);
    for i = 1:nP
        h(i) = plot(ax, OMEGA, squeeze(S(rr, i, :)), '-', 'Color', colours{i}, 'LineWidth', 2, ...
            'Marker', markers{i}, 'MarkerIndices', 1:mEvery:nW, 'MarkerSize', 6, ...
            'MarkerFaceColor', colours{i}, 'MarkerEdgeColor', 'w');
    end
    xline(ax, dip.p0, ':', 'heave-like dip', 'Color', ink, 'LabelVerticalAlignment', 'bottom', ...
        'FontSize', 8, 'HandleVisibility', 'off');
    set(ax, 'YScale', 'log', 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off');
    grid(ax, 'on'); set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
    xlim(ax, OMEGA([1 end])); ylim(ax, yl);
    xlabel(ax, '\omega (rad/s)');
    if rr == 1, ylabel(ax, '||P \partialH/\partialln p|| / ||H||'); end
    title(ax, repTitles{rr}, 'FontWeight', 'normal');
    if rr == 1, legend(ax, h, {'\beta', '\gamma', 'R'}, 'Location', 'best', 'Box', 'off'); end
end
saveFigure(fig1, [baseName, '_sensitivity.png']);

fig2 = figure('Color', 'w', 'Position', [80 80 1250 420]);
for rr = 1:nR
    ax = subplot(1, nR, rr); hold(ax, 'on');
    h = gobjects(1, nP);
    for i = 1:nP
        h(i) = plot(ax, OMEGA, squeeze(cum.jointStd(rr, i, :)), '-', 'Color', colours{i}, 'LineWidth', 2, ...
            'Marker', markers{i}, 'MarkerIndices', 1:mEvery:nW, 'MarkerSize', 6, ...
            'MarkerFaceColor', colours{i}, 'MarkerEdgeColor', 'w');
    end
    set(ax, 'YScale', 'log', 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off');
    grid(ax, 'on'); set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
    xlim(ax, OMEGA([1 end]));
    xlabel(ax, 'upper frequency \omega_{max} (rad/s)');
    if rr == 1, ylabel(ax, 'joint std of ln p (per unit relative noise density)'); end
    title(ax, repTitles{rr}, 'FontWeight', 'normal');
    if rr == 1, legend(ax, h, {'\beta', '\gamma', 'R'}, 'Location', 'best', 'Box', 'off'); end
end
saveFigure(fig2, [baseName, '_cumulative.png']);

fprintf('\nDone in %.1f min.\n', toc(tStart) / 60);
catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% =================================================================================================
%  Local functions
%% =================================================================================================
function J = stackJ(P, dH, normH)
% Rows [Re; Im] of (P dH(:,k,:)) / ||H(omega_k)||, stacked frequency by frequency.
[nS, nW, nP] = size(dH);
J = zeros(2 * nS * nW, nP);
for k = 1:nW
    blk = P * reshape(dH(:, k, :), nS, nP) / normH(k);
    J((k - 1) * 2 * nS + (1:2 * nS), :) = [real(blk); imag(blk)];
end
end

function w0 = locateDip(w, H, Xa)
% Minimum of |c0| (heave-like affine coefficient) by a parabola through the
% three smallest points; NaN if the minimum sits on the sweep edge.
c0 = abs(pinv(Xa(:, :)) * H);
c0 = c0(1, :);
[~, k] = min(c0);
if k == 1 || k == numel(w)
    w0 = NaN;
    return
end
a = polyfit(w(k - 1:k + 1), c0(k - 1:k + 1), 2);
w0 = -a(2) / (2 * a(1));
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);     % R2020a+
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end