%% runSensorFRFReconnaissance.m
% File 5 of the SSI bridge: acceleration FRFs of the four baseline sensors
% at the reference point p0, over the trusted EMM band. Its job is to show
% what sensor-level structure the EMM actually produces, so that the fixed
% pole-fit windows for fitLocalPoleModel.m (file 6) can be chosen from
% data rather than assumed.
%
% WHAT IT DOES
%   1  Dense sweep: computeSensorFRF over OMEGA_RANGE (acceleration; the
%      displacement FRF comes back in info at no extra cost).
%   2  Local refinement around every detected peak of ||H_a||, so each peak
%      is resolved well enough to measure its width.
%   3  Structure: peaks and dips of ||H_a|| and of each |H_j|, with
%      prominence and an EFFECTIVE half-power width. The width is only a
%      diagnostic: overlapping modes, antiresonances and the frequency-
%      dependent hydrodynamics mean it is NOT a damping estimate. File 6
%      measures damping properly from a shared-pole fit.
%   3b Rigid-body / non-rigid split. All four sensors lie on the x-axis
%      (x = r cos(theta)), so per frequency H_j = c_heave + c_pitch * x_j + E_j,
%      fitted by least squares. c_heave and c_pitch are the rigid heave and
%      pitch (small-angle) responses; E is everything else: flexure PLUS the
%      floe following the curvature of the incident wave. A preliminary
%      Octave run (34 points) found no resonance peaks in |H_j| or ||H||
%      (smooth, monotone magnitudes; relative phases growing ~linearly with
%      omega, i.e. dominated by wave propagation across the floe), so the
%      structure, if any, is expected in c_heave, c_pitch and E.
%   4  Proposed pole-fit windows (peak +- WINDOW_HALFWIDTHS * half-width,
%      merged where they overlap, clipped to the trusted band). PROPOSED
%      only: they are frozen in file 6 after looking at these plots.
%   5  Truncation spot-check: [50 10 10] vs TRUNC_CHECK at a few in-band
%      frequencies and a few just ABOVE the validated ceiling (information
%      for a possible later band extension; nothing here depends on it).
%   6  Figures and a .mat that files 6, 8 and 9 load.
%
% THE OLD EMM FEATURES are overlaid for comparison only (from
% anchorPredictionModel.mat: A20, A01, A11 = dry-mode-projection maxima,
% A21, A00 = near-zeros). They are projections onto dry modes, not sensor
% responses, so they are NOT expected to coincide with sensor peaks.
%
% SENSORS: Level 2C kappa-long configuration, candidate indices [2 6 14 26]
% on the level2C_finalCheck grid (rGrid = [0.3 0.5 0.7 0.9]*R, 8 angles
% from 0; index 1 = centre). Stored here as EXPLICIT coordinates so the
% saved file never depends on that indexing. A named baseline, not an
% optimum for this problem. All four lie on the incidence axis (theta = 0
% or pi), so only modes symmetric about that axis are observable, which
% is everything a single incident direction excites.
%
% COST: about 1 s per EMM solve in MATLAB; the default grid plus
% refinement and the truncation check is roughly 250-300 solves (~5 min).
%
% OUTPUT (in Inverse Mapping/Parameter Inversion/Results/)
%   sensorFRFRecon_<stamp>.mat / .txt / _frf.png / _rigid.png
%
% Lives in Inverse Mapping/Parameter Inversion/SSI Twin/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));      % computeSensorFRF.m (Inverse Mapping)
clearvars -except thisDir
close all, clc

%% ---- Configuration -----------------------------------------------------------
P0 = [4.6985e-5, 1.4548e-3, 0.3830];         % [beta gamma R], reference (not the physical floe)
R0 = P0(3);
SENSOR_IDX = [2 6 14 26];                    % level2C_finalCheck candidate indices (label only)
SENSORS = [0.3*R0, 0;                        % [r theta], non-dimensional r
           0.3*R0, pi;
           0.5*R0, pi;
           0.9*R0, 0];
SENSOR_NAMES = {'s2 (0.3R, 0)', 's6 (0.3R, \pi)', 's14 (0.5R, \pi)', 's26 (0.9R, 0)'};

DEPTH = 1.88; GRAVITY = 9.81;
OMEGA_RANGE = [3.0 8.5];                     % rad/s: alpha from 1.73 to 13.85 (validated [1.7 13.85])
N_OMEGA = 221;                               % dOmega = 0.025 rad/s
TRUNCATION = [50 10 10];

REFINE_PEAKS = true;
REFINE_POINTS = 15;                          % extra solves per detected ||H|| peak
REFINE_HALFSPAN = 2;                         % refine over peak +- 2 coarse spacings (locates the peak;
                                             % the coarse grid already resolves widths >~ 0.1 rad/s)
MIN_PROMINENCE_DB = 0.3;                     % peaks/dips smaller than this (in dB) are not reported
WINDOW_HALFWIDTHS = 3;                       % proposed window = peak +- 3 effective half-widths
MIN_WINDOW_HALFSPAN = 0.15;                  % rad/s, floor for narrow or width-less peaks

TRUNC_CHECK = [60 12 12];
TRUNC_OMEGA_IN = [4.0 5.5 7.0 8.3];          % in-band spot checks
TRUNC_OMEGA_ABOVE = [8.7 9.0 9.5];           % above the validated ceiling (information only)

USE_PARALLEL = false;

%% ---- Output ---------------------------------------------------------------------
resultsDir = fullfile(thisDir, '..', 'Results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['sensorFRFRecon_', stamp]);
diary([baseName, '.txt']);
try
tStart = tic;
frfArgs = {'Output', 'acceleration', 'Truncation', TRUNCATION, 'WaterDepth', DEPTH, ...
    'Gravity', GRAVITY, 'UseParallel', USE_PARALLEL};

fprintf('runSensorFRFReconnaissance  (%s)\n', stamp);
fprintf('p0 = [beta gamma R] = [%.6g %.6g %.6g]\n', P0);
fprintf('sensors (Level 2C kappa-long, candidate indices %s):\n', mat2str(SENSOR_IDX));
for j = 1:size(SENSORS, 1)
    fprintf('  %-16s r = %.4f (%.2f R), theta = %.4f\n', SENSOR_NAMES{j}, SENSORS(j, 1), ...
        SENSORS(j, 1) / R0, SENSORS(j, 2));
end
fprintf('band [%.2f %.2f] rad/s, %d points, truncation %s\n\n', OMEGA_RANGE, N_OMEGA, mat2str(TRUNCATION));

%% ---- Old EMM features (overlay only) ----------------------------------------------------
oldFeat = struct('omega', [], 'labels', {{}}, 'isMax', []);
modelFile = fullfile(thisDir, '..', '..', 'anchorPredictionModel.mat');
if isfile(modelFile)
    am = load(modelFile);
    oldFeat.omega = am.f0(:).';
    oldFeat.labels = am.labels;
    oldFeat.isMax = [true true true false false];     % A20 A01 A11 maxima; A21 A00 near-zeros
    fprintf('old EMM features at p0 (overlay only):');
    for k = 1:numel(oldFeat.omega)
        fprintf('  %s %.3f', oldFeat.labels{k}, oldFeat.omega(k));
    end
    fprintf('\n\n');
else
    fprintf('anchorPredictionModel.mat not found: no old-feature overlay.\n\n');
end

%% ---- 1  Dense sweep -----------------------------------------------------------------------------
fprintf('=== 1: dense sweep ===\n');
omega = linspace(OMEGA_RANGE(1), OMEGA_RANGE(2), N_OMEGA);
[Ha, info] = computeSensorFRF(omega, SENSORS, P0, frfArgs{:});
Hd = info.Hdisplacement;
fprintf('%d solves, %.1f s (%.2f s per solve)\n\n', N_OMEGA, sum(info.runtimeSec), mean(info.runtimeSec));

%% ---- 2  Refine around ||H_a|| peaks -------------------------------------------------------------
dW = omega(2) - omega(1);
if REFINE_PEAKS
    fprintf('=== 2: refinement around ||H_a|| peaks ===\n');
    pk = findLocalExtrema(vecnorm(Ha, 2, 1), 'max', MIN_PROMINENCE_DB);
    wNew = [];
    for q = 1:numel(pk.idx)
        c = omega(pk.idx(q));
        wNew = [wNew, linspace(c - REFINE_HALFSPAN * dW, c + REFINE_HALFSPAN * dW, REFINE_POINTS)]; %#ok<AGROW>
    end
    wNew = wNew(wNew >= OMEGA_RANGE(1) & wNew <= OMEGA_RANGE(2));
    wNew = setdiff(round(wNew * 1e10) / 1e10, round(omega * 1e10) / 1e10);   % drop duplicates of the coarse grid
    if ~isempty(wNew)
        [HaN, infoN] = computeSensorFRF(wNew, SENSORS, P0, frfArgs{:});
        [omega, ord] = sort([omega, wNew]);
        Ha = [Ha, HaN]; Ha = Ha(:, ord);
        Hd = [Hd, infoN.Hdisplacement]; Hd = Hd(:, ord);
        fprintf('%d peaks refined with %d extra solves (%.1f s)\n\n', numel(pk.idx), numel(wNew), sum(infoN.runtimeSec));
    else
        fprintf('no peaks to refine\n\n');
    end
end
nW = numel(omega);
alpha = DEPTH * omega.^2 / GRAVITY;

%% ---- 3  Structure: peaks, dips, effective widths --------------------------------------------------
fprintf('=== 3: response structure (acceleration) ===\n');
normA = vecnorm(Ha, 2, 1);
normD = vecnorm(Hd, 2, 1);
[~, refSensor] = max(mean(abs(Ha), 2));                 % largest average response: phase reference
struct3 = struct();
struct3.norm = describeCurve(omega, normA, MIN_PROMINENCE_DB);
struct3.normDisp = describeCurve(omega, normD, MIN_PROMINENCE_DB);
for j = 1:size(SENSORS, 1)
    struct3.sensor(j) = describeCurve(omega, abs(Ha(j, :)), MIN_PROMINENCE_DB); %#ok<SAGROW>
end
printCurve('||H_a|| (all sensors)', struct3.norm, oldFeat);
% relative phase slope: propagation check. The incident wave travels from
% +x towards -x (Montiel 2.2.2) under exp(+i w t), i.e. exp(i(k x + w t)), so a
% floe that simply follows it has phase(H_j / H_ref) ~ +k (x_j - x_ref).
xS = SENSORS(:, 1) .* cos(SENSORS(:, 2));              % non-dimensional x
relPh = unwrap(angle(Ha ./ Ha(refSensor, :)), [], 2);
kDeep = omega.^2 / GRAVITY * DEPTH;                    % non-dimensional deep-water wavenumber
fprintf('relative phase vs deep-water propagation k (x_j - x_ref), at omega = %.2f / %.2f / %.2f rad/s:\n', ...
    omega(1), omega(round(nW / 2)), omega(end));
for j = setdiff(1:size(SENSORS, 1), refSensor)
    kk = [1, round(nW / 2), nW];
    fprintf('  %-16s phase %s pi   propagation %s pi\n', SENSOR_NAMES{j}, mat2str(relPh(j, kk) / pi, 3), ...
        mat2str(kDeep(kk) * (xS(j) - xS(refSensor)) / pi, 3));
end
fprintf('(finite depth makes k larger than omega^2/g at low omega; the comparison is indicative only)\n\n');
printCurve('||H_eta|| (displacement, for comparison)', struct3.normDisp, oldFeat);
for j = 1:size(SENSORS, 1)
    printCurve(sprintf('|H_a| at %s', SENSOR_NAMES{j}), struct3.sensor(j), oldFeat);
end
fprintf('(width = effective half-power width of the peak; NaN = the -3 dB level is not reached\n');
fprintf(' on both sides before a higher value or the band edge, i.e. the peak is not isolated)\n\n');

%% ---- 3b  Rigid-body / non-rigid split ------------------------------------------------------------
fprintf('=== 3b: rigid heave + pitch vs non-rigid remainder (acceleration) ===\n');
Xr = [ones(size(xS)), xS];
if rank(Xr) < 2
    error('runSensorFRFReconnaissance:rigidFit', 'Sensors do not span two distinct x positions.');
end
Crig = Xr \ Ha;                                         % 2 x nW: [heave; pitch]
Enr = Ha - Xr * Crig;                                   % non-rigid remainder
fracNR = vecnorm(Enr, 2, 1) ./ normA;
struct3.heave = describeCurve(omega, abs(Crig(1, :)), MIN_PROMINENCE_DB);
struct3.pitch = describeCurve(omega, abs(Crig(2, :)), MIN_PROMINENCE_DB);
struct3.nonRigid = describeCurve(omega, vecnorm(Enr, 2, 1), MIN_PROMINENCE_DB);
printCurve('|heave|', struct3.heave, oldFeat);
printCurve('|pitch|', struct3.pitch, oldFeat);
printCurve('||non-rigid remainder||', struct3.nonRigid, oldFeat);
fprintf('non-rigid fraction ||E|| / ||H||: %.3f at %.2f, %.3f at %.2f, %.3f at %.2f rad/s (max %.3f at %.2f)\n', ...
    fracNR(1), omega(1), fracNR(round(nW / 2)), omega(round(nW / 2)), fracNR(end), omega(end), ...
    max(fracNR), omega(find(fracNR == max(fracNR), 1)));
fprintf('(E contains flexure AND quasi-static following of the incident wave curvature; it is not\n');
fprintf(' a pure elastic-mode signal)\n\n');

%% ---- 4  Proposed windows --------------------------------------------------------------------------
fprintf('=== 4: proposed pole-fit windows (freeze in file 6) ===\n');
% Peaks of ||H_a||, |heave|, |pitch| and ||E|| all propose windows.
P = struct('omega', [], 'width', []);
srcs = {'norm', 'heave', 'pitch', 'nonRigid'};
for q = 1:numel(srcs)
    P.omega = [P.omega, struct3.(srcs{q}).peaks.omega];
    P.width = [P.width, struct3.(srcs{q}).peaks.width];
end
win = zeros(0, 2);
for q = 1:numel(P.omega)
    hs = WINDOW_HALFWIDTHS * P.width(q) / 2;
    if ~isfinite(hs), hs = MIN_WINDOW_HALFSPAN; end
    hs = max(hs, MIN_WINDOW_HALFSPAN);
    win(end + 1, :) = [max(OMEGA_RANGE(1), P.omega(q) - hs), min(OMEGA_RANGE(2), P.omega(q) + hs)]; %#ok<SAGROW>
end
[win, members] = mergeIntervals(win);
for w = 1:size(win, 1)
    pts = sum(omega >= win(w, 1) & omega <= win(w, 2));
    edge = '';
    if win(w, 1) <= OMEGA_RANGE(1) + 1e-9 || win(w, 2) >= OMEGA_RANGE(2) - 1e-9
        edge = '  [touches the trusted-band edge]';
    end
    fprintf('  window %d: [%.3f %.3f] rad/s, %d grid points, contains peak(s) at %s%s\n', w, win(w, :), ...
        pts, mat2str(P.omega(members{w}), 4), edge);
end
if isempty(win)
    fprintf('  none: no peak above %.1f dB prominence in ||H_a||, |heave|, |pitch| or ||E||.\n', MIN_PROMINENCE_DB);
    fprintf('  A local shared-pole fit has nothing local to fit; file 6 must be redesigned (see header).\n');
end
fprintf('\n');

%% ---- 5  Truncation spot-check -----------------------------------------------------------------------
fprintf('=== 5: truncation spot-check, %s vs %s ===\n', mat2str(TRUNCATION), mat2str(TRUNC_CHECK));
wT = [TRUNC_OMEGA_IN, TRUNC_OMEGA_ABOVE];
wState = warning('off', 'computeSensorFRF:outsideValidatedAlpha');   % expected for the above-ceiling points
[HT1, iT1] = computeSensorFRF(wT, SENSORS, P0, frfArgs{:});
[HT2, ~] = computeSensorFRF(wT, SENSORS, P0, frfArgs{:}, 'Truncation', TRUNC_CHECK);
warning(wState);
truncCheck.omega = wT;
truncCheck.alpha = iT1.alpha;
truncCheck.relChange = vecnorm(HT2 - HT1, 2, 1) ./ vecnorm(HT1, 2, 1);
truncCheck.outsideValidated = iT1.outsideValidated;
truncCheck.truncations = [TRUNCATION; TRUNC_CHECK];
fprintf('  %8s %8s %12s %s\n', 'omega', 'alpha', 'rel change', '');
for k = 1:numel(wT)
    fprintf('  %8.3f %8.3f %12.2e %s\n', wT(k), truncCheck.alpha(k), truncCheck.relChange(k), ...
        repmat('above validated ceiling', 1, truncCheck.outsideValidated(k)));
end
fprintf('\n');

%% ---- 6  Save ------------------------------------------------------------------------------------------
recon = struct();
recon.p0 = P0;
recon.sensors = SENSORS;
recon.sensorIdx = SENSOR_IDX;
recon.sensorNames = SENSOR_NAMES;
recon.sensorGridNote = 'level2C_finalCheck grid: rGrid = [0.3 0.5 0.7 0.9]*R, thetaGrid = 8 angles from 0, index 1 = centre';
recon.omega = omega;
recon.alpha = alpha;
recon.Ha = Ha;
recon.Hdisp = Hd;
recon.truncation = TRUNCATION;
recon.depth = DEPTH; recon.gravity = GRAVITY;
recon.refSensor = refSensor;
recon.xSensors = xS;
recon.rigidCoeffs = Crig;                  % [heave; pitch] per frequency (acceleration)
recon.nonRigid = Enr;
recon.nonRigidFraction = fracNR;
recon.relPhase = relPh;
recon.structure = struct3;
recon.proposedWindows = win;
recon.windowMembers = members;
recon.truncCheck = truncCheck;
recon.oldFeatures = oldFeat;
recon.config = struct('OMEGA_RANGE', OMEGA_RANGE, 'N_OMEGA', N_OMEGA, 'REFINE_PEAKS', REFINE_PEAKS, ...
    'REFINE_POINTS', REFINE_POINTS, 'REFINE_HALFSPAN', REFINE_HALFSPAN, 'MIN_PROMINENCE_DB', MIN_PROMINENCE_DB, ...
    'WINDOW_HALFWIDTHS', WINDOW_HALFWIDTHS, 'MIN_WINDOW_HALFSPAN', MIN_WINDOW_HALFSPAN);
save([baseName, '.mat'], 'recon');
fprintf('Saved: %s.mat\n', baseName);

%% ---- 7  Figures ----------------------------------------------------------------------------------------
colours = {[42 120 214]/255, [235 104 52]/255, [27 175 122]/255, [237 161 0]/255};   % validated slots 1-4
ink = [0.35 0.35 0.33]; gridInk = [0.85 0.85 0.83]; refInk = [0.55 0.55 0.52]; winFill = [0.94 0.94 0.92];
markers = {'o', 's', '^', 'd'};
mEvery = max(1, round(nW / 12));

fig1 = figure('Color', 'w', 'Position', [60 60 1200 820]);
% (a) |H_a| per sensor
ax = subplot(2, 2, 1); hold(ax, 'on');
h = gobjects(1, size(SENSORS, 1));
for j = 1:size(SENSORS, 1)
    h(j) = plot(ax, omega, abs(Ha(j, :)), '-', 'Color', colours{j}, 'LineWidth', 2, 'Marker', markers{j}, ...
        'MarkerIndices', 1:mEvery:nW, 'MarkerSize', 6, 'MarkerFaceColor', colours{j}, 'MarkerEdgeColor', 'w');
end
set(ax, 'YScale', 'log'); styleAxes(ax, ink, gridInk);
decorate(ax, win, oldFeat, winFill, refInk);
ylabel(ax, '|H_a|  (m s^{-2} per m incident)');
title(ax, '(a) sensor acceleration response', 'FontWeight', 'normal');
legend(ax, h, SENSOR_NAMES, 'Location', 'best', 'Box', 'off');
% (b) norms, each normalised to its own max (one axis, no dual scale)
ax = subplot(2, 2, 2); hold(ax, 'on');
h2 = gobjects(1, 2);
h2(1) = plot(ax, omega, normA / max(normA), '-', 'Color', colours{1}, 'LineWidth', 2);
h2(2) = plot(ax, omega, normD / max(normD), '-', 'Color', colours{2}, 'LineWidth', 2);
plot(ax, struct3.norm.peaks.omega, struct3.norm.peaks.value / max(normA), 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', colours{1}, 'MarkerEdgeColor', 'w');
styleAxes(ax, ink, gridInk);
decorate(ax, win, oldFeat, winFill, refInk);
ylabel(ax, '||H|| / max ||H||');
title(ax, '(b) four-sensor norm: acceleration vs displacement', 'FontWeight', 'normal');
legend(ax, h2, {'acceleration', 'displacement'}, 'Location', 'best', 'Box', 'off');
% (c) normalised spatial vector magnitudes
ax = subplot(2, 2, 3); hold(ax, 'on');
V = Ha ./ normA;
for j = 1:size(SENSORS, 1)
    plot(ax, omega, abs(V(j, :)), '-', 'Color', colours{j}, 'LineWidth', 2, 'Marker', markers{j}, ...
        'MarkerIndices', 1:mEvery:nW, 'MarkerSize', 6, 'MarkerFaceColor', colours{j}, 'MarkerEdgeColor', 'w');
end
styleAxes(ax, ink, gridInk); ylim(ax, [0 1]);
decorate(ax, win, oldFeat, winFill, refInk);
xlabel(ax, '\omega (rad/s)'); ylabel(ax, '|H_j| / ||H||');
title(ax, '(c) spatial pattern (magnitude)', 'FontWeight', 'normal');
% (d) phase relative to the reference sensor
ax = subplot(2, 2, 4); hold(ax, 'on');
for j = setdiff(1:size(SENSORS, 1), refSensor)
    plot(ax, omega, unwrap(angle(Ha(j, :) ./ Ha(refSensor, :))) / pi, '-', 'Color', colours{j}, 'LineWidth', 2, ...
        'Marker', markers{j}, 'MarkerIndices', 1:mEvery:nW, 'MarkerSize', 6, 'MarkerFaceColor', colours{j}, ...
        'MarkerEdgeColor', 'w');
end
styleAxes(ax, ink, gridInk);
decorate(ax, win, oldFeat, winFill, refInk);
xlabel(ax, '\omega (rad/s)'); ylabel(ax, 'phase of H_j / H_{ref}  (\times\pi)');
title(ax, sprintf('(d) spatial pattern (phase), reference %s', SENSOR_NAMES{refSensor}), 'FontWeight', 'normal');
annotation(fig1, 'textbox', [0.01 0.005 0.98 0.03], 'String', ...
    'Shaded: proposed pole-fit windows. Dotted lines: old EMM features ("(zero)" = near-zero feature); overlay only, not expected to coincide.', ...
    'EdgeColor', 'none', 'Color', ink, 'FontSize', 9);
saveFigure(fig1, [baseName, '_frf.png']);

fig2 = figure('Color', 'w', 'Position', [80 80 1200 420]);
ax = subplot(1, 3, 1); hold(ax, 'on');
plot(ax, omega, abs(Crig(1, :)), '-', 'Color', colours{1}, 'LineWidth', 2);
styleAxes(ax, ink, gridInk); decorate(ax, win, oldFeat, winFill, refInk);
xlabel(ax, '\omega (rad/s)'); ylabel(ax, '|heave| (m s^{-2} per m incident)');
title(ax, '(a) rigid heave', 'FontWeight', 'normal');
ax = subplot(1, 3, 2); hold(ax, 'on');
plot(ax, omega, abs(Crig(2, :)), '-', 'Color', colours{2}, 'LineWidth', 2);
styleAxes(ax, ink, gridInk); decorate(ax, win, oldFeat, winFill, refInk);
xlabel(ax, '\omega (rad/s)'); ylabel(ax, '|pitch| (per unit non-dimensional x)');
title(ax, '(b) rigid pitch', 'FontWeight', 'normal');
ax = subplot(1, 3, 3); hold(ax, 'on');
plot(ax, omega, fracNR, '-', 'Color', colours{3}, 'LineWidth', 2);
styleAxes(ax, ink, gridInk); ylim(ax, [0 1]); decorate(ax, win, oldFeat, winFill, refInk);
xlabel(ax, '\omega (rad/s)'); ylabel(ax, '||E|| / ||H||');
title(ax, '(c) non-rigid fraction', 'FontWeight', 'normal');
saveFigure(fig2, [baseName, '_rigid.png']);

fprintf('\nDone in %.1f min.\n', toc(tStart) / 60);
catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ================================================================================================
%  Local functions
%% ================================================================================================
function ex = findLocalExtrema(y, kind, minPromDB)
% Interior local maxima ('max') or minima ('min') of a positive curve, with
% prominence in dB (20 log10). Toolbox-free (no findpeaks).
y = y(:).';
if strcmp(kind, 'min'), s = -1; else, s = 1; end
yd = s * 20 * log10(max(y, realmin));
n = numel(yd);
idx = []; prom = [];
for k = 2:n - 1
    if yd(k) > yd(k - 1) && yd(k) >= yd(k + 1)
        % prominence: drop to the lowest point before reaching a higher value (or the edge)
        L = k - 1; lowL = yd(k);
        while L >= 1 && yd(L) <= yd(k), lowL = min(lowL, yd(L)); L = L - 1; end
        Rr = k + 1; lowR = yd(k);
        while Rr <= n && yd(Rr) <= yd(k), lowR = min(lowR, yd(Rr)); Rr = Rr + 1; end
        pr = yd(k) - max(lowL, lowR);
        if pr >= minPromDB
            idx(end + 1) = k; prom(end + 1) = pr; %#ok<AGROW>
        end
    end
end
ex = struct('idx', idx, 'promDB', prom);
end

function d = describeCurve(omega, y, minPromDB)
% Peaks (with effective half-power width) and dips of a positive curve.
pk = findLocalExtrema(y, 'max', minPromDB);
dp = findLocalExtrema(y, 'min', minPromDB);
width = NaN(size(pk.idx));
for q = 1:numel(pk.idx)
    width(q) = halfPowerWidth(omega, y, pk.idx(q));
end
d.peaks = struct('omega', omega(pk.idx), 'value', y(pk.idx), 'promDB', pk.promDB, 'width', width, ...
    'widthRatio', width ./ (2 * omega(pk.idx)));
d.dips = struct('omega', omega(dp.idx), 'value', y(dp.idx), 'promDB', dp.promDB);
end

function w = halfPowerWidth(omega, y, k)
% Width between the -3 dB crossings either side of peak k, by linear
% interpolation. NaN if a side reaches a higher value or the band edge first.
lvl = y(k) / sqrt(2);
wL = NaN; wR = NaN;
for m = k - 1:-1:1
    if y(m) > y(k), break; end
    if y(m) <= lvl
        wL = omega(m) + (lvl - y(m)) * (omega(m + 1) - omega(m)) / (y(m + 1) - y(m));
        break;
    end
end
for m = k + 1:numel(y)
    if y(m) > y(k), break; end
    if y(m) <= lvl
        wR = omega(m - 1) + (y(m - 1) - lvl) * (omega(m) - omega(m - 1)) / (y(m - 1) - y(m));
        break;
    end
end
w = wR - wL;
end

function printCurve(name, d, oldFeat)
fprintf('%s\n', name);
if isempty(d.peaks.omega)
    fprintf('  no peaks above threshold\n');
end
for q = 1:numel(d.peaks.omega)
    fprintf('  peak %7.3f rad/s  prominence %5.1f dB  width %7.3f (ratio %.3f)%s\n', d.peaks.omega(q), ...
        d.peaks.promDB(q), d.peaks.width(q), d.peaks.widthRatio(q), nearestOld(d.peaks.omega(q), oldFeat));
end
for q = 1:numel(d.dips.omega)
    fprintf('  dip  %7.3f rad/s  prominence %5.1f dB%s\n', d.dips.omega(q), d.dips.promDB(q), ...
        nearestOld(d.dips.omega(q), oldFeat));
end
end

function s = nearestOld(w, oldFeat)
s = '';
if isempty(oldFeat.omega), return; end
[dmin, k] = min(abs(oldFeat.omega - w));
s = sprintf('   (nearest old feature %s at %.3f, %+.3f)', oldFeat.labels{k}, oldFeat.omega(k), w - oldFeat.omega(k));
if dmin > 0.5, s = ''; end
end

function [out, members] = mergeIntervals(win)
% Merge overlapping [lo hi] rows; members{w} lists the original rows merged.
out = zeros(0, 2); members = {};
if isempty(win), return; end
[~, ord] = sort(win(:, 1));
win = win(ord, :);
cur = win(1, :); mem = ord(1);
for q = 2:size(win, 1)
    if win(q, 1) <= cur(2)
        cur(2) = max(cur(2), win(q, 2)); mem(end + 1) = ord(q); %#ok<AGROW>
    else
        out(end + 1, :) = cur; members{end + 1} = mem; %#ok<AGROW>
        cur = win(q, :); mem = ord(q);
    end
end
out(end + 1, :) = cur; members{end + 1} = mem;
end

function decorate(ax, win, oldFeat, winFill, refInk)
% Call AFTER the data are plotted and the limits set: shading uses the
% current y-limits and is sent to the back so it never changes the limits.
yl = ylim(ax); xl = xlim(ax);
hf = [];
for w = 1:size(win, 1)
    hf(end + 1) = fill(ax, [win(w, 1) win(w, 2) win(w, 2) win(w, 1)], [yl(1) yl(1) yl(2) yl(2)], winFill, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off'); %#ok<AGROW>
end
if ~isempty(hf)                                         % send shading to the back (last child = bottom)
    ch = get(ax, 'Children');
    set(ax, 'Children', [ch(~ismember(ch, hf)); ch(ismember(ch, hf))]);
end
for k = 1:numel(oldFeat.omega)
    lab = oldFeat.labels{k};
    if ~oldFeat.isMax(k), lab = [lab, ' (zero)']; end
    plot(ax, oldFeat.omega(k) * [1 1], yl, ':', 'Color', refInk, 'LineWidth', 1, 'HandleVisibility', 'off');
    text(ax, oldFeat.omega(k), yl(1), [' ', lab], 'Rotation', 90, 'FontSize', 8, 'Color', refInk, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end
ylim(ax, yl); xlim(ax, xl);
end

function styleAxes(ax, ink, gridInk)
set(ax, 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off', 'Layer', 'top');
grid(ax, 'on'); set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
axis(ax, 'tight');
yl = ylim(ax);
if strcmp(get(ax, 'YScale'), 'log'), ylim(ax, yl .* [0.8 1.25]); end
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);     % R2020a+
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end