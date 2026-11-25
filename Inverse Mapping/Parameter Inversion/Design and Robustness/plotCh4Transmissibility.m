function files = plotCh4Transmissibility(results)
%PLOTCH4TRANSMISSIBILITY  Chapter 4 figure: parameter sensitivity of the EMM complex transmissibility.
%   plotCh4Transmissibility()          compute (or load the cached .mat) and draw
%   plotCh4Transmissibility(results)   redraw from a results struct already in memory
%
%   Two panels side by side, shared legend below:
%     (a) relative amplitude |T_j(omega)|
%     (b) relative phase arg T_j(omega) (unwrapped, degrees)
%   for one non-reference sensor of the reference layout, at
%     reference      (beta0, R0)       blue
%     stiffness      (1.2 beta0, R0)   orange
%     radius         (beta0, 1.2 R0)   green
%   EMM point transmissibilities T_j = H_aj / H_ar, no Welch smoothing, over the retained analysis band.
%
%   Choices fixed before looking at the curves (state them in the caption):
%     * plotted sensor: s14, the non-reference sensor of the reference layout farthest from the
%       reference s26 (s14 at (0.5R, 180 deg), s26 at (0.9R, 0 deg), separation 1.4R);
%     * perturbation: +20% in beta and in R, one at a time, gamma fixed;
%     * sensor positions held fixed in PHYSICAL coordinates when R changes, as in the inverse Jacobian.
%
%   The EMM is called through computeSensorFRF (Inverse Mapping/), wrapped in sensorAccelerationFRF()
%   at the bottom of this file. Three parameter sets x 102 frequencies = 306 EMM solves; the result is
%   cached, so later redraws are instant.
%
% Lives with the other plotting scripts; output goes to Results/ch4Transmissibility_<stamp>.

P = setupStudyPaths();                               % puts Reporting/ (sharcStyle etc.) on the path
S = sharcStyle();
files = {};

%% ---- configuration ---------------------------------------------------------------------------------------
cfg.H      = 1.88;                                   % water depth (m)
cfg.p0     = [4.6985e-5, 1.4548e-3, 0.3830];         % [beta0, gamma0, R0] (nondimensional)
cfg.scale  = 1.2;                                    % one-at-a-time perturbation factor
cfg.band   = [3.191, 8.283];                         % retained analysis band (rad/s)
cfg.dOmega = 0.05;                                   % plotting grid (rad/s); EMM point values, no interpolation
% reference layout, (r/R0, theta in degrees), theta = 0 on the +x axis (waves travel in -x)
cfg.sensorNames = {'s2', 's6', 's14', 's26'};
cfg.sensorPolar = [0.3 0; 0.3 180; 0.5 180; 0.9 0];
cfg.refName  = 's26';
cfg.plotName = 's14';

%% ---- compute or load -------------------------------------------------------------------------------------
if nargin < 1 || isempty(results)
    outDir = fullfile(P.results, 'ch4Transmissibility');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    cacheFile = fullfile(outDir, 'ch4Transmissibility_cache.mat');
    useCache = false;
    if exist(cacheFile, 'file')
        M0 = load(cacheFile);
        useCache = isfield(M0, 'cfg') && isequal(M0.cfg, cfg);       % stale if any setting changed
        if ~useCache, fprintf('  cache settings differ from cfg: recomputing\n'); end
    end
    if useCache
        results = M0.results;
        fprintf('  loaded cached transmissibilities from %s\n', cacheFile);
    else
        results = computeTransmissibilities(cfg);
        save(cacheFile, 'results', 'cfg');
        fprintf('  saved %s\n', cacheFile);
    end
else
    outDir = fullfile(P.results, 'ch4Transmissibility');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
end
run = struct('id', 'CH4', 'dir', outDir);

%% ---- figure ------------------------------------------------------------------------------------------------
w = results.omega;
[fig, ax] = sharcFigure('full', [1 2], 5.6);
hh = gobjects(1, 3);
for c = 1:3
    T = results.T(c, :);
    hh(c) = sharcSeries(ax(1), w, abs(T), c, 'Marker', 'none');
    sharcSeries(ax(2), w, rad2deg(unwrap(angle(T))), c, 'Marker', 'none', 'HandleVisibility', 'off');
end
for a = 1:2, xlim(ax(a), cfg.band); xlabel(ax(a), '\omega (rad/s)'); end
ylabel(ax(1), sprintf('|T_{%s}|', results.plotName));
ylabel(ax(2), sprintf('arg T_{%s} (deg)', results.plotName));
sharcPanel(ax(1), 'a', 'Relative amplitude');
sharcPanel(ax(2), 'b', 'Relative phase');
lab = {'Reference (\beta_0, R_0)', ...
       sprintf('Stiffer (%.1f\\beta_0, R_0)', cfg.scale), ...
       sprintf('Larger (\\beta_0, %.1fR_0)', cfg.scale)};
sharcSharedLegend(fig, ax, hh, lab, 3);
files = [files, saveSharcFigure(fig, run, 1, 'transmissibility_emm')];

%% ---- export for the Python thesis figure (plot_ch4_transmissibility.py) ----------------------------------
% columns: omega (rad/s), then Re and Im of T for the reference, stiffer and larger floes
csvFile = fullfile(outDir, 'ch4_transmissibility_emm.csv');
Tcsv = [w(:), real(results.T(1, :)).', imag(results.T(1, :)).', real(results.T(2, :)).', ...
    imag(results.T(2, :)).', real(results.T(3, :)).', imag(results.T(3, :)).'];
fid = fopen(csvFile, 'w');
fprintf(fid, 'omega,ReT_ref,ImT_ref,ReT_beta,ImT_beta,ReT_R,ImT_R\n');
fclose(fid);
writematrix(Tcsv, csvFile, 'WriteMode', 'append');
fprintf('  wrote %s\n', csvFile);
files{end + 1} = csvFile;

%% ---- printed check ---------------------------------------------------------------------------------------
dMag = max(abs(abs(results.T(2:3, :)) - abs(results.T(1, :))), [], 2);
dPh  = max(abs(rad2deg(unwrap(angle(results.T(2:3, :)), [], 2) - unwrap(angle(results.T(1, :))))), [], 2);
fprintf('\n  Ch4 transmissibility, sensor %s / ref %s, %d frequencies\n', results.plotName, ...
    results.refName, numel(w));
fprintf('    max |d|T||: stiffer %.3f, larger %.3f\n', dMag(1), dMag(2));
fprintf('    max |d arg T|: stiffer %.1f deg, larger %.1f deg\n', dPh(1), dPh(2));
end

%% ================================================================================================
function results = computeTransmissibilities(cfg)
% EMM point transmissibility of the plotted sensor for the three parameter sets.
omega = linspace(cfg.band(1), cfg.band(2), round(diff(cfg.band) / cfg.dOmega) + 1);   % both band edges included
R0phys = cfg.p0(3) * cfg.H;
% physical sensor coordinates from the REFERENCE radius, then held fixed for every parameter set
xy = [cfg.sensorPolar(:, 1) * R0phys .* cosd(cfg.sensorPolar(:, 2)), ...
      cfg.sensorPolar(:, 1) * R0phys .* sind(cfg.sensorPolar(:, 2))];
iRef  = find(strcmp(cfg.sensorNames, cfg.refName));
iPlot = find(strcmp(cfg.sensorNames, cfg.plotName));
pSet = [cfg.p0; cfg.p0 .* [cfg.scale 1 1]; cfg.p0 .* [1 1 cfg.scale]];
T = complex(zeros(3, numel(omega)));
for c = 1:3
    fprintf('  EMM set %d of 3: beta = %.4g, R = %.4g\n', c, pSet(c, 1), pSet(c, 3));
    Ha = sensorAccelerationFRF(pSet(c, :), xy([iPlot iRef], :), omega, cfg.H);   % 2 x nOmega
    T(c, :) = Ha(1, :) ./ Ha(2, :);
end
results = struct('omega', omega, 'T', T, 'pSet', pSet, 'plotName', cfg.plotName, ...
    'refName', cfg.refName, 'xy', xy);
end

function Ha = sensorAccelerationFRF(p, xy, omega, H)
% Complex vertical-acceleration FRF of the EMM disc at fixed PHYSICAL sensor positions.
%   p     = [beta, gamma, R], nondimensional
%   xy    = n_s x 2 physical coordinates (m), origin at the disc centre, theta = 0 on the incidence axis
%   omega = 1 x n_omega physical angular frequencies (rad/s)
%   H     = water depth (m)
% computeSensorFRF takes [r/H, theta(rad)] rows, so the physical positions are converted here and do not
% move when R changes. The common amplitude factor alpha(omega) of Eq. (3.10) cancels in T = Ha_j/Ha_r,
% so the plotted transmissibility does not depend on whether this version of computeSensorFRF applies it.
sensors = [hypot(xy(:, 1), xy(:, 2)) / H, atan2(xy(:, 2), xy(:, 1))];
Ha = computeSensorFRF(omega, sensors, p, 'Output', 'acceleration', 'Truncation', [50 10 10], ...
    'WaterDepth', H, 'UseParallel', false);
end