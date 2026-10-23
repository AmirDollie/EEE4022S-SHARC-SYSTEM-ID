function plotImuCalibration(matFile, outDir)
% plotImuCalibration  Appendix figure: static-hold magnitude error before and after
% the accelerometer calibration of runAccelCalibration.m.
%   plotImuCalibration(matFile)          saves to the current folder
%   plotImuCalibration(matFile, outDir)  saves to outDir (e.g. 'figs/appD')
%
%   matFile is the accelCalibration_<stamp>.mat written by runAccelCalibration.m
%   (struct 'out'). Nothing is refitted: the figure is drawn from the saved fit.
%
%   Layout: two stacked panels sharing the hold axis, (a) raw and (b) calibrated,
%   each with its own y-range so the calibrated residual is readable. Each hold is a
%   category, so values are drawn as stems from zero (no connecting lines), with the
%   three sensors offset side by side within each hold.
%   Output: imu_calibration.pdf (vector) / .png (preview). Include at natural width.

if nargin < 2, outDir = pwd; end
if ~exist(outDir, 'dir'), mkdir(outDir); end
S = load(matFile, 'out'); out = S.out;

% ---------------------------------------------------------------- style
FS = 9; FSL = 8.5;
ink  = [0.10 0.10 0.10]; grey = [0.55 0.55 0.55]; gridC = [0.88 0.88 0.88];
cols = [0.13 0.40 0.70;      % LSM6DS3TR-C
        0.85 0.45 0.10;      % LSM6DSV16X
        0.30 0.60 0.30];     % ICM-42688-P   (match these to the sensor colours used elsewhere)
figW = 16; figH = 8.5;       % cm, final size
dx = 0.22;                   % side-by-side offset of the three sensors within a hold

% ---------------------------------------------------------------- data (mg)
g0 = 9.80665; G = out.cfg.G;
nH = numel(out.holdNames); nS = numel(out.sensors);
errRaw = NaN(nH, nS); errCal = NaN(nH, nS);
for k = 1:nS
    A = out.Amean(:, :, k); c = out.cal{k};
    errRaw(:, k) = 1e3 * (sqrt(sum(A.^2, 2)) - G) / g0;
    u = c.holdsUsed;
    errCal(u, k) = 1e3 * (sqrt(sum(((A(u, :) - c.offset) * c.M.').^2, 2)) - G) / g0;
end
labels = prettyHold(out.holdNames);

% ---------------------------------------------------------------- figure
fig = figure('Units', 'centimeters', 'Position', [2 2 figW figH], 'Color', 'w');
pos = [0.08 0.58 0.70 0.36; 0.08 0.15 0.70 0.36];
data = {errRaw, errCal};
ttl  = {'(a) Raw', '(b) Calibrated'};
hl = gobjects(1, nS);
for p = 1:2
    ax = axes(fig, 'Units', 'normalized', 'Position', pos(p, :)); hold(ax, 'on');
    E = data{p};
    lim = max(abs(E(:)), [], 'omitnan'); lim = niceLim(1.15 * lim);
    for h = 1:nH                                   % light vertical bands separate the holds
        if mod(h, 2) == 0
            patch(ax, h + [-0.5 0.5 0.5 -0.5], lim * [-1 -1 1 1], [0.965 0.965 0.965], ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    plot(ax, [0.5 nH + 0.5], [0 0], '-', 'Color', grey, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    for k = 1:nS
        x = (1:nH).' + (k - 2) * dx;
        for h = 1:nH
            if isfinite(E(h, k))
                plot(ax, [x(h) x(h)], [0 E(h, k)], '-', 'Color', cols(k, :), 'LineWidth', 1.1, ...
                    'HandleVisibility', 'off');
            end
        end
        hl(k) = plot(ax, x, E(:, k), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', cols(k, :), ...
            'MarkerEdgeColor', cols(k, :), 'LineStyle', 'none');
    end
    xlim(ax, [0.5 nH + 0.5]); ylim(ax, lim * [-1 1]);
    set(ax, 'FontSize', FS, 'Box', 'off', 'XColor', ink, 'YColor', ink, 'TickDir', 'out', ...
        'XTick', 1:nH, 'YGrid', 'on', 'GridColor', gridC, 'GridAlpha', 1, 'Layer', 'top');
    if p == 1, set(ax, 'XTickLabel', {}); else, set(ax, 'XTickLabel', labels); end
    ylabel(ax, '|a| - g (mg)', 'FontSize', FS);
    title(ax, ttl{p}, 'FontSize', FS, 'FontWeight', 'normal', 'HorizontalAlignment', 'left', ...
        'Units', 'normalized', 'Position', [0 1.02 0]);
end
lg = legend(hl, out.sensors, 'FontSize', FSL, 'Box', 'off', 'Interpreter', 'none');
set(lg, 'Units', 'normalized', 'Position', [0.80 0.62 0.18 0.20]);

name = fullfile(outDir, 'imu_calibration');
if exist('exportgraphics', 'file') || exist('exportgraphics', 'builtin')
    exportgraphics(fig, [name '.pdf'], 'ContentType', 'vector');
    exportgraphics(fig, [name '.png'], 'Resolution', 300);
else
    set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [figW figH], 'PaperPosition', [0 0 figW figH]);
    print(fig, [name '.pdf'], '-dpdf');
    print(fig, [name '.png'], '-dpng', '-r300');
end
fprintf('Saved %s (.pdf, .png)\n', name);
end

% =====================================================================
function L = prettyHold(names)
% 'up' -> '+z up', 'side090' -> 'side 90°', 'tilt120' -> 'tilt 120°'
L = names;
for i = 1:numel(names)
    n = names{i};
    t = regexp(n, '^(side|tilt)(\d+)$', 'tokens', 'once');
    if strcmpi(n, 'up'),       L{i} = '+z up';
    elseif strcmpi(n, 'down'), L{i} = '-z up';
    elseif ~isempty(t),        L{i} = sprintf('%s %d%s', t{1}, str2double(t{2}), char(176));
    end
end
end

function v = niceLim(x)
% round an axis limit up to 1, 2 or 5 times a power of ten
if ~isfinite(x) || x <= 0, v = 1; return; end
p = 10^floor(log10(x)); m = x / p;
if m <= 1, v = p; elseif m <= 2, v = 2 * p; elseif m <= 5, v = 5 * p; else, v = 10 * p; end
end