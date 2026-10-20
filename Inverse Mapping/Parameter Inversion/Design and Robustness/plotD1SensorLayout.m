function files = plotD1SensorLayout(results, cfg, run)
%PLOTD1SENSORLAYOUT  Draw and save the D1 figures from saved results (no EMM, no Fisher recomputation).
%   plotD1SensorLayout()                     redraw from the newest Results/D1/D1_*.mat
%   plotD1SensorLayout('...\D1_<stamp>.mat') redraw from that run
%   plotD1SensorLayout(results, cfg, run)    called by runD1SensorLayout
%
%   Figure 1  metrics of the A-optimal design against n_s, one line per noise case (2 x 2):
%             (a) sigma_lnbeta with the step gains of the headline case, (b) sigma_lnR, (c) d_W, (d) kappa(F)
%   Figure 2  the A-optimal layouts (headline noise) drawn top-down on the floe, n_s = 2, 3, 4
%   Figure 3  how forgiving placement is: ECDF of sigma_lnbeta over all admissible layouts (best
%             reference per layout), one line per n_s, optima marked
%   Figure 4  n_s = 4 designs on the floe: Fisher A-optimal, Level 2C, Gramian lambda_min optimum,
%             each with its Fisher rank and efficiency
%   Figure 5  the 33-point candidate grid, every point labelled, Level 2C and its reference marked
%   Every floe panel shows the incident wave direction (arrow on the incidence axis). Under the EMM
%   convention (incident potential ~ e^{k0 x}, k0 = i k, time factor e^{i omega t}) the wave arrives
%   from the theta = 0 side and travels towards -x; set WAVE_FROM_DEG below if that ever changes.
%   It also prints two checks that use only the saved search: (i) mirror symmetry about the
%   incidence axis (a layout and its mirror image must have the same A), and (ii) direction
%   robustness: the efficiency of each optimum when the wave arrives from another direction
%   (equivalently, the layout rotated by 45, 90, ... degrees on the grid).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'D1', 'D1_*.mat'));
        if isempty(Lst), error('plotD1SensorLayout:none', 'No D1 results in %s', fullfile(P.results, 'D1')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M = load(matFile);
    results = M.results; cfg = M.cfg; %#ok<NASGU>
    run = struct('id', 'D1', 'dir', fileparts(matFile));
    fprintf('  redrawing D1 figures from %s\n', matFile);
end
S = sharcStyle();
G = results.geometry;
ns = results.ns; nK = numel(ns);
noise = results.noise; nN = numel(noise);
files = {};
bestA = @(c, k) results.best(c, k, 1).A;

%% ---- checks from the saved search ---------------------------------------------------------------------------
fprintf('\n  D1 checks from the saved search (headline %s):\n', results.headline);
for k = 1:nK
    Sx = results.search{1, k};
    L = double(Sx.layouts);
    rowM = rowsOf(L, mapIdx(L, @(a) mod(-a, 8)));
    rel = max(abs(Sx.A(rowM) ./ Sx.A - 1));
    fprintf('    n_s = %d: mirror symmetry about the incidence axis, max |A(mirror)/A - 1| = %.1e\n', ns(k), rel);
end
fprintf('    direction robustness of the A-optima (efficiency A_opt / A of the rotated layout):\n');
rotEta = NaN(nK, 8);
for k = 1:nK
    Sx = results.search{1, k};
    L = double(Sx.layouts);
    b = bestA(1, k);
    for q = 0:7
        Lr = sort(mapIdx(b.layout, @(a) mod(a + q, 8)), 2);
        r = rowsOf(L, Lr);
        rotEta(k, q + 1) = b.A / Sx.A(r);
    end
    fprintf('      n_s = %d %-18s rotated by 0:45:315 deg: %s\n', ns(k), nameOf(G, b.layout), ...
        sprintf('%.2f ', rotEta(k, :)));
end

%% ---- figure 1: metrics against n_s ----------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [2 2], [], 'vgap', 1.6);
hh = []; lab = {};
for c = 1:nN
    v = zeros(4, nK);
    for k = 1:nK
        b = bestA(c, k);
        v(:, k) = [100 * b.sB; 100 * b.sR; b.dW; b.kappa];
    end
    hh(end + 1) = sharcSensorSeries(ax(1), ns, v(1, :), noise{c}); %#ok<AGROW>
    for j = 2:4, sharcSensorSeries(ax(j), ns, v(j, :), noise{c}); end
    lab{end + 1} = S.sensorLabel(noise{c}); %#ok<AGROW>
end
yl = {'\sigma_{ln\beta} (%)', '\sigma_{lnR} (%)', 'Welch bias significance d_W', 'condition number \kappa(F)'};
tt = {'Stiffness uncertainty', 'Radius uncertainty', 'Welch bias', 'Conditioning'};
for j = 1:4
    xlim(ax(j), [ns(1) - 0.3, ns(end) + 0.3]);
    sharcTicks(ax(j), 'x', ns);
    xlabel(ax(j), 'number of sensors n_s');
    ylabel(ax(j), yl{j});
    sharcPanel(ax(j), char('a' + j - 1), tt{j});
end
yb = [results.best(:, :, 1)]; yb = arrayfun(@(x) 100 * x.A.sB, yb);
yTop = 1.25 * max(yb(:));                                                      % zoomed: the 1% bound lies far above
ylim(ax(1), [0, yTop]);
if yTop > 1, sharcThreshold(ax(1), 'y', 1, '', 'tick'); end                    % Identifiable: sigma <= 1%
g = results.count(1).gain;
sh = arrayfun(@(k) 100 * bestA(1, k).sB, 1:nK);
for k = 2:nK                                                                   % step gains, headline, below the line
    text(ax(1), mean(ns(k - 1:k)), mean(sh(k - 1:k)) - 0.06 * yTop, sprintf('-%.1f%%', 100 * g(k)), ...
        'FontName', S.font, 'FontSize', S.smallFontSize, 'Color', S.axis, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end
dws = arrayfun(@(x) x.A.dW, results.best(:, :, 1));
ylim(ax(3), [0, max(1.1, 1.15 * max(dws(:)))]);
sharcThreshold(ax(3), 'y', 0.5, '', 'tick'); sharcThreshold(ax(3), 'y', 1, '', 'tick');
srs = arrayfun(@(x) 100 * x.A.sR, results.best(:, :, 1)); ylim(ax(2), [0, 1.15 * max(srs(:))]);
kps = arrayfun(@(x) x.A.kappa, results.best(:, :, 1)); ylim(ax(4), [0, 1.15 * max(kps(:))]);
lab = strrep(lab, 'G3 assumption', 'design-stage stress level');             % no internal test codes in figures
sharcSharedLegend(fig, ax, hh, lab, nN);
files = [files, saveSharcFigure(fig, run, 1, 'sensor_count')];

%% ---- figure 2: A-optimal layouts on the floe --------------------------------------------------------------------
[fig, ax] = floeFigure(nK);
for k = 1:nK
    b = bestA(1, k);
    hl = plotFloeLayout(ax(k), G, b.layout, b.ref, S);
    sharcPanel(ax(k), char('a' + k - 1), sprintf('n_s = %d: \\sigma_{ln\\beta} = %.3f%%', ns(k), 100 * b.sB));
end
sharcSharedLegend(fig, ax, hl, {'candidate point', 'sensor', 'reference sensor', 'wave incidence axis'}, 4);
files = [files, saveSharcFigure(fig, run, 2, 'optimal_layouts')];

%% ---- figure 3: ECDF of sigma_lnbeta over all admissible layouts --------------------------------------------------
[fig, ax] = sharcFigure('single', [1 1], 7.0);
hh = []; lab = {};
xs = [];
for k = 1:nK
    Sx = results.search{1, k};
    ok = Sx.admissible(:, 1) & isfinite(Sx.sB);
    x = sort(100 * Sx.sB(ok)); y = (1:numel(x)).' / numel(x);
    xs = [xs; x]; %#ok<AGROW>
    hh(end + 1) = sharcSeries(ax, x, y, k, 'Marker', 'none'); %#ok<AGROW>
    sharcSeries(ax, x(1), y(1), k, 'LineStyle', 'none', 'HandleVisibility', 'off');   % the optimum
    lab{end + 1} = sprintf('n_s = %d (%d layouts)', ns(k), numel(x)); %#ok<AGROW>
end
set(ax, 'XScale', 'log');
xlim(ax, [10^floor(log10(min(xs))), 10^ceil(log10(prctileNR(xs, 99.5)))]);
ylim(ax, [0 1]);
sharcTicks(ax, 'x', [0.2 0.5 1 2 5 10 20 50]);
sharcTicks(ax, 'y', [0 0.25 0.5 0.75 1]);
xlabel(ax, '\sigma_{ln\beta} of the layout (%), best reference');
ylabel(ax, 'fraction of admissible layouts');
sharcThreshold(ax, 'x', 1, '1%', 'top');
sharcPanel(ax, 'a', sprintf('Placement sensitivity (%s)', S.sensorLabel(results.headline)));
sharcLegend(ax, hh, lab, 'southeast');
files = [files, saveSharcFigure(fig, run, 3, 'placement_ecdf')];

%% ---- figure 4: n_s = 4 designs compared ---------------------------------------------------------------------------
k4 = find(ns == numel(G.level2C));
if ~isempty(k4)
    designs = {};
    b = bestA(1, k4);
    nAdm = nnz(results.search{1, k4}.admissible(:, 1));
    designs(end + 1, :) = {'Fisher A-optimal', b.layout, b.ref, sprintf('rank 1 / %d, \\eta = 1', nAdm)};
    if isfield(results, 'level2C') && isfield(results.level2C, 'row') && ~isempty(results.level2C(1).row)
        l2 = results.level2C(1);
        designs(end + 1, :) = {'Level 2C (best reference)', G.level2C, l2.bestRef, sprintf('rank %d / %d, \\eta = %.2f', ...
            l2.rank_bestRef, l2.nAdmissible, l2.eff_bestRef)};
    end
    if isfield(results, 'gramian') && ~isempty(results.gramian)
        gm = results.gramian(k4); f = gm.fisher(1);
        Sx = results.search{1, k4};
        r = rowsOf(double(Sx.layouts), gm.tied(1, :));
        tie = ''; if size(gm.tied, 1) > 1, tie = sprintf(' (1 of %d tied)', size(gm.tied, 1)); end
        designs(end + 1, :) = {['Gramian \lambda_{min}' tie], gm.tied(1, :), double(Sx.refA(r)), ...
            sprintf('rank %d / %d, \\eta = %.2f', f.rank(1), f.nAdmissible, f.eta(1))};
    end
    nD = size(designs, 1);
    [fig, ax] = floeFigure(nD);
    for d = 1:nD
        hl = plotFloeLayout(ax(d), G, designs{d, 2}, designs{d, 3}, S);
        sharcPanel(ax(d), char('a' + d - 1), designs{d, 1});
        text(ax(d), 0, -1.22 * G.Rphys, designs{d, 4}, 'HorizontalAlignment', 'center', ...
            'FontName', S.font, 'FontSize', S.smallFontSize, 'Color', S.ink);
    end
    sharcSharedLegend(fig, ax, hl, {'candidate point', 'sensor', 'reference sensor', 'wave incidence axis'}, 4);
    files = [files, saveSharcFigure(fig, run, 4, 'design_comparison')];
end

%% ---- figure 5: the candidate grid, labelled ------------------------------------------------------------------------
[fig, ax] = sharcFigure('single', [1 1], 9.2, 'left', 0.2, 'right', 0.2, 'bottom', 1.0, 'top', 0.55);
hl = plotCandidateGrid(ax, G, S);
title(ax, sprintf('Candidate sensor positions (%d points)', size(G.xy, 1)), 'FontName', S.font, ...
    'FontSize', S.fontSize, 'FontWeight', 'normal', 'Color', S.ink);
refName = G.names{double(G.level2C(end))};
if isfield(G, 'level2CRef') && ~isempty(G.level2CRef), refName = G.names{double(G.level2CRef)}; end
lg = sharcLegend(ax, hl, {'candidate point', 'Level 2C sensor', sprintf('Level 2C reference (%s)', refName)}, ...
    'southoutside');
set(lg, 'Orientation', 'vertical', 'NumColumns', 2);                            % two rows: fits the single-column width
files = [files, saveSharcFigure(fig, run, 5, 'candidate_grid')];
end

%% ================================================================================================
function [fig, ax] = floeFigure(n)
% n square floe panels across the full text width (axes switched off; equal aspect)
[fig, ax] = sharcFigure('full', [1 n], 6.6, 'left', 0.25, 'right', 0.25, 'bottom', 0.45, 'top', 0.6, 'hgap', 0.35);
end

function h = plotFloeLayout(ax, G, layout, ref, S)
% top-down floe: boundary, 0.9 R candidate limit, all candidates, the layout and its reference
R = G.Rphys; t = linspace(0, 2 * pi, 361);
hold(ax, 'on');
plot(ax, R * cos(t), R * sin(t), '-', 'Color', S.ink, 'LineWidth', S.lineWidth, 'HandleVisibility', 'off');
plot(ax, G.maxRadiusFrac * R * cos(t), G.maxRadiusFrac * R * sin(t), ':', 'Color', S.threshold, ...
    'LineWidth', S.thinLineWidth, 'HandleVisibility', 'off');
h(4) = plot(ax, [-1.08 1.08] * R, [0 0], '--', 'Color', S.threshold, 'LineWidth', S.thinLineWidth);   % theta = 0 axis
h(1) = plot(ax, G.xy(:, 1), G.xy(:, 2), 'o', 'Color', S.threshold, 'MarkerSize', 3, 'MarkerFaceColor', 'w', ...
    'LineWidth', 0.6);
layout = double(layout(:).');
oth = layout(layout ~= ref);
h(2) = sharcSeries(ax, G.xy(oth, 1), G.xy(oth, 2), 1, 'LineStyle', 'none', 'MarkerSize', 8);
h(3) = sharcSeries(ax, G.xy(ref, 1), G.xy(ref, 2), 2, 'LineStyle', 'none', 'Marker', 'p', 'MarkerSize', 13);
onAxis = layout(abs(G.xy(layout, 2)) < 1e-9);         % points on the incidence axis, left to right
[~, o] = sort(G.xy(onAxis, 1)); onAxis = onAxis(o);
for i = layout
    a = atan2(G.xy(i, 2), G.xy(i, 1));
    rI = hypot(G.xy(i, 1), G.xy(i, 2)) / R;
    off = (0.13 + 0.09 * (rI > 0.8)) * R * [cos(a), sin(a)];   % off-axis: radially outwards, outer ring clear of the rim
    j = find(onAxis == i);
    if ~isempty(j)                                         % on the axis: alternate above / below
        off = [-0.07 * R * sign(G.xy(i, 1)) * (rI > 0.8), 0.13 * R * (2 * mod(j, 2) - 1)];
    end
    text(ax, G.xy(i, 1) + off(1), G.xy(i, 2) + off(2), G.names{i}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontName', S.font, 'FontSize', S.smallFontSize, 'Color', S.ink, ...
        'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 0.5);
end
drawIncidentWave(ax, R, S);
axis(ax, 'equal');
xlim(ax, [-1.15 1.62] * R); ylim(ax, [-1.3 1.32] * R);
axis(ax, 'off');
end

function h = plotCandidateGrid(ax, G, S)
% the full candidate grid: boundary, rings, spokes, every point labelled radially outwards
R = G.Rphys; t = linspace(0, 2 * pi, 361);
hold(ax, 'on');
plot(ax, R * cos(t), R * sin(t), '-', 'Color', S.ink, 'LineWidth', S.lineWidth, 'HandleVisibility', 'off');
rings = unique(round(1e6 * hypot(G.xy(:, 1), G.xy(:, 2)) / R) / 1e6);
rings = rings(rings > 0).';
for rr = rings                                                  % candidate rings
    plot(ax, rr * R * cos(t), rr * R * sin(t), ':', 'Color', S.threshold, 'LineWidth', S.thinLineWidth, ...
        'HandleVisibility', 'off');
end
for a = (0:7) * pi / 4                                          % spokes through the angular slots
    plot(ax, [0 max(rings)] * R * cos(a), [0 max(rings)] * R * sin(a), ':', 'Color', S.threshold, ...
        'LineWidth', S.thinLineWidth, 'HandleVisibility', 'off');
end
plot(ax, [-1.08 1.08] * R, [0 0], '--', 'Color', S.threshold, 'LineWidth', S.thinLineWidth, ...
    'HandleVisibility', 'off');                                  % incidence axis
h(1) = plot(ax, G.xy(:, 1), G.xy(:, 2), 'o', 'Color', S.axis, 'MarkerSize', 4, 'MarkerFaceColor', 'w', ...
    'LineWidth', 0.7);
l2 = double(G.level2C(:).');
ref = l2(end);
if isfield(G, 'level2CRef') && ~isempty(G.level2CRef), ref = double(G.level2CRef); end   % grid index
oth = l2(l2 ~= ref);
h(2) = sharcSeries(ax, G.xy(oth, 1), G.xy(oth, 2), 1, 'LineStyle', 'none', 'MarkerSize', 7);
h(3) = sharcSeries(ax, G.xy(ref, 1), G.xy(ref, 2), 2, 'LineStyle', 'none', 'Marker', 'p', 'MarkerSize', 11);
for i = 1:size(G.xy, 1)
    r = hypot(G.xy(i, 1), G.xy(i, 2));
    if r < 1e-9
        off = [0, 0.075 * R];                                    % centre: label above
    else
        a = atan2(G.xy(i, 2), G.xy(i, 1));
        off = 0.085 * R * [cos(a), sin(a)];                      % label radially outwards
        if abs(G.xy(i, 2)) < 1e-9, off = [0.0, 0.075 * R]; end   % on the incidence axis: label above
    end
    isL2 = any(l2 == i);
    text(ax, G.xy(i, 1) + off(1), G.xy(i, 2) + off(2), G.names{i}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontName', S.font, 'FontSize', S.smallFontSize - 1, ...
        'Color', S.ink, 'FontWeight', ternary(isL2, 'bold', 'normal'));
end
ringTxt = {'rings at r/R', strjoin(arrayfun(@(x) sprintf('%.1f', x), rings, 'UniformOutput', false), ', ')};
for q = 1:2                                                     % two single-line labels beside the floe
    text(ax, 1.36 * R, (-0.42 - 0.12 * (q - 1)) * R, ringTxt{q}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontName', S.font, 'FontSize', S.smallFontSize - 1, 'Color', S.axis);
end
drawIncidentWave(ax, R, S);
axis(ax, 'equal');
xlim(ax, [-1.12 1.62] * R); ylim(ax, [-1.12 1.12] * R);
axis(ax, 'off');
end

function drawIncidentWave(ax, R, S)
% arrow on the incidence axis, outside the floe, showing the direction of wave travel
WAVE_FROM_DEG = 0;                       % the wave arrives from this polar angle (EMM: theta = 0 side)
a = WAVE_FROM_DEG * pi / 180;
u = -[cos(a), sin(a)];                    % direction of travel (towards the floe centre)
p0 = 1.55 * R * [cos(a), sin(a)];         % tail
p1 = 1.17 * R * [cos(a), sin(a)];         % head tip
hw = 0.06 * R; hl = 0.11 * R;             % head half-width and length
n = [-u(2), u(1)];
base = p1 - hl * u;
plot(ax, [p0(1) base(1)], [p0(2) base(2)], '-', 'Color', S.ink, 'LineWidth', 1.4, 'HandleVisibility', 'off');
patch(ax, [p1(1), base(1) + hw * n(1), base(1) - hw * n(1)], [p1(2), base(2) + hw * n(2), base(2) - hw * n(2)], ...
    S.ink, 'EdgeColor', S.ink, 'HandleVisibility', 'off');
for c = [1.42 1.30]                       % two short crest lines behind the arrow
    q = c * R * [cos(a), sin(a)];
    plot(ax, q(1) + [-1 1] * 0.13 * R * n(1), q(2) + [-1 1] * 0.13 * R * n(2), '-', 'Color', S.threshold, ...
        'LineWidth', S.thinLineWidth, 'HandleVisibility', 'off');
end
m = 0.5 * (p0 + p1);
text(ax, m(1), m(2) + 0.2 * R, 'wave', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontName', S.font, 'FontSize', S.smallFontSize, 'Color', S.ink);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function Lm = mapIdx(L, f)
% apply an angular-slot map f (0..7 -> 0..7) to grid indices; the centre (index 1) is unchanged
Lm = L;
ring = L >= 2;
a = mod(L(ring) - 2, 8); rr = floor((L(ring) - 2) / 8);
Lm(ring) = 2 + 8 * rr + f(a);
end

function r = rowsOf(L, Lq)
% row of each query layout (any order within the row) in the layout table L (rows ascending)
Lq = sort(Lq, 2);
base = 34;
key = L * (base .^ (0:size(L, 2) - 1)).';
kq = Lq * (base .^ (0:size(Lq, 2) - 1)).';
[tf, r] = ismember(kq, key);
if ~all(tf), error('plotD1:row', 'Layout not found in the search table.'); end
end

function s = nameOf(G, L)
s = strjoin(G.names(double(L)), ' ');
end

function p = prctileNR(x, q)
x = sort(x(:)); p = x(max(1, min(numel(x), round(q / 100 * (numel(x) - 1)) + 1)));
end