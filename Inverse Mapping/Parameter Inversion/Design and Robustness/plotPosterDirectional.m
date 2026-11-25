function plotPosterDirectional(matFile)
%PLOTPOSTERDIRECTIONAL  Poster figure 'Direction matters': how little second-direction energy is too much.
%   plotPosterDirectional()                    newest Results/V2/V2_20*.mat
%   plotPosterDirectional('...\V2_<stamp>.mat') that run
%
%   From the saved V2 results (no records, no EMM, no solver):
%     left   schematic, top-down: the floe, the primary wave (from the left) and a secondary wave at an
%            example angle phi = 60 deg
%     right  for each evaluated angle phi, the secondary-wave energy fraction f_2 (%) at which
%              orange diamond  the stiffness bias becomes significant (d_sys = 0.5)
%              blue circle     the coherence drop becomes detectable
%            joined by a bar. Blue left of orange: coherence warns before the bias matters. An angle whose
%            coherence never crossed the threshold within the sweep is marked 'not detected'.
%     note   one line under the schematic stating whether coherence always warns first.
%
%   Printed at exactly 130 x 65 mm (fixed paper size, no cropping), so the font sizes are the printed sizes.
%   Lives next to plotDirectionalForcingStudy.m.

if nargin < 1 || isempty(matFile)
    P = setupStudyPaths();
    Lst = dir(fullfile(P.results, 'V2', 'V2_20*.mat'));
    if isempty(Lst), error('plotPosterDirectional:none', 'No V2 results in %s', fullfile(P.results, 'V2')); end
    [~, k] = max([Lst.datenum]);
    matFile = fullfile(Lst(k).folder, Lst(k).name);
end
M = load(matFile);
results = M.results; cfg = M.cfg;
fprintf('plotPosterDirectional: %s\n', matFile);

X  = results.crossings;
TH = cfg.V2.thresholds;
ang = results.angleDeg(:).'; nA = numel(ang);
fBias = NaN(1, nA); fDet = NaN(1, nA); detStatus = repmat({'missing'}, 1, nA);
for i = 1:nA
    onPhi = abs([X.phiDeg] - ang(i)) < 1e-6;
    xb = X(onPhi & strcmp({X.criterion}, 'd_sys') & abs([X.threshold] - TH.dSys(1)) < 1e-10);
    if ~isempty(xb) && xb(1).crossed, fBias(i) = xb(1).f2; end
    xd = X(onPhi & strcmp({X.criterion}, 'detection'));
    if ~isempty(xd)
        if xd(1).crossed, fDet(i) = xd(1).f2; detStatus{i} = 'detected';
        else, detStatus{i} = 'not detected'; end
    end
    fprintf('  phi %3g deg: significant bias at f2 = %.3g%%; coherence: %s %s\n', ang(i), 100 * fBias(i), ...
            detStatus{i}, ternary(isfinite(fDet(i)), sprintf('at f2 = %.3g%%', 100 * fDet(i)), ''));
end
warnsFirst = isfinite(fDet) & isfinite(fBias) & fDet < fBias;

%% ---- style (poster palette, printed sizes) --------------------------------------------------------------------
ink    = [ 30  42  50] / 255;
navy   = [ 11  37  69] / 255;
ocean  = [ 19  73 123] / 255;
accent = [235 104  52] / 255;
ice    = [234 243 249] / 255;
grey   = [0.60 0.64 0.68];
fontName = 'Times New Roman';
fsTick = 16; fsLabel = 17; fsLab = 16;
W = 13; H = 6.5;                                          % cm

fig = figure('Units', 'centimeters', 'Position', [2 2 W H], 'Color', 'w');
if exist('theme', 'file') || exist('theme', 'builtin'), theme(fig, 'light'); end

%% ---- left: schematic, one primary and one secondary wave ---------------------------------------------------------
axs = axes(fig, 'Units', 'normalized', 'Position', [0.00 0.16 0.33 0.84], 'Color', 'w');
hold(axs, 'on'); axis(axs, 'equal'); axis(axs, 'off');
tt = linspace(0, 2 * pi, 200);
fill(axs, cos(tt), sin(tt), ice, 'EdgeColor', ink, 'LineWidth', 1.4);
phiShow = 60;
arrowTo(axs, [-2.25 0], [-1.08 0], navy, 2.4);                           % primary, from the left
e2 = [cosd(180 - phiShow), sind(180 - phiShow)];
arrowTo(axs, 2.15 * e2, 1.08 * e2, accent, 2.0);                         % secondary, at phi
aa = linspace(180, 180 - phiShow, 30);
plot(axs, 1.45 * cosd(aa), 1.45 * sind(aa), '-', 'Color', grey, 'LineWidth', 1.1);
text(axs, 1.72 * cosd(180 - phiShow / 2), 1.72 * sind(180 - phiShow / 2), '\phi', 'FontName', fontName, ...
     'FontSize', fsLab + 2, 'Color', ink, 'HorizontalAlignment', 'center');
text(axs, -1.65, -0.22, 'primary', 'FontName', fontName, 'FontSize', fsLab, 'Color', navy, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(axs, -0.45, 2.10, 'secondary, f_2', 'FontName', fontName, 'FontSize', fsLab, 'Color', accent, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
xlim(axs, [-2.45 1.15]); ylim(axs, [-1.2 2.4]);
if all(warnsFirst(isfinite(fBias)))
    verdict = 'coherence warns first';
else
    verdict = 'coherence does not always warn first';
end
annotation(fig, 'textbox', [0.00 0.00 0.36 0.16], 'String', verdict, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontName', fontName, ...
    'FontSize', fsLab, 'FontAngle', 'italic', 'Color', ink);

%% ---- right: first-crossing energy fractions per angle -------------------------------------------------------------
ax = axes(fig, 'Units', 'normalized', 'Position', [0.40 0.22 0.57 0.60], 'Color', 'w');
hold(ax, 'on');
y = 1:nA;
for i = 1:nA
    if isfinite(fBias(i)) && isfinite(fDet(i))
        plot(ax, 100 * [fDet(i) fBias(i)], [y(i) y(i)], '-', 'Color', grey, 'LineWidth', 2.4);
    end
end
hB = plot(ax, 100 * fBias, y, 'd', 'MarkerSize', 11, 'MarkerFaceColor', accent, 'MarkerEdgeColor', 'w', ...
          'LineWidth', 1);
hD = plot(ax, 100 * fDet, y, 'o', 'MarkerSize', 10, 'MarkerFaceColor', ocean, 'MarkerEdgeColor', 'w', ...
          'LineWidth', 1);
set(ax, 'XScale', 'log', 'YDir', 'reverse');
allF = 100 * [fBias, fDet]; allF = allF(isfinite(allF) & allF > 0);
xl = [min(allF) / 2, max(allF) * 4];                       % room on the right for 'not detected'
xlim(ax, xl); ylim(ax, [0.4, nA + 0.6]);
for i = 1:nA                                               % explain an absent detection marker
    if ~strcmp(detStatus{i}, 'detected')
        txt = 'not detected'; if strcmp(detStatus{i}, 'missing'), txt = 'n/a'; end
        x0 = 100 * fBias(i); if ~isfinite(x0), x0 = xl(1); end
        text(ax, x0 * 1.35, y(i), txt, 'FontName', fontName, 'FontSize', fsTick - 2, 'Color', ocean, ...
             'FontAngle', 'italic', 'VerticalAlignment', 'middle');
    end
end
cand = [0.001 0.003 0.01 0.03 0.1 0.3 1 3];                % three or four readable ticks
tk = cand(cand >= xl(1) & cand <= xl(2));
if numel(tk) > 4, tk = tk(1:2:end); end
set(ax, 'XTick', tk, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), tk, 'UniformOutput', false));
set(ax, 'YTick', y, 'YTickLabel', arrayfun(@(a) sprintf('%g\\circ', a), ang, 'UniformOutput', false), ...
        'TickLabelInterpreter', 'tex', 'XMinorTick', 'off', 'YMinorTick', 'off', 'XGrid', 'off', ...
        'YGrid', 'off', 'XMinorGrid', 'off');
set(ax, 'FontName', fontName, 'FontSize', fsTick, 'LineWidth', 0.8, 'Box', 'on', 'TickDir', 'in', ...
        'XColor', ink, 'YColor', ink, 'Layer', 'top');
xlabel(ax, 'secondary-wave energy f_2 (%)', 'FontName', fontName, 'FontSize', fsLabel);
ylabel(ax, '\phi', 'FontName', fontName, 'FontSize', fsLabel + 1);
legend(ax, [hB hD], {'bias significant', 'coherence detects'}, 'Location', 'northoutside', ...
       'Orientation', 'horizontal', 'Box', 'off', 'FontName', fontName, 'FontSize', fsLab, ...
       'TextColor', ink, 'Color', 'none');

% fixed printed size, no cropping
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [W H], 'PaperPosition', [0 0 W H]);
print(fig, 'directional.pdf', '-dpdf', '-vector');
fprintf('saved directional.pdf (%g x %g cm)\n', W, H);
end

%% ================================================================================================
function arrowTo(a, p0, p1, col, lw)
% straight arrow from p0 to p1 with a filled head at p1
u = (p1 - p0) / norm(p1 - p0); n = [-u(2), u(1)];
hl = 0.28; hw = 0.15;
base = p1 - hl * u;
plot(a, [p0(1) base(1)], [p0(2) base(2)], '-', 'Color', col, 'LineWidth', lw);
P = [p1; base + hw * n; base - hw * n];
patch(a, P(:, 1), P(:, 2), col, 'EdgeColor', col);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end