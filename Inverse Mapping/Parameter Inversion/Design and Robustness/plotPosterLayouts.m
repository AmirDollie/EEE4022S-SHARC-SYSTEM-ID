function plotPosterLayouts(d3File, d1File)
%PLOTPOSTERLAYOUTS  Poster figure 'Few sensors suffice': the optimal layouts on the floe, top-down.
%   plotPosterLayouts()                          newest Results/D3/D3_20*.mat and Results/D1/D1_*.mat
%   plotPosterLayouts(d3File, d1File)            those runs
%
%   (a) Known incident heading:   the D1 A-optimal layout, n_s = 2, with sigma_lnbeta AT the known heading
%       taken from the D1 results (results.best(1, k, 1).A.sB, headline noise case). No fallback: the script
%       stops if that value is missing rather than quoting a worst-case bound.
%   (b) Unknown incident heading: the D3a robust (worst-case) optimum, n_s = 4, with its worst-case
%       sigma_lnbeta over the headings (heading estimated). ONE incident wave is drawn, with a double-headed
%       arc for its unknown heading psi, so it cannot be mistaken for several simultaneous waves (that is V2).
%   Every sensor is labelled (outer-ring labels inwards); the reference sensor, if saved, is purple.
%
%   Drawn mirrored left-right (x -> -x), so the wave arrives from the LEFT, as in the poster's hero figure;
%   the problem is symmetric under this reflection, so the layouts and numbers are unchanged.
%   Printed at exactly 130 x 65 mm (fixed paper size, no cropping), so the font sizes are the printed sizes.
%   Lives next to plotDirectionRobustLayoutStudy.m (Inverse Mapping/Parameter Inversion/Design and Robustness/).

P = setupStudyPaths();
if nargin < 1 || isempty(d3File), d3File = newest(fullfile(P.results, 'D3'), 'D3_20*.mat'); end
if nargin < 2 || isempty(d1File), d1File = newest(fullfile(P.results, 'D1'), 'D1_*.mat'); end
M3 = load(d3File); M1 = load(d1File);
fprintf('plotPosterLayouts:\n  D3 %s\n  D1 %s\n', d3File, d1File);

%% ---- D3: the two layouts and the worst-case sigma of the robust one ---------------------------------------------
R3 = M3.results;
A  = R3.d3a(1);
g  = R3.geom;
ns = R3.nsList;
pick = @(k, name) A.perNs(k).sel(strcmp({A.perNs(k).sel.name}, name));
k2 = find(ns == 2, 1); k4 = find(ns == 4, 1);
if isempty(k2) || isempty(k4), error('plotPosterLayouts:ns', 'n_s = 2 and 4 must both be in the D3 nsList'); end
d1sel = pick(k2, 'D1');
rw    = pick(k4, 'robustWorst');

%% ---- D1: sigma at the known heading for n_s = 2 (headline noise case) --------------------------------------------
R1 = M1.results;
kD1 = find(R1.ns == 2, 1);
if isempty(kD1), error('plotPosterLayouts:d1', 'n_s = 2 not found in the D1 results'); end
bD1 = R1.best(1, kD1, 1);
if ~isfield(bD1, 'A') || ~isfield(bD1.A, 'sB') || ~isfinite(bD1.A.sB)
    error('plotPosterLayouts:d1sigma', 'D1 sigma_lnbeta for n_s = 2 not found (results.best(1,k,1).A.sB)');
end
sKnown = 100 * bD1.A.sB;
if ~isequal(sort(double(bD1.A.layout(:))), sort(double(d1sel.layout(:))))
    warning('plotPosterLayouts:mismatch', ['The D1 n_s = 2 optimum (%s) differs from the D1 layout stored in D3 ' ...
            '(%s); the figure draws the D3 one.'], strjoin(g.names(bD1.A.layout), ' '), strjoin(g.names(d1sel.layout), ' '));
end
sWorst = 100 * rw.sBworstN;
fprintf('  (a) known heading,   n_s = 2: %s   sigma_lnbeta = %.3f%%  (D1, heading known)\n', ...
        strjoin(g.names(d1sel.layout), ' '), sKnown);
fprintf('  (b) unknown heading, n_s = 4: %s   worst-case sigma_lnbeta = %.3f%%  (D3a, heading estimated)\n', ...
        strjoin(g.names(rw.layout), ' '), sWorst);

%% ---- style (poster palette, printed sizes) --------------------------------------------------------------------
ink   = [ 30  42  50] / 255;
navy  = [ 11  37  69] / 255;
ocean = [ 19  73 123] / 255;
refc  = [107  70 193] / 255;
ice   = [234 243 249] / 255;
grey  = [0.72 0.76 0.79];
fontName = 'Times New Roman';
fsTitle = 20; fsNote = 17; fsLab = 15;
W = 13; H = 6.5;                                          % cm

fig = figure('Units', 'centimeters', 'Position', [2 2 W H], 'Color', 'w');
if exist('theme', 'file') || exist('theme', 'builtin'), theme(fig, 'light'); end
Rp = g.Rphys; tt = linspace(0, 2 * pi, 240);
xy = g.xy; xy(:, 1) = -xy(:, 1);                         % mirror: wave from the left
panels = {d1sel, '(a) known heading',   sprintf('2 sensors:  \\sigma_{ln\\beta} = %.3f%%', sKnown), 'known';
          rw,    '(b) unknown heading', sprintf('4 sensors:  worst \\sigma_{ln\\beta} = %.3f%%', sWorst), 'unknown'};
for p = 1:2
    sel = panels{p, 1};
    ax = axes(fig, 'Units', 'normalized', 'Position', [0.005 + 0.50 * (p - 1), 0.15, 0.49, 0.68], 'Color', 'w');
    hold(ax, 'on');
    fill(ax, Rp * cos(tt), Rp * sin(tt), ice, 'EdgeColor', ink, 'LineWidth', 1.5);           % floe
    plot(ax, xy(:, 1), xy(:, 2), '.', 'MarkerSize', 5, 'Color', grey);                       % candidate grid
    drawArrow(ax, [-1.62 0] * Rp, [-1.06 0] * Rp, navy, 2.0, 0.13 * Rp);                     % one incident wave
    if strcmp(panels{p, 4}, 'known')
        plot(ax, [-1.0 1.0] * Rp, [0 0], '--', 'Color', grey, 'LineWidth', 0.8);            % incidence axis
    else
        % heading uncertainty: a double-headed arc through the wave direction, labelled psi
        aa = linspace(125, 235, 60);
        rA = 1.38 * Rp;
        plot(ax, rA * cosd(aa(3:end-2)), rA * sind(aa(3:end-2)), '-', 'Color', navy, 'LineWidth', 1.4);
        arcHead(ax, rA, aa(end-2), aa(end), navy, 0.11 * Rp);
        arcHead(ax, rA, aa(3), aa(1), navy, 0.11 * Rp);
        text(ax, 1.30 * Rp * cosd(140) - 0.05 * Rp, 1.30 * Rp * sind(140) + 0.12 * Rp, '\psi ?', ...
             'FontName', fontName, 'FontSize', fsLab + 1, 'Color', navy, 'HorizontalAlignment', 'right');
    end
    lay = double(sel.layout(:).');
    ref = [];
    if isfield(sel, 'ref') && ~isempty(sel.ref), ref = double(sel.ref); end
    if p == 1 && isempty(ref) && isfield(bD1.A, 'ref') && ~isempty(bD1.A.ref)
        ref = double(bD1.A.ref);                              % D1 stores the best reference of its optimum
    end
    for i = lay
        c = ocean; if any(ref == i), c = refc; end
        plot(ax, xy(i, 1), xy(i, 2), 'o', 'MarkerSize', 11, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', ...
             'LineWidth', 1.4);
        r = hypot(xy(i, 1), xy(i, 2));
        if r < 1e-9 || abs(xy(i, 2)) < 1e-9
            off = [0, 0.19 * Rp];                                 % centre or incidence axis: above
        elseif r > 0.70 * Rp
            off = -0.20 * Rp * [xy(i, 1), xy(i, 2)] / r;          % outer ring: inwards, clear of the rim
        else
            off = 0.20 * Rp * [xy(i, 1), xy(i, 2)] / r;           % inner rings: outwards
        end
        nm = g.names{i}; if any(ref == i), nm = [nm ' (ref)']; end
        text(ax, xy(i, 1) + off(1), xy(i, 2) + off(2), nm, 'FontName', fontName, 'FontSize', fsLab, ...
             'FontWeight', 'bold', 'Color', c, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    axis(ax, 'equal'); axis(ax, 'off');
    xlim(ax, [-1.66 1.12] * Rp); ylim(ax, [-1.12 1.12] * Rp);
    annotation(fig, 'textbox', [0.005 + 0.50 * (p - 1), 0.84, 0.49, 0.15], 'String', panels{p, 2}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', fontName, 'FontSize', fsTitle, 'FontWeight', 'bold', 'Color', navy, 'Interpreter', 'tex');
    annotation(fig, 'textbox', [0.005 + 0.50 * (p - 1), 0.00, 0.49, 0.15], 'String', panels{p, 3}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', fontName, 'FontSize', fsNote, 'Color', ink, 'Interpreter', 'tex');
end

% fixed printed size, no cropping
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [W H], 'PaperPosition', [0 0 W H]);
print(fig, 'layouts.pdf', '-dpdf', '-vector');
fprintf('saved layouts.pdf (%g x %g cm)\n', W, H);
end

%% ================================================================================================
function f = newest(folder, pattern)
L = dir(fullfile(folder, pattern));
if isempty(L), error('plotPosterLayouts:none', 'No %s in %s', pattern, folder); end
[~, k] = max([L.datenum]); f = fullfile(L(k).folder, L(k).name);
end

function drawArrow(a, p0, p1, col, lw, hl)
% straight arrow p0 -> p1 with a filled head of length hl
u = (p1 - p0) / norm(p1 - p0); n = [-u(2), u(1)];
base = p1 - hl * u;
plot(a, [p0(1) base(1)], [p0(2) base(2)], '-', 'Color', col, 'LineWidth', lw);
P = [p1; base + 0.55 * hl * n; base - 0.55 * hl * n];
patch(a, P(:, 1), P(:, 2), col, 'EdgeColor', col);
end

function arcHead(a, r, angFrom, angTo, col, hl)
% filled arrowhead at the end of an arc of radius r, pointing from angFrom to angTo (degrees)
p1 = r * [cosd(angTo), sind(angTo)]; p0 = r * [cosd(angFrom), sind(angFrom)];
u = (p1 - p0) / norm(p1 - p0); n = [-u(2), u(1)];
base = p1 - hl * u;
P = [p1; base + 0.55 * hl * n; base - 0.55 * hl * n];
patch(a, P(:, 1), P(:, 2), col, 'EdgeColor', col);
end