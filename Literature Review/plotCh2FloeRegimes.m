function plotCh2FloeRegimes(outDir)
% plotCh2FloeRegimes  Chapter 2 schematic of rigid-body and flexural regimes.
%   plotCh2FloeRegimes()        saves to the current folder
%   plotCh2FloeRegimes(outDir)  saves to outDir (e.g. 'figs/ch2')
%
%   Four cross-sections of a floe in a regular wave, for combinations of
%   k*R_phys (floe size vs wavelength) and R_phys/ell_f (floe size vs flexural
%   length). Floe shapes come from a quasi-static 1D free-free beam on a
%   buoyant (Winkler) foundation driven by the wave surface:
%       ell_f^4 w'''' + w = eta(x),   w'' = w''' = 0 at x = +-R.
%   Inertia and the scattered wave are neglected: illustrative only, NOT the EMM.
%   Output: ch2_floe_regimes.pdf / .png
%
%   Written as a function file so local functions work in both MATLAB and Octave.

if nargin < 1, outDir = pwd; end
if ~exist(outDir, 'dir'), mkdir(outDir); end

% ---------------------------------------------------------------- style
FS = 9; FSL = 8.5;
ink = [0.10 0.10 0.10]; wave = [0.13 0.40 0.70]; grey = [0.55 0.55 0.55];
floeFill = [0.86 0.86 0.86];
surfI = [0.85 0.40 0.05];    % incident surface under the floe (orange, distinct from floe and wave)

% ---------------------------------------------------------------- cases
%          kR    R/ell_f   title
cases = { 0.4,   1.0,  '(a) Short floe: rides the wave';
          2.0,   0.5,  '(b) Stiff floe: heave and pitch';
          2.0,   2.0,  '(c) Comparable scales: bending set by D';
          2.0,  10.0,  '(d) Compliant floe: follows the surface'};
phi  = pi / 3;               % wave phase at the snapshot (gives both heave and pitch)
xr   = [-0.75 0.75];         % horizontal window in wavelengths
yr   = [-2.0 1.5];           % raise yr(2) to 1.7 if panel (b) clips at the top
tTop = 0.10; tBot = 0.20;    % drawn freeboard and draught (schematic, a = 1)

figW = 16; figH = 9;         % cm; set figH = 9.5 if the titles crowd the row above
fig = figure('Units', 'centimeters', 'Position', [2 2 figW figH], 'Color', 'w');
pos = [0.04 0.55 0.44 0.31; 0.54 0.55 0.44 0.31;
       0.04 0.12 0.44 0.31; 0.54 0.12 0.44 0.31];

for c = 1:4
    kR = cases{c, 1}; Rl = cases{c, 2};
    [xf, wf, R] = beamShape(kR, Rl, phi, 400);
    xw = linspace(xr(1), xr(2), 800);
    eta = cos(2 * pi * xw + phi);
    out = abs(xw) > R;
    etaO = eta; etaO(~out) = NaN;          % open-water surface
    etaI = eta; etaI(out) = NaN;           % incident surface beneath the floe

    ax = axes(fig, 'Units', 'normalized', 'Position', pos(c, :)); hold(ax, 'on');
    plot(ax, xr, [0 0], ':', 'Color', grey, 'LineWidth', 0.7);           % still water
    plot(ax, xw, etaO, '-', 'Color', wave, 'LineWidth', 1.4);
    patch(ax, [xf fliplr(xf)], [wf + tTop, fliplr(wf - tBot)], floeFill, ...
        'EdgeColor', ink, 'LineWidth', 0.9);
    % incident surface beneath the floe, drawn on top in a contrasting colour
    plot(ax, xw, etaI, '--', 'Color', surfI, 'LineWidth', 1.5);
    % floe length bracket
    yb = -1.45;
    plot(ax, [-R R], [yb yb], '-', 'Color', ink, 'LineWidth', 0.7);
    plot(ax, [-R -R], yb + [-0.1 0.1], '-', 'Color', ink, 'LineWidth', 0.7);
    plot(ax, [R R], yb + [-0.1 0.1], '-', 'Color', ink, 'LineWidth', 0.7);
    text(ax, 0, yb - 0.32, '$2R_{\rm phys}$', 'Interpreter', 'latex', ...
        'FontSize', FSL, 'Color', ink, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim(ax, xr); ylim(ax, yr);
    set(ax, 'YTick', [], 'FontSize', FS, 'Box', 'off', 'XColor', ink, 'YColor', 'none', ...
        'TickDir', 'out', 'XTick', -0.5:0.5:0.5);
    if c >= 3
        xlabel(ax, '$x/\lambda$', 'Interpreter', 'latex', 'FontSize', FS);
    else
        set(ax, 'XTickLabel', {});
    end
    % two-line title: heading, then the regime parameters in the same font, centred
    sub = sprintf('\\fontsize{%g}{\\itk}R_{phys} = %g,   R_{phys}/{\\it%s}_f = %g', ...
                  FSL, kR, char(8467), Rl);
    title(ax, {cases{c, 3}, sub}, 'Interpreter', 'tex', ...
          'FontSize', FS, 'FontWeight', 'normal', 'Color', ink);
end

name = fullfile(outDir, 'ch2_floe_regimes');
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
function [x, w, R] = beamShape(kR, Rl, phi, n)
% Quasi-static free-free beam on a Winkler foundation under eta = cos(kx + phi).
% Lengths in wavelengths (k = 2*pi). Returns the deflection w on x in [-R, R].
k = 2 * pi; R = kR / k; l = R / Rl;
x = linspace(-R, R, n);
A = exp(1i * phi);                        % complex surface amplitude
P = A / (1 + (k * l)^4);                  % particular (forced) response
mu = [1 + 1i, 1 - 1i, -1 + 1i, -1 - 1i] / (sqrt(2) * l);   % roots of l^4 mu^4 + 1 = 0
x0 = R * sign(real(mu));                  % anchor growing/decaying parts at the edge they grow toward
M = zeros(4); rhs = zeros(4, 1); row = 0;
for xb = [-R R]
    for m = [2 3]                          % zero bending moment and shear force
        row = row + 1;
        M(row, :) = mu.^m .* exp(mu .* (xb - x0));
        rhs(row) = -P * (1i * k)^m * exp(1i * k * xb);
    end
end
cc = M \ rhs;
w = P * exp(1i * k * x);
for i = 1:4
    w = w + cc(i) * exp(mu(i) * (x - x0(i)));
end
w = real(w);
end