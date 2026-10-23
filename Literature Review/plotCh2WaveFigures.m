% plotCh2WaveFigures  Chapter 2 background figures (no project dependencies).
%   Figure 1: a harmonic wave in space eta(x) at fixed t, and in time eta(t) at fixed x,
%             annotated with amplitude a, wave height 2a, wavelength lambda and period T.
%   Figure 2: JONSWAP spectra at the thesis reference H_s, all normalised to the same
%             variance (same H_s):
%               (a) gamma_J = 1, 3.3, 7 at omega_p = 5 rad/s   (peakedness)
%               (b) omega_p = 3.5, 5, 7.5 rad/s at gamma_J = 3.3 (peak frequency, R2 range)
%   Outputs: ch2_harmonic_wave.pdf/.png and ch2_jonswap_peakedness.pdf/.png in the
%   current folder. Copy them to figs/ch2/.
%   Include Figure 2 in LaTeX WITHOUT a width option, so the 9 pt text stays 9 pt.
1;

function exportFig(fig, name)
% vector PDF plus PNG preview; exportgraphics where available (MATLAB R2020a+)
if exist('exportgraphics', 'file') || exist('exportgraphics', 'builtin')
    exportgraphics(fig, [name '.pdf'], 'ContentType', 'vector');
    exportgraphics(fig, [name '.png'], 'Resolution', 300);
else
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, [name '.pdf'], '-dpdf');
    print(fig, [name '.png'], '-dpng', '-r300');
end
end

function dimArrow(ax, x0, y0, x1, y1, txt, tx, ty, ha, col, fs)
% double-headed dimension line with a label, in data coordinates
plot(ax, [x0 x1], [y0 y1], '-', 'Color', col, 'LineWidth', 0.8);
d = [x1 - x0, y1 - y0]; L = hypot(d(1), d(2)); u = d / L;
xl = xlim(ax); yl = ylim(ax);
s = [diff(xl), diff(yl)] * 0.018;                          % head size per axis
n = [-u(2), u(1)];
for e = [0 1]
    p = [x0 y0] + e * d; v = (1 - 2 * e) * u;              % head points inwards from each end
    q1 = p + v .* s * 1.6 + n .* s * 0.7; q2 = p + v .* s * 1.6 - n .* s * 0.7;
    patch(ax, [p(1) q1(1) q2(1)], [p(2) q1(2) q2(2)], col, 'EdgeColor', col);
end
text(ax, tx, ty, txt, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle', ...
    'FontSize', fs, 'Color', col, 'Interpreter', 'tex');
end

function S = jonswap(w, wp, gam, Hs)
% JONSWAP shape normalised so that 4 sqrt(m0) = Hs. m0 is integrated over a wide
% grid, not the plotted window, so curves whose tail leaves the window are not inflated.
shape = @(x) x.^-5 .* exp(-1.25 * (wp ./ x).^4) .* ...
    gam.^exp(-(x - wp).^2 ./ (2 * (0.07 * (x <= wp) + 0.09 * (x > wp)).^2 * wp^2));
wFull = linspace(0.2 * wp, 20 * wp, 20000);
S = shape(w) * (Hs / 4)^2 / trapz(wFull, shape(wFull));
end

% ---------------------------------------------------------------- style
FS = 10; FSL = 9;
ink = [0.10 0.10 0.10]; wave = [0.13 0.40 0.70]; grey = [0.55 0.55 0.55];
cols = [0.13 0.40 0.70; 0.85 0.45 0.10; 0.70 0.15 0.20];

% ================================================================ Figure 1
a = 1; lambda = 2 * pi; T = 2 * pi;                         % unit wave (shape only)
x = linspace(0, 2.25 * lambda, 600); t = linspace(0, 2.25 * T, 600);
fig = figure('Units', 'centimeters', 'Position', [2 2 16 5.6], 'Color', 'w');

ax1 = axes(fig, 'Units', 'normalized', 'Position', [0.08 0.18 0.40 0.70]);
hold(ax1, 'on');
plot(ax1, x, a * cos(2 * pi * x / lambda - pi / 2), '-', 'Color', wave, 'LineWidth', 1.6);
plot(ax1, x([1 end]), [0 0], ':', 'Color', grey, 'LineWidth', 0.8);   % still-water level
xlim(ax1, [x(1) x(end)]); ylim(ax1, [-1.75 1.75]);
xc = lambda / 4; xt = 3 * lambda / 4;                       % first crest and trough
dimArrow(ax1, xc, 1.25, xc + lambda, 1.25, '\lambda', xc + lambda / 2, 1.5, 'center', ink, FSL);
plot(ax1, [xc xc], [1 1.25], ':', 'Color', grey); plot(ax1, [xc + lambda xc + lambda], [1 1.25], ':', 'Color', grey);
dimArrow(ax1, xt, -1, xt, 0, 'a', xt + 0.2, -0.5, 'left', ink, FSL);
xh = xt + lambda + 0.35;                                    % wave height beside the second trough
plot(ax1, [xc + lambda xh + 0.2], [1 1], ':', 'Color', grey);
dimArrow(ax1, xh, -1, xh, 1, '2a', xh - 0.25, 0, 'right', ink, FSL);   % wave height (H is depth in this thesis)
set(ax1, 'XTick', [], 'YTick', 0, 'YTickLabel', {'0'}, 'FontSize', FS, 'Box', 'off', ...
    'XColor', ink, 'YColor', ink, 'TickDir', 'out');
xlabel(ax1, 'horizontal position x', 'FontSize', FS);
ylabel(ax1, '\eta(x, t_0)', 'FontSize', FS);
title(ax1, '(a) Snapshot in space', 'FontSize', FS, 'FontWeight', 'normal');

ax2 = axes(fig, 'Units', 'normalized', 'Position', [0.58 0.18 0.40 0.70]);
hold(ax2, 'on');
plot(ax2, t, a * cos(2 * pi * t / T - pi / 2), '-', 'Color', wave, 'LineWidth', 1.6);
plot(ax2, t([1 end]), [0 0], ':', 'Color', grey, 'LineWidth', 0.8);
xlim(ax2, [t(1) t(end)]); ylim(ax2, [-1.75 1.75]);
tc = T / 4;
dimArrow(ax2, tc, 1.25, tc + T, 1.25, 'T', tc + T / 2, 1.5, 'center', ink, FSL);
plot(ax2, [tc tc], [1 1.25], ':', 'Color', grey); plot(ax2, [tc + T tc + T], [1 1.25], ':', 'Color', grey);
set(ax2, 'XTick', [], 'YTick', 0, 'YTickLabel', {'0'}, 'FontSize', FS, 'Box', 'off', ...
    'XColor', ink, 'YColor', ink, 'TickDir', 'out');
xlabel(ax2, 'time t', 'FontSize', FS);
ylabel(ax2, '\eta(x_0, t)', 'FontSize', FS);
title(ax2, '(b) Record at a fixed point', 'FontSize', FS, 'FontWeight', 'normal');
exportFig(fig, 'ch2_harmonic_wave');

% ================================================================ Figure 2
FS2 = 9; FSL2 = 8.5;                                        % final size: include at natural width
Hs = 0.05; wp0 = 5.0; gam0 = 3.3;                           % thesis reference sea (Sec. 3.4)
gammas = [1 3.3 7];
wps    = [3.5 5.0 7.5];                                     % R2 peak-frequency range
colsW  = [0.55 0.74 0.90; 0.13 0.40 0.70; 0.04 0.17 0.38];  % light to dark: low to high omega_p
w = linspace(2, 12, 2000);

fig = figure('Units', 'centimeters', 'Position', [2 2 16 6.2], 'Color', 'w');

% (a) peakedness at the reference peak frequency
ax = axes(fig, 'Units', 'normalized', 'Position', [0.08 0.18 0.39 0.70]); hold(ax, 'on');
hl = zeros(1, 3); lab = cell(1, 3); Smax = 0;
for i = 1:3
    S = 1e4 * jonswap(w, wp0, gammas(i), Hs); Smax = max(Smax, max(S));
    hl(i) = plot(ax, w, S, '-', 'Color', cols(i, :), 'LineWidth', 1.5);
    lab{i} = sprintf('\\gamma_J = %g', gammas(i));
end
% two lines keep the legend narrow enough to stay right of the peak; char(8211) is an en dash
lab{gammas == 1} = sprintf('\\gamma_J = 1\n(Pierson%sMoskowitz)', char(8211));
plot(ax, [wp0 wp0], [0 Smax * 1.08], ':', 'Color', grey, 'LineWidth', 0.8);
text(ax, wp0 + 0.15, Smax * 1.04, '\omega_p', 'FontSize', FSL2, 'Color', ink);
xlim(ax, [w(1) w(end)]); ylim(ax, [0 Smax * 1.12]);
set(ax, 'FontSize', FS2, 'Box', 'off', 'XColor', ink, 'YColor', ink, 'TickDir', 'out');
xlabel(ax, 'angular frequency \omega (rad/s)', 'FontSize', FS2);
ylabel(ax, 'S_\eta(\omega)  (10^{-4} m^2 s/rad)', 'FontSize', FS2);
title(ax, sprintf('(a) Peakedness, \\omega_p = %g rad/s', wp0), 'FontSize', FS2, 'FontWeight', 'normal');
lg = legend(ax, hl, lab, 'Location', 'northeast', 'Box', 'off', 'FontSize', FSL2);
drawnow;                                                    % let the legend size itself first
lp = lg.Position;                                           % then push it flush with the right axis edge
lp(1) = ax.Position(1) + ax.Position(3) - lp(3);
lg.Position = lp;

% (b) peak frequency at the reference peakedness
ax = axes(fig, 'Units', 'normalized', 'Position', [0.59 0.18 0.39 0.70]); hold(ax, 'on');
hl = zeros(1, 3); lab = cell(1, 3); Smax = 0;
for i = 1:3
    S = 1e4 * jonswap(w, wps(i), gam0, Hs); Smax = max(Smax, max(S));
    hl(i) = plot(ax, w, S, '-', 'Color', colsW(i, :), 'LineWidth', 1.5);
    lab{i} = sprintf('\\omega_p = %g rad/s', wps(i));
end
xlim(ax, [w(1) w(end)]); ylim(ax, [0 Smax * 1.12]);
set(ax, 'FontSize', FS2, 'Box', 'off', 'XColor', ink, 'YColor', ink, 'TickDir', 'out');
xlabel(ax, 'angular frequency \omega (rad/s)', 'FontSize', FS2);
title(ax, sprintf('(b) Peak frequency, \\gamma_J = %g', gam0), 'FontSize', FS2, 'FontWeight', 'normal');
legend(ax, hl, lab, 'Location', 'northeast', 'Box', 'off', 'FontSize', FSL2);

exportFig(fig, 'ch2_jonswap_peakedness');
fprintf('Saved ch2_harmonic_wave and ch2_jonswap_peakedness (.pdf, .png)\n');