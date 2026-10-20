function plotTwinInversionVerification(matFile)
%PLOTTWININVERSIONVERIFICATION Thesis figure for Section 4.7: end-to-end inversion at floes A and B.
%   plotTwinInversionVerification()            uses the default result file (run 20260927_143438)
%   plotTwinInversionVerification(matFile)     uses the given transmissibilityTwinInversion_*.mat
%
%   Reads the saved 'out' struct of runTransmissibilityTwinInversion; no recomputation.
%   Reproduces the runner's own figure (same layout, labels and styling), with the
%   left-column legends fixed in the top-left corner, and exports it.
%   Writes fig_twin_inversion_verification.pdf (vector) and .png (300 dpi) next to the .mat file.

%% ---- Options --------------------------------------------------------------------------------

%% ---- Load -----------------------------------------------------------------------------------
if nargin < 1 || isempty(matFile)
    here = fileparts(mfilename('fullpath'));
    cands = {fullfile(here, '..', 'Results', 'transmissibilityTwinInversion_20260927_143438.mat'), ...
             fullfile(here, 'transmissibilityTwinInversion_20260927_143438.mat'), ...
             'transmissibilityTwinInversion_20260927_143438.mat'};
    matFile = '';
    for c = 1:numel(cands)
        if isfile(cands{c}), matFile = cands{c}; break; end
    end
    if isempty(matFile)
        [f, pth] = uigetfile('transmissibilityTwinInversion_*.mat', 'Select the twin-inversion result file');
        if isequal(f, 0), error('No result file selected.'); end
        matFile = fullfile(pth, f);
    end
end
S = load(matFile, 'out'); out = S.out;
res = out.truth; nT = numel(res); nVal = out.nVal;
fprintf('Loaded %s (%d truth points, %d validation records)\n', matFile, nT, nVal);
if all(arrayfun(@(r) isempty(r.nonlinear), res))
    warning('This file has no nonlinear inversions (N_NL = 0); use run 20260927_143438.');
end


%% ---- Figure (as in runTransmissibilityTwinInversion, section 4) ---------------------------------
colA = [42 120 214] / 255; colB = [235 104 52] / 255; colC = [120 120 120] / 255; ink = [0.35 0.35 0.33];
% Same size and proportions as the runner's original figure
fig = figure('Color', 'w', 'Position', [60 60 1250 420 * nT]);
FS = 10; FSL = 8;
ang = linspace(0, 2 * pi, 200); qq = linspace(0, 12, 200);
floeParams = {'1.2\beta_0, 0.97R_0', '0.8\beta_0, 1.03R_0'};
for it = 1:nT
    q = res(it);
    ax = subplot(nT, 2, 2 * it - 1); hold(ax, 'on');
    plot(ax, q.G3a.theta(1, :), q.G3a.theta(2, :), '.', 'Color', 1 - 0.5 * (1 - colA), 'MarkerSize', 7);
    E = sqrtm(inv(q.F)) * [cos(ang); sin(ang)];
    h1 = plot(ax, E(1, :), E(2, :), '-', 'Color', colA, 'LineWidth', 1.8);
    h2 = plot(ax, 2.448 * E(1, :), 2.448 * E(2, :), '--', 'Color', colA, 'LineWidth', 1.2);
    h3 = plot(ax, q.dthetaWelch(1), q.dthetaWelch(2), 'd', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
    hs = [h1 h2 h3]; lab = {'F^{-1} 1\sigma', '95% ellipse', 'Welch-expected bias'};
    for j = 1:numel(q.nonlinear)
        d = q.nonlinear(j).theta - q.thetaTrue;
        hn = plot(ax, d(1), d(2), 'p', 'Color', colB, 'MarkerFaceColor', colB, 'MarkerSize', 12);
        if j == 1, hs(end + 1) = hn; lab{end + 1} = 'nonlinear estimate'; end %#ok<AGROW>
    end
    set(ax, 'XColor', ink, 'YColor', ink, 'FontSize', FS, 'Box', 'off'); grid(ax, 'on');
    xlabel(ax, '\Delta ln \beta (from truth)'); ylabel(ax, '\Delta ln R (from truth)');
    title(ax, sprintf('(%c) Floe %c (%s): linearised estimates', 'a' + 2 * (it - 1), 'A' + (it - 1), floeParams{it}), ...
        'FontWeight', 'normal');
    legend(ax, hs, lab, 'Location', 'northwest', 'Box', 'off', 'FontSize', FSL);
    ax = subplot(nT, 2, 2 * it); hold(ax, 'on');
    plot(ax, qq, 1 - exp(-qq / 2), '-', 'Color', 'k', 'LineWidth', 1.2);
    plot(ax, sort(q.G3a.q), (1:nVal) / nVal, '-', 'Color', colA, 'LineWidth', 1.8);
    plot(ax, sort(q.G3b.q(:, 2)), (1:nVal) / nVal, ':', 'Color', colB, 'LineWidth', 1.8);
    plot(ax, sort(q.G3b.q(:, 4)), (1:nVal) / nVal, '-', 'Color', colB, 'LineWidth', 1.8);
    plot(ax, sort(q.G3c.q), (1:nVal) / nVal, '-', 'Color', colC, 'LineWidth', 1.2);
    xlim(ax, [0 12]); ylim(ax, [0 1]);
    set(ax, 'XColor', ink, 'YColor', ink, 'FontSize', FS, 'Box', 'off'); grid(ax, 'on');
    xlabel(ax, 'q = \Delta\theta^T C^{-1} \Delta\theta'); ylabel(ax, 'empirical CDF');
    title(ax, sprintf('(%c) Floe %c: calibration of reported uncertainty', 'b' + 2 * (it - 1), 'A' + (it - 1)), ...
        'FontWeight', 'normal');
    legend(ax, {'\chi^2_2 (ideal)', 'oracle covariance', 'jackknife (Hartlap)', 'pooled jackknife (Hartlap)', 'nominal covariance \Sigma(p_0)'}, ...
        'Location', 'southeast', 'Box', 'off', 'FontSize', FSL);
end

%% ---- Export ---------------------------------------------------------------------------------
[pth, ~] = fileparts(matFile); if isempty(pth), pth = pwd; end
base = fullfile(pth, 'fig_twin_inversion_verification');
drawnow;
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [base '.png'], 'Resolution', 300);
fprintf('Saved %s.pdf and .png\n', base);
end