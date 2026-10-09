function files = plotM1MonteCarloSpotChecks(results, cfg, run)
%PLOTM1MONTECARLOSPOTCHECKS  Draw and save the M1 figures from saved results (no records, no EMM, no solver).
%   plotM1MonteCarloSpotChecks()                     redraw from the newest Results/M1/M1_*.mat
%   plotM1MonteCarloSpotChecks('...\M1_<stamp>.mat') redraw from that run
%   plotM1MonteCarloSpotChecks(results, cfg, run)    from a results struct in memory
%
%   Figure 1  acceptance per case (ticks D1, A1, R2, R1 = the source study of the case): (a) empirical / predicted scatter of ln beta and ln R with +-1 SE, the
%             pre-registered band shaded; (b) 68% and 95% coverage of q with binomial +-1 SE, the 95% band
%             shaded, nominal levels dashed
%   Figure 2  calibration of the predicted uncertainty: empirical CDF of q = dtheta' F dtheta per case against the
%             nominal chi2(2) reference of G3 (q is centred on the truth, as G3a, so with the small Welch bias it is
%             not exactly central chi2(2): its expected mean is 2 + d_W^2, printed per panel)
%   Figure 3  the per-record linearised estimates (theta_hat - theta_true, %) per case against the analytic
%             1-sigma and 95% ellipses of F^-1, with the Welch-expected bias, the empirical mean and the
%             nonlinear estimates of record 1 (both starts)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'M1', 'M1_20*.mat'));
        if isempty(Lst), error('plotM1MonteCarloSpotChecks:none', 'No M1 results in %s', fullfile(P.results, 'M1')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M0 = load(matFile);
    results = M0.results; cfg = M0.cfg;
    run = struct('id', 'M1', 'dir', fileparts(matFile));
    fprintf('  redrawing M1 figures from %s\n', matFile);
end
S = sharcStyle();
files = {};
nC = numel(results.names);
lab = cellfun(@caseLabel, results.names, 'UniformOutput', false);
short = cellfun(@(c) strtok(c, '_'), results.names, 'UniformOutput', false);   % D1, A1, R2, R1 (axis ticks)
A = cfg.M1.accept;
chi = [-2 * log(1 - 0.6827), -2 * log(0.05)];
sm = results.summary; tg = results.targets; st = results.state;

%% ---- figure 1: acceptance ----------------------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.9, 'top', 0.75);
x = 1:nC; dx = 0.12;
xl = [0.4 nC + 0.6];
fillBand(ax(1), xl, A.empPred, S);
fillBand(ax(2), xl, A.coverage95, S);
r = cellfun(@(s) s.empPred(:), sm, 'UniformOutput', false); r = [r{:}];          % 2 x nC
se = cellfun(@(s) s.seEmpPred(:), sm, 'UniformOutput', false); se = [se{:}];
cv = cellfun(@(s) s.coverage(:), sm, 'UniformOutput', false); cv = [cv{:}];     % 2 x nC (68, 95)
n = cellfun(@(s) s.n, sm);
seCv = sqrt(cv .* (1 - cv) ./ [n; n]);
seCv(2, :) = cellfun(@(s) s.seCoverage95, sm);                                     % the rule's own SE (floored)
h = []; hl = {};
for j = 1:2
    xx = x + (2 * j - 3) * dx;
    errBars(ax(1), xx, r(j, :), se(j, :), j, S);
    sharcSeries(ax(1), xx, r(j, :), j, 'LineStyle', 'none');
end
for j = 1:2
    xx = x + (2 * j - 3) * dx;
    errBars(ax(2), xx, cv(3 - j, :), seCv(3 - j, :), 2 + j, S);
    sharcSeries(ax(2), xx, cv(3 - j, :), 2 + j, 'LineStyle', 'none');
end
h(end + 1) = sharcSeries(ax(end), NaN, NaN, 1, 'LineStyle', 'none'); hl{end + 1} = 'ln \beta (a)';
h(end + 1) = sharcSeries(ax(end), NaN, NaN, 2, 'LineStyle', 'none'); hl{end + 1} = 'ln R (a)';
h(end + 1) = sharcSeries(ax(end), NaN, NaN, 3, 'LineStyle', 'none'); hl{end + 1} = '95% coverage (b)';
h(end + 1) = sharcSeries(ax(end), NaN, NaN, 4, 'LineStyle', 'none'); hl{end + 1} = '68% coverage (b)';
for a = 1:2
    xlim(ax(a), xl);
    set(ax(a), 'XTick', x, 'XTickLabel', short, 'XTickLabelRotation', 0);
    xlabel(ax(a), 'M1 case (source study)');
end
ylim(ax(1), [0.75 1.30]); ylim(ax(2), [0.5 1.02]);
ylabel(ax(1), 'empirical / predicted \sigma'); ylabel(ax(2), 'coverage of q');
sharcPanel(ax(1), 'a', 'Scatter'); sharcPanel(ax(2), 'b', 'Coverage');
for k = 1:nC
    text(ax(1), x(k), 0.77, sprintf('n = %d', n(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.axis);
end
sharcThreshold(ax(1), 'y', 1, '', 'tick');
sharcThreshold(ax(2), 'y', 0.95, 'nominal 95%', 'left');
sharcThreshold(ax(2), 'y', 0.6827, 'nominal 68%', 'left');
sharcSharedLegend(fig, ax, h, hl, 4);
files = [files, saveSharcFigure(fig, run, 1, 'acceptance')];

%% ---- figure 2: q against chi2(2) ----------------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [2 2], 12.5, 'vgap', 1.45, 'hgap', 1.6, 'top', 0.75);
qq = linspace(0, 14, 300);
for k = 1:nC
    q = sort(st{k}.q(1:st{k}.nDone));
    sharcSeries(ax(k), qq, 1 - exp(-qq / 2), 7, 'Marker', 'none');
    sharcSeries(ax(k), [0, reshape([q; q], 1, []), 14], [0, reshape([(0:numel(q) - 1); (1:numel(q))] / numel(q), 1, []), 1], ...
        1, 'Marker', 'none');
    xlim(ax(k), [0 14]); ylim(ax(k), [0 1.02]);
    sharcTicks(ax(k), 'x', [0 2 4 6 8 10 12 14]);
    if k > nC - 2, xlabel(ax(k), 'q = \Delta\theta^T F \Delta\theta'); end
    if mod(k, 2) == 1, ylabel(ax(k), 'cumulative fraction'); end
    sharcPanel(ax(k), char('a' + k - 1), lab{k});
    note = {sprintf('n = %d, mean q %.2f', st{k}.nDone, sm{k}.meanQ), sprintf('(2 + d_W^2 = %.2f)', sm{k}.expectedQ), ...
        sprintf('coverage %.3f / %.3f', sm{k}.coverage)};
    for i = 1:3                                  % one text object per line (portable)
        text(ax(k), 13.6, 0.55 - 0.08 * i, note{i}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
            'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.ink);
    end
    sharcThreshold(ax(k), 'x', chi(1), '68%', 'bottom');
    sharcThreshold(ax(k), 'x', chi(2), '95%', 'bottom');
end
h = [sharcSeries(ax(end), NaN, NaN, 7, 'Marker', 'none'), sharcSeries(ax(end), NaN, NaN, 1, 'Marker', 'none')];
sharcSharedLegend(fig, ax, h, {'\chi^2_2 (nominal G3 reference)', 'records (linearised, analytic F)'}, 2);
files = [files, saveSharcFigure(fig, run, 2, 'q_distribution')];

%% ---- figure 3: estimates against the F^-1 ellipses ----------------------------------------------------------------
[fig, ax] = sharcFigure('full', [2 2], 14, 'vgap', 1.45, 'hgap', 1.75, 'left', 1.65, 'top', 0.75);
ang = linspace(0, 2 * pi, 200);
for k = 1:nC
    t = tg{k}; s = st{k};
    D = 100 * s.dtheta(:, 1:s.nDone);
    E = 100 * sqrtm(inv(t.F)) * [cos(ang); sin(ang)];
    sharcSeries(ax(k), D(1, :), D(2, :), 1, 'LineStyle', 'none', 'MarkerSize', 2.5, 'MarkerFaceColor', 'none');
    sharcSeries(ax(k), E(1, :), E(2, :), 7, 'Marker', 'none');
    sharcSeries(ax(k), sqrt(chi(2)) * E(1, :), sqrt(chi(2)) * E(2, :), 7, 'Marker', 'none', 'LineStyle', ':');
    sharcSeries(ax(k), 100 * t.dthetaW(1), 100 * t.dthetaW(2), 4, 'LineStyle', 'none', 'MarkerSize', 7);
    sharcSeries(ax(k), 100 * sm{k}.mean(1), 100 * sm{k}.mean(2), 3, 'LineStyle', 'none', 'MarkerSize', 7);
    for j = 1:numel(s.nl)
        d = 100 * (s.nl(j).theta - t.thetaTrue);
        sharcSeries(ax(k), d(1), d(2), 2, 'LineStyle', 'none', 'MarkerSize', 7 - 3 * (j - 1));
    end
    lim = 1.15 * sqrt(chi(2)) * max(abs(E), [], 2);
    lim = max(lim, 1.05 * max(abs(D), [], 2));
    xlim(ax(k), [-1 1] * lim(1)); ylim(ax(k), [-1 1] * lim(2));
    if k > nC - 2, xlabel(ax(k), '\Delta ln \beta (%)'); end
    if mod(k, 2) == 1, ylabel(ax(k), '\Delta ln R (%)'); end
    sharcPanel(ax(k), char('a' + k - 1), sprintf('%s (n = %d)', lab{k}, s.nDone));
end
h = [sharcSeries(ax(end), NaN, NaN, 1, 'LineStyle', 'none', 'MarkerSize', 2.5, 'MarkerFaceColor', 'none'), ...
    sharcSeries(ax(end), NaN, NaN, 7, 'Marker', 'none'), sharcSeries(ax(end), NaN, NaN, 7, 'Marker', 'none', 'LineStyle', ':'), ...
    sharcSeries(ax(end), NaN, NaN, 4, 'LineStyle', 'none', 'MarkerSize', 7), ...
    sharcSeries(ax(end), NaN, NaN, 3, 'LineStyle', 'none', 'MarkerSize', 7), ...
    sharcSeries(ax(end), NaN, NaN, 2, 'LineStyle', 'none', 'MarkerSize', 8)];
sharcSharedLegend(fig, ax, h, {'records (linearised)', 'F^{-1}, 1\sigma', 'F^{-1}, 95%', 'Welch-expected bias', ...
    'empirical mean', 'nonlinear, record 1 (both starts)'}, 3);
files = [files, saveSharcFigure(fig, run, 3, 'estimates')];

%% ---- printed check ------------------------------------------------------------------------------------------------
fprintf('\n  M1 figures: %d cases\n', nC);
for k = 1:nC
    fprintf('    %-18s n %3d  emp/pred [%.3f %.3f]  cov [%.3f %.3f]  NL-lin max %.3f sigma  %s\n', results.names{k}, ...
        sm{k}.n, sm{k}.empPred, sm{k}.coverage, max(reshape(abs([st{k}.nl.gapSigma]), [], 1)), ternary(sm{k}.pass, 'PASS', 'FAIL'));
end
end

%% ================================================================================================
function s = caseLabel(name)
switch name
    case 'D1_ns2', s = 'D1: n_s = 2';
    case 'A1_miz_R25_h1', s = 'A1: MIZ R 25 m';
    case 'A1_ring_beta_x2', s = 'A1: \beta \times 2';
    case 'R2_wp3.5_gJ7', s = 'R2: \omega_p 3.5, \gamma_J 7';
    case 'R2_wp6_gJ7', s = 'R2: \omega_p 6, \gamma_J 7';
    case 'R1_p0_lsm6dsv16x', s = 'R1: p_0, LSM6DSV16X';
    otherwise, s = strrep(name, '_', ' ');
end
end

function fillBand(ax, xl, band, S)
patch(ax, xl([1 2 2 1]), band([1 1 2 2]), S.fill, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function errBars(ax, x, y, e, slot, S)
for i = 1:numel(x)
    plot(ax, [x(i) x(i)], y(i) + [-1 1] * e(i), '-', 'Color', S.series(slot, :), 'LineWidth', S.thinLineWidth, ...
        'HandleVisibility', 'off');
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end