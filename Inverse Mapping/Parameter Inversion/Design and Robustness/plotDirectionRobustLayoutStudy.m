function files = plotDirectionRobustLayoutStudy(results, cfg, run)
%PLOTDIRECTIONROBUSTLAYOUTSTUDY  Draw and save the D3 figures from saved results (no EMM, no search).
%   plotDirectionRobustLayoutStudy()                       redraw from the newest Results/D3/D3_*.mat
%   plotDirectionRobustLayoutStudy('...\D3_<stamp>.mat')   redraw from that run
%   plotDirectionRobustLayoutStudy(results, cfg, run)      from a results struct in memory
%
%   Figure 1  layouts on the floe, n_s = 2, 3, 4 (headline noise): the D1 optimum (known heading 0), the D3a worst-case
%             optimum (unknown heading) and the D3b least-bad layout (smallest worst d_tot); the arrow is heading 0
%   Figure 2  (a) worst-case sigma_lnbeta over the 24 headings (heading estimated) against n_s for the D3a worst-case and
%             mean optima, the D1 optimum and Level 2C, with the D1 optimum with the heading known for contrast; the 1%
%             Identifiable limit; (b) efficiency against heading (heading known): best A at that heading / A of the
%             design, for the D1 optimum, the D3a worst-case optimum (n_s = 2) and Level 2C
%   Figure 3  (a) the heading-nuisance penalty sigma_lnbeta(psi estimated) / sigma_lnbeta(psi known) against heading;
%             (b) D3b: the least-bad worst-case d_tot per sea against n_s, with the d_tot = 0.5 gate and the number of
%             bias-feasible layouts
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'D3', 'D3_20*.mat'));
        if isempty(Lst), error('plotD3:none', 'No D3 results in %s', fullfile(P.results, 'D3')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M0 = load(matFile);
    results = M0.results; cfg = M0.cfg;
    run = struct('id', 'D3', 'dir', fileparts(matFile));
    fprintf('  redrawing D3 figures from %s\n', matFile);
end
S = sharcStyle();
files = {};
A = results.d3a(1); B = results.d3b; g = results.geom;
ns = results.nsList; nK = numel(ns);
h = results.headingsDeg;
pick = @(k, name) A.perNs(k).sel(strcmp({A.perNs(k).sel.name}, name));
nameOf = @(lay) strjoin(g.names(lay), '+');

%% ---- figure 1: layouts ----------------------------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 nK], [], 'hgap', 0.8, 'left', 0.6, 'right', 0.3, 'bottom', 0.6);
Rp = g.Rphys; tt = linspace(0, 2 * pi, 200);
for k = 1:nK
    a = ax(k);
    plot(a, Rp * cos(tt), Rp * sin(tt), '-', 'Color', S.axis, 'LineWidth', S.thinLineWidth, 'HandleVisibility', 'off');
    plot(a, g.xy(:, 1), g.xy(:, 2), '.', 'Color', S.grid * 0.8, 'MarkerSize', 6, 'HandleVisibility', 'off');
    d1 = pick(k, 'D1'); rw = pick(k, 'robustWorst'); lb = B(k).leastBad;
    sharcSeries(a, g.xy(d1.layout, 1), g.xy(d1.layout, 2), 2, 'LineStyle', 'none', 'MarkerSize', 9, 'MarkerFaceColor', 'none');
    sharcSeries(a, g.xy(rw.layout, 1), g.xy(rw.layout, 2), 1, 'LineStyle', 'none', 'MarkerSize', 6);
    sharcSeries(a, g.xy(lb.layout, 1), g.xy(lb.layout, 2), 3, 'LineStyle', 'none', 'MarkerSize', 4);
    quiver(a, -1.25 * Rp, 0, 0.35 * Rp, 0, 0, 'Color', S.axis, 'LineWidth', S.thinLineWidth, 'MaxHeadSize', 0.8, ...
        'HandleVisibility', 'off');
    axis(a, 'equal'); xlim(a, 1.3 * Rp * [-1 1]); ylim(a, 1.15 * Rp * [-1 1]);
    set(a, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    sharcPanel(a, char('a' + k - 1), sprintf('n_s = %d', ns(k)));
end
hh = [sharcSeries(ax(end), NaN, NaN, 2, 'LineStyle', 'none', 'MarkerSize', 9, 'MarkerFaceColor', 'none'), ...
    sharcSeries(ax(end), NaN, NaN, 1, 'LineStyle', 'none', 'MarkerSize', 6), ...
    sharcSeries(ax(end), NaN, NaN, 3, 'LineStyle', 'none', 'MarkerSize', 4)];
sharcSharedLegend(fig, ax, hh, {'D1 optimum (heading 0 known)', 'D3a worst-case optimum (heading unknown)', ...
    'D3b least-bad (smallest worst d_{tot})'}, 3);
files = [files, saveSharcFigure(fig, run, 1, 'layouts')];

%% ---- figure 2: sigma vs n_s, efficiency vs heading -------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.7, 'left', 1.75, 'top', 0.75);
series = {'robustWorst', 1, 'D3a worst-case optimum'; 'robustMean', 4, 'D3a mean optimum'; 'D1', 2, 'D1 optimum'};
hh = []; hl = {};
for s = 1:size(series, 1)
    v = arrayfun(@(k) 100 * pick(k, series{s, 1}).sBworstN, 1:nK);
    hh(end + 1) = sharcSeries(ax(1), ns, v, series{s, 2}); hl{end + 1} = [series{s, 3} ' (a)']; %#ok<AGROW>
end
vK = arrayfun(@(k) 100 * pick(k, 'D1').sBworstK, 1:nK);
hh(end + 1) = sharcSeries(ax(1), ns, vK, 2, 'LineStyle', ':', 'MarkerFaceColor', 'none'); hl{end + 1} = 'D1 optimum, heading known (a)';
k4 = find(arrayfun(@(k) any(strcmp({A.perNs(k).sel.name}, 'Level2C')), 1:nK), 1);
if ~isempty(k4)
    hh(end + 1) = sharcSeries(ax(1), ns(k4), 100 * pick(k4, 'Level2C').sBworstN, 5, 'LineStyle', 'none');
    hl{end + 1} = 'Level 2C (a)';
end
set(ax(1), 'YScale', 'log');
xlim(ax(1), [ns(1) - 0.4, ns(end) + 0.4]);
yy = vK;
for s = 1:size(series, 1), yy = [yy, arrayfun(@(k) 100 * pick(k, series{s, 1}).sBworstN, 1:nK)]; end %#ok<AGROW>
if ~isempty(k4), yy = [yy, 100 * pick(k4, 'Level2C').sBworstN]; end
yy = yy(isfinite(yy));
ylim(ax(1), [10^floor(log10(min(yy) * 0.8)), 10^ceil(log10(max([yy, 1.2])))]);
set(ax(1), 'XTick', ns);
xlabel(ax(1), 'number of sensors n_s'); ylabel(ax(1), 'worst-case \sigma_{ln\beta} over headings (%)');
sharcPanel(ax(1), 'a', 'Precision, heading unknown');
sharcThreshold(ax(1), 'y', 100 * cfg.class.identifiable.maxSigma, 'Identifiable limit', 'left');
% (b) efficiency vs heading, n_s = 2 and Level 2C
k2 = find(ns == 2, 1);
effs = {pick(k2, 'D1'), 2, 'D1 optimum, n_s = 2'; pick(k2, 'robustWorst'), 1, 'D3a worst-case optimum, n_s = 2'};
if ~isempty(k4), effs(end + 1, :) = {pick(k4, 'Level2C'), 5, 'Level 2C'}; end
for s = 1:size(effs, 1)
    hh(end + 1) = sharcSeries(ax(2), h, effs{s, 1}.perHeading.effKnown, effs{s, 2}, 'MarkerSize', 3.5); %#ok<AGROW>
    hl{end + 1} = [effs{s, 3} ' (b)']; %#ok<AGROW>
end
xlim(ax(2), [0 360]); ylim(ax(2), [0 1.05]);
sharcTicks(ax(2), 'x', 0:45:360);
xlabel(ax(2), 'wave heading \psi (deg)'); ylabel(ax(2), 'efficiency (best A / A), heading known');
sharcPanel(ax(2), 'b', 'Cost of placing for the wrong heading');
sharcSharedLegend(fig, ax, hh, hl, 3);
files = [files, saveSharcFigure(fig, run, 2, 'precision_efficiency')];

%% ---- figure 3: nuisance penalty, D3b bias -----------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.7, 'left', 1.75, 'top', 0.75);
hh = []; hl = {};
for s = 1:size(effs, 1)
    q = effs{s, 1}.perHeading;
    hh(end + 1) = sharcSeries(ax(1), h, q.sBn ./ q.sBk, effs{s, 2}, 'MarkerSize', 3.5); %#ok<AGROW>
    hl{end + 1} = [effs{s, 3} ' (a)']; %#ok<AGROW>
end
xlim(ax(1), [0 360]); sharcTicks(ax(1), 'x', 0:45:360);
yl = ylim(ax(1)); ylim(ax(1), [min(1, yl(1)), yl(2)]);
xlabel(ax(1), 'wave heading \psi (deg)'); ylabel(ax(1), '\sigma_{ln\beta}(\psi estimated) / \sigma_{ln\beta}(\psi known)');
sharcPanel(ax(1), 'a', 'Heading-nuisance penalty');
sharcThreshold(ax(1), 'y', 1, '', 'right');
% (b) D3b least-bad worst d_tot per sea
nS = numel(results.seaNames);
M = reshape([B.minWorstPerSea], nS, nK);           % nS x nK
for s = 1:nS
    hh(end + 1) = sharcSeries(ax(2), ns, M(s, :), s + (s >= 5), 'MarkerSize', 4); %#ok<AGROW>
    hl{end + 1} = [results.seaNames{s} ' (b)']; %#ok<AGROW>
end
set(ax(2), 'YScale', 'log');
xlim(ax(2), [ns(1) - 0.4, ns(end) + 0.4]); set(ax(2), 'XTick', ns);
ylim(ax(2), [10^floor(log10(min([M(:); cfg.D3.biasGate]) / 2)), 10^ceil(log10(max(M(:)) * 2))]);
xlabel(ax(2), 'number of sensors n_s'); ylabel(ax(2), 'least-bad worst-case d_{tot}');
sharcPanel(ax(2), 'b', 'D3b: can placement rescue the inverse?');
for k = 1:nK
    text(ax(2), ns(k), 10^(log10(max(M(:)) * 2) - 0.05), sprintf('%d feasible', B(k).nFeasible), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', S.smallFontSize, 'FontName', S.font, ...
        'Color', S.ink);
end
sharcThreshold(ax(2), 'y', cfg.D3.biasGate, sprintf('gate d_{tot} = %g', cfg.D3.biasGate), 'left');
sharcSharedLegend(fig, ax, hh, hl, 4);
files = [files, saveSharcFigure(fig, run, 3, 'penalty_bias')];

%% ---- printed check -------------------------------------------------------------------------------------------------------
fprintf('\n  D3 figures (%s):\n', A.noise);
for k = 1:nK
    fprintf('    n_s %d: D1 %s worst sigma %.3f%%; worst-case optimum %s %.3f%%; D3b least bad %s d_tot %.3g, %d feasible\n', ...
        ns(k), nameOf(pick(k, 'D1').layout), 100 * pick(k, 'D1').sBworstN, nameOf(pick(k, 'robustWorst').layout), ...
        100 * pick(k, 'robustWorst').sBworstN, nameOf(B(k).leastBad.layout), B(k).leastBad.dTotWorst, B(k).nFeasible);
end
end