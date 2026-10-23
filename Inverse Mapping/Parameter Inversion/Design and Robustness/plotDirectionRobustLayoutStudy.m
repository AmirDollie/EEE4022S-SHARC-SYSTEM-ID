function files = plotDirectionRobustLayoutStudy(results, cfg, run)
%PLOTDIRECTIONROBUSTLAYOUTSTUDY  Draw and save the D3 figures from saved results (no EMM, no search).
%   plotDirectionRobustLayoutStudy()                       redraw from the newest Results/D3/D3_*.mat
%   plotDirectionRobustLayoutStudy('...\D3_<stamp>.mat')   redraw from that run
%   plotDirectionRobustLayoutStudy(results, cfg, run)      from a results struct in memory
%
%   Figure 1  layouts on the floe (headline noise): rows the D1 optimum (heading 0 known), the D3a robust (worst-case)
%             optimum (heading unknown) and the D3b least-bad layout (smallest worst-case d_tot); columns n_s = 2, 3, 4;
%             each panel gives its worst-case sigma_lnbeta or d_tot. Arrows: D1 the known heading psi = 0;
%             D3a all 24 headings faded, the layout's worst-case heading dark; D3b all headings faded with a
%             schematic two-component sea (primary at 0, secondary grey at phi = 60 deg)
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

%% ---- figure 1: layouts (rows: design rule, columns: n_s; no overlapping markers) -----------------------------------
% Full text width, included at natural size. One colour per row (orange D1, blue D3a, green D3b), all
% markers filled; the floe outline, candidate points and incidence axis match the D1 figures. The wave
% arrow is drawn as in D1: arriving from the theta = 0 side and travelling in -x (the EMM convention;
% the earlier arrow on the left pointing +x showed the opposite direction).
[fig, ax] = sharcFigure('full', [3 nK], 12.5, 'hgap', 0.25, 'vgap', 0.55, 'left', 0.9, 'right', 0.15, ...
    'bottom', 0.2, 'top', 0.6);
Rp = g.Rphys; tt = linspace(0, 2 * pi, 200);
rowName = {'D1 optimum', 'D3a robust', 'D3b least-bad'};
rowSlot = [2 1 3];
fade = 0.80 * [1 1 1];      % faded heading arrows
second = 0.50 * [1 1 1];    % secondary wave component (D3b schematic)
phiShow = 60;               % schematic second-component angle, one of the D3b angles
% worst-case heading of each D3a robust layout (heading estimated), if the per-heading curve was saved
psiWorst = NaN(1, nK);
for k = 1:nK
    rw = pick(k, 'robustWorst');
    if isfield(rw, 'perHeading') && isfield(rw.perHeading, 'sBn')
        [~, iw] = max(rw.perHeading.sBn); psiWorst(k) = h(iw);
    end
end
for k = 1:nK
    d1 = pick(k, 'D1'); rw = pick(k, 'robustWorst'); lb = B(k).leastBad;
    lays = {d1.layout, rw.layout, lb.layout};
    note = {sprintf('worst \\sigma_{ln\\beta} = %.2f%%', 100 * d1.sBworstN), ...
        sprintf('worst \\sigma_{ln\\beta} = %.2f%%', 100 * rw.sBworstN), sprintf('worst d_{tot} = %.3g', lb.dTotWorst)};
    for r = 1:3
        a = ax((r - 1) * nK + k);
        plot(a, Rp * cos(tt), Rp * sin(tt), '-', 'Color', S.ink, 'LineWidth', S.thinLineWidth, 'HandleVisibility', 'off');
        plot(a, [-1.08 1.08] * Rp, [0 0], '--', 'Color', S.threshold, 'LineWidth', S.thinLineWidth, ...
            'HandleVisibility', 'off');                                              % incidence axis
        plot(a, g.xy(:, 1), g.xy(:, 2), 'o', 'Color', S.threshold, 'MarkerSize', 2.5, 'MarkerFaceColor', 'w', ...
            'LineWidth', 0.5, 'HandleVisibility', 'off');                            % candidate points
        sharcSeries(a, g.xy(lays{r}, 1), g.xy(lays{r}, 2), rowSlot(r), 'LineStyle', 'none', 'MarkerSize', 7);
        switch r
            case 1   % D1: one known heading
                drawWave(a, Rp, 0, 1.45, 1.10, S.ink, 1.2);
            case 2   % D3a: heading unknown; every evaluated heading faded, the layout's worst-case heading dark
                for j = 1:numel(h), drawWave(a, Rp, h(j), 1.24, 1.07, fade, 0.7, 0.5); end
                if ~isnan(psiWorst(k)), drawWave(a, Rp, psiWorst(k), 1.34, 1.07, S.ink, 1.2, 0.8); end
            case 3   % D3b: every heading faded; schematic two-component sea, primary and secondary at phi
                for j = 1:numel(h), drawWave(a, Rp, h(j), 1.24, 1.07, fade, 0.7, 0.5); end
                drawWave(a, Rp, 0, 1.34, 1.07, S.ink, 1.2, 0.8);
                drawWave(a, Rp, phiShow, 1.34, 1.07, second, 1.2, 0.8);
                text(a, 1.22 * Rp * cosd(phiShow / 2), 1.22 * Rp * sind(phiShow / 2), '\phi', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.ink);
        end
        axis(a, 'equal'); xlim(a, [-1.30 1.47] * Rp); ylim(a, [-1.58 1.36] * Rp);
        set(a, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
        text(a, 0, -1.37 * Rp, note{r}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.ink);
        if r == 1
            sharcPanel(a, char('a' + k - 1), sprintf('n_s = %d', ns(k)));
        end
        if k == 1
            text(a, -1.30 * Rp, 0, rowName{r}, 'Rotation', 90, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', S.fontSize, 'FontName', S.font, 'Color', S.ink);
        end
    end
end
files = [files, saveSharcFigure(fig, run, 1, 'layouts')];

%% ---- figure 2: sigma vs n_s, efficiency vs heading -------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.7, 'left', 1.75, 'top', 0.75);
series = {'robustWorst', 1, 'D3a robust optimum'; 'robustMean', 4, 'D3a mean-A optimum'; 'D1', 2, 'D1 optimum'};
hh = []; hl = {};
for s = 1:size(series, 1)
    v = arrayfun(@(k) 100 * pick(k, series{s, 1}).sBworstN, 1:nK);
    hh(end + 1) = sharcSeries(ax(1), ns, v, series{s, 2}); hl{end + 1} = [series{s, 3} ' (a)']; %#ok<AGROW>
end
vK = arrayfun(@(k) 100 * pick(k, 'D1').sBworstK, 1:nK);
hh(end + 1) = sharcSeries(ax(1), ns, vK, 2, 'LineStyle', ':', 'MarkerFaceColor', 'none'); hl{end + 1} = 'D1 optimum, heading known (a)';
k4 = find(arrayfun(@(k) any(strcmp({A.perNs(k).sel.name}, 'Level2C')), 1:nK), 1);
if ~isempty(k4)
    hh(end + 1) = sharcSeries(ax(1), ns(k4), 100 * pick(k4, 'Level2C').sBworstN, 3, 'LineStyle', 'none');
    hl{end + 1} = 'Reference configuration (a)';
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
effs = {pick(k2, 'D1'), 2, 'D1 optimum, n_s = 2'; pick(k2, 'robustWorst'), 1, 'D3a robust, n_s = 2'};
if ~isempty(k4), effs(end + 1, :) = {pick(k4, 'Level2C'), 3, 'Reference configuration'}; end
for s = 1:size(effs, 1)
    hh(end + 1) = sharcSeries(ax(2), h, effs{s, 1}.perHeading.effKnown, effs{s, 2}, 'MarkerSize', 3.5); %#ok<AGROW>
    hl{end + 1} = [effs{s, 3} ' (b)']; %#ok<AGROW>
end
xlim(ax(2), [0 360]); ylim(ax(2), [0 1.05]);
sharcTicks(ax(2), 'x', 0:90:360);
xlabel(ax(2), 'wave heading \psi (deg)'); ylabel(ax(2), 'efficiency (best A / A), heading known');
sharcPanel(ax(2), 'b', 'Efficiency versus heading (heading known)');
sharcSharedLegend(fig, ax, hh, hl, 3);
files = [files, saveSharcFigure(fig, run, 2, 'precision_efficiency')];

%% ---- figure 3: nuisance penalty, D3b bias -----------------------------------------------------------------------------
% One legend below both panels, without (a)/(b) suffixes: the three designs belong to (a), the seas to
% (b). In (b) the separation phi is the colour (blue 30, orange 60, green 90 deg), the energy share the
% fill (filled 10%, open 25%), and the spread sea is ink.
[fig, ax] = sharcFigure('full', [1 2], 7.0, 'hgap', 1.7, 'left', 1.75, 'top', 0.75);
hh = []; hl = {};
for s = 1:size(effs, 1)
    q = effs{s, 1}.perHeading;
    hh(end + 1) = sharcSeries(ax(1), h, q.sBn ./ q.sBk, effs{s, 2}, 'MarkerSize', 3.5); %#ok<AGROW>
    hl{end + 1} = effs{s, 3}; %#ok<AGROW>
end
xlim(ax(1), [0 360]); sharcTicks(ax(1), 'x', 0:90:360);
yl = ylim(ax(1)); ylim(ax(1), [min(1, yl(1)), yl(2)]);
xlabel(ax(1), 'wave heading \psi (deg)'); ylabel(ax(1), '\sigma_{ln\beta}(\psi estimated) / \sigma_{ln\beta}(\psi known)');
sharcPanel(ax(1), 'a', 'Heading-nuisance penalty');
sharcThreshold(ax(1), 'y', 1, '', 'right');
% (b) D3b least-bad worst d_tot per sea
nS = numel(results.seaNames);
M = reshape([B.minWorstPerSea], nS, nK);           % nS x nK
sea = cfg.D3.d3b;
phiList = [30 60 90];
for s = 1:nS
    if strcmp(sea(s).kind, 'spread')
        hh(end + 1) = sharcSeries(ax(2), ns, M(s, :), 1, 'Color', S.ink, 'MarkerFaceColor', S.ink, ...
            'MarkerEdgeColor', S.ink, 'MarkerSize', 4); %#ok<AGROW>
        hl{end + 1} = sprintf('cos^{2s} spread, s = %g', sea(s).s); %#ok<AGROW>
    else
        slot = find(phiList == sea(s).phiDeg, 1);
        args = {'MarkerSize', 4};
        if sea(s).f2 > 0.15, args = [args, {'MarkerFaceColor', 'w'}]; end %#ok<AGROW>
        hh(end + 1) = sharcSeries(ax(2), ns, M(s, :), slot, args{:}); %#ok<AGROW>
        hl{end + 1} = sprintf('\\phi = %d%s, f_2 = %d%%', sea(s).phiDeg, char(176), round(100 * sea(s).f2)); %#ok<AGROW>
    end
end
set(ax(2), 'YScale', 'log');
xlim(ax(2), [ns(1) - 0.4, ns(end) + 0.4]); set(ax(2), 'XTick', ns);
ylim(ax(2), [10^floor(log10(min([M(:); cfg.D3.biasGate]) / 2)), 10^ceil(log10(max(M(:)) * 2))]);
xlabel(ax(2), 'number of sensors n_s'); ylabel(ax(2), 'least-bad worst-case d_{tot}');
sharcPanel(ax(2), 'b', 'Directional seas');
for k = 1:nK
    text(ax(2), ns(k), 10^(log10(max(M(:)) * 2) - 0.05), sprintf('%d feasible', B(k).nFeasible), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', S.smallFontSize, 'FontName', S.font, ...
        'Color', S.ink);
end
sharcThreshold(ax(2), 'y', cfg.D3.biasGate, sprintf('gate d_{tot} = %g', cfg.D3.biasGate), 'left');
sharcSharedLegend(fig, ax, hh, hl, 5);
files = [files, saveSharcFigure(fig, run, 3, 'penalty_bias')];

%% ---- printed check -------------------------------------------------------------------------------------------------------
fprintf('\n  D3 figures (%s):\n', A.noise);
for k = 1:nK
    fprintf('    n_s %d: D1 %s worst sigma %.3f%%; worst-case optimum %s %.3f%%; D3b least bad %s d_tot %.3g, %d feasible\n', ...
        ns(k), nameOf(pick(k, 'D1').layout), 100 * pick(k, 'D1').sBworstN, nameOf(pick(k, 'robustWorst').layout), ...
        100 * pick(k, 'robustWorst').sBworstN, nameOf(B(k).leastBad.layout), B(k).leastBad.dTotWorst, B(k).nFeasible);
end
end

%% ================================================================================================
function drawWave(a, R, psiDeg, r0, r1, col, lw, sc)
% incident-wave arrow arriving from the psi side and pointing at the floe centre; psi = 0 reproduces the
% D1 arrow (arriving from theta = 0, travelling in -x). r0, r1: tail and tip radii in units of R;
% sc scales the arrowhead (default 1).
if nargin < 8, sc = 1; end
e = [cosd(psiDeg), sind(psiDeg)]; n = [-e(2), e(1)];
hw = 0.07 * R * sc; hl = 0.12 * R * sc;
tip = r1 * R * e; base = tip + hl * e; tail = r0 * R * e;
plot(a, [tail(1) base(1)], [tail(2) base(2)], '-', 'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
P = [tip; base + hw * n; base - hw * n];
patch(a, P(:, 1), P(:, 2), col, 'EdgeColor', col, 'HandleVisibility', 'off');
end