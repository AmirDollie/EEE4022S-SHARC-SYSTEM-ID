function files = plotA1FloePropertySweep(results, cfg, run)
%PLOTA1FLOEPROPERTYSWEEP  Draw and save the A1 figures from saved results (no EMM, no Fisher recomputation).
%   plotA1FloePropertySweep()                     redraw from the newest Results/A1/A1_*.mat
%   plotA1FloePropertySweep('...\A1_<stamp>.mat') redraw from that run
%   plotA1FloePropertySweep(results, cfg, run)    called by runA1FloePropertySweep
%
%   Colour = final class (reserved class colours); marker = scenario group. Open red markers are the
%   scenarios that failed the truncation metric gate but are Not identifiable at BOTH truncations
%   ("Not identifiable, class robust; sigma not converged"); their sigma is plotted but not quotable.
%
%   Figure 1  identifiability maps: (a) the parameter-space ring and corners in (beta/beta0, R/R0);
%             (b) every scenario in the non-dimensional plane (R/l_f, kR at the JONSWAP peak)
%   Figure 2  (a) sigma_lnbeta and (b) sigma_lnR against R/l_f, thresholds 1% and 5%
%   Figure 3  (a) Welch bias d_W against kR (thresholds 0.5, 1); (b) the feature-level truncation bias
%             of [50 10 10] against [70 15 15] against kR (model discrepancy for V1, NOT a gate)
%   Figure 4  the 12 physical classes: sigma_lnbeta and sigma_lnR per scenario (log scale)
%
%   Also prints and writes A1_reporting.csv: the reporting label of each scenario (the run's final
%   class, plus "class robust, sigma not converged" where that applies) and the R/l_f bracket of the
%   stiffness cliff.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'A1', 'A1_*.mat'));
        if isempty(Lst), error('plotA1FloePropertySweep:none', 'No A1 results in %s', fullfile(P.results, 'A1')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M = load(matFile);
    results = M.results; cfg = M.cfg; %#ok<NASGU>
    run = struct('id', 'A1', 'dir', fileparts(matFile));
    fprintf('  redrawing A1 figures from %s\n', matFile);
end
S = sharcStyle();
sc = results.scenarios;
n = numel(sc);
files = {};

%% ---- collect ------------------------------------------------------------------------------------------------
grpNames = {'ring', 'corner', 'basin', 'tank', 'miz'};
grpLabels = {'parameter ring', 'corner', 'basin p0', 'tank floe', 'MIZ floe (scaled)'};
grpMarker = {'o', 's', 'p', '^', 'd'};
D = struct('ok', false(1, n), 'grp', zeros(1, n), 'sB', NaN(1, n), 'sR', NaN(1, n), 'dW', NaN(1, n), ...
    'kR', NaN(1, n), 'RLf', NaN(1, n), 'dTr', NaN(1, n), 'code', NaN(1, n), 'robust', false(1, n), ...
    'bR', [sc.betaRatio], 'RR', [sc.RRatio]);
label = cell(1, n);
for i = 1:n
    D.grp(i) = find(strcmp(grpNames, sc(i).group));
    r = results.res{i};
    if isempty(r)
        label{i} = sprintf('no result (%s)', results.status{i});
        continue
    end
    D.ok(i) = true;
    D.sB(i) = r.sigmaBeta; D.sR(i) = r.sigmaR; D.dW(i) = r.dW; D.code(i) = r.finalCode;
    D.kR(i) = r.groups.kRpeak; D.RLf(i) = r.groups.RoverLf; D.dTr(i) = r.truncation.d;
    % class robust to truncation although sigma is not converged: Not identifiable at both truncations
    D.robust(i) = r.finalCode == 4 && r.modelAssessed && r.alphaValid && r.nodeGridValid && ...
        r.truncation.done && r.truncation.sameClass && strcmp(r.inverseClass, 'Not identifiable');
    label{i} = r.finalClass;
    if D.robust(i), label{i} = 'Not identifiable (class robust to truncation; sigma not converged)'; end
end
classCol = {S.class.Identifiable, S.class.Marginal, S.class.NotIdentifiable, S.class.NotIdentifiable};

%% ---- printed summary and reporting table --------------------------------------------------------------------
fprintf('\n  A1 reporting labels (final class of the run; robust-class note added where it applies):\n');
for i = 1:n
    fprintf('    %-26s %s\n', sc(i).name, label{i});
end
ok = D.ok & isfinite(D.RLf);
lostB = ok & D.sB > 0.05;  keptB = ok & D.sB <= 0.0125;
if any(lostB) && any(keptB)
    fprintf(['  stiffness cliff in R/l_f: largest R/l_f with sigma_lnbeta > 5%%: %.3g; smallest with <= 1.25%%: %.3g ' ...
        '(no scenario in between)\n'], max(D.RLf(lostB)), min(D.RLf(keptB)));
end
if isfield(run, 'dir') && ~isempty(run.dir)
    rows = cell(n, 9);
    for i = 1:n
        rows(i, :) = {sc(i).name, sc(i).group, D.RLf(i), D.kR(i), D.sB(i), D.sR(i), D.dW(i), D.dTr(i), label{i}};
    end
    files{end + 1} = writeStudyCsv(run, 'reporting', {'scenario', 'group', 'R_over_lf', 'kR_peak', 'sigma_lnbeta', ...
        'sigma_lnR', 'dW', 'trunc_feature_bias_d', 'reporting_label'}, rows);
end

%% ---- figure 1: identifiability maps ---------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.9);
par = find(ismember(D.grp, [1 2 3]));
for i = par
    drawPoint(ax(1), D.bR(i), D.RR(i), i);
end
set(ax(1), 'XScale', 'log', 'YScale', 'log');
xlim(ax(1), [0.18 6]); ylim(ax(1), [0.38 2.7]);
sharcTicks(ax(1), 'x', [0.25 0.5 1 2 4]); sharcTicks(ax(1), 'y', [0.5 0.75 1 1.5 2]);
for i = par
    if D.ok(i)
        text(ax(1), D.bR(i), D.RR(i) * 1.09, sprintf('%.2f%%', 100 * D.sB(i)), 'FontSize', S.smallFontSize - 1, ...
            'FontName', S.font, 'Color', S.axis, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    else
        text(ax(1), D.bR(i), D.RR(i) * 1.09, 'EMM failure', 'FontSize', S.smallFontSize - 1, ...
            'FontName', S.font, 'Color', S.axis, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end
xlabel(ax(1), 'stiffness \beta/\beta_0'); ylabel(ax(1), 'radius R/R_0');
sharcPanel(ax(1), 'a', 'Parameter space (\sigma_{ln\beta} labelled)');
for i = find(ok)
    drawPoint(ax(2), D.RLf(i), D.kR(i), i);
end
set(ax(2), 'XScale', 'log', 'YScale', 'log');
xlim(ax(2), [0.15 15]); ylim(ax(2), [0.3 5]);
sharcTicks(ax(2), 'x', [0.2 0.5 1 2 5 10]); sharcTicks(ax(2), 'y', [0.5 1 2 4]);
xlabel(ax(2), 'floe size / flexural length R/l_f'); ylabel(ax(2), 'kR at the spectral peak');
sharcPanel(ax(2), 'b', 'All scenarios, non-dimensional');
[h, l] = legendEntries();
sharcSharedLegend(fig, ax, h, l, 4);
files = [files, saveSharcFigure(fig, run, 1, 'identifiability_map')];

%% ---- figure 2: uncertainty against R / l_f ------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.9);
for i = find(ok)
    drawPoint(ax(1), D.RLf(i), 100 * D.sB(i), i);
    drawPoint(ax(2), D.RLf(i), 100 * D.sR(i), i);
end
for k = 1:2
    set(ax(k), 'XScale', 'log', 'YScale', 'log');
    xlim(ax(k), [0.15 15]);
    sharcTicks(ax(k), 'x', [0.2 0.5 1 2 5 10]);
    xlabel(ax(k), 'floe size / flexural length R/l_f');
end
ylim(ax(1), [0.1 1000]); sharcTicks(ax(1), 'y', [0.1 1 10 100 1000]);
ylim(ax(2), [0.02 5]); sharcTicks(ax(2), 'y', [0.02 0.1 1 5]);
ylabel(ax(1), '\sigma_{ln\beta} (%)'); ylabel(ax(2), '\sigma_{lnR} (%)');
sharcPanel(ax(1), 'a', 'Stiffness uncertainty'); sharcPanel(ax(2), 'b', 'Radius uncertainty');
for k = 1:2
    sharcThreshold(ax(k), 'y', 100 * cfg.class.identifiable.maxSigma, '1%', 'right');
    sharcThreshold(ax(k), 'y', 100 * cfg.class.marginal.maxSigma, '5%', 'right');
end
[h, l] = legendEntries();
sharcSharedLegend(fig, ax, h, l, 4);
files = [files, saveSharcFigure(fig, run, 2, 'uncertainty_vs_size')];

%% ---- figure 3: bias against kR ----------------------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.9);
for i = find(ok)
    drawPoint(ax(1), D.kR(i), D.dW(i), i);
    drawPoint(ax(2), D.kR(i), D.dTr(i), i);
end
for k = 1:2
    set(ax(k), 'XScale', 'log');
    xlim(ax(k), [0.3 5]); sharcTicks(ax(k), 'x', [0.5 1 2 4]);
    xlabel(ax(k), 'kR at the spectral peak');
end
ylim(ax(1), [0 1.1]); ylim(ax(2), [0 3.3]);
ylabel(ax(1), 'Welch bias significance d_W'); ylabel(ax(2), 'truncation bias d');
sharcPanel(ax(1), 'a', 'Welch bias (enters the class)'); sharcPanel(ax(2), 'b', 'Truncation [50] vs [70] (V1; not a gate)');
sharcThreshold(ax(1), 'y', cfg.class.identifiable.maxBiasDistance, 'd = 0.5', 'right');
sharcThreshold(ax(1), 'y', cfg.class.marginal.maxBiasDistance, 'd = 1', 'right');
sharcThreshold(ax(2), 'y', 0.5, 'd = 0.5', 'right');
[h, l] = legendEntries();
sharcSharedLegend(fig, ax, h, l, 4);
files = [files, saveSharcFigure(fig, run, 3, 'bias_vs_kR')];

%% ---- figure 4: the physical classes ----------------------------------------------------------------------------------
phys = find(ismember(D.grp, [3 4 5]));
[fig, ax] = sharcFigure('full', [1 2], 9.5, 'left', 3.6, 'hgap', 1.0);
ypos = numel(phys):-1:1;
names = strrep({sc(phys).name}, '_', '\_');
for j = 1:numel(phys)
    i = phys(j);
    if ~D.ok(i), continue, end
    drawPoint(ax(1), 100 * D.sB(i), ypos(j), i);
    drawPoint(ax(2), 100 * D.sR(i), ypos(j), i);
end
for k = 1:2
    set(ax(k), 'XScale', 'log', 'YTick', fliplr(ypos), 'YLim', [0.4 numel(phys) + 0.6]);
    set(ax(k), 'YGrid', 'on');
end
set(ax(1), 'YTickLabel', fliplr(names)); set(ax(2), 'YTickLabel', {});
xlim(ax(1), [0.1 1000]); sharcTicks(ax(1), 'x', [0.1 1 10 100 1000]);
xlim(ax(2), [0.02 5]); sharcTicks(ax(2), 'x', [0.02 0.1 1 5]);
xlabel(ax(1), '\sigma_{ln\beta} (%)'); xlabel(ax(2), '\sigma_{lnR} (%)');
sharcPanel(ax(1), 'a', 'Stiffness'); sharcPanel(ax(2), 'b', 'Radius');
for k = 1:2
    sharcThreshold(ax(k), 'x', 100 * cfg.class.identifiable.maxSigma, '1%', 'top');
end
[h, l] = legendEntries();
sharcSharedLegend(fig, ax, h, l, 4);
files = [files, saveSharcFigure(fig, run, 4, 'physical_classes')];

%% ================================================================================================
    function drawPoint(a, x, y, i)
        % one scenario: marker by group, colour by final class; open red = robust Not identifiable;
        % a scenario without a result is a grey cross
        mk = grpMarker{D.grp(i)};
        if ~D.ok(i)
            plot(a, x, y, 'x', 'Color', S.threshold, 'MarkerSize', S.markerSize + 2, 'LineWidth', S.lineWidth);
            return
        end
        c = classCol{D.code(i)};
        if D.code(i) == 4 && ~D.robust(i)
            c = S.threshold;                         % a genuine model-check failure: grey
        end
        face = c;
        if D.code(i) == 4, face = 'w'; end           % open marker: sigma not converged
        plot(a, x, y, 'LineStyle', 'none', 'Marker', mk, 'MarkerSize', S.markerSize + 1, ...
            'MarkerFaceColor', face, 'MarkerEdgeColor', c, 'LineWidth', S.thinLineWidth + 0.4);
    end

    function [h, l] = legendEntries()
        % class colours (filled squares), groups (black markers), and the two special symbols, drawn off-axis
        a = gca;
        h = []; l = {};
        cls = {'Identifiable', 'Marginal', 'Not identifiable'};
        for c = 1:3
            h(end + 1) = plot(a, NaN, NaN, 's', 'MarkerSize', S.markerSize + 1, 'MarkerFaceColor', classCol{c}, ...
                'MarkerEdgeColor', classCol{c}); %#ok<AGROW>
            l{end + 1} = cls{c}; %#ok<AGROW>
        end
        h(end + 1) = plot(a, NaN, NaN, 's', 'MarkerSize', S.markerSize + 1, 'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', S.class.NotIdentifiable, 'LineWidth', S.thinLineWidth + 0.4);
        l{end + 1} = 'Not identifiable, \sigma not converged';
        present = unique(D.grp);
        for g = present(:).'
            h(end + 1) = plot(a, NaN, NaN, grpMarker{g}, 'MarkerSize', S.markerSize + 1, 'MarkerFaceColor', S.ink, ...
                'MarkerEdgeColor', S.ink); %#ok<AGROW>
            l{end + 1} = grpLabels{g}; %#ok<AGROW>
        end
        if any(~D.ok)
            h(end + 1) = plot(a, NaN, NaN, 'x', 'Color', S.threshold, 'MarkerSize', S.markerSize + 2, ...
                'LineWidth', S.lineWidth);
            l{end + 1} = 'no result';
        end
    end
end