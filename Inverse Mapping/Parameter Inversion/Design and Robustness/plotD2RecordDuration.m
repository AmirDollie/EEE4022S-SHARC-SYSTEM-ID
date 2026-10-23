function files = plotD2RecordDuration(results, cfg, run)
%PLOTD2RECORDDURATION  Draw and save the D2 figures from saved results (no EMM, a few seconds).
%   plotD2RecordDuration()                      redraw from the newest Results/D2/D2_*.mat
%   plotD2RecordDuration('...\D2_<stamp>.mat')  redraw from that run
%   plotD2RecordDuration(results, cfg, run)     called by runD2RecordDuration at the end of a run
%
%   Figure 1 (headline noise case, 1 x 3): (a) sigma_lnbeta, (b) d_W, (c) RMSE_lnbeta against record
%     length, one line per Welch segment length L (admissible combinations only), log-log; y-ranges
%     fitted to the data ((a) and (c) shared), bias limits d = 0.5, 1 in (b); dashed vertical line at
%     the production record (109 min); shared legend below.
%   Figure 2 (all noise cases): the identifiable window [T_lo, T_hi] of every L as a horizontal bar,
%     one row per noise case, so "record at least X min with L = Y" can be read directly. A window still
%     open at the end of the search (T_hi > 600 min) ends in an arrow. Only L with at least one window
%     appear (L = 512 is never identifiable; say so in the caption).
%   L is coloured by categorical slot 1-4 in both figures (D2 does not colour noise cases).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'D2', 'D2_*.mat'));
        if isempty(Lst), error('plotD2RecordDuration:none', 'No D2 results in %s', fullfile(P.results, 'D2')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M = load(matFile);
    results = M.results; cfg = M.cfg;
    run = struct('id', 'D2', 'dir', fileparts(matFile));
    fprintf('  redrawing D2 figures from %s\n', matFile);
end
S = sharcStyle();
R = results.R; Ls = results.L; T = results.Treal; nL = numel(Ls);
noise = results.noise; nN = numel(noise);
h = find(strcmp(noise, results.headline));
files = {};

%% ---- figure 1: headline noise case against record length -----------------------------------------------
% one row of three panels across the text width, shared legend below. The y-ranges are fitted to the
% data: the 1% precision limit lies about an order of magnitude above every curve, so drawing it would
% leave most of panels (a) and (c) empty (the caption states this instead).
[fig, ax] = sharcFigure('full', [1 3], 6.0, 'hgap', 1.45);
hh = gobjectsOr(nL); lab = cell(1, nL);
for a = 1:nL
    ok = results.admit(a, :);
    hh(a) = sharcSeries(ax(1), T(ok), 100 * squeeze(R.sigB(h, a, ok)), a);
    sharcSeries(ax(2), T(ok), squeeze(R.dW(h, a, ok)), a);
    sharcSeries(ax(3), T(ok), 100 * squeeze(R.rmseB(h, a, ok)), a);
    lab{a} = sprintf('L = %d', Ls(a));
end
ylab = {'\sigma_{ln\beta} (%)', 'bias significance d_W', 'RMSE_{ln\beta} (%)'};
for k = 1:3
    set(ax(k), 'XScale', 'log', 'YScale', 'log');
    xlim(ax(k), [8 300]);
    xlabel(ax(k), 'record length (min)');
    ylabel(ax(k), ylab{k});
end
yS = 100 * [R.sigB(h, :, :), R.rmseB(h, :, :)]; yS = yS(isfinite(yS));
yLim = [0.8 * min(yS), 1.25 * max(yS)];                 % (a) and (c) share a range so they compare directly
ylim(ax(1), yLim); ylim(ax(3), yLim);
yD = R.dW(h, :, :); yD = yD(isfinite(yD));
dLim = [0.8 * min(yD), max(1.25, 1.15 * max(yD))];      % keep both bias limits (0.5 and 1) in view
ylim(ax(2), dLim);
for k = 1:3, sharcTicks(ax(k), 'x', [10 30 60 109 240]); end
sharcTicks(ax(1), 'y', logTicks(yLim)); sharcTicks(ax(3), 'y', logTicks(yLim));
sharcTicks(ax(2), 'y', logTicks(dLim));
if yLim(2) >= 1, sharcThreshold(ax(1), 'y', 1, '', 'tick'); sharcThreshold(ax(3), 'y', 1, '', 'tick'); end
sharcThreshold(ax(2), 'y', 0.5, '', 'tick'); sharcThreshold(ax(2), 'y', 1, '', 'tick');  % d = 0.5, 1
for k = 1:3, sharcThreshold(ax(k), 'x', 109.2); end                                % production record
sharcPanel(ax(1), 'a', 'Random error'); sharcPanel(ax(2), 'b', 'Welch bias'); sharcPanel(ax(3), 'c', 'Total error');
sharcSharedLegend(fig, ax, hh, lab, nL);
files = [files, saveSharcFigure(fig, run, 1, sprintf('duration_%s', results.headline))];

%% ---- figure 2: identifiable window per noise case and L ------------------------------------------------
[fig, ax] = sharcFigure('full', [1 1], 7.2, 'left', 2.9);
off = linspace(-0.27, 0.27, nL);
xMax = 600;
hasW = any(isfinite(results.TidentLo), 1);       % L with at least one identifiable window
hh = gobjectsOr(nnz(hasW)); k = 0;
for a = find(hasW)
    k = k + 1; hh(k) = sharcSeries(ax, NaN, NaN, a);  % legend handle
end
for n = 1:nN
    y = nN + 1 - n;                                  % first noise case at the top
    for a = 1:nL
        lo = results.TidentLo(n, a); hi = min(results.TidentHi(n, a), xMax);
        if isnan(lo), continue; end
        sharcSeries(ax, [lo hi], [y y] + off(a), a, 'Marker', 'none', 'LineWidth', 3.2, 'HandleVisibility', 'off');
        sharcSeries(ax, lo, y + off(a), a, 'LineStyle', 'none', 'HandleVisibility', 'off');   % start = shortest record
        if isinf(results.TidentHi(n, a))             % still identifiable at the end of the search
            sharcSeries(ax, xMax, y + off(a), a, 'LineStyle', 'none', 'Marker', '>', 'HandleVisibility', 'off');
        end
    end
end
set(ax, 'XScale', 'log');
xlim(ax, [3 xMax * 1.12]); ylim(ax, [0.4 nN + 0.6]);   % x room for the open-end arrows
sharcTicks(ax, 'x', [3 10 30 60 109 240 600]);
set(ax, 'YTick', 1:nN, 'YTickLabel', cellfun(S.sensorLabel, fliplr(noise), 'UniformOutput', false));
xlabel(ax, 'record length (min)');
sharcThreshold(ax, 'x', 109.2, 'production', 'top');
sharcPanel(ax, 'a', 'Identifiable record lengths by noise case (marker: shortest)');
sharcSharedLegend(fig, ax, hh, arrayfun(@(x) sprintf('L = %d', x), Ls(hasW), 'UniformOutput', false), nnz(hasW));
files = [files, saveSharcFigure(fig, run, 2, 'identifiable_window')];
end

%% ================================================================================================
function t = logTicks(lim)
c = [0.01 0.02 0.03 0.05 0.1 0.2 0.3 0.5 1 2 3 5 10 20 50 100];
t = c(c >= lim(1) & c <= lim(2));
end

function h = gobjectsOr(n)
if exist('gobjects', 'builtin') || exist('gobjects', 'file'), h = gobjects(1, n); else, h = zeros(1, n); end
end