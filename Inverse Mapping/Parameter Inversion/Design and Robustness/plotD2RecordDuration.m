function files = plotD2RecordDuration(results, cfg, run)
%PLOTD2RECORDDURATION  Draw and save the D2 figures from saved results (no EMM, a few seconds).
%   plotD2RecordDuration()                      redraw from the newest Results/D2/D2_*.mat
%   plotD2RecordDuration('...\D2_<stamp>.mat')  redraw from that run
%   plotD2RecordDuration(results, cfg, run)     called by runD2RecordDuration at the end of a run
%
%   Figure 1 (headline noise case, 2 x 2): (a) sigma_lnbeta, (b) d_W, (c) RMSE_lnbeta against record
%     length, one line per Welch segment length L (admissible combinations only), log-log; thresholds
%     as ticks (sigma = 1%; d = 0.5, 1); dashed vertical line at the production record (109 min);
%     legend in the fourth quadrant.
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
[fig, axAll] = sharcFigure('full', [2 2], [], 'vgap', 1.75);   % 2 x 2: three data panels, legend in the fourth
ax = axAll(1:3);
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
yLim = [10^floor(log10(min(yS))), 10^ceil(log10(max([yS(:); 1.5])))];
ylim(ax(1), yLim); ylim(ax(3), yLim);
yD = R.dW(h, :, :); yD = yD(isfinite(yD));
ylim(ax(2), [10^floor(log10(min(yD))), 10^ceil(log10(max([yD(:); 1.5])))]);
for k = 1:3, sharcTicks(ax(k), 'x', [10 30 60 109 240]); end
sharcTicks(ax(1), 'y', logTicks(yLim)); sharcTicks(ax(3), 'y', logTicks(yLim));
sharcTicks(ax(2), 'y', logTicks(get(ax(2), 'YLim')));
sharcThreshold(ax(1), 'y', 1, '', 'tick');                                         % sigma = 1%
sharcThreshold(ax(2), 'y', 0.5, '', 'tick'); sharcThreshold(ax(2), 'y', 1, '', 'tick');  % d = 0.5, 1
sharcThreshold(ax(3), 'y', 1, '', 'tick');
for k = 1:3, sharcThreshold(ax(k), 'x', 109.2); end                                % production record
sharcPanel(ax(1), 'a', 'Random error'); sharcPanel(ax(2), 'b', 'Welch bias'); sharcPanel(ax(3), 'c', 'Total error');
% legend in the empty fourth quadrant, so no panel loses height
set(axAll(4), 'Units', 'centimeters'); p4 = get(axAll(4), 'Position');
lg = sharcLegend(ax(3), hh, lab);
set(lg, 'Units', 'centimeters');
lp = get(lg, 'Position');
set(lg, 'Position', [p4(1) + 0.3, p4(2) + (p4(4) - lp(4)) / 2, lp(3), lp(4)]);
delete(axAll(4));
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
c = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100];
t = c(c >= lim(1) & c <= lim(2));
end

function h = gobjectsOr(n)
if exist('gobjects', 'builtin') || exist('gobjects', 'file'), h = gobjects(1, n); else, h = zeros(1, n); end
end