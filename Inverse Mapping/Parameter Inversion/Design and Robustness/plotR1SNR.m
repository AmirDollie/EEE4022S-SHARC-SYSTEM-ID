function files = plotR1SNR(results, cfg, run)
%PLOTR1SNR  Draw and save the R1 figure from saved results (no EMM, no Fisher recomputation).
%   plotR1SNR()                     redraw from the newest Results/R1/R1_*.mat
%   plotR1SNR('...\R1_<stamp>.mat') redraw from that run
%   plotR1SNR(results, cfg, run)    called by runR1SNR
%
%   Figure 1 (the SNR window), two stacked panels on one log ASD axis (ug/sqrt(Hz)):
%     (a) sigma_lnbeta and sigma_lnR (ink, labelled directly), 1% threshold
%     (b) Welch bias significance d_W, thresholds 0.5 and 1, refined class boundaries labelled
%   The class of each noise interval is a light band behind both panels, named at the top of (a).
%   The named noise cases are exact evaluations (results.marked), drawn in their sensor colours on both
%   curves and named, with their ASD, in the shared legend.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'R1', 'R1_*.mat'));
        if isempty(Lst), error('plotR1SNR:none', 'No R1 results in %s', fullfile(P.results, 'R1')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M0 = load(matFile);
    results = M0.results; cfg = M0.cfg;
    run = struct('id', 'R1', 'dir', fileparts(matFile));
    fprintf('  redrawing R1 figure from %s\n', matFile);
end
S = sharcStyle();
files = {};
asd = @(s) sqrt(2 * pi * s) / 9.81e-6;              % S_nn per rad/s -> ug/sqrt(Hz)
G = results.grid; Mk = results.marked; B = results.boundaries;
x = asd(results.Snn); xM = asd(results.SnnMarked);
names = results.names; nM = numel(names);
idn = cfg.class.identifiable; mar = cfg.class.marginal;
xl = [0.92 * x(1), 1.08 * x(end)];
iv = asd(results.intervals);                         % Identifiable intervals, ug/sqrt(Hz)

% Two panels on one ASD axis, shared legend below. The class of each noise interval is shown as a light
% band behind BOTH panels (red / amber / green), so no separate "where the cases fall" panel is needed;
% the band names sit at the top of (a), the two refined boundaries are dashed and labelled in (b).
% The noise cases are drawn on the curves in their sensor colours and named in the legend.
[fig, ax] = sharcFigure('full', [2 1], 10.5, 'vgap', 0.55, 'left', 1.6, 'right', 0.5, 'top', 0.75);
yl1 = [5e-3 2]; yl2 = [0 max(1.3, 1.08 * max(Mk.dW))];
ylim(ax(1), yl1); ylim(ax(2), yl2);
set(ax(1), 'YScale', 'log');
for a = 1:2, set(ax(a), 'XScale', 'log'); xlim(ax(a), xl); end

%% ---- class bands behind both panels ---------------------------------------------------------------------------
[segLo, segHi, segCls] = classSegments(results.Snn, G.class, B);
for j = 1:numel(segCls)
    c = classColour(S, segCls{j});
    for a = 1:2
        yl = ylim(ax(a));
        patch(ax(a), asd([segLo(j) segHi(j) segHi(j) segLo(j)]), yl([1 1 2 2]), c, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.12, 'HandleVisibility', 'off');
    end
    text(ax(1), asd(sqrt(segLo(j) * segHi(j))), yl1(2) / 1.25, segCls{j}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.ink);
end

%% ---- curves: neutral ink, labelled directly ---------------------------------------------------------------------
sharcSeries(ax(1), x, 100 * G.sB, 1, 'Marker', 'none', 'Color', S.ink, 'HandleVisibility', 'off');
sharcSeries(ax(1), x, 100 * G.sR, 1, 'Marker', 'none', 'Color', S.axis, 'LineWidth', S.thinLineWidth, ...
    'HandleVisibility', 'off');
sharcSeries(ax(2), x, G.dW, 1, 'Marker', 'none', 'Color', S.ink, 'HandleVisibility', 'off');
kL = round(0.82 * numel(x));                          % label both curves near the noisy end, above the line
text(ax(1), x(kL), 100 * G.sB(kL) * 1.6, '\sigma_{ln\beta}', 'FontSize', S.fontSize, 'FontName', S.font, ...
    'Color', S.ink, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax(1), x(kL), 100 * G.sR(kL) * 1.6, '\sigma_{lnR}', 'FontSize', S.fontSize, 'FontName', S.font, ...
    'Color', S.axis, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%% ---- named noise cases on the curves (legend entries) -----------------------------------------------------------
[~, o] = sort(results.SnnMarked);                    % legend from quietest to noisiest
hh = gobjects(1, nM); lab = cell(1, nM);
for r = 1:nM
    j = o(r); n = names{j};
    sharcSensorSeries(ax(1), xM(j), 100 * Mk.sB(j), n, 'LineStyle', 'none', 'HandleVisibility', 'off');
    hh(r) = sharcSensorSeries(ax(2), xM(j), Mk.dW(j), n, 'LineStyle', 'none');
    lab{r} = sprintf('%s (%.0f)', S.sensorLabel(n), xM(j));
end

%% ---- axes, labels -------------------------------------------------------------------------------------------------
xt = [30 50 100 200 300 500 1000 2000];
for a = 1:2, sharcTicks(ax(a), 'x', xt); end
set(ax(1), 'XTickLabel', {});
sharcTicks(ax(1), 'y', [0.01 0.1 1]);
ylabel(ax(1), 'standard deviation (%)'); ylabel(ax(2), 'bias significance d_W');
xlabel(ax(2), 'accelerometer noise ASD (\mug/\surdHz)');
sharcPanel(ax(1), 'a', 'Precision'); sharcPanel(ax(2), 'b', 'Welch bias');

%% ---- thresholds and boundaries (last: they freeze the limits) -------------------------------------------------
sharcThreshold(ax(1), 'y', 100 * idn.maxSigma, '1%', 'right');                  % 5% lies far above the curves
sharcThreshold(ax(2), 'y', idn.maxBiasDistance, 'd = 0.5', 'right');
sharcThreshold(ax(2), 'y', mar.maxBiasDistance, 'd = 1', 'right');
for b = 1:numel(B)
    xb = asd(B(b).Snn);
    sharcThreshold(ax(1), 'x', xb, '');
    sharcThreshold(ax(2), 'x', xb, sprintf('%.1f', xb), 'top');
end
sharcSharedLegend(fig, ax, hh, lab, 3);
files = [files, saveSharcFigure(fig, run, 1, 'snr_window')];

%% ---- printed check -------------------------------------------------------------------------------------------------
fprintf('\n  R1 figure: Identifiable for');
for r = 1:size(iv, 1), fprintf(' [%.1f, %.1f]', iv(r, 1), iv(r, 2)); end
fprintf(' ug/sqrt(Hz)\n');
for r = 1:nM
    j = o(r);
    fprintf('    %-24s %7.1f ug/rtHz  sig lnb %.3f%%  d_W %.3f  %s\n', S.sensorLabel(names{j}), xM(j), ...
        100 * Mk.sB(j), Mk.dW(j), Mk.class{j});
end
end

%% ================================================================================================
function [lo, hi, cls] = classSegments(Snn, cl, B)
% stretches of constant class along the grid; a change between grid points k and k+1 is placed at the
% refined boundary for that bracket, or at the geometric midpoint if none was refined
lo = Snn(1); hi = []; cls = cl(1);
for k = 1:numel(Snn) - 1
    if ~strcmp(cl{k}, cl{k + 1})
        e = sqrt(Snn(k) * Snn(k + 1));
        b = find(arrayfun(@(q) abs(B(q).bracket(1) - Snn(k)) <= 1e-12 * Snn(k), 1:numel(B)), 1);
        if ~isempty(b), e = B(b).Snn; end
        hi(end + 1) = e; lo(end + 1) = e; cls{end + 1} = cl{k + 1}; %#ok<AGROW>
    end
end
hi(end + 1) = Snn(end);
end

function s = criterionLabel(b)
% readable binding criterion for the annotation
if ~isempty(strfind(b.criterion, 'd_W'))
    s = sprintf('d_W = %g', b.threshold);
elseif ~isempty(strfind(b.criterion, 's_max')) || ~isempty(strfind(b.criterion, 'sigma'))
    s = sprintf('\\sigma_{ln\\beta} = %g%%', 100 * b.threshold);
else
    s = sprintf('%s = %g', b.criterion, b.threshold);
end
end

function c = classColour(S, name)
% Status colour of a class name ('Identifiable', 'Marginal', 'Not identifiable'). The field names in
% sharcStyle are 'NotIdentifiable' etc., so the name is mapped explicitly rather than by removing spaces.
switch lower(strtrim(name))
    case 'identifiable', c = S.class.Identifiable;
    case 'marginal', c = S.class.Marginal;
    otherwise, c = S.class.NotIdentifiable;
end
end