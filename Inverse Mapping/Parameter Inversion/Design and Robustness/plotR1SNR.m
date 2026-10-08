function files = plotR1SNR(results, cfg, run)
%PLOTR1SNR  Draw and save the R1 figure from saved results (no EMM, no Fisher recomputation).
%   plotR1SNR()                     redraw from the newest Results/R1/R1_*.mat
%   plotR1SNR('...\R1_<stamp>.mat') redraw from that run
%   plotR1SNR(results, cfg, run)    called by runR1SNR
%
%   Figure 1 (the SNR window), three stacked panels on one log noise axis: ASD in ug/sqrt(Hz) along the
%   bottom, S_nn per rad/s along the top.
%     (a) sigma_lnbeta and sigma_lnR, thresholds 1% and 5%
%     (b) Welch bias significance d_W, thresholds 0.5 and 1
%     (c) where each named noise case falls, over the class band of the sweep, with the regime of
%         each Marginal stretch (bias-limited or precision-limited) read from the binding metric
%   The Identifiable interval (refined boundaries, results.intervals) is shaded in (a) and (b); each
%   refined boundary is a dashed vertical line, labelled in the panel whose criterion binds there.
%   The named noise cases are exact evaluations (results.marked), not read off the grid; they are named
%   in (c), so the figure has no legend (the two curves in (a) are labelled directly).
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

[fig, ax] = sharcFigure('full', [3 1], 17, 'top', 1.35, 'vgap', 1.0, 'left', 1.6, 'right', 0.6, 'bottom', 1.35);
ylim(ax(1), [5e-3 20]); ylim(ax(2), [0 1.15]); ylim(ax(3), [0.4 nM + 1.5]);
set(ax(1), 'YScale', 'log');
for a = 1:3, set(ax(a), 'XScale', 'log'); xlim(ax(a), xl); end

%% ---- shading first, so the data sit on top ----------------------------------------------------------------------
for a = 1:2
    yl = ylim(ax(a));
    for r = 1:size(iv, 1)
        patch(ax(a), iv(r, [1 2 2 1]), yl([1 1 2 2]), S.fill, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
% (c) class band: one patch per stretch of constant class; edges at the refined boundaries
[segLo, segHi, segCls] = classSegments(results.Snn, G.class, B);
for j = 1:numel(segCls)
    c = classColour(S, segCls{j});
    patch(ax(3), asd([segLo(j) segHi(j) segHi(j) segLo(j)]), [0.4 0.4 nM + 1.5 nM + 1.5], c, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.16, 'HandleVisibility', 'off');
    xm = asd(sqrt(segLo(j) * segHi(j)));             % geometric centre of the stretch
    txt = segCls{j};
    if ~strcmp(segCls{j}, 'Identifiable')            % name the regime from the metric that binds there
        k = find(results.Snn >= segLo(j) & results.Snn <= segHi(j));
        if any(G.sMax(k) > idn.maxSigma) && ~any(G.dW(k) > idn.maxBiasDistance)
            txt = sprintf('%s: precision-limited', segCls{j});
        elseif any(G.dW(k) > idn.maxBiasDistance) && ~any(G.sMax(k) > idn.maxSigma)
            txt = sprintf('%s: bias-limited', segCls{j});
        end
    end
    text(ax(3), xm, nM + 1.0, txt, 'HorizontalAlignment', 'center', 'FontSize', S.smallFontSize, ...
        'FontName', S.font, 'Color', S.ink);
end

%% ---- (a) precision and (b) Welch bias -------------------------------------------------------------------------
sharcSeries(ax(1), x, 100 * G.sB, 6, 'Marker', 'none');
sharcSeries(ax(1), x, 100 * G.sR, 5, 'Marker', 'none');
sharcSeries(ax(2), x, G.dW, 6, 'Marker', 'none');
% curves labelled directly (no legend: the noise cases are named in (c))
k0 = 4;
text(ax(1), x(k0), 100 * G.sB(k0) * 1.9, '\sigma_{ln\beta}', 'FontSize', S.fontSize, 'FontName', S.font, ...
    'Color', S.series(6, :), 'HorizontalAlignment', 'center');
text(ax(1), x(k0), 100 * G.sR(k0) * 1.9, '\sigma_{lnR}', 'FontSize', S.fontSize, 'FontName', S.font, ...
    'Color', S.series(5, :), 'HorizontalAlignment', 'center');

%% ---- named noise cases: on the curves and one row each in (c) -------------------------------------------------
[~, o] = sort(results.SnnMarked);                    % rows from quietest (bottom) to noisiest (top)
for r = 1:nM
    j = o(r); n = names{j};
    sharcSensorSeries(ax(1), xM(j), 100 * Mk.sB(j), n, 'LineStyle', 'none', 'HandleVisibility', 'off');
    sharcSensorSeries(ax(2), xM(j), Mk.dW(j), n, 'LineStyle', 'none', 'HandleVisibility', 'off');
    sharcSensorSeries(ax(3), xM(j), r, n, 'LineStyle', 'none', 'HandleVisibility', 'off');
    lbl = sprintf('%s, %.0f', S.sensorLabel(n), xM(j));      % class: from the band behind it
    if xM(j) < sqrt(xl(1) * xl(2)) * 3          % label right of the marker, or left of it near the right edge
        text(ax(3), xM(j) * 1.08, r, lbl, 'VerticalAlignment', 'middle', 'FontSize', S.smallFontSize, ...
            'FontName', S.font, 'Color', S.ink);
    else
        text(ax(3), xM(j) / 1.08, r, lbl, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
            'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.ink);
    end
end

%% ---- axes, labels -------------------------------------------------------------------------------------------------
xt = [30 50 100 200 300 500 1000 2000];
for a = 1:3, sharcTicks(ax(a), 'x', xt); end
set(ax(1), 'XTickLabel', {}); set(ax(2), 'XTickLabel', {});
sharcTicks(ax(1), 'y', [0.01 0.1 1 10]);
set(ax(3), 'YTick', [], 'YColor', S.axis);
ylabel(ax(1), 'standard deviation (%)'); ylabel(ax(2), 'Welch bias significance d_W');
ylabel(ax(3), 'noise case');
xlabel(ax(3), 'accelerometer noise ASD (\mug/\surdHz)');
text(ax(1), 0.01, 0.97, '(a) Precision', 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', S.fontSize, 'FontName', S.font, 'Color', S.ink);    % (a) inside: the S_nn axis is above it
sharcPanel(ax(2), 'b', 'Welch bias'); sharcPanel(ax(3), 'c', 'Where the noise cases fall');

%% ---- thresholds and boundaries (last: they freeze the limits) -------------------------------------------------
sharcThreshold(ax(1), 'y', 100 * idn.maxSigma, '1%', 'right');
sharcThreshold(ax(1), 'y', 100 * mar.maxSigma, '5%', 'right');
sharcThreshold(ax(2), 'y', idn.maxBiasDistance, 'd = 0.5', 'right');
sharcThreshold(ax(2), 'y', mar.maxBiasDistance, 'd = 1', 'right');
for b = 1:numel(B)
    xb = asd(B(b).Snn);
    lbl = sprintf('%.1f \\mug/\\surdHz: %s', xb, criterionLabel(B(b)));
    isBias = ~isempty(strfind(B(b).criterion, 'd_W'));
    for a = 1:3
        if (a == 1 && ~isBias) || (a == 2 && isBias)
            sharcThreshold(ax(a), 'x', xb, lbl, 'bottom');
        else
            sharcThreshold(ax(a), 'x', xb, '');
        end
    end
end

% S_nn scale above (a): a second x axis on the final position of ax(1)
set(ax(1), 'Units', 'centimeters');
ax2 = axes('Parent', fig, 'Units', 'centimeters', 'Position', get(ax(1), 'Position'), 'Color', 'none', ...
    'XAxisLocation', 'top', 'YTick', [], 'XScale', 'log', 'XLim', xl, 'Box', 'off', ...
    'FontName', S.font, 'FontSize', S.fontSize, 'XColor', S.axis, 'YColor', 'none', 'LineWidth', S.axisLineWidth, ...
    'TickDir', 'out');
dec = -8:-4;
set(ax2, 'XTick', asd(10 .^ dec), 'XTickLabel', arrayfun(@(e) sprintf('10^{%d}', e), dec, 'UniformOutput', false));
try, set(ax2, 'XMinorTick', 'off'); catch, end
xlabel(ax2, 'noise PSD S_{nn} ((m/s^2)^2 per rad/s)');
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