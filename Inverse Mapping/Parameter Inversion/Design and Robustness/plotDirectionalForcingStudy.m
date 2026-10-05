function files = plotDirectionalForcingStudy(results, cfg, run)
%PLOTDIRECTIONALFORCINGSTUDY  Draw and save the V2 figures from saved results (no records, no EMM, no solver).
%   plotDirectionalForcingStudy()                          redraw from the newest Results/V2/V2_*.mat
%   plotDirectionalForcingStudy('...\V2_<stamp>.mat')      redraw from that run
%   plotDirectionalForcingStudy(results, cfg, run)         from a results struct in memory
%
%   Figure 1  inverse bias against the second-direction energy fraction f2, one series per angle phi, log-log:
%             (a) d_sys with the d = 0.5 and 1 limits; (b) |dbeta/beta| with the 1% limit. Registered points and the
%             bisected crossings; the nonlinear diagnostic as an open marker in (b)
%   Figure 2  detectability: (a) the standardised coherence drop D_max against f2 with the detection threshold;
%             (b) the first-crossing fractions per angle: detection against d_sys = 0.5, 1 and |dbeta/beta| = 1%.
%             Coherence warns first where the detection marker lies BELOW the d_sys = 0.5 marker
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'V2', 'V2_20*.mat'));
        if isempty(Lst), error('plotDirectionalForcingStudy:none', 'No V2 results in %s', fullfile(P.results, 'V2')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M0 = load(matFile);
    results = M0.results; cfg = M0.cfg;
    run = struct('id', 'V2', 'dir', fileparts(matFile));
    fprintf('  redrawing V2 figures from %s\n', matFile);
end
S = sharcStyle();
files = {};
g = [results.grid{:}];
X = results.crossings;
TH = cfg.V2.thresholds;
ang = results.angleDeg; nA = numel(ang);
angLab = arrayfun(@(a) sprintf('\\phi = %g\\circ', a), ang, 'UniformOutput', false);
fCross = [X([X.crossed]).f2];
fMin = 10^floor(log10(min([fCross, results.f2]) / 1.5));
xl = [fMin, 0.7];

%% ---- figure 1: bias ----------------------------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.7, 'left', 1.75, 'top', 0.75);
for i = 1:nA
    gi = g(abs([g.phiDeg] - ang(i)) < 1e-6);
    [x, y] = withCrossings(gi, 'dSys', X, ang(i), 'd_sys');
    sharcSeries(ax(1), x, y, i);
    [x, y] = withCrossings(gi, 'absBeta', X, ang(i), '|dbeta/beta|');
    sharcSeries(ax(2), x, y, i);
end
nl = results.nonlinear;
if isfield(nl, 'phiDeg')
    [~, i] = min(abs(ang - nl.phiDeg));
    sharcSeries(ax(2), nl.f2, abs(nl.betaPct), i, 'LineStyle', 'none', 'MarkerSize', 9, 'MarkerFaceColor', 'none');
end
yd = [g.dSys]; yb = abs([g.betaPct]);
if isfield(nl, 'phiDeg'), yb = [yb, abs(nl.betaPct)]; end
set(ax, 'XScale', 'log', 'YScale', 'log');
xlim(ax(1), xl); xlim(ax(2), xl);
ylim(ax(1), [10^floor(log10(min([TH.dSys / 3, yd(yd > 0)]))), 2 * max([yd, 1])]);
ylim(ax(2), [10^floor(log10(min([TH.betaPct / 10, yb(yb > 0)]))), 2 * max([yb, TH.betaPct])]);
for a = 1:2
    fracTicks(ax(a), 'x', xl);
    xlabel(ax(a), 'second-direction energy fraction f_2');
end
ylabel(ax(1), 'd_{sys}'); ylabel(ax(2), '|\delta\beta/\beta| (%)');
sharcPanel(ax(1), 'a', 'Bias significance'); sharcPanel(ax(2), 'b', 'Stiffness bias');
sharcThreshold(ax(1), 'y', TH.dSys(1), 'd = 0.5', 'left');
sharcThreshold(ax(1), 'y', TH.dSys(2), 'd = 1', 'left');
sharcThreshold(ax(2), 'y', TH.betaPct, sprintf('%g%%', TH.betaPct), 'left');
h = arrayfun(@(i) sharcSeries(ax(end), NaN, NaN, i), 1:nA);
hl = angLab;
if isfield(nl, 'phiDeg')
    h(end + 1) = sharcSeries(ax(end), NaN, NaN, i, 'LineStyle', 'none', 'MarkerSize', 9, 'MarkerFaceColor', 'none');
    hl{end + 1} = sprintf('nonlinear (b)');
end
sharcSharedLegend(fig, ax, h, hl, numel(hl));
files = [files, saveSharcFigure(fig, run, 1, 'bias')];

%% ---- figure 2: detectability -------------------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 2.0, 'left', 1.75, 'top', 0.75);
thrD = g(1).threshold;
for i = 1:nA
    gi = g(abs([g.phiDeg] - ang(i)) < 1e-6);
    [x, y] = withCrossings(gi, 'Dmax', X, ang(i), 'detection', thrD);
    sharcSeries(ax(1), x, max(y, eps), i);
end
yD = [g.Dmax];
set(ax(1), 'XScale', 'log', 'YScale', 'log');
xlim(ax(1), xl); ylim(ax(1), [10^floor(log10(min([thrD / 3, yD(yD > 0)]))), 2 * max([yD, thrD])]);
fracTicks(ax(1), 'x', xl);
xlabel(ax(1), 'second-direction energy fraction f_2'); ylabel(ax(1), 'coherence drop D_{max}');
sharcPanel(ax(1), 'a', 'Coherence drop');
sharcThreshold(ax(1), 'y', thrD, sprintf('detectable (%.2f)', thrD), 'left');
% (b) first crossings per angle
crit = {'detection', NaN, 5; 'd_sys', TH.dSys(1), 6; 'd_sys', TH.dSys(2), 7; '|dbeta/beta|', TH.betaPct, 8};
critLab = {'coherence detectable', 'd_{sys} = 0.5', 'd_{sys} = 1', sprintf('|\\delta\\beta/\\beta| = %g%%', TH.betaPct)};
dx = [-0.24 -0.08 0.08 0.24];
for c = 1:size(crit, 1)
    f = NaN(1, nA);
    for i = 1:nA
        x = X(abs([X.phiDeg] - ang(i)) < 1e-6 & strcmp({X.criterion}, crit{c, 1}) & ...
            (isnan(crit{c, 2}) | [X.threshold] == crit{c, 2}));
        if ~isempty(x) && x(1).crossed, f(i) = x(1).f2; end
    end
    sharcSeries(ax(2), (1:nA) + dx(c), f, crit{c, 3}, 'LineStyle', 'none');
end
set(ax(2), 'YScale', 'log');
xlim(ax(2), [0.5 nA + 0.5]); ylim(ax(2), [xl(1) / 10, xl(2)]);   % a decade below for the verdicts
set(ax(2), 'XTick', 1:nA, 'XTickLabel', arrayfun(@(a) sprintf('%g\\circ', a), ang, 'UniformOutput', false));
fracTicks(ax(2), 'y', xl);
xlabel(ax(2), 'second-direction angle \phi'); ylabel(ax(2), 'first-crossing f_2');
sharcPanel(ax(2), 'b', 'Which comes first');
for i = 1:nA
    w = results.warning(i);
    if w.early, s = 'warns'; else, s = 'no warning'; end
    text(ax(2), i, xl(1) / 10 * 1.3, s, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.ink);
end
h = [arrayfun(@(i) sharcSeries(ax(end), NaN, NaN, i), 1:nA), ...
    arrayfun(@(c) sharcSeries(ax(end), NaN, NaN, crit{c, 3}, 'LineStyle', 'none'), 1:size(crit, 1))];
sharcSharedLegend(fig, ax, h, [cellfun(@(s) [s ' (a)'], angLab, 'UniformOutput', false), ...
    cellfun(@(s) [s ' (b)'], critLab, 'UniformOutput', false)], 4);
files = [files, saveSharcFigure(fig, run, 2, 'detectability')];

%% ---- printed check -------------------------------------------------------------------------------------------------
fprintf('\n  V2 figures: %d angles x %d fractions; detection threshold %.3f\n', nA, numel(results.f2), thrD);
for i = 1:nA
    w = results.warning(i);
    fprintf('    phi %2g: f_detect %.4g, f(d_sys = 0.5) %.4g: %s\n', ang(i), w.fDetect, w.fBias, w.verdict);
end
if isfield(nl, 'phiDeg')
    fprintf('    nonlinear (%s): dbeta/beta %+.3f%% (linear %+.3f%%), d gap %.3f\n', nl.caseName, nl.betaPct, ...
        nl.linBetaPct, nl.dGap);
end
end

%% ================================================================================================
function [x, y] = withCrossings(gi, field, X, phiDeg, crit, thrOverride)
% registered points of one angle plus its bisected crossings (value = threshold there), sorted in f2
x = [gi.f2];
switch field
    case 'absBeta', y = abs([gi.betaPct]);
    otherwise, y = [gi.(field)];
end
xc = X(abs([X.phiDeg] - phiDeg) < 1e-6 & strcmp({X.criterion}, crit) & [X.crossed]);
for k = 1:numel(xc)
    x(end + 1) = xc(k).f2; %#ok<AGROW>
    if nargin >= 6, y(end + 1) = thrOverride; else, y(end + 1) = xc(k).threshold; end %#ok<AGROW>
end
[x, o] = sort(x); y = y(o);
end

function fracTicks(ax, axisName, lim)
% energy-fraction ticks on a log axis: decades as 10^{k} (short, so MATLAB does not rotate them) plus 0.1 and 0.5
k = ceil(log10(lim(1))):-1;
t = [10.^k, 0.5];
lab = [arrayfun(@(e) sprintf('10^{%d}', e), k(1:end - 1), 'UniformOutput', false), {'0.1', '0.5'}];
A = upper(axisName);
set(ax, [A 'Tick'], t, [A 'TickLabel'], lab, 'TickLabelInterpreter', 'tex');
try, set(ax, [A 'TickLabelRotation'], 0); catch, end
try, set(ax, [A 'MinorTick'], 'off'); catch, end
end