function files = plotV1ModelMismatch(results, cfg, run, ladder)
%PLOTV1MODELMISMATCH  Draw and save the V1 figures from saved results (no records, no EMM, no solver).
%   plotV1ModelMismatch()                            redraw from the newest Results/V1/V1_*.mat and the newest
%                                                    stage 0 ladder Results/V1/diagnostics/truncationLadder_*.mat
%   plotV1ModelMismatch('...\V1_<stamp>.mat')        redraw from that run (newest ladder)
%   plotV1ModelMismatch(results, cfg, run, ladder)   from structs in memory
%
%   Figure 1  truncation: (a) stage 0 ladder at p0, successive-step d against M with the d = 0.1 convergence
%             rule, plus the last step at the large-kR floes and the P check; (b) V1a per floe, [50 10 10]
%             inverse against the selected truth: d_W (matched), d_sys and d_tot, with the class change
%   Figure 2  bias significance: (a) d_sys and d_tot against the water-depth error; (b) d_sys against the
%             common radial sensor error, log-log, both signs, with the bisected tolerance crossings, a
%             slope-1 reference and the 1% beta crossings
%   Figure 3  parameter bias: dbeta/beta and dR/R (%) against (a) depth error and (b) radial error, with the
%             nonlinear inversion of the largest-d_sys case against its linear projection
%
%   d_W and the matched class per floe come from V1_truncation.csv of the same run (the .mat keeps them for
%   p0 only); without the CSV those two items are omitted.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'V1', 'V1_20*.mat'));
        if isempty(Lst), error('plotV1ModelMismatch:none', 'No V1 results in %s', fullfile(P.results, 'V1')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M0 = load(matFile);
    results = M0.results; cfg = M0.cfg;
    run = struct('id', 'V1', 'dir', fileparts(matFile));
    fprintf('  redrawing V1 figures from %s\n', matFile);
end
if nargin < 4 || isempty(ladder)
    ladder = newestLadder(fullfile(run.dir, 'diagnostics'), results.truncation.truth);
end
S = sharcStyle();
files = {};
tr = results.truncation;
nF = numel(tr.floes);
rowT = [tr.rows{:}];
mm = [results.mismatch{:}];
csv = readRunCsv(fullfile(run.dir, 'V1_truncation.csv'));
TH = cfg.V1.toleranceThresholds;
lab = @(t) sprintf('[%d %d %d]', t);

%% ---- figure 1: truncation ----------------------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.7, 'left', 1.75, 'top', 0.75);
% (a) ladder
Lv = ladder.levels;
M = arrayfun(@(l) l.tr(1), Lv); dStep = [Lv.d];
ok = ~isnan(dStep);
Mt = ladder.truth(1);
sharcSeries(ax(1), M(ok), dStep(ok), 1);
for i = 2:min(numel(ladder.dLastStep), 4)
    sharcSeries(ax(1), Mt + 2.2 * (i - 1), ladder.dLastStep(i), i, 'LineStyle', 'none');
end
sharcSeries(ax(1), Mt - 2.2, ladder.pCheck.d, 5, 'LineStyle', 'none');
xlim(ax(1), [min(M) - 8, max(M) + 14]);
ylim(ax(1), [0, max([dStep(ok), ladder.dLastStep, ladder.pCheck.d, ladder.rule]) * 1.25]);
sharcTicks(ax(1), 'x', M);
xlabel(ax(1), sprintf('M (with [P N] = [%d %d])', Lv(1).tr(2:3)));
ylabel(ax(1), 'step distance d');
sharcPanel(ax(1), 'a', 'Stage 0 truncation ladder');
for k = find(ok)
    text(ax(1), M(k), dStep(k), sprintf('  %.3f', dStep(k)), 'FontSize', S.smallFontSize, 'FontName', S.font, ...
        'Color', S.axis, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
sharcThreshold(ax(1), 'y', ladder.rule, sprintf('rule d = %g', ladder.rule), 'left');
sharcThreshold(ax(1), 'x', Mt, 'selected', 'top');
% (b) V1a per floe
x = 1:nF;
dSys = [rowT.dSys]; dTot = [rowT.dTot];
dW = NaN(1, nF); clsM = repmat({''}, 1, nF);
if ~isempty(csv)
    for i = 1:nF
        j = find(strcmp(csv.floe, tr.floes{i}), 1);
        if ~isempty(j), dW(i) = str2double(csv.d_W{j}); clsM{i} = csv.class_matched{j}; end
    end
end
if strcmp(tr.floes{1}, 'basin_p0') && isnan(dW(1))
    dW(1) = results.nominal.dW; clsM{1} = results.nominal.class;
end
dx = 0.2;
sharcSeries(ax(2), x - dx, dW, 6, 'LineStyle', 'none');
sharcSeries(ax(2), x, dSys, 7, 'LineStyle', 'none');
sharcSeries(ax(2), x + dx, dTot, 8, 'LineStyle', 'none');
yTop = 1.32 * max([dSys dTot dW 1]);
xlim(ax(2), [0.4 nF + 0.6]); ylim(ax(2), [0 yTop]);
set(ax(2), 'XTick', x, 'XTickLabel', cellfun(@floeLabel, tr.floes, 'UniformOutput', false));
xlabel(ax(2), sprintf('floe (inverse %s, truth %s)', lab(tr.inverse), lab(tr.truth)));
ylabel(ax(2), 'Fisher distance d');
sharcPanel(ax(2), 'b', 'V1a: truncation bias');
for i = 1:nF
    s = classShort(rowT(i).class);
    if ~isempty(clsM{i}), s = sprintf('%s \\rightarrow %s', classShort(clsM{i}), s); end
    text(ax(2), x(i), max([dSys(i) dTot(i) dW(i)]) + 0.05 * yTop, s, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.ink);
end
sharcThreshold(ax(2), 'y', TH.dSys(1), '', 'tick');
sharcThreshold(ax(2), 'y', TH.dSys(2), '', 'tick');
h = [sharcSeries(ax(end), NaN, NaN, 1), sharcSeries(ax(end), NaN, NaN, 2, 'LineStyle', 'none'), ...
    sharcSeries(ax(end), NaN, NaN, 3, 'LineStyle', 'none'), sharcSeries(ax(end), NaN, NaN, 4, 'LineStyle', 'none'), ...
    sharcSeries(ax(end), NaN, NaN, 5, 'LineStyle', 'none'), sharcSeries(ax(end), NaN, NaN, 6, 'LineStyle', 'none'), ...
    sharcSeries(ax(end), NaN, NaN, 7, 'LineStyle', 'none'), sharcSeries(ax(end), NaN, NaN, 8, 'LineStyle', 'none')];
hl = [{'p_0 ladder (a)'}, cellfun(@(f) [floeLabel(f) ', last step (a)'], tr.floes(2:min(nF, 4)), 'UniformOutput', false), ...
    {sprintf('p_0, P \\rightarrow %d (a)', ladder.pCheck.to(2)), 'd_W, matched model (b)', 'd_{sys} (b)', 'd_{tot} (b)'}];
sharcSharedLegend(fig, ax, h(1:numel(hl)), hl, 4);
files = [files, saveSharcFigure(fig, run, 1, 'truncation')];

%% ---- figure 2: bias significance against depth and radial error ---------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.7, 'top', 0.75);
% (a) depth, with the nominal point (d_sys = 0, d_tot = d_W)
dep = mm(strcmp({mm.kind}, 'depth'));
[e, o] = sort([0, 100 * [dep.value]]);
ys = [0, [dep.dSys]]; yt = [results.nominal.dW, [dep.dTot]];
sharcSeries(ax(1), e, ys(o), 7);
sharcSeries(ax(1), e, yt(o), 8);
xlim(ax(1), [min(e) - 0.6, max(e) + 0.6]); ylim(ax(1), [0 1.15 * max([TH.dSys ys yt])]);
sharcTicks(ax(1), 'x', e);
xlabel(ax(1), 'water-depth error \epsilon_H (%)'); ylabel(ax(1), 'Fisher distance d');
sharcPanel(ax(1), 'a', 'V1b: water depth');
sharcThreshold(ax(1), 'y', TH.dSys(1), 'd = 0.5', 'right');
sharcThreshold(ax(1), 'y', TH.dSys(2), 'd = 1', 'right');
% (b) radial, log-log, both signs
rad = mm(strcmp({mm.kind}, 'radial'));
tol = results.tolerances;
tR = tol(strcmp({tol.kind}, 'radial') & [tol.crossed] & isfinite([tol.x]));
isD = strcmp({tR.criterion}, 'd_sys'); isB = ~isD;
for sd = [-1 1]
    slot = 3.5 + sd / 2;                                         % 3 negative, 4 positive
    rs = rad(sign([rad.value]) == sd);
    td = tR(isD & [tR.side] == sd);
    xx = 1000 * abs([[td.x], rs.value]); yy = [[td.threshold], rs.dSys];
    [xx, o] = sort(xx);
    sharcSeries(ax(2), xx, yy(o), slot);
end
dRef = tR(isD);
if ~isempty(dRef)
    k = mean([dRef.threshold] ./ (1000 * abs([dRef.x])));      % local d per mm from the crossings
    xr = [0.01 40];
    sharcSeries(ax(2), xr, k * xr, 7, 'Marker', 'none', 'LineStyle', ':');
end
set(ax(2), 'XScale', 'log', 'YScale', 'log');
xlim(ax(2), [0.02 40]); ylim(ax(2), [0.2 2 * max([rad.dSys])]);
sharcTicks(ax(2), 'x', [0.05 0.1 0.3 1 3 10 20]);
sharcTicks(ax(2), 'y', [0.5 1 10 100]);
xlabel(ax(2), 'common radial sensor error |\Delta r| (mm)'); ylabel(ax(2), 'd_{sys}');
sharcPanel(ax(2), 'b', 'V1c: radial sensor position');
sharcThreshold(ax(2), 'y', TH.dSys(1), '', 'tick');
sharcThreshold(ax(2), 'y', TH.dSys(2), '', 'tick');
bx = 1000 * abs([tR(isB).x]);
if ~isempty(bx)
    sharcThreshold(ax(2), 'x', mean(bx), sprintf('|\\delta\\beta/\\beta| = %g%%', TH.betaPct), 'bottom');
end
h = [sharcSeries(ax(end), NaN, NaN, 7), sharcSeries(ax(end), NaN, NaN, 8), sharcSeries(ax(end), NaN, NaN, 3), ...
    sharcSeries(ax(end), NaN, NaN, 4), sharcSeries(ax(end), NaN, NaN, 7, 'Marker', 'none', 'LineStyle', ':')];
sharcSharedLegend(fig, ax, h, {'d_{sys} (a)', 'd_{tot} (a)', '\Delta r < 0 (b)', '\Delta r > 0 (b)', ...
    'slope 1 through the crossings (b)'}, 5);
files = [files, saveSharcFigure(fig, run, 2, 'significance')];

%% ---- figure 3: parameter bias ---------------------------------------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.7, 'left', 1.75, 'top', 0.75);
nl = results.nonlinear;
[nlKind, nlX] = strtok(nl.caseName);
nlX = str2double(nlX);
grp = {dep, rad}; sc = [100 1000];
xlab = {'water-depth error \epsilon_H (%)', 'common radial sensor error \Delta r (mm)'};
ttl = {'V1b: water depth', 'V1c: radial sensor position'};
kinds = {'depth', 'radial'};
for a = 1:2
    g = grp{a};
    [xx, o] = sort([0, sc(a) * [g.value]]);
    yb = [0, [g.betaPct]]; yr = [0, [g.RPct]];
    sharcSeries(ax(a), xx, yb(o), 1);
    sharcSeries(ax(a), xx, yr(o), 2);
    if strcmp(nlKind, kinds{a})
        sharcSeries(ax(a), sc(a) * nlX, nl.betaPct, 1, 'LineStyle', 'none', 'MarkerSize', 9, 'MarkerFaceColor', 'none');
        sharcSeries(ax(a), sc(a) * nlX, nl.RPct, 2, 'LineStyle', 'none', 'MarkerSize', 9, 'MarkerFaceColor', 'none');
    end
    yl = [min([yb yr nl.betaPct * strcmp(nlKind, kinds{a})]), max([yb yr nl.betaPct * strcmp(nlKind, kinds{a})])];
    pad = 0.12 * diff(yl);
    xlim(ax(a), [min(xx), max(xx)] + [-1 1] * 0.08 * diff([min(xx), max(xx)]));
    ylim(ax(a), yl + [-pad pad]);
    sharcTicks(ax(a), 'x', xx);
    xlabel(ax(a), xlab{a});
    if a == 1, ylabel(ax(a), 'parameter bias (%)'); end
    sharcPanel(ax(a), char('a' + a - 1), ttl{a});
    if TH.betaPct < max(abs(ylim(ax(a))))
        sharcThreshold(ax(a), 'y', TH.betaPct, '');                % unlabelled: named in the caption
        sharcThreshold(ax(a), 'y', -TH.betaPct, '');
    else
        text(ax(a), xx(1), yl(2) + 0.6 * pad, sprintf(' \\pm%g%% lies outside the axes', TH.betaPct), ...
            'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.axis, 'VerticalAlignment', 'middle');
    end
end
h = [sharcSeries(ax(end), NaN, NaN, 1), sharcSeries(ax(end), NaN, NaN, 2), ...
    sharcSeries(ax(end), NaN, NaN, 1, 'LineStyle', 'none', 'MarkerSize', 9, 'MarkerFaceColor', 'none'), ...
    sharcSeries(ax(end), NaN, NaN, 2, 'LineStyle', 'none', 'MarkerSize', 9, 'MarkerFaceColor', 'none')];
sharcSharedLegend(fig, ax, h, {'\delta\beta/\beta, linear', '\deltaR/R, linear', '\delta\beta/\beta, nonlinear', ...
    '\deltaR/R, nonlinear'}, 4);
files = [files, saveSharcFigure(fig, run, 3, 'parameter_bias')];

%% ---- printed check ---------------------------------------------------------------------------------------------------
fprintf('\n  V1 figures: ladder %s (selected %s, rule d < %g)\n', strjoin(arrayfun(@(l) lab(l.tr), Lv, ...
    'UniformOutput', false), ', '), lab(ladder.truth), ladder.rule);
for i = 1:nF
    fprintf('    %-24s d_W %.3f  d_sys %.3f  d_tot %.3f  %s -> %s\n', tr.floes{i}, dW(i), dSys(i), dTot(i), ...
        clsM{i}, rowT(i).class);
end
for k = 1:numel(tR)
    fprintf('    radial %s %-13s %-4g at %+.4f mm\n', ternary(tR(k).side < 0, '-', '+'), tR(k).criterion, ...
        tR(k).threshold, 1000 * tR(k).x);
end
fprintf('    nonlinear (%s): dbeta/beta %+.2f%% (linear %+.2f%%), d gap %.2f\n', nl.caseName, nl.betaPct, ...
    nl.linBetaPct, nl.dGap);
end

%% ================================================================================================
function ladder = newestLadder(dDir, truth)
Lst = dir(fullfile(dDir, 'truncationLadder_*.mat'));
[~, o] = sort([Lst.datenum], 'descend');
for k = o
    f = fullfile(Lst(k).folder, Lst(k).name);
    w = whos('-file', f);
    if any(strcmp({w.name}, 'ladder'))
        L0 = load(f, 'ladder');
        if isequal(L0.ladder.truth(:), truth(:))
            ladder = L0.ladder;
            fprintf('  stage 0 ladder from %s\n', f);
            return
        end
    end
end
error('plotV1ModelMismatch:ladder', 'No complete stage 0 ladder with truth [%d %d %d] in %s', truth, dDir);
end

function T = readRunCsv(file)
% minimal CSV reader (no quoted fields expected in V1_truncation.csv): struct of cellstr columns, [] if absent
T = [];
fid = fopen(file, 'r');
if fid < 0, fprintf('  (no %s: d_W and matched class per floe omitted)\n', file); return; end
hdr = strsplit(strtrim(fgetl(fid)), ',');
rows = {};
ln = fgetl(fid);
while ischar(ln)
    if ~isempty(strtrim(ln)), rows(end + 1, :) = strsplit(strtrim(ln), ','); end %#ok<AGROW>
    ln = fgetl(fid);
end
fclose(fid);
for j = 1:numel(hdr), T.(hdr{j}) = rows(:, j); end
end

function s = floeLabel(name)
switch name
    case 'basin_p0', s = 'p_0';
    case 'ring_R_x2', s = 'R\times2';
    case 'corner_beta_x0.5_R_x2', s = '\beta/2, R\times2';
    case 'corner_beta_x2_R_x2', s = '2\beta, R\times2';
    otherwise, s = strrep(name, '_', ' ');
end
end

function s = classShort(c)
switch lower(strtrim(c))
    case 'identifiable', s = 'I';
    case 'marginal', s = 'M';
    case {'not identifiable', 'notidentifiable'}, s = 'N';
    otherwise, s = c;
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end