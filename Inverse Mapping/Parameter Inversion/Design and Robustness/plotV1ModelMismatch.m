
function files = plotV1ModelMismatch(results, cfg, run, ladder)
%PLOTV1MODELMISMATCH  Thesis-ready V1 plots from saved results.
% No EMM solves or nonlinear inversions are rerun.
%
% Usage:
%   plotV1ModelMismatch()
%   plotV1ModelMismatch('...\V1_<stamp>.mat')
%   plotV1ModelMismatch(results, cfg, run, ladder)
%
% Figure 1:
%   (a) Truncation convergence at p0.
%   (b) Welch, model-mismatch and combined bias across floes.
%
% Figure 2:
%   (a) Water-depth sensitivity.
%   (b) Common radial sensor-position sensitivity.
%
% Figure 3:
%   (a) Depth-induced parameter bias.
%   (b) Radial-position-induced parameter bias.
%
% Output directories and naming are unchanged.
% The existing saveSharcFigure calls are retained.

%% ---- Load saved results ------------------------------------------------
if nargin < 1 || ischar(results) || isstring(results)

    P = setupStudyPaths();

    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'V1', 'V1_20*.mat'));

        if isempty(Lst)
            error('plotV1ModelMismatch:none', ...
                'No V1 results in %s', fullfile(P.results, 'V1'));
        end

        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end

    M0 = load(matFile);
    results = M0.results;
    cfg = M0.cfg;

    run = struct('id', 'V1', 'dir', fileparts(matFile));

    fprintf('  Redrawing V1 figures from %s\n', matFile);
end

if nargin < 4 || isempty(ladder)
    ladder = newestLadder( ...
        fullfile(run.dir, 'diagnostics'), ...
        results.truncation.truth);
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


%% =======================================================================
% FIGURE 1: TRUNCATION CONVERGENCE AND MODEL-MISMATCH BIAS
% ========================================================================
%
% Changes:
% - Simplified convergence ladder.
% - Removed crowded supplementary markers and point labels.
% - Logarithmic bias axis so all three bias measures remain visible.
% - Removed overlapping per-floe classification annotations.
% - Compact legends.
%
% The supplementary large-kR and P convergence checks remain
% in the saved diagnostic results but are not plotted here.

[fig, ax] = sharcFigure('full', [1 2], 9.3, ...
    'hgap', 1.8, ...
    'left', 1.45, ...
    'right', 0.65, ...
    'top', 0.70, ...
    'bottom', 1.35);


%% Panel (a): truncation convergence

Lv = ladder.levels;

M = arrayfun(@(l) l.tr(1), Lv);
dStep = [Lv.d];

ok = isfinite(dStep);
Mt = ladder.truth(1);

sharcSeries(ax(1), M(ok), dStep(ok), 1);

xlim(ax(1), [min(M(ok))-9, Mt+10]);
ylim(ax(1), [0, 1.3*max(dStep(ok))]);

sharcTicks(ax(1), 'x', M(ok));

xlabel(ax(1), 'Vertical truncation order M');
ylabel(ax(1), 'Successive-step distance d');

sharcPanel(ax(1), 'a', 'Truncation convergence');

% Registered convergence threshold.
sharcThreshold(ax(1), 'y', ladder.rule, '', 'tick');

% State explicitly that the selected reference is unconverged.
lastStep = dStep(find(ok, 1, 'last'));

text(ax(1), 0.97, 0.89, ...
    sprintf('Final step: d = %.3f > %.1f', ...
        lastStep, ladder.rule), ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontName', S.font, ...
    'FontSize', S.smallFontSize, ...
    'Color', S.axis);


%% Panel (b): truncation mismatch across four floes

x = 1:nF;

dSys = [rowT.dSys];
dTot = [rowT.dTot];

dW = NaN(1, nF);
clsM = repmat({''}, 1, nF);

% Retrieve matched-model Welch bias from the run CSV.
if ~isempty(csv)
    for i = 1:nF
        j = find(strcmp(csv.floe, tr.floes{i}), 1);

        if ~isempty(j)
            dW(i) = str2double(csv.d_W{j});
            clsM{i} = csv.class_matched{j};
        end
    end
end

% Fallback for the nominal floe.
if strcmp(tr.floes{1}, 'basin_p0') && isnan(dW(1))
    dW(1) = results.nominal.dW;
    clsM{1} = results.nominal.class;
end

% Slight horizontal offsets distinguish the three quantities.
dx = 0.16;

hW = sharcSeries(ax(2), ...
    x-dx, dW, 6, 'LineStyle', 'none');

hSys = sharcSeries(ax(2), ...
    x, dSys, 7, 'LineStyle', 'none');

hTot = sharcSeries(ax(2), ...
    x+dx, dTot, 8, 'LineStyle', 'none');

% Logarithmic axis exposes the matched-model Welch bias.
set(ax(2), 'YScale', 'log');

xlim(ax(2), [0.50, nF+0.50]);
ylim(ax(2), [0.20, 55]);

set(ax(2), ...
    'YTick', [0.2 0.5 1 3 10 30], ...
    'XTick', x, ...
    'XTickLabel', { ...
        'p_0', ...
        '2R_0', ...
        '0.5\beta_0, 2R_0', ...
        '2\beta_0, 2R_0'}, ...
    'XTickLabelRotation', 0);

xlabel(ax(2), 'Test floe');
ylabel(ax(2), 'Fisher distance d (log scale)');

sharcPanel(ax(2), 'b', 'Model-mismatch significance');

% Identifiable and Marginal bias limits.
for thresh = TH.dSys(:)'
    yline(ax(2), thresh, '--', ...
        'Color', [0.52 0.52 0.52], ...
        'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end

legend(ax(2), [hW hSys hTot], ...
    {'Matched d_W', ...
     'Mismatch d_{sys}', ...
     'Combined d_{tot}'}, ...
    'Location', 'northwest', ...
    'Box', 'off', ...
    'FontName', S.font, ...
    'FontSize', S.smallFontSize);

% Original export destination and filename retained.
files = [files, saveSharcFigure(fig, run, 1, 'truncation')];


%% =======================================================================
% FIGURE 2: DEPTH AND RADIAL-POSITION BIAS SIGNIFICANCE
% ========================================================================

[fig, ax] = sharcFigure('full', [1 2], 9.3, ...
    'hgap', 1.8, ...
    'left', 1.45, ...
    'right', 0.65, ...
    'top', 0.70, ...
    'bottom', 1.30);


%% Panel (a): water-depth sensitivity

dep = mm(strcmp({mm.kind}, 'depth'));

[e, o] = sort([0, 100*[dep.value]]);

ys = [0, [dep.dSys]];
yt = [results.nominal.dW, [dep.dTot]];

h1 = sharcSeries(ax(1), e, ys(o), 7);
h2 = sharcSeries(ax(1), e, yt(o), 8);

xlim(ax(1), [-5.6, 5.6]);
ylim(ax(1), [0, 1.25]);

sharcTicks(ax(1), 'x', [-5 -2 0 2 5]);

xlabel(ax(1), 'Water-depth error \epsilon_H (%)');
ylabel(ax(1), 'Fisher distance d');

sharcPanel(ax(1), 'a', 'Water-depth error');

for thresh = TH.dSys(:)'
    yline(ax(1), thresh, '--', ...
        'Color', [0.52 0.52 0.52], ...
        'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end

legend(ax(1), [h1 h2], ...
    {'Isolated d_{sys}', 'Combined d_{tot}'}, ...
    'Location', 'southwest', ...
    'Box', 'off', ...
    'FontName', S.font, ...
    'FontSize', S.smallFontSize);


%% Panel (b): radial sensor-position sensitivity

rad = mm(strcmp({mm.kind}, 'radial'));

tol = results.tolerances;

tR = tol( ...
    strcmp({tol.kind}, 'radial') & ...
    [tol.crossed] & ...
    isfinite([tol.x]));

isD = strcmp({tR.criterion}, 'd_sys');
isB = ~isD;

hSign = gobjects(1, 2);

% Plot inward and outward shifts separately.
for ii = 1:2

    signs = [-1 1];
    sd = signs(ii);

    slot = 3.5 + sd/2;

    rs = rad(sign([rad.value]) == sd);
    td = tR(isD & [tR.side] == sd);

    xx = 1000*abs([[td.x], rs.value]);
    yy = [[td.threshold], rs.dSys];

    [xx, o] = sort(xx);

    hSign(ii) = sharcSeries( ...
        ax(2), xx, yy(o), slot);
end

% Unit-slope reference through the small-error crossings.
dRef = tR(isD);

if ~isempty(dRef)

    k = mean( ...
        [dRef.threshold] ./ ...
        (1000*abs([dRef.x])));

    xr = [0.008 30];

    sharcSeries(ax(2), xr, k*xr, 7, ...
        'Marker', 'none', ...
        'LineStyle', ':', ...
        'HandleVisibility', 'off');
end

set(ax(2), ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'XTick', [0.01 0.02 0.05 0.1 0.5 1 5 20], ...
    'YTick', [0.5 1 10 100 1000]);

% Include the 8.9 and 17.9 micrometre crossings.
xlim(ax(2), [0.006, 35]);

ylim(ax(2), ...
    [0.30, max(1500, 1.3*max([rad.dSys]))]);

xlabel(ax(2), 'Common radial error |\Delta r| (mm)');
ylabel(ax(2), 'Isolated d_{sys}');

sharcPanel(ax(2), 'b', 'Sensor-radius error');

for thresh = TH.dSys(:)'
    yline(ax(2), thresh, '--', ...
        'Color', [0.52 0.52 0.52], ...
        'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end

% Mark the physical 1% stiffness-error crossing.
bx = 1000*abs([tR(isB).x]);

if ~isempty(bx)
    xline(ax(2), mean(bx), ':', ...
        'Color', [0.35 0.35 0.35], ...
        'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
end

legend(ax(2), hSign, ...
    {'Inward shift', 'Outward shift'}, ...
    'Location', 'northwest', ...
    'Box', 'off', ...
    'FontName', S.font, ...
    'FontSize', S.smallFontSize);

% Original export destination and filename retained.
files = [files, saveSharcFigure(fig, run, 2, 'significance')];


%% =======================================================================
% FIGURE 3: PARAMETER BIAS AND NONLINEAR CHECK
% ========================================================================

[fig, ax] = sharcFigure('full', [1 2], 9.3, ...
    'hgap', 1.8, ...
    'left', 1.45, ...
    'right', 0.65, ...
    'top', 0.70, ...
    'bottom', 1.30);

nl = results.nonlinear;

[nlKind, nlX] = strtok(nl.caseName);
nlX = str2double(nlX);

grp = {dep, rad};
sc = [100 1000];

xlab = { ...
    'Water-depth error \epsilon_H (%)', ...
    'Radial error \Delta r (mm)'};

ttl = { ...
    'Depth-induced parameter bias', ...
    'Radial position-induced bias'};

kinds = {'depth', 'radial'};

hBeta = [];
hRadius = [];
hNonB = [];
hNonR = [];

for a = 1:2

    g = grp{a};

    [xx, o] = sort([0, sc(a)*[g.value]]);

    yb = [0, [g.betaPct]];
    yr = [0, [g.RPct]];

    hBeta = sharcSeries(ax(a), xx, yb(o), 1);
    hRadius = sharcSeries(ax(a), xx, yr(o), 2);

    yy = [yb yr];

    % Overlay the nonlinear result where available.
    if strcmp(nlKind, kinds{a})

        hNonB = sharcSeries( ...
            ax(a), sc(a)*nlX, nl.betaPct, 1, ...
            'LineStyle', 'none', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', 'none');

        hNonR = sharcSeries( ...
            ax(a), sc(a)*nlX, nl.RPct, 2, ...
            'LineStyle', 'none', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', 'none');

        yy = [yy nl.betaPct nl.RPct];
    end

    pad = max(0.10*(max(yy)-min(yy)), 0.002);

    xlim(ax(a), ...
        [min(xx), max(xx)] + ...
        [-1 1]*0.08*(max(xx)-min(xx)));

    ylim(ax(a), ...
        [min(yy)-pad, max(yy)+pad]);

    sharcTicks(ax(a), 'x', xx);

    xlabel(ax(a), xlab{a});

    if a == 1
        ylabel(ax(a), 'Parameter bias (%)');
    end

    sharcPanel(ax(a), char('a'+a-1), ttl{a});
end

% Compact legend in the second panel.
if ~isempty(hNonB) && ~isempty(hNonR)

    legend(ax(2), ...
        [hBeta hRadius hNonB hNonR], ...
        {'\delta\beta/\beta (linear)', ...
         '\delta R/R (linear)', ...
         '\delta\beta/\beta (nonlinear)', ...
         '\delta R/R (nonlinear)'}, ...
        'Location', 'northeast', ...
        'Box', 'off', ...
        'FontName', S.font, ...
        'FontSize', S.smallFontSize);

else

    legend(ax(2), [hBeta hRadius], ...
        {'\delta\beta/\beta', '\delta R/R'}, ...
        'Location', 'northeast', ...
        'Box', 'off');
end

% Original export destination and filename retained.
files = [files, saveSharcFigure(fig, run, 3, 'parameter_bias')];


%% ---- Diagnostic printout ----------------------------------------------

fprintf('\n  V1 figures: ladder %s (selected %s, rule d < %g)\n', ...
    strjoin(arrayfun(@(l) lab(l.tr), Lv, ...
    'UniformOutput', false), ', '), ...
    lab(ladder.truth), ladder.rule);

for i = 1:nF
    fprintf(['    %-24s d_W %.3f  d_sys %.3f  ' ...
             'd_tot %.3f  %s -> %s\n'], ...
        tr.floes{i}, dW(i), dSys(i), dTot(i), ...
        clsM{i}, rowT(i).class);
end

for k = 1:numel(tR)
    fprintf('    radial %s %-13s %-4g at %+.4f mm\n', ...
        ternary(tR(k).side < 0, '-', '+'), ...
        tR(k).criterion, ...
        tR(k).threshold, ...
        1000*tR(k).x);
end

fprintf(['    nonlinear (%s): dbeta/beta %+.2f%% ' ...
         '(linear %+.2f%%), d gap %.2f\n'], ...
    nl.caseName, nl.betaPct, nl.linBetaPct, nl.dGap);

end


%% =======================================================================
% HELPER FUNCTIONS
% ========================================================================

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

            fprintf('  Stage 0 ladder from %s\n', f);

            return
        end
    end
end

error('plotV1ModelMismatch:ladder', ...
    'No complete Stage 0 ladder with truth [%d %d %d] in %s', ...
    truth, dDir);

end


function T = readRunCsv(file)
% Minimal CSV reader for V1_truncation.csv.
% No quoted fields are expected.

T = [];

fid = fopen(file, 'r');

if fid < 0
    fprintf('  (No %s: d_W and matched classes omitted)\n', file);
    return;
end

hdr = strsplit(strtrim(fgetl(fid)), ',');

rows = {};

ln = fgetl(fid);

while ischar(ln)

    if ~isempty(strtrim(ln))
        rows(end+1, :) = strsplit(strtrim(ln), ','); %#ok<AGROW>
    end

    ln = fgetl(fid);
end

fclose(fid);

for j = 1:numel(hdr)
    T.(hdr{j}) = rows(:, j);
end

end


function s = ternary(c, a, b)

if c
    s = a;
else
    s = b;
end

end
