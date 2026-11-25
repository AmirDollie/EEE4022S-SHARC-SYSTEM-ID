function demoM1Live(varargin)
%DEMOM1LIVE  Live demonstration of the output-only inverse at p0: M1 in miniature.
%
%   demoM1Live                     % R1_p0_lsm6dsv16x, the measured-noise operating point
%   demoM1Live('Case', 'D1_ns2')   % any M1 case name
%   demoM1Live('M1File', 'C:\...\M1_20261009_143730.mat')   % stored run for the nonlinear replay
%
%   WHAT IT SHOWS. Each click generates one independent 109-minute digital-twin record (same seed schedule as
%   the thesis M1 run, so record 1 is the record M1 solved nonlinearly), estimates the noise-corrected Welch
%   transmissibilities and computes the linearised parameter estimate with the analytic J, Sigma_k and F,
%   exactly as M1 did (m1Tools.estimateRecord, g3EstimatorCore unchanged). The estimate is added to the
%   parameter plane over the Fisher ellipses that were drawn before any record existed.
%     top left      control panel (buttons)
%     left          floe, the 33 candidate positions (faint), the layout, reference and incident direction
%     centre        parameter plane: estimates, predicted 68/95% ellipses about truth + predicted Welch bias
%     right         legend and running check: scatter / prediction, 95% coverage, mean offset vs predicted bias
%     bottom left   this record's transmissibility (markers) against the EMM model (line), magnitude and
%                   phase (phase unwrapped on the model's branch)
%     bottom right  replay of the STORED nonlinear Gauss-Newton solves of record 1 (from beta x 0.8 and x 1.2):
%                   Fisher distance to the linearised estimate per iteration. A live nonlinear solve takes
%                   10 to 15 min, so it is replayed from the M1 results file, not recomputed.
%
%   BUTTONS. Next record | +10 records | Replay nonlinear | Reset (clears records and the replay).
%   SETUP COST. evaluateScenario on the cached field: seconds. Run it before visitors arrive.
%
%   Place in Inverse Mapping/Parameter Inversion/Design and Robustness/ (uses setupStudyPaths,
%   defineStudyScenarios, resolveM1Scenarios, m1Tools). Changes no file and writes nothing.

ip = inputParser;
ip.addParameter('Case', 'R1_p0_lsm6dsv16x', @ischar);
ip.addParameter('M1File', '', @ischar);
ip.addParameter('Sensor', [], @(x) isempty(x) || isscalar(x));   % which non-reference sensor to plot (default: last)
ip.parse(varargin{:});
opt = ip.Results;

%% ---- setup (do this before the demo starts) ------------------------------------------------------------------------
fprintf('Setting up %s ...\n', opt.Case);
tSetup = tic;
P = setupStudyPaths();
cfg = defineStudyScenarios();
T = m1Tools();
cases = resolveM1Scenarios(cfg, P);
qc = find(strcmp({cases.name}, opt.Case), 1);
if isempty(qc), error('demoM1Live:case', 'Unknown case %s. Options: %s', opt.Case, strjoin({cases.name}, ', ')); end
ctx = T.context(cases(qc), cfg);
t = ctx.tgt;
fprintf('  analytic target: sigma [%.4f%% %.4f%%], d_W %.3f, %s (%.1f s)\n', 100 * t.sigma, t.dW, t.class, toc(tSetup));

stored = loadStored(opt.M1File, P, opt.Case);
if isempty(stored)
    fprintf('  no stored M1 run found: the nonlinear replay is disabled (pass ''M1File'').\n');
else
    fprintf('  stored M1 run: %s (%d records, %d nonlinear solves)\n', stored.file, stored.st.nDone, numel(stored.st.nl));
end

m = t.pFeat / 2;                                   % non-reference sensors
jPlot = opt.Sensor; if isempty(jPlot), jPlot = m; end
f0 = reshape(t.f0, t.pFeat, []);                   % model features on all bins, [Re; Im] per frequency
Tmod = f0(jPlot, :) + 1i * f0(jPlot + m, :);
w = t.omega(:).';
chi95 = -2 * log(0.05); chi68 = -2 * log(1 - 0.6827);

%% ---- figure ---------------------------------------------------------------------------------------------------------
navy = [0.07 0.29 0.48]; sky = [0.36 0.61 0.78]; ice = [0.92 0.95 0.98]; orange = [0.92 0.41 0.20]; grey = [0.45 0.45 0.45];
fig = figure('Name', 'Output-only inverse: live demonstration', 'Color', 'w', 'Units', 'normalized', ...
    'Position', [0.02 0.05 0.96 0.86], 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', ...
    'DefaultAxesFontSize', 10.5, 'DefaultTextFontSize', 10.5);
annotation(fig, 'textbox', [0.02 0.945 0.96 0.05], 'String', ...
    'Recovering floe rigidity \beta and radius R from accelerations alone', 'EdgeColor', 'none', ...
    'FontSize', 18, 'FontWeight', 'bold', 'Color', navy, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
annotation(fig, 'textbox', [0.02 0.910 0.96 0.032], 'String', sprintf(['Case %s:  %.0f-minute records,  ' ...
    'Welch segment L = %d,  LSM6DSV16X-class sensor noise.  Each record is new; the analytic prediction is fixed.'], ...
    strrep(t.name, '_', '\_'), t.N * t.dt / 60, t.L), 'EdgeColor', 'none', 'FontSize', 11.5, 'Color', grey, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% ---- controls (top left) ----
pnl = uipanel(fig, 'Units', 'normalized', 'Position', [0.02 0.725 0.25 0.150], 'BackgroundColor', ice, ...
    'BorderType', 'line', 'HighlightColor', sky, 'Title', 'Controls', 'FontSize', 11, 'FontWeight', 'bold', ...
    'ForegroundColor', navy);
mkBtn(pnl, 'Next record', [0.04 0.55 0.44 0.38], @(~, ~) addRecords(1), 'Generate one new record and invert it', true);
mkBtn(pnl, '+10 records', [0.52 0.55 0.44 0.38], @(~, ~) addRecords(10), 'Generate and invert ten records', false);
mkBtn(pnl, 'Replay nonlinear', [0.04 0.08 0.44 0.38], @(~, ~) replayNL(), ...
    'Replay the stored full nonlinear solves of record 1', false);
mkBtn(pnl, 'Reset', [0.52 0.08 0.44 0.38], @(~, ~) resetDemo(), 'Clear all records', false);

% ---- floe and sensors (left) ----
axF = axes(fig, 'Position', [0.02 0.385 0.25 0.31]); hold(axF, 'on'); axis(axF, 'equal'); axis(axF, 'off');
R = t.p(3); th = linspace(0, 2 * pi, 200);
fill(axF, R * cos(th), R * sin(th), ice, 'EdgeColor', sky, 'LineWidth', 1.8);
G = cfg.layout.gridFrac;
plot(axF, G(:, 1) * R .* cos(G(:, 2)), G(:, 1) * R .* sin(G(:, 2)), '.', 'Color', [0.75 0.80 0.85], 'MarkerSize', 8);
Sp = t.sensors(t.layout, :);                       % [r theta], nondimensional r, theta in rad
S = [Sp(:, 1) .* cos(Sp(:, 2)), Sp(:, 1) .* sin(Sp(:, 2))];
names = sensorNames(t, cfg);
others = setdiff(1:numel(t.layout), t.ref);
plot(axF, S(:, 1), S(:, 2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', navy, 'MarkerEdgeColor', 'w');
plot(axF, S(t.ref, 1), S(t.ref, 2), 's', 'MarkerSize', 16, 'LineWidth', 2.2, 'MarkerEdgeColor', orange);
plot(axF, S(others(jPlot), 1), S(others(jPlot), 2), 'o', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', sky);
[~, ord] = sort(S(:, 1));                          % alternate above / below along the wave axis
for i = 1:numel(ord)
    k = ord(i); up = mod(i, 2) == 1;
    text(axF, S(k, 1), S(k, 2) + (2 * up - 1) * 0.30 * R, names{k}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', navy, 'FontWeight', 'bold', ...
        'BackgroundColor', [1 1 1], 'Margin', 1);
end
quiver(axF, 1.80 * R, 0, -0.45 * R, 0, 0, 'LineWidth', 2.2, 'MaxHeadSize', 1, 'Color', grey);
text(axF, 1.57 * R, 0.22 * R, 'waves', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', grey);
text(axF, -1.10 * R, -1.40 * R, {'\color[rgb]{0.92,0.41,0.20}\bforange square\rm\color{black}: reference sensor', ...
    sprintf('\\color[rgb]{0.36,0.61,0.78}\\bfblue ring\\rm\\color{black}: sensor %s, plotted below', names{others(jPlot)})}, ...
    'FontSize', 9.5, 'VerticalAlignment', 'top');
xlim(axF, [-1.25 1.90] * R); ylim(axF, [-1.95 1.25] * R);
title(axF, 'Floe and sensor layout', 'FontSize', 11.5, 'Color', navy);

% ---- parameter plane (centre) ----
axE = axes(fig, 'Position', [0.33 0.44 0.44 0.435]); hold(axE, 'on'); grid(axE, 'on'); box(axE, 'on');
set(axE, 'GridAlpha', 0.12, 'Layer', 'top');
mu = 100 * t.dthetaW(:);
C = 1e4 * t.C;
[e95x, e95y] = ellipse(mu, C, chi95); [e68x, e68y] = ellipse(mu, C, chi68);
hE95 = fill(axE, e95x, e95y, ice, 'EdgeColor', sky, 'LineWidth', 1.6);
hE68 = plot(axE, e68x, e68y, '--', 'Color', sky, 'LineWidth', 1.3);
hPts = scatter(axE, NaN, NaN, 26, navy, 'filled', 'MarkerFaceAlpha', 0.5);
hLast = plot(axE, NaN, NaN, 'o', 'MarkerSize', 13, 'LineWidth', 2.2, 'Color', orange);
hTrue = plot(axE, 0, 0, 'k+', 'MarkerSize', 16, 'LineWidth', 2.2);
hBias = plot(axE, mu(1), mu(2), 'd', 'MarkerSize', 10, 'MarkerFaceColor', orange, 'MarkerEdgeColor', 'k');
hNL = plot(axE, NaN, NaN, 'p', 'MarkerSize', 18, 'MarkerFaceColor', [1 0.85 0.2], 'MarkerEdgeColor', 'k');
lim = 1.35 * sqrt(chi95 * diag(C)).' + abs(mu.');
xlim(axE, [-1 1] * lim(1)); ylim(axE, [-1 1] * lim(2));
xlabel(axE, 'estimated \beta  relative to the true floe (%)');
ylabel(axE, 'estimated R  relative to the true floe (%)');
title(axE, 'Parameter estimates: one dot per independent record', 'FontSize', 12.5, 'Color', navy);
legend(axE, [hE95 hE68 hPts hLast hTrue hBias hNL], {'predicted 95% region', 'predicted 68% region', ...
    'estimate from one record', 'latest record', 'true floe', 'predicted Welch bias', 'nonlinear solve (stored)'}, ...
    'Position', [0.795 0.665 0.185 0.21], 'FontSize', 9.5, 'Box', 'off');

% ---- running statistics (right) ----
annotation(fig, 'textbox', [0.795 0.625 0.185 0.03], 'String', 'Running check', 'EdgeColor', 'none', ...
    'FontSize', 11.5, 'FontWeight', 'bold', 'Color', navy, 'Margin', 0);
hStat = annotation(fig, 'textbox', [0.795 0.365 0.185 0.255], 'String', '', 'EdgeColor', sky, ...
    'BackgroundColor', ice, 'FontSize', 10, 'VerticalAlignment', 'top', 'Margin', 7, 'Interpreter', 'tex');

% ---- this record's transmissibility (bottom left) ----
phM = rad2deg(unwrap(angle(Tmod)));
axM = axes(fig, 'Position', [0.05 0.075 0.22 0.21]); hold(axM, 'on'); grid(axM, 'on'); box(axM, 'on');
hMagD = plot(axM, NaN, NaN, 'o', 'MarkerSize', 4, 'Color', orange, 'LineWidth', 1);
plot(axM, w, abs(Tmod), '-', 'Color', navy, 'LineWidth', 2);
ylabel(axM, '|T|'); xlabel(axM, '\omega (rad/s)'); xlim(axM, [min(w) max(w)]);
axP = axes(fig, 'Position', [0.32 0.075 0.22 0.21]); hold(axP, 'on'); grid(axP, 'on'); box(axP, 'on');
hPhD = plot(axP, NaN, NaN, 'o', 'MarkerSize', 4, 'Color', orange, 'LineWidth', 1);
plot(axP, w, phM, '-', 'Color', navy, 'LineWidth', 2);
ylabel(axP, 'arg T (deg)'); xlabel(axP, '\omega (rad/s)'); xlim(axP, [min(w) max(w)]);
legend(axP, {'this record', 'EMM model'}, 'Location', 'southwest', 'FontSize', 9, 'Box', 'off');
annotation(fig, 'textbox', [0.05 0.300 0.49 0.03], 'String', sprintf(['Transmissibility of %s relative to %s ' ...
    '(the unknown wave cancels in this ratio)'], names{others(jPlot)}, names{t.ref}), 'EdgeColor', 'none', ...
    'FontSize', 11.5, 'FontWeight', 'bold', 'Color', navy, 'Margin', 0, 'VerticalAlignment', 'middle');

% ---- nonlinear replay (bottom right) ----
axN = axes(fig, 'Position', [0.62 0.075 0.36 0.21]);
resetNLAxes();

D = NaN(2, 0); Q = NaN(1, 0); nRec = 0; tLast = NaN;
updateStats();
fprintf('Ready (%.1f s). Click "Next record".\n', toc(tSetup));

%% ---- callbacks ------------------------------------------------------------------------------------------------------
    function addRecords(n)
        for k = 1:n
            nRec = nRec + 1;
            seed = T.seed(cfg, qc, nRec);
            tr = tic;
            e = T.estimateRecord(ctx, seed);
            tLast = toc(tr);
            D(:, end + 1) = e.dtheta; Q(end + 1) = e.q; %#ok<AGROW>
            Tz = e.z(jPlot, :) + 1i * e.z(jPlot + m, :);
            phD = rad2deg(angle(Tz));
            phD = phD + 360 * round((phM - phD) / 360);          % same branch as the unwrapped model
            set(hMagD, 'XData', w(e.use), 'YData', abs(Tz(e.use)));
            set(hPhD, 'XData', w(e.use), 'YData', phD(e.use));
            set(hPts, 'XData', 100 * D(1, :), 'YData', 100 * D(2, :));
            set(hLast, 'XData', 100 * e.dtheta(1), 'YData', 100 * e.dtheta(2));
            updateStats();
            if nRec == 1 && ~isempty(stored) && stored.st.nDone >= 1
                gap = max(abs(e.dtheta - stored.st.dtheta(:, 1)));
                fprintf('  record 1 vs the thesis M1 run: max difference %.1e (reproduces from its seed)\n', gap);
            end
            drawnow;
        end
    end

    function updateStats()
        if nRec == 0
            str = {'No records yet.', '', 'Press \bfNext record\rm.'};
        else
            str = {sprintf('Records: \\bf%d\\rm   (%.2f s each)', nRec, tLast), ...
                sprintf('Last: \\beta %+.3f%%, R %+.4f%%', 100 * D(1, end), 100 * D(2, end)), ''};
            if nRec >= 5
                r = std(D, 0, 2) ./ t.sigma(:);
                str = [str, {sprintf('Scatter/pred: \\beta %.2f, R %.2f', r(1), r(2)), ...
                    sprintf('In 95%% ellipse: \\bf%.0f%%\\rm', 100 * mean(Q <= chi95)), ...
                    sprintf('Mean \\beta: %+.3f%% (pred %+.3f%%)', 100 * mean(D(1, :)), mu(1))}];
            else
                str = [str, {'(statistics after 5 records)'}];
            end
        end
        set(hStat, 'String', str);
    end

    function resetNLAxes()
        cla(axN); hold(axN, 'on'); grid(axN, 'on'); box(axN, 'on');
        set(axN, 'YScale', 'log', 'YLim', [1e-4 3e3], 'YTick', 10.^(-4:3), 'XLim', [0 5], 'XTick', 0:5, ...
            'YMinorGrid', 'off', 'GridAlpha', 0.12);
        yline(axN, 1, ':', '1 \sigma', 'Color', grey, 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.2);
        xlabel(axN, 'Gauss-Newton iteration'); ylabel(axN, 'distance (\sigma)');
        title(axN, 'Full nonlinear solve of record 1 (stored; about 10 min each live)', 'FontSize', 11.5, 'Color', navy);
    end

    function replayNL()
        if isempty(stored) || isempty(stored.st.nl)
            set(hStat, 'String', {'No stored nonlinear solve', 'was found (pass ''M1File'').'}); return
        end
        resetNLAxes();
        cols = [navy; orange];
        F = stored.t.F; hs = gobjects(0); leg = {};
        for k = 1:numel(stored.st.nl)
            nl = stored.st.nl(k);
            H = nl.history(1:2, :);
            d = max(sqrt(sum((F * (H - nl.thetaLin)) .* (H - nl.thetaLin), 1)), 1e-4);
            c = cols(1 + mod(k - 1, 2), :);
            h = plot(axN, NaN, NaN, '-o', 'Color', c, 'LineWidth', 2, 'MarkerFaceColor', c, 'MarkerSize', 6);
            for i = 1:numel(d)
                set(h, 'XData', 0:i - 1, 'YData', d(1:i)); drawnow; pause(0.5);
            end
            hs(end + 1) = h; %#ok<AGROW>
            leg{end + 1} = sprintf('%.1f \\beta start: %.0f \\sigma \\rightarrow %.0e \\sigma', ...
                1 + nl.betaOffset, d(1), d(end)); %#ok<AGROW>
            dth = nl.theta - stored.t.thetaTrue(:);
            set(hNL, 'XData', 100 * dth(1), 'YData', 100 * dth(2));
        end
        legend(axN, hs, leg, 'Location', 'northeast', 'FontSize', 9.5, 'Box', 'off');
    end

    function resetDemo()
        D = NaN(2, 0); Q = NaN(1, 0); nRec = 0; tLast = NaN;
        set([hPts hLast hNL hMagD hPhD], 'XData', NaN, 'YData', NaN);
        resetNLAxes();
        updateStats();
    end
end

%% ================================================================================================
function [x, y] = ellipse(mu, C, c)
th = linspace(0, 2 * pi, 200);
L = chol(C, 'lower');
p = mu(:) + sqrt(c) * L * [cos(th); sin(th)];
x = p(1, :); y = p(2, :);
end

function mkBtn(parent, label, pos, cb, tip, primary)
b = uicontrol(parent, 'Style', 'pushbutton', 'String', label, 'Units', 'normalized', 'Position', pos, ...
    'FontSize', 12, 'Callback', cb, 'TooltipString', tip, 'BackgroundColor', [1 1 1]);
if primary, set(b, 'FontWeight', 'bold', 'ForegroundColor', [0.07 0.29 0.48]); end
end

function names = sensorNames(t, cfg)
% grid names (s1..s33) where the layout indexes the 33-point grid, else the Level 2C names
n = size(t.sensors, 1);
if n == size(cfg.layout.gridFrac, 1)
    idx = t.layout;
elseif n == numel(cfg.layout.level2C.idx)
    idx = cfg.layout.level2C.idx(t.layout);
else
    idx = t.layout;
end
names = arrayfun(@(i) sprintf('s%d', i), idx, 'UniformOutput', false);
end

function s = loadStored(file, P, caseName)
% the stored M1 run: the given file, else the newest Results/M1/*/M1_*.mat
s = [];
if isempty(file)
    hits = dir(fullfile(P.results, 'M1', '**', 'M1_*.mat'));
    hits = hits(~contains({hits.name}, 'ckpt'));
    if isempty(hits), return; end
    [~, i] = max([hits.datenum]);
    file = fullfile(hits(i).folder, hits(i).name);
end
if ~exist(file, 'file'), return; end
L = load(file, 'results');
if ~isfield(L, 'results'), return; end
R = L.results;
if isfield(R, 'names'), names = cellstr(R.names); else, names = {R.cases.name}; end
i = find(strcmp(names, caseName), 1);
if isempty(i), return; end
s = struct('file', file, 't', R.targets{i}, 'st', R.state{i});
end