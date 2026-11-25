function floeFlexDemo(cacheFile)
%FLOEFLEXDEMO  Interactive animation of a floe flexing under incident waves, from the EMM forward model.
%
%   floeFlexDemo                      % uses floeFlexCache.mat next to this file (buildFloeFlexCache first)
%
%   LEFT   3D floe surface, coloured by bending curvature (Laplacian of the deflection: zero for rigid heave
%          and pitch, so colour = bending only), vertical exaggeration shown; sensors s2, s6, s14, s26 marked.
%          Cross-section along the wave axis: incident wave (grey) and floe deflection (navy), true scale.
%          Sensor traces: vertical acceleration at the reference s26 and at s14, as the IMUs would record it.
%   RIGHT  Floe: beta and R sliders (snap to the cached grid) and physical presets.
%          Sensor layout: Fisher-optimal pair, Gramian-optimal pair, a mirrored (no-information) pair,
%          Level 2C and the direction-robust cross, with sigma_lnbeta, sigma_lnR, d_W and class from the
%          D1 engine at p0 (buildFloeFlexLayouts), and a bar chart comparing them.
%          Sea: monochromatic (omega) or JONSWAP (omega_p, gamma_J), with the spectrum.
%          Information: physical size, flexural length, R/l_f, k R, regime, peak deflection.
%          Play / pause and speed.
%
%   Nothing is solved live: every floe was solved by buildFloeFlexCache at production truncation; a sea is
%   only a re-weighting of the cached frequency responses (linear superposition), so changing it is instant.
%   Run buildFloeFlexCache and buildFloeFlexLayouts once beforehand; on the day only this file runs.
%   Lives in Forward Model/Flex Demo/.

here = fileparts(mfilename('fullpath'));
if nargin < 1, cacheFile = fullfile(here, 'floeFlexCache.mat'); end
if ~exist(cacheFile, 'file'), error('floeFlexDemo:cache', 'No cache: run buildFloeFlexCache first.'); end
S = load(cacheFile, 'cache'); cache = S.cache;
if any(cellfun(@isempty, cache.pat)), error('floeFlexDemo:cache', 'Cache incomplete: rerun buildFloeFlexCache.'); end
layFile = fullfile(here, 'floeFlexLayouts.mat');
hasLay = exist(layFile, 'file') == 2;
if hasLay
    S2 = load(layFile, 'layouts'); LL = S2.layouts.L; gridFrac = S2.layouts.gridFrac;
else                                                  % sensors still shown, no uncertainty panel
    warning('floeFlexDemo:layouts', 'No floeFlexLayouts.mat: run buildFloeFlexLayouts for the layout study.');
    [TTg, RRg0] = meshgrid((0:7) * pi / 4, [0.3 0.5 0.7 0.9]);
    gridFrac = [0 0; reshape(RRg0.', [], 1), reshape(TTg.', [], 1)];
    LL = struct('name', 'Level 2C (thesis reference)', 'short', 'Level 2C', 'idx', [2 6 14 26], 'ref', 26, ...
        'sigma', [NaN NaN], 'dW', NaN, 'class', '', 'singular', false, 'C', NaN(2), 'A', NaN, ...
        'note', 'Run buildFloeFlexLayouts for the uncertainties.', 'recordFactor', NaN);
end
iLay0 = find(startsWith({LL.name}, 'Level 2C'), 1); if isempty(iLay0), iLay0 = 1; end
drDir = fullfile(here, '..', '..', 'Inverse Mapping', 'Parameter Inversion', 'Design and Robustness');

% ---- constants and colours ----
navy = [0.07 0.29 0.48]; sky = [0.36 0.61 0.78]; ice = [0.92 0.95 0.98]; orange = [0.92 0.41 0.20];
grey = [0.45 0.45 0.45];
H = cache.H; g = cache.g;
wN = cache.omega;                                      % cached nodes
wF = wN(1):0.02:wN(end);                               % fine grid for JONSWAP superposition (repeat 314 s)
rng(7); epsF = 2 * pi * rand(1, numel(wF));            % fixed random phases
Hs = 0.05; aMono = 0.025;                              % reference sea H_s; monochromatic amplitude (m)
iP0 = find(strcmp({cache.floes.kind}, 'grid') & [cache.floes.iB] == 3 & [cache.floes.iR] == 3, 1);
traceWin = 12;                                         % s of sensor trace shown

% ---- state ----
st = struct('busy', false, 'warned', false, 'floe', iP0, 'layout', iLay0, 'mode', 'JONSWAP', 'w', 5, 'wp', 5, 'gJ', 3.3, 't', 0, 'speed', 1, ...
    'playing', true, 'lastTic', tic);
cur = struct();                                        % matrices for the current floe and sea
trace = struct('t', [], 'a', zeros(2, 0));

%% ---- window -----------------------------------------------------------------------------------------------------
fig = uifigure('Name', 'Floe flexing under incident waves', 'Color', 'w', 'Units', 'normalized', ...
    'Position', [0.02 0.05 0.96 0.86]);
main = uigridlayout(fig, [2 2], 'ColumnWidth', {'1.55x', '1x'}, 'RowHeight', {34, '1x'}, ...
    'Padding', [12 6 12 6], 'ColumnSpacing', 14, 'RowSpacing', 4, 'BackgroundColor', 'w');
hdr = uilabel(main, 'Text', 'How a floe bends under waves, and why that decides whether its stiffness can be measured', ...
    'FontSize', 17, 'FontWeight', 'bold', 'FontColor', navy);
hdr.Layout.Row = 1; hdr.Layout.Column = [1 2];

% Greek letters as characters (plain labels, one font throughout)
cB = char(946); cG = char(947); cW = char(969); cL = char(8467); c0 = char(8320);

% left column: [3D floe | information], cross-section, sensor traces
left = uigridlayout(main, [3 1], 'RowHeight', {'1x', 165, 135}, 'Padding', [0 0 0 0], ...
    'RowSpacing', 4, 'BackgroundColor', 'w');
left.Layout.Row = 2; left.Layout.Column = 1;
top = uigridlayout(left, [1 2], 'ColumnWidth', {'1.65x', '1x'}, 'Padding', [0 0 0 0], 'ColumnSpacing', 10, ...
    'BackgroundColor', 'w');
top.Layout.Row = 1;
ax3 = uiaxes(top); ax3.Layout.Column = 1;
axX = uiaxes(left); axX.Layout.Row = 2;
axT = uiaxes(left); axT.Layout.Row = 3;

% right column: floe, sea, playback
right = uigridlayout(main, [4 1], 'RowHeight', {'fit', 'fit', 'fit', 'fit'}, 'Padding', [0 0 0 0], ...
    'RowSpacing', 5, 'BackgroundColor', 'w');
right.Layout.Row = 2; right.Layout.Column = 2;

% floe panel
pF = uipanel(right, 'Title', 'Floe', 'FontWeight', 'bold', 'FontSize', 13, 'ForegroundColor', navy, ...
    'BackgroundColor', ice);
gF = uigridlayout(pF, [4 2], 'ColumnWidth', {100, '1x'}, 'RowHeight', {38, 38, 'fit', 'fit'}, ...
    'BackgroundColor', ice, 'RowSpacing', 3, 'Padding', [8 4 8 4]);
uilabel(gF, 'Text', ['Rigidity  ' cB ' / ' cB c0], 'FontSize', 13);
sB = uislider(gF, 'Limits', [1 5], 'MajorTicks', 1:5, 'MajorTickLabels', {'0.25', '0.5', '1', '2', '4'}, ...
    'MinorTicks', [], 'Value', 3, 'FontSize', 11);
uilabel(gF, 'Text', ['Radius  R / R' c0], 'FontSize', 13);
sR = uislider(gF, 'Limits', [1 5], 'MajorTicks', 1:5, 'MajorTickLabels', {'0.5', '0.75', '1', '1.5', '2'}, ...
    'MinorTicks', [], 'Value', 3, 'FontSize', 11);
lp = uilabel(gF, 'Text', 'Presets', 'FontSize', 13); lp.Layout.Row = 3; lp.Layout.Column = 1;
gP = uigridlayout(gF, [1 5], 'Padding', [0 0 0 0], 'RowSpacing', 4, 'ColumnSpacing', 4, 'BackgroundColor', ice, ...
    'RowHeight', {26});
gP.Layout.Row = 3; gP.Layout.Column = 2;
presetNames = {'p0 basin disc', '1 m PP disc, 5 mm', 'SHARC model floe', 'MIZ floe, R = 10 m', 'MIZ floe, R = 50 m'};
presetLabels = {'p0 disc', '1 m PP', 'SHARC', 'MIZ 10 m', 'MIZ 50 m'};
for k = 1:numel(presetNames)
    uibutton(gP, 'Text', presetLabels{k}, 'FontSize', 10, 'ButtonPushedFcn', @(~, ~) choosePreset(presetNames{k}));
end
lFloe = uilabel(gF, 'Text', '', 'FontSize', 11, 'FontWeight', 'bold', 'FontColor', navy);
lFloe.Layout.Row = 4; lFloe.Layout.Column = [1 2];

% sensor layout panel
pL = uipanel(right, 'Title', 'Sensor layout (D1 engine, at p0)', 'FontWeight', 'bold', 'FontSize', 13, ...
    'ForegroundColor', navy, 'BackgroundColor', ice);
gL = uigridlayout(pL, [3 2], 'ColumnWidth', {100, '1x'}, 'RowHeight', {26, 'fit', 'fit'}, ...
    'BackgroundColor', ice, 'RowSpacing', 3, 'Padding', [8 4 8 4]);
uilabel(gL, 'Text', 'Layout', 'FontSize', 13);
dLay = uidropdown(gL, 'Items', {LL.name}, 'ItemsData', 1:numel(LL), 'Value', iLay0, 'FontSize', 12, ...
    'ValueChangedFcn', @(src, ~) selectLayout(src.Value));
lLay = uilabel(gL, 'Text', '', 'FontSize', 11, 'WordWrap', 'on');
lLay.Layout.Row = 2; lLay.Layout.Column = [1 2];
gT = uigridlayout(gL, [numel(LL) + 1, 4], 'ColumnWidth', {'1.5x', '1x', '1x', '1.2x'}, ...
    'RowHeight', repmat({'fit'}, 1, numel(LL) + 1), 'Padding', [0 0 0 0], 'RowSpacing', 0, ...
    'BackgroundColor', 'w');
gT.Layout.Row = 3; gT.Layout.Column = [1 2];
hdrT = {'Layout', [cB ' uncertainty'], 'Record vs Fisher', 'Class'};
for k = 1:4
    uilabel(gT, 'Text', hdrT{k}, 'FontSize', 10, 'FontWeight', 'bold', 'FontColor', grey);
end
tblLab = gobjects(numel(LL), 4);
for k = 1:numel(LL)
    q = LL(k);
    if q.singular || ~isfinite(q.sigma(1))
        v = {q.short, 'no information', [char(8734) ' (none)'], 'Not identifiable'};
    else
        v = {q.short, sprintf('%.3f %%', 100 * q.sigma(1)), sprintf('%.2g%s', q.recordFactor, char(215)), q.class};
    end
    for j = 1:4, tblLab(k, j) = uilabel(gT, 'Text', v{j}, 'FontSize', 11); end
end

% sea panel
pS = uipanel(right, 'Title', 'Incident waves', 'FontWeight', 'bold', 'FontSize', 13, 'ForegroundColor', navy, ...
    'BackgroundColor', ice);
gS = uigridlayout(pS, [3 2], 'ColumnWidth', {100, '1x'}, 'RowHeight', {26, 38, 38}, ...
    'BackgroundColor', ice, 'RowSpacing', 3, 'Padding', [8 4 8 4]);
uilabel(gS, 'Text', 'Sea', 'FontSize', 13);
dMode = uidropdown(gS, 'Items', {'JONSWAP (irregular sea)', 'Monochromatic (one frequency)'}, 'FontSize', 12);
lS1 = uilabel(gS, 'Text', '', 'FontSize', 13, 'WordWrap', 'on');
s1 = uislider(gS, 'FontSize', 11);
lS2 = uilabel(gS, 'Text', '', 'FontSize', 13, 'WordWrap', 'on');
s2 = uislider(gS, 'FontSize', 11);

% information panel
side = uigridlayout(top, [2 1], 'RowHeight', {'fit', '1x'}, 'Padding', [0 0 0 0], 'RowSpacing', 6, ...
    'BackgroundColor', 'w');
side.Layout.Column = 2;
pI = uipanel(side, 'Title', 'What the floe is doing', 'FontWeight', 'bold', 'FontSize', 12, ...
    'ForegroundColor', navy, 'BackgroundColor', ice);
gI = uigridlayout(pI, [5 2], 'ColumnWidth', {92, '1x'}, 'RowHeight', repmat({'fit'}, 1, 5), ...
    'BackgroundColor', ice, 'RowSpacing', 2, 'Padding', [8 4 8 4]);
infoNames = {'Size', ['R / ' cL 'f'], 'Wave', 'Deflection', 'Verdict'};
infoVal = gobjects(1, 5);
pSp = uipanel(side, 'Title', 'Sea spectrum', 'FontWeight', 'bold', 'FontSize', 12, 'ForegroundColor', navy, ...
    'BackgroundColor', 'w');
pSp.Layout.Row = 2;
gSp = uigridlayout(pSp, [1 1], 'Padding', [4 4 8 4], 'BackgroundColor', 'w');
axS = uiaxes(gSp);
for k = 1:5
    uilabel(gI, 'Text', infoNames{k}, 'FontSize', 11, 'FontColor', grey, 'WordWrap', 'on');
    infoVal(k) = uilabel(gI, 'Text', '', 'FontSize', 11, 'WordWrap', 'on');
end
infoVal(5).FontWeight = 'bold';
infoVal(2).FontWeight = 'bold';

% playback
gB = uigridlayout(right, [1 5], 'ColumnWidth', {100, 50, 80, '1x', 170}, 'Padding', [0 0 0 0], ...
    'BackgroundColor', 'w');
bPlay = uibutton(gB, 'Text', 'Pause', 'FontSize', 13, 'FontWeight', 'bold', 'ButtonPushedFcn', @(~, ~) togglePlay());
uilabel(gB, 'Text', 'Speed', 'FontSize', 12, 'HorizontalAlignment', 'right');
dSpeed = uidropdown(gB, 'Items', {'0.25x', '0.5x', '1x'}, 'Value', '1x', 'FontSize', 12, ...
    'ValueChangedFcn', @(src, ~) setSpeed(src.Value));
uilabel(gB, 'Text', '');
uibutton(gB, 'Text', 'Open inversion demo', 'FontSize', 12, 'ButtonPushedFcn', @(~, ~) openM1(), ...
    'Tooltip', 'Every floe solved with the EMM at production truncation [50 10 10]');

%% ---- static plot elements ---------------------------------------------------------------------------------------
cmap = divergingMap(256, navy, orange);
hold(ax3, 'on');
hSurf = surf(ax3, zeros(2), zeros(2), zeros(2), zeros(2), 'EdgeColor', 'none', 'FaceColor', 'interp', ...
    'FaceLighting', 'gouraud', 'AmbientStrength', 0.6);
hRim = plot3(ax3, NaN, NaN, NaN, '-', 'Color', navy, 'LineWidth', 1.5);
hCand = plot3(ax3, NaN, NaN, NaN, '.', 'Color', [0.55 0.6 0.65], 'MarkerSize', 9);
hSens = plot3(ax3, NaN, NaN, NaN, 'o', 'MarkerSize', 7, 'MarkerFaceColor', navy, 'MarkerEdgeColor', 'w');
hRef = plot3(ax3, NaN, NaN, NaN, 's', 'MarkerSize', 11, 'LineWidth', 2, 'MarkerEdgeColor', orange);
hArrow = quiver3(ax3, 0, 0, 0, 0, 0, 0, 0, 'Color', grey, 'LineWidth', 2.5, 'MaxHeadSize', 0.8);
hTxt = text(ax3, 0, 0, 0, 'waves', 'Color', grey, 'FontSize', 12, 'HorizontalAlignment', 'center');
colormap(ax3, cmap); light(ax3, 'Position', [-1 -2 3]);
view(ax3, -32, 28); axis(ax3, 'off');
hTitle3 = title(ax3, '', 'FontSize', 12, 'Color', navy);
for axh = [ax3 axX axT axS]
    axh.Toolbar.Visible = 'off'; disableDefaultInteractivity(axh);
end

hold(axX, 'on'); grid(axX, 'on'); box(axX, 'on');
hInc = plot(axX, NaN, NaN, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.3);
hDefl = plot(axX, NaN, NaN, '-', 'Color', navy, 'LineWidth', 3);
hSensX = plot(axX, NaN, NaN, 'o', 'MarkerSize', 6, 'MarkerFaceColor', navy, 'MarkerEdgeColor', 'w');
xlabel(axX, 'position along the wave direction (m)'); ylabel(axX, 'elevation (mm)');
title(axX, ['Cross-section along the wave axis (true scale):  \color[rgb]{0.55,0.55,0.55}incident wave' ...
    '\color[rgb]{0.07,0.29,0.48}  floe'], 'FontSize', 11, 'Color', navy);

hold(axT, 'on'); grid(axT, 'on'); box(axT, 'on');
hTr1 = plot(axT, NaN, NaN, '-', 'Color', orange, 'LineWidth', 1.4);
hTr2 = plot(axT, NaN, NaN, '--', 'Color', navy, 'LineWidth', 1.4);
xlabel(axT, 'time (s)'); ylabel(axT, 'accel. (m/s^2)');
hTitleT = title(axT, '', 'FontSize', 11, 'Color', navy);

hold(axS, 'on'); box(axS, 'on');
xlabel(axS, '\omega (rad/s)'); ylabel(axS, 'S(\omega)');

%% ---- callbacks wired, first draw -------------------------------------------------------------------------------
sB.ValueChangedFcn = @(src, ~) chooseGrid(round(src.Value), round(sR.Value));
sR.ValueChangedFcn = @(src, ~) chooseGrid(round(sB.Value), round(src.Value));
dMode.ValueChangedFcn = @(src, ~) setMode(src.Value);
s1.ValueChangedFcn = @(~, ~) seaChanged();
s2.ValueChangedFcn = @(~, ~) seaChanged();
setMode(dMode.Value);
selectFloe(iP0);
selectLayout(iLay0);

tmr = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.04, 'BusyMode', 'drop', 'TimerFcn', @(~, ~) frame());
fig.CloseRequestFcn = @(~, ~) closeApp();
st.lastTic = tic;
start(tmr);

%% =================================================================================================================
%  callbacks
%% =================================================================================================================
    function chooseGrid(iB, iR)
        sB.Value = iB; sR.Value = iR;
        k = find(strcmp({cache.floes.kind}, 'grid') & [cache.floes.iB] == iB & [cache.floes.iR] == iR, 1);
        selectFloe(k);
    end

    function choosePreset(name)
        if strcmp(name, 'p0 basin disc'), chooseGrid(3, 3); return; end
        selectFloe(find(strcmp({cache.floes.name}, name), 1));
    end

    function setMode(v)
        st.mode = ternary(startsWith(v, 'JONSWAP'), 'JONSWAP', 'Mono');
        if strcmp(st.mode, 'JONSWAP')
            setSlider(s1, [3.5 7.5], st.wp, 3.5:1:7.5); s1.Enable = 'on';
            setSlider(s2, [1 7], st.gJ, 1:7); s2.Enable = 'on';
            lS1.Text = ['Peak ' cW 'p (rad/s)']; lS2.Text = ['Peakedness ' cG 'J'];
        else
            setSlider(s1, [wN(1) wN(end)], st.w, 3:1:8); s1.Enable = 'on';
            set(s2, 'Enable', 'off');
            lS1.Text = ['Frequency ' cW ' (rad/s)']; lS2.Text = '';
        end
        seaChanged();
    end

    function seaChanged()
        if strcmp(st.mode, 'JONSWAP'), st.wp = s1.Value; st.gJ = s2.Value; else, st.w = s1.Value; end
        if isfield(cur, 'Pz'), rebuild(); end
    end

    function selectFloe(k)
        st.floe = k;
        f = cache.floes(k);
        isGrid = strcmp(f.kind, 'grid');
        if isGrid, lFloe.Text = sprintf('Showing: %s', f.name);
        else, lFloe.Text = sprintf('Showing preset: %s  (sliders not in use)', f.name); end
        set([sB sR], 'Enable', 'on');
        rebuild();
    end

    function rebuild()
        % matrices for the current floe and sea, then scales, then static geometry. The timer is held off
        % (st.busy) while cur is rebuilt, so it never draws a half-built state.
        st.busy = true;
        f = cache.floes(st.floe); pat = cache.pat{st.floe};
        R = f.p(3); nT = numel(pat.theta); nR = numel(pat.r);
        if strcmp(st.mode, 'JONSWAP')
            w = wF; dw = wF(2) - wF(1);
            Sw = jonswap(wF, st.wp, st.gJ, Hs);
            amp = sqrt(2 * Sw * dw);
            ph = epsF;
        else
            w = st.w; amp = aMono; ph = 0; Sw = [];
        end
        cur.w = w; cur.c = amp .* exp(1i * ph);
        cur.Pz = interpW(reshape(double(pat.eta), nT * nR, []), pat.omega, w);
        cur.Pl = interpW(reshape(double(pat.lap), nT * nR, []), pat.omega, w);
        cur.nT = nT; cur.nR = nR;
        % the 33 candidate points lie exactly on the polar grid (r steps of R/20, theta steps of 5 deg)
        iT = mod(round(gridFrac(:, 2) / (2 * pi) * (nT - 1)), nT - 1) + 1;
        iRr = round(gridFrac(:, 1) * (nR - 1)) + 1;
        cur.cand = sub2ind([nT nR], iT, iRr);
        cur.cx = gridFrac(:, 1) * R .* cos(gridFrac(:, 2)) * H;
        cur.cy = gridFrac(:, 1) * R .* sin(gridFrac(:, 2)) * H;
        % cross-section: theta = pi ray (reversed) then theta = 0 ray, and the incident wave over +-2R
        iPi = round((nT - 1) / 2) + 1;
        idxPi = sub2ind([nT nR], repmat(iPi, 1, nR), nR:-1:1);
        idx0 = sub2ind([nT nR], ones(1, nR - 1), 2:nR);
        cur.xIdx = [idxPi idx0];
        cur.xFloe = [-fliplr(pat.r) pat.r(2:end)] * H;
        xInc = linspace(-2 * R, 2 * R, 241);
        k0 = interp1(pat.omega, pat.k0, w, 'linear');
        cur.xInc = xInc * H;
        cur.Pinc = -exp(xInc(:) * k0(:).');                % per unit amplitude, solver convention (see test 2)
        % sensors: physical acceleration = -omega^2 * deflection, at every candidate point
        cur.Acand = cur.Pz(cur.cand, :) .* (-(w(:).'.^2));
        % scales from a sample of times
        ts = linspace(0, 60, 61);
        E = exp(1i * w(:) * ts);
        Z = real(cur.Pz * (cur.c(:) .* E));
        cur.zMax = max(abs(Z), [], 'all');
        Zp = real(interpW(reshape(double(cache.pat{iP0}.lap), [], numel(cache.pat{iP0}.omega)), ...
            cache.pat{iP0}.omega, w) * (cur.c(:) .* E));
        cur.lapRef = quantileAbs(Zp(:), 0.99) / H^2;        % colour scale: p0 under the same sea
        cur.exag = niceExag(0.18 * R * H / max(cur.zMax, eps));
        A = real(cur.Acand * (cur.c(:) .* E)); cur.aMax = max(abs(A), [], 'all');
        % geometry
        [TT, RRg] = ndgrid(pat.theta, pat.r);
        cur.X = RRg .* cos(TT) * H; cur.Y = RRg .* sin(TT) * H;
        Rm = R * H;
        set(hSurf, 'XData', cur.X, 'YData', cur.Y, 'ZData', zeros(size(cur.X)), 'CData', zeros(size(cur.X)));
        caxis(ax3, [-1 1] * cur.lapRef);
        zl = 0.3 * Rm;
        set(ax3, 'XLim', [-1.05 1.7] * Rm, 'YLim', [-1.05 1.05] * Rm, 'ZLim', [-zl zl], 'DataAspectRatio', [1 1 1]);
        set(hArrow, 'XData', 1.65 * Rm, 'YData', 0, 'ZData', 0, 'UData', -0.45 * Rm, 'VData', 0, 'WData', 0);
        set(hTxt, 'Position', [1.45 * Rm, 0.22 * Rm, 0]);
        hTitle3.String = {f.name, sprintf(['\\rm\\fontsize{10}\\color[rgb]{0.3,0.3,0.3}vertical \\times%s;  colour = bending: ' ...
            '\\color[rgb]{0.07,0.29,0.48}concave\\color[rgb]{0.3,0.3,0.3} / white flat / ' ...
            '\\color[rgb]{0.92,0.41,0.20}convex'], num2str(cur.exag))};
        set(axX, 'XLim', [-2 2] * Rm, 'YLim', [-1.15 1.15] * 1e3 * max([cur.zMax, max(abs(cur.c)) * ...
            ternary(strcmp(st.mode, 'JONSWAP'), 3, 1)]));
        trace = struct('t', [], 'a', zeros(2, 0));
        set(axT, 'YLim', [-1.15 1.15] * max(cur.aMax, eps));
        drawSpectrum(Sw);
        updateInfo(f, pat);
        updateLayoutText();
        st.busy = false;
        frame(true);
    end

    function frame(force)
        if nargin < 1, force = false; end
        if ~isvalid(fig) || st.busy || ~isfield(cur, 'Pz'), return; end
        try
            drawFrame(force);
        catch err                                     % never let one bad frame stop the timer
            if ~st.warned, warning('floeFlexDemo:frame', '%s', err.message); st.warned = true; end
        end
    end

    function drawFrame(force)
        dtw = toc(st.lastTic); st.lastTic = tic;
        if st.playing, st.t = st.t + st.speed * min(dtw, 0.2);
        elseif ~force, return; end
        e = cur.c(:) .* exp(1i * cur.w(:) * st.t);
        z = reshape(real(cur.Pz * e), cur.nT, cur.nR);
        lp = reshape(real(cur.Pl * e), cur.nT, cur.nR) / H^2;
        zc = z(cur.cand);
        lay = LL(st.layout); s = lay.idx; rf = lay.ref; ot = trOther(lay);
        lift = 0.01 * cur.X(1, end);
        set(hSurf, 'ZData', z * cur.exag, 'CData', lp);
        set(hRim, 'XData', cur.X(:, end), 'YData', cur.Y(:, end), 'ZData', z(:, end) * cur.exag);
        set(hCand, 'XData', cur.cx, 'YData', cur.cy, 'ZData', zc * cur.exag + lift);
        set(hSens, 'XData', cur.cx(s), 'YData', cur.cy(s), 'ZData', zc(s) * cur.exag + 2 * lift);
        set(hRef, 'XData', cur.cx(rf), 'YData', cur.cy(rf), 'ZData', zc(rf) * cur.exag + 2 * lift);
        % cross-section, true scale (mm)
        set(hInc, 'XData', cur.xInc, 'YData', 1e3 * real(cur.Pinc * e));
        set(hDefl, 'XData', cur.xFloe, 'YData', 1e3 * z(cur.xIdx));
        onAxis = s(abs(cur.cy(s)) < 1e-9);
        set(hSensX, 'XData', cur.cx(onAxis), 'YData', 1e3 * zc(onAxis));
        % traces: the reference and one other sensor of the layout
        a = real(cur.Acand([rf ot], :) * e);
        trace.t(end + 1) = st.t; trace.a(:, end + 1) = a;
        keep = trace.t > st.t - traceWin; trace.t = trace.t(keep); trace.a = trace.a(:, keep);
        set(hTr1, 'XData', trace.t, 'YData', trace.a(1, :));
        set(hTr2, 'XData', trace.t, 'YData', trace.a(2, :));
        set(axT, 'XLim', [max(0, st.t - traceWin), max(traceWin, st.t)]);
        drawnow limitrate;
    end

    function drawSpectrum(Sw)
        cla(axS);
        wb = [cache.omega(1) cache.omega(end)];
        if strcmp(st.mode, 'JONSWAP')
            ww = linspace(2.5, 9, 400); SS = 1e4 * jonswap(ww, st.wp, st.gJ, Hs);   % cm^2 s
            yl = 1.15 * max(SS);
            fill(axS, [wb(1) wb(2) wb(2) wb(1)], [0 0 yl yl], ice, 'EdgeColor', 'none');
            plot(axS, ww, SS, '-', 'Color', navy, 'LineWidth', 2);
            ylabel(axS, 'S(\omega) (cm^2 s)');
            pSp.Title = sprintf('Sea spectrum: JONSWAP, Hs %.0f cm, peak period %.2f s', 100 * Hs, 2 * pi / st.wp);
        else
            yl = 1.3 * 100 * aMono;
            fill(axS, [wb(1) wb(2) wb(2) wb(1)], [0 0 yl yl], ice, 'EdgeColor', 'none');
            plot(axS, [st.w st.w], [0 100 * aMono], '-', 'Color', navy, 'LineWidth', 3);
            plot(axS, st.w, 100 * aMono, 'o', 'MarkerFaceColor', navy, 'MarkerEdgeColor', navy);
            ylabel(axS, 'amplitude (cm)');
            pSp.Title = sprintf('Sea: one frequency, period %.2f s, amplitude %.1f cm', 2 * pi / st.w, 100 * aMono);
        end
        set(axS, 'XLim', [2.5 9], 'YLim', [0 yl]);
        text(axS, mean(wb), 0.92 * yl, 'model band', 'Color', grey, 'HorizontalAlignment', 'center', 'FontSize', 9);
        xlabel(axS, '\omega (rad/s)');
    end

    function updateInfo(f, pat)
        R = f.p(3); beta = f.p(1); lf = beta^0.25;
        wRef = ternary(strcmp(st.mode, 'JONSWAP'), st.wp, st.w);
        k = abs(interp1(pat.omega, pat.k0, min(max(wRef, pat.omega(1)), pat.omega(end))));
        Hd = f.Hdisp; ts = sqrt(Hd / H);               % Froude time scale (1 for basin floes)
        sc = ''; if Hd ~= H, sc = '  (full scale)'; end
        infoVal(1).Text = sprintf('R %s,  %sf %s%s', lenStr(R * Hd), cL, lenStr(lf * Hd), sc);
        infoVal(2).Text = sprintf('%.2f', R / lf);
        infoVal(3).Text = sprintf('%s long, %.1f s;  kR %.2f', lenStr(2 * pi / k * Hd), 2 * pi / wRef * ts, k * R);
        infoVal(4).Text = sprintf('peak %.0f mm', 1e3 * cur.zMax);
        [txt, col] = regime(R / lf, k * R);
        infoVal(5).Text = txt; infoVal(5).FontColor = col;
    end

    function selectLayout(k)
        st.layout = k; dLay.Value = k;
        lay = LL(k);
        hTitleT.String = sprintf(['What the accelerometers record:  \\color[rgb]{0.92,0.41,0.20}s%d (reference)' ...
            '\\color[rgb]{0.07,0.29,0.48}  s%d (dashed)'], lay.ref, trOther(lay));
        trace = struct('t', [], 'a', zeros(2, 0));
        updateLayoutText();
        drawLayoutBars();
        frame(true);
    end

    function updateLayoutText()
        if ~hasLay, lLay.Text = LL(1).note; return; end
        lay = LL(st.layout);
        names = strjoin(arrayfun(@(i) sprintf('s%d', i), lay.idx, 'UniformOutput', false), ', ');
        if lay.singular
            t1 = sprintf('%s (reference s%d):  no information, uncertainty infinite (Fisher matrix singular).', ...
                names, lay.ref);
        else
            t1 = sprintf(['%s (reference s%d):  ' cB ' %.3f%%,  R %.4f%%,  d_W %.2f,  %s.'], names, lay.ref, ...
                100 * lay.sigma(1), 100 * lay.sigma(2), lay.dW, lay.class);
            if lay.recordFactor > 1.05
                t1 = sprintf('%s  Needs %.1f%s the Fisher pair''s record for the same precision.', t1, ...
                    lay.recordFactor, char(215));
            elseif lay.recordFactor < 0.95
                t1 = sprintf('%s  Reaches the Fisher pair''s precision with %.0f%% of its record.', t1, ...
                    100 * lay.recordFactor);
            end
        end
        if st.floe ~= iP0, t1 = [t1 '  (Computed for p0, not the floe shown.)']; end
        lLay.Text = t1;
        lLay.Tooltip = lay.note;
    end

    function drawLayoutBars()
        % highlight the selected layout's row in the comparison table
        for i = 1:size(tblLab, 1)
            sel = i == st.layout;
            set(tblLab(i, :), 'FontWeight', ternary(sel, 'bold', 'normal'), 'FontColor', ternary(sel, orange, [0 0 0]));
        end
    end

    function openM1()
        try
            addpath(drDir); demoM1Live;
        catch err
            uialert(fig, err.message, 'Inversion demo');
        end
    end

    function togglePlay()
        st.playing = ~st.playing; st.lastTic = tic;
        bPlay.Text = ternary(st.playing, 'Pause', 'Play');
    end

    function setSpeed(v)
        st.speed = str2double(erase(v, 'x'));
    end

    function closeApp()
        try, stop(tmr); delete(tmr); catch, end
        delete(fig);
    end
end

%% =================================================================================================================
%  helpers
%% =================================================================================================================
function Y = interpW(X, wN, w)
% complex responses (rows) at cached nodes wN, linearly interpolated to w
if isscalar(w)
    Y = interp1(wN(:), X.', w, 'linear').';
    Y = Y(:);
else
    Y = interp1(wN(:), X.', w(:), 'linear').';
end
end

function S = jonswap(w, wp, gJ, Hs)
% JONSWAP elevation spectrum (one-sided, per rad/s), scaled to significant wave height Hs
s = 0.07 + 0.02 * (w > wp);
S0 = w.^-5 .* exp(-1.25 * (wp ./ w).^4) .* gJ.^exp(-(w - wp).^2 ./ (2 * s.^2 * wp^2));
ww = linspace(0.3 * wp, 6 * wp, 4000);
s2 = 0.07 + 0.02 * (ww > wp);
m0 = trapz(ww, ww.^-5 .* exp(-1.25 * (wp ./ ww).^4) .* gJ.^exp(-(ww - wp).^2 ./ (2 * s2.^2 * wp^2)));
S = S0 * (Hs / 4)^2 / m0;
end

function [txt, col] = regime(q, kR)
if q < 0.96
    txt = 'Nearly rigid: stiffness not identifiable';
    col = [0.75 0.2 0.15];
elseif q < 1.64
    txt = 'Transitional (A1 bracket 0.96 to 1.64)';
    col = [0.85 0.5 0.1];
elseif q > 8
    txt = 'Flexes strongly: Welch bias limits it';
    col = [0.85 0.5 0.1];
elseif kR < 1
    txt = 'Flexes, but small vs wavelength';
    col = [0.85 0.5 0.1];
else
    txt = 'Flexes: stiffness identifiable';
    col = [0.13 0.5 0.25];
end
end

function o = trOther(lay)
% the sensor traced against the reference: the last non-reference sensor of the layout
oth = lay.idx(lay.idx ~= lay.ref); o = oth(end);
end

function s = factorStr(f)
if ~isfinite(f), s = 'an infinite multiple of'; elseif abs(f - 1) < 0.05, s = 'the same as'; ...
else, s = sprintf('%.1f%s', f, char(215)); end
end

function setSlider(s, lim, v, ticks)
% change limits and value without ever leaving the slider in an invalid state
s.Limits = [min(s.Limits(1), lim(1)), max(s.Limits(2), lim(2))];
s.Value = min(max(v, lim(1)), lim(2));
s.Limits = lim; s.MajorTicks = ticks; s.MinorTicks = [];
end

function q = quantileAbs(x, p)
x = sort(abs(x(:)));
q = x(max(1, round(p * numel(x))));
end

function e = niceExag(x)
if x <= 1, e = 1; return; end
p = 10^floor(log10(x)); m = x / p;
if m < 2, m = 1; elseif m < 5, m = 2; else, m = 5; end
e = m * p;
end

function s = lenStr(L)
if L >= 1, s = sprintf('%.3g m', L); elseif L >= 0.01, s = sprintf('%.3g cm', 100 * L); else, s = sprintf('%.3g mm', 1e3 * L); end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function C = divergingMap(n, cNeg, cPos)
h = floor(n / 2); w = [1 1 1];
C = [interp1([0 1], [cNeg; w], linspace(0, 1, h)); interp1([0 1], [w; cPos], linspace(0, 1, n - h))];
end