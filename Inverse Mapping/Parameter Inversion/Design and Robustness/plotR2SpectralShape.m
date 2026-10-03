function files = plotR2SpectralShape(results, cfg, run)
%PLOTR2SPECTRALSHAPE  Draw and save the R2 figures from saved results (no EMM, no Fisher recomputation).
%   plotR2SpectralShape()                     redraw from the newest Results/R2/R2_*.mat
%   plotR2SpectralShape('...\R2_<stamp>.mat') redraw from that run
%   plotR2SpectralShape(results, cfg, run)    called by runR2SpectralShape
%
%   Figure 1  production L = 2048, the 15 unimodal seas: (a) sigma_lnbeta and (b) d_W against the peak
%             frequency, one line per gamma_J. The bimodal sea has no single peak frequency, so it is
%             not placed on this axis; its values are annotated.
%   Figure 2  the L comparison over all 16 seas: (a) d_W and (b) sigma_lnbeta at L = 1024, 2048, 4096,
%             coloured with the D2 slots (slot = position of L in cfg.D2.L), so L looks the same as in D2
%   Figure 3  mechanism at L = 2048 for four seas: (a) incident PSD on the bins (all seas have the same
%             in-band variance), (b) share of the diagonal Fisher information F_beta,beta per accepted bin
%             (not a decomposition of sigma_lnbeta, which also depends on F_betaR and F_RR), (c) local
%             feature-bias significance sqrt(b_W,k' Sigma_k^-1 b_W,k) per bin (a diagnostic, not a
%             decomposition of d_W: the global bias is projected through F and contributions can cancel)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'R2', 'R2_*.mat'));
        if isempty(Lst), error('plotR2SpectralShape:none', 'No R2 results in %s', fullfile(P.results, 'R2')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M0 = load(matFile);
    results = M0.results; cfg = M0.cfg;
    run = struct('id', 'R2', 'dir', fileparts(matFile));
    fprintf('  redrawing R2 figures from %s\n', matFile);
end
S = sharcStyle();
seas = results.seas; M = results.metrics; Ls = results.L; iL = results.iL; mech = results.mechanism;
nQ = numel(seas); nL = numel(Ls);
files = {};
isUni = strcmp({seas.type}, 'jonswap');
wp1 = arrayfun(@(x) x.omegaP(1), seas); gj1 = arrayfun(@(x) x.gammaJ(1), seas);
gJs = cfg.R2.gammaJ; wps = cfg.R2.omegaP;
iBi = find(~isUni, 1);
B = cfg.R2.bimodal;
gjLabel = @(g) sprintf('\\gamma_J = %g', g);
nsum = @(x) sum(x(~isnan(x)));                 % NaN = rejected bin (portable, no 'omitnan')

%% ---- figure 1: production L against peak frequency -----------------------------------------------------------
[fig, ax] = sharcFigure('full', [1 2], [], 'hgap', 1.9);
h = []; lab = {};
for g = 1:numel(gJs)
    q = find(isUni & gj1 == gJs(g));
    [~, o] = sort(wp1(q)); q = q(o);
    sharcSeries(ax(1), wp1(q), 100 * M.sigma(q, iL, 1), g);
    sharcSeries(ax(2), wp1(q), M.dW(q, iL), g);
    h(end + 1) = sharcSeries(ax(end), NaN, NaN, g); %#ok<AGROW> legend proxy in the legend's own axes
    lab{end + 1} = gjLabel(gJs(g)); %#ok<AGROW>
end
for k = 1:2
    xlim(ax(k), [3.2 7.8]); sharcTicks(ax(k), 'x', wps);
    xlabel(ax(k), 'peak frequency \omega_p (rad/s)');
end
ylim(ax(1), [0 1.1]); ylim(ax(2), [0 0.6]);
ylabel(ax(1), '\sigma_{ln\beta} (%)'); ylabel(ax(2), 'Welch bias significance d_W');
sharcPanel(ax(1), 'a', sprintf('Stiffness uncertainty (L = %d)', Ls(iL)));
sharcPanel(ax(2), 'b', sprintf('Welch bias (L = %d)', Ls(iL)));
sharcThreshold(ax(1), 'y', 100 * cfg.class.identifiable.maxSigma, '1%', 'right');
sharcThreshold(ax(2), 'y', cfg.class.identifiable.maxBiasDistance, 'd = 0.5', 'right');
if ~isempty(iBi)                                  % the bimodal sea: annotated, not placed on the omega_p axis
    note = {sprintf('bimodal (%g + %g rad/s): %.3f%%', B.omegaP, 100 * M.sigma(iBi, iL, 1)), ...
        sprintf('bimodal (%g + %g rad/s): d_W = %.3f', B.omegaP, M.dW(iBi, iL))};
    for a = 1:2
        text(ax(a), 0.97, 0.70, note{a}, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontSize', S.smallFontSize, 'FontName', S.font, 'Color', S.axis);
    end
end
sharcSharedLegend(fig, ax, h, lab, numel(lab));
files = [files, saveSharcFigure(fig, run, 1, 'shape_production')];

%% ---- figure 2: the L comparison over all seas ---------------------------------------------------------------------
% x: seas grouped by gamma_J (5 peak frequencies each), then the bimodal sea
order = []; xpos = []; tick = {}; x0 = 0;
for g = 1:numel(gJs)
    q = find(isUni & gj1 == gJs(g)); [~, o] = sort(wp1(q)); q = q(o);
    order = [order, q]; xpos = [xpos, x0 + (1:numel(q))]; %#ok<AGROW>
    tick = [tick, arrayfun(@(v) sprintf('%g', v), wp1(q), 'UniformOutput', false)]; %#ok<AGROW>
    x0 = x0 + numel(q) + 1;
end
if ~isempty(iBi), order = [order, iBi]; xpos = [xpos, x0 + 1]; tick = [tick, {'bimodal'}]; end
[fig, ax] = sharcFigure('full', [2 1], 12, 'vgap', 1.7, 'bottom', 1.5);       % stacked: 16 tick labels need the full width
lSlot = arrayfun(@(L) find(cfg.D2.L == L, 1), Ls);   % the D2 colour of each L
h = []; lab = {};
for k = 1:nL
    for g = 0:numel(gJs)                         % one segment per gamma_J group, no line across groups
        if g == 0, sel = 1:numel(order); else, sel = find(isUni(order) & gj1(order) == gJs(g)); end
        if g == 0
            sharcSeries(ax(1), xpos, M.dW(order, k), lSlot(k), 'LineStyle', 'none');          % all 16 markers
            sharcSeries(ax(2), xpos, 100 * M.sigma(order, k, 1), lSlot(k), 'LineStyle', 'none');
            hh = sharcSeries(ax(end), NaN, NaN, lSlot(k));
        else                                                                                % lines within a group
            sharcSeries(ax(1), xpos(sel), M.dW(order(sel), k), lSlot(k), 'Marker', 'none');
            sharcSeries(ax(2), xpos(sel), 100 * M.sigma(order(sel), k, 1), lSlot(k), 'Marker', 'none');
        end
    end
    h(end + 1) = hh; lab{end + 1} = sprintf('L = %d (%d segments)', Ls(k), M.K(k)); %#ok<AGROW>
end
for a = 1:2
    xlim(ax(a), [0.3 max(xpos) + 0.7]);
    set(ax(a), 'XTick', xpos, 'XTickLabel', tick, 'XTickLabelRotation', 0);   % MATLAB would tilt them otherwise
end
xlabel(ax(2), '\omega_p (rad/s), in groups \gamma_J = 1 | 3.3 | 7, then the bimodal sea');
for g = 1:numel(gJs)                              % group labels above the first panel
    sel = isUni(order) & gj1(order) == gJs(g);
    text(ax(1), mean(xpos(sel)), 1.22, gjLabel(gJs(g)), 'HorizontalAlignment', 'center', 'FontSize', S.smallFontSize, ...
        'FontName', S.font, 'Color', S.axis);
end
ylim(ax(1), [0 1.3]); ylim(ax(2), [0 1.1]);
ylabel(ax(1), 'Welch bias significance d_W'); ylabel(ax(2), '\sigma_{ln\beta} (%)');
sharcPanel(ax(1), 'a', 'Welch bias against L'); sharcPanel(ax(2), 'b', 'Stiffness uncertainty against L');
sharcThreshold(ax(1), 'y', cfg.class.identifiable.maxBiasDistance, 'd = 0.5', 'right');
sharcThreshold(ax(1), 'y', cfg.class.marginal.maxBiasDistance, 'd = 1', 'right');
sharcThreshold(ax(2), 'y', 100 * cfg.class.identifiable.maxSigma, '1%', 'right');
sharcSharedLegend(fig, ax, h, lab, nL);
files = [files, saveSharcFigure(fig, run, 2, 'L_comparison')];

%% ---- figure 3: mechanism --------------------------------------------------------------------------------------------
pick = [find(isUni & wp1 == 3.5 & gj1 == 7), results.iRef, find(isUni & wp1 == 7.5 & gj1 == 7), iBi];
w = mech.omega;
[fig, ax] = sharcFigure('full', [1 3], [], 'hgap', 1.75, 'left', 1.6, 'top', 0.75);
h = []; lab = {};
for j = 1:numel(pick)
    q = pick(j);
    sharcSeries(ax(1), w, 1e4 * mech.S(q, :), j, 'Marker', 'none');
    h(end + 1) = sharcSeries(ax(end), NaN, NaN, j, 'Marker', 'none'); %#ok<AGROW>
    sharcSeries(ax(2), w, 100 * mech.info(q, :, 1) / nsum(mech.info(q, :, 1)), j, 'Marker', 'none');
    sharcSeries(ax(3), w, mech.biasStd(q, :), j, 'Marker', 'none');
    lab{end + 1} = seaLabel(seas(q)); %#ok<AGROW>
end
for a = 1:3
    xlim(ax(a), [3 8.5]); sharcTicks(ax(a), 'x', [3 4 5 6 7 8]);
    xlabel(ax(a), '\omega (rad/s)');
end
ylabel(ax(1), 'S_\eta (10^{-4} m^2 s/rad)'); ylabel(ax(2), 'F_{\beta\beta} share (%)');
ylabel(ax(3), 'local bias');
sharcPanel(ax(1), 'a', 'Incident spectrum'); sharcPanel(ax(2), 'b', 'Stiffness information');
sharcPanel(ax(3), 'c', 'Welch feature bias');
sharcSharedLegend(fig, ax, h, lab, numel(lab));
files = [files, saveSharcFigure(fig, run, 3, 'mechanism')];

%% ---- printed check --------------------------------------------------------------------------------------------------
share = @(q, lo, hi) nsum(mech.info(q, w >= lo & w <= hi, 1)) / nsum(mech.info(q, :, 1));
fprintf('\n  share of F_beta,beta above 6 rad/s at L = %d (accepted bins):\n', Ls(iL));
for q = pick
    fprintf('    %-22s %.0f%%\n', seas(q).name, 100 * share(q, 6, Inf));
end
end

function s = seaLabel(x)
if strcmp(x.type, 'jonswap')
    s = sprintf('\\omega_p = %g, \\gamma_J = %g', x.omegaP, x.gammaJ);
else
    s = sprintf('bimodal %g + %g', x.omegaP);
end
end