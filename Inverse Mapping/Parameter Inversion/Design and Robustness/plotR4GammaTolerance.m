function files = plotR4GammaTolerance(results, cfg, run)
%PLOTR4GAMMATOLERANCE  Draw and save the R4 figure from saved results (no EMM, a few seconds).
%   plotR4GammaTolerance()                      redraw from the newest Results/R4/R4_*.mat
%   plotR4GammaTolerance('...\R4_20261001_183822.mat')   redraw from that run
%   plotR4GammaTolerance(results, cfg, run)     called by runR4GammaTolerance at the end of a run
%
%   Writes Results/R4/figures/R4_fig1_gamma_tolerance.pdf/.png (replacing the previous version).
%   (a) induced beta and R error against the gamma error, headline noise case, with the exact-EMM
%       check points; linear epsilon axis, legend in the empty lower-right quadrant.
%   (b) bias significance d_sys against |eps_gamma| on log-log axes, one line per noise case. d_sys is
%       proportional to |ln(1 + eps)|, so the lines are straight and parallel and the crossing of the
%       d = 0.5 and d = 1 thresholds (the tolerances) can be read directly. At each |eps| the larger
%       of the two signs is plotted (the conservative side, as the tolerances are defined).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        L = dir(fullfile(P.results, 'R4', 'R4_*.mat'));
        if isempty(L), error('plotR4GammaTolerance:none', 'No R4 results in %s', fullfile(P.results, 'R4')); end
        [~, k] = max([L.datenum]);
        matFile = fullfile(L(k).folder, L(k).name);
    else
        matFile = char(results);
    end
    M = load(matFile);
    results = M.results; cfg = M.cfg;
    run = struct('id', 'R4', 'dir', fileparts(matFile));
    fprintf('  redrawing R4 figure from %s\n', matFile);
end
S = sharcStyle();
res = results.res;
epsGrid = results.epsGrid(:).';
hIdx = find(strcmp({res.noise}, results.headline));
r = res(hIdx);
ePct = 100 * epsGrid;

[fig, ax] = sharcFigure('full', [1 2]);

%% (a) parameter error, headline noise case
% slots 6 and 5 are not used by any noise case, so these lines cannot be confused with panel (b)
h1 = sharcSeries(ax(1), ePct, 100 * r.dBeta, 6);
h2 = sharcSeries(ax(1), ePct, 100 * r.dR, 5);
hh = [h1 h2];
lab = {'\delta\beta/\beta (linear)', '\deltaR/R (linear)'};
nl = results.nonlinear;
if ~isempty(nl) && ~isempty(nl(1).dthetaNL)
    eN = 100 * [nl.eps];
    dN = [nl.dthetaNL];
    h3 = plot(ax(1), eN, 100 * (exp(dN(1, :)) - 1), 'o', 'Color', S.ink, 'MarkerSize', 7, 'LineWidth', 0.9);
    plot(ax(1), eN, 100 * (exp(dN(2, :)) - 1), 'o', 'Color', S.ink, 'MarkerSize', 7, 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');                  % check points: open black circles, not a series
    hh(end + 1) = h3; lab{end + 1} = 'exact-EMM check';
end
xlabel(ax(1), 'error in assumed \gamma, \epsilon_\gamma (%)');
ylabel(ax(1), 'induced parameter error (%)');
yMax = max(1.5, 1.15 * max(abs(100 * r.dBeta)));
xlim(ax(1), [-21 21]); ylim(ax(1), [-1 1] * yMax);
sharcTicks(ax(1), 'x', [-20 -10 0 10 20]);
sharcTicks(ax(1), 'y', [-2 -1 0 1 2]);
sharcThreshold(ax(1), 'y', 1, '1%', 'left'); sharcThreshold(ax(1), 'y', -1, '', 'left');
sharcLegend(ax(1), hh, lab, 'southeast');
sharcPanel(ax(1), 'a', sprintf('Parameter error (%s)', S.sensorLabel(results.headline)));

%% (b) bias significance against |eps|, every noise case, log-log
aE = unique(abs(epsGrid));
hb = []; lb = {};
for n = 1:numel(res)
    d = zeros(size(aE));
    for k = 1:numel(aE)
        d(k) = max(res(n).dSys(abs(abs(epsGrid) - aE(k)) < 1e-12));   % worse of the two signs
    end
    hb(end + 1) = sharcSensorSeries(ax(2), 100 * aE, d, res(n).noise); %#ok<AGROW>
    lb{end + 1} = S.sensorLabel(res(n).noise); %#ok<AGROW>
end
set(ax(2), 'XScale', 'log', 'YScale', 'log');
xlabel(ax(2), 'error in assumed \gamma, |\epsilon_\gamma| (%)');
ylabel(ax(2), 'bias significance d_{sys}');
xlim(ax(2), [0.4 25]);
allD = [res.dSys];
ylim(ax(2), [10^floor(log10(min(allD(allD > 0)))) 10^ceil(log10(max(allD)))]);
sharcTicks(ax(2), 'x', [0.5 1 2 5 10 20]);
sharcTicks(ax(2), 'y', [0.01 0.1 1 10]);
sharcThreshold(ax(2), 'y', 0.5, '', 'tick'); sharcThreshold(ax(2), 'y', 1, '', 'tick');   % d = 0.5, 1 (caption)
sharcPanel(ax(2), 'b', 'Bias significance by noise case');

sharcSharedLegend(fig, ax, hb, lb, 3);
files = saveSharcFigure(fig, run, 1, 'gamma_tolerance');
end