function files = plotR3NoiseMisspecification(results, cfg, run)
%PLOTR3NOISEMISSPECIFICATION  Draw and save the R3 figure from saved results (no EMM, seconds).
%   plotR3NoiseMisspecification()                     redraw from the newest Results/R3/R3_*.mat
%   plotR3NoiseMisspecification('...\R3_<stamp>.mat') redraw from that run
%   plotR3NoiseMisspecification(results, cfg, run)    called by runR3NoiseMisspecification
%
%   (a) |delta beta / beta| and (b) bias significance d_sys against the error eps in the assumed
%   reference noise floor, one line per noise case (log y): a smooth curve on the fine tolerance grid
%   with markers at the pre-registered eps values. eps = -100% is "no correction"; dashed
%   vertical line there. Thresholds as ticks: 1% in (a); d = 0.5 and 1 in (b). Points where the
%   corrected denominator would reach zero (d_sys = Inf) are left out.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 1 || ischar(results) || isstring(results)
    P = setupStudyPaths();
    if nargin < 1 || isempty(results)
        Lst = dir(fullfile(P.results, 'R3', 'R3_*.mat'));
        if isempty(Lst), error('plotR3NoiseMisspecification:none', 'No R3 results in %s', fullfile(P.results, 'R3')); end
        [~, k] = max([Lst.datenum]);
        matFile = fullfile(Lst(k).folder, Lst(k).name);
    else
        matFile = char(results);
    end
    M = load(matFile);
    results = M.results; cfg = M.cfg; %#ok<NASGU>
    run = struct('id', 'R3', 'dir', fileparts(matFile));
    fprintf('  redrawing R3 figure from %s\n', matFile);
end
S = sharcStyle();
res = results.res;
ePct = 100 * results.epsGrid;

[fig, ax] = sharcFigure('full', [1 2], [], 'left', 1.75, 'hgap', 1.9);   % room for the log tick labels
hh = []; lab = {};
for n = 1:numel(res)
    y1 = 100 * abs(res(n).dBeta); y1(~isfinite(res(n).dSys) | y1 <= 0) = NaN;
    y2 = res(n).dSys;             y2(~isfinite(y2) | y2 <= 0) = NaN;
    if isfield(res, 'epsFine') && ~isempty(res(n).epsFine)
        % smooth curve on the fine grid, markers at the pre-registered eps values
        eF = 100 * res(n).epsFine;
        f1 = 100 * res(n).dBetaFine; f1(~isfinite(f1) | f1 <= 0) = NaN;
        f2 = res(n).dSysFine;        f2(~isfinite(f2) | f2 <= 0) = NaN;
        sharcSensorSeries(ax(1), eF, f1, res(n).noise, 'Marker', 'none', 'HandleVisibility', 'off');
        sharcSensorSeries(ax(2), eF, f2, res(n).noise, 'Marker', 'none', 'HandleVisibility', 'off');
        sharcSensorSeries(ax(1), ePct, y1, res(n).noise, 'LineStyle', 'none', 'HandleVisibility', 'off');
        sharcSensorSeries(ax(2), ePct, y2, res(n).noise, 'LineStyle', 'none', 'HandleVisibility', 'off');
        hh(end + 1) = sharcSensorSeries(ax(1), NaN, NaN, res(n).noise); %#ok<AGROW>   legend: line + marker
    else
        hh(end + 1) = sharcSensorSeries(ax(1), ePct, y1, res(n).noise); %#ok<AGROW>
        sharcSensorSeries(ax(2), ePct, y2, res(n).noise);
    end
    lab{end + 1} = S.sensorLabel(res(n).noise); %#ok<AGROW>
end
set(ax(1), 'YScale', 'log'); set(ax(2), 'YScale', 'log');
a1 = 100 * abs([res.dBeta]); a1 = a1(isfinite(a1) & a1 > 0); a1 = a1(:);
a2 = [res.dSys]; a2 = a2(isfinite(a2) & a2 > 0); a2 = a2(:);
% lower limits from the pre-registered points only (the fine curves dip towards zero at eps = 0 and are
% simply clipped there); upper limits also cover the fine curves, ignoring near-singular spikes
lo1 = 10^floor(log10(min(a1))); lo2 = 10^floor(log10(min(a2)));
if isfield(res, 'dBetaFine')
    f1 = 100 * [res.dBetaFine]; f2 = [res.dSysFine];
    f1 = f1(isfinite(f1) & f1 > 0 & f1 <= 100 * max(a1));
    f2 = f2(isfinite(f2) & f2 > 0 & f2 <= 10 * max(a2));
    a1 = [a1(:); f1(:)];
    a2 = [a2(:); f2(:)];
end
ylim(ax(1), [lo1, 10^ceil(log10(max([a1(:); 1.5])))]);
ylim(ax(2), [lo2, 10^ceil(log10(max([a2(:); 1.5])))]);
sharcTicks(ax(1), 'y', decades(get(ax(1), 'YLim')));
sharcTicks(ax(2), 'y', decades(get(ax(2), 'YLim')));
for k = 1:2
    xlim(ax(k), [-108 58]);
    sharcTicks(ax(k), 'x', [-100 -50 -30 -10 10 30 50]);
    xlabel(ax(k), 'error in assumed noise floor, \epsilon (%)');
end
ylabel(ax(1), '|\delta\beta/\beta| (%)');
ylabel(ax(2), 'bias significance d_{sys}');
sharcThreshold(ax(1), 'y', 1, '', 'tick');                                      % |dbeta/beta| = 1%
sharcThreshold(ax(2), 'y', 0.5, '', 'tick'); sharcThreshold(ax(2), 'y', 1, '', 'tick');   % d = 0.5, 1
sharcThreshold(ax(1), 'x', -100, 'no correction', 'top');
sharcThreshold(ax(2), 'x', -100, 'no correction', 'top');
sharcPanel(ax(1), 'a', 'Stiffness error'); sharcPanel(ax(2), 'b', 'Bias significance by noise case');
sharcSharedLegend(fig, ax, hh, lab, 3);
files = saveSharcFigure(fig, run, 1, 'noise_misspecification');
end

function t = decades(lim)
t = 10 .^ (ceil(log10(lim(1))):floor(log10(lim(2))));
if any(abs(t - 1) < 1e-12), t = unique([t 0.5]); end   % keep 0.5 visible where 1 is
t = t(t >= lim(1) & t <= lim(2));
end