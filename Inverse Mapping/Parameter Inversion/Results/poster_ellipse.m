%PLOTPOSTERELLIPSE  Poster figure: predicted Fisher ellipse against simulated estimates.
%
%   Redraws panel (a) of the transmissibility twin-inversion check (truth A,
%   1.2 beta0 and 0.97 R0) for the poster, and writes ellipse.pdf at 130 x 65 mm,
%   the size the poster reserves for it, so the font sizes below are the printed sizes.
%
%   Cleaned up compared with the original panel:
%     - axes in percent (100 x Delta ln, i.e. relative error), so no x10^-3 exponent
%     - only the 1-sigma and 95% ellipses, the 200 estimates and the Welch bias
%       (the origin is the truth, since the axes are errors from it)
%     - ellipses labelled directly; a two-entry legend in the empty upper-left corner
%     - poster colours and Times text
%
%   Input: the .mat file saved by runTransmissibilityTwinInversion (run 20260927_143438).
%   Its field names are not hard-coded: the script finds the covariance and the 200
%   estimates by matching the values printed in the run log. If it cannot, it lists
%   every numeric array in the file so the right one can be named in FORCE below.

clear; clc;

%% ---- settings ---------------------------------------------------------------------
matFile = fullfile('..', 'Results', 'transmissibilityTwinInversion_20260927_143438.mat');
outFile = 'ellipse.pdf';
figSize = [13 6.5];                    % cm: width x height on the poster

% values printed in the run log for truth A
sigLog   = [0.00498 0.00095];          % predicted sigma(theta) = sqrt(diag(F^-1))
biasLog  = [0.00108 0.00013];          % Welch-expected bias
meanLog  = [0.00114 0.00015];          % mean offset of the 200 oracle estimates
ratioLog = [1.11 1.10];                % empirical / predicted scatter
nEst     = 200;

% If the automatic search fails, name the arrays here, e.g. FORCE.C = 'S.res(1).Cpred'
FORCE.C   = '';                        % 2x2 covariance F^-1 (or F itself)
FORCE.est = '';                        % 200x2 estimate offsets from truth

% poster palette
navy    = [ 11  37  69] / 255;
ocean   = [ 19  73 123] / 255;
glacial = [ 91 155 200] / 255;
accent  = [235 104  52] / 255;
fontName = 'Times New Roman';
fsTick = 16; fsLabel = 18; fsLegend = 15;   % pt, as printed

%% ---- load and find the arrays -------------------------------------------------------
S = load(matFile);
items = collectNumeric(S, 'S', {});

C = [];
if ~isempty(FORCE.C)
    C = eval(FORCE.C);
    if any(abs(sqrt(abs(diag(C)))'./sigLog - 1) > 0.5), C = inv(C); end   % F given, not F^-1
end
for k = 1:numel(items)
    if ~isempty(C), break; end
    v = items{k}.value;
    if ~isequal(size(v), [2 2]) || ~isreal(v) || norm(v - v.') > 1e-9 * norm(v), continue; end
    if all(abs(sqrt(abs(diag(v)))'./sigLog - 1) < 0.03)
        C = v; fprintf('covariance: %s\n', items{k}.name);
    elseif det(v) > 0 && all(abs(sqrt(abs(diag(inv(v))))'./sigLog - 1) < 0.03)
        C = inv(v); fprintf('Fisher matrix (inverted): %s\n', items{k}.name);
    end
end

est = [];
if ~isempty(FORCE.est), est = eval(FORCE.est); if size(est, 2) ~= 2, est = est.'; end, end
bestScore = inf;
for k = 1:numel(items)
    if ~isempty(FORCE.est), break; end
    v = items{k}.value;
    if isequal(size(v), [2 nEst]), v = v.'; end
    if ~isequal(size(v), [nEst 2]) || ~isreal(v), continue; end
    score = norm((mean(v) - meanLog) ./ sigLog) + norm(std(v) ./ sigLog - ratioLog);
    if score < bestScore, bestScore = score; est = v; estName = items{k}.name; end
end
if ~isempty(est) && isempty(FORCE.est)
    if bestScore > 0.5, est = []; else, fprintf('estimates: %s (match score %.2f)\n', estName, bestScore); end
end

if isempty(C) || isempty(est)
    fprintf('\nCould not identify the arrays automatically. Numeric arrays in the file:\n');
    for k = 1:numel(items)
        fprintf('  %-60s %s\n', items{k}.name, mat2str(size(items{k}.value)));
    end
    error(['Set FORCE.C and/or FORCE.est at the top of the script to the right names ' ...
           '(the covariance is 2x2; the estimates are %dx2).'], nEst);
end

% sanity check against the log
fprintf('predicted sigma  [%.5f %.5f]   (log [%.5f %.5f])\n', sqrt(diag(C)), sigLog);
fprintf('estimate mean    [%+.5f %+.5f] (log [%+.5f %+.5f])\n', mean(est), meanLog);
fprintf('scatter / sigma  [%.2f %.2f]     (log [%.2f %.2f])\n', std(est)./sqrt(diag(C))', ratioLog);

%% ---- plot (percent: 100 x Delta ln = relative error) ---------------------------------
P  = 100 * est;
Cp = 1e4 * C;
b  = 100 * biasLog;

fig = figure('Units', 'centimeters', 'Position', [2 2 figSize], 'Color', 'w');
% MATLAB R2025a+ opens figures in the dark theme when the desktop is dark; force light
if exist('theme', 'file') || exist('theme', 'builtin'), theme(fig, 'light'); end
ax  = axes(fig, 'Color', 'w'); hold(ax, 'on');

hS = scatter(ax, P(:,1), P(:,2), 12, glacial, 'filled', ...
             'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
[x1, y1] = ellipseXY(Cp, 1);                    % 1 sigma (Mahalanobis radius 1)
[x2, y2] = ellipseXY(Cp, sqrt(-2*log(0.05)));   % 95% for two parameters (chi2_2)
h1 = plot(ax, x1, y1, '-',  'Color', ocean, 'LineWidth', 2.4);
h2 = plot(ax, x2, y2, '--', 'Color', ocean, 'LineWidth', 1.6);
hB = plot(ax, b(1), b(2), 'd', 'MarkerSize', 10, 'MarkerFaceColor', accent, ...
          'MarkerEdgeColor', 'w', 'LineWidth', 0.8);

% limits: the 95% ellipse plus a margin, with room on the left for the legend
xr = [min(x2) max(x2)]; yr = [min(y2) max(y2)];
xlim(ax, xr + [-0.12 0.05] * diff(xr));
ylim(ax, yr + [-0.08 0.38] * diff(yr));

% direct labels at the top of each ellipse
[~, i1] = max(y1); [~, i2] = max(y2);
text(ax, x1(i1), y1(i1), '1\sigma', 'Color', ocean, 'FontName', fontName, ...
     'FontSize', fsLegend, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
     'BackgroundColor', 'w', 'Margin', 0.5);
text(ax, x2(i2), y2(i2), '95%', 'Color', ocean, 'FontName', fontName, ...
     'FontSize', fsLegend, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
     'BackgroundColor', 'w', 'Margin', 0.5);

set(ax, 'FontName', fontName, 'FontSize', fsTick, 'LineWidth', 0.8, 'Box', 'on', ...
        'TickDir', 'in', 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, ...
        'XColor', navy, 'YColor', navy, 'Layer', 'top');
xlabel(ax, 'stiffness error (%)', 'FontName', fontName, 'FontSize', fsLabel);
ylabel(ax, 'radius error (%)',    'FontName', fontName, 'FontSize', fsLabel);

lg = legend(ax, [hS hB], {'simulated estimates', 'predicted bias'}, ...
            'Location', 'northwest', 'Box', 'off', 'FontName', fontName, ...
            'FontSize', fsLegend, 'TextColor', navy, 'Color', 'none');
lg.ItemTokenSize = [18 10];

%% ---- export at the printed size -----------------------------------------------------
exportgraphics(fig, outFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('saved %s\n', outFile);

%% ---- local functions ----------------------------------------------------------------
function [x, y] = ellipseXY(C, r)
% Points on the ellipse  v' * inv(C) * v = r^2.
    t = linspace(0, 2*pi, 400);
    [V, D] = eig((C + C.') / 2);
    xy = V * sqrt(D) * [cos(t); sin(t)] * r;
    x = xy(1, :); y = xy(2, :);
end

function items = collectNumeric(x, name, items)
% Every numeric array inside x (structs, struct arrays and cells, searched recursively).
    if isnumeric(x) && ~isempty(x)
        items{end+1} = struct('name', name, 'value', double(x));
    elseif isstruct(x)
        f = fieldnames(x);
        for e = 1:numel(x)
            for i = 1:numel(f)
                if numel(x) > 1, nm = sprintf('%s(%d).%s', name, e, f{i});
                else,            nm = sprintf('%s.%s', name, f{i}); end
                items = collectNumeric(x(e).(f{i}), nm, items);
            end
        end
    elseif iscell(x)
        for e = 1:numel(x)
            items = collectNumeric(x{e}, sprintf('%s{%d}', name, e), items);
        end
    end
end