function lg = sharcSharedLegend(fig, ax, handles, labels, nCols)
%SHARCSHAREDLEGEND  One legend for several panels, in a band below them, so no legend covers data.
%   lg = sharcSharedLegend(fig, ax, handles, labels)          up to 6 entries per row
%   lg = sharcSharedLegend(fig, ax, handles, labels, 3)       3 columns
%
%   The legend is attached to the LAST axes (an in-panel legend on another axes is kept).
%   The figure grows downwards by a legend band; the panels keep their size and alignment. Call it after the
%   axis labels are set and before saveSharcFigure.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
s = sharcStyle();
n = numel(labels);
if nargin < 5 || isempty(nCols), nCols = min(n, 6); end
nRows = ceil(n / nCols);
% grow the figure downwards by a band and move every panel up by the same amount (centimetres)
set(fig, 'Units', 'centimeters');
figPos = get(fig, 'Position');
bandCm = 0.42 * nRows + 0.3;
for k = 1:numel(ax)
    set(ax(k), 'Units', 'centimeters');
    p = get(ax(k), 'Position');
    set(ax(k), 'Position', [p(1), p(2) + bandCm, p(3), p(4)]);
end
H = figPos(4) + bandCm;
set(fig, 'Position', [figPos(1:3) H], 'PaperSize', [figPos(3) H], 'PaperPosition', [0 0 figPos(3) H]);
lg = legend(ax(end), handles, labels, 'Box', 'off', 'FontSize', s.smallFontSize, 'FontName', s.font, ...
    'TextColor', s.ink);
try
    set(lg, 'NumColumns', nCols);
catch
    set(lg, 'Orientation', 'horizontal');
end
set(lg, 'Units', 'centimeters');
lp = get(lg, 'Position');
set(lg, 'Position', [(figPos(3) - lp(3)) / 2, 0.08, lp(3), bandCm - 0.1]);
end
