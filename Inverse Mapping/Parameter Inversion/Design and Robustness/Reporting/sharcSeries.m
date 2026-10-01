function h = sharcSeries(ax, x, y, slot, varargin)
%SHARCSERIES  Plot one data series in the standard style for categorical slot k.
%   h = sharcSeries(ax, x, y, 1)                         slot 1: colour, marker, line width, size
%   h = sharcSeries(ax, x, y, 2, 'DisplayName', 'R')     any plot property may follow and overrides
%   h = sharcSeries(ax, x, y, 3, 'Marker', 'none')       a plain line
%   h = sharcSeries(ax, x, y, 3, 'LineStyle', 'none')    markers only
%
%   Every data series in a study figure should go through this function (or sharcSensorSeries for
%   the noise cases), so the style lives only in sharcStyle. Lines are solid by default: dashed lines
%   are reserved for thresholds (sharcThreshold), and identity is carried by colour AND marker.
%   Slots are used in fixed order (1 for the first entity, 2 for the second, ...), never cycled.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
s = sharcStyle();
if slot < 1 || slot > size(s.series, 1)
    error('sharcSeries:slot', 'slot must be 1..%d', size(s.series, 1));
end
c = s.series(slot, :);
h = plot(ax, x, y, 'LineStyle', '-', 'Color', c, 'LineWidth', s.lineWidth, 'Marker', s.markers{slot}, ...
    'MarkerSize', s.markerSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, varargin{:});
end
