function lg = sharcLegend(ax, handles, labels, location)
%SHARCLEGEND  Standard legend: small font, no box, entries in the order given.
%   lg = sharcLegend(ax, handles, labels)
%   lg = sharcLegend(ax, handles, labels, 'northwest')     (default 'best')
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
s = sharcStyle();
if nargin < 4, location = 'best'; end
lg = legend(ax, handles, labels, 'Location', location, 'Box', 'off', 'FontSize', s.smallFontSize, ...
    'FontName', s.font, 'TextColor', s.ink);
end
