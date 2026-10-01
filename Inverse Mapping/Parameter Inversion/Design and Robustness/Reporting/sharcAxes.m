function sharcAxes(ax)
%SHARCAXES  Apply the standard axis style: recessive axes and hairline grid, outward ticks, no box.
%   sharcAxes(ax)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
s = sharcStyle();
hold(ax, 'on');
set(ax, 'FontName', s.font, 'FontSize', s.fontSize, 'XColor', s.axis, 'YColor', s.axis, ...
    'LineWidth', s.axisLineWidth, 'Box', 'off', 'TickDir', 'out', 'Color', 'w');
grid(ax, 'on');
try
    set(ax, 'GridColor', s.grid, 'GridAlpha', 1, 'MinorGridColor', s.grid, 'MinorGridAlpha', 0.6);
catch
    % Octave: grid colour properties differ; the default light grid is acceptable for tests
end
end
