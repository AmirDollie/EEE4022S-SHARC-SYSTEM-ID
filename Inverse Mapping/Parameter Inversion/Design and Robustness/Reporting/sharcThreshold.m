function h = sharcThreshold(ax, axisName, value, label, where)
%SHARCTHRESHOLD  Dashed reference line with a small label, in the reserved threshold colour.
%   CALL IT LAST on that axes: after plotting the data and setting xlim/ylim (and the log scale). The
%   line spans the axis limits at the time of the call and those limits are then frozen.
%   sharcThreshold(ax, 'y', 0.5, 'd = 0.5')              horizontal line at y = 0.5
%   sharcThreshold(ax, 'x', 109, 'production', 'top')     vertical line at x = 109
%   where: y lines 'right' (default, inside at the right end) | 'left' | 'outside' (right margin)
%          | 'tick' (no text; the value is added as a labelled y tick: use where the data crowd the line);
%          x lines 'top' (default) | 'bottom'
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
s = sharcStyle();
if nargin < 4, label = ''; end
xl = xlim(ax); yl = ylim(ax);
if strcmpi(axisName, 'y')
    h = plot(ax, xl, [value value], '--', 'Color', s.threshold, 'LineWidth', s.thinLineWidth, ...
        'HandleVisibility', 'off');
    if nargin < 5, where = 'right'; end
    if strcmp(where, 'tick')
        % no text: the value becomes a labelled tick, and the caption names the line
        t = unique([get(ax, 'YTick'), value]);
        sharcTicks(ax, 'y', t);
    elseif ~isempty(label)
        % 'right' / 'left': inside the axes at that end, just above the line;
        % 'outside': in the right margin, beside the line (needs sharcStyle.margin.right room)
        switch where
            case 'left',    x = xl(1); ha = 'left';  va = 'bottom'; lbl = [' ' label];
            case 'outside', x = xl(2); ha = 'left';  va = 'middle'; lbl = [' ' label];
            otherwise,      x = xl(2); ha = 'right'; va = 'bottom'; lbl = [label ' '];
        end
        text(ax, x, value, lbl, 'Color', s.axis, 'FontSize', s.smallFontSize, 'FontName', s.font, ...
            'HorizontalAlignment', ha, 'VerticalAlignment', va, 'Clipping', 'off');
    end
else
    h = plot(ax, [value value], yl, '--', 'Color', s.threshold, 'LineWidth', s.thinLineWidth, ...
        'HandleVisibility', 'off');
    if nargin < 5, where = 'top'; end
    if ~isempty(label)
        if strcmp(where, 'bottom'), y = yl(1); va = 'bottom'; else, y = yl(2); va = 'top'; end
        text(ax, value, y, [' ' label], 'Color', s.axis, 'FontSize', s.smallFontSize, 'FontName', s.font, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', va, 'Rotation', 0);
    end
end
xlim(ax, xl); ylim(ax, yl);
end
