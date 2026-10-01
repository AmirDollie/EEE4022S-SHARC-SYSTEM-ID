function sharcPanel(ax, letter, titleText)
%SHARCPANEL  Panel label "(a)" and optional short title, left-aligned above the axes.
%   sharcPanel(ax, 'a')
%   sharcPanel(ax, 'b', 'Bias significance')
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
s = sharcStyle();
if nargin < 3, titleText = ''; end
str = sprintf('(%s)', letter);
if ~isempty(titleText), str = sprintf('(%s) %s', letter, titleText); end
t = title(ax, str, 'FontName', s.font, 'FontSize', s.fontSize, 'FontWeight', 'normal', 'Color', s.ink);
try
    set(ax, 'TitleHorizontalAlignment', 'left');
catch
    set(t, 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02 0]);
end
end
