function [fig, ax] = sharcFigure(width, nPanels, heightCm, varargin)
%SHARCFIGURE  New figure in the standard thesis style, sized in centimetres, with styled axes.
%   [fig, ax] = sharcFigure('single')            one panel, half text width
%   [fig, ax] = sharcFigure('full', [1 3])        one row of three panels, full text width
%   [fig, ax] = sharcFigure('full', [2 2], 12)    2 x 2 panels, 12 cm tall
%   [fig, ax] = sharcFigure('single', [1 1], [], 'left', 2.6)   override a margin (cm) for this figure
%
%   width     'single' | 'full' | a number in cm
%   nPanels   [rows cols] (default [1 1]); ax is returned row-major, ax(k) is panel k
%   heightCm  default: panel height = panel width x style aspect, plus the fixed margins
%   Axes are placed with fixed centimetre margins (sharcStyle.margin), so layouts are identical
%   across figures.
%
%   Axes come pre-styled (sharcAxes) and with hold on. Label panels with sharcPanel(ax(k), 'a').
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.

s = sharcStyle();
for i = 1:2:numel(varargin), s.margin.(varargin{i}) = varargin{i + 1}; end
if nargin < 1 || isempty(width), width = 'single'; end
if nargin < 2 || isempty(nPanels), nPanels = [1 1]; end
if ischar(width), W = s.width.(width); else, W = width; end
if nargin < 3 || isempty(heightCm)
    pw0 = (W - s.margin.left - s.margin.right - (nPanels(2) - 1) * s.margin.hgap) / nPanels(2);
    heightCm = nPanels(1) * pw0 * s.aspect + s.margin.bottom + s.margin.top + (nPanels(1) - 1) * s.margin.vgap;
end
vis = 'on';
if exist('OCTAVE_VERSION', 'builtin') && isempty(getenv('DISPLAY')), vis = 'off'; end   % headless tests
fig = figure('Color', 'w', 'Visible', vis, 'Units', 'centimeters', 'Position', [2 2 W heightCm], ...
    'PaperUnits', 'centimeters', 'PaperSize', [W heightCm], 'PaperPosition', [0 0 W heightCm], ...
    'InvertHardcopy', 'off');
ax = gobjects0(prod(nPanels));
m = s.margin;
nr = nPanels(1); nc = nPanels(2);
pw = (W - m.left - m.right - (nc - 1) * m.hgap) / nc;              % panel size, cm
ph = (heightCm - m.bottom - m.top - (nr - 1) * m.vgap) / nr;
for k = 1:prod(nPanels)
    r = ceil(k / nc); c = k - (r - 1) * nc;
    x = m.left + (c - 1) * (pw + m.hgap);
    y = heightCm - m.top - r * ph - (r - 1) * m.vgap;
    ax(k) = axes('Parent', fig, 'Units', 'centimeters', 'Position', [x y pw ph]);
    sharcAxes(ax(k));
end
end

function h = gobjects0(n)
if exist('gobjects', 'builtin') || exist('gobjects', 'file')
    h = gobjects(1, n);
else
    h = zeros(1, n);                                   % Octave
end
end
