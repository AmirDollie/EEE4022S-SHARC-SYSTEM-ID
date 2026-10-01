function sharcTicks(ax, axisName, ticks, fmt)
%SHARCTICKS  Set ticks with plain-number labels (no 10^n or 2x10^1 notation on log axes).
%   sharcTicks(ax, 'x', [10 30 60 109 240])
%   sharcTicks(ax, 'y', [0.01 0.1 0.5 1 10], '%g')
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
if nargin < 4, fmt = '%g'; end
ticks = sort(ticks(:).');
labels = arrayfun(@(v) sprintf(fmt, v), ticks, 'UniformOutput', false);
A = upper(axisName);
set(ax, [A 'Tick'], ticks, [A 'TickLabel'], labels);
try, set(ax, [A 'MinorTick'], 'off'); catch, end
end
