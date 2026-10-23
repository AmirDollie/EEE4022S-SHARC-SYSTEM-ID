function h = sharcSensorSeries(ax, x, y, sensorName, varargin)
%SHARCSENSORSERIES  Plot one series for a sensor noise case in that sensor's fixed colour and style.
%   h = sharcSensorSeries(ax, x, y, 'lsm6dsv16x')
%   h = sharcSensorSeries(ax, x, y, 'lsm6ds3trc', 'Marker', 'none')     overrides as in sharcSeries
%
%   sensorName is one of cfg.noise.names. The colour, line style, marker fill and legend label
%   (DisplayName) come from sharcStyle, so a sensor looks the same in every figure (R1, R3, R4, D2, ...).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
s = sharcStyle();
k = find(strcmp(s.sensorNames, sensorName), 1);
if isempty(k)
    error('sharcSensorSeries:name', 'Unknown sensor ''%s''; known: %s', sensorName, strjoin(s.sensorNames, ', '));
end
extra = {'DisplayName', s.sensorLabels{k}, 'LineStyle', s.sensorLineStyles{k}};
if s.sensorOpenMarker(k), extra = [extra, {'MarkerFaceColor', 'w'}]; end
h = sharcSeries(ax, x, y, s.sensorSlot(k), extra{:}, varargin{:});
end