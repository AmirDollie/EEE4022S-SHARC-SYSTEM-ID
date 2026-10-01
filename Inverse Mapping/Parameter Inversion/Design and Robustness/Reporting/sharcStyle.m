function s = sharcStyle()
%SHARCSTYLE  The single figure style for every Design-and-Robustness study (thesis figures).
%   s = sharcStyle()
%
%   One place for colours, markers, line styles, fonts and sizes, so every figure in the results
%   chapter looks the same. Change a value here, rerun the study scripts, and every figure follows.
%
%   COLOUR RULES
%   - s.series(k, :) is categorical slot k, used in FIXED order (slot 1 for the first entity, slot 2
%     for the second, ...), never cycled. The same entity keeps the same slot in every figure: use
%     s.sensorColour(name) for the noise cases so a sensor is always the same colour.
%     The palette passes a colour-vision-deficiency check for adjacent slots, but three slots are
%     below 3:1 contrast on white, so every series ALSO gets its own marker and a legend entry:
%     identity is never carried by colour alone (also needed for greyscale printing).
%   - s.class.(...) are reserved for the identifiability classes and never used for a series.
%   - Thresholds (sigma = 1%, d = 0.5, d = 1, ...) are dashed lines in s.threshold with a text label.
%
%   DATA SERIES: always plot with sharcSeries(ax, x, y, slot) or sharcSensorSeries(ax, x, y, name),
%   never with plot(..., 'Color', ...) directly, so this file is the only place the style is set.
%
%   SIZES (cm): s.width.single (half text width, one panel), s.width.full (full text width).
%   Export: saveSharcFigure writes a vector PDF (for LaTeX) and a 300 dpi PNG (for quick viewing).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.

hex = @(h) [hex2dec(h(2:3)), hex2dec(h(4:5)), hex2dec(h(6:7))] / 255;
slots = {'#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948'};
s.series = cell2mat(cellfun(hex, slots(:), 'UniformOutput', false));
s.markers = {'o', 's', '^', 'd', 'v', '>', '<', 'p'};
s.lineStyles = {'-', '--', '-.', ':', '-', '--', '-.', ':'};

s.ink = hex('#0b0b0b');            % titles, labels, annotations
s.axis = hex('#52514e');           % axis lines and tick labels
s.grid = hex('#e1e0d9');           % hairline grid
s.threshold = hex('#898781');      % dashed threshold lines
s.fill = hex('#f0efec');           % shaded regions (e.g. inadmissible combinations)

% identifiability classes (status colours, reserved)
s.class.Identifiable = hex('#0ca30c');
s.class.Marginal = hex('#fab219');
s.class.NotIdentifiable = hex('#d03b3b');

% fixed sensor-to-slot mapping (cfg.noise.names order), so a sensor is the same colour everywhere
s.sensorNames = {'datasheet', 'fused3', 'lsm6dsv16x', 'g3', 'icm42688p', 'lsm6ds3trc'};
s.sensorLabels = {'Datasheet grade', 'Three IMUs fused', 'LSM6DSV16X', 'G3 assumption', 'ICM-42688-P', ...
    'LSM6DS3TR-C'};
s.sensorSlot = [7 3 1 8 2 4];
s.sensorColour = @(name) s.series(s.sensorSlot(strcmp(s.sensorNames, name)), :);
s.sensorMarker = @(name) s.markers{s.sensorSlot(strcmp(s.sensorNames, name))};
s.sensorLabel = @(name) s.sensorLabels{strcmp(s.sensorNames, name)};

s.font = 'Helvetica';               % sans-serif figures throughout (decided; do not change mid-thesis)
s.fontSize = 9;                    % axis labels, ticks
s.smallFontSize = 8;               % legends, annotations
s.lineWidth = 1.3;                 % data lines
s.thinLineWidth = 0.8;             % thresholds, references
s.axisLineWidth = 0.6;
s.markerSize = 5;

s.width.single = 8.4;              % cm
s.width.full = 16.5;               % cm
s.aspect = 0.68;                   % default height / width for one panel row
% fixed margins (cm) so every figure has the same layout in MATLAB and Octave
s.margin = struct('left', 1.45, 'right', 0.45, 'bottom', 1.05, 'top', 0.55, 'hgap', 1.6, 'vgap', 1.25);
s.dpi = 300;
end
