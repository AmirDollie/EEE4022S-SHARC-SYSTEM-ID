%% runAccelCalibration.m
% H1a: static accelerometer calibration of the three IMUs from the hold logs in Data/Raw/Calibration
% (imu_cal_NN_<orientation>.log, one PuTTY log per static hold, 10 Hz).
%
% Per hold: trim the ends, check that the board was still (isStaticHold: magnitude, vector, drift and
% gyro criteria), average. Per sensor: PRIMARY fit a_cal = diag(s) (a_raw - b) with ||a_cal|| = g
% (calibrateAccelSixPosition 'diag', 6 parameters) whatever the number of holds; the full model with
% non-orthogonality is run only as a DIAGNOSTIC when >= 12 holds exist. Plus the two-point flip
% estimate on z from the 'up' and 'down' holds. Robustness: leave-one-hold-out refits (>= 8 holds).
%
% Headline for the thesis: the z-axis gain of each sensor, the RELATIVE z gain between the raw
% sensors, and the calibration's own scale uncertainty against the 0.42% tolerance. This is NOT a
% measured post-calibration gain match: that claim belongs to H1b (co-located transmissibility).
% Offsets are reported but do not enter the transmissibility (0 Hz, removed by detrending).
%
% OUTPUT (Results/): accelCalibration_<stamp>.txt / .mat / .png
%
% Lives in IMU Characterisation/Processing/.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
clearvars -except thisDir ACCEL_CAL_CONFIG
close all

%% ---- Configuration ------------------------------------------------------------------------------
cfg.DATA_DIR = fullfile(thisDir, '..', 'Data', 'Raw', 'Calibration');
cfg.PATTERN = 'imu_cal_*.log';
cfg.EXCLUDE = 'shake';                % dynamic logs are not static holds
cfg.RESULTS_DIR = fullfile(thisDir, '..', 'Results');
cfg.G = 9.796;                        % m/s^2, Cape Town
cfg.TRIM_S = 10;                      % s removed from each end of a hold
cfg.STILL = {'MagStd', 0.05, 'VecStd', 0.08, 'GyroRms', 0.01, 'AvgLoss', 1e-4, 'Drift', 0.02};   % isStaticHold
cfg.FULL_DIAGNOSTIC_MIN_HOLDS = 12;
cfg.TOL_REL_GAIN = 0.0042;            % measurement-tolerance study
cfg.SHOW = true;
if exist('ACCEL_CAL_CONFIG', 'var') && isstruct(ACCEL_CAL_CONFIG)
    fn = fieldnames(ACCEL_CAL_CONFIG);
    for i = 1:numel(fn), cfg.(fn{i}) = ACCEL_CAL_CONFIG.(fn{i}); end
end
if ~exist(cfg.RESULTS_DIR, 'dir'), mkdir(cfg.RESULTS_DIR); end
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(cfg.RESULTS_DIR, ['accelCalibration_', stamp]);

%% ---- 1  Hold means ------------------------------------------------------------------------------
files = dir(fullfile(cfg.DATA_DIR, cfg.PATTERN));
files = files(cellfun(@isempty, strfind({files.name}, cfg.EXCLUDE)));
[~, o] = sort({files.name}); files = files(o);
nH = numel(files);
if nH == 0, error('runAccelCalibration:noFiles', 'No %s in %s.', cfg.PATTERN, cfg.DATA_DIR); end
holdNames = regexprep({files.name}, '^imu_cal_\d+_|\.log$', '');
names = {};
Amean = []; Astd = []; nUsed = []; still = [];
for h = 1:nH
    rec = loadImuRecord(fullfile(cfg.DATA_DIR, files(h).name), 'Verbose', false);
    if h == 1
        names = {rec.sensors.name}; nS = numel(names);
        Amean = NaN(nH, 3, nS); Astd = NaN(nH, nS); nUsed = zeros(nH, nS); still = false(nH, nS);
        stillInfo = cell(nH, nS);
    end
    for k = 1:nS
        idx = find(strcmp({rec.sensors.name}, names{k}), 1);
        if isempty(idx), continue; end
        s = rec.sensors(idx);
        sel = s.t >= s.t(1) + cfg.TRIM_S & s.t <= s.t(end) - cfg.TRIM_S;
        a = s.acc(sel, :);
        Amean(h, :, k) = mean(a, 1);
        Astd(h, k) = std(sqrt(sum(a.^2, 2)));
        nUsed(h, k) = nnz(sel);
        [still(h, k), stillInfo{h, k}] = isStaticHold(a, s.gyr(sel, :), cfg.STILL{:});
    end
end

%% ---- 2  Fits --------------------------------------------------------------------------------------
iUp = find(strcmpi(holdNames, 'up'), 1); iDown = find(strcmpi(holdNames, 'down'), 1);
model = 'diag';                                           % primary, whatever the number of holds
runFull = nH >= cfg.FULL_DIAGNOSTIC_MIN_HOLDS;
cal = cell(1, nS); calFull = cell(1, nS);
for k = 1:nS
    ok = still(:, k) & all(isfinite(Amean(:, :, k)), 2);
    opts = {'Model', model};
    if ~isempty(iUp) && ~isempty(iDown) && ok(iUp) && ok(iDown)
        opts = [opts, {'TwoPointAxis', [nnz(ok(1:iUp)) nnz(ok(1:iDown)) 3]}]; %#ok<AGROW>
    end
    cal{k} = calibrateAccelSixPosition(Amean(ok, :, k), cfg.G, opts{:});
    cal{k}.holdsUsed = find(ok).';
    if runFull && nnz(ok) >= cfg.FULL_DIAGNOSTIC_MIN_HOLDS
        calFull{k} = calibrateAccelSixPosition(Amean(ok, :, k), cfg.G, 'Model', 'full');
    end
end

%% ---- 3  Relative z gain between sensors ----------------------------------------------------------
gz = cellfun(@(c) 1 / c.scale(3), cal);                  % reading per unit true acceleration, z
gzRaw = NaN(1, nS);
for k = 1:nS, if isfield(cal{k}, 'twoPoint'), gzRaw(k) = cal{k}.twoPoint.gain; end, end
pairs = nchoosek(1:nS, 2);
relFit = zeros(size(pairs, 1), 1); relTwo = relFit; relScaleUncertainty = relFit;
for p = 1:size(pairs, 1)
    i = pairs(p, 1); j = pairs(p, 2);
    relFit(p) = gz(i) / gz(j) - 1;                       % uncalibrated mismatch, from the fit
    relTwo(p) = gzRaw(i) / gzRaw(j) - 1;                 % uncalibrated mismatch, two-point
    % calibration scale uncertainty of the pair (NOT a measured post-calibration match, not a 95%
    % bound): the larger of the leave-one-out half-range and the formal hold-noise floor, per sensor
    % (max ignores NaN, which is what the half-range is when leave-one-out is unavailable)
    ui = max([cal{i}.looSpreadScale(3), cal{i}.formalStdScale(3)]) / cal{i}.scale(3);
    uj = max([cal{j}.looSpreadScale(3), cal{j}.formalStdScale(3)]) / cal{j}.scale(3);
    relScaleUncertainty(p) = hypot(ui, uj);
end

%% ---- 4  Summary ------------------------------------------------------------------------------------
L = {};
L{end + 1} = sprintf('runAccelCalibration  %s', stamp);
L{end + 1} = sprintf('%d holds from %s; g = %.3f m/s^2; trim %g s; primary model ''%s''%s', nH, cfg.DATA_DIR, ...
    cfg.G, cfg.TRIM_S, model, ternary(runFull, ', full model as diagnostic', ''));
L{end + 1} = sprintf(['stillness (exclude): std|a| <= %g m/s^2, vector std <= %g m/s^2, gyro RMS <= %g rad/s, ' ...
    'averaging loss <= %g; (warn): drift <= %g m/s^2'], cfg.STILL{2:2:end});
L{end + 1} = '';
L{end + 1} = 'HOLDS (mean raw acceleration, m/s^2; std of |a| in the hold)';
for h = 1:nH
    row = sprintf('  %-10s', holdNames{h});
    for k = 1:nS
        row = [row, sprintf(' | %-11s [%7.3f %7.3f %7.3f] |a| %6.3f sd %.3f%s', names{k}, Amean(h, :, k), ...
            norm(Amean(h, :, k)), Astd(h, k), holdFlag(still(h, k), stillInfo{h, k}))]; %#ok<AGROW>
    end
    L{end + 1} = row; %#ok<SAGROW>
end
L{end + 1} = '';
L{end + 1} = 'PRIMARY FIT  a_cal = diag(s) (a_raw - b), ||a_cal|| = g';
L{end + 1} = sprintf('  %-12s %5s %4s %10s %10s %10s %12s %12s %12s %9s %8s', 'sensor', 'holds', 'dof', 'gainErr x', ...
    'gainErr y', 'gainErr z', 'offset x', 'offset y', 'offset z', 'rms res', 'cond');
L{end + 1} = sprintf('  %-12s %5s %4s %10s %10s %10s %12s %12s %12s %9s %8s', '', '', '', '(%)', '(%)', '(%)', '(mg)', ...
    '(mg)', '(mg)', '(mg)', '');
for k = 1:nS
    c = cal{k};
    L{end + 1} = sprintf('  %-12s %5d %4d %10.3f %10.3f %10.3f %12.1f %12.1f %12.1f %9.2f %8.1f', names{k}, ...
        numel(c.holdsUsed), c.dof, 100 * c.gainError, 1e3 * c.offset / 9.80665, 1e3 * c.rmsResidual / 9.80665, c.cond); %#ok<SAGROW>
    L{end + 1} = sprintf('  %-12s %5s %4s %10.3f %10.3f %10.3f %12.1f %12.1f %12.1f   (leave-one-out half-range)', '', '', ...
        '', 100 * c.looSpreadScale ./ c.scale, 1e3 * c.looSpreadOffset / 9.80665); %#ok<SAGROW>
    if ~c.loo.available
        L{end + 1} = sprintf('  %-12s leave-one-out unavailable: %d holds leave no residual in a 6-parameter refit (need >= 8)', ...
            '', numel(c.holdsUsed)); %#ok<SAGROW>
    elseif ~isempty(c.loo.essentialHolds)
        L{end + 1} = sprintf('  %-12s leave-one-out: %d of %d refits determined; without hold(s) %s the fit is undetermined', ...
            '', c.loo.nDetermined, numel(c.holdsUsed), strjoin(holdNames(c.holdsUsed(c.loo.essentialHolds)), ', ')); %#ok<SAGROW>
    end
end
L{end + 1} = '';
L{end + 1} = 'TWO-POINT FLIP, z axis (up / down holds)';
for k = 1:nS
    if isfield(cal{k}, 'twoPoint')
        L{end + 1} = sprintf('  %-12s gain error %+7.3f %%, offset %+7.1f mg   (fit: %+7.3f %%)', names{k}, ...
            100 * cal{k}.twoPoint.gainError, 1e3 * cal{k}.twoPoint.offset / 9.80665, 100 * cal{k}.gainError(3)); %#ok<SAGROW>
    end
end
L{end + 1} = '';
if runFull
    L{end + 1} = 'FULL-MODEL DIAGNOSTIC (upper-triangular M; non-orthogonality in mrad)';
    for k = 1:nS
        if isempty(calFull{k}), continue; end
        c = calFull{k};
        L{end + 1} = sprintf('  %-12s dof %d, rms res %.2f mg (diag: %.2f mg), off-diagonal [%.2f %.2f %.2f] mrad', names{k}, ...
            c.dof, 1e3 * c.rmsResidual / 9.80665, 1e3 * cal{k}.rmsResidual / 9.80665, 1e3 * [c.M(1, 2) / c.M(2, 2), ...
            c.M(1, 3) / c.M(3, 3), c.M(2, 3) / c.M(3, 3)]); %#ok<SAGROW>
    end
    L{end + 1} = '';
end
L{end + 1} = sprintf('RELATIVE z GAIN BETWEEN SENSORS (tolerance %.2f %%)', 100 * cfg.TOL_REL_GAIN);
L{end + 1} = sprintf('  %-26s %14s %14s %24s', 'pair', 'raw (fit)', 'raw (2-point)', 'calibration scale unc.');
for p = 1:size(pairs, 1)
    L{end + 1} = sprintf('  %-26s %+13.3f%% %+13.3f%% %21.3f%%   %s', [names{pairs(p, 1)}, '/', names{pairs(p, 2)}], ...
        100 * relFit(p), 100 * relTwo(p), 100 * relScaleUncertainty(p), verdict(relScaleUncertainty(p), cfg.TOL_REL_GAIN, ...
        thinGeometry(cal{pairs(p, 1)}) || thinGeometry(cal{pairs(p, 2)}))); %#ok<SAGROW>
end
L{end + 1} = '';
L{end + 1} = 'NOTES';
L{end + 1} = ' - gain error = reading / true - 1 (positive reads high). Offsets do not affect the transmissibility.';
L{end + 1} = ' - "calibration scale unc." is the uncertainty of the calibration itself (leave-one-out half-range or';
L{end + 1} = '   hold-noise floor), not a 95% bound and not a measured post-calibration match: H1b makes that claim.';
L{end + 1} = ' - a sphere fit fixes axis lengths, not directions: misalignment to the board frame is not estimated.';
L{end + 1} = ' - dof = holds - 6 parameters. With dof = 1 the residual RMS is weak evidence and leave-one-out is';
L{end + 1} = '   unavailable (needs >= 8 holds), so the uncertainty is the hold-noise floor only. Add tilted holds';
L{end + 1} = '   (about 45 deg, three roll angles) to test the model; the full model needs >= 12 holds.';
L{end + 1} = ' - the two-point z estimate assumes the up/down holds were aligned with z; a tilt of theta biases';
L{end + 1} = '   it by about (1 - cos theta), 0.06% at 2 deg. The fit does not need aligned holds.';
summary = strjoin(L, sprintf('\n'));
fprintf('%s\n', summary);
fid = fopen([baseName, '.txt'], 'w'); fprintf(fid, '%s\n', summary); fclose(fid);
out = struct('cfg', cfg, 'files', {{files.name}}, 'holdNames', {holdNames}, 'sensors', {names}, ...
    'Amean', Amean, 'Astd', Astd, 'nUsed', nUsed, 'still', still, 'cal', {cal}, 'pairs', pairs, ...
    'relFit', relFit, 'relTwoPoint', relTwo, 'relScaleUncertainty', relScaleUncertainty, 'calFull', {calFull}, ...
    'stillInfo', {stillInfo}, 'summary', summary);
save([baseName, '.mat'], 'out', '-v7');

%% ---- 5  Figure: |a| per hold before and after calibration ---------------------------------------
fig = figure('Color', 'w', 'Visible', ternary(cfg.SHOW, 'on', 'off'), 'Position', [60 60 1100 380]);
ax = axes('Parent', fig); hold(ax, 'on');
mk = {'o', 's', '^'};
for k = 1:nS
    raw = sqrt(sum(Amean(:, :, k).^2, 2));
    c = cal{k}; calN = NaN(nH, 1);
    calN(c.holdsUsed) = sqrt(sum(((Amean(c.holdsUsed, :, k) - c.offset) * c.M.').^2, 2));
    plot(ax, 1:nH, 1e3 * (raw - cfg.G) / 9.80665, ['--', mk{k}], 'DisplayName', [names{k}, ' raw']);
    plot(ax, 1:nH, 1e3 * (calN - cfg.G) / 9.80665, ['-', mk{k}], 'LineWidth', 1.4, 'DisplayName', [names{k}, ' calibrated']);
end
set(ax, 'XTick', 1:nH, 'XTickLabel', holdNames); grid(ax, 'on');
ylabel(ax, '|a| - g  (mg)'); title(ax, 'Static holds: magnitude error before and after calibration');
legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none');
try
    exportgraphics(fig, [baseName, '.png'], 'Resolution', 150);
catch
    print(fig, [baseName, '.png'], '-dpng', '-r150');
end
if ~cfg.SHOW, close(fig); end
fprintf('Saved: %s.txt / .mat / .png\n', baseName);

%% ================================================================================================
function s = verdict(u, tol, thin)
% The uncertainty is only quotable if the hold geometry supports it.
if thin
    s = 'noise floor only: too few holds to test the model';
elseif u <= tol
    s = 'calibration scale uncertainty below tolerance';
else
    s = 'calibration scale uncertainty EXCEEDS tolerance';
end
end

function s = holdFlag(still, info)
if ~still
    s = [' EXCLUDED (', strjoin(info.fails, ','), ')'];
elseif ~isempty(info.warnings)
    s = sprintf(' warn: %s %.3f m/s^2, kept (averaging loss %.1e)', strjoin(info.warnings, ','), info.drift, info.avgLoss);
else
    s = '';
end
end

function tf = thinGeometry(c)
tf = ~c.loo.available || ~isempty(c.loo.essentialHolds);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end