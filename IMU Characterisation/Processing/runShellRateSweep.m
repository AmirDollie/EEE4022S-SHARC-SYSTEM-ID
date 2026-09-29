%% runShellRateSweep.m
% Shell-rate sweep: which commanded `imu stream on <rate>` does the Zephyr DIAGNOSTIC SHELL
% sustain, with three IMUs printed per tick, without appreciable message loss?
%
% This measures the shell/print path only. It says nothing about the IMUs' own output
% data rate (a firmware/ODR question): the result is to be worded as
%   "the Zephyr diagnostic shell sustained about X Hz for three simultaneously printed IMUs".
%
% INPUT  all logs in Data/Raw whose names end in _<rate>Hz.log (newest file per rate), or
%        an explicit FILES list. Each is analysed by reportImuRecord (figures and text per
%        log saved to Results/).
% OUTPUT Results/shellRateSweep_<stamp>.txt / .mat / .png:
%        one row per (commanded rate, sensor): rows, received rate f_recv, f_recv / f_cmd, tick-rate estimate, median / p99 / max
%        interval, timing residual, gap-like intervals, malformed lines, Zephyr drops during
%        streaming; and a verdict.
%
% VERDICT: two separate questions (a rate passes only if every sensor passes; edit below)
%   LOSS-FREE  Zephyr drop notices during streaming == 0, malformed lines <=
%              MALFORMED_MAX_FRAC of all attempted IMU lines in the log, and p99(dt) <= P99_MAX_PERIODS / f_cmd (no stalls)
%   ON-RATE    f_recv / f_cmd >= RATIO_MIN, with f_recv = (N-1)/span the lines actually received
%              per second (the tick-rate estimate f_tick fills inferred slots and overstates the
%              delivered rate once lines are lost; it is reported but not used)
%   Both also require every EXPECTED_SENSORS entry to be present, and the row counts to be
%   balanced (fewest / most >= ROW_BALANCE_MIN): a sensor that vanishes or is depleted
%   would otherwise leave the remaining rows passing on their own.
%   A sleep-based loop loses a fixed overhead per tick, so f_recv / f_cmd falls with rate even
%   with no loss at all; the effective rate is known from the timestamps, so a loss-free rate
%   below its command is still usable. The recommendation is the highest LOSS-FREE rate,
%   quoted with its effective rate.
% The inferred missing slots are reported but NOT used: without a sample counter a slow loop
% iteration and a lost one look the same.
%
% Lives in IMU Characterisation/Processing/.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
clearvars -except thisDir
close all

%% ---- Configuration -----------------------------------------------------------------------
DATA_DIR = fullfile(thisDir, '..', 'Data', 'Raw');
RESULTS_DIR = fullfile(thisDir, '..', 'Results');
FILES = {};                          % {} = auto: newest *_<rate>Hz.log per rate in DATA_DIR
RATIO_MIN = 0.95;
P99_MAX_PERIODS = 1.5;
MALFORMED_MAX_FRAC = 0.005;
ROW_BALANCE_MIN = 0.98;              % fewest rows / most rows across the sensors of one log
EXPECTED_SENSORS = {'LSM6DS3TR-C', 'LSM6DSV16X', 'ICM-42688-P'};

%% ---- Find the logs --------------------------------------------------------------------------
if isempty(FILES)
    d = dir(fullfile(DATA_DIR, '*Hz.log'));
    rates = NaN(1, numel(d));
    for i = 1:numel(d)
        tk = regexp(d(i).name, '_(\d+(\.\d+)?)Hz\.log$', 'tokens', 'once');
        if ~isempty(tk), rates(i) = str2double(tk{1}); end
    end
    keep = ~isnan(rates); d = d(keep); rates = rates(keep);
    [ur, ~, g] = unique(rates);
    FILES = cell(1, numel(ur));
    for r = 1:numel(ur)
        cand = find(g == r);
        [~, newest] = max([d(cand).datenum]);
        FILES{r} = fullfile(DATA_DIR, d(cand(newest)).name);
    end
end
if isempty(FILES), error('runShellRateSweep:noLogs', 'No *_<rate>Hz.log files in %s.', DATA_DIR); end
if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(RESULTS_DIR, ['shellRateSweep_', stamp]);

%% ---- Analyse each log ------------------------------------------------------------------------------
reps = cell(1, numel(FILES)); T = struct([]);
for i = 1:numel(FILES)
    fprintf('\n##### %s\n', FILES{i});
    reps{i} = reportImuRecord(FILES{i}, 'SaveDir', RESULTS_DIR, 'Show', false);
    T = [T, reps{i}.rows]; %#ok<AGROW>
end
[~, ord] = sort([T.fCmd]); T = T(ord);

%% ---- Table and verdict ----------------------------------------------------------------------------------
L = {};
L{end + 1} = sprintf('runShellRateSweep  (%s)', stamp);
L{end + 1} = 'Zephyr DIAGNOSTIC SHELL path, three IMUs printed per tick; timestamps are logging-path times.';
L{end + 1} = sprintf(['Criteria: f_recv/f_cmd >= %.2f, p99(dt) <= %.2f commanded periods, 0 drops during streaming, ' ...
    'malformed <= %.1f%% of attempted IMU lines, all %d sensors present, rows balanced >= %.2f'], ...
    RATIO_MIN, P99_MAX_PERIODS, 100 * MALFORMED_MAX_FRAC, numel(EXPECTED_SENSORS), ROW_BALANCE_MIN);
L{end + 1} = '';
L{end + 1} = sprintf('%5s %-12s %6s %9s %7s %9s %8s %8s %8s %7s %6s %6s %5s %5s  %s', 'cmd', 'sensor', 'rows', 'f_recv', ...
    'ratio', 'f_tick', 'med dt', 'p99 dt', 'max dt', 'resid', 'gapLk', 'slots', 'malf', 'drops', 'loss-free / on-rate');
L{end + 1} = sprintf('%5s %-12s %6s %9s %7s %9s %8s %8s %8s %7s', '(Hz)', '', '', '(Hz)', '', '(Hz, est)', '(ms)', '(ms)', '(ms)', '(ms)');
lossFree = false(size(T)); onRate = false(size(T));
for k = 1:numel(T)
    r = T(k);
    malformedFrac = r.malformed / (r.totalParsedRows + r.malformed);   % per log: all sensors' attempted lines
    lossFree(k) = r.dropsDuringStream == 0 && malformedFrac <= MALFORMED_MAX_FRAC ...
        && r.p99Dt <= 1e3 * P99_MAX_PERIODS / r.fCmd;
    onRate(k) = r.ratio >= RATIO_MIN;
    L{end + 1} = sprintf('%5g %-12s %6d %9.4f %7.3f %9.4f %8.1f %8.1f %8.1f %7.2f %6d %6d %5d %5d  %s / %s', r.fCmd, r.sensor, ...
        r.rows, r.receivedRate, r.ratio, r.tickRateEst, r.medDt, r.p99Dt, r.maxDt, r.residMs, r.gapLike, r.inferredSlots, r.malformed, ...
        r.dropsDuringStream, yesno(lossFree(k)), yesno(onRate(k))); %#ok<SAGROW>
end
L{end + 1} = '(malf and drops are per log, repeated on each sensor row; gapLk / slots are inferred, not proven losses)';
ratesTested = unique([T.fCmd]);
ok = false(size(ratesTested)); okRate = ok; allSensorsPresent = ok; balanced = ok; balance = NaN(size(ratesTested));
L{end + 1} = '';
for i = 1:numel(ratesTested)
    f = ratesTested(i);
    idx = [T.fCmd] == f;
    namesHere = {T(idx).sensor};
    allSensorsPresent(i) = numel(namesHere) == numel(EXPECTED_SENSORS) && all(ismember(EXPECTED_SENSORS, namesHere));
    balance(i) = min([T(idx).rows]) / max([T(idx).rows]);
    balanced(i) = balance(i) >= ROW_BALANCE_MIN;
    ok(i) = allSensorsPresent(i) && balanced(i) && all(lossFree(idx));
    okRate(i) = allSensorsPresent(i) && all(onRate(idx));
    L{end + 1} = sprintf('  %5g Hz: sensors %-3s  rows balanced %-3s (%.3f)  loss-free %-3s  on-rate %-3s  (received %.3f Hz)', ...
        f, yesno(allSensorsPresent(i)), yesno(balanced(i)), balance(i), yesno(ok(i)), yesno(okRate(i)), ...
        mean([T(idx).receivedRate])); %#ok<SAGROW>
end
if any(ok)
    best = max(ratesTested(ok));
    if any(~ok & ratesTested < best)
        L{end + 1} = sprintf('NOTE: a lower rate failed while %g Hz passed: check those logs (touched board? typing?).', best);
    end
    L{end + 1} = sprintf(['VERDICT: the Zephyr diagnostic shell sustained %g Hz commanded (received %.3f Hz per sensor) for ' ...
        'three simultaneously printed IMUs without appreciable shell-path loss under the stated criteria.'], best, ...
        mean([T([T.fCmd] == best).receivedRate]));
    if ~okRate(ratesTested == best)
        L{end + 1} = sprintf('  (below its command by %.1f%%: use the measured received rate, not the commanded one)', ...
            100 * (1 - mean([T([T.fCmd] == best).ratio])));
    end
else
    L{end + 1} = 'VERDICT: no tested rate was loss-free.';
end
L{end + 1} = 'This is the shell/print path limit, not the IMU output data rate.';
summary = strjoin(L, sprintf('\n'));
fprintf('\n%s\n', summary);
fid = fopen([baseName, '.txt'], 'w'); fprintf(fid, '%s\n', summary); fclose(fid);
out = struct('files', {FILES}, 'table', T, 'lossFree', lossFree, 'onRate', onRate, 'ratesTested', ratesTested, ...
    'rateLossFree', ok, 'rateOnRate', okRate, 'allSensorsPresent', allSensorsPresent, 'rowBalance', balance, ...
    'criteria', struct('RATIO_MIN', RATIO_MIN, 'P99_MAX_PERIODS', P99_MAX_PERIODS, 'MALFORMED_MAX_FRAC', MALFORMED_MAX_FRAC, ...
    'ROW_BALANCE_MIN', ROW_BALANCE_MIN, 'EXPECTED_SENSORS', {EXPECTED_SENSORS}), ...
    'summary', summary);
save([baseName, '.mat'], 'out');

%% ---- Figure -------------------------------------------------------------------------------------------------
names = unique({T.sensor}, 'stable');
mk = {'o', 's', '^'};
cols = {[0.85 0.33 0.10], [0.93 0.69 0.13], [0.49 0.18 0.56], [0.30 0.75 0.93], [0.47 0.67 0.19]};
fig = figure('Color', 'w', 'Position', [60 60 1150 420]);
ax1 = subplot(1, 2, 1); hold(ax1, 'on');
lim = [0, 1.1 * max(ratesTested)];
plot(ax1, lim, lim, 'k:');
for n = 1:numel(names)
    s = T(strcmp({T.sensor}, names{n}));
    plot(ax1, [s.fCmd], [s.receivedRate], ['-', mk{1 + mod(n - 1, 3)}], 'Color', cols{1 + mod(n - 1, 5)}, 'LineWidth', 1.4);
end
xlabel(ax1, 'commanded rate (Hz)'); ylabel(ax1, 'received rate per sensor (Hz)'); grid(ax1, 'on');
legend(ax1, [{'ideal'}, names], 'Location', 'northwest', 'Interpreter', 'none');
title(ax1, 'Received vs commanded rate (shell path)');
ax2 = subplot(1, 2, 2); hold(ax2, 'on'); h2 = zeros(1, numel(names));
for n = 1:numel(names)
    s = T(strcmp({T.sensor}, names{n}));
    c = cols{1 + mod(n - 1, 5)};
    h2(n) = plot(ax2, [s.fCmd], [s.medDt] .* [s.fCmd] / 1e3, ['-', mk{1 + mod(n - 1, 3)}], 'Color', c, 'LineWidth', 1.4); %#ok<SAGROW>
    plot(ax2, [s.fCmd], [s.p99Dt] .* [s.fCmd] / 1e3, ['--', mk{1 + mod(n - 1, 3)}], 'Color', c, 'LineWidth', 1.0);
end
plot(ax2, lim, [1 1], 'k:'); plot(ax2, lim, P99_MAX_PERIODS * [1 1], 'r:');
legend(ax2, h2, names, 'Location', 'northwest', 'Interpreter', 'none');
xlabel(ax2, 'commanded rate (Hz)'); ylabel(ax2, 'interval / commanded period'); grid(ax2, 'on');
title(ax2, 'Median (solid) and 99th-percentile (dashed) interval');
try
    exportgraphics(fig, [baseName, '.png'], 'Resolution', 150);
catch
    print(fig, [baseName, '.png'], '-dpng', '-r150');
end
fprintf('Saved: %s.txt / .mat / .png\n', baseName);

%% ================================================================================================
function s = yesno(tf)
if tf, s = 'yes'; else, s = 'no'; end
end