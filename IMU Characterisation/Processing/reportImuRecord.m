function rep = reportImuRecord(file, varargin)
%REPORTIMURECORD  One-command first look at a SHARC IMU stream log.
%   rep = reportImuRecord(file)
%   rep = reportImuRecord(file, 'NominalRate', 20, 'SaveDir', 'Results', 'Show', true)
%
%   Runs, for one raw PuTTY log:
%     1  loadImuRecord: parser summary (rows per sensor, malformed lines, Zephyr drop notices)
%     2  checkImuTimestamps for every sensor: received rate (N-1)/span, the tick rate that
%        checkImuTimestamps estimates after filling inferred slots, median / 99th-percentile / max
%        interval, timing residual. Without a sample counter, long intervals are only
%        'gap-like intervals' and 'inferred missing nominal slots': a lost iteration and a slow
%        loop iteration look the same. Zephyr drop notices during streaming are the direct
%        evidence of logging loss. The timestamps are
%        LOGGING-PATH times (MCU uptime when the log line was created) until the firmware
%        confirms they are acquisition times; they are labelled as such everywhere.
%     3  per-sensor channel statistics: mean and std of ax ay az gx gy gz, mean |a|
%     4  pairwise comparison of the co-located sensors: rows matched tick by tick (nearest
%        time within half a nominal period), correlation of the de-meaned signals per axis.
%        High correlation = common motion of the bench, not independent sensor noise.
%     5  figures (time series per sensor; dt histogram) and a text summary, saved to SaveDir
%
%   Name-value options
%     'NominalRate'  commanded rate (Hz). Default: parsed from the 'imu stream on <rate>'
%                    command in the log, else from a '_<rate>Hz' tag in the file name.
%     'SaveDir'      output folder (default: ../Results next to the Data folder, i.e.
%                    IMU Characterisation/Results; created if missing). '' = do not save.
%     'Show'         show figures (default true)
%
%   Output rep: file, nominalRate, rec (from loadImuRecord), timing(k) (from
%   checkImuTimestamps), rows(k) (one-line timing summary per sensor, used by
%   runShellRateSweep), stats(k), corr (pairs, axes, R per axis), summary (text).
%
% Lives in IMU Characterisation/Processing/.

opt = struct('NominalRate', [], 'SaveDir', [], 'Show', true);
for i = 1:2:numel(varargin)
    if ~isfield(opt, varargin{i}), error('reportImuRecord:args', 'Unknown option ''%s''.', varargin{i}); end
    opt.(varargin{i}) = varargin{i + 1};
end
[fdir, fname] = fileparts(file);
if isempty(opt.SaveDir)
    opt.SaveDir = fullfile(fdir, '..', '..', 'Results');
end

rec = loadImuRecord(file, 'Verbose', false);
S = rec.sensors;
if isempty(S), error('reportImuRecord:empty', 'No IMU data lines found in %s.', file); end

% ---- nominal rate ------------------------------------------------------------------------
fs = opt.NominalRate; src = 'option';
if isempty(fs)
    for c = numel(rec.meta.streamCommands):-1:1
        tk = regexp(rec.meta.streamCommands{c}, 'imu stream on\s+(\d+(\.\d+)?)', 'tokens', 'once');
        if ~isempty(tk), fs = str2double(tk{1}); src = 'log command'; break; end
    end
end
if isempty(fs)
    tk = regexp(fname, '_(\d+(\.\d+)?)Hz', 'tokens', 'once');
    if ~isempty(tk), fs = str2double(tk{1}); src = 'file name'; end
end
if isempty(fs)
    fs = 1 / median(diff(S(1).t)); src = 'median interval (no command or tag found)';
end

% ---- timing and statistics -------------------------------------------------------------------
L = {};
L{end + 1} = sprintf('reportImuRecord  %s', file);
L{end + 1} = sprintf('commanded rate %.4g Hz (from %s); timestamps are LOGGING-PATH times (apparent resolution below)', fs, src);
L{end + 1} = sprintf('%d lines, %d malformed IMU lines skipped, Zephyr messages dropped: %d during streaming (%d in total incl. boot)', ...
    rec.meta.nLines, rec.meta.nMalformed, rec.meta.droppedDuringStream, rec.meta.droppedMessages);
L{end + 1} = '';
L{end + 1} = sprintf('%-12s %6s %8s %9s %8s %9s %8s %8s %8s %8s %6s %6s %7s', 'sensor', 'rows', 'span s', 'f_recv', ...
    'vs cmd%', 'f_tick', 'med dt', 'p99 dt', 'max dt', 'resid', 'gapLk', 'slots', 'res ms');
L{end + 1} = sprintf('%-12s %6s %8s %9s %8s %9s %8s %8s %8s %8s %6s %6s %7s', '', '', '', '(Hz)', '', '(Hz, est)', '(ms)', ...
    '(ms)', '(ms)', '(ms)', '', '', '');
timing = struct([]); stats = struct([]); rows = struct([]);
for k = 1:numel(S)
    o = checkImuTimestamps(S(k).t, 'NominalRate', fs, 'Verbose', false);
    timing = [timing, o]; %#ok<AGROW>
    dts = 1e3 * diff(S(k).t);
    fRecv = (S(k).nRows - 1) / (S(k).t(end) - S(k).t(1));          % lines actually received per second
    rw = struct('sensor', S(k).name, 'fCmd', fs, 'rows', S(k).nRows, 'span', S(k).t(end) - S(k).t(1), ...
        'receivedRate', fRecv, 'tickRateEst', o.rate, 'ratio', fRecv / fs, 'medDt', median(dts), 'p99Dt', pctl(dts, 99), 'maxDt', max(dts), ...
        'residMs', 1e3 * o.jitterStd, 'gapLike', o.nGaps, 'inferredSlots', o.nMissing, ...
        'resolutionMs', 1e3 * o.apparentResolution, 'malformed', rec.meta.nMalformed, ...
        'totalParsedRows', sum([S.nRows]), ...
        'dropsDuringStream', rec.meta.droppedDuringStream);
    rows = [rows, rw]; %#ok<AGROW>
    L{end + 1} = sprintf('%-12s %6d %8.1f %9.4f %+8.2f %9.4f %8.1f %8.1f %8.1f %8.2f %6d %6d %7.3g', S(k).name, S(k).nRows, ...
        rw.span, rw.receivedRate, 100 * (rw.ratio - 1), rw.tickRateEst, rw.medDt, rw.p99Dt, rw.maxDt, rw.residMs, ...
        rw.gapLike, rw.inferredSlots, rw.resolutionMs); %#ok<AGROW>
    X = [S(k).acc, S(k).gyr];
    st = struct('name', S(k).name, 'mean', mean(X, 1), 'std', std(X, 0, 1), ...
        'aMag', mean(sqrt(sum(S(k).acc.^2, 2))));
    stats = [stats, st]; %#ok<AGROW>
end
L{end + 1} = '(f_recv = lines received per second, (N-1)/span; f_tick = underlying tick rate estimated by checkImuTimestamps';
L{end + 1} = ' after filling inferred slots, meaningful only while losses are rare; gapLk = gap-like intervals, slots = inferred';
L{end + 1} = ' missing nominal slots: no sample counter, so NOT proven losses)';
L{end + 1} = '';
L{end + 1} = sprintf('%-12s %8s | %-26s | %-26s | %-26s | %-26s', 'sensor', '|a|', 'mean a (m/s^2)', 'std a (m/s^2)', ...
    'mean g (rad/s)', 'std g (rad/s)');
for k = 1:numel(stats)
    s = stats(k);
    L{end + 1} = sprintf('%-12s %8.4f | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f | %8.5f %8.5f %8.5f', ...
        s.name, s.aMag, s.mean(1:3), s.std(1:3), s.mean(4:6), s.std(4:6)); %#ok<AGROW>
end
L{end + 1} = '(gyro std at or below the 0.001 rad/s print step is limited by the printout, not the sensor)';

% ---- pairwise correlation of co-located sensors ------------------------------------------------------
axesNames = {'ax', 'ay', 'az', 'gx', 'gy', 'gz'};
pairs = nchoosek(1:numel(S), 2);
R = NaN(size(pairs, 1), 6); nMatch = zeros(size(pairs, 1), 1);
for p = 1:size(pairs, 1)
    i = pairs(p, 1); j = pairs(p, 2);
    [ii, jj] = matchTicks(S(i).t, S(j).t, 0.5 / fs);
    nMatch(p) = numel(ii);
    if nMatch(p) < 10, continue; end
    Xi = [S(i).acc(ii, :), S(i).gyr(ii, :)]; Xj = [S(j).acc(jj, :), S(j).gyr(jj, :)];
    for c = 1:6
        a = Xi(:, c) - mean(Xi(:, c)); b = Xj(:, c) - mean(Xj(:, c));
        R(p, c) = (a' * b) / sqrt((a' * a) * (b' * b));
    end
end
L{end + 1} = '';
L{end + 1} = sprintf('PAIRWISE CORRELATION (matched ticks; ~1 = common bench motion, ~0 = independent noise)');
L{end + 1} = sprintf('%-26s %6s %s', 'pair', 'ticks', sprintf('%7s', axesNames{:}));
for p = 1:size(pairs, 1)
    L{end + 1} = sprintf('%-26s %6d %s', [S(pairs(p, 1)).name, ' / ', S(pairs(p, 2)).name], nMatch(p), ...
        sprintf('%7.2f', R(p, :))); %#ok<AGROW>
end
L{end + 1} = '(correlation is also limited where the printout quantises the signal, e.g. gyro z)';
summary = strjoin(L, sprintf('\n'));
fprintf('%s\n', summary);

rep = struct('file', file, 'nominalRate', fs, 'rateSource', src, 'rec', rec, 'timing', {timing}, 'rows', {rows}, 'stats', {stats}, ...
    'corr', struct('pairs', pairs, 'axes', {axesNames}, 'R', R, 'nMatched', nMatch), 'summary', summary);

% ---- figures -----------------------------------------------------------------------------------------------
vis = 'on'; if ~opt.Show, vis = 'off'; end
f1 = figure('Color', 'w', 'Visible', vis, 'Position', [50 50 1300 250 * numel(S)]);
t0 = min(arrayfun(@(s) s.t(1), S));
for k = 1:numel(S)
    ax = subplot(numel(S), 2, 2 * k - 1);
    plot(ax, S(k).t - t0, S(k).acc - mean(S(k).acc, 1)); grid(ax, 'on');
    ylabel(ax, 'a - mean (m/s^2)'); title(ax, sprintf('%s accelerometer (de-meaned)', S(k).name), 'Interpreter', 'none');
    if k == 1, legend(ax, 'x', 'y', 'z'); end
    ax = subplot(numel(S), 2, 2 * k);
    plot(ax, S(k).t - t0, S(k).gyr - mean(S(k).gyr, 1)); grid(ax, 'on');
    ylabel(ax, 'g - mean (rad/s)'); title(ax, sprintf('%s gyroscope (de-meaned)', S(k).name), 'Interpreter', 'none');
end
xlabel(ax, 'logging-path time (s)');
f2 = figure('Color', 'w', 'Visible', vis, 'Position', [80 80 700 400]); hold on;
allDt = cell2mat(arrayfun(@(s) 1e3 * diff(s.t), S(:), 'UniformOutput', false));
edges = floor(min([allDt; 1e3 / fs])) - 1.5:1:ceil(max(allDt)) + 1.5;
cmax = 1;
for k = 1:numel(S)
    cnt = histc(1e3 * diff(S(k).t), edges); cnt = max(cnt, 0.5);     % floor so the log axis shows empty bins
    stairs(edges, cnt, 'LineWidth', 1.4); cmax = max(cmax, max(cnt));
end
set(gca, 'YScale', 'log'); ylim([0.5, 2 * cmax]); xlim(edges([1 end]));
plot(1e3 / fs * [1 1], [0.5, 2 * cmax], 'k--'); grid on;           % log counts: rare long intervals stay visible
xlabel('\Delta t between lines (ms, logging-path time)'); ylabel('count');
legend([{S.name}, {'commanded'}], 'Interpreter', 'none'); title(sprintf('%s: interval distribution', fname), 'Interpreter', 'none');

% ---- save -------------------------------------------------------------------------------------------------------
if ~isempty(opt.SaveDir)
    if ~exist(opt.SaveDir, 'dir'), mkdir(opt.SaveDir); end
    base = fullfile(opt.SaveDir, ['report_', fname]);
    fid = fopen([base, '.txt'], 'w'); fprintf(fid, '%s\n', summary); fclose(fid);
    save([base, '.mat'], 'rep');
    saveFig(f1, [base, '_series.png']); saveFig(f2, [base, '_dt.png']);
    fprintf('Saved: %s.txt / .mat / _series.png / _dt.png\n', base);
end
end

function p = pctl(x, q)
x = sort(x(:)); p = x(max(1, min(numel(x), ceil(q / 100 * numel(x)))));
end

function [ii, jj] = matchTicks(ti, tj, tol)
% nearest-neighbour match of tj to ti within tol (both sorted)
k = interp1(tj, 1:numel(tj), ti, 'nearest', 'extrap');
k = max(1, min(numel(tj), round(k)));
ok = abs(tj(k) - ti) <= tol;
ii = find(ok); jj = k(ok);
[jj, u] = unique(jj, 'stable'); ii = ii(u);
end

function saveFig(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 150);
catch
    print(fig, file, '-dpng', '-r150');
end
end