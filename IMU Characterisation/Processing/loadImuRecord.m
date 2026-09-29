function rec = loadImuRecord(file, varargin)
%LOADIMURECORD  Parse a SHARC board (Zephyr shell) IMU stream log into per-sensor arrays.
%   rec = loadImuRecord(file)
%   rec = loadImuRecord(file, 'Verbose', false)
%
%   INPUT FORMAT (PuTTY "All session output" log of `imu stream on <rate>`): Zephyr log
%   lines, ANSI colour codes and shell prompts interleaved, e.g.
%     uart:~$ [00:05:24.390,000] <inf> imu: ICM-42688-P  a   1.086  -0.086   9.849  g  -0.000  -0.005   0.000
%   timestamp [hh:mm:ss.mmm,uuu] = MCU uptime when the log message was created (NOT the IMU
%   sample time; 1 ms resolution in the first log); a = accelerometer (m/s^2); g = gyroscope
%   (rad/s, Zephyr SI convention: confirm with the firmware).
%
%   ROBUSTNESS: shell echo and log output can interleave inside a line (seen in the first log:
%   '1.0801.7438', ',001.74300', 'L<ESC>[DSM6DSV16X'). Only lines matching the strict pattern
%   (3-digit ms and us fields, known sensor name, signed decimal numbers) are kept; everything
%   else is counted in rec.meta. Zephyr's '--- N messages dropped ---' notices are counted too.
%   There is no sample counter, so a missing line cannot be told from a slow loop iteration.
%
%   OUTPUT rec.sensors(k): name, t (s, device uptime), acc (N x 3), gyr (N x 3), nRows
%          rec.meta: file, nLines, nMalformed, malformed (first 20 lines), droppedMessages (all),
%                    droppedDuringStream (notices after the first IMU data line: the ones that matter),
%                    (sum of Zephyr drop notices), notices (non-data log lines: boot, warnings,
%                    errors), streamCommands (lines containing 'imu stream')
%   Pass rec.sensors(k).t to checkImuTimestamps with the commanded rate.
%
% Lives in IMU Characterisation/Processing/.

verbose = true;
for i = 1:2:numel(varargin)
    if strcmpi(varargin{i}, 'Verbose'), verbose = varargin{i + 1}; end
end
txt = fileread(file);
txt = regexprep(txt, '\x1B\[[0-9;]*[A-Za-z]', '');           % ANSI escape sequences
lines = regexp(txt, '\r\n|\n|\r', 'split');
lines = lines(~cellfun(@isempty, strtrim(lines)));

names = {'LSM6DS3TR-C', 'LSM6DSV16X', 'ICM-42688-P'};
num = '(-?\d+\.\d+)';
pat = ['\[(\d+):(\d\d):(\d\d)\.(\d{3}),(\d{3})\]\s*<inf> imu: (\S+)\s+a\s+', num, '\s+', num, '\s+', num, ...
    '\s+g\s+', num, '\s+', num, '\s+', num, '\s*$'];
data = cell(1, numel(names)); for k = 1:numel(names), data{k} = zeros(0, 7); end
malformed = {}; notices = {}; cmds = {}; dropped = 0; droppedStream = 0; seenData = false;
for i = 1:numel(lines)
    L = lines{i};
    tok = regexp(L, pat, 'tokens', 'once');
    if ~isempty(tok)
        k = find(strcmp(names, tok{6}), 1);
        if isempty(k), malformed{end + 1} = L; continue; end %#ok<AGROW>
        v = str2double(tok([1:5, 7:12])); v = v(:).';
        t = v(1) * 3600 + v(2) * 60 + v(3) + v(4) / 1e3 + v(5) / 1e6;
        data{k}(end + 1, :) = [t, v(6:11)]; %#ok<AGROW>
        seenData = true;
        continue
    end
    dm = regexp(L, '---\s*(\d+) messages dropped', 'tokens', 'once');
    if ~isempty(dm)
        dropped = dropped + str2double(dm{1});
        if seenData, droppedStream = droppedStream + str2double(dm{1}); end   % after streaming began
    end
    if ~isempty(strfind(L, 'imu stream')), cmds{end + 1} = strtrim(L); end %#ok<AGROW>
    if ~isempty(regexp(L, '<inf> imu: \S+\s+a\s', 'once'))
        malformed{end + 1} = L; %#ok<AGROW>
    elseif ~isempty(regexp(L, '<(inf|wrn|err)>', 'once'))
        notices{end + 1} = strtrim(regexprep(L, '^(uart:~\$\s*)+', '')); %#ok<AGROW>
    end
end

sensors = struct('name', {}, 't', {}, 'acc', {}, 'gyr', {}, 'nRows', {});
for k = 1:numel(names)
    D = data{k};
    if isempty(D), continue; end
    [~, ord] = sort(D(:, 1)); D = D(ord, :);                     % arrival order is time order here
    sensors(end + 1) = struct('name', names{k}, 't', D(:, 1), 'acc', D(:, 2:4), 'gyr', D(:, 5:7), ...
        'nRows', size(D, 1)); %#ok<AGROW>
end
rec.sensors = sensors;
rec.meta = struct('file', file, 'nLines', numel(lines), 'nMalformed', numel(malformed), ...
    'malformed', {malformed(1:min(20, end))}, 'droppedMessages', dropped, 'droppedDuringStream', droppedStream, 'notices', {notices}, ...
    'streamCommands', {cmds});

if verbose
    fprintf('loadImuRecord: %s\n  %d lines, %d malformed IMU lines skipped, Zephyr messages dropped: %d in total, %d during streaming\n', ...
        file, numel(lines), numel(malformed), dropped, droppedStream);
    for k = 1:numel(sensors)
        s = sensors(k);
        fprintf('  %-12s %6d rows, t %.3f to %.3f s (%.1f s), median dt %.1f ms, |a| mean %.4f m/s^2\n', s.name, s.nRows, ...
            s.t(1), s.t(end), s.t(end) - s.t(1), 1e3 * median(diff(s.t)), mean(sqrt(sum(s.acc.^2, 2))));
    end
end
end