function rec = loadImuRecord(file, varargin)
%LOADIMURECORD  Parse a SHARC board (Zephyr shell) IMU stream log into per-sensor arrays.
%   rec = loadImuRecord(file)
%   rec = loadImuRecord(file, 'Verbose', false)
%
%   INPUT FORMAT (PuTTY "All session output" log of `imu stream on <rate>`): Zephyr log
%   lines, ANSI colour codes and shell prompts interleaved, e.g.
%     uart:~$ [00:05:24.390,000] <inf> imu: ICM-42688-P  a   1.086  -0.086   9.849  g  -0.000  -0.005   0.000
%   timestamp [hh:mm:ss.mmm,uuu] = MCU uptime when the log message was created (NOT the IMU
%   sample time; 1 ms resolution in the first logs); a = accelerometer (m/s^2); g = gyroscope
%   (rad/s, Zephyr SI convention: confirm with the firmware).
%
%   ROBUSTNESS: shell echo and log output can interleave inside a line (seen in the first log:
%   '1.0801.7438', ',001.74300', 'L<ESC>[DSM6DSV16X'). Only lines matching the strict pattern
%   (3-digit ms and us fields, known sensor name, signed decimal numbers) are kept; everything
%   else is counted in rec.meta. Zephyr's '--- N messages dropped ---' notices are counted too.
%   There is no sample counter, so a missing line cannot be told from a slow loop iteration.
%
%   ORDER AND REBOOTS: rows are kept in ARRIVAL order (not sorted by time). A board reboot
%   restarts the uptime clock, so a record that spans one contains a backward time jump; these
%   are counted per sensor in rec.meta.timeResets. Sorting would interleave the sessions, so
%   callers that need one monotonic session must split at the resets (runStaticNoiseCharacterisation
%   keeps the longest session).
%
%   OUTPUT rec.sensors(k): name, t (s, device uptime), acc (N x 3), gyr (N x 3), nRows
%          rec.meta: file, nLines, nMalformed, malformed (first 20 lines), droppedMessages (all),
%                    droppedDuringStream (notices after the first IMU data line: the ones that matter),
%                    timeResets (backward time jumps per sensor), notices (non-data log lines:
%                    boot, warnings, errors), streamCommands (lines containing 'imu stream')
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
lines = lines(:);

names = {'LSM6DS3TR-C', 'LSM6DSV16X', 'ICM-42688-P'};
num = '(-?\d+\.\d+)';
pat = ['\[(\d+):(\d\d):(\d\d)\.(\d{3}),(\d{3})\]\s*<inf> imu: (\S+)\s+a\s+', num, '\s+', num, '\s+', num, ...
    '\s+g\s+', num, '\s+', num, '\s+', num, '\s*$'];

% ---- data lines, vectorised -------------------------------------------------------------------
tok = regexp(lines, pat, 'tokens', 'once');
isData = ~cellfun(@isempty, tok);
T = cellfun(@(c) c(:).', tok(isData), 'UniformOutput', false);  % row cells (orientation differs by platform)
T = vertcat(T{:});                                            % nData x 12 cellstr
if isempty(T), T = cell(0, 12); end
known = ismember(T(:, 6), names);
dataIdx = find(isData);
badName = dataIdx(~known);
T = T(known, :); dataIdx = dataIdx(known);
nD = size(T, 1);
if nD > 0
    numStr = T(:, [1:5, 7:12]).';
    V = sscanf(sprintf('%s ', numStr{:}), '%f');
    V = reshape(V, 11, nD).';
    tt = V(:, 1) * 3600 + V(:, 2) * 60 + V(:, 3) + V(:, 4) / 1e3 + V(:, 5) / 1e6;
else
    V = zeros(0, 11); tt = zeros(0, 1);
end
firstData = min([dataIdx; numel(lines) + 1]);

% ---- everything else -------------------------------------------------------------------------------
other = find(~isData);
L = lines(other);
dm = regexp(L, '---\s*(\d+) messages dropped', 'tokens', 'once');
hasDrop = ~cellfun(@isempty, dm);
dropN = zeros(numel(L), 1);
if any(hasDrop), dropN(hasDrop) = cellfun(@(c) str2double(c{1}), dm(hasDrop)); end
dropped = sum(dropN);
droppedStream = sum(dropN(other > firstData));
cmds = strtrim(L(~cellfun(@isempty, strfind(L, 'imu stream'))));
isImuLike = ~cellfun(@isempty, regexp(L, '<inf> imu: \S+\s+a\s', 'once'));
malformed = [lines(badName); L(isImuLike)];
isNotice = ~isImuLike & ~cellfun(@isempty, regexp(L, '<(inf|wrn|err)>', 'once'));
notices = strtrim(regexprep(L(isNotice), '^(uart:~\$\s*)+', ''));

% ---- per sensor, arrival order -------------------------------------------------------------------------
sensors = struct('name', {}, 't', {}, 'acc', {}, 'gyr', {}, 'nRows', {});
resets = zeros(1, 0);
for k = 1:numel(names)
    sel = strcmp(T(:, 6), names{k});
    if ~any(sel), continue; end
    sensors(end + 1) = struct('name', names{k}, 't', tt(sel), 'acc', V(sel, 6:8), 'gyr', V(sel, 9:11), ...
        'nRows', nnz(sel)); %#ok<AGROW>
    resets(end + 1) = nnz(diff(tt(sel)) < -1); %#ok<AGROW>
end
rec.sensors = sensors;
rec.meta = struct('file', file, 'nLines', numel(lines), 'nMalformed', numel(malformed), ...
    'malformed', {malformed(1:min(20, end)).'}, 'droppedMessages', dropped, 'droppedDuringStream', droppedStream, ...
    'timeResets', resets, 'notices', {notices.'}, 'streamCommands', {cmds.'});

if verbose
    fprintf('loadImuRecord: %s\n  %d lines, %d malformed IMU lines skipped, Zephyr messages dropped: %d in total, %d during streaming\n', ...
        file, numel(lines), numel(malformed), dropped, droppedStream);
    for k = 1:numel(sensors)
        s = sensors(k);
        fprintf('  %-12s %6d rows, t %.3f to %.3f s, median dt %.1f ms, |a| mean %.4f m/s^2, time resets %d\n', s.name, ...
            s.nRows, s.t(1), s.t(end), 1e3 * median(diff(s.t)), mean(sqrt(sum(s.acc.^2, 2))), resets(k));
    end
end
end