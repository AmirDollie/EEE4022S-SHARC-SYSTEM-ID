function rec = loadImuRecord(file, varargin)
%LOADIMURECORD  Parse a SHARC board (Zephyr shell) IMU stream log into per-sensor arrays.
%   rec = loadImuRecord(file)
%   rec = loadImuRecord(file, 'Verbose', false)
%   rec = loadImuRecord(file, 'BlockBytes', 16e6)
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
%   MEMORY: the file is read in blocks of whole lines ('BlockBytes', default 16e6), so a long
%   overnight log (hundreds of MB) is parsed with bounded memory; results do not depend on the
%   block size.
%
%   OUTPUT rec.sensors(k): name, t (s, device uptime), acc (N x 3), gyr (N x 3), nRows
%          rec.meta: file, nLines, nMalformed, malformed (first 20 lines), droppedMessages (all),
%                    droppedDuringStream (notices after the first IMU data line: the ones that matter),
%                    timeResets (backward time jumps per sensor), notices (non-data log lines:
%                    boot, warnings, errors), streamCommands (lines containing 'imu stream')
%   Pass rec.sensors(k).t to checkImuTimestamps with the commanded rate.
%
% Lives in IMU Characterisation/Processing/.

verbose = true; blockBytes = 16e6;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'verbose', verbose = varargin{i + 1};
        case 'blockbytes', blockBytes = varargin{i + 1};
        otherwise, error('loadImuRecord:args', 'Unknown option ''%s''.', varargin{i});
    end
end

names = {'LSM6DS3TR-C', 'LSM6DSV16X', 'ICM-42688-P'};
num = '(-?\d+\.\d+)';
pat = ['\[(\d+):(\d\d):(\d\d)\.(\d{3}),(\d{3})\]\s*<inf> imu: (\S+)\s+a\s+', num, '\s+', num, '\s+', num, ...
    '\s+g\s+', num, '\s+', num, '\s+', num, '\s*$'];

% ---- read in blocks of whole lines (bounded memory for long records) ------------------------------
fid = fopen(file, 'r');
if fid < 0, error('loadImuRecord:open', 'Cannot open %s.', file); end
cleaner = onCleanup(@() fclose(fid));
Vc = {}; Nc = {};                               % per-block numeric rows and sensor index
malformed = {}; notices = {}; cmds = {};
dropped = 0; droppedStream = 0; seenData = false; nLines = 0; nMalformed = 0;
carry = '';
while true
    chunk = fread(fid, blockBytes, '*char').';
    atEnd = numel(chunk) < blockBytes;
    txt = [carry, chunk];
    if ~atEnd
        cut = find(txt == char(10), 1, 'last');
        if isempty(cut)                          % no newline in this block: keep reading
            carry = txt; continue
        end
        carry = txt(cut + 1:end); txt = txt(1:cut);
    else
        carry = '';
    end
    [Vb, Nb, info] = parseBlock(txt, pat, names, seenData);
    Vc{end + 1} = Vb; Nc{end + 1} = Nb; %#ok<AGROW>
    nLines = nLines + info.nLines;
    nMalformed = nMalformed + numel(info.malformed);
    if numel(malformed) < 20, malformed = [malformed; info.malformed(1:min(end, 20 - numel(malformed)))]; end %#ok<AGROW>
    notices = [notices; info.notices]; cmds = [cmds; info.cmds]; %#ok<AGROW>
    dropped = dropped + info.dropped; droppedStream = droppedStream + info.droppedStream;
    seenData = seenData || info.anyData;
    if atEnd, break; end
end
V = vertcat(Vc{:}); sid = vertcat(Nc{:});
if isempty(V), V = zeros(0, 11); sid = zeros(0, 1); end
tt = V(:, 1) * 3600 + V(:, 2) * 60 + V(:, 3) + V(:, 4) / 1e3 + V(:, 5) / 1e6;

% ---- per sensor, arrival order -------------------------------------------------------------------------
sensors = struct('name', {}, 't', {}, 'acc', {}, 'gyr', {}, 'nRows', {});
resets = zeros(1, 0);
for k = 1:numel(names)
    sel = sid == k;
    if ~any(sel), continue; end
    sensors(end + 1) = struct('name', names{k}, 't', tt(sel), 'acc', V(sel, 6:8), 'gyr', V(sel, 9:11), ...
        'nRows', nnz(sel)); %#ok<AGROW>
    resets(end + 1) = nnz(diff(tt(sel)) < -1); %#ok<AGROW>
end
rec.sensors = sensors;
rec.meta = struct('file', file, 'nLines', nLines, 'nMalformed', nMalformed, ...
    'malformed', {malformed.'}, 'droppedMessages', dropped, 'droppedDuringStream', droppedStream, ...
    'timeResets', resets, 'notices', {notices.'}, 'streamCommands', {cmds.'});

if verbose
    fprintf('loadImuRecord: %s\n  %d lines, %d malformed IMU lines skipped, Zephyr messages dropped: %d in total, %d during streaming\n', ...
        file, nLines, nMalformed, dropped, droppedStream);
    for k = 1:numel(sensors)
        s = sensors(k);
        fprintf('  %-12s %6d rows, t %.3f to %.3f s, median dt %.1f ms, |a| mean %.4f m/s^2, time resets %d\n', s.name, ...
            s.nRows, s.t(1), s.t(end), 1e3 * median(diff(s.t)), mean(sqrt(sum(s.acc.^2, 2))), resets(k));
    end
end
end

function [V, sid, info] = parseBlock(txt, pat, names, seenBefore)
% Parse one block of whole lines. Numeric rows V (n x 11: h m s ms us ax ay az gx gy gz) and
% sensor index sid, in arrival order; counts and notices for everything else.
txt = regexprep(txt, '\x1B\[[0-9;]*[A-Za-z]', '');            % ANSI escape sequences
lines = regexp(txt, '\r\n|\n|\r', 'split');
lines = lines(~cellfun(@isempty, strtrim(lines)));
lines = lines(:);
tok = regexp(lines, pat, 'tokens', 'once');
isData = ~cellfun(@isempty, tok);
T = cellfun(@(c) c(:).', tok(isData), 'UniformOutput', false);  % row cells (orientation differs by platform)
T = vertcat(T{:});                                            % nData x 12 cellstr
if isempty(T), T = cell(0, 12); end
[known, sid] = ismember(T(:, 6), names);
dataIdx = find(isData);
badName = dataIdx(~known);
T = T(known, :); sid = sid(known); dataIdx = dataIdx(known);
nD = size(T, 1);
if nD > 0
    numStr = T(:, [1:5, 7:12]).';
    V = reshape(sscanf(sprintf('%s ', numStr{:}), '%f'), 11, nD).';
else
    V = zeros(0, 11);
end
firstData = min([dataIdx; numel(lines) + 1]);
if seenBefore, firstData = 0; end                             % streaming began in an earlier block
other = find(~isData);
L = lines(other);
dm = regexp(L, '---\s*(\d+) messages dropped', 'tokens', 'once');
hasDrop = ~cellfun(@isempty, dm);
dropN = zeros(numel(L), 1);
if any(hasDrop), dropN(hasDrop) = cellfun(@(c) str2double(c{1}), dm(hasDrop)); end
isImuLike = ~cellfun(@isempty, regexp(L, '<inf> imu: \S+\s+a\s', 'once'));
isNotice = ~isImuLike & ~cellfun(@isempty, regexp(L, '<(inf|wrn|err)>', 'once'));
info = struct('nLines', numel(lines), 'dropped', sum(dropN), 'droppedStream', sum(dropN(other > firstData)), ...
    'cmds', {strtrim(L(~cellfun(@isempty, strfind(L, 'imu stream'))))}, ...
    'malformed', {[lines(badName); L(isImuLike)]}, ...
    'notices', {strtrim(regexprep(L(isNotice), '^(uart:~\$\s*)+', ''))}, 'anyData', nD > 0);
sid = sid(:);
end