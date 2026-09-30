%% runStaticNoiseCharacterisation.m
% Static noise characterisation of one long, stationary 10 Hz shell record (all three IMUs).
% A thin driver around the validated functions; it adds no new estimators:
%
%   loadImuRecord -> longest reboot-free session -> checkImuTimestamps (per sensor)
%   -> alignImuSensors -> LONGEST CLEAN CONTIGUOUS BLOCK -> computeNoisePSD (+ coherence)
%   -> allanDeviation -> bias drift -> summary, .mat and figures
%
% CLEAN BLOCK RULE. After alignment the record is split wherever a tick interval exceeds
% GAP_FACTOR x the median interval or time does not increase; the longest piece (by duration)
% is analysed. Missing ticks are NEVER deleted and the two sides joined: that would hide
% elapsed time from the Allan deviation and the PSD. All blocks are listed in the summary.
%
% REBOOTS. A reboot restarts the uptime clock; each sensor's rows are split at backward time
% jumps and only the longest session is kept (loadImuRecord keeps arrival order for this).
%
% OUTPUT (IMU Characterisation/Results): staticNoise_<logname>.txt / .mat /
%   _asd.png / _coherence.png / _adev_acc.png / _adev_gyr.png
%
% Datasheet noise densities (one-sided, compare with the 'white' fields) are in DATASHEET_ACC /
% DATASHEET_GYR; NaN where not yet confirmed from the datasheet.
%
% Overrides (e.g. from the tester): set a struct STATIC_NOISE_CONFIG with any of the
% configuration names below before running the script.
%
% Lives in IMU Characterisation/Processing/.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
clearvars -except thisDir STATIC_NOISE_CONFIG
close all

%% ---- Configuration -----------------------------------------------------------------------------
cfg.LOG_FILE = '';                   % '' = newest *long*.log in DATA_DIR
cfg.DATA_DIR = fullfile(thisDir, '..', 'Data', 'Raw');
cfg.RESULTS_DIR = fullfile(thisDir, '..', 'Results');
cfg.COMMANDED_RATE = 10;             % Hz, for checkImuTimestamps
cfg.SEGMENT_LENGTH = 8192;           % Welch segment (samples); reduced automatically for short blocks
cfg.GAP_FACTOR = 1.5;                % block break: dt > GAP_FACTOR * median(dt)
cfg.MIN_BLOCK_S = 600;               % warn if the longest clean block is shorter
cfg.DRIFT_WINDOW_S = 3600;           % bias drift = mean(last window) - mean(first window)
cfg.SHOW = true;
% datasheet noise densities, one-sided; rows = sensors (LSM6DS3TR-C, LSM6DSV16X, ICM-42688-P), columns = x y z
cfg.DATASHEET_ACC = [90 90 90; 60 60 60; 65 65 70];          % ug/sqrt(Hz) (ST: high-performance mode)
cfg.DATASHEET_GYR = [NaN NaN NaN; NaN NaN NaN; 2.8 2.8 2.8]; % mdps/sqrt(Hz); NaN = not yet confirmed
if exist('STATIC_NOISE_CONFIG', 'var')
    fn = fieldnames(STATIC_NOISE_CONFIG);
    for i = 1:numel(fn), cfg.(fn{i}) = STATIC_NOISE_CONFIG.(fn{i}); end
end
G0 = 9.80665;

%% ---- 1  Load, session, timing -----------------------------------------------------------------------
if isempty(cfg.LOG_FILE)
    d = dir(fullfile(cfg.DATA_DIR, '*long*.log'));
    if isempty(d), error('runStaticNoiseCharacterisation:noLog', 'No *long*.log in %s; set cfg.LOG_FILE.', cfg.DATA_DIR); end
    [~, newest] = max([d.datenum]);
    cfg.LOG_FILE = fullfile(cfg.DATA_DIR, d(newest).name);
end
[~, logName] = fileparts(cfg.LOG_FILE);
if ~exist(cfg.RESULTS_DIR, 'dir'), mkdir(cfg.RESULTS_DIR); end
baseName = fullfile(cfg.RESULTS_DIR, ['staticNoise_', logName]);
Ltxt = {};
Ltxt{end + 1} = sprintf('runStaticNoiseCharacterisation  %s', datestr(now, 'yyyy-mm-dd HH:MM')); %#ok<TNOW1,DATST>
Ltxt{end + 1} = sprintf('log: %s', cfg.LOG_FILE);

tic;
rec = loadImuRecord(cfg.LOG_FILE, 'Verbose', false);
Ltxt{end + 1} = sprintf('parsed %d lines in %.0f s: %d malformed IMU lines, Zephyr drops %d during streaming (%d total), time resets %s', ...
    rec.meta.nLines, toc, rec.meta.nMalformed, rec.meta.droppedDuringStream, rec.meta.droppedMessages, ...
    mat2str(rec.meta.timeResets));
[rec, sess] = longestSession(rec);
Ltxt{end + 1} = sprintf('sessions (reboots split the uptime clock): %d; kept session %d, %.2f h', sess.n, sess.kept, sess.spanH);

Ltxt{end + 1} = '';
Ltxt{end + 1} = 'TIMING (logging-path timestamps; gap-like intervals are inferred, not proven losses)';
Ltxt{end + 1} = sprintf('  %-12s %8s %9s %9s %8s %8s %6s %9s', 'sensor', 'rows', 'f_recv', 'f_tick', 'resid', 'max dt', 'gapLk', 'drift');
Ltxt{end + 1} = sprintf('  %-12s %8s %9s %9s %8s %8s %6s %9s', '', '', '(Hz)', '(Hz)', '(ms)', '(ms)', '', '(ppm)');
timing = struct([]);
for k = 1:numel(rec.sensors)
    s = rec.sensors(k);
    o = checkImuTimestamps(s.t, 'NominalRate', cfg.COMMANDED_RATE, 'Verbose', false);
    timing = [timing, o]; %#ok<AGROW>
    Ltxt{end + 1} = sprintf('  %-12s %8d %9.4f %9.4f %8.2f %8.0f %6d %+9.1f', s.name, s.nRows, (s.nRows - 1) / (s.t(end) - s.t(1)), ...
        o.rate, 1e3 * o.jitterStd, 1e3 * max(diff(s.t)), o.nGaps, o.driftPpm); %#ok<SAGROW>
end

%% ---- 2  Align and select the longest clean block --------------------------------------------------------
A = alignImuSensors(rec);
[blk, blocks] = longestBlock(A.t, cfg.GAP_FACTOR);
nB = numel(blk); tB = A.t(blk);
fs = (nB - 1) / (tB(end) - tB(1));
Ltxt{end + 1} = '';
Ltxt{end + 1} = sprintf('ALIGNMENT: %d common ticks (%d dropped where a sensor was missing)', numel(A.t), A.nDropped);
Ltxt{end + 1} = sprintf('CLEAN BLOCKS (break where dt > %.1f x median or time does not increase): %d', cfg.GAP_FACTOR, size(blocks, 1));
for b = 1:size(blocks, 1)
    Ltxt{end + 1} = sprintf('  block %d: ticks %d-%d, t %.1f-%.1f s, %.2f h%s', b, blocks(b, 1), blocks(b, 2), ...
        A.t(blocks(b, 1)), A.t(blocks(b, 2)), (A.t(blocks(b, 2)) - A.t(blocks(b, 1))) / 3600, ...
        ternary(blocks(b, 1) == blk(1), '   <- analysed', '')); %#ok<SAGROW>
end
Ltxt{end + 1} = sprintf('analysed block: %d ticks, %.2f h, fs = %.4f Hz (received rate in the block)', nB, (tB(end) - tB(1)) / 3600, fs);
if tB(end) - tB(1) < cfg.MIN_BLOCK_S
    Ltxt{end + 1} = sprintf('WARNING: longest clean block is shorter than %.0f s', cfg.MIN_BLOCK_S);
end
L = min(cfg.SEGMENT_LENGTH, 2^floor(log2(max(8, nB / 4))));

%% ---- 3  PSD, coherence, Allan deviation, drift ---------------------------------------------------------------
kinds = {'acc', 'gyr'}; unitsK = {'m/s^2', 'rad/s'}; axN = 'xyz';
psd = struct(); adev = struct(); res = struct();
for q = 1:2
    for ax = 1:3
        X = squeeze(A.(kinds{q})(blk, ax, :));
        p = computeNoisePSD(X, fs, 'Time', tB, 'SegmentLength', L, 'ChannelNames', A.names, 'Units', unitsK{q});
        a = allanDeviation(X, fs, 'Time', tB, 'ChannelNames', A.names, 'Units', unitsK{q});
        psd.(kinds{q}){ax} = p; adev.(kinds{q}){ax} = a;
        wide = p.f >= 0.5 & p.f <= 0.45 * fs;
        nw = round(cfg.DRIFT_WINDOW_S * fs);
        nw = min(nw, floor(nB / 3));
        drift = mean(X(end - nw + 1:end, :), 1) - mean(X(1:nw, :), 1);
        res.(kinds{q})(ax) = struct('whitePSD', sqrt(median(p.Pxx(wide, :), 1)), 'whiteADEV', [a.terms.white], ...
            'basinASD', sqrt(p.bands.basin.meanPSD), 'mizASD', sqrt(p.bands.miz.meanPSD), 'mizBins', p.bands.miz.nBins, ...
            'basinSnnPerRadS', p.bands.basin.meanPSD / (2 * pi), 'cohBasin', p.bands.basin.meanCoherence, ...
            'cohMiz', p.bands.miz.meanCoherence, 'cohAbove95', mean(p.coherence(p.f > 0, :) > p.coherence95Threshold, 1), ...
            'cohBasinAbove95', mean(p.coherence(inBand(p, 'basin'), :) > p.coherence95Threshold, 1), ...
            'cohMizAbove95', mean(p.coherence(inBand(p, 'miz'), :) > p.coherence95Threshold, 1), ...
            'coh95', p.coherence95Threshold, 'plateau', [a.terms.plateau], 'flicker', [a.terms.flicker], ...
            'tauAtMin', [a.terms.tauAtMin], 'randomWalk', [a.terms.randomWalk], 'drift', drift, 'driftWindowS', nw / fs);
    end
end

%% ---- 4  Summary ------------------------------------------------------------------------------------------------
sc = {1e6 / G0, 1e3 * 180 / pi};                               % to ug/sqrt(Hz) and mdps/sqrt(Hz)
scN = {'ug/rtHz', 'mdps/rtHz'};
dsh = {cfg.DATASHEET_ACC, cfg.DATASHEET_GYR};
p1 = psd.acc{1};
Ltxt{end + 1} = '';
Ltxt{end + 1} = sprintf('Welch: L = %d (%.0f s), %d segments, df = %.4g Hz, dof %.0f, 95%% CI x[%.2f %.2f]; coherence 95%% threshold %.3f', ...
    L, L / fs, p1.nSegments, p1.df, p1.dofEffective, p1.ci95, p1.coherence95Threshold);
Ltxt{end + 1} = sprintf('Allan: trusted tau up to %.0f s (edf >= 4)', max(adev.acc{1}.tau(adev.acc{1}.edf >= 4)));
for q = 1:2
    Ltxt{end + 1} = '';
    Ltxt{end + 1} = sprintf('%s  (white: one-sided density in %s; PSD = median over 0.5 Hz..0.45 fs, ADEV = slope -1/2 fit)', ...
        upper(ternary(q == 1, 'accelerometer', 'gyroscope')), scN{q});
    Ltxt{end + 1} = sprintf('  %-3s %-12s %9s %9s %8s %9s %9s %11s %6s %6s %6s %6s %10s %9s %10s %11s', 'ax', 'sensor', 'white PSD', ...
        'white AD', 'x sheet', 'basin', 'MIZ', 'S_nn basin', 'cohB', '>95B', 'cohM', '>95M', 'plateau', 'tauMin', 'K (RW)', 'drift');
    Ltxt{end + 1} = sprintf('  %-3s %-12s %9s %9s %8s %9s %9s %11s %6s %6s %6s %6s %10s %9s %10s %11s', '', '', '', '', '', '(same)', ...
        '(same)', '(per rad/s)', 'mean', '(%)', 'mean', '(%)', ['(', unitsK{q}, ')'], '(s)', '', ['(', unitsK{q}, ')']);
    for ax = 1:3
        r = res.(kinds{q})(ax);
        for k = 1:numel(A.names)
            pk = pairsOf(k, p1.pairs);
            Ltxt{end + 1} = sprintf('  %-3s %-12s %9.1f %9.1f %8.2f %9.1f %9.1f %11.3g %6.2f %6.1f %6.2f %6.1f %10.3g %9.0f %10.3g %+11.2e', ...
                axN(ax), A.names{k}, sc{q} * r.whitePSD(k), sc{q} * r.whiteADEV(k), sc{q} * r.whitePSD(k) / dsh{q}(k, ax), ...
                sc{q} * r.basinASD(k), sc{q} * r.mizASD(k), r.basinSnnPerRadS(k), mean(r.cohBasin(pk)), ...
                100 * mean(r.cohBasinAbove95(pk)), mean(r.cohMiz(pk)), 100 * mean(r.cohMizAbove95(pk)), r.plateau(k), ...
                r.tauAtMin(k), r.randomWalk(k), r.drift(k)); %#ok<SAGROW>
        end
    end
end
Ltxt{end + 1} = '';
Ltxt{end + 1} = 'NOTES';
Ltxt{end + 1} = ' - x sheet = PSD white level / datasheet density; the IMUs'' internal ODR and filter are unknown (aliasing hypothesis).';
Ltxt{end + 1} = [' - cohB / cohM = mean coherence with the other two sensors in the basin / MIZ band (descriptive only);'];
Ltxt{end + 1} = ['   >95B / >95M = % of band bins above the single-bin 95% threshold: about 5% is expected for'];
Ltxt{end + 1} = ['   independent sensors, clearly more means common motion in that band.'];
Ltxt{end + 1} = [' - gyro: the 0.001 rad/s print step is comparable to the gyro noise, so gyro levels are print-limited in'];
Ltxt{end + 1} = ['   EITHER direction (rounding adds ~1.3e-4 rad/s/rtHz when the noise spans several steps, but hides noise'];
Ltxt{end + 1} = ['   when it spans less than one: a synthetic 5.7 mdps/rtHz gyro reads ~4 mdps/rtHz). Not quotable yet.'];
Ltxt{end + 1} = ' - plateau is NaN where no flat Allan region was found; K is NaN where no +1/2 region was found in the trusted taus.';
Ltxt{end + 1} = sprintf(' - drift = mean of the last %.0f s minus the first %.0f s of the block (thermal and bias drift together).', ...
    res.acc(1).driftWindowS, res.acc(1).driftWindowS);
summary = strjoin(Ltxt, sprintf('\n'));
fprintf('%s\n', summary);
fid = fopen([baseName, '.txt'], 'w'); fprintf(fid, '%s\n', summary); fclose(fid);
out = struct('cfg', cfg, 'meta', rec.meta, 'session', sess, 'timing', timing, 'blocks', blocks, 'block', blk([1 end]), ...
    'fs', fs, 'names', {A.names}, 'psd', psd, 'adev', adev, 'results', res, 'summary', summary);
save([baseName, '.mat'], 'out', '-v7');

%% ---- 5  Figures --------------------------------------------------------------------------------------------------
% Each figure is drawn into explicit axes (never via figure(h) / gcf, which makes hidden figures visible
% and lets a click or a closed window redirect or delete them), saved at once, then closed if hidden.
% A failed export is a warning: the .txt and .mat above are the results, figures must never abort the run.
vis = ternary(cfg.SHOW, 'on', 'off');
nSaved = 0;
fig = newFig(vis);                                            % 1  accelerometer ASD
for ax = 1:3
    p = psd.acc{ax}; k = p.f > 0; h = panel(fig, ax);
    shade(h, p, [min(p.ASD(k, :), [], 1), max(p.ASD(k, :), [], 1)] .* [0.5 0.5 0.5 2 2 2]);
    loglog(h, p.f(k), p.ASD(k, :) * sc{1}, 'LineWidth', 1.2);
    set(h, 'XScale', 'log', 'YScale', 'log'); grid(h, 'on'); xlim(h, p.f([2 end]));
    xlabel(h, 'f (Hz)'); ylabel(h, 'ASD (\mug/\surdHz)'); title(h, sprintf('a_%s', axN(ax)));
    if ax == 1, legend(h, A.names, 'Interpreter', 'none', 'Location', 'southwest'); end
end
nSaved = nSaved + saveFig(fig, [baseName, '_asd.png'], cfg.SHOW);
fig = newFig(vis);                                            % 2  accelerometer coherence
for ax = 1:3
    p = psd.acc{ax}; k = p.f > 0; h = panel(fig, ax);
    shade(h, p, [0 1]);
    plot(h, p.f(k), p.coherence(k, :), 'LineWidth', 1.0);
    plot(h, p.f([2 end]), p.coherenceNoiseFloor * [1 1], 'k:', p.f([2 end]), p.coherence95Threshold * [1 1], 'k--');
    set(h, 'XScale', 'log'); ylim(h, [0 1]); grid(h, 'on'); xlim(h, p.f([2 end]));
    xlabel(h, 'f (Hz)'); ylabel(h, '\gamma^2'); title(h, sprintf('a_%s coherence', axN(ax)));
    if ax == 1
        lab = arrayfun(@(i) [A.names{p.pairs(i, 1)}, '/', A.names{p.pairs(i, 2)}], 1:size(p.pairs, 1), 'UniformOutput', false);
        legend(h, [lab, {'expected', '95%'}], 'Interpreter', 'none', 'Location', 'northeast');
    end
end
nSaved = nSaved + saveFig(fig, [baseName, '_coherence.png'], cfg.SHOW);
yLab = {'ADEV (\mug)', 'ADEV (mdps)'}; pre = 'ag'; suf = {'_adev_acc.png', '_adev_gyr.png'};
for q = 1:2                                                   % 3, 4  Allan deviation
    fig = newFig(vis);
    for ax = 1:3
        a = adev.(kinds{q}){ax}; h = panel(fig, ax);
        loglog(h, a.tau, a.adev * sc{q}, 'LineWidth', 1.3);
        set(h, 'XScale', 'log', 'YScale', 'log'); grid(h, 'on');
        xlabel(h, '\tau (s)'); ylabel(h, yLab{q}); title(h, sprintf('%s_%s', pre(q), axN(ax)));
        if ax == 1, legend(h, A.names, 'Interpreter', 'none', 'Location', 'southwest'); end
    end
    nSaved = nSaved + saveFig(fig, [baseName, suf{q}], cfg.SHOW);
end
fprintf('Saved: %s.txt / .mat and %d of 4 figures (_asd, _coherence, _adev_acc, _adev_gyr .png)\n', baseName, nSaved);

%% ================================================================================================
function [rec, sess] = longestSession(rec)
% Split each sensor at backward time jumps (> 1 s: a reboot) and keep the session, common to
% all sensors, with the longest span.
S = rec.sensors; nS = numel(S);
sid = cell(1, nS); nSess = zeros(1, nS);
for k = 1:nS
    sid{k} = 1 + cumsum([0; diff(S(k).t) < -1]);
    nSess(k) = sid{k}(end);
end
n = min(nSess);
span = zeros(1, n);
for s = 1:n                                                   % span common to ALL sensors
    tStart = zeros(1, nS); tStop = zeros(1, nS);
    for k = 1:nS
        tk = S(k).t(sid{k} == s);
        tStart(k) = tk(1); tStop(k) = tk(end);
    end
    span(s) = min(tStop) - max(tStart);
end
[~, kept] = max(span);
for k = 1:nS
    sel = sid{k} == kept;
    S(k).t = S(k).t(sel); S(k).acc = S(k).acc(sel, :); S(k).gyr = S(k).gyr(sel, :); S(k).nRows = nnz(sel);
end
rec.sensors = S;
sess = struct('n', n, 'kept', kept, 'spanH', span(kept) / 3600, 'spansH', span / 3600);
end

function [idx, blocks] = longestBlock(t, factor)
% Contiguous blocks between breaks (dt > factor * median(dt), or dt <= 0); longest by duration.
dt = diff(t);
brk = find(dt > factor * median(dt) | dt <= 0);
starts = [1; brk + 1]; ends = [brk; numel(t)];
blocks = [starts, ends];
dur = t(ends) - t(starts);
[~, b] = max(dur);
idx = (starts(b):ends(b)).';
end

function in = inBand(p, name)
r = p.bands.(name).range;
in = p.f >= r(1) & p.f <= r(2) & p.f > 0;
end

function k = pairsOf(c, pairs)
k = find(any(pairs == c, 2));
end

function shade(h, p, yl)
if numel(yl) > 2, yl = [min(yl(1:3)), max(yl(4:6))] * 1e6 / 9.80665; end
cols = [1.0 0.92 0.85; 0.85 0.92 1.0];
bn = {'miz', 'basin'};
for b = 1:2
    r = p.bands.(bn{b}).range;
    patch(h, [r(1) r(2) r(2) r(1)], [yl(1) yl(1) yl(2) yl(2)], cols(b, :), 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function fig = newFig(vis)
fig = figure('Color', 'w', 'Visible', vis, 'Position', [40 40 1400 420]);
end

function h = panel(fig, ax)
% axes ax of a 1 x 3 row, parented explicitly (no dependence on the current figure)
h = axes('Parent', fig, 'Position', [0.05 + (ax - 1) * 0.325, 0.14, 0.27, 0.76]);
hold(h, 'on');
end

function ok = saveFig(fig, file, keepOpen)
% Export PNG; fall back to print, then saveas; warn (never error) if all fail. Close hidden figures.
ok = false;
if ~ishghandle(fig)
    warning('runStaticNoiseCharacterisation:figure', 'Figure for %s was closed before saving; skipped.', file);
    return
end
if exist(file, 'file') == 2, delete(file); end               % no stale PNG can pass for a new one
drawnow;
attempts = {@() exportgraphics(fig, file, 'Resolution', 150), @() print(fig, file, '-dpng', '-r150'), ...
    @() saveas(fig, file)};
msg = '';
for i = 1:numel(attempts)
    if ~ishghandle(fig), break; end
    try
        attempts{i}(); ok = exist(file, 'file') == 2;
        if ok, break; end
    catch err
        msg = err.message;
    end
end
if ~ok
    warning('runStaticNoiseCharacterisation:figure', 'Could not save %s (%s). Results .txt/.mat are unaffected.', file, msg);
end
if ~keepOpen && ishghandle(fig), close(fig); end
end