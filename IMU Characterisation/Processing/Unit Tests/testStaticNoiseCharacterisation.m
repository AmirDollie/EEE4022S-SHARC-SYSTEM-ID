%% testStaticNoiseCharacterisation.m
% End-to-end test of runStaticNoiseCharacterisation.m on a synthetic Zephyr-format log:
% 2.5 h at ~9.75 Hz, three IMUs with known white noise, a 60 s outage at 1.5 h, then a reboot
% (uptime restarts) followed by 10 min of data. The driver must keep the first session, split
% it at the outage, analyse the 1.5 h side (not the 1 h side, not the post-reboot session), and
% recover the known white levels through both the PSD and the Allan routes.
%
% Lives in IMU Characterisation/Processing/Unit-Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
nPass = 0; nFail = 0;
work = fullfile(tempdir, 'staticNoiseTest');
if ~exist(fullfile(work, 'Raw'), 'dir'), mkdir(fullfile(work, 'Raw')); end
logFile = fullfile(work, 'Raw', 'imu_synthetic_10Hz_long.log');

% ---- write the synthetic log ------------------------------------------------------------------------
dtTick = 0.1025; fs = 1 / dtTick;
Nacc = [8.5e-3 2.3e-3 3.0e-3];                  % one-sided white, m/s^2/sqrt(Hz), like the real sensors
bias = [0.10 0.02 10.18; 0.38 0.12 10.52; 0.35 -0.06 9.92];
names = {'LSM6DS3TR-C  ', 'LSM6DSV16X   ', 'ICM-42688-P  '};
rng(51);
fid = fopen(logFile, 'w');
fprintf(fid, '=~=~=~=~= PuTTY log synthetic =~=~=~=~=\r\n');
fprintf(fid, 'uart:~$ *** Booting Zephyr OS build synthetic ***\r\n');
fprintf(fid, 'uart:~$ imu stream on 10\r\n');
sessions = {[0, 1.5 * 3600], [1.5 * 3600 + 60, 2.5 * 3600], [0, 600]};   % [start end] uptime (s); 3rd after reboot
t0 = 196.1;
for sidx = 1:3
    if sidx == 3
        fprintf(fid, 'uart:~$ *** Booting Zephyr OS build synthetic ***\r\n');
        fprintf(fid, 'uart:~$ imu stream on 10\r\n');
        t0 = 5.0;
    end
    ticks = (ceil(sessions{sidx}(1) / dtTick):floor(sessions{sidx}(2) / dtTick)).';
    tt = t0 + ticks * dtTick;
    M = zeros(numel(tt), 3 * 10);
    for k = 1:3
        ts = floor((tt + (k - 1) * 1e-3) * 1e3) / 1e3;                 % 1 ms print resolution
        a = bias(k, :) + Nacc(k) * sqrt(fs / 2) * randn(numel(tt), 3);
        g = [0.02 -0.06 0.01] + 1e-4 * sqrt(fs / 2) * randn(numel(tt), 3);
        hh = floor(ts / 3600); mm = floor(mod(ts, 3600) / 60); ss = mod(ts, 60);
        M(:, (k - 1) * 10 + (1:10)) = [hh, mm, floor(ss), round((ss - floor(ss)) * 1e3), a, g];
    end
    fmt = '';
    for k = 1:3
        fmt = [fmt, char(27), '[1;32muart:~$ ', char(27), '[m[%02d:%02d:%02d.%03d,000] ', char(27), '[0m<inf> imu: ', ...
            names{k}, 'a %7.3f %7.3f %7.3f  g %7.3f %7.3f %7.3f', char(27), '[0m\r\n']; %#ok<AGROW>
    end
    fprintf(fid, fmt, M.');
end
fclose(fid);

% ---- run the driver ----------------------------------------------------------------------------------
runDriver(fullfile(thisDir, '..', 'runStaticNoiseCharacterisation.m'), struct('LOG_FILE', logFile, ...
    'RESULTS_DIR', fullfile(work, 'Results'), 'SHOW', false, 'MIN_BLOCK_S', 600));
S = load(fullfile(work, 'Results', 'staticNoise_imu_synthetic_10Hz_long.mat')); out = S.out;

% ---- checks ----------------------------------------------------------------------------------------------
[nPass, nFail] = check(isequal(out.meta.timeResets, [1 1 1]), 'one reboot detected in every sensor', nPass, nFail);
[nPass, nFail] = check(out.session.n == 2 && out.session.kept == 1 && abs(out.session.spanH - 2.5) < 0.01, ...
    sprintf('kept the pre-reboot session (%.2f h of %d sessions)', out.session.spanH, out.session.n), nPass, nFail);
[nPass, nFail] = check(size(out.blocks, 1) == 2, sprintf('outage splits the session into 2 blocks (found %d)', size(out.blocks, 1)), ...
    nPass, nFail);
tA = out.psd.acc{1}.N / out.fs / 3600;
[nPass, nFail] = check(abs(tA - 1.5) < 0.01, sprintf('analysed the longer 1.5 h side (%.3f h), not joined across the gap', tA), ...
    nPass, nFail);
[nPass, nFail] = check(abs(out.fs / fs - 1) < 1e-3, sprintf('fs from the block %.4f Hz (truth %.4f)', out.fs, fs), nPass, nFail);
r = out.results.acc;
wP = reshape([r.whitePSD], 3, 3); wA = reshape([r.whiteADEV], 3, 3);           % sensor x axis
eP = max(max(abs(wP ./ Nacc(:) - 1))); eA = max(max(abs(wA ./ Nacc(:) - 1)));
[nPass, nFail] = check(eP < 0.05, sprintf('white level via PSD within 5%% on all 9 channels (max err %.1f%%)', 100 * eP), nPass, nFail);
[nPass, nFail] = check(eA < 0.10, sprintf('white level via Allan within 10%% on all 9 channels (max err %.1f%%)', 100 * eA), nPass, nFail);
ex = [r.cohBasinAbove95, r.cohMizAbove95];                  % pooled over axes, pairs and both bands
[nPass, nFail] = check(mean(ex) < 0.15, sprintf(['independent sensors: %.1f%% of band bins above the 95%% threshold ' ...
    '(about 5%% expected by chance)'], 100 * mean(ex)), nPass, nFail);
[nPass, nFail] = check(exist(fullfile(work, 'Results', 'staticNoise_imu_synthetic_10Hz_long_adev_acc.png'), 'file') == 2 && ...
    exist(fullfile(work, 'Results', 'staticNoise_imu_synthetic_10Hz_long.txt'), 'file') == 2, 'summary and figures written', nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testStaticNoiseCharacterisation:failed', '%d test(s) failed.', nFail); end

function runDriver(script, STATIC_NOISE_CONFIG) %#ok<INUSD>
% run the driver script in this function's own workspace (it clears variables)
run(script);
end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end