function [still, m] = isStaticHold(acc, gyr, varargin)
%ISSTATICHOLD  Is a calibration hold usable, and did the sensor move during it?
%   [still, m] = isStaticHold(acc, gyr)
%   [still, m] = isStaticHold(acc, gyr, Name, Value, ...)
%
%   acc  N x 3 (m/s^2), gyr N x 3 (rad/s) or [] if unavailable.
%   A slow rotation keeps |a| = g, so a magnitude test alone cannot see it. The calibration uses only
%   the hold's MEAN vector, so what matters is whether motion corrupts that mean:
%
%   EXCLUSION criteria (hold unusable, still = false)
%     MagStd    std of |a|                                        (m/s^2)  vibration, bumps
%     VecStd    norm of the per-axis std of a                     (m/s^2)  rocking
%     GyroRms   RMS of |omega - mean omega|                       (rad/s)  rotation, if gyro present
%     AvgLoss   mean(|a|) - |mean(a)|, relative to mean(|a|)      (-)      how much averaging over a
%               changing direction shortens the mean vector; this is the error the sphere fit sees
%   WARNING criterion (reported, hold kept)
%     Drift     |mean(last quarter) - mean(first quarter)|        (m/s^2)  slow rotation or settling
%               (0.02 m/s^2 ~ 0.1 deg). A slow 0.3 deg drift shortens the mean by ~1e-6 relative:
%               harmless for the fit, but it means the hold was not mechanically settled.
%
%   Defaults sit above the measured 10 Hz noise of all three IMUs (LSM6DS3TR-C per-axis std about
%   0.02 m/s^2); the gyro limit is loose because the shell prints 0.001 rad/s.
%   m returns the measured values, m.fails (exclusion) and m.warnings.
%
% Lives in IMU Characterisation/Processing/.

opt = struct('MagStd', 0.05, 'VecStd', 0.08, 'GyroRms', 0.01, 'AvgLoss', 1e-4, 'Drift', 0.02);
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i + 1}; end
n = size(acc, 1); q = max(1, floor(n / 4));
mag = sqrt(sum(acc.^2, 2));
m.magStd = std(mag);
m.vecStd = norm(std(acc, 0, 1));
m.avgLoss = (mean(mag) - norm(mean(acc, 1))) / mean(mag);
m.drift = norm(mean(acc(end - q + 1:end, :), 1) - mean(acc(1:q, :), 1));
m.gyroRms = NaN;
if ~isempty(gyr)
    w = gyr - mean(gyr, 1);
    m.gyroRms = sqrt(mean(sum(w.^2, 2)));
end
fails = {};
if m.magStd > opt.MagStd, fails{end + 1} = 'magnitude'; end
if m.vecStd > opt.VecStd, fails{end + 1} = 'vector'; end
if m.avgLoss > opt.AvgLoss, fails{end + 1} = 'averaging'; end
if ~isnan(m.gyroRms) && m.gyroRms > opt.GyroRms, fails{end + 1} = 'gyro'; end
warns = {};
if m.drift > opt.Drift, warns{end + 1} = 'drift'; end
m.fails = fails;
m.warnings = warns;
still = isempty(fails);
end