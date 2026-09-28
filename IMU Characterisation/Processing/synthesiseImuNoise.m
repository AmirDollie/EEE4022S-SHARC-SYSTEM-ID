function rec = synthesiseImuNoise(duration, fs, varargin)
%SYNTHESISEIMUNOISE  Synthetic IMU record with known timing and noise properties.
%   rec = synthesiseImuNoise(duration, fs, Name, Value, ...)
%
%   Test data for the IMU processing pipeline: every property the pipeline estimates
%   (actual rate, clock errors, jitter, gaps, duplicates, counter wraps, noise density,
%   flicker and random-walk coefficients) is set here and returned in rec.truth.
%
%   TIMING MODEL (three clocks)
%     true sample times   t_n = n / (fs (1 + OdrPpm 1e-6)),  n = 0 .. N-1   (IMU oscillator)
%     device timestamps   tau_n = TsOffset + (1 + TsPpm 1e-6) t_n + jitter, quantised to
%                         TsResolution and wrapped modulo TsWrap                (MCU clock)
%     sample counter      n modulo 2^CounterBits
%     host arrival times  (1 + HostPpm 1e-6) t_n + latency, latency >= 0 with exponential
%                         spread and occasional bursts                          (PC clock)
%   The pipeline can only see tau (and the host times): the apparent sample period in
%   device time is (1 + TsPpm 1e-6) / (fs (1 + OdrPpm 1e-6)) (rec.truth.periodDevice).
%
%   NOISE MODEL (each channel = constant + white + flicker + random walk)
%     white          ONE-SIDED amplitude density N (units / sqrt(Hz), the datasheet convention and
%                    the transmissibility estimator's): variance N^2 fs / 2, one-sided PSD N^2,
%                    Allan deviation N / sqrt(2 tau). (IEEE 952 quotes N with the two-sided PSD,
%                    giving N / sqrt(tau): check which convention a datasheet uses.)
%     flicker        one-sided PSD B^2 / (2 pi f), B a COEFFICIENT of this generator; the
%                    Allan plateau this implies (sqrt(ln 2 / pi) B = 0.470 B for this one-sided
%                    form, since sigma_A^2 = 2 ln2 h for S = h / f) is to be confirmed by the
%                    allanDeviation test before B is compared with any datasheet "bias instability"
%     random walk    increments K sqrt(dt);  Allan deviation K sqrt(tau / 3)
%   Accelerometer [ax ay az] (m/s^2, gravity on z) and gyroscope [gx gy gz] (rad/s) get
%   their own [whiteDensity, flickerCoefficient, randomWalkCoefficient] triples.
%
%   DEFECTS (applied after generation, in the order drop -> duplicate)
%     'Drop'       k x 2 array of [firstIndex count]: samples removed (counter jumps)
%     'Duplicate'  vector of indices repeated once (same counter, same data)
%
%   Name-value options (defaults in brackets)
%     'OdrPpm' [0]  'TsPpm' [0]  'TsDriftPpm' [0] (device-clock rate change start to end)  'TsOffset' [0]  'TsJitter' [0] (s, std)  'TsResolution' [1e-6]
%     'TsWrap' [Inf] (s)  'CounterBits' [16]  'HostPpm' [0]  'HostLatency' [2e-3] (s, mean)
%     'HostBurstProb' [0.01]  'HostBurstDelay' [20e-3]
%     'Accel' [N B K] [2e-3 5e-4 1e-5]   'Gyro' [N B K] [1e-4 2e-5 5e-7]   (white, flicker, random walk)
%     'Gravity' [9.81]  'Drop' [zeros(0,2)]  'Duplicate' [[]]  'Seed' [1]
%
%   Output rec: tau (device time, wrapped), counter (wrapped), host, acc (N x 3), gyr (N x 3),
%   index (true sample index of each row), truth (all settings plus periodDevice and the
%   lists of dropped and duplicated true indices).
%
% Lives in IMU Characterisation/Processing/.

opt = struct('OdrPpm', 0, 'TsPpm', 0, 'TsDriftPpm', 0, 'TsOffset', 0, 'TsJitter', 0, 'TsResolution', 1e-6, 'TsWrap', Inf, ...
    'CounterBits', 16, 'HostPpm', 0, 'HostLatency', 2e-3, 'HostBurstProb', 0.01, 'HostBurstDelay', 20e-3, ...
    'Accel', [2e-3 5e-4 1e-5], 'Gyro', [1e-4 2e-5 5e-7], 'Gravity', 9.81, 'Drop', zeros(0, 2), ...
    'Duplicate', [], 'Seed', 1);
opt = parseOptions(opt, varargin);
rng(opt.Seed);

N = round(duration * fs);
n = (0:N - 1).';
t = n / (fs * (1 + 1e-6 * opt.OdrPpm));
dtTrue = 1 / (fs * (1 + 1e-6 * opt.OdrPpm));

acc = zeros(N, 3); gyr = zeros(N, 3);
for c = 1:3
    acc(:, c) = channelNoise(N, dtTrue, opt.Accel);
    gyr(:, c) = channelNoise(N, dtTrue, opt.Gyro);
end
acc(:, 3) = acc(:, 3) + opt.Gravity;

tau = opt.TsOffset + (1 + 1e-6 * opt.TsPpm) * t + 0.5e-6 * opt.TsDriftPpm * t.^2 / max(t(end), eps) ...
    + opt.TsJitter * randn(N, 1);
tau = opt.TsResolution * round(tau / opt.TsResolution);
lat = opt.HostLatency * (0.5 + 0.5 * exprnd1(N));
burst = rand(N, 1) < opt.HostBurstProb;
lat(burst) = lat(burst) + opt.HostBurstDelay * rand(nnz(burst), 1);
host = (1 + 1e-6 * opt.HostPpm) * t + lat;
host = max(host, cummax1(host));                   % a serial stream arrives in order

keep = true(N, 1);
for r = 1:size(opt.Drop, 1)
    keep(opt.Drop(r, 1):min(N, opt.Drop(r, 1) + opt.Drop(r, 2) - 1)) = false;
end
idx = find(keep);
dupIdx = intersect(opt.Duplicate(:), idx);
rows = sort([idx; dupIdx]);                        % a duplicate follows its original

rec.index = rows - 1;
rec.tau = tau(rows);
if isfinite(opt.TsWrap), rec.tau = mod(rec.tau, opt.TsWrap); end
rec.counter = mod(rows - 1, 2^opt.CounterBits);
rec.host = host(rows);
rec.acc = acc(rows, :);
rec.gyr = gyr(rows, :);
rec.fs = fs;
truth = opt;
truth.duration = duration; truth.N = N; truth.fs = fs;
truth.periodDevice = (1 + 1e-6 * opt.TsPpm) * dtTrue;
truth.rateDevice = 1 / truth.periodDevice;
truth.rateErrorPpm = 1e6 * (truth.rateDevice / fs - 1);   % actual rate in device time vs requested fs
truth.droppedIndex = find(~keep) - 1;
truth.duplicatedIndex = dupIdx - 1;
truth.hostVsDevicePpm = 1e6 * ((1 + 1e-6 * opt.HostPpm) / (1 + 1e-6 * opt.TsPpm) - 1);
rec.truth = truth;
end

function x = channelNoise(N, dt, p)
% white + flicker (1/f) + random walk, one channel
fs = 1 / dt;
x = p(1) * sqrt(fs / 2) * randn(N, 1);             % one-sided density p(1): variance p(1)^2 fs / 2
if p(2) > 0
    x = x + flicker(N, fs, p(2));
end
if p(3) > 0
    x = x + cumsum(p(3) * sqrt(dt) * randn(N, 1));
end
end

function y = flicker(N, fs, B)
% one-sided PSD B^2 / (2 pi f) by spectral shaping of white noise
M = 2^nextpow2(2 * N);
X = fft(randn(M, 1));
f = (0:M - 1).' * fs / M;
f(f > fs / 2) = fs - f(f > fs / 2);                % symmetric frequency magnitude
H = zeros(M, 1);
H(f > 0) = sqrt(B^2 * fs / (4 * pi) ./ f(f > 0));
y = real(ifft(X .* H));
y = y(1:N) - mean(y(1:N));
end

function e = exprnd1(N)
e = -log(rand(N, 1));
end

function c = cummax1(x)
c = x;
for i = 2:numel(x), if c(i) < c(i - 1), c(i) = c(i - 1); end, end
end

function opt = parseOptions(opt, args)
if mod(numel(args), 2) ~= 0, error('synthesiseImuNoise:args', 'Name-value pairs expected.'); end
for i = 1:2:numel(args)
    if ~isfield(opt, args{i}), error('synthesiseImuNoise:args', 'Unknown option ''%s''.', args{i}); end
    opt.(args{i}) = args{i + 1};
end
end