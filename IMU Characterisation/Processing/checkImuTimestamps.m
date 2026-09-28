function out = checkImuTimestamps(tau, varargin)
%CHECKIMUTIMESTAMPS  Timing audit of one IMU record: rate, jitter, drift, gaps, duplicates.
%   out = checkImuTimestamps(tau, 'NominalRate', fs, Name, Value, ...)
%
%   tau   device timestamps (s), one per received sample, in arrival order (may wrap).
%
%   Name-value options
%     'NominalRate'  requested output rate fs (Hz). Required.
%     'Counter'      sample counter per row (optional; strongly preferred: without it,
%                    gaps are inferred from the timestamps alone)
%     'CounterBits'  counter width in bits (default 16)
%     'TsWrap'       timestamp wrap period (s), e.g. 2^32 * 1e-6 for a 32-bit microsecond
%                    counter (default Inf; negative jumps are then reported, not unwrapped)
%     'Host'         PC arrival times (s) per row (optional): fits host vs device clock
%     'Windows'      number of windows for the windowed rate / drift estimate (default 10)
%     'Verbose'      print a summary (default true)
%
%   METHOD
%     1  unwrap timestamps and counter; per-row step s_i in true sample index
%        (counter difference, corrected by whole wraps using the timestamp difference;
%        without a counter, round(dtau / T) with T refined from a first fit):
%        s = 0 duplicate, s = 1 normal, s > 1 gap of s - 1 missing samples
%     2  least-squares fit tau_n = a + b n over the reconstructed index n (duplicates removed):
%        actual rate 1 / b, rate error in ppm against fs, timing residual e_n (jitter plus
%        any clock wander)
%     3  windowed slopes and a quadratic fit: how much the rate changes over the record
%     4  host clock: host = c + d tau, relative rate (d - 1) in ppm, latency spread
%   All rates are in the DEVICE clock's units: an absolute check needs an independent
%   reference (GNSS PPS or hours of host time), see the week plan.
%
%   Output struct fields: nRows, nSamples (reconstructed span), nDuplicates, duplicateRows,
%   nGaps, nMissing, gaps (table rows [rowBefore, firstMissingIndex, missingCount, tau]),
%   nBackwardJumps, nWraps, index (reconstructed n per row), period, rate, rateErrorPpm,
%   offset, residual (per kept row), jitterStd, jitterMaxAbs, jitterP99, apparentResolution
%   (smallest spacing between distinct timestamp increments: it estimates the timestamp
%   quantum only when jitter or a rate offset makes the increments visit adjacent levels, may
%   return a multiple of it, and is NaN when every increment is identical; take the true
%   unit from the firmware or register map),
%   windowRatePpm, driftPpm (quadratic: rate change end minus start), elapsed,
%   nominalSamplesFromDeviceTime (elapsed device time x fs + 1: a rough check only, it differs
%   from nSamples by the rate error), and host (struct) when Host is given.
%
% Lives in IMU Characterisation/Processing/.

opt = struct('NominalRate', [], 'Counter', [], 'CounterBits', 16, 'TsWrap', Inf, 'Host', [], ...
    'Windows', 10, 'Verbose', true);
if mod(numel(varargin), 2) ~= 0, error('checkImuTimestamps:args', 'Name-value pairs expected.'); end
for i = 1:2:numel(varargin)
    if ~isfield(opt, varargin{i}), error('checkImuTimestamps:args', 'Unknown option ''%s''.', varargin{i}); end
    opt.(varargin{i}) = varargin{i + 1};
end
if isempty(opt.NominalRate), error('checkImuTimestamps:args', 'NominalRate is required.'); end
tau = double(tau(:)); nR = numel(tau);
T0 = 1 / opt.NominalRate;

% ---- 1  unwrap -----------------------------------------------------------------------------
d = diff(tau);
nWraps = 0;
if isfinite(opt.TsWrap)
    w = d < -opt.TsWrap / 2;
    nWraps = nnz(w);
    tau = tau + opt.TsWrap * [0; cumsum(w)];
    d = diff(tau);
end
nBack = nnz(d < 0);

if ~isempty(opt.Counter)
    c = double(opt.Counter(:));
    M = 2^opt.CounterBits;
    dc = mod(diff(c), M);
    steps = dc + M * round((d / T0 - dc) / M);            % whole counter wraps from the timestamps
    usedCounter = true;
else
    usedCounter = false;
    Tfit = T0;
    for pass = 1:2                                          % refine T so large gaps round correctly
        steps = round(d / Tfit);
        n = [0; cumsum(max(steps, 0))];
        keepF = [true; steps > 0];
        p = polyfit(n(keepF), tau(keepF), 1);
        Tfit = p(1);
    end
end
steps = max(steps, 0);                                      % backward jumps are counted separately
n = [0; cumsum(steps)];
dupRows = find(steps == 0) + 1;
keep = true(nR, 1); keep(dupRows) = false;
gapRows = find(steps > 1);
gaps = [gapRows, n(gapRows) + 1, steps(gapRows) - 1, tau(gapRows)];

% ---- 2  linear clock fit ---------------------------------------------------------------------
nk = n(keep); tk = tau(keep);
nc = nk - mean(nk);
b = (nc' * (tk - mean(tk))) / (nc' * nc);
a = mean(tk) - b * mean(nk);
r = tk - (a + b * nk);
du = unique(round(diff(tk) / 1e-9) * 1e-9);
res = min(diff(du));
if isempty(res), res = NaN; end                              % all increments identical: undetermined

% ---- 3  drift ------------------------------------------------------------------------------------
W = max(1, min(opt.Windows, floor(numel(nk) / 20)));
edges = round(linspace(1, numel(nk) + 1, W + 1));
winPpm = zeros(1, W);
for k = 1:W
    ii = edges(k):edges(k + 1) - 1;
    pk = polyfit(nk(ii) - nk(ii(1)), tk(ii), 1);
    winPpm(k) = 1e6 * (1 / (pk(1) * opt.NominalRate) - 1);
end
q = polyfit(nc / max(abs(nc)), tk, 2);                      % tau = q1 x^2 + q2 x + q3, x in [-1, 1]
x0 = nc(1) / max(abs(nc)); x1 = nc(end) / max(abs(nc));      % actual first and last samples
slopeStart = q(2) + 2 * q(1) * x0; slopeEnd = q(2) + 2 * q(1) * x1;
driftPpm = 1e6 * (slopeStart / slopeEnd - 1);               % rate(end) / rate(start) - 1

out = struct('nRows', nR, 'usedCounter', usedCounter, 'nSamples', n(end) + 1, ...
    'nDuplicates', numel(dupRows), 'duplicateRows', dupRows, 'nGaps', numel(gapRows), ...
    'nMissing', sum(steps(gapRows) - 1), 'gaps', gaps, 'nBackwardJumps', nBack, 'nWraps', nWraps, ...
    'index', n, 'keep', keep, 'period', b, 'rate', 1 / b, 'rateErrorPpm', 1e6 * (1 / (b * opt.NominalRate) - 1), ...
    'offset', a, 'residual', r, 'jitterStd', std(r), 'jitterMaxAbs', max(abs(r)), ...
    'jitterP99', prctile1(abs(r), 99), 'apparentResolution', res, 'windowRatePpm', winPpm, 'driftPpm', driftPpm, ...
    'elapsed', tk(end) - tk(1), 'nominalSamplesFromDeviceTime', round((tk(end) - tk(1)) * opt.NominalRate) + 1);

% ---- 4  host clock -------------------------------------------------------------------------------
if ~isempty(opt.Host)
    h = double(opt.Host(:)); hk = h(keep);
    tcen = tk - mean(tk);
    dh = (tcen' * (hk - mean(hk))) / (tcen' * tcen);
    lat = hk - (mean(hk) + dh * tcen);
    lat = lat - min(lat);                                   % latency above the fastest arrival
    out.host = struct('relativePpm', 1e6 * (dh - 1), 'latencyMedian', median(lat), ...
        'latencyP99', prctile1(lat, 99), 'latencyMax', max(lat), 'nonMonotonic', nnz(diff(h) < 0));
end

if opt.Verbose, report(out, opt); end
end

function p = prctile1(x, q)
x = sort(x(:));
if isempty(x), p = NaN; return; end
p = x(max(1, min(numel(x), ceil(q / 100 * numel(x)))));
end

function report(o, opt)
fprintf('checkImuTimestamps: %d rows, nominal %.6g Hz, %s\n', o.nRows, opt.NominalRate, ...
    ternary(o.usedCounter, 'gaps from the sample counter', 'NO counter: gaps inferred from timestamps'));
fprintf('  elapsed %.3f s (device clock); nominal count from elapsed time %d, reconstructed span %d, received %d\n', o.elapsed, ...
    o.nominalSamplesFromDeviceTime, o.nSamples, o.nRows);
fprintf('  duplicates %d; gaps %d (%d samples missing); backward jumps %d; timestamp wraps %d\n', ...
    o.nDuplicates, o.nGaps, o.nMissing, o.nBackwardJumps, o.nWraps);
fprintf('  actual rate %.6f Hz (%+.3f ppm vs nominal, in device-clock units)\n', o.rate, o.rateErrorPpm);
fprintf('  timing residual after the linear fit: std %.3g us, 99%% %.3g us, max %.3g us; apparent resolution %.3g us\n', ...
    1e6 * o.jitterStd, 1e6 * o.jitterP99, 1e6 * o.jitterMaxAbs, 1e6 * o.apparentResolution);
fprintf('  windowed rate error (ppm): %s; quadratic drift over record %+.3f ppm\n', ...
    mat2str(round(o.windowRatePpm * 1000) / 1000), o.driftPpm);
if isfield(o, 'host')
    fprintf('  host vs device clock %+.2f ppm; host latency median %.3g ms, 99%% %.3g ms, max %.3g ms\n', ...
        o.host.relativePpm, 1e3 * o.host.latencyMedian, 1e3 * o.host.latencyP99, 1e3 * o.host.latencyMax);
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end