function A = alignImuSensors(rec, varargin)
%ALIGNIMUSENSORS  Match the co-located sensors of one record tick by tick.
%   A = alignImuSensors(rec)                 rec from loadImuRecord
%   A = alignImuSensors(rec, 'Tolerance', 0.02, 'Reference', 3)
%
%   The three IMUs are read in the same loop tick, about 1 ms apart. For every row of the
%   reference sensor, the nearest row of each other sensor within Tolerance (s; default 0.2 of
%   the median interval) is taken; ticks where any sensor is missing are dropped. Returns
%     A.t      n x 1 reference times (s, logging path)
%     A.acc    n x 3 x S, A.gyr n x 3 x S   (sample x axis x sensor)
%     A.names  sensor names, A.index n x S row indices into rec.sensors(k), A.nDropped
%   e.g. the z-accelerations of all sensors: squeeze(A.acc(:, 3, :))   (n x S)
%
% Lives in IMU Characterisation/Processing/.

S = rec.sensors; nS = numel(S);
ref = nS; tol = [];
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'tolerance', tol = varargin{i + 1};
        case 'reference', ref = varargin{i + 1};
        otherwise, error('alignImuSensors:args', 'Unknown option ''%s''.', varargin{i});
    end
end
tr = S(ref).t;
if isempty(tol), tol = 0.2 * median(diff(tr)); end
idx = zeros(numel(tr), nS); ok = true(numel(tr), 1);
for k = 1:nS
    if k == ref, idx(:, k) = (1:numel(tr)).'; continue; end
    j = interp1(S(k).t, 1:numel(S(k).t), tr, 'nearest', 'extrap');
    j = max(1, min(numel(S(k).t), round(j)));
    idx(:, k) = j;
    ok = ok & abs(S(k).t(j) - tr) <= tol;
end
idx = idx(ok, :);
for k = 1:nS                                                  % one-to-one: drop re-used rows
    [~, u] = unique(idx(:, k), 'stable');
    keep = false(size(idx, 1), 1); keep(u) = true; idx = idx(keep, :);
end
n = size(idx, 1);
A.t = tr(idx(:, ref));
A.acc = zeros(n, 3, nS); A.gyr = zeros(n, 3, nS);
for k = 1:nS
    A.acc(:, :, k) = S(k).acc(idx(:, k), :);
    A.gyr(:, :, k) = S(k).gyr(idx(:, k), :);
end
A.names = {S.name}; A.index = idx; A.nDropped = numel(tr) - n;
end