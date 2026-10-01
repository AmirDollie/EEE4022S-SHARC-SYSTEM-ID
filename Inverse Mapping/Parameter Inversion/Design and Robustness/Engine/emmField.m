function fld = emmField(p, sensors, varargin)
%EMMFIELD  EMM acceleration FRF at fixed sensor points, and its derivatives in ln-parameters.
%   fld = emmField(p, sensors)
%   fld = emmField(p, sensors, 'ThetaIdx', [1 3], 'NodeSpacing', 0.2, 'CacheDir', d, ...)
%
%   p        [beta gamma R], EMM scaling.
%   sensors  nS x 2 [r theta], NON-DIMENSIONAL r (r / H), FIXED physical points. They are NOT
%            moved when R is perturbed: dH/dlnR is the response change of the same points.
%            Any number of points (e.g. the whole 33-point candidate grid): one EMM solve per
%            node frequency serves every point.
%
%   Method (no new physics): the node FRF of synthesiseTwinRecords (EMM solved on a node grid
%   'NodeSpacing' apart over the validated band, H/(-omega^2) spline-interpolated, cached on
%   disk) at p and at p .* exp(+-h e_i) for each i in ThetaIdx. Derivatives are central
%   differences of the interpolated FRF, so they are evaluated at any frequency in the band.
%
%   fld.H(w)       nS x nW complex acceleration FRF at p
%   fld.dH(w)      nS x nW x nTheta, dH / d ln p(ThetaIdx(i))
%   fld.twin       prepared twin at p (fld.twin.Hacc == fld.H), for the Welch-expected bias
%   fld.p, sensors, thetaIdx, h, nodeSpacing, band, cacheStatus (1 + 2 nTheta entries)
%
%   Options: 'ThetaIdx' [1 3] (which ln-parameters to differentiate; evaluateScenario passes the
%   estimated ones first, then any nuisance ones, e.g. [1 3 2] for R4), 'Step' 1e-3,
%   'NodeSpacing' 0.2 (rad/s), 'FRFOptions' {} (e.g. {'Truncation', [70 15 15]}), 'Band' []
%   (validated band), 'CacheDir' '' (no cache), 'CheckPoints' 3 at p, 0 at the perturbed points,
%   'FRF' [] or a handle @(omega, p) -> nS x nW acceleration FRF that REPLACES the EMM (analytic
%   test fields; not cached). Everything downstream is identical.
%   'Verbose' false; true prints one line per node FRF (1 + 2 nTheta of them): which point, whether
%   it was loaded from the cache or solved, and how long it took.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Engine/.

opt = struct('ThetaIdx', [1 3], 'Step', 1e-3, 'NodeSpacing', 0.2, 'FRFOptions', {{}}, 'Band', [], ...
    'CacheDir', '', 'CheckPoints', 3, 'FRF', [], 'Verbose', false);
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i + 1}; end
p = p(:).';
if ~(numel(p) == 3 && all(p > 0) && all(isfinite(p))), error('emmField:p', 'p must be [beta gamma R] > 0.'); end
if any(sensors(:, 1) > p(3) * (1 - 1e-9)) || any(sensors(:, 1) < 0)
    error('emmField:sensorOutsideFloe', 'All sensor radii must lie in [0, R).');
end
if isempty(opt.FRF)
    spec = @(pp, cp) struct('p', pp, 'sensors', sensors, 'nodeSpacing', opt.NodeSpacing, ...
        'frfOptions', {opt.FRFOptions}, 'band', opt.Band, 'checkPoints', cp);
    cacheArgs = {};
    if ~isempty(opt.CacheDir), cacheArgs = {'CacheDir', opt.CacheDir}; end
else
    band = opt.Band; if isempty(band), band = [2.9784 8.5012]; end
    spec = @(pp, cp) struct('p', pp, 'sensors', sensors, 'nodeSpacing', opt.NodeSpacing, ...
        'band', band, 'checkPoints', cp, 'frf', @(w) opt.FRF(w, pp));
    cacheArgs = {};
end

nT = numel(opt.ThetaIdx);
nTot = 1 + 2 * nT;
pn = {'beta', 'gamma', 'R'};
tc = tic;
twin0 = synthesiseTwinRecords(spec(p, opt.CheckPoints), cacheArgs{:});
report(opt.Verbose, 1, nTot, 'nominal', twin0, toc(tc));
twP = cell(1, nT); twM = cell(1, nT);
status = cell(1, nTot); status{1} = twin0.cacheStatus;
for i = 1:nT
    e = zeros(1, 3); e(opt.ThetaIdx(i)) = opt.Step;
    tc = tic;
    twP{i} = synthesiseTwinRecords(spec(p .* exp(e), 0), cacheArgs{:});
    report(opt.Verbose, 2 * i, nTot, sprintf('ln %s + h', pn{opt.ThetaIdx(i)}), twP{i}, toc(tc));
    tc = tic;
    twM{i} = synthesiseTwinRecords(spec(p .* exp(-e), 0), cacheArgs{:});
    report(opt.Verbose, 2 * i + 1, nTot, sprintf('ln %s - h', pn{opt.ThetaIdx(i)}), twM{i}, toc(tc));
    status{2 * i} = twP{i}.cacheStatus; status{2 * i + 1} = twM{i}.cacheStatus;
end
h = opt.Step;
fld.H = @(w) twin0.Hacc(w);
fld.dH = @(w) derivs(w, twP, twM, h, size(sensors, 1));
fld.twin = twin0;
fld.p = p; fld.sensors = sensors; fld.thetaIdx = opt.ThetaIdx; fld.h = h;
fld.nodeSpacing = opt.NodeSpacing; fld.band = twin0.band; fld.cacheStatus = status;
end

function report(verbose, i, n, what, tw, sec)
if ~verbose, return; end
fprintf('      emmField %d/%d  %-12s %-8s %3d nodes  %6.1f s\n', i, n, what, tw.cacheStatus, numel(tw.omegaNodes), sec);
if exist('OCTAVE_VERSION', 'builtin'), fflush(stdout); else, drawnow; end   % show the line immediately
end

function D = derivs(w, twP, twM, h, nS)
nT = numel(twP);
D = complex(zeros(nS, numel(w), nT));
for i = 1:nT
    D(:, :, i) = (twP{i}.Hacc(w) - twM{i}.Hacc(w)) / (2 * h);
end
end