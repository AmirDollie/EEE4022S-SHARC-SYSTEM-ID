function T = m1Tools()
%M1TOOLS  The M1 machinery around the G3 code, shared by runM1MonteCarloSpotChecks and
%   testM1MonteCarloSpotChecks so the tester exercises exactly the code the run uses.
%   T = m1Tools();   then e.g.  ctx = T.context(cases(q), cfg)
%
%   No new estimator: records come from synthesiseTwinRecords on the engine's field twin, the estimator,
%   information matrix, gradient, objective, Jacobian and Gauss-Newton solver are G3's own
%   (g3EstimatorCore, verbatim). This file only wires them to the analytic (E1) prediction of a scenario.
%
%   WHAT IS TESTED. For one scenario the engine predicts, for one record, F, sigma = sqrt(diag F^-1) and
%   the Welch-expected bias dtheta_W. M1 generates records and applies the frozen estimator:
%     z        corrected features of the record (welchSegments + zFromSegments, G3)
%     dtheta   = F^-1 J' Sigma^-1 (z - f0) on the accepted bins: the linearised estimate minus the truth
%              (G3a, with the ANALYTIC J, Sigma_k and F of evaluateScenario)
%     q        = dtheta' F dtheta (G3a: centred on the truth; E[q] = 2 + d_W^2)
%   and compares the scatter with sigma, the coverage of q with chi2(2), and the mean with dtheta_W.
%
%   FUNCTIONS
%   ctx = T.context(case, cfg)         evaluateScenario on the case -> the analytic target (ctx.tgt, saved)
%                                      and everything needed to make and estimate records (ctx, not saved)
%   ok  = T.compareSource(tgt, source, cfg)   does the target reproduce the source study? (struct)
%   s   = T.seed(cfg, caseIndex, record)      the record seed (noise seed = s + 1e6, synthesiseTwinRecords)
%   e   = T.estimateRecord(ctx, seed)         one record: z, dtheta (linearised), q, rejected bins
%   e   = T.linearisedFromZ(tgt, G, z)        the linearised step alone, from features z (p x nF)
%   [b, what] = T.isBorderline(empPred, cov95, n, cfg)   the pre-registered extension rule
%   st  = T.newState(ctx, cfg)                empty Monte Carlo / nonlinear state for a case
%   st  = T.mcBatch(ctx, st, nTarget, ck)     records st.nDone+1 .. nTarget; checkpoint every ck.every
%                                             to ck.file (ck.file '' = none); ck.verbose prints progress
%   sm  = T.summarise(st, tgt, cfg)           emp/pred, mean offset, coverage, borderline, pass flags
%   fun = T.makeModel(ctx)                    the nonlinear model theta -> features on ALL bins
%                                             (emmField at the trial point, same node grid, truncation,
%                                             nu and band as the target; memoised)
%   nl  = T.nonlinear(ctx, st, betaOffset, fun, cfg)   G3 solver from beta x (1 + offset), R true
%   key = T.checkpointKey(ctx, cfg)           everything a checkpoint's validity depends on
%   [st, ok, msg] = T.loadCheckpoint(file, key)
%   T.saveCheckpoint(file, st)
%   [clean, files] = T.cleanTree(porcelainText)   uncommitted .m files in `git status --porcelain` output
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
T.context = @context;
T.compareSource = @compareSource;
T.seed = @seedOf;
T.estimateRecord = @estimateRecord;
T.linearisedFromZ = @linearisedFromZ;
T.isBorderline = @isBorderline;
T.newState = @newState;
T.mcBatch = @mcBatch;
T.summarise = @summarise;
T.makeModel = @makeModel;
T.nonlinear = @nonlinear;
T.checkpointKey = @checkpointKey;
T.loadCheckpoint = @loadCheckpoint;
T.saveCheckpoint = @saveCheckpoint;
T.cleanTree = @cleanTree;
T.chi2q = [-2 * log(1 - 0.6827), -2 * log(0.05)];       % chi2(2) quantiles, as G3
T.version = 'm1Tools 2026-10-03';
end

%% ================================================================================================
%  Analytic target and record context
%% ================================================================================================
function ctx = context(cs, cfg)
out = evaluateScenario(cs.scn, cfg);
scn = out.scn;
nL = numel(scn.layout); others = setdiff(1:nL, scn.ref); p = 2 * numel(others);
nF = numel(out.omega);
valid = out.validMask(:).';
rows = reshape((find(valid) - 1) * p + (1:p).', [], 1);
tgt = struct('name', cs.name, 'index', cs.index, 'thetaTrue', log(scn.p([1 3])).', 'p', scn.p, ...
    'sigma', out.sigma(:), 'corr', out.corr, 'condF', out.condF, 'F', out.F, 'C', out.covTheta, ...
    'dthetaW', out.dthetaW(:), 'dW', out.dW, 'sMax', out.sMax, 'class', out.class, 'singular', out.singular, ...
    'omega', out.omega, 'kBins', out.kBins, 'validMask', valid, 'nUsedBins', out.nUsedBins, ...
    'invalidFraction', out.invalidFraction, 'f0', out.f0, 'J', out.Jall, 'blocks', out.blocksAll, ...
    'rows', rows, 'pFeat', p, 'nF', nF, 'layout', scn.layout, 'ref', scn.ref, 'sensors', scn.sensors, ...
    'Snn', scn.Snn, 'N', scn.N, 'dt', scn.dt, 'L', scn.L, 'overlap', scn.overlap, 'bandUse', scn.bandUse, ...
    'nodeSpacing', scn.nodeSpacing, 'frfOptions', {scn.frfOptions}, 'seaText', seaText(scn.sea), 'seconds', out.seconds);
if numel(scn.Snn) > 1, error('m1Tools:Snn', 'M1 assumes one noise PSD for every sensor.'); end
ctx = struct('name', cs.name, 'index', cs.index, 'tgt', tgt, 'scn', scn, 'twin', out.field.twin, ...
    'others', others, 'pBlk', p, 'specArgs', {seaArgs(scn.sea)}, 'noiseStd', sqrt(pi * scn.Snn / scn.dt), ...
    'G', g3EstimatorCore());
end

function r = compareSource(tgt, src, cfg)
r.relSigma = max(abs(tgt.sigma(:).' ./ src.sigma(:).' - 1));
r.dDW = abs(tgt.dW - src.dW);
r.sameClass = strcmp(tgt.class, src.class);
r.ok = r.relSigma <= cfg.M1.sourceRelTol && r.dDW <= cfg.M1.sourceRelTol * max(1, src.dW) && r.sameClass;
end

function s = seedOf(cfg, caseIndex, record)
s = cfg.M1.seedBase + cfg.M1.seedStride * caseIndex + record;
end

%% ================================================================================================
%  One record
%% ================================================================================================
function e = estimateRecord(ctx, seed)
t = ctx.tgt; G = ctx.G;
Y = synthesiseTwinRecords(ctx.twin, t.N, t.dt, 'Seed', seed, ctx.specArgs{:}, 'NoiseStd', ctx.noiseStd);
Y = Y(t.layout, :);
Pw = G.welchSegments(Y, t.dt, t.L, t.overlap, t.kBins, t.ref, ctx.others);
z = G.zFromSegments(Pw, 1:Pw.K, t.Snn);                       % p x nF, NaN = bin rejected by the estimator
e = linearisedFromZ(t, G, z);
e.seed = seed;
end

function e = linearisedFromZ(t, G, z)
% linearised estimate from corrected features z (p x nF) on the analytic accepted bins this record kept
use = t.validMask & all(isfinite(z), 1);
nBad = nnz(t.validMask & ~use);
rows = reshape((find(use) - 1) * t.pFeat + (1:t.pFeat).', [], 1);
r = z(:, use) - reshape(t.f0(rows), t.pFeat, []);
if nBad == 0
    F = t.F;
else
    F = G.infoMatrix(t.J(rows, :), t.blocks(:, :, use));     % this record lost bins: its own information
end
dth = F \ G.gradTerm(t.J(rows, :), t.blocks(:, :, use), r);
e = struct('seed', NaN, 'z', z, 'dtheta', dth, 'q', dth.' * t.F * dth, 'nBad', nBad, 'use', use);
end

%% ================================================================================================
%  Monte Carlo state and batches
%% ================================================================================================
function st = newState(ctx, cfg)
nMax = cfg.M1.nRecordsExtended;
st = struct('name', ctx.name, 'index', ctx.index, 'key', checkpointKey(ctx, cfg), 'nDone', 0, ...
    'nTarget', cfg.M1.nRecords, 'extended', false, 'seeds', NaN(1, nMax), 'dtheta', NaN(2, nMax), ...
    'q', NaN(1, nMax), 'nBad', NaN(1, nMax), 'seconds', NaN(1, nMax), 'nl', struct([]), ...
    'commits', {{}}, 'decision', '');
end

function st = mcBatch(ctx, st, nTarget, ck)
if nargin < 4 || isempty(ck), ck = struct(); end
if ~isfield(ck, 'file'), ck.file = ''; end
if ~isfield(ck, 'every'), ck.every = 10; end
if ~isfield(ck, 'verbose'), ck.verbose = false; end
if ~isfield(ck, 'cfg'), error('m1Tools:ck', 'ck.cfg is required (the seed schedule).'); end
t0 = tic; n0 = st.nDone;
for r = st.nDone + 1:nTarget
    tr = tic;
    s = seedOf(ck.cfg, ctx.index, r);
    e = estimateRecord(ctx, s);
    st.seeds(r) = s; st.dtheta(:, r) = e.dtheta; st.q(r) = e.q; st.nBad(r) = e.nBad; st.seconds(r) = toc(tr);
    st.nDone = r;
    if ck.verbose && r == n0 + 1
        fprintf('      record %d: %.1f s; %d records to go, about %.0f min\n', r, st.seconds(r), nTarget - r, ...
            (nTarget - r) * st.seconds(r) / 60);
        flushOut();
    end
    if ~isempty(ck.file) && (mod(r, ck.every) == 0 || r == nTarget)
        saveCheckpoint(ck.file, st);
        if ck.verbose, fprintf('      %d / %d records (%.0f s), checkpoint saved\n', r, nTarget, toc(t0)); flushOut(); end
    end
end
end

function sm = summarise(st, tgt, cfg)
n = st.nDone;
D = st.dtheta(:, 1:n); q = st.q(1:n);
A = cfg.M1.accept; chi = [-2 * log(1 - 0.6827), -2 * log(0.05)];
sm.n = n;
sm.sigmaPred = tgt.sigma(:);
sm.sigmaEmp = std(D, 0, 2);                              % as G3a
sm.empPred = sm.sigmaEmp ./ sm.sigmaPred;
sm.mean = mean(D, 2); sm.se = sm.sigmaEmp / sqrt(n);
sm.dthetaW = tgt.dthetaW(:);
sm.offsetMatchesWelch = all(abs(sm.mean - sm.dthetaW) < 3 * sm.se);    % G3.3 rule, diagnostic only
sm.meanQ = mean(q); sm.expectedQ = 2 + tgt.dW^2;
sm.coverage = [mean(q <= chi(1)), mean(q <= chi(2))];
sm.nRecordsWithRejectedBins = nnz(st.nBad(1:n) > 0);
% acceptance (pre-registered values, cfg.M1.accept)
sm.passScatter = all(sm.empPred >= A.empPred(1) & sm.empPred <= A.empPred(2));
sm.passCoverage = sm.coverage(2) >= A.coverage95(1) && sm.coverage(2) <= A.coverage95(2);
% borderline: the 1-SE interval of a statistic contains one of its acceptance limits
[sm.borderline, sm.borderlineWhat, sm.seEmpPred, sm.seCoverage95] = isBorderline(sm.empPred, sm.coverage(2), n, cfg);
% nonlinear part (filled once the inversions exist)
sm.nNonlinear = numel(st.nl);
if isempty(st.nl)
    sm.passNonlinear = NaN;
else
    sm.passNonlinear = all([st.nl.pass]);
end
sm.pass = sm.passScatter && sm.passCoverage && isequal(sm.passNonlinear, true);
end

function [b, what, seR, seC] = isBorderline(empPred, cov95, n, cfg)
% the 1-standard-error interval of emp/pred (se = r / sqrt(2(n - 1))) or of the 95% coverage (binomial,
% p(1 - p) floored at 0.95 x 0.05) contains an acceptance limit (cfg.M1.borderline, logged 2026-10-03)
A = cfg.M1.accept; Bd = cfg.M1.borderline;
seR = empPred(:) / sqrt(2 * (n - 1));
seC = sqrt(max(cov95 * (1 - cov95), Bd.coverageVarianceFloor) / n);
nearR = any(abs(empPred(:) - A.empPred(:).') < Bd.nSE * seR, 2);
nearC = any(abs(cov95 - A.coverage95) < Bd.nSE * seC);
b = any(nearR) || nearC;
what = strjoin([repmat({'emp/pred'}, 1, any(nearR)), repmat({'coverage'}, 1, nearC)], ', ');
end

%% ================================================================================================
%  Nonlinear inversion (G3 solver)
%% ================================================================================================
function fun = makeModel(ctx)
% features on ALL bins at theta = [ln beta; ln R]; gamma, sensors (fixed physical points), node grid,
% truncation / nu (frfOptions), band, layout and reference as the target. Memoised by theta.
t = ctx.tgt; scn = ctx.scn;
memo = containers.Map('KeyType', 'char', 'ValueType', 'any');
sens = scn.sensors(scn.layout, :);
fun = @(th) modelAt(th, t, scn, sens, memo);
end

function f = modelAt(th, t, scn, sens, memo)
key = sprintf('%.15e_', th);
if isKey(memo, key), f = memo(key); return; end
p = t.p; p([1 3]) = exp(th(:).');
args = {'ThetaIdx', [], 'NodeSpacing', scn.nodeSpacing, 'FRFOptions', scn.frfOptions, 'CheckPoints', 0};
if ~isempty(scn.frf), args = [args, {'FRF', scn.frf}]; end
fld = emmField(p, sens, args{:});
[~, f] = transmissibilityFromField(fld.H(t.omega), [], 1:size(sens, 1), scn.ref);
memo(key) = f; %#ok<NASGU>   containers.Map is a handle: the entry persists
end

function nl = nonlinear(ctx, st, betaOffset, fun, cfg)
t = ctx.tgt; G = ctx.G; S = cfg.M1.solver;
r1 = cfg.M1.nonlinearRecord;
if st.nDone < r1 || st.seeds(r1) ~= seedOf(cfg, ctx.index, r1)
    error('m1Tools:nonlinearRecord', 'Record %d of %s is not in the state.', r1, ctx.name);
end
e = estimateRecord(ctx, st.seeds(r1));                       % regenerated from its seed (bitwise the same)
if any(abs(e.dtheta - st.dtheta(:, r1)) > 1e-12 * max(1, abs(st.dtheta(:, r1))))
    error('m1Tools:reproduce', 'Record %d of %s does not reproduce from its seed.', r1, ctx.name);
end
rows = reshape((find(e.use) - 1) * t.pFeat + (1:t.pFeat).', [], 1);
fr = @(th) subsetRows(fun(th), rows);
z = e.z(:, e.use);
blocks = t.blocks(:, :, e.use);
thLin = t.thetaTrue + e.dtheta;
th0 = t.thetaTrue + [log(1 + betaOffset); 0];
tt = tic;
J0 = G.jacobianCD(fr, th0, S.hStep);
g = G.gaussNewton(fr, z, blocks, th0, J0, t.sigma, S.hStep, S.maxIt, S.stepSigma, S.tolSigma);
nl = struct('betaOffset', betaOffset, 'record', r1, 'seed', st.seeds(r1), 'thetaStart', th0, 'theta', g.theta, ...
    'thetaLin', thLin, 'converged', g.converged, 'it', g.it, 'jRecomputes', g.jRecomputes, 'evals', g.evals + 5, ...
    'chi2', g.obj, 'dof', numel(z) - 2, 'chi2dof', g.obj / (numel(z) - 2), 'seconds', toc(tt), ...
    'gapSigma', (g.theta - thLin) ./ t.sigma, 'errSigma', (g.theta - t.thetaTrue) ./ t.sigma, ...
    'history', g.history, 'pass', false);
nl.pass = g.converged && all(abs(nl.gapSigma) < cfg.M1.accept.nonlinVsLinSigma);
end

function v = subsetRows(f, rows)
v = f(rows);
end

%% ================================================================================================
%  Checkpoints and provenance
%% ================================================================================================
function key = checkpointKey(ctx, cfg)
t = ctx.tgt;
key = struct('version', 'm1Tools 2026-10-03', 'name', t.name, 'index', t.index, 'p', t.p, 'sensors', t.sensors, ...
    'layout', t.layout, 'ref', t.ref, 'Snn', t.Snn, 'sea', t.seaText, 'N', t.N, 'dt', t.dt, 'L', t.L, ...
    'overlap', t.overlap, 'bandUse', t.bandUse, 'nodeSpacing', t.nodeSpacing, 'frfOptions', optText(t.frfOptions), ...
    'seedBase', cfg.M1.seedBase, 'seedStride', cfg.M1.seedStride, 'nRecords', cfg.M1.nRecords, ...
    'nRecordsExtended', cfg.M1.nRecordsExtended, 'nonlinearRecord', cfg.M1.nonlinearRecord, ...
    'nonlinearBetaOffsets', cfg.M1.nonlinearBetaOffsets, 'solver', cfg.M1.solver, 'accept', cfg.M1.accept, ...
    'borderline', cfg.M1.borderline, 'sigma', t.sigma(:).', 'dW', t.dW);
end

function [st, ok, msg] = loadCheckpoint(file, key)
st = []; ok = false; msg = 'no checkpoint';
if ~exist(file, 'file'), return; end
S = load(file, 'st');
if ~isfield(S, 'st') || ~isequal(S.st.key, key)
    msg = 'stale checkpoint (scenario, seeds, M1 settings or analytic target differ): ignored';
    return
end
st = S.st; ok = true; msg = sprintf('resumed: %d records, %d nonlinear runs', st.nDone, numel(st.nl));
end

function saveCheckpoint(file, st) %#ok<INUSD>
d = fileparts(file);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
save([file '.tmp'], 'st', '-v7');
movefile([file '.tmp'], file, 'f');               % never leave a half-written checkpoint
end

function [clean, files] = cleanTree(txt)
% uncommitted (modified, added, deleted, renamed or untracked) MATLAB files in `git status --porcelain`
lines = regexp(txt, '\r?\n', 'split');
files = {};
for i = 1:numel(lines)
    s = lines{i};
    if numel(s) < 4, continue; end
    f = strtrim(s(4:end));
    if numel(f) > 1 && f(1) == '"' && f(end) == '"', f = f(2:end - 1); end
    a = strfind(f, ' -> '); if ~isempty(a), f = f(a(end) + 4:end); end
    if numel(f) >= 2 && strcmpi(f(end - 1:end), '.m'), files{end + 1} = f; end %#ok<AGROW>
end
clean = isempty(files);
end

%% ================================================================================================
function a = seaArgs(sea)
% exactly as evaluateScenario passes the sea to synthesiseTwinRecords
if isa(sea, 'function_handle')
    a = {'Spectrum', sea};
else
    a = {'Spectrum', 'jonswap', 'Hs', sea.Hs, 'PeakFrequency', sea.omegaP, 'PeakEnhancement', sea.gammaJ};
end
end

function s = seaText(sea)
if isa(sea, 'function_handle'), s = ['handle ' func2str(sea)];
else, s = sprintf('jonswap Hs %.15g omegaP %.15g gammaJ %.15g', sea.Hs, sea.omegaP, sea.gammaJ); end
end

function s = optText(c)
s = '';
for i = 1:numel(c)
    if ischar(c{i}), s = [s c{i} ' ']; else, s = [s mat2str(c{i}, 15) ' ']; end %#ok<AGROW>
end
end

function flushOut()
if exist('OCTAVE_VERSION', 'builtin'), fflush(stdout); else, drawnow; end
end