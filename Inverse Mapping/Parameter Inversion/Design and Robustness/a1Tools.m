function A = a1Tools()
%A1TOOLS  The numerical core of study A1 (floe-property and domain sweep), shared by
%   runA1FloePropertySweep and testA1FloePropertySweep so both exercise exactly the same code.
%   A = a1Tools();   then e.g.  sc = A.buildScenarios(cfg)
%
%   A1 asks: under ONE fixed sensing arrangement (Level 2C fractions of each floe's true radius,
%   reference s26, LSM6DSV16X noise, setting F), which floes can the frozen inverse identify?
%   Every scenario gets two answers, kept apart:
%     1. Can the forward model be trusted here?  alpha inside the validated window; the 0.2 rad/s node
%        grid against 0.1 rad/s (a FEATURE error b projected through the scenario's own Fisher
%        information, d = sqrt(dth' F dth), dth = F^-1 J' Sigma^-1 b, d < 0.1); and truncation
%        [50 10 10] against [70 15 15] as METRIC convergence (sigma, d_W and class from a full
%        evaluation at the higher truncation; cfg.A1.truncationGate). The truncation feature bias d is
%        reported alongside as model-discrepancy information, not as a gate.
%     2. Only then: the inverse class from the frozen E1 chain (evaluateScenario, unchanged).
%   A failed model check is reported as 'Model check failed', never as 'Not identifiable'. If F is
%   singular the checks cannot be projected: modelAssessed = false, modelValid = false, and the final
%   class is the inverse class (Not identifiable), never a model pass that was not tested.
%
%   FREQUENCY AND SCALING RULES (cfg.A1, pre-registered 30 Sep, MIZ rule amended 2 Oct before the run)
%   * Finite depth (ring, corners, basin, tank floes): the basin itself, H = 1.88 m, the basin band
%     cfg.band.basin (alpha 1.94 to 13.2).
%   * Deep water (MIZ): screenFloeRegime gives the effective depth H_eff (min(k, kappa) H_eff = pi at
%     12.5 s) and the EMM parameters at H_eff. The EMM is non-dimensional (it sees only alpha, beta,
%     gamma, R), so the case is evaluated FROUDE-SCALED to the basin: lengths x 1.88 / H_eff, times x
%     sqrt(1.88 / H_eff). The 6 to 12.5 s band maps to alpha = H_eff omega^2 / g and then to basin
%     omega = sqrt(alpha g / 1.88), intersected with the basin band (which excludes the synthesis
%     tapers), so periods slightly shorter than 6 s may be cut; the covered periods are reported.
%     Sea state, record length and noise are the basin ones in scaled units: the result is the
%     identifiability of the non-dimensional MIZ case under the basin-equivalent experiment. The
%     equivalent full-scale experiment (Hs, Tp, duration, sample interval, noise PSD) is reported.
%   * Poisson ratio: passed to the EMM ('Nu') whenever it differs from 0.3 (polypropylene 0.40). The
%     free-edge conditions depend on nu, not only on D. nu = 0.3 passes nothing, so p0 shares the
%     cache with every earlier study.
%
%   FUNCTIONS
%   sc   = A.buildScenarios(cfg)        24 scenarios (8 ring, 4 corners, 12 physical), all converted
%   s    = A.physicalToEMM(phys, cfg)   one physical class -> p, band, H_eff, scale (used by the above)
%   band = A.mizBand(H_eff, cfg)        basin-equivalent band [rad/s] and the periods it covers
%   pos  = A.placeSensors(sc, cfg)      Level 2C fractions x the TRUE R, [r theta] non-dimensional
%   scn  = A.buildScn(sc, cfg, opts)    the evaluateScenario input (sensors fixed at the true R)
%   g    = A.dimensionlessGroups(sc, omega, cfg)   alpha, kR (min / peak / max), R / l_f
%   res  = A.evaluate(sc, cfg, opts)    production metrics + both model checks + groups (one scenario)
%   q    = A.projectBias(out, b)        d and dtheta of a feature error b with out's own F
%   f    = A.features(sc, scn, omega, nodeSpacing, frfOptions, opts)   stacked T at the bins
%   key  = A.checkpointKey(sc, cfg)     everything a checkpoint's validity depends on
%   file = A.checkpointPath(dir, sc)
%   [res, ok] = A.loadCheckpoint(file, key)   ok = false if missing or stale (key mismatch)
%   A.saveCheckpoint(file, res, key)
%   c    = A.classCode(name)            1 Identifiable, 2 Marginal, 3 Not identifiable, 4 Model check failed
%
%   opts (A.evaluate / A.buildScn): struct with cacheDir (''), verbose (false), nodeSpacing
%   (cfg.emm.nodeSpacing), runChecks (true), frf ([] = EMM; an analytic handle for tests).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

% screenFloeRegime (MIZ effective depth) lives in Mission Physics/, which setupStudyPaths does not add
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(root, 'Mission Physics'));

A.buildScenarios = @buildScenarios;
A.physicalToEMM = @physicalToEMM;
A.mizBand = @mizBand;
A.placeSensors = @placeSensors;
A.buildScn = @buildScn;
A.dimensionlessGroups = @dimensionlessGroups;
A.evaluate = @evaluateA1;
A.projectBias = @projectBias;
A.features = @features;
A.checkpointKey = @checkpointKey;
A.checkpointPath = @checkpointPath;
A.loadCheckpoint = @loadCheckpoint;
A.saveCheckpoint = @saveCheckpoint;
A.classCode = @classCode;
A.version = 'a1Tools 2026-10-02b';
end

%% ================================================================================================
%  Scenarios
%% ================================================================================================
function sc = buildScenarios(cfg)
p0 = cfg.p0.vec; Hb = cfg.physics.Hbasin; g = cfg.physics.g; rhoW = cfg.physics.rhoWbasin;
sc = [];
% (a) ring and corners: basin, p0 gamma, nu = 0.3
for i = 1:numel(cfg.A1.param)
    q = cfg.A1.param(i);
    s = blank();
    s.name = q.name; s.group = q.group; s.p = q.p;
    s.betaRatio = q.p(1) / p0(1); s.RRatio = q.p(3) / p0(3); s.gammaRatio = q.p(2) / p0(2);
    s.nu = cfg.physics.nuEMM; s.deepWater = false; s.Hphys = Hb; s.rhoW = rhoW;
    s.Rphys = q.p(3) * Hb; s.D = q.p(1) * rhoW * g * Hb^4; s.m = q.p(2) * Hb * rhoW;
    s.lf = (s.D / (rhoW * g))^(1/4);
    s.band = cfg.A1.bandRule.finiteDepthBand; s.scale = 1;
    s.note = 'parameter-space case (basin, nu = 0.3)';
    sc = [sc, s]; %#ok<AGROW>
end
% (b) physical classes
for i = 1:numel(cfg.A1.physical)
    sc = [sc, physicalToEMM(cfg.A1.physical(i), cfg)]; %#ok<AGROW>
end
for i = 1:numel(sc)
    sc(i).index = i;
    sc(i).frfOptions = {};
    if abs(sc(i).nu - 0.3) > 1e-12, sc(i).frfOptions = {'Nu', sc(i).nu}; end
    sc(i).bandPeriods = sort(2 * pi ./ (sc(i).band / sqrt(sc(i).scale)), 'ascend');   % full-scale periods, s
end
end

function s = physicalToEMM(phys, cfg)
p0 = cfg.p0.vec; Hb = cfg.physics.Hbasin; g = cfg.physics.g;
s = blank();
s.name = phys.name; s.group = phys.class; s.provisional = phys.provisional; s.note = phys.note;
s.deepWater = phys.deepWater; s.Rphys = phys.Rphys; s.h = phys.h; s.E = phys.E; s.nu = phys.nu;
s.rhoI = phys.rhoI; s.rhoW = phys.rhoW; s.D = phys.D; s.m = phys.m; s.lf = phys.lf;
if ~phys.deepWater
    s.p = phys.p; s.Hphys = phys.H; s.scale = 1;
    s.band = cfg.A1.bandRule.finiteDepthBand;
    if s.m / s.rhoW >= s.Hphys, error('a1Tools:draught', '%s: draught exceeds the depth.', s.name); end
else
    floe = struct('R', phys.Rphys, 'h', phys.h, 'E', phys.E, 'nu', phys.nu, 'rhoIce', phys.rhoI);
    scr = screenFloeRegime(floe, cfg.A1.bandRule.mizPeriods, 'Depth', Inf, 'RhoWater', phys.rhoW, ...
        'Gravity', g, 'AlphaRange', cfg.emm.alphaValidated);
    s.Hphys = scr.emm.H;                                  % effective depth, m
    s.p = [scr.emm.beta, scr.emm.gamma, scr.emm.Rnd];
    s.scale = s.Hphys / Hb;                               % full-scale length / basin length
    [s.band, s.alphaFullScale] = mizBand(s.Hphys, cfg);
    s.screen = struct('note', scr.emm.note, 'alpha', scr.emm.alpha, 'T', scr.T, 'kR', scr.kR, ...
        'kappaRatio', scr.kappaRatio, 'RoverLf', scr.RoverLf);
end
s.betaRatio = s.p(1) / p0(1); s.RRatio = s.p(3) / p0(3); s.gammaRatio = s.p(2) / p0(2);
end

function [band, alphaFull] = mizBand(Heff, cfg)
% the full-scale period band, as alpha at the effective depth, as basin rad/s, cut to the basin band
g = cfg.physics.g; Hb = cfg.physics.Hbasin;
T = sort(cfg.A1.bandRule.mizPeriods, 'descend');          % [12.5 6] s -> increasing omega
alphaFull = Heff * (2 * pi ./ T).^2 / g;
wb = sqrt(alphaFull * g / Hb);
band = [max(wb(1), cfg.band.basin(1)), min(wb(2), cfg.band.basin(2))];
if band(1) >= band(2), error('a1Tools:mizBand', 'The mapped MIZ band does not overlap the basin band.'); end
end

function s = blank()
s = struct('index', 0, 'name', '', 'group', '', 'provisional', false, 'note', '', 'p', [NaN NaN NaN], ...
    'betaRatio', NaN, 'RRatio', NaN, 'gammaRatio', NaN, 'deepWater', false, 'Hphys', NaN, 'scale', 1, ...
    'Rphys', NaN, 'h', NaN, 'E', NaN, 'nu', 0.3, 'rhoI', NaN, 'rhoW', NaN, 'D', NaN, 'm', NaN, 'lf', NaN, ...
    'band', [NaN NaN], 'bandPeriods', [NaN NaN], 'alphaFullScale', [NaN NaN], 'screen', [], 'frfOptions', {{}});
end

%% ================================================================================================
%  Sensors, scenario struct, groups
%% ================================================================================================
function pos = placeSensors(sc, cfg)
% Level 2C fractions of the TRUE radius, converted ONCE; emmField then keeps these points fixed while
% R is perturbed (dH/dlnR is the response change at the same physical points)
pos = cfg.A1.layoutFrac .* [sc.p(3) 1];
end

function scn = buildScn(sc, cfg, opts)
opts = fillOpts(opts, cfg);
scn = struct('p', sc.p, 'sensors', placeSensors(sc, cfg), 'layout', 1:size(cfg.A1.layoutFrac, 1), ...
    'ref', cfg.layout.level2C.ref, 'Snn', cfg.noise.(cfg.A1.noiseCase), 'sea', cfg.sea, ...
    'bandUse', sc.band, 'nodeSpacing', opts.nodeSpacing, 'frfOptions', {sc.frfOptions}, ...
    'cacheDir', opts.cacheDir, 'verbose', opts.verbose, 'frf', opts.frf);
end

function o = fillOpts(o, cfg)
if nargin < 1 || isempty(o), o = struct(); end
d = struct('cacheDir', '', 'verbose', false, 'nodeSpacing', cfg.emm.nodeSpacing, 'runChecks', true, 'frf', []);
fn = fieldnames(d);
for i = 1:numel(fn), if ~isfield(o, fn{i}), o.(fn{i}) = d.(fn{i}); end, end
end

function G = dimensionlessGroups(sc, omega, cfg)
% everything in the EMM scaling (lengths / H): alpha = H omega^2 / g at the BASIN depth, which equals
% the full-scale alpha for a Froude-scaled MIZ case; kH from alpha = kH tanh(kH); kR = kH R;
% l_f / H = beta^(1/4), so R / l_f = R beta^(-1/4)
Hb = cfg.physics.Hbasin; g = cfg.physics.g;
alpha = Hb * omega.^2 / g;
kH = arrayfun(@openWaterKH, alpha);
alphaPk = Hb * cfg.sea.omegaP^2 / g;
G.alphaMin = min(alpha); G.alphaMax = max(alpha);
G.alphaValid = G.alphaMin >= cfg.emm.alphaValidated(1) && G.alphaMax <= cfg.emm.alphaValidated(2);
G.kR = kH * sc.p(3);
G.kRmin = min(G.kR); G.kRmax = max(G.kR);
G.kRpeak = openWaterKH(alphaPk) * sc.p(3);
G.RoverLf = sc.p(3) * sc.p(1)^(-1/4);
G.lfPhys = sc.p(1)^(1/4) * sc.Hphys;                     % m, full scale
end

function x = openWaterKH(alpha)
% x tanh x = alpha (finite-depth dispersion in the EMM scaling)
x = max(alpha, sqrt(alpha));
for it = 1:60
    t = tanh(x); f = x * t - alpha; fp = t + x * (1 - t^2);
    dx = f / fp; x = x - dx;
    if abs(dx) < 1e-14 * x, break; end
end
end

%% ================================================================================================
%  One scenario: production metrics, model checks, groups
%% ================================================================================================
function res = evaluateA1(sc, cfg, opts)
opts = fillOpts(opts, cfg);
t0 = tic;
res = struct('name', sc.name, 'scenario', sc, 'status', 'ok', 'message', '');

% ---- production (frozen E1 chain, unchanged) ----
scn = buildScn(sc, cfg, opts);
out = evaluateScenario(scn, cfg);
gridUsed = opts.nodeSpacing;

% ---- check 1: node grid (production spacing vs the record grid, 0.1 rad/s) ----
% a check counts as passed only if it was run and passed, or if cfg deliberately does not require it
chkGrid = struct('done', false, 'd', NaN, 'relT', NaN, 'valid', ~(opts.runChecks && cfg.A1.requireNodeGridCheck), ...
    'fallback', false, 'spacing', gridUsed);
if opts.runChecks && cfg.A1.requireNodeGridCheck && ~out.singular
    fFine = features(sc, scn, out.omega, cfg.emm.nodeSpacingCheck, sc.frfOptions, opts);
    [chkGrid.d, chkGrid.relT] = checkMetrics(out, fFine);
    chkGrid.done = true;
    chkGrid.valid = chkGrid.d < cfg.emm.maxInterpBiasDistance;
    if ~chkGrid.valid
        % pre-registered fallback: evaluate on the 0.1 rad/s grid itself (no finer grid to check it against)
        if opts.verbose, fprintf('      node-grid d = %.3g >= %.2g: falling back to %.2g rad/s\n', chkGrid.d, ...
                cfg.emm.maxInterpBiasDistance, cfg.emm.nodeSpacingCheck); end
        scn.nodeSpacing = cfg.emm.nodeSpacingCheck;
        out = evaluateScenario(scn, cfg);
        gridUsed = scn.nodeSpacing;
        chkGrid.fallback = true; chkGrid.spacing = gridUsed;
        chkGrid.valid = true;                    % resolved by the fallback rule (unverified beyond 0.1)
    end
end

% ---- check 2: truncation (nominal vs higher, at the grid actually used) ----
% GATE (amended 2026-10-02 after the p0 diagnostic): do the A1 METRICS converge? A full evaluation at the
% check truncation must give sigma within cfg.A1.truncationGate.maxRelSigma, d_W within maxDeltaDW and
% the same inverse class. The feature-level bias d of the production model against the check model is
% still computed and reported (model-discrepancy information for V1), but it is not a gate: it measures
% how wrong an inversion of real data with this model would be, not whether A1's metrics are reliable.
chkTr = struct('done', false, 'd', NaN, 'relT', NaN, 'relSigma', NaN, 'deltaDW', NaN, 'classHi', '', ...
    'sameClass', false, 'valid', ~(opts.runChecks && cfg.A1.requireTruncationCheck), ...
    'truncation', cfg.emm.truncationCheck);
if opts.runChecks && cfg.A1.requireTruncationCheck && ~out.singular
    scnHi = scn;
    scnHi.frfOptions = [sc.frfOptions, {'Truncation', cfg.emm.truncationCheck}];
    outHi = evaluateScenario(scnHi, cfg);
    [chkTr.d, chkTr.relT] = checkMetrics(out, outHi.f0);
    chkTr.classHi = outHi.class;
    chkTr.sameClass = strcmp(outHi.class, out.class);
    if outHi.singular
        chkTr.relSigma = Inf; chkTr.deltaDW = Inf;
    else
        chkTr.relSigma = max(abs(outHi.sigma(:) ./ out.sigma(:) - 1));
        chkTr.deltaDW = abs(outHi.dW - out.dW);
    end
    chkTr.done = true;
    G = cfg.A1.truncationGate;
    chkTr.valid = chkTr.relSigma < G.maxRelSigma && chkTr.deltaDW < G.maxDeltaDW && chkTr.sameClass;
end

% ---- groups and alpha ----
G = dimensionlessGroups(sc, out.omega, cfg);

% ---- classification: model validity first ----
% A singular F cannot project the grid or truncation discrepancy into parameter space, so the model
% checks are NOT ASSESSABLE (modelAssessed = false, modelValid = false). The inverse conclusion still
% stands: F singular means Not identifiable, and that is the final class. It is never recorded as a
% scenario whose forward model passed its checks.
inverseClass = out.class;
why = {};
if ~G.alphaValid, why{end + 1} = sprintf('alpha [%.2f %.2f] outside validated', G.alphaMin, G.alphaMax); end
if out.singular
    modelAssessed = false;
    modelValid = false;
    finalClass = inverseClass;
    why{end + 1} = 'model checks not assessable: F singular';
else
    modelAssessed = true;
    modelValid = G.alphaValid && chkGrid.valid && chkTr.valid;
    if modelValid, finalClass = inverseClass; else, finalClass = cfg.A1.modelCheckClass; end
    if ~chkGrid.valid, why{end + 1} = sprintf('node grid d = %.3g', chkGrid.d); end
    if ~chkTr.valid, why{end + 1} = sprintf('truncation: sigma %.2g rel, d_W %.2g abs, class %s', chkTr.relSigma, ...
            chkTr.deltaDW, chkTr.classHi); end
end

% ---- results (no function handles: checkpoints stay small and portable) ----
res.sigma = out.sigma(:).'; res.sigmaBeta = out.sigma(1); res.sigmaR = out.sigma(2);
res.rho = out.corr; res.condF = out.condF; res.singular = out.singular; res.F = out.F;
res.dthetaW = out.dthetaW(:).'; res.dW = out.dW; res.dTot = out.dTot; res.sMax = out.sMax;
res.inverseClass = inverseClass; res.modelValid = modelValid; res.finalClass = finalClass;
res.finalCode = classCode(finalClass); res.inverseCode = classCode(inverseClass);
res.modelReason = strjoin(why, '; '); res.modelAssessed = modelAssessed;
res.alphaValid = G.alphaValid; res.nodeGridValid = chkGrid.valid; res.truncationValid = chkTr.valid;
res.nodeGrid = chkGrid; res.truncation = chkTr; res.nodeSpacingUsed = gridUsed;
res.groups = G;
res.omega = out.omega; res.nEff = out.nEff; res.perBinInfo = out.perBinInfo;
res.minCoherence = out.minCoherence; res.invalidFraction = out.invalidFraction; res.nUsedBins = out.nUsedBins;
res.sensors = scn.sensors; res.Snn = scn.Snn; res.frfOptions = sc.frfOptions;
res.equivalent = equivalentExperiment(sc, scn, cfg);
res.seconds = toc(t0);
end

function [d, relT] = checkMetrics(out, fRef)
% feature error of the production model against a better-resolved one, projected with the
% production F: b = (reference features) - (production features), data minus model
b = fRef - out.f0;
q = projectBias(out, b);
d = q.d;
v = out.validMask;
m = size(out.blocksAll, 1) / 2;                % transmissibilities in the layout
Tp = unstackTransmissibility(out.f0, m); Tr = unstackTransmissibility(fRef, m);
relT = max(max(abs(Tr(:, v) - Tp(:, v)) ./ abs(Tr(:, v))));
end

function q = projectBias(out, b)
% the E1 engine's projector on the bins the estimator accepts (identical to evaluateScenario's)
if out.singular, q = struct('d', NaN, 'dtheta', [NaN; NaN]); return; end
v = out.validMask;
pBlk = size(out.blocksAll, 1);
rows = reshape((find(v) - 1) * pBlk + (1:pBlk).', [], 1);
r = fisherFromBlocks(out.Jall(rows, :), out.blocksAll(:, :, v), b(rows));
q = struct('d', r.d, 'dtheta', r.dtheta);
end

function f = features(sc, scn, omega, nodeSpacing, frfOptions, opts) %#ok<INUSL>
% stacked transmissibility features of the layout at the bins, from the EMM at p only (no derivatives)
if nargin < 6 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'cacheDir'), opts.cacheDir = ''; end          % callers may pass a partial opts struct
if ~isfield(opts, 'verbose'), opts.verbose = false; end
if ~isfield(opts, 'frf'), opts.frf = []; end
fld = emmField(scn.p, scn.sensors, 'ThetaIdx', [], 'NodeSpacing', nodeSpacing, 'FRFOptions', frfOptions, ...
    'CacheDir', opts.cacheDir, 'Verbose', opts.verbose, 'FRF', opts.frf);
[~, f] = transmissibilityFromField(fld.H(omega), [], scn.layout, scn.ref);
end

function e = equivalentExperiment(sc, scn, cfg)
% the full-scale experiment the scaled evaluation represents (identity for the basin and tank cases)
s = sc.scale; ts = sqrt(s);
e = struct('scale', s, 'timeScale', ts, 'Hs', cfg.sea.Hs * s, 'Tp', 2 * pi / cfg.sea.omegaP * ts, ...
    'recordMin', cfg.acq.Nref * cfg.acq.dt * ts / 60, 'dt', cfg.acq.dt * ts, 'Snn', scn.Snn * ts, ...
    'periods', sc.bandPeriods, 'note', '');
if sc.deepWater
    e.note = sprintf(['Froude-scaled: lengths x %.1f, times x %.2f; noise PSD per rad/s x %.2f (a real ' ...
        'LSM6DSV16X at full scale is %.1f times quieter than this case assumes)'], s, ts, ts, ts);
end
end

%% ================================================================================================
%  Checkpoints and class codes
%% ================================================================================================
function key = checkpointKey(sc, cfg)
key = struct('version', 'a1Tools 2026-10-02b', 'name', sc.name, 'p', sc.p, 'band', sc.band, ...
    'frfOptions', {sc.frfOptions}, 'sensors', placeSensors(sc, cfg), 'Snn', cfg.noise.(cfg.A1.noiseCase), ...
    'sea', cfg.sea, 'N', cfg.acq.Nref, 'L', cfg.welch.L, 'nodeSpacing', cfg.emm.nodeSpacing, ...
    'nodeSpacingCheck', cfg.emm.nodeSpacingCheck, 'truncation', cfg.emm.truncation, ...
    'truncationCheck', cfg.emm.truncationCheck, 'maxGrid', cfg.emm.maxInterpBiasDistance, ...
    'truncGate', cfg.A1.truncationGate, 'class', cfg.class);
end

function file = checkpointPath(dirName, sc)
file = fullfile(dirName, sprintf('A1_%02d_%s.mat', sc.index, sc.name));
end

function [res, ok] = loadCheckpoint(file, key)
res = []; ok = false;
if ~exist(file, 'file'), return; end
S = load(file);
if isfield(S, 'key') && isequaln(S.key, key) && isfield(S, 'res')
    res = S.res; ok = true;
end
end

function saveCheckpoint(file, res, key) %#ok<INUSL>
d = fileparts(file);
if ~exist(d, 'dir'), mkdir(d); end
saved = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<NASGU>
save(file, 'res', 'key', 'saved', '-v7');
end

function c = classCode(name)
switch name
    case 'Identifiable', c = 1;
    case 'Marginal', c = 2;
    case 'Not identifiable', c = 3;
    case 'Model check failed', c = 4;
    otherwise, c = NaN;
end
end