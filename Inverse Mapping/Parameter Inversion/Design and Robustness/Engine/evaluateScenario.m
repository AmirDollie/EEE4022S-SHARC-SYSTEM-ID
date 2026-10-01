function out = evaluateScenario(scn, cfg)
%EVALUATESCENARIO  E1: one scenario -> inverse-quality metrics, analytically (no Monte Carlo).
%   out = evaluateScenario(scn)
%   out = evaluateScenario(scn, cfg)        cfg = defineStudyScenarios() supplies defaults and
%                                           the classification thresholds
%
%   CHAIN   scenario -> H, dH (emmField) -> bins (welchBins) -> T, J (transmissibilityFromField)
%           -> Sigma_k (analyticTransmissibilityCovariance on the Welch-smoothed spectral matrix)
%           -> keep only the bins the estimator accepts -> F, sigma, rho, kappa (fisherFromBlocks)
%           -> Welch bias and any systematic bias projected as vectors -> d_tot -> class
%
%   ESTIMATED vs NUISANCE parameters. The inverse estimates theta = [ln beta, ln R] with gamma
%   KNOWN (frozen pipeline). 'estIdx' (default [1 3]) are the ONLY columns that enter F.
%   'nuisanceIdx' (e.g. 2 for gamma, R4) are differentiated too, but returned separately as
%   out.Jnuisance and never estimated. Sign convention for systematic errors: bias = data features
%   minus model features at the truth. If the inverse fixes gamma_inv = gamma_true (1 + eps), the
%   model is evaluated at the wrong gamma, so the feature error is
%       bSys = f(gamma_true) - f(gamma_inv) ~ -Jnuisance(:, k) * log(1 + eps)
%   (pass it as scn.sysBias, or as a handle that builds it).
%
%   REJECTED BINS. Bins failing the estimator's denominator rule (MinDenominatorFraction) carry no
%   information in the real inverse, so they are removed from F and from every bias projection.
%   out.invalidFraction reports how many were removed; out.Jall / blocksAll / validMask keep the
%   full set for diagnostics.
%
%   SCN fields (defaults from cfg in brackets)
%     p          [beta gamma R]                                  [cfg.p0.vec]
%     sensors    nP x 2 [r theta], NON-DIMENSIONAL, fixed points [Level 2C at p(3)]
%     layout     indices into sensors used                         [1:nP]
%     ref        reference, index INTO layout                      [cfg.layout.level2C.ref]
%     estIdx     estimated ln-parameters (enter F)                   [1 3]
%     nuisanceIdx  differentiated but NOT estimated (e.g. 2 = gamma)  []
%     Snn        scalar or per-layout-sensor noise PSD (per rad/s)  [cfg.noise.(default)]
%     sea        struct Hs, omegaP, gammaJ, OR spectrum handle @(w) [cfg.sea]
%     N, dt, L, overlap, binStep, bandEst, bandUse                   [cfg.acq / cfg.welch / cfg.band]
%     nodeSpacing, frfOptions, cacheDir                             [cfg.emm.nodeSpacing, {}, '']
%     frf        [] (EMM) or handle @(omega, p) analytic FRF (tests)  []
%     verbose    print one progress line per node FRF                false
%     sysBias    OPTIONAL: deliberate systematic feature error, a feature vector on the bins, or a
%                handle @(T, omega, ctx) -> feature vector (ctx carries H, Ssmoothed, Snn, cv ...)
%     computeWelchBias  true
%
%   OUT: sigma (ln beta, ln R), corr, condF, singular, F, covTheta, dthetaW, dW, dthetaSys, dSys,
%        dthetaTot, dTot, sMax, class, minCoherence (accepted bins), invalidFraction, nUsedBins,
%        degenerateBins, omega, nEff, K, nDof, validMask, blocksAll, Jall, Jnuisance, f0, bW,
%        bSys, perBinInfo, field (emmField struct), scn (resolved), seconds
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Engine/.

if nargin < 2 || isempty(cfg), cfg = defineStudyScenarios(); end
t0 = tic;
scn = resolve(scn, cfg);

% ---- field and bins ----
nEst = numel(scn.estIdx);
fld = emmField(scn.p, scn.sensors, 'ThetaIdx', [scn.estIdx(:).', scn.nuisanceIdx(:).'], 'NodeSpacing', ...
    scn.nodeSpacing, 'FRFOptions', scn.frfOptions, 'CacheDir', scn.cacheDir, 'Step', scn.step, 'FRF', scn.frf, ...
    'Verbose', scn.verbose);
[omega, kBins] = welchBins(scn.L, scn.dt, scn.bandEst, scn.binStep, scn.bandUse);
H = fld.H(omega); dH = fld.dH(omega);
[T0, f0, Jboth] = transmissibilityFromField(H, dH, scn.layout, scn.ref);
J = Jboth(:, 1:nEst);                        % estimated parameters only
Jnuis = Jboth(:, nEst + 1:end);              % nuisance sensitivities (never enter F)

% ---- record statistics: spectrum, taper, Welch kernel ----
specArgs = seaArgs(scn.sea);
[~, ~, info] = synthesiseTwinRecords(fld.twin, scn.N, scn.dt, 'Seed', 1, specArgs{:});
[nEff, K] = welchNEff(scn.N, scn.L, scn.overlap);
Ssm = welchSmoothedSpectrum(omega, info, scn.L, scn.layout);
Snn = scn.Snn; if isscalar(Snn), Snn = Snn * ones(numel(scn.layout), 1); end
cv = analyticTransmissibilityCovariance(Ssm, [], Snn, scn.ref, nEff);

% ---- biases ----
others = setdiff(1:numel(scn.layout), scn.ref);
bias = zeros(numel(f0), 0); bW = zeros(numel(f0), 1);
if scn.computeWelchBias
    infoL = info; infoL.H = info.H(scn.layout, :);
    TW = welchExpectedTransmissibility(omega, infoL, scn.L, scn.ref, others, 0);
    bW = stackTransmissibility(TW) - f0;
end
bSys = zeros(numel(f0), 1); hasSys = false;
if ~isempty(scn.sysBias)
    hasSys = true;
    if isa(scn.sysBias, 'function_handle')
        ctx = struct('H', H(scn.layout, :), 'Ssmoothed', Ssm, 'Snn', Snn, 'cv', cv, 'info', info, 'scn', scn);
        bSys = scn.sysBias(T0, omega, ctx);
    else
        bSys = scn.sysBias(:);
    end
end
bias = [bW, bSys, bW + bSys];

% ---- keep only the bins the estimator accepts ----
pBlk = size(cv.blocks, 1);
rows = reshape((find(cv.valid) - 1) * pBlk + (1:pBlk).', [], 1);
q = fisherFromBlocks(J(rows, :), cv.blocks(:, :, cv.valid), bias(rows, :));

% ---- metrics ----
out = struct();
out.sigma = q.sigma; out.corr = q.corr; out.condF = q.cond; out.singular = q.singular;
out.F = q.F; out.covTheta = q.covTheta;
out.dthetaW = q.dtheta(:, 1); out.dW = q.d(1);
out.dthetaSys = q.dtheta(:, 2); out.dSys = q.d(2);
if ~hasSys, out.dSys = 0; end
out.dthetaTot = q.dtheta(:, 3); out.dTot = q.d(3);
out.sMax = max(q.sigma(1:min(2, end)));
out.class = classify(out, cfg);
out.minCoherence = min(min(cv.coherence(:, cv.valid)));
if ~any(cv.valid), out.minCoherence = NaN; end
out.invalidFraction = mean(~cv.valid);
out.nUsedBins = q.nUsed; out.degenerateBins = q.degenerateBins;
out.omega = omega; out.kBins = kBins; out.nEff = nEff; out.K = K; out.nDof = cv.nDof;
out.validMask = cv.valid; out.blocksAll = cv.blocks; out.Jall = J; out.Jnuisance = Jnuis;
out.f0 = f0; out.T0 = T0; out.bW = bW; out.bSys = bSys;
out.perBinInfo = q.perBin; out.coherence = cv.coherence;
out.field = fld; out.scn = scn; out.seconds = toc(t0);
end

%% ================================================================================================
function scn = resolve(scn, cfg)
if isfield(scn, 'thetaIdx')
    error('evaluateScenario:thetaIdx', ['scn.thetaIdx was replaced by estIdx (estimated, enter F) and ' ...
        'nuisanceIdx (differentiated only). gamma must never be estimated.']);
end
d = struct('p', cfg.p0.vec, 'sensors', [], 'layout', [], 'ref', cfg.layout.level2C.ref, 'estIdx', [1 3], ...
    'nuisanceIdx', [], 'frf', [], 'verbose', false, ...
    'Snn', cfg.noise.(cfg.noise.default), 'sea', cfg.sea, 'N', cfg.acq.Nref, 'dt', cfg.acq.dt, ...
    'L', cfg.welch.L, 'overlap', cfg.welch.overlap, 'binStep', cfg.welch.binStep, 'bandEst', cfg.band.estimate, ...
    'bandUse', cfg.band.basin, 'nodeSpacing', cfg.emm.nodeSpacing, 'frfOptions', {{}}, 'cacheDir', '', ...
    'step', cfg.emm.jacobianStep, 'sysBias', [], 'computeWelchBias', true);
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(scn, fn{i}) || (isempty(scn.(fn{i})) && ~any(strcmp(fn{i}, {'sysBias', 'frfOptions', 'cacheDir', ...
            'nuisanceIdx', 'frf'})))
        scn.(fn{i}) = d.(fn{i});
    end
end
scn.p = scn.p(:).';
scn.estIdx = scn.estIdx(:).'; scn.nuisanceIdx = scn.nuisanceIdx(:).';
if ~isequal(scn.estIdx, [1 3])               % order matters: sigma(1) = ln beta, sigma(2) = ln R
    error('evaluateScenario:estIdx', 'The frozen inverse estimates only [ln beta, ln R]; estIdx must be exactly [1 3].');
end
if ~(isempty(scn.nuisanceIdx) || isequal(scn.nuisanceIdx, 2))
    error('evaluateScenario:nuisanceIdx', 'nuisanceIdx must be [] or 2 (gamma) and must not overlap estIdx.');
end
if isempty(scn.sensors), scn.sensors = cfg.layout.level2C.frac .* [scn.p(3) 1]; end
if isempty(scn.layout), scn.layout = 1:size(scn.sensors, 1); end
end

function a = seaArgs(sea)
if isa(sea, 'function_handle')
    a = {'Spectrum', sea};
else
    a = {'Spectrum', 'jonswap', 'Hs', sea.Hs, 'PeakFrequency', sea.omegaP, 'PeakEnhancement', sea.gammaJ};
end
end

function c = classify(m, cfg)
I = cfg.class.identifiable; M = cfg.class.marginal;
if m.singular || ~all(isfinite(m.sigma))
    c = cfg.class.names{3};
    return
end
d = m.dTot;
if m.sMax <= I.maxSigma && d <= I.maxBiasDistance && m.condF <= I.maxCondF
    c = cfg.class.names{1};
elseif m.sMax <= M.maxSigma && d <= M.maxBiasDistance && m.condF <= M.maxCondF
    c = cfg.class.names{2};
else
    c = cfg.class.names{3};
end
end