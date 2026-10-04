function T = v1Tools()
%V1TOOLS  The numerical core of study V1 (forward-model mismatch), shared by runV1ModelMismatch and
%   testV1ModelMismatch so the tester exercises exactly the code the run uses.
%   T = v1Tools();   then e.g.  m = T.truthModel(base, 'depth', 0.02, cfg)
%
%   QUESTION. M1 showed the estimator behaves as predicted when data and inverse share the EMM. V1 breaks that:
%   the data come from a perturbed TRUTH model, the inverse keeps the frozen nominal model, and the result is
%   the parameter bias the mismatch causes.
%
%   DEFINITION (Welch expectation minus Welch expectation, so the known Welch bias is not counted twice)
%       b_sys   = f_W(truth) - f_W(nominal)                 feature discrepancy on the retained bins
%       dth_sys = F^-1 J' Sigma^-1 b_sys,  d_sys = sqrt(dth_sys' F dth_sys)     with the NOMINAL J, Sigma, F
%       dth_tot = dth_W + dth_sys (bias vectors added, the programme rule), d_tot, and the class from d_tot
%   f_W is the Welch-expected transmissibility of the frozen estimator (welchExpectedTransmissibility on the
%   record statistics of synthesiseTwinRecords), computed by ONE function for truth and nominal alike.
%
%   TRUTH MODELS (base = the nominal scenario; the inverse is never changed)
%     'none'        the nominal model itself (zero discrepancy by construction)
%     'truncation'  value [M P N]: the EMM at a higher truncation (frfOptions 'Truncation')
%     'depth'       value eps: true depth H0 (1 + eps), PHYSICALLY consistent: the floe (D, m, R_phys), the
%                   sensors' physical positions and the physical-frequency sea are fixed, so in the EMM scaling
%                   beta ~ H^-4, gamma ~ H^-1, R ~ H^-1, r_i ~ H^-1, alpha = H omega^2 / g (frfOptions
%                   'WaterDepth'). The synthesis band and taper are held at the NOMINAL physical band, so the
%                   excitation is identical; at +5% the EMM is then evaluated slightly outside its validated
%                   alpha window (reported by T.alphaRange, logged as a caveat).
%     'radial'      value dr (m): every sensor with r > 0 moved radially by dr (common signed error), angles
%                   unchanged; the inverse keeps the nominal positions
%   Every truth model has the nominal physical band (the validated band at H0; the depth model gets it explicitly),
%   so the node grid and taper are the same.
%
%   FUNCTIONS
%   base = T.base(out)                          nominal p, sensors, frfOptions and physical band of a scenario
%   m    = T.truthModel(base, kind, value, cfg) truth model struct: p, sensors, frfOptions, band, H, kind, value
%   [fW, f, info] = T.welchFeatures(m, out, cacheDir)   Welch-expected (and exact) features on out's bins
%   q    = T.project(out, b)                    d and dtheta of a feature error b (the engine projector, accepted bins)
%   r    = T.evaluate(base, out, fW0, kind, value, cfg, cacheDir)   one mismatch case: d_sys, dbeta/beta, dR/R,
%                                               d_tot, class, alpha range
%   c    = T.classify(out, dTot, cfg)           the pre-registered class with the total bias
%   a    = T.alphaRange(H, omegaRetained, band, cfg)   alpha min / max (retained bins and synthesis support)
%   nl   = T.nonlinear(out, zW, thetaLin, cfg)  G3 solver on the noiseless Welch features zW from theta_true
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
T.base = @base;
T.truthModel = @truthModel;
T.welchFeatures = @welchFeatures;
T.project = @project;
T.evaluate = @evaluateCase;
T.classify = @classify;
T.alphaRange = @alphaRange;
T.nonlinear = @nonlinear;
T.version = 'v1Tools 2026-10-04';
end

%% ================================================================================================
function b = base(out)
scn = out.scn;
b = struct('p', scn.p, 'sensors', scn.sensors, 'frfOptions', {scn.frfOptions}, 'band', out.field.band, ...
    'nodeSpacing', scn.nodeSpacing, 'H', depthOf(scn.frfOptions));
end

function m = truthModel(b, kind, value, cfg)
m = struct('p', b.p, 'sensors', b.sensors, 'frfOptions', {b.frfOptions}, 'band', [], ...
    'nodeSpacing', b.nodeSpacing, 'H', b.H, 'kind', kind, 'value', value);
switch kind
    case 'none'
    case 'truncation'
        m.frfOptions = [stripOpt(b.frfOptions, 'Truncation'), {'Truncation', value}];
    case 'depth'
        s = 1 / (1 + value);                                   % H0 / H_true
        m.H = b.H * (1 + value);
        m.p = b.p .* [s^4, s, s];                              % beta ~ H^-4, gamma ~ H^-1, R ~ H^-1
        m.sensors = [b.sensors(:, 1) * s, b.sensors(:, 2)];    % physical positions fixed
        m.frfOptions = [stripOpt(b.frfOptions, 'WaterDepth'), {'WaterDepth', m.H}];
        m.band = b.band;                                       % the NOMINAL physical band: same sea and taper
    case 'radial'
        r = b.sensors(:, 1); moved = r > 0;
        r(moved) = r(moved) + value / b.H;                     % dr in metres, EMM length scale H
        m.sensors = [r, b.sensors(:, 2)];
    otherwise
        error('v1Tools:kind', 'Unknown mismatch kind ''%s''.', kind);
end
if any(m.sensors(:, 1) >= m.p(3)), error('v1Tools:offFloe', 'A perturbed sensor lies off the floe.'); end
if nargin >= 4 && ~isempty(cfg), m.alpha = cfg.emm.alphaValidated; end
end

function [fW, f, info] = welchFeatures(m, out, cacheDir)
% Welch-expected and exact features of model m on out's bins, exactly as evaluateScenario builds them
scn = out.scn;
% band [] = the validated band of the model's own depth, which is the nominal band for every case except depth
% (where truthModel sets the nominal band explicitly); default CheckPoints, so the fields of evaluateScenario and
% of diagnoseV1Truncation are reused from the cache
fld = emmField(m.p, m.sensors, 'ThetaIdx', [], 'NodeSpacing', m.nodeSpacing, 'FRFOptions', m.frfOptions, ...
    'Band', m.band, 'CacheDir', cacheDir, 'FRF', scn.frf);
omega = out.omega;
[~, f] = transmissibilityFromField(fld.H(omega), [], scn.layout, scn.ref);
sa = seaArgs(scn.sea);
[~, ~, info] = synthesiseTwinRecords(fld.twin, scn.N, scn.dt, 'Seed', 1, sa{:});
infoL = info; infoL.H = info.H(scn.layout, :);
others = setdiff(1:numel(scn.layout), scn.ref);
TW = welchExpectedTransmissibility(omega, infoL, scn.L, scn.ref, others, 0);
fW = stackTransmissibility(TW);
fW = fW(:); f = f(:);
end

function q = project(out, b)
% the E1 engine's projector on the bins the estimator accepts (identical to evaluateScenario's)
v = out.validMask;
pBlk = size(out.blocksAll, 1);
rows = reshape((find(v) - 1) * pBlk + (1:pBlk).', [], 1);
r = fisherFromBlocks(out.Jall(rows, :), out.blocksAll(:, :, v), b(rows));
q = struct('d', r.d, 'dtheta', r.dtheta(:));
end

function r = evaluateCase(b, out, fW0, kind, value, cfg, cacheDir)
m = truthModel(b, kind, value, cfg);
fW = welchFeatures(m, out, cacheDir);
bSys = fW - fW0;
qs = project(out, bSys);
qt = project(out, out.bW + bSys);
r = struct('kind', kind, 'value', value, 'dSys', qs.d, 'dthetaSys', qs.dtheta, ...
    'betaPct', 100 * (exp(qs.dtheta(1)) - 1), 'RPct', 100 * (exp(qs.dtheta(2)) - 1), ...
    'dTot', qt.d, 'dthetaTot', qt.dtheta, 'class', classify(out, qt.d, cfg), 'fW', fW, 'bSys', bSys, ...
    'alpha', alphaRange(m.H, out.omega(out.validMask), bandOr(m.band, b.band), cfg), 'p', m.p, 'sensors', m.sensors, ...
    'frfOptions', {m.frfOptions});
end

function c = classify(out, dTot, cfg)
% the pre-registered classification (as evaluateScenario), with the total bias distance
I = cfg.class.identifiable; M = cfg.class.marginal;
if out.singular || ~all(isfinite(out.sigma)), c = cfg.class.names{3}; return; end
if out.sMax <= I.maxSigma && dTot <= I.maxBiasDistance && out.condF <= I.maxCondF
    c = cfg.class.names{1};
elseif out.sMax <= M.maxSigma && dTot <= M.maxBiasDistance && out.condF <= M.maxCondF
    c = cfg.class.names{2};
else
    c = cfg.class.names{3};
end
end

function a = alphaRange(H, omegaRetained, band, cfg)
g = cfg.physics.g; v = cfg.emm.alphaValidated;
a.H = H;
a.minRetained = H * min(omegaRetained)^2 / g; a.maxRetained = H * max(omegaRetained)^2 / g;
a.minSupport = H * band(1)^2 / g; a.maxSupport = H * band(2)^2 / g;
a.retainedOutside = a.minRetained < v(1) || a.maxRetained > v(2);
a.extrapolated = a.minSupport < v(1) * (1 - 1e-6) || a.maxSupport > v(2) * (1 + 1e-6);
end

%% ================================================================================================
function nl = nonlinear(out, zW, thetaLin, cfg)
% G3 safeguarded Gauss-Newton (verbatim, g3EstimatorCore) on noiseless Welch features zW, nominal model,
% nominal Sigma, from theta_true. Compared with the TOTAL linear prediction theta_true + dth_W + dth_sys.
G = g3EstimatorCore(); S = cfg.V1.nonlinear.solver;
scn = out.scn; v = out.validMask; pBlk = size(out.blocksAll, 1);
rows = reshape((find(v) - 1) * pBlk + (1:pBlk).', [], 1);
memo = containers.Map('KeyType', 'char', 'ValueType', 'any');
sens = scn.sensors(scn.layout, :);
fun = @(th) modelAt(th, scn, sens, out.omega, rows, memo);
th0 = log(scn.p([1 3])).';
tt = tic;
J0 = G.jacobianCD(fun, th0, S.hStep);
z = reshape(zW(rows), pBlk, []);
g = G.gaussNewton(fun, z, out.blocksAll(:, :, v), th0, J0, out.sigma(:), S.hStep, S.maxIt, S.stepSigma, S.tolSigma);
gap = g.theta - thetaLin;
nl = struct('theta', g.theta, 'thetaTrue', th0, 'thetaLin', thetaLin, 'converged', g.converged, 'it', g.it, ...
    'evals', g.evals + 5, 'jRecomputes', g.jRecomputes, 'chi2', g.obj, 'seconds', toc(tt), ...
    'gapSigma', gap ./ out.sigma(:), 'dGap', sqrt(gap.' * out.F * gap), ...
    'betaPct', 100 * (exp(g.theta(1) - th0(1)) - 1), 'RPct', 100 * (exp(g.theta(2) - th0(2)) - 1), ...
    'linBetaPct', 100 * (exp(thetaLin(1) - th0(1)) - 1), 'linRPct', 100 * (exp(thetaLin(2) - th0(2)) - 1));
end

function f = modelAt(th, scn, sens, omega, rows, memo)
key = sprintf('%.15e_', th);
if isKey(memo, key), f = memo(key); return; end
p = scn.p; p([1 3]) = exp(th(:).');
args = {'ThetaIdx', [], 'NodeSpacing', scn.nodeSpacing, 'FRFOptions', scn.frfOptions, 'CheckPoints', 0};
if ~isempty(scn.frf), args = [args, {'FRF', scn.frf}]; end
fld = emmField(p, sens, args{:});
[~, f] = transmissibilityFromField(fld.H(omega), [], 1:size(sens, 1), scn.ref);
f = f(:); f = f(rows);
memo(key) = f; %#ok<NASGU>   containers.Map is a handle
end

%% ================================================================================================
function b = bandOr(b, d)
if isempty(b), b = d; end
end

function H = depthOf(opts)
H = 1.88;                                                   % computeSensorFRF default WaterDepth
k = find(strcmpi(opts(1:2:end), 'WaterDepth'), 1);
if ~isempty(k), H = opts{2 * k}; end
end

function o = stripOpt(opts, name)
o = opts;
k = find(strcmpi(o(1:2:end), name));
for i = fliplr(k), o(2 * i - 1:2 * i) = []; end
end

function a = seaArgs(sea)
% exactly as evaluateScenario passes the sea to synthesiseTwinRecords
if isa(sea, 'function_handle')
    a = {'Spectrum', sea};
else
    a = {'Spectrum', 'jonswap', 'Hs', sea.Hs, 'PeakFrequency', sea.omegaP, 'PeakEnhancement', sea.gammaJ};
end
end