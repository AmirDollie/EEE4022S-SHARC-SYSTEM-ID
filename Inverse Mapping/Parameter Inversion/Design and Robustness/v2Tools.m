function T = v2Tools()
%V2TOOLS  The numerical core of study V2 (directional forcing), shared by runDirectionalForcingStudy and
%   testDirectionalForcingStudy so the tester exercises exactly the code the run uses. D3 reuses headingPoints,
%   components, spreadWeights and context.
%   T = v2Tools();   then e.g.  ctx = T.context(out0, deg2rad([0 15 30 60 90]), cfg, cacheDir)
%
%   QUESTION. The frozen inverse assumes ONE incident direction (rank-one spectral matrix, S_jr / S_rr = H_j / H_r).
%   With a second independent component the matrix is rank two and the estimated transmissibility is a mixture.
%   How large is the parameter bias, and does the coherence drop reveal the failure before the bias matters?
%
%   DEFINITIONS (as V1: Welch expectation minus Welch expectation)
%       b_dir   = f_W(directional) - f_W(single)        both = stack(cv.T) of the engine's covariance function on
%                                                       directionalSpectralMatrix, so b_dir = 0 exactly at f2 = 0
%       dth_dir = F^-1 J' Sigma^-1 b_dir, d_sys         NOMINAL J, Sigma, F, nominal accepted bins (v1Tools.project)
%       dth_tot = dth_W + dth_dir, d_tot, class         vector sum (programme rule)
%   COHERENCE. gamma^2_jr = |S_jr|^2 / (S_jj S_rr) of the noisy spectral matrix, single-direction baseline vs the
%   directional sea. Fisher z = atanh(|gamma|). Per pair j the expected standardised drop over the nominal accepted bins
%       D_j = mean_k(z_single - z_dir) / sd_j
%   where sd_j is the standard deviation of the bin-mean Fisher z of the REAL estimator under the single-direction null,
%   calibrated by Monte Carlo (T.nullCoherence: records from the twin, extractTransmissibility coherence), so any
%   correlation between the retained bins is included. Without the calibration sd_j = 1/sqrt(2 nDof n_bins) (bins
%   independent, nDof = N_eff - 1; tests only). Detectable when max_j D_j >= z(1 - alpha/m) + z(power): the 80%-power
%   detectability threshold of a one-sided 5% test, Bonferroni over the m pairs (the rejection threshold itself is
%   z(1 - alpha/m)). The class uses the NOMINAL sigma (the frozen inverse's uncertainty ruler) with the directional
%   total bias; sigmaDir reports the scatter the estimator would actually have under the directional sea.
%
%   FUNCTIONS
%   pts  = T.headingPoints(sensors, headings)      the sensors rotated by -phi for each heading, stacked;
%                                                  pts.rows(h, :) are the rows of heading h (pts.sensors order)
%   c    = T.components(pts, layout, headings, fractions)   component struct array for directionalSpectralMatrix
%   c    = T.twoComponents(pts, layout, phi, f2)   primary at 0 with 1 - f2, secondary at phi with f2 (f2 = 0: one)
%   [h, w] = T.spreadWeights(meanHeading, s, nDir) discretised cos^2s((theta - mean) / 2) spread, weights sum 1 (D3)
%   ctx  = T.context(out0, headings, cfg, cacheDir, frfFor)   one field at all rotated points, record statistics,
%                                                  and the single-direction baseline through the directional path
%   r    = T.evaluate(ctx, out0, phi, f2, cfg)     one (phi, f2) case: bias, class, coherence statistic, sigma under
%                                                  the directional sea (diagnostic)
%   nc   = T.nullCoherence(out0, cfg, nRec, seedBase)   Monte Carlo null of the bin-mean Fisher z (per pair): sd, mean
%                                                  offset, the analytic sd and the effective number of independent bins
%   st   = T.coherenceStats(coh0, coh, valid, nDof, cfgCoh, sdNull)   the detection statistic (sdNull [] = analytic)
%   L    = T.linearity(ctx, out0, phi, f2s, cfg)   d_sys / f2 at small f2 (local regime check: constant if b ~ f2)
%   z    = T.normInv(p)                            standard normal quantile (no toolbox)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
T.headingPoints = @headingPoints;
T.components = @components;
T.twoComponents = @twoComponents;
T.spreadWeights = @spreadWeights;
T.context = @context;
T.evaluate = @evaluateCase;
T.coherenceStats = @coherenceStats;
T.nullCoherence = @nullCoherence;
T.linearity = @linearity;
T.normInv = @normInv;
T.version = 'v2Tools 2026-10-05b';
end

%% ================================================================================================
function pts = headingPoints(sensors, headings)
% H_phi(r, theta) = H_0(r, theta - phi): heading h samples the nominal field at the sensors rotated by -phi
headings = reshape(headings, 1, []);
nS = size(sensors, 1); nH = numel(headings);
S = zeros(nS * nH, 2);
rows = zeros(nH, nS);
for h = 1:nH
    k = (h - 1) * nS + (1:nS);
    S(k, :) = [sensors(:, 1), sensors(:, 2) - headings(h)];
    rows(h, :) = k;
end
pts = struct('sensors', S, 'rows', rows, 'headings', headings, 'nSensors', nS);
end

function c = components(pts, layout, headings, fractions)
headings = reshape(headings, 1, []); fractions = reshape(fractions, 1, []);
if numel(headings) ~= numel(fractions), error('v2Tools:components', 'One fraction per heading.'); end
c = struct('rows', {}, 'fraction', {}, 'heading', {});
for i = 1:numel(headings)
    h = find(abs(angle(exp(1i * (pts.headings - headings(i))))) < 1e-12, 1);
    if isempty(h), error('v2Tools:heading', 'Heading %.6g rad is not among the sampled headings.', headings(i)); end
    c(end + 1) = struct('rows', pts.rows(h, layout), 'fraction', fractions(i), 'heading', headings(i)); %#ok<AGROW>
end
end

function c = twoComponents(pts, layout, phi, f2)
if f2 == 0
    c = components(pts, layout, 0, 1);
else
    c = components(pts, layout, [0 phi], [1 - f2, f2]);
end
end

function [h, w] = spreadWeights(meanHeading, s, nDir)
% cos^2s((theta - mean) / 2) on nDir equally spaced headings over the full circle, normalised to unit energy
h = meanHeading + (0:nDir - 1) * 2 * pi / nDir;
h = angle(exp(1i * h));
d = angle(exp(1i * (h - meanHeading)));                       % wrapped difference: symmetric at any mean heading
w = cos(d / 2).^(2 * s);
w = w / sum(w);
end

function ctx = context(out0, headings, cfg, cacheDir, frfFor)
% one field (no derivatives) at the layout sensors rotated for every heading; same p, model, band, record
if nargin < 5, frfFor = []; end
scn = out0.scn;
if ~isempty(scn.frf) && isempty(frfFor)
    error('v2Tools:frf', 'An analytic scenario needs frfFor (a factory @(sensors) -> @(w, p) FRF).');
end
headings = reshape(headings, 1, []);
if ~any(abs(headings) < 1e-12), headings = [0, headings]; end   % the single-direction baseline is always sampled
pts = headingPoints(scn.sensors(scn.layout, :), headings);
args = {'ThetaIdx', [], 'NodeSpacing', scn.nodeSpacing, 'FRFOptions', scn.frfOptions, 'CacheDir', cacheDir};
if ~isempty(frfFor), args = [args, {'FRF', frfFor(pts.sensors)}]; end
tc = tic;
fld = emmField(scn.p, pts.sensors, args{:});
sa = seaArgs(scn.sea);
[~, ~, info] = synthesiseTwinRecords(fld.twin, scn.N, scn.dt, 'Seed', 1, sa{:});
lay = 1:pts.nSensors;                                         % pts already holds the layout sensors in layout order
Snn = scn.Snn; if isscalar(Snn), Snn = Snn * ones(numel(lay), 1); end
nEff = welchNEff(scn.N, scn.L, scn.overlap);
src = struct('omega', out0.omega, 'info', info, 'L', scn.L);
ctx = struct('pts', pts, 'layout', lay, 'ref', scn.ref, 'src', src, 'Snn', Snn(:), 'nEff', nEff, ...
    'valid', out0.validMask, 'fieldSeconds', toc(tc), 'cacheStatus', {fld.cacheStatus});
% single-direction baseline through the SAME path (directionalSpectralMatrix -> engine covariance)
D0 = directionalSpectralMatrix(src, components(pts, lay, 0, 1), Snn);
cv0 = analyticTransmissibilityCovariance(D0.Ssig, [], Snn, scn.ref, nEff);
ctx.fW0 = reshape(stackTransmissibility(cv0.T), [], 1);
ctx.coh0 = cv0.coherence;
ctx.blocks0 = cv0.blocks;
ctx.valid0 = cv0.valid;
ctx.nDof = cv0.nDof;
end

function r = evaluateCase(ctx, out0, phi, f2, cfg)
T1 = v1Tools();
comp = twoComponents(ctx.pts, ctx.layout, phi, f2);
D = directionalSpectralMatrix(ctx.src, comp, ctx.Snn);
cv = analyticTransmissibilityCovariance(D.Ssig, [], ctx.Snn, ctx.ref, ctx.nEff);
fW = reshape(stackTransmissibility(cv.T), [], 1);
b = fW - ctx.fW0;
qs = T1.project(out0, b);
qt = T1.project(out0, out0.bW(:) + b);
sdNull = []; if isfield(ctx, 'sdNull'), sdNull = ctx.sdNull; end
st = coherenceStats(ctx.coh0, cv.coherence, ctx.valid, ctx.nDof, cfg.V2.coherence, sdNull);
% diagnostic: the scatter the estimator would actually have under this sea (nominal J and accepted bins)
v = ctx.valid; pBlk = size(cv.blocks, 1);
rows = reshape((find(v) - 1) * pBlk + (1:pBlk).', [], 1);
qd = fisherFromBlocks(out0.Jall(rows, :), cv.blocks(:, :, v), zeros(numel(rows), 1));
r = struct('phiDeg', round(phi * 180 / pi * 1e9) / 1e9, 'f2', f2, 'dSys', qs.d, 'dthetaSys', qs.dtheta, ...
    'betaPct', 100 * (exp(qs.dtheta(1)) - 1), 'RPct', 100 * (exp(qs.dtheta(2)) - 1), ...
    'dTot', qt.d, 'dthetaTot', qt.dtheta, 'class', T1.classify(out0, qt.d, cfg), ...
    'D', st.D, 'Dmax', st.Dmax, 'threshold', st.threshold, 'detect', st.detect, 'perBinFrac', st.perBinFrac, ...
    'meanCoh0', st.meanCoh0, 'meanCoh', st.meanCoh, 'minCoh', st.minCoh, ...
    'sigmaDir', qd.sigma(:), 'sigmaRatio', qd.sigma(:) ./ out0.sigma(:), ...
    'nValidChanged', sum(cv.valid ~= ctx.valid0), 'fW', fW, 'b', b);
end

function st = coherenceStats(coh0, coh, valid, nDof, c, sdNull)
if nargin < 6, sdNull = []; end
m = size(coh0, 1);
nB = sum(valid);
z0 = atanh(sqrt(min(coh0(:, valid), 1 - 1e-15)));
z1 = atanh(sqrt(min(coh(:, valid), 1 - 1e-15)));
dz = z0 - z1;                                           % m x nB, positive = coherence lost
if isempty(sdNull), sdNull = ones(m, 1) / sqrt(2 * nDof * nB); end   % independent bins (uncalibrated)
D = mean(dz, 2) ./ sdNull(:);
thr = normInv(1 - c.alpha / m) + normInv(c.power);
st = struct('D', D, 'Dmax', max(D), 'threshold', thr, 'detect', max(D) >= thr, ...
    'perBinFrac', mean(dz * sqrt(2 * nDof) > normInv(c.perBinOneSided), 2), ...
    'meanCoh0', mean(coh0(:, valid), 2), 'meanCoh', mean(coh(:, valid), 2), 'minCoh', min(coh(:, valid), [], 2), ...
    'nBins', nB, 'm', m, 'sdNull', sdNull(:));
end

function nc = nullCoherence(out0, cfg, nRec, seedBase)
% bin-mean Fisher z of the REAL estimator's coherence (extractTransmissibility: periodic Hann, L, 50% overlap, every
% second bin, segment mean removed) on records of the single-direction twin with the nominal sea and noise, minus the
% analytic expectation z0: its spread is the null scatter that D is measured against (bins NOT assumed independent)
scn = out0.scn; v = out0.validMask;
lay = scn.layout; Snn = scn.Snn;
noiseStd = sqrt(pi * Snn / scn.dt);
sa = seaArgs(scn.sea);
z0 = atanh(sqrt(out0.coherence(:, v)));
m = size(z0, 1); nB = sum(v);
Z = zeros(m, nRec);
tc = tic;
for r = 1:nRec
    Y = synthesiseTwinRecords(out0.field.twin, scn.N, scn.dt, 'Seed', seedBase + r, sa{:}, 'NoiseStd', noiseStd);
    Y = Y(lay, :);
    ex = extractTransmissibility(Y, scn.dt, 'Reference', scn.ref, 'SegmentLength', scn.L, 'Overlap', scn.overlap, ...
        'Band', scn.bandEst, 'NoisePSD', Snn(min(end, scn.ref)));
    [hit, at] = ismember(out0.kBins, ex.binIndex);
    if ~all(hit), error('v2Tools:bins', 'The estimator does not return every retained bin.'); end
    c = ex.coherence(:, at);
    Z(:, r) = mean(atanh(sqrt(min(c(:, v), 1 - 1e-15))) - z0, 2);
end
nDof = out0.nEff - 1;
sdA = 1 / sqrt(2 * nDof * nB);
sd = std(Z, 0, 2);
nc = struct('nRec', nRec, 'seedBase', seedBase, 'sd', sd, 'mean', mean(Z, 2), 'sdAnalytic', sdA, ...
    'varianceRatio', (sd / sdA).^2, 'nBins', nB, 'nBinsEff', 1 ./ (2 * nDof * sd.^2), ...
    'sdOfSd', sd / sqrt(2 * (nRec - 1)), 'seconds', toc(tc));
end

function L = linearity(ctx, out0, phi, f2s, cfg)
% d_sys / f2 and |dbeta/beta| / f2 at small f2: in the local regime b_dir ~ f2, so both are constant
d = zeros(size(f2s)); b = d;
for i = 1:numel(f2s)
    r = evaluateCase(ctx, out0, phi, f2s(i), cfg);
    d(i) = r.dSys / f2s(i); b(i) = abs(r.betaPct) / f2s(i);
end
L = struct('phiDeg', round(phi * 180 / pi * 1e9) / 1e9, 'f2', f2s, 'dPerF2', d, 'betaPerF2', b, ...
    'spread', (max(d) - min(d)) / min(d));
end

function z = normInv(p)
z = sqrt(2) * erfinv(2 * p - 1);
end

function a = seaArgs(sea)
% exactly as evaluateScenario passes the sea to synthesiseTwinRecords
if isa(sea, 'function_handle')
    a = {'Spectrum', sea};
else
    a = {'Spectrum', 'jonswap', 'Hs', sea.Hs, 'PeakFrequency', sea.omegaP, 'PeakEnhancement', sea.gammaJ};
end
end