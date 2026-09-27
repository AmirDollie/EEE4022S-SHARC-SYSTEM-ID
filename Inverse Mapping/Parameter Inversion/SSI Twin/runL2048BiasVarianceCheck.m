%% runL2048BiasVarianceCheck.m
% Follow-up to U1 (runTransmissibilityCovariance): U1 found that the
% production setting (L = 1024, every second bin, 42 frequencies) has a
% correct covariance model but a Welch-resolution parameter bias of about
% one sigma (d_bias 0.96 Monte Carlo, 0.91 predicted). Doubling L from 512 to
% 1024 cut d_pred by about 4, so this script tests L = 2048 as a mitigation.
%
% NO NEW EMM SOLVES. The L = 2048 grid has half the L = 1024 spacing, so
%   every 4th L = 2048 bin  = the 42 production frequencies (J, f0 from U1)
%   every 2nd L = 2048 bin  = the 84 L1024-full frequencies (JFull, f0Full from U1)
% Both are verified by ismember against the saved U1 grids, not assumed.
%
% RECORDS: the same 200 records as U1 (same twin, N, JONSWAP, noise, seeds
% SEED0 + r), regenerated and checked: the L = 1024 production sigma
% recomputed here must reproduce the saved U1 value.
%
% CANDIDATES (all corrected, blkdiag(Sigma_k), same band)
%   B    L = 1024, every 2nd bin, 42 bins   (U1 production, reference)
%   E    L = 2048, every 4th bin, 42 bins   (same frequencies as B)
%   F    L = 2048, every 2nd bin, 84 bins   (uses all bins the lag structure
%        allows if L = 2048 adjacent-bin correlation is again confined to lag 1;
%        keeps the information that E throws away)
%
% PER CANDIDATE: nEff; median mean|corr| between RETAINED neighbours (floor
% sqrt(2/(pi N_MC))); sigma(ln beta), sigma(ln R) from all draws; held-out
% emp/pred (2-fold cross-fitting, within-fold scatter, prediction from the
% same half-sample covariance, as in U1); dtheta_bias and d_bias from the
% Monte Carlo bias and from the Welch-expected bias (known JONSWAP, so a
% diagnostic of this synthetic experiment, not an output-only correction).
%
% A CANDIDATE IS ACCEPTABLE IF
%   d_pred < 0.3 and d_MC < 0.3, emp/pred in [0.8, 1.25] for both parameters,
%   retained lag-1 correlation <= LAG_OK.
% Recommendation: the acceptable candidate with the smallest sigma(ln beta).
%
% OUTPUT (in Inverse Mapping/Parameter Inversion/Results/)
%   L2048BiasVarianceCheck_<stamp>.mat / .txt
%   out.cand(i) holds omega, f0, J, Sigma blocks and bias for the U2 choice.
%
% Lives in Inverse Mapping/Parameter Inversion/SSI Twin/.

thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..', '..', '..');
addpath(fullfile(thisDir, '..', '..'));
addpath(fullfile(rootDir, 'Transmissibility'));
addpath(fullfile(rootDir, 'SSI'));
addpath(fullfile(rootDir, 'Forward Model'));
addpath(fullfile(rootDir, 'Forward Model', 'Animation'));
clearvars -except thisDir rootDir
close all, clc

%% ---- Configuration -----------------------------------------------------------
U1_FILE = '';                                                % '' = newest transmissibilityCovariance_*.mat
SEED0 = 7000;                                                % must match runTransmissibilityCovariance
BAND_EST = [3 8.5]; OVERLAP = 0.5;
HS = 0.05; WP = 5.0; GAMMA_J = 3.3;                          % must match U1
D_BIAS_OK = 0.3; EMP_RANGE = [0.8 1.25]; LAG_OK = 0.10;

%% ---- Output ----------------------------------------------------------------------
resultsDir = fullfile(thisDir, '..', 'Results');
cacheDir = fullfile(resultsDir, 'twinCache');
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['L2048BiasVarianceCheck_', stamp]);
diary([baseName, '.txt']);
try

if isempty(U1_FILE)
    d = dir(fullfile(resultsDir, 'transmissibilityCovariance_*.mat'));
    if isempty(d), error('runL2048BiasVarianceCheck:noU1', 'No U1 result in %s.', resultsDir); end
    [~, newest] = max([d.datenum]);
    U1_FILE = fullfile(resultsDir, d(newest).name);
end
U = load(U1_FILE); U = U.out;
P0 = U.p0; SENSORS = U.sensors; REF = U.reference; N = U.N; DT = U.dt; N_MC = U.nMC;
nS = size(SENSORS, 1); m = nS - 1; p = 2 * m; others = setdiff(1:nS, REF);
fprintf('runL2048BiasVarianceCheck  (%s)\n', stamp);
fprintf('U1 input: %s\n', U1_FILE);
fprintf('N = %d, N_MC = %d, noise std %.4g, S_nn = %.3e\n\n', N, N_MC, U.prod.noiseStd, U.prod.Snn);

%% ---- 1  Twin and grids ------------------------------------------------------------------------
twin = synthesiseTwinRecords(struct('p', P0, 'sensors', SENSORS, 'nodeSpacing', 0.1), 'CacheDir', cacheDir);
specJ = {'Spectrum', 'jonswap', 'Hs', HS, 'PeakFrequency', WP, 'PeakEnhancement', GAMMA_J};
[~, ~, iJ] = synthesiseTwinRecords(twin, N, DT, 'Seed', 1, specJ{:});
assert(abs(iJ.signalStd(REF) * U.noiseRel - U.prod.noiseStd) < 1e-12 * U.prod.noiseStd, ...
    'Twin / spectrum differ from U1 (noise level does not reproduce)');
Snn = U.prod.Snn; noiseStd = U.prod.noiseStd;
opt = @(L) {'Reference', REF, 'SegmentLength', L, 'Overlap', OVERLAP, 'Band', BAND_EST, 'NoisePSD', Snn};
o1024 = opt(1024); o2048 = opt(2048);
e0 = extractTransmissibility(randn(nS, 1024), DT, o1024{1:8}); w1024 = e0.omega;
e0 = extractTransmissibility(randn(nS, 2048), DT, o2048{1:8}); w2048 = e0.omega;
key = @(w) round(w(:) * 1e9);
[okB, idxB] = ismember(key(U.omegaDec), key(w1024));
[okE, idxE] = ismember(key(U.omegaDec), key(w2048));
[okF, idxF] = ismember(key(U.omegaFull), key(w2048));
assert(all(okB) && all(okE) && all(okF), 'Grids do not align: L = 2048 bins do not contain the U1 frequencies');
assert(all(diff(idxE) == 4) && all(diff(idxF) == 2), 'L = 2048 retained bins are not every 4th / every 2nd');
fprintf('grid check: the %d production frequencies are every 4th L = 2048 bin, the %d L1024-full\n', ...
    numel(U.omegaDec), numel(U.omegaFull));
fprintf('frequencies every 2nd L = 2048 bin (verified); U1 Jacobians reused, no EMM solves\n\n');

cand = struct('name', {'B  L1024 every 2nd (U1 production)', 'E  L2048 every 4th', 'F  L2048 every 2nd'}, ...
    'L', {1024, 2048, 2048}, 'idx', {idxB(:).', idxE(:).', idxF(:).'});
cand(1).omega = U.omegaDec; cand(1).f0 = U.prod.f0; cand(1).J = U.prod.J;
cand(2).omega = U.omegaDec; cand(2).f0 = U.prod.f0; cand(2).J = U.prod.J;
cand(3).omega = U.omegaFull; cand(3).f0 = U.f0Full; cand(3).J = U.JFull;

%% ---- 2  Same records, two segment lengths ---------------------------------------------------------------
t0 = tic;
Z1 = NaN(p, numel(w1024), N_MC); Z2 = NaN(p, numel(w2048), N_MC);
for r = 1:N_MC
    Y = synthesiseTwinRecords(twin, N, DT, 'Seed', SEED0 + r, specJ{:}, 'NoiseStd', noiseStd);
    ea = extractTransmissibility(Y, DT, o1024{:});
    eb = extractTransmissibility(Y, DT, o2048{:});
    Z1(:, :, r) = reshape(stackTransmissibility(ea.T), p, []);
    Z2(:, :, r) = reshape(stackTransmissibility(eb.T), p, []);
end
nEff = [ea.nEffective, eb.nEffective, eb.nEffective];
fprintf('Monte Carlo: %d records, nEff %.1f (L = 1024), %.1f (L = 2048), %.0f s\n', N_MC, nEff(1:2), toc(t0));
zp1 = welchExpectedZ(w1024, iJ, 1024, DT, REF, others, 0);
zp2 = welchExpectedZ(w2048, iJ, 2048, DT, REF, others, 0);
cand(1).Z = Z1(:, cand(1).idx, :); cand(1).zPred = zp1(:, cand(1).idx);
cand(2).Z = Z2(:, cand(2).idx, :); cand(2).zPred = zp2(:, cand(2).idx);
cand(3).Z = Z2(:, cand(3).idx, :); cand(3).zPred = zp2(:, cand(3).idx);

%% ---- 3  Analysis ----------------------------------------------------------------------------------------------
floorRho = sqrt(2 / (pi * N_MC));
for c = 1:numel(cand)
    Zc = cand(c).Z;
    okDraw = squeeze(~any(any(isnan(Zc), 1), 2));
    Zc = Zc(:, :, okDraw); cand(c).nDropped = sum(~okDraw);
    D = reshape(Zc, [], size(Zc, 3)) - cand(c).f0;
    cand(c).b = mean(D, 2);
    cand(c).bPred = cand(c).zPred(:) - cand(c).f0;
    cand(c).Sigma = buildSigma(Zc);
    cand(c).res = infoAnalysis(cand(c).J, cand(c).Sigma, cand(c).b, cand(c).bPred);
    [cand(c).sigmaEmp, cand(c).sigmaPredCV] = crossFit(Zc, D, cand(c).J);
    cand(c).lag1 = lag1Corr(Zc);
    cand(c).nEff = nEff(c);
    cand(c).blocks = blocksOf(Zc);
    q = cand(c).res;
    ep = cand(c).sigmaEmp ./ cand(c).sigmaPredCV;
    cand(c).ok = q.dBias < D_BIAS_OK && q.dBiasPred < D_BIAS_OK && all(ep >= EMP_RANGE(1) & ep <= EMP_RANGE(2)) ...
        && cand(c).lag1 <= LAG_OK;
    cand(c).Z = [];                                                    % not saved (size)
end
repro = max(abs(cand(1).res.sigma - U.treatments(2).res.sigma) ./ U.treatments(2).res.sigma);
fprintf('record check: L = 1024 production sigma reproduces U1 to %.1e (relative)\n\n', repro);
if repro > 1e-6
    warning('runL2048BiasVarianceCheck:notSameRecords', ...
        'Production sigma differs from U1 by %.2e: records are not identical (check SEED0).', repro);
end

fprintf('BIAS-VARIANCE CHECK (one record of %.0f min; lag-1 floor %.2f)\n', N * DT / 60, floorRho);
fprintf('%-36s %3s %6s %5s | %8s %8s %6s | %7s %7s | %9s %9s | %6s %6s | %s\n', 'candidate', 'nF', 'nEff', 'lag1', ...
    'sig lnB', 'sig lnR', 'corr', 'emp/pB', 'emp/pR', 'dth_B', 'dth_R', 'd_MC', 'd_pred', 'ok');
for c = 1:numel(cand)
    q = cand(c).res; ep = cand(c).sigmaEmp ./ cand(c).sigmaPredCV;
    fprintf('%-36s %3d %6.1f %5.2f | %8.5f %8.5f %6.3f | %7.2f %7.2f | %+9.5f %+9.5f | %6.2f %6.2f | %s\n', cand(c).name, ...
        numel(cand(c).omega), cand(c).nEff, cand(c).lag1, q.sigma, q.corr, ep, q.dthetaPred, q.dBias, q.dBiasPred, ...
        flagStr(cand(c).ok));
end
fprintf('(dth = Welch-expected parameter bias in ln units; d_MC noise floor sqrt(2/N_MC) = %.2f)\n\n', sqrt(2 / N_MC));
fprintf('In percent (1 sigma, 100(exp(sigma)-1)) and bias (100(exp(dth)-1)):\n');
for c = 1:numel(cand)
    q = cand(c).res;
    fprintf('  %-36s beta %.2f%% (bias %+.2f%%), R %.3f%% (bias %+.3f%%)\n', cand(c).name, 100 * (exp(q.sigma(1)) - 1), ...
        100 * (exp(q.dthetaPred(1)) - 1), 100 * (exp(q.sigma(2)) - 1), 100 * (exp(q.dthetaPred(2)) - 1));
end
rB = cand(1).res; rE = cand(2).res;
fprintf('\nL 1024 -> 2048 on the same production frequencies: d_pred x %.2f, sigma x [%.2f %.2f]\n', ...
    rE.dBiasPred / rB.dBiasPred, rE.sigma ./ rB.sigma);
okList = find([cand.ok]);
if isempty(okList)
    fprintf('DECISION: no candidate meets all criteria; U2 needs an explicit treatment of the Welch bias.\n');
    choice = NaN;
else
    [~, i] = min(arrayfun(@(c) cand(c).res.sigma(1), okList));
    choice = okList(i);
    fprintf('DECISION: use %s for U2 (acceptable, smallest sigma(ln beta)).\n', cand(choice).name);
end

%% ---- 4  Save --------------------------------------------------------------------------------------------------
out = struct('u1File', U1_FILE, 'p0', P0, 'sensors', SENSORS, 'reference', REF, 'theta0', U.theta0, 'N', N, ...
    'dt', DT, 'nMC', N_MC, 'overlap', OVERLAP, 'Snn', Snn, 'noiseStd', noiseStd, 'cand', cand, 'choice', choice, ...
    'reproducibility', repro, 'criteria', struct('dBias', D_BIAS_OK, 'empRange', EMP_RANGE, 'lag1', LAG_OK));
save([baseName, '.mat'], 'out');
fprintf('Saved: %s.mat\n', baseName);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ================================================================================================
%  Local functions
%% ================================================================================================
function [sigEmp, sigPredCV] = crossFit(Z, D, J)
% 2-fold cross-fitting: covariance from one half, linearised estimates on the
% other; within-fold scatter against the same half-sample F_h^-1 (as in U1).
nD = size(D, 2); half = {1:floor(nD / 2), floor(nD / 2) + 1:nD};
varCV = zeros(2, 1); ss = zeros(2, 1); df = 0;
for h = 1:2
    fitIdx = half{h}; appIdx = half{3 - h};
    L = chol(buildSigma(Z(:, :, fitIdx)), 'lower');
    Jw = L \ J; Fh = Jw.' * Jw;
    th = Fh \ (Jw.' * (L \ D(:, appIdx)));
    varCV = varCV + numel(appIdx) * diag(inv(Fh));
    th = th - mean(th, 2);
    ss = ss + sum(th.^2, 2); df = df + numel(appIdx) - 1;
end
sigEmp = sqrt(ss / df);
sigPredCV = sqrt(varCV / nD);
end

function rho = lag1Corr(Z)
% Median over retained neighbours of mean |corr| between matching components.
[p, nF, ~] = size(Z);
r = NaN(1, nF - 1);
for k = 1:nF - 1
    A = reshape(Z(:, k, :), p, []); B = reshape(Z(:, k + 1, :), p, []);
    A = A - mean(A, 2); B = B - mean(B, 2);
    r(k) = mean(abs(sum(A .* B, 2) ./ sqrt(sum(A.^2, 2) .* sum(B.^2, 2))));
end
rho = median(r);
end

function B = blocksOf(Z)
[p, nF, ~] = size(Z);
B = zeros(p, p, nF);
for k = 1:nF
    B(:, :, k) = cov(reshape(Z(:, k, :), p, []).');
end
end

function S = buildSigma(Z)
% Block-diagonal covariance from the per-bin sample covariances.
[p, nF, ~] = size(Z);
S = zeros(p * nF);
blocks = blocksOf(Z);
for k = 1:nF
    ix = (k - 1) * p + (1:p);
    S(ix, ix) = blocks(:, :, k);
end
end

function q = infoAnalysis(J, S, b, bPred)
L = chol(S, 'lower');
Jw = L \ J;
F = Jw.' * Jw;
C = inv(F);
q.F = F; q.covTheta = C;
q.sigma = sqrt(diag(C));
q.corr = C(1, 2) / (q.sigma(1) * q.sigma(2));
q.dtheta = F \ (Jw.' * (L \ b));
q.dBias = sqrt(q.dtheta.' * F * q.dtheta);
q.dthetaPred = F \ (Jw.' * (L \ bPred));
q.dBiasPred = sqrt(q.dthetaPred.' * F * q.dthetaPred);
end

function z = welchExpectedZ(wBins, info, L, dt, ref, others, SnnDen)
% Welch-expected transmissibility (as in runTransmissibilityEstimatorTest).
wq = info.omega;
G = info.inputPSD .* info.taper.^2;
Hq = info.H;
U = 3 * L / 8;
K = abs(hannDTFT((wBins(:) - wq) * dt, L)).^2 / (info.N * U);
Srr = K * (G .* abs(Hq(ref, :)).^2).' + SnnDen;
Sjr = K * (G .* Hq(others, :) .* conj(Hq(ref, :))).';
T = (Sjr ./ Srr).';
z = reshape(stackTransmissibility(T), 2 * numel(others), numel(wBins));
end

function W = hannDTFT(nu, L)
W = 0.5 * dirichletSum(nu, L) - 0.25 * dirichletSum(nu - 2 * pi / L, L) - 0.25 * dirichletSum(nu + 2 * pi / L, L);
end

function D = dirichletSum(x, L)
den = 1 - exp(-1i * x);
D = (1 - exp(-1i * x * L)) ./ den;
small = abs(den) < 1e-9;
D(small) = L;
end

function s = flagStr(ok)
if ok, s = 'PASS'; else, s = 'FLAG'; end
end