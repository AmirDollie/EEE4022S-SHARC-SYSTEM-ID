%% runTransmissibilityCovariance.m
% U1: finite-record information analysis. Converts the real Monte Carlo
% statistics of the output-only transmissibility features into predicted
% uncertainty and bias of theta = [ln beta; ln R] (gamma fixed):
%
%   Monte Carlo T_hat -> Sigma_T, b -> J_T -> F = J' Sigma^-1 J -> sigma(theta), bias(theta)
%
% DATA: one fresh production Monte Carlo (S3 twin -> S2 estimator):
%   N = 65536, JONSWAP, 10% sensor-noise stress case, reference-noise
%   correction ON, N_MC = 200. Every record is analysed twice, with
%   L = 1024 (production) and L = 512 (benchmark), so all treatments see the
%   same records. The corrected L = 1024 case was not part of Saturday's sweep.
%
% FREQUENCY SETS (all inside the ramp-free band BAND_USE, which excludes the
% S3 taper ramps; a twin artefact, not physics)
%   L1024-full : every L = 1024 Welch bin in BAND_USE
%   L1024-dec  : every second bin of the L = 1024 grid on [3, 8.5] (Saturday's
%                phase: first bin 3.007 rad/s), then trimmed: the PRODUCTION grid
%                (42 bins, 3.252 to 8.284 rad/s)
%   L512       : every L = 512 bin in BAND_USE (these are the even-k bins of
%                the L = 1024 grid, so one Jacobian serves every set)
%
% COVARIANCE TREATMENTS
%   A  L512, all bins, blkdiag(Sigma_k)               optimistic benchmark (Friday style)
%   B  L1024-dec, blkdiag(Sigma_k)                    PRODUCTION
%   C  L1024-full, block tridiagonal (Sigma_k, C_k,k+1) complete benchmark
%   D  L1024-full, blkdiag(Sigma_k)                   deliberately wrong (ignores lag 1)
%   C is assembled from the per-bin blocks and the adjacent 12 x 12 pair
%   covariances; if it is not positive definite the off-diagonal blocks are
%   shrunk (C_k,k+1 -> a C_k,k+1, largest a in 1, 0.95, ...) and a is reported.
%
% PER TREATMENT
%   F = J' Sigma^-1 J (Cholesky whitening), Cov(theta) ~ F^-1,
%   sigma_pred = sqrt(diag F^-1), corr(ln beta, ln R), cond(F)
%   sigma_emp  = scatter of the LINEARISED estimator
%                theta_hat = theta0 + F^-1 J' Sigma^-1 (z - f0)
%                applied to every Monte Carlo draw, with 2-fold CROSS-FITTING
%                (Sigma estimated on one half of the draws, applied to the
%                other), so weights and draws are independent. sigma_emp /
%                sigma_pred,CV ~ 1 means the covariance model is right; > 1 means
%                it overstates the information (expected for D).
%   bias projection, for the Monte Carlo bias b = mean(z) - f0 (f0 = S1 EMM
%   features at theta0, so interpolation error is included) and for the
%   deterministic Welch-expected bias b_pred (Saturday's Eq. for T_pred with
%   the correction on, no Monte Carlo noise):
%     dtheta_bias = F^-1 J' Sigma^-1 b,   ratio dtheta_bias ./ sigma_pred,
%     d_bias = sqrt(dtheta_bias' F dtheta_bias)
%   The Monte Carlo d_bias has a noise floor: with no true bias,
%   E[d_bias^2] = 2 / N_MC (floor sqrt(2/N_MC) = 0.10 for 200 draws).
%   b_pred uses the twin's KNOWN JONSWAP spectrum: it is the true large-sample
%   resolution bias of this synthetic experiment, not a correction available to
%   the output-only field problem (that needs S_u: Scenario W or P).
%   If d_bias ~ 1: a shorter record would only hide the offset inside a wider
%   confidence region; the remedies are larger L, a Welch observation operator
%   built from information about S_u, or an explicit systematic-error term.
%
% DECISIONS PRINTED
%   U1.1 decimation cost: sigma_pred(B) / sigma_pred(C) <= 1.2 for both parameters
%   U1.2 covariance model: sigma_emp / sigma_pred,CV within [0.8, 1.25] for B and C, where
%        sigma_pred,CV = sqrt(sum_h n_h diag(F_h^-1) / n) uses the SAME half-sample covariance
%        that weighted the held-out estimates (apples to apples); the reported sigma
%        itself uses all draws
%   C is a reliable benchmark only if its off-diagonal shrinkage a >= 0.8 (all draws
%   and both folds); otherwise U1.1 flags and C must be inspected first
%   U1.3 parameter bias: d_bias(B) < 0.3 (Monte Carlo and predicted)
%   U1.4 Jacobian step: J with h and h/2 agree; sigma changes < 1%
%
% EMM COST: S1 on the L1024-full grid (about 83 frequencies) at theta0 and
% theta0 +- h e_i (5 x 83 solves) plus the step check (4 x 42): about 10 min in
% MATLAB, cached in Results/twinCache (instant afterwards). Set USE_PARALLEL to
% use parfor inside computeSensorFRF.
%
% OUTPUT (in Inverse Mapping/Parameter Inversion/Results/)
%   transmissibilityCovariance_<stamp>.mat / .txt / _ellipses.png
%   out.prod holds the production grid, f0, J, Sigma_k and bias for U2.
%
% Lives in Inverse Mapping/Parameter Inversion/SSI Twin/.

thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..', '..', '..');
addpath(fullfile(thisDir, '..', '..'));                        % transmissibilityFeatureVector, computeSensorFRF
addpath(fullfile(rootDir, 'Transmissibility'));
addpath(fullfile(rootDir, 'SSI'));
addpath(fullfile(rootDir, 'Forward Model'));
addpath(fullfile(rootDir, 'Forward Model', 'Animation'));
clearvars -except thisDir rootDir
close all, clc

%% ---- Configuration -----------------------------------------------------------
P0 = [4.6985e-5, 1.4548e-3, 0.3830]; R0 = P0(3);
SENSORS = [0.3*R0, 0; 0.3*R0, pi; 0.5*R0, pi; 0.9*R0, 0];    % Level 2C: s2, s6, s14, s26
REF = 4;
THETA_IDX = [1 3];                                           % [beta R]; gamma fixed
DT = 0.1; OVERLAP = 0.5; N = 65536; N_MC = 200;
L_PROD = 1024; L_BENCH = 512;
BAND_EST = [3 8.5];                                          % estimator band (Saturday)
BAND_USE = [3.178 8.301];                                    % excludes the S3 taper ramps
HS = 0.05; WP = 5.0; GAMMA_J = 3.3; NOISE_REL = 0.10;
H_STEP = 1e-3;                                               % central-difference step in ln-parameters
STEP_CHECK = true;                                           % also J with H_STEP/2 on the production grid
USE_PARALLEL = false;
SEED0 = 7000;
D_BIAS_OK = 0.3; DEC_COST_OK = 1.2; EMP_RANGE = [0.8 1.25]; SHRINK_OK = 0.8;

%% ---- Output ----------------------------------------------------------------------
resultsDir = fullfile(thisDir, '..', 'Results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
cacheDir = fullfile(resultsDir, 'twinCache');
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['transmissibilityCovariance_', stamp]);
diary([baseName, '.txt']);
try

nS = size(SENSORS, 1); m = nS - 1; others = setdiff(1:nS, REF);
theta0 = log(P0(THETA_IDX)).';
frfOpt = {};
if USE_PARALLEL, frfOpt = {'UseParallel', true}; end
fprintf('runTransmissibilityCovariance  (%s)\n', stamp);
fprintf('p0 = [%.4e %.4e %.4f]; theta = [ln beta, ln R]; reference s26\n', P0);
fprintf('N = %d (%.0f s), N_MC = %d, JONSWAP Hs %.2f wp %.1f, noise %.0f%%, correction ON\n\n', ...
    N, N * DT, N_MC, HS, WP, 100 * NOISE_REL);

%% ---- 1  Twin, noise level, frequency sets ------------------------------------------------
twin = synthesiseTwinRecords(struct('p', P0, 'sensors', SENSORS, 'nodeSpacing', 0.1), 'CacheDir', cacheDir);
specJ = {'Spectrum', 'jonswap', 'Hs', HS, 'PeakFrequency', WP, 'PeakEnhancement', GAMMA_J};
[~, ~, iJ] = synthesiseTwinRecords(twin, N, DT, 'Seed', 1, specJ{:});
noiseStd = NOISE_REL * iJ.signalStd(REF);
Snn = noiseStd^2 * DT / pi;
fprintf('twin %s; noise std %.4g m/s^2, S_nn = %.3e\n', twin.cacheStatus, noiseStd, Snn);

opt1024 = {'Reference', REF, 'SegmentLength', L_PROD, 'Overlap', OVERLAP, 'Band', BAND_EST, 'NoisePSD', Snn};
opt512 = {'Reference', REF, 'SegmentLength', L_BENCH, 'Overlap', OVERLAP, 'Band', BAND_EST, 'NoisePSD', Snn};
e1 = extractTransmissibility(randn(nS, L_PROD), DT, opt1024{1:8});
e5 = extractTransmissibility(randn(nS, L_PROD), DT, opt512{1:8});
w1024 = e1.omega; w512 = e5.omega;
inUse = @(w) w >= BAND_USE(1) & w <= BAND_USE(2);
sets.full = find(inUse(w1024));                               % indices into the L = 1024 estimator grid
dec = 1:2:numel(w1024);
sets.dec = dec(inUse(w1024(dec)));
sets.b512 = find(inUse(w512));                                % indices into the L = 512 estimator grid
wFull = w1024(sets.full); wDec = w1024(sets.dec); w512u = w512(sets.b512);
% L = 512 bins must be L = 1024 bins (even k), so one Jacobian serves all sets
[isSub, loc512] = ismember(round(w512u * 1e9), round(wFull * 1e9));
assert(all(isSub), 'L = 512 bins are not a subset of the L = 1024 grid');
[~, locDec] = ismember(round(wDec * 1e9), round(wFull * 1e9));
fprintf('frequency sets in [%.3f %.3f] rad/s: L1024-full %d bins, L1024-dec %d bins (%.3f to %.3f), L512 %d bins\n\n', ...
    BAND_USE, numel(wFull), numel(wDec), wDec([1 end]), numel(w512u));

%% ---- 2  Jacobian through S1 (cached) ---------------------------------------------------------
t0 = tic;
[f0, J, jstat] = cachedJacobian(cacheDir, theta0, wFull, SENSORS, P0, THETA_IDX, REF, H_STEP, frfOpt);
fprintf('Jacobian on %d frequencies: %s (%.0f s)\n', numel(wFull), jstat, toc(t0));
p = 2 * m;                                                   % features per frequency
Jb = reshape(J, p, numel(wFull), 2);                         % p x nF x 2
f0b = reshape(f0, p, numel(wFull));
grab = @(idx) deal(reshape(f0b(:, idx), [], 1), reshape(Jb(:, idx, :), [], 2));
[f0Dec, JDec] = grab(locDec);
[f0Full, JFull] = grab(1:numel(wFull));
[f0512, J512] = grab(loc512);
if STEP_CHECK
    t0 = tic;
    [~, Jh2] = cachedJacobian(cacheDir, theta0, wDec, SENSORS, P0, THETA_IDX, REF, H_STEP / 2, frfOpt);
    stepDiff = norm(Jh2 - JDec, 'fro') / norm(JDec, 'fro');
    fprintf('step check (h = %.0e vs %.0e, production grid): relative difference %.2e (%.0f s)\n', ...
        H_STEP, H_STEP / 2, stepDiff, toc(t0));
else
    Jh2 = []; stepDiff = NaN;
end
fprintf('\n');

%% ---- 3  Production Monte Carlo ------------------------------------------------------------------
t0 = tic;
Z1024 = NaN(p, numel(w1024), N_MC); Z512 = NaN(p, numel(w512), N_MC);
for r = 1:N_MC
    Y = synthesiseTwinRecords(twin, N, DT, 'Seed', SEED0 + r, specJ{:}, 'NoiseStd', noiseStd);
    ea = extractTransmissibility(Y, DT, opt1024{:});
    eb = extractTransmissibility(Y, DT, opt512{:});
    Z1024(:, :, r) = reshape(stackTransmissibility(ea.T), p, []);
    Z512(:, :, r) = reshape(stackTransmissibility(eb.T), p, []);
end
fprintf('Monte Carlo: %d records, L = %d (nEff %.1f) and L = %d (nEff %.1f), %.0f s\n', N_MC, L_PROD, ...
    ea.nEffective, L_BENCH, eb.nEffective, toc(t0));
ZFull = Z1024(:, sets.full, :); ZDec = Z1024(:, sets.dec, :); Z512u = Z512(:, sets.b512, :);
nanDraws = [sum(any(any(isnan(ZFull), 1), 2)), sum(any(any(isnan(Z512u), 1), 2))];
fprintf('draws with an invalid (NaN) bin: L1024 %d, L512 %d (dropped)\n', nanDraws);
okFull = squeeze(~any(any(isnan(ZFull), 1), 2)); ok512 = squeeze(~any(any(isnan(Z512u), 1), 2));
ZFull = ZFull(:, :, okFull); ZDec = ZDec(:, :, okFull); Z512u = Z512u(:, :, ok512);

% Welch-expected (deterministic) bias, correction on: no S_nn in the denominator
zp1024 = welchExpectedZ(w1024, iJ, L_PROD, DT, REF, others, 0);
zp512 = welchExpectedZ(w512, iJ, L_BENCH, DT, REF, others, 0);
bPredFull = reshape(zp1024(:, sets.full), [], 1) - f0Full;
bPredDec = reshape(zp1024(:, sets.dec), [], 1) - f0Dec;
bPred512 = reshape(zp512(:, sets.b512), [], 1) - f0512;

%% ---- 4  Treatments --------------------------------------------------------------------------------
tr = struct('name', {'A  L512 blkdiag', 'B  L1024-dec blkdiag (production)', 'C  L1024-full lag-1 banded', ...
    'D  L1024-full blkdiag (ignores lag 1)'}, 'key', {'A', 'B', 'C', 'D'});
tr(1).Z = Z512u; tr(1).f0 = f0512; tr(1).J = J512; tr(1).bPred = bPred512; tr(1).band = false;
tr(2).Z = ZDec;  tr(2).f0 = f0Dec; tr(2).J = JDec; tr(2).bPred = bPredDec; tr(2).band = false;
tr(3).Z = ZFull; tr(3).f0 = f0Full; tr(3).J = JFull; tr(3).bPred = bPredFull; tr(3).band = true;
tr(4).Z = ZFull; tr(4).f0 = f0Full; tr(4).J = JFull; tr(4).bPred = bPredFull; tr(4).band = false;
for t = 1:numel(tr)
    [tr(t).Sigma, tr(t).shrink] = buildSigma(tr(t).Z, tr(t).band);
    D = reshape(tr(t).Z, [], size(tr(t).Z, 3)) - tr(t).f0;          % dim x nDraws
    tr(t).b = mean(D, 2);
    tr(t).res = infoAnalysis(tr(t).J, tr(t).Sigma, tr(t).b, tr(t).bPred);
    % 2-fold cross-fitted empirical scatter of the linearised estimator
    nD = size(D, 2); half = {1:floor(nD / 2), floor(nD / 2) + 1:nD};
    thetaHat = zeros(2, nD); varCV = zeros(2, 1); ssEmp = zeros(2, 1); dfEmp = 0; tr(t).shrinkFold = [1 1];
    for h = 1:2
        fitIdx = half{h}; appIdx = half{3 - h};
        [Sh, tr(t).shrinkFold(h)] = buildSigma(tr(t).Z(:, :, fitIdx), tr(t).band);
        Lh = chol(Sh, 'lower');
        Jw = Lh \ tr(t).J;
        Fh = Jw.' * Jw;
        thetaHat(:, appIdx) = Fh \ (Jw.' * (Lh \ D(:, appIdx)));
        varCV = varCV + numel(appIdx) * diag(inv(Fh));               % what this fold's estimator claims
        Xh = thetaHat(:, appIdx) - mean(thetaHat(:, appIdx), 2);     % within-fold scatter only: the two folds'
        ssEmp = ssEmp + sum(Xh.^2, 2);                               % weights give slightly different mean shifts
        dfEmp = dfEmp + numel(appIdx) - 1;
    end
    tr(t).thetaHat = thetaHat;                                        % deviations from theta0 (plotted raw)
    tr(t).sigmaEmp = sqrt(ssEmp / dfEmp);                             % pooled within-fold std
    tr(t).sigmaPredCV = sqrt(varCV / nD);                             % same weights as the held-out estimates
    tr(t).nFreq = size(tr(t).Z, 2); tr(t).dim = numel(tr(t).f0);
end

%% ---- 5  Results -----------------------------------------------------------------------------------------
fprintf('\nU1 RESULTS (theta = [ln beta, ln R]; sigma for ONE record of %.0f min)\n', N * DT / 60);
fprintf('%-38s %4s %4s | %9s %9s %7s %9s | %9s %9s | %8s %8s %7s | %7s\n', 'treatment', 'nF', 'dim', ...
    'sig lnB', 'sig lnR', 'corr', 'cond F', 'emp/pred B', 'emp/pred R', 'bias/sig B', 'bias/sig R', 'd_bias', 'd_pred');
for t = 1:numel(tr)
    q = tr(t).res;
    fprintf('%-38s %4d %4d | %9.4f %9.4f %7.3f %9.2e | %9.2f %9.2f | %8.2f %8.2f %7.2f | %7.2f\n', tr(t).name, ...
        tr(t).nFreq, tr(t).dim, q.sigma, q.corr, q.cond, tr(t).sigmaEmp ./ tr(t).sigmaPredCV, q.ratio, q.dBias, q.dBiasPred);
end
fprintf(['(sig: from the covariance of all draws, the final uncertainty estimate; emp/pred: held-out scatter of the\n' ...
    ' linearised estimator over the sigma claimed by the SAME half-sample covariance that weighted it, F_h^-1;\n' ...
    ' expected slightly above 1 even for a correct model (half-sample weights are noisy and F_h is optimistic);\n' ...
    ' bias/sig and d_bias from the Monte Carlo bias, noise floor of d_bias sqrt(2/N_MC) = %.2f; d_pred from the Welch-expected bias)\n'], sqrt(2 / N_MC));
C_RELIABLE = min([tr(3).shrink, tr(3).shrinkFold]) >= SHRINK_OK;
fprintf('*** C shrinkage a (lag-1 blocks x a for positive definiteness): all draws %.2f, folds %.2f %.2f -> %s\n', ...
    tr(3).shrink, tr(3).shrinkFold, ternary(C_RELIABLE, 'C usable as benchmark', ...
    sprintf('C NOT a reliable benchmark (a < %.2f): inspect before concluding on decimation', SHRINK_OK)));
fprintf('\n');

checks = struct();
costB = tr(2).res.sigma ./ tr(3).res.sigma;
checks.U11 = all(costB <= DEC_COST_OK) && C_RELIABLE; checks.decimationCost = costB; checks.CReliable = C_RELIABLE;
fprintf('U1.1 decimation cost sigma(B)/sigma(C): ln beta %.2f, ln R %.2f (limit %.2f)%s   %s\n', costB, DEC_COST_OK, ...
    ternary(C_RELIABLE, '', ', C unreliable'), flagStr(checks.U11));
empB = tr(2).sigmaEmp ./ tr(2).sigmaPredCV; empC = tr(3).sigmaEmp ./ tr(3).sigmaPredCV;
checks.U12 = all([empB; empC] >= EMP_RANGE(1) & [empB; empC] <= EMP_RANGE(2));
fprintf('U1.2 covariance model emp/pred: B %.2f %.2f, C %.2f %.2f (D, which ignores lag 1: %.2f %.2f)   %s\n', ...
    empB, empC, tr(4).sigmaEmp ./ tr(4).sigmaPredCV, flagStr(checks.U12));
checks.U13 = tr(2).res.dBias < D_BIAS_OK && tr(2).res.dBiasPred < D_BIAS_OK;
fprintf('U1.3 production parameter bias: d_bias %.2f (Monte Carlo), %.2f (Welch-expected), limit %.2f   %s\n', ...
    tr(2).res.dBias, tr(2).res.dBiasPred, D_BIAS_OK, flagStr(checks.U13));
fprintf('     dtheta_bias (Monte Carlo) = [%+.4f %+.4f], (Welch-expected) = [%+.4f %+.4f] in [ln beta, ln R]\n', ...
    tr(2).res.dtheta, tr(2).res.dthetaPred);
if STEP_CHECK
    resH2 = infoAnalysis(Jh2, tr(2).Sigma, tr(2).b, tr(2).bPred);
    sigChange = max(abs(resH2.sigma ./ tr(2).res.sigma - 1));
    checks.U14 = sigChange < 0.01; checks.stepDiff = stepDiff;
    fprintf('U1.4 Jacobian step: |J_h - J_h/2| / |J_h| = %.2e, max sigma change %.2e   %s\n', stepDiff, sigChange, flagStr(checks.U14));
end
fprintf('\nIn physical terms (production B): beta known to a factor exp(sigma) = %.2f, R to %.3f (1 sigma)\n\n', ...
    exp(tr(2).res.sigma));

%% ---- 6  Save ------------------------------------------------------------------------------------------------
trSave = rmfield(tr, {'Z'});
prodGrid = struct('omega', wDec, 'estimatorIdx', sets.dec, 'segmentLength', L_PROD, 'overlap', OVERLAP, ...
    'band', BAND_USE, 'f0', f0Dec, 'J', JDec, 'SigmaBlocks', blocksOf(ZDec), 'bias', tr(2).b, ...
    'biasPred', bPredDec, 'F', tr(2).res.F, 'Snn', Snn, 'noiseStd', noiseStd);
out = struct('p0', P0, 'sensors', SENSORS, 'reference', REF, 'theta0', theta0, 'N', N, 'dt', DT, 'nMC', N_MC, ...
    'jonswap', struct('Hs', HS, 'wp', WP, 'gammaJ', GAMMA_J), 'noiseRel', NOISE_REL, 'hStep', H_STEP, ...
    'omegaFull', wFull, 'omegaDec', wDec, 'omega512', w512u, 'f0Full', f0Full, 'JFull', JFull, 'JDecHalfStep', Jh2, ...
    'treatments', trSave, 'checks', checks, 'prod', prodGrid, 'ZFull', ZFull, 'Z512', Z512u, ...
    'jacobianStatus', jstat, 'twinCacheFile', twin.cacheFile);
save([baseName, '.mat'], 'out');
fprintf('Saved: %s.mat\n', baseName);

%% ---- 7  Figure --------------------------------------------------------------------------------------------------
cols = {[120 120 120] / 255, [42 120 214] / 255, [27 175 122] / 255, [235 104 52] / 255};
ink = [0.35 0.35 0.33]; gridInk = [0.85 0.85 0.83];
fig = figure('Color', 'w', 'Position', [60 60 1250 470]);
ax = subplot(1, 2, 1); hold(ax, 'on'); h = []; lab = {};
ang = linspace(0, 2 * pi, 200);
for t = 1:numel(tr)
    E = sqrtm(tr(t).res.covTheta) * [cos(ang); sin(ang)];
    h(end + 1) = plot(ax, E(1, :), E(2, :), '-', 'Color', cols{t}, 'LineWidth', 1.8); %#ok<AGROW>
    lab{end + 1} = [tr(t).key, ' predicted 1\sigma']; %#ok<AGROW>
end
for t = [2 3]
    plot(ax, tr(t).thetaHat(1, :), tr(t).thetaHat(2, :), '.', 'Color', 1 - 0.45 * (1 - cols{t}), 'MarkerSize', 6);
end
h(end + 1) = plot(ax, tr(2).res.dtheta(1), tr(2).res.dtheta(2), 'p', 'Color', cols{2}, 'MarkerFaceColor', cols{2}, 'MarkerSize', 11);
lab{end + 1} = 'B bias (Monte Carlo)';
h(end + 1) = plot(ax, tr(2).res.dthetaPred(1), tr(2).res.dthetaPred(2), 'd', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
lab{end + 1} = 'B bias (Welch-expected)';
set(ax, 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off'); grid(ax, 'on');
set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
xlabel(ax, '\Delta ln \beta'); ylabel(ax, '\Delta ln R');
title(ax, '(a) predicted 1\sigma ellipses, linearised Monte Carlo (dots: B, C)', 'FontWeight', 'normal');
legend(ax, h, lab, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
ax = subplot(1, 2, 2); hold(ax, 'on');
contrib = perFrequencyInfo(JDec, blocksOf(ZDec));           % 2 x nF: diagonal of J_k' Sigma_k^-1 J_k
plot(ax, wDec, contrib(1, :) / sum(contrib(1, :)), '-o', 'Color', cols{2}, 'LineWidth', 1.6, 'MarkerSize', 4, 'MarkerFaceColor', cols{2});
plot(ax, wDec, contrib(2, :) / sum(contrib(2, :)), '-s', 'Color', cols{4}, 'LineWidth', 1.6, 'MarkerSize', 4, 'MarkerFaceColor', cols{4});
set(ax, 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off'); grid(ax, 'on');
set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
xlabel(ax, '\omega (rad/s)'); ylabel(ax, 'share of F_{ii}');
title(ax, '(b) production grid: where the information comes from', 'FontWeight', 'normal');
legend(ax, {'ln \beta', 'ln R'}, 'Location', 'best', 'Box', 'off');
saveFigure(fig, [baseName, '_ellipses.png']);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ================================================================================================
%  Local functions
%% ================================================================================================
function [f0, J, status] = cachedJacobian(cacheDir, theta0, w, sensors, P0, idx, ref, h, frfOpt)
% Central-difference Jacobian of the S1 features in theta, cached on disk by key.
% physical EMM options are part of the key; purely computational ones are not
compOpts = {'useparallel', 'verbose'};
phys = {};
for a = 1:2:numel(frfOpt) - 1
    if ~any(strcmpi(frfOpt{a}, compOpts)), phys = [phys, frfOpt(a:a + 1)]; end %#ok<AGROW>
end
key = struct('theta0', theta0(:).', 'omega', w(:).', 'sensors', sensors, 'P0', P0, 'idx', idx, 'ref', ref, 'h', h, ...
    'frfOptions', {phys});
if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
d = dir(fullfile(cacheDir, 'transJacobian_*.mat'));
for i = 1:numel(d)
    S = load(fullfile(cacheDir, d(i).name), 'key');
    if isfield(S, 'key') && keyMatch(S.key, key)
        S = load(fullfile(cacheDir, d(i).name));
        f0 = S.f0; J = S.J; status = ['loaded ', d(i).name];
        return
    end
end
F = @(th) transmissibilityFeatureVector(th, w, sensors, 'ThetaIdx', idx, 'PBase', P0, 'Reference', ref, ...
    'FRFOptions', frfOpt);
f0 = F(theta0);
J = zeros(numel(f0), numel(theta0));
for k = 1:numel(theta0)
    e = zeros(size(theta0)); e(k) = h;
    J(:, k) = (F(theta0 + e) - F(theta0 - e)) / (2 * h);
end
n = 1; file = fullfile(cacheDir, sprintf('transJacobian_%s_%02d.mat', datestr(now, 'yyyymmdd_HHMMSS'), n)); %#ok<TNOW1,DATST>
while exist(file, 'file'), n = n + 1; file = fullfile(cacheDir, sprintf('transJacobian_%s_%02d.mat', datestr(now, 'yyyymmdd_HHMMSS'), n)); end %#ok<TNOW1,DATST>
save(file, 'key', 'f0', 'J', '-v7');
status = sprintf('computed (%d EMM solves), cached', numel(w) * (1 + 2 * numel(theta0)));
end

function tf = keyMatch(a, b)
tf = false;
try
    if ~isequal(size(a.omega), size(b.omega)) || any(abs(a.omega - b.omega) > 1e-9), return; end
    if any(abs(a.theta0 - b.theta0) > 1e-12) || any(abs(a.P0 - b.P0) > 1e-12 * abs(b.P0)), return; end
    if ~isequal(size(a.sensors), size(b.sensors)) || any(abs(a.sensors(:) - b.sensors(:)) > 1e-12), return; end
    if ~isfield(a, 'frfOptions') || ~isequal(a.frfOptions, b.frfOptions), return; end
    tf = isequal(a.idx, b.idx) && isequal(a.ref, b.ref) && a.h == b.h;
catch
    tf = false;
end
end

function B = blocksOf(Z)
% Per-bin p x p sample covariances of draws Z (p x nF x nR).
[p, nF, ~] = size(Z);
B = zeros(p, p, nF);
for k = 1:nF
    B(:, :, k) = cov(squeeze(Z(:, k, :)).');
end
end

function [S, a] = buildSigma(Z, banded)
% Block-diagonal covariance, or block tridiagonal with lag-1 cross-covariances
% (shrunk towards block diagonal until positive definite).
[p, nF, ~] = size(Z);
Bd = zeros(p * nF);
blocks = blocksOf(Z);
for k = 1:nF
    ix = (k - 1) * p + (1:p);
    Bd(ix, ix) = blocks(:, :, k);
end
a = 1;
S = Bd;
if ~banded, return; end
Off = zeros(p * nF);
for k = 1:nF - 1
    X = [squeeze(Z(:, k, :)); squeeze(Z(:, k + 1, :))];
    C = cov(X.');
    ix = (k - 1) * p + (1:p); iy = k * p + (1:p);
    Off(ix, iy) = C(1:p, p + 1:end);
    Off(iy, ix) = C(1:p, p + 1:end).';
end
for a = 1:-0.05:0
    S = Bd + a * Off;
    [~, flag] = chol(S);
    if flag == 0, return; end
end
end

function q = infoAnalysis(J, S, b, bPred)
L = chol(S, 'lower');
Jw = L \ J;
F = Jw.' * Jw;
C = inv(F);
q.F = F; q.covTheta = C;
q.sigma = sqrt(diag(C));
q.corr = C(1, 2) / prod(q.sigma);
q.cond = cond(F);
q.dtheta = F \ (Jw.' * (L \ b));
q.ratio = q.dtheta ./ q.sigma;
q.dBias = sqrt(q.dtheta.' * F * q.dtheta);
q.dthetaPred = F \ (Jw.' * (L \ bPred));
q.dBiasPred = sqrt(q.dthetaPred.' * F * q.dthetaPred);
end

function c = perFrequencyInfo(J, blocks)
[p, ~, nF] = size(blocks);
c = zeros(2, nF);
for k = 1:nF
    Jk = J((k - 1) * p + (1:p), :);
    Fk = Jk.' * (blocks(:, :, k) \ Jk);
    c(:, k) = diag(Fk);
end
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

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function s = flagStr(ok)
if ok, s = 'PASS'; else, s = 'FLAG'; end
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end