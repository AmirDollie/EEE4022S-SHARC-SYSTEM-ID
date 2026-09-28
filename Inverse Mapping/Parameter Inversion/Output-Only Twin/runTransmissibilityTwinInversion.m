%% runTransmissibilityTwinInversion.m
% U2 / Gate G3: close the inverse loop. Finite stochastic twin records at
% truth points away from p0 are turned into corrected transmissibility
% features (production setting F) and inverted for theta = [ln beta; ln R]
% (gamma fixed), starting from p0. The question:
%
%   Do finite-record output-only transmissibilities recover (beta, R) with
%   errors consistent with the predicted covariance F^-1?
%
% PRODUCTION SETTING F (frozen after runL2048BiasVarianceCheck)
%   N = 65536, dt = 0.1 s, L = 2048, 50% overlap, reference-noise correction on,
%   every second L = 2048 bin in [3.178, 8.301] rad/s: 84 frequencies, 504 features,
%   Sigma_T = blkdiag(Sigma_k), model side the point ratio T_j = H_j / H_r.
%
% TRUTH POINTS (gamma = gamma0; sensors stay at their p0 positions, all inside
% R >= 0.97 R0):  A = (1.2 beta0, 0.97 R0),  B = (0.8 beta0, 1.03 R0).
% Sea state and ABSOLUTE sensor noise are the same at every point (noise std fixed
% from p0, a physical sensor property).
%
% MODEL EVALUATOR (nonlinear solver): the EMM on a node grid MODEL_SPACING apart
% (synthesiseTwinRecords node cache), spline-interpolated to the 84 frequencies,
% memoised by theta. At each truth point it is checked against direct EMM
% features (S1): the interpolation error is projected into parameter space,
% d_interp = |dtheta_interp|_F, and must stay below D_INTERP_OK.
%
% ENSEMBLES (disjoint seeds per truth point)
%   calibration  N_CAL records  -> Sigma(theta_true) (oracle), Monte Carlo bias
%   validation   N_VAL records  -> linearised coverage tests G3a, G3b, G3c
%   nonlinear    N_NL records   -> full nonlinear inversions from p0
%   plus N_CAL calibration records at p0 for G3c.
%
% TESTS (linearised: theta_hat = theta_true + C J' Sigma^-1 (z - f_true), J = J(theta_true))
%   G3a  oracle Sigma(theta_true): Cov(theta_hat) vs F^-1, mean offset vs the
%        Welch-expected bias dtheta_W (not only vs zero), coverage of the F^-1
%        ellipses: q = (theta_hat - theta_true)' F (theta_hat - theta_true),
%        E[q] = 2 + d_W^2, P(q <= 5.991) ~ 0.95 when unbiased
%   G3b  record-derived Sigma: each validation record estimates its own per-bin
%        6 x 6 covariance by a BLOCK jackknife over its Welch segments (blocks of
%        JK_BLOCK consecutive segments, delete one block), inverts with it and
%        reports its own C_i = (J' Sigma_i^-1 J)^-1; q_i uses C_i. Reported raw and
%        with the Hartlap precision factor (B - p - 2)/(B - 1) (B blocks, p = 6),
%        which removes the known optimism of an inverted sample covariance.
%        A 6 x 6 covariance from about 30 blocks is noisy, so a POOLED variant
%        averages each bin's jackknife covariance with its POOL_HALFWIDTH
%        neighbours on each side (the covariance varies smoothly with omega),
%        with the per-bin Hartlap factor for the pooled degrees of freedom.
%        This is the deployable uncertainty: no knowledge of S_u or theta.
%   G3c  Sigma(theta0) from p0 records used at the truth point (sensitivity only)
% NONLINEAR: safeguarded modified Gauss-Newton with the oracle Sigma, from p0;
%   J reused, recomputed when a step fails to reduce the objective even after
%   halving, when an accepted step exceeds STEP_SIGMA predicted sigmas, or after
%   3 iterations without convergence. Converged when every step component is
%   below TOL_SIGMA predicted sigmas. Compared with the linearised estimate from
%   the same record (agreement means the linear theory applies there).
%
% GATE G3 (per truth point; PASS / FLAG)
%   G3.1 nonlinear: converged, |theta_NL - theta_true| < 3 sigma, |theta_NL - theta_lin| < 0.2 sigma
%   G3.2 G3a: sigma_emp / sigma_pred in [0.85, 1.20], 95% coverage in [0.90, 0.99]
%   G3.3 G3a: mean offset within 3 standard errors of dtheta_W
%   G3.4 G3b (pooled, Hartlap): 95% coverage >= 0.90 (raw variants reported alongside)
%   plus the interpolation check d_interp < D_INTERP_OK
%
% RUNTIME (MATLAB): per truth point about 60 EMM solves (twin) + 84 (direct S1) +
% 116 (J at the truth) + about 300-600 per nonlinear inversion; records about
% 3 min in total. Roughly 30-40 min for two truth points and one nonlinear
% record each. Twins are cached in Results/twinCache.
%
% OUTPUT (in Inverse Mapping/Parameter Inversion/Results/)
%   transmissibilityTwinInversion_<stamp>.mat / .txt / _g3.png
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
P0 = [4.6985e-5, 1.4548e-3, 0.3830]; R0 = P0(3);
SENSORS = [0.3*R0, 0; 0.3*R0, pi; 0.5*R0, pi; 0.9*R0, 0];
REF = 4; THETA_IDX = [1 3];
TRUTH = [1.2 * P0(1), P0(2), 0.97 * R0;                     % A
         0.8 * P0(1), P0(2), 1.03 * R0];                    % B
TRUTH_NAMES = {'A (1.2 beta0, 0.97 R0)', 'B (0.8 beta0, 1.03 R0)'};
DT = 0.1; N = 65536; L = 2048; OVERLAP = 0.5;
BAND_EST = [3 8.5]; BAND_USE = [3.178 8.301];
HS = 0.05; WP = 5.0; GAMMA_J = 3.3; NOISE_REL = 0.10;
N_CAL = 200; N_VAL = 200; N_NL = 1; % CHANGE N_NL = 1 FOR ACTUAL RUN!!!
JK_BLOCK = 2;                                               % Welch segments per jackknife block
POOL_HALFWIDTH = 2;                                         % G3b pooled: +- bins averaged
MODEL_SPACING = 0.2;                                        % EMM node spacing of the solver's model (rad/s)
H_STEP = 1e-3;
MAX_IT = 20; STEP_SIGMA = 3; TOL_SIGMA = 0.01;
D_INTERP_OK = 0.1;
SEED_CAL = 20000; SEED_VAL = 30000; SEED_NL = 40000; SEED_P0 = 50000;

%% ---- Output ----------------------------------------------------------------------
resultsDir = fullfile(thisDir, '..', 'Results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
cacheDir = fullfile(resultsDir, 'twinCache');
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['transmissibilityTwinInversion_', stamp]);
diary([baseName, '.txt']);
try

nS = size(SENSORS, 1); m = nS - 1; p = 2 * m; others = setdiff(1:nS, REF);
theta0 = log(P0(THETA_IDX)).';
specJ = {'Spectrum', 'jonswap', 'Hs', HS, 'PeakFrequency', WP, 'PeakEnhancement', GAMMA_J};
chi95 = -2 * log(0.05); chi68 = -2 * log(1 - 0.6827);       % chi2(2) quantiles
fprintf('runTransmissibilityTwinInversion  (%s)\n', stamp);
fprintf('setting F: N = %d, L = %d, every 2nd bin in [%.3f %.3f]; N_CAL %d, N_VAL %d, N_NL %d; jackknife blocks of %d segments\n\n', ...
    N, L, BAND_USE, N_CAL, N_VAL, N_NL, JK_BLOCK);

%% ---- 1  p0 twin, noise level, frequency set -------------------------------------------------------
twin0 = synthesiseTwinRecords(struct('p', P0, 'sensors', SENSORS, 'nodeSpacing', 0.1), 'CacheDir', cacheDir);
[Y0, ~, info0] = synthesiseTwinRecords(twin0, N, DT, 'Seed', SEED_P0, specJ{:});
noiseStd = NOISE_REL * info0.signalStd(REF);
Snn = noiseStd^2 * DT / pi;
eOpt = {'Reference', REF, 'SegmentLength', L, 'Overlap', OVERLAP, 'Band', BAND_EST, 'NoisePSD', Snn};
e0 = extractTransmissibility(Y0, DT, eOpt{:});
keep = mod(e0.binIndex, 2) == 0 & e0.omega >= BAND_USE(1) & e0.omega <= BAND_USE(2);
kBins = e0.binIndex(keep); omega = e0.omega(keep); nF = numel(omega);
% own segment-level Welch must reproduce extractTransmissibility exactly
P = welchSegments(Y0, DT, L, OVERLAP, kBins, REF, others);
zChk = zFromSegments(P, 1:P.K, Snn);
zRef = reshape(stackTransmissibility(e0.T(:, keep)), p, nF);
assert(max(abs(zChk(:) - zRef(:))) < 1e-10 * max(abs(zRef(:))), 'segment-level Welch does not reproduce extractTransmissibility');
fprintf('noise std %.4g m/s^2 (S_nn %.3e); %d frequencies %.3f to %.3f rad/s; %d Welch segments, %d jackknife blocks\n', ...
    noiseStd, Snn, nF, omega([1 end]), P.K, floor(P.K / JK_BLOCK));
fprintf('segment-level Welch reproduces extractTransmissibility (verified)\n');

% calibration at p0 (for G3c)
Zp0 = runEnsemble(twin0, N_CAL, SEED_P0 + 1, N, DT, specJ, noiseStd, L, OVERLAP, kBins, REF, others, Snn);
blocksP0 = blockCov(Zp0);

% model evaluator and J at theta0 (shared by all nonlinear inversions)
modelBand = [max(twin0.band(1), omega(1) - MODEL_SPACING), min(twin0.band(2), omega(end) + MODEL_SPACING)];
clear modelFeatures                                          % reset the memo of the model evaluator
modelFun = @(th) modelFeatures(th, P0, THETA_IDX, SENSORS, omega, REF, others, MODEL_SPACING, modelBand);
t0 = tic;
J0 = jacobianCD(modelFun, theta0, H_STEP);
fprintf('J(theta0) through the node-interpolated model: %.0f s\n\n', toc(t0));

%% ---- 2  Truth points ------------------------------------------------------------------------------------
nT = size(TRUTH, 1);
res = struct([]);
for it = 1:nT
    pT = TRUTH(it, :);
    thT = log(pT(THETA_IDX)).';
    fprintf('=== truth %s: theta_true - theta0 = [%+.4f %+.4f] ===\n', TRUTH_NAMES{it}, thT - theta0);
    t1 = tic;
    twinT = synthesiseTwinRecords(struct('p', pT, 'sensors', SENSORS, 'nodeSpacing', 0.1), 'CacheDir', cacheDir);
    [~, ~, infoT] = synthesiseTwinRecords(twinT, N, DT, 'Seed', 1, specJ{:});
    fTrue = reshape(transmissibilityFeatureVector(thT, omega, SENSORS, 'ThetaIdx', THETA_IDX, 'PBase', P0, ...
        'Reference', REF), p, nF);                                     % direct EMM (S1)
    JT = jacobianCD(modelFun, thT, H_STEP);
    fprintf('  twin %s (interp check T %.1e); direct S1 and J(theta_true): %.0f s\n', twinT.cacheStatus, ...
        twinT.interpCheck.maxRelErrT, toc(t1));

    % calibration: oracle Sigma(theta_true) and Monte Carlo bias
    Zc = runEnsemble(twinT, N_CAL, SEED_CAL + 1000 * it, N, DT, specJ, noiseStd, L, OVERLAP, kBins, REF, others, Snn);
    blocksT = blockCov(Zc);
    FT = infoMatrix(JT, blocksT);
    CT = inv(FT); sigT = sqrt(diag(CT));
    bMC = mean(Zc, 3) - fTrue;
    bW = welchExpectedZ(omega, infoT, L, DT, REF, others, 0) - fTrue;
    dthW = CT * gradTerm(JT, blocksT, bW); dthMC = CT * gradTerm(JT, blocksT, bMC);
    % interpolation check of the solver's model, in parameter space
    fInterp = reshape(modelFun(thT), p, nF);
    dthI = CT * gradTerm(JT, blocksT, fInterp - fTrue);
    dInterp = sqrt(dthI.' * FT * dthI);
    relI = max(abs(fInterp(:) - fTrue(:))) / max(abs(fTrue(:)));
    fprintf('  sigma(theta) = [%.5f %.5f]; Welch-expected bias [%+.5f %+.5f] (d %.2f), Monte Carlo [%+.5f %+.5f]\n', ...
        sigT, dthW, sqrt(dthW.' * FT * dthW), dthMC);
    fprintf('  model interpolation (spacing %.2f): max rel error %.1e, d_interp %.3f  %s\n', MODEL_SPACING, relI, dInterp, ...
        flagStr(dInterp < D_INTERP_OK));

    % validation: G3a / G3b / G3c
    [Zv, Sjk, nBlk] = runEnsemble(twinT, N_VAL, SEED_VAL + 1000 * it, N, DT, specJ, noiseStd, L, OVERLAP, kBins, REF, ...
        others, Snn, JK_BLOCK);
    hart = (nBlk - p - 2) / (nBlk - 1);
    thA = zeros(2, N_VAL); thC = zeros(2, N_VAL); qA = zeros(1, N_VAL); qC = zeros(1, N_VAL);
    Fp0 = infoMatrix(JT, blocksP0); Cp0 = inv(Fp0);
    vNames = {'raw', 'raw Hartlap', 'pooled', 'pooled Hartlap'};
    thBv = zeros(2, N_VAL, 4); qBv = zeros(N_VAL, 4); sigBv = zeros(2, N_VAL, 4);
    for i = 1:N_VAL
        r = Zv(:, :, i) - fTrue;
        thA(:, i) = CT * gradTerm(JT, blocksT, r);
        qA(i) = thA(:, i).' * FT * thA(:, i);
        thC(:, i) = Cp0 * gradTerm(JT, blocksP0, r);
        qC(i) = thC(:, i).' * Fp0 * thC(:, i);
        raw = Sjk(:, :, :, i);
        [pooled, hPool] = poolBlocks(raw, POOL_HALFWIDTH, nBlk, p);
        variants = {raw, raw / hart, pooled, pooled ./ reshape(hPool, 1, 1, [])};   % Sigma / h = Hartlap precision
        for v = 1:4
            Fi = infoMatrix(JT, variants{v});
            thBv(:, i, v) = Fi \ gradTerm(JT, variants{v}, r);
            qBv(i, v) = thBv(:, i, v).' * Fi * thBv(:, i, v);
            sigBv(:, i, v) = sqrt(diag(inv(Fi)));
        end
    end
    covBv = [mean(qBv <= chi68, 1); mean(qBv <= chi95, 1)];                % 2 x 4
    sA = std(thA, 0, 2); mA = mean(thA, 2); seA = sA / sqrt(N_VAL);
    covA = [mean(qA <= chi68), mean(qA <= chi95)];
    covC = [mean(qC <= chi68), mean(qC <= chi95)];
    fprintf('  G3a oracle:     emp/pred [%.2f %.2f]; mean offset [%+.5f %+.5f] +- [%.5f %.5f] (Welch-expected [%+.5f %+.5f]);\n', ...
        sA ./ sigT, mA, seA, dthW);
    fprintf('                  mean q %.2f (2 + d_W^2 = %.2f); coverage 68%% %.3f, 95%% %.3f\n', mean(qA), ...
        2 + dthW.' * FT * dthW, covA);
    fprintf('  G3b jackknife:  %d blocks of %d segments (Hartlap %.3f); pooled over +-%d bins (interior Hartlap %.3f)\n', ...
        nBlk, JK_BLOCK, hart, POOL_HALFWIDTH, ((2 * POOL_HALFWIDTH + 1) * (nBlk - 1) + 1 - p - 2) / ((2 * POOL_HALFWIDTH + 1) * (nBlk - 1)));
    for v = 1:4
        fprintf('     %-15s coverage 68%% %.3f 95%% %.3f; scatter / median reported sigma [%.2f %.2f]; scatter / oracle scatter [%.2f %.2f]\n', ...
            vNames{v}, covBv(:, v), std(thBv(:, :, v), 0, 2) ./ median(sigBv(:, :, v), 2), std(thBv(:, :, v), 0, 2) ./ sA);
    end
    fprintf('  G3c Sigma(p0):  emp/claimed [%.2f %.2f]; coverage 68%% %.3f, 95%% %.3f\n', std(thC, 0, 2) ./ sqrt(diag(Cp0)), covC);

    % nonlinear inversions from p0 (oracle Sigma)
    nl = struct([]);
    for j = 1:N_NL
        [Yn] = synthesiseTwinRecords(twinT, N, DT, 'Seed', SEED_NL + 1000 * it + j, specJ{:}, 'NoiseStd', noiseStd);
        Pn = welchSegments(Yn, DT, L, OVERLAP, kBins, REF, others);
        zn = zFromSegments(Pn, 1:Pn.K, Snn);
        thLin = thT + CT * gradTerm(JT, blocksT, zn - fTrue);
        t2 = tic;
        g = gaussNewton(modelFun, zn, blocksT, theta0, J0, sigT, H_STEP, MAX_IT, STEP_SIGMA, TOL_SIGMA);
        nl(j).theta = g.theta; nl(j).thetaLin = thLin; nl(j).it = g.it; nl(j).jRecomputes = g.jRecomputes;
        nl(j).converged = g.converged; nl(j).chi2 = g.obj; nl(j).dof = numel(zn) - 2; nl(j).seconds = toc(t2);
        nl(j).errSigma = (g.theta - thT) ./ sigT; nl(j).linGapSigma = (g.theta - thLin) ./ sigT;
        nl(j).evals = g.evals; nl(j).history = g.history;
        fprintf('  nonlinear %d: %s in %d it (%d J recomputes, %d model evals, %.0f s); error [%+.2f %+.2f] sigma;\n', j, ...
            ternary(g.converged, 'converged', 'NOT converged'), g.it, g.jRecomputes, g.evals, nl(j).seconds, nl(j).errSigma);
        fprintf('               vs linearised [%+.3f %+.3f] sigma; chi2/dof %.3f; beta %+.3f%%, R %+.4f%% from truth\n', ...
            nl(j).linGapSigma, g.obj / nl(j).dof, 100 * (exp(g.theta(1) - thT(1)) - 1), 100 * (exp(g.theta(2) - thT(2)) - 1));
    end

    % gate
    gate.interp = dInterp < D_INTERP_OK;
    if isempty(nl)
        gate.g31 = NaN;                                                   % N_NL = 0: not tested
    else
        eS = [nl.errSigma]; lG = [nl.linGapSigma];
        gate.g31 = all([nl.converged]) && all(abs(eS(:)) < 3) && all(abs(lG(:)) < 0.2);
    end
    gate.g32 = all(sA ./ sigT >= 0.85 & sA ./ sigT <= 1.20) && covA(2) >= 0.90 && covA(2) <= 0.99;
    gate.g33 = all(abs(mA - dthW) < 3 * seA);
    gate.g34 = covBv(2, 4) >= 0.90;
    fprintf('  GATE: interp %s | G3.1 nonlinear %s | G3.2 oracle coverage %s | G3.3 offset = Welch bias %s | G3.4 jackknife %s\n\n', ...
        flagStr(gate.interp), flagStr(gate.g31), flagStr(gate.g32), flagStr(gate.g33), flagStr(gate.g34));

    res(it).name = TRUTH_NAMES{it}; res(it).p = pT; res(it).thetaTrue = thT; res(it).fTrue = fTrue; res(it).J = JT;
    res(it).blocks = blocksT; res(it).F = FT; res(it).sigma = sigT; res(it).dthetaWelch = dthW; res(it).dthetaMC = dthMC;
    res(it).dInterp = dInterp; res(it).relInterp = relI;
    res(it).G3a = struct('theta', thA, 'sigmaEmp', sA, 'mean', mA, 'se', seA, 'q', qA, 'coverage', covA);
    res(it).G3b = struct('variants', {vNames}, 'theta', thBv, 'q', qBv, 'coverage', covBv, 'sigmaReported', sigBv, ...
        'nBlocks', nBlk, 'hartlap', hart, 'poolHalfwidth', POOL_HALFWIDTH);
    res(it).G3c = struct('theta', thC, 'q', qC, 'coverage', covC, 'sigmaClaimed', sqrt(diag(Cp0)));
    res(it).nonlinear = nl; res(it).gate = gate;
end

%% ---- 3  Save ------------------------------------------------------------------------------------------------
out = struct('p0', P0, 'sensors', SENSORS, 'reference', REF, 'theta0', theta0, 'omega', omega, 'kBins', kBins, ...
    'N', N, 'dt', DT, 'L', L, 'overlap', OVERLAP, 'noiseStd', noiseStd, 'Snn', Snn, 'modelSpacing', MODEL_SPACING, ...
    'jkBlock', JK_BLOCK, 'nCal', N_CAL, 'nVal', N_VAL, 'nNL', N_NL, 'J0', J0, 'blocksP0', blocksP0, 'truth', res);
save([baseName, '.mat'], 'out');
fprintf('Saved: %s.mat\n', baseName);

%% ---- 4  Figure --------------------------------------------------------------------------------------------------
colA = [42 120 214] / 255; colB = [235 104 52] / 255; colC = [120 120 120] / 255; ink = [0.35 0.35 0.33];
fig = figure('Color', 'w', 'Position', [60 60 1250 420 * nT]);
ang = linspace(0, 2 * pi, 200); qq = linspace(0, 12, 200);
for it = 1:nT
    q = res(it);
    ax = subplot(nT, 2, 2 * it - 1); hold(ax, 'on');
    plot(ax, q.G3a.theta(1, :), q.G3a.theta(2, :), '.', 'Color', 1 - 0.5 * (1 - colA), 'MarkerSize', 7);
    E = sqrtm(inv(q.F)) * [cos(ang); sin(ang)];
    h1 = plot(ax, E(1, :), E(2, :), '-', 'Color', colA, 'LineWidth', 1.8);
    h2 = plot(ax, 2.448 * E(1, :), 2.448 * E(2, :), '--', 'Color', colA, 'LineWidth', 1.2);
    h3 = plot(ax, q.dthetaWelch(1), q.dthetaWelch(2), 'd', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
    hs = [h1 h2 h3]; lab = {'F^{-1} 1\sigma', '95% ellipse', 'Welch-expected bias'};
    for j = 1:numel(q.nonlinear)
        d = q.nonlinear(j).theta - q.thetaTrue;
        hn = plot(ax, d(1), d(2), 'p', 'Color', colB, 'MarkerFaceColor', colB, 'MarkerSize', 12);
        if j == 1, hs(end + 1) = hn; lab{end + 1} = 'nonlinear estimate'; end %#ok<AGROW,SAGROW>
    end
    set(ax, 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off'); grid(ax, 'on');
    xlabel(ax, '\Delta ln \beta (from truth)'); ylabel(ax, '\Delta ln R (from truth)');
    title(ax, sprintf('(%c) truth %s: G3a linearised validation', 'a' + 2 * (it - 1), q.name), 'FontWeight', 'normal');
    legend(ax, hs, lab, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
    ax = subplot(nT, 2, 2 * it); hold(ax, 'on');
    plot(ax, qq, 1 - exp(-qq / 2), '-', 'Color', 'k', 'LineWidth', 1.2);
    plot(ax, sort(q.G3a.q), (1:N_VAL) / N_VAL, '-', 'Color', colA, 'LineWidth', 1.8);
    plot(ax, sort(q.G3b.q(:, 2)), (1:N_VAL) / N_VAL, ':', 'Color', colB, 'LineWidth', 1.8);
    plot(ax, sort(q.G3b.q(:, 4)), (1:N_VAL) / N_VAL, '-', 'Color', colB, 'LineWidth', 1.8);
    plot(ax, sort(q.G3c.q), (1:N_VAL) / N_VAL, '-', 'Color', colC, 'LineWidth', 1.2);
    xlim(ax, [0 12]); ylim(ax, [0 1]);
    set(ax, 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off'); grid(ax, 'on');
    xlabel(ax, 'q = \Delta\theta^T C^{-1} \Delta\theta'); ylabel(ax, 'empirical CDF');
    title(ax, sprintf('(%c) calibration of reported uncertainty', 'b' + 2 * (it - 1)), 'FontWeight', 'normal');
    legend(ax, {'\chi^2_2 (ideal)', 'G3a oracle', 'G3b jackknife (Hartlap)', 'G3b pooled jackknife (Hartlap)', 'G3c \Sigma(p_0)'}, ...
        'Location', 'southeast', 'Box', 'off', 'FontSize', 8);
end
saveFigure(fig, [baseName, '_g3.png']);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ================================================================================================
%  Local functions
%% ================================================================================================
function P = welchSegments(Y, dt, L, ov, kBins, ref, others)
% Per-segment scaled cross-spectra, identical conventions to extractTransmissibility.
n = size(Y, 2);
w = 0.5 * (1 - cos(2 * pi * (0:L - 1).' / L));
U = sum(w.^2);
hop = max(1, round(L * (1 - ov)));
starts = 1:hop:(n - L + 1);
K = numel(starts);
scale = 2 * dt / (U * 2 * pi);
nF = numel(kBins);
P.Sjr = complex(zeros(numel(others), nF, K)); P.Srr = zeros(1, nF, K); P.K = K;
for s = 1:K
    seg = Y(:, starts(s):starts(s) + L - 1);
    seg = seg - mean(seg, 2);
    X = fft(seg .* w.', [], 2);
    X = X(:, kBins + 1);
    P.Sjr(:, :, s) = scale * X(others, :) .* conj(X(ref, :));
    P.Srr(1, :, s) = scale * abs(X(ref, :)).^2;
end
end

function z = zFromSegments(P, segIdx, Snn)
% Corrected transmissibility features (p x nF) from a subset of segments.
Sjr = mean(P.Sjr(:, :, segIdx), 3); Srr = mean(P.Srr(:, :, segIdx), 3);
den = Srr - Snn;
T = Sjr ./ den;
T(:, ~(den > 0.1 * Srr)) = complex(NaN, NaN);
z = reshape(stackTransmissibility(T), 2 * size(T, 1), []);
end

function [blocks, nBlk] = jackknifeCov(P, blockSize, Snn)
% Delete-one-block jackknife covariance of z per bin (p x p x nF).
nBlk = floor(P.K / blockSize);
edges = [(0:nBlk - 1) * blockSize, P.K];                 % last block absorbs the remainder
zb = [];
for b = 1:nBlk
    keepSeg = setdiff(1:P.K, edges(b) + 1:edges(b + 1));
    zb(:, :, b) = zFromSegments(P, keepSeg, Snn); %#ok<AGROW>
end
[p, nF, ~] = size(zb);
blocks = zeros(p, p, nF);
for k = 1:nF
    X = reshape(zb(:, k, :), p, nBlk);
    X = X - mean(X, 2);
    blocks(:, :, k) = (nBlk - 1) / nBlk * (X * X.');
end
end

function [Z, Sjk, nBlk] = runEnsemble(twin, nR, seed0, N, dt, specJ, noiseStd, L, ov, kBins, ref, others, Snn, jkBlock)
% Features of nR records (p x nF x nR); optionally the per-record jackknife covariance.
wantJK = nargin > 13;
Z = []; Sjk = []; nBlk = NaN;
for r = 1:nR
    Y = synthesiseTwinRecords(twin, N, dt, 'Seed', seed0 + r, specJ{:}, 'NoiseStd', noiseStd);
    P = welchSegments(Y, dt, L, ov, kBins, ref, others);
    z = zFromSegments(P, 1:P.K, Snn);
    if any(isnan(z(:)))
        error('runTransmissibilityTwinInversion:invalidBin', 'Record %d has an invalid corrected bin.', r);
    end
    Z(:, :, r) = z; %#ok<AGROW>
    if wantJK
        [Sjk(:, :, :, r), nBlk] = jackknifeCov(P, jkBlock, Snn); %#ok<AGROW>
    end
end
end

function [P, h] = poolBlocks(B, w, nBlk, p)
% Average each bin's covariance with up to w neighbours per side; per-bin
% Hartlap factor for the pooled degrees of freedom (nBins * (nBlk - 1)).
nF = size(B, 3);
P = zeros(size(B)); h = zeros(1, nF);
for k = 1:nF
    idx = max(1, k - w):min(nF, k + w);
    P(:, :, k) = mean(B(:, :, idx), 3);
    dof = numel(idx) * (nBlk - 1);
    h(k) = (dof - p - 1) / dof;
end
end

function B = blockCov(Z)
[p, nF, nR] = size(Z);
B = zeros(p, p, nF);
for k = 1:nF
    B(:, :, k) = cov(reshape(Z(:, k, :), p, nR).');
end
end

function F = infoMatrix(J, blocks)
[p, ~, nF] = size(blocks);
F = zeros(size(J, 2));
for k = 1:nF
    Jk = J((k - 1) * p + (1:p), :);
    F = F + Jk.' * (blocks(:, :, k) \ Jk);
end
end

function g = gradTerm(J, blocks, r)
% J' Sigma^-1 r for blkdiag Sigma; r is p x nF or a stacked vector.
[p, ~, nF] = size(blocks);
r = reshape(r, p, nF);
g = zeros(size(J, 2), 1);
for k = 1:nF
    Jk = J((k - 1) * p + (1:p), :);
    g = g + Jk.' * (blocks(:, :, k) \ r(:, k));
end
end

function v = objective(blocks, r)
[p, ~, nF] = size(blocks);
r = reshape(r, p, nF);
v = 0;
for k = 1:nF
    v = v + r(:, k).' * (blocks(:, :, k) \ r(:, k));
end
end

function f = modelFeatures(theta, P0, idx, sensors, omega, ref, others, spacing, band)
% Node-interpolated EMM transmissibility features (stacked), memoised by theta
% (the key also carries the grid and node settings, so a changed setup never hits).
persistent keys vals
key = sprintf('%.12e_', theta, P0, spacing, band, omega(1), omega(end), numel(omega));
if ~isempty(keys)
    hit = find(strcmp(keys, key), 1);
    if ~isempty(hit), f = vals{hit}; return; end
end
p = P0; p(idx) = exp(theta(:).');
tw = synthesiseTwinRecords(struct('p', p, 'sensors', sensors, 'nodeSpacing', spacing, 'band', band, 'checkPoints', 0));
H = tw.Hacc(omega);
f = stackTransmissibility(H(others, :) ./ H(ref, :));
keys{end + 1} = key; vals{end + 1} = f;
end

function J = jacobianCD(fun, theta, h)
f1 = fun(theta);
J = zeros(numel(f1), numel(theta));
for k = 1:numel(theta)
    e = zeros(size(theta)); e(k) = h;
    J(:, k) = (fun(theta + e) - fun(theta - e)) / (2 * h);
end
end

function g = gaussNewton(fun, z, blocks, theta, J, sigma, h, maxIt, stepSigma, tolSigma)
% Safeguarded modified Gauss-Newton (J reused, recomputed when needed).
evals = 0;
f = fun(theta); evals = evals + 1;
r = z(:) - f; obj = objective(blocks, r);
g.jRecomputes = 0; g.converged = false; g.history = [theta; obj];
sinceJ = 0; jFresh = true;                                  % J(theta0) was computed at the start point
for it = 1:maxIt
    F = infoMatrix(J, blocks);
    step = F \ gradTerm(J, blocks, r);
    accepted = false;
    for a = [1 0.5 0.25 0.125]
        thNew = theta + a * step;
        fNew = fun(thNew); evals = evals + 1;
        rNew = z(:) - fNew; objNew = objective(blocks, rNew);
        if objNew < obj, accepted = true; break; end
    end
    if ~accepted
        if jFresh
            break                                            % J is current and still no descent: at the minimum
        end
        J = jacobianCD(fun, theta, h); evals = evals + 4; g.jRecomputes = g.jRecomputes + 1; sinceJ = 0; jFresh = true;
        continue
    end
    theta = thNew; r = rNew; obj = objNew; sinceJ = sinceJ + 1; jFresh = false;
    g.history(:, end + 1) = [theta; obj];
    if all(abs(a * step) ./ sigma < tolSigma)
        g.converged = true; break
    end
    if any(abs(a * step) ./ sigma > stepSigma) || sinceJ >= 3
        J = jacobianCD(fun, theta, h); evals = evals + 4; g.jRecomputes = g.jRecomputes + 1; sinceJ = 0; jFresh = true;
    end
end
g.theta = theta; g.obj = obj; g.it = it; g.evals = evals;
end

function z = welchExpectedZ(wBins, info, L, dt, ref, others, SnnDen)
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
if isnan(ok), s = 'n/a'; elseif ok, s = 'PASS'; else, s = 'FLAG'; end
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end