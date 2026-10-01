%% testAnalyticEngine.m
% E1 validation gate. The analytic engine must reproduce the trusted Monte Carlo results of U1
% (transmissibilityCovariance_20260927_124013), the L2048 check (L2048BiasVarianceCheck_20260927_130816)
% and G3 before any study may use it.
%
% QUICK part (cache only, no EMM solve, about a minute): covariance algebra against a direct complex
% Gaussian simulation; Sigma_k, sigma and d_W at settings F (L 2048), B (L 1024) and A (L 512) against
% the cached Monte Carlo, using the cached p0 twin FRF and the cached Jacobians; covariance STRUCTURE
% (every element, standardised by its Monte Carlo sampling error); invariances; and evaluateScenario
% end to end on an analytic test field (gamma never estimated, rejected bins removed from F,
% a non-identifiable field returns a clean 'Not identifiable').
% FULL part (EMM solves through evaluateScenario, the production chain; set FULL_POINTS):
%   {'p0'}             today's E1 exit: F(p0) end to end              (~140 EMM solves at 0.2 rad/s)
%   {'p0', 'A', 'B'}   Friday's gate: adds the oracle sigma at truths A and B
% EMM solves are cached in Results/twinCache, so a rerun is fast.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

FULL_POINTS = {'p0', 'A', 'B'};  % {} quick only; {'p0'} today; {'p0', 'A', 'B'} for the gate

thisDir = fileparts(mfilename('fullpath'));
drRoot = fullfile(thisDir, '..');
piRoot = fullfile(drRoot, '..');
root = fullfile(piRoot, '..', '..');
addpath(drRoot); addpath(fullfile(drRoot, 'Engine'));
addpath(fullfile(root, 'Transmissibility')); addpath(fullfile(root, 'Inverse Mapping'));
addpath(fullfile(root, 'Forward Model')); addpath(fullfile(root, 'Forward Model', 'Animation')); addpath(fullfile(root, 'SSI'));
resDir = fullfile(piRoot, 'Results'); cacheDir = fullfile(resDir, 'twinCache');
cfg = defineStudyScenarios();
nPass = 0; nFail = 0;

%% ---- 1  covariance algebra against a direct complex Gaussian simulation ------------------------------
rng(11);
Hs4 = [1 + 0.5i; -0.7 + 0.9i; 0.4 - 1.1i; 1.3 + 0.2i]; n = 400; R = 20000;
Snn4 = [0.05; 0.08; 0.03; 0.2];
cv = analyticTransmissibilityCovariance(Hs4, 1, Snn4, 4, n);
Z = zeros(6, R);
for r = 1:R
    u = (randn(1, n) + 1i * randn(1, n)) / sqrt(2);
    e = (randn(4, n) + 1i * randn(4, n)) / sqrt(2) .* sqrt(Snn4);
    X = Hs4 * u + e; S = X * X' / n;
    T = S(1:3, 4) / (real(S(4, 4)) - Snn4(4)); Z(:, r) = [real(T); imag(T)];
end
Cmc = cov(Z.');
rd = diag(Cmc) ./ diag(cv.blocks); fr = norm(Cmc - cv.blocks, 'fro') / norm(cv.blocks, 'fro');
[nPass, nFail] = check(all(abs(rd - 1) < 0.04) && fr < 0.04, sprintf(['covariance algebra vs %d simulated Wishart draws ' ...
    '(n = %d, noisy reference): variance ratios %.3f-%.3f, Frobenius %.3f'], R, n, min(rd), max(rd), fr), nPass, nFail);

%% ---- 2  cached Monte Carlo: settings F, B, A --------------------------------------------------------
U1 = load(fullfile(resDir, 'transmissibilityCovariance_20260927_124013.mat')); U1 = U1.out;
L2 = load(fullfile(resDir, 'L2048BiasVarianceCheck_20260927_130816.mat')); L2 = L2.out;
twin = synthesiseTwinRecords(struct('p', U1.p0, 'sensors', U1.sensors, 'nodeSpacing', 0.1), 'CacheDir', cacheDir);
[~, ~, info] = synthesiseTwinRecords(twin, U1.N, U1.dt, 'Seed', 1, 'Spectrum', 'jonswap', 'Hs', 0.05, ...
    'PeakFrequency', 5, 'PeakEnhancement', 3.3);
Snn = L2.Snn; ref = 4;
[nEffF, KF] = welchNEff(U1.N, 2048, 0.5);
[nPass, nFail] = check(abs(nEffF - L2.cand(3).nEff) < 1e-9, sprintf('Welch nEff at L = 2048: %.3f (%d segments), as the estimator', ...
    nEffF, KF), nPass, nFail);
cases = struct('name', {'F (L 2048, production)', 'B (L 1024)', 'A (L 512)'}, 'L', {2048, 1024, 512}, ...
    'omega', {L2.cand(3).omega, U1.omegaDec, U1.omega512}, 'J', {L2.cand(3).J, U1.treatments(2).J, U1.treatments(1).J}, ...
    'f0', {L2.cand(3).f0, U1.treatments(2).f0, U1.treatments(1).f0}, ...
    'sigTarget', {cfg.E1.sigmaF(:), U1.treatments(2).res.sigma, U1.treatments(1).res.sigma}, ...
    'dTarget', {cfg.E1.dPred.F2048, cfg.E1.dPred.B1024, cfg.E1.dPred.A512});
for c = 1:numel(cases)
    w = cases(c).omega;
    nE = welchNEff(U1.N, cases(c).L, 0.5);
    Ssm = welchSmoothedSpectrum(w, info, cases(c).L, 1:4);
    cvc = analyticTransmissibilityCovariance(Ssm, [], Snn, ref, nE);
    TW = welchExpectedTransmissibility(w, info, cases(c).L, ref, 1:3, 0);
    q = fisherFromBlocks(cases(c).J, cvc.blocks, stackTransmissibility(TW) - cases(c).f0);
    eS = q.sigma ./ cases(c).sigTarget - 1; eD = q.d / cases(c).dTarget - 1;
    [nPass, nFail] = check(all(abs(eS) <= cfg.E1.tolSigma) && abs(eD) <= cfg.E1.tolD, sprintf(['%s: sigma [%.5f %.6f] vs ' ...
        '[%.5f %.6f] (%+.1f%%, %+.1f%%); d_W %.3f vs %.2f (%+.0f%%)'], cases(c).name, q.sigma, cases(c).sigTarget, ...
        100 * eS, q.d, cases(c).dTarget, 100 * eD), nPass, nFail);
    if c == 1, cvF = cvc; end
end

% covariance STRUCTURE at F: every element of every block, standardised by the sampling error of a
% 200-draw covariance, Var(S_ij) = (S_ij^2 + S_ii S_jj) / (n - 1). Monte Carlo noise alone gives z with
% rms 1 and about 4.6% beyond |z| = 2; off-diagonal (correlation) elements are reported separately.
B = L2.cand(3).blocks; nB = size(B, 3); nMC = L2.nMC; zd = []; zo = [];
for k = 1:nB
    A = cvF.blocks(:, :, k); M = B(:, :, k);
    for i = 1:6
        for j = i:6
            z = (M(i, j) - A(i, j)) / sqrt((A(i, j)^2 + A(i, i) * A(j, j)) / (nMC - 1));
            if i == j, zd(end + 1) = z; else, zo(end + 1) = z; end %#ok<AGROW>
        end
    end
end
rz = @(z) sqrt(mean(z.^2)); fz = @(z) mean(abs(z) > 2);
ok = rz(zd) <= 1.3 && rz(zo) <= 1.3 && fz(zd) <= 0.10 && fz(zo) <= 0.10 && abs(mean(zo)) < 0.2;
[nPass, nFail] = check(ok, sprintf(['Sigma_k structure at F, %d bins: variances z rms %.2f (%.1f%% beyond 2), ' ...
    'covariances z rms %.2f (%.1f%% beyond 2), mean %+.2f; pure Monte Carlo noise: rms 1, 4.6%%'], nB, rz(zd), ...
    100 * fz(zd), rz(zo), 100 * fz(zo), mean(zo)), nPass, nFail);

%% ---- 3  invariances ------------------------------------------------------------------------------------
wF = L2.cand(3).omega; Hb = twin.Hacc(wF);
dHs = complex(randn(4, numel(wF), 2), randn(4, numel(wF), 2));       % any derivative field will do
Ssm = welchSmoothedSpectrum(wF, info, 2048, 1:4);
lay2 = [3 1 2 4];                                                     % same sensors, permuted
[~, ~, J1] = transmissibilityFromField(Hb, dHs, 1:4, 4);
[~, ~, J2] = transmissibilityFromField(Hb, dHs, lay2, 4);
c1 = analyticTransmissibilityCovariance(Ssm, [], Snn, 4, nEffF);
c2 = analyticTransmissibilityCovariance(Ssm(lay2, lay2, :), [], Snn, 4, nEffF);
q1 = fisherFromBlocks(J1, c1.blocks); q2 = fisherFromBlocks(J2, c2.blocks);
eP = norm(q1.F - q2.F) / norm(q1.F);
[nPass, nFail] = check(eP < 1e-10, sprintf('permuting the non-reference sensors leaves F unchanged (rel. diff %.1e)', eP), nPass, nFail);
c3 = analyticTransmissibilityCovariance(Ssm, [], Snn, 4, 2 * nEffF);
q3 = fisherFromBlocks(J1, c3.blocks);
ratio = q3.F ./ q1.F; want = (2 * nEffF - 1) / (nEffF - 1);
[nPass, nFail] = check(max(abs(ratio(:) / want - 1)) < 1e-10, sprintf('F scales with the Welch degrees of freedom (x%.3f)', want), ...
    nPass, nFail);

%% ---- 4  evaluateScenario on an analytic test field (fast, no EMM) -----------------------------------
Sx = cfg.layout.level2C.frac .* [cfg.p0.R 1];
frf = @(w, p) testField(w, p, Sx, cfg.p0.vec, true);
base = struct('sensors', Sx, 'frf', frf, 'Snn', 1e-6);
m0 = evaluateScenario(base, cfg);
mN = evaluateScenario(setfield(base, 'nuisanceIdx', 2), cfg); %#ok<SFLD>
eF = norm(mN.F - m0.F) / norm(m0.F);
ok = isequal(size(mN.F), [2 2]) && eF < 1e-12 && size(mN.Jnuisance, 2) == 1 && norm(mN.Jnuisance) > 0;
[nPass, nFail] = check(ok, sprintf(['gamma as nuisance: F stays 2 x 2 and unchanged (rel. diff %.1e), d f/d ln gamma ' ...
    'returned separately (norm %.2g)'], eF, norm(mN.Jnuisance)), nPass, nFail);
threw = false;
try, evaluateScenario(setfield(base, 'thetaIdx', [1 3 2]), cfg); catch, threw = true; end %#ok<SFLD>
[nPass, nFail] = check(threw, 'the old thetaIdx field is refused, so gamma cannot be estimated by accident', nPass, nFail);
threw = false;
try, evaluateScenario(setfield(base, 'estIdx', [1 2 3]), cfg); catch, threw = true; end %#ok<SFLD>
threw2 = false;
try, evaluateScenario(setfield(base, 'estIdx', [3 1]), cfg); catch, threw2 = true; end %#ok<SFLD>
[nPass, nFail] = check(threw && threw2, 'estIdx locked to [1 3]: [1 2 3] (gamma into F) and [3 1] (swapped sigma) refused', ...
    nPass, nFail);

mR = evaluateScenario(setfield(base, 'Snn', 3e-3), cfg); %#ok<SFLD>
v = mR.validMask; pB = size(mR.blocksAll, 1);
rows = reshape((find(v) - 1) * pB + (1:pB).', [], 1);
qv = fisherFromBlocks(mR.Jall(rows, :), mR.blocksAll(:, :, v));
qa = fisherFromBlocks(mR.Jall, mR.blocksAll);
dF = qa.F - qv.F;                                 % information the rejected bins would have added
ok = mR.invalidFraction > 0 && mR.invalidFraction < 1 && norm(mR.F - qv.F) / norm(qv.F) < 1e-12 && ...
    min(eig((dF + dF.') / 2)) >= -1e-12 * norm(qv.F) && mR.nUsedBins == nnz(v);
[nPass, nFail] = check(ok, sprintf(['rejected bins: %.0f%% of bins fail the estimator''s rule; F equals the manual ' ...
    'valid-bin F, uses %d bins, and the removed bins would only have added information (%.2g%% of F, low SNR)'], ...
    100 * mR.invalidFraction, mR.nUsedBins, 100 * norm(dF) / norm(qv.F)), nPass, nFail);

mS = evaluateScenario(setfield(base, 'frf', @(w, p) testField(w, p, Sx, cfg.p0.vec, false)), cfg); %#ok<SFLD>
ok = mS.singular && all(isinf(mS.sigma)) && strcmp(mS.class, 'Not identifiable');
[nPass, nFail] = check(ok, sprintf(['field with no beta dependence: singular F handled, sigma = [%g %g], class ''%s'' ' ...
    '(no error, no regularisation)'], mS.sigma, mS.class), nPass, nFail);

%% ---- 5  FULL: production chain through evaluateScenario (EMM) -------------------------------------------
for i = 1:numel(FULL_POINTS)
    name = FULL_POINTS{i};
    switch name
        case 'p0', p = cfg.p0.vec; tgt = cfg.E1.sigmaF(:); dT = cfg.E1.dPred.F2048;
        case 'A',  p = cfg.E1.truthA; tgt = cfg.E1.sigmaTruthA(:); dT = NaN;
        case 'B',  p = cfg.E1.truthB; tgt = cfg.E1.sigmaTruthB(:); dT = NaN;
    end
    scn = struct('p', p, 'sensors', cfg.layout.level2C.frac .* [cfg.p0.R 1], 'Snn', Snn, 'cacheDir', cacheDir, ...
        'verbose', true);
    fprintf(['  ... [%d/%d] evaluateScenario at %s: 5 EMM field evaluations of about 28 node solves each (''solved'' takes ' ...
        'minutes, ''loaded'' is instant)\n'], i, numel(FULL_POINTS), name);
    m = evaluateScenario(scn, cfg);
    eS = m.sigma ./ tgt - 1;
    ok = all(abs(eS) <= cfg.E1.tolSigma);
    msg = sprintf('evaluateScenario at %s: sigma [%.5f %.6f] vs [%.5f %.6f] (%+.1f%%, %+.1f%%), corr %.3f, kappa %.2e, %s', ...
        name, m.sigma, tgt, 100 * eS, m.corr, m.condF, m.class);
    if ~isnan(dT)
        ok = ok && abs(m.dW / dT - 1) <= cfg.E1.tolD;
        msg = [msg, sprintf(', d_W %.3f vs %.2f', m.dW, dT)]; %#ok<AGROW>
    end
    [nPass, nFail] = check(ok, sprintf('%s (%.0f s)', msg, m.seconds), nPass, nFail);
end

fprintf('\n%d passed, %d failed%s\n', nPass, nFail, ternary(isempty(FULL_POINTS), ...
    ' (quick part only: set FULL_POINTS for the EMM chain)', ''));
if nFail > 0, error('testAnalyticEngine:failed', '%d check(s) failed.', nFail); end

function H = testField(w, p, S, p0, betaDependent)
% Smooth analytic acceleration FRF with known dependence on beta, gamma, R (test only).
b = log(p(1) / p0(1)) * betaDependent; gm = log(p(2) / p0(2)); r = p(3) / p0(3);
w = w(:).'; H = zeros(size(S, 1), numel(w));
for j = 1:size(S, 1)
    x = S(j, 1) / p0(3); th = S(j, 2);
    H(j, :) = -w.^2 .* (1 + 0.6 * x * exp(1i * (w * (0.4 + 0.3 * x) * (1 + 0.2 * b) / r + th)) + ...
        0.05 * gm * cos(th) * w / 8);
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end