%% testDirectionalForcingStudy.m
% Pre-flight checks for study V2 and for the shared forcing model directionalSpectralMatrix (D3 uses it). All must
% pass before runDirectionalForcingStudy is run. Checks 1 to 11 and 13 use an analytic field (check 9 also simulates
% 100 twin records, about 10 to 30 s); check 12 needs one EMM field at the rotated Level 2C points (29 nodes, ~40 s).
%
%   1  cfg.V2: registered angles and fractions unchanged; the 2026-10-05 clarification is logged
%   2  one component at heading 0 (Welch form) is the engine: the smoothed matrix equals welchSmoothedSpectrum, and
%      through the engine covariance f_W = evaluateScenario f0 + b_W, coherence and Sigma blocks equal the engine's
%   3  f2 = 0 through the two-component path ([1 0] fractions) reproduces the single direction exactly
%   4  two components with the SAME heading (0.7 / 0.3) reproduce the single direction (splitting energy is a no-op)
%   5  Hermitian; signal matrix positive semidefinite; with noise positive definite (every bin, every case)
%   6  energy: fractions must sum to 1 (0.6 + 0.6 is refused); a sensor at the floe centre sees the same S_rr for
%      every heading and fraction (total incident energy preserved)
%   7  rotation: +phi and -phi give the same matrix at the on-axis Level 2C layout (mirror symmetry); a 90 deg
%      heading on the 33-point grid equals the heading-0 field at the grid rotated by two 45 deg steps
%   8  noiseless single source (point form): coherence 1 and rank one; a second direction gives rank two and
%      coherence < 1
%   9  Fisher z: Monte Carlo of the sample coherence from n = 60 complex Gaussian segments gives var(atanh|g|) within
%      10% of 1/(2 (n - 1)); the 80%-power detectability threshold is z(1 - 0.05/3) + z(0.8) = 2.97 for m = 3 pairs;
%      the END-TO-END null (v2Tools.nullCoherence: cfg.V2.coherence.nullRecordsTest twin records through
%      extractTransmissibility) runs, matches every retained bin, and its sd of the bin-mean Fisher z is within a factor
%      2 in variance of the independent-bin value (machinery check; the run uses the calibrated sd whatever it is)
%  10  v2Tools.evaluate at f2 = 0: b_dir = 0, d_sys = 0, D = 0, not detected, class and d_tot of the engine; at
%      f2 > 0 the bias is non-zero and identical for +phi and -phi
%  11  spread weights (D3): unit sum, symmetric about the mean (also at 170 deg, across the branch cut),
%      concentrating as s grows
%  12  EMM: the rotated-point field at p0 reproduces evaluateScenario (f0 + b_W, coherence, Sigma) to 1e-10 and the
%      mirror symmetry of check 7 holds on the EMM to 1e-9
%  13  local regime (analytic): d_sys / f2 constant to 1% over cfg.V2.linearityF2 at every angle (b_dir ~ f2)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

if ~exist('RUN_EMM', 'var'), RUN_EMM = true; end   % set RUN_EMM = false before calling to skip check 12
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
P = setupStudyPaths();
cfg = defineStudyScenarios();
T = v2Tools(); V = cfg.V2;
nPass = 0; nFail = 0;
rel = @(a, b) max(abs(a(:) - b(:))) / max(abs(b(:)));
p0 = cfg.p0.vec;
frfFor = @(sens) @(w, pp) fieldDir(w, pp, p0, sens);
L2 = cfg.layout.level2C.frac .* [p0(3) 1];

%% ---- 1  configuration ----------------------------------------------------------------------------------------------
ok = isequal(V.angleDeg, [15 30 60 90]) && isequal(V.secondEnergyFraction, [0.05 0.10 0.25 0.50]) && ...
    V.primaryHeadingDeg == 0 && strcmp(V.energyConvention, 'fixedTotal') && ...
    any(arrayfun(@(c) ~isempty(strfind(c.field, 'V2.primaryHeadingDeg')), cfg.meta.changes));
[nPass, nFail] = check(ok, sprintf('1  cfg.V2: angles %s deg, f2 %s, primary %g deg, %s; clarification logged', ...
    mat2str(V.angleDeg), mat2str(V.secondEnergyFraction), V.primaryHeadingDeg, V.energyConvention), nPass, nFail);

%% ---- analytic nominal scenario and context ----------------------------------------------------------------------
oA = evaluateScenario(struct('Snn', cfg.noise.lsm6dsv16x, 'frf', frfFor(L2)), cfg);
headings = deg2radP([0 V.angleDeg -V.mirrorCheckDeg]);
cA = T.context(oA, headings, cfg, '', frfFor);

%% ---- 2  one component = the engine -----------------------------------------------------------------------------------
c0 = T.components(cA.pts, cA.layout, 0, 1);
D0 = directionalSpectralMatrix(cA.src, c0, cA.Snn);
Seng = welchSmoothedSpectrum(cA.src.omega, cA.src.info, cA.src.L, cA.pts.rows(1, :));
e2a = rel(D0.Ssig, Seng);
e2b = rel(cA.fW0, oA.f0(:) + oA.bW(:));
e2c = rel(cA.coh0, oA.coherence); e2d = rel(cA.blocks0, oA.blocksAll);
ok = e2a < 1e-14 && e2b < 1e-10 && e2c < 1e-10 && e2d < 1e-10 && isequal(cA.valid0, oA.validMask);
[nPass, nFail] = check(ok, sprintf(['2  one component = engine: smoothed S %.1e, f_W vs f0 + b_W %.1e, coherence %.1e, ' ...
    'Sigma blocks %.1e, accepted bins equal %d'], e2a, e2b, e2c, e2d, isequal(cA.valid0, oA.validMask)), nPass, nFail);

%% ---- 3  f2 = 0 through the two-component path ---------------------------------------------------------------------
c3 = T.components(cA.pts, cA.layout, [0 deg2radP(30)], [1 0]);
D3 = directionalSpectralMatrix(cA.src, c3, cA.Snn);
ok = isequal(D3.Ssig, D0.Ssig) && D3.nComponents == 2;
[nPass, nFail] = check(ok, sprintf('3  fractions [1 0]: identical to the single direction (max diff %.1e)', ...
    max(abs(D3.Ssig(:) - D0.Ssig(:)))), nPass, nFail);

%% ---- 4  same heading, split energy ---------------------------------------------------------------------------------
c4 = T.components(cA.pts, cA.layout, [0 0], [0.7 0.3]);
D4 = directionalSpectralMatrix(cA.src, c4, cA.Snn);
e4 = rel(D4.Ssig, D0.Ssig);
[nPass, nFail] = check(e4 < 1e-14, sprintf('4  0.7 + 0.3 at the same heading = single direction (%.1e)', e4), nPass, nFail);

%% ---- 5  Hermitian, PSD, PD with noise -----------------------------------------------------------------------------
worstH = 0; worstNeg = 0; minPD = Inf;
for phi = deg2radP(V.angleDeg)
    for f2 = V.secondEnergyFraction
        D = directionalSpectralMatrix(cA.src, T.twoComponents(cA.pts, cA.layout, phi, f2), cA.Snn);
        for k = 1:size(D.Ssig, 3)
            A = D.Ssig(:, :, k); B = D.S(:, :, k);
            worstH = max(worstH, norm(A - A', 'fro') / norm(A, 'fro'));
            ev = eig((A + A') / 2); worstNeg = max(worstNeg, -min(ev) / max(ev));
            minPD = min(minPD, min(eig((B + B') / 2)));
        end
    end
end
ok = worstH < 1e-14 && worstNeg < 1e-10 && minPD > 0;
[nPass, nFail] = check(ok, sprintf(['5  16 cases x %d bins: Hermitian (%.1e), signal PSD (most negative eigenvalue ' ...
    '%.1e of the largest), with noise PD (min eigenvalue %.2e)'], size(D0.Ssig, 3), worstH, worstNeg, minPD), nPass, nFail);

%% ---- 6  energy --------------------------------------------------------------------------------------------------------
refused = false;
try
    directionalSpectralMatrix(cA.src, T.components(cA.pts, cA.layout, [0 deg2radP(30)], [0.6 0.6]), cA.Snn);
catch err
    refused = strcmp(err.identifier, 'directionalSpectralMatrix:energy');
end
centre = [0 0; L2];
c6 = T.context(evaluateScenario(struct('Snn', cfg.noise.lsm6dsv16x, 'frf', frfFor(centre), 'sensors', centre, ...
    'ref', 5), cfg), headings, cfg, '', frfFor);
srr = squeeze(real(directionalSpectralMatrix(c6.src, T.components(c6.pts, c6.layout, 0, 1), c6.Snn).Ssig(1, 1, :)));
e6 = 0;
for phi = deg2radP(V.angleDeg)
    for f2 = V.secondEnergyFraction
        D = directionalSpectralMatrix(c6.src, T.twoComponents(c6.pts, c6.layout, phi, f2), c6.Snn);
        e6 = max(e6, rel(squeeze(real(D.Ssig(1, 1, :))), srr));
    end
end
ok = refused && e6 < 1e-13;
[nPass, nFail] = check(ok, sprintf(['6  energy: unit sum enforced (%s); centre sensor S_00 unchanged over 16 cases ' ...
    '(%.1e)'], tf(refused), e6), nPass, nFail);

%% ---- 7  rotation ------------------------------------------------------------------------------------------------------
phiM = deg2radP(V.mirrorCheckDeg);
Dp = directionalSpectralMatrix(cA.src, T.twoComponents(cA.pts, cA.layout, phiM, 0.25), cA.Snn);
Dm = directionalSpectralMatrix(cA.src, T.twoComponents(cA.pts, cA.layout, -phiM, 0.25), cA.Snn);
e7a = rel(Dp.Ssig, Dm.Ssig);
grid = cfg.layout.gridFrac .* [p0(3) 1];
pg = T.headingPoints(grid, deg2radP([0 90]));
Hg = fieldDir(cA.src.omega, p0, p0, pg.sensors);
perm = zeros(1, size(grid, 1));                                % grid point at theta - 90 deg
for i = 1:size(grid, 1)
    if grid(i, 1) == 0, perm(i) = i; continue; end           % the centre is its own rotation
    tgt = [grid(i, 1), grid(i, 2) - pi / 2];
    perm(i) = find(abs(grid(:, 1) - tgt(1)) < 1e-12 & abs(angle(exp(1i * (grid(:, 2) - tgt(2))))) < 1e-12, 1);
end
e7b = rel(Hg(pg.rows(2, :), :), Hg(pg.rows(1, perm), :));
ok = e7a < 1e-12 && e7b < 1e-12;
[nPass, nFail] = check(ok, sprintf(['7  rotation: +%g and -%g deg equal at Level 2C (%.1e); 90 deg on the grid = ' ...
    'heading 0 at the grid permuted by two 45 deg steps (%.1e)'], V.mirrorCheckDeg, V.mirrorCheckDeg, e7a, e7b), ...
    nPass, nFail);

%% ---- 8  rank and coherence, point form ------------------------------------------------------------------------------
w8 = cA.src.omega; H8 = fieldDir(w8, p0, p0, cA.pts.sensors); Su8 = ones(1, numel(w8));
src8 = struct('H', H8, 'Su', Su8);
P1 = directionalSpectralMatrix(src8, T.components(cA.pts, cA.layout, 0, 1), 0);
P2 = directionalSpectralMatrix(src8, T.twoComponents(cA.pts, cA.layout, deg2radP(60), 0.25), 0);
cv1 = analyticTransmissibilityCovariance(P1.Ssig, [], 0, cA.ref, cA.nEff);
cv2 = analyticTransmissibilityCovariance(P2.Ssig, [], 0, cA.ref, cA.nEff);
r1 = 0; r2 = Inf; r3 = 0;
for k = 1:numel(w8)
    e1 = sort(real(eig(P1.Ssig(:, :, k))), 'descend'); e2 = sort(real(eig(P2.Ssig(:, :, k))), 'descend');
    r1 = max(r1, e1(2) / e1(1)); r2 = min(r2, e2(2) / e2(1)); r3 = max(r3, e2(3) / e2(1));
end
ok = max(abs(cv1.coherence(:) - 1)) < 1e-12 && r1 < 1e-12 && r2 > 1e-6 && r3 < 1e-12 && max(cv2.coherence(:)) < 1;
[nPass, nFail] = check(ok, sprintf(['8  noiseless point form: single coherence 1 (%.1e), lambda2/lambda1 %.1e; two ' ...
    'directions: lambda2/lambda1 >= %.1e, lambda3/lambda1 %.1e, max coherence %.6f'], max(abs(cv1.coherence(:) - 1)), ...
    r1, r2, r3, max(cv2.coherence(:))), nPass, nFail);

%% ---- 9  Fisher z variance and the threshold -------------------------------------------------------------------------
n9 = 60; nMC = 4000; g9 = 0.85;                                % true |coherence|
Cs = [1, g9; g9, 1]; Lc = chol(Cs, 'lower');
zz = zeros(nMC, 1); rng(9, 'twister');
for i = 1:nMC
    X = Lc * (randn(2, n9) + 1i * randn(2, n9)) / sqrt(2);
    S9 = X * X' / n9;
    zz(i) = atanh(abs(S9(1, 2)) / sqrt(real(S9(1, 1)) * real(S9(2, 2))));
end
vr = var(zz) * 2 * (n9 - 1);
thr = T.normInv(1 - V.coherence.alpha / 3) + T.normInv(V.coherence.power);
st0 = T.coherenceStats(cA.coh0, cA.coh0, cA.valid0, cA.nDof, V.coherence);
fprintf('  ... end-to-end coherence null (%d analytic twin records)\n', V.coherence.nullRecordsTest);
nc = T.nullCoherence(oA, cfg, V.coherence.nullRecordsTest, V.coherence.nullSeedBase);
stC = T.coherenceStats(cA.coh0, cA.coh0 * 0.999, cA.valid0, cA.nDof, V.coherence, nc.sd);
dzC = mean(atanh(sqrt(cA.coh0(:, cA.valid0))) - atanh(sqrt(0.999 * cA.coh0(:, cA.valid0))), 2);
ok = abs(vr - 1) < 0.10 && abs(thr - 2.9700) < 5e-4 && all(st0.D == 0) && ~st0.detect && ...
    all(isfinite(nc.sd)) && all(nc.varianceRatio > 0.5 & nc.varianceRatio < 2) && max(abs(stC.D - dzC ./ nc.sd)) < 1e-12;
[nPass, nFail] = check(ok, sprintf(['9  Fisher z: var x 2(n - 1) = %.3f (n = %d, %d draws); 80%%-power threshold %.4f ' ...
    '(m = 3); D = 0 for identical coherence; null MC (%d records, %.0f s): variance ratio to independent bins %s, ' ...
    'effective bins %s of %d, mean offset %s sd'], vr, n9, nMC, thr, nc.nRec, nc.seconds, mat2str(nc.varianceRatio.', 3), ...
    mat2str(round(nc.nBinsEff.'), 4), nc.nBins, mat2str((nc.mean ./ nc.sd).', 3)), nPass, nFail);

%% ---- 10  v2Tools.evaluate ------------------------------------------------------------------------------------------
r0 = T.evaluate(cA, oA, deg2radP(30), 0, cfg);
rp = T.evaluate(cA, oA, deg2radP(V.mirrorCheckDeg), 0.25, cfg);
rm = T.evaluate(cA, oA, -deg2radP(V.mirrorCheckDeg), 0.25, cfg);
ok = all(r0.b == 0) && r0.dSys == 0 && all(r0.D == 0) && ~r0.detect && strcmp(r0.class, oA.class) && ...
    abs(r0.dTot - oA.dW) < 1e-12 && rp.dSys > 0 && abs(rp.dSys - rm.dSys) < 1e-10 * rp.dSys && ...
    all(abs(r0.sigmaRatio - 1) < 1e-10);
[nPass, nFail] = check(ok, sprintf(['10 evaluate: f2 = 0 gives b = 0, d_sys 0, D 0, d_tot = d_W, class %s, sigma ratio 1; ' ...
    'analytic f2 = 0.25 at +-%g deg: d_sys %.3f / %.3f, D_max %.2f'], r0.class, V.mirrorCheckDeg, rp.dSys, rm.dSys, ...
    rp.Dmax), nPass, nFail);

%% ---- 11  spread weights ---------------------------------------------------------------------------------------------
[h1, w1] = T.spreadWeights(0, 2, 24); [~, w2] = T.spreadWeights(0, 20, 24);
sym = max(abs(w1 - w1(mod(-(0:23), 24) + 1)));
mu = deg2radP(170); [h3, w3] = T.spreadWeights(mu, 2, 24);
sym3 = max(abs(w3 - w3(mod(-(0:23), 24) + 1)));                % mirror about the mean, index k <-> -k
peak3 = abs(angle(exp(1i * (h3(w3 == max(w3)) - mu))));
ok = abs(sum(w1) - 1) < 1e-14 && sym < 1e-14 && w2(1) > w1(1) && abs(h1(1)) < 1e-15 && ...
    abs(sum(w3) - 1) < 1e-14 && sym3 < 1e-14 && all(peak3 < 1e-12) && max(abs(sort(w3) - sort(w1))) < 1e-14;
[nPass, nFail] = check(ok, sprintf(['11 spread cos^2s (24 headings): sum %.15f, symmetry %.1e, weight at the mean %.3f ' ...
    '(s = 2) -> %.3f (s = 20); mean 170 deg: symmetry %.1e, peak at the mean, same weights as mean 0'], sum(w1), sym, ...
    w1(1), w2(1), sym3), nPass, nFail);

%% ---- 13  local regime (analytic) ----------------------------------------------------------------------------------
sp = zeros(1, numel(V.angleDeg));
for i = 1:numel(V.angleDeg)
    Li = T.linearity(cA, oA, deg2radP(V.angleDeg(i)), V.linearityF2, cfg);
    sp(i) = Li.spread;
end
[nPass, nFail] = check(all(sp < 0.01), sprintf(['13 local regime: d_sys / f2 spread over f2 = %s is %s %% (phi = %s); ' ...
    'b_dir ~ f2'], mat2str(V.linearityF2), mat2str(100 * sp, 3), mat2str(V.angleDeg)), nPass, nFail);

%% ---- 12  EMM ------------------------------------------------------------------------------------------------------------
if RUN_EMM
    fprintf('  ... EMM: nominal p0 (cached) and one field at the rotated Level 2C points (about 40 s)\n');
    oE = evaluateScenario(struct('cacheDir', P.cache), cfg);
    cE = T.context(oE, headings, cfg, P.cache);
    e12a = rel(cE.fW0, oE.f0(:) + oE.bW(:)); e12b = rel(cE.coh0, oE.coherence); e12c = rel(cE.blocks0, oE.blocksAll);
    Ep = directionalSpectralMatrix(cE.src, T.twoComponents(cE.pts, cE.layout, phiM, 0.25), cE.Snn);
    Em = directionalSpectralMatrix(cE.src, T.twoComponents(cE.pts, cE.layout, -phiM, 0.25), cE.Snn);
    e12d = rel(Ep.Ssig, Em.Ssig);
    ok = e12a < 1e-10 && e12b < 1e-10 && e12c < 1e-10 && e12d < 1e-9 && isequal(cE.valid0, oE.validMask);
    [nPass, nFail] = check(ok, sprintf(['12 EMM p0 (%.0f s): f_W vs f0 + b_W %.1e, coherence %.1e, Sigma %.1e; mirror ' ...
        '+-%g deg %.1e'], cE.fieldSeconds, e12a, e12b, e12c, V.mirrorCheckDeg, e12d), nPass, nFail);
else
    fprintf('  SKIP 12 (RUN_EMM = false)\n');
end

fprintf('\n  testDirectionalForcingStudy: %d passed, %d failed (13 checks)\n', nPass, nFail);
if nFail > 0, error('testDirectionalForcingStudy:failed', '%d check(s) failed.', nFail); end

%% ================================================================================================
function [nPass, nFail] = check(ok, msg, nPass, nFail)
if ok
    nPass = nPass + 1; fprintf('  PASS %s\n', msg);
else
    nFail = nFail + 1; fprintf('  FAIL %s\n', msg);
end
end

function s = tf(b)
if b, s = 'ok'; else, s = 'NO'; end
end

function r = deg2radP(d)
r = d * pi / 180;
end

function H = fieldDir(w, pp, p0, sens)
% smooth analytic test field of a circular floe, mirror symmetric in theta (incidence along +x), depending on beta
% and R; the centre (r = 0) is rotation invariant
w = w(:).'; b = (pp(1) / p0(1))^0.35; x = sens(:, 1) / pp(3); c = cos(sens(:, 2)); c2 = cos(2 * sens(:, 2));
H = -(ones(size(sens, 1), 1) * w.^2) .* (1 + 0.6 * (x.^2 .* c) * ones(1, numel(w)) .* exp(1i * (0.9 * b * x .* c) * w) ...
    + 0.25 * b * (x * w) + 0.2 * (x.^2 .* c2) * ones(1, numel(w)) .* exp(1i * (0.5 * x) * w));
end