%% testDirectionRobustLayoutStudy.m
% Pre-flight checks for study D3. All must pass before runDirectionRobustLayoutStudy is run. Checks 1 to 12 and 14 use an
% analytic field (a few minutes in Octave); check 13 builds the production EMM field (10 grid copies, 330 points, 5 solves; a few
% minutes the first time, then cached for the run).
%
%   1  cfg.D3: headings 7.5 + 15k (24), dpsi 0.5 deg, the seven D3b seas of the programme, gate 0.5; clarification logged
%   2  geometry: every P_k is a permutation, P_k = P_1^k, P_8 = identity; admissible layouts = D1's for n_s = 2, 3, 4
%   3  heading 0, known: d3Tools reproduces d1LayoutTools (engine path) for D1's and Level 2C's layouts and
%      references, and the n_s = 2 search has the same optimum as D1's search on the same field
%   4  the heading nuisance never improves (beta, R): A_nuisance >= A_known for every n_s = 2 layout at delta = 7.5
%   5  with dT/dpsi set to zero the nuisance metrics equal the known-heading metrics
%   6  dT/dpsi by central differences is step-stable: +-0.5 deg and +-0.25 deg give the nuisance A within 0.1%
%   7  rotation: a grid copy sampled directly at 52.5 deg gives exactly the heading table's 52.5 deg column
%      (layout P_1(L) at 7.5 deg) for every n_s = 2 layout
%   8  full-circle invariance: mean and worst A over the 24 headings are identical for L and its 45 deg rotation
%   9  D3b with f2 = 0: b_dir = 0, Sigma = the single-direction Sigma, so A and d_tot equal the known-heading A, d_W
%  10  D3b Sigma: the sub-block used for a layout equals the engine covariance computed directly on that layout's
%      directional spectral matrix (directionalSpectralMatrix restricted to the layout)
%  11  D3b nominal ruler at Level 2C, heading 0, reference s26: d_sys and d_tot equal V2's (v2Tools.evaluate) for
%      phi = 30 deg, f2 = 0.25 (cross-check of the two implementations)
%  12  spread sea: 24 components, unit energy, every heading on the sampled copies; the centre point's S_00 equals the
%      single-direction value (total energy preserved)
%  14  D3b reference choice: per layout, the feasible-best reference has d_tot <= gate and A >= the precision-optimal
%      A; the least-biased d_tot <= the precision-optimal d_tot; a feasible reference exists iff d_tot,min <= gate
%  13  EMM: the D3 field's heading-0 copy reproduces D1's full field (H, dH, smoothed S) to 1e-12, and the heading-0
%      known-heading n_s = 2 optimum is D1's (s6 + s18, reference s18)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

if ~exist('RUN_EMM', 'var'), RUN_EMM = true; end   % set RUN_EMM = false before calling to skip check 13
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
P = setupStudyPaths();
cfg = defineStudyScenarios();
T = d3Tools(); D1 = d1LayoutTools(); V = cfg.D3;
nPass = 0; nFail = 0;
rel = @(a, b) max(abs(a(:) - b(:))) / max(abs(b(:)));
p0 = cfg.p0.vec;
frfFor = @(sens) @(w, pp) fieldDir(w, pp, p0, sens);
Snn = cfg.noise.lsm6dsv16x;

%% ---- 1  configuration ----------------------------------------------------------------------------------------------
sd = V.d3b;
ok = numel(V.headingsDeg) == 24 && max(abs(V.headingsDeg - (7.5 + 15 * (0:23)))) == 0 && V.headingStepDeg == 0.5 && ...
    numel(sd) == 7 && isequal([sd(1:6).phiDeg], [30 30 60 60 90 90]) && isequal([sd(1:6).f2], [0.1 0.25 0.1 0.25 0.1 0.25]) && ...
    strcmp(sd(7).kind, 'spread') && V.biasGate == 0.5 && ...
    any(arrayfun(@(c) strcmp(c.field, 'D3.* (added)'), cfg.meta.changes)) && ...
    any(arrayfun(@(c) ~isempty(strfind(c.field, 'D3.headingsDeg (AMENDED')), cfg.meta.changes));
[nPass, nFail] = check(ok, sprintf(['1  cfg.D3: 24 headings 7.5:15:352.5 deg, dpsi %g deg, seas 30/60/90 deg x 10/25%% + ' ...
    'spread s %g, gate %.1f, ruler %s; clarification and heading amendment logged'], V.headingStepDeg, sd(7).s, V.biasGate, V.d3bRuler), ...
    nPass, nFail);

%% ---- 2  geometry and layouts ----------------------------------------------------------------------------------------
geom = T.geometry(cfg, p0);
pm = geom.perm; okP = true;
for k = 0:7
    okP = okP && isequal(sort(pm(k + 1, :)), 1:33);
    if k >= 1, okP = okP && isequal(pm(k + 1, :), pm(2, pm(k, :))); end
end
okP = okP && isequal(pm(2, pm(8, :)), 1:33);
nD3 = zeros(1, 3); nD1 = nD3;
for k = 1:3
    Lk = T.layoutSet(geom, k + 1, cfg); nD3(k) = Lk.n;
    [~, okD1] = D1.enumerateLayouts(geom, k + 1, cfg, cfg.D1.minSpacing); nD1(k) = nnz(okD1);
end
[nPass, nFail] = check(okP && isequal(nD3, nD1), sprintf(['2  geometry: P_k permutations, P_k = P_1^k, P_8 = I (%s); ' ...
    'admissible layouts %s = D1 %s'], tf(okP), mat2str(nD3), mat2str(nD1)), nPass, nFail);

%% ---- analytic field (heading 0, the three deltas, +-0.5 and +-0.25 deg, 52.5 deg with its +-0.5 deg copies, 30 deg) -
fprintf('  ... analytic field at the grid copies\n');
ang = T.copyAngles(cfg, [7.5 + [-0.25 0.25], 52.5 + [-0.5 0 0.5], 30]);
F = T.buildField(cfg, p0, geom, ang, '', frfFor, false);
L2 = T.layoutSet(geom, 2, cfg);

%% ---- 3  heading 0 = D1 ----------------------------------------------------------------------------------------------
full = D1.buildFullField(cfg, p0, geom, '', false, frfFor(geom.sensorsND));
pre1 = D1.precompute(full, Snn, cfg);
pre0 = T.precomputeA(F, 0, Snn, false);
cases = {[6 18], 18; sort(geom.level2C), geom.level2C(cfg.layout.level2C.ref); [2 14 26], 2};
e3 = 0;
for i = 1:size(cases, 1)
    q1 = D1.evaluateLayout(pre1, cases{i, 1}, cases{i, 2}, cfg, 'direct');
    [Fm, G] = T.layoutFisher(pre0(cases{i, 2}), cases{i, 1}, cases{i, 2}, {'bW'});
    q3 = T.metricsA(Fm, G, cfg);
    e3 = max([e3, abs(q3.known.A / q1.A - 1), max(abs(q3.known.sigma ./ q1.sigma - 1)), abs(q3.known.d - q1.dW)]);
end
S1 = D1.searchLayouts(pre1, double(L2.rows), cfg, 'D1');
[~, i1] = min(S1.A); S0 = T.searchA(pre0, L2, cfg, false);
[~, i3] = min(S0.Ak);
ok = e3 < 1e-10 && i1 == i3 && double(S1.refA(i1)) == S0.refK(i3) && abs(S0.Ak(i3) / S1.A(i1) - 1) < 1e-10;
[nPass, nFail] = check(ok, sprintf(['3  heading 0 known = D1: A, sigma, d_W within %.1e for 3 layouts; n_s = 2 optimum ' ...
    '%s ref %s in both'], e3, strjoin(geom.names(L2.rows(i3, :)), '+'), geom.names{S0.refK(i3)}), nPass, nFail);

%% ---- 4  nuisance never improves -----------------------------------------------------------------------------------------
preA = T.precomputeA(F, 7.5, Snn, true);
SA = T.searchA(preA, L2, cfg, true);
worst4 = min(SA.An ./ SA.Ak);
[nPass, nFail] = check(worst4 >= 1 - 1e-12, sprintf(['4  nuisance: min A_nuisance / A_known over %d layouts = %.6f ' ...
    '(median %.3f)'], L2.n, worst4, median(SA.An ./ SA.Ak)), nPass, nFail);

%% ---- 5  zero heading derivative -------------------------------------------------------------------------------------------
lay = [6 18]; r = 18;
P5 = preA(r); P5.dT(:, :, 3) = 0;
[Fm, G] = T.layoutFisher(P5, lay, r, {'bW'}); q5 = T.metricsA(Fm, G, cfg);
ok = abs(q5.nuis.A / q5.known.A - 1) < 1e-14 && max(abs(q5.nuis.sigma - q5.known.sigma)) == 0;
[nPass, nFail] = check(ok, sprintf('5  dT/dpsi = 0: nuisance A / known A - 1 = %.1e', q5.nuis.A / q5.known.A - 1), nPass, nFail);

%% ---- 6  step halving --------------------------------------------------------------------------------------------------------
Fh = F; Fh.dpsi = V.headingStepCheckDeg * pi / 180;
preH = T.precomputeA(Fh, 7.5, Snn, true);
e6 = 0;
for lay6 = {[6 18], sort(geom.level2C), [10 22], [3 15 27]}
    l = lay6{1}; r = l(end);
    [Fa, Ga] = T.layoutFisher(preA(r), l, r, {'bW'}); [Fb, Gb] = T.layoutFisher(preH(r), l, r, {'bW'});
    qa = T.metricsA(Fa, Ga, cfg); qb = T.metricsA(Fb, Gb, cfg);
    e6 = max(e6, abs(qa.nuis.A / qb.nuis.A - 1));
end
[nPass, nFail] = check(e6 < 1e-3, sprintf('6  dpsi %g deg vs %g deg: nuisance A differs by at most %.2e', ...
    V.headingStepDeg, V.headingStepCheckDeg, e6), nPass, nFail);

%% ---- 7  rotation: direct 52.5 deg copy = table column ------------------------------------------------------------------
pre52 = T.precomputeA(F, 52.5, Snn, true);
S52 = T.searchA(pre52, L2, cfg, true);
perD = struct('deltaDeg', 7.5, 'S', SA);
M = T.headingTable(perD, 'An', L2, F, 52.5); Mk = T.headingTable(perD, 'Ak', L2, F, 52.5);
e7 = max(rel(M, S52.An), rel(Mk, S52.Ak));
[nPass, nFail] = check(e7 < 1e-10, sprintf(['7  rotation: direct 52.5 deg copy vs P_1 lookup at 7.5 deg, %d layouts: ' ...
    '%.1e'], L2.n, e7), nPass, nFail);

%% ---- 8  full-circle invariance -----------------------------------------------------------------------------------------------
pd = struct('deltaDeg', num2cell([7.5 22.5 37.5]), 'S', []);
pd(1).S = SA; pd(2).S = T.searchA(T.precomputeA(F, 22.5, Snn, true), L2, cfg, true);
pd(3).S = T.searchA(T.precomputeA(F, 37.5, Snn, true), L2, cfg, true);
An = T.headingTable(pd, 'An', L2, F, V.headingsDeg);
mA = mean(An, 2); wA = max(An, [], 2);
e8 = max(max(abs(mA(L2.rot(:, 2)) - mA) ./ mA), max(abs(wA(L2.rot(:, 2)) - wA) ./ wA));
[nPass, nFail] = check(e8 < 1e-12, sprintf('8  full circle: mean and worst A of L and P_1(L) agree to %.1e', e8), nPass, nFail);

%% ---- 9  D3b with f2 = 0 -----------------------------------------------------------------------------------------------------------
sea0 = struct('kind', 'two', 'phiDeg', 30, 'f2', 0, 's', NaN, 'nDir', NaN);
preB0 = T.precomputeB(F, 7.5, sea0, Snn, 'directional');
SB0 = T.searchB(preB0, L2, cfg);
bmax = max(arrayfun(@(p) max(abs(p.bD(:))) / max(abs(p.bW(:)) + 1), preB0));   % rounding only (two code paths)
ok = bmax < 1e-12 && rel(SB0.A, SA.Ak) < 1e-10 && max(abs(SB0.dTot - SA.dWk)) < 1e-10 && max(SB0.dSys) < 1e-8;
[nPass, nFail] = check(ok, sprintf(['9  D3b f2 = 0: b_dir max %.1e (rounding); A vs known A %.1e; d_tot vs d_W %.1e; ' ...
    'd_sys %.1e'], bmax, ...
    rel(SB0.A, SA.Ak), max(abs(SB0.dTot - SA.dWk)), max(SB0.dSys)), nPass, nFail);

%% ---- 10  D3b Sigma sub-block = engine on the layout -----------------------------------------------------------------------
seaB = struct('kind', 'two', 'phiDeg', 60, 'f2', 0.25, 's', NaN, 'nDir', NaN);
preB = T.precomputeB(F, 7.5, seaB, Snn, 'directional');
lay = [3 15 27]; r = 15;
comps = T.seaComponents(F, 7.5, seaB);
for i = 1:numel(comps), comps(i).rows = comps(i).rows(lay); end
Dl = directionalSpectralMatrix(struct('omega', F.omega, 'info', F.info, 'L', F.L), comps, Snn);
cvl = analyticTransmissibilityCovariance(Dl.Ssig, [], Snn, find(lay == r), F.nEff);
o = lay(lay ~= r); pos = o - (o > r); idx = [pos(:); 32 + pos(:)];
e10 = rel(preB(r).B(idx, idx, :), cvl.blocks(:, :, cvl.valid));
ok = e10 < 1e-12 && isequal(preB(r).valid, cvl.valid);
[nPass, nFail] = check(ok, sprintf('10 D3b Sigma: layout sub-block vs engine on the layout''s own S: %.1e, same accepted bins %d', ...
    e10, isequal(preB(r).valid, cvl.valid)), nPass, nFail);

%% ---- 11  V2 cross-check ------------------------------------------------------------------------------------------------------
L2C = cfg.layout.level2C.frac .* [p0(3) 1];
oA = evaluateScenario(struct('Snn', Snn, 'frf', frfFor(L2C)), cfg);
T2 = v2Tools(); c2 = T2.context(oA, [0 30] * pi / 180, cfg, '', frfFor);
rV = T2.evaluate(c2, oA, 30 * pi / 180, 0.25, cfg);
sea11 = struct('kind', 'two', 'phiDeg', 30, 'f2', 0.25, 's', NaN, 'nDir', NaN);
preN = T.precomputeB(F, 0, sea11, Snn, 'nominal');
ref = geom.level2C(cfg.layout.level2C.ref);
[Fm, G] = T.layoutFisher(preN(ref), sort(geom.level2C), ref, {'bW', 'bD'}); qB = T.metricsB(Fm, G, cfg);
e11 = max(abs(qB.dSys / rV.dSys - 1), abs(qB.dTot / rV.dTot - 1));
[nPass, nFail] = check(e11 < 1e-8, sprintf(['11 V2 cross-check (Level 2C, s26, 30 deg, 0.25, nominal ruler): d_sys %.4f vs %.4f, ' ...
    'd_tot %.4f vs %.4f'], qB.dSys, rV.dSys, qB.dTot, rV.dTot), nPass, nFail);

%% ---- 12  spread sea -----------------------------------------------------------------------------------------------------------
cS = T.seaComponents(F, 7.5, V.d3b(7));
src = struct('omega', F.omega, 'info', F.info, 'L', F.L);
DS = directionalSpectralMatrix(src, cS, Snn);
D1c = directionalSpectralMatrix(src, T.seaComponents(F, 7.5, sea0), Snn);
e12 = rel(squeeze(real(DS.Ssig(1, 1, :))), squeeze(real(D1c.Ssig(1, 1, :))));
ok = numel(cS) == V.d3b(7).nDir && abs(sum([cS.fraction]) - 1) < 1e-14 && e12 < 1e-12;
[nPass, nFail] = check(ok, sprintf('12 spread: %d components, energy %.15f, centre S_00 vs single %.1e', numel(cS), ...
    sum([cS.fraction]), e12), nPass, nFail);

%% ---- 14  D3b reference choice -----------------------------------------------------------------------------------------
S14 = T.searchB(preB, T.layoutSet(geom, 3, cfg), cfg);
fin = isfinite(S14.Af);
ok = all(S14.dTotf(fin) <= V.biasGate) && all(S14.Af(fin) >= S14.A(fin) * (1 - 1e-12)) && ...
    all(S14.dTotMin <= S14.dTot + 1e-12 | isnan(S14.dTot)) && isequal(fin, S14.dTotMin <= V.biasGate);
[nPass, nFail] = check(ok, sprintf(['14 D3b references (n_s = 3, 60 deg, 0.25): %d layouts with a gate-passing reference; ' ...
    'precision-optimal reference passes for %d; min d_tot %.3g (precision-optimal %.3g)'], nnz(fin), ...
    nnz(S14.dTot <= V.biasGate), min(S14.dTotMin), min(S14.dTot)), nPass, nFail);

%% ---- 13  EMM ----------------------------------------------------------------------------------------------------------------------
if RUN_EMM
    fprintf('  ... EMM: the production D3 field (10 copies, cached for the run) and D1''s full field (cached)\n');
    FE = T.buildField(cfg, p0, geom, T.copyAngles(cfg), P.cache, [], true);
    fullE = D1.buildFullField(cfg, p0, geom, P.cache, false);
    r0 = 1:33;
    e13 = [rel(FE.H(r0, :), fullE.H), rel(FE.dH(r0, :, :), fullE.dH), ...
        rel(welchSmoothedSpectrum(FE.omega, FE.info, FE.L, r0), fullE.Ssm)];
    S0E = T.searchA(T.precomputeA(FE, 0, Snn, false), L2, cfg, false);
    [~, iE] = min(S0E.Ak);
    ok = all(e13 < 1e-12) && isequal(L2.rows(iE, :), V.d1Gate.layout) && S0E.refK(iE) == V.d1Gate.ref;
    [nPass, nFail] = check(ok, sprintf(['13 EMM (%.0f s): heading-0 copy vs D1 field H %.1e, dH %.1e, S %.1e; n_s = 2 ' ...
        'optimum %s ref %s (D1: s6+s18 ref s18)'], FE.seconds, e13, strjoin(geom.names(L2.rows(iE, :)), '+'), ...
        geom.names{S0E.refK(iE)}), nPass, nFail);
else
    fprintf('  SKIP 13 (RUN_EMM = false)\n');
end

fprintf('\n  testDirectionRobustLayoutStudy: %d passed, %d failed (14 checks)\n', nPass, nFail);
if nFail > 0, error('testDirectionRobustLayoutStudy:failed', '%d check(s) failed.', nFail); end

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

function H = fieldDir(w, pp, p0, sens)
% smooth analytic test field of a circular floe (as testDirectionalForcingStudy), mirror symmetric about the incidence
w = w(:).'; b = (pp(1) / p0(1))^0.35; x = sens(:, 1) / pp(3); c = cos(sens(:, 2)); c2 = cos(2 * sens(:, 2));
H = -(ones(size(sens, 1), 1) * w.^2) .* (1 + 0.6 * (x.^2 .* c) * ones(1, numel(w)) .* exp(1i * (0.9 * b * x .* c) * w) ...
    + 0.25 * b * (x * w) + 0.2 * (x.^2 .* c2) * ones(1, numel(w)) .* exp(1i * (0.5 * x) * w));
end