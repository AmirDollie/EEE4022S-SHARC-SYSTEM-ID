function T = d3Tools()
%D3TOOLS  The numerical core of study D3 (direction-robust sensor count and layout), shared by
%   runDirectionRobustLayoutStudy and testDirectionRobustLayoutStudy so the tester exercises exactly the code the run uses.
%   T = d3Tools();   then e.g.  geom = T.geometry(cfg, p)
%
%   TWO DIFFERENT PROBLEMS (never blended)
%   D3a  unknown but SINGLE heading: T_j = H_j / H_r still holds (model correct); the heading psi is an extra unknown.
%        theta = [ln beta, ln R, psi]; the (beta, R) information after estimating psi is the Schur complement
%        F_eff = F_tt - F_tp F_pp^-1 F_pt (F_pp = 0 exactly: F_eff = F_tt). Scored PER HEADING (A = tr F_eff^-1), then
%        averaged / worst case over headings; never A of a heading-averaged F.
%   D3b  multi-directional sea: the frozen single-direction inverse (primary heading known) is WRONG. Per layout and
%        condition: Sigma from the rank > 1 spectral matrix (directionalSpectralMatrix, engine covariance), F and the
%        projection weights under that Sigma, b_dir = f_W(dir) - f_W(single), d_tot from b_W + b_dir.
%
%   HEADINGS BY SYMMETRY. Heading psi samples the nominal field at the points rotated by -psi. The 33-point grid has
%   eight-fold symmetry, so heading 45k + delta of layout L equals heading delta of layout P_k(L), where P_k maps each grid
%   point to the grid point at angle theta - 45k. Only the grid copies at delta (and delta +- dpsi for d/dpsi) are
%   sampled; every other heading is an index permutation.
%
%   FUNCTIONS
%   geom = T.geometry(cfg, p)                  D1 geometry plus perm (8 x 33): perm(k + 1, i) = P_k(i)
%   L    = T.layoutSet(geom, ns, cfg)          admissible layouts (D1 rule), keys, rot (nL x 8): row of P_k(layout)
%   a    = T.copyAngles(cfg, extra)            grid-copy angles (deg): 0, the heading offsets delta, delta +- dpsi, extra
%   F    = T.buildField(cfg, p, geom, angles, cacheDir, frfFor, verbose)   ONE emmField at all copies (ln beta, ln R
%                                              derivatives), the retained bins, the record statistics, N_eff
%   r    = T.headingRows(F, hDeg)              rows of the stacked field for the 33 grid points at absolute heading h
%   pre  = T.precomputeA(F, deltaDeg, Snn, withPsi)       per reference: Sigma blocks, dT/d[lnb lnR (psi)], b_W
%   c    = T.seaComponents(F, deltaDeg, sea)   directionalSpectralMatrix components of a D3b sea, primary heading delta
%   pre  = T.precomputeB(F, deltaDeg, sea, Snn, ruler)    per reference under that sea: Sigma, dT/d[lnb lnR], b_W, b_dir
%   [Fm, G] = T.layoutFisher(P, layout, ref, biasFields) the layout's F and J' Sigma^-1 b (sub-block of pre(ref))
%   q    = T.metricsA(F3, g, cfg)              known-heading and heading-nuisance metrics (sigma, A, d_W, kappa, class)
%   q    = T.metricsB(F2, G2, cfg)             D3b metrics (A, sigma, d_sys, d_tot, kappa, class)
%   S    = T.searchA(pre, L, cfg, withPsi, label)   every layout x reference; best reference per layout
%   S    = T.searchB(pre, L, cfg, label)       A-best, feasible-best (d_tot <= gate) and least-biased reference
%   M    = T.headingTable(perDelta, field, L, F, hDeg)  nL x nH table of a search field over absolute headings
%   c    = T.classOf(sigma, d, kappa, singular, cfg)    1 Identifiable, 2 Marginal, 3 Not identifiable
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
T.geometry = @geometry;
T.layoutSet = @layoutSet;
T.copyAngles = @copyAngles;
T.buildField = @buildField;
T.headingRows = @headingRows;
T.precomputeA = @precomputeA;
T.seaComponents = @seaComponents;
T.precomputeB = @precomputeB;
T.layoutFisher = @layoutFisher;
T.metricsA = @metricsA;
T.metricsB = @metricsB;
T.searchA = @searchA;
T.searchB = @searchB;
T.headingTable = @headingTable;
T.classOf = @classOf;
T.version = 'd3Tools 2026-10-05';
end

%% ================================================================================================
function geom = geometry(cfg, p)
D1 = d1LayoutTools();
geom = D1.buildGeometry(cfg, p);
nC = size(geom.frac, 1);
ring = zeros(nC, 1); ang = zeros(nC, 1);
rings = unique(geom.frac(geom.frac(:, 1) > 0, 1));
for i = 1:nC
    if geom.frac(i, 1) > 0
        ring(i) = find(abs(rings - geom.frac(i, 1)) < 1e-12);
        ang(i) = mod(round(geom.frac(i, 2) / (pi / 4)), 8);
    end
end
perm = zeros(8, nC);
for k = 0:7
    for i = 1:nC
        if ring(i) == 0, perm(k + 1, i) = i; continue; end
        perm(k + 1, i) = find(ring == ring(i) & ang == mod(ang(i) - k, 8));
    end
end
geom.perm = perm; geom.ring = ring; geom.angIdx = ang;
end

function L = layoutSet(geom, ns, cfg)
D1 = d1LayoutTools();
[lay, ok] = D1.enumerateLayouts(geom, ns, cfg, cfg.D3.minSpacing);
rows = lay(ok, :);
w = 34.^(0:ns - 1).';
key = rows * w;
rot = zeros(size(rows, 1), 8);
for k = 0:7
    R = sort(reshape(geom.perm(k + 1, rows(:)), size(rows)), 2);
    [tf, loc] = ismember(R * w, key);
    if ~all(tf), error('d3:closure', 'The admissible set is not closed under the 45 deg rotation.'); end
    rot(:, k + 1) = loc;
end
L = struct('ns', ns, 'rows', rows, 'key', key, 'rot', rot, 'nAll', size(lay, 1), 'n', size(rows, 1));
end

function a = copyAngles(cfg, extra)
if nargin < 2, extra = []; end
d = unique(round(mod(cfg.D3.headingsDeg, 45) * 1e9) / 1e9);
s = cfg.D3.headingStepDeg;
a = unique([0, d, d + s, d - s, extra(:).']);
end

%% ================================================================================================
function F = buildField(cfg, p, geom, angles, cacheDir, frfFor, verbose)
if nargin < 6, frfFor = []; end
if nargin < 7, verbose = true; end
nC = size(geom.sensorsND, 1);
pts = zeros(nC * numel(angles), 2);
for i = 1:numel(angles)
    pts((i - 1) * nC + (1:nC), :) = [geom.sensorsND(:, 1), geom.sensorsND(:, 2) - angles(i) * pi / 180];
end
args = {'ThetaIdx', [1 3], 'NodeSpacing', cfg.emm.nodeSpacing, 'CacheDir', cacheDir, 'Step', cfg.emm.jacobianStep, ...
    'Verbose', verbose};
if ~isempty(frfFor), args = [args, {'FRF', frfFor(pts)}]; end
tc = tic;
fld = emmField(p, pts, args{:});
[omega, kBins] = welchBins(cfg.welch.L, cfg.acq.dt, cfg.band.estimate, cfg.welch.binStep, cfg.band.basin);
sea = cfg.sea;
% record statistics (fine grid, spectrum, taper) from a two-sensor view of the prepared twin, so no 330-channel
% record is synthesised; then the FRF of EVERY point on that fine grid (info.H = twin.Hacc(info.omega), as
% synthesiseTwinRecords builds it)
tw2 = fld.twin; tw2.sensors = fld.twin.sensors(1:2, :); tw2.Hacc = @(w) firstRows(fld.twin.Hacc(w));
[~, ~, info] = synthesiseTwinRecords(tw2, cfg.acq.Nref, cfg.acq.dt, 'Seed', 1, 'Spectrum', 'jonswap', ...
    'Hs', sea.Hs, 'PeakFrequency', sea.omegaP, 'PeakEnhancement', sea.gammaJ);
info.H = fld.twin.Hacc(info.omega);
[nEff, K] = welchNEff(cfg.acq.Nref, cfg.welch.L, cfg.welch.overlap);
F = struct('angles', angles, 'nC', nC, 'pts', pts, 'omega', omega, 'kBins', kBins, 'H', fld.H(omega), ...
    'dH', fld.dH(omega), 'info', info, 'nEff', nEff, 'K', K, 'L', cfg.welch.L, 'perm', geom.perm, ...
    'dpsi', cfg.D3.headingStepDeg * pi / 180, 'seconds', toc(tc), 'cacheStatus', {fld.cacheStatus});
end

function H = firstRows(H)
H = H(1:2, :);
end

function r = rowsOf(F, aDeg)
i = find(abs(F.angles - aDeg) < 1e-9, 1);
if isempty(i), error('d3:angle', 'Grid copy at %.6g deg was not sampled.', aDeg); end
r = (i - 1) * F.nC + (1:F.nC);
end

function r = headingRows(F, hDeg)
h = mod(hDeg, 360);
d = mod(h, 45); if d > 45 - 1e-9, d = 0; end
k = mod(round((h - d) / 45), 8);
base = rowsOf(F, round(d * 1e9) / 1e9);
r = base(F.perm(k + 1, :));
end

%% ================================================================================================
function pre = precomputeA(F, deltaDeg, Snn, withPsi)
rows0 = rowsOf(F, deltaDeg);
H = F.H(rows0, :); dH = F.dH(rows0, :, :);
if withPsi
    dHp = (F.H(rowsOf(F, deltaDeg + F.dpsi * 180 / pi), :) - F.H(rowsOf(F, deltaDeg - F.dpsi * 180 / pi), :)) / (2 * F.dpsi);
    dH = cat(3, dH, dHp);
end
Ssm = welchSmoothedSpectrum(F.omega, F.info, F.L, rows0);
infoC = F.info; infoC.H = F.info.H(rows0, :);
pre = repmat(struct('ref', 0, 'valid', [], 'B', [], 'dT', [], 'bW', [], 'fullPD', true, 'nValid', 0), 1, F.nC);
for r = 1:F.nC
    others = setdiff(1:F.nC, r);
    cv = analyticTransmissibilityCovariance(Ssm, [], Snn, r, F.nEff);
    TW = welchExpectedTransmissibility(F.omega, infoC, F.L, r, others, 0);
    pre(r) = fill(r, cv.blocks, cv.valid, H, dH, others, TW - H(others, :) ./ H(r, :), []);
end
end

function c = seaComponents(F, deltaDeg, sea)
switch sea.kind
    case 'two'
        h = [deltaDeg, deltaDeg + sea.phiDeg]; w = [1 - sea.f2, sea.f2];
    case 'spread'
        V2 = v2Tools();
        [hr, w] = V2.spreadWeights(deltaDeg * pi / 180, sea.s, sea.nDir);
        h = round(mod(hr * 180 / pi, 360) * 1e9) / 1e9;
    otherwise
        error('d3:sea', 'Unknown D3b sea kind ''%s''.', sea.kind);
end
c = struct('rows', {}, 'fraction', {}, 'heading', {});
for i = 1:numel(h)
    c(end + 1) = struct('rows', headingRows(F, h(i)), 'fraction', w(i), 'heading', h(i) * pi / 180); %#ok<AGROW>
end
end

function pre = precomputeB(F, deltaDeg, sea, Snn, ruler)
rows0 = rowsOf(F, deltaDeg);
H = F.H(rows0, :); dH = F.dH(rows0, :, :);
src = struct('omega', F.omega, 'info', F.info, 'L', F.L);
D = directionalSpectralMatrix(src, seaComponents(F, deltaDeg, sea), Snn);
Ssm = welchSmoothedSpectrum(F.omega, F.info, F.L, rows0);
infoC = F.info; infoC.H = F.info.H(rows0, :);
pre = repmat(struct('ref', 0, 'valid', [], 'B', [], 'dT', [], 'bW', [], 'fullPD', true, 'nValid', 0, 'bD', []), ...
    1, F.nC);
for r = 1:F.nC
    others = setdiff(1:F.nC, r);
    cvD = analyticTransmissibilityCovariance(D.Ssig, [], Snn, r, F.nEff);
    TW = welchExpectedTransmissibility(F.omega, infoC, F.L, r, others, 0);
    if strcmp(ruler, 'directional')
        B = cvD.blocks; v = cvD.valid;
    else                                               % 'nominal': the V2 ruler (tester cross-check)
        cvS = analyticTransmissibilityCovariance(Ssm, [], Snn, r, F.nEff);
        B = cvS.blocks; v = cvS.valid;
    end
    pre(r) = fill(r, B, v, H, dH, others, TW - H(others, :) ./ H(r, :), cvD.T - TW);
end
end

function s = fill(r, blocks, v, H, dH, others, bW, bD)
Hr = H(r, :); Ho = H(others, :);
nT = size(dH, 3);
dT = complex(zeros(numel(others), size(H, 2), nT));
for i = 1:nT
    dT(:, :, i) = (dH(others, :, i) .* Hr - Ho .* dH(r, :, i)) ./ Hr.^2;
end
B = blocks(:, :, v);
pd = true;
for k = 1:size(B, 3)
    [~, flag] = chol(B(:, :, k)); if flag ~= 0, pd = false; break; end
end
s = struct('ref', r, 'valid', v, 'B', B, 'dT', dT(:, v, :), 'bW', bW(:, v), 'fullPD', pd, 'nValid', nnz(v));
if ~isempty(bD), s.bD = bD(:, v); end
end

%% ================================================================================================
function [Fm, G] = layoutFisher(P, layout, ref, biasFields)
% F = sum_k J_k' Sigma_k^-1 J_k and G = sum_k J_k' Sigma_k^-1 [b_1 b_2 ...] for the layout (sub-block of P = pre(ref))
o = layout(layout ~= ref);
pos = o - (o > ref); m = numel(pos); n32 = size(P.dT, 1);
idx = [pos(:); n32 + pos(:)];
Bs = P.B(idx, idx, :);
nT = size(P.dT, 3);
Jp = [permute(real(P.dT(pos, :, :)), [1 3 2]); permute(imag(P.dT(pos, :, :)), [1 3 2])];   % 2m x nT x Kv
nb = numel(biasFields);
bp = zeros(2 * m, nb, size(Bs, 3));
for j = 1:nb
    b = P.(biasFields{j});
    bp(:, j, :) = [permute(real(b(pos, :)), [1 3 2]); permute(imag(b(pos, :)), [1 3 2])];
end
Kv = size(Bs, 3);
if Kv == 0, Fm = zeros(nT); G = zeros(nT, nb); return; end
if exist('pagemldivide', 'builtin') || exist('pagemldivide', 'file')
    X = pagemldivide(Bs, [Jp, bp]);
    JT = pagetranspose(Jp);
    Y = sum(pagemtimes(JT, X), 3);
else
    Y = zeros(nT, nT + nb);
    for k = 1:Kv
        Y = Y + Jp(:, :, k).' * (Bs(:, :, k) \ [Jp(:, :, k), bp(:, :, k)]);
    end
end
Fm = (Y(:, 1:nT) + Y(:, 1:nT).') / 2;
G = Y(:, nT + 1:end);
end

function q = metricsA(F3, g, cfg)
% known heading: the (beta, R) block; heading nuisance: the Schur complement (F_pp = 0 exactly: no penalty)
Fk = F3(1:2, 1:2);
q.known = core(Fk, g(1:2), [], cfg);
if size(F3, 1) < 3
    q.nuis = q.known; return;
end
if F3(3, 3) > 0
    Fe = Fk - F3(1:2, 3) * F3(3, 1:2) / F3(3, 3);
    Fe = (Fe + Fe.') / 2;
else
    Fe = Fk;
end
dth = [];
if F3(3, 3) > 0 && isfinite(rcond(F3)) && rcond(F3) >= 1e-14
    dth = F3 \ g; dth = dth(1:2);                    % bias of (beta, R) with psi estimated
end
q.nuis = core(Fe, g(1:2), dth, cfg);                % dth [] (psi unidentifiable): the known-heading bias
end

function q = metricsB(F2, G2, cfg)
c = core(F2, sum(G2, 2), [], cfg);                   % b_W + b_dir: d_tot
q = c;
q.dTot = c.d;
q.dSys = NaN; q.dthetaSys = [NaN; NaN];
if ~c.singular
    dd = F2 \ G2(:, end);
    q.dSys = sqrt(max(dd.' * F2 * dd, 0));
    q.dthetaSys = dd;
end
end

function q = core(Fm, g, dth, cfg)
% metrics of an information matrix Fm with the bias dth (or Fm \ g when dth is []); no solve if Fm is singular
sing = ~all(isfinite(Fm(:))) || rcond(Fm) < 1e-12;
if sing
    q = struct('A', Inf, 'sigma', [Inf; Inf], 'd', NaN, 'kappa', Inf, 'singular', true, 'dtheta', [NaN; NaN]);
else
    if isempty(dth), dth = Fm \ g; end
    C = Fm \ eye(2);
    q = struct('A', trace(C), 'sigma', sqrt(diag(C)), 'd', sqrt(max(dth.' * Fm * dth, 0)), 'kappa', cond(Fm), ...
        'singular', false, 'dtheta', dth);
end
q.classCode = classOf(q.sigma, q.d, q.kappa, q.singular, cfg);
end

function c = classOf(sigma, d, kappa, singular, cfg)
I = cfg.class.identifiable; M = cfg.class.marginal;
if singular || any(~isfinite(sigma)), c = 3; return; end
s = max(sigma);
if s <= I.maxSigma && d <= I.maxBiasDistance && kappa <= I.maxCondF, c = 1;
elseif s <= M.maxSigma && d <= M.maxBiasDistance && kappa <= M.maxCondF, c = 2;
else, c = 3;
end
end

%% ================================================================================================
function S = searchA(pre, L, cfg, withPsi, label)
if nargin < 5, label = ''; end
nL = L.n; ns = L.ns;
nanv = NaN(nL, 1);
S = struct('An', Inf(nL, 1), 'sBn', nanv, 'sRn', nanv, 'dWn', nanv, 'kapN', nanv, 'clsN', 3 * ones(nL, 1), ...
    'refN', zeros(nL, 1), 'Ak', Inf(nL, 1), 'sBk', nanv, 'sRk', nanv, 'dWk', nanv, 'kapK', nanv, ...
    'clsK', 3 * ones(nL, 1), 'refK', zeros(nL, 1), 'withPsi', withPsi);
t0 = tic; nextPct = 25;
for i = 1:nL
    lay = L.rows(i, :);
    for j = 1:ns
        r = lay(j);
        [Fm, G] = layoutFisher(pre(r), lay, r, {'bW'});
        q = metricsA(Fm, G, cfg);
        if q.known.A < S.Ak(i) || S.refK(i) == 0
            S.Ak(i) = q.known.A; S.sBk(i) = q.known.sigma(1); S.sRk(i) = q.known.sigma(2); S.dWk(i) = q.known.d;
            S.kapK(i) = q.known.kappa; S.clsK(i) = q.known.classCode; S.refK(i) = r;
        end
        if q.nuis.A < S.An(i) || S.refN(i) == 0
            S.An(i) = q.nuis.A; S.sBn(i) = q.nuis.sigma(1); S.sRn(i) = q.nuis.sigma(2); S.dWn(i) = q.nuis.d;
            S.kapN(i) = q.nuis.kappa; S.clsN(i) = q.nuis.classCode; S.refN(i) = r;
        end
    end
    if 100 * i / nL >= nextPct
        fprintf('      %s n_s = %d: %3.0f%% of %d layouts (%.0f s)\n', label, ns, 100 * i / nL, nL, toc(t0));
        if exist('OCTAVE_VERSION', 'builtin'), fflush(stdout); else, drawnow; end
        nextPct = nextPct + 25;
    end
end
S.seconds = toc(t0);
end

function S = searchB(pre, L, cfg, label)
% per layout and condition, three reference choices (the reference is a processing choice):
%   A-best        the precision-optimal reference (what an analyst choosing by sigma would use): A, sB, dTot, ...
%   feasible-best the smallest A among references with d_tot <= cfg.D3.biasGate (Af = Inf if none): the gate test
%   least-biased  the smallest d_tot over references (dTotMin): the least-bad diagnostic
% Choosing the reference with knowledge of the bias is optimistic, so "no feasible layout" is a robust conclusion.
if nargin < 4, label = ''; end
nL = L.n; ns = L.ns; gate = cfg.D3.biasGate;
nanv = NaN(nL, 1);
S = struct('A', Inf(nL, 1), 'sB', nanv, 'sR', nanv, 'dSys', nanv, 'dTot', nanv, 'kap', nanv, 'cls', 3 * ones(nL, 1), ...
    'ref', zeros(nL, 1), 'Af', Inf(nL, 1), 'sBf', Inf(nL, 1), 'dTotf', nanv, 'clsf', 3 * ones(nL, 1), 'reff', zeros(nL, 1), ...
    'dTotMin', Inf(nL, 1), 'refMinD', zeros(nL, 1));
t0 = tic;
for i = 1:nL
    lay = L.rows(i, :);
    for j = 1:ns
        r = lay(j);
        [Fm, G] = layoutFisher(pre(r), lay, r, {'bW', 'bD'});
        q = metricsB(Fm, G, cfg);
        if q.A < S.A(i) || S.ref(i) == 0
            S.A(i) = q.A; S.sB(i) = q.sigma(1); S.sR(i) = q.sigma(2); S.dSys(i) = q.dSys; S.dTot(i) = q.dTot;
            S.kap(i) = q.kappa; S.cls(i) = q.classCode; S.ref(i) = r;
        end
        if ~isnan(q.dTot)
            if q.dTot < S.dTotMin(i), S.dTotMin(i) = q.dTot; S.refMinD(i) = r; end
            if q.dTot <= gate && q.A < S.Af(i)
                S.Af(i) = q.A; S.sBf(i) = q.sigma(1); S.dTotf(i) = q.dTot; S.clsf(i) = q.classCode; S.reff(i) = r;
            end
        end
    end
end
S.seconds = toc(t0);
if ~isempty(label), fprintf('      %s n_s = %d: %d layouts (%.0f s)\n', label, ns, nL, S.seconds); end
end

function M = headingTable(perDelta, field, L, F, hDeg)
% perDelta: struct array with fields deltaDeg and S (a search result); M(:, j) = field of each layout at heading hDeg(j)
M = NaN(L.n, numel(hDeg));
for j = 1:numel(hDeg)
    h = mod(hDeg(j), 360); d = mod(h, 45); k = mod(round((h - d) / 45), 8);
    c = find(abs([perDelta.deltaDeg] - d) < 1e-9, 1);
    if isempty(c), error('d3:heading', 'No search at delta = %.6g deg.', d); end
    v = perDelta(c).S.(field);
    M(:, j) = v(L.rot(:, k + 1));
end
end