function T = d1LayoutTools()
%D1LAYOUTTOOLS  The numerical core of study D1 (sensor count, layout and reference), shared by
%   runD1SensorLayout and testD1SensorLayout so both exercise exactly the same code.
%   T = d1LayoutTools();   then e.g.  geom = T.buildGeometry(cfg, p)
%
%   DESIGN. The expensive physics is computed ONCE for all 33 candidate points (one emmField: nominal,
%   ln beta +- h, ln R +- h, sensors fixed at physical positions). For each of the 33 possible
%   references the transmissibilities to the other 32 points, their Jacobian, the Welch bias and the
%   FULL 64 x 64 real covariance block per bin (shared-reference correlation kept) are precomputed.
%   A layout with reference r is then only a sub-selection of rows/columns of pre(r): the covariance
%   of T_j for j in the layout depends only on S_jl, S_jr and S_rr, so the sub-block of the 32-sensor
%   block is EXACTLY the covariance of the smaller layout. No EMM solve happens inside the search.
%
%   FUNCTIONS
%   geom = T.buildGeometry(cfg, p)
%       33-point candidate grid at the true R: geom.frac [r/R theta], geom.sensorsND (EMM scaling,
%       r/H), geom.xy (m), geom.rPhys (m), geom.Rphys, geom.dist (33 x 33, m), geom.names {'s1'..},
%       geom.level2C (grid indices of the Level 2C sensors, found by matching coordinates).
%   [layouts, ok] = T.enumerateLayouts(geom, ns, cfg, spacing)
%       all nchoosek(33, ns) layouts (rows, ascending grid index) and the admissibility mask:
%       radius <= maxRadiusFrac R, buoy wholly on the floe (r + d/2 <= R), centre spacing >= spacing.
%   full = T.buildFullField(cfg, p, geom, cacheDir, verbose)
%       emmField at all 33 points, the retained bins, H, dH, the twin record info, N_eff and the
%       Welch-smoothed 33 x 33 signal spectral matrix (as evaluateScenario builds them).
%   pre = T.precompute(full, Snn, cfg)
%       33-element struct array, one per reference (see DESIGN), accepted bins only.
%   q = T.evaluateLayout(pre, layout, ref, cfg, mode)
%       metrics of one layout (grid indices, any order) with reference ref (a GRID index in layout).
%       mode 'fast' (paged solves, default) or 'direct' (fisherFromBlocks on the extracted J, blocks
%       and bias: the E1 engine path itself). q: F, sigma [ln beta; ln R], A = tr(F^-1), D = log det F,
%       E = lambda_min(F), kappa, dW, class, classCode (1 Identifiable, 2 Marginal, 3 Not), singular.
%       Singular F (rcond < 1e-12): A = Inf, D = -Inf, E = 0, Not identifiable; never regularised.
%   S = T.searchLayouts(pre, layouts, cfg, label)
%       every layout x every reference in it: S.Aall, sBall, sRall, dWall, kappaAll, classAll (nL x ns,
%       column = position of the reference in the layout), and per layout the A-best reference with
%       its sigma, dW, kappa, class, plus the D-best and E-best values and references (each criterion
%       picks its own reference).
%   G = T.gramianLayouts(cfg, nsList)
%       the Level 2C Gramian lambda_min optimum per n_s, recomputed exactly as level2C_finalCheck
%       (two known modes at 2.70 and 6.00 rad/s, dt 0.02 s, N = 4001 instants, exhaustive search),
%       with all numerically tied optima.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
T.buildGeometry = @buildGeometry;
T.enumerateLayouts = @enumerateLayouts;
T.buildFullField = @buildFullField;
T.precompute = @precompute;
T.evaluateLayout = @evaluateLayout;
T.searchLayouts = @searchLayouts;
T.gramianLayouts = @gramianLayouts;
T.classCode = @classCode;
end

%% ================================================================================================
function geom = buildGeometry(cfg, p)
frac = cfg.layout.gridFrac;
nC = size(frac, 1);
H = cfg.physics.Hbasin;
geom.frac = frac;
geom.sensorsND = [frac(:, 1) * p(3), frac(:, 2)];        % r/H, fixed physical points at the true R
geom.Rphys = p(3) * H;
geom.rPhys = frac(:, 1) * geom.Rphys;
geom.xy = [geom.rPhys .* cos(frac(:, 2)), geom.rPhys .* sin(frac(:, 2))];
d = zeros(nC);
for i = 1:nC
    d(:, i) = sqrt(sum((geom.xy - geom.xy(i, :)).^2, 2));
end
geom.dist = d;
geom.names = arrayfun(@(i) sprintf('s%d', i), 1:nC, 'UniformOutput', false);
L2 = cfg.layout.level2C.frac;
geom.level2C = zeros(1, size(L2, 1));
for k = 1:size(L2, 1)
    m = find(abs(frac(:, 1) - L2(k, 1)) < 1e-12 & abs(angDiff(frac(:, 2), L2(k, 2))) < 1e-12);
    if numel(m) ~= 1, error('d1:level2C', 'Level 2C sensor %d not found uniquely on the grid.', k); end
    geom.level2C(k) = m;
end
geom.level2CRef = geom.level2C(cfg.layout.level2C.ref);
end

function a = angDiff(a, b)
a = mod(a - b + pi, 2 * pi) - pi;
end

%% ================================================================================================
function [layouts, ok] = enumerateLayouts(geom, ns, cfg, spacing)
nC = size(geom.xy, 1);
layouts = nchoosek(1:nC, ns);
pointOk = geom.frac(:, 1) <= cfg.D1.maxRadiusFrac + 1e-12 & ...
    geom.rPhys + cfg.D1.edgeClearance <= geom.Rphys + 1e-12;
ok = all(pointOk(layouts), 2);
for a = 1:ns - 1
    for b = a + 1:ns
        ok = ok & geom.dist(sub2ind([nC nC], layouts(:, a), layouts(:, b))) >= spacing - 1e-12;
    end
end
end

%% ================================================================================================
function full = buildFullField(cfg, p, geom, cacheDir, verbose, frf)
% frf (optional): analytic FRF handle @(omega, p) -> 33 x nW, replacing the EMM (tests only)
if nargin < 5, verbose = true; end
if nargin < 6, frf = []; end
fld = emmField(p, geom.sensorsND, 'ThetaIdx', [1 3], 'NodeSpacing', cfg.emm.nodeSpacing, ...
    'CacheDir', cacheDir, 'Step', cfg.emm.jacobianStep, 'Verbose', verbose, 'FRF', frf);
[omega, kBins] = welchBins(cfg.welch.L, cfg.acq.dt, cfg.band.estimate, cfg.welch.binStep, cfg.band.basin);
sea = cfg.sea;
[~, ~, info] = synthesiseTwinRecords(fld.twin, cfg.acq.Nref, cfg.acq.dt, 'Seed', 1, 'Spectrum', 'jonswap', ...
    'Hs', sea.Hs, 'PeakFrequency', sea.omegaP, 'PeakEnhancement', sea.gammaJ);
[nEff, K] = welchNEff(cfg.acq.Nref, cfg.welch.L, cfg.welch.overlap);
nC = size(geom.sensorsND, 1);
full = struct('field', fld, 'omega', omega, 'kBins', kBins, 'H', fld.H(omega), 'dH', fld.dH(omega), ...
    'info', info, 'nEff', nEff, 'K', K, 'L', cfg.welch.L, ...
    'Ssm', welchSmoothedSpectrum(omega, info, cfg.welch.L, 1:nC), 'nC', nC);
end

%% ================================================================================================
function pre = precompute(full, Snn, cfg) %#ok<INUSD>
nC = full.nC;
pre = repmat(struct('ref', 0, 'others', [], 'valid', [], 'B', [], 'dT', [], 'bW', [], 'fullPD', true, ...
    'nValid', 0), 1, nC);
for r = 1:nC
    others = setdiff(1:nC, r);
    cv = analyticTransmissibilityCovariance(full.Ssm, [], Snn, r, full.nEff);
    v = cv.valid;
    Hr = full.H(r, :); Ho = full.H(others, :);
    dT = complex(zeros(numel(others), numel(full.omega), 2));
    for i = 1:2
        dT(:, :, i) = (full.dH(others, :, i) .* Hr - Ho .* full.dH(r, :, i)) ./ Hr.^2;
    end
    TW = welchExpectedTransmissibility(full.omega, full.info, full.L, r, others, 0);
    bW = TW - Ho ./ Hr;                                % Welch-expected minus exact, per transmissibility
    B = cv.blocks(:, :, v);
    pd = true;
    for k = 1:size(B, 3)
        [~, flag] = chol(B(:, :, k)); if flag ~= 0, pd = false; break; end
    end
    pre(r).ref = r; pre(r).others = others; pre(r).valid = v; pre(r).nValid = nnz(v);
    pre(r).B = B; pre(r).dT = dT(:, v, :); pre(r).bW = bW(:, v); pre(r).fullPD = pd;
end
end

%% ================================================================================================
function q = evaluateLayout(pre, layout, ref, cfg, mode)
if nargin < 5 || isempty(mode), mode = 'fast'; end
P = pre(ref);
o = layout(layout ~= ref);
if numel(o) ~= numel(layout) - 1, error('d1:ref', 'The reference must be one of the layout sensors.'); end
pos = o - (o > ref);                                % position of each sensor in P.others
m = numel(pos); n32 = numel(P.others);
idx = [pos(:); n32 + pos(:)];
Bs = P.B(idx, idx, :);
Jp = [permute(real(P.dT(pos, :, :)), [1 3 2]); permute(imag(P.dT(pos, :, :)), [1 3 2])];   % 2m x 2 x Kv
bp = [permute(real(P.bW(pos, :)), [1 3 2]); permute(imag(P.bW(pos, :)), [1 3 2])];         % 2m x 1 x Kv
Kv = size(Bs, 3);
useFast = strcmp(mode, 'fast') && P.fullPD && Kv > 0;
if useFast
    if exist('pagemldivide', 'builtin') || exist('pagemldivide', 'file')
        X = pagemldivide(Bs, [Jp, bp]);
        JT = pagetranspose(Jp);
        F = sum(pagemtimes(JT, X(:, 1:2, :)), 3);
        g = sum(pagemtimes(JT, X(:, 3, :)), 3);
    else                                            % Octave / older MATLAB
        F = zeros(2); g = zeros(2, 1);
        for k = 1:Kv
            X = Bs(:, :, k) \ [Jp(:, :, k), bp(:, :, k)];
            F = F + Jp(:, :, k).' * X(:, 1:2);
            g = g + Jp(:, :, k).' * X(:, 3);
        end
    end
    F = (F + F.') / 2;
    singular = Kv == 0 || ~all(isfinite(F(:))) || rcond(F) < 1e-12;
    nUsed = Kv;
else                                                % the E1 engine path itself
    Jf = reshape(permute(Jp, [1 3 2]), 2 * m * Kv, 2);
    bf = reshape(bp, 2 * m * Kv, 1);
    qq = fisherFromBlocks(Jf, Bs, bf);
    F = qq.F; singular = qq.singular; nUsed = qq.nUsed;
    if ~singular, g = F * qq.dtheta; end
end
q = struct('F', F, 'singular', singular, 'nUsed', nUsed, 'layout', layout, 'ref', ref);
if singular
    q.sigma = [Inf; Inf]; q.A = Inf; q.D = -Inf; q.E = 0; q.kappa = Inf; q.dW = NaN; q.corr = NaN;
else
    C = F \ eye(2);
    q.sigma = sqrt(diag(C)); q.A = trace(C); q.D = log(det(F)); q.E = min(eig(F)); q.kappa = cond(F);
    q.corr = C(1, 2) / sqrt(C(1, 1) * C(2, 2));
    dth = F \ g;
    q.dW = sqrt(dth.' * F * dth);
end
q.classCode = classCode(q, cfg);
q.class = cfg.class.names{q.classCode};
end

function c = classCode(q, cfg)
I = cfg.class.identifiable; M = cfg.class.marginal;
if q.singular || ~all(isfinite(q.sigma)), c = 3; return; end
s = max(q.sigma);
if s <= I.maxSigma && q.dW <= I.maxBiasDistance && q.kappa <= I.maxCondF, c = 1;
elseif s <= M.maxSigma && q.dW <= M.maxBiasDistance && q.kappa <= M.maxCondF, c = 2;
else, c = 3;
end
end

%% ================================================================================================
function S = searchLayouts(pre, layouts, cfg, label)
if nargin < 4, label = ''; end
[nL, ns] = size(layouts);
S.layouts = uint8(layouts);
S.Aall = Inf(nL, ns); S.sBall = NaN(nL, ns); S.sRall = NaN(nL, ns); S.dWall = NaN(nL, ns);
S.kappaAll = NaN(nL, ns); S.classAll = zeros(nL, ns, 'uint8');
nanv = NaN(nL, 1);
S.refA = zeros(nL, 1, 'uint8'); S.A = Inf(nL, 1); S.sB = nanv; S.sR = nanv; S.dW = nanv; S.kappa = nanv;
S.corr = nanv; S.classCode = zeros(nL, 1, 'uint8');
S.refD = zeros(nL, 1, 'uint8'); S.D = -Inf(nL, 1); S.refE = zeros(nL, 1, 'uint8'); S.E = zeros(nL, 1);
t0 = tic; nextPct = 10;
for i = 1:nL
    lay = layouts(i, :);
    for k = 1:ns
        q = evaluateLayout(pre, lay, lay(k), cfg, 'fast');
        S.Aall(i, k) = q.A; S.sBall(i, k) = q.sigma(1); S.sRall(i, k) = q.sigma(2); S.dWall(i, k) = q.dW;
        S.kappaAll(i, k) = q.kappa; S.classAll(i, k) = q.classCode;
        if q.A < S.A(i) || (k == 1 && isinf(q.A))
            S.A(i) = q.A; S.refA(i) = lay(k); S.sB(i) = q.sigma(1); S.sR(i) = q.sigma(2);
            S.dW(i) = q.dW; S.kappa(i) = q.kappa; S.corr(i) = q.corr; S.classCode(i) = q.classCode;
        end
        % k == 1 initialises, so a layout singular for every reference (D = -Inf, E = 0 throughout)
        % still gets a valid reference rather than 0
        if k == 1 || q.D > S.D(i), S.D(i) = q.D; S.refD(i) = lay(k); end
        if k == 1 || q.E > S.E(i), S.E(i) = q.E; S.refE(i) = lay(k); end
    end
    if 100 * i / nL >= nextPct
        fprintf('    %s n_s = %d: %3.0f%% of %d layouts (%.0f s)\n', label, ns, 100 * i / nL, nL, toc(t0));
        if exist('OCTAVE_VERSION', 'builtin'), fflush(stdout); else, drawnow; end
        nextPct = nextPct + 10;
    end
end
S.seconds = toc(t0);
end

%% ================================================================================================
function G = gramianLayouts(cfg, nsList)
% exactly the Level 2C setup of level2C_finalCheck.m (Check 2: lambda_min at N = 4001)
p = cfg.p0.vec; H = cfg.physics.Hbasin; g = cfg.physics.g;
beta = p(1); gamma = p(2); R = p(3); nu = 0.3; M = 50; P = 10; Ntr = 10; dt = 0.02;
omega1 = 2.70; omega2 = 6.00;
A = buildKnownOscillatorA(omega1, omega2, dt);
b1 = precomputeDeflectionData(H * omega1^2 / g, beta, gamma, R, nu, M, P, Ntr);
b2 = precomputeDeflectionData(H * omega2^2 / g, beta, gamma, R, nu, M, P, Ntr);
loc = [cfg.layout.gridFrac(:, 1) * R, cfg.layout.gridFrac(:, 2)];
nC = size(loc, 1);
phi1 = zeros(nC, 1); phi2 = phi1;
for c = 1:nC
    phi1(c) = evaluateDeflection(b1, loc(c, 1), loc(c, 2));
    phi2(c) = evaluateDeflection(b2, loc(c, 1), loc(c, 2));
end
phi1 = phi1 / norm(phi1); phi2 = phi2 / norm(phi2);
rows = zeros(nC, 4);
for c = 1:nC, rows(c, :) = buildCandidateCrows(phi1(c), phi2(c)); end
Nrec = numel(0:dt:80);                               % 4001 measurement instants
lambdaMinFn = @(Wo) gramianMetrics(Wo).lambdaMin;
G = struct('ns', num2cell(nsList), 'layout', [], 'score', [], 'tied', [], 'N', Nrec);
for k = 1:numel(nsList)
    [idx, best, allS] = exhaustiveSensorSelection(rows, A, Nrec, lambdaMinFn, nsList(k));
    combos = nchoosek(1:nC, nsList(k));
    tol = 1e-10 * max(1, abs(best));
    G(k).layout = sort(idx); G(k).score = best;
    G(k).tied = combos(abs(allS - best) <= tol, :);
end
end