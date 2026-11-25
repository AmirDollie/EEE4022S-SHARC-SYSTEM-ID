%% diag_F1_Mo_ladder.m  --  F1: open-water tower to very large Mo, then a clean extrapolation in M
%
% CONTEXT
%   diag_F1_floor.m showed the ~1e-4 "floor" is a smooth beat in the OPEN-WATER tower (period ~1/gamma
%   ~ 690 modes at p0), while the PLATE tower converges cleanly at p ~ 1.74 up to M = 500. The two
%   towers couple only through the P+1 Gegenbauer coefficients, so they can be truncated separately.
%
% WHY LARGE Mo IS CHEAP
%   The open-water modes are orthogonal (D0 diagonal), so MIn and M0 are elementwise and the whole
%   open-water tower enters the solve through one (P+1)x(P+1) sum,
%       G0(Mo) = sum_{m=0..Mo} D0g(m,:).' * D0g(m,:) / (d0_m E0'_m),
%   plus m = 0 terms for the incident wave. A cumulative sum gives G0 at EVERY Mo in one pass, so the
%   cost is linear in Mo and no (Mo x Mo) matrix is ever formed. The plate side (dense, ~M^3) is built
%   once with the production functions. Open-water Bessel ratios use the exponentially scaled
%   besselk/besseli (no overflow at k_m R ~ 5e4), and open-water roots come from one bracket per root.
%   All of this lives in this file; nothing in Forward Model is modified.
%
% WHAT IT RUNS  (reference floe p0, P = 10, N = 10, sensors s2, s6, s14 against s26)
%   Part 0  Regression: the reduced solver against the production pipeline at Mo = M = 50 and 200
%           (8.5 rad/s, same roots). Must agree to ~1e-12. Also: fast open-water roots vs dispersionRoots
%           at Mo = 500, and the d0 identity (tan k = -alpha/k at a root) against direct tan.
%   Part 1  Mo sweep 100..MoMax (dense to 2000, log-spaced beyond) at plate M = 100 and 200, at 8.5 and
%           5.0 rad/s. Error relative to Mo = MoMax, its envelope over one beat period (1/gamma), the
%           envelope decay rate, and the estimated open-water error left at MoMax.
%   Part 2  Plate M ladder 50, 100, 200, 400 at Mo = MoMax: rates p from three levels, two Richardson
%           extrapolates and their difference (the extrapolation uncertainty), and the error of each
%           truncation (production M = 50 included) against the extrapolate.
%
% RUNTIME  A few minutes (quick mode under a minute).
% WRITES   results/diag_F1_Mo_ladder.log, .mat, .png
% HOW      In MATLAB from the repo root:   run('diagnostics/diag_F1_Mo_ladder.m')
%          Quick smoke test:               MOLAD_QUICK = true; run('diagnostics/diag_F1_Mo_ladder.m')
%          (clear MOLAD_QUICK before the full run.)

%% ---- setup --------------------------------------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
root = fileparts(here);
addpath(fullfile(root, 'Forward Model'));
addpath(fullfile(root, 'Forward Model', 'Animation'));
resDir = fullfile(root, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
logFile = fullfile(resDir, 'diag_F1_Mo_ladder.log');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile); diary('on');
cleanupDiary = onCleanup(@() diary('off'));

QUICK = exist('MOLAD_QUICK', 'var') && MOLAD_QUICK;
fprintf('==== diag_F1_Mo_ladder  (%s)%s ====\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ternary(QUICK, '  [QUICK MODE]', ''));

bet0 = 4.6985e-5; gam0 = 1.4548e-3; R = 0.3830; nu = 0.3; depth = 1.88; g = 9.81;
sensRR = [0.3 0; 0.3 pi; 0.5 pi; 0.9 0];
sensR = sensRR(:, 1).' * R; sensTh = sensRR(:, 2).';
N = 10; P = 10;

omegaList = [8.5 5.0];
MoMax = 40000;
MoGrid = unique([100:10:2000, round(logspace(log10(2000), log10(MoMax), 400))]);
Mfix = [100 200];
Mladder = [50 100 200 400];
regM = [50 200];
if QUICK
    omegaList = 8.5; MoMax = 5000;
    MoGrid = unique([100:50:1000, round(logspace(3, log10(MoMax), 60))]);
    Mfix = 100; Mladder = [50 100 200]; regM = 50;
end
period = 1 / gam0;
fprintf('P = %d, N = %d, MoMax = %d (%d grid points), plate M ladder %s, beat period 1/gamma = %.0f modes\n\n', ...
    P, N, MoMax, numel(MoGrid), mat2str(Mladder), period);

%% ---- Part 0: regression of the reduced solver ---------------------------------------------------------------------
fprintf('PART 0  Reduced solver vs production pipeline (8.5 rad/s, identical roots, d0 from tan as in production)\n');
alpha = depth * 8.5^2 / g;
Mr = max(regM);
kapR = dispersionRoots(alpha, bet0, gam0, Mr);
kOWR = dispersionRoots(alpha, 0, 0, Mr);
for M = regM
    [Tp, bap] = productionSolve(alpha, bet0, gam0, R, nu, P, N, kapR(1:M + 3), kOWR(1:M + 1), sensR, sensTh);
    PL = plateSide(alpha, bet0, gam0, R, nu, P, N, kapR(1:M + 3), sensR, sensTh);
    O = openSide(alpha, gam0, R, P, N, kOWR(1:M + 1), M, 'tan');
    [Tr, baR] = combine(PL, O, N, 1);
    fprintf('  Mo = M = %3d:  |T_red - T_prod|/|T| = %.2e   |b/a_red - b/a_prod|/|b/a| = %.2e\n', ...
        M, norm(Tr - Tp) / norm(Tp), norm(baR - bap) / norm(bap));
    if M == Mr
        Oid = openSide(alpha, gam0, R, P, N, kOWR(1:M + 1), M, 'identity');
        Tid = combine(PL, Oid, N, 1);
        fprintf('  d0 from the identity tan k = -alpha/k instead of tan(k), Mo = %d: |dT|/|T| = %.2e\n', M, norm(Tid - Tr) / norm(Tr));
    end
end
kChk = 500; if QUICK, kChk = 200; end
tc = tic; kProd = dispersionRoots(alpha, 0, 0, kChk); tProd = toc(tc);
tc = tic; kFast = openWaterRoots(alpha, kChk); tFast = toc(tc);
fprintf('  Fast open-water roots vs dispersionRoots, Mo = %d: max rel diff %.2e  (%.2f s vs %.2f s)\n\n', ...
    kChk, max(abs(kFast - kProd) ./ abs(kProd)), tFast, tProd);

%% ---- Parts 1 and 2, per frequency -----------------------------------------------------------------------------------
res = struct('omega', {}, 'MoGrid', {}, 'Tmo', {}, 'Bmo', {}, 'errT', {}, 'envT', {}, 'qT', {}, 'tailT', {}, ...
    'Mladder', {}, 'TM', {}, 'pM', {}, 'Tinf', {}, 'unc', {}, 'errM', {});
showMo = [100 200 500 1000 2000 5000 10000 20000];
showMo = showMo(showMo < MoMax);
for w = 1:numel(omegaList)
    omega = omegaList(w); alpha = depth * omega^2 / g;
    fprintf('==================== omega = %.1f rad/s (alpha = %.4f) ====================\n', omega, alpha);
    tc = tic;
    kOW = openWaterRoots(alpha, MoMax);
    kap = dispersionRoots(alpha, bet0, gam0, max([Mladder Mfix]));
    O = openSide(alpha, gam0, R, P, N, kOW, MoGrid, 'identity');
    fprintf('Open-water roots and tower sums to Mo = %d: %.1f s\n', MoMax, toc(tc));

    % ---- Part 1: Mo sweep at fixed plate M
    fprintf('\nPART 1  Mo sweep: e(Mo) = |T(Mo) - T(MoMax)| / |T|, envelope = max of e over [Mo, Mo + 1/gamma]\n');
    Tmo = cell(1, numel(Mfix)); Bmo = cell(1, numel(Mfix)); errT = cell(1, numel(Mfix)); envT = cell(1, numel(Mfix));
    for i = 1:numel(Mfix)
        PL = plateSide(alpha, bet0, gam0, R, nu, P, N, kap(1:Mfix(i) + 3), sensR, sensTh);
        [Tmo{i}, Bmo{i}] = combine(PL, O, N, 1:numel(MoGrid));
        errT{i} = vecnorm(Tmo{i} - Tmo{i}(:, end), 2, 1) / norm(Tmo{i}(:, end));
        envT{i} = envelope(MoGrid, errT{i}, period);
    end
    fprintf('  %8s', 'Mo'); for i = 1:numel(Mfix), fprintf('   e (M=%3d)  env (M=%3d)', Mfix(i), Mfix(i)); end; fprintf('\n');
    for Mo = showMo
        [~, s] = min(abs(MoGrid - Mo));
        fprintf('  %8d', MoGrid(s)); for i = 1:numel(Mfix), fprintf('   %10.2e  %11.2e', errT{i}(s), envT{i}(s)); end; fprintf('\n');
    end
    % envelope decay fitted where the envelope is not biased by the MoMax reference
    sel = MoGrid >= 2000 & MoGrid <= MoMax / 4;
    if nnz(sel) >= 5
        c = polyfit(log(MoGrid(sel)), log(envT{end}(sel)), 1);
        qT = -c(1); tailT = exp(polyval(c, log(MoMax)));
        fprintf('  Envelope decay (M = %d, Mo %d..%d): env ~ Mo^-%.2f; extrapolated open-water error left at Mo = %d: %.2e\n', ...
            Mfix(end), 2000, MoMax / 4, qT, MoMax, tailT);
    else
        qT = NaN; tailT = NaN;
        fprintf('  (too few points for an envelope fit in quick mode)\n');
    end
    if numel(Mfix) > 1
        fprintf('  Plate-M independence: max over Mo of |e(M=%d) - e(M=%d)| = %.2e (small => truncations separable)\n', ...
            Mfix(1), Mfix(2), max(abs(errT{1} - errT{2})));
    end

    % ---- Part 2: plate M ladder at Mo = MoMax
    fprintf('\nPART 2  Plate M ladder at Mo = %d\n', MoMax);
    TM = complex(zeros(3, numel(Mladder)));
    for i = 1:numel(Mladder)
        PL = plateSide(alpha, bet0, gam0, R, nu, P, N, kap(1:Mladder(i) + 3), sensR, sensTh);
        TM(:, i) = combine(PL, O, N, numel(MoGrid));
    end
    sM = vecnorm(diff(TM, 1, 2), 2, 1) / norm(TM(:, end));
    pM = log2(sM(1:end - 1) ./ sM(2:end));
    fprintf('  %12s %12s %8s\n', 'M -> 2M', 's_T', 'p');
    for i = 1:numel(sM)
        if i <= numel(pM), fprintf('  %5d->%-6d %12.3e %8.2f\n', Mladder(i), Mladder(i + 1), sM(i), pM(i));
        else, fprintf('  %5d->%-6d %12.3e\n', Mladder(i), Mladder(i + 1), sM(i)); end
    end
    Tinf = NaN(3, 1); unc = NaN; errM = NaN(1, numel(Mladder));
    if numel(Mladder) >= 4
        % Richardson: (100, 200) with p from (50,100,200), and (200, 400) with p from (100,200,400)
        TinfA = TM(:, 3) + (TM(:, 3) - TM(:, 2)) / (2^pM(1) - 1);
        Tinf  = TM(:, 4) + (TM(:, 4) - TM(:, 3)) / (2^pM(2) - 1);
        unc = norm(Tinf - TinfA) / norm(Tinf);
        errM = vecnorm(TM - Tinf, 2, 1) / norm(Tinf);
        fprintf('  Richardson extrapolate T_inf from (200, 400), p = %.2f; uncertainty |T_inf(200,400) - T_inf(100,200)|/|T| = %.2e\n', pM(2), unc);
        fprintf('  Error against T_inf:'); fprintf('   M = %d: %.2e', [Mladder; errM]); fprintf('\n');
        for target = [1e-5 1e-6]
            if errM(end) <= target, fprintf('  M for |dT|/|T| <= %.0e: <= %d\n', target, Mladder(end));
            else, fprintf('  M for |dT|/|T| <= %.0e: about %.0f (power law from M = 400)\n', target, Mladder(end) * (errM(end) / target)^(1 / pM(2))); end
        end
    else
        fprintf('  (Richardson needs four levels; skipped in quick mode)\n');
    end
    fprintf('\n');
    res(w) = struct('omega', omega, 'MoGrid', MoGrid, 'Tmo', {Tmo}, 'Bmo', {Bmo}, 'errT', {errT}, 'envT', {envT}, ...
        'qT', qT, 'tailT', tailT, 'Mladder', Mladder, 'TM', TM, 'pM', pM, 'Tinf', Tinf, 'unc', unc, 'errM', errM);
end

%% ---- figure -------------------------------------------------------------------------------------------------------
cols = [0.85 0.33 0.10; 0.00 0.45 0.74];
fig = figure('Visible', 'off', 'Position', [100 100 1100 420]);
subplot(1, 2, 1); hold on;
for w = 1:numel(res)
    plot(res(w).MoGrid, res(w).errT{end}, '-', 'Color', [cols(w, :) 0.35]);
    plot(res(w).MoGrid, res(w).envT{end}, '-', 'Color', cols(w, :), 'LineWidth', 1.5);
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; box on;
xlabel('open-water modes Mo'); ylabel('|T(Mo) - T(Mo_{max})| / |T|');
title(sprintf('(a) open-water tower, plate M = %d (thin: error, thick: envelope)', Mfix(end)));
lab = cell(1, 2 * numel(res));
for w = 1:numel(res)
    lab{2 * w - 1} = sprintf('%.1f rad/s', res(w).omega);
    lab{2 * w} = sprintf('%.1f rad/s envelope', res(w).omega);
end
legend(lab, 'Location', 'southwest');
subplot(1, 2, 2); hold on;
for w = 1:numel(res)
    if all(isfinite(res(w).errM)), plot(res(w).Mladder, res(w).errM, '-o', 'Color', cols(w, :), 'LineWidth', 1.2); end
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; box on;
xlabel('plate modes M'); ylabel('|T(M) - T_\infty| / |T|');
title(sprintf('(b) plate tower at Mo = %d, against the Richardson extrapolate', MoMax));
legend(arrayfun(@(r) sprintf('%.1f rad/s', r.omega), res, 'UniformOutput', false), 'Location', 'southwest');
pngFile = fullfile(resDir, 'diag_F1_Mo_ladder.png');
exportgraphics(fig, pngFile, 'Resolution', 150);
close(fig);

meta = struct('beta', bet0, 'gamma', gam0, 'R', R, 'nu', nu, 'depth', depth, 'P', P, 'N', N, 'MoMax', MoMax, ...
    'Mfix', Mfix, 'sensorsRR', sensRR, 'quick', QUICK, 'when', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
save(fullfile(resDir, 'diag_F1_Mo_ladder.mat'), 'meta', 'res', '-v7');
fprintf('Saved %s\n      %s\n==== done ====\n', fullfile(resDir, 'diag_F1_Mo_ladder.mat'), pngFile);
diary('off');

%% ==== local functions ================================================================================================
function PL = plateSide(alpha, beta, gam, R, nu, P, N, kappa, sensR, sensTh)
% Everything that depends only on the plate tower, per angular order n:
%   B{n} = Dg.' * En_R * Mmat              (plate part of the (P+1)x(P+1) interface system)
%   H{n} = sensor weights x [Ln*Mmat; Mmat]  so that eta(sensors) = sum_n H{n} * U_n
    V = bennettsCoefficients(alpha, beta, gam, kappa);
    D  = innerProductMatrix(kappa(3:end), gam);
    Dg = gegenbauerCrossMatrix(kappa(3:end), gam, P, 1 - gam);
    nModes = numel(kappa);
    dphi = zeros(nModes, 1);
    for idx = 1:nModes, dphi(idx) = edgeDerivative(kappa(idx), gam); end
    PL.B = cell(1, 2 * N + 1); PL.H = cell(1, 2 * N + 1);
    for nIdx = 1:(2 * N + 1)
        n = nIdx - N - 1;
        Ln = edgeMap(n, alpha, beta, gam, R, nu, kappa);
        Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);
        Mmat = (D * Enprime) \ Dg;
        En_R = radialMatrix(n, R, R, kappa, V, Ln);
        PL.B{nIdx} = Dg.' * En_R * Mmat;
        Q = [Ln * Mmat; Mmat];                                         % Afull = Q * U
        Ih = ihat(n, kappa(:), sensR(:).', R);                         % nModes x nSensors
        Wt = (dphi .* Ih).' .* exp(1i * n * sensTh(:)) / (1i * alpha);  % nSensors x nModes
        PL.H{nIdx} = Wt * Q;
    end
end

function O = openSide(alpha, gam, R, P, N, k, MoSel, d0mode)
% Everything that depends only on the open-water tower k = [k0; k1; ...; kMo], with the cumulative
% interface sum G0 extracted at the requested truncations MoSel (numbers of evanescent modes).
    k = k(:); p1 = P + 1;
    D0g = gegenbauerVec(k, gam, P, 1);
    switch d0mode
        case 'tan'                                     % as innerProductMatrix(k, 0), diagonal
            d0 = 0.5 * (tan(k) ./ k + 1 ./ cos(k).^2);
        case 'identity'                                % at a root, k tan k = -alpha exactly
            t = -alpha ./ k;
            d0 = 0.5 * (t ./ k + 1 + t.^2);
    end
    x = k * R; x0 = x(1); k0 = k(1);
    rows = MoSel(:) + 1; nS = numel(rows);
    [jj, ll] = ndgrid(1:p1, 1:p1);
    O.G0 = cell(1, 2 * N + 1); O.rhs = cell(1, 2 * N + 1); O.m0row = cell(1, 2 * N + 1);
    O.mIn1 = zeros(1, 2 * N + 1); O.a = zeros(1, 2 * N + 1); O.baf = zeros(1, 2 * N + 1);
    p0 = imag(k0);
    for nIdx = 1:(2 * N + 1)
        n = nIdx - N - 1;
        E0p = k .* (-besselk(n - 1, x, 1) ./ besselk(n, x, 1) - n ./ x);   % d/dr Khat_n at r = R
        wgt = 1 ./ (d0 .* E0p);
        C = cumsum(D0g(:, jj(:)) .* D0g(:, ll(:)) .* wgt, 1);              % (Mo+1) x (P+1)^2
        O.G0{nIdx} = reshape(C(rows, :).', p1, p1, nS);
        EIp1 = k0 * (besseli(n - 1, x0, 1) / besseli(n, x0, 1) - n / x0);  % d/dr Ihat_n at r = R, m = 0
        mIn1 = -EIp1 / E0p(1);
        a = besseli(n, k0 * R) / (1i * alpha);                              % incidentAmplitude, m = 0
        O.rhs{nIdx} = D0g(1, :).' * (1 + mIn1) * a;
        O.m0row{nIdx} = D0g(1, :) * wgt(1);
        O.mIn1(nIdx) = mIn1; O.a(nIdx) = a;
        O.baf(nIdx) = besselj(n, p0 * R) / (besselh(n, 2, p0 * R) * a);
    end
end

function [T, ba] = combine(PL, O, N, sIdx)
% Solve the (P+1)x(P+1) interface system for each angular order and each selected truncation.
    nS = numel(sIdx);
    T = complex(zeros(3, nS)); ba = complex(zeros(2 * N + 1, nS));
    for q = 1:nS
        s = sIdx(q);
        e = zeros(size(PL.H{1}, 1), 1);
        for nIdx = 1:(2 * N + 1)
            U = (PL.B{nIdx} - O.G0{nIdx}(:, :, s)) \ O.rhs{nIdx};
            e = e + PL.H{nIdx} * U;
            ba(nIdx, q) = (O.mIn1(nIdx) * O.a(nIdx) + O.m0row{nIdx} * U) * O.baf(nIdx);
        end
        T(:, q) = e(1:3) / e(4);
    end
end

function [T, ba] = productionSolve(alpha, beta, gam, R, nu, P, N, kappa, k, sensR, sensTh)
% The production pipeline, unchanged (as precomputeDeflectionData.m), for the Part 0 regression.
    M = numel(k) - 1;
    V = bennettsCoefficients(alpha, beta, gam, kappa);
    k0 = k(1); p0 = imag(k0);
    D0  = innerProductMatrix(k, 0);
    D   = innerProductMatrix(kappa(3:end), gam);
    D0g = gegenbauerCrossMatrix(k, gam, P, 1);
    Dg  = gegenbauerCrossMatrix(kappa(3:end), gam, P, 1 - gam);
    nModes = numel(kappa);
    dphiFull = zeros(nModes, 1);
    for idx = 1:nModes, dphiFull(idx) = edgeDerivative(kappa(idx), gam); end
    AfullMat = complex(zeros(nModes, 2 * N + 1)); ba = complex(zeros(2 * N + 1, 1));
    for nIdx = 1:(2 * N + 1)
        n = nIdx - N - 1;
        Ln = edgeMap(n, alpha, beta, gam, R, nu, kappa);
        [EIprime, E0prime] = openWaterRadialDerivatives(n, R, k);
        Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);
        [MIn, M0, Mmat] = matchingMatrices(EIprime, E0prime, Enprime, D0, D, D0g, Dg);
        En_R = radialMatrix(n, R, R, kappa, V, Ln);
        [~, Sn, Sn0] = scatteringMatrices(D0g, Dg, En_R, MIn, M0, Mmat);
        AIn = incidentAmplitude(alpha, n, k0, R, M);
        An = Sn * AIn; Aen = Ln * An;
        AfullMat(:, nIdx) = [Aen(1); Aen(2); An];
        A0 = Sn0 * AIn;
        ba(nIdx) = (A0(1) / besselh(n, 2, p0 * R)) / (AIn(1) / besselj(n, p0 * R));
    end
    data = struct('alpha', alpha, 'R', R, 'N', N, 'kappa', kappa, 'dphiFull', dphiFull, 'AfullMat', AfullMat);
    e = evaluateDeflection(data, sensR, sensTh); e = e(:);
    T = e(1:3) / e(4);
end

function Dg = gegenbauerVec(r, gam, P, normScale)
% Vectorised gegenbauerCrossMatrix (same closed form, delta = 1/6).
    delta = 1 / 6; G = builtin('gamma', delta);
    r = r(:); x = r * (1 - gam);
    Dg = zeros(numel(r), P + 1);
    for j = 0:P
        Dg(:, j + 1) = G * (2 ./ x).^delta * (2 * j + delta) * (-1)^j .* besselj(2 * j + delta, x) ./ cos(r * normScale);
    end
end

function k = openWaterRoots(alpha, Mo)
% [k0; k1; ...; kMo] as dispersionRoots(alpha, 0, 0, Mo). k_m lies in ((m - 1/2) pi, m pi), where
% f(k) = k sin k + alpha cos k changes sign. fzero for m <= 50; vectorised Newton from
% k ~ m pi - atan(alpha/(m pi)) beyond, checked against the bracket (fzero for any that fail).
    k0 = dispersionRoots(alpha, 0, 0, 1); k0 = k0(1);
    f  = @(x) x .* sin(x) + alpha * cos(x);
    fp = @(x) (1 - alpha) * sin(x) + x .* cos(x);
    m = (1:Mo).';
    km = zeros(Mo, 1);
    mSmall = min(50, Mo);
    for i = 1:mSmall, km(i) = fzero(f, [(i - 0.5) * pi, i * pi]); end
    if Mo > mSmall
        mm = m(mSmall + 1:end);
        x = mm * pi - atan(alpha ./ (mm * pi));
        for it = 1:8, x = x - f(x) ./ fp(x); end
        bad = ~(x > (mm - 0.5) * pi & x < mm * pi) | ~isfinite(x);
        for i = find(bad).', x(i) = fzero(f, [(mm(i) - 0.5) * pi, mm(i) * pi]); end
        km(mSmall + 1:end) = x;
    end
    k = [k0; km];
end

function I = ihat(n, kap, r, R)
% I_n(kappa r) / I_n(kappa R), overflow-safe (exponentially scaled besseli). kap column, r row.
    I = besseli(n, kap .* r, 1) ./ besseli(n, kap * R, 1) .* exp(abs(real(kap .* r)) - abs(real(kap * R)));
end

function env = envelope(Mo, e, period)
% Running maximum of e over the window [Mo, Mo + period].
    env = zeros(size(e));
    for i = 1:numel(Mo), env(i) = max(e(Mo >= Mo(i) & Mo <= Mo(i) + period)); end
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
