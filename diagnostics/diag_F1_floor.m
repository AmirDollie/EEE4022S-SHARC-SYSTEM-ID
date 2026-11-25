%% diag_F1_floor.m  --  F1 follow-up: what sets the ~1e-4 floor in the M ladder?
%
% CONTEXT
%   diag_F1_convergence.m (MATLAB, 2026-11-25) showed T, U and b/a converging at p ~ 1.7 up to M = 100,
%   then the doubling steps stop shrinking (~1e-4 relative from M = 200 to 400), at every radius and with
%   or without the Bennetts reduction. Because b/a and U plateau too, the floor is in the matching
%   solution, not in the deflection reconstruction. This script tests three explanations.
%
%   H1  Two-tower beat. Open-water roots k_m ~ m*pi and plate roots kappa_m ~ m*pi/(1-gamma). The
%       Gegenbauer projections of the OPEN-WATER modes are besselj(., k_m (1-gamma)), whose phase drifts
%       by ~pi*gamma*m against the plate tower, so the open-water tail is modulated with a period of
%       about 1/gamma ~ 690 modes. That would make convergence look stalled for M of order 1/gamma.
%   H2  Precision floor: an ill-conditioned step amplifies rounding error to ~1e-4.
%   H3  Genuine slow convergence that oscillates, which a coarse doubling ladder cannot resolve.
%
% WHAT IT RUNS  (reference floe p0, omega = 8.5 rad/s, the worst frequency in the ladder; P = 10, N = 10;
%               reference layout s2, s6, s14 against s26)
%   Part 1  Dense M sweep (M = 100:10:500, open-water and plate towers equal) at gamma0/2, gamma0, 2*gamma0
%           with beta and R held fixed. Reports the 10-mode step s(M) of T and b/a, its local slope, and a
%           jaggedness ratio J = median|2nd difference| / median|1st difference|.
%   Part 2  Tower split at gamma0. The towers couple only through the P+1 Gegenbauer coefficients, so
%           they can be truncated separately: (a) open-water count Mo swept with the plate tower fixed at
%           Mmax; (b) plate count M swept with the open-water tower fixed at Mmax.
%   Part 3  Precision probe at M = Mo = 400: all roots perturbed by a relative 1e-13 (three draws), plus
%           the smallest rcond of the three linear solves per angular order.
%
% HOW TO READ IT
%   H1 (beat):        J small (smooth); the slope flattens EARLIER at 2*gamma0 and LATER at gamma0/2; the
%                     compensated curves in figure panel (b) line up when plotted against gamma*M; the
%                     floor sits in the Mo sweep (2a), not the M sweep (2b).
%   H2 (precision):   J ~ 1 or more (jagged); Part 3 change ~1e-5 or more (amplification > 1e8); eps/rcond
%                     close to the floor; insensitive to gamma.
%   H3 (oscillating): J small, but no consistent shift with gamma.
%
% RUNTIME  About 5 minutes (quick mode under 1 minute).
% WRITES   results/diag_F1_floor.log, results/diag_F1_floor.mat, results/diag_F1_floor.png
% HOW      In MATLAB from the repo root:   run('diagnostics/diag_F1_floor.m')
%          Quick smoke test:               FLOOR_QUICK = true; run('diagnostics/diag_F1_floor.m')
%          (clear FLOOR_QUICK before the full run: the desktop workspace keeps it.)
% Nothing in Forward Model is modified.

%% ---- setup --------------------------------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
root = fileparts(here);
addpath(fullfile(root, 'Forward Model'));
addpath(fullfile(root, 'Forward Model', 'Animation'));
resDir = fullfile(root, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
logFile = fullfile(resDir, 'diag_F1_floor.log');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile); diary('on');
cleanupDiary = onCleanup(@() diary('off'));

QUICK = exist('FLOOR_QUICK', 'var') && FLOOR_QUICK;
fprintf('==== diag_F1_floor  (%s)%s ====\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ternary(QUICK, '  [QUICK MODE]', ''));

% reference floe p0 and the worst frequency of the ladder
bet0 = 4.6985e-5; gam0 = 1.4548e-3; R = 0.3830; nu = 0.3;
depth = 1.88; g = 9.81;
omega = 8.5; alpha = depth * omega^2 / g;
sensRR = [0.3 0; 0.3 pi; 0.5 pi; 0.9 0];          % s2, s6, s14, s26  [r/R, theta]
sensR = sensRR(:, 1).' * R; sensTh = sensRR(:, 2).';
refIdx = 4;
N = 10; P = 10;

Mgrid = 100:10:500;
perturbM = 400;
if QUICK, Mgrid = 100:50:300; perturbM = 200; end
Mmax = max(Mgrid);
gamList = gam0 * [0.5 1 2];
gamNames = {'gamma0/2', 'gamma0', '2*gamma0'};
pertSize = 1e-13; nPert = 3;
pRef = 1.7;                                        % clean-regime rate from diag_F1_convergence (M <= 100)

kRmax = (Mmax + 0.5) * pi / (1 - max(gamList)) * R;
fprintf('omega = %.1f rad/s, alpha = %.4f, P = %d, N = %d, M grid %d:%d:%d\n', omega, alpha, P, N, ...
    Mgrid(1), Mgrid(2) - Mgrid(1), Mgrid(end));
fprintf('Largest tower argument kappa_M R = %.1f (overflow near 700): %s\n\n', kRmax, ternary(kRmax < 650, 'OK', 'TOO LARGE'));
if kRmax >= 650, error('diag_F1_floor:overflow', 'Mmax too large for unscaled Bessel functions at this R.'); end

%% ---- roots -----------------------------------------------------------------------------------------------------
tc = tic;
kOW = dispersionRoots(alpha, 0, 0, Mmax);                       % open water: independent of gamma
kapG = cell(1, numel(gamList));
for gi = 1:numel(gamList)
    kapG{gi} = dispersionRoots(alpha, bet0, gamList(gi), Mmax);
end
fprintf('Root searches at M = %d for open water and three draughts: %.1f s\n\n', Mmax, toc(tc));

nM = numel(Mgrid);
Mmid = Mgrid(1:end - 1);                                        % s(M) = |f(M+dM) - f(M)| / |f|

%% ---- Part 1: dense M sweep at three draughts ---------------------------------------------------------------------
fprintf('PART 1  Dense M sweep (Mo = M), 10-mode steps s(M) = |f(M+10) - f(M)| / |f(Mmax)|\n');
P1 = struct('gam', {}, 'T', {}, 'ba', {}, 'eta', {}, 'rc', {}, 'sT', {}, 'sB', {}, 'sec', {});
for gi = 1:numel(gamList)
    gam = gamList(gi);
    T = complex(zeros(3, nM)); ba = complex(zeros(2 * N + 1, nM)); eta = complex(zeros(4, nM)); rc = zeros(1, nM);
    tc = tic;
    for i = 1:nM
        M = Mgrid(i);
        out = solveOrders(alpha, bet0, gam, R, nu, P, N, kapG{gi}(1:M + 3), kOW(1:M + 1));
        e = evaluateDeflection(out.data, sensR, sensTh); e = e(:);
        eta(:, i) = e; T(:, i) = e(1:3) / e(refIdx); ba(:, i) = out.ba; rc(i) = out.rc;
    end
    P1(gi) = struct('gam', gam, 'T', T, 'ba', ba, 'eta', eta, 'rc', rc, ...
        'sT', stepSeries(T), 'sB', stepSeries(ba), 'sec', toc(tc));
    fprintf('  %-9s (gamma = %.4e, 1/gamma = %4.0f modes)  %.1f s\n', gamNames{gi}, gam, 1 / gam, P1(gi).sec);
end

showM = Mmid(ismember(Mmid, [100 150 200 250 300 350 400 450 490]));
fprintf('\n  s_T(M):\n  %8s', 'M'); fprintf('%12s', gamNames{:}); fprintf('\n');
for M = showM
    i = find(Mmid == M);
    fprintf('  %8d', M); fprintf('%12.3e', arrayfun(@(r) r.sT(i), P1)); fprintf('\n');
end
fprintf('\n  s_b/a(M):\n  %8s', 'M'); fprintf('%12s', gamNames{:}); fprintf('\n');
for M = showM
    i = find(Mmid == M);
    fprintf('  %8d', M); fprintf('%12.3e', arrayfun(@(r) r.sB(i), P1)); fprintf('\n');
end

fprintf('\n  Local slope of log s_T against log M (a clean power law E ~ M^-p gives about -(p+1) = %.1f),\n', -(pRef + 1));
fprintf('  jaggedness J (smooth << 1, noise ~ 1 or more) and eps/rcond (rough precision floor):\n');
fprintf('  %-9s %14s %14s %8s %8s %14s\n', '', 'slope 100-200', 'slope 300-500', 'J_T', 'J_b/a', 'eps/min rcond');
for gi = 1:numel(gamList)
    fprintf('  %-9s %14.2f %14.2f %8.2f %8.2f %14.2e\n', gamNames{gi}, localSlope(Mmid, P1(gi).sT, 100, 200), ...
        localSlope(Mmid, P1(gi).sT, 300, 500), jaggedness(P1(gi).T), jaggedness(P1(gi).ba), eps / min(P1(gi).rc));
end
fprintf('  Reading: H1 predicts the slope flattens earliest at 2*gamma0 and stays steep longest at gamma0/2.\n\n');

%% ---- Part 2: which tower carries the floor ---------------------------------------------------------------------------
fprintf('PART 2  Tower split at gamma0 (towers couple only through the P+1 = %d Gegenbauer coefficients)\n', P + 1);
gi0 = 2;
Ta = complex(zeros(3, nM)); Ba = complex(zeros(2 * N + 1, nM));
Tb = complex(zeros(3, nM)); Bb = complex(zeros(2 * N + 1, nM));
tc = tic;
for i = 1:nM
    M = Mgrid(i);
    outA = solveOrders(alpha, bet0, gam0, R, nu, P, N, kapG{gi0}(1:Mmax + 3), kOW(1:M + 1));   % sweep Mo
    e = evaluateDeflection(outA.data, sensR, sensTh); e = e(:);
    Ta(:, i) = e(1:3) / e(refIdx); Ba(:, i) = outA.ba;
    outB = solveOrders(alpha, bet0, gam0, R, nu, P, N, kapG{gi0}(1:M + 3), kOW(1:Mmax + 1));   % sweep M
    e = evaluateDeflection(outB.data, sensR, sensTh); e = e(:);
    Tb(:, i) = e(1:3) / e(refIdx); Bb(:, i) = outB.ba;
end
P2 = struct('sTa', stepSeries(Ta), 'sBa', stepSeries(Ba), 'sTb', stepSeries(Tb), 'sBb', stepSeries(Bb), ...
    'Ta', Ta, 'Ba', Ba, 'Tb', Tb, 'Bb', Bb, 'sec', toc(tc));
fprintf('  (%.1f s)\n  %8s %16s %16s %16s %16s\n', P2.sec, 'M', 's_T (a: Mo)', 's_T (b: M)', 's_b/a (a: Mo)', 's_b/a (b: M)');
for M = showM
    i = find(Mmid == M);
    fprintf('  %8d %16.3e %16.3e %16.3e %16.3e\n', M, P2.sTa(i), P2.sTb(i), P2.sBa(i), P2.sBb(i));
end
fprintf('  %8s %16.2f %16.2f %16.2f %16.2f   (slope 300-500)\n', '', localSlope(Mmid, P2.sTa, 300, 500), ...
    localSlope(Mmid, P2.sTb, 300, 500), localSlope(Mmid, P2.sBa, 300, 500), localSlope(Mmid, P2.sBb, 300, 500));
fprintf('  %8s %16.2f %16.2f %16.2f %16.2f   (jaggedness J)\n', '', jaggedness(Ta), jaggedness(Tb), jaggedness(Ba), jaggedness(Bb));
fprintf('  Reading: H1 predicts the floor in (a), the open-water tower, and clean decay in (b).\n\n');

%% ---- Part 3: precision probe ------------------------------------------------------------------------------------------
fprintf('PART 3  Precision probe at M = Mo = %d, gamma0: all roots scaled by (1 + %.0e * randn), %d draws\n', perturbM, pertSize, nPert);
kap0 = kapG{gi0}(1:perturbM + 3); k0v = kOW(1:perturbM + 1);
base = solveOrders(alpha, bet0, gam0, R, nu, P, N, kap0, k0v);
e = evaluateDeflection(base.data, sensR, sensTh); e = e(:); Tbase = e(1:3) / e(refIdx);
dT = zeros(1, nPert); dB = zeros(1, nPert);
for d = 1:nPert
    rng(d);
    kapP = kap0 .* (1 + pertSize * randn(size(kap0)));
    kP = k0v .* (1 + pertSize * randn(size(k0v)));
    pert = solveOrders(alpha, bet0, gam0, R, nu, P, N, kapP, kP);
    e = evaluateDeflection(pert.data, sensR, sensTh); e = e(:);
    dT(d) = norm(e(1:3) / e(refIdx) - Tbase) / norm(Tbase);
    dB(d) = norm(pert.ba - base.ba) / norm(base.ba);
end
fprintf('  relative change in T:   '); fprintf('%10.2e', dT); fprintf('   (amplification up to %.1e)\n', max(dT) / pertSize);
fprintf('  relative change in b/a: '); fprintf('%10.2e', dB); fprintf('   (amplification up to %.1e)\n', max(dB) / pertSize);
fprintf('  smallest rcond over the three solves and all orders: %.2e  (eps/rcond = %.2e)\n', base.rc, eps / base.rc);
fprintf('  Reading: a change near 1e-5 or more (amplification > 1e8) means a precision floor (H2).\n\n');

%% ---- figure -----------------------------------------------------------------------------------------------------------
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
fig = figure('Visible', 'off', 'Position', [100 100 1350 400]);
subplot(1, 3, 1); hold on;
for gi = 1:numel(gamList), plot(Mmid, P1(gi).sT, '-o', 'Color', cols(gi, :), 'MarkerSize', 3); end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; box on;
xlabel('M'); ylabel('s_T(M) = |T(M+10) - T(M)| / |T|'); title('(a) 10-mode step in T');
legend(gamNames, 'Location', 'southwest');
subplot(1, 3, 2); hold on;
for gi = 1:numel(gamList)
    plot(gamList(gi) * Mmid, P1(gi).sT .* (Mmid / 100).^(pRef + 1), '-o', 'Color', cols(gi, :), 'MarkerSize', 3);
end
set(gca, 'YScale', 'log'); grid on; box on;
xlabel('\gamma M'); ylabel(sprintf('s_T (M/100)^{%.1f}', pRef + 1));
title('(b) compensated step vs \gammaM (flat = clean power law)');
subplot(1, 3, 3); hold on;
plot(Mmid, P1(gi0).sT, '-o', 'Color', [0.4 0.4 0.4], 'MarkerSize', 3);
plot(Mmid, P2.sTa, '-s', 'Color', cols(1, :), 'MarkerSize', 3);
plot(Mmid, P2.sTb, '-^', 'Color', cols(2, :), 'MarkerSize', 3);
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; box on;
xlabel('M (swept tower)'); ylabel('s_T(M)'); title('(c) which tower: \gamma_0');
legend({'both towers', 'open water (Mo) swept', 'plate (M) swept'}, 'Location', 'southwest');
pngFile = fullfile(resDir, 'diag_F1_floor.png');
exportgraphics(fig, pngFile, 'Resolution', 150);
close(fig);

%% ---- save -------------------------------------------------------------------------------------------------------------
meta = struct('beta', bet0, 'gamma0', gam0, 'gamList', gamList, 'R', R, 'nu', nu, 'depth', depth, 'omega', omega, ...
    'alpha', alpha, 'P', P, 'N', N, 'Mgrid', Mgrid, 'sensorsRR', sensRR, 'ref', refIdx, 'perturbM', perturbM, ...
    'pertSize', pertSize, 'quick', QUICK, 'when', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
P3 = struct('dT', dT, 'dB', dB, 'rc', base.rc);
save(fullfile(resDir, 'diag_F1_floor.mat'), 'meta', 'P1', 'P2', 'P3', '-v7');
fprintf('Saved %s\n      %s\n', fullfile(resDir, 'diag_F1_floor.mat'), pngFile);
fprintf('==== done ====\n');
diary('off');

%% ==== local functions =================================================================================================
function out = solveOrders(alpha, beta, gam, R, nu, P, N, kappa, k)
% Production matching pipeline (as precomputeDeflectionData.m and diag_F1_convergence.m), with the
% open-water tower k (Mo+1 roots) and the plate tower kappa (M+3 roots) allowed to differ in length.
% Also returns b/a of the propagating open-water mode and the smallest rcond of the three solves.
    Mo = numel(k) - 1;
    V = bennettsCoefficients(alpha, beta, gam, kappa);
    k0 = k(1);
    D0  = innerProductMatrix(k, 0);
    D   = innerProductMatrix(kappa(3:end), gam);
    D0g = gegenbauerCrossMatrix(k, gam, P, 1);
    Dg  = gegenbauerCrossMatrix(kappa(3:end), gam, P, 1 - gam);
    nModes = numel(kappa);
    dphiFull = zeros(nModes, 1);
    for idx = 1:nModes, dphiFull(idx) = edgeDerivative(kappa(idx), gam); end
    AfullMat = complex(zeros(nModes, 2 * N + 1));
    ba = complex(zeros(2 * N + 1, 1));
    rc = Inf;
    p0 = imag(k0);
    for nIdx = 1:(2 * N + 1)
        n = nIdx - N - 1;
        Ln = edgeMap(n, alpha, beta, gam, R, nu, kappa);
        [EIprime, E0prime] = openWaterRadialDerivatives(n, R, k);
        Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);
        [MIn, M0, Mmat] = matchingMatrices(EIprime, E0prime, Enprime, D0, D, D0g, Dg);
        En_R = radialMatrix(n, R, R, kappa, V, Ln);
        [~, Sn, Sn0] = scatteringMatrices(D0g, Dg, En_R, MIn, M0, Mmat);
        AIn = incidentAmplitude(alpha, n, k0, R, Mo);
        An = Sn * AIn;
        Aen = Ln * An;
        AfullMat(:, nIdx) = [Aen(1); Aen(2); An];
        A0 = Sn0 * AIn;
        ba(nIdx) = (A0(1) / besselh(n, 2, p0 * R)) / (AIn(1) / besselj(n, p0 * R));
        rc = min([rc, rcond(D * Enprime), rcond(D0 * E0prime), rcond(Dg.' * En_R * Mmat - D0g.' * M0)]);
    end
    out.data = struct('alpha', alpha, 'R', R, 'N', N, 'kappa', kappa, 'dphiFull', dphiFull, 'AfullMat', AfullMat);
    out.ba = ba;
    out.rc = rc;
end

function s = stepSeries(F)
% Relative step between neighbouring grid points, normalised by the value at the finest point.
    s = zeros(1, size(F, 2) - 1);
    for i = 1:numel(s), s(i) = norm(F(:, i + 1) - F(:, i)); end
    s = s / norm(F(:, end));
end

function J = jaggedness(F)
% median |second difference| / median |first difference| along the M grid.
% A smooth power-law decay on a 10-mode grid gives J well below 1; rounding noise gives J ~ 1 to 2.
    d1 = vecnorm(diff(F, 1, 2), 2, 1);
    d2 = vecnorm(diff(F, 2, 2), 2, 1);
    J = median(d2) / median(d1);
end

function sl = localSlope(Mm, s, Mlo, Mhi)
% Least-squares slope of log s against log M over Mlo <= M <= Mhi (NaN if fewer than 3 points).
    sel = Mm >= Mlo & Mm <= Mhi & s > 0;
    if nnz(sel) < 3, sl = NaN; return; end
    c = polyfit(log(Mm(sel)), log(s(sel)), 1);
    sl = c(1);
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
