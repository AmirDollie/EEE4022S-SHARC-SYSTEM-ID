%% diag_F1_convergence.m  --  F1 forward-model convergence diagnostic (paper plan, Week 1, Wed)
%
% WHAT IT RUNS
%   Read-only use of the production Forward Model (dispersionRoots, bennettsCoefficients, edgeMap,
%   matchingMatrices, scatteringMatrices, evaluateDeflection, ...). Nothing in the repo is modified.
%   Reference floe p0 (beta = 4.6985e-5, gamma = 1.4548e-3, R = 0.3830, nu = 0.3), depth 1.88 m,
%   reference layout s2, s6, s14, s26 (reference sensor s26), incidence along theta = 0,
%   six frequencies across the validated band (3.0 to 8.5 rad/s).
%
%   Part A  Energy conservation. For an inviscid, lossless floe each angular order n must satisfy
%           |1 + 2 b_n / a_n| = 1, where a_n J_n(p r) + b_n H2_n(p r) is the propagating open-water
%           mode (a_n incident, b_n scattered). Also checks the n <-> -n symmetry b_n/a_n = b_-n/a_-n.
%           A residual that does not fall with M is evidence of a defect, not of slow convergence.
%   Part B  Truncation ladders.
%           M ladder at P = 10, N = 10:  M = 25, 50, 100, 200, 400  (doubling, so the rate can be fitted)
%           P ladder at M = 400, N = 10: P = 5, 10, 15, 20, 30
%           Recorded per configuration and frequency: sensor deflections, the three transmissibilities
%           T_j = eta_j / eta_s26, the Gegenbauer interface coefficients U, the first tower amplitudes,
%           and the far-field ratios b_n / a_n. Reported as RELATIVE norms (not Fisher units), with the
%           algebraic rate p fitted from the doubling ladder: |f(2M) - f(M)| ~ M^-p.
%   Part C  Bennetts reduction. Relative L2 error of phi_{-j}(z) ~ sum_m V(m,j) phi_m(z) on [-1, -gamma]
%           against M, and the mismatch of its z-derivative at z = -gamma. If Part B's rate matches
%           Part C's, the reduction is a candidate bottleneck.
%   Part D  Same M ladder with the damped modes projected EXACTLY (closed-form <phi_i, phi_-j> and
%           Gegenbauer projections of phi_-j) instead of through the truncated Bennetts series. If the
%           rate improves, the reduction is the bottleneck; if not, the slow convergence is intrinsic to
%           the truncated tower projection and extrapolation (or larger M) is the remedy.
%   Part E  Radial probe: deflection at r/R = 0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95 (theta = 0) on the M ladder,
%           per point. Tests whether convergence depends on sensor radius, and whether the change from
%           M = 200 to 400 keeps falling (Octave showed a ~1e-4 plateau there that conditioning, tan
%           precision and P do not explain).
%
% RUNTIME  Roughly 10 to 25 minutes on a laptop (dominated by the M = 400 root searches, which are done
%          once per frequency and re-used by every configuration).
% WRITES   results/diag_F1_convergence.mat  and  results/diag_F1_convergence.log  (diary of everything
%          printed), in the repo root. Delete or rename them to rerun cleanly.
% HOW      From the repo root:   run diagnostics/diag_F1_convergence.m
%          or from any shell:     matlab -batch "run('diagnostics/diag_F1_convergence.m')"

%% ---- setup ------------------------------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
root = fileparts(here);
addpath(fullfile(root, 'Forward Model'));
addpath(fullfile(root, 'Forward Model', 'Animation'));
resDir = fullfile(root, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
logFile = fullfile(resDir, 'diag_F1_convergence.log');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile); diary('on');
cleanupDiary = onCleanup(@() diary('off'));

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
QUICK = exist('F1_QUICK', 'var') && F1_QUICK;   % set F1_QUICK = true before running for a smoke test

fprintf('==== diag_F1_convergence  (%s)%s ====\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
    ternary(QUICK, '  [QUICK MODE]', ''));

% reference floe p0 and physical scaling
bet = 4.6985e-5; gam = 1.4548e-3; R = 0.3830; nu = 0.3;
depth = 1.88; g = 9.81;
omega = [3.0 4.0 5.0 6.0 7.0 8.5];
if QUICK, omega = [4.0 7.0]; end
alphaW = depth * omega.^2 / g;
sensorNames = {'s2', 's6', 's14', 's26'};
sensRR = [0.3 0; 0.3 pi; 0.5 pi; 0.9 0];          % [r/R, theta]
sensR = sensRR(:, 1).' * R; sensTh = sensRR(:, 2).';
refIdx = 4;                                        % s26
probeRR = [0 0.3 0.5 0.7 0.8 0.9 0.95]; probeR = probeRR * R; probeTh = zeros(size(probeR));
N = 10;

Mladder = [25 50 100 200 400];  Pfix = 10;
Pladder = [5 10 15 20 30];      Mfix = 400;
if QUICK, Mladder = [25 50 100]; Mfix = 100; Pladder = [5 10 20]; end
Mmax = max([Mladder Mfix]);

% Bessel overflow guard: besseli/besselk lose range near |x| ~ 700
kRmax = (Mmax + 0.5) * pi / (1 - gam) * R;
fprintf('Largest tower argument kappa_M R = %.1f (overflow near 700): %s\n\n', kRmax, ...
    ternary(kRmax < 650, 'OK', 'TOO LARGE - lower Mmax'));
if kRmax >= 650, error('diag_F1:overflow', 'Mmax too large for unscaled Bessel functions at this R.'); end

%% ---- roots once per frequency at Mmax -------------------------------------------------------------------------
nW = numel(omega);
rts = struct('kappa', {}, 'k', {});
fprintf('Root searches at M = %d\n', Mmax);
for w = 1:nW
    tc = tic;
    rts(w).kappa = dispersionRoots(alphaW(w), bet, gam, Mmax);   % [kappa_-2; kappa_-1; kappa_0; ...; kappa_M]
    rts(w).k     = dispersionRoots(alphaW(w), 0, 0, Mmax);       % [k_0; ...; k_M]
    fprintf('  omega = %4.1f rad/s  alpha = %7.4f   %.1f s\n', omega(w), alphaW(w), toc(tc));
end
fprintf('\n');

%% ---- configurations ------------------------------------------------------------------------------------------
cfgs = [ [Mladder(:), repmat(Pfix, numel(Mladder), 1)]; [repmat(Mfix, numel(Pladder), 1), Pladder(:)] ];
cfgs = unique(cfgs, 'rows', 'stable');
nC = size(cfgs, 1);
res = struct('M', {}, 'P', {}, 'eta', {}, 'T', {}, 'U', {}, 'Atower', {}, 'ba', {}, 'eRes', {}, 'symRes', {}, 'sec', {}, 'probe', {});

for c = 1:nC
    M = cfgs(c, 1); P = cfgs(c, 2);
    eta = complex(zeros(numel(sensR), nW));
    U = complex(zeros(P + 1, 2 * N + 1, nW));
    Atw = complex(zeros(min(M + 1, 11), 2 * N + 1, nW));
    ba = complex(zeros(2 * N + 1, nW));
    probe = complex(zeros(numel(probeR), nW));
    tc = tic;
    for w = 1:nW
        kap = rts(w).kappa(1:M + 3);
        kk  = rts(w).k(1:M + 1);
        out = solveAllOrders(alphaW(w), bet, gam, R, nu, M, P, N, kap, kk);
        e = evaluateDeflection(out.data, sensR, sensTh);
        eta(:, w) = e(:);
        ep = evaluateDeflection(out.data, probeR, probeTh); probe(:, w) = ep(:);
        U(:, :, w) = out.U;
        Atw(:, :, w) = out.Atower(1:size(Atw, 1), :);
        ba(:, w) = out.ba(:);
    end
    T = eta(1:3, :) ./ repmat(eta(refIdx, :), 3, 1);
    eRes = abs(abs(1 + 2 * ba) - 1);                    % (2N+1) x nW
    symRes = abs(ba(1:N, :) - flipud(ba(N + 2:end, :))) ./ max(abs(ba(1:N, :)), eps);
    res(c) = struct('M', M, 'P', P, 'eta', eta, 'T', T, 'U', U, 'Atower', Atw, 'ba', ba, ...
        'eRes', eRes, 'symRes', symRes, 'sec', toc(tc), 'probe', probe);
    fprintf('config [M P N] = [%3d %2d %2d]  %6.1f s   max energy residual %.2e   max n/-n asymmetry %.2e\n', ...
        M, P, N, res(c).sec, max(eRes(:)), max(symRes(:)));
end
fprintf('\n');

%% ---- Part A: energy ---------------------------------------------------------------------------------------------
fprintf('PART A  Energy conservation |1 + 2 b_n/a_n| - 1 (max over n), per frequency\n');
fprintf('  %-12s', '[M P]'); fprintf('%11.1f', omega); fprintf('   (rad/s)\n');
for c = 1:nC
    fprintf('  [%3d %2d]    ', res(c).M, res(c).P); fprintf('%11.2e', max(res(c).eRes, [], 1)); fprintf('\n');
end
fprintf('  Reading: a defect-free lossless solver gives residuals at round-off or falling with M.\n\n');

%% ---- Part B: ladders ----------------------------------------------------------------------------------------------
rel = @(a, b) norm(a(:) - b(:)) / norm(b(:));
iM = find(cfgs(:, 2) == Pfix & ismember(cfgs(:, 1), Mladder));
[~, o] = sort(cfgs(iM, 1)); iM = iM(o);
iP = find(cfgs(:, 1) == Mfix & ismember(cfgs(:, 2), Pladder));
[~, o] = sort(cfgs(iP, 2)); iP = iP(o);

fprintf('PART B1  M ladder at P = %d: relative change of each quantity against the finest M (%d)\n', Pfix, cfgs(iM(end), 1));
fprintf('  %6s %12s %12s %12s %12s %12s\n', 'M', 'T (all)', 'eta (all)', '|eta_s26|', 'U', 'b/a');
ref = res(iM(end));
for i = iM(:).'
    fprintf('  %6d %12.3e %12.3e %12.3e %12.3e %12.3e\n', res(i).M, rel(res(i).T, ref.T), rel(res(i).eta, ref.eta), ...
        rel(abs(res(i).eta(refIdx, :)), abs(ref.eta(refIdx, :))), rel(res(i).U, ref.U), rel(res(i).ba, ref.ba));
end
fprintf('\n  Successive doubling steps s(M) = |f(2M) - f(M)| / |f|, and fitted rate p = log2(s(M)/s(2M)):\n');
fprintf('  %12s %12s %8s %12s %8s %12s %8s\n', 'M -> 2M', 's_T', 'p_T', 's_U', 'p_U', 's_ba', 'p_ba');
sT = []; sU = []; sB = [];
for j = 1:numel(iM) - 1
    a = res(iM(j)); b = res(iM(j + 1));
    sT(j) = rel(a.T, b.T); sU(j) = rel(a.U, b.U); sB(j) = rel(a.ba, b.ba); %#ok<SAGROW>
end
for j = 1:numel(sT)
    if j < numel(sT)
        pT = log2(sT(j) / sT(j + 1)); pU = log2(sU(j) / sU(j + 1)); pB = log2(sB(j) / sB(j + 1));
        fprintf('  %5d->%-6d %12.3e %8.2f %12.3e %8.2f %12.3e %8.2f\n', res(iM(j)).M, res(iM(j + 1)).M, sT(j), pT, sU(j), pU, sB(j), pB);
    else
        fprintf('  %5d->%-6d %12.3e %8s %12.3e %8s %12.3e %8s\n', res(iM(j)).M, res(iM(j + 1)).M, sT(j), '', sU(j), '', sB(j), '');
    end
end
% per-frequency, per-transmissibility relative change, finest step
a = res(iM(end - 1)); b = res(iM(end));
fprintf('\n  Last M step, |dT|/|T| per transmissibility and frequency:\n');
fprintf('  %-10s', ''); fprintf('%11.1f', omega); fprintf('\n');
for j = 1:3
    fprintf('  %-10s', [sensorNames{j} '/s26']);
    fprintf('%11.2e', abs(a.T(j, :) - b.T(j, :)) ./ abs(b.T(j, :))); fprintf('\n');
end
fprintf('  %-10s', '|eta_s26|'); fprintf('%11.2e', abs(abs(a.eta(refIdx, :)) - abs(b.eta(refIdx, :))) ./ abs(b.eta(refIdx, :))); fprintf('\n');

% Richardson extrapolation of T using the fitted rate of the last two steps
if numel(sT) >= 2
    pFit = log2(sT(end - 1) / sT(end));
    Tlast = res(iM(end)).T; Tprev = res(iM(end - 1)).T;
    Tinf = Tlast + (Tlast - Tprev) / (2^pFit - 1);
    fprintf('\n  Richardson (p = %.2f): |T_inf - T_M%d| / |T| = %.3e  (estimated remaining truncation error at M = %d)\n', ...
        pFit, res(iM(end)).M, rel(Tlast, Tinf), res(iM(end)).M);
    fprintf('  Production M = 50 against T_inf: %.3e\n', rel(res(iM(cfgs(iM, 1) == 50)).T, Tinf));
else
    Tinf = []; pFit = NaN;
end

fprintf('\nPART B2  P ladder at M = %d: relative change against the finest P (%d)\n', Mfix, cfgs(iP(end), 2));
fprintf('  %6s %12s %12s %12s\n', 'P', 'T (all)', 'eta (all)', 'b/a');
refP = res(iP(end));
for i = iP(:).'
    fprintf('  %6d %12.3e %12.3e %12.3e\n', res(i).P, rel(res(i).T, refP.T), rel(res(i).eta, refP.eta), rel(res(i).ba, refP.ba));
end
fprintf('  Reading: with the correct corner weight, T should settle within a few P once M is large.\n\n');

%% ---- Part C: Bennetts reduction ------------------------------------------------------------------------------------
fprintf('PART C  Bennetts reduction phi_{-j}(z) ~ sum_m V(m,j) phi_m(z), at omega = %.1f rad/s\n', omega(min(3, nW)));
wC = min(3, nW);
z = linspace(-1, -gam, 4001);
fprintf('  %6s %14s %14s %16s %16s\n', 'M', 'L2 rel j=-2', 'L2 rel j=-1', 'dphi(-gam) j=-2', 'dphi(-gam) j=-1');
Cerr = zeros(numel(Mladder), 4);
for i = 1:numel(Mladder)
    M = Mladder(i);
    kap = rts(wC).kappa(1:M + 3);
    V = bennettsCoefficients(alphaW(wC), bet, gam, kap);
    tower = kap(3:end);
    Phi = zeros(numel(tower), numel(z));
    for m = 1:numel(tower), Phi(m, :) = verticalEigenfunction(z, tower(m), gam); end
    dTower = edgeDerivative(tower, gam);
    for j = 1:2
        damped = kap(j);                              % column j of V <-> kappa(j)  (kappa_-2, kappa_-1)
        exact = verticalEigenfunction(z, damped, gam);
        approx = V(:, j).' * Phi;
        Cerr(i, j) = sqrt(trapz(z, abs(approx - exact).^2) / trapz(z, abs(exact).^2));
        dEx = edgeDerivative(damped, gam);
        Cerr(i, 2 + j) = abs(V(:, j).' * dTower(:) - dEx) / abs(dEx);
    end
    fprintf('  %6d %14.3e %14.3e %16.3e %16.3e\n', M, Cerr(i, :));
end
if numel(Mladder) >= 3
    pC = log2(Cerr(end - 1, 1:2) ./ Cerr(end, 1:2));
    fprintf('  fitted L2 rate p (last doubling): %.2f, %.2f\n', pC);
end
fprintf('  Reading: the derivative column is expected NOT to converge (reduction is integration-valid only).\n');
fprintf('  Compare the L2 rate with p_T in Part B1.\n\n');

%% ---- Part D: exact damped-mode projection ----------------------------------------------------------------------
fprintf('PART D  M ladder at P = %d with the damped modes projected exactly (no Bennetts series)\n', Pfix);
resX = struct('M', {}, 'T', {}, 'eta', {}, 'eRes', {});
for i = 1:numel(Mladder)
    M = Mladder(i);
    etaX = complex(zeros(numel(sensR), nW)); eX = zeros(2 * N + 1, nW);
    for w = 1:nW
        out = solveAllOrdersExact(alphaW(w), bet, gam, R, nu, M, Pfix, N, rts(w).kappa(1:M + 3), rts(w).k(1:M + 1));
        e = evaluateDeflection(out.data, sensR, sensTh); etaX(:, w) = e(:);
        eX(:, w) = abs(abs(1 + 2 * out.ba) - 1);
    end
    resX(i) = struct('M', M, 'T', etaX(1:3, :) ./ repmat(etaX(refIdx, :), 3, 1), 'eta', etaX, 'eRes', eX);
end
fprintf('  max energy residual (exact variant): %.2e\n', max(arrayfun(@(r) max(r.eRes(:)), resX)));
fprintf('  %12s %12s %8s %12s %8s\n', 'M -> 2M', 's_T exact', 'p', 's_T reduced', 'p');
sX = zeros(1, numel(resX) - 1);
for j = 1:numel(resX) - 1, sX(j) = rel(resX(j).T, resX(j + 1).T); end
for j = 1:numel(sX)
    if j < numel(sX)
        fprintf('  %5d->%-6d %12.3e %8.2f %12.3e %8.2f\n', resX(j).M, resX(j + 1).M, sX(j), log2(sX(j) / sX(j + 1)), sT(j), log2(sT(j) / sT(j + 1)));
    else
        fprintf('  %5d->%-6d %12.3e %8s %12.3e\n', resX(j).M, resX(j + 1).M, sX(j), '', sT(j));
    end
end
fprintf('  reduced vs exact at the finest M: %.3e (relative, T)\n\n', rel(res(iM(end)).T, resX(end).T));

%% ---- Part E: radial probe ---------------------------------------------------------------------------------------
fprintf('PART E  Radial probe on the M ladder (P = %d): max over frequency of |eta(2M) - eta(M)| / |eta(2M)|\n', Pfix);
fprintf('  %12s', 'r/R'); fprintf('%9.2f', probeRR); fprintf('\n');
sE = zeros(numel(iM) - 1, numel(probeRR));
for j = 1:numel(iM) - 1
    a = res(iM(j)).probe; b = res(iM(j + 1)).probe;
    sE(j, :) = max(abs(b - a) ./ abs(b), [], 2).';
    fprintf('  %5d->%-6d', res(iM(j)).M, res(iM(j + 1)).M); fprintf('%9.1e', sE(j, :)); fprintf('\n');
end
fprintf('  rate p per point:\n');
for j = 1:size(sE, 1) - 1
    fprintf('  %12s', sprintf('%d/%d', res(iM(j)).M, res(iM(j + 1)).M)); fprintf('%9.2f', log2(sE(j, :) ./ sE(j + 1, :))); fprintf('\n');
end
fprintf('  Reading: p should stay roughly constant down each column. A negative p in the last row means the\n');
fprintf('  changes stopped shrinking (a plateau), and Richardson extrapolation is NOT safe.\n\n');

%% ---- save ------------------------------------------------------------------------------------------------------
meta = struct('beta', bet, 'gamma', gam, 'R', R, 'nu', nu, 'depth', depth, 'omega', omega, 'alpha', alphaW, ...
    'sensorsRR', sensRR, 'sensorNames', {sensorNames}, 'ref', refIdx, 'N', N, 'Mladder', Mladder, 'Pfix', Pfix, ...
    'Pladder', Pladder, 'Mfix', Mfix, 'quick', QUICK, 'platform', ternary(isOctave, 'Octave', 'MATLAB'), ...
    'when', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
save(fullfile(resDir, 'diag_F1_convergence.mat'), 'meta', 'res', 'resX', 'Cerr', 'Tinf', 'pFit', 'cfgs', 'sE', 'probeRR', '-v7');
fprintf('Saved %s\n', fullfile(resDir, 'diag_F1_convergence.mat'));
fprintf('==== done ====\n');
diary('off');

%% ==== local functions ==============================================================================================
function out = solveAllOrders(alpha, beta, gam, R, nu, M, P, N, kappa, k)
% Same pipeline as precomputeDeflectionData.m (production), with the roots passed in and the
% interface coefficients U and the open-water scattered amplitudes kept for the diagnostics.
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
    U = complex(zeros(P + 1, 2 * N + 1));
    ba = complex(zeros(2 * N + 1, 1));
    p0 = imag(k0);                                % k0 = 1i * p0, p0 > 0
    for nIdx = 1:(2 * N + 1)
        n = nIdx - N - 1;
        Ln = edgeMap(n, alpha, beta, gam, R, nu, kappa);
        [EIprime, E0prime] = openWaterRadialDerivatives(n, R, k);
        Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);
        [MIn, M0, Mmat] = matchingMatrices(EIprime, E0prime, Enprime, D0, D, D0g, Dg);
        En_R = radialMatrix(n, R, R, kappa, V, Ln);
        [MU, Sn, Sn0] = scatteringMatrices(D0g, Dg, En_R, MIn, M0, Mmat);
        AIn = incidentAmplitude(alpha, n, k0, R, M);
        An = Sn * AIn;
        Aen = Ln * An;
        AfullMat(:, nIdx) = [Aen(1); Aen(2); An];
        U(:, nIdx) = MU * AIn;
        A0 = Sn0 * AIn;
        % propagating open-water mode: AIn(1) I_n(k0 r)/I_n(k0 R) + A0(1) K_n(k0 r)/K_n(k0 R)
        %   = a J_n(p0 r) + b H2_n(p0 r),  a = AIn(1)/J_n(p0 R),  b = A0(1)/H2_n(p0 R)
        a = AIn(1) / besselj(n, p0 * R);
        b = A0(1) / besselh(n, 2, p0 * R);
        ba(nIdx) = b / a;
    end
    out.data = struct('alpha', alpha, 'R', R, 'N', N, 'kappa', kappa, 'dphiFull', dphiFull, 'AfullMat', AfullMat);
    out.U = U;
    out.Atower = AfullMat(3:end, :);
    out.ba = ba;
end

function out = solveAllOrdersExact(alpha, beta, gam, R, nu, M, P, N, kappa, k)
% As solveAllOrders, but the damped modes are projected EXACTLY instead of through the truncated
% Bennetts series phi_{-j} ~ sum V(:,j) phi_m:
%   velocity rows:  D*V    -> C   = <phi_i, phi_{-j}>          (closed form, innerProductMatrix)
%   pressure rows:  Dg'*V  -> Dgd' (Gegenbauer projections of phi_{-j}, gegenbauerCrossMatrix)
    k0 = k(1);
    tower = kappa(3:end); damped = kappa(1:2);
    D0  = innerProductMatrix(k, 0);
    Dall = innerProductMatrix([tower; damped], gam);
    D   = Dall(1:end-2, 1:end-2);
    C   = Dall(1:end-2, end-1:end);
    D0g = gegenbauerCrossMatrix(k, gam, P, 1);
    Dg  = gegenbauerCrossMatrix(tower, gam, P, 1 - gam);
    Dgd = gegenbauerCrossMatrix(damped, gam, P, 1 - gam);
    nModes = numel(kappa);
    dphiFull = zeros(nModes, 1);
    for idx = 1:nModes, dphiFull(idx) = edgeDerivative(kappa(idx), gam); end
    AfullMat = complex(zeros(nModes, 2*N+1)); U = complex(zeros(P+1, 2*N+1)); ba = complex(zeros(2*N+1,1));
    p0 = imag(k0); M1 = numel(tower);
    for nIdx = 1:(2*N+1)
        n = nIdx - N - 1;
        Ln = edgeMap(n, alpha, beta, gam, R, nu, kappa);
        [EIprime, E0prime] = openWaterRadialDerivatives(n, R, k);
        Ep  = zeros(M1,1); for m = 1:M1, Ep(m) = scaledBesselDerivative('I', n, tower(m), R, R); end
        Epd = zeros(2,1);  for j = 1:2,  Epd(j) = scaledBesselDerivative('I', n, damped(j), R, R); end
        MIn = -(E0prime \ EIprime);
        M0  = (D0*E0prime) \ D0g;
        Avel = D*diag(Ep) + C*diag(Epd)*Ln;           % exact plate-side velocity projection
        Mmat = Avel \ Dg;
        Bp = Dg.' + Dgd.'*Ln;                          % exact plate-side pressure projection (Ihat(R)=1)
        Aeq = Bp*Mmat - D0g.'*M0;
        MU = Aeq \ (D0g.'*(eye(M1) + MIn));
        Sn = Mmat*MU; Sn0 = MIn + M0*MU;
        AIn = incidentAmplitude(alpha, n, k0, R, M);
        An = Sn*AIn; Aen = Ln*An;
        AfullMat(:, nIdx) = [Aen(1); Aen(2); An];
        U(:, nIdx) = MU*AIn;
        A0 = Sn0*AIn;
        ba(nIdx) = (A0(1)/besselh(n,2,p0*R)) / (AIn(1)/besselj(n,p0*R));
    end
    out.data = struct('alpha',alpha,'R',R,'N',N,'kappa',kappa,'dphiFull',dphiFull,'AfullMat',AfullMat);
    out.U = U; out.Atower = AfullMat(3:end,:); out.ba = ba;
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
