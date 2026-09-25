%% level4A_broadbandKnownOscillator.m
% LEVEL 4A: does broadband stochastic excitation, synthesised in the
% frequency domain (stochasticSensorSynthesis.m), let the existing SSI
% chain recover the DAMPED poles of a known system -- or does SSI model
% the finite set of synthesis tones as undamped oscillators instead?
%
% This is the gate in front of the EMM bridge (runBroadbandSSIBridgeTest).
% No EMM solves: the "floe" is an analytical modal model with known
% (omega_n, zeta, mode shapes), observed as ACCELERATION.
%
% ROUTES (same system, same record length N*DT, same SSI settings)
%   synthGaussian    : stochasticSensorSynthesis, circular complex Gaussian
%                      input, band BAND with a cosine taper, PLUS white
%                      measurement noise. THE ROUTE UNDER TEST.
%   synthRandomPhase : same, but deterministic amplitudes and random phases
%                      (the older forced-tone style), plus noise. Comparison.
%   timeDomain       : direct time-domain simulation of the same modal
%                      system driven by sampled white noise (exact ZOH
%                      discretisation, burn-in discarded), plus the same
%                      measurement noise. CONTROL.
%   timeBandpassed   : the time-domain output, band-passed in the FFT domain
%                      with EXACTLY the synthesis band and taper, plus the
%                      same noise. Same spectrum as synthGaussian but a
%                      genuinely continuous (non-tone) input. SECOND CONTROL.
%   synthNoiseless   : synthGaussian WITHOUT measurement noise. Expected to
%                      FAIL (see below); kept to document why noise is needed.
%   synthGaussianCont: synthGaussian built from the CONTINUOUS FRF H_a(i w)
%                      instead of the exact ZOH FRF. Shows the size of the
%                      sampling-model difference (the EMM FRF is continuous).
%
% SAMPLING MODEL (revised after the first MATLAB run, 20260925_104125)
% The synthesis routes use the EXACT discrete-time FRF of the same ZOH
% system the time-domain routes simulate,
%     H_d(w) = C (e^{i w DT} I - A_d)^{-1} B_d + D,
% so synthGaussian and timeBandpassed describe the same sampled system
% and differ only in the input mechanism (periodic tone set vs a
% genuinely non-periodic white sequence, band-passed by FFT). In the first
% run the synthesis used the continuous H_a(i w), so the Gate A difference
% (+0.169 for twoClose mode 1) mixed synthesis and sampling-model effects.
% Minor residual: FFT band-passing a non-periodic record wraps a few
% correlation times (~1/(zeta*w_n) = 3-5 s) at the record ends.
%
% READING THE ROUTES AS A DECOMPOSITION (at equal record length)
%   timeDomain  vs timeBandpassed : effect of BAND LIMITATION (+ noise floor)
%   timeBandpassed vs synthGaussian: effect of the SYNTHESIS itself (tones,
%                                    circularity, Gaussian bin amplitudes)
% The gate is about the second difference; the first is a property of any
% band-limited twin and is reported, not gated.
%
% WHY MEASUREMENT NOISE IS PART OF THE ROUTE UNDER TEST (found while
% building this script, then designed in)
% With DT = 0.1 s the Nyquist frequency is 31.4 rad/s, but the synthesis
% excites only BAND = [3 9]: the sampled output PSD is EXACTLY zero on
% most of [0, pi/DT]. A stationary process whose spectrum vanishes on a
% set of positive measure has integral(log S) = -Inf, so by the Szego-
% Kolmogorov formula its one-step prediction error is ZERO: it is
% deterministic (Paley-Wiener). SSI's stochastic model has no innovation
% to work with and fits the (perfectly predictable) signal with
% near-undamped poles. Development checks, single mode, N = 4096:
%   full-band synthesis (0 to Nyquist)         -> zeta recovered (0.059 vs 0.050)
%   band [3 9], with or without taper         -> zeta ~ 0.001, no mode found
%   time-domain white noise, ideally band-
%   passed to [3 9] after simulation          -> fails identically
%   band [3 9] + white noise at 1% of std     -> zeta recovered (6/6)
% With 30 realisations at N = 4096 (1% noise), median zeta_hat/zeta was
% 0.92 (synthGaussian), 0.96 (synthRandomPhase), 0.98 (timeDomain) and
% 0.76 (timeBandpassed); IQRs are ~0.3 wide, hence N_REAL = 30 below.
% So the failure is BAND LIMITATION, not the tone discretisation, and any
% broadband floor removes it. The added white noise is best read as a
% CONTROLLED INNOVATION FLOOR: in real data that floor comes from sensor
% noise, excitation outside the EMM band, unmodelled dynamics and
% environmental disturbance together, and it will be replaced by the
% measured IMU noise spectrum later. NOISE_FRACTION sets it; Stage 4
% sweeps it (including zero). Stage 5 sweeps the BAND WIDTH, which says
% how wide the EMM FRF sweep must be for the band-limited twin to stop
% biasing damping (the EMM truncation is validated only to alpha ~ 13.85,
% i.e. ~8.5 rad/s, so this sets the cost of widening it).
%
% WHY THE CONTROL MATTERS
% A synthesised record is exactly one period of a set of tones spaced
% dOmega = 2*pi/(N*DT) apart, so refining the grid and lengthening the
% record are the same thing. The time-domain route has the SAME length
% but a continuous input spectrum. Reading the two together:
%   both fail at this N         -> record too short for SSI (not synthesis)
%   time-domain passes, synth
%   fails at the same N         -> the tone discretisation is the cause
%   both pass                   -> synthesis is adequate at this N_BW
% N_BW = 2*zeta*omega_n / dOmega = number of synthesis bins across the
% full half-power bandwidth. No transition value is assumed; this script
% measures it.
%
% DIFFERENCES BETWEEN ROUTES, acknowledged: the synthesis routes are band
% limited to BAND (as the EMM FRF will be), the time-domain route is
% broadband up to Nyquist and includes the acceleration feedthrough term.
% Neither changes the system poles.
%
% CASES
%   single   : omega_n = 6.0, zeta = 0.05, 3 sensors.
%   twoClose : omega_n = [5.5 6.3], zeta = [0.04 0.06]. The half-power
%              bands (0.44 and 0.76 rad/s wide) nearly touch: a harder,
%              EMM-like case.
% Mode shapes are real (classical damping), so the reference-sensor
% ratios q_j = Phi_j / Phi_ref are real; SSI's complex shapes are compared
% through the same ratio, which removes SSI's arbitrary complex scaling.
%
% POLE SELECTION (per realisation)
% Existing chain: buildHankelMatrix -> hankelProjection -> sweepModelOrders
% -> classifyPoleStability -> buildModalBranches -> computeBranchPersistence.
% Each PERSISTENT branch is assigned to the nearest true mode (if its
% median frequency is within ASSIGN_TOL of it); for each mode the branch
% with the most class-3 transitions wins (ties: longest, then closest).
% Estimates are medians over that branch's class-3 orders (all its orders
% if it has none). Persistent in-band branches not assigned are counted
% as spurious; persistent in-band branches with damping below
% UNDAMPED_FRACTION * min(zeta_true) are counted as tone-like.
%
% SPEED: sweepModelOrders is handed R21 (from hankelProjection's LQ factor)
% instead of P = R21*Q1'. Since Q1 has orthonormal columns, P*P' = R21*R21',
% so U and S of the SVD (all extractStateSpace uses) are identical, but the
% SVD is of an (l*i) x (r*i) matrix instead of (l*i) x j. Stage 0 checks
% this equivalence numerically before anything else runs.
%
% MEASUREMENT NOISE: independent white Gaussian per sensor and sample,
% std = NOISE_FRACTION * std(noise-free record, all channels pooled), i.e.
% the same absolute noise level on every sensor (same IMU type), drawn
% from its own seed.
%
% REFERENCE SENSOR per mode = the sensor with the largest |true shape|
% component, frozen before identification, so q_j = Phi_j / Phi_ref never
% divides by a small number and the shape error is relative to 1.
%
% OUTPUT (in SSI/Validation Experiments/Results/)
%   level4A_<stamp>.mat  : config + every per-realisation estimate
%   level4A_<stamp>.txt  : console log
%   level4A_<stamp>_damping.png, level4A_<stamp>_stabilisation.png,
%   level4A_<stamp>_noise.png
%
% Lives in SSI/Validation Experiments/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));                 % SSI chain + stochasticSensorSynthesis
clearvars -except thisDir
close all, clc

%% ---- Configuration -------------------------------------------------------
DT        = 0.1;                                   % s  (Nyquist 31.4 rad/s)
N_LIST    = [128 256 512 1024 2048 4096 8192];     % record length = N*DT = 2*pi/dOmega
N_REAL    = 30;                                    % realisations per (case, route, N); median IQR ~0.3
ROUTES    = {'synthGaussian', 'synthRandomPhase', 'timeDomain', 'timeBandpassed', 'synthNoiseless', ...
             'synthGaussianCont'};
NOISE_FRACTION = 0.01;                             % measurement noise std / signal std
NOISE_LIST = [0 1e-3 3e-3 1e-2 3e-2 1e-1];         % Stage 4 sweep (synthGaussian, N = N_LIST(end-1))
BAND_LIST = {[3 9], [2 10], [1 12], [0.5 15], [0 Inf]};   % Stage 5 sweep; [0 Inf] = everything below Nyquist
BAND      = [3 9];                                 % rad/s, synthesis band
TAPER_FRACTION = 0.1;
SU        = 1;                                     % one-sided input PSD per rad/s
BURN_IN   = 2000;                                  % samples discarded (time-domain route)

I_BLOCK   = 20;                                    % Hankel block rows (2 s)
N_MAX     = 20;                                    % highest model order
FREQ_TOL  = 0.01;  DAMP_TOL = 0.10;  MAC_TOL = 0.95;
PERSIST_W = 3;     PERSIST_K = 2;
ASSIGN_TOL = 0.10;                                 % relative: branch -> true mode
UNDAMPED_FRACTION = 0.2;

% pass criteria, applied to the MEDIAN over realisations
PASS.detectMin   = 0.9;                            % fraction of realisations mode found
PASS.freqBiasMax = 0.005;                          % |median w_hat / w - 1|
PASS.dampBiasMax = 0.15;                           % |median z_hat / z - 1|
PASS.shapeErrMax = 0.15;                           % median of max_j |q_hat_j - q_j| (q_ref = 1)
PASS.synthVsBandpassedMax = 0.15;                  % |median z ratio, synthGaussian - timeBandpassed|

cases(1).name = 'single';
cases(1).wn   = 6.0;
cases(1).zeta = 0.05;
cases(1).Phi  = [1; -0.6; 0.3];
cases(2).name = 'twoClose';
cases(2).wn   = [5.5 6.3];
cases(2).zeta = [0.04 0.06];
cases(2).Phi  = [1 0.3; 0.5 -1; -0.4 0.8];

%% ---- Output ------------------------------------------------------------------
resultsDir = fullfile(thisDir, 'Results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['level4A_', stamp]);
diary([baseName, '.txt']);
try   % close the diary even if the run errors
warnState = warning('off', 'modalParameters:realEigenvalue');
tStart = tic;

fprintf('level4A_broadbandKnownOscillator  (%s)\n', stamp);
fprintf('DT = %.3g s, band [%.2g %.2g] rad/s, taper %.2g, i = %d, nMax = %d, %d realisations\n', ...
    DT, BAND, TAPER_FRACTION, I_BLOCK, N_MAX, N_REAL);
fprintf('tolerances: freq %.3g, damp %.3g, MAC %.3g, persistence %d of %d\n', ...
    FREQ_TOL, DAMP_TOL, MAC_TOL, PERSIST_K, PERSIST_W);
fprintf('measurement noise: %.3g x signal std (all routes except synthNoiseless)\n\n', NOISE_FRACTION);

%% ---- Stage 0: R21 shortcut reproduces the P-based sweep ----------------------------
fprintf('=== Stage 0: sweepModelOrders on R21 vs on P ===\n');
[yChk, ~] = synthesiseCase(cases(2), 'synthGaussian', 2048, DT, BAND, TAPER_FRACTION, SU, BURN_IN, NOISE_FRACTION, 7);
[Yp, Yf] = buildHankelMatrix(yChk, I_BLOCK);
[Pfull, ~, Rlq] = hankelProjection(Yp, Yf);
R21 = Rlq(size(Yp, 1) + 1:end, 1:size(Yp, 1));
nChk = 12;
resP = sweepModelOrders(Pfull, size(yChk, 1), DT, nChk);
resR = sweepModelOrders(R21, size(yChk, 1), DT, nChk);
worst = 0;
for n = 1:nChk
    assert(numel(resP(n).frequency) == numel(resR(n).frequency), 'Stage 0: pole count differs at order %d', n);
    worst = max([worst; abs(sort(resP(n).frequency) - sort(resR(n).frequency)) ./ sort(resP(n).frequency)]);
end
assert(worst < 1e-8, 'Stage 0: R21 shortcut changes the poles (rel %.2e)', worst);
fprintf('PASS: identical pole counts at orders 1-%d, max relative frequency difference %.1e\n\n', nChk, worst);
clear yChk Yp Yf Pfull Rlq R21 resP resR

fprintf('=== Stage 0b: ZOH FRF matches the simulated sampled system ===\n');
[Ad0, Bd0, Cc0, Dc0] = modalStateSpaceZOH(cases(2), DT);
w0 = 6.1; nT = 3000; x0 = zeros(size(Ad0, 1), 1); yT0 = zeros(size(Cc0, 1), nT);
for n = 1:nT                                          % u_n = cos(w0 n DT), run to steady state
    un = cos(w0 * (n - 1) * DT);
    yT0(:, n) = Cc0 * x0 + Dc0 * un;
    x0 = Ad0 * x0 + Bd0 * un;
end
tt = (0:nT - 1) * DT;
yF0 = real(modalAccelerationFRFzoh(cases(2), w0, DT) * exp(1i * w0 * tt));
tail = nT - 499:nT;
dZ = max(max(abs(yT0(:, tail) - yF0(:, tail)))) / max(max(abs(yF0(:, tail))));
yC0 = real(modalAccelerationFRF(cases(2), w0) * exp(1i * w0 * tt));
dC = max(max(abs(yT0(:, tail) - yC0(:, tail)))) / max(max(abs(yF0(:, tail))));
assert(dZ < 1e-8, 'Stage 0b: ZOH FRF does not match the simulation (rel %.2e)', dZ);
fprintf('PASS: ZOH FRF vs simulation %.1e; continuous FRF vs simulation %.2f (sampling-model gap at %.1f rad/s)\n\n', ...
    dZ, dC, w0);
clear Ad0 Bd0 Cc0 Dc0 x0 yT0 yF0 yC0 tt tail

%% ---- Stage 1: sweep -----------------------------------------------------------------
nC = numel(cases); nRt = numel(ROUTES); nN = numel(N_LIST); mMax = 2;
est.w        = NaN(nC, nRt, nN, N_REAL, mMax);
est.z        = NaN(nC, nRt, nN, N_REAL, mMax);
est.qErr     = NaN(nC, nRt, nN, N_REAL, mMax);
est.detected = false(nC, nRt, nN, N_REAL, mMax);
est.nSpurious = zeros(nC, nRt, nN, N_REAL);
est.nToneLike = zeros(nC, nRt, nN, N_REAL);
stabExample = struct('caseName', {}, 'route', {}, 'N', {}, 'results', {}, 'branches', {}, ...
    'persistence', {}, 'binOmega', {});

tol = struct('freq', FREQ_TOL, 'damp', DAMP_TOL, 'mac', MAC_TOL, 'W', PERSIST_W, 'K', PERSIST_K, ...
    'assign', ASSIGN_TOL, 'undamped', UNDAMPED_FRACTION, 'band', BAND);

for c = 1:nC
    cs = cases(c); nM = numel(cs.wn);
    for iN = 1:nN
        N = N_LIST(iN);
        dW = 2 * pi / (N * DT);
        for iR = 1:nRt
            for r = 1:N_REAL
                seed = 100000 * c + 1000 * iN + r;       % same seed across routes
                [y, binOmega] = synthesiseCase(cs, ROUTES{iR}, N, DT, BAND, TAPER_FRACTION, SU, BURN_IN, ...
                    NOISE_FRACTION, seed);
                [results, branches, persistence] = runSSIChain(y, DT, I_BLOCK, N_MAX, tol);
                e = assignModes(results, branches, persistence, cs, tol);
                est.w(c, iR, iN, r, 1:nM) = e.w;
                est.z(c, iR, iN, r, 1:nM) = e.z;
                est.qErr(c, iR, iN, r, 1:nM) = e.qErr;
                est.detected(c, iR, iN, r, 1:nM) = e.detected;
                est.nSpurious(c, iR, iN, r) = e.nSpurious;
                est.nToneLike(c, iR, iN, r) = e.nToneLike;
                if r == 1 && c == 2 && any(N == [N_LIST(1), N_LIST(ceil(nN / 2)), N_LIST(nN)])
                    stabExample(end + 1) = struct('caseName', cs.name, 'route', ROUTES{iR}, 'N', N, ... %#ok<SAGROW>
                        'results', results, 'branches', branches, 'persistence', persistence, ...
                        'binOmega', binOmega);
                end
            end
        end
        fprintf('  %-9s N = %5d (T = %6.1f s, dOmega = %.4f)  done  [%.0f s elapsed]\n', ...
            cs.name, N, N * DT, dW, toc(tStart));
    end
end

%% ---- Stage 2: summary tables ----------------------------------------------------------
blank = struct('nBW', NaN(1, nN), 'det', NaN(1, nN), 'wBias', NaN(1, nN), 'zMed', NaN(1, nN), ...
    'zIQR', NaN(nN, 2), 'qMed', NaN(1, nN), 'pass', false(1, nN));
summ = repmat(struct('route', {repmat({blank}, 1, nRt)}), nC, mMax);
passTab = false(nC, nRt, nN);
fprintf('\n=== Stage 2: medians over %d realisations ===\n', N_REAL);
fprintf('det = detection rate, w = median w_hat/w - 1, z = median z_hat/z (IQR), q = median shape error,\n');
fprintf('spur/tone = mean persistent in-band spurious / tone-like branches per realisation\n');
for c = 1:nC
    cs = cases(c); nM = numel(cs.wn);
    for m = 1:nM
        fprintf('\n--- case %s, mode %d: omega_n = %.2f, zeta = %.3f ---\n', cs.name, m, cs.wn(m), cs.zeta(m));
        fprintf('%-17s %6s %7s %6s %10s %17s %8s %6s %6s %5s\n', 'route', 'N', 'N_BW', 'det', ...
            'w bias', 'z ratio (IQR)', 'q err', 'spur', 'tone', 'pass');
        for iR = 1:nRt
            for iN = 1:nN
                N = N_LIST(iN);
                nBW = 2 * cs.zeta(m) * cs.wn(m) * N * DT / (2 * pi);
                det = mean(squeeze(est.detected(c, iR, iN, :, m)));
                wr = squeeze(est.w(c, iR, iN, :, m)) / cs.wn(m) - 1;
                zr = squeeze(est.z(c, iR, iN, :, m)) / cs.zeta(m);
                qe = squeeze(est.qErr(c, iR, iN, :, m));
                wB = nanMedian(wr); zM = nanMedian(zr); qM = nanMedian(qe);
                zQ = nanQuantiles(zr, [0.25 0.75]);
                ok = det >= PASS.detectMin && abs(wB) <= PASS.freqBiasMax && ...
                     abs(zM - 1) <= PASS.dampBiasMax && qM <= PASS.shapeErrMax;
                if m == 1
                    passTab(c, iR, iN) = ok;
                else
                    passTab(c, iR, iN) = passTab(c, iR, iN) && ok;
                end
                summ(c, m).route{iR}.nBW(iN) = nBW;
                summ(c, m).route{iR}.det(iN) = det;
                summ(c, m).route{iR}.wBias(iN) = wB;
                summ(c, m).route{iR}.zMed(iN) = zM;
                summ(c, m).route{iR}.zIQR(iN, :) = zQ;
                summ(c, m).route{iR}.qMed(iN) = qM;
                summ(c, m).route{iR}.pass(iN) = ok;
                fprintf('%-17s %6d %7.1f %6.2f %+10.2e %7.3f (%.2f-%.2f) %8.3f %6.2f %6.2f %5s\n', ...
                    ROUTES{iR}, N, nBW, det, wB, zM, zQ(1), zQ(2), qM, ...
                    mean(est.nSpurious(c, iR, iN, :)), mean(est.nToneLike(c, iR, iN, :)), yesNo(ok));
            end
        end
    end
end

%% ---- Stage 3: verdict ------------------------------------------------------------------
fprintf('\n=== Stage 3: verdict ===\n');
iG = find(strcmp(ROUTES, 'synthGaussian'));
iT = find(strcmp(ROUTES, 'timeDomain'));
iZ = find(strcmp(ROUTES, 'synthNoiseless'));
iB = find(strcmp(ROUTES, 'timeBandpassed'));
controlOK = all(passTab(:, iT, nN));
bandLimitedOK = all(all(passTab(:, iG, nN - 1:nN)));
for c = 1:nC
    nBWmin = min(2 * cases(c).zeta .* cases(c).wn) * N_LIST * DT / (2 * pi);   % limiting mode
    fprintf('case %s:\n', cases(c).name);
    for iR = 1:nRt
        p = squeeze(passTab(c, iR, :)).';
        firstStable = find(fliplr(cumprod(fliplr(double(p)))), 1);                     % from here on, all pass
        if ~any(p)
            fprintf('  %-17s passes at no N\n', ROUTES{iR});
        else
            fprintf('  %-17s passes individually at N = %s', ROUTES{iR}, mat2str(N_LIST(p)));
            if isempty(firstStable)
                fprintf('; no N from which every longer record passes\n');
            else
                fprintf('; every N >= %d passes (N_BW >= %.1f, narrowest mode)\n', ...
                    N_LIST(firstStable), nBWmin(firstStable));
            end
        end
    end
    if any(squeeze(passTab(c, iZ, :)))
        fprintf('  NOTE: synthNoiseless passed somewhere -- unexpected given the band-limitation argument; inspect.\n');
    end
    discret = find(squeeze(passTab(c, iT, :)) & ~squeeze(passTab(c, iG, :)));
    if ~isempty(discret)
        fprintf('  time-domain passes but synthGaussian fails at N = %s (band limitation and/or synthesis;\n', ...
            mat2str(N_LIST(discret)));
        fprintf('  Gate A below separates the two)\n');
    end
end
fprintf('\nGate A, synthesis fidelity: median zeta ratio, synthGaussian minus timeBandpassed\n');
fprintf('(same sampled system, band, taper and noise; only the input mechanism differs)\n');
faithfulOK = true;
for c = 1:nC
    for m = 1:numel(cases(c).wn)
        d = summ(c, m).route{iG}.zMed - summ(c, m).route{iB}.zMed;
        zG = est.z(c, iG, nN - 1:nN, :, m); zB = est.z(c, iB, nN - 1:nN, :, m);
        dPool = nanMedian(zG(:) / cases(c).zeta(m)) - nanMedian(zB(:) / cases(c).zeta(m));
        okM = abs(dPool) <= PASS.synthVsBandpassedMax;
        faithfulOK = faithfulOK && okM;
        fprintf('  %-9s mode %d: per N %s | pooled two longest %+6.3f %s\n', cases(c).name, m, ...
            sprintf('%+6.3f ', d), dPool, yesNo(okM));
    end
end
fprintf('  (per-N columns are N = %s; gate uses the pooled value, tolerance %.2f)\n', ...
    mat2str(N_LIST), PASS.synthVsBandpassedMax);
iC = find(strcmp(ROUTES, 'synthGaussianCont'));
fprintf('Sampling model: median zeta ratio, synthGaussianCont (continuous FRF) minus synthGaussian (ZOH FRF)\n');
for c = 1:nC
    for m = 1:numel(cases(c).wn)
        fprintf('  %-9s mode %d: per N %s\n', cases(c).name, m, ...
            sprintf('%+6.3f ', summ(c, m).route{iC}.zMed - summ(c, m).route{iG}.zMed));
    end
end
fprintf('Gate B, band-limited recovery vs truth (synthGaussian, two longest N): %s\n', yesNo(bandLimitedOK));
fprintf('Control (time domain, longest N): %s\n', yesNo(controlOK));
gateOK = controlOK && faithfulOK;

if ~controlOK
    fprintf(['\nRESULT: CONTROL FAILED. The time-domain route does not pass at the longest record, so the\n' ...
        'SSI settings / pole selection are the problem, not the synthesis. Diagnose before using level4A.\n']);
elseif gateOK && bandLimitedOK
    fprintf(['\nRESULT: GATE PASSED. The synthesis is faithful (Gate A) and, with measurement noise, SSI recovers\n' ...
        'the damped poles and shape ratios from the band-limited record (Gate B). Use the N_BW transition\n' ...
        'above to choose dOmega for the EMM bridge (runBroadbandSSIBridgeTest.m).\n']);
elseif gateOK
    fprintf(['\nRESULT: GATE PASSED WITH A CAVEAT. The synthesis is faithful (Gate A: it behaves like a genuinely\n' ...
        'continuous input with the same band), but band limitation itself costs accuracy (Gate B fails):\n' ...
        'compare timeDomain with timeBandpassed above for zeta bias, spurious in-band branches and close-mode\n' ...
        'shape error. The EMM twin is band-limited by the truncation ceiling, so the bridge test inherits\n' ...
        'this cost; it is a property of the band-limited twin, not of the synthesis. Proceed, and report it.\n']);
else
    fprintf(['\nRESULT: GATE FAILED. The control passes but Gaussian synthesis does not (vs truth, or vs the\n' ...
        'band-passed control) at the longest records. Stop: diagnose the synthesis / SSI relationship\n' ...
        'before any EMM bridge test. If synthGaussian tracks timeBandpassed but both miss the truth, the\n' ...
        'cause is band limitation, not synthesis: widen BAND or raise the noise floor (Stage 4).\n']);
end

%% ---- Stage 4: measurement-noise sweep (synthGaussian) --------------------------------------
N4 = N_LIST(max(1, nN - 1));
nNz = numel(NOISE_LIST);
est4.w = NaN(nC, nNz, N_REAL, mMax); est4.z = est4.w; est4.qErr = est4.w;
est4.detected = false(nC, nNz, N_REAL, mMax);
fprintf('\n=== Stage 4: noise sweep, synthGaussian, N = %d ===\n', N4);
for c = 1:nC
    cs = cases(c); nM = numel(cs.wn);
    for iz = 1:nNz
        for r = 1:N_REAL
            seed = 900000 + 100000 * c + 1000 * iz + r;
            y = synthesiseCase(cs, 'synthGaussian', N4, DT, BAND, TAPER_FRACTION, SU, BURN_IN, NOISE_LIST(iz), seed);
            [results, branches, persistence] = runSSIChain(y, DT, I_BLOCK, N_MAX, tol);
            e = assignModes(results, branches, persistence, cs, tol);
            est4.w(c, iz, r, 1:nM) = e.w; est4.z(c, iz, r, 1:nM) = e.z;
            est4.qErr(c, iz, r, 1:nM) = e.qErr; est4.detected(c, iz, r, 1:nM) = e.detected;
        end
    end
    for m = 1:nM
        fprintf('case %s, mode %d (omega_n = %.2f, zeta = %.3f)\n', cs.name, m, cs.wn(m), cs.zeta(m));
        fprintf('  %8s %6s %10s %17s %8s\n', 'noise', 'det', 'w bias', 'z ratio (IQR)', 'q err');
        for iz = 1:nNz
            zr = squeeze(est4.z(c, iz, :, m)) / cs.zeta(m);
            zQ = nanQuantiles(zr, [0.25 0.75]);
            fprintf('  %8.0e %6.2f %+10.2e %7.3f (%.2f-%.2f) %8.3f\n', NOISE_LIST(iz), ...
                mean(squeeze(est4.detected(c, iz, :, m))), nanMedian(squeeze(est4.w(c, iz, :, m)) / cs.wn(m) - 1), ...
                nanMedian(zr), zQ(1), zQ(2), nanMedian(squeeze(est4.qErr(c, iz, :, m))));
        end
    end
end

%% ---- Stage 5: band-width sweep (synthGaussian, ZOH FRF) ------------------------------------
nBd = numel(BAND_LIST);
est5.w = NaN(nC, nBd, N_REAL, mMax); est5.z = est5.w; est5.qErr = est5.w;
est5.detected = false(nC, nBd, N_REAL, mMax);
est5.nSpurious = zeros(nC, nBd, N_REAL);
fprintf('\n=== Stage 5: band-width sweep, synthGaussian, N = %d, noise %.3g ===\n', N4, NOISE_FRACTION);
dW4 = 2 * pi / (N4 * DT);
bandsUsed = zeros(nBd, 2);
for ib = 1:nBd
    bnd = BAND_LIST{ib};
    bandsUsed(ib, :) = [max(bnd(1), dW4), min(bnd(2), pi / DT - dW4)];
end
for c = 1:nC
    cs = cases(c); nM = numel(cs.wn);
    for ib = 1:nBd
        tolB = tol; tolB.band = BAND;                  % spurious count stays over the ORIGINAL band
        for r = 1:N_REAL
            seed = 1900000 + 100000 * c + 1000 * ib + r;
            y = synthesiseCase(cs, 'synthGaussian', N4, DT, bandsUsed(ib, :), TAPER_FRACTION, SU, BURN_IN, ...
                NOISE_FRACTION, seed);
            [results, branches, persistence] = runSSIChain(y, DT, I_BLOCK, N_MAX, tolB);
            e = assignModes(results, branches, persistence, cs, tolB);
            est5.w(c, ib, r, 1:nM) = e.w; est5.z(c, ib, r, 1:nM) = e.z;
            est5.qErr(c, ib, r, 1:nM) = e.qErr; est5.detected(c, ib, r, 1:nM) = e.detected;
            est5.nSpurious(c, ib, r) = e.nSpurious;
        end
    end
    for m = 1:nM
        fprintf('case %s, mode %d (omega_n = %.2f, zeta = %.3f)\n', cs.name, m, cs.wn(m), cs.zeta(m));
        fprintf('  %-16s %6s %10s %17s %8s %6s\n', 'band (rad/s)', 'det', 'w bias', 'z ratio (IQR)', 'q err', 'spur');
        for ib = 1:nBd
            zr = squeeze(est5.z(c, ib, :, m)) / cs.zeta(m);
            zQ = nanQuantiles(zr, [0.25 0.75]);
            fprintf('  [%5.2f %6.2f]  %6.2f %+10.2e %7.3f (%.2f-%.2f) %8.3f %6.2f\n', bandsUsed(ib, :), ...
                mean(squeeze(est5.detected(c, ib, :, m))), nanMedian(squeeze(est5.w(c, ib, :, m)) / cs.wn(m) - 1), ...
                nanMedian(zr), zQ(1), zQ(2), nanMedian(squeeze(est5.qErr(c, ib, :, m))), ...
                mean(est5.nSpurious(c, ib, :)));
        end
    end
end
fprintf('(spur = persistent unassigned branches inside the ORIGINAL band [%g %g])\n', BAND);

%% ---- Save ----------------------------------------------------------------------------------
config = struct('DT', DT, 'N_LIST', N_LIST, 'N_REAL', N_REAL, 'ROUTES', {ROUTES}, 'BAND', BAND, ...
    'TAPER_FRACTION', TAPER_FRACTION, 'SU', SU, 'BURN_IN', BURN_IN, 'I_BLOCK', I_BLOCK, 'N_MAX', N_MAX, ...
    'NOISE_FRACTION', NOISE_FRACTION, 'NOISE_LIST', NOISE_LIST, 'N_STAGE4', N4, ...
    'BAND_LIST', {BAND_LIST}, 'bandsUsed', bandsUsed, 'frfModel', 'zoh', ...
    'tol', tol, 'PASS', PASS, 'cases', cases);
save([baseName, '.mat'], 'config', 'est', 'est4', 'est5', 'summ', 'passTab', 'controlOK', 'faithfulOK', ...
    'bandLimitedOK', 'gateOK', 'stabExample');
fprintf('\nSaved: %s.mat\n', baseName);

%% ---- Figures ---------------------------------------------------------------------------------
colours = {[42 120 214]/255, [235 104 52]/255, [27 175 122]/255, [237 161 0]/255, [232 123 164]/255, ...
    [0 131 0]/255};                                                % validated categorical slots 1-6
ink = [0.35 0.35 0.33]; gridInk = [0.85 0.85 0.83]; refInk = [0.55 0.55 0.52];
markers = {'o', 's', '^', 'd', 'v', 'p'};

% (a) damping ratio vs N_BW, one panel per (case, mode)
panels = [1 1; 2 1; 2 2];
fig1 = figure('Color', 'w', 'Position', [100 100 1150 360]);
for p = 1:size(panels, 1)
    c = panels(p, 1); m = panels(p, 2);
    ax = subplot(1, size(panels, 1), p); hold(ax, 'on');
    fill(ax, [1e-1 1e4 1e4 1e-1], 1 + PASS.dampBiasMax * [-1 -1 1 1], [0.94 0.94 0.92], 'EdgeColor', 'none');
    plot(ax, [1e-1 1e4], [1 1], ':', 'Color', refInk, 'LineWidth', 1);
    hLines = gobjects(1, nRt);
    for iR = 1:nRt
        s = summ(c, m).route{iR};
        s.zMed(isnan(s.zMed)) = 0.02;                          % undetected: plotted on the floor
        s.zIQR(isnan(s.zIQR)) = 0.02;
        x = s.nBW * (1 + 0.04 * (iR - 3.5));                   % small offset so bars don't overlap
        errorbar(ax, x, s.zMed, s.zMed - s.zIQR(:, 1).', s.zIQR(:, 2).' - s.zMed, ...
            'LineStyle', 'none', 'Color', colours{iR}, 'LineWidth', 1, 'CapSize', 0);
        hLines(iR) = plot(ax, x, s.zMed, '-', 'Color', colours{iR}, 'LineWidth', 2, ...
            'Marker', markers{iR}, 'MarkerSize', 7, 'MarkerFaceColor', colours{iR}, 'MarkerEdgeColor', 'w');
    end
    set(ax, 'XScale', 'log', 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off');
    grid(ax, 'on'); set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
    xlim(ax, [0.5 150]); ylim(ax, [0 2]);
    text(ax, 0.6, 0.08, 'points on the floor = mode not found', 'FontSize', 8, 'Color', ink);
    xlabel(ax, 'N_{BW} = 2\zeta\omega_n / \Delta\omega  (bins across half-power band)');
    if p == 1, ylabel(ax, 'median \zeta_{SSI} / \zeta_{true}  (IQR bars)'); end
    title(ax, sprintf('%s, \\omega_n = %.1f, \\zeta = %.2f', cases(c).name, cases(c).wn(m), cases(c).zeta(m)), ...
        'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);
    if p == size(panels, 1)
        legend(ax, hLines, {'synthesis, Gaussian', 'synthesis, random phase', 'time domain (control)', ...
            'time domain, band-passed', 'synthesis, no noise', 'synthesis, continuous FRF'}, ...
            'Location', 'southeast', 'Box', 'off');
    end
end
saveFigure(fig1, [baseName, '_damping.png']);

% (b) stabilisation diagrams, twoClose, Gaussian synthesis vs control, short/medium/long record
fig2 = figure('Color', 'w', 'Position', [100 520 1150 600]);
exRoutes = {'synthGaussian', 'timeDomain'};
Ns = unique([stabExample.N]);
for a = 1:numel(exRoutes)
    for b = 1:numel(Ns)
        k = find(strcmp({stabExample.route}, exRoutes{a}) & [stabExample.N] == Ns(b), 1);
        ax = subplot(numel(exRoutes), numel(Ns), (a - 1) * numel(Ns) + b); hold(ax, 'on');
        ex = stabExample(k);
        if strcmp(exRoutes{a}, 'synthGaussian') && Ns(b) <= 256
            for wb = ex.binOmega                               % synthesis tones, visible when coarse
                plot(ax, [wb wb], [0 N_MAX + 1], '-', 'Color', [0.93 0.93 0.91], 'LineWidth', 1);
            end
        end
        for wn = cases(2).wn
            plot(ax, [wn wn], [0 N_MAX + 1], ':', 'Color', refInk, 'LineWidth', 1);
        end
        for n = 1:N_MAX
            f = ex.results(n).frequency;
            plot(ax, f, n * ones(size(f)), '.', 'Color', [0.7 0.7 0.68], 'MarkerSize', 8);
        end
        for bb = 1:numel(ex.branches)
            if ex.persistence(bb).isPersistent
                cls3 = ex.branches(bb).class == 3;
                plot(ax, ex.branches(bb).frequency(cls3), ex.branches(bb).order(cls3), 'o', ...
                    'MarkerSize', 5, 'MarkerFaceColor', colours{1 + 2 * (a == 2)}, 'MarkerEdgeColor', 'w');
            end
        end
        xlim(ax, BAND); ylim(ax, [0 N_MAX + 1]);
        set(ax, 'XColor', ink, 'YColor', ink, 'FontSize', 9, 'Box', 'off');
        grid(ax, 'on'); set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
        nBWex = 2 * min(cases(2).zeta .* cases(2).wn) * Ns(b) * DT / (2 * pi);
        title(ax, sprintf('%s, N = %d (N_{BW} = %.1f)', exRoutes{a}, Ns(b), nBWex), ...
            'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);
        if b == 1, ylabel(ax, 'model order'); end
        if a == numel(exRoutes), xlabel(ax, '\omega (rad/s)'); end
    end
end
saveFigure(fig2, [baseName, '_stabilisation.png']);

% (c) noise sweep
fig3 = figure('Color', 'w', 'Position', [100 100 760 360]);
ax = axes(fig3); hold(ax, 'on');
xN = NOISE_LIST; xN(xN == 0) = 3e-4;                               % zero noise plotted at the left edge
fill(ax, [1e-4 1 1 1e-4], 1 + PASS.dampBiasMax * [-1 -1 1 1], [0.94 0.94 0.92], 'EdgeColor', 'none');
plot(ax, [1e-4 1], [1 1], ':', 'Color', refInk, 'LineWidth', 1);
hN = gobjects(1, 3); labs = cell(1, 3); p = 0;
for c = 1:nC
    for m = 1:numel(cases(c).wn)
        p = p + 1;
        zm = arrayfun(@(iz) nanMedian(squeeze(est4.z(c, iz, :, m)) / cases(c).zeta(m)), 1:nNz);
        zm(isnan(zm)) = 0.02;
        hN(p) = plot(ax, xN, zm, '-', 'Color', colours{p}, 'LineWidth', 2, 'Marker', markers{p}, ...
            'MarkerSize', 7, 'MarkerFaceColor', colours{p}, 'MarkerEdgeColor', 'w');
        labs{p} = sprintf('%s, \\omega_n = %.1f', cases(c).name, cases(c).wn(m));
    end
end
set(ax, 'XScale', 'log', 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off');
grid(ax, 'on'); set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
xlim(ax, [2e-4 0.2]); ylim(ax, [0 2]);
xlabel(ax, 'measurement noise std / signal std  (leftmost point = no noise)');
ylabel(ax, 'median \zeta_{SSI} / \zeta_{true}');
title(ax, sprintf('Gaussian synthesis, band [%g %g] rad/s, N = %d', BAND, N4), 'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);
legend(ax, hN, labs, 'Location', 'southeast', 'Box', 'off');
saveFigure(fig3, [baseName, '_noise.png']);

fprintf('Done in %.1f min.\n', toc(tStart) / 60);
warning(warnState);
catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ==================================================================================================
%  Local functions
%% ==================================================================================================
function [y, binOmega] = synthesiseCase(cs, route, N, dt, band, taperFrac, Su, burnIn, noiseFrac, seed)
% Acceleration response of the modal model for one route, plus white
% measurement noise (std = noiseFrac * pooled signal std; none for
% synthNoiseless). binOmega = synthesis frequencies (empty for time domain).
binOmega = [];
switch route
    case {'synthGaussian', 'synthRandomPhase', 'synthNoiseless', 'synthGaussianCont'}
        dW = 2 * pi / (N * dt);
        k = ceil(band(1) / dW - 1e-9):floor(band(2) / dW + 1e-9);
        binOmega = k * dW;
        if strcmp(route, 'synthGaussianCont')
            Ha = modalAccelerationFRF(cs, binOmega);
        else
            Ha = modalAccelerationFRFzoh(cs, binOmega, dt);
        end
        mode = 'gaussian';
        if strcmp(route, 'synthRandomPhase'), mode = 'randomPhase'; end
        y = stochasticSensorSynthesis(binOmega, Ha, dt, 'Seed', seed, 'InputPSD', Su, ...
            'Taper', 'cosine', 'TaperFraction', taperFrac, 'AmplitudeMode', mode);
    case 'timeDomain'
        y = simulateModalTimeDomain(cs, dt, N, Su, burnIn, seed);
    case 'timeBandpassed'
        % Same band and taper as the synthesis routes: the taper is taken
        % from stochasticSensorSynthesis itself so the two cannot drift apart.
        y = simulateModalTimeDomain(cs, dt, N, Su, burnIn, seed);
        dW = 2 * pi / (N * dt);
        k = ceil(band(1) / dW - 1e-9):floor(band(2) / dW + 1e-9);
        [~, ~, tinfo] = stochasticSensorSynthesis(k * dW, ones(1, numel(k)), dt, 'Taper', 'cosine', ...
            'TaperFraction', taperFrac, 'InputCoefficients', zeros(numel(k), 1));
        mask = zeros(1, N);
        mask(k + 1) = tinfo.taper;
        mask(N - k + 1) = tinfo.taper;
        y = real(ifft(fft(y, [], 2) .* mask, [], 2));
    otherwise
        error('level4A:badRoute', 'Unknown route %s', route);
end
if ~strcmp(route, 'synthNoiseless') && noiseFrac > 0
    y = y + noiseFrac * std(y(:), 1) * drawNormal(size(y), seed + 5000000);
end
end

function Ha = modalAccelerationFRF(cs, w)
% H_a(w) = sum_m Phi(:,m) * (-w^2) / (wn_m^2 - w^2 + 2i zeta_m wn_m w), for e^{+iwt}.
Ha = zeros(size(cs.Phi, 1), numel(w));
for m = 1:numel(cs.wn)
    Ha = Ha + cs.Phi(:, m) * (-w.^2 ./ (cs.wn(m)^2 - w.^2 + 2i * cs.zeta(m) * cs.wn(m) * w));
end
end

function [Ad, Bd, Cc, Dc] = modalStateSpaceZOH(cs, dt)
% q_m'' + 2 zeta_m wn_m q_m' + wn_m^2 q_m = u(t);  y = Phi * q''.
% Exact ZOH discretisation (state x = [q; q']).
m = numel(cs.wn);
K = -diag(cs.wn.^2); Cd = -diag(2 * cs.zeta .* cs.wn);
Ac = [zeros(m), eye(m); K, Cd];
Bc = [zeros(m, 1); ones(m, 1)];
Cc = cs.Phi * [K, Cd];
Dc = cs.Phi * ones(m, 1);
Md = expm([Ac, Bc; zeros(1, 2 * m + 1)] * dt);
Ad = Md(1:2 * m, 1:2 * m); Bd = Md(1:2 * m, end);
end

function Hd = modalAccelerationFRFzoh(cs, w, dt)
% Exact FRF of the sampled (ZOH) system: H_d = C (e^{iw dt} I - A_d)^{-1} B_d + D.
% Same poles as the continuous model; numerator differs (ZOH hold, delay).
[Ad, Bd, Cc, Dc] = modalStateSpaceZOH(cs, dt);
nx = size(Ad, 1);
Hd = zeros(size(Cc, 1), numel(w));
for q = 1:numel(w)
    Hd(:, q) = Cc * ((exp(1i * w(q) * dt) * eye(nx) - Ad) \ Bd) + Dc;
end
end

function y = simulateModalTimeDomain(cs, dt, N, Su, burnIn, seed)
% Time-domain simulation of the ZOH system (modalStateSpaceZOH); u is
% sampled white noise with one-sided PSD Su per rad/s over [0, pi/dt]
% -> variance pi*Su/dt.
[Ad, Bd, Cc, Dc] = modalStateSpaceZOH(cs, dt);
m = numel(cs.wn);
u = sqrt(pi * Su / dt) * drawNormal([1, N + burnIn], seed);
x = zeros(2 * m, 1);
y = zeros(size(cs.Phi, 1), N);
for n = 1:N + burnIn
    if n > burnIn
        y(:, n - burnIn) = Cc * x + Dc * u(n);
    end
    x = Ad * x + Bd * u(n);
end
end

function [results, branches, persistence] = runSSIChain(y, dt, iBlock, nMax, tol)
l = size(y, 1);
[Yp, Yf] = buildHankelMatrix(y, iBlock);
[~, ~, Rlq] = hankelProjection(Yp, Yf);
ri = size(Yp, 1);
R21 = Rlq(ri + 1:end, 1:ri);                 % same U, S as P (see header)
results = sweepModelOrders(R21, l, dt, nMax);
matches = cell(nMax, 1);
for n = 2:nMax
    matches{n} = classifyPoleStability(results(n - 1), results(n), tol.freq, tol.damp, tol.mac);
end
branches = buildModalBranches(results, matches);
persistence = computeBranchPersistence(branches, tol.W, tol.K);
end

function e = assignModes(results, branches, persistence, cs, tol)
nM = numel(cs.wn);
e.w = NaN(1, nM); e.z = NaN(1, nM); e.qErr = NaN(1, nM); e.detected = false(1, nM);
e.nSpurious = 0; e.nToneLike = 0;
nB = numel(branches);
if nB == 0, return; end
isP = [persistence.isPersistent];
medF = arrayfun(@(b) median(b.frequency), branches);
medZ = arrayfun(@(b) median(b.damping), branches);
inBand = medF >= tol.band(1) & medF <= tol.band(2);
chosen = zeros(1, nM);
for m = 1:nM
    [dist, nearest] = min(abs(medF(:) - cs.wn(:).') ./ cs.wn(:).', [], 2);   % nearest true mode per branch
    cand = find(isP(:) & nearest == m & dist <= tol.assign);
    if isempty(cand), continue; end
    score = [[persistence(cand).bestWindowCount].', arrayfun(@(b) numel(branches(b).order), cand(:)), -dist(cand)];
    [~, order] = sortrows(score, [-1 -2 -3]);
    b = cand(order(1));
    chosen(m) = b;
    br = branches(b);
    use = br.class == 3;
    if ~any(use), use = true(size(br.order)); end
    e.w(m) = median(br.frequency(use));
    e.z(m) = median(br.damping(use));
    [~, ref] = max(abs(cs.Phi(:, m)));             % frozen per mode from the truth
    qTrue = cs.Phi(:, m) / cs.Phi(ref, m);
    ords = br.order(use); fr = br.frequency(use);
    Q = NaN(numel(qTrue), numel(ords));
    for t = 1:numel(ords)
        idx = find(results(ords(t)).frequency == fr(t), 1);
        phi = results(ords(t)).modeShapes(:, idx);
        Q(:, t) = phi / phi(ref);
    end
    qHat = nanMedianRows(real(Q)) + 1i * nanMedianRows(imag(Q));
    e.qErr(m) = max(abs(qHat - qTrue));
    e.detected(m) = true;
end
others = setdiff(find(isP & inBand), chosen(chosen > 0));
e.nSpurious = numel(others);
e.nToneLike = sum(isP & inBand & medZ < tol.undamped * min(cs.zeta));
end

function g = drawNormal(sz, seed)
if exist('OCTAVE_VERSION', 'builtin') > 0
    saved = randn('state'); randn('state', seed); %#ok<RAND>
    g = randn(sz); randn('state', saved);         %#ok<RAND>
else
    g = randn(RandStream('mt19937ar', 'Seed', seed), sz);
end
end

function v = nanMedian(x)
x = x(~isnan(x));
if isempty(x), v = NaN; else, v = median(x); end
end

function v = nanMedianRows(X)
v = NaN(size(X, 1), 1);
for j = 1:size(X, 1), v(j) = nanMedian(X(j, :)); end
end

function q = nanQuantiles(x, p)
% Toolbox-free quantiles (linear interpolation, as prctile's default).
x = sort(x(~isnan(x(:))));
if isempty(x), q = NaN(size(p)); return; end
n = numel(x);
pos = 1 + (n - 1) * p(:);
lo = floor(pos); hi = ceil(pos);
q = x(lo) + (pos - lo) .* (x(hi) - x(lo));
q = reshape(q, size(p));
end

function s = yesNo(tf)
if tf, s = 'yes'; else, s = 'no'; end
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);     % R2020a+
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end