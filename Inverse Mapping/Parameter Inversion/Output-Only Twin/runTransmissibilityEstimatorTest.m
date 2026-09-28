%% runTransmissibilityEstimatorTest.m
% S4 / Gate G2: finite-record statistics of the output-only transmissibility
% estimator. Many S3 twin records (synthesiseTwinRecords) go through the S2
% estimator (extractTransmissibility); per Welch bin we measure the bias
% and the 6 x 6 covariance of
%     z_k = [Re T_1; Re T_2; Re T_3; Im T_1; Im T_2; Im T_3](omega_k),
% i.e. exactly the frequency-major blocks of the S1 feature vector, as a
% function of record length, input colour and sensor noise.
%
% CASES (per record length; common random numbers: the same wave
% realisation xi for white and JONSWAP, and the same noise draw)
%   1 W-clean     white input, no sensor noise
%   2 J-clean     JONSWAP input, no sensor noise
%   3 W-noise     white input, sensor noise, no correction
%   4 W-noise-c   same records as 3, reference-noise correction (NoisePSD)
%   5 J-noise     JONSWAP input, the SAME absolute sensor noise
%   6 J-noise-c   same records as 5, corrected
% White and JONSWAP are matched in IN-BAND input variance (same HsInBand),
% so the comparison isolates spectral colour. The absolute noise std is set
% once (NOISE_REL x expected JONSWAP signal std of the reference channel)
% and used for every noisy case: fixed physical sensor noise.
%
% TRUTH AND BIAS MODEL
%   zTrue  interpolated twin FRF at the Welch bins (what the records follow)
%   zS1    S1 EMM features at the Welch bins (CHECK_S1); zS1 - zTrue is the
%          interpolation error, printed once
%   zPred  the Welch-EXPECTED transmissibility: ratio of the expected cross-
%          and auto-spectra of a Hann-windowed segment,
%            T_pred,j(w_k) = sum_q K(w_k - w_q) G_q H_j H_r^*
%                            / (sum_q K(w_k - w_q) G_q |H_r|^2 + S_nn),
%          K = |W_Hann|^2 / (N U) (normalised spectral window), G_q the
%          tapered input PSD on the record's fine grid, S_nn only for the
%          uncorrected noisy cases. It contains the two deterministic biases:
%          RESOLUTION bias (H_j/H_r varies across the window main lobe, and
%          the input colour weights that variation) and REFERENCE-NOISE bias.
%          What remains (mean - zPred) is the ratio-estimator bias, O(1/nEff).
%
% STATISTICS per case, record length and bin (valid realisations only)
%   bias b_k = mean(z) - zTrue, relative |b_k| / |zTrue_k|
%   residual bias mean(z) - zPred (what the bias model does not explain)
%   Sigma_k  sample covariance (6 x 6); relative std sqrt(tr Sigma_k)/|zTrue_k|
%   Hotelling T2_k = n b' Sigma^-1 b, fraction of bins above 16.81 (the
%            ASYMPTOTIC chi2(6) 99% point, not the exact finite-sample
%            Hotelling threshold; bins may have as few as 20 valid draws):
%            a diagnostic significance check, against zTrue and against zPred
%   validFrac, mean coherence, Bendat-Piersol ratio (analytical varTdiag /
%   Monte Carlo E|dT_j|^2), inter-frequency correlation mean|corr(z_k, z_{k+l})|,
%   l = 1..3 (the U1 blkdiag(Sigma_k) assumption), with its Monte Carlo
%   floor sqrt(2 / (pi N_MC)) for truly uncorrelated bins
%
% GATE G2 CHECKS (PASS / FLAG; FLAG = read the numbers, not a crash)
%   G2.1  bias model: against zPred at most 10% of bins show a significant
%         residual bias (diagnostic Hotelling check) in every case at the longest record,
%         i.e. the deterministic bias is explained to Monte Carlo resolution
%         (against zTrue the fraction is typically near 1)
%   G2.2  covariance ~ 1/nEff: log-log slope of median(tr Sigma/|z|^2) vs
%         nEff within [-1.15, -0.85], every case
%   G2.3  correction: in bins with predicted noise bias > 1%, the NOISE part
%         of the bias (distance of the mean from the noise-free prediction)
%         is reduced at least 3-fold by the correction (W and J). The
%         resolution part is untouched by design, so the total is not used.
%   G2.4  Bendat-Piersol diagnostic within a factor 2 of Monte Carlo (noisy)
%   G2.5  colour: in noise-dominated bins (noisy variance > 3 x clean
%         variance for W and J) the slope of log(var_J/var_W) vs
%         log(G_rJ/G_rW) is near -1 (variance ~ S_nn/G_r)
%   G2.6  inter-frequency correlation at lag 1 vs Monte Carlo floor [INFO]
%   G2.7  segment-length trade-off at the longest record, SEG_SWEEP [INFO]
%   G2.8  resolution-bias correction without knowing S_u: residual if the
%         model side is Welch-smoothed with WHITE weighting while the sea is
%         JONSWAP, |zPred_J - zPred_W| / |z| [INFO, for U2]
%
% RUNTIME (MATLAB): twin 57 EMM solves (about 1 min, cached in
% Results/twinCache, instant afterwards), CHECK_S1 45 solves (about 45 s),
% Monte Carlo 3 lengths x N_MC x 4 records (a few minutes).
%
% OUTPUT (in Inverse Mapping/Parameter Inversion/Results/)
%   transmissibilityEstimatorTest_<stamp>.mat / .txt / _stats.png / _snr.png
%   out.stats{iN, c}.Sigma (6 x 6 x nF), .bias, .zPred, ... are the inputs for U1;
%   out.Z{iN, c} (6 x nF x N_MC raw draws) allows a full or banded
%   cross-frequency covariance; out.Zsw / omegaSw / nEffSw keep the L sweep.
%
% Lives in Inverse Mapping/Parameter Inversion/SSI Twin/.

thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..', '..', '..');
addpath(fullfile(thisDir, '..', '..'));                        % Inverse Mapping: transmissibilityFeatureVector, computeSensorFRF
addpath(fullfile(rootDir, 'Transmissibility'));               % S2, S3, stack helpers
addpath(fullfile(rootDir, 'SSI'));
addpath(fullfile(rootDir, 'Forward Model'));
addpath(fullfile(rootDir, 'Forward Model', 'Animation'));
clearvars -except thisDir rootDir
close all, clc

%% ---- Configuration -----------------------------------------------------------
P0 = [4.6985e-5, 1.4548e-3, 0.3830]; R0 = P0(3);
SENSORS = [0.3*R0, 0; 0.3*R0, pi; 0.5*R0, pi; 0.9*R0, 0];    % Level 2C: s2, s6, s14, s26
REF = 4;                                                     % s26
DT = 0.1; SEG = 512; OVERLAP = 0.5; BAND = [3 8.5];          % 45 Welch bins
N_LIST = [4096 16384 65536];                                 % 410 s, 27 min, 109 min records (multiples of SEG)
N_MC = 200;
HS = 0.05; WP = 5.0; GAMMA_J = 3.3;                          % JONSWAP (basin scale)
NOISE_REL = 0.10;                                            % noise std / expected JONSWAP s26 signal std: a controlled
                                                             % 10% measurement-noise STRESS TEST, not an IMU specification
SEG_SWEEP = [256 512 1024];                                  % G2.7, longest record only
TWIN_SPACING = 0.1;
CHECK_S1 = true;
SEED0 = 4000;

%% ---- Output ----------------------------------------------------------------------
resultsDir = fullfile(thisDir, '..', 'Results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
cacheDir = fullfile(resultsDir, 'twinCache');
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['transmissibilityEstimatorTest_', stamp]);
diary([baseName, '.txt']);
try

caseNames = {'W-clean', 'J-clean', 'W-noise', 'W-noise-c', 'J-noise', 'J-noise-c'};
caseInput = [1 2 1 1 2 2];                                   % 1 white, 2 JONSWAP
caseSnnDen = [0 0 1 0 1 0];                                  % S_nn in the expected denominator
nC = numel(caseNames); nL = numel(N_LIST); m = size(SENSORS, 1) - 1; nS = m + 1;
others = setdiff(1:nS, REF);
fprintf('runTransmissibilityEstimatorTest  (%s)\n', stamp);
fprintf('p0 = [%.4e %.4e %.4f], reference s26, dt %.2f s, L %d, overlap %.2f, band [%.2f %.2f] rad/s\n', ...
    P0, DT, SEG, OVERLAP, BAND);
fprintf('record lengths %s samples, N_MC = %d\n', mat2str(N_LIST), N_MC);
assert(all(mod(N_LIST, max(SEG_SWEEP)) == 0), 'N_LIST entries must be multiples of every segment length');
selfCheckHannKernel();
fprintf('Hann spectral-window kernel: closed form and normalisation verified\n\n');

%% ---- 0  Twin (cached node FRF) ----------------------------------------------------------
t0 = tic;
twin = synthesiseTwinRecords(struct('p', P0, 'sensors', SENSORS, 'nodeSpacing', TWIN_SPACING), ...
    'CacheDir', cacheDir, 'Verbose', true);
fprintf('twin: %s, %d nodes in [%.4f %.4f] rad/s, interpolation check H %.1e, T %.1e (%.1f s)\n', ...
    twin.cacheStatus, numel(twin.omegaNodes), twin.band, twin.interpCheck.maxRelErrH, ...
    twin.interpCheck.maxRelErrT, toc(t0));

%% ---- 1  Calibration: matched in-band variance, fixed absolute noise --------------------
specJ = {'Spectrum', 'jonswap', 'Hs', HS, 'PeakFrequency', WP, 'PeakEnhancement', GAMMA_J};
[~, ~, iJ] = synthesiseTwinRecords(twin, N_LIST(end), DT, 'Seed', 1, specJ{:});
whiteLevel = (iJ.HsInBand / 4)^2 / (sum(iJ.taper.^2) * iJ.dOmega);
specW = {'Spectrum', 'white', 'WhiteLevel', whiteLevel};
[~, ~, iW] = synthesiseTwinRecords(twin, N_LIST(end), DT, 'Seed', 1, specW{:});
noiseStd = NOISE_REL * iJ.signalStd(REF);
Snn = noiseStd^2 * DT / pi;
specs = {specW, specJ};
fprintf('JONSWAP Hs %.3f m (in band %.4f m), wp %.2f rad/s; white level %.3e m^2 s/rad (in band %.4f m)\n', ...
    HS, iJ.HsInBand, WP, whiteLevel, iW.HsInBand);
fprintf('expected signal std (m/s^2)  JONSWAP %s   white %s\n', mat2str(iJ.signalStd, 3), mat2str(iW.signalStd, 3));
fprintf('sensor noise std %.4g m/s^2 on every channel (%.0f%% of JONSWAP s26 std), S_nn = %.3e\n\n', ...
    noiseStd, 100 * NOISE_REL, Snn);

%% ---- 2  Welch grid, truth, S1 check --------------------------------------------------------------
estOpt = {'Reference', REF, 'SegmentLength', SEG, 'Overlap', OVERLAP, 'Band', BAND};
est0 = extractTransmissibility(randn(nS, SEG), DT, estOpt{:});
omega = est0.omega; nF = numel(omega);
zTrue = trueZ(twin, omega, REF, others);
zNorm = sqrt(sum(zTrue.^2, 1));
GrW = interp1(iW.omega, iW.expectedPSD(REF, :), omega);    % expected s26 signal PSD at the Welch bins
GrJ = interp1(iJ.omega, iJ.expectedPSD(REF, :), omega);
fprintf('%d Welch bins in [%.3f %.3f] rad/s\n', nF, omega([1 end]));
if CHECK_S1
    t1 = tic;
    fS1 = transmissibilityFeatureVector(log(P0([1 3])), omega, SENSORS, 'PBase', P0, 'Reference', REF);
    zS1 = reshape(fS1, 2 * m, nF);
    dS1 = sqrt(sum((zS1 - zTrue).^2, 1)) ./ zNorm;
    fprintf('S1 EMM features vs interpolated truth: max relative difference %.2e (median %.2e)  (%.0f s)\n', ...
        max(dS1), median(dS1), toc(t1));
else
    zS1 = []; dS1 = [];
end
fprintf('reference SNR G_r/S_nn over the band: white %.1f to %.1f dB, JONSWAP %.1f to %.1f dB\n\n', ...
    10 * log10(min(GrW / Snn)), 10 * log10(max(GrW / Snn)), 10 * log10(min(GrJ / Snn)), 10 * log10(max(GrJ / Snn)));

%% ---- 3  Expected (Welch-smoothed) transmissibility per case and record length -------------------------
zPred = cell(nL, nC);
for iN = 1:nL
    infos = cell(1, 2);
    for s = 1:2
        [~, ~, infos{s}] = synthesiseTwinRecords(twin, N_LIST(iN), DT, 'Seed', 1, specs{s}{:});
    end
    for c = 1:nC
        zPred{iN, c} = welchExpectedZ(omega, infos{caseInput(c)}, SEG, DT, REF, others, caseSnnDen(c) * Snn);
    end
end

%% ---- 4  Monte Carlo ----------------------------------------------------------------------------------
Z = cell(nL, nC); VT = cell(nL, nC); COH = cell(nL, nC); nEff = zeros(1, nL); recSec = zeros(1, nL);
nSw = numel(SEG_SWEEP); swCases = [1 2 5];                    % W-clean, J-clean, J-noise
Zsw = cell(nSw, numel(swCases)); omegaSw = cell(1, nSw); nEffSw = zeros(1, nSw);
for iN = 1:nL
    N = N_LIST(iN);
    for c = 1:nC
        Z{iN, c} = NaN(2 * m, nF, N_MC); VT{iN, c} = zeros(m, nF); COH{iN, c} = zeros(m, nF);
    end
    tN = tic;
    for r = 1:N_MC
        seed = SEED0 + 1e5 * iN + r;
        Ywc = synthesiseTwinRecords(twin, N, DT, 'Seed', seed, specW{:});
        Yjc = synthesiseTwinRecords(twin, N, DT, 'Seed', seed, specJ{:});
        Ywn = synthesiseTwinRecords(twin, N, DT, 'Seed', seed, specW{:}, 'NoiseStd', noiseStd);
        Yjn = synthesiseTwinRecords(twin, N, DT, 'Seed', seed, specJ{:}, 'NoiseStd', noiseStd);
        e = cell(1, nC);
        e{1} = extractTransmissibility(Ywc, DT, estOpt{:});
        e{2} = extractTransmissibility(Yjc, DT, estOpt{:});
        e{3} = extractTransmissibility(Ywn, DT, estOpt{:});
        e{4} = extractTransmissibility(Ywn, DT, estOpt{:}, 'NoisePSD', Snn);
        e{5} = extractTransmissibility(Yjn, DT, estOpt{:});
        e{6} = extractTransmissibility(Yjn, DT, estOpt{:}, 'NoisePSD', Snn);
        for c = 1:nC
            Z{iN, c}(:, :, r) = reshape(stackTransmissibility(e{c}.T), 2 * m, nF);
            VT{iN, c} = VT{iN, c} + e{c}.varTdiag / N_MC;
            COH{iN, c} = COH{iN, c} + e{c}.coherence / N_MC;
        end
        if iN == nL                                           % G2.7 segment-length sweep
            recs = {Ywc, Yjc, Yjn};
            for iL = 1:nSw
                for i = 1:numel(swCases)
                    es = extractTransmissibility(recs{i}, DT, 'Reference', REF, 'SegmentLength', SEG_SWEEP(iL), ...
                        'Overlap', OVERLAP, 'Band', BAND);
                    if r == 1 && i == 1
                        omegaSw{iL} = es.omega; nEffSw(iL) = es.nEffective;
                        Zsw{iL, 1} = NaN(2 * m, numel(es.omega), N_MC); Zsw{iL, 2} = Zsw{iL, 1}; Zsw{iL, 3} = Zsw{iL, 1};
                    end
                    Zsw{iL, i}(:, :, r) = reshape(stackTransmissibility(es.T), 2 * m, numel(es.omega));
                end
            end
        end
    end
    nEff(iN) = e{1}.nEffective;
    recSec(iN) = toc(tN);
    fprintf('N = %6d (%.0f s record, %d segments, nEff %.1f): %d x 4 records in %.1f s\n', ...
        N, N * DT, e{1}.nSegments, nEff(iN), N_MC, recSec(iN));
end
fprintf('\n');

%% ---- 5  Statistics -------------------------------------------------------------------------------------
stats = cell(nL, nC);
for iN = 1:nL
    for c = 1:nC
        s = binStats(Z{iN, c}, zTrue, zPred{iN, c}, m);
        s.bpRatio = VT{iN, c} ./ s.varTj;                                   % analytical / Monte Carlo
        s.coherence = COH{iN, c};
        s.zPred = zPred{iN, c};
        stats{iN, c} = s;
    end
end
predNoiseBias = Snn ./ ([GrW; GrJ] + Snn);                                 % 2 x nF, reference-noise part
lagFloor = sqrt(2 / (pi * N_MC));

%% ---- 6  Summary tables ---------------------------------------------------------------------------------
med = @(x) median(x(isfinite(x)));
fprintf('SUMMARY (median over bins; rel = relative to |zTrue_k|)\n');
fprintf('%-10s %7s %10s %10s %10s %10s %9s %8s %8s %8s %8s %8s\n', 'case', 'N', 'rel bias', 'resid', 'rel std', ...
    'std*sqrtN', 'bias/std', 'sig.b', 'sig.res', 'valid', 'BP/MC', 'lag1');
for c = 1:nC
    for iN = 1:nL
        s = stats{iN, c};
        fprintf('%-10s %7d %10.2e %10.2e %10.2e %10.3f %9.2f %8.2f %8.2f %8.3f %8.2f %8.2f\n', caseNames{c}, N_LIST(iN), ...
            med(s.relBias), med(s.relResid), med(s.relStd), med(s.relStd) * sqrt(nEff(iN)), med(s.relBias ./ s.relStd), ...
            mean(s.T2 > 16.81), mean(s.T2res > 16.81), mean(s.validFrac), med(s.bpRatio), med(s.lagCorr(1, :)));
    end
end
fprintf(['(resid = |mean - zPred|/|z|: bias NOT explained by resolution + noise; std*sqrtN constant if var ~ 1/nEff;\n' ...
    ' bias/std = median over bins of |bias|/std for ONE record; sig.b / sig.res = fraction of bins with\n' ...
    ' Hotelling T2 > 16.81 against zTrue / zPred; lag1 floor for uncorrelated bins = %.2f)\n\n'], lagFloor);

iN = nL;
fprintf('PER-BIN DETAIL, N = %d, JONSWAP (every 4th bin)\n', N_LIST(iN));
fprintf('%8s %8s %10s %10s %10s %10s %10s %8s %9s %8s %7s\n', 'omega', 'SNR dB', 'bias clean', 'pred clean', ...
    'bias unc', 'pred unc', 'bias corr', 'valid-c', 'std noisy', 'BP/MC', 'coh');
pr = @(c, k) norm(stats{iN, c}.zPred(:, k) - zTrue(:, k)) / zNorm(k);
for k = 1:4:nF
    fprintf('%8.3f %8.1f %10.2e %10.2e %10.2e %10.2e %10.2e %8.2f %9.2e %8.2f %7.3f\n', omega(k), ...
        10 * log10(GrJ(k) / Snn), stats{iN, 2}.relBias(k), pr(2, k), stats{iN, 5}.relBias(k), pr(5, k), ...
        stats{iN, 6}.relBias(k), stats{iN, 6}.validFrac(k), stats{iN, 5}.relStd(k), ...
        med(stats{iN, 5}.bpRatio(:, k)), mean(stats{iN, 5}.coherence(:, k)));
end
fprintf('\n');

%% ---- 7  Segment-length trade-off (G2.7) ---------------------------------------------------------------------
[~, ~, iWL] = synthesiseTwinRecords(twin, N_LIST(nL), DT, 'Seed', 1, specW{:});
[~, ~, iJL] = synthesiseTwinRecords(twin, N_LIST(nL), DT, 'Seed', 1, specJ{:});
infoSw = {iWL, iJL, iJL}; snnSw = [0 0 Snn];
sweep = struct('L', num2cell(SEG_SWEEP));
fprintf('SEGMENT-LENGTH TRADE-OFF, N = %d (median over bins)\n', N_LIST(nL));
fprintf('%6s %5s %6s   %-26s %-26s %-26s\n', 'L', 'bins', 'nEff', 'W-clean bias / pred / std', ...
    'J-clean bias / pred / std', 'J-noise bias / pred / std');
for iL = 1:nSw
    w = omegaSw{iL};
    zt = trueZ(twin, w, REF, others);
    row = '';
    for i = 1:numel(swCases)
        zp = welchExpectedZ(w, infoSw{i}, SEG_SWEEP(iL), DT, REF, others, snnSw(i));
        si = binStats(Zsw{iL, i}, zt, zp, m);
        pb = sqrt(sum((zp - zt).^2, 1)) ./ sqrt(sum(zt.^2, 1));
        sweep(iL).bias(i) = med(si.relBias); sweep(iL).pred(i) = med(pb); sweep(iL).std(i) = med(si.relStd);
        row = [row, sprintf('%8.2e %8.2e %8.2e  ', sweep(iL).bias(i), sweep(iL).pred(i), sweep(iL).std(i))]; %#ok<AGROW>
    end
    fprintf('%6d %5d %6.1f   %s\n', SEG_SWEEP(iL), numel(w), nEffSw(iL), row);
    sweep(iL).omega = w; sweep(iL).nEff = nEffSw(iL);
end
fprintf('(std is for one record; bias does not shrink with record length, std does)\n\n');

%% ---- 8  Gate G2 checks ------------------------------------------------------------------------------------
checks = struct();
fprintf('GATE G2 CHECKS\n');
% G2.1 bias model
sigRes = zeros(1, nC); sigTrue = zeros(1, nC);
for c = 1:nC
    sigRes(c) = mean(stats{nL, c}.T2res > 16.81); sigTrue(c) = mean(stats{nL, c}.T2 > 16.81);
end
checks.G21sigResid = sigRes; checks.G21sigTrue = sigTrue;
checks.G21 = all(sigRes <= 0.1);
fprintf('  G2.1 fraction of bins with significant bias, vs zTrue -> vs zPred (N = %d): %s   %s\n', N_LIST(nL), ...
    strjoin(arrayfun(@(c) sprintf('%s %.2f->%.2f', caseNames{c}, sigTrue(c), sigRes(c)), 1:nC, 'UniformOutput', false), ', '), ...
    flagStr(checks.G21));
% G2.2 1/nEff
slopes = NaN(1, nC);
for c = 1:nC
    v = zeros(1, nL);
    for iN = 1:nL, v(iN) = med(stats{iN, c}.relStd.^2); end
    q = polyfit(log(nEff), log(v), 1); slopes(c) = q(1);
end
checks.G22slopes = slopes;
checks.G22 = all(slopes > -1.15 & slopes < -0.85);
fprintf('  G2.2 slope of log var vs log nEff: %s   %s\n', ...
    strjoin(arrayfun(@(c) sprintf('%s %.2f', caseNames{c}, slopes(c)), 1:nC, 'UniformOutput', false), ', '), flagStr(checks.G22));
% G2.3 correction
mW = predNoiseBias(1, :) > 1e-2; mJ = predNoiseBias(2, :) > 1e-2;
% noise part of the bias: mean(z) - noise-free prediction (zPred of the corrected case, S_nn not in it)
nb = @(c, cFree, mask) med(sqrt(sum((stats{nL, c}.bias(:, mask) + zTrue(:, mask) - zPred{nL, cFree}(:, mask)).^2, 1)) ./ zNorm(mask));
cW = [nb(3, 4, mW), nb(4, 4, mW)]; cJ = [nb(5, 6, mJ), nb(6, 6, mJ)];
okW = ~any(mW) || cW(2) < cW(1) / 3; okJ = ~any(mJ) || cJ(2) < cJ(1) / 3;
checks.G23 = okW && okJ; checks.G23bins = [sum(mW) sum(mJ)]; checks.G23noiseBias = [cW; cJ];
fprintf('  G2.3 bins with predicted noise bias > 1%%: white %d, JONSWAP %d; median NOISE part of bias unc -> corr: W %.2e -> %.2e, J %.2e -> %.2e   %s\n', ...
    sum(mW), sum(mJ), cW, cJ, flagStr(checks.G23));
% G2.4 Bendat-Piersol
bp = [med(stats{nL, 3}.bpRatio), med(stats{nL, 5}.bpRatio)];
checks.G24 = all(bp > 0.5 & bp < 2); checks.G24ratio = bp;
fprintf('  G2.4 Bendat-Piersol / Monte Carlo variance (noisy, uncorrected): white %.2f, JONSWAP %.2f   %s\n', bp, flagStr(checks.G24));
% G2.5 colour, noise-dominated bins only
vWc = stats{nL, 1}.varTj; vJc = stats{nL, 2}.varTj; vW = stats{nL, 3}.varTj; vJ = stats{nL, 5}.varTj;
x = repmat(log(GrJ ./ GrW), m, 1); y = log(vJ ./ vW);
nd = vW > 3 * vWc & vJ > 3 * vJc & isfinite(x) & isfinite(y);
if sum(nd(:)) >= 5 && numel(unique(round(x(nd) * 1e6))) >= 3
    q = polyfit(x(nd), y(nd), 1); checks.G25slope = q(1);
    checks.G25 = q(1) > -1.25 && q(1) < -0.75;
    fprintf('  G2.5 colour, %d noise-dominated (bin, channel) pairs: slope of log(var_J/var_W) vs log(G_rJ/G_rW) = %.2f (theory -1)   %s\n', ...
        sum(nd(:)), q(1), flagStr(checks.G25));
else
    checks.G25slope = NaN; checks.G25 = NaN;
    fprintf('  G2.5 colour: only %d noise-dominated pairs, slope not fitted   INFO\n', sum(nd(:)));
end
% G2.6 inter-frequency correlation
lag = zeros(3, 4); lagCase = [1 2 3 5];
for l = 1:3
    for i = 1:4, lag(l, i) = med(stats{nL, lagCase(i)}.lagCorr(l, :)); end
end
checks.G26lag = lag; checks.G26floor = lagFloor;
fprintf('  G2.6 median mean|corr(z_k, z_{k+l})| (W-clean, J-clean, W-noise, J-noise), floor %.2f:\n', lagFloor);
for l = 1:3, fprintf('         lag %d: %.2f %.2f %.2f %.2f\n', l, lag(l, :)); end
fprintf('       INFO: lag-1 well above the floor means blkdiag(Sigma_k) overstates information in U1\n');
fprintf('       (use every other bin, or a banded covariance).\n');
fprintf('  G2.7 segment-length trade-off: see table above   INFO\n');
% G2.8 white-weighted smoothing when the sea is JONSWAP
mis = sqrt(sum((zPred{nL, 2} - zPred{nL, 1}).^2, 1)) ./ zNorm;
checks.G28mismatch = mis;
fprintf(['  G2.8 resolution bias left if the model is Welch-smoothed with WHITE weighting but the sea is JONSWAP:\n' ...
    '       median %.2e, max %.2e (vs raw JONSWAP resolution bias median %.2e)   INFO\n\n'], med(mis), max(mis), ...
    med(sqrt(sum((zPred{nL, 2} - zTrue).^2, 1)) ./ zNorm));

%% ---- 9  Save --------------------------------------------------------------------------------------------------
out = struct('p0', P0, 'sensors', SENSORS, 'reference', REF, 'dt', DT, 'segmentLength', SEG, ...
    'overlap', OVERLAP, 'band', BAND, 'NList', N_LIST, 'nEff', nEff, 'nMC', N_MC, 'omega', omega, ...
    'zTrue', zTrue, 'zS1', zS1, 'dS1', dS1, 'caseNames', {caseNames}, 'stats', {stats}, ...
    'whiteLevel', whiteLevel, 'jonswap', struct('Hs', HS, 'wp', WP, 'gammaJ', GAMMA_J, 'HsInBand', iJ.HsInBand), ...
    'noiseStd', noiseStd, 'noiseRel', NOISE_REL, 'Snn', Snn, 'GrW', GrW, 'GrJ', GrJ, ...
    'predNoiseBias', predNoiseBias, 'sweep', sweep, 'checks', checks, ...
    'twinInterpCheck', twin.interpCheck, 'twinCacheFile', twin.cacheFile);
% raw Monte Carlo draws (about 8 MB): Z{iN, c} is 6 x nF x N_MC, NaN where the
% corrected estimate is invalid. Kept so U1 can build the full or banded
% cross-frequency covariance without rerunning S4.
out.Z = Z;
out.Zsw = Zsw; out.omegaSw = omegaSw; out.nEffSw = nEffSw; out.segSweep = SEG_SWEEP; out.sweepCases = caseNames(swCases);
save([baseName, '.mat'], 'out');
fprintf('Saved: %s.mat\n', baseName);

%% ---- 10  Figures ----------------------------------------------------------------------------------------------
colW = [42 120 214] / 255; colJ = [235 104 52] / 255; ink = [0.35 0.35 0.33]; gridInk = [0.85 0.85 0.83];
cc = {colW, colJ, colW, colW, colJ, colJ};
ls = {':', ':', '-', '--', '-', '--'};
fig = figure('Color', 'w', 'Position', [60 60 1200 820]);
ax = subplot(2, 2, 1); hold(ax, 'on'); h = []; lab = {};
mk = {'o', 'o', '', '', 's', '^'};
for c = [1 2 5 6]
    h(end + 1) = plot(ax, omega, stats{nL, c}.relBias, mk{c}, 'Color', cc{c}, 'MarkerSize', 4, 'MarkerFaceColor', cc{c}); %#ok<AGROW>
    lab{end + 1} = [caseNames{c}, ' measured']; %#ok<AGROW>
    pb = sqrt(sum((stats{nL, c}.zPred - zTrue).^2, 1)) ./ zNorm;
    h(end + 1) = plot(ax, omega, pb, ls{c}, 'Color', cc{c}, 'LineWidth', 1.6); %#ok<AGROW>
    lab{end + 1} = [caseNames{c}, ' predicted']; %#ok<AGROW>
end
styleAx(ax, ink, gridInk, 'relative bias |b_k| / |z_k|', sprintf('(a) bias vs Welch-expected prediction, N = %d', N_LIST(nL)));
legend(ax, h, lab, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
ax = subplot(2, 2, 2); hold(ax, 'on');
for c = 1:nC
    plot(ax, omega, stats{nL, c}.relStd, ls{c}, 'Color', cc{c}, 'LineWidth', 1.8);
end
styleAx(ax, ink, gridInk, 'relative std  sqrt(tr \Sigma_k) / |z_k|', sprintf('(b) spread of one record, N = %d', N_LIST(nL)));
legend(ax, caseNames, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
ax = subplot(2, 2, 3); hold(ax, 'on'); h = []; lab = {};
shade = linspace(0.35, 1, nL);
for iN = 1:nL
    h(end + 1) = plot(ax, omega, stats{iN, 3}.relStd * sqrt(nEff(iN)), '-', 'Color', 1 - shade(iN) * (1 - colW), 'LineWidth', 1.8); %#ok<AGROW>
    lab{end + 1} = sprintf('W-noise N = %d', N_LIST(iN)); %#ok<AGROW>
    h(end + 1) = plot(ax, omega, stats{iN, 5}.relStd * sqrt(nEff(iN)), '-', 'Color', 1 - shade(iN) * (1 - colJ), 'LineWidth', 1.8); %#ok<AGROW>
    lab{end + 1} = sprintf('J-noise N = %d', N_LIST(iN)); %#ok<AGROW>
end
styleAx(ax, ink, gridInk, 'relative std \times sqrt(n_{eff})', '(c) 1/n_{eff} scaling: curves collapse if it holds');
legend(ax, h, lab, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
ax = subplot(2, 2, 4); hold(ax, 'on'); h = [];
for i = 1:numel(lagCase)
    c = lagCase(i);
    h(end + 1) = plot(ax, omega(1:end - 1), stats{nL, c}.lagCorr(1, :), ls{c}, 'Color', cc{c}, 'LineWidth', 1.8); %#ok<AGROW>
end
plot(ax, omega([1 end]), lagFloor * [1 1], '-', 'Color', ink, 'LineWidth', 1);
ylim(ax, [0 1]);
styleAx(ax, ink, gridInk, 'mean |corr(z_k, z_{k+1})|', '(d) adjacent-bin correlation (grey: Monte Carlo floor)', true);
legend(ax, h, caseNames(lagCase), 'Location', 'best', 'Box', 'off', 'FontSize', 8);
saveFigure(fig, [baseName, '_stats.png']);

fig2 = figure('Color', 'w', 'Position', [80 80 1500 420]);
ax = subplot(1, 3, 1); hold(ax, 'on');
plot(ax, omega, 10 * log10(GrW / Snn), '-', 'Color', colW, 'LineWidth', 1.8);
plot(ax, omega, 10 * log10(GrJ / Snn), '-', 'Color', colJ, 'LineWidth', 1.8);
styleAx(ax, ink, gridInk, 'reference SNR  G_r / S_{nn}  (dB)', '(a) reference-channel SNR', true);
legend(ax, {'white', 'JONSWAP'}, 'Location', 'best', 'Box', 'off');
ax = subplot(1, 3, 2); hold(ax, 'on');
plot(ax, omega, stats{nL, 6}.validFrac, '-', 'Color', colJ, 'LineWidth', 1.8);
plot(ax, omega, stats{nL, 4}.validFrac, '--', 'Color', colW, 'LineWidth', 1.8);
ylim(ax, [0 1.05]);
styleAx(ax, ink, gridInk, 'fraction of valid corrected estimates', sprintf('(b) correction validity, N = %d', N_LIST(nL)), true);
legend(ax, {'J-noise-c', 'W-noise-c'}, 'Location', 'best', 'Box', 'off');
ax = subplot(1, 3, 3); hold(ax, 'on'); h = []; lab = {};
swCol = {colW, colJ, colJ}; swLs = {':', ':', '-'}; swNames = caseNames(swCases);
for i = 1:numel(swCases)
    h(end + 1) = plot(ax, SEG_SWEEP, arrayfun(@(sw) sw.bias(i), sweep), [swLs{i}, 'o'], 'Color', swCol{i}, 'LineWidth', 1.6); %#ok<AGROW>
    lab{end + 1} = [swNames{i}, ' bias']; %#ok<AGROW>
    h(end + 1) = plot(ax, SEG_SWEEP, arrayfun(@(sw) sw.std(i), sweep), [swLs{i}, 's'], 'Color', 1 - 0.5 * (1 - swCol{i}), 'LineWidth', 1.6); %#ok<AGROW>
    lab{end + 1} = [swNames{i}, ' std']; %#ok<AGROW>
end
set(ax, 'XScale', 'log');
styleAx(ax, ink, gridInk, 'median relative bias / std', sprintf('(c) segment length trade-off, N = %d', N_LIST(nL)));
xlabel(ax, 'segment length L (samples)');
legend(ax, h, lab, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
saveFigure(fig2, [baseName, '_snr.png']);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ================================================================================================
%  Local functions
%% ================================================================================================
function z = trueZ(twin, w, ref, others)
H = twin.Hacc(w);
z = reshape(stackTransmissibility(H(others, :) ./ H(ref, :)), 2 * numel(others), numel(w));
end

function z = welchExpectedZ(wBins, info, L, dt, ref, others, SnnDen)
% Welch-expected transmissibility (ratio of expected spectra) for records
% synthesised on the fine grid info.omega with tapered input PSD, estimated
% with a periodic-Hann window of L samples. SnnDen adds white reference noise.
wq = info.omega;
G = info.inputPSD .* info.taper.^2;                        % 1 x nq
Hq = info.H;                                               % nS x nq
U = 3 * L / 8;                                             % sum of periodic Hann squared
K = abs(hannDTFT((wBins(:) - wq) * dt, L)).^2 / (info.N * U);   % nF x nq, normalised spectral window
Srr = K * (G .* abs(Hq(ref, :)).^2).' + SnnDen;            % nF x 1
Sjr = K * (G .* Hq(others, :) .* conj(Hq(ref, :))).';       % nF x m
T = (Sjr ./ Srr).';                                        % m x nF
z = reshape(stackTransmissibility(T), 2 * numel(others), numel(wBins));
end

function W = hannDTFT(nu, L)
% DTFT of the periodic Hann window w_n = 0.5 - 0.5 cos(2 pi n / L), n = 0..L-1,
% at normalised frequency nu (rad/sample): W = 0.5 D(nu) - 0.25 D(nu - 2pi/L) - 0.25 D(nu + 2pi/L),
% D(x) = sum_n exp(-i x n) = (1 - exp(-i x L)) / (1 - exp(-i x)).
W = 0.5 * dirichletSum(nu, L) - 0.25 * dirichletSum(nu - 2 * pi / L, L) - 0.25 * dirichletSum(nu + 2 * pi / L, L);
end

function D = dirichletSum(x, L)
den = 1 - exp(-1i * x);
D = (1 - exp(-1i * x * L)) ./ den;
small = abs(den) < 1e-9;
D(small) = L;
end

function selfCheckHannKernel()
% Closed form vs direct sum, and normalisation sum_d |W(2 pi d/N)|^2 = N U.
L = 64; n = 0:L - 1; w = 0.5 - 0.5 * cos(2 * pi * n / L);
nu = [0, 0.013, 2 * pi / L, 0.5, 1.7, -0.9];
Wd = exp(-1i * nu(:) * n) * w(:);
assert(max(abs(hannDTFT(nu(:), L) - Wd)) < 1e-9 * L, 'Hann DTFT closed form wrong');
N = 1024; d = 0:N - 1;
assert(abs(sum(abs(hannDTFT(2 * pi * d / N, L)).^2) / (N * 3 * L / 8) - 1) < 1e-10, 'Hann kernel normalisation wrong');
end

function s = binStats(Zc, zTrue, zPred, m)
% Zc: 2m x nF x nR draws (NaN where invalid). Per-bin statistics over the
% realisations whose 2m components are all finite.
[p, nF, nR] = size(Zc);
s = struct('bias', NaN(p, nF), 'resid', NaN(p, nF), 'Sigma', NaN(p, p, nF), 'relBias', NaN(1, nF), ...
    'relResid', NaN(1, nF), 'relStd', NaN(1, nF), 'T2', NaN(1, nF), 'T2res', NaN(1, nF), ...
    'validFrac', zeros(1, nF), 'nValid', zeros(1, nF), 'varTj', NaN(m, nF), 'lagCorr', NaN(3, nF - 1));
zn = sqrt(sum(zTrue.^2, 1));
for k = 1:nF
    X = reshape(Zc(:, k, :), p, nR);
    ok = all(isfinite(X), 1);
    s.nValid(k) = sum(ok); s.validFrac(k) = mean(ok);
    if sum(ok) < 20, continue; end
    X = X(:, ok);
    mu = mean(X, 2);
    C = cov(X.');
    Ci = pinv(C);
    b = mu - zTrue(:, k); br = mu - zPred(:, k);
    s.bias(:, k) = b; s.resid(:, k) = br; s.Sigma(:, :, k) = C;
    s.relBias(k) = norm(b) / zn(k); s.relResid(k) = norm(br) / zn(k);
    s.relStd(k) = sqrt(trace(C)) / zn(k);
    s.T2(k) = sum(ok) * (b.' * Ci * b); s.T2res(k) = sum(ok) * (br.' * Ci * br);
    s.varTj(:, k) = diag(C(1:m, 1:m)) + diag(C(m + 1:end, m + 1:end));   % E|dT_j|^2 about the mean
end
for l = 1:3
    for k = 1:nF - l
        A = reshape(Zc(:, k, :), p, nR); B = reshape(Zc(:, k + l, :), p, nR);
        ok = all(isfinite(A), 1) & all(isfinite(B), 1);
        if sum(ok) < 20, continue; end
        A = A(:, ok) - mean(A(:, ok), 2); B = B(:, ok) - mean(B(:, ok), 2);
        r = sum(A .* B, 2) ./ sqrt(sum(A.^2, 2) .* sum(B.^2, 2));
        s.lagCorr(l, k) = mean(abs(r));
    end
end
end

function styleAx(ax, ink, gridInk, yl, ttl, linearY)
if nargin < 6 || ~linearY, set(ax, 'YScale', 'log'); end
set(ax, 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off');
grid(ax, 'on'); set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
xlabel(ax, '\omega (rad/s)'); ylabel(ax, yl);
title(ax, ttl, 'FontWeight', 'normal');
end

function s = flagStr(ok)
if ok, s = 'PASS'; else, s = 'FLAG'; end
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);     % R2020a+
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end