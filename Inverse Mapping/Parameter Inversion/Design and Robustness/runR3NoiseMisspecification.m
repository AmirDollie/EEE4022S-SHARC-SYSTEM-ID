%% runR3NoiseMisspecification.m
% STUDY R3 (SHOULD): how accurately must the reference sensor's noise floor be known?
%
% QUESTION. The estimator removes the reference noise from the denominator,
%     T_hat_j = S_jr / (S_rr - S_nn,assumed),
% so an error in the assumed noise floor biases every transmissibility in the same direction.
% All spectra are the window-smoothed (Welch-expected) ones, so with T_W = S_jr,W / S_rr,sig,W the
% additional bias relative to the correctly corrected Welch expectation is b = T_W Delta / (S_rr,sig,W - Delta)
% (exact in expectation for a fixed bin set) and the total feature error is b_W + b.
% Delta = S_nn,assumed - S_nn,true. For a white error Delta = eps S_nn: eps = -1 is "no correction at all".
%
% METHOD (pre-registered, cfg.R3). eps in cfg.R3.epsNoise at every noise case in cfg.R3.noiseCases, at p0 with Level 2C,
% bias projected by the E1 engine (dtheta_sys, d_sys, combined with the Welch bias as d_tot and the
% class). Tolerances by first crossing on the fine eps grid cfg.R3.toleranceSearch / toleranceN:
% |delta beta/beta| = 1%, d_sys = 0.5, 1. The bias is the EXACT expected one (window-smoothed reference
% spectrum, full ratio), not the first-order T eps S_nn / (S_rr - S_nn); the latter is its small-eps limit.
% Coloured case (cfg.R3.colouredCase, a SECONDARY check reported separately from the tolerance table): the true reference noise has the MEASURED basin-band shape of the
% LSM6DS3TR-C z axis (16.9 h record; power law fitted over the band, scaled to the measured band mean),
% while the estimator subtracts the white band mean. The same is reported for the other two measured
% sensors (additional, not pre-registered). The covariance Sigma_k uses the white level (the shape
% changes it by a few percent; the bias is the quantity of interest).
% Consistency check: the bias built here is passed through evaluateScenario's own sysBias path for
% one case and must give the same d_sys.
%
% DELIVERABLE. How accurately the noise floor must be known, per sensor type.
% CAVEAT. The real estimator rejects bins where its denominator falls below 0.1 S_rr; the analysis keeps
% the bins accepted under the TRUE noise and reports how many would flip under the assumed one.
%
% COST. No new EMM solves; seconds.
% OUTPUT. Results/R3/ (log, .mat, summary, CSV tables, figures/R3_fig1_noise_misspecification via
% plotR3NoiseMisspecification, which can also redraw the figure later from the .mat).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();
run = startStudy('R3', 'Noise-floor mis-specification', P.results);

try
noiseNames = cfg.R3.noiseCases;
nN = numel(noiseNames);
epsGrid = cfg.R3.epsNoise(:).';
epsFine = linspace(cfg.R3.toleranceSearch(1), cfg.R3.toleranceSearch(2), cfg.R3.toleranceN);
headline = cfg.R3.headlineNoise;
hIdx = find(strcmp(noiseNames, headline));
minDenFrac = 0.1;                                   % estimator rule: denominator > 0.1 S_rr

%% ---- measured noise shapes (z axis, basin band) ---------------------------------------------------------
nf = dir(fullfile(P.imuResults, 'staticNoise_*_long.mat'));
if isempty(nf), error('runR3:noiseFile', 'Measured noise record (staticNoise_*_long.mat) not found.'); end
[~, kNew] = max([nf.datenum]);                     % newest measured record
nz = load(fullfile(P.imuResults, nf(kNew).name)); nz = nz.out.psd.acc{3};   % acc{3} = z axis
measuredMap = struct('lsm6ds3trc', 'LSM6DS3TR-C', 'lsm6dsv16x', 'LSM6DSV16X', 'icm42688p', 'ICM-42688-P');
primaryColoured = regexp(cfg.R3.colouredCase, '^[a-z0-9]+', 'match', 'once');   % pre-registered case
colouredNames = [{primaryColoured}, setdiff({'lsm6ds3trc', 'lsm6dsv16x', 'icm42688p'}, {primaryColoured}, 'stable')];
if ~isfield(measuredMap, primaryColoured)
    error('runR3:colouredCase', 'cfg.R3.colouredCase ''%s'' is not a measured sensor.', cfg.R3.colouredCase);
end
band = cfg.band.basin;
shape = struct();
for c = 1:numel(colouredNames)
    k = find(strcmp(nz.channelNames, measuredMap.(colouredNames{c})));
    in = nz.omega >= band(1) & nz.omega <= band(2);
    w = nz.omega(in); Pm = nz.PxxPerRadS(in, k);
    pf = polyfit(log(w), log(Pm), 1);                 % power law over the basin band
    fit = @(om) exp(polyval(pf, log(om)));
    scale = cfg.noise.(colouredNames{c}) / mean(fit(w));   % match the band mean used everywhere else
    shape.(colouredNames{c}) = struct('exponent', pf(1), 'logCoeffs', pf, 'scale', scale, ...
        'measuredMean', mean(Pm), 'maxDevFromFit', max(abs(movmeanSafe(Pm, 21) ./ fit(w) - 1)));
    fprintf('Measured %s z noise over the basin band: PSD ~ omega^%.2f, mean %.3g (cfg %.3g), smoothed data within %.1f%% of the fit\n', ...
        measuredMap.(colouredNames{c}), pf(1), mean(Pm), cfg.noise.(colouredNames{c}), 100 * shape.(colouredNames{c}).maxDevFromFit);
end

%% ---- sweep ---------------------------------------------------------------------------------------------
res = struct('noise', noiseNames);
for n = 1:nN
    Snn = cfg.noise.(noiseNames{n});
    m = evaluateScenario(struct('Snn', Snn, 'cacheDir', P.cache), cfg);
    [Srr, Tm] = referenceSignalPSD(m, cfg);           % smoothed reference signal PSD and Welch-expected T_W
    pBlk = size(m.blocksAll, 1);
    rows = reshape((find(m.validMask) - 1) * pBlk + (1:pBlk).', [], 1);
    Jv = m.Jall(rows, :); Bv = m.blocksAll(:, :, m.validMask);
    biasOf = @(Delta) stackTransmissibility(Tm .* (Delta ./ (Srr - Delta)));
    % pre-registered grid
    nE = numel(epsGrid);
    Bmat = zeros(numel(m.f0), nE); flips = zeros(1, nE); blown = false(1, nE);
    for k = 1:nE
        Delta = epsGrid(k) * Snn * ones(size(Srr));
        den = Srr - Delta;
        blown(k) = any(den(m.validMask) <= 0);
        flips(k) = mean(den(m.validMask) < minDenFrac * (Srr(m.validMask) + Snn));
        Bmat(:, k) = biasOf(Delta);
    end
    q = fisherFromBlocks(Jv, Bv, [Bmat(rows, :), Bmat(rows, :) + m.bW(rows)]);
    res(n).Snn = Snn; res(n).sigma = m.sigma(:).'; res(n).dW = m.dW; res(n).condF = m.condF;
    res(n).dthetaSys = q.dtheta(:, 1:nE);
    res(n).dSys = q.d(1:nE); res(n).dTot = q.d(nE + (1:nE));
    res(n).dSys(blown) = Inf; res(n).dTot(blown) = Inf;
    res(n).dBeta = exp(res(n).dthetaSys(1, :)) - 1; res(n).dR = exp(res(n).dthetaSys(2, :)) - 1;
    res(n).flipFraction = flips;
    res(n).class = arrayfun(@(k) classOf(max(m.sigma), res(n).dTot(k), m.condF, m.singular, cfg), 1:nE, ...
        'UniformOutput', false);
    res(n).minSrrOverSnn = min(Srr(m.validMask)) / Snn;
    % fine grid for tolerances (the projection is linear in b, b is nonlinear in eps)
    Bf = zeros(numel(rows), numel(epsFine));
    for k = 1:numel(epsFine)
        bk = biasOf(epsFine(k) * Snn * ones(size(Srr))); Bf(:, k) = bk(rows);
    end
    qf = fisherFromBlocks(Jv, Bv, Bf);
    okF = arrayfun(@(e) all(Srr(m.validMask) - e * Snn > 0), epsFine);
    dF = qf.d; dF(~okF) = Inf;
    bF = abs(exp(qf.dtheta(1, :)) - 1); bF(~okF) = Inf;
    res(n).epsFine = epsFine; res(n).dSysFine = dF; res(n).dBetaFine = bF;   % for the figure
    res(n).tolBeta1 = firstCrossing(epsFine, bF, 0.01);
    res(n).tolD05 = firstCrossing(epsFine, dF, 0.5);
    res(n).tolD1 = firstCrossing(epsFine, dF, 1.0);
    % coloured noise (measured sensors only)
    res(n).coloured = [];
    if isfield(shape, noiseNames{n})
        sh = shape.(noiseNames{n});
        Strue = sh.scale * exp(polyval(sh.logCoeffs, log(m.omega(:).')));
        Delta = Snn - Strue;                          % assumed white mean minus true coloured
        bc = biasOf(Delta);
        qc = fisherFromBlocks(Jv, Bv, [bc(rows), bc(rows) + m.bW(rows)]);
        res(n).coloured = struct('exponent', shape.(noiseNames{n}).exponent, 'dthetaSys', qc.dtheta(:, 1).', ...
            'dSys', qc.d(1), 'dTot', qc.d(2), 'dBeta', exp(qc.dtheta(1, 1)) - 1, 'maxRelDelta', max(abs(Delta)) / Snn);
    end
    fprintf('  %-11s S_nn %.3g, min S_rr/S_nn %.0f: tolerance |eps| for 1%% beta %s, d 0.5 %s, d 1 %s; no correction (eps = -1): d_sys %.3f, dbeta %+.3f%%\n', ...
        noiseNames{n}, Snn, res(n).minSrrOverSnn, tolStr(res(n).tolBeta1), tolStr(res(n).tolD05), tolStr(res(n).tolD1), ...
        res(n).dSys(epsGrid == -1), 100 * res(n).dBeta(epsGrid == -1));
    if ~isempty(res(n).coloured)
        fprintf('              coloured (measured shape, omega^%.2f, max |Delta|/S_nn %.2f): d_sys %.4f, dbeta %+.4f%%\n', ...
            res(n).coloured.exponent, res(n).coloured.maxRelDelta, res(n).coloured.dSys, 100 * res(n).coloured.dBeta);
    end
end

%% ---- consistency with the engine's own sysBias path -----------------------------------------------------
Snn = cfg.noise.(headline); e0 = 0.3;
mC = evaluateScenario(struct('Snn', Snn, 'cacheDir', P.cache, 'sysBias', ...
    @(T, omega, ctx) noiseMisspecBias(ctx, e0 * Snn)), cfg);   % T_W rebuilt from the engine's own smoothed spectra
dHere = res(hIdx).dSys(epsGrid == e0);
fprintf('\nConsistency: d_sys at eps = +30%% (%s) here %.6f, through evaluateScenario %.6f (rel. diff %.1e)\n', ...
    headline, dHere, mC.dSys, abs(dHere / mC.dSys - 1));

%% ---- tables --------------------------------------------------------------------------------------------
rowsT = {};
for n = 1:nN
    for k = 1:numel(epsGrid)
        rowsT(end + 1, :) = {noiseNames{n}, epsGrid(k), 100 * res(n).dBeta(k), 100 * res(n).dR(k), res(n).dSys(k), ...
            res(n).dTot(k), res(n).flipFraction(k), res(n).class{k}}; %#ok<AGROW>
    end
end
writeStudyCsv(run, 'sweep', {'noise_case', 'eps_noise', 'dbeta_pct', 'dR_pct', 'd_sys', 'd_tot', ...
    'bins_flipped_fraction', 'class'}, rowsT);
rowsT = cell(nN, 7);
for n = 1:nN
    rowsT(n, :) = {noiseNames{n}, res(n).Snn, res(n).minSrrOverSnn, 100 * res(n).tolBeta1, 100 * res(n).tolD05, ...
        100 * res(n).tolD1, res(n).dSys(epsGrid == -1)};
end
writeStudyCsv(run, 'tolerance', {'noise_case', 'Snn', 'min_Srr_over_Snn', 'tol_beta1pct_pct', 'tol_d05_pct', ...
    'tol_d1_pct', 'd_sys_no_correction'}, rowsT);
rowsT = {};
for n = 1:nN
    if ~isempty(res(n).coloured)
        c = res(n).coloured;
        rowsT(end + 1, :) = {noiseNames{n}, c.exponent, c.maxRelDelta, 100 * c.dBeta, c.dSys, c.dTot}; %#ok<AGROW>
    end
end
writeStudyCsv(run, 'coloured', {'noise_case', 'psd_exponent', 'max_rel_Delta', 'dbeta_pct', 'd_sys', 'd_tot'}, rowsT);

%% ---- summary -------------------------------------------------------------------------------------------
L = {};
L{end + 1} = sprintf('Setting: p0, Level 2C (ref s26), N = %d, L = %d, JONSWAP Hs %.2f m.', cfg.acq.Nref, cfg.welch.L, cfg.sea.Hs);
L{end + 1} = sprintf('Consistency with evaluateScenario sysBias path: rel. diff %.1e.', abs(dHere / mC.dSys - 1));
L{end + 1} = '';
L{end + 1} = 'DELIVERABLE: required accuracy of the reference noise floor (first crossing, smaller |eps| of the two signs;';
L{end + 1} = '             eps = -100% is no correction; "not reached" = not reached within the search range)';
L{end + 1} = sprintf('  %-12s %9s %10s %12s %12s %12s %16s', 'noise case', 'S_nn', 'minSNR', '|db/b|<=1%', 'd_sys<=0.5', ...
    'd_sys<=1', 'd_sys no corr.');
for n = 1:nN
    L{end + 1} = sprintf('  %-12s %9.3g %10.0f %12s %12s %12s %16.3f', noiseNames{n}, res(n).Snn, res(n).minSrrOverSnn, ...
        tolStr(res(n).tolBeta1), tolStr(res(n).tolD05), tolStr(res(n).tolD1), res(n).dSys(epsGrid == -1)); %#ok<AGROW>
end
L{end + 1} = '  minSNR = smallest S_rr,signal / S_nn over the accepted bins.';
L{end + 1} = '';
L{end + 1} = 'Coloured noise (true = measured basin-band shape, assumed = white band mean):';
for n = 1:nN
    if ~isempty(res(n).coloured)
        c = res(n).coloured;
        L{end + 1} = sprintf('  %-12s PSD ~ omega^%.2f, max |Delta|/S_nn %.2f: dbeta %+.4f%%, d_sys %.4f, d_tot %.3f', ...
            noiseNames{n}, c.exponent, c.maxRelDelta, 100 * c.dBeta, c.dSys, c.dTot); %#ok<AGROW>
    end
end
results = struct('epsGrid', epsGrid, 'res', res, 'headline', headline, 'shape', shape);
try
    plotR3NoiseMisspecification(results, cfg, run);   % figure (also redrawable later from the .mat)
catch errPlot                                         % a figure problem must not lose the results
    warning('runR3:plot', 'Figure failed (%s); results are still saved. Redraw with plotR3NoiseMisspecification.', ...
        errPlot.message);
end
finishStudy(run, results, cfg, L);
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================
function [Srr, TW] = referenceSignalPSD(m, cfg) %#ok<INUSD>
% window-smoothed SIGNAL auto-spectrum of the reference at the bins and the Welch-expected
% transmissibility T_W = S_jr,W / S_rr,sig,W from the same smoothed matrix, built as evaluateScenario does
scn = m.scn;
sea = scn.sea;
if isa(sea, 'function_handle'), a = {'Spectrum', sea};
else, a = {'Spectrum', 'jonswap', 'Hs', sea.Hs, 'PeakFrequency', sea.omegaP, 'PeakEnhancement', sea.gammaJ};
end
[~, ~, info] = synthesiseTwinRecords(m.field.twin, scn.N, scn.dt, 'Seed', 1, a{:});
Ssm = welchSmoothedSpectrum(m.omega, info, scn.L, scn.layout);
[Srr, TW] = welchExpectedT(Ssm, scn.ref);
end

function [Srr, TW] = welchExpectedT(Ssm, ref)
Srr = real(squeeze(Ssm(ref, ref, :))).';
other = setdiff(1:size(Ssm, 1), ref);               % same order as stackTransmissibility / the engine
Sjr = reshape(Ssm(other, ref, :), numel(other), []);
TW = Sjr ./ Srr;
end

function b = noiseMisspecBias(ctx, Delta)
% R3 bias relative to the correctly corrected Welch expectation: T_W Delta / (S_rr,sig,W - Delta)
[Srr, TW] = welchExpectedT(ctx.Ssmoothed, ctx.scn.ref);
b = stackTransmissibility(TW .* (Delta ./ (Srr - Delta)));
end

function y = movmeanSafe(x, k)
% centred moving average (Octave and MATLAB)
x = x(:); h = floor(k / 2); y = zeros(size(x));
for i = 1:numel(x)
    y(i) = mean(x(max(1, i - h):min(numel(x), i + h)));
end
end

function e = firstCrossing(epsGrid, metric, level)
e = NaN; c = [];
ip = find(epsGrid > 0 & metric >= level, 1, 'first');
if ~isempty(ip), c(end + 1) = interpCross(epsGrid, metric, ip, level); end
in = find(epsGrid < 0 & metric >= level, 1, 'last');
if ~isempty(in), c(end + 1) = interpCross(epsGrid, metric, in, level); end
if ~isempty(c), e = min(abs(c)); end
end

function x = interpCross(e, m, i, level)
j = round(i - sign(e(i)));
if j < 1 || j > numel(e) || m(j) >= level || ~isfinite(m(i)), x = e(i); return; end
x = e(j) + (level - m(j)) * (e(i) - e(j)) / (m(i) - m(j));
end

function s = tolStr(e)
if isnan(e), s = 'not reached'; else, s = sprintf('%.1f%%', 100 * e); end
end

function c = classOf(sMax, dTot, condF, singular, cfg)
I = cfg.class.identifiable; M = cfg.class.marginal;
if singular || ~isfinite(dTot), c = cfg.class.names{3}; return; end
if sMax <= I.maxSigma && dTot <= I.maxBiasDistance && condF <= I.maxCondF, c = cfg.class.names{1};
elseif sMax <= M.maxSigma && dTot <= M.maxBiasDistance && condF <= M.maxCondF, c = cfg.class.names{2};
else, c = cfg.class.names{3};
end
end