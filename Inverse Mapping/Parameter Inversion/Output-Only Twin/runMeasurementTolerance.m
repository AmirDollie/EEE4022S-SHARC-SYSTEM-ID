%% runMeasurementTolerance.m
% Measurement-chain tolerances from the frozen output-only inversion (U2).
% How large may a timing offset, a clock-rate mismatch or a relative
% scale-factor error between IMUs be before it biases (beta, R)?
%
% ERROR MODELS (sensor k perturbed, all others perfect; features on the U2 grid)
%   timing   sensor k records y_k(t - dt_k(t)), dt_k(t) = dt0 + eps * (t - t_sync)
%            -> each Welch segment m picks up exp(-i w dt_k(t_m)); the Welch average
%            of the cross-spectrum gives
%              T_j -> T_j * c_j(w),  c_j(w) = mean_m exp(-i w [dt_j(t_m) - dt_r(t_m)])
%            (t_m = segment centres of the actual N, L, overlap). A constant offset is
%            a pure phase ramp; a rate mismatch adds a mean offset eps*(tbar - t_sync)
%            AND a magnitude loss (the phase rotates across segments, roughly
%            sinc(w eps T_rec / 2)). Two synchronisation conventions for eps:
%            'start' (clocks agree at the record start) and 'centre' (at mid-record:
%            the mean offset vanishes, only the rotation loss remains).
%   scale    T_j -> T_j (1 + s_j) / (1 + s_r): a real factor at every frequency.
%   An error on a non-reference sensor j changes only T_j; an error on the reference
%   (s26) changes all three transmissibilities.
%
% PROJECTION (the U2 truth points A and B, their J, Sigma blocks and F)
%   b = z(T_perturbed) - z(T_true)          (exact in feature space)
%   dtheta = F^-1 J' Sigma^-1 b,  d = sqrt(dtheta' F dtheta)   (linear in theta)
%   Reported: d (in units of one-record sigma, so specific to U2's 10% noise and
%   109-minute record) and dtheta in percent, 100(exp(dtheta) - 1), which does not
%   depend on an overall noise scaling. Tolerance = largest error in the contiguous safe
%   interval starting from zero (first crossing of d = 0.1, 0.3; the rate-mismatch response
%   can become non-monotone, so later dips are deliberately ignored),
%   worst case over the four sensors, interpolated on the log sweep.
%
% VERIFICATION (VERIFY = true, no EMM solves): noise-free S3 records at truth A are
% warped in time (spline resampling of one sensor onto t - dt(t)) or scaled, pushed
% through the S2 estimator, and the feature change (warped minus unwarped, same
% records) is compared with the analytic b at a few test points.
%
% Also reported: tolerance by parameter bias (BIAS_LEVELS, percent), the number to carry
% to the hardware spec. It is independent of an overall covariance scaling (unlike d) but
% still depends on the band and its weighting, and, for rate mismatch, on the record length
% (dt = eps t grows with T_rec). Every error is swept with both signs and the worse sign kept;
% the +/- asymmetry is printed (small errors are linear, so it should be small).
%
% INPUT: newest transmissibilityTwinInversion_*.mat in ../Results (U2).
% OUTPUT (../Results): measurementTolerance_<stamp>.mat / .txt / .png
%
% Lives in Inverse Mapping/Parameter Inversion/Output-Only Twin/.

thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..', '..', '..');
addpath(fullfile(thisDir, '..', '..'));
addpath(fullfile(rootDir, 'Transmissibility'));
addpath(fullfile(rootDir, 'SSI'));
addpath(fullfile(rootDir, 'Forward Model'));
addpath(fullfile(rootDir, 'Forward Model', 'Animation'));
clearvars -except thisDir rootDir
close all, clc

%% ---- Configuration -----------------------------------------------------------
U2_FILE = '';                                   % '' = newest transmissibilityTwinInversion_*.mat
DT0_SWEEP = logspace(-6, -1, 51);               % constant offset: 1 us .. 100 ms
EPS_SWEEP = logspace(-10, -3, 71);              % rate mismatch: 1e-4 .. 1000 ppm
SCALE_SWEEP = logspace(-5, -1, 41);             % relative scale error: 0.001 % .. 10 %
D_LEVELS = [0.1 0.3];
BIAS_LEVELS = [0.5 1 2];                        % percent, max(|beta bias|, |R bias|)
VERIFY = true; N_VERIFY = 5;
VERIFY_POINTS = struct('kind', {'dt0', 'eps', 'scale'}, 'value', {5e-3, 1e-5, 1e-2}, 'sensor', {1, 1, 4});
HS = 0.05; WP = 5.0; GAMMA_J = 3.3;             % sea state used by U2

%% ---- Output ----------------------------------------------------------------------
resultsDir = fullfile(thisDir, '..', 'Results');
cacheDir = fullfile(resultsDir, 'twinCache');
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['measurementTolerance_', stamp]);
diary([baseName, '.txt']);
try

if isempty(U2_FILE)
    d = dir(fullfile(resultsDir, 'transmissibilityTwinInversion_*.mat'));
    if isempty(d), error('runMeasurementTolerance:noU2', 'No U2 result in %s.', resultsDir); end
    [~, newest] = max([d.datenum]);
    U2_FILE = fullfile(resultsDir, d(newest).name);
end
U = load(U2_FILE); U = U.out;
omega = U.omega(:).'; nF = numel(omega); N = U.N; dt = U.dt; L = U.L; ov = U.overlap;
nS = size(U.sensors, 1); ref = U.reference; others = setdiff(1:nS, ref); m = numel(others); p = 2 * m;
sensorNames = {'s2', 's6', 's14', 's26 (ref)'};
hop = max(1, round(L * (1 - ov)));
tSeg = ((1:hop:(N - L + 1)) - 1 + L / 2) * dt;        % segment centre times (s)
Trec = N * dt;
fprintf('runMeasurementTolerance  (%s)\nU2 input: %s\n', stamp, U2_FILE);
fprintf('%d frequencies %.3f to %.3f rad/s; record %.0f s, %d Welch segments (L = %d)\n\n', nF, omega([1 end]), ...
    Trec, numel(tSeg), L);

%% ---- 1  Analytic sweeps ---------------------------------------------------------------------------------
kinds = {'dt0', 'epsStart', 'epsCentre', 'scale'};
sweeps = {DT0_SWEEP, EPS_SWEEP, EPS_SWEEP, SCALE_SWEEP};
units = {'s', '', '', ''};
nT = numel(U.truth);
R = struct();
for it = 1:nT
    tr = U.truth(it);
    Ttrue = unstackTransmissibility(tr.fTrue(:), m);                 % m x nF complex
    for q = 1:numel(kinds)
        x = sweeps{q};
        dS = zeros(nS, numel(x), 2); thS = zeros(2, nS, numel(x), 2);      % last index: +x, -x
        for sg = 1:2
            for k = 1:nS
                for i = 1:numel(x)
                    Tp = perturbT(Ttrue, omega, kinds{q}, (3 - 2 * sg) * x(i), k, ref, others, tSeg, Trec);
                    b = reshape(stackTransmissibility(Tp), p, nF) - tr.fTrue;
                    [dS(k, i, sg), thS(:, k, i, sg)] = project(tr.J, tr.blocks, tr.F, b);
                end
            end
        end
        pctS = 100 * squeeze(max(abs(exp(thS) - 1), [], 1));                 % nS x nx x 2
        pctS = reshape(pctS, nS, numel(x), 2);
        [dAll, sgD] = max(dS, [], 3);                                        % worse sign in d
        [pctAll, sgP] = max(pctS, [], 3);                                    % worse sign in percent bias
        thD = zeros(2, nS, numel(x)); thP = thD;
        for k = 1:nS
            for i = 1:numel(x)
                thD(:, k, i) = thS(:, k, i, sgD(k, i)); thP(:, k, i) = thS(:, k, i, sgP(k, i));
            end
        end
        R(it).(kinds{q}).x = x; R(it).(kinds{q}).d = dAll; R(it).(kinds{q}).dtheta = thD;
        R(it).(kinds{q}).pct = pctAll; R(it).(kinds{q}).dthetaPct = thP;
        R(it).(kinds{q}).dSigned = dS; R(it).(kinds{q}).dthetaSigned = thS;
        R(it).(kinds{q}).asymmetry = abs(dS(:, :, 1) - dS(:, :, 2)) ./ max(dAll, realmin);   % relative +/- difference
    end
    R(it).name = tr.name;
end

%% ---- 2  Tolerances ------------------------------------------------------------------------------------------
labels = {'constant offset dt0', 'rate mismatch eps (synced at record start)', ...
    'rate mismatch eps (synced at record centre)', 'relative scale factor s_k - s_r'};
fmt = {@(v) sprintf('%.3g ms', 1e3 * v), @(v) sprintf('%.3g ppm', 1e6 * v), @(v) sprintf('%.3g ppm', 1e6 * v), ...
    @(v) sprintf('%.3g %%', 100 * v)};
tol = struct();
fprintf('TOLERANCES (worst sensor and sign; d in units of one U2 record sigma; beta/R bias at the tolerance)\n');
for q = 1:numel(kinds)
    fprintf('  %s\n', labels{q});
    for it = 1:nT
        S = R(it).(kinds{q});
        for lv = 1:numel(D_LEVELS)
            [xt, kWorst, below] = toleranceAt(S.x, S.d, D_LEVELS(lv));
            if below
                fprintf('    truth %-24s d < %.1f: below the sweep, < %s (extend the sweep)\n', R(it).name, ...
                    D_LEVELS(lv), fmt{q}(xt));
                continue
            end
            if isnan(xt)
                fprintf('    truth %s: d < %.1f over the whole sweep\n', R(it).name, D_LEVELS(lv));
                continue
            end
            thAt = interpDtheta(S.x, squeeze(S.dtheta(:, kWorst, :)), xt);
            fprintf('    truth %-24s d < %.1f: %-12s (worst: %-9s)  beta %+.3f %%, R %+.4f %%\n', R(it).name, D_LEVELS(lv), ...
                fmt{q}(xt), sensorNames{kWorst}, 100 * (exp(thAt(1)) - 1), 100 * (exp(thAt(2)) - 1));
            tol(it).(kinds{q})(lv) = struct('level', D_LEVELS(lv), 'value', xt, 'sensor', kWorst, 'dtheta', thAt);
        end
    end
end
fprintf('\n(d is linear in the error while it is small: the parameter bias is linearised about the truth)\n');
fprintf('SIGN ASYMMETRY at the d = %.1f tolerance, |d(+x) - d(-x)| / d (worst sensor):\n', D_LEVELS(1));
for q = 1:numel(kinds)
    for it = 1:nT
        S = R(it).(kinds{q});
        [xt, kWorst, below] = toleranceAt(S.x, S.d, D_LEVELS(1));
        if isnan(xt) || below, continue; end
        a = interp1(log(S.x), S.asymmetry(kWorst, :), log(xt));
        fprintf('  %-45s truth %-24s %.2e\n', labels{q}, R(it).name, a);
    end
end
fprintf('\n');

fprintf('TOLERANCES BY PARAMETER BIAS (worst sensor and sign; max(|beta bias|, |R bias|) in %%)\n');
fprintf('(independent of an overall noise scaling, unlike d; still depends on the band and its weighting,\n');
fprintf(' and for rate mismatch explicitly on the record length)\n');
tolPct = struct();
for q = 1:numel(kinds)
    fprintf('  %s\n', labels{q});
    for it = 1:nT
        S = R(it).(kinds{q});
        pct = S.pct;
        for lv = 1:numel(BIAS_LEVELS)
            [xt, kWorst, below] = toleranceAt(S.x, max(pct, realmin), BIAS_LEVELS(lv));
            if below
                fprintf('    truth %-24s bias < %.1f %%: below the sweep, < %s (extend the sweep)\n', R(it).name, ...
                    BIAS_LEVELS(lv), fmt{q}(xt));
                continue
            end
            if isnan(xt)
                fprintf('    truth %-24s bias < %.1f %%: over the whole sweep\n', R(it).name, BIAS_LEVELS(lv));
                continue
            end
            thAt = interpDtheta(S.x, squeeze(S.dthetaPct(:, kWorst, :)), xt);
            fprintf('    truth %-24s bias < %.1f %%: %-12s (worst: %-9s)  beta %+.3f %%, R %+.4f %%\n', R(it).name, ...
                BIAS_LEVELS(lv), fmt{q}(xt), sensorNames{kWorst}, 100 * (exp(thAt(1)) - 1), 100 * (exp(thAt(2)) - 1));
            tolPct(it).(kinds{q})(lv) = struct('levelPct', BIAS_LEVELS(lv), 'value', xt, 'sensor', kWorst, 'dtheta', thAt);
        end
    end
end
fprintf('\n');

%% ---- 3  Verification with warped S3 records -------------------------------------------------------------------
ver = struct([]);
if VERIFY
    twin = synthesiseTwinRecords(struct('p', U.truth(1).p, 'sensors', U.sensors, 'nodeSpacing', 0.1), 'CacheDir', cacheDir);
    specJ = {'Spectrum', 'jonswap', 'Hs', HS, 'PeakFrequency', WP, 'PeakEnhancement', GAMMA_J};
    eOpt = {'Reference', ref, 'SegmentLength', L, 'Overlap', ov, 'Band', [3 8.5]};
    tr = U.truth(1); Ttrue = unstackTransmissibility(tr.fTrue(:), m);
    t = (0:N - 1) * dt;
    fprintf('VERIFICATION at truth %s (%d noise-free records, twin %s)\n', tr.name, N_VERIFY, twin.cacheStatus);
    for v = 1:numel(VERIFY_POINTS)
        vp = VERIFY_POINTS(v);
        kind = vp.kind; if strcmp(kind, 'eps'), kind = 'epsStart'; end
        bSim = zeros(p, nF);
        for r = 1:N_VERIFY
            Y = synthesiseTwinRecords(twin, N, dt, 'Seed', 90000 + r, specJ{:});
            Yp = perturbRecord(Y, t, kind, vp.value, vp.sensor);
            z0 = featuresOf(Y, dt, eOpt, U.kBins, p);
            z1 = featuresOf(Yp, dt, eOpt, U.kBins, p);
            bSim = bSim + (z1 - z0) / N_VERIFY;
        end
        Tp = perturbT(Ttrue, omega, kind, vp.value, vp.sensor, ref, others, tSeg, Trec);
        bAna = reshape(stackTransmissibility(Tp), p, nF) - tr.fTrue;
        [dS, thS] = project(tr.J, tr.blocks, tr.F, bSim);
        [dA, thA] = project(tr.J, tr.blocks, tr.F, bAna);
        relB = norm(bSim(:) - bAna(:)) / norm(bAna(:));
        fprintf('  %-8s %-10g on %-9s: d simulated %.3f, analytic %.3f; dtheta sim [%+.2e %+.2e] ana [%+.2e %+.2e]; |b_sim-b_ana|/|b_ana| %.2f\n', ...
            vp.kind, vp.value, sensorNames{vp.sensor}, dS, dA, thS, thA, relB);
        ver(v).kind = vp.kind; ver(v).value = vp.value; ver(v).sensor = vp.sensor; %#ok<SAGROW>
        ver(v).dSim = dS; ver(v).dAna = dA; ver(v).thetaSim = thS; ver(v).thetaAna = thA; ver(v).relB = relB; %#ok<SAGROW>
    end
    fprintf('  (the analytic form weights Welch segments equally, a record weights them by their power:\n');
    fprintf('   agreement to a few percent in d is expected, and exact for dt0 and scale)\n\n');
end

%% ---- 4  Save and figure -----------------------------------------------------------------------------------------
out = struct('u2File', U2_FILE, 'omega', omega, 'N', N, 'dt', dt, 'L', L, 'overlap', ov, 'tSeg', tSeg, ...
    'sweeps', R, 'tolerances', tol, 'dLevels', D_LEVELS, 'tolerancesPct', tolPct, 'biasLevels', BIAS_LEVELS, 'verification', ver, 'sensorNames', {sensorNames});
save([baseName, '.mat'], 'out');
fprintf('Saved: %s.mat\n', baseName);

cols = {[42 120 214] / 255, [235 104 52] / 255, [27 175 122] / 255, [120 120 120] / 255};
xl = {'constant offset \delta t_0 (s)', 'rate mismatch \epsilon (synced at start)', ...
    'rate mismatch \epsilon (synced at centre)', 'relative scale error s_k - s_r'};
fig = figure('Color', 'w', 'Position', [60 60 1400 380]);
for q = 1:numel(kinds)
    ax = subplot(1, numel(kinds), q); hold(ax, 'on');
    h = [];
    for k = 1:nS
        h(end + 1) = plot(ax, R(1).(kinds{q}).x, R(1).(kinds{q}).d(k, :), '-', 'Color', cols{k}, 'LineWidth', 1.6); %#ok<AGROW>
        if nT > 1
            plot(ax, R(2).(kinds{q}).x, R(2).(kinds{q}).d(k, :), '--', 'Color', cols{k}, 'LineWidth', 1.2);
        end
    end
    for lv = 1:numel(D_LEVELS)
        plot(ax, R(1).(kinds{q}).x([1 end]), D_LEVELS(lv) * [1 1], ':', 'Color', 'k');
    end
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 9, 'Box', 'off'); grid(ax, 'on');
    ylim(ax, [1e-3 1e2]);
    xlabel(ax, xl{q}); ylabel(ax, 'd_{bias}');
    if q == 1, legend(ax, h, sensorNames, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8); end
end
sgtitleSafe(fig, 'Measurement tolerances: solid truth A, dashed truth B; dotted d = 0.1, 0.3');
saveFigure(fig, [baseName, '.png']);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ================================================================================================
%  Local functions
%% ================================================================================================
function Tp = perturbT(T, omega, kind, x, k, ref, others, tSeg, Trec)
% Welch-expected transmissibility when sensor k carries one error of size x.
w = omega(:).';
switch kind
    case 'dt0'
        dtk = x * ones(size(tSeg));
    case 'epsStart'
        dtk = x * tSeg;
    case 'epsCentre'
        dtk = x * (tSeg - Trec / 2);
    case 'scale'
        g = 1 + x;
end
if strcmp(kind, 'scale')
    fac = g * ones(size(w));                                          % gain of sensor k
    if k == ref
        Tp = T ./ fac;
    else
        Tp = T; Tp(others == k, :) = T(others == k, :) .* fac;
    end
    return
end
c = mean(exp(-1i * w(:) * dtk(:).'), 2).';                          % 1 x nF: factor on sensor k's spectrum
if k == ref
    Tp = T .* conj(c);                                                % exact Welch form: S_jr -> S_jr * conj(c_r)
else
    Tp = T; Tp(others == k, :) = T(others == k, :) .* c;
end
end

function Yp = perturbRecord(Y, t, kind, x, k)
% Apply one measurement error to sensor k of a record (for the verification).
Yp = Y;
switch kind
    case 'dt0'
        tk = t - x;
    case 'epsStart'
        tk = t - x * t;
    case 'scale'
        Yp(k, :) = (1 + x) * Y(k, :);
        return
end
Yp(k, :) = interp1(t, Y(k, :), tk, 'spline', 0);
end

function z = featuresOf(Y, dt, eOpt, kBins, p)
e = extractTransmissibility(Y, dt, eOpt{:});
[~, loc] = ismember(kBins, e.binIndex);
z = reshape(stackTransmissibility(e.T(:, loc)), p, numel(kBins));
end

function [d, dtheta] = project(J, blocks, F, b)
[p, ~, nF] = size(blocks);
b = reshape(b, p, nF);
g = zeros(size(J, 2), 1);
for k = 1:nF
    Jk = J((k - 1) * p + (1:p), :);
    g = g + Jk.' * (blocks(:, :, k) \ b(:, k));
end
dtheta = F \ g;
d = sqrt(dtheta.' * F * dtheta);
end

function [xt, kWorst, below] = toleranceAt(x, d, level)
% Largest x in the contiguous safe interval from zero (first crossing, log-interpolated),
% worst case over sensors. below = true if some sensor already exceeds the level at the
% first sweep point: the tolerance is then only known to be < x(1) (extend the sweep).
xt = Inf; kWorst = NaN; below = false;
for k = 1:size(d, 1)
    i = find(d(k, :) >= level, 1);
    if isempty(i), continue; end
    if i == 1
        below = true; xt = x(1); kWorst = NaN;
        return
    end
    f = (log(level) - log(d(k, i - 1))) / (log(d(k, i)) - log(d(k, i - 1)));
    xk = exp(log(x(i - 1)) + f * (log(x(i)) - log(x(i - 1))));
    if xk < xt, xt = xk; kWorst = k; end
end
if isinf(xt), xt = NaN; end
end

function th = interpDtheta(x, dth, xt)
th = [interp1(log(x), dth(1, :), log(xt)); interp1(log(x), dth(2, :), log(xt))];
end

function sgtitleSafe(fig, txt)
try
    sgtitle(fig, txt, 'FontSize', 10);
catch
    annotation(fig, 'textbox', [0 0.95 1 0.05], 'String', txt, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end