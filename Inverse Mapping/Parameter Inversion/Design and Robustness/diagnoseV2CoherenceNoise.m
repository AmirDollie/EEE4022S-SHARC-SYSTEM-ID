%% diagnoseV2CoherenceNoise.m
% V2 post-run DIAGNOSTIC (not a registered criterion): how well must the sensor noise PSD be known for the coherence
% early warning to mean anything?
%
% WHY. V2's detection statistic compares the observed coherence with the coherence PREDICTED by the single-direction
% model, and at p0 that prediction is set mostly by the noise floor (mean pair coherence ~0.99, so 1 - gamma^2 ~ 0.01
% is largely S_nn / S_sig). If the assumed S_nn is off by eps, the predicted baseline moves by about eps / 2 in Fisher z
% per bin, while the 80%-power detection threshold is only ~0.03 in bin-mean z. A first-order estimate puts the noise
% error that mimics (eps < 0) or masks (eps > 0) a detection at about 6%, whereas R3 only required +-30.6% for the
% inverse itself. This script measures it exactly on the run's own field and null calibration.
%
% WHAT. Truth: single-direction (or directional) sea with the nominal noise. Analyst's baseline: the single-direction
% model with S_nn,assumed = S_nn (1 + eps). D_j = mean_k(z_baseline - z_expected) / sd_j with the run's Monte Carlo
% sd_j (results.nullCoherence.sd), exactly as v2Tools.coherenceStats.
%   (a) no second direction: the false detection D_max(eps); the eps at which it reaches the 5% rejection threshold
%       (false alarms become likely) and the 80%-power threshold, per sign
%   (b) f_detect(eps) per angle for eps in cfg-free list below, against the run's f(d_sys = 0.5): does the early
%       warning at 60 and 90 deg survive a realistic noise-floor error?
% COST. Cached fields only (seconds). OUTPUT. Printed; Results/V2/diagnostics/V2noise_<stamp>.mat and
% V2_noise_falseAlarm.csv, V2_noise_fdetect.csv.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();
T = v2Tools(); V = cfg.V2;
Lst = dir(fullfile(P.results, 'V2', 'V2_20*.mat'));
if isempty(Lst), error('V2noise:none', 'Run runDirectionalForcingStudy first.'); end
[~, k] = max([Lst.datenum]); R = load(fullfile(Lst(k).folder, Lst(k).name)); res = R.results;
dDir = fullfile(P.results, 'V2', 'diagnostics'); if ~exist(dDir, 'dir'), mkdir(dDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS');
[~, gc] = system(sprintf('git -C "%s" rev-parse --short HEAD', P.root));
fprintf('\n==== V2 diagnostic: coherence early warning vs noise-PSD error (%s, code %s, results %s) ====\n', stamp, ...
    strtrim(gc), Lst(k).name);

out0 = evaluateScenario(struct('cacheDir', P.cache), cfg);
ctx = T.context(out0, [V.primaryHeadingDeg, V.angleDeg] * pi / 180, cfg, P.cache);
sd = res.nullCoherence.sd(:);
m = size(ctx.coh0, 1);
thr80 = T.normInv(1 - V.coherence.alpha / m) + T.normInv(V.coherence.power);
thr05 = T.normInv(1 - V.coherence.alpha / m);
D0 = directionalSpectralMatrix(ctx.src, T.components(ctx.pts, ctx.layout, 0, 1), ctx.Snn);
cohA = @(e) getfield(analyticTransmissibilityCovariance(D0.Ssig, [], ctx.Snn * (1 + e), ctx.ref, ctx.nEff), 'coherence');
cohDir = @(phi, f2) getfield(analyticTransmissibilityCovariance(directionalSpectralMatrix(ctx.src, ...
    T.twoComponents(ctx.pts, ctx.layout, phi, f2), ctx.Snn).Ssig, [], ctx.Snn, ctx.ref, ctx.nEff), 'coherence');
Dmax = @(cohBase, cohObs) max(getfield(T.coherenceStats(cohBase, cohObs, ctx.valid, ctx.nDof, V.coherence, sd), 'D'));

%% ---- (a) false detection with no second direction -----------------------------------------------------------------
epsA = [-0.30 -0.10 -0.05 -0.02 -0.01 0.01 0.02 0.05 0.10 0.30];
fprintf('\n(a) no second direction: D_max when the baseline assumes S_nn (1 + eps)  [5%% rejection %.3f, 80%% power %.3f]\n', ...
    thr05, thr80);
FA = zeros(size(epsA));
for i = 1:numel(epsA)
    FA(i) = Dmax(cohA(epsA(i)), ctx.coh0);
    fprintf('  eps %+5.0f%%  D_max %+8.3f  %s\n', 100 * epsA(i), FA(i), verdictA(FA(i), thr05, thr80));
end
fprintf('  noise error at which a FALSE detection appears (underestimated S_nn):\n');
epsCross = NaN(1, 2); thrs = [thr05 thr80];
for c = 1:2
    f = @(u) Dmax(cohA(-10^u), ctx.coh0);                      % eps = -10^u, u in log10 |eps|
    try
        [u, ~] = r1Boundary(f, thrs(c), [-4 log10(0.5)], 0.005);
        epsCross(c) = -10^u;
        fprintf('    D_max = %.3f at eps = %.2f%%\n', thrs(c), 100 * epsCross(c));
    catch
        fprintf('    D_max = %.3f not reached for |eps| in [1e-4, 0.5]\n', thrs(c));
    end
end

%% ---- (b) does the early warning survive? ---------------------------------------------------------------------------
epsB = [-0.10 -0.05 0 0.05 0.10];
fprintf('\n(b) f_detect per angle when the baseline assumes S_nn (1 + eps); f(d_sys = 0.5) from the run\n');
X = res.crossings;
FD = NaN(numel(V.angleDeg), numel(epsB)); FB = NaN(numel(V.angleDeg), 1);
for a = 1:numel(V.angleDeg)
    xb = X(abs([X.phiDeg] - V.angleDeg(a)) < 1e-6 & strcmp({X.criterion}, 'd_sys') & [X.threshold] == V.thresholds.dSys(1));
    FB(a) = xb.f2;
    line = sprintf('  phi %2g  f(d_sys = 0.5) %.3g  f_detect:', V.angleDeg(a), FB(a));
    for i = 1:numel(epsB)
        cb = cohA(epsB(i));
        met = @(u) Dmax(cb, cohDir(V.angleDeg(a) * pi / 180, 10^u)) - thr80;
        if met(-7) >= 0
            FD(a, i) = 0; s = 'false alarm';
        elseif met(log10(0.5)) < 0
            FD(a, i) = Inf; s = 'never';
        else
            u = r1Boundary(met, 0, [-7 log10(0.5)], 0.005); FD(a, i) = 10^u; s = sprintf('%.3g', FD(a, i));
        end
        if FD(a, i) > 0 && isfinite(FD(a, i)), s = [s, ternary(FD(a, i) < FB(a), ' (warns)', ' (late)')]; end
        line = [line, sprintf('  [%+3.0f%%] %s', 100 * epsB(i), s)]; %#ok<AGROW>
    end
    fprintf('%s\n', line);
end

save(fullfile(dDir, sprintf('V2noise_%s.mat', stamp)), 'epsA', 'FA', 'epsCross', 'epsB', 'FD', 'FB', 'thr05', 'thr80', ...
    'sd', '-v7');
run = struct('id', 'V2', 'dir', dDir);
writeStudyCsv(run, 'noise_falseAlarm', {'eps_pct', 'D_max'}, num2cell([100 * epsA(:), FA(:)]));
rowsB = {};
for a = 1:numel(V.angleDeg)
    for i = 1:numel(epsB)
        rowsB(end + 1, :) = {V.angleDeg(a), 100 * epsB(i), FD(a, i), FB(a)}; %#ok<SAGROW>
    end
end
writeStudyCsv(run, 'noise_fdetect', {'phi_deg', 'eps_pct', 'f_detect', 'f_dsys_0p5'}, rowsB);
fprintf('\n  saved %s\n', fullfile(dDir, sprintf('V2noise_%s.mat', stamp)));

function s = verdictA(D, t05, t80)
if D >= t80, s = 'FALSE DETECTION (>= 80%-power threshold)';
elseif D >= t05, s = 'false alarm likely (>= 5% rejection threshold)';
elseif D <= -t05, s = 'masks real drops';
else, s = '';
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end