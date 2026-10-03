%% runR1SNR.m
% STUDY R1 (MUST): at what signal-to-noise ratio does the inverse stop working, and where do the real
% sensors sit on that axis?
%
% ONE AXIS. The system is linear and the signal spectrum scales as Hs^2, so only Hs^2 / S_nn matters:
% halving the noise PSD is the same as raising Hs by sqrt(2). R1 sweeps S_nn at the fixed reference sea
% (JONSWAP Hs 0.05 m, omega_p 5 rad/s, gamma_J 3.3) and every boundary is also given as the Hs at which
% the LSM6DSV16X would reach it.
%
% WHY THERE CAN BE TWO BOUNDARIES. Lower noise raises F, but the Welch feature bias b_W does not depend on
% the noise, so d_W = sqrt(dth' F dth) grows as the noise falls. The Identifiable class can therefore be
% a WINDOW: limited by precision at the noisy end and by Welch bias at the quiet end (cfg.R1.boundaryRule:
% every transition is found, none is assumed).
%
% METHOD (pre-registered cfg.R1; clarification logged 2026-10-03 before the run)
%   1. 41 log-spaced S_nn from 1e-8 to 1e-4 through evaluateScenario (p0, Level 2C, reference s26,
%      N = 65536, L = 2048, cached field): sigma, rho, kappa, d_W, invalid fraction, min coherence, class.
%   2. Every named noise case (cfg.R1.marked) evaluated exactly, not read off the grid.
%   3. Every class change between neighbouring grid points is refined with r1Boundary on the criterion
%      that changes there (s_max against 1% / 5%, d_W against 0.5 / 1, kappa against 1e8), to
%      cfg.R1.refineTolDecades in log10(S_nn).
%   4. Identifiable interval(s); each boundary as S_nn, ASD (ug/sqrt(Hz)) and equivalent Hs for
%      cfg.R1.equivalentHsSensor; sensor / G3 ratios in PSD and ASD.
%
% RUN ORDER. testR1SNR must pass first. COST: about 41 + 6 + 10 per boundary evaluations, a minute or two.
% OUTPUT. Results/R1/ (log, .mat, summary, R1_sweep.csv, R1_sensor_points.csv, R1_boundaries.csv; figures via
% plotR1SNR once it exists).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();
run = startStudy('R1', 'SNR axis and measured noise', P.results);

try
Snn = cfg.R1.Snn(:).'; nS = numel(Snn); xg = log10(Snn);
names = cfg.R1.marked(:).'; nM = numel(names);
I = cfg.class.identifiable; Mg = cfg.class.marginal;
ev = @(s) evaluateScenario(struct('Snn', s, 'sea', cfg.R1.sea, 'cacheDir', P.cache), cfg);
asdUg = @(s) sqrt(2 * pi * s) / 9.81e-6;                 % one-sided per rad/s -> ug/sqrt(Hz)
fprintf('Step 1: %d S_nn from %.3g to %.3g (%.1f to %.0f ug/sqrt(Hz)); sea Hs %.3g m, omega_p %.3g, gamma_J %.3g\n', ...
    nS, Snn(1), Snn(end), asdUg(Snn(1)), asdUg(Snn(end)), cfg.R1.sea.Hs, cfg.R1.sea.omegaP, cfg.R1.sea.gammaJ);

%% ---- 1  sweep -------------------------------------------------------------------------------------------------
G = sweepMetrics(nS);
tic;
for k = 1:nS
    G = putMetrics(G, k, ev(Snn(k)), cfg);
end
fprintf('  sweep done in %.0f s\n', toc);
fprintf('  %10s %9s %9s %9s %7s %8s %8s  %s\n', 'S_nn', 'ug/rtHz', 'sig lnb', 'sig lnR', 'd_W', 'kappa', 'invalid', 'class');
for k = 1:nS
    fprintf('  %10.3g %9.1f %8.3f%% %8.3f%% %7.3f %8.3g %7.1f%%  %s\n', Snn(k), asdUg(Snn(k)), 100 * G.sB(k), ...
        100 * G.sR(k), G.dW(k), G.kappa(k), 100 * G.invalid(k), G.class{k});
end

%% ---- 2  named noise cases -------------------------------------------------------------------------------------
Mk = sweepMetrics(nM); SnnM = cellfun(@(n) cfg.noise.(n), names);
for j = 1:nM
    Mk = putMetrics(Mk, j, ev(SnnM(j)), cfg);
end
SnnG3 = cfg.noise.g3;

%% ---- 3  transitions -------------------------------------------------------------------------------------------
crit = struct('name', {'s_max', 's_max', 'd_W', 'd_W', 'kappa'}, 'field', {'sMax', 'sMax', 'dW', 'dW', 'kappa'}, ...
    'threshold', {I.maxSigma, Mg.maxSigma, I.maxBiasDistance, Mg.maxBiasDistance, I.maxCondF});
B = struct('x', {}, 'Snn', {}, 'from', {}, 'to', {}, 'criterion', {}, 'threshold', {}, 'nEval', {}, 'bracket', {});
cache = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k = find(G.code(1:end - 1) ~= G.code(2:end))
    found = false;
    for c = 1:numel(crit)
        v = G.(crit(c).field)([k k + 1]); t = crit(c).threshold;
        if sign(v(1) - t) ~= sign(v(2) - t)
            f = @(x) metricAt(x, crit(c).field, ev, cache, cfg);
            [xb, rinfo] = r1Boundary(f, t, xg([k k + 1]), cfg.R1.refineTolDecades);
            B(end + 1) = struct('x', xb, 'Snn', 10^xb, 'from', G.class{k}, 'to', G.class{k + 1}, ...
                'criterion', crit(c).name, 'threshold', t, 'nEval', rinfo.nEval, 'bracket', Snn([k k + 1])); %#ok<SAGROW>
            found = true;
        end
    end
    if ~found
        B(end + 1) = struct('x', mean(xg([k k + 1])), 'Snn', 10^mean(xg([k k + 1])), 'from', G.class{k}, ...
            'to', G.class{k + 1}, 'criterion', 'undetermined (no single criterion crosses)', 'threshold', NaN, ...
            'nEval', 0, 'bracket', Snn([k k + 1])); %#ok<SAGROW>
    end
end
SnnRef = cfg.noise.(cfg.R1.equivalentHsSensor);
HsEq = @(s) cfg.R1.sea.Hs * sqrt(SnnRef / s);         % Hs at which the reference sensor has this SNR
fprintf('\nStep 3: %d class transition(s) on [%.3g, %.3g]\n', numel(B), Snn(1), Snn(end));
for b = 1:numel(B)
    fprintf('  S_nn = %.4g (%.1f ug/sqrt(Hz), Hs_eq %s = %.4f m): %s -> %s, binding %s = %.4g  (%d evaluations)\n', ...
        B(b).Snn, asdUg(B(b).Snn), cfg.R1.equivalentHsSensor, HsEq(B(b).Snn), B(b).from, B(b).to, B(b).criterion, ...
        B(b).threshold, B(b).nEval);
end
% Identifiable interval(s): runs of class 1 on the grid, ends replaced by refined boundaries
runs = {}; inRun = false;
for k = 1:nS
    if G.code(k) == 1 && ~inRun, lo = k; inRun = true; end
    if inRun && (k == nS || G.code(k + 1) ~= 1)
        runs{end + 1} = [lo k]; inRun = false; %#ok<SAGROW>
    end
end
intervals = zeros(numel(runs), 2);
for r = 1:numel(runs)
    a = runs{r}(1); z = runs{r}(2);
    intervals(r, :) = [edgeOf(B, Snn, a, 'lo'), edgeOf(B, Snn, z, 'hi')];
end

%% ---- 4  tables ---------------------------------------------------------------------------------------------------
rows = arrayfun(@(k) {Snn(k), asdUg(Snn(k)), HsEq(Snn(k)), G.sB(k), G.sR(k), G.rho(k), G.kappa(k), G.dW(k), ...
    G.invalid(k), G.nUsed(k), G.minCoh(k), G.class{k}}, 1:nS, 'UniformOutput', false);
writeStudyCsv(run, 'sweep', {'Snn_per_rad_s', 'asd_ug_rtHz', 'Hs_equivalent_m', 'sigma_lnbeta', 'sigma_lnR', 'rho_betaR', ...
    'condF', 'dW', 'invalid_fraction', 'bins_used', 'min_coherence', 'class'}, vertcat(rows{:}));
rows = arrayfun(@(j) {names{j}, cfg.noise.labels{strcmp(cfg.noise.names, names{j})}, SnnM(j), asdUg(SnnM(j)), ...
    SnnM(j) / SnnG3, sqrt(SnnM(j) / SnnG3), HsEq(SnnM(j)), Mk.sB(j), Mk.sR(j), Mk.dW(j), Mk.kappa(j), Mk.invalid(j), ...
    Mk.class{j}}, 1:nM, 'UniformOutput', false);
writeStudyCsv(run, 'sensor_points', {'case', 'label', 'Snn_per_rad_s', 'asd_ug_rtHz', 'psd_over_G3', 'asd_over_G3', ...
    'Hs_equivalent_m', 'sigma_lnbeta', 'sigma_lnR', 'dW', 'condF', 'invalid_fraction', 'class'}, vertcat(rows{:}));
if ~isempty(B)
    rows = arrayfun(@(b) {B(b).Snn, asdUg(B(b).Snn), HsEq(B(b).Snn), B(b).from, B(b).to, B(b).criterion, B(b).threshold, ...
        B(b).bracket(1), B(b).bracket(2)}, 1:numel(B), 'UniformOutput', false);
    writeStudyCsv(run, 'boundaries', {'Snn_per_rad_s', 'asd_ug_rtHz', 'Hs_equivalent_m', 'from_class', 'to_class', ...
        'criterion', 'threshold', 'bracket_lo', 'bracket_hi'}, vertcat(rows{:}));
end

%% ---- 5  summary --------------------------------------------------------------------------------------------------
L = {};
L{end + 1} = sprintf(['p0, Level 2C (reference s26), N = %d, L = %d, JONSWAP Hs %.3g m, omega_p %.3g rad/s; S_nn swept ' ...
    '%.3g to %.3g per rad/s (%.1f to %.0f ug/sqrt(Hz)). SNR ~ Hs^2 / S_nn.'], cfg.acq.Nref, cfg.welch.L, cfg.R1.sea.Hs, ...
    cfg.R1.sea.omegaP, Snn(1), Snn(end), asdUg(Snn(1)), asdUg(Snn(end)));
L{end + 1} = '';
if isempty(intervals)
    L{end + 1} = 'No Identifiable interval on the swept range.';
end
for r = 1:size(intervals, 1)
    L{end + 1} = sprintf(['Identifiable for S_nn in [%.4g, %.4g] (%.1f to %.1f ug/sqrt(Hz)); for the %s at this sea, ' ...
        'Hs from %.4f to %.4f m.'], intervals(r, :), asdUg(intervals(r, 1)), asdUg(intervals(r, 2)), ...
        cfg.R1.equivalentHsSensor, HsEq(intervals(r, 2)), HsEq(intervals(r, 1))); %#ok<SAGROW>
end
L{end + 1} = 'Class transitions (low to high noise):';
for b = 1:numel(B)
    L{end + 1} = sprintf('  S_nn = %.4g (%.1f ug/sqrt(Hz)): %s -> %s, binding %s = %.4g; Hs_eq %.4f m', B(b).Snn, ...
        asdUg(B(b).Snn), B(b).from, B(b).to, B(b).criterion, B(b).threshold, HsEq(B(b).Snn)); %#ok<SAGROW>
end
L{end + 1} = '';
L{end + 1} = sprintf('%-34s %10s %9s %8s %8s %9s %9s %7s  %s', 'measured / assumed noise', 'S_nn', 'ug/rtHz', ...
    'PSD/G3', 'ASD/G3', 'sig lnb', 'sig lnR', 'd_W', 'class');
for j = 1:nM
    L{end + 1} = sprintf('%-34s %10.3g %9.1f %8.3f %8.3f %8.3f%% %8.3f%% %7.3f  %s', ...
        cfg.noise.labels{strcmp(cfg.noise.names, names{j})}, SnnM(j), asdUg(SnnM(j)), SnnM(j) / SnnG3, ...
        sqrt(SnnM(j) / SnnG3), 100 * Mk.sB(j), 100 * Mk.sR(j), Mk.dW(j), Mk.class{j}); %#ok<SAGROW>
end

results = struct('Snn', Snn, 'grid', G, 'names', {names}, 'SnnMarked', SnnM, 'marked', Mk, 'boundaries', B, ...
    'intervals', intervals, 'SnnG3', SnnG3, 'equivalentHsSensor', cfg.R1.equivalentHsSensor, 'SnnRef', SnnRef);
if exist('plotR1SNR', 'file') == 2
    try
        plotR1SNR(results, cfg, run);
    catch perr
        fprintf('  plotting failed (results are still saved): %s\n', perr.message);
    end
end
finishStudy(run, results, cfg, L);
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================
function G = sweepMetrics(n)
G = struct('sB', NaN(1, n), 'sR', NaN(1, n), 'rho', NaN(1, n), 'kappa', NaN(1, n), 'dW', NaN(1, n), ...
    'sMax', NaN(1, n), 'invalid', NaN(1, n), 'nUsed', NaN(1, n), 'minCoh', NaN(1, n), 'code', NaN(1, n));
G.class = cell(1, n);
end

function G = putMetrics(G, k, out, cfg)
G.sB(k) = out.sigma(1); G.sR(k) = out.sigma(2); G.rho(k) = out.corr; G.kappa(k) = out.condF; G.dW(k) = out.dW;
G.sMax(k) = out.sMax; G.invalid(k) = out.invalidFraction; G.nUsed(k) = out.nUsedBins; G.minCoh(k) = out.minCoherence;
G.class{k} = out.class; G.code(k) = find(strcmp(cfg.class.names, out.class));
end

function v = metricAt(x, field, ev, cache, cfg)
% one criterion at S_nn = 10^x; evaluations are cached so two criteria in one bracket share them
if isKey(cache, x)
    G = cache(x);
else
    G = putMetrics(sweepMetrics(1), 1, ev(10^x), cfg);
    cache(x) = G; %#ok<NASGU>   containers.Map is a handle: the entry persists
end
v = G.(field);
end

function e = edgeOf(B, Snn, k, side)
% refined boundary at the edge of an Identifiable run, or the sweep end if the run reaches it
if strcmp(side, 'lo'), e = Snn(1); else, e = Snn(end); end
if strcmp(side, 'lo') && k > 1
    b = find(arrayfun(@(q) B(q).bracket(2) == Snn(k), 1:numel(B)), 1);
    if ~isempty(b), e = B(b).Snn; end
elseif strcmp(side, 'hi') && k < numel(Snn)
    b = find(arrayfun(@(q) B(q).bracket(1) == Snn(k), 1:numel(B)), 1);
    if ~isempty(b), e = B(b).Snn; end
end
end