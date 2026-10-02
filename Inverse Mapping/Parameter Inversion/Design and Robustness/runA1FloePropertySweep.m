%% runA1FloePropertySweep.m
% STUDY A1 (MUST): for which floes does the transmissibility inverse work?
%
% QUESTION. Under ONE fixed sensing arrangement (Level 2C fractions of each floe's true radius,
% reference s26, LSM6DSV16X noise, setting F), which floes are identifiable? The sensing is NOT
% re-optimised per floe, so a change in performance comes from the floe, not from the layout (D3 is
% where the layout is revisited).
%
% SCENARIOS (pre-registered, cfg.A1): 8 ring cases (beta x 0.25..4 at R0; R x 0.5..2 at beta0), 4 corners
% ({0.5, 2}^2), 12 physical classes (basin p0; four polypropylene discs; the SHARC model floe at three
% provisional moduli; four MIZ floes). 24 in total.
%
% METHOD, per scenario (a1Tools.evaluate, the same code the tester checks)
%   1. Physical -> EMM scaling with ONE depth throughout. Finite depth: the basin, H = 1.88 m, band
%      cfg.band.basin. MIZ: effective depth from screenFloeRegime, evaluated Froude-scaled to the basin
%      (AMENDMENT 2026-10-02, cfg.meta.changes); the equivalent full-scale experiment is reported.
%   2. Production metrics from the frozen E1 chain (evaluateScenario, unchanged): sigma, rho, kappa(F),
%      Welch bias d_W, inverse class.
%   3. Forward-model and numerical checks:
%        node grid   0.2 vs 0.1 rad/s: feature discrepancy projected with the scenario's own F, d < 0.1;
%                    on failure the scenario is re-evaluated on 0.1 rad/s
%        truncation  [50 10 10] vs [70 15 15] at the grid actually used: convergence of the A1 METRICS,
%                    sigma within 1%, d_W within 0.02 and the same class (amended 2 Oct after the p0
%                    diagnostic); the feature-bias d is retained for V1, not gated
%        alpha       retained bins inside the validated EMM window [1.7, 13.85]
%      Any failure: final class 'Model check failed' (the inverse class is still recorded). Singular F:
%      checks not assessable, final class Not identifiable (modelAssessed = false).
%   4. Dimensionless groups: kR (band minimum, JONSWAP peak, maximum) and R / l_f.
%
% RESUMABLE. Each scenario is saved to Results/A1/checkpoints/ as soon as it finishes, with a key of
% everything it depends on; a later run loads a matching checkpoint instead of recomputing, and recomputes
% a stale one. A scenario that throws is logged and skipped (not checkpointed, so a rerun retries it);
% the batch continues.
%
% RUN ORDER. testA1FloePropertySweep must pass first (it also warms the EMM cache for p0). Then
%   PREFLIGHT_ONLY = true   prints the 24 scenarios, bands, groups and the time estimate; no EMM
%   ONLY = {'modelFloe_E3.0GPa', 'miz_R25_h1'}  real-EMM smoke test of the nu = 0.4 and Froude-scaled MIZ
%                           routes, which the tester exercises only with analytic fields
%   then the full batch overnight (ONLY = {}).
%
% COST. Per scenario 5 node-FRF sets at 0.2 rad/s (29 nodes, ~30 s each), one at 0.1 rad/s (57 nodes) and
% 5 at [70 15 15] (~40 s each): roughly 7 min, so about 3 h for the batch; the EMM cache
% (Results/twinCache) is reused across runs.
%
% OUTPUT. Results/A1/ (log, .mat, summary, A1_scenarios.csv, A1_validity.csv, A1_equivalent.csv;
% figures via plotA1FloePropertySweep once it exists).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

PREFLIGHT_ONLY = false;   % true: scenario table, bands and groups only (no EMM solve)
ONLY = {'modelFloe_E3.0GPa', 'miz_R25_h1'};                % {} = all 24; or a cell of scenario names for a trial run
FORCE = false;            % true: ignore existing checkpoints and recompute
VERBOSE = true;           % one line per node-FRF set inside each scenario

P = setupStudyPaths();
cfg = defineStudyScenarios();
A = a1Tools();
run = startStudy('A1', 'Floe-property and domain sweep', P.results);
ckDir = fullfile(run.dir, 'checkpoints');

try
%% ---- 1  scenarios ------------------------------------------------------------------------------------
sc = A.buildScenarios(cfg);
groups = {sc.group};
fprintf('Step 1: %d scenarios: %d ring, %d corner, %d physical (%d basin, %d tank, %d MIZ)\n', numel(sc), ...
    nnz(strcmp(groups, 'ring')), nnz(strcmp(groups, 'corner')), nnz(~ismember(groups, {'ring', 'corner'})), ...
    nnz(strcmp(groups, 'basin')), nnz(strcmp(groups, 'tank')), nnz(strcmp(groups, 'miz')));
fprintf('  sensing: Level 2C fractions of each true R, reference s26, noise %s (%.3g), N = %d, L = %d\n', ...
    cfg.A1.noiseCase, cfg.noise.(cfg.A1.noiseCase), cfg.acq.Nref, cfg.welch.L);
fprintf(['  checks: node grid %.1f vs %.1f rad/s (feature d < %.2g); truncation %s vs %s (sigma within %.0f%%, ' ...
    'd_W within %.2g, same class; feature bias reported, not gated)\n\n'], cfg.emm.nodeSpacing, cfg.emm.nodeSpacingCheck, ...
    cfg.emm.maxInterpBiasDistance, mat2str(cfg.emm.truncation), mat2str(cfg.emm.truncationCheck), ...
    100 * cfg.A1.truncationGate.maxRelSigma, cfg.A1.truncationGate.maxDeltaDW);

%% ---- 2  preflight table (no EMM) ---------------------------------------------------------------------
[omegaAll] = deal(cell(1, numel(sc)));
fprintf('%-3s %-26s %-7s %9s %9s %9s %7s %15s %14s %7s %7s %6s\n', '#', 'scenario', 'group', 'beta/b0', ...
    'R/R0', 'gamma/g0', 'nu', 'band (rad/s)', 'periods (s)', 'kR pk', 'R/l_f', 'scale');
for i = 1:numel(sc)
    w = welchBins(cfg.welch.L, cfg.acq.dt, cfg.band.estimate, cfg.welch.binStep, sc(i).band);
    omegaAll{i} = w;
    G = A.dimensionlessGroups(sc(i), w, cfg);
    prov = ''; if sc(i).provisional, prov = ' (prov.)'; end
    fprintf('%-3d %-26s %-7s %9.3g %9.3g %9.3g %7.2f %7.3f-%-7.3f %6.2f-%-7.2f %7.3f %7.3f %6.1f%s\n', i, sc(i).name, ...
        sc(i).group, sc(i).betaRatio, sc(i).RRatio, sc(i).gammaRatio, sc(i).nu, sc(i).band, sc(i).bandPeriods, ...
        G.kRpeak, G.RoverLf, sc(i).scale, prov);
    if ~G.alphaValid
        fprintf('     ! alpha [%.2f %.2f] outside the validated window\n', G.alphaMin, G.alphaMax);
    end
end
fprintf('\n');
if PREFLIGHT_ONLY
    fprintf('PREFLIGHT_ONLY: stopping before any EMM solve.\n');
    diary('off');
    return
end

%% ---- 3  the batch ------------------------------------------------------------------------------------
sel = 1:numel(sc);
if ~isempty(ONLY)
    sel = find(ismember({sc.name}, ONLY));
    missing = setdiff(ONLY, {sc.name});
    if ~isempty(missing), error('runA1:only', 'Unknown scenario name(s): %s', strjoin(missing, ', ')); end
    fprintf('ONLY: %s\n\n', strjoin({sc(sel).name}, ', '));
end
opts = struct('cacheDir', P.cache, 'verbose', VERBOSE);
res = cell(1, numel(sc)); status = repmat({'not run'}, 1, numel(sc)); msg = repmat({''}, 1, numel(sc));
tBatch = tic;
for n = 1:numel(sel)
    i = sel(n);
    key = A.checkpointKey(sc(i), cfg);
    file = A.checkpointPath(ckDir, sc(i));
    fprintf('[%02d/%02d] %s  (%s)\n', n, numel(sel), sc(i).name, datestr(now, 'HH:MM:SS'));
    if ~FORCE
        [r, ok] = A.loadCheckpoint(file, key);
        if ok
            res{i} = r; status{i} = 'loaded';
            fprintf('  loaded checkpoint: %s\n', classLine(r));
            continue
        elseif exist(file, 'file')
            fprintf('  checkpoint is stale (configuration changed): recomputing\n');
        end
    end
    try
        r = A.evaluate(sc(i), cfg, opts);
        A.saveCheckpoint(file, r, key);
        res{i} = r; status{i} = 'computed';
        if ~r.modelAssessed, fprintf('  model checks ... NOT ASSESSABLE (F singular)\n'); end
        if r.nodeGrid.done
            fprintf('  node grid ...... %s (d = %.3g, max |dT|/|T| = %.2e%s)\n', passStr(r.nodeGridValid), ...
                r.nodeGrid.d, r.nodeGrid.relT, ternary(r.nodeGrid.fallback, ', FALLBACK to 0.1 rad/s', ''));
        else
            fprintf('  node grid ...... not run\n');
        end
        if r.truncation.done
            fprintf('  truncation ..... %s (sigma %.1e rel, d_W %.1e abs, class %s; feature bias d = %.3g, max |dT|/|T| = %.2e)\n', ...
                passStr(r.truncationValid), r.truncation.relSigma, r.truncation.deltaDW, r.truncation.classHi, ...
                r.truncation.d, r.truncation.relT);
        else
            fprintf('  truncation ..... not run\n');
        end
        fprintf('  alpha .......... %s ([%.2f %.2f])\n', passStr(r.alphaValid), r.groups.alphaMin, r.groups.alphaMax);
        fprintf('  E1 metrics ..... %s\n', classLine(r));
        fprintf('  saved (%.0f s; batch %.1f min)\n', r.seconds, toc(tBatch) / 60);
    catch err
        status{i} = 'FAILED'; msg{i} = err.message;
        fprintf('  FAILED: %s\n', err.message);
        if ~isempty(err.stack), fprintf('    at %s line %d\n', err.stack(1).name, err.stack(1).line); end
    end
end
fprintf('\nBatch finished in %.1f min: %d computed, %d loaded, %d failed, %d not run\n\n', toc(tBatch) / 60, ...
    nnz(strcmp(status, 'computed')), nnz(strcmp(status, 'loaded')), nnz(strcmp(status, 'FAILED')), ...
    nnz(strcmp(status, 'not run')));

%% ---- 4  gather ---------------------------------------------------------------------------------------
done = find(~cellfun(@isempty, res));
tab = struct('name', {sc.name}, 'group', {sc.group}, 'status', status, 'message', msg);
for i = 1:numel(sc)
    r = res{i};
    if isempty(r)
        tab(i).finalClass = status{i}; tab(i).finalCode = NaN; tab(i).inverseCode = NaN;
        continue
    end
    tab(i).finalClass = r.finalClass; tab(i).finalCode = r.finalCode; tab(i).inverseCode = r.inverseCode;
end

%% ---- 5  tables -----------------------------------------------------------------------------------------
rows = {}; vrows = {}; erows = {};
for i = done
    r = res{i}; s = sc(i); G = r.groups;
    rows(end + 1, :) = {s.name, s.group, s.provisional, s.betaRatio, s.RRatio, s.gammaRatio, s.p(1), s.p(2), s.p(3), ...
        s.Rphys, s.h, s.E, s.nu, s.Hphys, s.lf, G.kRmin, G.kRpeak, G.kRmax, G.RoverLf, ...
        r.sigmaBeta, r.sigmaR, r.rho, r.condF, r.dW, r.nUsedBins, r.minCoherence, r.inverseClass, r.modelAssessed, ...
        r.modelValid, r.finalClass}; %#ok<AGROW>
    vrows(end + 1, :) = {s.name, s.band(1), s.band(2), G.alphaMin, G.alphaMax, r.alphaValid, r.nodeGrid.d, ...
        r.nodeGrid.relT, r.nodeGrid.fallback, r.nodeSpacingUsed, r.nodeGridValid, r.truncation.relSigma, r.truncation.deltaDW, ...
        r.truncation.classHi, r.truncation.d, r.truncation.relT, ...
        r.truncationValid, r.modelAssessed, r.modelValid, r.modelReason, r.seconds}; %#ok<AGROW>
    e = r.equivalent;
    erows(end + 1, :) = {s.name, e.scale, e.timeScale, e.Hs, e.Tp, e.recordMin, e.dt, e.Snn, e.periods(1), ...
        e.periods(2), e.note}; %#ok<AGROW>
end
if ~isempty(done)
    writeStudyCsv(run, 'scenarios', {'scenario', 'group', 'provisional', 'beta_over_beta0', 'R_over_R0', ...
        'gamma_over_gamma0', 'beta', 'gamma', 'R', 'Rphys_m', 'h_m', 'E_Pa', 'nu', 'depth_m', 'lf_m', 'kR_min', ...
        'kR_peak', 'kR_max', 'R_over_lf', 'sigma_lnbeta', 'sigma_lnR', 'rho_betaR', 'condF', 'dW', 'bins_used', ...
        'min_coherence', 'inverse_class', 'model_assessed', 'model_valid', 'final_class'}, rows);
    writeStudyCsv(run, 'validity', {'scenario', 'band_lo_rad_s', 'band_hi_rad_s', 'alpha_min', 'alpha_max', ...
        'alpha_valid', 'd_node_grid', 'relT_node_grid', 'grid_fallback', 'node_spacing_used', 'node_grid_valid', ...
        'trunc_rel_sigma', 'trunc_delta_dW', 'trunc_class', 'trunc_feature_bias_d', 'trunc_relT', 'truncation_valid', 'model_assessed', 'model_valid', 'reason', 'seconds'}, vrows);
    writeStudyCsv(run, 'equivalent', {'scenario', 'length_scale', 'time_scale', 'Hs_m', 'Tp_s', 'record_min', ...
        'dt_s', 'Snn_per_rad_s', 'period_lo_s', 'period_hi_s', 'note'}, erows);
end

%% ---- 6  summary ------------------------------------------------------------------------------------------
L = {};
L{end + 1} = sprintf(['Sensing fixed: Level 2C fractions of each true R, reference s26, %s noise, N = %d, L = %d. ' ...
    'Model checks: node grid feature d < %.2g; truncation metrics converged (sigma %.0f%%, d_W %.2g, class); ' ...
    'alpha in [%.2f, %.2f]. Truncation feature bias reported for V1, not gated.'], cfg.A1.noiseCase, cfg.acq.Nref, ...
    cfg.welch.L, cfg.emm.maxInterpBiasDistance, 100 * cfg.A1.truncationGate.maxRelSigma, cfg.A1.truncationGate.maxDeltaDW, ...
    cfg.emm.alphaValidated);
L{end + 1} = sprintf('Scenarios with results: %d of %d (failed: %d).', numel(done), numel(sc), nnz(strcmp(status, 'FAILED')));
names = [cfg.class.names, {cfg.A1.modelCheckClass}];
codes = [tab.finalCode];
L{end + 1} = ['Final classes: ' strjoin(arrayfun(@(c) sprintf('%s %d', names{c}, nnz(codes == c)), 1:4, ...
    'UniformOutput', false), ', ') '.'];
L{end + 1} = '';
L{end + 1} = sprintf('%-26s %-7s %8s %8s %8s %8s %7s %6s %8s %8s  %s', 'scenario', 'group', 'b/b0', 'R/R0', ...
    'sig lnb', 'sig lnR', 'd_W', 'kR pk', 'R/l_f', 'kappa', 'final class (inverse class if different)');
for i = 1:numel(sc)
    r = res{i};
    if isempty(r)
        L{end + 1} = sprintf('%-26s %-7s  %s %s', sc(i).name, sc(i).group, status{i}, msg{i}); %#ok<SAGROW>
        continue
    end
    extra = '';
    if ~r.modelValid && r.modelAssessed, extra = sprintf(' (inverse: %s; %s)', r.inverseClass, r.modelReason); end
    if ~r.modelAssessed, extra = sprintf(' (%s)', r.modelReason); end
    if sc(i).provisional, extra = [extra ' [provisional]']; end %#ok<AGROW>
    L{end + 1} = sprintf('%-26s %-7s %8.3g %8.3g %7.3f%% %7.3f%% %7.3f %6.2f %8.3f %8.2g  %s%s', sc(i).name, ...
        sc(i).group, sc(i).betaRatio, sc(i).RRatio, 100 * r.sigmaBeta, 100 * r.sigmaR, r.dW, r.groups.kRpeak, ...
        r.groups.RoverLf, r.condF, r.finalClass, extra); %#ok<SAGROW>
end
mf = find(~cellfun(@isempty, res) & strncmp({sc.name}, 'modelFloe', 9));
if ~isempty(mf)
    L{end + 1} = '';
    L{end + 1} = ['Recorded hypothesis (model floe near rigid in the basin band; beta weak, R may remain): ' ...
        strjoin(arrayfun(@(i) sprintf('%s sigma_lnbeta / sigma_lnR = %.3g', sc(i).name, ...
        res{i}.sigmaBeta / res{i}.sigmaR), mf, 'UniformOutput', false), '; ') '.'];
end
mz = find(~cellfun(@isempty, res) & strcmp({sc.group}, 'miz'));
if ~isempty(mz)
    L{end + 1} = '';
    L{end + 1} = 'MIZ cases are Froude-scaled basin-equivalent evaluations (regime and domain screening). Equivalent full scale:';
    for i = mz
        e = res{i}.equivalent;
        L{end + 1} = sprintf('  %-22s depth %.0f m, Hs %.2f m, Tp %.1f s, record %.0f min, dt %.2f s, periods %.1f-%.1f s', ...
            sc(i).name, sc(i).Hphys, e.Hs, e.Tp, e.recordMin, e.dt, e.periods); %#ok<SAGROW>
    end
end

results = struct('scenarios', sc, 'res', {res}, 'table', tab, 'status', {status}, 'messages', {msg}, ...
    'omega', {omegaAll}, 'only', {ONLY}, 'tools', A.version);
if exist('plotA1FloePropertySweep', 'file') == 2 && ~isempty(done)
    try
        plotA1FloePropertySweep(results, cfg, run);
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
function s = classLine(r)
s = sprintf('sigma = [%.3f%% %.3f%%], d_W = %.3f, kappa = %.3g, %s', 100 * r.sigmaBeta, 100 * r.sigmaR, r.dW, ...
    r.condF, r.finalClass);
if ~r.modelAssessed, s = sprintf('%s (%s)', s, r.modelReason);
elseif ~r.modelValid, s = sprintf('%s (inverse: %s)', s, r.inverseClass); end
end

function s = passStr(tf)
if tf, s = 'PASS'; else, s = 'FAIL'; end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end