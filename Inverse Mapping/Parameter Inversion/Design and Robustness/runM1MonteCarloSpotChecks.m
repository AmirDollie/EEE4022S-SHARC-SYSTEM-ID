%% runM1MonteCarloSpotChecks.m
% M1: Monte Carlo and nonlinear spot checks of the ANALYTIC (E1) predictions away from p0.
%
% QUESTION. D1, A1, R1 and R2 rest on analytic F, sigma and Welch bias rather than on simulated records. G3
% showed at p0 (and A, B) that real finite records follow those predictions. Do they still, in the four
% deliberately awkward cases fixed by the pre-registered rules (cfg.M1)?
%   1 D1_ns2            the minimal two-sensor design (s6 + reference s18): one transmissibility per bin
%   2 A1_miz_R25_h1     the A1 case nearest the Identifiable/Marginal boundary (Froude-scaled MIZ floe)
%   3 R2_wp3.5_gJ7      the R2 stress sea (narrow, at the lower band edge)
%   4 R1_p0_lsm6dsv16x  p0 at the measured LSM6DSV16X noise (the measured-noise twin rerun)
% Truth and inverse share the EMM: M1 tests the estimator and the Fisher/bias machinery, not the physics (V1).
%
% METHOD (pre-registered: programme rev. 5 and cfg.M1; clarification logged 2026-10-03 before the run)
%   Stage 1  resolve the cases (resolveM1Scenarios) and evaluate each with the frozen engine; each must
%            reproduce its source study (sigma, d_W, class; cfg.M1.sourceRelTol). Any mismatch STOPS the
%            run before a record is generated.
%   Stage 2  cfg.M1.nRecords records per case (synthesiseTwinRecords on the engine's field twin, seeds
%            cfg.M1.seedBase + cfg.M1.seedStride * case + record), the G3 estimator (g3EstimatorCore,
%            verbatim), the linearised estimate with the analytic J, Sigma_k and F. Scatter: emp/pred in
%            cfg.M1.accept.empPred. Mean offset against the Welch-expected bias: diagnostic (G3.3 rule).
%   Stage 3  coverage of q = dtheta' F dtheta (G3a, centred on the truth) at 68 and 95%; accept the 95%
%            coverage in cfg.M1.accept.coverage95. If a statistic's 1-SE interval contains an acceptance
%            limit (borderline), records 101-200 are added once and the case is judged on all 200.
%   Stage 4  record 1 inverted with the G3 safeguarded Gauss-Newton from beta x 0.8 and x 1.2 (R true):
%            converged and |theta_NL - theta_lin| / sigma < 0.2 per component (G3.1).
%   A case passes if Stages 2, 3 and 4 all pass. A failure is reported, not tuned away.
%
% RESUMABLE. Each case has a checkpoint in Results/M1/checkpoints/ (records every cfg.M1.checkpointEvery and
% every nonlinear run). Rerunning the script continues from it; a checkpoint whose scenario, seeds or
% analytic target differ is ignored. FORCE = true starts every selected case afresh.
%
% PROVENANCE. The run refuses to start with uncommitted .m files in the repository (cfg.M1.requireCleanTree),
% so the code commit in the log is the code that ran. ALLOW_DIRTY = true overrides it (logged loudly).
%
% RUN ORDER. testM1MonteCarloSpotChecks must pass first. COST: Stage 1 seconds (cached fields); Stages 2-3
% about 1 s per record (400 to 800 records); Stage 4 about 300-600 EMM solves per inversion (G3), so roughly
% 10-15 min each and 1.5-2 h for all eight. The log prints a projection after the first record and the first
% inversion.
%
% OUTPUT. Results/M1/ (log, .mat with every per-record estimate, summary, M1_analytic.csv, M1_summary.csv,
% M1_records.csv, M1_nonlinear.csv; checkpoints/). Figures: plotM1MonteCarloSpotChecks (after the run).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

ONLY = {};              % e.g. {'R1_p0_lsm6dsv16x'} to run some cases; {} = all four
FORCE = false;          % true: ignore checkpoints of the selected cases
ALLOW_DIRTY = false;    % true: run with uncommitted .m files (provenance is then NOT the logged commit)

P = setupStudyPaths();
cfg = defineStudyScenarios();
T = m1Tools();
run = startStudy('M1', 'Monte Carlo and nonlinear spot checks', P.results);
ckDir = fullfile(run.dir, 'checkpoints');

try
%% ---- 0  provenance -------------------------------------------------------------------------------------------------
[st0, gitOut] = system(sprintf('git -C "%s" status --porcelain', P.root));
if st0 ~= 0
    fprintf('  WARNING: git status failed (%s); the code state cannot be checked.\n', strtrim(gitOut));
else
    [clean, dirty] = T.cleanTree(gitOut);
    if ~clean
        fprintf('  uncommitted .m files:\n'); fprintf('    %s\n', dirty{:});
        if cfg.M1.requireCleanTree && ~ALLOW_DIRTY
            error('M1:dirtyTree', ['%d uncommitted .m file(s): commit first so that the logged commit %s is the ' ...
                'code that runs (or set ALLOW_DIRTY = true).'], numel(dirty), run.commit);
        end
        fprintf('  *** ALLOW_DIRTY: this run is NOT reproducible from commit %s ***\n', run.commit);
    else
        fprintf('  working tree clean (.m files): the code that runs is commit %s\n', run.commit);
    end
end

%% ---- 1  scenarios and analytic targets ---------------------------------------------------------------------------
cases = resolveM1Scenarios(cfg, P);
sel = 1:numel(cases);
if ~isempty(ONLY), sel = find(ismember({cases.name}, ONLY)); end
if isempty(sel), error('M1:only', 'ONLY matches no case.'); end
fprintf('\nStage 1: %d case(s); analytic targets from the frozen engine, checked against their source studies\n', numel(sel));
ctx = cell(1, numel(cases)); cmp = cell(1, numel(cases));
allOk = true;
for q = sel
    fprintf('  %d %-18s %s\n', q, cases(q).name, cases(q).note);
    ctx{q} = T.context(cases(q), cfg);
    t = ctx{q}.tgt; s = cases(q).source;
    cmp{q} = T.compareSource(t, s, cfg);
    fprintf(['      analytic sigma [%.4f%% %.4f%%], d_W %.4f, kappa %.3g, %s, %d of %d bins (%.0f s)\n' ...
        '      source %s: sigma [%.4f%% %.4f%%], d_W %.4f, %s -> rel. sigma %.1e, d_W %.1e  %s\n'], ...
        100 * t.sigma, t.dW, t.condF, t.class, t.nUsedBins, t.nF, t.seconds, cases(q).study, 100 * s.sigma, s.dW, ...
        s.class, cmp{q}.relSigma, cmp{q}.dDW, okStr(cmp{q}.ok));
    allOk = allOk && cmp{q}.ok;
end
if ~allOk
    error('M1:target', 'An analytic target does not reproduce its source study: no records generated.');
end
fprintf('  all targets reproduce their source studies\n');

%% ---- 2-3  Monte Carlo and coverage ----------------------------------------------------------------------------------
fprintf('\nStage 2-3: %d records per case (%d if borderline); linearised estimate with the analytic J, Sigma, F\n', ...
    cfg.M1.nRecords, cfg.M1.nRecordsExtended);
state = cell(1, numel(cases)); sm = cell(1, numel(cases));
for q = sel
    file = fullfile(ckDir, sprintf('M1_ckpt_%s.mat', cases(q).name));
    st = [];
    if ~FORCE, [st, ok, msg] = T.loadCheckpoint(file, T.checkpointKey(ctx{q}, cfg)); else, ok = false; msg = 'FORCE'; end
    if ~ok, st = T.newState(ctx{q}, cfg); end
    st.commits{end + 1} = run.commit;
    if numel(unique(st.commits)) > 1
        fprintf('  *** %s: checkpoint written under commit(s) %s, resumed under %s ***\n', cases(q).name, ...
            strjoin(unique(st.commits(1:end - 1)), ', '), run.commit);
    end
    fprintf('  %d %s: %s\n', q, cases(q).name, msg);
    ck = struct('file', file, 'every', cfg.M1.checkpointEvery, 'verbose', true, 'cfg', cfg);
    st = T.mcBatch(ctx{q}, st, st.nTarget, ck);
    r = T.summarise(st, ctx{q}.tgt, cfg);
    if ~st.extended && st.nTarget == cfg.M1.nRecords && r.borderline
        st.extended = true; st.nTarget = cfg.M1.nRecordsExtended;
        st.decision = sprintf('borderline at %d records (%s): extended to %d', r.n, r.borderlineWhat, st.nTarget);
        fprintf('      %s\n', st.decision);
        T.saveCheckpoint(file, st);
        st = T.mcBatch(ctx{q}, st, st.nTarget, ck);
        r = T.summarise(st, ctx{q}.tgt, cfg);
    elseif isempty(st.decision)
        st.decision = sprintf('not borderline at %d records', r.n);
    end
    T.saveCheckpoint(file, st);
    state{q} = st; sm{q} = r;
    printMC(r, ctx{q}.tgt, cfg);
end

%% ---- 4  nonlinear spot checks --------------------------------------------------------------------------------------
offs = cfg.M1.nonlinearBetaOffsets;
fprintf('\nStage 4: record %d of each case inverted from beta x %s (R true), G3 solver and settings\n', ...
    cfg.M1.nonlinearRecord, strjoin(arrayfun(@(o) sprintf('%.2f', 1 + o), offs, 'UniformOutput', false), ' and '));
nDoneNL = 0; tNL = 0; nTotNL = numel(sel) * numel(offs);
for q = sel
    st = state{q};
    file = fullfile(ckDir, sprintf('M1_ckpt_%s.mat', cases(q).name));
    fun = T.makeModel(ctx{q});
    for o = offs
        if ~isempty(st.nl) && any(abs([st.nl.betaOffset] - o) < 1e-12)
            fprintf('  %s beta x %.2f: from checkpoint\n', cases(q).name, 1 + o);
            nTotNL = nTotNL - 1;
            continue
        end
        nl = T.nonlinear(ctx{q}, st, o, fun, cfg);
        if isempty(st.nl), st.nl = nl; else, st.nl(end + 1) = nl; end
        T.saveCheckpoint(file, st);
        nDoneNL = nDoneNL + 1; tNL = tNL + nl.seconds;
        fprintf(['  %s beta x %.2f: %s in %d it (%d J recomputes, %d model evals, %.0f s); chi2/dof %.3f;\n' ...
            '      NL - lin [%+.3f %+.3f] sigma, NL - truth [%+.2f %+.2f] sigma  %s\n'], cases(q).name, 1 + o, ...
            ternary(nl.converged, 'converged', 'NOT converged'), nl.it, nl.jRecomputes, nl.evals, nl.seconds, ...
            nl.chi2dof, nl.gapSigma, nl.errSigma, okStr(nl.pass));
        if nDoneNL == 1 && nTotNL > 1
            fprintf('      about %.0f min for the remaining %d inversion(s)\n', (nTotNL - 1) * tNL / 60, nTotNL - 1);
        end
        flushOut();
    end
    state{q} = st;
    sm{q} = T.summarise(st, ctx{q}.tgt, cfg);
end

%% ---- 5  verdict, tables ---------------------------------------------------------------------------------------------
fprintf('\nVerdict (emp/pred in [%.2f, %.2f], 95%% coverage in [%.2f, %.2f], nonlinear within %.1f sigma of linearised)\n', ...
    cfg.M1.accept.empPred, cfg.M1.accept.coverage95, cfg.M1.accept.nonlinVsLinSigma);
for q = sel
    r = sm{q};
    fprintf('  %-18s n %3d  emp/pred [%.3f %.3f] %s  cov95 %.3f %s  nonlinear %s  =>  %s\n', cases(q).name, r.n, ...
        r.empPred, okStr(r.passScatter), r.coverage(2), okStr(r.passCoverage), okStr(r.passNonlinear), ...
        ternary(r.pass, 'PASS', 'FAIL'));
end

rowsA = {}; rowsS = {}; rowsR = {}; rowsN = {};
for q = sel
    t = ctx{q}.tgt; s = cases(q).source; r = sm{q}; st = state{q};
    rowsA(end + 1, :) = {cases(q).name, cases(q).study, 100 * t.sigma(1), 100 * t.sigma(2), t.dW, t.condF, t.class, ...
        t.nUsedBins, 100 * s.sigma(1), 100 * s.sigma(2), s.dW, cmp{q}.relSigma, cmp{q}.dDW, cmp{q}.ok}; %#ok<SAGROW>
    rowsS(end + 1, :) = {cases(q).name, r.n, 100 * r.sigmaPred(1), 100 * r.sigmaEmp(1), r.empPred(1), ...
        100 * r.sigmaPred(2), 100 * r.sigmaEmp(2), r.empPred(2), 100 * r.dthetaW(1), 100 * r.mean(1), 100 * r.se(1), ...
        100 * r.dthetaW(2), 100 * r.mean(2), 100 * r.se(2), r.offsetMatchesWelch, r.meanQ, r.expectedQ, ...
        r.coverage(1), r.coverage(2), r.nRecordsWithRejectedBins, r.borderline, st.decision, r.passScatter, ...
        r.passCoverage, r.passNonlinear, r.pass}; %#ok<SAGROW>
    for k = 1:st.nDone
        d = st.dtheta(:, k);
        rowsR(end + 1, :) = {cases(q).name, k, st.seeds(k), 100 * (exp(d(1)) - 1), 100 * (exp(d(2)) - 1), d(1), d(2), ...
            st.q(k), t.nUsedBins - st.nBad(k), st.nBad(k)}; %#ok<SAGROW>
    end
    for j = 1:numel(st.nl)
        n = st.nl(j);
        rowsN(end + 1, :) = {cases(q).name, 1 + n.betaOffset, n.record, n.seed, n.converged, n.it, n.jRecomputes, ...
            n.evals, n.chi2dof, n.gapSigma(1), n.gapSigma(2), n.errSigma(1), n.errSigma(2), n.seconds, n.pass}; %#ok<SAGROW>
    end
end
writeStudyCsv(run, 'analytic', {'case', 'source_study', 'sigma_lnbeta_pct', 'sigma_lnR_pct', 'dW', 'kappa', 'class', ...
    'bins_used', 'source_sigma_lnbeta_pct', 'source_sigma_lnR_pct', 'source_dW', 'rel_diff_sigma', 'abs_diff_dW', ...
    'reproduced'}, rowsA);
writeStudyCsv(run, 'summary', {'case', 'n_records', 'pred_sigma_lnbeta_pct', 'emp_sigma_lnbeta_pct', 'emp_pred_beta', ...
    'pred_sigma_lnR_pct', 'emp_sigma_lnR_pct', 'emp_pred_R', 'welch_bias_lnbeta_pct', 'mean_offset_lnbeta_pct', ...
    'se_lnbeta_pct', 'welch_bias_lnR_pct', 'mean_offset_lnR_pct', 'se_lnR_pct', 'offset_matches_welch', 'mean_q', ...
    'expected_q', 'coverage68', 'coverage95', 'records_with_rejected_bins', 'borderline', 'decision', 'pass_scatter', ...
    'pass_coverage', 'pass_nonlinear', 'pass'}, rowsS);
writeStudyCsv(run, 'records', {'case', 'record', 'seed', 'beta_err_pct', 'R_err_pct', 'dlnbeta', 'dlnR', 'q', ...
    'bins_used', 'bins_rejected'}, rowsR);
writeStudyCsv(run, 'nonlinear', {'case', 'beta_start_factor', 'record', 'seed', 'converged', 'iterations', ...
    'J_recomputes', 'model_evals', 'chi2_dof', 'gap_lnbeta_sigma', 'gap_lnR_sigma', 'err_lnbeta_sigma', ...
    'err_lnR_sigma', 'seconds', 'pass'}, rowsN);

%% ---- 6  summary and save --------------------------------------------------------------------------------------------
L = {};
L{end + 1} = sprintf('M1: Monte Carlo and nonlinear spot checks of the analytic (E1) predictions');
L{end + 1} = sprintf('run %s, %s, code commit %s', run.stamp, run.platform, run.commit);
L{end + 1} = '';
L{end + 1} = sprintf(['Acceptance (pre-registered): emp/pred in [%.2f, %.2f]; 95%% coverage of q in [%.2f, %.2f]; nonlinear ' ...
    'converged and within %.1f sigma of the linearised estimate (both starts).'], cfg.M1.accept.empPred, ...
    cfg.M1.accept.coverage95, cfg.M1.accept.nonlinVsLinSigma);
for q = sel
    r = sm{q}; t = ctx{q}.tgt;
    L{end + 1} = ''; %#ok<SAGROW>
    L{end + 1} = sprintf('%s (%s): %s', cases(q).name, cases(q).rule, ternary(r.pass, 'PASS', 'FAIL')); %#ok<SAGROW>
    L{end + 1} = sprintf('  predicted sigma [%.4f%% %.4f%%], d_W %.3f, %s; %d records (%s)', 100 * t.sigma, t.dW, ...
        t.class, r.n, state{q}.decision); %#ok<SAGROW>
    L{end + 1} = sprintf('  emp/pred [%.3f +- %.3f, %.3f +- %.3f]; coverage 68%% %.3f, 95%% %.3f; mean q %.2f (2 + d_W^2 = %.2f)', ...
        [r.empPred(:), r.seEmpPred(:)].', r.coverage, r.meanQ, r.expectedQ); %#ok<SAGROW>
    L{end + 1} = sprintf('  mean offset [%+.4f%% %+.4f%%] +- [%.4f%% %.4f%%] vs Welch-expected [%+.4f%% %+.4f%%]: %s', ...
        100 * r.mean, 100 * r.se, 100 * r.dthetaW, ternary(r.offsetMatchesWelch, 'consistent (3 SE)', ...
        'NOT within 3 SE (diagnostic)')); %#ok<SAGROW>
    for j = 1:numel(state{q}.nl)
        n = state{q}.nl(j);
        L{end + 1} = sprintf('  nonlinear from beta x %.2f: %s, %d it, NL - lin [%+.3f %+.3f] sigma, chi2/dof %.3f', ...
            1 + n.betaOffset, ternary(n.converged, 'converged', 'NOT converged'), n.it, n.gapSigma, n.chi2dof); %#ok<SAGROW>
    end
end
results = struct('cases', rmfield(cases, 'scn'), 'targets', {cellfun(@(c) c.tgt, ctx(sel), 'UniformOutput', false)}, ...
    'compare', {cmp(sel)}, 'state', {state(sel)}, 'summary', {sm(sel)}, 'selected', sel, ...
    'names', {{cases(sel).name}}, 'tools', T.version);
finishStudy(run, results, cfg, L);
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================
function printMC(r, t, cfg)
fprintf(['      %d records: emp/pred [%.3f +- %.3f, %.3f +- %.3f]  %s; coverage 68%% %.3f, 95%% %.3f  %s; mean q %.2f ' ...
    '(2 + d_W^2 = %.2f)\n      mean offset [%+.4f%% %+.4f%%] +- [%.4f%% %.4f%%], Welch-expected [%+.4f%% %+.4f%%]: %s; ' ...
    'records with rejected bins %d\n'], r.n, [r.empPred(:), r.seEmpPred(:)].', okStr(r.passScatter), r.coverage, ...
    okStr(r.passCoverage), r.meanQ, r.expectedQ, 100 * r.mean, 100 * r.se, 100 * t.dthetaW, ...
    ternary(r.offsetMatchesWelch, 'consistent', 'NOT within 3 SE'), r.nRecordsWithRejectedBins);
if r.borderline, fprintf('      borderline: %s\n', r.borderlineWhat); end
if cfg.M1.nRecords < 2, fprintf('      (too few records for statistics)\n'); end
end

function s = okStr(ok)
if isnan(ok), s = 'n/a'; elseif ok, s = 'PASS'; else, s = 'FAIL'; end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function flushOut()
if exist('OCTAVE_VERSION', 'builtin'), fflush(stdout); else, drawnow; end
end