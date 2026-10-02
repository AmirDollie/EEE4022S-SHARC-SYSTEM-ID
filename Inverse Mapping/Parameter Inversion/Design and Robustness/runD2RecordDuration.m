%% runD2RecordDuration.m
% STUDY D2 (SHOULD): how long must a record be, and which Welch segment length L should be used?
%
% QUESTION. For each Welch segment length, how does record duration trade random uncertainty against
% finite-resolution bias, and what is the shortest record that satisfies the pre-registered
% identifiability criteria? At fixed L the random error falls as sigma ~ 1/sqrt(N_eff - 1) while the
% Welch bias dtheta_W does not change, so its significance d_W grows as sqrt(N_eff - 1): a longer record
% makes the same bias more significant. Short L gives many segments but strong window smoothing
% (large bias); long L resolves the response but needs a long record for enough segments.
% "Precise" (sigma small) is therefore not the same as "identifiable" (sigma, d_W and kappa all pass).
%
% METHOD (pre-registered, cfg.D2). Durations cfg.D2.durationMin x L in cfg.D2.L x noise cases
% cfg.D2.noiseCases, at p0 with the Level 2C layout, through the E1 engine (evaluateScenario). Only the
% record and the Welch settings change: the floe, layout, reference, sea state and EMM field are fixed
% (no new EMM solves). A combination is admitted only with K >= cfg.D2.minSegments Welch segments;
% the others are reported as inadmissible and not evaluated. The 109 min point uses the production
% record (N = 65536), so it reproduces E1 exactly. Also reported: RMSE = sqrt(sigma^2 + dtheta_W^2).
% Between grid points the engine's exact scaling at fixed L (sigma sqrt(N_eff - 1) and
% d_W / sqrt(N_eff - 1) constant, kappa unchanged) is used on a 0.1 min grid; it is checked on the grid.
%
% DELIVERABLES (per noise case; headline cfg.D2.headlineNoise, the LSM6DSV16X)
%   1. shortest record with sigma_lnbeta <= 1% (precision only), per L;
%   2. shortest record classified Identifiable (sigma <= 1%, d_W <= 0.5, kappa <= 1e8), over all L,
%      and the identifiable window [T_lo, T_hi] of every L (d_W grows with T, so the window can close);
%   3. the preferred L at that duration (earliest T_lo; ties go to the widest window);
%   4. a verdict on the production choice (109 min, L = 2048): class there, and how it compares with (2).
%   (d_tot = d_W here: no systematic error is applied in D2.)
%
% COST. No new EMM solves (the p0 field is cached); about 140 engine evaluations, a few minutes.
%
% OUTPUT. Results/D2/ (log, .mat, summary, CSV tables, figures D2_fig1, D2_fig2 via plotD2RecordDuration).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();
run = startStudy('D2', 'Record duration and Welch segment length', P.results);

try
Tmin = cfg.D2.durationMin(:).';
Ls = cfg.D2.L(:).';
noiseNames = cfg.D2.noiseCases;
nT = numel(Tmin); nL = numel(Ls); nN = numel(noiseNames);
dt = cfg.acq.dt;
Nof = @(T) round(T * 60 / dt);
Nrec = arrayfun(Nof, Tmin);
Nrec(Tmin == 109) = cfg.acq.Nref;                  % production record, 109.2 min
Treal = Nrec * dt / 60;
headline = cfg.D2.headlineNoise;
I = cfg.class.identifiable;

%% ---- 1  admissibility (noise independent) --------------------------------------------------------------
Kseg = zeros(nL, nT); admit = false(nL, nT);
for a = 1:nL
    for t = 1:nT
        [~, Kseg(a, t)] = welchNEff(Nrec(t), Ls(a), cfg.welch.overlap);
        admit(a, t) = Kseg(a, t) >= cfg.D2.minSegments;
    end
end
fprintf('Step 1: admissible (K >= %d) combinations: %d of %d\n', cfg.D2.minSegments, nnz(admit), numel(admit));
fprintf('  %-8s %s\n', 'T (min)', sprintf('%7g', Tmin));
for a = 1:nL, fprintf('  L %-6d K %s\n', Ls(a), sprintf('%7d', Kseg(a, :))); end

%% ---- 2  sweep ------------------------------------------------------------------------------------------
nanA = NaN(nN, nL, nT);
R = struct('sigB', nanA, 'sigR', nanA, 'dthW', NaN(nN, nL, nT, 2), 'dW', nanA, 'rmseB', nanA, 'rmseR', nanA, ...
    'nEff', nanA, 'nBins', nanA, 'nUsed', nanA, 'invalid', nanA, 'condF', nanA, 'class', {repmat({'inadmissible (K < 8)'}, nN, nL, nT)});
fprintf('\nStep 2: sweep, %d noise cases x %d admissible (L, T)\n', nN, nnz(admit));
for n = 1:nN
    for a = 1:nL
        for t = 1:nT
            if ~admit(a, t), continue; end
            scn = struct('Snn', cfg.noise.(noiseNames{n}), 'L', Ls(a), 'N', Nrec(t), 'cacheDir', P.cache);
            m = evaluateScenario(scn, cfg);
            R.sigB(n, a, t) = m.sigma(1); R.sigR(n, a, t) = m.sigma(2);
            R.dthW(n, a, t, :) = m.dthetaW(:);
            R.dW(n, a, t) = m.dW;
            R.rmseB(n, a, t) = sqrt(m.sigma(1)^2 + m.dthetaW(1)^2);
            R.rmseR(n, a, t) = sqrt(m.sigma(2)^2 + m.dthetaW(2)^2);
            R.nEff(n, a, t) = m.nEff; R.nBins(n, a, t) = numel(m.omega); R.nUsed(n, a, t) = m.nUsedBins;
            R.invalid(n, a, t) = m.invalidFraction; R.condF(n, a, t) = m.condF;
            R.class{n, a, t} = m.class;
        end
    end
    a0 = Ls == 2048; t0 = Tmin == 109;
    fprintf('  %-11s done (109 min, L = 2048: sigma_lnbeta %.4f, d_W %.3f, %s)\n', noiseNames{n}, ...
        R.sigB(n, a0, t0), R.dW(n, a0, t0), R.class{n, a0, t0});
end

%% ---- 3  scaling check at fixed L ----------------------------------------------------------------------
cS = max(R.sigB, R.sigR) .* sqrt(R.nEff - 1);       % s_max sqrt(N_eff - 1)
cD = R.dW ./ sqrt(R.nEff - 1);
devS = 0; devD = 0;
for n = 1:nN
    for a = 1:nL
        v = squeeze(cS(n, a, :)); v = v(isfinite(v)); if numel(v) > 1, devS = max(devS, max(abs(v / median(v) - 1))); end
        v = squeeze(cD(n, a, :)); v = v(isfinite(v)); if numel(v) > 1, devD = max(devD, max(abs(v / median(v) - 1))); end
    end
end
fprintf('\nStep 3: scaling check: s_max sqrt(N_eff - 1) varies by at most %.2f%%, d_W / sqrt(N_eff - 1) by %.2f%%\n', ...
    100 * devS, 100 * devD);

%% ---- 4  deliverables on a fine duration grid ----------------------------------------------------------
Tfine = 1:0.1:600;                                  % minutes
nEffF = zeros(nL, numel(Tfine)); Kf = zeros(nL, numel(Tfine));
for a = 1:nL
    for i = 1:numel(Tfine)
        Ni = Nof(Tfine(i));
        if Ni < Ls(a)                               % shorter than one segment: no estimate at all
            nEffF(a, i) = NaN; Kf(a, i) = 0;
        else
            [nEffF(a, i), Kf(a, i)] = welchNEff(Ni, Ls(a), cfg.welch.overlap);
        end
    end
end
Tsig = NaN(nN, nL); TsigKbound = false(nN, nL);    % (1) sigma_lnbeta <= 1% (precision only)
Tlo = NaN(nN, nL); Thi = NaN(nN, nL);              % (2) identifiable window per L
for n = 1:nN
    for a = 1:nL
        okK = Kf(a, :) >= cfg.D2.minSegments;
        if ~any(isfinite(cS(n, a, :))), continue; end
        cb = median(R.sigB(n, a, isfinite(R.sigB(n, a, :))) .* sqrt(R.nEff(n, a, isfinite(R.sigB(n, a, :))) - 1));
        cs = median(cS(n, a, isfinite(cS(n, a, :))));
        cd = median(cD(n, a, isfinite(cD(n, a, :))));
        kv = squeeze(R.condF(n, a, :)); kv = kv(isfinite(kv));
        if isempty(kv), continue; end
        kap = max(kv);                              % kappa does not change with N at fixed L
        % scaled metrics only where the combination is admissible (K >= 8, so N_eff - 1 > 0)
        sB = NaN(size(Tfine)); sM = sB; dW = sB;
        nd = nEffF(a, okK) - 1;
        sB(okK) = cb ./ sqrt(nd);
        sM(okK) = cs ./ sqrt(nd);
        dW(okK) = cd .* sqrt(nd);
        i1 = find(okK & sB <= 0.01, 1, 'first');
        if ~isempty(i1)
            Tsig(n, a) = Tfine(i1);
            firstAdmissible = find(okK, 1, 'first');   % the shortest admissible record already meets 1%:
            TsigKbound(n, a) = ~isempty(firstAdmissible) && i1 == firstAdmissible;   % set by K >= 8, not sigma
        end
        ident = okK & sM <= I.maxSigma & dW <= I.maxBiasDistance & kap <= I.maxCondF;
        i2 = find(ident, 1, 'first');
        if ~isempty(i2)
            Tlo(n, a) = Tfine(i2);
            i3 = find(~ident(i2:end), 1, 'first');
            if isempty(i3), Thi(n, a) = Inf; else, Thi(n, a) = Tfine(i2 + i3 - 2); end
        end
    end
end
% shortest identifiable record and the preferred L there
Tbest = NaN(nN, 1); Lbest = NaN(nN, 1);
for n = 1:nN
    if all(isnan(Tlo(n, :))), continue; end
    Tbest(n) = min(Tlo(n, :));
    cand = find(Tlo(n, :) <= Tbest(n) + 0.05);
    [~, k] = max(Thi(n, cand));                    % ties: widest identifiable window
    Lbest(n) = Ls(cand(k));
end
% RMSE-optimal L at each grid duration (ln beta)
Lopt = NaN(nN, nT);
for n = 1:nN
    for t = 1:nT
        v = R.rmseB(n, :, t);
        if any(isfinite(v)), [~, a] = min(v); Lopt(n, t) = Ls(a); end
    end
end
% production verdict
a0 = find(Ls == 2048); t0 = find(Tmin == 109);
verdict = cell(nN, 1);
for n = 1:nN
    c = R.class{n, a0, t0};
    if strcmp(c, cfg.class.names{1})
        verdict{n} = sprintf('Identifiable; 109 min is %.1fx the shortest identifiable record (%s min, L = %d)', ...
            109.2 / Tbest(n), num2str(Tbest(n)), Lbest(n));
    else
        why = {};
        if max(R.sigB(n, a0, t0), R.sigR(n, a0, t0)) > I.maxSigma, why{end + 1} = 'sigma > 1%'; end %#ok<AGROW>
        if R.dW(n, a0, t0) > I.maxBiasDistance, why{end + 1} = sprintf('d_W %.2f > 0.5', R.dW(n, a0, t0)); end %#ok<AGROW>
        alt = '';
        okL = find(Tlo(n, :) <= 109.2 & Thi(n, :) >= 109.2);
        if ~isempty(okL), alt = sprintf('; at 109 min use L = %s', strjoin(arrayfun(@num2str, Ls(okL), 'UniformOutput', false), ' or ')); end
        verdict{n} = sprintf('%s (%s)%s', c, strjoin(why, ', '), alt);
    end
end
% the old local rule d ~ 0.26 sqrt(N/65536) (2048/L)^2 (G3 noise), against the engine
gI = find(strcmp(noiseNames, 'g3'));
ruleRatio = NaN(nL, nT);
for a = 1:nL
    for t = 1:nT
        if admit(a, t), ruleRatio(a, t) = 0.26 * sqrt(Nrec(t) / 65536) * (2048 / Ls(a))^2 / R.dW(gI, a, t); end
    end
end

%% ---- 5  tables ------------------------------------------------------------------------------------------
rowsT = {};
for n = 1:nN
    for a = 1:nL
        for t = 1:nT
            rowsT(end + 1, :) = {noiseNames{n}, Ls(a), Tmin(t), Nrec(t), Kseg(a, t), R.nEff(n, a, t), ...
                100 * R.sigB(n, a, t), 100 * R.sigR(n, a, t), R.dW(n, a, t), 100 * R.rmseB(n, a, t), ...
                100 * R.rmseR(n, a, t), R.condF(n, a, t), R.nBins(n, a, t), R.nUsed(n, a, t), R.invalid(n, a, t), ...
                R.class{n, a, t}}; %#ok<AGROW>
        end
    end
end
writeStudyCsv(run, 'sweep', {'noise_case', 'L', 'T_min', 'N', 'K', 'N_eff', 'sigma_lnbeta_pct', 'sigma_lnR_pct', ...
    'd_W', 'rmse_lnbeta_pct', 'rmse_lnR_pct', 'kappa_F', 'bins_candidate', 'bins_used_in_F', 'invalid_fraction', ...
    'class'}, rowsT);
rowsT = cell(nN, 4 + 3 * nL);
for n = 1:nN
    rowsT(n, :) = [{noiseNames{n}, Tbest(n), Lbest(n)}, num2cell(Tsig(n, :)), num2cell(Tlo(n, :)), ...
        num2cell(Thi(n, :)), {verdict{n}}];
end
writeStudyCsv(run, 'min_duration', [{'noise_case', 'T_identifiable_min', 'L_preferred'}, ...
    arrayfun(@(L) sprintf('T_sigma1pct_L%d', L), Ls, 'UniformOutput', false), ...
    arrayfun(@(L) sprintf('T_ident_lo_L%d', L), Ls, 'UniformOutput', false), ...
    arrayfun(@(L) sprintf('T_ident_hi_L%d', L), Ls, 'UniformOutput', false), {'production_verdict'}], rowsT);
rowsT = cell(nN, 1 + nT);
for n = 1:nN, rowsT(n, :) = [{noiseNames{n}}, num2cell(Lopt(n, :))]; end
writeStudyCsv(run, 'rmse_optimal_L', [{'noise_case'}, arrayfun(@(T) sprintf('T%gmin', T), Tmin, ...
    'UniformOutput', false)], rowsT);

%% ---- 6  summary ----------------------------------------------------------------------------------------
Lhdr = arrayfun(@(x) sprintf('L=%d', x), Ls, 'UniformOutput', false);
Thdr = arrayfun(@(x) sprintf('%gmin', x), Tmin, 'UniformOutput', false);
fmtT = @(v) sprintf('%10s', tStr(v));
L = {};
L{end + 1} = sprintf('Setting: p0, Level 2C (ref s26), JONSWAP Hs %.2f m, omega_p %.1f rad/s; admitted if K >= %d.', ...
    cfg.sea.Hs, cfg.sea.omegaP, cfg.D2.minSegments);
L{end + 1} = sprintf('Engine scaling at fixed L: s_max sqrt(N_eff - 1) within %.2f%%, d_W / sqrt(N_eff - 1) within %.2f%%.', ...
    100 * devS, 100 * devD);
L{end + 1} = '';
L{end + 1} = 'HEADLINE: shortest record classified Identifiable, preferred L, and the production verdict';
for n = 1:nN
    L{end + 1} = sprintf('  %-12s %s min with L = %s | 109 min, L = 2048: %s', noiseNames{n}, tStr(Tbest(n)), ...
        num2str(Lbest(n)), verdict{n}); %#ok<AGROW>
end
L{end + 1} = '';
L{end + 1} = '(1) Shortest record (min) with sigma_lnbeta <= 1% (precision only; * = set by the K >= 8 rule)';
L{end + 1} = sprintf('  %-12s %s', 'noise case', sprintf('%10s', Lhdr{:}));
for n = 1:nN
    s = '';
    for a = 1:nL
        v = tStr(Tsig(n, a)); if TsigKbound(n, a), v = [v '*']; end
        s = [s sprintf('%10s', v)]; %#ok<AGROW>
    end
    L{end + 1} = sprintf('  %-12s %s', noiseNames{n}, s); %#ok<AGROW>
end
L{end + 1} = '';
L{end + 1} = '(2) Identifiable window [T_lo, T_hi] (min) per L: sigma <= 1%, d_W <= 0.5, kappa <= 1e8';
L{end + 1} = sprintf('  %-12s %s', 'noise case', sprintf('%16s', Lhdr{:}));
for n = 1:nN
    s = '';
    for a = 1:nL
        if isnan(Tlo(n, a)), w = 'never'; else, w = sprintf('[%s, %s]', tStr(Tlo(n, a)), tStr(Thi(n, a))); end
        s = [s sprintf('%16s', w)]; %#ok<AGROW>
    end
    L{end + 1} = sprintf('  %-12s %s', noiseNames{n}, s); %#ok<AGROW>
end
L{end + 1} = '  (windows close because d_W grows with record length at fixed L; searched to 600 min)';
L{end + 1} = '';
L{end + 1} = '(3) RMSE-optimal L (ln beta) at each grid duration';
L{end + 1} = sprintf('  %-12s %s', 'noise case', sprintf('%8s', Thdr{:}));
for n = 1:nN, L{end + 1} = sprintf('  %-12s %s', noiseNames{n}, sprintf('%8.0f', Lopt(n, :))); end %#ok<AGROW>
L{end + 1} = '';
rr = ruleRatio(isfinite(ruleRatio));
L{end + 1} = sprintf('Old local rule d ~ 0.26 sqrt(N/65536)(2048/L)^2 vs engine (G3 noise): ratio %.2f to %.2f.', ...
    min(rr), max(rr));
results = struct('Tmin', Tmin, 'Treal', Treal, 'L', Ls, 'N', Nrec, 'K', Kseg, 'admit', admit, ...
    'noise', {noiseNames}, 'headline', headline, 'R', R, 'Tfine', Tfine, 'nEffFine', nEffF, 'Kfine', Kf, ...
    'Tsigma1pct', Tsig, 'TsigmaKbound', TsigKbound, 'TidentLo', Tlo, 'TidentHi', Thi, 'Tbest', Tbest, ...
    'Lbest', Lbest, 'verdict', {verdict}, 'Lopt', Lopt, 'ruleRatio', ruleRatio, 'scalingDev', [devS devD]);
plotD2RecordDuration(results, cfg, run);            % figures (also redrawable later from the .mat)
finishStudy(run, results, cfg, L);
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================
function s = tStr(v)
if isnan(v), s = '-'; elseif isinf(v), s = '>600'; else, s = sprintf('%.1f', v); end
end