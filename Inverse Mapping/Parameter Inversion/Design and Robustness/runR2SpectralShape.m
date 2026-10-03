%% runR2SpectralShape.m
% STUDY R2 (MUST): does the SHAPE of the sea spectrum change what the inverse can do?
%
% QUESTION. With the floe (p0), layout (Level 2C, reference s26), noise (LSM6DSV16X), record (N = 65536)
% and the in-band signal level all fixed, how do the peak frequency, the peakedness and a bimodal sea
% change sigma, d_W and the class, and is L = 2048 still adequate for narrow and bimodal seas?
%
% WHY THE SHAPE CAN MATTER AT ALL. For one incident component S_jr = H_j H_r* S_uu and S_rr = |H_r|^2 S_uu,
% so the exact transmissibility T_j = H_j / H_r does not depend on the sea: f0 and J are identical for
% every spectrum (checked here). But (i) the covariance Sigma_k depends on the per-bin SNR, so the sea
% decides WHICH bins carry the information, and (ii) the Welch estimate smooths S_jr and S_rr with the
% spectral window BEFORE the ratio, so where T varies across the window the result depends on how S_uu
% weights neighbouring frequencies: the Welch bias b_W depends on the sea.
%
% METHOD (pre-registered, cfg.R2; clarifications logged 2026-10-03 before the run)
%   1. 16 seas from buildR2Seas: JONSWAP omega_p in cfg.R2.omegaP x gamma_J in cfg.R2.gammaJ, plus the
%      bimodal sea; every sea normalised to the reference sea's elevation variance in the retained band.
%   2. Each sea at each L in cfg.R2.checkL through evaluateScenario (the frozen E1 chain, cached p0 field).
%   3. Mechanism (L = 2048): per-bin Fisher information and per-bin standardised Welch bias, the sea's PSD
%      on the bins, and the invariance of f0 and J across seas.
%   4. L verdict per sea: L = 2048 adequate if Identifiable at 2048 (cfg.R2.adequacyRule); otherwise
%      whether 1024 or 4096 is. Also the L with the smallest RMSE in ln beta, sqrt(sigma^2 + dtheta_W^2).
%   5. M1 case (cfg.R2.m1Rule, fixed before the run): among the gamma_J = 7 seas and the bimodal sea, the
%      largest q = max(s_max / 0.01, d_W / 0.5, kappa / 1e8) at L = 2048.
%
% RUN ORDER. testR2SpectralShape must pass first.
% COST. 48 analytic evaluations on the cached p0 field: minutes, no EMM solve.
% OUTPUT. Results/R2/ (log, .mat, summary, R2_sweep.csv, R2_L_comparison.csv, R2_M1_selection.csv;
% figures via plotR2SpectralShape once it exists).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();
[seas, mkScn] = buildR2Seas(cfg);
run = startStudy('R2', 'Spectral shape at fixed in-band variance', P.results);

try
nQ = numel(seas); Ls = cfg.R2.checkL(:).'; nL = numel(Ls);
iL = find(Ls == cfg.welch.L, 1);
if isempty(iL), error('runR2:L', 'cfg.R2.checkL must contain the production L = %d.', cfg.welch.L); end
iRef = find([seas.isReference]);
fprintf(['Step 1: %d seas (%d JONSWAP + bimodal), in-band variance %.4g m^2 in [%.3f, %.3f] rad/s ' ...
    '(spread %.1e relative); reference %s; L = %s; noise %s\n'], nQ, nQ - 1, seas(iRef).Vinband, ...
    cfg.R2.normalisationBand, max(abs([seas.Vinband] / seas(iRef).Vinband - 1)), seas(iRef).name, mat2str(Ls), ...
    cfg.R2.noiseCase);
fprintf('  %-22s %8s %8s %9s\n', 'sea', 'omega_p', 'gamma_J', 'Hs (m)');
for q = 1:nQ
    if strcmp(seas(q).type, 'jonswap')
        fprintf('  %-22s %8.3g %8.3g %9.4f\n', seas(q).name, seas(q).omegaP, seas(q).gammaJ, seas(q).Hs);
    else
        fprintf('  %-22s %8s %8s %9s  (components %s rad/s, gamma_J %s, fractions %s)\n', seas(q).name, '-', '-', ...
            '-', mat2str(seas(q).omegaP), mat2str(seas(q).gammaJ), mat2str(cfg.R2.bimodal.energyFraction));
    end
end

%% ---- 2  the sweep -------------------------------------------------------------------------------------------
fprintf('\nStep 2: %d evaluations\n', nQ * nL);
M = struct('sigma', NaN(nQ, nL, 2), 'corr', NaN(nQ, nL), 'condF', NaN(nQ, nL), 'dW', NaN(nQ, nL), ...
    'dthetaW', NaN(nQ, nL, 2), 'sMax', NaN(nQ, nL), 'code', NaN(nQ, nL), 'nBins', NaN(nQ, nL), ...
    'nUsed', NaN(nQ, nL), 'invalid', NaN(nQ, nL), 'minCoh', NaN(nQ, nL), 'K', NaN(1, nL), 'nEff', NaN(1, nL));
M.class = cell(nQ, nL);
mech = struct('omega', [], 'S', [], 'info', [], 'biasStd', [], 'f0diff', NaN(1, nQ), 'Jdiff', NaN(1, nQ), ...
    'bWdiff', NaN(1, nQ));
store = cell(1, nQ);
tAll = tic;
for q = 1:nQ
    for k = 1:nL
        out = evaluateScenario(mkScn(seas(q), Ls(k), P.cache), cfg);
        M.sigma(q, k, :) = out.sigma(1:2); M.corr(q, k) = out.corr; M.condF(q, k) = out.condF;
        M.dW(q, k) = out.dW; M.dthetaW(q, k, :) = out.dthetaW(1:2); M.sMax(q, k) = out.sMax;
        M.class{q, k} = out.class; M.code(q, k) = find(strcmp(cfg.class.names, out.class));
        M.nBins(q, k) = numel(out.omega); M.nUsed(q, k) = out.nUsedBins; M.invalid(q, k) = out.invalidFraction;
        M.minCoh(q, k) = out.minCoherence; M.K(k) = out.K; M.nEff(k) = out.nEff;
        if k == iL
            if q == 1
                mech.omega = out.omega; nF = numel(out.omega);
                mech.S = NaN(nQ, nF); mech.info = NaN(nQ, nF, 2); mech.biasStd = NaN(nQ, nF);
            end
            mech.S(q, :) = seas(q).S(out.omega);
            vIdx = find(out.validMask);              % perBinInfo covers the accepted bins only
            mech.info(q, vIdx, :) = permute(out.perBinInfo, [3 2 1]);
            mech.biasStd(q, :) = perBinBias(out);
            store{q} = struct('f0', out.f0, 'J', out.Jall, 'bW', out.bW);
        end
    end
    fprintf('  %-22s L=2048: sigma = [%.3f%% %.3f%%], d_W = %.3f, kappa = %.3g, %s\n', seas(q).name, ...
        100 * M.sigma(q, iL, 1), 100 * M.sigma(q, iL, 2), M.dW(q, iL), M.condF(q, iL), M.class{q, iL});
end
for q = 1:nQ                                    % invariance of the exact features and Jacobian across seas
    mech.f0diff(q) = max(abs(store{q}.f0 - store{iRef}.f0));
    mech.Jdiff(q) = max(max(abs(store{q}.J - store{iRef}.J)));
    mech.bWdiff(q) = norm(store{q}.bW - store{iRef}.bW) / norm(store{iRef}.bW);
end
fprintf('  done in %.0f s\n', toc(tAll));
fprintf(['\nStep 3: exact T is sea-independent: max |f0 - f0_ref| = %.1e, max |J - J_ref| = %.1e over all seas; ' ...
    'the Welch bias is not: ||b_W - b_W,ref|| / ||b_W,ref|| from %.2f to %.3g\n'], max(mech.f0diff), max(mech.Jdiff), ...
    min(mech.bWdiff(setdiff(1:nQ, iRef))), max(mech.bWdiff));

%% ---- 4  L verdict -----------------------------------------------------------------------------------------------
rmseB = sqrt(M.sigma(:, :, 1).^2 + M.dthetaW(:, :, 1).^2);
[~, bestL] = min(rmseB, [], 2);
adequate = M.code(:, iL) == 1;
fprintf('\nStep 4: L verdict (adequate = Identifiable at L = %d)\n', cfg.welch.L);
fprintf('  %-22s %s   %s\n', 'sea', strjoin(arrayfun(@(L) sprintf('%-16s', sprintf('L=%d', L)), Ls, 'UniformOutput', false), ''), ...
    'smallest RMSE(ln beta)');
for q = 1:nQ
    cells = arrayfun(@(k) sprintf('%-16s', sprintf('%s d%.2f', abbrev(M.class{q, k}), M.dW(q, k))), 1:nL, ...
        'UniformOutput', false);
    fprintf('  %-22s %s   L=%d\n', seas(q).name, strjoin(cells, ''), Ls(bestL(q)));
end

%% ---- 5  M1 case -----------------------------------------------------------------------------------------------
I = cfg.class.identifiable;
isUni = strcmp({seas.type}, 'jonswap');
gj1 = arrayfun(@(x) x.gammaJ(1), seas);          % per sea (the bimodal sea has two components)
cand = find((isUni & gj1 == 7) | ~isUni);
qScore = max([M.sMax(cand, iL) / I.maxSigma, M.dW(cand, iL) / I.maxBiasDistance, M.condF(cand, iL) / I.maxCondF], [], 2);
[~, w] = max(qScore); m1 = cand(w);
fprintf('\nStep 5: M1 case (%s): %s, q = %.3f\n', cfg.R2.m1Rule, seas(m1).name, qScore(w));

%% ---- 6  tables ------------------------------------------------------------------------------------------------
rows = {};
for q = 1:nQ
    for k = 1:nL
        rows(end + 1, :) = {seas(q).name, seas(q).type, wpStr(seas(q)), gjStr(seas(q)), seas(q).Hs, seas(q).Vinband, ...
            Ls(k), M.K(k), M.nBins(q, k), M.nUsed(q, k), M.sigma(q, k, 1), M.sigma(q, k, 2), M.corr(q, k), ...
            M.condF(q, k), M.dW(q, k), M.dthetaW(q, k, 1), rmseB(q, k), M.minCoh(q, k), M.class{q, k}}; %#ok<SAGROW>
    end
end
writeStudyCsv(run, 'sweep', {'sea', 'type', 'omega_p', 'gamma_J', 'Hs_m', 'inband_variance_m2', 'L', 'segments', ...
    'bins', 'bins_used', 'sigma_lnbeta', 'sigma_lnR', 'rho_betaR', 'condF', 'dW', 'dlnbeta_W', 'rmse_lnbeta', ...
    'min_coherence', 'class'}, rows);
rows = {};
for q = 1:nQ
    r = {seas(q).name};
    for k = 1:nL, r = [r, {M.class{q, k}, M.sigma(q, k, 1), M.dW(q, k)}]; end %#ok<AGROW>
    rows(end + 1, :) = [r, {adequate(q), Ls(bestL(q))}]; %#ok<SAGROW>
end
hdr = {'sea'};
for k = 1:nL, hdr = [hdr, {sprintf('class_L%d', Ls(k)), sprintf('sigma_lnbeta_L%d', Ls(k)), sprintf('dW_L%d', Ls(k))}]; end %#ok<AGROW>
writeStudyCsv(run, 'L_comparison', [hdr, {'L2048_adequate', 'best_L_rmse'}], rows);
rows = arrayfun(@(j) {seas(cand(j)).name, M.sMax(cand(j), iL), M.dW(cand(j), iL), M.condF(cand(j), iL), qScore(j), ...
    cand(j) == m1}, 1:numel(cand), 'UniformOutput', false);
writeStudyCsv(run, 'M1_selection', {'sea', 'sMax', 'dW', 'condF', 'q', 'selected'}, vertcat(rows{:}));

%% ---- 7  summary ----------------------------------------------------------------------------------------------------
L = {};
L{end + 1} = sprintf(['p0, Level 2C (reference s26), %s noise, N = %d; every sea normalised to the reference in-band ' ...
    'variance %.3g m^2 over [%.3f, %.3f] rad/s.'], cfg.R2.noiseCase, cfg.acq.Nref, seas(iRef).Vinband, cfg.R2.normalisationBand);
L{end + 1} = sprintf(['Exact transmissibility is sea-independent (max |df0| = %.1e, max |dJ| = %.1e); the Welch bias ' ...
    'is not: ||b_W - b_W,ref|| / ||b_W,ref|| from %.2f to %.3g.'], max(mech.f0diff), max(mech.Jdiff), ...
    min(mech.bWdiff(setdiff(1:nQ, iRef))), max(mech.bWdiff));
L{end + 1} = '';
L{end + 1} = sprintf('%-22s %10s %10s %8s %8s  %-16s %-16s %-16s', 'sea (L = 2048)', 'sig lnb', 'sig lnR', 'd_W', ...
    'kappa', sprintf('L=%d', Ls(1)), sprintf('L=%d', Ls(2)), sprintf('L=%d', Ls(min(3, nL))));
for q = 1:nQ
    L{end + 1} = sprintf('%-22s %9.3f%% %9.3f%% %8.3f %8.3g  %-16s %-16s %-16s', seas(q).name, 100 * M.sigma(q, iL, 1), ...
        100 * M.sigma(q, iL, 2), M.dW(q, iL), M.condF(q, iL), M.class{q, 1}, M.class{q, min(2, nL)}, M.class{q, min(3, nL)}); %#ok<SAGROW>
end
L{end + 1} = '';
L{end + 1} = sprintf('L = %d adequate (Identifiable) for %d of %d seas.', cfg.welch.L, nnz(adequate), nQ);
bad = find(~adequate);
for q = bad(:).'
    alt = Ls(M.code(q, :) == 1);
    if isempty(alt), s = 'no L in the set is Identifiable'; else, s = sprintf('Identifiable at L = %s', mat2str(alt)); end
    L{end + 1} = sprintf('  %s: %s at L = %d (sigma %.3f%%, d_W %.3f); %s.', seas(q).name, M.class{q, iL}, ...
        cfg.welch.L, 100 * M.sigma(q, iL, 1), M.dW(q, iL), s); %#ok<SAGROW>
end
[~, iBest] = min(M.sigma(:, iL, 1)); [~, iWorst] = max(M.sigma(:, iL, 1));
L{end + 1} = sprintf('sigma_lnbeta at L = 2048 ranges from %.3f%% (%s) to %.3f%% (%s); reference %.3f%%.', ...
    100 * M.sigma(iBest, iL, 1), seas(iBest).name, 100 * M.sigma(iWorst, iL, 1), seas(iWorst).name, 100 * M.sigma(iRef, iL, 1));
L{end + 1} = sprintf('M1 case (%s): %s (q = %.3f).', cfg.R2.m1Rule, seas(m1).name, qScore(w));

seaTable = rmfield(seas, {'sea', 'S'});          % handles are not saved; buildR2Seas rebuilds them
results = struct('seas', seaTable, 'L', Ls, 'iL', iL, 'iRef', iRef, 'metrics', M, 'mechanism', mech, ...
    'rmseBeta', rmseB, 'bestL', Ls(bestL), 'adequate2048', adequate, 'm1', struct('candidates', {{seas(cand).name}}, ...
    'q', qScore, 'selected', seas(m1).name));
if exist('plotR2SpectralShape', 'file') == 2
    try
        plotR2SpectralShape(results, cfg, run);
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
function s = perBinBias(out)
% standardised Welch bias per bin, sqrt(b_k' Sigma_k^-1 b_k) (NaN for rejected bins)
pB = size(out.blocksAll, 1); nF = numel(out.omega);
s = NaN(1, nF);
for k = find(out.validMask(:).')
    rows = (k - 1) * pB + (1:pB);
    [Lk, flag] = chol(out.blocksAll(:, :, k), 'lower');
    if flag == 0, s(k) = norm(Lk \ out.bW(rows)); end
end
end

function a = abbrev(c)
switch c
    case 'Identifiable', a = 'I';
    case 'Marginal', a = 'M';
    otherwise, a = 'N';
end
end

function s = wpStr(x)
s = strjoin(arrayfun(@(v) sprintf('%.4g', v), x.omegaP, 'UniformOutput', false), '+');
end

function s = gjStr(x)
s = strjoin(arrayfun(@(v) sprintf('%.3g', v), x.gammaJ, 'UniformOutput', false), '+');
end