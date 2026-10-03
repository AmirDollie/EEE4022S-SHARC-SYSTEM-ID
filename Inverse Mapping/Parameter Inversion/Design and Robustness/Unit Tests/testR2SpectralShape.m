%% testR2SpectralShape.m
% Pre-flight checks for study R2. Must pass (all) before runR2SpectralShape is run. The seas come from
% buildR2Seas, the same code the runner uses; checks 6-10 use the cached p0 field (no EMM solve).
%
%   1   15 JONSWAP + 1 bimodal = 16 unique seas; every (omega_p, gamma_J) once
%   2   the reference sea exists once, with Hs exactly cfg.sea.Hs (scale 1)
%   3   every sea has the reference retained-band variance
%   4   the bimodal components carry 40% / 60% of that variance
%   5   spectra finite and non-negative; every JONSWAP peaks at its omega_p
%   6   exact T is sea-independent: f0 and J identical for two very different seas
%   7   the reference sea at L = 2048 reproduces the frozen evaluateScenario default exactly
%   8   the Welch bias DOES depend on the sea (the mechanism R2 measures)
%   9   L = 1024, 2048, 4096: R2 path = direct evaluateScenario with that L; segments >= cfg.D2.minSegments
%   10  only the sea and L change: p, sensors, layout, reference, noise, N, dt, band, node grid, truncation
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
P = setupStudyPaths();
cfg = defineStudyScenarios();
nPass = 0; nFail = 0;
[seas, mkScn] = buildR2Seas(cfg);
nQ = numel(seas);
uni = strcmp({seas.type}, 'jonswap');
wp1 = arrayfun(@(x) x.omegaP(1), seas); gj1 = arrayfun(@(x) x.gammaJ(1), seas);   % bimodal: first component

%% ---- 1-5  construction (no EMM) --------------------------------------------------------------------------
combos = [wp1(uni).', gj1(uni).'];
[A, B] = meshgrid(cfg.R2.omegaP, cfg.R2.gammaJ);
want = sortrows([A(:), B(:)]);
ok = nQ == 16 && nnz(uni) == 15 && nnz(strcmp({seas.type}, 'bimodal')) == 1 && numel(unique({seas.name})) == 16 && ...
    isequal(sortrows(combos), want);
[nPass, nFail] = check(ok, sprintf('1  %d seas: %d JONSWAP (every omega_p x gamma_J once) + %d bimodal, unique names', ...
    nQ, nnz(uni), nnz(~uni)), nPass, nFail);

iRef = find([seas.isReference]);
ok = numel(iRef) == 1 && seas(iRef).Hs == cfg.sea.Hs && seas(iRef).scale == 1 && ...
    seas(iRef).omegaP == cfg.sea.omegaP && seas(iRef).gammaJ == cfg.sea.gammaJ && isequal(seas(iRef).sea, cfg.sea);
[nPass, nFail] = check(ok, sprintf('2  reference sea %s once, Hs = %.4g m exactly, passed as the E1 default struct', ...
    seas(iRef).name, seas(iRef).Hs), nPass, nFail);

V = [seas.Vinband];
relV = max(abs(V / V(iRef) - 1));
[nPass, nFail] = check(relV < 1e-9, sprintf('3  retained-band variance %.4g m^2 for every sea (max relative spread %.1e); Hs from %.4f to %.4f m', ...
    V(iRef), relV, min([seas(uni).Hs]), max([seas(uni).Hs])), nPass, nFail);

b = seas(~uni);
fr = b.componentVinband / sum(b.componentVinband);
ok = max(abs(fr - cfg.R2.bimodal.energyFraction)) < 1e-9 && abs(sum(b.componentVinband) / V(iRef) - 1) < 1e-9;
[nPass, nFail] = check(ok, sprintf('4  bimodal components carry %.4f / %.4f of the in-band variance (configured %s)', fr, ...
    mat2str(cfg.R2.bimodal.energyFraction)), nPass, nFail);

wq = linspace(0.5, 15, 30001);
ok = true; worst = 0;
for q = 1:nQ
    Sq = seas(q).S(wq);
    ok = ok && all(isfinite(Sq)) && all(Sq >= 0);
    if uni(q)
        [~, im] = max(Sq);
        worst = max(worst, abs(wq(im) / seas(q).omegaP - 1));
    end
end
ok = ok && worst < 1e-3;
[nPass, nFail] = check(ok, sprintf('5  all spectra finite and non-negative; JONSWAP peaks within %.1e of omega_p', worst), nPass, nFail);

%% ---- 6-8  the E1 chain ------------------------------------------------------------------------------------
fprintf('  ... evaluations on the cached p0 field\n');
qa = find(uni & wp1 == 3.5 & gj1 == 7);
qb = find(uni & wp1 == 7.5 & gj1 == 1);
oa = evaluateScenario(mkScn(seas(qa), cfg.welch.L, P.cache), cfg);
ob = evaluateScenario(mkScn(seas(qb), cfg.welch.L, P.cache), cfg);
ok = isequal(oa.f0, ob.f0) && isequal(oa.Jall, ob.Jall) && isequal(oa.omega, ob.omega);
[nPass, nFail] = check(ok, sprintf('6  exact features and Jacobian identical for %s and %s (max |df0| %.1e, max |dJ| %.1e)', ...
    seas(qa).name, seas(qb).name, max(abs(oa.f0 - ob.f0)), max(max(abs(oa.Jall - ob.Jall)))), nPass, nFail);

Snn = cfg.noise.(cfg.R2.noiseCase);
oR = evaluateScenario(mkScn(seas(iRef), cfg.welch.L, P.cache), cfg);
e1 = evaluateScenario(struct('Snn', Snn, 'cacheDir', P.cache), cfg);
ok = isequal(oR.sigma, e1.sigma) && isequal(oR.F, e1.F) && oR.dW == e1.dW && oR.condF == e1.condF && ...
    strcmp(oR.class, e1.class) && isequal(oR.omega, e1.omega) && oR.nEff == e1.nEff;
[nPass, nFail] = check(ok, sprintf(['7  reference sea at L = %d = evaluateScenario default: sigma [%.4f%% %.4f%%], d_W %.4f, ' ...
    'kappa %.4g, %s, %d bins, N_eff %.3f (bitwise)'], cfg.welch.L, 100 * oR.sigma, oR.dW, oR.condF, oR.class, ...
    numel(oR.omega), oR.nEff), nPass, nFail);

dA = norm(oa.bW - oR.bW) / norm(oR.bW); dB = norm(ob.bW - oR.bW) / norm(oR.bW);
ok = dA > 1e-3 && dB > 1e-3 && ~isequal(oa.F, ob.F);
[nPass, nFail] = check(ok, sprintf(['8  the Welch bias depends on the sea: ||b_W - b_W,ref|| / ||b_W,ref|| = %.3f (%s), %.3f (%s); ' ...
    'F differs too (per-bin SNR weighting)'], dA, seas(qa).name, dB, seas(qb).name), nPass, nFail);

%% ---- 9  L bookkeeping ------------------------------------------------------------------------------------------
ok = true; txt = {};
for L = cfg.R2.checkL
    o1 = evaluateScenario(mkScn(seas(iRef), L, P.cache), cfg);
    o2 = evaluateScenario(struct('Snn', Snn, 'L', L, 'cacheDir', P.cache), cfg);
    Kexp = floor(2 * cfg.acq.Nref / L) - 1;
    ok = ok && isequal(o1.sigma, o2.sigma) && isequal(o1.F, o2.F) && o1.dW == o2.dW && o1.K == Kexp && ...
        Kexp >= cfg.D2.minSegments;
    txt{end + 1} = sprintf('L=%d: K=%d, %d bins, sigma_lnbeta %.3f%%, d_W %.3f', L, o1.K, numel(o1.omega), ...
        100 * o1.sigma(1), o1.dW); %#ok<SAGROW>
end
[nPass, nFail] = check(ok, ['9  R2 path = direct evaluateScenario at each L; ' strjoin(txt, '; ')], nPass, nFail);

%% ---- 10  nothing else changes ------------------------------------------------------------------------------------------
fields = {'p', 'sensors', 'layout', 'ref', 'estIdx', 'nuisanceIdx', 'Snn', 'N', 'dt', 'overlap', 'binStep', 'bandEst', ...
    'bandUse', 'nodeSpacing', 'frfOptions', 'step'};
ok = true; diffs = {};
for q = [qa qb iRef]
    s = oa.scn; if q == qb, s = ob.scn; elseif q == iRef, s = oR.scn; end
    for f = fields
        if ~isequal(s.(f{1}), e1.scn.(f{1})), ok = false; diffs{end + 1} = f{1}; end %#ok<SAGROW>
    end
end
ok = ok && isequal(fieldnames(mkScn(seas(1), 1024, '')).', {'Snn', 'sea', 'L', 'cacheDir'});
msg = 'only the sea and L differ from the E1 default (p0, Level 2C, s26, noise, N, dt, bands, node grid, truncation)';
if ~ok, msg = [msg ': DIFFERS in ' strjoin(unique(diffs), ', ')]; end
[nPass, nFail] = check(ok, ['10 ' msg], nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testR2SpectralShape:failed', '%d check(s) failed: do not run R2.', nFail); end

%% ================================================================================================
function [nPass, nFail] = check(ok, msg, nPass, nFail)
if ok
    fprintf('PASS  %s\n', msg); nPass = nPass + 1;
else
    fprintf('FAIL  %s\n', msg); nFail = nFail + 1;
end
end