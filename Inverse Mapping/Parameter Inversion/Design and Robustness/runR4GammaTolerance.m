%% runR4GammaTolerance.m
% STUDY R4 (SHOULD): how accurately must the KNOWN gamma (draught) be known?
%
% QUESTION. The inverse fixes gamma at an assumed value gamma_inv = gamma_true (1 + eps). The fit then
% absorbs the error through beta and R. How large can eps be before the induced bias matters?
%
% METHOD (pre-registered, cfg.R4).
%   Linear: the feature error at the truth is b = f(gamma_true) - f(gamma_inv) ~ -J_gamma ln(1 + eps)
%   (J_gamma from emmField, never estimated), projected by the E1 engine:
%       dtheta_gamma = F^-1 J' Sigma^-1 b,   d_sys = sqrt(dtheta' F dtheta),
%   and combined with the Welch bias as d_tot. Grid cfg.R4.epsGamma, noise cases cfg.R4.noiseCases.
%   Nonlinear check at cfg.R4.nonlinearCheck: the EXACT feature error (EMM solved at gamma_inv) and a
%   fixed-Jacobian (chord) Gauss-Newton fit (the G3 solver's update with J held at the truth) of the noiseless model,
%   iterated until the step is below 0.01 sigma. Agreement with the linear dtheta shows the linear
%   tolerance is trustworthy.
%   Tolerances by first crossing (programme, Section 3.2): |delta beta / beta| = 1%, d_sys = 0.5, d_sys = 1,
%   found on the fine eps grid cfg.R4.toleranceSearch / toleranceN with the same linear formula.
%
% DELIVERABLE. Required gamma (draught) accuracy for |delta beta/beta| <= 1% and for d_sys <= 0.5, per
% noise case; gamma = rho_i h / (rho_w H), so a relative error in gamma is the same relative error in the
% draught (or in h, with rho_i known). Caveat: beta ~ h^3, so 1% in h is 3% in an independent beta
% reference.
%
% COST. Two extra EMM field evaluations for J_gamma (about 1 minute in MATLAB), plus about 3 per
% nonlinear case (about 6 minutes). Cached in Results/twinCache, so a rerun is fast.
%
% OUTPUT. Results/R4/ (log, .mat, summary, CSV tables, figures/R4_fig1_gamma_tolerance.pdf/.png).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

RUN_NONLINEAR = true;          % false skips the exact-EMM check (testing only)

P = setupStudyPaths();
cfg = defineStudyScenarios();
run = startStudy('R4', 'Uncertainty in the known gamma (draught)', P.results);

try
p0 = cfg.p0.vec;
sensors = cfg.layout.level2C.frac .* [p0(3) 1];      % fixed physical points at the truth
base = struct('p', p0, 'sensors', sensors, 'nuisanceIdx', 2, 'cacheDir', P.cache, 'verbose', true);
noiseNames = cfg.R4.noiseCases;
nN = numel(noiseNames);
epsGrid = cfg.R4.epsGamma(:).';
epsFine = linspace(cfg.R4.toleranceSearch(1), cfg.R4.toleranceSearch(2), cfg.R4.toleranceN);
headline = cfg.R4.headlineNoise;

%% ---- 1  J_gamma at p0 --------------------------------------------------------
fprintf('Step 1: field, Jacobian and gamma sensitivity at p0 (7 field evaluations)\n');
m0 = evaluateScenario(setfield(base, 'Snn', cfg.noise.(headline)), cfg); %#ok<SFLD>
Jg = m0.Jnuisance;                                    % d f / d ln gamma, all bins
pBlk = size(m0.blocksAll, 1);
validRows = @(m) reshape((find(m.validMask) - 1) * pBlk + (1:pBlk).', [], 1);

%% ---- 2  linear sweep, every noise case ---------------------------------------------------------------
fprintf('\nStep 2: linear projection for %d eps values x %d noise cases\n', numel(epsGrid), nN);
res = struct('noise', noiseNames, 'Snn', [], 'sigma', [], 'slope', [], 'dthetaSys', [], 'dSys', [], ...
    'dTot', [], 'class', [], 'dBeta', [], 'dR', [], 'tolBeta1', [], 'tolD05', [], 'tolD1', [], 'dW', []);
for n = 1:nN
    Snn = cfg.noise.(noiseNames{n});
    m = evaluateScenario(setfield(setfield(base, 'Snn', Snn), 'verbose', false), cfg); %#ok<SFLD>
    rows = validRows(m);
    % unit sensitivity: dtheta per unit ln(1 + eps), and the grid (sys alone and with the Welch bias)
    Bg = -Jg * log(1 + epsGrid);
    q = fisherFromBlocks(m.Jall(rows, :), m.blocksAll(:, :, m.validMask), [-Jg(rows), Bg(rows, :), ...
        Bg(rows, :) + m.bW(rows)]);
    nE = numel(epsGrid);
    res(n).Snn = Snn;
    res(n).sigma = m.sigma(:).';
    res(n).dW = m.dW;
    res(n).slope = q.dtheta(:, 1).';                 % d theta / d ln(1 + eps)
    res(n).dthetaSys = q.dtheta(:, 1 + (1:nE));
    res(n).dSys = q.d(1 + (1:nE));
    res(n).dTot = q.d(1 + nE + (1:nE));
    res(n).dBeta = exp(res(n).dthetaSys(1, :)) - 1;
    res(n).dR = exp(res(n).dthetaSys(2, :)) - 1;
    res(n).class = arrayfun(@(k) classOf(max(m.sigma), res(n).dTot(k), m.condF, m.singular, cfg), 1:nE, ...
        'UniformOutput', false);
    % tolerances by first crossing on the fine grid (the projection is linear in ln(1 + eps))
    dthF = res(n).slope(:) * log(1 + epsFine);
    Fm = m.F;
    dF = sqrt(sum(dthF .* (Fm * dthF), 1));
    bF = abs(exp(dthF(1, :)) - 1);
    res(n).tolBeta1 = firstCrossing(epsFine, bF, 0.01);
    res(n).tolD05 = firstCrossing(epsFine, dF, 0.5);
    res(n).tolD1 = firstCrossing(epsFine, dF, 1.0);
    fprintf('  %-11s sigma [%.4f %.5f], d_W %.3f; per unit ln gamma error: d ln beta %+.4f, d ln R %+.5f; ', ...
        noiseNames{n}, m.sigma, m.dW, res(n).slope);
    fprintf('tolerance |eps|: 1%% beta %s, d 0.5 %s, d 1 %s\n', tolStr(res(n).tolBeta1), tolStr(res(n).tolD05), ...
        tolStr(res(n).tolD1));
end
hIdx = find(strcmp(noiseNames, headline));

%% ---- 3  nonlinear check: exact EMM at gamma_inv, chord Gauss-Newton -----------------------------------
nl = struct('eps', num2cell(cfg.R4.nonlinearCheck(:).'), 'dthetaLin', [], 'dthetaNL', [], 'iters', [], ...
    'diffOverSigma', [], 'relDiff', [], 'exactVsLinearFeature', []);
if RUN_NONLINEAR
    fprintf('\nStep 3: nonlinear model check at eps = %s (exact EMM at gamma_inv; fixed-Jacobian Gauss-Newton)\n', ...
        mat2str(cfg.R4.nonlinearCheck));
    mH = evaluateScenario(setfield(setfield(base, 'Snn', cfg.noise.(headline)), 'verbose', false), cfg); %#ok<SFLD>
    rows = validRows(mH);
    Jv = mH.Jall(rows, :); Bv = mH.blocksAll(:, :, mH.validMask);
    sig = mH.sigma(:);
    fTrue = mH.f0;                                    % noiseless data at the truth (exact T, no Welch bias)
    theta0 = log(p0([1 3])).';
    for k = 1:numel(nl)
        e = nl(k).eps;
        gInv = p0(2) * (1 + e);
        fModel = @(th) modelFeatures([exp(th(1)), gInv, exp(th(2))], sensors, mH.omega, mH.scn, P.cache);
        % exact feature error at the truth versus the linear one
        bExact = fTrue - fModel(theta0);
        bLin = -Jg * log(1 + e);
        nl(k).exactVsLinearFeature = norm(bExact(rows) - bLin(rows)) / norm(bLin(rows));
        qL = fisherFromBlocks(Jv, Bv, bLin(rows));
        nl(k).dthetaLin = qL.dtheta(:);
        th = theta0; it = 0;
        while it < 6
            it = it + 1;
            r = fTrue - fModel(th);
            q = fisherFromBlocks(Jv, Bv, r(rows));
            th = th + q.dtheta(:);
            if all(abs(q.dtheta(:)) < 0.01 * sig), break; end
        end
        nl(k).dthetaNL = th - theta0;
        nl(k).iters = it;
        nl(k).diffOverSigma = (nl(k).dthetaNL - nl(k).dthetaLin) ./ sig;
        nl(k).relDiff = nl(k).dthetaNL ./ nl(k).dthetaLin - 1;      % nonlinear relative to linear
        fprintf(['  eps %+5.0f%%: linear d ln beta %+.5f, d ln R %+.6f | nonlinear %+.5f, %+.6f (%d iterations) | ' ...
            'difference [%+.3f %+.3f] sigma, [%+.1f%% %+.1f%%] relative; exact vs linear feature error %.1f%%\n'], ...
            100 * e, nl(k).dthetaLin, nl(k).dthetaNL, it, nl(k).diffOverSigma, 100 * nl(k).relDiff, ...
            100 * nl(k).exactVsLinearFeature);
    end
end

%% ---- 4  tables ------------------------------------------------------------------------------------------
rowsT = {};
for n = 1:nN
    for k = 1:numel(epsGrid)
        rowsT(end + 1, :) = {noiseNames{n}, epsGrid(k), 100 * res(n).dBeta(k), 100 * res(n).dR(k), ...
            res(n).dSys(k), res(n).dTot(k), res(n).class{k}}; %#ok<AGROW>
    end
end
writeStudyCsv(run, 'sweep', {'noise_case', 'eps_gamma', 'dbeta_pct', 'dR_pct', 'd_sys', 'd_tot', 'class'}, rowsT);
rowsT = cell(nN, 8);
for n = 1:nN
    rowsT(n, :) = {noiseNames{n}, res(n).Snn, res(n).sigma(1), res(n).slope(1), res(n).slope(2), ...
        100 * res(n).tolBeta1, 100 * res(n).tolD05, 100 * res(n).tolD1};
end
writeStudyCsv(run, 'tolerance', {'noise_case', 'Snn', 'sigma_lnbeta', 'dlnbeta_per_lngamma', 'dlnR_per_lngamma', ...
    'tol_beta1pct_pct', 'tol_d05_pct', 'tol_d1_pct'}, rowsT);

%% ---- 5  summary ----------------------------------------------------------------------------------------
L = {};
L{end + 1} = sprintf('Setting: p0, Level 2C (ref s26), N = %d, L = %d, JONSWAP Hs %.2f m, omega_p %.1f rad/s.', ...
    m0.scn.N, m0.scn.L, cfg.sea.Hs, cfg.sea.omegaP);
L{end + 1} = sprintf('Sensitivity per unit ln(gamma) error: d ln beta %+.4f, d ln R %+.5f (%s noise).', ...
    res(hIdx).slope, headline);
L{end + 1} = '(earlier FRF-based estimate for comparison: 0.057 for beta, 0.008 for R; magnitudes)';
L{end + 1} = '';
L{end + 1} = 'DELIVERABLE: required accuracy of gamma (= relative accuracy of the draught) by first crossing';
L{end + 1} = sprintf('  %-12s %10s %12s %12s %12s', 'noise case', 'sigma_lnb', '|db/b|<=1%', 'd_sys<=0.5', 'd_sys<=1');
for n = 1:nN
    L{end + 1} = sprintf('  %-12s %9.2f%% %12s %12s %12s', noiseNames{n}, 100 * res(n).sigma(1), ...
        tolStr(res(n).tolBeta1), tolStr(res(n).tolD05), tolStr(res(n).tolD1)); %#ok<AGROW>
end
L{end + 1} = 'Tolerances are the smaller |eps| of the two signs. gamma ~ draught ~ h (rho_i known): the same';
L{end + 1} = 'relative accuracy applies to the measured draught; note beta ~ h^3 for an independent reference.';
draft = 1e3 * p0(2) * cfg.physics.Hbasin;                % draught at p0, mm
L{end + 1} = sprintf('In absolute terms at p0 (draught gamma H = %.2f mm): |db/b| <= 1%% needs +-%.2f mm (%s);', ...
    draft, draft * res(hIdx).tolBeta1, headline);
L{end + 1} = sprintf('d_sys <= 0.5 needs +-%.3f mm. These are p0 values: other floes have other sensitivities (A1).', ...
    draft * res(hIdx).tolD05);
if RUN_NONLINEAR
    L{end + 1} = '';
    L{end + 1} = 'Nonlinear model check (exact EMM at gamma_inv, fixed-Jacobian Gauss-Newton, headline noise):';
    for k = 1:numel(nl)
        L{end + 1} = sprintf('  eps %+4.0f%%: dbeta/beta %+.3f%% (linear %+.3f%%), nonlinear/linear - 1 = [%+.1f%% %+.1f%%], %d iterations', ...
            100 * nl(k).eps, 100 * (exp(nl(k).dthetaNL(1)) - 1), 100 * (exp(nl(k).dthetaLin(1)) - 1), ...
            100 * nl(k).relDiff, nl(k).iters); %#ok<AGROW>
    end
end
results = struct('epsGrid', epsGrid, 'res', res, 'nonlinear', nl, 'headline', headline, 'Jgamma', Jg, ...
    'omega', m0.omega);
plotR4GammaTolerance(results, cfg, run);         % figure (also redrawable later from the .mat)
finishStudy(run, results, cfg, L);
catch err
    diary('off');
    rethrow(err);
end

%% ================================================================================================
function f = modelFeatures(p, sensors, omega, scn, cacheDir)
% noiseless transmissibility features of the frozen model at p (no derivatives)
fld = emmField(p, sensors, 'ThetaIdx', [], 'NodeSpacing', scn.nodeSpacing, 'CacheDir', cacheDir, ...
    'FRFOptions', scn.frfOptions, 'Verbose', true);
[~, f] = transmissibilityFromField(fld.H(omega), [], scn.layout, scn.ref);
end

function e = firstCrossing(epsGrid, metric, level)
% smallest |eps| at which metric reaches level, over both signs (NaN if never within the grid)
e = NaN;
pos = epsGrid > 0; neg = epsGrid < 0;
c = [];
ip = find(pos & metric >= level, 1, 'first');
if ~isempty(ip), c(end + 1) = interpCross(epsGrid, metric, ip, level); end
in = find(neg & metric >= level, 1, 'last');
if ~isempty(in), c(end + 1) = interpCross(epsGrid, metric, in, level); end
if ~isempty(c), e = min(abs(c)); end
end

function x = interpCross(e, m, i, level)
% linear interpolation towards eps = 0 from the first grid point beyond the level
j = i - sign(e(i)) ;                                  % neighbour closer to zero
j = round(j);
if j < 1 || j > numel(e) || m(j) >= level, x = e(i); return; end
x = e(j) + (level - m(j)) * (e(i) - e(j)) / (m(i) - m(j));
end

function s = tolStr(e)
% NaN: never reached within cfg.R4.toleranceSearch
if isnan(e), s = 'beyond search'; else, s = sprintf('%.2f%%', 100 * e); end
end

function c = classOf(sMax, dTot, condF, singular, cfg)
I = cfg.class.identifiable; M = cfg.class.marginal;
if singular, c = cfg.class.names{3}; return; end
if sMax <= I.maxSigma && dTot <= I.maxBiasDistance && condF <= I.maxCondF, c = cfg.class.names{1};
elseif sMax <= M.maxSigma && dTot <= M.maxBiasDistance && condF <= M.maxCondF, c = cfg.class.names{2};
else, c = cfg.class.names{3};
end
end