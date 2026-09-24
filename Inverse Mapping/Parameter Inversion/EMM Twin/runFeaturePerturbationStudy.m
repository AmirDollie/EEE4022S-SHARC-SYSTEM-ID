%% runFeaturePerturbationStudy.m
% Controlled feature-perturbation / inverse-sensitivity study around the
% solved noise-free EMM twin (runEMMTwinInversion.m).
%
% Question: if the five observed feature frequencies carry relative errors,
% how large are the resulting errors in (beta, gamma, R)?
%
% PARTS
%   A. Load the solved twin; form J_rel = d r / d ln x at the solution from
%      the Jacobian lsqnonlin already returned (no EMM calls); SVD.
%   B. Linear Monte-Carlo sweep (no EMM calls). NOISE MODEL (baseline, stated
%      explicitly): eta_i iid N(0, sigma_f^2), the SAME relative standard
%      deviation on all five features,  f_obs = f_true .* (1 + eta).
%      The residual at the truth is then d = -eta ./ (1 + eta) and the local
%      parameter error is  dlnx = -pinv(J_rel) * d.
%      Cross-checked against the first-order linear-Gaussian benchmark
%      RMS(dlnx_j) = sigma_f * sqrt(diag(pinv(J_rel)*pinv(J_rel)'))_j
%      (first-order because the samples use the exact d = -eta./(1+eta)).
%      The 'outside bounds' fractions are LINEARISED prediction rates under
%      this noise model, not nonlinear inversion failure probabilities.
%   C. Deterministic worst-direction demonstration (no EMM calls): residual
%      d = s*a*u3 at the truth gives exactly dlnx = -s*(a/sigma3)*v3. Not a
%      new result -- the SVD made visible.
%   D. OPTIONAL nonlinear verification (EXPENSIVE, off by default): full
%      solveInverse runs on f_obs = f_true ./ (1 + s*a*u3) (residual at the
%      truth exactly s*a*u3), warm-started at x_true, to see where the linear
%      picture breaks down. A solution touching a bound is reported as
%      'bound-limited', never as a recovered value.
%
% Linear predictions are in ln x; they are mapped back with
% x = x_true .* exp(dlnx) before comparing with the bounds, and relative
% parameter errors are exp(dlnx) - 1.
%
% OUTPUT (in Parameter Inversion/Results/)
%   FeaturePerturbation_<timestamp>.mat / .txt / .png / .fig
%
% Lives in Inverse Mapping/Parameter Inversion/EMM Twin/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));         % solveInverse.m, inverseResidual.m
addpath(fullfile(thisDir, '..', '..'));   % modalFeatureVector.m, predictFeatureAnchors.m

clearvars -except thisDir
close all, clc

%% ============================ CONFIGURATION ==============================
resultsDir = fullfile(thisDir, '..', 'Results');
TWIN_FILE  = fullfile(resultsDir, 'EMMTwin_noiseFree_20260924_094431.mat');  % start 1
TWIN_START = 1;                         % which solved start's Jacobian to use

NOISE_LEVELS = [1e-5 3e-5 1e-4 3e-4 1e-3 3e-3 1e-2];   % sigma_f (relative)
N_MC = 1e4;                             % Monte-Carlo samples per level
RNG_SEED = 1;                           % reproducible sampling

U3_AMPLITUDES = NOISE_LEVELS;           % part C amplitudes a

RUN_NONLINEAR = true;                   % part D: ~20-45 min per case
% [a, sign]: residual at truth = sign*a*u3. First the small-perturbation
% regime (1e-4, linear picture expected to hold), then the transition
% (3e-4; the -1 case is predicted at gamma ~ 0.857, just above lb = 0.85).
% The 1e-3 pair is predicted outside the box on both sides (bound-limited
% demonstration); add it later if wanted:  1e-3 +1;  1e-3 -1
NONLINEAR_CASES = [1e-4  +1
                   1e-4  -1
                   3e-4  +1
                   3e-4  -1];
NL_MAX_FEVALS = 70;                     % ~10 iterations
NL_MAX_ITERS = 12;
BOUND_TOL = 1e-4;                       % fraction of box width counted as "on a bound"

modelFile = fullfile(thisDir, '..', '..', 'anchorPredictionModel.mat');
EXPECTED_ERROR_IDS = {'modalFeatureVector:invalidFeatureFit'};
paramNames = {'beta', 'gamma', 'R'};
%% =========================================================================

stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['FeaturePerturbation_', stamp]);
diary([baseName, '.txt']);
try   % closes the log even if the run errors

%% A. Twin, Jacobian, SVD
T = load(TWIN_FILE);
twin = T.twin;
res = T.results{TWIN_START};
assert(~isempty(res) && strcmp(res.status, 'solved'), 'runFeaturePerturbationStudy:notSolved', ...
    'Start %d in %s is not a solved start.', TWIN_START, TWIN_FILE);
xTrue = twin.xTrue; fTrue = twin.fObs(:); pRef = twin.pRef; lb = twin.lb; ub = twin.ub;
m = numel(fTrue); n = numel(xTrue);

% res.jacobian = dr/dx at xHat (rows already divided by fTrue);
% d/d(ln x_j) = x_j * d/dx_j.
Jrel = res.jacobian * diag(res.xHat);
[U, S, V] = svd(Jrel, 'econ');
sigma = diag(S);
if V(2, end) > 0, U(:, end) = -U(:, end); V(:, end) = -V(:, end); end   % sign convention: v3(gamma) < 0
u3 = U(:, end); v3 = V(:, end); sigma3 = sigma(end);
P = pinv(Jrel);
Cunit = P * P.';                        % first-order covariance per unit sigma_f^2 (= inv(J'J), without squaring kappa)
rmsUnit = sqrt(diag(Cunit));            % analytic RMS(dlnx_j) / sigma_f

fprintf('runFeaturePerturbationStudy  (%s)\n', stamp);
fprintf('twin: %s, start %d  (xTrue = %s, xHat = %s)\n', TWIN_FILE, TWIN_START, ...
    mat2str(xTrue.', 4), mat2str(res.xHat.', 8));
fprintf('bounds (scaled): lb = %s   ub = %s\n\n', mat2str(lb.', 4), mat2str(ub.', 4));
fprintf('A. J_rel = dr/dln x at the solution\n');
fprintf('   sigma = [%s]   kappa = %.0f\n', sprintf(' %.4g', sigma), sigma(1) / sigma3);
fprintf('   v3 (weak parameter direction, beta gamma R) = [%s]\n', sprintf(' %+.4f', v3));
fprintf('   u3 (weak feature direction)                 = [%s]\n', sprintf(' %+.4f', u3));
fprintf('   analytic RMS(dlnx)/sigma_f (beta gamma R)    = [%s]\n\n', sprintf(' %.4g', rmsUnit));

%% B. Linear Monte-Carlo sweep
rng(RNG_SEED);
nL = numel(NOISE_LEVELS);
mc = struct('sigmaF', num2cell(NOISE_LEVELS), 'rmsDlnx', [], 'rmsAnalytic', [], 'mcOverAnalytic', [], ...
    'medianAbsRelErr', [], 'p95AbsRelErr', [], 'fracOutParam', [], 'fracOutAny', []);
fprintf('B. Linear Monte Carlo, N = %d per level, eta iid N(0, sigma_f^2) on all %d features\n', N_MC, m);
fprintf('   %-9s | %-30s | %-22s | %-30s | %-24s | %s\n', 'sigma_f', 'RMS error in ln x (b g R)', ...
    'MC/analytic (b g R)', '95th pct |rel. err| (b g R)', 'outside bounds (b g R)', 'any');
for L = 1:nL
    s = NOISE_LEVELS(L);
    eta = s * randn(m, N_MC);
    d = -eta ./ (1 + eta);
    dLogX = -P * d;                                   % n x N
    relErr = exp(dLogX) - 1;
    xLin = xTrue .* exp(dLogX);
    outP = xLin < lb | xLin > ub;

    mc(L).rmsDlnx = sqrt(mean(dLogX.^2, 2));
    mc(L).rmsAnalytic = s * rmsUnit;
    mc(L).mcOverAnalytic = mc(L).rmsDlnx ./ mc(L).rmsAnalytic;
    mc(L).medianAbsRelErr = median(abs(relErr), 2);
    mc(L).p95AbsRelErr = pct95(abs(relErr));
    mc(L).fracOutParam = mean(outP, 2);
    mc(L).fracOutAny = mean(any(outP, 1));
    fprintf('   %-9.0e | %-30s | %-22s | %-30s | %-24s | %.3f\n', s, ...
        sprintf('%.2e ', mc(L).rmsDlnx), sprintf('%.3f ', mc(L).mcOverAnalytic), ...
        sprintf('%.2e ', mc(L).p95AbsRelErr), sprintf('%.3f ', mc(L).fracOutParam), mc(L).fracOutAny);
end
worstRatio = max(max(abs([mc.mcOverAnalytic] - 1)));
fprintf('   MC vs first-order benchmark RMS: worst deviation %.1f%% (sampling error ~%.1f%%; the top level also carries the\n', ...
    100 * worstRatio, 100 / sqrt(2 * N_MC));
fprintf('   small eta/(1+eta) nonlinearity). RMS error ratio gamma/beta = %.0f, gamma/R = %.0f.\n\n', ...
    rmsUnit(2) / rmsUnit(1), rmsUnit(2) / rmsUnit(3));

%% C. Worst-direction demonstration
fprintf('C. Residual d = s*a*u3 at the truth -> dlnx = -s*(a/sigma3)*v3  (linear; x = xTrue.*exp(dlnx))\n');
fprintf('   %-8s %4s | %-32s | %-26s | %s\n', 'a', 's', 'predicted x (b g R)', 'rel. error (b g R)', 'inside bounds?');
u3tab = struct('a', {}, 'sign', {}, 'dLogX', {}, 'xLin', {}, 'relErr', {}, 'inBounds', {});
for a = U3_AMPLITUDES
    for sgn = [+1, -1]
        dLogX = -P * (sgn * a * u3);
        assert(norm(dLogX - (-sgn * a / sigma3 * v3)) <= 1e-9 * norm(dLogX), ...
            'runFeaturePerturbationStudy:svdMismatch', 'pinv and SVD predictions disagree');
        xLin = xTrue .* exp(dLogX);
        inB = all(xLin >= lb & xLin <= ub);
        u3tab(end+1) = struct('a', a, 'sign', sgn, 'dLogX', dLogX, 'xLin', xLin, ...
            'relErr', exp(dLogX) - 1, 'inBounds', inB); %#ok<AGROW>
        fprintf('   %-8.0e %+4d | %-32s | %-26s | %s\n', a, sgn, sprintf('%.5f ', xLin), ...
            sprintf('%+.2e ', exp(dLogX) - 1), mat2str(inB));
    end
end
fprintf('\n');

%% Figure: RMS relative parameter error vs feature noise (linear MC)
fig = figure('Color', 'w', 'Position', [100 100 620 440]);
ax = axes(fig); hold(ax, 'on');
colours = {[42 120 214]/255, [235 104 52]/255, [27 175 122]/255};   % validated categorical slots 1-3 (#2a78d6 #eb6834 #1baf7a)
markers = {'o', 's', '^'};
R = [mc.rmsDlnx];                                   % n x nL
for j = 1:n
    plot(ax, NOISE_LEVELS, R(j, :), '-', 'Color', colours{j}, 'LineWidth', 2, 'Marker', markers{j}, ...
        'MarkerSize', 8, 'MarkerFaceColor', colours{j}, 'MarkerEdgeColor', 'w', 'DisplayName', paramNames{j});
end
plot(ax, NOISE_LEVELS, NOISE_LEVELS, ':', 'Color', [0.55 0.55 0.52], 'LineWidth', 1, ...
    'DisplayName', 'error = \sigma_f');
set(ax, 'XScale', 'log', 'YScale', 'log', 'Box', 'off', 'TickDir', 'out', ...
    'XColor', [0.35 0.35 0.33], 'YColor', [0.35 0.35 0.33], 'FontSize', 10);
grid(ax, 'on'); set(ax, 'GridColor', [0.85 0.85 0.83], 'GridAlpha', 1);
xlabel(ax, 'relative feature noise \sigma_f');
ylabel(ax, 'RMS error in ln(parameter)');
title(ax, sprintf('Linearised parameter error vs feature noise (\\kappa = %.0f)', sigma(1) / sigma3), ...
    'FontWeight', 'normal');
lgd = legend(ax, [paramNames, {'error = \sigma_f'}]);
set(lgd, 'Location', 'northwest'); legend(ax, 'boxoff');
saveas(fig, [baseName, '.png']);
savefig(fig, [baseName, '.fig']);

study = struct('stamp', stamp, 'twinFile', TWIN_FILE, 'twinStart', TWIN_START, 'xTrue', xTrue, ...
    'xHat', res.xHat, 'fTrue', fTrue, 'lb', lb, 'ub', ub, 'Jrel', Jrel, 'sigma', sigma, 'U', U, 'V', V, ...
    'rmsUnit', rmsUnit, 'noiseLevels', NOISE_LEVELS, 'nMC', N_MC, 'rngSeed', RNG_SEED, ...
    'noiseModel', 'eta_i iid N(0,sigma_f^2), same relative sd on all features; f_obs = f_true.*(1+eta)');
nonlinear = struct('a', {}, 'sign', {}, 'status', {}, 'boundLimited', {}, 'xHat', {}, ...
    'dLogXnl', {}, 'dLogXlin', {}, 'gammaRatio', {}, 'resnorm', {}, 'funcCount', {}, 'runtimeMin', {}, 'result', {});
save([baseName, '.mat'], 'study', 'mc', 'u3tab', 'nonlinear');

%% D. Optional nonlinear verification
if RUN_NONLINEAR
    model = load(modelFile);
    assert(isequal(model.p0, twin.model.p0) && isequal(model.Js, twin.model.Js), ...
        'runFeaturePerturbationStudy:modelChanged', ...
        'anchorPredictionModel.mat differs from the one the twin was run with.');
    fwd = @(p) modalFeatureVector(p(1), p(2), p(3), 'Strict', true, ...
        'AnchorMode', 'predicted', 'AnchorModel', modelFile);
    fprintf('D. Nonlinear verification (warm start at xTrue), %d case(s)\n', size(NONLINEAR_CASES, 1));
    for c = 1:size(NONLINEAR_CASES, 1)
        a = NONLINEAR_CASES(c, 1); sgn = NONLINEAR_CASES(c, 2);
        fObsPert = fTrue ./ (1 + sgn * a * u3);       % residual at the truth = sgn*a*u3 exactly
        dLogXlin = -sgn * a / sigma3 * v3;
        fprintf('\n---- case %d: a = %.0e, sign %+d; linear prediction x = %s ----\n', ...
            c, a, sgn, mat2str((xTrue .* exp(dLogXlin)).', 6));
        r = solveInverse(fwd, fObsPert, pRef, xTrue, lb, ub, ...
            'ExpectedErrorIDs', EXPECTED_ERROR_IDS, 'MaxFunctionEvaluations', NL_MAX_FEVALS, ...
            'MaxIterations', NL_MAX_ITERS, 'Display', 'iter');
        r = rmfield(r, 'solverOptions');
        w = ub - lb;
        onBound = (r.xHat - lb) ./ w < BOUND_TOL | (ub - r.xHat) ./ w < BOUND_TOL;
        boundLimited = ~strcmp(r.status, 'leftValidRegion') && any(onBound);   % on a bound, whatever the stop reason
        dLogXnl = log(r.xHat ./ xTrue);
        e = struct('a', a, 'sign', sgn, 'status', r.status, 'boundLimited', boundLimited, 'xHat', r.xHat, ...
            'dLogXnl', dLogXnl, 'dLogXlin', dLogXlin, 'gammaRatio', dLogXnl(2) / dLogXlin(2), ...
            'resnorm', r.resnorm, 'funcCount', r.funcCount, 'runtimeMin', r.runtimeSec / 60, 'result', r);
        nonlinear(end+1) = e; %#ok<AGROW>
        if boundLimited
            verdict = sprintf('BOUND-LIMITED (%s on a bound; status %s): not a recovered value', ...
                strjoin(paramNames(onBound), ', '), r.status);
        elseif ~strcmp(r.status, 'solved')
            verdict = sprintf('%s (interior, not converged): not a recovered value', r.status);
        else
            verdict = sprintf('interior; gamma error nonlinear/linear = %.3f', e.gammaRatio);
        end
        fprintf('  xHat = %s   resnorm %.2e   %d calls, %.1f min\n  -> %s\n', mat2str(r.xHat.', 6), ...
            r.resnorm, r.funcCount, r.runtimeSec / 60, verdict);
        save([baseName, '.mat'], 'study', 'mc', 'u3tab', 'nonlinear');   % after every case
    end
else
    fprintf('D. Nonlinear verification skipped (RUN_NONLINEAR = false).\n');
end

fprintf('\nSaved: %s.mat / .png / .fig\n', baseName);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ------------------------------------------------------------------------
function q = pct95(A)
% 95th percentile of each row (nearest-rank; no Statistics Toolbox needed).
As = sort(A, 2);
q = As(:, max(1, ceil(0.95 * size(As, 2))));
end