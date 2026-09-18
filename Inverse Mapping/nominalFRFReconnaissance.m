% Step 1 of the identifiability gameplan: sweep broadly at the nominal
% (beta0, gamma0, R0) using a deliberately broad search interval, without
% assuming resonance locations, to identify resonance peaks VISIBLE IN
% THE NOMINAL FORCED RESPONSE before anything is tracked or perturbed.
% Note: a forced FRF cannot guarantee every natural mode of the floe is
% excited or observable at the chosen evaluation points -- this finds
% what is visible here, not a provably exhaustive list of every mode.
%
% METHOD: a 400-point broad reconnaissance sweep, followed by LOCAL
% high-resolution refinement of every candidate the broad sweep already
% detected. This is a peak-LOCATION refinement check, not a global
% grid-convergence check -- it confirms the coarse grid did not
% MISLOCATE a detected peak, but cannot discover a genuinely narrow
% resonance that fell entirely between two coarse samples and was never
% detected as a candidate in the first place.
%
% SENSOR LOCATIONS: not subject to the 4-sensor deployment constraint --
% multi-radius, staggered-angle grid to guard against angular and radial
% nodes.
%
% METRIC: the sampled spatial response norm ||Phi(alpha)||_2.
%
% PEAK DETECTION: no prominence filter on the coarse sweep -- every local
% maximum is detected and reported. RANKING uses the coarse, full-sweep
% prominence (promsCoarse) throughout, since a local refinement window's
% own prominence is only comparable to other features WITHIN that same
% window, not across different candidates' differently-shaped windows.

%% Add paths BEFORE starting the pool, so new workers inherit them
% directly -- no asynchronous broadcast, no race condition.
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', 'Forward Model', 'Animation'));
addpath(fullfile(thisDir, '..', 'SSI'));

if isempty(gcp('nocreate'))
    parpool;
else
    % A pool already existed before this script's addpath calls ran, so
    % its workers may not have the updated path. Broadcast and WAIT for
    % it to finish before any parfor below can start.
    currentPath = path;
    futures = parfevalOnAll(@() addpath(currentPath), 0);
    wait(futures);
end

% Do not use clear all here: it would clear the parallel pool handle,
% forcing a fresh pool to spin up. Use clearvars if variables need
% clearing between runs.
close all, clc

H = 1.88; nu = 0.3;
g = 9.81;
beta0 = 4.6985e-5; gamma0 = 1.4548e-3; R0 = 0.3830;
M = 50; P = 10; N = 10;

%% Multi-radius, staggered-angle evaluation grid
rFrac = [0.3, 0.6, 0.9];
thetaRing = linspace(0, 2*pi, 9); thetaRing(end) = [];

evalLocations = [0, 0];
for ri = 1:length(rFrac)
    stagger = (ri-1) * (pi/16);
    for th = thetaRing
        evalLocations(end+1,:) = [rFrac(ri)*R0, mod(th + stagger, 2*pi)]; %#ok<AGROW>
    end
end
nEval = size(evalLocations,1);
fprintf('Reconnaissance evaluation grid: %d points (centre + %d radii x %d angles, staggered)\n', ...
    nEval, length(rFrac), length(thetaRing));

%% Coarse sweep (400 points), parallelised
fprintf('\n=== Coarse sweep (400 points) ===\n');
[alphaCoarse, responseNormCoarse] = sweepResponseNorm(400, evalLocations, beta0, gamma0, R0, nu, M, P, N);

%% Peak detection on the coarse grid: NO prominence filter
[pksCoarse, locsCoarse, ~, promsCoarse] = findpeaks(responseNormCoarse, alphaCoarse);
fprintf('Coarse grid: %d local maxima detected at alpha = %s\n', length(locsCoarse), mat2str(locsCoarse, 4));

%% Local refinement around each coarse candidate, parallelised per window.
% Preallocated with NaN so a failed refinement at one index cannot shift
% or corrupt the mapping between coarse candidate k and its refined
% counterpart.
fprintf('\n=== Local refinement of each coarse-detected candidate ===\n');
refineHalfWidth = 0.3;
refinePointsPerWindow = 60;

locsFine  = nan(size(locsCoarse));
pksFine   = nan(size(locsCoarse));
promsFine = nan(size(locsCoarse));   % LOCAL-window prominence, diagnostic only -- NOT used for ranking

for k = 1:length(locsCoarse)
    windowLo = max(0.5, locsCoarse(k) - refineHalfWidth);
    windowHi = min(20, locsCoarse(k) + refineHalfWidth);
    alphaWindow = linspace(windowLo, windowHi, refinePointsPerWindow);

    responseWindow = zeros(size(alphaWindow));
    parfor a = 1:length(alphaWindow)
        binData = precomputeDeflectionData(alphaWindow(a), beta0, gamma0, R0, nu, M, P, N);
        Phi = zeros(nEval, 1);
        for e = 1:nEval
            Phi(e) = evaluateDeflection(binData, evalLocations(e,1), evalLocations(e,2));
        end
        responseWindow(a) = norm(Phi);
    end

    [pk, loc, ~, prom] = findpeaks(responseWindow, alphaWindow);
    if ~isempty(pk)
        % Select the peak NEAREST the original coarse location, not the
        % largest in the window -- prevents branch-jumping onto a
        % neighbouring resonance if two peaks sit close together.
        [~, best] = min(abs(loc - locsCoarse(k)));
        locsFine(k) = loc(best);
        pksFine(k) = pk(best);
        promsFine(k) = prom(best);
    end
    fprintf('  refined candidate %d/%d (coarse alpha=%.3f)\n', k, length(locsCoarse), locsCoarse(k));
end

fprintf('\nCoarse locations:  %s\n', mat2str(locsCoarse, 4));
fprintf('Refined locations: %s\n', mat2str(locsFine, 4));

%% Plot: coarse curve, with coarse and refined peak locations both marked
figure;
plot(alphaCoarse, responseNormCoarse, 'LineWidth', 1.2); hold on;
plot(locsCoarse, pksCoarse, 'x', 'MarkerSize', 8, 'LineWidth', 1.2);
plot(locsFine, pksFine, 'ro', 'MarkerSize', 8, 'LineWidth', 1.2);
xlabel('\alpha (non-dimensional frequency)');
ylabel('||\Phi(\alpha)||_2 (sampled spatial response norm)');
title('Nominal FRF reconnaissance: coarse curve, coarse peaks (x), refined peaks (o)');
legend('coarse sweep (400 pts)', 'coarse-detected peaks', 'refined peaks');
grid on;

%% Peak-location refinement check (NOT a global grid-convergence claim),
% with an explicit guard against the zero-peaks-detected case
validRefined = ~isnan(locsFine);
if ~isempty(locsFine) && all(validRefined) && max(abs(locsFine - locsCoarse)) < 0.05
    fprintf('\nPeak locations confirmed stable under local refinement (within 0.05 in alpha).\n');
    fprintf('This confirms the coarse grid did not MISLOCATE any detected peak. It does NOT\n');
    fprintf('confirm the coarse grid detected every resonance that genuinely exists.\n');
elseif isempty(locsFine)
    fprintf('\nNo local maxima were detected at all in the coarse sweep -- either this floe\n');
    fprintf('genuinely has no resonance visible in this evaluation grid/alpha range, or something\n');
    fprintf('upstream (evaluation points, alpha range, or the forward model call itself) needs\n');
    fprintf('checking before trusting that conclusion.\n');
else
    fprintf('\nAt least one candidate shifted under refinement, or failed to refine at all --\n');
    fprintf('inspect the plot and per-candidate output above before proceeding.\n');
end

%% Report every coarse-detected candidate, ranked by the COARSE, full-
% sweep prominence (the only prominence value comparable across all
% candidates), with physical frequency and truncation-trust warning
fprintf('\nAll detected local maxima, ranked by coarse-sweep prominence:\n');
[~, rankOrder] = sort(promsCoarse, 'descend');
for idx = 1:length(rankOrder)
    k = rankOrder(idx);
    alphaToReport = locsFine(k);
    if isnan(alphaToReport)
        alphaToReport = locsCoarse(k);   % fall back to coarse if refinement failed
    end
    omega_k = sqrt(alphaToReport*g/H);
    trustNote = '';
    if alphaToReport > 13.85
        trustNote = '  [ABOVE validated truncation ceiling alpha~13.85 -- re-check M,P,N before trusting]';
    end
    fprintf('  alpha=%.4f  ->  omega=%.4f rad/s  (coarse prominence=%.4f)%s\n', ...
        alphaToReport, omega_k, promsCoarse(k), trustNote);
end

%% Diagnostic: does removing the known 1/alpha kinematic prefactor reveal
% resonant structure hidden underneath it? This directly tests eq. (A.2)'s
% own explicit 1/alpha scaling as the likely cause of the featureless
% monotonic decay, using data ALREADY computed -- no new forward model
% calls needed.
compensatedResponse = alphaCoarse(:) .* responseNormCoarse(:);   % undo the 1/alpha factor

figure;
plot(alphaCoarse, compensatedResponse, 'LineWidth', 1.2);
xlabel('\alpha (non-dimensional frequency)');
ylabel('\alpha \times ||\Phi(\alpha)||_2 (1/\alpha kinematic factor removed)');
title('Reconnaissance with kinematic prefactor compensated out');
grid on;

[pksComp, locsComp, ~, promsComp] = findpeaks(compensatedResponse, alphaCoarse);
fprintf('Compensated sweep: %d local maxima detected at alpha = %s\n', length(locsComp), mat2str(locsComp, 4));

%% Diagnostic: phase behavior, not amplitude, is the unambiguous test for
% genuine resonant structure. The 1/alpha kinematic factor contributes a
% FIXED, CONSTANT -90 degree phase shift, never a phase SWING -- so any
% observed swing in arg(eta) as alpha varies is unambiguous evidence of a
% real pole crossing, completely independent of the amplitude-scaling
% question raised (and corrected) above.
%
% Requires COMPLEX eta values, one representative sensor (not the
% ||Phi||_2 norm, which discards phase entirely).

referenceLocation = evalLocations(1,:);   % e.g. the centre point already in the grid

etaPhase = zeros(size(alphaCoarse));
for a = 1:length(alphaCoarse)
    binData = precomputeDeflectionData(alphaCoarse(a), beta0, gamma0, R0, nu, M, P, N);
    etaPhase(a) = angle(evaluateDeflection(binData, referenceLocation(1), referenceLocation(2)));
end

figure;
plot(alphaCoarse, unwrap(etaPhase)*180/pi, 'LineWidth', 1.2);
xlabel('\alpha (non-dimensional frequency)');
ylabel('arg(\eta) (degrees, unwrapped)');
title('Phase diagnostic: a genuine resonance shows a ~180 deg swing; the 1/\alpha term alone cannot');
grid on;

%% Free diagnostic (step 1 of 4): does removing the known 1/alpha
% kinematic factor reveal any underlying structure at all? Not treated
% as a resonance detector on its own -- purely a check for whether that
% explicit scaling term is obscuring something worth investigating
% further. Costs nothing: uses data already computed.
compensatedResponse = alphaCoarse(:) .* responseNormCoarse(:);

figure;
plot(alphaCoarse, compensatedResponse, 'LineWidth', 1.2);
grid on
xlabel('\alpha');
ylabel('\alpha \times ||\Phi(\alpha)||_2 (kinematic factor removed)');
title('Diagnostic: is any structure hidden under the 1/\alpha term?');

%% Check what's actually available inside binData before assuming a new
% expensive sweep is needed
binDataSample = precomputeDeflectionData(2.0, beta0, gamma0, R0, nu, M, P, N);
disp(fieldnames(binDataSample))

%% Confirm what AfullMat actually is before trusting it as "the" system
% matrix -- check its size/structure against what's expected from the
% EMM matching system, rather than assuming based on the field name alone.
binDataSample = precomputeDeflectionData(2.0, beta0, gamma0, R0, nu, M, P, N);
fprintf('Size of AfullMat: %s\n', mat2str(size(binDataSample.AfullMat)));
fprintf('Is it square? %s\n', mat2str(size(binDataSample.AfullMat,1) == size(binDataSample.AfullMat,2)));
fprintf('Condition number at alpha=2.0: %.4e\n', cond(binDataSample.AfullMat));
fprintf('Smallest singular value at alpha=2.0: %.4e\n', min(svd(binDataSample.AfullMat)));
%% Local function, placed at the end of the script per MATLAB compatibility
function [alphaGrid, responseNorm] = sweepResponseNorm(nPoints, evalLocations, beta0, gamma0, R0, nu, M, P, N)
    alphaGrid = linspace(0.5, 20, nPoints);
    nEval = size(evalLocations,1);
    responseNorm = zeros(nPoints, 1);

    ticStart = tic;
    parfor a = 1:nPoints
        binData = precomputeDeflectionData(alphaGrid(a), beta0, gamma0, R0, nu, M, P, N);
        Phi = zeros(nEval, 1);
        for e = 1:nEval
            Phi(e) = evaluateDeflection(binData, evalLocations(e,1), evalLocations(e,2));
        end
        responseNorm(a) = norm(Phi);
    end
    fprintf('  Sweep of %d points completed in %.1fs\n', nPoints, toc(ticStart));
end
