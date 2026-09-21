%% modalFeatureRefinement.m
% Local, targeted refinement of the handful of modal-response features
% spotted in the 35-point modalRAOReconnaissance.m sweep: three broad
% maxima and two suspected modal-coefficient zero/sign-change crossings
% (one visually clear, one candidate right at the edge of the validated
% truncation range). This is NOT another global sweep -- five small
% local windows only, at a much finer omega step than the coarse sweep.
%
% For each window this script:
%   (1) resweeps that omega range at a finer, tiered step
%   (2) extracts the feature location (argmax omega, or zero-crossing
%       omega + the phase change either side of it)
%   (3) ALSO keeps the full locally-sampled complex A_{n,j}(omega_k)
% -- because whether the eventual Jacobian should differentiate a
% feature LOCATION or a sampled COMPLEX VALUE at fixed omega is exactly
% the open question this script's output is meant to help answer: a
% broad, flat maximum is a poor finite-difference observable if its
% location isn't stable under refinement, and Re/Im samples at fixed
% omega are usually much smoother functions of the parameters.
%
% Extra check, A_{0,0} window only: that candidate zero sits right next
% to the validated truncation ceiling (alpha~13.85 for M=50,P=10,N=10 --
% see nominalFRFReconnaissance.m), so this window is ALSO re-solved at a
% bumped (M,P,N) -- a TRUNCATION-convergence check, not an M-only check,
% since all three orders move together -- to see whether the
% zero-crossing location is a real feature or a truncation-order
% artifact. Both truncations are evaluated on the EXACT SAME omega grid:
% comparing a result on a different grid to the base result would
% conflate grid quantisation with a genuine physical shift (an earlier
% version of this check did exactly that -- Delta-omega=0.1 for the
% check grid vs 0.05 for the base grid means the check grid cannot even
% land on the base grid's own zero location, so any "shift" it reports
% is partly just which grid point happened to be nearest). Toggle
% CHECK_TRUNCATION_CONVERGENCE below to skip this and save time.
%
% SELECTIVE RERUN: if modalFeatureRefinementResults.mat from a previous
% run already exists next to this script, a window is only actually
% re-solved if (a) its label is listed in WINDOWS_TO_RUN below, or (b)
% its omegaLo/omegaHi/dOmega changed since that previous run (caught
% automatically, so a spec edit can never silently reuse stale data even
% if you forget to add it to WINDOWS_TO_RUN) -- everything else is
% loaded from the previous save. Set WINDOWS_TO_RUN = {} to force
% recomputing every window from scratch.
%
% Depends on modalRAOReconnaissance.m having been run already in this
% folder: loads its saved parameters and the coarse curve from
% modalRAOReconnaissanceResults.mat, so nothing here can silently drift
% out of sync with the coarse sweep.
%
% *** RUNTIME NOTE *** each point costs one full EMM solve at the real
% nominal truncation orders (same cost as modalRAOReconnaissance.m's own
% sweep points). Only windows actually being (re)run pay this cost --
% see SELECTIVE RERUN above.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', 'Forward Model', 'Animation'));

clearvars -except thisDir
close all, clc

% Runs whenever true, regardless of whether the A_{0,0} window itself was
% freshly computed or reused from cache above (it's an independent,
% separate re-solve) -- set false too on a run where you're only
% touching other windows and don't need to redo it.
CHECK_TRUNCATION_CONVERGENCE = true;

% Which windows to force-recompute this run (by label). Everything else
% is reused from a previous modalFeatureRefinementResults.mat if present
% and its spec hasn't changed -- see SELECTIVE RERUN in the header.
% Cleanup-run default: only the two windows that actually needed fixing
% (A_{2,0} was widened; A_{0,0}'s truncation-convergence grid was fixed).
WINDOWS_TO_RUN = {'A_{2,0}', 'A_{0,0}'};

%% Load the coarse sweep's parameters and results -- single source of
% truth for (beta0,gamma0,R0,nu,M,P,N,H,g,modeList), so this script
% cannot silently use different physics than the sweep it is refining.
coarseFile = fullfile(thisDir, 'modalRAOReconnaissanceResults.mat');
if ~isfile(coarseFile)
    error('modalFeatureRefinement:missingCoarseSweep', ...
        'modalRAOReconnaissanceResults.mat not found next to this script -- run modalRAOReconnaissance.m first.');
end
S = load(coarseFile);
beta0 = S.beta0; gamma0 = S.gamma0; R0 = S.R0; nu = S.nu;
M = S.M; P = S.P; N = S.N; H = S.H; g = S.g;
modeList = S.modeList;
NTHETA = 90;

%% Refinement windows. 'label' is the plain A_{n,j} LaTeX tag (for axis
% labels/titles); 'row' indexes modeList; 'kind' is 'max' or 'zero'.
% Tiered dOmega: coarser for the three broad maxima (what's being tested
% there is whether the argmax location is STABLE, not sub-percent
% precision), fine for the two zero-crossings (location IS the point).
% A_{2,0} widened from the original 4.5-5.5 (Delta-omega=0.15): the
% first refinement pass showed |A_{2,0}| still rising at omega=5.5 (its
% argmax landed on the window edge), so the true maximum lies further
% right. Widened to 5.3-6.5 at a finer Delta-omega=0.1, per the alpha=10
% bridge result (omega~7.22) already showing |A_{2,0}| declining again
% by then, so the true maximum is bracketed between the old edge and
% there. If this run's argmax STILL lands on the omega=6.5 edge, widen
% again towards ~7.2 and rerun just this window.
featureSpecs = struct( ...
    'label',  {'A_{2,0}', 'A_{0,1}', 'A_{1,1}', 'A_{2,1}', 'A_{0,0}'}, ...
    'row',    {3,         6,         7,         8,         1}, ...
    'kind',   {'max',     'max',     'max',     'zero',    'zero'}, ...
    'omegaLo',{5.3,       5.2,       6.0,       4.5,       7.3}, ...
    'omegaHi',{6.5,       6.2,       7.2,       5.5,       8.236}, ...
    'dOmega', {0.1,       0.15,      0.15,      0.05,      0.05});
nWindows = numel(featureSpecs);

%% Load any previous refinement results, for the selective-rerun logic
prevResultsFile = fullfile(thisDir, 'modalFeatureRefinementResults.mat');
havePrev = isfile(prevResultsFile);
if havePrev
    Prev = load(prevResultsFile);
    fprintf('Found a previous modalFeatureRefinementResults.mat -- windows not in\n');
    fprintf('WINDOWS_TO_RUN (and whose spec is unchanged) will be reused from it.\n\n');
end

%% Refine each window (or reuse a previous result for it)
refined = repmat(struct('label','', 'row',0, 'kind','', 'omega',[], 'A',[], ...
    'featureOmega',NaN, 'featureValue',NaN, 'atEdge',false, 'phaseFlipDeg',NaN), 1, nWindows);

for w = 1:nWindows
    spec = featureSpecs(w);

    % Decide whether this window can be reused as-is: only if a previous
    % result for this exact label+kind exists AND its omegaLo/omegaHi/
    % dOmega are identical to the current spec -- a spec edit (e.g.
    % widening a window) is therefore caught automatically even if the
    % label was left out of WINDOWS_TO_RUN by mistake.
    canReuse = false;
    if havePrev && ~ismember(spec.label, WINDOWS_TO_RUN)
        prevSpecIdx = find(strcmp({Prev.featureSpecs.label}, spec.label) & ...
                            strcmp({Prev.featureSpecs.kind}, spec.kind), 1);
        if ~isempty(prevSpecIdx)
            pSpec = Prev.featureSpecs(prevSpecIdx);
            canReuse = (pSpec.omegaLo == spec.omegaLo) && (pSpec.omegaHi == spec.omegaHi) && (pSpec.dOmega == spec.dOmega);
        end
    end

    if canReuse
        refined(w) = Prev.refined(prevSpecIdx);
        fprintf('=== Window %d/%d: %s (%s) -- REUSED from previous run (spec unchanged) ===\n', ...
            w, nWindows, spec.label, spec.kind);
        fprintf('  (previous result: feature omega = %.4f, |A| there = %.4e)\n', ...
            refined(w).featureOmega, refined(w).featureValue);
        continue;
    end

    omegaWin = spec.omegaLo:spec.dOmega:spec.omegaHi;
    if omegaWin(end) < spec.omegaHi - 1e-9
        omegaWin(end+1) = spec.omegaHi; %#ok<AGROW>
    end
    nW = numel(omegaWin);
    Awin = complex(zeros(1, nW));

    fprintf('=== Window %d/%d: %s (%s), omega = %.3f..%.3f, dOmega=%.3f, %d points ===\n', ...
        w, nWindows, spec.label, spec.kind, spec.omegaLo, spec.omegaHi, spec.dOmega, nW);
    ticW = tic;
    for k = 1:nW
        alpha = H*omegaWin(k)^2/g;
        data = precomputeDeflectionData(alpha, beta0, gamma0, R0, nu, M, P, N);
        wFun = @(r,theta) evaluateDeflection(data, r, theta);
        Ak = projectEMMOntoPlateModes(wFun, modeList(spec.row,:), R0, nu, 'NTheta', NTHETA);
        Awin(k) = Ak;
        fprintf('  [%2d/%2d] omega=%.3f done (%.1fs elapsed)\n', k, nW, omegaWin(k), toc(ticW));
    end

    refined(w).label = spec.label;
    refined(w).row = spec.row;
    refined(w).kind = spec.kind;
    refined(w).omega = omegaWin;
    refined(w).A = Awin;

    if strcmp(spec.kind, 'max')
        [maxVal, idx] = max(abs(Awin));
        refined(w).featureOmega = omegaWin(idx);
        refined(w).featureValue = maxVal;
        refined(w).atEdge = (idx == 1) || (idx == nW);
        fprintf('  argmax |%s| = %.4e at omega = %.4f\n', spec.label, maxVal, omegaWin(idx));
    else
        [minVal, idx] = min(abs(Awin));
        refined(w).featureOmega = omegaWin(idx);
        refined(w).featureValue = minVal;
        refined(w).atEdge = (idx == 1) || (idx == nW);
        if idx > 1 && idx < nW
            refined(w).phaseFlipDeg = abs(mod(angle(Awin(idx+1)) - angle(Awin(idx-1)) + pi, 2*pi) - pi) * 180/pi;
        end
        fprintf('  min |%s| = %.4e at omega = %.4f (phase change either side ~ %.1f deg)\n', ...
            spec.label, minVal, omegaWin(idx), refined(w).phaseFlipDeg);
    end
    if refined(w).atEdge
        fprintf('  WARNING: %s feature sits at the window EDGE -- window is too narrow, widen omegaLo/omegaHi and rerun.\n', spec.label);
    end
end

%% Extra check on the A_{0,0} candidate zero: re-solve at a bumped
% truncation order ON THE EXACT SAME OMEGA GRID as the base window (not
% a separately-chosen grid -- see header note on why that confounds grid
% quantisation with a genuine shift) and compare where the minimum
% falls.
if CHECK_TRUNCATION_CONVERGENCE
    a00Idx = find(strcmp({featureSpecs.label}, 'A_{0,0}') & strcmp({featureSpecs.kind}, 'zero'), 1);
    if ~isempty(a00Idx)
        spec = featureSpecs(a00Idx);
        Mbig = M + 20; Pbig = P + 5; Nbig = N + 5;
        omegaCheck = refined(a00Idx).omega;   % SAME grid as the base result, on purpose
        nC = numel(omegaCheck);
        Acheck = complex(zeros(1, nC));

        fprintf('\n=== Truncation-convergence check on %s candidate zero ===\n', spec.label);
        fprintf('  base:   M=%d,P=%d,N=%d\n', M, P, N);
        fprintf('  bumped: M=%d,P=%d,N=%d\n', Mbig, Pbig, Nbig);
        fprintf('  (same omega grid as the base window, %d points, so any difference is\n', nC);
        fprintf('  physical, not grid quantisation)\n');
        ticC = tic;
        for k = 1:nC
            alpha = H*omegaCheck(k)^2/g;
            data = precomputeDeflectionData(alpha, beta0, gamma0, R0, nu, Mbig, Pbig, Nbig);
            wFun = @(r,theta) evaluateDeflection(data, r, theta);
            Ak = projectEMMOntoPlateModes(wFun, modeList(spec.row,:), R0, nu, 'NTheta', NTHETA);
            Acheck(k) = Ak;
            fprintf('  [%2d/%2d] omega=%.3f done (%.1fs elapsed)\n', k, nC, omegaCheck(k), toc(ticC));
        end
        [~, idxC] = min(abs(Acheck));
        omegaZeroBig = omegaCheck(idxC);
        omegaZeroBase = refined(a00Idx).featureOmega;
        shiftFrac = abs(omegaZeroBig - omegaZeroBase) / spec.dOmega;

        fprintf('  zero-crossing, base truncation:   omega = %.4f\n', omegaZeroBase);
        fprintf('  zero-crossing, bumped truncation: omega = %.4f\n', omegaZeroBig);
        if shiftFrac > 1
            fprintf('  WARNING: shifted by more than one grid step under truncation-convergence --\n');
            fprintf('  treat as a likely truncation-order artifact, not a physical feature, until (M,P,N) is validated further out.\n');
        else
            fprintf('  Stable under truncation-convergence (shift <= 1 grid step, on the SAME grid --\n');
            fprintf('  not a quantisation artifact). Looks genuine.\n');
        end
    end
end

%% One figure per mode: refined curve, coarse 35-point samples overlaid
% (where they fall inside the window), feature location marked. Kept as
% separate figures with their own vertical scale -- exactly what the
% combined j=0 plot from the coarse sweep could not show for A3,0/A4,0.
for w = 1:nWindows
    figure('Name', sprintf('%s refinement', refined(w).label));
    plot(refined(w).omega, abs(refined(w).A), 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
    hold on;
    % Build the legend labels dynamically -- do NOT assume the coarse
    % overlay is present: a window can legitimately contain none of the
    % original 35 sweep points, and a fixed 3-label legend() call then
    % silently mislabels whichever lines DO exist (caught by testing
    % this script, not obvious from reading it).
    legendEntries = {'refined (fine grid)'};
    inWin = S.omegaGrid >= refined(w).omega(1) & S.omegaGrid <= refined(w).omega(end);
    if any(inWin)
        plot(S.omegaGrid(inWin), abs(S.Amodal(refined(w).row, inWin)), 'kx', 'MarkerSize', 8, 'LineWidth', 1.2);
        legendEntries{end+1} = 'coarse 35-pt sweep';
    end
    yl = ylim;
    plot([refined(w).featureOmega, refined(w).featureOmega], yl, 'r--', 'LineWidth', 1.0);
    legendEntries{end+1} = 'feature location';
    ylim(yl);
    xlabel('\omega (rad/s)');
    ylabel(sprintf('|%s|', refined(w).label));
    title(sprintf('%s refinement: %s(%s) at \\omega = %.4f', refined(w).kind, refined(w).kind, refined(w).label, refined(w).featureOmega));
    legend(legendEntries, 'Location', 'best');
    grid on;
end

%% Save: feature locations + the full locally-sampled complex curves, so
% modalFeatureVector.m can use EITHER representation without rerunning
% this refinement.
resultsFile = fullfile(thisDir, 'modalFeatureRefinementResults.mat');
save(resultsFile, 'refined', 'featureSpecs', 'beta0', 'gamma0', 'R0', 'nu', 'M', 'P', 'N', 'H', 'g', 'modeList');
fprintf('\nRefinement complete. Saved to:\n  %s\n', resultsFile);

fprintf('\nFeature summary:\n');
for w = 1:nWindows
    if strcmp(refined(w).kind, 'max')
        kindLabel = 'omega_max';
    else
        kindLabel = 'omega_zero';
    end
    edgeFlag = '';
    if refined(w).atEdge
        edgeFlag = '  [AT WINDOW EDGE -- widen window]';
    end
    fprintf('  %-10s  %s = %.4f   |A| there = %.4e%s\n', refined(w).label, kindLabel, ...
        refined(w).featureOmega, refined(w).featureValue, edgeFlag);
end