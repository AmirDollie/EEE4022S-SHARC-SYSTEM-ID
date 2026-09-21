%% modalRAOReconnaissance.m
% Sweeps the nominal EMM forward model over a physical frequency grid and
% projects the resulting deflection field onto Montiel's 8 dominant dry
% circular-plate modes at each frequency, producing projected modal-
% response curves |A_{n,j}(omega)| and phase arg(A_{n,j}(omega)).
%
% NOT yet "Montiel RAOs": these are raw projection coefficients. Montiel's
% published RAOs are normalized by incident wave amplitude, with a
% separate normalization for pitch -- that normalization hasn't been
% applied here, so treat |A_{n,j}| as a reconnaissance signal (is there
% frequency-localised modal structure?), not yet as directly comparable
% to published RAO magnitudes.
%
% Deliberately kept simple: this generates trustworthy modal-response
% curves ONLY. It does NOT run findpeaks and does NOT attempt any
% inversion -- that is later work, once these curves exist to look at.
%
% Pipeline:
%   (beta0, gamma0, R0) -> {omega_k} -> {alpha_k}
%     -> EMM solve at each alpha_k   (precomputeDeflectionData/evaluateDeflection)
%     -> project onto the 8 dry modes (projectEMMOntoPlateModes)
%     -> A_{n,j}(omega_k)  ->  plots + saved complex data
%
% Depends on:
%   Forward Model             :  circularPlateMode.m, projectEMMOntoPlateModes.m
%   Forward Model/Animation   :  precomputeDeflectionData.m, evaluateDeflection.m
% (all four already live flat in Forward Model/ and Forward Model/Animation/
% in this repo, so the two addpath lines below cover everything -- no
% separate path needed for the dry-mode basis/projector.)
%
% *** RUNTIME NOTE *** each sweep point costs one full EMM solve
% (precomputeDeflectionData, at M=50,P=10,N=10) plus one 8-mode
% projection. The single-alpha bridge check that validated this pipeline
% (see modalRAOReconnaissanceTester.m, Test 1) took real time per point
% at these truncation orders -- with 35 points, do not expect this to
% finish in seconds. If it is uncomfortably slow on your machine, drop
% nFreq below to something like 8-10 for a first look at the curve shape
% before committing to the full sweep.
%
% See modalRAOReconnaissanceTester.m for validation of the SWEEP PLUMBING
% itself. This script assumes the EMM forward model and the dry
% basis/projector are already validated elsewhere (deflectionTester.m,
% circularPlateModeTester.m, projectEMMOntoPlateModesTester.m).

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', 'Forward Model', 'Animation'));

% NOTE: clearvars, not clear all -- thisDir is needed again below to
% build resultsFile, so it must survive the clear. `save` with a bare
% filename writes to MATLAB's current working directory, which is not
% necessarily this script's folder (e.g. if launched from the repo
% root) -- building an explicit path avoids the tester silently looking
% in the wrong place for this file.
clearvars -except thisDir
close all, clc

%% Nominal physical parameters (matches nominalFRFReconnaissance.m)
H = 1.88; g = 9.81; nu = 0.3;
beta0 = 4.6985e-5; gamma0 = 1.4548e-3; R0 = 0.3830;
M = 50; P = 10; N = 10;

%% Montiel's 8 dominant modes (Montiel.2013_2.pdf Sec. 3: An,0 n=0-4, An,1 n=0-2)
modeList = [0 0; 1 0; 2 0; 3 0; 4 0; 0 1; 1 1; 2 1];
nModes = size(modeList, 1);

%% Frequency grid, defined in physical omega (not alpha), since SSI will
% eventually hand back frequencies in physical units, not alpha. Bounds
% chosen to stay inside the truncation ceiling already validated for
% M=50,P=10,N=10 (alpha ~13.85, see nominalFRFReconnaissance.m header
% comment) with a margin, and starting away from alpha=0 for the same
% reason as that script's own coarse-sweep lower bound.
alphaMin = 0.5; alphaMax = 13;
omegaMin = sqrt(alphaMin*g/H);
omegaMax = sqrt(alphaMax*g/H);
nFreq = 35;
omegaGrid = linspace(omegaMin, omegaMax, nFreq);
alphaGrid = H * omegaGrid.^2 / g;

fprintf('Sweeping %d frequencies: omega = %.3f .. %.3f rad/s  (alpha = %.3f .. %.3f)\n', ...
    nFreq, omegaGrid(1), omegaGrid(end), alphaGrid(1), alphaGrid(end));

%% Main sweep: one EMM solve + one modal projection per frequency
NTHETA = 90;   % validated vs. NTheta=720 in the single-alpha bridge check

Amodal = complex(zeros(nModes, nFreq));
ticStart = tic;
for k = 1:nFreq
    alpha = alphaGrid(k);

    % expensive EMM solve, ONCE per frequency
    data = precomputeDeflectionData(alpha, beta0, gamma0, R0, nu, M, P, N);

    % cheap, elementwise-vectorized spatial evaluator
    wFun = @(r, theta) evaluateDeflection(data, r, theta);

    % modal decomposition onto the 8 dry modes
    Amodal(:, k) = projectEMMOntoPlateModes(wFun, modeList, R0, nu, 'NTheta', NTHETA);

    fprintf('  [%2d/%2d] omega=%.3f alpha=%.3f done (%.1fs elapsed)\n', ...
        k, nFreq, omegaGrid(k), alpha, toc(ticStart));
end

%% Figure 1: rigid/low-order family |A_{n,0}|, n=0..4 (modeList rows 1-5)
figure;
hold on;
legendLabels1 = cell(5,1);
for row = 1:5
    plot(omegaGrid, abs(Amodal(row,:)), 'LineWidth', 1.2);
    legendLabels1{row} = sprintf('|A_{%d,0}|', modeList(row,1));
end
xlabel('\omega (rad/s)');
ylabel('|A_{n,0}|');
title('Projected modal response: rigid/low-order family (j=0)');
legend(legendLabels1);
grid on;

%% Figure 2: j=1 family |A_{n,1}|, n=0..2 (modeList rows 6-8) -- kept as a
% SEPARATE figure rather than overlaid on Figure 1, so a large heave/tilt
% curve cannot visually hide a small flexural feature (the exact problem
% this reconnaissance is trying to escape).
figure;
hold on;
legendLabels2 = cell(3,1);
for row = 6:8
    plot(omegaGrid, abs(Amodal(row,:)), 'LineWidth', 1.2);
    legendLabels2{row-5} = sprintf('|A_{%d,1}|', modeList(row,1));
end
xlabel('\omega (rad/s)');
ylabel('|A_{n,1}|');
title('Projected modal response: j=1 family');
legend(legendLabels2);
grid on;

%% Figure 3: phase diagnostic, one subplot per mode. Phase behaviour
% around any amplitude feature may help distinguish a genuine
% resonance-like response from a random numerical bump.
figure('Name', 'Projected modal response phase (unwrapped)');
for row = 1:nModes
    subplot(4,2,row);
    plot(omegaGrid, unwrap(angle(Amodal(row,:)))*180/pi, 'LineWidth', 1.0);
    title(sprintf('arg(A_{%d,%d})', modeList(row,1), modeList(row,2)));
    xlabel('\omega (rad/s)'); ylabel('deg');
    grid on;
end

%% Save everything: later feature-selection work operates on this file
% without rerunning the EMM sweep. Saved next to THIS script (not to
% pwd), so modalRAOReconnaissanceTester.m's diagnostic section finds it
% regardless of where the sweep was launched from.
resultsFile = fullfile(thisDir, 'modalRAOReconnaissanceResults.mat');
save(resultsFile, ...
    'omegaGrid', 'alphaGrid', 'Amodal', 'modeList', ...
    'beta0', 'gamma0', 'R0', 'nu', 'M', 'P', 'N', 'H', 'g');

fprintf('\nSweep complete. Saved to:\n  %s\n', resultsFile);
fprintf('Look at the plots: is there clear frequency-localised structure in any\n');
fprintf('individual structural mode? (No peak detection run here, deliberately.)\n');