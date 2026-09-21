%% modalRAOReconnaissanceTester.m
% Tests the SWEEP PLUMBING behind modalRAOReconnaissance.m -- not any
% physical hypothesis about where a resonance "should" be. Five checks,
% each with a hard assert, plus one non-fatal diagnostic at the end.
%
% COST NOTE -- deliberately mixed parameter regimes:
%   Test 1 (bridge regression) is pinned to a SPECIFIC previously-obtained
%   result and therefore MUST use the real nominal physical parameters
%   (M=50,P=10,N=10) -- there is no substitute. That is the slow part of
%   this file (one EMM solve at those truncation orders).
%
%   Tests 2-5 check generic numerical/structural properties of the
%   projection pipeline (quadrature convergence, cos/sin symmetry, basis-
%   set independence, determinism) that do not depend on which physical
%   floe parameters are used. They deliberately reuse the cheap sanity
%   parameters already validated elsewhere in this repo
%   (deflectionTester.m: alpha=10, beta=1e-2, gamma=0.1, R=1, M=10, P=6,
%   N=10) so this file runs in seconds, not minutes. This is a scope
%   choice: if you want these re-checked at the real nominal parameters
%   too, swap the parameter block in each test -- nothing else changes.
%
% This tester lives in Inverse Mapping/Unit-Tests/, so '..','..' reaches
% test_repo's top level. circularPlateMode.m and projectEMMOntoPlateModes.m
% live flat in Forward Model/ alongside deflection.m etc., so the first
% addpath line below already covers them -- no separate path needed.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));

% NOTE: clearvars, not clear all -- thisDir is needed again at the very
% end (the diagnostic section looks for modalRAOReconnaissanceResults.mat
% relative to it), so it must survive the clear.
clearvars -except thisDir
close all, clc

nu = 0.3;
modeList8 = [0 0; 1 0; 2 0; 3 0; 4 0; 0 1; 1 1; 2 1];
NTHETA = 90;

%% ==================================================================
%% Test 1: bridge regression against the previously-obtained result,
%% at alpha=10 with the REAL nominal physical parameters.
%% ==================================================================
fprintf('=== Test 1: bridge regression (nominal physical parameters, alpha=10) ===\n');
alphaNom = 10;
beta0 = 4.6985e-5; gamma0 = 1.4548e-3; R0 = 0.3830;
M0 = 50; P0 = 10; N0 = 10;

dataNom = precomputeDeflectionData(alphaNom, beta0, gamma0, R0, nu, M0, P0, N0);
wFunNom = @(r,theta) evaluateDeflection(dataNom, r, theta);
Anew = projectEMMOntoPlateModes(wFunNom, modeList8, R0, nu, 'NTheta', NTHETA);

% Reference values: recomputed independently with this exact pipeline
% (precomputeDeflectionData -> evaluateDeflection -> projectEMMOntoPlateModes,
% NTheta=90) at these exact parameters, and cross-checked against the two
% values quoted from the earlier single-alpha bridge run (A_{0,0},
% A_{1,0}) -- both matched to the precision quoted there.
Aref = [ -5.441518995200536e-03 + -2.531314313159397e-03i;   % n=0,j=0
          1.123214486612836e-02 + -2.811665238715328e-02i;   % n=1,j=0
          2.974695833021309e-02 +  7.469014640443048e-03i;   % n=2,j=0
         -3.311955803111737e-03 +  1.234308519051019e-02i;   % n=3,j=0
         -3.560405655852490e-03 + -3.984977105824656e-04i;   % n=4,j=0
         -1.660336716430219e-02 + -7.723641319033946e-03i;   % n=0,j=1
          2.464857292222929e-03 + -6.170107000521650e-03i;   % n=1,j=1
          1.139484613025274e-03 +  2.861074790494173e-04i ]; % n=2,j=1

fprintf('  n  j       Re(A)         Im(A)       Re(Aref)      Im(Aref)\n');
for k = 1:8
    fprintf('  %d  %d   %10.6f   %10.6f   %10.6f   %10.6f\n', ...
        modeList8(k,1), modeList8(k,2), real(Anew(k)), imag(Anew(k)), real(Aref(k)), imag(Aref(k)));
end

REGRESSION_TOL = 1e-6;
relErr = norm(Anew - Aref) / norm(Aref);
fprintf('  ||Anew - Aref|| / ||Aref|| = %.3e\n', relErr);
assert(relErr < REGRESSION_TOL, 'Test1:regression', ...
    'modal coefficients drifted from the reference bridge result by relative %.3e (tol %.3e) -- something in precomputeDeflectionData/evaluateDeflection/projectEMMOntoPlateModes changed', ...
    relErr, REGRESSION_TOL);
fprintf('  PASS\n');

%% ==================================================================
%% Shared cheap-parameter setup for Tests 2-5
%% ==================================================================
betaC = 1e-2; gammaC = 0.1; RC = 1; MC = 10; PC = 6; NC = 10;
alphaTest = [1, 7, 13];   % low/middle/high WITHIN the actual sweep's alpha
                           % range (0.5..13 in modalRAOReconnaissance.m) --
                           % cheap-parameter regime only, see COST NOTE above

dataC = cell(size(alphaTest));
wFunC = cell(size(alphaTest));
for a = 1:numel(alphaTest)
    dataC{a} = precomputeDeflectionData(alphaTest(a), betaC, gammaC, RC, nu, MC, PC, NC);
    wFunC{a} = @(r,theta) evaluateDeflection(dataC{a}, r, theta);
end

%% ==================================================================
%% Test 2: angular-quadrature convergence, NTheta = 45/90/180, at each
%% of the 3 representative alphas above.
%% ==================================================================
fprintf('\n=== Test 2: angular-quadrature convergence (NTheta=45/90/180) ===\n');
QUAD_TOL = 1e-6;
for a = 1:numel(alphaTest)
    A45  = projectEMMOntoPlateModes(wFunC{a}, modeList8, RC, nu, 'NTheta', 45);
    A90  = projectEMMOntoPlateModes(wFunC{a}, modeList8, RC, nu, 'NTheta', 90);
    A180 = projectEMMOntoPlateModes(wFunC{a}, modeList8, RC, nu, 'NTheta', 180);
    rel90  = norm(A90  - A180) / norm(A180);
    rel45  = norm(A45  - A180) / norm(A180);
    fprintf('  alpha=%5.2f:  ||A90-A180||/||A180||=%.3e   ||A45-A180||/||A180||=%.3e\n', ...
        alphaTest(a), rel90, rel45);
    assert(rel90 < QUAD_TOL, 'Test2:quadratureConvergence', ...
        'alpha=%.2f: NTheta=90 vs 180 relative difference %.3e exceeds tolerance %.3e', ...
        alphaTest(a), rel90, QUAD_TOL);
end
fprintf('  PASS (NTheta=90 justified at all 3 test alphas)\n');

%% ==================================================================
%% Test 3: cosine/sine symmetry, at the same 3 alphas, for the n>0
%% modes only (sin is undefined/meaningless at n=0).
%% ==================================================================
fprintf('\n=== Test 3: cosine/sine symmetry (n>0 modes) ===\n');
SYMMETRY_TOL = 1e-6;
modeListNpos = modeList8(modeList8(:,1) > 0, :);
for a = 1:numel(alphaTest)
    Acos = projectEMMOntoPlateModes(wFunC{a}, modeListNpos, RC, nu, 'NTheta', NTHETA, 'AngularKind', 'cos');
    Asin = projectEMMOntoPlateModes(wFunC{a}, modeListNpos, RC, nu, 'NTheta', NTHETA, 'AngularKind', 'sin');
    ratio = norm(Asin) / norm(Acos);
    fprintf('  alpha=%5.2f:  ||Asin||/||Acos|| = %.3e\n', alphaTest(a), ratio);
    assert(ratio < SYMMETRY_TOL, 'Test3:symmetry', ...
        'alpha=%.2f: sin/cos amplitude ratio %.3e exceeds tolerance %.3e -- assumed theta=0 symmetry may not hold', ...
        alphaTest(a), ratio, SYMMETRY_TOL);
end
fprintf('  PASS (response is confirmed cos-family-only, as assumed)\n');

%% ==================================================================
%% Test 4: coefficient independence from the requested mode set --
%% A_{2,0} must not depend on whether it is requested alone, as part of
%% the 8 Montiel modes, or as part of a larger 13-mode basis.
%% ==================================================================
fprintf('\n=== Test 4: coefficient independence from requested mode set ===\n');
INDEP_TOL = 1e-8;
aMid = 2;   % index into alphaTest/dataC/wFunC for alpha=10
modeListBig = [modeList8; 5 0; 3 1; 4 1; 0 2; 1 2];

Asingle = projectEMMOntoPlateModes(wFunC{aMid}, [2 0], RC, nu, 'NTheta', NTHETA);
A8      = projectEMMOntoPlateModes(wFunC{aMid}, modeList8, RC, nu, 'NTheta', NTHETA);
Abig    = projectEMMOntoPlateModes(wFunC{aMid}, modeListBig, RC, nu, 'NTheta', NTHETA);

A20_single = Asingle(1);
A20_eight  = A8(3);          % row 3 of modeList8 is [2 0]
A20_big    = Abig(3);        % same row index, since modeListBig extends modeList8

fprintf('  A_{2,0} alone      = %+.8f %+.8fi\n', real(A20_single), imag(A20_single));
fprintf('  A_{2,0} in 8-set   = %+.8f %+.8fi\n', real(A20_eight),  imag(A20_eight));
fprintf('  A_{2,0} in 13-set  = %+.8f %+.8fi\n', real(A20_big),    imag(A20_big));

relDiff81  = abs(A20_eight - A20_single) / abs(A20_single);
relDiffBig = abs(A20_big   - A20_single) / abs(A20_single);
fprintf('  relative diff (8-set vs alone)  = %.3e\n', relDiff81);
fprintf('  relative diff (13-set vs alone) = %.3e\n', relDiffBig);
assert(relDiff81 < INDEP_TOL && relDiffBig < INDEP_TOL, 'Test4:basisIndependence', ...
    'A_{2,0} changed with the requested mode set (%.3e, %.3e) beyond tolerance %.3e -- projection is coupling coefficients through basis size', ...
    relDiff81, relDiffBig, INDEP_TOL);
fprintf('  PASS\n');

%% ==================================================================
%% Test 5: sweep repeatability -- the SAME alpha solved twice, from
%% scratch (fresh precomputeDeflectionData both times, not a reused
%% `data` struct), must reproduce the same modal coefficients. This is
%% the determinism check: no hidden state, no randomness anywhere in
%% the pipeline.
%% ==================================================================
fprintf('\n=== Test 5: sweep repeatability (fresh solve, same alpha, twice) ===\n');
REPEAT_TOL = 1e-10;
alphaRepeat = 10;
dataRep1 = precomputeDeflectionData(alphaRepeat, betaC, gammaC, RC, nu, MC, PC, NC);
A_rep1 = projectEMMOntoPlateModes(@(r,theta) evaluateDeflection(dataRep1, r, theta), modeList8, RC, nu, 'NTheta', NTHETA);

dataRep2 = precomputeDeflectionData(alphaRepeat, betaC, gammaC, RC, nu, MC, PC, NC);
A_rep2 = projectEMMOntoPlateModes(@(r,theta) evaluateDeflection(dataRep2, r, theta), modeList8, RC, nu, 'NTheta', NTHETA);

repeatRelErr = norm(A_rep1 - A_rep2) / norm(A_rep1);
fprintf('  ||A_rep1 - A_rep2|| / ||A_rep1|| = %.3e\n', repeatRelErr);
assert(repeatRelErr < REPEAT_TOL, 'Test5:repeatability', ...
    'two fresh solves at the same alpha disagreed by relative %.3e (tol %.3e) -- look for hidden state', ...
    repeatRelErr, REPEAT_TOL);
fprintf('  PASS\n');

fprintf('\nAll 5 sweep-plumbing tests passed.\n');

%% ==================================================================
%% Diagnostic (NOT a hard assertion): frame-to-frame continuity of the
%% ACTUAL sweep output. Only runs if modalRAOReconnaissance.m has
%% already been run in this folder and left its results file behind --
%% this deliberately does not redo a 35-point nominal-parameter sweep
%% here just for a diagnostic plot. A genuine narrow response feature or
%% phase transition can legitimately cause a locally large jump; this is
%% a flag to go look, not a failure.
%% ==================================================================
fprintf('\n=== Diagnostic: frame-to-frame continuity of the real sweep (if available) ===\n');
resultsFile = fullfile(thisDir, '..', 'modalRAOReconnaissanceResults.mat');
if isfile(resultsFile)
    S = load(resultsFile);
    nF = size(S.Amodal, 2);
    jumpNorms = zeros(nF-1, 1);
    for k = 1:nF-1
        jumpNorms(k) = norm(S.Amodal(:,k+1) - S.Amodal(:,k));
    end
    figure;
    plot(S.omegaGrid(1:end-1), jumpNorms, 'o-', 'LineWidth', 1.0);
    xlabel('\omega (rad/s)'); ylabel('||A(k+1) - A(k)||');
    title('Diagnostic: frame-to-frame modal-vector jump across the real sweep');
    grid on;

    medJump = median(jumpNorms);
    flagged = find(jumpNorms > 10*medJump);
    if isempty(flagged)
        fprintf('  No jump exceeds 10x the median step (median=%.3e). Sweep looks continuous.\n', medJump);
    else
        fprintf('  %d step(s) exceed 10x the median jump (median=%.3e) -- inspect these before trusting\n', numel(flagged), medJump);
        fprintf('  that they are genuine features rather than dispersionRoots.m branch-switching:\n');
        for f = flagged'
            fprintf('    omega=%.3f -> %.3f  (jump=%.3e)\n', S.omegaGrid(f), S.omegaGrid(f+1), jumpNorms(f));
        end
    end
else
    fprintf('  modalRAOReconnaissanceResults.mat not found next to Inverse Mapping/ --\n');
    fprintf('  run modalRAOReconnaissance.m first, then rerun this tester for this diagnostic.\n');
end