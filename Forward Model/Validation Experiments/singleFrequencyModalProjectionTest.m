%% Single-frequency bridge: real EMM field -> dry modal coefficients -> reconstruction
% This is the "bridge" milestone between the validated dry-mode machinery
% (circularPlateMode.m + circularPlateModeTester.m, projectEMMOntoPlateModes.m
% + projectEMMOntoPlateModesTester.m, both validated on synthetic fields
% only) and the real wet EMM output from Forward Model/deflection.m. It
% does NOT sweep frequency -- that is modalRAOReconnaissance.m, later.
%
% Pipeline exercised here:
%   choose one alpha -> precomputeDeflectionData (solve EMM ONCE)
%     -> evaluateDeflection (cheap repeated evaluation, wraps into wFun)
%     -> projectEMMOntoPlateModes (real EMM field -> A_{n,j})
%     -> reconstruct w_recon from A_{n,j}
%     -> compare w_recon against the real EMM field on an INDEPENDENT
%        grid, for Montiel's 8 dominant modes AND for a nested sequence
%        of larger mode sets, checking the reconstruction error decreases
%        systematically as more modes are retained.
%
% *** NOTE ON PARAMETERS -- READ BEFORE TRUSTING THE NUMBERS BELOW ***
% This uses alpha=10, beta=1e-2, gamma=0.1, R=1, nu=0.3, M=10, P=6, N=10
% -- the EXACT parameter set already used and shown to converge in your
% own Forward Model/Unit-Tests/deflectionTester.m (its N-convergence
% check). This was a deliberate choice: I do NOT have your physical
% floe/water-depth values or the alpha=H*omega^2/g conversion constants
% (H, g, and which omega=2.7 rad/s case) in this session, and guessing
% them would risk silently validating against the wrong physics. Swap in
% your real (alpha, beta, gamma, R, nu, M, P, N) at the top of this
% script -- nothing else needs to change, since the rest of the pipeline
% is parameter-agnostic. Every number reported when this was last run
% (see comments inline below) is for the deflectionTester.m sanity
% parameters, NOT your physical floe.
%
% RUNTIME NOTE: the real EMM field is much more expensive to evaluate
% than the synthetic fields used in projectEMMOntoPlateModesTester.m.
% precomputeDeflectionData is the expensive part (Steps 0-6, run ONCE,
% ~1-2s here); each mode projection after that costs roughly 1.5-3s at
% NTheta=90 (see below), scaling with the number of modes requested. The
% full script (33 mode-projections: 8 dominant + 10 + 15 + 24 nested,
% with overlap) took on the order of a couple of minutes in Octave --
% expect MATLAB to be faster, but don't be surprised if this is not an
% instant script.
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', 'Animation')); % precomputeDeflectionData.m, evaluateDeflection.m live here, not in Forward Model/ itself
clear all, close all, clc %#ok<CLALL>

alpha = 10; beta = 1e-2; gamma = 0.1; M = 10; P = 6; N = 10;
R = 1; nu = 0.3;

%% Step 0: sanity cross-check -- precompute+evaluate must equal deflection.m
% deflection.m and precomputeDeflectionData.m/evaluateDeflection.m are
% two SEPARATE implementations of the same pipeline in this repo (the
% latter split into a one-time solve + cheap repeated evaluation, for
% exactly the performance reason noted above). Before trusting either as
% the "real EMM field", confirm they agree -- if they didn't, that would
% mean one of them has a bug, and nothing downstream would be
% trustworthy. This is NOT validating the EMM physics itself (that's a
% separate, much larger validation this session hasn't done -- see the
% scope note in projectEMMOntoPlateModesTester.m), only that these two
% code paths compute the same thing.
fprintf('=== Step 0: precompute+evaluate vs. direct deflection.m ===\n');
CROSSCHECK_TOL = 1e-10;
rCheck  = [0.1000, 0.5000, 0.9000, 0.9900, 0.0001];
thCheck = [0.2000, 1.0472, 2.5000, 4.2000, 5.9000];
etaDirect = zeros(size(rCheck));
for k = 1:numel(rCheck)
    etaDirect(k) = deflection(rCheck(k), thCheck(k), alpha, beta, gamma, R, nu, M, P, N);
end

fprintf('  solving EMM once via precomputeDeflectionData...\n');
data = precomputeDeflectionData(alpha, beta, gamma, R, nu, M, P, N);
etaFast = evaluateDeflection(data, rCheck, thCheck);

crossCheckMaxDiff = max(abs(etaDirect(:) - etaFast(:)));
fprintf('  max |deflection.m - evaluateDeflection(precomputed)| = %.3e (last run: exactly 0)\n', crossCheckMaxDiff);
assert(crossCheckMaxDiff < CROSSCHECK_TOL, 'Step0:crossCheck', ...
    'precomputeDeflectionData/evaluateDeflection disagree with deflection.m by %.3e -- do not trust results below until this is resolved', ...
    crossCheckMaxDiff);

% wFun: the cheap, elementwise-vectorized field handle projectEMMOntoPlateModes
% expects. evaluateDeflection is ALREADY vectorized this way (besseli(n,km*r)
% and exp(1i*n*theta) are both elementwise over matching-size r/theta
% arrays) -- no arrayfun wrapper needed, confirmed by the cross-check above
% using vector rCheck/thCheck directly.
wFun = @(r, theta) evaluateDeflection(data, r, theta);

% NTheta=90 instead of projectEMMOntoPlateModes' default 720: for this
% expensive real field, 720 is needless cost, not needless accuracy --
% confirmed separately that NTheta=90 vs NTheta=720 change a projected
% amplitude by ~1e-18 (the periodic trapezoidal rule is spectrally
% accurate regardless; see projectEMMOntoPlateModes.m's docstring). Use
% NTheta=90 throughout this script for that reason.
NTHETA = 90;

%% Step 1-2: project the real EMM field onto Montiel's 8 dominant modes
fprintf('\n=== Step 1-2: project onto Montiel''s 8 dominant modes ===\n');
% An,0 for n=0..4, An,1 for n=0..2 -- see Montiel.2013_2.pdf, Figure 1 /
% Sec. 3: "We consider a subset of eight modal amplitudes ... These are
% An,0, n=0-4, and An,1, n=0-2."
modeList = [0 0; 1 0; 2 0; 3 0; 4 0; 0 1; 1 1; 2 1];
A = projectEMMOntoPlateModes(wFun, modeList, R, nu, 'NTheta', NTHETA);

fprintf('  n  j       Re(A)         Im(A)        |A|        angle (rad)\n');
for k = 1:size(modeList,1)
    fprintf('  %d  %d   %10.6f   %10.6f   %10.6f   %8.4f\n', ...
        modeList(k,1), modeList(k,2), real(A(k)), imag(A(k)), abs(A(k)), angle(A(k)));
end
% Last run (alpha=10, deflectionTester.m parameters): A_{3,0} was the
% largest-magnitude flexible term (|A|~0.024), noticeably bigger than the
% rigid-body heave/tilt terms A_{0,0}, A_{1,0} (|A|~0.0025 each) -- i.e.
% at this alpha, flexural response already dominates rigid-body motion.
% Do not assume this generalises to other alpha -- that is exactly what
% modalRAOReconnaissance.m (frequency sweep, later) is for.

%% Step 3: reconstruct and check reconstruction error on an INDEPENDENT grid
% "Independent" = not the projector's own internal 90-point theta grid,
% and not any grid the projector's adaptive radial quadrature happened to
% pick -- a separate, fixed (r,theta) grid built fresh here.
fprintf('\n=== Step 3: reconstruction error (8-mode set, independent grid) ===\n');
modeCache = buildModeCache(0, 6, 0, 3, R, nu); % covers every (n,j) used below
rGrid  = linspace(1e-4, R-1e-4, 30);
thGrid = linspace(0, 2*pi*(59/60), 60); % 60 points over one period, no duplicate endpoint
[Rg, Tg] = ndgrid(rGrid, thGrid);
wTrue = wFun(Rg, Tg);
totalEnergy = sum(sum(abs(wTrue).^2 .* Rg));

wRecon8 = reconstructFromCache(Rg, Tg, modeCache, modeList, A);
relErr8 = sqrt(sum(sum(abs(wTrue - wRecon8).^2 .* Rg)) / totalEnergy);
fprintf('  relative L2 reconstruction error, 8 dominant modes = %.4f (last run: ~0.027)\n', relErr8);
fprintf('  (per the earlier review: this does NOT need to be tiny with only 8 modes --\n');
fprintf('   the real question is Step 4 below: does it improve as more modes are added?)\n');

%% Step 4: basis-convergence check -- does reconstruction improve with more modes?
% A NESTED sequence (each set is a strict superset of the previous one):
%   L0: n=0..4, j=0..1  (10 modes)
%   L1: n=0..4, j=0..2  (15 modes, adds j=2 at each n)
%   L2: n=0..5, j=0..3  (24 modes, adds n=5 and j=3)
% Nesting matters: for an orthonormal basis, enlarging the retained
% subspace can only decrease (or leave unchanged) the best achievable
% least-squares error against a FIXED field -- so relErr(L1) <= relErr(L0)
% and relErr(L2) <= relErr(L1) are not just hoped for, they are guaranteed
% by construction as long as (a) the sets really do nest and (b) the same
% fixed grid is used to measure error at every level (both true here).
fprintf('\n=== Step 4: nested basis-convergence check ===\n');
levelSpecs = { 0,4, 0,1;
               0,4, 0,2;
               0,5, 0,3 };
MONOTONIC_TOL = 1e-4; % forgiving vs. the exact guarantee above
relErrs = zeros(size(levelSpecs,1), 1);
for L = 1:size(levelSpecs,1)
    nLo = levelSpecs{L,1}; nHi = levelSpecs{L,2};
    jLo = levelSpecs{L,3}; jHi = levelSpecs{L,4};
    ml = [];
    for n = nLo:nHi
        for j = jLo:jHi
            ml = [ml; n j]; %#ok<AGROW>
        end
    end
    Al = projectEMMOntoPlateModes(wFun, ml, R, nu, 'NTheta', NTHETA);
    wReconL = reconstructFromCache(Rg, Tg, modeCache, ml, Al);
    relErrs(L) = sqrt(sum(sum(abs(wTrue - wReconL).^2 .* Rg)) / totalEnergy);
    fprintf('  n=%d..%d, j=%d..%d (%d modes): relative L2 error = %.4f\n', ...
        nLo, nHi, jLo, jHi, size(ml,1), relErrs(L));
end
% Last run: 0.0246 -> 0.0226 -> 0.0111 -- systematic decrease, consistent
% with the dry basis correctly representing more and more of the real
% EMM field as it grows, exactly the qualitative result the earlier
% review was looking for (not any particular target numbers).

fprintf('\n  checking monotonic non-increase:\n');
for L = 2:numel(relErrs)
    ok = relErrs(L) <= relErrs(L-1) + MONOTONIC_TOL;
    fprintf('    level %d (%.4f) <= level %d (%.4f) + %.0e ?  %d\n', L, relErrs(L), L-1, relErrs(L-1), MONOTONIC_TOL, ok);
    assert(ok, 'Step4:monotonicConvergence', ...
        'reconstruction error increased from level %d (%.4f) to level %d (%.4f) -- either a nesting mistake above, or a real problem with the projection/reconstruction', ...
        L-1, relErrs(L-1), L, relErrs(L));
end

fprintf('\nBridge test complete: the real EMM field is being represented in the dry basis,\n');
fprintf('and reconstruction improves systematically as more modes are retained.\n');
fprintf('Because the incident wave is symmetric about the propagation axis, this used the\n');
fprintf('cos(n*theta) family only (AngularKind default) -- consistent with Montiel_Thesis_6.pdf''s\n');
fprintf('statement that the disc response has no roll/antisymmetric component.\n');

%% Local functions, placed at the end per MATLAB compatibility (matches
% the convention used elsewhere in this repo, e.g. nominalFRFReconnaissance.m)

function cache = buildModeCache(nLo, nHi, jLo, jHi, R, nu)
    cache = struct();
    for n = nLo:nHi
        for j = jLo:jHi
            key = sprintf('n%d_j%d', n, j);
            cache.(key) = circularPlateMode(n, j, R, nu);
        end
    end
end

function w = reconstructFromCache(r, theta, cache, modeList, A)
    w = zeros(size(r));
    for k = 1:size(modeList,1)
        n = modeList(k,1); j = modeList(k,2);
        key = sprintf('n%d_j%d', n, j);
        m = cache.(key);
        w = w + A(k) * m.radial(r) .* cos(n*theta);
    end
end