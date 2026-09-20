%% Validation of the EMM-field-to-dry-mode projection (projectEMMOntoPlateModes.m)
% Step 2/4 of the dry-mode inverse-mapping plan (circularPlateMode.m ->
% circularPlateModeTester.m -> projectEMMOntoPlateModes.m -> THIS TESTER
% -> modalRAOReconnaissance.m).
%
% IMPORTANT SCOPE NOTE (read before trusting this too far): both tests
% below use SYNTHETIC fields, not real EMM output from Forward
% Model/deflection.m. That is deliberate, per the earlier review: the
% first target for this stage is "verify that the modal coefficients
% reconstruct [a] spatial field, and that increasing the number of
% retained dry modes reduces the reconstruction error" -- a property of
% the PROJECTION MATH, which doesn't require a real hydroelastic field to
% check. Wiring in deflection.m (once its own dependency chain is
% validated) is future work, not this file.
%
% Checks: (1) exact recovery -- build a field as a KNOWN finite
% combination of dry modes with known complex coefficients, confirm the
% projector recovers those exact coefficients, confirm projecting onto
% modes NOT in that combination gives ~0 (no cross-leakage), and confirm
% reconstructing from the recovered coefficients reproduces the field
% pointwise, off the projector's own internal quadrature grid; (2)
% reconstruction error is monotonically non-increasing as a NESTED
% sequence of larger mode sets is used on a field that is NOT built from
% the basis (an off-centre complex Gaussian bump) -- a fundamental
% property of projecting onto an orthonormal basis (a larger subspace
% can only get you a better or equal least-squares fit), checked here
% against an INDEPENDENT quadrature method (integral2) rather than the
% trapezoidal-in-theta method used inside the function under test, so a
% shared quadrature bug could not silently pass both sides.
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc %#ok<CLALL>

R = 2.3; nu = 0.3;

%% Test 1: exact recovery of a known finite modal combination
% A deliberately mixed set: one rigid-body mode (n=0,j=0), a rigid tilt
% (n=1,j=0), two flexible modes, and a higher (n,j) pair -- complex
% coefficients throughout, since a real hydroelastic field is complex
% (phase varies across the disc).
fprintf('=== Test 1: exact recovery of a known modal combination ===\n');
trueModes = [0 0; 1 0; 2 0; 0 1; 3 1];
trueA = [1.5+0.3i; -0.8i; 2.1; 0.6-1.2i; -0.4+0.4i];

nTrue = size(trueModes, 1);
trueModeStructs = cell(nTrue, 1);
for k = 1:nTrue
    trueModeStructs{k} = circularPlateMode(trueModes(k,1), trueModes(k,2), R, nu);
end
wFun = @(r, theta) syntheticField(r, theta, trueModeStructs, trueModes, trueA);

RECOVER_TOL = 1e-8; % forgiving vs. the ~1e-16 actually observed
Arecovered = projectEMMOntoPlateModes(wFun, trueModes, R, nu);
recoverRelErr = norm(Arecovered - trueA) / norm(trueA);
fprintf('  true A vs. recovered A:\n');
for k = 1:nTrue
    fprintf('    n=%d j=%d  true=%8.4f%+8.4fi  recovered=%8.4f%+8.4fi\n', ...
        trueModes(k,1), trueModes(k,2), real(trueA(k)), imag(trueA(k)), real(Arecovered(k)), imag(Arecovered(k)));
end
fprintf('  ||recovered - true|| / ||true|| = %.2e\n', recoverRelErr);
assert(recoverRelErr < RECOVER_TOL, 'Test1:recovery', ...
    'recovered amplitudes differ from true amplitudes by relative %.2e, exceeds tolerance %.2e', recoverRelErr, RECOVER_TOL);

% Test 1b: modes NOT present in the combination must project to ~0 --
% confirms no cross-leakage between modes (i.e. orthogonality is
% actually being exploited correctly here, not just assumed).
fprintf('\n=== Test 1b: decoy modes (absent from the field) project to ~0 ===\n');
DECOY_TOL = 1e-8;
decoyModes = [4 0; 2 1; 1 2];
Adecoy = projectEMMOntoPlateModes(wFun, decoyModes, R, nu);
for k = 1:size(decoyModes,1)
    fprintf('  n=%d j=%d  |A|=%.2e\n', decoyModes(k,1), decoyModes(k,2), abs(Adecoy(k)));
end
maxDecoy = max(abs(Adecoy));
fprintf('  max |decoy amplitude| = %.2e\n', maxDecoy);
assert(maxDecoy < DECOY_TOL, 'Test1b:decoyLeakage', ...
    'a mode absent from the field projected to %.2e, exceeds tolerance %.2e -- possible cross-leakage', maxDecoy, DECOY_TOL);

% Test 1c: reconstructing from the recovered coefficients reproduces the
% field POINTWISE, at test points deliberately off the projector's own
% internal 720-point theta grid and off any obviously special radius --
% checks the reconstruction as an actual spatial field, not just that
% the coefficient numbers happened to match.
fprintf('\n=== Test 1c: pointwise reconstruction matches the true field ===\n');
RECON_TOL = 1e-8;
rTest  = [0.30 1.10 1.90 2.20 0.01];
thTest = [0.37 1.90 3.30 5.10 0.05]; % not multiples of 2*pi/720
wTrue = syntheticField(rTest, thTest, trueModeStructs, trueModes, trueA);
wRecon = zeros(size(rTest));
for k = 1:nTrue
    wRecon = wRecon + Arecovered(k) * trueModeStructs{k}.radial(rTest) .* cos(trueModes(k,1)*thTest);
end
reconMaxErr = max(abs(wTrue - wRecon));
fprintf('  max pointwise |true - reconstructed| = %.2e\n', reconMaxErr);
assert(reconMaxErr < RECON_TOL, 'Test1c:pointwiseReconstruction', ...
    'pointwise reconstruction error %.2e exceeds tolerance %.2e', reconMaxErr, RECON_TOL);

%% Test 2: reconstruction error decreases monotonically as more modes are retained
% Target field is an off-centre, complex-amplitude Gaussian bump -- NOT
% constructed from circularPlateMode at all, so this genuinely tests
% approximation quality, not just self-consistency of Test 1's
% construction. An off-centre bump is a deliberately hard case: it has
% angular content spread across every n (breaking the disc's rotational
% symmetry), so this is not expected to converge to zero error quickly
% with only low-order modes -- what matters, and what is actually
% asserted, is that the error never goes UP as the retained mode set
% grows (a basic property of least-squares projection onto a nested
% sequence of subspaces of an orthonormal basis), not that it reaches
% some particular small value.
fprintf('\n=== Test 2: reconstruction error vs. number of retained modes ===\n');
x0 = 0.3*R; y0 = 0.15*R; sigma = 0.6*R;
bumpAmp = 1.2 - 0.5i;
wFunBump = @(r, theta) bumpAmp .* exp( -(((r.*cos(theta))-x0).^2 + ((r.*sin(theta))-y0).^2) ./ (2*sigma^2) );

Lmax = 5;
modeCache = cell(Lmax+1, Lmax+1);
for n = 0:Lmax
    for j = 0:Lmax
        modeCache{n+1, j+1} = circularPlateMode(n, j, R, nu);
    end
end

QUAD_OPTS = {'AbsTol', 1e-12, 'RelTol', 1e-10};
totalEnergy = integral2(@(r,theta) abs(wFunBump(r,theta)).^2 .* r, 0, R, 0, 2*pi, QUAD_OPTS{:});
fprintf('  total field energy (integral2, independent of the projector''s own quadrature) = %.6f\n', totalEnergy);

relErrs = zeros(Lmax+1, 1);
MONOTONIC_TOL = 1e-6; % forgiving vs. the essentially-exact non-increase expected
for L = 0:Lmax
    modeList = zeros((L+1)^2, 2);
    idx = 0;
    for n = 0:L
        for j = 0:L
            idx = idx + 1;
            modeList(idx,:) = [n j];
        end
    end
    Al = projectEMMOntoPlateModes(wFunBump, modeList, R, nu);
    reconFun = @(r, theta) reconstructFromCache(r, theta, modeCache, modeList, Al);
    errEnergy = integral2(@(r,theta) abs(wFunBump(r,theta) - reconFun(r,theta)).^2 .* r, 0, R, 0, 2*pi, QUAD_OPTS{:});
    relErrs(L+1) = sqrt(errEnergy / totalEnergy);
    fprintf('  L=%d  (%d modes, n=0..%d, j=0..%d)  relative L2 reconstruction error = %.6e\n', ...
        L, size(modeList,1), L, L, relErrs(L+1));
end

fprintf('  checking monotonic non-increase as L grows:\n');
for L = 1:Lmax
    ok = relErrs(L+1) <= relErrs(L) + MONOTONIC_TOL;
    fprintf('    err(L=%d)=%.3e <= err(L=%d)=%.3e + %.0e ?  %d\n', L, relErrs(L+1), L-1, relErrs(L), MONOTONIC_TOL, ok);
    assert(ok, 'Test2:monotonicError', ...
        'reconstruction error INCREASED going from L=%d (%.3e) to L=%d (%.3e) -- retaining more modes should never make the least-squares fit worse', ...
        L-1, relErrs(L), L, relErrs(L+1));
end
assert(relErrs(end) < relErrs(1), 'Test2:overallImprovement', ...
    'error at Lmax (%.3e) is not meaningfully smaller than at L=0 (%.3e) -- projection does not appear to be resolving any structure', ...
    relErrs(end), relErrs(1));

fprintf('\nAll tests passed.\n');

%% Local functions, placed at the end per MATLAB compatibility (matches
% the convention used elsewhere in this repo, e.g. nominalFRFReconnaissance.m)

function w = syntheticField(r, theta, modeStructs, modes, A)
    w = zeros(size(r));
    for k = 1:numel(modeStructs)
        w = w + A(k) * modeStructs{k}.radial(r) .* cos(modes(k,1)*theta);
    end
end

function w = reconstructFromCache(r, theta, modeCache, modeList, A)
    w = zeros(size(r));
    for k = 1:size(modeList,1)
        n = modeList(k,1); j = modeList(k,2);
        m = modeCache{n+1, j+1};
        w = w + A(k) * m.radial(r) .* cos(n*theta);
    end
end