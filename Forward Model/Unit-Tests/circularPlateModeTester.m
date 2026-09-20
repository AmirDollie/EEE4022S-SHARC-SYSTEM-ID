%% Validation of the dry, free-edge circular-plate eigenmode basis
% (circularPlateMode.m). This is the basis Meylan & Squire (1996) use
% directly and Montiel (2013) projects the wet EMM response onto --
% Step 1/4 of the dry-mode inverse-mapping plan (circularPlateMode.m ->
% this tester -> projectEMMOntoPlateModes.m -> modalRAOReconnaissance.m).
%
% Checks: (1) both free-edge conditions (moment, shear) vanish at the
% roots found; (2) modal orthogonality across distinct j at fixed n,
% INCLUDING the rigid-body modes against their own flexible modes --
% this is the check that actually caught two real bugs during
% development, see the note before Test 2 -- AND that every mode is
% correctly normalised (diagonal of the same inner product = 1, not just
% the off-diagonal cross terms = 0); (3) rigid-body modes are identified
% separately from the flexible root search, with lambda=0, and the next
% index at that n is a genuine nonzero root, not another lambda=0;
% (4) eigenvalues are stable as root-search resolution (NBrackets)
% increases, now over the actual low-order basis this project uses
% (n=0..4, j=0,1), not just one spot-checked (n,j); (5) the numerical
% root ordering at each n produces exactly j interior nodal circles,
% confirming the index j really does mean "number of nodal circles" and
% not just "j-th root found"; (6) an EXTERNAL benchmark against natural
% frequencies published in Montiel's own PhD thesis (independent of this
% codebase, computed with different code) -- see the note before Test 6.
%
% All numeric checks now assert (tolerances deliberately forgiving, not
% machine-epsilon-tight -- see each test) so a regression fails loudly
% instead of only showing up as a number in the console someone has to
% notice.
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc %#ok<CLALL>

R = 2.3; nu = 0.3;
TOL = 1e-10; % forgiving vs. the ~1e-14 to 1e-16 actually observed -- tight
             % enough to catch a real regression, loose enough not to
             % nuisance-fail on machine/toolbox floating-point noise

%% Test 1: free-edge BC residuals vanish at every found root
% Independent recompute of the moment (BC1) and shear (BC2) conditions,
% NOT calling circularPlateMode's own internal edgeConditionValues --
% a residual check that reused the same (possibly buggy) formula the
% root-finder used would be checking the root-finder against itself.
fprintf('=== Test 1: free-edge residuals ===\n');
cases = [0 1; 0 2; 1 1; 1 2; 2 0; 2 1; 3 0; 4 0];
for k = 1:size(cases,1)
    n = cases(k,1); j = cases(k,2);
    m = circularPlateMode(n, j, R, nu);
    [bc1J, bc1I, bc2J, bc2I] = localEdgeConditionValues(m.lambda, n, nu);
    bc1 = bc1J + m.C*bc1I;
    bc2 = bc2J + m.C*bc2I;
    scale = max([abs(bc1J), abs(bc1I), abs(bc2J), abs(bc2I), 1]);
    relResidMoment = abs(bc1)/scale;
    relResidShear = abs(bc2)/scale;
    fprintf('  n=%d j=%d  lambda=%.6f  relresid_moment=%.2e  relresid_shear=%.2e  (both should be ~1e-14 or smaller)\n', ...
        n, j, m.lambda, relResidMoment, relResidShear);
    assert(relResidMoment < TOL, 'Test1:momentResidual', ...
        'n=%d j=%d: moment (BC1) residual %.2e exceeds tolerance %.2e', n, j, relResidMoment, TOL);
    assert(relResidShear < TOL, 'Test1:shearResidual', ...
        'n=%d j=%d: shear (BC2) residual %.2e exceeds tolerance %.2e', n, j, relResidShear, TOL);
end

%% Test 2: modal orthogonality AND normalisation, including rigid-body
% IMPORTANT: this is the check that matters most in this file. Two
% independent algebra bugs made it into circularPlateMode.m during
% development (a dropped R''/r term in the shear BC, then a wrong
% coefficient in the Bessel third-derivative formula) and BOTH produced
% a root-finder that was internally self-consistent -- Test 1 above
% would have passed at machine precision for both bugs -- but was
% solving the WRONG equation. Only this orthogonality check caught
% them: the first bug showed up as ~5e-3 cross-orthogonality at n=0,
% the second (rarer, n-dependent, invisible at n=0) as ~5e-3 between
% the n=1 rigid tilt mode and its own flexible modes specifically.
% After both fixes this comes out at machine precision (~1e-16)
% everywhere tested, including rigid-vs-flexible pairs. If this test
% ever regresses to anything above TOL, do not trust Test 1 passing
% as evidence the formulas are right -- it isn't.
%
% The diagonal (self) inner products are checked here too: circularPlateMode
% claims int_0^2pi int_0^R |w|^2 r dr dtheta = 1 for every mode, but that
% was previously only exercised implicitly (the off-diagonal cross terms
% would still vanish even if every mode's normalisation constant A were
% wrong by the same wrong factor -- orthogonality alone cannot catch a
% scaling error common to all modes at a given n).
fprintf('\n=== Test 2: modal orthogonality + normalisation (fixed n, distinct j) ===\n');
maxCrossIP = 0;
maxNormErr = 0;
for n = [0 1 2 3]
    m0 = circularPlateMode(n, 0, R, nu);
    m1 = circularPlateMode(n, 1, R, nu);
    m2 = circularPlateMode(n, 2, R, nu);

    ip01 = radialInnerProduct(m0, m1, R);
    ip02 = radialInnerProduct(m0, m2, R);
    ip12 = radialInnerProduct(m1, m2, R);
    maxCrossIP = max([maxCrossIP, abs(ip01), abs(ip02), abs(ip12)]);

    if n == 0, angularFactor = 2*pi; else, angularFactor = pi; end
    norm0 = angularFactor * radialInnerProduct(m0, m0, R);
    norm1 = angularFactor * radialInnerProduct(m1, m1, R);
    norm2 = angularFactor * radialInnerProduct(m2, m2, R);
    maxNormErr = max([maxNormErr, abs(norm0-1), abs(norm1-1), abs(norm2-1)]);

    fprintf('  n=%d: |<j0,j1>|=%.2e  |<j0,j2>|=%.2e  |<j1,j2>|=%.2e  norms=[%.12f %.12f %.12f]  (rigid-body involved: j0 for n=0,1)\n', ...
        n, abs(ip01), abs(ip02), abs(ip12), norm0, norm1, norm2);
end
fprintf('  max |cross inner product| across all pairs tested: %.2e  (should be ~1e-14 or smaller)\n', maxCrossIP);
fprintf('  max |normalisation - 1| across all modes tested:   %.2e  (should be ~1e-14 or smaller)\n', maxNormErr);
assert(maxCrossIP < TOL, 'Test2:orthogonality', ...
    'max cross inner product %.2e exceeds tolerance %.2e', maxCrossIP, TOL);
assert(maxNormErr < TOL, 'Test2:normalisation', ...
    'max normalisation error %.2e exceeds tolerance %.2e', maxNormErr, TOL);

%% Test 3: rigid-body modes separated correctly from the flexible search
fprintf('\n=== Test 3: rigid-body modes ===\n');
heave = circularPlateMode(0, 0, R, nu);
tilt  = circularPlateMode(1, 0, R, nu);
fprintf('  n=0,j=0 (heave): isRigidBody=%d, lambda=%.6f (expect 1, 0)\n', heave.isRigidBody, heave.lambda);
fprintf('  n=1,j=0 (tilt):  isRigidBody=%d, lambda=%.6f (expect 1, 0)\n', tilt.isRigidBody, tilt.lambda);
assert(heave.isRigidBody && heave.lambda == 0, 'Test3:heave', 'n=0,j=0 is not flagged/valued as rigid-body');
assert(tilt.isRigidBody && tilt.lambda == 0, 'Test3:tilt', 'n=1,j=0 is not flagged/valued as rigid-body');

heave1 = circularPlateMode(0, 1, R, nu);
tilt1  = circularPlateMode(1, 1, R, nu);
fprintf('  n=0,j=1 (next):  isRigidBody=%d, lambda=%.6f (expect 0, a genuine nonzero root, not another 0)\n', ...
    heave1.isRigidBody, heave1.lambda);
fprintf('  n=1,j=1 (next):  isRigidBody=%d, lambda=%.6f (expect 0, a genuine nonzero root, not another 0)\n', ...
    tilt1.isRigidBody, tilt1.lambda);
assert(~heave1.isRigidBody && heave1.lambda > 0, 'Test3:heave1', 'n=0,j=1 should be a genuine flexible root, not another rigid-body zero');
assert(~tilt1.isRigidBody && tilt1.lambda > 0, 'Test3:tilt1', 'n=1,j=1 should be a genuine flexible root, not another rigid-body zero');

%% Test 4: eigenvalue stability as root-search resolution increases
% Not a convergence-to-truth claim -- just confirms a coarse bracket
% scan does not mislocate a root relative to a much finer one, same
% caution as dispersionTester.m in Forward Model/Unit-Tests/. Expanded
% from the single (n,j)=(2,1) spot-check to the actual low-order basis
% (n=0..4, j=0,1) this project intends to use for the modal projection --
% the dry solve is cheap enough there is no reason to only check one pair.
fprintf('\n=== Test 4: eigenvalue stability vs root-search resolution (n=0..4, j=0,1) ===\n');
RES_TOL = 1e-8; % forgiving vs. the ~1e-15 to 1e-16 actually observed
maxResDiff = 0;
for n = 0:4
    for j = 0:1
        mCoarse = circularPlateMode(n, j, R, nu, 'NBrackets', 500);
        mFine   = circularPlateMode(n, j, R, nu, 'NBrackets', 8000);
        d = abs(mCoarse.lambda - mFine.lambda);
        maxResDiff = max(maxResDiff, d);
        fprintf('  n=%d j=%d: lambda(NBrackets=500)=%.10f, lambda(NBrackets=8000)=%.10f, diff=%.2e\n', ...
            n, j, mCoarse.lambda, mFine.lambda, d);
        assert(d < RES_TOL, 'Test4:resolutionStability', ...
            'n=%d j=%d: coarse/fine lambda differ by %.2e, exceeds tolerance %.2e', n, j, d, RES_TOL);
    end
end
fprintf('  max diff across all (n,j) tested: %.2e (should be ~0)\n', maxResDiff);

%% Test 5: nodal-circle count matches the index j
% Montiel (2013) and Meylan & Squire (1996) both interpret n and j as the
% number of nodal DIAMETERS and nodal CIRCLES of the mode shape,
% respectively -- see e.g. Montiel_Thesis_6.pdf, Sec. 6.1: "the angular
% and radial mode indexes (n and j, respectively) represent the number of
% nodal diameters and nodal circles of the mode shapes". This test counts
% actual sign changes of the radial part R_{n,j}(r) on the open interval
% (0,R) and checks it equals j exactly. This is the check the earlier
% review specifically asked for: Tests 1-4 confirm the eigenproblem was
% solved correctly, but say nothing about whether the ROOT ORDERING
% (i.e. which root gets labelled j=0, j=1, j=2, ...) matches the physical
% convention the rest of the project (and any RAO labelling, A_{n,j})
% depends on. A root-finder that silently returned roots out of order
% would pass every test above and still be wrong here.
fprintf('\n=== Test 5: nodal-circle count vs index j ===\n');
nodalCases = [0 0; 1 0; 0 1; 0 2; 1 1; 1 2; 2 0; 2 1; 3 0; 3 1; 4 0];
rGrid = linspace(1e-6, R - 1e-9, 20000); % avoid the exact endpoints (r=0
                                          % is a removable singularity of
                                          % J_n/I_n for n>0; r=R is the
                                          % free edge, not necessarily 0)
for k = 1:size(nodalCases,1)
    n = nodalCases(k,1); j = nodalCases(k,2);
    m = circularPlateMode(n, j, R, nu);
    vals = m.radial(rGrid);
    s = sign(vals);
    s(s == 0) = 1; % guard against landing exactly on a zero
    nodalCircles = sum(abs(diff(s)) == 2);
    fprintf('  n=%d j=%d: nodal circles counted = %d  (expect %d)\n', n, j, nodalCircles, j);
    assert(nodalCircles == j, 'Test5:nodalCircleCount', ...
        'n=%d j=%d: counted %d interior nodal circles, expected exactly j=%d', n, j, nodalCircles, j);
end

%% Test 6: external benchmark against Montiel's PhD thesis (independent source)
% Everything above is INTERNAL consistency (Tests 1,4,5: the root-finder
% agrees with itself and with the (n,j) labelling convention) or a
% property that even a wrong-but-self-adjoint eigenproblem could
% accidentally satisfy (Test 2: orthogonality). None of that proves the
% governing equations themselves are the physically correct free-edge
% plate equations -> only comparing against a number this codebase did
% NOT produce can do that.
%
% Montiel_Thesis_9_conclusion_and_appendix.pdf, Table F.2, gives measured
% natural frequencies f_{n,j} (Hz) of the free-edge circular PVC disks
% used in Montiel's wave-tank experiments, computed from this exact same
% dry free-edge eigenvalue problem with Montiel's own (different) code --
% they are NOT a hand solution, they come from solving det B_n(lambda)=0
% independently and converting lambda to a physical frequency via the
% disk's real material properties, which are given in
% Montiel_Thesis_4.pdf, Table 4.2.
% 
%  Using the h=3mm disk and the
% EXPERIMENTALLY MEASURED properties column of Table 4.2 (R=0.72m,
% E=838 MPa, rho=623 kg/m^3, nu=0.3 -- NOT the manufacturer
% "specifications" column, E=1300 MPa/rho=700 kg/m^3, which was also
% tried and gives a uniform ~17.5% offset across every mode, consistent
% with Montiel's own text explaining measured vs. spec'd material
% properties disagree substantially for this PVC):
%
%   f_{n,j} = lambda_{n,j}^2 / (2*pi) * sqrt(D / (rho*h*R^4)),   D = E*h^3/(12*(1-nu^2))
%
% Checked against all 18 free-edge (n,j) pairs tabulated (n=0..4, j up to
% 3) -- NOT the adjacent f^(eb) column in the same table, which is a
% different (edge-beam) boundary condition circularPlateMode.m does not
% implement. Max relative error found: 0.12% (n=2,j=0), most modes agree
% to <0.03% -- well within what 2-decimal-place rounding in the
% published table alone would explain.
fprintf('\n=== Test 6: external benchmark vs. Montiel thesis Table F.2 (h=3mm disk) ===\n');
R_disk = 0.72; nu_disk = 0.3;
h_disk = 0.003; E_disk = 838e6; rho_disk = 623; % Table 4.2, "experiments" column
D_disk = E_disk*h_disk^3 / (12*(1-nu_disk^2));
rhoh_disk = rho_disk*h_disk;
% columns: n, j, f_{n,j} (Hz), free-edge model, Table F.2
benchmarkRows = [
    2 0  1.73;  0 1  2.91;  3 0  4.02;  1 1  6.62;  4 0  7.06;  2 1 11.40;
    0 2 12.43;  3 1 17.13;  1 2 19.33;  4 1 23.77;  2 2 27.27;  0 3 28.37;
    3 2 36.19;  1 3 38.45;  4 2 46.04;  2 3 49.56;  3 3 61.64;  4 3 74.68 ];
BENCH_TOL_PCT = 1.0; % forgiving vs. the ~0.12% max actually observed --
                      % the table itself is only published to 2 decimal
                      % places, so sub-percent agreement is the practical
                      % ceiling, not a sign of a looser check
maxBenchRelErr = 0;
for k = 1:size(benchmarkRows,1)
    n = benchmarkRows(k,1); j = benchmarkRows(k,2); fTable = benchmarkRows(k,3);
    m = circularPlateMode(n, j, R_disk, nu_disk);
    omega = m.lambda^2 * sqrt(D_disk / (rhoh_disk*R_disk^4));
    fMine = omega / (2*pi);
    relErrPct = abs(fMine - fTable) / fTable * 100;
    maxBenchRelErr = max(maxBenchRelErr, relErrPct);
    fprintf('  n=%d j=%d: f_mine=%.4f Hz  f_table=%.4f Hz  relerr=%.2f%%\n', n, j, fMine, fTable, relErrPct);
    assert(relErrPct < BENCH_TOL_PCT, 'Test6:externalBenchmark', ...
        'n=%d j=%d: %.2f%% relative error vs. Montiel Table F.2 exceeds tolerance %.1f%%', n, j, relErrPct, BENCH_TOL_PCT);
end
fprintf('  max relative error across all 18 benchmarked modes: %.2f%%\n', maxBenchRelErr);

fprintf('\nAll tests passed.\n');

%% Local functions, placed at the end per MATLAB compatibility (matches
% the convention used elsewhere in this repo, e.g. nominalFRFReconnaissance.m)

function ip = radialInnerProduct(mA, mB, R)
    ip = integral(@(r) mA.radial(r).*mB.radial(r).*r, 0, R);
end

function [bc1J, bc1I, bc2J, bc2I] = localEdgeConditionValues(lambda, n, nu)
    [Jv, Jp, Jpp, Jppp] = localBesselJderivs(n, lambda);
    [Iv, Ip, Ipp, Ippp] = localBesselIderivs(n, lambda);
    bc1J = lambda^2*Jpp + nu*lambda*Jp - nu*n^2*Jv;
    bc1I = lambda^2*Ipp + nu*lambda*Ip - nu*n^2*Iv;
    bc2J = lambda^3*Jppp + lambda^2*Jpp - lambda*(1+n^2*(2-nu))*Jp + n^2*(3-nu)*Jv;
    bc2I = lambda^3*Ippp + lambda^2*Ipp - lambda*(1+n^2*(2-nu))*Ip + n^2*(3-nu)*Iv;
end

function [Jv, Jp, Jpp, Jppp] = localBesselJderivs(n, x)
    Jv  = besselj(n, x);
    Jp  = besselj(n-1, x) - (n/x)*Jv;
    Jpp = -Jp/x - (1 - n^2/x^2)*Jv;
    Jppp = Jp*((n^2+2)/x^2 - 1) + Jv*(1/x - 3*n^2/x^3);
end

function [Iv, Ip, Ipp, Ippp] = localBesselIderivs(n, x)
    Iv  = besseli(n, x);
    Ip  = besseli(n-1, x) - (n/x)*Iv;
    Ipp = -Ip/x + (1 + n^2/x^2)*Iv;
    Ippp = Ip*((n^2+2)/x^2 + 1) + Iv*(-1/x - 3*n^2/x^3);
end