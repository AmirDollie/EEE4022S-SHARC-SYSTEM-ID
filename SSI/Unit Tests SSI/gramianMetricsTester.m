%% gramianMetricsTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

%% Test 1: hand-computable diagonal case (eigenvalues are just the diagonal)
fprintf('=== Test 1: diagonal Wo, all metrics hand-checkable ===\n');
Wo1 = diag([4, 9, 25]);
m1 = gramianMetrics(Wo1);

assert(abs(m1.trace - 38) < 1e-12, 'Test 1 failed: trace should be 4+9+25=38.');
assert(abs(m1.lambdaMin - 4) < 1e-12, 'Test 1 failed: lambdaMin should be 4.');
assert(abs(m1.kappa - 25/4) < 1e-12, 'Test 1 failed: kappa should be 25/4=6.25.');
assert(abs(m1.logDet - log(4*9*25)) < 1e-10, 'Test 1 failed: logDet should be log(900).');
assert(isequal(m1.eigenvalues, [4;9;25]), 'Test 1 failed: eigenvalues should be sorted [4;9;25].');
fprintf('PASS: trace=%.4f, lambdaMin=%.4f, kappa=%.4f, logDet=%.4f\n', m1.trace, m1.lambdaMin, m1.kappa, m1.logDet);

%% Test 2: rank-deficient Wo (one EXACTLY zero eigenvalue) -> logDet=-Inf, kappa=Inf
fprintf('\n=== Test 2: rank-deficient Wo -> -Inf / Inf, not NaN ===\n');
Wo2 = diag([0, 5, 10]);
m2 = gramianMetrics(Wo2);

assert(m2.logDet == -Inf, 'Test 2 failed: logDet should be exactly -Inf.');
assert(m2.kappa == Inf, 'Test 2 failed: kappa should be exactly Inf.');
assert(m2.lambdaMin == 0, 'Test 2 failed: lambdaMin should be exactly 0.');
assert(abs(m2.trace - 15) < 1e-12, 'Test 2 failed: trace should still be finite (15).');
fprintf('PASS: logDet=%.4f, kappa=%.4f, lambdaMin=%.4f, trace=%.4f\n', m2.logDet, m2.kappa, m2.lambdaMin, m2.trace);

%% Test 3: tiny NEGATIVE eigenvalue from roundoff is clipped to 0, not propagated
fprintf('\n=== Test 3: tiny negative eigenvalue (roundoff) is clipped, not propagated ===\n');
Wo3 = diag([-1e-15, 5, 10]);
m3 = gramianMetrics(Wo3);
assert(m3.lambdaMin == 0, 'Test 3 failed: tiny negative eigenvalue should be clipped to exactly 0.');
assert(m3.logDet == -Inf, 'Test 3 failed: clipping to 0 should still correctly give -Inf for logDet.');
fprintf('PASS: tiny negative eigenvalue (-1e-15) correctly clipped to 0\n');

%% Test 4: non-diagonal case, cross-checked against MATLAB's own det()/eig() directly
fprintf('\n=== Test 4: non-diagonal case, cross-checked against det()/eig() directly ===\n');
M = [2 1 0; 1 3 1; 0 1 4];
Wo4 = M'*M;
m4 = gramianMetrics(Wo4);

directLogDet = log(det(Wo4));
directTrace = trace(Wo4);
directEig = sort(eig(Wo4), 'ascend');

assert(abs(m4.logDet - directLogDet) < 1e-8, 'Test 4 failed: logDet mismatch vs direct log(det(Wo)).');
assert(abs(m4.trace - directTrace) < 1e-10, 'Test 4 failed: trace mismatch vs direct trace(Wo).');
assert(max(abs(m4.eigenvalues - directEig)) < 1e-10, 'Test 4 failed: eigenvalues mismatch vs direct eig(Wo).');
fprintf('PASS: logDet, trace, and eigenvalues all match MATLAB''s direct calls\n');

%% Test 5: apply to a REAL Gramian from this project's own pipeline
fprintf('\n=== Test 5: applied to a real finiteObservabilityGramian output ===\n');
omega1 = 2.70; omega2 = 6.00; dt = 0.02;
A = buildKnownOscillatorA(omega1, omega2, dt);
C = [buildCandidateCrows(1+0.5i, 0.3-0.2i); buildCandidateCrows(-0.2+0.9i, 0.7-0.1i)];
Wo5 = finiteObservabilityGramian(A, C, 500);
m5 = gramianMetrics(Wo5);
fprintf('  trace=%.4f, lambdaMin=%.6f, kappa=%.4f, logDet=%.4f\n', m5.trace, m5.lambdaMin, m5.kappa, m5.logDet);
assert(all(isfinite([m5.trace, m5.kappa, m5.logDet])), 'Test 5 failed: expected finite metrics for this well-observed 2-sensor case.');
fprintf('PASS: all metrics finite and sensible for a real 2-sensor Gramian\n');

%% Test 6 (NEW): tiny POSITIVE eigenvalue is numerically rank-deficient,
% not silently treated as a valid, if small, observable direction
fprintf('\n=== Test 6: tiny positive eigenvalue treated as numerically rank-deficient ===\n');
WoTiny = diag([1e-18, 5, 10]);
mTiny = gramianMetrics(WoTiny);
assert(mTiny.lambdaMin == 0, 'Test 6 failed: 1e-18 should be treated as numerically zero.');
assert(mTiny.logDet == -Inf, 'Test 6 failed: logDet should be -Inf for this case.');
assert(mTiny.kappa == Inf, 'Test 6 failed: kappa should be Inf for this case.');
fprintf('PASS: eigenvalue=1e-18 correctly treated as numerically rank-deficient (not finite-but-tiny)\n');

%% Test 7 (NEW): a genuinely indefinite input must error, not be silently repaired
fprintf('\n=== Test 7: genuinely indefinite input must NOT be silently repaired ===\n');
WoBad = diag([-1, 5, 10]);
didError = false;
try
    gramianMetrics(WoBad);
catch ME
    didError = strcmp(ME.identifier, 'gramianMetrics:notPSD');
end
assert(didError, 'Test 7 failed: a genuinely indefinite matrix (eigenvalue=-1) should raise gramianMetrics:notPSD.');
fprintf('PASS: genuinely indefinite input (eigenvalue=-1) correctly raised an error instead of being silently clipped\n');

fprintf('\nAll 7 tests passed.\n');