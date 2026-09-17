%% finiteObservabilityGramianTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

%% Test 1: hand-computable case, N=1 (trivial: Wo = C'*C exactly)
fprintf('=== Test 1: N=1 reduces to C''*C exactly ===\n');
A_test = [0 1; -1 0];   % arbitrary, doesn't matter for N=1 since A^0=I
C_test = [1 2];
Wo1 = finiteObservabilityGramian(A_test, C_test, 1);
expected1 = C_test' * C_test;
assert(max(abs(Wo1(:) - expected1(:))) < 1e-14, 'Test 1 failed: N=1 case should equal C''*C exactly.');
fprintf('PASS: Wo(N=1) matches C''*C exactly.\n');

%% Test 2: hand-computable case, N=2, checked by explicit hand expansion
fprintf('\n=== Test 2: N=2, checked against explicit hand expansion ===\n');
Wo2 = finiteObservabilityGramian(A_test, C_test, 2);
expected2 = (C_test)'*(C_test) + (C_test*A_test)'*(C_test*A_test);
assert(max(abs(Wo2(:) - expected2(:))) < 1e-12, 'Test 2 failed: N=2 case does not match hand expansion.');
fprintf('PASS: Wo(N=2) matches explicit hand expansion C''C + (CA)''(CA).\n');

%% Test 3: symmetry and positive semi-definiteness (must ALWAYS hold, any A, any N)
fprintf('\n=== Test 3: Wo is always symmetric and positive semi-definite ===\n');
Wo3 = finiteObservabilityGramian(A_test, C_test, 10);
asymmetry = max(abs(Wo3(:) - reshape(Wo3', [], 1)));
eigWo3 = eig(Wo3);
assert(asymmetry < 1e-12, 'Test 3 failed: Wo should be exactly symmetric.');
assert(all(eigWo3 > -1e-10), 'Test 3 failed: Wo should be positive semi-definite (all eigenvalues >= 0).');
fprintf('PASS: Wo is symmetric (asymmetry=%.3e) and PSD (min eigenvalue=%.3e)\n', asymmetry, min(eigWo3));

%% Test 4: THE KEY CHECK -> for the undamped oscillator A, Wo genuinely
% does NOT converge as N grows (it keeps growing, roughly linearly with
% N, since each cycle adds equal observable energy). This is a
% deliberate demonstration, not a bug: it confirms directly why raw Wo
% magnitude must never be compared across different N, only rankings
% derived from it at a FIXED, consistently-chosen N.
fprintf('\n=== Test 4: Wo genuinely does NOT converge for an undamped system ===\n');
omega1 = 2.70; omega2 = 6.00; dt = 0.02;
A_undamped = buildKnownOscillatorA(omega1, omega2, dt);
C_undamped = buildCandidateCrows(1+0.5i, 0.3-0.2i);   % arbitrary single sensor

Ns = [10, 50, 200, 1000, 5000];
traces = zeros(size(Ns));
for k = 1:length(Ns)
    Wo_k = finiteObservabilityGramian(A_undamped, C_undamped, Ns(k));
    traces(k) = trace(Wo_k);
end

fprintf('  N       trace(Wo)    trace(Wo)/N\n');
for k = 1:length(Ns)
    fprintf('  %5d   %10.4f   %10.6f\n', Ns(k), traces(k), traces(k)/Ns(k));
end

% If Wo were converging, trace(Wo)/N would shrink toward 0 as N grows.
% Instead it should stay roughly CONSTANT, confirming trace(Wo) grows
% linearly with N rather than settling to a fixed value.
ratioSpread = max(traces./Ns) - min(traces./Ns);
assert(ratioSpread < 0.01, ...
    'Test 4: expected trace(Wo)/N to stay roughly constant (linear growth), got spread=%.4f -- check whether this still reflects the expected undamped behaviour.');
assert(traces(end) > 10*traces(1), ...
    'Test 4 failed: trace(Wo) should grow substantially with N for an undamped system, not plateau.');
fprintf('\nPASS: trace(Wo)/N stays roughly constant across N (spread=%.4f) -- confirms trace(Wo)\n', ratioSpread);
fprintf('grows linearly with N rather than converging, exactly as expected for an undamped\n');
fprintf('system. Raw Wo magnitude is therefore meaningless in isolation; only rankings at a\n');
fprintf('fixed, consistently-chosen N are meaningful (see the orchestration script''s N-sensitivity check).\n');

fprintf('\nAll 4 tests passed.\n');