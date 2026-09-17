%% buildKnownOscillatorATester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

omega1 = 2.70; omega2 = 6.00; dt = 0.02;
A = buildKnownOscillatorA(omega1, omega2, dt);

%% Test 1: shape and block-diagonal structure
fprintf('=== Test 1: shape and structure ===\n');
assert(isequal(size(A), [4 4]), 'Test 1 failed: A must be 4x4.');
assert(max(abs(A(1:2,3:4)), [], 'all') < 1e-14, 'Test 1 failed: off-diagonal block (1:2,3:4) should be exactly zero.');
assert(max(abs(A(3:4,1:2)), [], 'all') < 1e-14, 'Test 1 failed: off-diagonal block (3:4,1:2) should be exactly zero.');
assert(max(abs(imag(A(:)))) == 0, 'Test 1 failed: A must be purely real.');
fprintf('PASS: A is 4x4, real, block-diagonal.\n');

%% Test 2: eigenvalues are EXACTLY the two known frequencies, on the unit circle
fprintf('\n=== Test 2: eigenvalues match e^{+-i*omega*dt} exactly ===\n');
eigA = eig(A);
expected = [exp(1i*omega1*dt); exp(-1i*omega1*dt); exp(1i*omega2*dt); exp(-1i*omega2*dt)];

% match each computed eigenvalue to its closest expected value (order from
% eig() is not guaranteed to match the order we listed expected in)
matched = false(4,1);
maxErr = 0;
for k = 1:4
    [err, idx] = min(abs(eigA(k) - expected));
    maxErr = max(maxErr, err);
    matched(idx) = true;
end
assert(all(matched), 'Test 2 failed: not all four expected eigenvalues were matched.');
assert(maxErr < 1e-12, 'Test 2 failed: eigenvalue error too large (%.3e).', maxErr);
fprintf('PASS: all 4 eigenvalues match e^{+-i*omega*dt} to within %.3e\n', maxErr);

%% Test 3: |eigenvalue| = 1 exactly (undamped, zeta=0, matching Level 1)
fprintf('\n=== Test 3: eigenvalues lie exactly on the unit circle (zeta=0) ===\n');
magErr = max(abs(abs(eigA) - 1));
assert(magErr < 1e-14, 'Test 3 failed: eigenvalue magnitudes should be exactly 1 (undamped).');
fprintf('PASS: |eigenvalue|=1 to within %.3e (undamped, matching Level 1''s zeta~0 result)\n', magErr);

%% Test 4: the actual GENERATED SIGNAL matches cos(omega*t) directly,
% not just that the eigenvalues look right in isolation -- this is the
% test that would have caught the earlier C-mapping-style bug if it were
% present here instead.
fprintf('\n=== Test 4: simulated state trajectory matches cos/sin(omega*t) directly ===\n');
nSteps = 200;
x = zeros(4, nSteps);
x(:,1) = [1; 0; 1; 0];   % t=0: cos(0)=1, sin(0)=0, for BOTH frequencies
for k = 1:nSteps-1
    x(:,k+1) = A*x(:,k);
end

tVec = (0:nSteps-1)*dt;
expected_cos1 = cos(omega1*tVec);
expected_sin1 = sin(omega1*tVec);
expected_cos2 = cos(omega2*tVec);
expected_sin2 = sin(omega2*tVec);

err1 = max(abs(x(1,:) - expected_cos1));
err2 = max(abs(x(2,:) - expected_sin1));
err3 = max(abs(x(3,:) - expected_cos2));
err4 = max(abs(x(4,:) - expected_sin2));

assert(err1 < 1e-9, 'Test 4 failed: state 1 does not match cos(omega1*t), error=%.3e', err1);
assert(err2 < 1e-9, 'Test 4 failed: state 2 does not match sin(omega1*t), error=%.3e', err2);
assert(err3 < 1e-9, 'Test 4 failed: state 3 does not match cos(omega2*t), error=%.3e', err3);
assert(err4 < 1e-9, 'Test 4 failed: state 4 does not match sin(omega2*t), error=%.3e', err4);
fprintf('PASS: simulated trajectory matches cos/sin(omega*t) exactly for both frequencies\n');
fprintf('  (max errors: %.3e, %.3e, %.3e, %.3e over %d steps)\n', err1, err2, err3, err4, nSteps);

fprintf('\nAll 4 tests passed.\n');