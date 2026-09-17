%% buildCandidateCrowsTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

%% Test 1: known algebraic case, hand-computable
fprintf('=== Test 1: hand-computable case ===\n');
phi1 = 2 + 3i;
phi2 = -1 + 0.5i;
C_j = buildCandidateCrows(phi1, phi2);
expected = [2, -3, -1, -0.5];
assert(isequal(size(C_j), [1 4]), 'Test 1 failed: wrong shape.');
assert(max(abs(C_j - expected)) < 1e-14, 'Test 1 failed: C_j does not match hand-computed expected value.');
fprintf('PASS: C_j = %s (matches hand-computed [2, -3, -1, -0.5])\n', mat2str(C_j));

%% Test 2: the real end-to-end check -> does y=C*x actually reproduce
% Re{phi*e^{i*omega*t}}, computed independently via complex arithmetic?

% This is the test that would catch a state-ordering or sign error that
% Test 1 alone (a single hand-picked numeric case) might not expose if
% the error happened to cancel out for that particular phi.
fprintf('\n=== Test 2: y=C*x matches Re{phi*e^{i*omega*t}} independently ===\n');
omega1 = 2.70; omega2 = 6.00; dt = 0.02;
phi1 = 0.8 - 1.2i;
phi2 = -0.3 + 0.6i;

A = buildKnownOscillatorA(omega1, omega2, dt);
C_j = buildCandidateCrows(phi1, phi2);

nSteps = 300;
x = zeros(4, nSteps);
x(:,1) = [1; 0; 1; 0];
for k = 1:nSteps-1
    x(:,k+1) = A*x(:,k);
end
y = C_j * x;   % 1 x nSteps

tVec = (0:nSteps-1)*dt;
yExpected = real(phi1*exp(1i*omega1*tVec) + phi2*exp(1i*omega2*tVec));

err = max(abs(y - yExpected));
assert(err < 1e-9, 'Test 2 failed: y=C*x does not match the independently-computed physical signal, error=%.3e', err);
fprintf('PASS: y=C*x matches Re{phi1*e^{i*omega1*t}+phi2*e^{i*omega2*t}} to within %.3e\n', err);

%% Test 3: multiple candidate rows stack correctly into a valid C matrix
fprintf('\n=== Test 3: multiple candidate rows stack into a valid multi-sensor C ===\n');
sensorPhis = [1+0.5i, 0.2-0.1i;
              -0.4+0.3i, 0.9+0.2i;
              0.1-0.7i, -0.5-0.5i];   % 3 sensors x 2 modes
C_multi = zeros(3,4);
for s = 1:3
    C_multi(s,:) = buildCandidateCrows(sensorPhis(s,1), sensorPhis(s,2));
end
assert(isequal(size(C_multi), [3 4]), 'Test 3 failed: wrong shape for multi-sensor C.');
y_multi = C_multi * x;
for s = 1:3
    yExpected_s = real(sensorPhis(s,1)*exp(1i*omega1*tVec) + sensorPhis(s,2)*exp(1i*omega2*tVec));
    errS = max(abs(y_multi(s,:) - yExpected_s));
    assert(errS < 1e-9, 'Test 3 failed: sensor %d row does not match its expected signal, error=%.3e', s, errS);
end
fprintf('PASS: all 3 stacked sensor rows independently match their expected signals\n');

fprintf('\nAll 3 tests passed.\n');