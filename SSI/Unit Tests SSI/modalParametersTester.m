%% modalParametersTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

%% Test 1: known oscillator, check (omega, zeta) recovered exactly
fprintf('=== Test 1: exact recovery from a known noise-free system ===\n');
omega0_true = 1.3; zeta0_true = 0.05; dt = 0.05;

Ac_true = [0 1; -omega0_true^2 -2*zeta0_true*omega0_true];
A_true = expm(Ac_true*dt);
C_true = [1 0];

nSamples = 2000;
x = zeros(2, nSamples);
x(:,1) = [1; 0];
for k = 1:nSamples-1
    x(:,k+1) = A_true*x(:,k);
end
y = C_true*x;

i = 20;
[Yp_ref, Yf, ~] = buildHankelMatrix(y, i);
P = hankelProjection(Yp_ref, Yf);
[A_id, C_id, ~] = extractStateSpace(P, 1, 2);

[omega_id, zeta_id, Phi_id, ~] = modalParameters(A_id, C_id, dt);

fprintf('True:       omega0=%.6f rad/s, zeta0=%.6f\n', omega0_true, zeta0_true);
fprintf('Identified: omega =%.6f rad/s, zeta =%.6f (mode 1)\n', omega_id(1), zeta_id(1));
fprintf('            omega =%.6f rad/s, zeta =%.6f (mode 2)\n', omega_id(2), zeta_id(2));

%% Test 2: independent algebraic check, no SSI at all
fprintf('\n=== Test 2: direct algebraic check, bypassing SSI entirely ===\n');
omega_a = 2.0; zeta_a = 0.03;
omega_b = 5.0; zeta_b = 0.10;
dt2 = 0.02;

lambda_c_a = -zeta_a*omega_a + 1i*omega_a*sqrt(1-zeta_a^2);
lambda_c_b = -zeta_b*omega_b + 1i*omega_b*sqrt(1-zeta_b^2);
lambda_d_a = exp(lambda_c_a*dt2);
lambda_d_b = exp(lambda_c_b*dt2);

% Real 2x2 block per conjugate pair: r*e^{+/-i*theta} as eigenvalues of
% [r*cos(theta) -r*sin(theta); r*sin(theta) r*cos(theta)]. This is the
% correct way to build a REAL matrix with prescribed complex-conjugate
% eigenvalues -- using Psi=I (as before) forces real eigenvalues and
% cannot represent a genuine complex pair at all.
r_a = abs(lambda_d_a); theta_a = angle(lambda_d_a);
r_b = abs(lambda_d_b); theta_b = angle(lambda_d_b);

block_a = r_a * [cos(theta_a) -sin(theta_a); sin(theta_a) cos(theta_a)];
block_b = r_b * [cos(theta_b) -sin(theta_b); sin(theta_b) cos(theta_b)];

A_known = blkdiag(block_a, block_b);
C_known = [1 0 0 0; 0 0 1 0];

[omega_check, zeta_check, ~, ~] = modalParameters(A_known, C_known, dt2);

[omega_sorted, idx] = sort(omega_check);
zeta_sorted = zeta_check(idx);

fprintf('Built with: omega=[%.4f, %.4f], zeta=[%.4f, %.4f]\n', omega_a, omega_b, zeta_a, zeta_b);
fprintf('Recovered omega (sorted, should show %.4f and %.4f twice each): %s\n', ...
    omega_a, omega_b, mat2str(omega_sorted, 6));
fprintf('Recovered zeta, correctly paired to the sorted omega above: %s\n', mat2str(zeta_sorted, 6));

%% Test 3: real-eigenvalue warning actually fires
fprintf('\n=== Test 3: warning check for a purely real (non-oscillatory) pole ===\n');
A_real = 0.9;   % scalar, single real eigenvalue, no oscillation
C_real = 1;
lastwarn('');
[~, ~, ~, ~] = modalParameters(A_real, C_real, 0.05);
[warnMsg, warnId] = lastwarn;
fprintf('Warning fired: %d (id: %s)\n', ~isempty(warnMsg), warnId);