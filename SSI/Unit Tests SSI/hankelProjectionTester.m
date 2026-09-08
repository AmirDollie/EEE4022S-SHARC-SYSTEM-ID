%% hankelProjectionTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

%% Build a small test case via buildHankelMatrix, reusing the earlier example
y = [1 2 3 4 5 6 7 8;
     10 20 30 40 50 60 70 80];
i = 2;
[Yp_ref, Yf, ~] = buildHankelMatrix(y, i);

%% Test 1: fast (QR-based) projection vs. the direct textbook formula
[P_fast, ~, ~] = hankelProjection(Yp_ref, Yf);

P_direct = Yf * Yp_ref' * pinv(Yp_ref * Yp_ref') * Yp_ref;

fprintf('=== Test 1: QR-based projection vs. direct formula ===\n');
fprintf('max|P_fast - P_direct| = %.3e (should be ~0)\n', ...
    max(abs(P_fast(:) - P_direct(:))));

%% Test 2: R is genuinely lower triangular (sanity check on the QR step itself)
r_i = size(Yp_ref,1);
H = [Yp_ref; Yf];
[Qfull, Rfull] = qr(H', 0);
R = Rfull';
upperPart = triu(R, 1);   % strictly-upper entries; should all be ~0
fprintf('\n=== Test 2: R is lower triangular ===\n');
fprintf('max|strictly-upper entries of R| = %.3e (should be ~0)\n', max(abs(upperPart(:))));

%% Test 3: H reconstructs correctly from R and Q (did the QR step itself work?)
Q1_full = Qfull(:, 1:r_i);
Q2_full = Qfull(:, r_i+1:end);
Q_check = [Q1_full, Q2_full];
H_reconstructed = R * Q_check';
fprintf('\n=== Test 3: H reconstructs from R*Q'' ===\n');
fprintf('max|H - R*Q''| = %.3e (should be ~0)\n', max(abs(H(:) - H_reconstructed(:))));

%% Test 4: projection is idempotent (a real projection matrix satisfies P*P^T's
%  column space = P's column space; a cheap proxy check: projecting Yf onto
%  Yp_ref twice should give the same result as projecting once)
fprintf('\n=== Test 4: projection is stable under re-projection ===\n');
[P_fast2, ~, ~] = hankelProjection(Yp_ref, P_fast);
fprintf('max|P_fast - P_fast2| = %.3e (should be ~0, P is already in the target row space)\n', ...
    max(abs(P_fast(:) - P_fast2(:))));