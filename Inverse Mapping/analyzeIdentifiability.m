%% analyzeIdentifiability.m
clear; clc;

%% Load saved Jacobians
load('Jacobian_eps1e-3.mat', 'J1', 'info1');
load('Jacobian_eps1e-4.mat', 'J2', 'info2');
load('Jacobian_eps1e-5.mat', 'J3', 'info3');

%% Reference parameters
beta0  = 4.6985e-5;
gamma0 = 1.4548e-3;
R0     = 0.3830;

p0 = [beta0; gamma0; R0];

%% Check that baseline feature vectors are consistent
rel12_f0 = norm(info2.f0 - info1.f0) / norm(info2.f0);
rel23_f0 = norm(info3.f0 - info2.f0) / norm(info3.f0);

fprintf('Baseline f0 consistency:\n');
fprintf('  eps1 vs eps2: %.3e\n', rel12_f0);
fprintf('  eps2 vs eps3: %.3e\n\n', rel23_f0);

f0 = info2.f0;

%% Dimensionless scaling
Dp = diag(p0);
Df = diag(f0);

Js1 = Df \ (J1 * Dp);
Js2 = Df \ (J2 * Dp);
Js3 = Df \ (J3 * Dp);

%% SVDs
[U1,S1,V1] = svd(Js1, 'econ');
[U2,S2,V2] = svd(Js2, 'econ');
[U3,S3,V3] = svd(Js3, 'econ');

s1 = diag(S1);
s2 = diag(S2);
s3 = diag(S3);

fprintf('Scaled singular values:\n');
fprintf('  eps=1e-3: %s\n', mat2str(s1.', 6));
fprintf('  eps=1e-4: %s\n', mat2str(s2.', 6));
fprintf('  eps=1e-5: %s\n\n', mat2str(s3.', 6));

fprintf('Scaled condition numbers:\n');
fprintf('  eps=1e-3: %.6g\n', s1(1)/s1(end));
fprintf('  eps=1e-4: %.6g\n', s2(1)/s2(end));
fprintf('  eps=1e-5: %.6g\n\n', s3(1)/s3(end));

%% Weakest parameter direction
fprintf('Weakest right singular vectors:\n');
fprintf('  eps=1e-3: %s\n', mat2str(V1(:,3).', 6));
fprintf('  eps=1e-4: %s\n', mat2str(V2(:,3).', 6));
fprintf('  eps=1e-5: %s\n\n', mat2str(V3(:,3).', 6));

%% Sign-invariant agreement of weakest directions
c12 = abs(V1(:,3)' * V2(:,3));
c23 = abs(V2(:,3)' * V3(:,3));
c13 = abs(V1(:,3)' * V3(:,3));

fprintf('Weak-direction alignment |v_i^T v_j|:\n');
fprintf('  1e-3 vs 1e-4: %.8f\n', c12);
fprintf('  1e-4 vs 1e-5: %.8f\n', c23);
fprintf('  1e-3 vs 1e-5: %.8f\n\n', c13);

%% Working result
Js = Js2;
U  = U2;
S  = S2;
V  = V2;

disp('Working scaled Jacobian Js (epsilon = 1e-4):')
disp(Js)