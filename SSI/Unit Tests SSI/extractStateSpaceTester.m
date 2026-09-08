% Need access to the full Forward Model 
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'JONSWAP'));
clear all, close all, clc

%% Test 1: noise-free sanity check with a hand-built, exactly-known system
% A simple 2-state oscillator, no noise at all -> if SSI can't recover
% THIS exactly, nothing downstream can be trusted at all.
fprintf('=== Test 1: noise-free known-system recovery ===\n');
omega0 = 1.3; zeta0 = 0.05; dt = 0.05;
Ac_true = [0 1; -omega0^2 -2*zeta0*omega0];
A_true = expm(Ac_true*dt);
C_true = [1 0];   % observe position only

nSamples = 2000;
x = zeros(2, nSamples);
x(:,1) = [1; 0];
for k = 1:nSamples-1
    x(:,k+1) = A_true*x(:,k);
end
y = C_true*x;   % 1 x nSamples, noise-free

i = 20;
[Yp_ref, Yf, ~] = buildHankelMatrix(y, i);
P = hankelProjection(Yp_ref, Yf);
[A_id, C_id, sv] = extractStateSpace(P, 1, 2);

lam_true = eig(A_true);
lam_id = eig(A_id);
fprintf('True discrete eigenvalues:      %s\n', mat2str(lam_true, 6));
fprintf('Identified discrete eigenvalues: %s\n', mat2str(sort(lam_id), 6));

%% Test 2: real twin test against the validated forward model
fprintf('\n=== Test 2: twin test against Forward Model synthetic data ===\n');
Hs = 0.1; w0 = 5.5; gammaJONSWAP = 3.3;
wMin = 5.49; wMax = 5.51; numBins = 1;   % near-single-frequency, as before

H = 1.88; beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
M = 50; P_trunc = 10; N = 10;

specData = precomputeSpectralData(Hs, w0, gammaJONSWAP, wMin, wMax, numBins, ...
    H, beta, gamma, R, nu, M, P_trunc, N);

omega_true = specData.omega(1);
fprintf('True forcing angular frequency: %.6f rad/s\n', omega_true);

rTest = 0.7*R; thetaTest = pi/4;
dt2 = 0.05;
tVec = 0:dt2:200;
zeta_ts = evaluateSpectralDeflection(specData, rTest, thetaTest, tVec);

i2 = 30;
[Yp_ref2, Yf2, ~] = buildHankelMatrix(zeta_ts, i2);
P2 = hankelProjection(Yp_ref2, Yf2);
[A_id2, ~, sv2] = extractStateSpace(P2, 1, 2);

lam_id2 = eig(A_id2);
identified_omega = abs(log(lam_id2(1))/dt2);
fprintf('Identified angular frequency (from poles): %.6f rad/s\n', identified_omega);
fprintf('Relative error: %.4f%%\n', 100*abs(identified_omega-omega_true)/omega_true);