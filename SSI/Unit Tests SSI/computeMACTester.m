%% 
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

%% Test 1: MAC of a shape against itself must be exactly 1
fprintf('=== Test 1: self-MAC = 1 ===\n');
phi = [1+2i; -3; 0.5i; 4-1i];
m = computeMAC(phi, phi);
fprintf('MAC(phi,phi) = %.10f (should be exactly 1)\n', m);

%% Test 2: MAC of orthogonal shapes must be exactly 0
fprintf('\n=== Test 2: orthogonal shapes give MAC = 0 ===\n');
phi_a = [1; 0; 0];
phi_b = [0; 1; 0];
m2 = computeMAC(phi_a, phi_b);
fprintf('MAC(orthogonal) = %.10f (should be exactly 0)\n', m2);

%% Test 3: MAC is invariant to arbitrary complex scaling
% This is the whole point of MAC -- output-only mode shapes are only
% known up to scale, so scaling one shape by any complex number should
% not change the MAC value at all.
fprintf('\n=== Test 3: scale invariance ===\n');
phi_c = [2-1i; 0.3; -1.5i; 3+3i];
scaleFactor = 5*exp(1i*0.7);   % arbitrary complex scale
m3_unscaled = computeMAC(phi, phi_c);
m3_scaled = computeMAC(phi, scaleFactor*phi_c);
fprintf('MAC unscaled = %.10f, MAC scaled = %.10f (should match exactly)\n', ...
    m3_unscaled, m3_scaled);

%% Test 4: explicitly confirm the squared-numerator bug is NOT present
% This is the specific error found in a third-party implementation
% during this project: computing |phi_a^H*phi_b| instead of the correct
% |phi_a^H*phi_b|^2. Confirm our version gives a DIFFERENT (correct)
% answer from that broken formula on a case where they visibly diverge.
fprintf('\n=== Test 4: confirm squared-numerator bug is absent ===\n');
phi_d = [1; 0.5; 0.2];
phi_e = [0.9; 0.6; 0.1];   % similar but not identical -- MAC should be
                            % well below 1, making the bug's effect visible
mac_correct = computeMAC(phi_d, phi_e);
mac_buggy = abs(phi_d'*phi_e) / sqrt(abs((phi_d'*phi_d)*(phi_e'*phi_e)));  % the broken (unsquared) version
fprintf('Correct MAC (squared numerator)   = %.6f\n', mac_correct);
fprintf('Buggy MAC (unsquared, sqrt(true)) = %.6f\n', mac_buggy);
fprintf('sqrt(mac_correct) = %.6f (should equal the buggy value above, confirming the relationship)\n', ...
    sqrt(mac_correct));

%% Test 5: matrix inputs return a full na x nb matrix correctly
fprintf('\n=== Test 5: matrix-input shape check ===\n');
PhiMat_a = [phi, phi_c];      % 4x2
PhiMat_b = [phi, phi_c];      % 4x2, reuse for a clean square case
Mfull = computeMAC(PhiMat_a, PhiMat_b);
fprintf('Output size: %dx%d (expect 2x2)\n', size(Mfull,1), size(Mfull,2));
fprintf('Diagonal (self-MAC, should both be 1): [%.6f, %.6f]\n', Mfull(1,1), Mfull(2,2));