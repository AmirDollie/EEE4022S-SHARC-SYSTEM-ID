%% buildHankelMatrixTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

%% Tiny hand-computable case: l=2 channels, n=6 samples, i=2
y = [1 2 3 4 5 6;
     10 20 30 40 50 60];
i = 2;

%% Test 1: default refIdx (all channels), check shapes
[Yp_ref, Yf, H] = buildHankelMatrix(y, i);
fprintf('=== Test 1: shapes (default refIdx = all channels) ===\n');
l = 2; r = 2; j = size(y,2) - 2*i + 1;   % j = 6-4+1 = 3
fprintf('Yp_ref: got %dx%d, expected %dx%d\n', size(Yp_ref,1), size(Yp_ref,2), r*i, j);
fprintf('Yf:     got %dx%d, expected %dx%d\n', size(Yf,1), size(Yf,2), l*i, j);
fprintf('H:      got %dx%d, expected %dx%d\n', size(H,1), size(H,2), (r+l)*i, j);

%% Test 2: hand-computed values, unscaled, to check the actual entries
% With i=2, j=3: Yp_ref block rows are y_ref shifted by 0 and 1 columns
% (columns 1:3 and 2:4 of y_ref); Yf block rows are shifted by 2 and 3
% (columns 3:5 and 4:6).
scale = 1/sqrt(j);
expected_Yp_ref_unscaled = [1 2 3; 10 20 30; 2 3 4; 20 30 40];   % blocks k=0,1, both channels
expected_Yf_unscaled     = [3 4 5; 30 40 50; 4 5 6; 40 50 60];   % blocks k=2,3, both channels

fprintf('\n=== Test 2: hand-computed entries (unscaled) ===\n');
fprintf('max|Yp_ref/scale - expected| = %.3e (should be ~0)\n', ...
    max(abs(Yp_ref(:)/scale - expected_Yp_ref_unscaled(:))));
fprintf('max|Yf/scale - expected|     = %.3e (should be ~0)\n', ...
    max(abs(Yf(:)/scale - expected_Yf_unscaled(:))));

%% Test 3: refIdx actually restricts which channels feed Yp_ref
[Yp_ref_ch1, ~, ~] = buildHankelMatrix(y, i, 1);   % only channel 1 as reference
fprintf('\n=== Test 3: refIdx=1 restricts Yp_ref to channel 1 only ===\n');
fprintf('Yp_ref size with refIdx=1: %dx%d (expect %dx%d, r=1)\n', ...
    size(Yp_ref_ch1,1), size(Yp_ref_ch1,2), 1*i, j);
expected_ch1_unscaled = [1 2 3; 2 3 4];
fprintf('max|Yp_ref_ch1/scale - expected| = %.3e (should be ~0)\n', ...
    max(abs(Yp_ref_ch1(:)/scale - expected_ch1_unscaled(:))));

%% Test 4: input validation actually fires
fprintf('\n=== Test 4: error handling ===\n');
try
    buildHankelMatrix(y, 5);   % i=5 needs n>=10, but n=6
    fprintf('ERROR: should have thrown for too-few-samples, but did not\n');
catch err
    fprintf('Correctly threw: %s\n', err.identifier);
end
try
    buildHankelMatrix(y, i, 3);   % refIdx=3 invalid, only 2 channels
    fprintf('ERROR: should have thrown for bad refIdx, but did not\n');
catch err
    fprintf('Correctly threw: %s\n', err.identifier);
end