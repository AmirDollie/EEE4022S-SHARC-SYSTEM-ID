%% classifyPoleStabilityTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

freqTol = 0.01; dampTol = 0.05; macTol = 0.99;

%% Test 1: perfect match (identical results) should all be class 3
fprintf('=== Test 1: identical consecutive orders -> all class 3 ===\n');
prevR.frequency = [1.2; 3.7];
prevR.damping = [0.02; 0.04];
prevR.modeShapes = [1 0.3; 0.5 1];

currR = prevR;

m1 = classifyPoleStability(prevR, currR, freqTol, dampTol, macTol);
assert(isequal([m1.class]', [3;3]), 'Test 1 failed: expected classes [3;3]');
assert(isequal([m1.previousIndex]', [1;2]), 'Test 1 failed: expected previousIndex [1;2]');
fprintf('PASS: Classes = %s, previousIndex = %s\n', mat2str([m1.class]'), mat2str([m1.previousIndex]'));

%% Test 2: one pole drifts in frequency beyond tolerance -> class 0
fprintf('\n=== Test 2: frequency drift beyond tolerance -> unmatched ===\n');
currR2 = prevR;
currR2.frequency(1) = prevR.frequency(1) * 1.05;
m2 = classifyPoleStability(prevR, currR2, freqTol, dampTol, macTol);
assert(isequal([m2.class]', [0;3]), 'Test 2 failed: expected classes [0;3]');
fprintf('PASS: Classes = %s\n', mat2str([m2.class]'));

%% Test 3: frequency+damping OK, MAC fails -> class 2 (not just "<3")
fprintf('\n=== Test 3: MAC fails, damping unchanged -> class 2 exactly ===\n');
currR3 = prevR;
currR3.modeShapes(:,1) = [0.1; -1];   % MAC(prev,curr) for pole 1 ~ 0.0016, fails macTol
m3 = classifyPoleStability(prevR, currR3, freqTol, dampTol, macTol);
assert(isequal([m3.class]', [2;3]), 'Test 3 failed: expected classes [2;3]');
fprintf('PASS: Classes = %s\n', mat2str([m3.class]'));

%% Test 4: one-to-one assignment -- the BETTER-MAC candidate must win,
%  not merely "no double claim"
fprintf('\n=== Test 4: one-to-one assignment, best-MAC candidate wins ===\n');
prevR4.frequency = [2.0];
prevR4.damping = [0.03];
prevR4.modeShapes = [1; 0.5];

currR4.frequency = [2.001; 2.002];
currR4.damping = [0.0301; 0.0299];
currR4.modeShapes = [1 1; 0.5 0.9];   % pole 1's shape is closer to prevR4's

m4 = classifyPoleStability(prevR4, currR4, freqTol, dampTol, macTol);
prevIdx = [m4.previousIndex]';
assert(sum(prevIdx==1) <= 1, 'Test 4 failed: prev pole 1 claimed by more than one current pole');
assert(prevIdx(1)==1 && prevIdx(2)==0, 'Test 4 failed: expected the better-MAC pole (1) to win, got previousIndex=%s', mat2str(prevIdx));
fprintf('PASS: previousIndex = %s (pole 1, the better MAC match, correctly won)\n', mat2str(prevIdx));

%% Test 5: empty previous order (n=1 case) -> all class 0
fprintf('\n=== Test 5: empty previous results -> all class 0 ===\n');
emptyPrev.frequency = [];
emptyPrev.damping = [];
emptyPrev.modeShapes = [];
m5 = classifyPoleStability(emptyPrev, prevR, freqTol, dampTol, macTol);
assert(isequal([m5.class]', [0;0]), 'Test 5 failed: expected classes [0;0]');
fprintf('PASS: Classes = %s\n', mat2str([m5.class]'));

%% Test 6: frequency OK, damping fails, MAC also fails -> class 1
fprintf('\n=== Test 6: frequency OK, damping fails -> class 1 ===\n');
currR6 = prevR;
currR6.damping(1) = 0.10;
currR6.modeShapes(:,1) = [0.1; -1];
m6 = classifyPoleStability(prevR, currR6, freqTol, dampTol, macTol);
assert(isequal([m6.class]', [1;3]), 'Test 6 failed: expected classes [1;3]');
fprintf('PASS: Classes = %s\n', mat2str([m6.class]'));

fprintf('\nAll 6 tests passed.\n');