%% sweepModelOrdersTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

%% Build a known 2-mode synthetic system (correct-order regime)
% Two independent damped oscillators, summed into ONE output channel --
% a genuinely multi-mode case, unlike the single-oscillator tests used
% for extractStateSpace/modalParameters individually.
fprintf('=== Test 1: correct-order regime, 2 known modes ===\n');
omega1 = 1.2; zeta1 = 0.02;
omega2 = 3.7; zeta2 = 0.04;
dt = 0.03;

Ac1 = [0 1; -omega1^2 -2*zeta1*omega1];
Ac2 = [0 1; -omega2^2 -2*zeta2*omega2];
A_true = blkdiag(expm(Ac1*dt), expm(Ac2*dt));
C_true = [1 0 1 0];   % single sensor, sees both modes summed

nSamples = 4000;
x = zeros(4, nSamples);
x(:,1) = [1; 0; 0.5; 0];
for k = 1:nSamples-1
    x(:,k+1) = A_true*x(:,k);
end
y = C_true*x;

i = 30;
[Yp_ref, Yf, ~] = buildHankelMatrix(y, i);
P = hankelProjection(Yp_ref, Yf);

nMax = 8;
results = sweepModelOrders(P, 1, dt, nMax);

%% Test 1a: check deduplication -- results(4) should have at most 2 modes
% (true order is 4, i.e. 2 conjugate pairs -> 2 distinct modes after dedup)
fprintf('At order n=4 (true order): found %d distinct mode(s) (expect 2)\n', ...
    length(results(4).frequency));
fprintf('Frequencies found: %s (true: [%.4f, %.4f])\n', ...
    mat2str(sort(results(4).frequency), 6), omega1, omega2);

%% Test 1b: confirm no duplicate (near-identical) frequencies survived
% i.e. dedup actually removed the conjugate partner, not just coincidence
fprintf('\n=== Test 1b: confirm no duplicate frequencies at any order ===\n');
anyDuplicates = false;
for n = 1:nMax
    f = sort(results(n).frequency);
    if length(f) > 1
        gaps = diff(f);
        if any(gaps < 1e-6)
            anyDuplicates = true;
            fprintf('  order %d: found near-duplicate frequencies! %s\n', n, mat2str(f));
        end
    end
end
if ~anyDuplicates
    fprintf('No duplicate frequencies found at any order (dedup working correctly)\n');
end

%% Test 2: struct array shapes are self-consistent at every order
fprintf('\n=== Test 2: internal shape consistency at every order ===\n');
allConsistent = true;
for n = 1:nMax
    k = length(results(n).frequency);
    okDamping = length(results(n).damping) == k;
    okShapes  = size(results(n).modeShapes, 2) == k;
    okPoles   = length(results(n).poles) == k;
    if ~(okDamping && okShapes && okPoles)
        allConsistent = false;
        fprintf('  order %d: MISMATCH (freq=%d, damp=%d, shapes=%d, poles=%d)\n', ...
            n, k, length(results(n).damping), size(results(n).modeShapes,2), length(results(n).poles));
    end
end
if allConsistent
    fprintf('All %d orders have internally consistent field lengths\n', nMax);
end

%% Test 3: number of distinct modes never exceeds floor(n/2)
% A real system of order n can have at most floor(n/2) conjugate pairs
fprintf('\n=== Test 3: mode count never exceeds floor(n/2) ===\n');
boundOK = true;
for n = 1:nMax
    k = length(results(n).frequency);
    maxPossible = floor(n/2);
    if k > maxPossible
        boundOK = false;
        fprintf('  order %d: found %d modes, but max possible is %d!\n', n, k, maxPossible);
    end
end
if boundOK
    fprintf('Mode count respects the floor(n/2) bound at every order\n');
end