%% computeBranchPersistenceTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

windowSize = 3; minStableInWindow = 2;   % keep at the ORIGINAL values -- do not
                                          % lower windowSize as a fix; see the
                                          % diagnostic discussion below

%% Test 1: one non-class-3 transition inside an otherwise class-3-stable
% run should NOT break persistence (the non-class-3 transition could be
% class 0, 1, or 2 -- none of them count toward N_3, but a single one
% among several class-3 neighbours should still leave a qualifying window)
fprintf('=== Test 1: one non-class-3 transition should NOT break persistence ===\n');
branches(1).order = [1 2 3 4 5 6];
branches(1).class = [NaN 3 3 0 3 3];   % one non-class-3 transition at position 4

p1 = computeBranchPersistence(branches(1), windowSize, minStableInWindow);
assert(p1.isPersistent == true, 'Test 1 failed: a single non-class-3 transition should not prevent persistence.');
fprintf('PASS: isPersistent = %d, best window starts at order %d with count %d\n', ...
    p1.isPersistent, p1.bestWindowStart, p1.bestWindowCount);

%% Test 2: a genuinely unstable branch (mostly class 0) should NOT be persistent
fprintf('\n=== Test 2: mostly-unstable branch should NOT be persistent ===\n');
branches2(1).order = [1 2 3 4 5 6];
branches2(1).class = [NaN 0 3 0 0 3];   % only isolated single hits, never 2-in-a-row-of-3

p2 = computeBranchPersistence(branches2(1), windowSize, minStableInWindow);
assert(p2.isPersistent == false, 'Test 2 failed: sparse class-3 hits should not count as persistent.');
fprintf('PASS: isPersistent = %d (correctly not persistent)\n', p2.isPersistent);

%% Test 3: too-short branch (fewer transitions than windowSize) -> not persistent, no error
fprintf('\n=== Test 3: branch shorter than windowSize ===\n');
branches3(1).order = [1 2];
branches3(1).class = [NaN 3];   % only 1 transition, windowSize=3

p3 = computeBranchPersistence(branches3(1), windowSize, minStableInWindow);
assert(p3.isPersistent == false, 'Test 3 failed: too-short branch should not be persistent.');
assert(isnan(p3.bestWindowStart), 'Test 3 failed: bestWindowStart should be NaN for too-short branch.');
fprintf('PASS: isPersistent = %d, bestWindowStart = %s (correctly NaN)\n', p3.isPersistent, mat2str(p3.bestWindowStart));

%% Test 4: perfectly stable branch -> best window count should equal windowSize exactly
fprintf('\n=== Test 4: fully class-3 branch, best window should be maximal ===\n');
branches4(1).order = [1 2 3 4 5];
branches4(1).class = [NaN 3 3 3 3];

p4 = computeBranchPersistence(branches4(1), windowSize, minStableInWindow);
assert(p4.bestWindowCount == windowSize, 'Test 4 failed: fully stable branch should have bestWindowCount == windowSize.');
fprintf('PASS: bestWindowCount = %d (equals windowSize=%d)\n', p4.bestWindowCount, windowSize);

%% Test 5: apply to the actual Pass 1 noisy-data branches -- exploratory,
% not pass/fail, to see how the sliding window behaves on the case that
% broke the strict-consecutiveness check
fprintf('\n=== Test 5: sliding window applied to the noisy Pass 1 branches (exploratory) ===\n');
freqTol = 0.01; dampTol = 0.05; macTol = 0.99; dt = 0.03;
omega_true = [1.2; 2.7; 4.1]; zeta_true = [0.02; 0.015; 0.03];
Ac = cell(3,1);
for q = 1:3
    Ac{q} = [0 1; -omega_true(q)^2 -2*zeta_true(q)*omega_true(q)];
end
A_true = blkdiag(expm(Ac{1}*dt), expm(Ac{2}*dt), expm(Ac{3}*dt));
C_true = [1 0 0.8 0 0.5 0];
nSamples = 6000;
x = zeros(6, nSamples); x(:,1) = [1;0; 0.7;0; 0.4;0];
for k = 1:nSamples-1, x(:,k+1) = A_true*x(:,k); end
rng(42);
noiseLevel = 0.01 * std(C_true*x);
y_noisy = C_true*x + noiseLevel*randn(1, nSamples);

i = 40; nMax = 14;
[Yp_ref, Yf, ~] = buildHankelMatrix(y_noisy, i);
P = hankelProjection(Yp_ref, Yf);
r5 = sweepModelOrders(P, 1, dt, nMax);
m5 = cell(nMax,1);
for n = 2:nMax
    m5{n} = classifyPoleStability(r5(n-1), r5(n), freqTol, dampTol, macTol);
end
branches5 = buildModalBranches(r5, m5);
p5 = computeBranchPersistence(branches5, windowSize, minStableInWindow);

for q = 1:3
    fprintf('  f=%.4f rad/s:\n', omega_true(q));
    for b = 1:length(branches5)
        if any(abs(branches5(b).frequency - omega_true(q)) / omega_true(q) < freqTol*3)
            fprintf('    branch %d: isPersistent=%d, bestWindowCount=%d/%d, spans orders [%d,%d]\n', ...
                b, p5(b).isPersistent, p5(b).bestWindowCount, windowSize, ...
                branches5(b).order(1), branches5(b).order(end));
            fprintf('    classes: %s\n', mat2str(branches5(b).class));
        end
    end
end

%% Diagnostic: freqTol sweep on the SAME noisy dataset, with raw per-order
% frequencies printed near each true mode, to distinguish "no compatible
% pole found" (cause 1: SSI estimate itself moved too far) from
% "compatible pole found but rejected/degraded" (causes 2-4: matching or
% extraction issue). This does NOT modify computeBranchPersistence,
% classifyPoleStability, or buildModalBranches -- it only re-runs the
% EXISTING, already-validated chain at different freqTol values to see
% where fragmentation actually comes from, per the discussion that
% lowering windowSize would treat the symptom, not the cause.
fprintf('\n=== Diagnostic: freqTol sweep on noisy data (NOT a fix, an experiment) ===\n');
for freqTolTest = [0.01, 0.02, 0.05]
    fprintf('\n--- freqTol = %.2f ---\n', freqTolTest);
    m5b = cell(nMax,1);
    for n = 2:nMax
        m5b{n} = classifyPoleStability(r5(n-1), r5(n), freqTolTest, dampTol, macTol);
    end
    branches5b = buildModalBranches(r5, m5b);

    for q = 1:3
        fprintf('  f=%.4f rad/s:\n', omega_true(q));
        for n = 6:nMax
            fVals = r5(n).frequency;
            nearby = fVals(abs(fVals - omega_true(q))/omega_true(q) < 0.1);   % generous, just to SEE what's near
            fprintf('    n=%2d: nearby frequencies = %s\n', n, mat2str(nearby, 4));
        end
        for b = 1:length(branches5b)
            if any(abs(branches5b(b).frequency - omega_true(q)) / omega_true(q) < freqTolTest*3)
                fprintf('    branch %d: spans [%d,%d], classes=%s\n', ...
                    b, branches5b(b).order(1), branches5b(b).order(end), mat2str(branches5b(b).class));
            end
        end
    end
end