%% experimentsAB_freqTolValidation.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

freqTol = 0.05; dampTol = 0.05; macTol = 0.99;
windowSize = 3; minStableInWindow = 2;   % ORIGINAL values, unchanged
dt = 0.03;

%% Sanity check: confirm MAC is trivially 1 for single-channel mode shapes
fprintf('=== Sanity check: MAC triviality at l=1 ===\n');
scalarPairs = {[1+0i, 5-3i], [0.001, 1000], [2i, -7]};
for k = 1:length(scalarPairs)
    a = scalarPairs{k}(1); b = scalarPairs{k}(2);
    fprintf('  MAC(%s, %s) = %.10f\n', num2str(a), num2str(b), computeMAC(a,b));
end
fprintf('(All exactly 1 -- confirms MAC provides NO discrimination with only 1 sensor channel.\n');
fprintf(' Every prior single-channel test has effectively been frequency+damping matching only.)\n');

%% Experiment A: confirm persistence at freqTol=0.05, ORIGINAL windowSize/minStableInWindow
fprintf('\n=== Experiment A: persistence at freqTol=0.05 (original W=3,k=2) ===\n');
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
rA = sweepModelOrders(P, 1, dt, nMax);
mA = cell(nMax,1);
for n = 2:nMax
    mA{n} = classifyPoleStability(rA(n-1), rA(n), freqTol, dampTol, macTol);
end
branchesA = buildModalBranches(rA, mA);
pA = computeBranchPersistence(branchesA, windowSize, minStableInWindow);

for q = 1:3
    foundPersistent = false;
    for b = 1:length(branchesA)
        if any(abs(branchesA(b).frequency - omega_true(q)) / omega_true(q) < freqTol) && pA(b).isPersistent
            foundPersistent = true;
            fprintf('  f=%.4f rad/s: PERSISTENT (branch %d, spans [%d,%d], bestWindowCount=%d/%d)\n', ...
                omega_true(q), b, branchesA(b).order(1), branchesA(b).order(end), ...
                pA(b).bestWindowCount, windowSize);
        end
    end
    if ~foundPersistent
        fprintf('  f=%.4f rad/s: NOT identified as persistent\n', omega_true(q));
    end
end

%% Experiment B1: two CLOSE true modes, SINGLE channel (l=1) -- tests
% frequency+damping discrimination ONLY, since MAC is trivial here
fprintf('\n=== Experiment B1: two close modes (2.70, 2.80 rad/s), single channel ===\n');
omega_B = [2.70; 2.80]; zeta_B = [0.015; 0.04];   % distinct damping, close frequency
fprintf('Frequency separation: %.2f%% (relative to 2.70)\n', 100*(omega_B(2)-omega_B(1))/omega_B(1));

AcB1 = [0 1; -omega_B(1)^2 -2*zeta_B(1)*omega_B(1)];
AcB2 = [0 1; -omega_B(2)^2 -2*zeta_B(2)*omega_B(2)];
A_trueB = blkdiag(expm(AcB1*dt), expm(AcB2*dt));
C_trueB1 = [1 0 1 0];   % single channel, sees both modes equally -- MAC cannot help here

xB = zeros(4, nSamples); xB(:,1) = [1;0; 1;0];
for k = 1:nSamples-1, xB(:,k+1) = A_trueB*xB(:,k); end
rng(7);
noiseLevelB = 0.01 * std(C_trueB1*xB);
yB1 = C_trueB1*xB + noiseLevelB*randn(1, nSamples);

[Yp_refB1, YfB1, ~] = buildHankelMatrix(yB1, i);
PB1 = hankelProjection(Yp_refB1, YfB1);
rB1 = sweepModelOrders(PB1, 1, dt, nMax);
mB1 = cell(nMax,1);
for n = 2:nMax
    mB1{n} = classifyPoleStability(rB1(n-1), rB1(n), freqTol, dampTol, macTol);
end
branchesB1 = buildModalBranches(rB1, mB1);

fprintf('Branches found near either true frequency:\n');
mergedDetected = false;
for b = 1:length(branchesB1)
    nearBoth = any(abs(branchesB1(b).frequency - omega_B(1))/omega_B(1) < freqTol) && ...
               any(abs(branchesB1(b).frequency - omega_B(2))/omega_B(2) < freqTol);
    if nearBoth
        mergedDetected = true;
        fprintf('  *** branch %d spans BOTH true frequencies -- possible wrongful merge! freqs=%s\n', ...
            b, mat2str(branchesB1(b).frequency, 4));
    elseif any(abs(branchesB1(b).frequency - omega_B(1))/omega_B(1) < freqTol) || ...
           any(abs(branchesB1(b).frequency - omega_B(2))/omega_B(2) < freqTol)
        fprintf('  branch %d: spans [%d,%d], freqs=%s, damping=%s\n', ...
            b, branchesB1(b).order(1), branchesB1(b).order(end), ...
            mat2str(branchesB1(b).frequency,4), mat2str(branchesB1(b).damping,4));
    end
end
if ~mergedDetected
    fprintf('No single branch spans both true frequencies -- damping difference alone kept them separate here.\n');
end

%% Experiment B2: SAME two close modes, but TWO channels with genuinely
% different mode-shape weighting -- this is where MAC can actually help
fprintf('\n=== Experiment B2: same two close modes, 2-channel sensor (MAC now meaningful) ===\n');
C_trueB2 = [1 0 0.2 0;    % channel 1: mostly mode 1
            0.2 0 1 0];   % channel 2: mostly mode 2

yB2 = C_trueB2*xB + noiseLevelB*randn(2, nSamples);

[Yp_refB2, YfB2, ~] = buildHankelMatrix(yB2, i);
PB2 = hankelProjection(Yp_refB2, YfB2);
rB2 = sweepModelOrders(PB2, 2, dt, nMax);
mB2 = cell(nMax,1);
for n = 2:nMax
    mB2{n} = classifyPoleStability(rB2(n-1), rB2(n), freqTol, dampTol, macTol);
end
branchesB2 = buildModalBranches(rB2, mB2);

fprintf('Branches found near either true frequency:\n');
for b = 1:length(branchesB2)
    if any(abs(branchesB2(b).frequency - omega_B(1))/omega_B(1) < freqTol) || ...
       any(abs(branchesB2(b).frequency - omega_B(2))/omega_B(2) < freqTol)
        fprintf('  branch %d: spans [%d,%d], freqs=%s\n', ...
            b, branchesB2(b).order(1), branchesB2(b).order(end), mat2str(branchesB2(b).frequency,4));
    end
end