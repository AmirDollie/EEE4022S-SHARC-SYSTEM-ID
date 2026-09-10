%% pass1IntegrationTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

freqTol = 0.01; dampTol = 0.05; macTol = 0.99;
dt = 0.03;

%% Build a known 3-mode system, true order = 6
omega_true = [1.2; 2.7; 4.1];
zeta_true  = [0.02; 0.015; 0.03];

Ac = cell(3,1);
for q = 1:3
    Ac{q} = [0 1; -omega_true(q)^2 -2*zeta_true(q)*omega_true(q)];
end
A_true = blkdiag(expm(Ac{1}*dt), expm(Ac{2}*dt), expm(Ac{3}*dt));
C_true = [1 0 0.8 0 0.5 0];   % one sensor, all three modes summed with different weights

nSamples = 6000;
x = zeros(6, nSamples);
x(:,1) = [1;0; 0.7;0; 0.4;0];
for k = 1:nSamples-1
    x(:,k+1) = A_true*x(:,k);
end

rng(42);   % fixed seed, so this test is reproducible
noiseLevel = 0;
%noiseLevel = 0.01 * std(C_true*x);
y_clean = C_true*x;
y_noisy = y_clean + noiseLevel*randn(size(y_clean));

i = 40;

%% Helper: run full sweep + classify, return persistence info per candidate frequency
function report = runSweepAndClassify(y, i, dt, nMax, freqTol, dampTol, macTol)
    [Yp_ref, Yf, ~] = buildHankelMatrix(y, i);
    P = hankelProjection(Yp_ref, Yf);
    results = sweepModelOrders(P, 1, dt, nMax);

    matches = cell(nMax,1);
    for n = 2:nMax
        matches{n} = classifyPoleStability(results(n-1), results(n), freqTol, dampTol, macTol);
    end
    report.results = results;
    report.matches = matches;
end

%% Test A: correct-order-and-above regime, noisy data, sweep up to nMax=14
fprintf('=== Test A: sweep with true order=6, 3 known modes, noisy data ===\n');
nMax_A = 14;
repA = runSweepAndClassify(y_noisy, i, dt, nMax_A, freqTol, dampTol, macTol);

fprintf('Order : frequencies found (rad/s)\n');
for n = 1:nMax_A
    fprintf('  n=%2d: %s\n', n, mat2str(sort(repA.results(n).frequency), 4));
end

fprintf('\nClass-3 links found at each order transition:\n');
for n = 2:nMax_A
    classes = [repA.matches{n}.class];
    fprintf('  n=%2d->%2d: classes = %s\n', n-1, n, mat2str(classes));
end

%% Check: does a class-3 chain exist near each true frequency, persisting
%  for at least 3 consecutive order transitions, somewhere in the sweep?
fprintf('\n=== Checking for persistent (>=3 consecutive) class-3 chains near true frequencies ===\n');
for q = 1:3
    found = false;
    for nStart = 2:(nMax_A-2)
        chainOK = true;
        for n = nStart:(nStart+2)
            classes = [repA.matches{n}.class];
            freqs = repA.results(n).frequency;
            hit = any(classes==3 & abs(freqs - omega_true(q))/omega_true(q) < freqTol*3);
            if ~hit, chainOK = false; break; end
        end
        if chainOK, found = true; break; end
    end
    fprintf('  f=%.4f rad/s: persistent class-3 chain found = %d\n', omega_true(q), found);
end

%% Test B: under-ordering (nMax deliberately below true order 6)
fprintf('\n=== Test B: under-ordering, nMax=4 (< true order 6) ===\n');
repB = runSweepAndClassify(y_noisy, i, dt, 4, freqTol, dampTol, macTol);
for n = 1:4
    fprintf('  n=%d: %d mode(s) found: %s\n', n, length(repB.results(n).frequency), ...
        mat2str(sort(repB.results(n).frequency), 4));
end
fprintf('(Expect: fewer than 3 distinct modes recoverable; system cannot represent all 3 true modes)\n');

%% Test C: over-ordering (nMax well above true order 6)
fprintf('\n=== Test C: over-ordering, nMax=20 (>> true order 6) ===\n');
nMax_C = 20;
repC = runSweepAndClassify(y_noisy, i, dt, nMax_C, freqTol, dampTol, macTol);
fprintf('Order : #modes found : #class-3 matches (of those)\n');
for n = 2:nMax_C
    nModes = length(repC.results(n).frequency);
    nClass3 = sum([repC.matches{n}.class]==3);
    fprintf('  n=%2d: %d modes, %d class-3\n', n, nModes, nClass3);
end
fprintf('(Expect: extra poles beyond the true 3 modes should NOT consistently form long stable chains)\n');