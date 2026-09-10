%% buildModalBranchesTester.m
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

%% Test 1: hand-constructed 3-order case with a KNOWN branch structure
results(1).frequency = 2.0;  results(1).damping = 0.02;
results(2).frequency = 2.01; results(2).damping = 0.021;
results(3).frequency = [2.02, 5.0]; results(3).damping = [0.022, 0.03];

matches = cell(3,1);
matches{2} = struct('previousIndex', 1, 'frequencyStable', true, ...
                     'dampingStable', true, 'macStable', true, 'class', 3);
matches{3} = [ ...
    struct('previousIndex', 1, 'frequencyStable', true, 'dampingStable', true, 'macStable', true, 'class', 3), ...
    struct('previousIndex', 0, 'frequencyStable', false, 'dampingStable', false, 'macStable', false, 'class', 0) ...
];

branches = buildModalBranches(results, matches);

fprintf('=== Test 1: known branch structure ===\n');
assert(length(branches) == 2, 'Test 1 failed: expected exactly 2 branches.');
assert(isequal(branches(1).order, [1 2 3]), 'Test 1 failed: branch 1 orders incorrect.');
assert(isequal(branches(2).order, 3), 'Test 1 failed: branch 2 should only appear at order 3.');
assert(isequaln(branches(1).class, [NaN 3 3]), 'Test 1 failed: branch 1 classes incorrect.');
assert(max(abs(branches(1).frequency - [2.0 2.01 2.02])) < 1e-10, 'Test 1 failed: branch 1 frequencies incorrect.');
assert(abs(branches(2).frequency - 5.0) < 1e-10, 'Test 1 failed: branch 2 frequency incorrect.');
fprintf('PASS: 2 branches, orders/classes/frequencies all correct.\n');

%% Test 2: branches from the validated noise-free 3-mode sweep
fprintf('\n=== Test 2: branches from noise-free 3-mode sweep ===\n');
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
y = C_true*x;   % noise-free, matching the validated Pass 1 integration case

i = 40; nMax = 14;
[Yp_ref, Yf, ~] = buildHankelMatrix(y, i);
P = hankelProjection(Yp_ref, Yf);
r2 = sweepModelOrders(P, 1, dt, nMax);
m2 = cell(nMax,1);
for n = 2:nMax
    m2{n} = classifyPoleStability(r2(n-1), r2(n), freqTol, dampTol, macTol);
end
branches2 = buildModalBranches(r2, m2);

for q = 1:3
    found = false;
    for b = 1:length(branches2)
        if length(branches2(b).order) >= 8 && ...
           any(abs(branches2(b).frequency - omega_true(q)) / omega_true(q) < 0.01)
            found = true;
            break
        end
    end
    assert(found, 'Test 2 failed: no long (>=8 order) branch found near %.2f rad/s.', omega_true(q));
    fprintf('PASS: long branch found near %.4f rad/s (branch %d, spans orders [%d,%d])\n', ...
        omega_true(q), b, branches2(b).order(1), branches2(b).order(end));
end