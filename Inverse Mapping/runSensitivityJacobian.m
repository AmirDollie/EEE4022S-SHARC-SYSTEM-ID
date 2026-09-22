%% runSensitivityJacobian.m
clear; clc;

% Reference parameter set
beta  = 4.6985e-5;
gamma = 1.4548e-3;
R     = 0.3830;

% Relative finite-difference step
epsilon = 1e-3;

fprintf('Running sensitivity Jacobian at:\n');
fprintf('  beta    = %.6e\n', beta);
fprintf('  gamma   = %.6e\n', gamma);
fprintf('  R       = %.6f\n', R);
fprintf('  epsilon = %.1e\n\n', epsilon);

tic
[J, info] = computeSensitivityJacobian(beta, gamma, R, epsilon);
elapsed = toc;

fprintf('\nCompleted in %.1f s\n\n', elapsed);

disp('J =')
disp(J)

disp('deltaP =')
disp(info.deltaP)

disp('fPlus - fMinus =')
disp(info.fPlus - info.fMinus)
J1 = J; info1 = info;

save('Jacobian_eps1e-3.mat', 'J1', 'info1');


%% Save the first one:
beta0  = 4.6985e-5;
gamma0 = 1.4548e-3;
R0     = 0.3830;

tic
[J1, info1] = computeSensitivityJacobian(beta0, gamma0, R0, 1e-3);
toc

save('Jacobian_eps1e-3.mat', 'J1', 'info1');
%% Next steps: Remaining epsilon sweep + convergence check
beta0 = 4.6985e-5; gamma0 = 1.4548e-3; R0 = 0.3830;

tic; [J2, info2] = computeSensitivityJacobian(beta0, gamma0, R0, 1e-4); toc
save('Jacobian_eps1e-4.mat', 'J2', 'info2');

tic; [J3, info3] = computeSensitivityJacobian(beta0, gamma0, R0, 1e-5); toc
save('Jacobian_eps1e-5.mat', 'J3', 'info3');

fprintf('\nColumn-wise relative change between step sizes:\n');
paramNames = {'beta', 'gamma', 'R'};
for j = 1:3
    e12 = norm(J2(:,j) - J1(:,j)) / norm(J2(:,j));
    e23 = norm(J3(:,j) - J2(:,j)) / norm(J3(:,j));
    fprintf('  %-6s  1e-3->1e-4 = %.3e   1e-4->1e-5 = %.3e\n', paramNames{j}, e12, e23);
end