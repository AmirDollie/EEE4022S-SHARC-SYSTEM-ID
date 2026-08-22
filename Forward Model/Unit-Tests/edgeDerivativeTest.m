% Numerical test for the edgeDerivative.m

xi = 2.5832 + 1.3285i;  % complex, to stress-test that too
gamma = 0.1;
h = 1e-6;

%Checks central difference vs my closed form soln
fd = (verticalEigenfunction(-gamma+h, xi, gamma) - verticalEigenfunction(-gamma-h, xi, gamma)) / (2*h);
closedForm = edgeDerivative(xi, gamma);

disp(fd)
disp(closedForm)
a = abs(fd - closedForm);
fprintf('The difference between closed form and central diff: %E\n', a)   % should be tiny, e.g. < 1e-6
