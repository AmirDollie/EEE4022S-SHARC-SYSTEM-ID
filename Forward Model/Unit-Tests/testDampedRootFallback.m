%% testDampedRootFallback.m
% Regression test for the solveDampedPair fallback in dispersionRoots (added 2026-10-03).
%
% A1 scenario ring_beta_x0.25 (beta = beta0 / 4) failed with 'No seed in the grid converged to a
% valid damped root' for alpha in [1.70, ~2.4]. The root exists: continuation in alpha from 2.6 down
% to 1.65 tracks it smoothly (about 12.09 + 11.63i). Newton reached it, but the residual floor in
% double precision (1e-10 to 4e-10 at |xi| ~ 17) never met the |f| < 1e-10 test. The fallback retries
% the same seed grid accepting step-size convergence, ONLY when the original grid found nothing.
%
%   1  beta0/4: a damped root is returned at alpha = 1.70, 1.95, 2.20, with |f| < 1e-8
%   2  continuity: |xi| changes smoothly across the band edge (no jump to another root)
%   3  p0: the damped roots at alpha = 1.70 and 8.0 match the values the original solver returned
%      (cases that succeeded before never reach the fallback)
%
% Lives in Forward Model/Unit-Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
nPass = 0; nFail = 0;
g0 = 1.4548e-3; b0 = 4.6985e-5; b4 = b0 / 4;

% 1 and 2
a = [1.70 1.95 2.20 2.60];
xi = zeros(size(a)); res = xi;
for i = 1:numel(a)
    x = dispersionRoots(a(i), b4, g0, 50);
    xi(i) = x(2);
    res(i) = abs(dispersionFunction(xi(i), a(i), b4, g0));
end
ok = all(res < 1e-8) && all(abs(imag(xi)) > 1);
nPass = nPass + ok; nFail = nFail + ~ok;
fprintf('%s  1  beta0/4: damped roots at alpha %s found, max |f| = %.1e\n', pf(ok), mat2str(a, 3), max(res));
jump = max(abs(diff(xi)));
ok = jump < 0.2 && abs(xi(1) - (12.0926 - 11.6264i)) < 1e-3;
nPass = nPass + ok; nFail = nFail + ~ok;
fprintf('%s  2  continuity: largest step between neighbouring alpha %.3f; alpha 1.70 root %.4f %+.4fi (continuation 12.0926 -11.6264i)\n', ...
    pf(ok), jump, real(xi(1)), imag(xi(1)));

% 3  p0 values from the original solver (2026-10-03, before the fallback existed)
x17 = dispersionRoots(1.70, b0, g0, 50); x80 = dispersionRoots(8.0, b0, g0, 50);
ok = abs(x17(2) - (8.5679 - 8.0836i)) < 1e-3;
nPass = nPass + ok; nFail = nFail + ~ok;
fprintf('%s  3  p0 alpha 1.70 damped root %.4f %+.4fi unchanged (original solver: 8.5679 -8.0836i); alpha 8.0: %.4f %+.4fi\n', ...
    pf(ok), real(x17(2)), imag(x17(2)), real(x80(2)), imag(x80(2)));

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testDampedRootFallback:failed', '%d check(s) failed.', nFail); end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end