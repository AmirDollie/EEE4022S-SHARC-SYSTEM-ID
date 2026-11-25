%% precomputeDeflectionDataSplitTester.m   (F1, 2026-11-25)
% precomputeDeflectionDataSplit (open-water tower truncated separately) against precomputeDeflectionData.
%   Test 1  openWaterRoots(alpha, 500) vs dispersionRoots(alpha, 0, 0, 500)       pass: max rel diff <= 1e-14
%   Test 2  Mo = M = 50: AfullMat and sensor deflection vs the production path,
%           p0 / 2R0 / a thick test plate, three frequencies                        pass: max rel diff <= 1e-11
%   Test 3  Mo > M runs, stays finite, and the energy check on the far field is not available here, so
%           instead: deflection changes smoothly as Mo grows (Mo = 1000 vs 2000)      pass: rel change <= 1e-3
% Run from anywhere: run('Forward Model/Unit-Tests/precomputeDeflectionDataSplitTester.m')

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', 'Animation'));

depth = 1.88; g = 9.81; nu = 0.3; P = 10; N = 10;
floes = struct('name', {'p0', '2R0', 'thick'}, 'beta', {4.6985e-5, 4.6985e-5, 1e-2}, ...
    'gam', {1.4548e-3, 1.4548e-3, 0.1}, 'R', {0.3830, 0.7660, 1.0});
omega = [3.0 5.0 8.5]; alphaW = depth * omega.^2 / g;
sens = [0.3 0; 0.3 pi; 0.5 pi; 0.9 0];
allPass = true;

% Test 1
e1 = 0;
for w = 1:numel(alphaW)
    kp = dispersionRoots(alphaW(w), 0, 0, 500);
    kf = openWaterRoots(alphaW(w), 500);
    e1 = max(e1, max(abs(kf - kp) ./ abs(kp)));
end
ok = e1 <= 1e-14; allPass = allPass && ok;
fprintf('Test 1  openWaterRoots vs dispersionRoots (Mo = 500): max rel diff %.2e  %s\n', e1, pf(ok));

% Test 2
M = 50; e2A = 0; e2E = 0;
for f = 1:numel(floes)
    fl = floes(f);
    for w = 1:numel(alphaW)
        dP = precomputeDeflectionData(alphaW(w), fl.beta, fl.gam, fl.R, nu, M, P, N);
        dS = precomputeDeflectionDataSplit(alphaW(w), fl.beta, fl.gam, fl.R, nu, M, P, N, M);
        e2A = max(e2A, norm(dS.AfullMat - dP.AfullMat, 'fro') / norm(dP.AfullMat, 'fro'));
        eP = evaluateDeflection(dP, sens(:, 1).' * fl.R, sens(:, 2).');
        eS = evaluateDeflection(dS, sens(:, 1).' * fl.R, sens(:, 2).');
        e2E = max(e2E, max(abs(eS - eP) ./ abs(eP)));
    end
end
ok = max(e2A, e2E) <= 1e-11; allPass = allPass && ok;
fprintf('Test 2  Split (Mo = M = %d) vs production: AfullMat %.2e, sensor deflection %.2e  %s\n', M, e2A, e2E, pf(ok));

% Test 3
fl = floes(1); a = alphaW(3);
d1 = precomputeDeflectionDataSplit(a, fl.beta, fl.gam, fl.R, nu, M, P, N, 1000);
d2 = precomputeDeflectionDataSplit(a, fl.beta, fl.gam, fl.R, nu, M, P, N, 2000);
e1k = evaluateDeflection(d1, sens(:, 1).' * fl.R, sens(:, 2).');
e2k = evaluateDeflection(d2, sens(:, 1).' * fl.R, sens(:, 2).');
e3 = max(abs(e2k - e1k) ./ abs(e2k));
ok = all(isfinite(e2k)) && e3 <= 1e-3 && isequal(d2.truncation, [M P N 2000]); allPass = allPass && ok;
fprintf('Test 3  Mo = 1000 vs 2000 (p0, 8.5 rad/s, M = %d): rel change %.2e, finite %d  %s\n', M, e3, all(isfinite(e2k)), pf(ok));

fprintf('OVERALL: %s\n', pf(allPass));

function s = pf(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end
