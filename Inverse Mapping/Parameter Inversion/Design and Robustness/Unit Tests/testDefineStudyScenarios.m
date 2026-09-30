%% testDefineStudyScenarios.m
% Consistency checks on the pre-registered configuration (no EMM solves). Prints the A1 physical
% scenarios in the EMM's scaling so the table can go straight into the notes.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
cfg = defineStudyScenarios();
nPass = 0; nFail = 0;

% 1  p0 round trip through the physical conversion
b0 = cfg.A1.physical(strcmp({cfg.A1.physical.name}, 'basin_p0'));
[nPass, nFail] = check(max(abs(b0.p ./ cfg.p0.vec - 1)) < 1e-12, 'basin_p0 converts back to p0 exactly', nPass, nFail);

% 1b independent check of the physical conversion (toEMM special-cases basin_p0, so test another)
q = cfg.A1.physical(strcmp({cfg.A1.physical.name}, 'pp3mm_R072'));
E = 1.5e9; h = 0.003; nu = 0.40; rhoI = 910; Hb = 1.88; rhoW = 1000; gg = 9.81;
D = E * h^3 / (12 * (1 - nu^2));                     % 4.018 N m
pHand = [D / (rhoW * gg * Hb^4), rhoI * h / (rhoW * Hb), 0.72 / Hb];
ok = max(abs(q.p ./ pHand - 1)) < 1e-12 && abs(D - 4.018) < 1e-3 && abs(q.p(1) / cfg.p0.beta - 0.698) < 1e-3;
[nPass, nFail] = check(ok, sprintf('pp3mm_R072 by hand: D = %.3f N m, beta = %.3f beta0, gamma = %.3g, R = %.4f', ...
    D, pHand(1) / cfg.p0.beta, pHand(2), pHand(3)), nPass, nFail);

% 2  candidate grid and Level 2C
G = cfg.layout.gridFrac;
ok = size(G, 1) == 33 && isequal(G(cfg.layout.level2C.idx, :), cfg.layout.level2C.frac);
[nPass, nFail] = check(ok, '33-point grid; indices [2 6 14 26] are the Level 2C sensors (s2 s6 s14 s26)', nPass, nFail);
SENSORS = [0.3 * cfg.p0.R, 0; 0.3 * cfg.p0.R, pi; 0.5 * cfg.p0.R, pi; 0.9 * cfg.p0.R, 0];  % G3 script
ok = max(max(abs(cfg.layout.level2C.frac .* [cfg.p0.R 1] - SENSORS))) < 1e-15;
[nPass, nFail] = check(ok, 'Level 2C fractions reproduce the G3 SENSORS matrix at p0', nPass, nFail);

% 3  basin band inside the validated alpha range
alpha = cfg.physics.Hbasin * cfg.band.basin.^2 / cfg.physics.g;
ok = alpha(1) >= cfg.emm.alphaValidated(1) && alpha(2) <= cfg.emm.alphaValidated(2);
[nPass, nFail] = check(ok, sprintf('basin band alpha [%.3f %.3f] inside [%.2f %.2f]', alpha, cfg.emm.alphaValidated), ...
    nPass, nFail);

% 4  class thresholds nested
I = cfg.class.identifiable; M = cfg.class.marginal;
ok = I.maxSigma < M.maxSigma && I.maxBiasDistance < M.maxBiasDistance && I.maxCondF <= M.maxCondF;
[nPass, nFail] = check(ok, 'Identifiable thresholds strictly inside Marginal', nPass, nFail);

% 5  noise: G3 value reproduces the twin's noise std (NOISE_REL x reference signal std)
sg3 = cfg.noise.sigma(strcmp(cfg.noise.names, 'g3'));
[nPass, nFail] = check(abs(sg3 - 6.67e-3) < 0.01e-3, sprintf('G3 noise sigma %.4g m/s^2 (twin 6.67e-3)', sg3), ...
    nPass, nFail);

% 6  D2 segment rule
K = @(Tmin, L) floor(2 * round(Tmin * 60 / cfg.acq.dt) / L) - 1;
ok = K(10, 2048) < cfg.D2.minSegments && K(10, 512) >= cfg.D2.minSegments && K(109, 2048) >= cfg.D2.minSegments;
[nPass, nFail] = check(ok, sprintf('K>=8 rule: 10 min L=2048 K=%d rejected, L=512 K=%d kept, 109 min L=2048 K=%d kept', ...
    K(10, 2048), K(10, 512), K(109, 2048)), nPass, nFail);

% 7  A1 scenario set
P = cfg.A1.physical; fin = ~[P.deepWater];
pp = vertcat(P(fin).p);
ok = all(isfinite(pp(:))) && all(pp(:) > 0) && all(pp(:, 2) < 1) && all(isnan([P(~fin).p]));
nTot = numel(cfg.A1.param) + numel(P);
ok = ok && nTot >= 20 && nTot <= 25;
[nPass, nFail] = check(ok, sprintf('A1: %d scenarios (%d parameter-space, %d physical); finite-depth p valid, deep water deferred', ...
    nTot, numel(cfg.A1.param), numel(P)), nPass, nFail);

% 8  D1 physical admissibility of Level 2C at p0
xy = @(rt, R) [rt(:, 1) .* R .* cos(rt(:, 2)), rt(:, 1) .* R .* sin(rt(:, 2))];
X = xy(cfg.layout.level2C.frac, cfg.p0.Rphys);
dmin = min(pdist2local(X));
ok = dmin >= cfg.D1.minSpacing && max(cfg.layout.level2C.frac(:, 1)) * cfg.p0.Rphys <= cfg.p0.Rphys - cfg.D1.edgeClearance;
[nPass, nFail] = check(ok, sprintf('Level 2C admissible at p0 with a %.0f cm buoy (min spacing %.3f m)', ...
    100 * cfg.D1.buoyDiameter, dmin), nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);

% ---- informative: A1 physical scenarios in EMM scaling --------------------------------------
fprintf('\nA1 physical scenarios (finite depth in EMM scaling; deep water via screenFloeRegime later)\n');
fprintf('  %-22s %6s %7s %8s %10s %9s %9s %7s %8s\n', 'name', 'class', 'R (m)', 'h (mm)', 'beta/beta0', 'gamma', 'R/H', ...
    'R/l_f', 'prov.');
for i = 1:numel(P)
    s = P(i);
    fprintf('  %-22s %6s %7.3f %8.1f %10.3g %9.3g %9.4f %7.3f %8s\n', s.name, s.class, s.Rphys, 1e3 * s.h, ...
        s.p(1) / cfg.p0.beta, s.p(2), s.p(3), s.RoverLf, yn(s.provisional));
end
% D1: how many n_s = 4 layouts the buoy constraint removes at p0
Xg = xy(cfg.layout.gridFrac, cfg.p0.Rphys);
okR = cfg.layout.gridFrac(:, 1) <= cfg.D1.maxRadiusFrac;
C = nchoosek(find(okR).', 4); nOk = 0;
for i = 1:size(C, 1), nOk = nOk + (min(pdist2local(Xg(C(i, :), :))) >= cfg.D1.minSpacing); end
fprintf('\nD1 at p0: %d of %d four-sensor layouts admissible with a %.0f cm buoy\n', nOk, size(C, 1), 100 * cfg.D1.buoyDiameter);

if nFail > 0, error('testDefineStudyScenarios:failed', '%d check(s) failed.', nFail); end

function d = pdist2local(X)
n = size(X, 1); d = zeros(1, n * (n - 1) / 2); c = 0;
for i = 1:n - 1
    for j = i + 1:n
        c = c + 1; d(c) = norm(X(i, :) - X(j, :));
    end
end
end

function s = yn(tf)
if tf, s = 'yes'; else, s = 'no'; end
end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end