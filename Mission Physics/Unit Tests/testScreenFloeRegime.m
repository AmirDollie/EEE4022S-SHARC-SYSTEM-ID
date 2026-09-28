%% testScreenFloeRegime.m
% Unit tests for waveNumber.m and screenFloeRegime.m (Mission Physics).
% Checks the dispersion solvers, the scaling laws, the p0 round trip into the EMM's
% nondimensional parameters, and the rigidity sensitivity against a finite difference.
%
% Lives in Mission Physics/Unit-Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
nPass = 0; nFail = 0;
g = 9.81;

% 1  deep water is omega^2 / g exactly
w = linspace(0.3, 9, 40);
[nPass, nFail] = check(max(abs(waveNumber(w, Inf) - w.^2 / g)) == 0, 'deep water k = omega^2/g', nPass, nFail);

% 2  finite-depth residual
for H = [0.3 1.88 5 50]
    k = waveNumber(w, H);
    res = max(abs(g * k .* tanh(k * H) - w.^2) ./ w.^2);
    [nPass, nFail] = check(res < 1e-12, sprintf('finite-depth residual < 1e-12 at H = %g m (%.1e)', H, res), nPass, nFail);
end

% 3  large depth tends to deep water; shallow limit k ~ omega / sqrt(g H)
k = waveNumber(w, 1e4);
[nPass, nFail] = check(max(abs(k ./ (w.^2 / g) - 1)) < 1e-12, 'H = 1e4 m matches deep water', nPass, nFail);
Hs = 0.05; ws = 0.05;
[nPass, nFail] = check(abs(waveNumber(ws, Hs) / (ws / sqrt(g * Hs)) - 1) < 1e-3, 'shallow limit k = omega/sqrt(gH)', nPass, nFail);

% 4  scaling laws
ice = struct('R', 20, 'h', 1, 'E', 5e9, 'nu', 0.3);
T = [6 10 15];
a = screenFloeRegime(ice, T);
b = screenFloeRegime(setfield(ice, 'R', 40), T); %#ok<SFLD>
[nPass, nFail] = check(max(abs(b.kR ./ a.kR - 2)) < 1e-12 && abs(b.RoverLf / a.RoverLf - 2) < 1e-12 ...
    && max(abs(b.kLf - a.kLf)) < 1e-12, 'doubling R doubles kR and R/l_f, leaves k l_f', nPass, nFail);
c = screenFloeRegime(setfield(ice, 'h', 2), T); %#ok<SFLD>
[nPass, nFail] = check(abs(c.flexuralLength / a.flexuralLength - 2^(3/4)) < 1e-12, 'l_f ~ h^(3/4)', nPass, nFail);
e = screenFloeRegime(setfield(ice, 'E', 4 * ice.E), T); %#ok<SFLD>
[nPass, nFail] = check(abs(e.flexuralLength / a.flexuralLength - 4^(1/4)) < 1e-12, 'l_f ~ E^(1/4)', nPass, nFail);
ref = (5e9 * 1 / (12 * 0.91) / (1025 * g))^(1/4);
[nPass, nFail] = check(abs(a.flexuralLength / ref - 1) < 1e-12, 'l_f formula', nPass, nFail);

% 5  ice-covered root: residual, zero-plate limit, stiffening, sensitivity by finite difference
for H = [Inf 30]
    s = screenFloeRegime(ice, T, 'Depth', H);
    hw = H - s.draught;
    if isinf(hw), tf = ones(size(s.kappa)); else, tf = tanh(s.kappa * hw); end
    res = max(abs((s.D * s.kappa.^4 + 1025 * g - s.m * s.omega.^2) .* s.kappa .* tf - 1025 * s.omega.^2) ./ (1025 * s.omega.^2));
    [nPass, nFail] = check(res < 1e-10, sprintf('ice dispersion residual (H = %g) %.1e', H, res), nPass, nFail);
    % stiffness always lengthens the wave (dln kappa/dln D < 0); where flexure dominates (k l_f > 0.5)
    % kappa < k, while at long periods the plate's mass loading can make kappa slightly > k
    fl = s.kLf > 0.5;
    [nPass, nFail] = check(all(s.dlnKappa_dlnD < 0) && all(s.kappaRatio(fl) < 1), ...
        sprintf('stiffness lengthens the wave; kappa < k where k l_f > 0.5 (H = %g)', H), nPass, nFail);
    sp = screenFloeRegime(struct('R', 20, 'D', s.D * 1.001, 'm', s.m), T, 'Depth', H);
    sm = screenFloeRegime(struct('R', 20, 'D', s.D / 1.001, 'm', s.m), T, 'Depth', H);
    fd = (log(sp.kappa) - log(sm.kappa)) / (2 * log(1.001));
    [nPass, nFail] = check(max(abs(fd - s.dlnKappa_dlnD)) < 1e-6, ...
        sprintf('dln kappa/dln D matches finite difference (H = %g), max err %.1e', H, max(abs(fd - s.dlnKappa_dlnD))), nPass, nFail);
end
z = screenFloeRegime(struct('R', 20, 'D', 1e-12, 'm', 1e-9), T);
[nPass, nFail] = check(max(abs(z.kappaRatio - 1)) < 1e-9, 'vanishing plate gives kappa = k', nPass, nFail);

% 6  p0 round trip: (beta, gamma, R) -> (D, m, R) -> EMM scaling at H = 1.88 recovers p0
p0 = [4.6985e-5, 1.4548e-3, 0.3830]; H0 = 1.88; rho0 = 1000;
basin = struct('R', p0(3) * H0, 'D', p0(1) * rho0 * g * H0^4, 'm', p0(2) * H0 * rho0);
s = screenFloeRegime(basin, 2 * pi ./ [2.979 8.501], 'Depth', H0, 'RhoWater', rho0);   % U2 band edges
[nPass, nFail] = check(abs(s.emm.beta / p0(1) - 1) < 1e-12 && abs(s.emm.gamma / p0(2) - 1) < 1e-12 ...
    && abs(s.emm.Rnd / p0(3) - 1) < 1e-12, 'p0 round trip (beta, gamma, R)', nPass, nFail);
[nPass, nFail] = check(max(abs(s.emm.alpha - [1.7 13.85]) ./ [1.7 13.85]) < 2e-3 && all(s.emm.alphaInRange), ...
    sprintf('p0 band maps to alpha [%.3f %.3f]', s.emm.alpha), nPass, nFail);

% 7  deep-water effective depth puts min(k, kappa) H = pi at the longest period
s = screenFloeRegime(ice, T);
[nPass, nFail] = check(abs(min(s.k(end), s.kappa(end)) * s.emm.H - pi) < 1e-12, ...
    sprintf('effective depth %.2f m, alpha %.2f to %.2f', s.emm.H, s.emm.alpha([end 1])), nPass, nFail);

% 8  worked example from the plan: T = 10 s deep water, lambda about 156 m
s = screenFloeRegime(struct('R', 25, 'h', 1, 'E', 6e9, 'nu', 0.3), 10);
[nPass, nFail] = check(abs(s.lambda - 156.1) < 0.5 && abs(s.kR - 1) < 0.01 && abs(s.flexuralLength - 15.3) < 0.2, ...
    sprintf('plan example: lambda %.1f m, kR %.3f, l_f %.2f m', s.lambda, s.kR, s.flexuralLength), nPass, nFail);

% 9  parity with the production EMM dispersion code (Forward Model/dispersionRoots.m):
%    the travelling ice-covered root is xi0 = i kappa H in its scaling
fmDir = fullfile(thisDir, '..', '..', 'Forward Model');
if exist(fullfile(fmDir, 'dispersionRoots.m'), 'file')
    addpath(fmDir);
    alphas = [1.7 4 8 13.85];
    s = screenFloeRegime(basin, 2 * pi ./ sqrt(alphas * g / H0), 'Depth', H0, 'RhoWater', rho0);
    err = zeros(size(alphas));
    for i = 1:numel(alphas)
        xi = dispersionRoots(alphas(i), p0(1), p0(2), 3);
        err(i) = abs(imag(xi(3)) - s.kappa(i) * H0) / (s.kappa(i) * H0);
    end
    [nPass, nFail] = check(max(err) < 1e-8, sprintf('kappa H matches dispersionRoots travelling root (max rel err %.1e)', max(err)), ...
        nPass, nFail);
else
    fprintf('  SKIP  parity with Forward Model/dispersionRoots.m (not found at %s)\n', fmDir);
end

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testScreenFloeRegime:failed', '%d test(s) failed.', nFail); end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end