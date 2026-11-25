function testFloeFlexPatterns()
%TESTFLOEFLEXPATTERNS  Checks floeFlexPatterns before the cache is built. Run once; every line should say PASS.
%   1  the vectorised deflection equals alpha * evaluateDeflection (the existing, verified code)
%   2  long-wave limit: |alpha eta| -> 1 at the centre (thesis Eq. E.1), and its phase against the
%      incident wave as the demo draws it (-e^{k0 x}; see note)
%   3  the analytic Laplacian matches a finite-difference Laplacian of the deflection
%   4  the preset conversions reproduce the A1 table (R_phys / l_f)
%   Lives in Forward Model/Flex Demo/.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); addpath(fullfile(here, '..', 'Animation'));
p0 = [4.6985e-5, 1.4548e-3, 0.3830]; H = 1.88; g = 9.81;

% 1 against evaluateDeflection
w = 5.0;
pat = floeFlexPatterns(p0, 0.3, w, 'nr', 9, 'ntheta', 13);
data = precomputeDeflectionData(H * w^2 / g, p0(1), p0(2), p0(3), 0.3, 50, 10, 10);
[RR, TT] = meshgrid(pat.r, pat.theta);
ref = pat.alpha * evaluateDeflection(data, RR, TT);
e1 = max(abs(double(pat.eta(:, :, 1)) - ref), [], 'all') / max(abs(ref), [], 'all');
refS = pat.alpha * evaluateDeflection(data, pat.sensorFrac(:, 1) * p0(3), pat.sensorFrac(:, 2));
e1s = max(abs(pat.etaS(:, 1) - refS)) / max(abs(refS));
report(1, sprintf('grid vs alpha*evaluateDeflection: rel. error %.1e; sensors %.1e', e1, e1s), max(e1, e1s) < 1e-5);

% 2 long-wave limit (alpha = 0.05)
wL = sqrt(0.05 * g / H);
pl = floeFlexPatterns(p0, 0.3, wL, 'nr', 5, 'ntheta', 5);
c = double(pl.eta(1, 1, 1));                       % r = 0
report(2, sprintf('long-wave |alpha eta(0)| = %.4f (expect about 1); alpha eta(0) = %.4f %+.4fi', abs(c), ...
    real(c), imag(c)), abs(abs(c) - 1) < 0.01);
% the demo draws the incident elevation per unit amplitude as -e^{k0 x} (from d_z phi = i alpha eta with the
% solver's incident potential (1/(i alpha)) e^{k0 x}); at x = 0 that is -1, so a floe following the surface
% must give alpha eta(0) close to -1
report(2, sprintf('phase: alpha eta(0) / incident(0) = %.4f %+.4fi (expect about 1)', real(-c), imag(-c)), ...
    abs(-c - 1) < 0.02);

% 3 Laplacian against finite differences on a fine grid
pf = floeFlexPatterns(p0, 0.3, w, 'nr', 161, 'ntheta', 361);
E = double(pf.eta(:, :, 1)); Lan = double(pf.lap(:, :, 1));
dr = pf.r(2) - pf.r(1); dt = pf.theta(2) - pf.theta(1);
ir = 40:120; it = 2:size(E, 1) - 1;                % interior points away from the centre and the edge
Err = (E(it, ir + 1) - 2 * E(it, ir) + E(it, ir - 1)) / dr^2;
Er = (E(it, ir + 1) - E(it, ir - 1)) / (2 * dr);
Ett = (E(it + 1, ir) - 2 * E(it, ir) + E(it - 1, ir)) / dt^2;
rr = pf.r(ir);
Lfd = Err + Er ./ rr + Ett ./ rr.^2;
e3 = max(abs(Lfd - Lan(it, ir)), [], 'all') / max(abs(Lan(it, ir)), [], 'all');
report(3, sprintf('analytic vs finite-difference Laplacian: rel. difference %.1e', e3), e3 < 1e-2);

% 4 presets against the A1 table
tbl = {'1 m PP disc, 5 mm', 2.40; 'SHARC model floe', 0.22; 'MIZ floe, R = 10 m', 0.65; 'MIZ floe, R = 50 m', 3.27};
conv = @(E, h, Rp, nu, rhoI, rhoW, Hd) [E * h^3 / (12 * (1 - nu^2)) / (rhoW * 9.81 * Hd^4), Rp / Hd];
v = {conv(1.5e9, 0.005, 0.5, 0.4, 910, 1000, 1.88), conv(1.5e9, 0.04, 0.22, 0.4, 910, 1000, 1.88), ...
    conv(6e9, 1, 10, 0.3, 917, 1025, 122), conv(6e9, 1, 50, 0.3, 917, 1025, 122)};
for k = 1:4
    q = v{k}(2) / v{k}(1)^0.25;
    report(4, sprintf('%-20s R/l_f = %.3f (A1 table %.2f)', tbl{k, 1}, q, tbl{k, 2}), abs(q - tbl{k, 2}) < 0.02);
end
fprintf('Timing: one floe at 29 frequencies takes about %.0f s (from test 1).\n', 29 * pat.seconds);
end

function report(n, msg, ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
fprintf('  [%d] %s  %s\n', n, s, msg);
end