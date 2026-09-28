function out = screenFloeRegime(floe, T, varargin)
%SCREENFLOEREGIME  Dimensional screening of a floe-wave case before any EMM run.
%   out = screenFloeRegime(floe, T, Name, Value, ...)
%
%   floe  struct with R (radius, m) and EITHER
%           h (m), E (Pa), nu, rhoIce (kg/m^3, default 917)       physical plate
%         OR
%           D (N m), m (kg/m^2)                                   rigidity and mass per area
%         (the second form is how the basin reference p0 enters: D = beta rho g H^4,
%          m = gamma H rho, see runFloeRegimeScreen.m)
%   T     wave periods (s), any vector.
%
%   Name-value options
%     'Depth'        water depth H (m); Inf = deep water (default Inf)
%     'RhoWater'     1025 (sea water; use 1000 for a fresh-water tank)
%     'Gravity'      9.81
%     'AlphaRange'   validated EMM range of alpha = H omega^2 / g, default [1.7 13.85]
%
%   Output struct (per-period quantities are row vectors the size of T)
%     inputs         R, h, E, nu, D, m, draught d = m / rhoWater, H, rhoWater, g
%     flexuralLength l_f = (D / (rhoWater g))^(1/4)
%     RoverLf        R / l_f                                  (independent of T)
%     omega, k, lambda, kH, kR, kLf                           open-water wave
%     kappa, kappaR  ice-covered travelling wavenumber, from
%                      (D kappa^4 + rhoWater g - m omega^2) kappa tanh(kappa (H - d)) = rhoWater omega^2
%     kappaRatio     kappa / k: how much the plate changes the wave (1 = none)
%     dlnKappa_dlnD  infinite-plate dispersion sensitivity to stiffness (0 = the rigidity does
%                    not change the flexural-gravity wavenumber); a cheap screen only, not the
%                    finite-floe transmissibility sensitivity dT/d ln beta
%     emm            the case in the EMM's scaling (lengths / H):
%                      H (tank) or an effective depth for deep water, chosen so that
%                      min(k, kappa) H = pi at the longest period (depth no longer matters),
%                      alpha, beta = (l_f / H)^4, gamma = d / H, Rnd = R / H (for deep water these
%                      depend on the chosen effective depth: numerical regime, not similarity),
%                      alphaInRange (per period), note
%
%   This is screening only: it does not run the EMM, compute information or decide
%   identifiability. It tells you which cases are worth sending to the EMM, and whether
%   they sit inside its validated alpha range.
%
% Lives in Mission Physics/.

opt = struct('Depth', Inf, 'RhoWater', 1025, 'Gravity', 9.81, 'AlphaRange', [1.7 13.85]);
if mod(numel(varargin), 2) ~= 0, error('screenFloeRegime:args', 'Name-value pairs expected.'); end
for i = 1:2:numel(varargin)
    name = varargin{i};
    if ~isfield(opt, name), error('screenFloeRegime:args', 'Unknown option ''%s''.', name); end
    opt.(name) = varargin{i + 1};
end
H = opt.Depth; rho = opt.RhoWater; g = opt.Gravity;

% ---- plate -------------------------------------------------------------------------------
if isfield(floe, 'D')
    D = floe.D; m = floe.m;
    h = NaN; E = NaN; nu = NaN;
    if isfield(floe, 'h'), h = floe.h; end
else
    h = floe.h; E = floe.E; nu = floe.nu;
    rhoIce = 917; if isfield(floe, 'rhoIce'), rhoIce = floe.rhoIce; end
    D = E * h^3 / (12 * (1 - nu^2));
    m = rhoIce * h;
end
R = floe.R;
d = m / rho;                                          % draught (floating equilibrium)
if ~isinf(H) && d >= H, error('screenFloeRegime:draught', 'Draught %.3g m exceeds the depth %.3g m.', d, H); end
Lf = (D / (rho * g))^(1/4);

% ---- open water --------------------------------------------------------------------------
T = T(:).';
omega = 2 * pi ./ T;
k = waveNumber(omega, H, g);
lambda = 2 * pi ./ k;

% ---- ice-covered travelling root and its rigidity sensitivity ------------------------------
kappa = zeros(size(omega)); sens = kappa;
for i = 1:numel(omega)
    [kappa(i), sens(i)] = iceRoot(omega(i), D, m, rho, g, H - d, k(i));
end

% ---- EMM scaling ---------------------------------------------------------------------------
if isinf(H)
    [~, iLong] = max(T);
    Hemm = pi / min(k(iLong), kappa(iLong));
    note = sprintf('deep water: effective depth %.4g m (min(k, kappa) H = pi at T = %.3g s)', Hemm, T(iLong));
else
    Hemm = H;
    note = sprintf('finite depth H = %.4g m', H);
end
alpha = Hemm * omega.^2 / g;
emm = struct('H', Hemm, 'alpha', alpha, 'beta', (Lf / Hemm)^4, 'gamma', d / Hemm, 'Rnd', R / Hemm, ...
    'alphaInRange', alpha >= opt.AlphaRange(1) & alpha <= opt.AlphaRange(2), 'alphaRange', opt.AlphaRange, ...
    'note', note);

out = struct('R', R, 'h', h, 'E', E, 'nu', nu, 'D', D, 'm', m, 'draught', d, 'H', H, 'rhoWater', rho, 'g', g, ...
    'flexuralLength', Lf, 'RoverLf', R / Lf, 'T', T, 'omega', omega, 'k', k, 'lambda', lambda, ...
    'kH', k * H, 'kR', k * R, 'kLf', k * Lf, 'kappa', kappa, 'kappaR', kappa * R, 'kappaRatio', kappa ./ k, ...
    'dlnKappa_dlnD', sens, 'emm', emm);
end

function [kap, sens] = iceRoot(w, D, m, rho, g, hw, k0)
% Travelling root of (D kap^4 + rho g - m w^2) kap tanh(kap hw) = rho w^2 (hw = Inf: deep).
% f is increasing in kap > 0 while m w^2 < rho g + D kap^4, which holds for any floating plate
% at these periods (checked below); safeguarded Newton on a bracket.
if m * w^2 >= rho * g
    error('screenFloeRegime:mass', 'm omega^2 >= rho g (T = %.3g s): plate inertia exceeds buoyancy.', 2 * pi / w);
end
if isinf(hw)
    th = @(x) ones(size(x)); dth = @(x) zeros(size(x));
else
    th = @(x) tanh(x * hw); dth = @(x) hw * (1 - tanh(x * hw).^2);
end
f = @(x) (D * x.^4 + rho * g - m * w^2) .* x .* th(x) - rho * w^2;
lo = 0; hi = max(k0, eps);
while f(hi) < 0, lo = hi; hi = 2 * hi; end
kap = min(max(k0, lo), hi);
for it = 1:200
    A = D * kap^4 + rho * g - m * w^2;
    t = th(kap);
    fp = (5 * D * kap^4 + rho * g - m * w^2) * t + A * kap * dth(kap);
    fv = A * kap * t - rho * w^2;
    if fv > 0, hi = kap; else, lo = kap; end
    kn = kap - fv / fp;
    if ~(kn > lo && kn < hi), kn = 0.5 * (lo + hi); end
    if abs(kn - kap) <= 1e-15 * kap, kap = kn; break; end
    kap = kn;
end
A = D * kap^4 + rho * g - m * w^2; t = th(kap);
fp = (5 * D * kap^4 + rho * g - m * w^2) * t + A * kap * dth(kap);
sens = -(D * kap^5 * t) / (kap * fp);                  % d ln kappa / d ln D (implicit function theorem)

end