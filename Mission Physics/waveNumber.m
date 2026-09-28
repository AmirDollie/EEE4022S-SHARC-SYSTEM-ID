function k = waveNumber(omega, H, g)
%WAVENUMBER  Open-water wavenumber from omega^2 = g k tanh(k H).
%   k = waveNumber(omega, H)      H = Inf (or omitted) gives deep water, k = omega^2 / g.
%   k = waveNumber(omega, H, g)   gravity (default 9.81).
%   omega may be any array; H a positive scalar. Newton iteration from Eckart's
%   approximation, relative residual < 1e-13.
%
% Lives in Mission Physics/.

if nargin < 2 || isempty(H), H = Inf; end
if nargin < 3 || isempty(g), g = 9.81; end
kDeep = omega.^2 / g;
if isinf(H)
    k = kDeep;
    return
end
k = kDeep ./ sqrt(tanh(kDeep * H));                  % Eckart's approximation (within about 5%)
k(omega == 0) = 0;
for it = 1:50
    th = tanh(k * H);
    f = g * k .* th - omega.^2;
    fp = g * (th + k * H .* (1 - th.^2));
    dk = f ./ fp; dk(omega == 0) = 0;
    k = k - dk;
    if all(abs(dk(:)) <= 1e-15 * max(k(:), realmin)), break; end
end
end