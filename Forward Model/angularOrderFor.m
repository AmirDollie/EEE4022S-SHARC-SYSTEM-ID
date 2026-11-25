function N = angularOrderFor(alphaMax, R, margin, Nmin)
%ANGULARORDERFOR  Angular truncation N for a floe of non-dimensional radius R.   (F1, 2026-11-25)
%   N = angularOrderFor(alphaMax, R)              margin 10, floor 10
%   N = angularOrderFor(alphaMax, R, margin, Nmin)
%
%   The incident plane wave's angular orders n are significant up to |n| ~ p0 R, where p0 is the
%   open-water propagating wavenumber (p0 tanh p0 = alpha). Beyond that the Fourier-Bessel content of
%   the wave decays rapidly, so
%       N = max(Nmin, ceil(p0(alphaMax) * R) + margin),
%   evaluated at the TOP of the band, where p0 R is largest.
%
%   Evidence (diag_F1_plate.m, 8.5 rad/s, error in T relative to N = 20):
%       p0 R =  5.3 (p0):   N = 10 -> 8e-6,   N = 15 -> 5e-10
%       p0 R = 10.6 (2R0):  N = 10 -> 2.2e-2, N = 15 -> 1.3e-4
%   so a margin of ~5 gives ~1e-5 to 1e-4 and a margin of 10 is comfortably converged. Check any
%   production choice against N + 5.
%
%   Use ONE N per floe (from the nominal R and the band's top frequency) for the nominal field AND
%   its perturbed-parameter fields, so finite-difference derivatives never straddle a change of N.
    if nargin < 3 || isempty(margin), margin = 10; end
    if nargin < 4 || isempty(Nmin), Nmin = 10; end
    k0 = dispersionRoots(alphaMax, 0, 0, 1);
    p0 = imag(k0(1));
    N = max(Nmin, ceil(p0 * R) + margin);
end
