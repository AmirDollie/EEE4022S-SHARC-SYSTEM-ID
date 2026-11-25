function pat = floeFlexPatterns(p, nu, omega, varargin)
%FLOEFLEXPATTERNS  Complex deflection and curvature patterns of one floe, per unit incident wave amplitude,
%   at each frequency in OMEGA, on a polar grid and at the sensor points. Used by buildFloeFlexCache and
%   floeFlexDemo.
%
%   pat = floeFlexPatterns([beta gamma R], nu, omega)
%   pat = floeFlexPatterns(..., 'Depth', 1.88, 'Truncation', [50 10 10], 'nr', 21, 'ntheta', 73, ...
%                          'SensorFrac', [0.3 0; 0.3 pi; 0.5 pi; 0.9 0])
%
%   PHYSICS. For each omega the EMM is solved by the existing precomputeDeflectionData (unchanged). The
%   deflection is then evaluated exactly as evaluateDeflection does, but vectorised over the grid, and
%   MULTIPLIED BY ALPHA: evaluateDeflection returns the response to an incident wave of amplitude 1/alpha
%   (solver normalisation, thesis Appendix D.2), so alpha * eta is the response per unit incident amplitude
%   (the displacement RAO of Eq. 3.14; the same correction as computeSensorFRF).
%
%   CURVATURE. Each term of the interior expansion, I_n(kappa_m r) e^{i n theta}, satisfies the modified
%   Helmholtz equation (Appendix D.2), so its Laplacian is exactly kappa_m^2 times itself. The Laplacian of
%   the deflection is therefore the same sum with an extra factor kappa_m^2: no finite differences. Rigid
%   heave and pitch (constant and linear in x) have zero Laplacian, so this field shows bending only.
%
%   OUTPUT (nondimensional, lengths / depth H)
%     pat.omega, pat.alpha     1 x nW
%     pat.k0                   1 x nW   open-water propagating root (incident wave e^{k0 x}, solver convention)
%     pat.r (1 x nr), pat.theta (ntheta x 1)
%     pat.eta                  ntheta x nr x nW   alpha * eta_EMM        (deflection per unit incident amplitude)
%     pat.lap                  ntheta x nr x nW   alpha * Lap eta_EMM    (nondimensional curvature, per unit amplitude)
%     pat.etaS                 nS x nW            at the sensor points
%     pat.sensorFrac, pat.p, pat.nu, pat.depth, pat.truncation, pat.seconds
%
%   Lives in Forward Model/Flex Demo/. Needs Forward Model/ and Forward Model/Animation/ on the path.

ip = inputParser;
ip.addParameter('Depth', 1.88);
ip.addParameter('Truncation', [50 10 10]);
ip.addParameter('nr', 21);
ip.addParameter('ntheta', 73);
ip.addParameter('SensorFrac', [0.3 0; 0.3 pi; 0.5 pi; 0.9 0]);    % Level 2C: s2, s6, s14, s26
ip.parse(varargin{:});
o = ip.Results;

beta = p(1); gamma = p(2); R = p(3);
M = o.Truncation(1); P = o.Truncation(2); N = o.Truncation(3);
g = 9.81;
omega = omega(:).';
nW = numel(omega);
r = linspace(0, R, o.nr);
theta = linspace(0, 2 * pi, o.ntheta).';
rs = o.SensorFrac(:, 1).' * R; ts = o.SensorFrac(:, 2).';
nS = numel(rs);

eta = complex(zeros(o.ntheta, o.nr, nW, 'single'));
lap = eta;
etaS = complex(zeros(nS, nW));
alpha = o.Depth * omega.^2 / g;
k0 = complex(zeros(1, nW));
t0 = tic;
for w = 1:nW
    a = alpha(w);
    data = precomputeDeflectionData(a, beta, gamma, R, nu, M, P, N);
    kOW = dispersionRoots(a, 0, 0, M);
    k0(w) = kOW(1);
    kap = data.kappa(:);
    E = zeros(o.ntheta, o.nr); L = E; Es = zeros(1, nS);
    for nIdx = 1:(2 * N + 1)
        n = nIdx - N - 1;
        c = data.dphiFull(:) .* data.AfullMat(:, nIdx);                 % nModes x 1
        den = besseli(n, kap * R);                                      % nModes x 1
        Ihat = besseli(n, kap * r) ./ den;                              % nModes x nr
        IhatS = besseli(n, kap * rs) ./ den;                            % nModes x nS
        f = c.' * Ihat;                                                 % 1 x nr
        fl = (c .* kap.^2).' * Ihat;                                    % 1 x nr  (Laplacian)
        ang = exp(1i * n * theta);
        E = E + ang * f;
        L = L + ang * fl;
        Es = Es + (c.' * IhatS) .* exp(1i * n * ts);
    end
    % evaluateDeflection divides by (1i * alpha); times alpha for unit incident amplitude
    eta(:, :, w) = single(E / 1i);
    lap(:, :, w) = single(L / 1i);
    etaS(:, w) = Es.' / 1i;
end
pat = struct('omega', omega, 'alpha', alpha, 'k0', k0, 'r', r, 'theta', theta, 'eta', eta, 'lap', lap, ...
    'etaS', etaS, 'sensorFrac', o.SensorFrac, 'p', p(:).', 'nu', nu, 'depth', o.Depth, ...
    'truncation', o.Truncation, 'seconds', toc(t0));
end