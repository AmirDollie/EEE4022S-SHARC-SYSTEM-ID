function cv = analyticTransmissibilityCovariance(H, Su, Snn, ref, nEff, varargin)
%ANALYTICTRANSMISSIBILITYCOVARIANCE  First-order covariance of the corrected Welch transmissibility.
%   cv = analyticTransmissibilityCovariance(H, Su, Snn, ref, nEff)
%   cv = analyticTransmissibilityCovariance(Ssig, [], Snn, ref, nEff)   signal spectral matrix given
%   cv = analyticTransmissibilityCovariance(..., 'MinDenominatorFraction', 0.1, 'FiniteSample', true)
%
%   H     nS x nF complex acceleration FRF at the bins, layout sensors only (rows in layout order)
%         OR nS x nS x nF signal spectral matrix (e.g. welchSmoothedSpectrum, which includes the
%         window's resolution decoherence; use this form, the point form is optimistic at short L)
%   Su    1 x nF incident elevation spectrum at the bins (one-sided, per rad/s, tapered if the
%         record is tapered); [] when H is a spectral matrix
%   Snn   scalar or nS x 1 sensor noise PSD (one-sided, per rad/s): sigma^2 dt / pi
%   ref   reference row of H
%   nEff  effective number of Welch averages (extractTransmissibility's nEffective)
%
%   MODEL. Spectral matrix S = H H^H Su + diag(Snn) (or Ssig + diag(Snn)). The Welch estimate is
%   treated as complex Wishart with n = nEff - 1 degrees of freedom ('FiniteSample', true: a
%   finite-sample degrees-of-freedom correction, validated against direct complex Gaussian
%   simulation and the cached Monte Carlo; false uses nEff):
%     E[dS_ab conj(dS_cd)] = S_ac S_db / nEff,   E[dS_ab dS_cd] = S_ad S_cb / nEff.
%   The estimator T_j = S_jr / (S_rr - Snn_r) (the project's corrected estimator), linearised:
%     dT_j = (dS_jr - T_j dS_rr) / Den,   Den = S_rr - Snn_r,
%   gives the complex covariance C and pseudo-covariance P of the m transmissibilities:
%     C_jl = S_rr (S_jl - conj(T_l) S_jr - T_j S_rl + T_j conj(T_l) S_rr) / (nEff Den^2)
%     P_jl = (S_jr - T_j S_rr)(S_lr - T_l S_rr) / (nEff Den^2)   (nonzero only through Snn_r)
%   and the REAL covariance of [Re T_1..T_m, Im T_1..T_m] (stackTransmissibility order):
%     E[x x'] = Re(C + P)/2,  E[y y'] = Re(C - P)/2,  E[x y'] = Im(P - C)/2.
%   Shared-reference correlation between different T_j is kept in full (not diagonal).
%
%   OUTPUT cv: blocks (2m x 2m x nF), C, P (m x m x nF), T (m x nF), Srr, Den (1 x nF),
%   coherence (m x nF, |S_jr|^2 / (S_jj S_rr)), valid (1 x nF, Den > MinDenominatorFraction S_rr,
%   the estimator's rejection rule), others.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Engine/.

opt = struct('MinDenominatorFraction', 0.1, 'FiniteSample', true);
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i + 1}; end
isMatrix = ndims(H) == 3 || (size(H, 1) == size(H, 2) && isempty(Su));
if isMatrix
    nS = size(H, 1); nF = size(H, 3);
else
    [nS, nF] = size(H);
end
nDof = nEff - opt.FiniteSample;
if isscalar(Snn), Snn = Snn * ones(nS, 1); end
Snn = Snn(:);
others = setdiff(1:nS, ref); m = numel(others);
blocks = zeros(2 * m, 2 * m, nF);
Cs = complex(zeros(m, m, nF)); Ps = Cs;
T = complex(zeros(m, nF)); Srr = zeros(1, nF); Den = Srr; coh = zeros(m, nF);
for k = 1:nF
    if isMatrix
        S = H(:, :, k) + diag(Snn);
    else
        h = H(:, k);
        S = (h * h') * Su(k) + diag(Snn);
    end
    srr = real(S(ref, ref));
    den = srr - Snn(ref);
    sjr = S(others, ref);                         % m x 1, S_jr = E[Y_j conj(Y_r)]
    Sjl = S(others, others);
    t = sjr / den;
    Ck = srr * (Sjl - sjr * t' - t * sjr' + srr * (t * t')) / (nDof * den^2);
    Ck = (Ck + Ck') / 2;
    u = sjr - t * srr;
    Pk = (u * u.') / (nDof * den^2);
    Pk = (Pk + Pk.') / 2;
    xx = real(Ck + Pk) / 2; yy = real(Ck - Pk) / 2; xy = imag(Pk - Ck) / 2;
    B = [xx, xy; xy.', yy];
    blocks(:, :, k) = (B + B.') / 2;
    Cs(:, :, k) = Ck; Ps(:, :, k) = Pk;
    T(:, k) = t; Srr(k) = srr; Den(k) = den;
    coh(:, k) = abs(sjr).^2 ./ (real(diag(Sjl)) * srr);
end
cv = struct('blocks', blocks, 'C', Cs, 'P', Ps, 'T', T, 'Srr', Srr, 'Den', Den, 'coherence', coh, ...
    'valid', Den > opt.MinDenominatorFraction * Srr, 'others', others, 'nEff', nEff, 'nDof', nDof);
end