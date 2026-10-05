function D = directionalSpectralMatrix(src, components, Snn, varargin)
%DIRECTIONALSPECTRALMATRIX  Sensor spectral matrix under several independent incident wave components.
%   Shared by V2 (directional forcing) and D3 (direction-robust layout), so the bias V2 reports is exactly the
%   bias D3 includes in its layout ranking.
%   D = directionalSpectralMatrix(src, components, Snn)
%   D = directionalSpectralMatrix(src, components, Snn, 'RequireUnitEnergy', false)
%
%   MODEL. y = sum_c H_c u_c + n, the components mutually independent and independent of the sensor noise:
%       S(w) = sum_c f_c H_c H_c^H S_u(w) + diag(S_nn)
%   S_u is the TOTAL incident spectrum and f_c the energy fraction of component c (sum f_c = 1 by default: the total
%   incident energy is fixed, so directionality is not confounded with SNR). One component gives the rank-one
%   matrix of the frozen inverse; two or more give rank > 1, and S_jr / S_rr is no longer H_j / H_r.
%   For a circular floe the response to a component from heading phi is the nominal (heading 0) field at the
%   sensor positions rotated by -phi:  H_phi(r, theta) = H_0(r, theta - phi). Each component is therefore a set of
%   ROWS of one field sampled at the rotated points (v2Tools.headingPoints builds them).
%
%   src, one of
%     Welch form (what the estimator sees; default): struct with fields
%        omega   retained bins (rad/s)
%        info    synthesiseTwinRecords INFO of a twin whose sensors include every rotated point
%        L       Welch segment length
%        S_sig = sum_c f_c welchSmoothedSpectrum(omega, info, L, rows_c): the engine's own kernel, so one
%        component at heading 0 reproduces evaluateScenario's smoothed spectral matrix exactly
%     Point form: struct with fields H (nP x nB complex FRF at the bins) and Su (1 x nB incident spectrum)
%        S_sig = sum_c f_c H(rows_c, :) H(rows_c, :)^H Su  (no window smoothing: tests and quick looks)
%   components  struct array with fields
%        rows      1 x nS, rows of info.H (or src.H) giving the layout sensors under that component, layout order
%        fraction  energy fraction f_c >= 0
%        heading   (optional) heading in rad, carried through for reporting
%   Snn         scalar or nS x 1 sensor noise PSD (one-sided, per rad/s), added on the diagonal of D.S
%
%   D.Ssig      nS x nS x nB signal spectral matrix (Hermitian, positive semidefinite)
%   D.S         D.Ssig + diag(Snn)
%   D.fractions, D.headings, D.form ('welch' | 'point'), D.nComponents
%
%   DOWNSTREAM (exactly the engine's path, no new estimator):
%     cv = analyticTransmissibilityCovariance(D.Ssig, [], Snn, ref, nEff)
%   gives the covariance of the corrected transmissibility under this sea (cv.blocks), its Welch expectation
%   cv.T = S_jr / (S_rr - S_nn,r) = Ssig_jr / Ssig_rr, the expected magnitude-squared coherence cv.coherence and
%   the estimator's bin acceptance cv.valid.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

opt = struct('RequireUnitEnergy', true);
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i + 1}; end
if isempty(components), error('directionalSpectralMatrix:none', 'At least one component is required.'); end
nC = numel(components);
fr = reshape([components.fraction], 1, []);
if any(~isfinite(fr)) || any(fr < 0)
    error('directionalSpectralMatrix:fraction', 'Energy fractions must be finite and non-negative.');
end
if opt.RequireUnitEnergy && abs(sum(fr) - 1) > 1e-12
    error('directionalSpectralMatrix:energy', 'Energy fractions sum to %.15g, not 1 (total energy is fixed).', sum(fr));
end
nS = numel(components(1).rows);
for c = 2:nC
    if numel(components(c).rows) ~= nS
        error('directionalSpectralMatrix:rows', 'Every component must map the same number of sensors.');
    end
end
hd = NaN(1, nC);
if isfield(components, 'heading')
    for c = 1:nC, if ~isempty(components(c).heading), hd(c) = components(c).heading; end, end
end

if isfield(src, 'info')
    form = 'welch';
    nB = numel(src.omega);
    Ssig = complex(zeros(nS, nS, nB));
    for c = 1:nC
        if fr(c) == 0, continue; end
        Ssig = Ssig + fr(c) * welchSmoothedSpectrum(src.omega, src.info, src.L, components(c).rows);
    end
elseif isfield(src, 'H') && isfield(src, 'Su')
    form = 'point';
    nB = size(src.H, 2);
    Ssig = complex(zeros(nS, nS, nB));
    su = reshape(src.Su, 1, []);
    for c = 1:nC
        if fr(c) == 0, continue; end
        Hc = src.H(components(c).rows, :);
        for k = 1:nB
            h = Hc(:, k);
            Ssig(:, :, k) = Ssig(:, :, k) + fr(c) * su(k) * (h * h');
        end
    end
    for k = 1:nB                                                 % exact Hermitian symmetry (the Welch form is
        Ssig(:, :, k) = (Ssig(:, :, k) + Ssig(:, :, k)') / 2;    % Hermitian by construction in welchSmoothedSpectrum)
    end
else
    error('directionalSpectralMatrix:src', 'src needs fields omega, info, L (Welch form) or H, Su (point form).');
end
if isscalar(Snn), Snn = Snn * ones(nS, 1); end
S = Ssig + repmat(diag(Snn(:)), [1 1 nB]);
D = struct('Ssig', Ssig, 'S', S, 'fractions', fr, 'headings', hd, 'form', form, 'nComponents', nC);
end