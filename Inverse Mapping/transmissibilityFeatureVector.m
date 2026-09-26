function [f, out] = transmissibilityFeatureVector(theta, omega, sensors, varargin)
%TRANSMISSIBILITYFEATUREVECTOR Model-side transmissibility features from the
%EMM: (theta, omega) -> H(omega; theta) -> T = H_j / H_r -> stacked f_T.
%
%   [F, OUT] = TRANSMISSIBILITYFEATUREVECTOR(THETA, OMEGA, SENSORS)
%   [F, OUT] = TRANSMISSIBILITYFEATUREVECTOR(..., Name, Value)
%
%   THETA   - LOG of the free parameters, in the order of 'ThetaIdx'
%             (default [1 3]: theta = [ln beta; ln R], gamma fixed).
%   OMEGA   - frequencies (rad/s) at which the features are evaluated.
%             Pass est.omega from extractTransmissibility so the model is
%             evaluated exactly where the data-side estimate exists.
%   SENSORS - nS x 2 [r theta] (non-dimensional r), as computeSensorFRF.
%
%   F       - stacked real feature vector, stackTransmissibility layout
%             (frequency-major: per frequency [Re T_1..T_m; Im T_1..T_m]).
%   OUT     - struct: p (full [beta gamma R] used), theta, thetaIdx, omega,
%             reference, others, T (m x nF), H (nS x nF acceleration FRF),
%             frfInfo (from computeSensorFRF).
%
%   PARAMETERS: p = PBase; p(ThetaIdx) = exp(theta). The fixed parameters
%   (gamma, by default) come from 'PBase'. ThetaIdx = 1:3 gives the
%   three-parameter map, so nothing here assumes gamma is fixed.
%
%   T is a ratio of two acceleration FRFs, so it is identical for the
%   displacement FRFs (the -omega^2 cancels); acceleration is used for
%   consistency with the rest of the pipeline.
%
%   NAME-VALUE OPTIONS
%     'ThetaIdx'    indices into [beta gamma R] (default [1 3])
%     'PBase'       full [beta gamma R] supplying the fixed parameters
%                   (default p0 = [4.6985e-5 1.4548e-3 0.3830])
%     'Reference'   reference sensor r (default 4, s26); must match the
%                   'Reference' used in extractTransmissibility
%     'FRFOptions'  cell of extra name-value pairs for computeSensorFRF
%                   (e.g. {'Truncation', [60 12 12]}); default {}
%     'MinRefFraction'  error if |H_r| < this fraction of max_j |H_j| at any
%                   frequency (default 1e-6): a near-zero reference makes T
%                   meaningless. Must be a finite scalar in [0, 1).
%
%   ERROR IDENTIFIERS (transmissibilityFeatureVector:<id>)
%     badTheta, badOmega, badSensors, badOption, referenceNearZero,
%     nonFiniteFeature
%   computeSensorFRF errors (including sensorOutsideFloe when R shrinks
%   below a sensor radius) propagate unchanged.
%
%   Lives in Inverse Mapping/. Needs Transmissibility/ (stackTransmissibility)
%   on the path. Tester:
%   Inverse Mapping/Inverse Mapping Unit Tests/transmissibilityFeatureVectorTester.m

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Transmissibility'));

%% ---- Options ---------------------------------------------------------------
opt = struct('ThetaIdx', [1 3], 'PBase', [4.6985e-5, 1.4548e-3, 0.3830], 'Reference', 4, ...
    'FRFOptions', {{}}, 'MinRefFraction', 1e-6);
if mod(numel(varargin), 2) ~= 0
    error('transmissibilityFeatureVector:badOption', 'Options must be name-value pairs.');
end
names = fieldnames(opt);
for a = 1:2:numel(varargin)
    key = varargin{a};
    if isa(key, 'string'), key = char(key); end
    hit = ischar(key) & strcmpi(key, names);
    if ~any(hit)
        error('transmissibilityFeatureVector:badOption', 'Unknown option.');
    end
    opt.(names{hit}) = varargin{a + 1};
end
idx = opt.ThetaIdx;
if ~(isnumeric(idx) && isvector(idx) && all(ismember(idx, 1:3)) && numel(unique(idx)) == numel(idx))
    error('transmissibilityFeatureVector:badOption', 'ThetaIdx must be distinct indices in 1..3.');
end
pb = opt.PBase;
if ~(isnumeric(pb) && isreal(pb) && numel(pb) == 3 && all(isfinite(pb)) && all(pb > 0))
    error('transmissibilityFeatureVector:badOption', 'PBase must be three positive finite values [beta gamma R].');
end
if ~iscell(opt.FRFOptions)
    error('transmissibilityFeatureVector:badOption', 'FRFOptions must be a cell of name-value pairs.');
end
mrf = opt.MinRefFraction;
if ~(isnumeric(mrf) && isscalar(mrf) && isreal(mrf) && isfinite(mrf) && mrf >= 0 && mrf < 1)
    error('transmissibilityFeatureVector:badOption', 'MinRefFraction must be a finite scalar in [0, 1).');
end

%% ---- Inputs ---------------------------------------------------------------------
if ~(isnumeric(theta) && isreal(theta) && numel(theta) == numel(idx) && all(isfinite(theta)))
    error('transmissibilityFeatureVector:badTheta', ...
        'THETA must hold %d finite real log-parameters (one per ThetaIdx entry).', numel(idx));
end
if ~(isnumeric(omega) && isreal(omega) && isvector(omega) && ~isempty(omega) && all(isfinite(omega)) && all(omega > 0))
    error('transmissibilityFeatureVector:badOmega', 'OMEGA must be a nonempty vector of positive frequencies.');
end
omega = omega(:).';
if ~(isnumeric(sensors) && ismatrix(sensors) && size(sensors, 2) == 2 && size(sensors, 1) >= 2)
    error('transmissibilityFeatureVector:badSensors', 'SENSORS must be nS x 2 with nS >= 2.');
end
nS = size(sensors, 1);
r = opt.Reference;
if ~(isnumeric(r) && isscalar(r) && r == round(r) && r >= 1 && r <= nS)
    error('transmissibilityFeatureVector:badOption', 'Reference must be a sensor index in 1..%d.', nS);
end

p = pb(:).';
p(idx) = exp(theta(:).');
if ~all(isfinite(p))
    error('transmissibilityFeatureVector:badTheta', 'Exponentiating THETA produced non-finite parameters.');
end

%% ---- EMM -> H -> T -> f ---------------------------------------------------------------
[H, frfInfo] = computeSensorFRF(omega, sensors, p, 'Output', 'acceleration', opt.FRFOptions{:});
others = setdiff(1:nS, r);
if any(abs(H(r, :)) < mrf * max(abs(H), [], 1))
    bad = omega(abs(H(r, :)) < mrf * max(abs(H), [], 1));
    error('transmissibilityFeatureVector:referenceNearZero', ...
        'Reference sensor response is near zero at omega = %s rad/s; choose another reference.', mat2str(bad, 5));
end
T = H(others, :) ./ H(r, :);
f = stackTransmissibility(T);
if ~all(isfinite(f))
    error('transmissibilityFeatureVector:nonFiniteFeature', 'Non-finite transmissibility feature.');
end

if nargout > 1
    out = struct('p', p, 'theta', theta(:).', 'thetaIdx', idx, 'omega', omega, 'reference', r, ...
        'others', others, 'T', T, 'H', H, 'frfInfo', frfInfo, 'layout', ...
        'frequency-major: per frequency [Re T_1..T_m; Im T_1..T_m]');
end
end