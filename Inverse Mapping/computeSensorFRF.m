function [H, info] = computeSensorFRF(omega, sensors, p, varargin)
%COMPUTESENSORFRF EMM frequency response at fixed sensor locations.
%
%   [H, INFO] = COMPUTESENSORFRF(OMEGA, SENSORS, P)
%   [H, INFO] = COMPUTESENSORFRF(..., Name, Value)
%
%   OMEGA   - 1 x nW (or nW x 1) PHYSICAL angular frequencies (rad/s), > 0.
%   SENSORS - nS x 2, each row [r, theta]: NON-DIMENSIONAL radius (r/H,
%             same convention as evaluateDeflection, 0 <= r <= R) and
%             angle (rad; theta = 0 is the incidence axis).
%   P       - [beta, gamma, R] (non-dimensional, as in the Forward Model).
%
%   H       - nS x nW complex frequency response, row order = SENSORS order.
%             'displacement': H_eta = eta, the plate deflection per unit
%                 incident wave amplitude (dimensionless RAO: metres of
%                 deflection per metre of incident amplitude).
%             'acceleration': H_a = -omega.^2 .* H_eta, in (m/s^2) per metre
%                 of incident amplitude (omega physical).
%   INFO    - struct: omega, alpha, sensors, p, output, truncation, depth,
%             gravity, nu, Hdisplacement (always), validatedAlphaRange,
%             outsideValidated (1 x nW logical), runtimeSec (1 x nW).
%
%   METHOD (no new physics): for each frequency,
%       alpha = depth * omega^2 / g,
%       data  = precomputeDeflectionData(alpha, beta, gamma, R, nu, M, P, N),
%       eta   = evaluateDeflection(data, r, theta)   (all sensors in one call).
%   One EMM solve per frequency; the per-sensor evaluation is cheap.
%
%   TIME CONVENTION: the Forward Model follows Montiel (2.9)-(2.10),
%   zeta(t) = Re{eta * exp(+i*omega*t)}, the same convention as
%   evaluateSpectralDeflection.m and stochasticSensorSynthesis.m, so H can
%   be passed to stochasticSensorSynthesis directly. For a passive system
%   under exp(+i w t) a local pole fit must return zeta > 0; that is the
%   end-to-end check of this convention (fitLocalPoleModel.m).
%
%   VALIDATED RANGE: truncation [M P N] = [50 10 10] was validated for
%   alpha in [1.7, 13.85] (about 2.9 to 8.5 rad/s at depth 1.88 m).
%   Frequencies outside 'ValidatedAlphaRange' are computed but flagged in
%   info.outsideValidated, and ONE warning
%   (computeSensorFRF:outsideValidatedAlpha) lists them. If you raise the
%   truncation to check convergence above the ceiling, the flag still
%   fires: convergence is only established by comparing truncations.
%
%   NAME-VALUE OPTIONS
%   'Output'              'acceleration' (default) or 'displacement'
%   'Truncation'          [M P N], default [50 10 10]
%   'WaterDepth'          H in metres, default 1.88
%   'Gravity'             default 9.81
%   'Nu'                  Poisson ratio, default 0.3
%   'ValidatedAlphaRange' [aMin aMax], default [1.7 13.85]
%   'UseParallel'         false (default) or true (parfor over frequencies;
%                         falls back to serial without the toolbox)
%   'Verbose'             false (default) or true (per-frequency progress)
%
%   ERROR IDENTIFIERS (computeSensorFRF:<id>)
%     badOmega, badSensors, sensorOutsideFloe, badP, badOption,
%     nonFiniteResponse
%   Errors raised inside the Forward Model propagate unchanged.
%
%   Lives in Inverse Mapping/. Tester:
%   Inverse Mapping/Inverse Mapping Unit Tests/computeSensorFRFTester.m

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', 'Forward Model', 'Animation'));

%% ---- Inputs --------------------------------------------------------------
if ~(isnumeric(omega) && isreal(omega) && isvector(omega) && ~isempty(omega) ...
        && all(isfinite(omega)) && all(omega > 0))
    error('computeSensorFRF:badOmega', 'OMEGA must be a nonempty real vector of positive finite frequencies (rad/s).');
end
omega = omega(:).';
nW = numel(omega);

if ~(isnumeric(p) && isreal(p) && numel(p) == 3 && all(isfinite(p)) && p(1) > 0 && p(2) >= 0 && p(3) > 0)
    error('computeSensorFRF:badP', 'P must be [beta gamma R] with beta > 0, gamma >= 0, R > 0.');
end
p = p(:).';
beta = p(1); gamma = p(2); R = p(3);

if ~(isnumeric(sensors) && isreal(sensors) && ismatrix(sensors) && size(sensors, 2) == 2 ...
        && size(sensors, 1) >= 1 && all(isfinite(sensors(:))) && all(sensors(:, 1) >= 0))
    error('computeSensorFRF:badSensors', 'SENSORS must be an nS x 2 real matrix of [r theta] rows with r >= 0.');
end
if any(sensors(:, 1) > R * (1 + 1e-12))
    error('computeSensorFRF:sensorOutsideFloe', ...
        'Sensor radius %.6g exceeds the floe radius R = %.6g (evaluation outside the plate is not meaningful).', ...
        max(sensors(:, 1)), R);
end
nS = size(sensors, 1);

opt = struct('Output', 'acceleration', 'Truncation', [50 10 10], 'WaterDepth', 1.88, 'Gravity', 9.81, ...
    'Nu', 0.3, 'ValidatedAlphaRange', [1.7 13.85], 'UseParallel', false, 'Verbose', false);
if mod(numel(varargin), 2) ~= 0
    error('computeSensorFRF:badOption', 'Options must be name-value pairs.');
end
names = fieldnames(opt);
for a = 1:2:numel(varargin)
    key = varargin{a};
    if isa(key, 'string'), key = char(key); end
    if ~ischar(key)
        error('computeSensorFRF:badOption', 'Option names must be text.');
    end
    hit = strcmpi(key, names);
    if ~any(hit)
        error('computeSensorFRF:badOption', 'Unknown option ''%s''.', key);
    end
    opt.(names{hit}) = varargin{a + 1};
end

out = opt.Output;
if isa(out, 'string'), out = char(out); end
if ~(ischar(out) && any(strcmpi(out, {'acceleration', 'displacement'})))
    error('computeSensorFRF:badOption', 'Output must be ''acceleration'' or ''displacement''.');
end
out = lower(out);
tr = opt.Truncation;
if ~(isnumeric(tr) && numel(tr) == 3 && all(isfinite(tr)) && all(tr >= 1) && all(tr == round(tr)))
    error('computeSensorFRF:badOption', 'Truncation must be [M P N], positive integers.');
end
tr = tr(:).';
scalarPos = @(v) isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v > 0;
if ~scalarPos(opt.WaterDepth) || ~scalarPos(opt.Gravity)
    error('computeSensorFRF:badOption', 'WaterDepth and Gravity must be positive finite scalars.');
end
if ~(isnumeric(opt.Nu) && isscalar(opt.Nu) && isfinite(opt.Nu) && opt.Nu > -1 && opt.Nu < 0.5)
    error('computeSensorFRF:badOption', 'Nu must be a scalar in (-1, 0.5).');
end
vr = opt.ValidatedAlphaRange;
if ~(isnumeric(vr) && numel(vr) == 2 && all(isfinite(vr)) && vr(1) >= 0 && vr(2) > vr(1))
    error('computeSensorFRF:badOption', 'ValidatedAlphaRange must be [aMin aMax] with 0 <= aMin < aMax.');
end
vr = vr(:).';
if ~(islogical(opt.UseParallel) || isnumeric(opt.UseParallel)) || ~isscalar(opt.UseParallel)
    error('computeSensorFRF:badOption', 'UseParallel must be true or false.');
end
if ~(islogical(opt.Verbose) || isnumeric(opt.Verbose)) || ~isscalar(opt.Verbose)
    error('computeSensorFRF:badOption', 'Verbose must be true or false.');
end

%% ---- Frequency loop ------------------------------------------------------------
depth = opt.WaterDepth; g = opt.Gravity; nu = opt.Nu;
alpha = depth * omega.^2 / g;
outsideValidated = alpha < vr(1) | alpha > vr(2);

r = sensors(:, 1).'; th = sensors(:, 2).';
M = tr(1); Pt = tr(2); Nt = tr(3);
Heta = complex(zeros(nS, nW));
runtimeSec = zeros(1, nW);

% Parallel only in MATLAB with the Parallel Computing Toolbox. Workers must
% already have the Forward Model on their path (addpath BEFORE parpool, as
% in nominalFRFReconnaissance.m), otherwise they cannot find the EMM.
useParallel = logical(opt.UseParallel) && exist('OCTAVE_VERSION', 'builtin') == 0 && ~isempty(ver('parallel'));
if useParallel
    parfor k = 1:nW
        tk = tic;
        data = precomputeDeflectionData(alpha(k), beta, gamma, R, nu, M, Pt, Nt);
        Heta(:, k) = reshape(evaluateDeflection(data, r, th), [], 1);
        runtimeSec(k) = toc(tk);
    end
else
    for k = 1:nW
        tk = tic;
        data = precomputeDeflectionData(alpha(k), beta, gamma, R, nu, M, Pt, Nt);
        Heta(:, k) = reshape(evaluateDeflection(data, r, th), [], 1);
        runtimeSec(k) = toc(tk);
        if opt.Verbose
            fprintf('  computeSensorFRF: %3d/%3d  omega = %7.4f rad/s  alpha = %7.4f  %.2f s%s\n', ...
                k, nW, omega(k), alpha(k), runtimeSec(k), repmat(' [outside validated alpha]', 1, outsideValidated(k)));
        end
    end
end

if ~all(isfinite(Heta(:)))
    bad = find(any(~isfinite(Heta), 1));
    error('computeSensorFRF:nonFiniteResponse', ...
        'Non-finite response at omega = %s rad/s.', mat2str(omega(bad), 6));
end

switch out
    case 'displacement'
        H = Heta;
    case 'acceleration'
        H = -(omega.^2) .* Heta;
end

if any(outsideValidated)
    % Message built with sprintf first: MATLAB's warning() rejects non-scalar
    % format arguments, so the text is passed through a single '%s'.
    msg = sprintf(['%d of %d frequencies lie outside the validated alpha range [%.4g %.4g] (omega = %s rad/s). ' ...
        'Computed with truncation [%d %d %d]; check convergence before trusting them.'], ...
        sum(outsideValidated), nW, vr(1), vr(2), mat2str(omega(outsideValidated), 5), tr(1), tr(2), tr(3));
    warning('computeSensorFRF:outsideValidatedAlpha', '%s', msg);
end

if nargout > 1
    info = struct('omega', omega, 'alpha', alpha, 'sensors', sensors, 'p', p, 'output', out, ...
        'truncation', tr, 'depth', depth, 'gravity', g, 'nu', nu, 'Hdisplacement', Heta, ...
        'validatedAlphaRange', vr, 'outsideValidated', outsideValidated, 'runtimeSec', runtimeSec, ...
        'usedParallel', useParallel);
end
end