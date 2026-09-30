function cal = calibrateAccelSixPosition(A, g, varargin)
%CALIBRATEACCELSIXPOSITION  Gravity-sphere calibration of one 3-axis accelerometer from static holds.
%   cal = calibrateAccelSixPosition(A, g)
%   cal = calibrateAccelSixPosition(A, g, 'Model', 'diag', 'Weights', w)
%
%   A     nH x 3 mean raw readings (m/s^2), one row per static hold, any orientations
%         (the name is kept from the plan; it takes N >= 7 holds, not six fixed faces).
%   g     local gravity (m/s^2), 9.796 in Cape Town.
%
%   MODEL   a_cal = M (a_raw - b),  fitted so that ||a_cal|| = g in every hold.
%     'diag'  M = diag(s): per-axis scale and offset, 6 parameters (default; needs >= 7 holds,
%             >= 8 for leave-one-out)
%     'full'  M upper triangular: adds 3 non-orthogonality terms, 9 parameters (needs >= 10 holds;
%             use as a diagnostic only, with 12 or more well-spread holds)
%   Each hold gives ONE equation (the norm), so the number of holds must exceed the number of
%   parameters for the residual to mean anything. The sphere fit cannot see a common rotation
%   of the axes (it only fixes lengths), so misalignment to the board frame is not estimated.
%
%   'Weights'  per-hold weights (default equal).
%   'TwoPointAxis', [iUp iDown axis]  also return the classic flip estimate for one axis
%              (scale = (aUp - aDown) / 2g, offset = (aUp + aDown) / 2), e.g. [1 2 3] for z.
%
%   OUTPUT cal: model, scale (1x3), offset (1x3), M (3x3), residual (nH x 1, m/s^2), rmsResidual,
%          dof, cond (Jacobian condition number), formalStdScale / formalStdOffset (from
%          'HoldNoise', the std of a hold mean, default 2e-4 m/s^2: a floor only), loo
%          (leave-one-hold-out: available (nH >= np + 2), scale, offset, determined, cond, nDetermined, essentialHolds =
%          holds whose removal leaves the fit undetermined, cond > 'MaxCond' default 1e3),
%          looSpreadScale / looSpreadOffset (half-range over the determined refits),
%          twoPoint (if requested), iterations, converged.
%   GAIN: the sensor reads 1/scale(j) times the true acceleration on axis j, so its gain error is
%   cal.gainError(j) = 1/scale(j) - 1 (positive = reads high).
%
% Lives in IMU Characterisation/Processing/.

opt = struct('Model', 'diag', 'Weights', [], 'TwoPointAxis', [], 'MaxIter', 50, 'Tol', 1e-12, ...
    'HoldNoise', 2e-4, 'MaxCond', 1e3);   % MaxCond: numerical diagnostic in these units, not a physical threshold
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i + 1}; end
if size(A, 2) ~= 3, error('calibrateAccelSixPosition:A', 'A must be nH x 3.'); end
nH = size(A, 1);
w = opt.Weights; if isempty(w), w = ones(nH, 1); end
w = w(:);

[x, J, it, conv] = fitModel(A, g, w, opt.Model, opt.MaxIter, opt.Tol);
if isempty(x), error('calibrateAccelSixPosition:holds', 'Too few holds for model ''%s''.', opt.Model); end
[M, b] = unpack(x, opt.Model);
r = residuals(A, M, b, g);
np = numel(x);

cal.model = opt.Model;
cal.M = M;
cal.scale = diag(M).';
cal.offset = b(:).';
cal.gainError = 1 ./ cal.scale - 1;             % sensor reading / true - 1
cal.residual = r;
cal.rmsResidual = sqrt(mean(r.^2));
cal.dof = nH - np;
cal.cond = cond(J);
cal.iterations = it;
cal.converged = conv;

% formal uncertainty from the noise of the hold means only (a floor, not the whole error budget)
Jw = J .* sqrt(w);
if rank(Jw) < np, error('calibrateAccelSixPosition:rank', 'Fit Jacobian is rank deficient.'); end
C = opt.HoldNoise^2 * ((Jw.' * Jw) \ eye(np));
cal.formalStdScale = sqrt(diag(C(1:3, 1:3))).';
if strcmp(opt.Model, 'full'), cal.formalStdScale = sqrt(diag(C([1 4 6], [1 4 6]))).'; end
cal.formalStdOffset = sqrt(diag(C(end - 2:end, end - 2:end))).';

% leave-one-hold-out robustness. Available only when every refit still has more holds than
% parameters (nH >= np + 2). A refit that is ill-conditioned or does not converge is recorded as
% not determined: it shows which holds the calibration cannot do without.
cal.loo = struct('available', nH - 1 > np, 'scale', NaN(nH, 3), 'offset', NaN(nH, 3), ...
    'determined', false(nH, 1), 'cond', NaN(nH, 1));
if cal.loo.available
    for h = 1:nH
        keep = setdiff(1:nH, h);
        [xh, Jh, ~, convh] = fitModel(A(keep, :), g, w(keep), opt.Model, opt.MaxIter, opt.Tol);
        if isempty(xh), continue; end
        cal.loo.cond(h) = cond(Jh);
        if ~convh || cal.loo.cond(h) > opt.MaxCond || any(~isfinite(xh)), continue; end
        [Mh, bh] = unpack(xh, opt.Model);
        cal.loo.scale(h, :) = diag(Mh).';
        cal.loo.offset(h, :) = bh(:).';
        cal.loo.determined(h) = true;
    end
end
cal.loo.nDetermined = nnz(cal.loo.determined);
cal.loo.essentialHolds = [];
if cal.loo.available, cal.loo.essentialHolds = find(~cal.loo.determined).'; end
d = cal.loo.determined;
if any(d)
    cal.looSpreadScale = (max(cal.loo.scale(d, :), [], 1) - min(cal.loo.scale(d, :), [], 1)) / 2;
    cal.looSpreadOffset = (max(cal.loo.offset(d, :), [], 1) - min(cal.loo.offset(d, :), [], 1)) / 2;
else
    cal.looSpreadScale = NaN(1, 3); cal.looSpreadOffset = NaN(1, 3);
end

if ~isempty(opt.TwoPointAxis)
    iu = opt.TwoPointAxis(1); id = opt.TwoPointAxis(2); ax = opt.TwoPointAxis(3);
    sc = (A(iu, ax) - A(id, ax)) / (2 * g);     % reading per unit true acceleration
    cal.twoPoint = struct('axis', ax, 'gain', sc, 'gainError', sc - 1, 'offset', (A(iu, ax) + A(id, ax)) / 2);
end
end

%% ================================================================================================
function [x, J, it, conv] = fitModel(A, g, w, model, maxIter, tol)
switch model
    case 'diag', np = 6;
    case 'full', np = 9;
    otherwise, error('calibrateAccelSixPosition:model', 'Model must be ''diag'' or ''full''.');
end
x = []; J = []; it = 0; conv = false;
if size(A, 1) <= np, return; end                 % need more holds than parameters (residual information)
x = [ones(3, 1); zeros(3, 1)];                   % [s; b] or [M upper (6); b]
if strcmp(model, 'full'), x = [1; 0; 0; 1; 0; 1; zeros(3, 1)]; end
sw = sqrt(w);
for it = 1:maxIter
    [M, b] = unpack(x, model);
    [r, J] = residuals(A, M, b, g, model);
    Jw = J .* sw; rw = r .* sw;
    dx = -(Jw \ rw);
    x = x + dx;
    if norm(dx) < tol * max(1, norm(x)), conv = true; break; end
end
[M, b] = unpack(x, model);
[~, J] = residuals(A, M, b, g, model);
end

function [M, b] = unpack(x, model)
if strcmp(model, 'diag')
    M = diag(x(1:3)); b = x(4:6);
else
    M = [x(1) x(2) x(3); 0 x(4) x(5); 0 0 x(6)]; b = x(7:9);
end
end

function [r, J] = residuals(A, M, b, g, model)
nH = size(A, 1);
D = A - b(:).';                                  % nH x 3
U = D * M.';                                     % calibrated vectors
n = sqrt(sum(U.^2, 2));
r = n - g;
if nargout < 2, return; end
if nargin < 5, model = 'diag'; end
dB = -(U * M) ./ n;                              % d n / d b = -(M' u) / n
if strcmp(model, 'diag')
    dS = U .* D ./ n;                            % d n / d s_j = u_j d_j / n
    J = [dS, dB];
else
    pairs = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];      % (row, col) of the upper-triangular entries
    dM = zeros(nH, 6);
    for k = 1:6
        dM(:, k) = U(:, pairs(k, 1)) .* D(:, pairs(k, 2)) ./ n;
    end
    J = [dM, dB];
end
end