%% runOutputOnlyInformationBound.m
% File 5d of the SSI bridge: how much (ln beta, ln R) information can ANY
% output-only method (SSI poles included) extract, given what is assumed
% about the unknown wave (input) spectrum? gamma fixed.
%
% WHY: with a single common input, y = H u, the output cross-spectral
% matrix is S_y(omega) = H H^H S_u. At each frequency this is fully
% described by
%     T_j = H_j / H_r   (j ~= r, transmissibilities: spatial structure)
%     l   = ln(|H_r|^2 S_u)   (log auto-spectrum level of the reference)
% The phase of H_r relative to the incident wave is gone for good, so the
% complex-FRF reference of file 5b is unreachable output-only. Anything
% derived from the outputs (poles, damping, shapes) is a function of
% {T, l} and cannot carry more information than they do. The auto-
% spectrum level only helps to the extent that S_u is known, because
% S_u and |H_r|^2 enter l as a product.
%
% SCOPE (keep the claim rigorous): an output-only SECOND-ORDER information
% bound, under the single-common-input / single-incident-field model and
% the SAME local noise surrogate as files 5b/5c (independent real FRF
% components with variance ||H(omega)||^2/dOmega per unit relative noise
% density). The real finite-record covariance of estimated cross-spectra
% (e.g. log-spectrum variance set by the number of averages) comes later,
% in G2.
%
% METHOD
%   Per frequency, one transformation g: [Re H; Im H] (8) -> [Re T; Im T; l]
%   (7), linearised: dz = G_k [Re dH; Im dH] + e_l * d ln S_u(omega_k), with
%       dT = (dH_j - T_j dH_r) / H_r,    dl = 2 Re(dH_r / H_r) + d ln S_u.
%   T and l share the noisy reference channel, so their covariance is
%   propagated JOINTLY: Sigma_k = s_k^2 G_k G_k' (never treated as
%   independent, which would double-count H_r and allow lambda > 1).
%   Unknown spectrum parameters eta (ln S_u = B eta) are removed by the
%   Schur complement, computed stably: whiten the full Jacobian, project
%   the theta columns onto the orthogonal complement of the nuisance
%   columns, J_eff = (I - Q Q') W J_theta. Then featureInformation compares
%   J_eff with the complex-FRF reference.
%
% VARIANTS (nested nuisance spaces, so information can only decrease downwards)
%   B   T + l, S_u exactly known (e.g. measured by a separate wave buoy)
%   C   T + l, S_u known up to one overall scale
%   D2, D4, D6  T + l, ln S_u in a smooth orthonormal-polynomial basis of
%       dimension 2, 4, 6 over normalised frequency (unknown smooth sea state)
%   J   T + l, S_u a JONSWAP spectrum with unknown (Hs, omega_p, gamma_J)
%       (local: derivative columns at a nominal sea state; not nested with D)
%   E   T + l, ln S_u arbitrary at every frequency (one nuisance per bin)
%   A   T only; E must reproduce A exactly (internal check), and A must
%       reproduce file 5c
%   and the complex FRF reference (lambda = 1).
%
% CHECKS printed: every lambda <= 1; E == A; PSD ordering
%   FRF >= B >= C >= D2 >= D4 >= D6 >= E.
%
% NO EMM SOLVES (uses the file 5b .mat). Runs in seconds.
%
% OUTPUT (in Inverse Mapping/Parameter Inversion/Results/)
%   outputOnlyBound_<stamp>.mat / .txt / _cumulative.png
%
% Lives in Inverse Mapping/Parameter Inversion/SSI Twin/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));      % featureInformation.m (Inverse Mapping)
clearvars -except thisDir
close all, clc

%% ---- Configuration -----------------------------------------------------------
SENS_FILE = '';                              % '' = newest sensorFRFSensitivity_*.mat in Results
THETA_IDX = [1 3];                           % [beta R]; gamma fixed
REF = 4;                                     % reference sensor for T and l (s26); T-information is
                                             % reference-invariant, l depends on it
POLY_DIMS = [2 4 6];                         % smooth ln S_u bases
JONSWAP_NOMINAL = struct('Hs', 0.1, 'wp', 5.5, 'gammaJ', 3.3);   % basin-scale nominal sea state
REF_CHECK_5C = [3.985 0.784];                % file 5c, complex T, joint std [beta R]

%% ---- Output ----------------------------------------------------------------------
resultsDir = fullfile(thisDir, '..', 'Results');
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['outputOnlyBound_', stamp]);
diary([baseName, '.txt']);
try

if isempty(SENS_FILE)
    d = dir(fullfile(resultsDir, 'sensorFRFSensitivity_*.mat'));
    if isempty(d)
        error('runOutputOnlyInformationBound:noInput', 'No sensorFRFSensitivity_*.mat in %s (run file 5b first).', resultsDir);
    end
    [~, newest] = max([d.datenum]);
    SENS_FILE = fullfile(resultsDir, d(newest).name);
end
S = load(SENS_FILE); sens = S.sens;
H0 = sens.H0; dH = sens.dH(:, :, THETA_IDX, end);
omega = sens.omega; nW = numel(omega); nS = size(H0, 1);
dW = omega(2) - omega(1);
normH = vecnorm(H0, 2, 1);
sensorNames = {'s2 (0.3R,0)', 's6 (0.3R,pi)', 's14 (0.5R,pi)', 's26 (0.9R,0)'};

fprintf('runOutputOnlyInformationBound  (%s)\n', stamp);
fprintf('input: %s\n', SENS_FILE);
fprintf('%d frequencies in [%.2f %.2f] rad/s; theta = [ln beta, ln R], gamma fixed; reference sensor %s\n\n', ...
    nW, omega(1), omega(end), sensorNames{REF});

%% ---- 1  Per-frequency transformation, covariance and parameter Jacobian ------------------
% z_k = [Re T (m); Im T (m); l] ; rows of G_k act on [Re dH; Im dH].
m = nS - 1; others = setdiff(1:nS, REF);
nz = 2 * m + 1; lRow = nz;
G = zeros(nz, 2 * nS, nW); Jz = zeros(nz, 2, nW); s2 = zeros(1, nW);
for k = 1:nW
    h = H0(:, k); hr = h(REF);
    T = h(others) / hr;
    A = zeros(m, nS); A(:, others) = eye(m) / hr; A(:, REF) = -T / hr;   % dT = A dH (complex)
    a = zeros(1, nS); a(REF) = 2 / hr;                                   % dl = Re(a dH)
    G(:, :, k) = [real(A), -imag(A); imag(A), real(A); real(a), -imag(a)];
    b = reshape(dH(:, k, :), nS, 2);
    Jz(:, :, k) = G(:, :, k) * [real(b); imag(b)];
    s2(k) = normH(k)^2 / dW;
end

% complex-FRF reference (identical to files 5b/5c)
Jr = zeros(2 * nS * nW, 2); Sr = zeros(2 * nS * nW, 1);
for k = 1:nW
    b = reshape(dH(:, k, :), nS, 2);
    rows = (k - 1) * 2 * nS + (1:2 * nS);
    Jr(rows, :) = [real(b); imag(b)];
    Sr(rows) = s2(k);
end
Fr = Jr.' * (Jr ./ Sr);
fprintf('reference (complex FRF) joint std: beta %.4f  R %.4f\n\n', sqrt(diag(inv(Fr))));

%% ---- 2  Nuisance bases for ln S_u ------------------------------------------------------------
x = 2 * (omega - omega(1)) / (omega(end) - omega(1)) - 1;      % normalised frequency in [-1, 1]
polyBasis = @(dim) orthBasis(x(:) .^ (0:dim - 1));
jonswapBasis = jonswapLogDerivatives(omega, JONSWAP_NOMINAL);   % nW x 3: d ln S / d ln(Hs, wp, gammaJ)

variants = struct('key', {}, 'label', {}, 'rows', {}, 'B', {});
variants(end + 1) = struct('key', 'A',  'label', 'T only',                       'rows', 1:2 * m, 'B', []);
variants(end + 1) = struct('key', 'B',  'label', 'T + l, S_u known',             'rows', 1:nz,   'B', []);
variants(end + 1) = struct('key', 'C',  'label', 'T + l, S_u up to scale',       'rows', 1:nz,   'B', ones(nW, 1));
for dd = POLY_DIMS
    variants(end + 1) = struct('key', sprintf('D%d', dd), 'label', sprintf('T + l, smooth ln S_u (dim %d)', dd), ...
        'rows', 1:nz, 'B', polyBasis(dd)); %#ok<SAGROW>
end
variants(end + 1) = struct('key', 'J',  'label', 'T + l, JONSWAP (Hs, wp, gJ) unknown', 'rows', 1:nz, 'B', jonswapBasis);
variants(end + 1) = struct('key', 'E',  'label', 'T + l, S_u arbitrary per bin', 'rows', 1:nz,   'B', eye(nW));

%% ---- 3  Information per variant ---------------------------------------------------------------
fprintf('=== 3: output-only information on theta = [ln beta, ln R] (vs complex FRF) ===\n');
fprintf('  %-4s %-36s %8s %8s  %-16s %8s %8s %7s %7s\n', 'key', 'variant', 'lambda1', 'lambda2', ...
    'weak v (std u.)', 'std b', 'std R', 'x FRF', 'x FRF');
res = struct('key', {}, 'label', {}, 'F', {}, 'out', {});
for v = 1:numel(variants)
    Jeff = effectiveJacobian(1:nW, variants(v), G, Jz, s2, lRow);
    o = featureInformation(Jeff, 1, Jr, Sr, 'ParamNames', {'lnBeta', 'lnR'}, 'WarnInconsistent', false);
    res(v).key = variants(v).key; res(v).label = variants(v).label; res(v).F = o.Ff; res(v).out = o;
    fprintf('  %-4s %-36s %8.4f %8.4f  %-16s %8.3f %8.3f %7.2f %7.2f%s\n', variants(v).key, variants(v).label, ...
        o.lambda, mat2str(o.Vscaled(:, 2).', 3), o.stdF, o.stdRatio, repmat('  [lambda > 1!]', 1, ~o.consistent));
end
fprintf('  %-4s %-36s %8.4f %8.4f  %-16s %8.3f %8.3f\n', 'F', 'complex FRF (unreachable output-only)', 1, 1, ...
    '(reference)', sqrt(diag(inv(Fr))));
fprintf('\n');

%% ---- 4  Checks ------------------------------------------------------------------------------------
fprintf('=== 4: checks ===\n');
key = @(k) find(strcmp({res.key}, k), 1);
FA = res(key('A')).F; FE = res(key('E')).F;
dAE = norm(FE - FA) / norm(FA);
fprintf('  arbitrary-S_u (E) vs T only (A): relative difference %.2e  %s\n', dAE, passFail(dAE < 1e-8));
stdA = res(key('A')).out.stdF.';
dC5 = max(abs(stdA - REF_CHECK_5C) ./ REF_CHECK_5C);
fprintf('  T only vs file 5c joint std [%.4f %.4f] vs [%.4f %.4f]: %s\n', stdA, REF_CHECK_5C, passFail(dC5 < 2e-3));
allCons = all(arrayfun(@(r) r.out.consistent, res));
fprintf('  all lambda <= 1: %s\n', passFail(allCons));
chain = [{'B', 'C'}, arrayfun(@(dd) sprintf('D%d', dd), POLY_DIMS, 'UniformOutput', false), {'E'}];
prevF = Fr; prevName = 'FRF'; orderOK = true;
for c = 1:numel(chain)
    Fc = res(key(chain{c})).F;
    mineig = min(eig((prevF - Fc + (prevF - Fc).') / 2)) / norm(prevF);
    ok = mineig > -1e-10;
    orderOK = orderOK && ok;
    fprintf('  %s >= %s (min eig of difference / ||F||: %+.2e)  %s\n', prevName, chain{c}, mineig, passFail(ok));
    prevF = Fc; prevName = chain{c};
end
checks = struct('dAE', dAE, 'dC5', dC5, 'allConsistent', allCons, 'orderOK', orderOK);
fprintf('\n');

%% ---- 5  Cumulative joint std vs upper frequency (selected variants) -----------------------------
fprintf('=== 5: cumulative joint std vs upper frequency ===\n');
sel = {'A', 'B', 'C', 'D4', 'J'};
cum = struct('key', [{'FRF'}, sel], 'std', []);
Fk = zeros(2, 2, nW);
for k = 1:nW
    rows = (k - 1) * 2 * nS + (1:2 * nS);
    Fk(:, :, k) = Jr(rows, :).' * (Jr(rows, :) ./ Sr(rows));
end
cum(1).std = NaN(2, nW); Fc = zeros(2);
for k = 1:nW
    Fc = Fc + Fk(:, :, k);
    if rcond(Fc) > 1e-14, cum(1).std(:, k) = sqrt(diag(inv(Fc))); end
end
for c = 1:numel(sel)
    vv = variants(strcmp({variants.key}, sel{c}));
    cum(c + 1).std = NaN(2, nW);
    for k = 1:nW
        vk = vv;
        if ~isempty(vv.B), vk.B = vv.B(1:k, :); end
        Jeff = effectiveJacobian(1:k, vk, G, Jz, s2, lRow);
        Fk2 = Jeff.' * Jeff;
        if rcond(Fk2) > 1e-14, cum(c + 1).std(:, k) = sqrt(diag(inv(Fk2))); end
    end
end
kk = unique(round(linspace(1, nW, 7)));
fprintf('  %-9s %s\n', 'omega_max', sprintf('%-22s', cum.key));
for k = kk
    row = '';
    for c = 1:numel(cum)
        row = [row, sprintf('b %6.3g R %6.3g   ', cum(c).std(:, k))]; %#ok<AGROW>
    end
    fprintf('  %9.2f %s\n', omega(k), row);
end
fprintf('\n');

%% ---- 6  Save -------------------------------------------------------------------------------------
out = struct('sensFile', SENS_FILE, 'omega', omega, 'thetaIdx', THETA_IDX, 'refSensor', REF, ...
    'variants', variants, 'result', res, 'Fr', Fr, 'checks', checks, 'cumulative', cum, ...
    'jonswapNominal', JONSWAP_NOMINAL);
save([baseName, '.mat'], 'out');
fprintf('Saved: %s.mat\n', baseName);

%% ---- 7  Figure ---------------------------------------------------------------------------------------
colours = {[42 120 214]/255, [235 104 52]/255, [27 175 122]/255, [237 161 0]/255, [232 123 164]/255, [0 131 0]/255};
markers = {'o', 's', '^', 'd', 'v', 'p'};
labels = {'complex FRF', 'T only', 'T + l, S_u known', 'T + l, S_u up to scale', 'T + l, smooth S_u (dim 4)', ...
    'T + l, JONSWAP unknown'};
ink = [0.35 0.35 0.33]; gridInk = [0.85 0.85 0.83];
mEvery = max(1, round(nW / 10));
fig = figure('Color', 'w', 'Position', [60 60 1150 450]);
pNames = {'\beta', 'R'};
for pp = 1:2
    ax = subplot(1, 2, pp); hold(ax, 'on');
    h = gobjects(1, numel(cum));
    for c = 1:numel(cum)
        h(c) = plot(ax, omega, cum(c).std(pp, :), '-', 'Color', colours{c}, 'LineWidth', 2, ...
            'Marker', markers{c}, 'MarkerIndices', 1:mEvery:nW, 'MarkerSize', 6, ...
            'MarkerFaceColor', colours{c}, 'MarkerEdgeColor', 'w');
    end
    set(ax, 'YScale', 'log', 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off');
    grid(ax, 'on'); set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
    xlim(ax, omega([1 end]));
    xlabel(ax, 'upper frequency \omega_{max} (rad/s)');
    ylabel(ax, sprintf('joint std of ln %s (per unit relative noise density)', pNames{pp}));
    title(ax, sprintf('(%c) ln %s, \\gamma fixed: output-only bounds', 'a' + pp - 1, pNames{pp}), 'FontWeight', 'normal');
    if pp == 1, legend(ax, h, labels, 'Location', 'best', 'Box', 'off'); end
end
saveFigure(fig, [baseName, '_cumulative.png']);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ================================================================================================
%  Local functions
%% ================================================================================================
function Jeff = effectiveJacobian(kIdx, variant, G, Jz, s2, lRow)
% Whitened theta-Jacobian with the nuisance (ln S_u) directions projected out.
rows = variant.rows;
hasNuis = ~isempty(variant.B) && any(rows == lRow);
nr = numel(rows); nk = numel(kIdx);
Jt = zeros(nr * nk, size(Jz, 2));
if hasNuis
    q = size(variant.B, 2);
    Jn = zeros(nr * nk, q);
end
for i = 1:nk
    k = kIdx(i);
    Gk = G(rows, :, k);
    L = chol(s2(k) * (Gk * Gk.'), 'lower');          % joint covariance of the kept rows
    blk = (i - 1) * nr + (1:nr);
    Jt(blk, :) = L \ Jz(rows, :, k);
    if hasNuis
        e = double(rows(:) == lRow);                  % ln S_u enters the l row only
        Jn(blk, :) = L \ (e * variant.B(i, :));
    end
end
Jeff = Jt;
if hasNuis
    [Q, Rq] = qr(Jn, 0);
    dR = abs(diag(Rq));
    rk = sum(dR > 1e-10 * max(dR));
    Q = Q(:, 1:rk);
    Jeff = Jt - Q * (Q.' * Jt);                       % orthogonal complement of the nuisance columns
end
end

function Q = orthBasis(V)
[Q, ~] = qr(V, 0);
end

function D = jonswapLogDerivatives(omega, p)
% Columns: d ln S / d ln Hs, d ln S / d ln wp, d ln S / d ln gammaJ at the
% nominal JONSWAP (central differences in the log-parameters).
f = @(Hs, wp, gJ) log(jonswap(omega, Hs, wp, gJ));
h = 1e-4;
D = zeros(numel(omega), 3);
D(:, 1) = (f(p.Hs * exp(h), p.wp, p.gammaJ) - f(p.Hs * exp(-h), p.wp, p.gammaJ)).' / (2 * h);
D(:, 2) = (f(p.Hs, p.wp * exp(h), p.gammaJ) - f(p.Hs, p.wp * exp(-h), p.gammaJ)).' / (2 * h);
D(:, 3) = (f(p.Hs, p.wp, p.gammaJ * exp(h)) - f(p.Hs, p.wp, p.gammaJ * exp(-h))).' / (2 * h);
end

function S = jonswap(w, Hs, wp, gJ)
% Standard JONSWAP shape (Pierson-Moskowitz tail, peak enhancement gJ),
% scaled so that Hs = 4 sqrt(m0) holds over the full spectrum
% (the scale only matters through ln Hs, which is a pure level shift).
sig = 0.07 * ones(size(w)); sig(w > wp) = 0.09;
r = exp(-(w - wp).^2 ./ (2 * sig.^2 * wp^2));
S = w.^-5 .* exp(-1.25 * (wp ./ w).^4) .* gJ.^r;
S = S * (Hs / 4)^2 / (wp^-4 / 5);             % approximate m0 normalisation (level only)
end

function s = passFail(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);     % R2020a+
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end