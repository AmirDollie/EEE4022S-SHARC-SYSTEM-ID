%% runTransmissibilityInformation.m
% File 5c of the SSI bridge: how much of the complex-FRF information on
% theta = [ln beta, ln R] (gamma fixed) do sensor TRANSMISSIBILITIES keep?
%
%   T_j(omega) = H_j(omega) / H_r(omega),   j ~= r (reference sensor r).
%
% WHY TRANSMISSIBILITIES: under the present model every sensor sees the
% same single incident field, y_j = H_j u, so the output cross-spectral
% ratio S_jr / S_rr = H_j H_r^* S_u / (|H_r|^2 S_u) = H_j / H_r: the unknown
% incident spectrum cancels. They can be estimated OUTPUT-ONLY, need no
% resonances and no model fit. "Input-independent" holds UNDER THE COMMON-
% INPUT / SINGLE-INCIDENT-FIELD ASSUMPTION; multi-directional seas and
% reference-sensor noise break the exact cancellation (a G2 question).
% What they discard is the common factor H_r(omega): absolute response
% level and phase relative to the incident wave.
%
% NO EMM SOLVES: uses H0 and the ln-parameter derivatives saved by
% runSensorFRFSensitivity.m (file 5b).
%
% NOISE MODEL (identical to file 5b, so the comparison is fair)
%   Reference observable: [Re H_j(omega_k); Im H_j(omega_k)], independent
%   real components with variance ||H(omega_k)||^2 / dOmega per unit
%   relative FRF noise density sigma_n^2. This reproduces file 5b's
%   information matrix exactly (checked in section 1).
%   Transmissibility noise: first-order propagation of that same noise
%   through T = H_j / H_r (dT_j = (dH_j - T_j dH_r) / H_r), so the
%   features are a local linear function of the reference and
%   featureInformation's lambda must lie in [0, 1].
%
% VARIANTS (all with consistent propagated covariance)
%   T      : complex transmissibilities [Re T; Im T]
%   |T|    : log-magnitude only, Re(dT / T)
%   arg T  : phase only, Im(dT / T)
%   for every choice of reference sensor r (primary: PRIMARY_REF).
%   The reference choice cannot matter to first order: ratios to any
%   other reference are an invertible function of these (T' = T / T_r'),
%   and log-magnitudes / phases change by an invertible linear map, with
%   the covariance transforming consistently. Section 2 checks this
%   numerically and prints the primary reference only.
%
% OUTPUTS: generalised eigenvalues lambda (fraction of FRF information
% kept per direction of theta) and eigenvectors, joint stds, cumulative
% joint std versus upper frequency; also the 3-parameter case for the
% primary reference, for the record.
%
% OUTPUT (in Inverse Mapping/Parameter Inversion/Results/)
%   transmissibilityInfo_<stamp>.mat / .txt / _cumulative.png
%
% Lives in Inverse Mapping/Parameter Inversion/SSI Twin/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));      % featureInformation.m (Inverse Mapping)
clearvars -except thisDir
close all, clc

%% ---- Configuration -----------------------------------------------------------
SENS_FILE = '';                              % '' = newest sensorFRFSensitivity_*.mat in Results
THETA_IDX = [1 3];                           % [beta R]; gamma fixed
PRIMARY_REF = 4;                             % s26 (0.9R, 0): largest response
REF_CHECK_2P = [2.3154 0.2627];              % file 5b, gamma fixed (sanity check)
REF_CHECK_3P = [3.21 39.3 0.406];            % file 5b, three parameters

%% ---- Output ----------------------------------------------------------------------
resultsDir = fullfile(thisDir, '..', 'Results');
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['transmissibilityInfo_', stamp]);
diary([baseName, '.txt']);
try

if isempty(SENS_FILE)
    d = dir(fullfile(resultsDir, 'sensorFRFSensitivity_*.mat'));
    if isempty(d)
        error('runTransmissibilityInformation:noInput', 'No sensorFRFSensitivity_*.mat in %s (run file 5b first).', resultsDir);
    end
    [~, newest] = max([d.datenum]);
    SENS_FILE = fullfile(resultsDir, d(newest).name);
end
S = load(SENS_FILE); sens = S.sens;
H0 = sens.H0;                                % nS x nW
dH = sens.dH(:, :, :, end);                  % primary step
omega = sens.omega; nW = numel(omega); nS = size(H0, 1);
dW = omega(2) - omega(1);
normH = vecnorm(H0, 2, 1);
sensorNames = {'s2 (0.3R,0)', 's6 (0.3R,pi)', 's14 (0.5R,pi)', 's26 (0.9R,0)'};

fprintf('runTransmissibilityInformation  (%s)\n', stamp);
fprintf('input: %s\n', SENS_FILE);
fprintf('%d frequencies in [%.2f %.2f] rad/s, %d sensors, step eps %.0e; theta = [ln %s, ln %s], gamma fixed\n\n', ...
    nW, omega(1), omega(end), nS, sens.epsList(end), sens.paramNames{THETA_IDX});

%% ---- 1  Reference (complex FRF) and consistency with file 5b ------------------------
fprintf('=== 1: reference information (complex FRF) ===\n');
[Jr2, Sr2] = frfReference(dH(:, :, THETA_IDX), normH, dW);
[Jr3, Sr3] = frfReference(dH, normH, dW);
Fr2 = Jr2.' * (Jr2 ./ Sr2); Fr3 = Jr3.' * (Jr3 ./ Sr3);
std2 = sqrt(diag(inv(Fr2))).'; std3 = sqrt(diag(inv(Fr3))).';
fprintf('  2-parameter joint std: beta %.4f  R %.4f   (file 5b: %.4f %.4f)\n', std2, REF_CHECK_2P);
fprintf('  3-parameter joint std: beta %.3g  gamma %.3g  R %.3g   (file 5b: %.3g %.3g %.3g)\n', std3, REF_CHECK_3P);
if max(abs(std2 - REF_CHECK_2P) ./ REF_CHECK_2P) > 1e-3
    warning('runTransmissibilityInformation:refMismatch', '%s', ...
        'Reference does not reproduce the file 5b values; check the input file.');
end
fprintf('\n');

%% ---- 2  Transmissibility variants for every reference sensor ---------------------------
fprintf('=== 2: information kept, theta = [ln beta, ln R] ===\n');
fprintf('(lambda = fraction of FRF information kept along the corresponding direction)\n');
variants = {'T', 'logMag', 'phase'};
varLabels = {'T (complex)', '|T| only', 'arg T only'};
res = struct();
for r = 1:nS
    for v = 1:numel(variants)
        [Jf, Sf] = transmissibilityFeatures(H0, dH(:, :, THETA_IDX), normH, dW, r, variants{v});
        res(r, v).out = featureInformation(Jf, Sf, Jr2, Sr2, 'ParamNames', {'lnBeta', 'lnR'}); %#ok<SAGROW>
    end
end
refSpread = 0;
for v = 1:numel(variants)
    L = cell2mat(arrayfun(@(r) res(r, v).out.lambda, 1:nS, 'UniformOutput', false));
    refSpread = max(refSpread, max(max(abs(L - L(:, PRIMARY_REF)))) / max(L(:)));
end
fprintf('  reference-sensor invariance: max relative spread of lambda over the %d choices = %.1e\n', nS, refSpread);
fprintf('  %-12s %9s %9s  %-18s %9s %9s %9s %9s\n', 'variant', 'lambda1', 'lambda2', 'weak v (std units)', ...
    'std beta', 'std R', 'x FRF', 'x FRF');
for v = 1:numel(variants)
    o = res(PRIMARY_REF, v).out;
    fprintf('  %-12s %9.4f %9.4f  %-18s %9.3f %9.3f %9.2f %9.2f%s\n', varLabels{v}, o.lambda, ...
        mat2str(o.Vscaled(:, 2).', 3), o.stdF, o.stdRatio, repmat('  [inconsistent!]', 1, ~o.consistent));
end
fprintf('  %-12s %9.4f %9.4f  %-18s %9.3f %9.3f\n', 'complex FRF', 1, 1, '(reference)', std2);
fprintf('  ("x FRF" = std ratio to the complex FRF; weak v = direction of least retained information, written\n');
fprintf('   in units of the reference stds [beta R], so it shows which parameter actually loses information)\n');
fprintf('\n');

%% ---- 3  Three-parameter case, primary reference (for the record) ----------------------
fprintf('=== 3: three parameters, reference %s ===\n', sensorNames{PRIMARY_REF});
[Jf3, Sf3] = transmissibilityFeatures(H0, dH, normH, dW, PRIMARY_REF, 'T');
o3 = featureInformation(Jf3, Sf3, Jr3, Sr3, 'ParamNames', {'lnBeta', 'lnGamma', 'lnR'});
fprintf('  lambda %s\n', mat2str(o3.lambda.', 4));
fprintf('  eigenvectors in reference-std units (columns, [beta gamma R]):\n'); disp(round(o3.Vscaled * 1000) / 1000);
fprintf('  joint std  T: beta %.3g gamma %.3g R %.3g   FRF: beta %.3g gamma %.3g R %.3g\n\n', o3.stdF, o3.stdR);

%% ---- 4  Cumulative joint std vs upper frequency (primary reference) --------------------
fprintf('=== 4: cumulative joint std vs upper frequency, reference %s ===\n', sensorNames{PRIMARY_REF});
curves = struct('label', {'complex FRF', 'T (complex)', '|T| only', 'arg T only'}, 'std', []);
Fk = perFrequencyInfo(Jr2, Sr2, 2 * nS);
curves(1).std = cumulativeStd(Fk);
for v = 1:numel(variants)
    [Jf, Sf] = transmissibilityFeatures(H0, dH(:, :, THETA_IDX), normH, dW, PRIMARY_REF, variants{v});
    blk = size(Jf, 1) / nW;
    curves(v + 1).std = cumulativeStd(perFrequencyInfoFull(Jf, Sf, blk));
end
kk = unique(round(linspace(1, nW, 7)));
fprintf('  %-9s %s\n', 'omega_max', sprintf('%-24s', curves.label));
for k = kk
    row = '';
    for c = 1:numel(curves)
        row = [row, sprintf('b %7.3g  R %7.3g    ', curves(c).std(:, k))]; %#ok<AGROW>
    end
    fprintf('  %9.2f %s\n', omega(k), row);
end
fprintf('\n');

%% ---- 5  Save ---------------------------------------------------------------------------------
out = struct('sensFile', SENS_FILE, 'omega', omega, 'thetaIdx', THETA_IDX, 'primaryRef', PRIMARY_REF, ...
    'variants', {variants}, 'variantLabels', {varLabels}, 'sensorNames', {sensorNames}, ...
    'refStd2', std2, 'refStd3', std3, 'result', res, 'threeParam', o3, 'cumulative', curves);
save([baseName, '.mat'], 'out');
fprintf('Saved: %s.mat\n', baseName);

%% ---- 6  Figure -----------------------------------------------------------------------------
colours = {[42 120 214]/255, [235 104 52]/255, [27 175 122]/255, [237 161 0]/255};   % validated slots 1-4
markers = {'o', 's', '^', 'd'};
ink = [0.35 0.35 0.33]; gridInk = [0.85 0.85 0.83];
mEvery = max(1, round(nW / 10));
fig = figure('Color', 'w', 'Position', [60 60 1100 430]);
pNames = {'\beta', 'R'};
for pp = 1:2
    ax = subplot(1, 2, pp); hold(ax, 'on');
    h = gobjects(1, numel(curves));
    for c = 1:numel(curves)
        h(c) = plot(ax, omega, curves(c).std(pp, :), '-', 'Color', colours{c}, 'LineWidth', 2, ...
            'Marker', markers{c}, 'MarkerIndices', 1:mEvery:nW, 'MarkerSize', 6, ...
            'MarkerFaceColor', colours{c}, 'MarkerEdgeColor', 'w');
    end
    set(ax, 'YScale', 'log', 'XColor', ink, 'YColor', ink, 'FontSize', 10, 'Box', 'off');
    grid(ax, 'on'); set(ax, 'GridColor', gridInk, 'GridAlpha', 1);
    xlim(ax, omega([1 end]));
    xlabel(ax, 'upper frequency \omega_{max} (rad/s)');
    ylabel(ax, sprintf('joint std of ln %s (per unit relative noise density)', pNames{pp}));
    title(ax, sprintf('(%c) ln %s, \\gamma fixed, reference %s', 'a' + pp - 1, pNames{pp}, sensorNames{PRIMARY_REF}), ...
        'FontWeight', 'normal');
    if pp == 1, legend(ax, h, {curves.label}, 'Location', 'best', 'Box', 'off'); end
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
function [J, Sv] = frfReference(dH, normH, dW)
% Rows [Re dH(:,k,:); Im dH(:,k,:)] per frequency; variance ||H_k||^2 / dW
% for every real component (file 5b's noise model).
[nS, nW, nP] = size(dH);
J = zeros(2 * nS * nW, nP); Sv = zeros(2 * nS * nW, 1);
for k = 1:nW
    b = reshape(dH(:, k, :), nS, nP);
    rows = (k - 1) * 2 * nS + (1:2 * nS);
    J(rows, :) = [real(b); imag(b)];
    Sv(rows) = normH(k)^2 / dW;
end
end

function [J, Sig] = transmissibilityFeatures(H0, dH, normH, dW, r, variant)
% First-order transmissibility features and their propagated covariance.
[nS, nW, nP] = size(dH);
others = setdiff(1:nS, r);
m = numel(others);
switch variant
    case 'T',      rowsKeep = 1:2 * m;
    case 'logMag', rowsKeep = 1:m;
    case 'phase',  rowsKeep = m + 1:2 * m;
end
blk = numel(rowsKeep);
J = zeros(blk * nW, nP);
Sig = zeros(blk * nW);
for k = 1:nW
    h = H0(:, k);
    T = h(others) / h(r);
    A = zeros(m, nS);                                  % complex: dT = A dH
    A(:, others) = eye(m) / h(r);
    A(:, r) = -T / h(r);
    if ~strcmp(variant, 'T')
        A = diag(1 ./ T) * A;                          % d ln T = dT / T
    end
    Ar = [real(A), -imag(A); imag(A), real(A)];        % real form on [Re dH; Im dH]
    Ar = Ar(rowsKeep, :);
    b = reshape(dH(:, k, :), nS, nP);
    rows = (k - 1) * blk + (1:blk);
    J(rows, :) = Ar * [real(b); imag(b)];
    Sig(rows, rows) = (normH(k)^2 / dW) * (Ar * Ar.');
end
end

function Fk = perFrequencyInfo(J, Sv, blk)
% Per-frequency information contributions for a diagonal covariance.
nW = size(J, 1) / blk; nP = size(J, 2);
Fk = zeros(nP, nP, nW);
for k = 1:nW
    rows = (k - 1) * blk + (1:blk);
    Fk(:, :, k) = J(rows, :).' * (J(rows, :) ./ Sv(rows));
end
end

function Fk = perFrequencyInfoFull(J, Sig, blk)
% Per-frequency information contributions for a block-diagonal covariance.
nW = size(J, 1) / blk; nP = size(J, 2);
Fk = zeros(nP, nP, nW);
for k = 1:nW
    rows = (k - 1) * blk + (1:blk);
    Jk = J(rows, :);
    Fk(:, :, k) = Jk.' * (Sig(rows, rows) \ Jk);
end
end

function st = cumulativeStd(Fk)
nP = size(Fk, 1); nW = size(Fk, 3);
st = NaN(nP, nW);
F = zeros(nP);
for k = 1:nW
    F = F + Fk(:, :, k);
    if rcond(F) > 1e-14
        st(:, k) = sqrt(diag(inv(F)));
    end
end
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);     % R2020a+
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end