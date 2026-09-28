%% runFloeRegimeScreen.m
% Track B screening table: the basin reference p0, a candidate tank floe (B1) and
% representative Antarctic MIZ floes (B2), through screenFloeRegime.m. No EMM runs.
%
% For each case: flexural length l_f, R/l_f, and over the case's period band kR, k l_f,
% the ice-covered wavenumber ratio kappa/k and the rigidity sensitivity d ln kappa / d ln D,
% plus the case expressed in the EMM's scaling (alpha, beta, gamma, R/H) and whether its
% alpha band lies inside the validated range [1.7 13.85].
%
% READING IT: p0 is the only case where beta is known to be identifiable (U2). The
% infinite-plate dispersion response (kappa / k, d ln kappa / d ln D) is organised mainly by
% k l_f; finite-floe size effects enter through kR and R / l_f (the latter also carries the
% stiffness). A case overlapping p0 in k l_f with comparable |d ln kappa / d ln D| occupies a
% similar flexural wave-scale regime and is a candidate for an EMM run; finite-floe
% similarity still needs kR and R / l_f to match. The EMM-scaled beta of a deep-water case
% depends on the chosen effective depth, so it is a numerical-regime note, not a similarity
% measure. None of these numbers decides identifiability (dT/d ln beta and F = J' Sigma^-1 J,
% Stage 3).
%
% THE TANK CASE IS ILLUSTRATIVE (polypropylene sheet, 1 m flume): replace it with the real
% tank depth, wave-maker range and a realisable floe once known.
%
% OUTPUT (Mission Physics/Results): floeRegimeScreen_<stamp>.mat / .txt / .png
%
% Lives in Mission Physics/.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
clearvars -except thisDir
close all, clc

%% ---- Cases -------------------------------------------------------------------------------------------
g = 9.81;
P0 = [4.6985e-5, 1.4548e-3, 0.3830]; H0 = 1.88; RHO0 = 1000;          % basin reference (EMM scaling)
nT = 60;
cases = struct('name', {}, 'floe', {}, 'T', {}, 'opts', {});
cases(end + 1) = struct('name', 'basin p0', ...
    'floe', struct('R', P0(3) * H0, 'D', P0(1) * RHO0 * g * H0^4, 'm', P0(2) * H0 * RHO0), ...
    'T', 2 * pi ./ linspace(2.979, 8.501, nT), 'opts', {{'Depth', H0, 'RhoWater', RHO0}});
cases(end + 1) = struct('name', 'tank (illustrative)', ...
    'floe', struct('R', 0.5, 'h', 0.010, 'E', 1.6e9, 'nu', 0.40, 'rhoIce', 905), ...
    'T', linspace(0.6, 2.0, nT), 'opts', {{'Depth', 1.0, 'RhoWater', 1000}});
MIZ_ICE = struct('h', 1.0, 'E', 6e9, 'nu', 0.3, 'rhoIce', 917);
MIZ_T = linspace(6, 15, nT);
for R = [10 25 50]
    fl = MIZ_ICE; fl.R = R;
    cases(end + 1) = struct('name', sprintf('MIZ R = %g m, h = 1 m', R), 'floe', fl, 'T', MIZ_T, ...
        'opts', {{'Depth', Inf}}); %#ok<SAGROW>
end
fl = struct('R', 25, 'h', 0.5, 'E', 3e9, 'nu', 0.3, 'rhoIce', 917);
cases(end + 1) = struct('name', 'MIZ R = 25 m, h = 0.5 m, E = 3 GPa', 'floe', fl, 'T', MIZ_T, 'opts', {{'Depth', Inf}});

%% ---- Output ----------------------------------------------------------------------------------------------
resultsDir = fullfile(thisDir, 'Results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['floeRegimeScreen_', stamp]);
diary([baseName, '.txt']);
try

fprintf('runFloeRegimeScreen  (%s)\n', stamp);
fprintf('Screening only: no EMM runs, no identifiability verdicts.\n\n');
res = cell(1, numel(cases));
for c = 1:numel(cases)
    s = screenFloeRegime(cases(c).floe, cases(c).T, cases(c).opts{:});
    res{c} = s;
    fprintf('=== %s ===\n', cases(c).name);
    if isnan(s.h)
        fprintf('  R = %.4g m, D = %.4g N m, m = %.4g kg/m^2 (draught %.3g mm)\n', s.R, s.D, s.m, 1e3 * s.draught);
    else
        fprintf('  R = %.4g m, h = %.4g m, E = %.3g GPa, nu = %.2f -> D = %.4g N m, draught %.3g m\n', ...
            s.R, s.h, s.E / 1e9, s.nu, s.D, s.draught);
    end
    fprintf('  l_f = %.4g m, R/l_f = %.3f;  %s\n', s.flexuralLength, s.RoverLf, s.emm.note);
    fprintf('  %7s %8s %9s %7s %7s %8s %9s %8s %7s\n', 'T (s)', 'lambda', 'kR', 'k l_f', 'kH', 'kap/k', 'dlnk/dlnD', 'alpha', 'in EMM');
    idx = unique(round(linspace(1, numel(s.T), 5)));
    [~, ord] = sort(s.T(idx)); idx = idx(ord);
    for i = idx
        fprintf('  %7.3f %8.3g %9.3f %7.3f %7.3g %8.3f %9.4f %8.3f %7s\n', s.T(i), s.lambda(i), s.kR(i), s.kLf(i), ...
            s.kH(i), s.kappaRatio(i), s.dlnKappa_dlnD(i), s.emm.alpha(i), yesno(s.emm.alphaInRange(i)));
    end
    fprintf('  bands: kR %.3g-%.3g, k l_f %.3g-%.3g, |dln kappa/dln D| %.3g-%.3g; EMM: alpha %.3g-%.3g (%d%% in range), beta %.3g, gamma %.3g, R/H %.3g\n\n', ...
        min(s.kR), max(s.kR), min(s.kLf), max(s.kLf), min(abs(s.dlnKappa_dlnD)), max(abs(s.dlnKappa_dlnD)), ...
        min(s.emm.alpha), max(s.emm.alpha), round(100 * mean(s.emm.alphaInRange)), s.emm.beta, s.emm.gamma, s.emm.Rnd);
end

% overlap with the validated reference
b = res{1};
fprintf('OVERLAP WITH p0 (fraction of each case''s log-range inside p0''s log-range; sampling-independent)\n');
fprintf('  %-36s %8s %8s %12s %10s\n', 'case', 'k l_f', 'kR', '|dlnk/dlnD|', 'R/l_f');
for c = 1:numel(cases)
    s = res{c};
    fprintf('  %-36s %8.2f %8.2f %12.2f %10.2f\n', cases(c).name, logOverlap(s.kLf, b.kLf), logOverlap(s.kR, b.kR), ...
        logOverlap(abs(s.dlnKappa_dlnD), abs(b.dlnKappa_dlnD)), s.RoverLf / b.RoverLf);
end
fprintf('  (last column: R/l_f relative to p0''s %.2f)\n\n', b.RoverLf);

out = struct('cases', cases, 'results', {res}, 'stamp', stamp);
save([baseName, '.mat'], 'out');
fprintf('Saved: %s.mat\n', baseName);

%% ---- Figure ----------------------------------------------------------------------------------------------
cols = {[0 0 0], [120 120 120] / 255, [42 120 214] / 255, [235 104 52] / 255, [27 175 122] / 255, [160 90 190] / 255};
fig = figure('Color', 'w', 'Position', [60 60 1350 400]);
yq = {@(s) s.kR, @(s) s.kappaRatio, @(s) abs(s.dlnKappa_dlnD)};
yl = {'kR', '\kappa / k', '|d ln \kappa / d ln D|'};
for q = 1:3
    ax = subplot(1, 3, q); hold(ax, 'on');
    h = zeros(1, numel(cases));
    for c = 1:numel(cases)
        s = res{c};
        lw = 1.6 + 1.2 * (c == 1);
        h(c) = plot(ax, s.kLf, yq{q}(s), '-', 'Color', cols{c}, 'LineWidth', lw);
    end
    set(ax, 'XScale', 'log', 'FontSize', 9, 'Box', 'off'); grid(ax, 'on');
    if q ~= 2, set(ax, 'YScale', 'log'); end
    xlabel(ax, '$k\ell_f$', 'Interpreter', 'latex'); ylabel(ax, yl{q});
    if q == 1, legend(ax, h, {cases.name}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8); end
end
saveFigure(fig, [baseName, '.png']);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');

%% ================================================================================================
function s = yesno(tf)
if tf, s = 'yes'; else, s = 'no'; end
end

function f = logOverlap(x, ref)
% Fraction of the case's log-range [min x, max x] that lies inside p0's log-range.
% Independent of how the period band was sampled (p0 is uniform in omega, MIZ in T).
lo = max(min(x), min(ref)); hi = min(max(x), max(ref));
f = max(0, log(hi) - log(lo)) / (log(max(x)) - log(min(x)));
end

function saveFigure(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 200);
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end