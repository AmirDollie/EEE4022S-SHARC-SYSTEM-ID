%% diag_F1_plate.m  --  F1 (Thursday): plate-tower rate out to M = 800, angular truncation N, hard floes
%
% CONTEXT
%   diag_F1_Mo_ladder.m: with the open-water tower at Mo = 40000 the plate tower converges as a clean power
%   law, but p was still drifting (8.5 rad/s: 1.69 -> 1.76; 5.0 rad/s: 1.49 -> 1.54), so the Richardson
%   uncertainty was ~60% of the remaining error. This run uses the new split solver
%   (precomputeDeflectionDataSplit, open water at Mo, plate at M) and the overflow-safe Bessel functions.
%
%   Part 1  p0 at 8.5 and 5.0 rad/s: M = 50 ... 800 at Mo = 40000. Three rates; does p settle? Richardson
%           extrapolates from (100,200), (200,400), (400,800) and the spread between the last two.
%   Part 2  Angular truncation: N = 10, 15, 20 at p0 and 2R0, 8.5 rad/s (M = 200, Mo = 20000). The review
%           flagged p0 R ~ 10.6 > N = 10 for 2R0 at the top of the band.
%   Part 3  Hard floes at 8.5 rad/s (Mo = 20000): 2R0, beta0/4 and the (beta0/2, 2R0) corner. M ladder
%           50 ... 400: rates, the production M = 50 error, and the Richardson spread.
%
% All errors are relative norms of the three transmissibilities T = eta(s2, s6, s14) / eta(s26), sensors at
% fixed fractions of each floe's radius. Fisher units (d_trunc) come next, through the V1a pipeline.
%
% RUNTIME  About 5-10 minutes (the M = 800 solves dominate); quick mode about a minute.
% WRITES   results/diag_F1_plate.log, .mat, .png
% HOW      run('diagnostics/diag_F1_plate.m')      quick:  PLATE_QUICK = true; run(...)   (clear it after)

%% ---- setup --------------------------------------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
root = fileparts(here);
addpath(fullfile(root, 'Forward Model'));
addpath(fullfile(root, 'Forward Model', 'Animation'));
resDir = fullfile(root, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
logFile = fullfile(resDir, 'diag_F1_plate.log');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile); diary('on');
cleanupDiary = onCleanup(@() diary('off'));

QUICK = exist('PLATE_QUICK', 'var') && PLATE_QUICK;
fprintf('==== diag_F1_plate  (%s)%s ====\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ternary(QUICK, '  [QUICK MODE]', ''));

bet0 = 4.6985e-5; gam0 = 1.4548e-3; R0 = 0.3830; nu = 0.3; depth = 1.88; g = 9.81; P = 10;
sens = [0.3 0; 0.3 pi; 0.5 pi; 0.9 0];
aOf = @(w) depth * w^2 / g;

M1 = [50 100 200 400 800]; Mo1 = 40000; om1 = [8.5 5.0];
Nlist = [10 15 20]; M2 = 200; Mo2 = 20000;
M3 = [50 100 200 400]; Mo3 = 20000;
if QUICK
    M1 = [50 100 200]; Mo1 = 5000; om1 = 8.5; Nlist = [10 15]; M2 = 50; Mo2 = 2000; M3 = [50 100]; Mo3 = 2000;
end
floe = @(name, b, gm, R) struct('name', name, 'beta', b, 'gam', gm, 'R', R);
fl0 = floe('p0', bet0, gam0, R0);

%% ---- Part 1: p0 out to M = 800 --------------------------------------------------------------------------------------
fprintf('PART 1  p0, plate M ladder %s at Mo = %d\n', mat2str(M1), Mo1);
P1 = struct('omega', {}, 'M', {}, 'T', {}, 'sec', {}, 's', {}, 'p', {}, 'Tinf', {}, 'spread', {}, 'err', {});
for w = 1:numel(om1)
    T = complex(zeros(3, numel(M1))); sec = zeros(1, numel(M1));
    for i = 1:numel(M1)
        tc = tic;
        T(:, i) = transmissibility(aOf(om1(w)), fl0, nu, M1(i), P, 10, Mo1, sens);
        sec(i) = toc(tc);
    end
    [s, p, Tinf, spread, err] = ladderStats(T);
    P1(w) = struct('omega', om1(w), 'M', M1, 'T', T, 'sec', sec, 's', s, 'p', p, 'Tinf', Tinf, 'spread', spread, 'err', err);
    fprintf('\n  omega = %.1f rad/s\n  %10s %12s %8s %10s\n', om1(w), 'M -> 2M', 's_T', 'p', 'time 2M');
    for i = 1:numel(s)
        fprintf('  %4d->%-5d %12.3e %8s %9.1fs\n', M1(i), M1(i + 1), s(i), num2strOr(p, i), sec(i + 1));
    end
    if ~isempty(spread)
        fprintf('  Richardson extrapolates: successive spreads'); fprintf(' %.2e', spread);
        fprintf('  (last = uncertainty of T_inf)\n  Error against T_inf:'); fprintf('  M = %d: %.2e', [M1; err]); fprintf('\n');
    end
end
fprintf('  Reading: p settling (last two within ~0.05) and a shrinking spread mean the extrapolate can be trusted.\n\n');

%% ---- Part 2: angular truncation N ------------------------------------------------------------------------------------
fprintf('PART 2  Angular truncation at 8.5 rad/s (M = %d, Mo = %d): |T(N) - T(N_max)| / |T|\n', M2, Mo2);
fls2 = [fl0, floe('2R0', bet0, gam0, 2 * R0)];
P2 = struct('floe', {}, 'N', {}, 'T', {}, 'err', {});
a = aOf(8.5); p0R = imag(openWaterRoots(a, 0)) ;
for f = 1:numel(fls2)
    T = complex(zeros(3, numel(Nlist)));
    for i = 1:numel(Nlist), T(:, i) = transmissibility(a, fls2(f), nu, M2, P, Nlist(i), Mo2, sens); end
    err = vecnorm(T - T(:, end), 2, 1) / norm(T(:, end));
    P2(f) = struct('floe', fls2(f).name, 'N', Nlist, 'T', T, 'err', err);
    fprintf('  %-4s (p0 R = %5.2f):', fls2(f).name, p0R * fls2(f).R); fprintf('   N = %d: %.2e', [Nlist; err]); fprintf('\n');
end
fprintf('  Reading: the N = 10 entry is the error of the production angular truncation.\n\n');

%% ---- Part 3: hard floes -------------------------------------------------------------------------------------------
fprintf('PART 3  Hard floes at 8.5 rad/s, plate M ladder %s at Mo = %d\n', mat2str(M3), Mo3);
fls3 = [floe('2R0', bet0, gam0, 2 * R0), floe('beta0/4', bet0 / 4, gam0, R0), floe('corner', bet0 / 2, gam0, 2 * R0)];
P3 = struct('floe', {}, 'M', {}, 'T', {}, 's', {}, 'p', {}, 'Tinf', {}, 'spread', {}, 'err', {});
for f = 1:numel(fls3)
    T = complex(zeros(3, numel(M3)));
    for i = 1:numel(M3), T(:, i) = transmissibility(a, fls3(f), nu, M3(i), P, 10, Mo3, sens); end
    [s, p, Tinf, spread, err] = ladderStats(T);
    P3(f) = struct('floe', fls3(f).name, 'M', M3, 'T', T, 's', s, 'p', p, 'Tinf', Tinf, 'spread', spread, 'err', err);
    fprintf('  %-8s steps', fls3(f).name); fprintf(' %.2e', s); fprintf('   p'); fprintf(' %.2f', p);
    if ~isempty(spread)
        fprintf('   spread %.2e   err(M=50) %.2e   err(M=%d) %.2e', spread(end), err(1), M3(end), err(end));
    end
    fprintf('\n');
end
fprintf('\n');

%% ---- figure ---------------------------------------------------------------------------------------------------------
fig = figure('Visible', 'off', 'Position', [100 100 1100 420]);
subplot(1, 2, 1); hold on;
for w = 1:numel(P1), plot(P1(w).M(1:end - 1), P1(w).s, '-o', 'LineWidth', 1.2); end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; box on;
xlabel('M'); ylabel('|T(2M) - T(M)| / |T|'); title(sprintf('(a) p0, plate steps at Mo = %d', Mo1));
legend(arrayfun(@(r) sprintf('%.1f rad/s', r.omega), P1, 'UniformOutput', false), 'Location', 'southwest');
subplot(1, 2, 2); hold on;
for f = 1:numel(P3), plot(P3(f).M(1:end - 1), P3(f).s, '-o', 'LineWidth', 1.2); end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; box on;
xlabel('M'); ylabel('|T(2M) - T(M)| / |T|'); title(sprintf('(b) hard floes, 8.5 rad/s, Mo = %d', Mo3));
legend({P3.floe}, 'Location', 'southwest');
pngFile = fullfile(resDir, 'diag_F1_plate.png');
exportgraphics(fig, pngFile, 'Resolution', 150);
close(fig);

save(fullfile(resDir, 'diag_F1_plate.mat'), 'P1', 'P2', 'P3', 'QUICK', '-v7');
fprintf('Saved %s\n      %s\n==== done ====\n', fullfile(resDir, 'diag_F1_plate.mat'), pngFile);
diary('off');

%% ==== local functions ================================================================================================
function T = transmissibility(alpha, fl, nu, M, P, N, Mo, sens)
    data = precomputeDeflectionDataSplit(alpha, fl.beta, fl.gam, fl.R, nu, M, P, N, Mo);
    e = evaluateDeflection(data, sens(:, 1).' * fl.R, sens(:, 2).');
    e = e(:);
    T = e(1:3) / e(4);
end

function [s, p, Tinf, spread, err] = ladderStats(T)
% Doubling ladder: steps, rates from consecutive steps, Richardson extrapolates from each pair of levels
% (using the rate from the three levels ending there), their successive spread, and errors against the last.
    nL = size(T, 2);
    s = vecnorm(diff(T, 1, 2), 2, 1) / norm(T(:, end));
    p = log2(s(1:end - 1) ./ s(2:end));
    ext = complex(zeros(3, numel(p)));
    for j = 1:numel(p)                                   % levels j, j+1, j+2 -> extrapolate from (j+1, j+2)
        ext(:, j) = T(:, j + 2) + (T(:, j + 2) - T(:, j + 1)) / (2^p(j) - 1);
    end
    if size(ext, 2) >= 2
        spread = vecnorm(diff(ext, 1, 2), 2, 1) / norm(ext(:, end));
        Tinf = ext(:, end);
        err = vecnorm(T - Tinf, 2, 1) / norm(Tinf);
    else
        spread = []; Tinf = NaN(3, 1); err = NaN(1, nL);
    end
end

function s = num2strOr(p, i)
    if i <= numel(p), s = sprintf('%.2f', p(i)); else, s = ''; end
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
