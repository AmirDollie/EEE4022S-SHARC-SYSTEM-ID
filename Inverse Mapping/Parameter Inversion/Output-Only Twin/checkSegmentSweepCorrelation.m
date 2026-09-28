%% checkSegmentSweepCorrelation.m
% Post-processing of the newest S4 result (runTransmissibilityEstimatorTest):
% inter-frequency correlation of the transmissibility features for every
% segment length in the sweep, on the full Welch grid and after keeping
% every second bin. Decides the U1/U2 data representation (L, decimation).
%
% Uses out.Zsw{iL, i} (6 x nF x N_MC raw draws) for the sweep cases
% (W-clean, J-clean, J-noise). J-noise is UNCORRECTED: the corrected case
% was not part of the sweep, but at L = 512 the correction left the lag
% structure unchanged (lag 1: 0.38 both), so J-noise is the proxy.
%
% Printed per L and case: median over bins of mean|corr(z_k, z_{k+l})|,
% l = 1..3, full grid and decimated grid (lag 1 decimated = lag 2 full),
% against the Monte Carlo floor sqrt(2 / (pi N_MC)). Also the median
% relative bias and std on the decimated grid, and the feature band edges.
% No EMM, runs in seconds.
%
% Lives in Inverse Mapping/Parameter Inversion/SSI Twin/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..', '..', 'Transmissibility'));
resultsDir = fullfile(thisDir, '..', 'Results');
d = dir(fullfile(resultsDir, 'transmissibilityEstimatorTest_*.mat'));
if isempty(d), error('checkSegmentSweepCorrelation:noInput', 'No S4 result in %s.', resultsDir); end
[~, newest] = max([d.datenum]);
S = load(fullfile(resultsDir, d(newest).name)); out = S.out;
if ~isfield(out, 'Zsw'), error('checkSegmentSweepCorrelation:noRaw', 'This S4 result has no raw sweep draws (Zsw).'); end
nMC = size(out.Zsw{1, 1}, 3);
floorRho = sqrt(2 / (pi * nMC));
fprintf('checkSegmentSweepCorrelation: %s (N = %d, N_MC = %d, floor %.2f)\n\n', d(newest).name, out.NList(end), nMC, floorRho);
fprintf('%6s %-8s %5s | %-17s | %-17s | %s\n', 'L', 'case', 'bins', 'full lag 1  2  3', 'decim. lag 1  2', 'decimated: bins, med rel std');

res = struct();
for iL = 1:numel(out.segSweep)
    for i = 1:numel(out.sweepCases)
        Zc = out.Zsw{iL, i};
        full = lagCorr(Zc, 3);
        dec = lagCorr(Zc(:, 1:2:end, :), 2);
        zd = Zc(:, 1:2:end, :);
        mu = mean(zd, 3);
        sd = sqrt(sum(var(zd, 0, 3), 1)) ./ sqrt(sum(mu.^2, 1));
        fprintf('%6d %-8s %5d | %5.2f %5.2f %5.2f | %5.2f %5.2f     | %3d, %.2e\n', out.segSweep(iL), out.sweepCases{i}, ...
            size(Zc, 2), full, dec, size(zd, 2), median(sd));
        res(iL, i).full = full; res(iL, i).dec = dec; res(iL, i).relStdDec = sd; %#ok<SAGROW>
    end
end
w = out.omegaSw{end};
fprintf('\nL = %d grid: %d bins, spacing %.4f rad/s, first/last %.3f / %.3f rad/s\n', out.segSweep(end), numel(w), ...
    w(2) - w(1), w(1), w(end));
TWIN_BAND = [2.9784 8.5012]; TAPER_WIDTH = 0.2;               % S3 defaults used by S4
ramp = [TWIN_BAND(1), TWIN_BAND(1) + TAPER_WIDTH; TWIN_BAND(2) - TAPER_WIDTH, TWIN_BAND(2)];
fprintf('Synthesis taper ramps (twin artefact): [%.3f %.3f] and [%.3f %.3f] rad/s\n', ramp(1, :), ramp(2, :));
wd = w(1:2:end);
fprintf('Decimated L = %d bins inside the ramps: %d of %d (exclude them from the inversion band)\n', ...
    out.segSweep(end), sum(wd < ramp(1, 2) | wd > ramp(2, 1)), numel(wd));

function rho = lagCorr(Zc, maxLag)
[p, nF, ~] = size(Zc);
rho = NaN(1, maxLag);
for l = 1:maxLag
    r = NaN(1, nF - l);
    for k = 1:nF - l
        A = squeeze(Zc(:, k, :)); B = squeeze(Zc(:, k + l, :));
        ok = all(isfinite(A), 1) & all(isfinite(B), 1);
        if sum(ok) < 20, continue; end
        A = A(:, ok) - mean(A(:, ok), 2); B = B(:, ok) - mean(B(:, ok), 2);
        r(k) = mean(abs(sum(A .* B, 2) ./ sqrt(sum(A.^2, 2) .* sum(B.^2, 2))));
    end
    rho(l) = median(r(isfinite(r)));
end
if p < 1, rho = []; end
end