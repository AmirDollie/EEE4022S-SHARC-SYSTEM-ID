%% checkAlphaScaling.m
% Impact of the incident-amplitude fix (computeSensorFRF, 2026-10-08) on the E1 headline numbers.
%
% The Forward Model returns eta for an incident wave of amplitude 1/alpha, so every cached
% acceleration FRF is low by a factor alpha(omega) = depth omega^2 / g (1.7 to 13.85 in band).
% This script evaluates the production chain (evaluateScenario) at p0, Level 2C, for every noise
% case in cfg.noise, twice:
%   OLD  the FRF the studies used (cached EMM, unscaled)
%   NEW  the same FRF multiplied by alpha(omega), passed in as an FRF handle
% The correction is an exact per-frequency scalar, so NEW needs no new EMM solves: it reuses the
% cached nodes. Nothing is written to twinCache by the NEW path.
%
% RUN THIS BEFORE patching computeSensorFRF.m or clearing twinCache. It refuses to run if the
% patched computeSensorFRF is already on the path (OLD would then not be the old FRF).
%
% Runtime: if the p0 Level 2C FRF and its four Jacobian neighbours are already in twinCache (E1
% full run), about a minute. Otherwise the first noise case solves about 150 EMM frequencies at
% [50 10 10] (several minutes); every later case and every NEW evaluation is cache-only.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();

% ---- guard: must run against the UNPATCHED computeSensorFRF ----
ws = warning('off', 'computeSensorFRF:outsideValidatedAlpha');
[~, gi] = computeSensorFRF(6.0, [0 0], cfg.p0.vec, 'Truncation', [14 4 3]);
warning(ws);
if isfield(gi, 'HdeflectionRaw')
    error('checkAlphaScaling:patched', ['computeSensorFRF is already patched. Restore the old version ' ...
        '(git stash / checkout) and keep the old twinCache, run this, then apply the patch.']);
end

g = cfg.physics.g; depth = cfg.physics.Hbasin;
alphaOf = @(w) depth * (w(:).') .^ 2 / g;
p0 = cfg.p0.vec;
sens = cfg.layout.level2C.frac .* [p0(3) 1];
frfNew = @(w, pp) correctedFRF(w, pp, p0, sens, cfg, P.cache, alphaOf);

stamp = datestr(now, 'yyyymmdd_HHMMSS');
outDir = fullfile(P.results, 'E1'); if ~exist(outDir, 'dir'), mkdir(outDir); end
logFile = fullfile(outDir, sprintf('alphaScalingCheck_%s.txt', stamp));
diary(logFile);
fprintf('checkAlphaScaling  %s\n', stamp);
fprintf('p0 = [%.4e %.4e %.4f], Level 2C, ref s26, L = %d, N = %d, sea Hs %.3f m, wp %.2f rad/s\n', ...
    p0, cfg.welch.L, cfg.acq.Nref, cfg.sea.Hs, cfg.sea.omegaP);
fprintf('alpha over the band: %.2f to %.2f (amplitude factor applied in NEW)\n\n', ...
    alphaOf(2.9784), alphaOf(8.5012));
fprintf('%-14s %-4s %10s %10s %8s %8s %8s %8s  %s\n', 'noise', 'FRF', 'sig lnB %', 'sig lnR %', ...
    'dW', 'dTot', 'minCoh', 'invalid', 'class');

names = cfg.noise.names;
res = struct('noise', {}, 'old', {}, 'new', {});
for i = 1:numel(names)
    nm = names{i};
    tc = tic;
    old = evaluateScenario(struct('Snn', cfg.noise.(nm), 'cacheDir', P.cache), cfg);
    new = evaluateScenario(struct('Snn', cfg.noise.(nm), 'frf', frfNew), cfg);
    row(nm, 'old', old); row('', 'new', new);
    fprintf('%-14s %-4s %10.3f %10.3f   (NEW/OLD sigma ratio; %.0f s)\n\n', '', 'ratio', ...
        new.sigma(1) / old.sigma(1), new.sigma(2) / old.sigma(2), toc(tc));
    res(end + 1) = struct('noise', nm, 'old', strip(old), 'new', strip(new)); %#ok<AGROW>
end
save(fullfile(outDir, sprintf('alphaScalingCheck_%s.mat', stamp)), 'res', 'stamp');
diary off;
fprintf('Saved %s\n', logFile);

% ------------------------------------------------------------------------------------------------
function H = correctedFRF(w, pp, p0, sens, cfg, cacheDir, alphaOf)
% Old cached node FRF at pp (same key emmField uses: checkPoints 3 at p0, 0 at the neighbours),
% times alpha(omega). Frequencies are clamped to the node band (the handle band in emmField is
% rounded to 4 decimals and can sit 1e-5 rad/s outside the validated band).
if isequal(pp, p0), cp = 3; else, cp = 0; end
spec = struct('p', pp, 'sensors', sens, 'nodeSpacing', cfg.emm.nodeSpacing, 'frfOptions', {{}}, ...
    'band', [], 'checkPoints', cp);
tw = synthesiseTwinRecords(spec, 'CacheDir', cacheDir);
w = w(:).';
w = min(max(w, tw.omegaNodes(1)), tw.omegaNodes(end));
H = alphaOf(w) .* tw.Hacc(w);
end

function row(nm, tag, o)
fprintf('%-14s %-4s %10.3f %10.3f %8.3f %8.3f %8.3f %8.3f  %s\n', nm, tag, 100 * o.sigma(1), ...
    100 * o.sigma(2), o.dW, o.dTot, o.minCoherence, o.invalidFraction, o.class);
end

function s = strip(o)
s = struct('sigma', o.sigma, 'dW', o.dW, 'dTot', o.dTot, 'minCoherence', o.minCoherence, ...
    'invalidFraction', o.invalidFraction, 'class', o.class, 'condF', o.condF);
end