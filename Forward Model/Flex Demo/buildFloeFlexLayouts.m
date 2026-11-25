function buildFloeFlexLayouts()
%BUILDFLOEFLEXLAYOUTS  Precompute the sensor-layout comparison shown in floeFlexDemo (run once).
%
%   Uses the D1 machinery unchanged (d1LayoutTools: the 33-point emmField from the twin cache, the
%   per-reference precompute, evaluateLayout), at p0, LSM6DSV16X noise, the production record and Welch
%   settings, heading 0 (waves along the x axis). Each layout uses its A-best reference, as D1 does.
%
%   LAYOUTS
%     Fisher-optimal pair     cfg.M1.D1.layout (the D1 n_s = 2 optimum, s6 + s18)
%     Gramian-optimal pair    recomputed by d1LayoutTools.gramianLayouts (Level 2C lambda_min optimum)
%     Mirrored pair           s20 + s24 (0.7 R at +-90 deg): symmetric about the wave axis, so both
%                             sensors move identically, T = 1 at every frequency, F singular
%     Level 2C                the thesis reference layout (s2 s6 s14 s26, reference s26)
%     Direction-robust cross  s2 s6 s20 s24 (D3a; its registered worst case over 24 headings is quoted
%                             from the thesis, the value here is at heading 0)
%
%   Writes floeFlexLayouts.mat next to this file. Takes about a minute with the D1 field cached
%   (Results/twinCache); without the cache emmField recomputes it (several minutes).
%   Lives in Forward Model/Flex Demo/.

here = fileparts(mfilename('fullpath'));
dr = fullfile(here, '..', '..', 'Inverse Mapping', 'Parameter Inversion', 'Design and Robustness');
addpath(dr);
P = setupStudyPaths();
cfg = defineStudyScenarios();
T = d1LayoutTools();
p0 = cfg.p0.vec;
t0 = tic;

geom = T.buildGeometry(cfg, p0);
full = T.buildFullField(cfg, p0, geom, P.cache, false);
Snn = cfg.noise.(cfg.noise.default);
pre = T.precompute(full, Snn, cfg);
G = T.gramianLayouts(cfg, 2);

defs = { ...   % name, grid indices, note, short label, reference ([] = A-best over the layout, as D1)
    'Fisher-optimal pair', cfg.M1.D1.layout, ...
    'D1 optimum: two sensors on the wave axis, chosen for parameter information', 'Fisher pair', cfg.M1.D1.ref; ...
    'Gramian-optimal pair', G(1).layout, ...
    'Modal-observability optimum: rim sensors, where the motion is largest', 'Gramian pair', []; ...
    'Mirrored pair (poor)', [20 24], ...
    'Symmetric about the wave axis: both sensors move identically, so the ratio carries no information', ...
    'Mirrored pair', []; ...
    'Level 2C (thesis reference)', geom.level2C, ...
    'Four sensors on the wave axis; the layout used throughout the thesis', 'Level 2C', geom.level2CRef; ...
    'Direction-robust cross', [2 6 20 24], ...
    'D3a: robust to an unknown heading (worst case over 24 headings 0.069%)', 'Cross (4)', []};

L = struct('name', {}, 'short', {}, 'idx', {}, 'ref', {}, 'sigma', {}, 'dW', {}, 'class', {}, 'singular', {}, ...
    'C', {}, 'A', {}, 'note', {});
for k = 1:size(defs, 1)
    lay = sort(defs{k, 2});
    best = [];
    refs = defs{k, 5}; if isempty(refs), refs = lay; end   % thesis reference where it has one, else A-best
    for r = refs
        q = T.evaluateLayout(pre, lay, r, cfg, 'direct');
        if isempty(best) || q.A < best.A, best = q; end
    end
    best.noInfo = best.singular || max(best.sigma) > 0.1;   % > 10%: numerically singular (mirrored pair)
    if best.noInfo, C = NaN(2); else, C = best.F \ eye(2); end
    L(end + 1) = struct('name', defs{k, 1}, 'short', defs{k, 4}, 'idx', lay, 'ref', best.ref, 'sigma', best.sigma(:).', ...
        'dW', best.dW, 'class', best.class, 'singular', best.noInfo, 'C', C, 'A', best.A, ...
        'note', defs{k, 3}); %#ok<AGROW>
    if best.noInfo
        fprintf('  %-28s %-18s ref s%-2d  no information (F numerically singular)\n', defs{k, 1}, ...
            strjoin(arrayfun(@(i) sprintf('s%d', i), lay, 'UniformOutput', false), ' '), best.ref);
    else
        fprintf('  %-28s %-18s ref s%-2d  sigma [%8.4f%% %8.4f%%]  dW %6.3f  %s\n', defs{k, 1}, ...
            strjoin(arrayfun(@(i) sprintf('s%d', i), lay, 'UniformOutput', false), ' '), best.ref, ...
            100 * best.sigma, best.dW, best.class);
    end
end
if numel(G(1).tied) > 0 && size(G(1).tied, 1) > 1
    fprintf('  note: %d tied Gramian optima; the first is shown.\n', size(G(1).tied, 1));
end
Aopt = L(1).A;
for k = 1:numel(L), L(k).recordFactor = ternary(L(k).singular, Inf, L(k).A / Aopt); end   % record length for the Fisher pair's precision

layouts = struct('L', {L}, 'gridFrac', cfg.layout.gridFrac, 'p0', p0, 'noise', cfg.noise.default, ...
    'recordMin', cfg.acq.Nref * cfg.acq.dt / 60, 'segment', cfg.welch.L, 'built', datestr(now)); %#ok<TNOW1,DATST>
save(fullfile(here, 'floeFlexLayouts.mat'), 'layouts');
fprintf('Saved floeFlexLayouts.mat (%.0f s).\n', toc(t0));
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end