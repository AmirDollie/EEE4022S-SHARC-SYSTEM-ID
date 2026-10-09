function cases = resolveM1Scenarios(cfg, P, readSources)
%RESOLVEM1SCENARIOS  The four M1 cases, rebuilt by the code of the studies that chose them.
%   cases = resolveM1Scenarios(cfg, P)          scenarios AND the source-study values (reads Results/)
%   cases = resolveM1Scenarios(cfg, P, false)   scenarios only (no results files needed)
%
%   Provenance only: no Monte Carlo and no evaluation. Each scenario is built by its source study's own
%   construction code, never reconstructed by hand, so the runner and the tester resolve them identically:
%     1 D1_ns2            D1 n_s = 2 optimum: grid points cfg.M1.D1.layout, reference cfg.M1.D1.ref, on the
%                         D1 33-point grid at p0 (d1LayoutTools.buildGeometry; the cached D1 field), LSM6DSV16X
%     2 A1_<A1scenario>   a1Tools.buildScenarios + buildScn: the A1 scenario cfg.M1.A1scenario exactly as A1
%                         built it (its band, nu, Level 2C at its true R, A1 noise case)
%     3 R2_<R2sea>        buildR2Seas: the R2-normalised sea cfg.M1.R2sea (fixed in-band variance), L = cfg.welch.L
%     4 R1_p0_lsm6dsv16x  the E1 default (p0, Level 2C, s26) at LSM6DSV16X with the R1 sea
%
%   cases(q): index, name, rule (cfg.M1.rules{q}), study, scn (the evaluateScenario input), note, and
%   source: file (newest source-study .mat), sigma [ln beta, ln R], dW, class, plus study-specific
%   identity fields (D1 layout/ref; A1 node spacing used; R2 selected sea). The source values are what
%   the M1 analytic target must reproduce (cfg.M1.sourceRelTol) before any record is generated.
%   An identity mismatch with the source study (e.g. D1's optimum is not s6 + s18) is an error.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
if nargin < 3, readSources = true; end
names = cfg.M1.caseNames;
cases = struct('index', num2cell(1:numel(names)), 'name', names, 'rule', cfg.M1.rules, 'study', '', ...
    'scn', [], 'note', '', 'source', []);

%% ---- 1  D1 n_s = 2 optimum -------------------------------------------------------------------------------------
D = cfg.M1.D1;
T = d1LayoutTools();
geom = T.buildGeometry(cfg, cfg.p0.vec);
lay = sort(D.layout);
cases(1).study = 'D1';
cases(1).scn = struct('p', cfg.p0.vec, 'sensors', geom.sensorsND, 'layout', lay, 'ref', find(lay == D.ref), ...
    'Snn', cfg.noise.(D.noiseCase), 'cacheDir', P.cache);
cases(1).note = sprintf('%s with reference %s on the D1 grid (33 points; the cached D1 field), %s', ...
    strjoin(geom.names(lay), ' + '), geom.names{D.ref}, D.noiseCase);

%% ---- 2  A1 boundary case ---------------------------------------------------------------------------------------
A = a1Tools();
sc = A.buildScenarios(cfg);
i = find(strcmp({sc.name}, cfg.M1.A1scenario));
if numel(i) ~= 1, error('resolveM1Scenarios:A1', 'A1 scenario %s not found exactly once.', cfg.M1.A1scenario); end
cases(2).study = 'A1';
cases(2).scn = A.buildScn(sc(i), cfg, struct('cacheDir', P.cache));
cases(2).note = sprintf('%s (%s), p = [%.4g %.4g %.4g], band %.3f to %.3f rad/s, scale %.4g', sc(i).name, ...
    sc(i).group, sc(i).p, sc(i).band, sc(i).scale);

%% ---- 3  R2 stress sea ------------------------------------------------------------------------------------------
[seas, mkScn] = buildR2Seas(cfg);
j = find(strcmp({seas.name}, cfg.M1.R2sea));
if numel(j) ~= 1, error('resolveM1Scenarios:R2', 'R2 sea %s not found exactly once.', cfg.M1.R2sea); end
cases(3).study = 'R2';
cases(3).scn = mkScn(seas(j), cfg.welch.L, P.cache);
cases(3).note = sprintf('%s: JONSWAP omega_p %g, gamma_J %g, Hs %.5f m (fixed in-band variance), L = %d', ...
    seas(j).name, seas(j).omegaP, seas(j).gammaJ, seas(j).Hs, cfg.welch.L);

%% ---- 4  p0 at the measured LSM6DSV16X noise ---------------------------------------------------------------------
cases(4).study = 'R1';
cases(4).scn = struct('Snn', cfg.noise.(cfg.M1.R1noiseCase), 'sea', cfg.R1.sea, 'cacheDir', P.cache);
cases(4).note = sprintf('E1 default (p0, Level 2C, reference s26) at %s, R1 sea', cfg.M1.R1noiseCase);

if ~readSources, return; end

%% ---- source-study values ---------------------------------------------------------------------------------------
% D1: best A-optimal admissible layout at n_s = 2, headline noise, primary spacing
[S, f] = newest(P, 'D1');
c = find(strcmp(S.results.noise, D.noiseCase), 1); k = find(S.results.ns == D.ns, 1);
b = S.results.best(c, k, D.spacingIndex).A;
if ~isequal(sort(b.layout), lay) || b.ref ~= D.ref
    error('resolveM1Scenarios:D1identity', 'D1 optimum in %s is [%s] ref %d, not the pre-registered [%s] ref %d.', ...
        f, num2str(b.layout), b.ref, num2str(lay), D.ref);
end
cases(1).source = struct('file', f, 'sigma', [b.sB b.sR], 'dW', b.dW, 'class', cfg.class.names{b.classCode}, ...
    'layout', b.layout, 'ref', b.ref);

% A1: the scenario's production metrics, and the node grid A1 actually used
[S, f] = newest(P, 'A1');
r = S.results.res{strcmp({S.results.scenarios.name}, cfg.M1.A1scenario)};
if r.nodeSpacingUsed ~= cases(2).scn.nodeSpacing
    error('resolveM1Scenarios:A1grid', 'A1 evaluated %s on a %.2g rad/s node grid, M1 resolves %.2g.', ...
        cfg.M1.A1scenario, r.nodeSpacingUsed, cases(2).scn.nodeSpacing);
end
cases(2).source = struct('file', f, 'sigma', r.sigma(1:2), 'dW', r.dW, 'class', r.finalClass, ...
    'nodeSpacingUsed', r.nodeSpacingUsed);

% R2: the selected M1 sea at L = 2048
[S, f] = newest(P, 'R2');
if ~strcmp(S.results.m1.selected, cfg.M1.R2sea)
    error('resolveM1Scenarios:R2identity', 'R2 (%s) selected %s, not the pre-registered %s.', f, ...
        S.results.m1.selected, cfg.M1.R2sea);
end
q = find(strcmp({S.results.seas.name}, cfg.M1.R2sea)); iL = S.results.iL;
cases(3).source = struct('file', f, 'sigma', reshape(S.results.metrics.sigma(q, iL, 1:2), 1, 2), ...
    'dW', S.results.metrics.dW(q, iL), 'class', S.results.metrics.class{q, iL}, 'selected', S.results.m1.selected);

% R1: the exact evaluation of the named noise case
[S, f] = newest(P, 'R1');
m = find(strcmp(S.results.names, cfg.M1.R1noiseCase));
cases(4).source = struct('file', f, 'sigma', [S.results.marked.sB(m) S.results.marked.sR(m)], ...
    'dW', S.results.marked.dW(m), 'class', S.results.marked.class{m});
end

%% ================================================================================================
function [S, file] = newest(P, id)
L = dir(fullfile(P.results, id, sprintf('%s_20*.mat', id)));
if isempty(L), error('resolveM1Scenarios:noSource', 'No %s results in %s: M1 needs the source study.', id, ...
        fullfile(P.results, id)); end
[~, k] = max([L.datenum]);
file = fullfile(L(k).folder, L(k).name);
S = load(file, 'results');
end