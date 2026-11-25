%% diagnoseTruncationF1.m  --  F1 (Gate G1): truncation error of candidate production models in Fisher units
%
% WHY. The thesis V1a compared [50 10 10] with [210 20 20] and found d_sys = 10-29, with no converged reference.
% The F1 diagnostics (diagnostics/diag_F1_*.m, 25 Nov 2026) separated the error into three independent parts:
%   open-water tower Mo   staircase, period 1/gamma, envelope ~ Mo^-1.4    -> Mo = 20000-40000, nearly free
%   plate tower M         clean power law, p ~ 1.5-1.8, settled by M = 800  -> M = 400-800 + Richardson
%   angular orders N      N = 10 leaves 2% in T at R = 2 R0                -> N per floe (angularOrderFor)
% and computeSensorFRF now takes 'Truncation' [M P N Mo]. This script turns that into the numbers Gate G1 needs:
% for each candidate truncation, the feature error against a converged reference, projected through the
% scenario's own F, J, Sigma (the V1 ruler, as diagnoseV1Truncation):
%       d = sqrt(dth' F dth),   dth = F^-1 J' Sigma^-1 b,   b = f(candidate) - f_inf
% plus dbeta/beta, dR/R and d_tot with the Welch bias (out.bW + b).
%
% REFERENCE f_inf. Plate ladder [M 10 N* Mo_ref], M = 100, 200, 400, 800, with N* = angularOrderFor(alpha_top, R)
% and Mo_ref = 40000. Rates p from three consecutive levels; Richardson extrapolates from (200, 400) and (400, 800);
% f_inf is the latter and the d of the difference between the two is its UNCERTAINTY. Exact (point) features, as
% diagnoseV1Truncation's preview; Welch-minus-Welch differs only at the window-smoothing level.
%
% CANDIDATES (per floe):  [50 10 10] (thesis production, = out.f0), [50 10 N*], [50 10 N* 20000],
%   [200 10 N* 20000], [400 10 N* 20000], [800 10 N* 40000]
% CHECKS:  N  -> [400 10 N*+5 20000] vs [400 10 N* 20000];   Mo -> [400 10 N* 40000] vs [400 10 N* 20000]
%
% FLOES. basin_p0 and cfg.V1.truncationFloes (the large-kR A1 floes, R = 2 R0, where N matters).
% RULER CAVEAT. F, J, Sigma are the scenario's own at [50 10 10]; for the R = 2 R0 floes that model carries the N
% error, so the ruler is approximate there (as it was in V1). Distances of a few tenths are what matter.
%
% COST. One nominal field (29 nodes, no derivatives) per truncation per floe; M = 800 takes ~6-13 s per node, so
% about 10 min per floe, about 40 min in all. Every field is cached in Results/twinCache, so a rerun after an
% interruption repeats only what was running.  Quick mode (F1D_QUICK = true): p0 only, M <= 200, small Mo.
% OUTPUT. Printed, plus Results/F1/truncationF1_<stamp>.log and .mat.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.   (F1, 2026-11-25)

P = setupStudyPaths();
cfg = defineStudyScenarios();
A = a1Tools();
QUICK = exist('F1D_QUICK', 'var') && F1D_QUICK;
dDir = fullfile(P.results, 'F1');
if ~exist(dDir, 'dir'), mkdir(dDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS');
logFile = fullfile(dDir, sprintf('truncationF1_%s.log', stamp));
diary(logFile); diary('on');
fprintf('\n==== F1: truncation error of candidate models in Fisher units (%s)%s ====\n', stamp, ...
    ternaryStr(QUICK, '  [QUICK MODE]', ''));

try
aTop = cfg.emm.alphaValidated(2);
Mlad = [100 200 400 800]; MoRef = 40000; MoStd = 20000;
names = [{'basin_p0'}, cfg.V1.truncationFloes(:).'];
if QUICK, Mlad = [50 100 200]; MoRef = 5000; MoStd = 2000; names = {'basin_p0'}; end
% optional overrides (set in the workspace before running; clear them afterwards), e.g. to tighten the reference
% at p0 with one more level:   F1D_MLAD = [200 400 800 1600]; F1D_FLOES = {'basin_p0'};
if exist('F1D_MLAD', 'var'), Mlad = F1D_MLAD; end
if exist('F1D_FLOES', 'var'), names = F1D_FLOES; end
fprintf('  reference ladder M = %s at Mo = %d; candidates use Mo = %d; N* = angularOrderFor(alpha = %.2f, R)\n', ...
    mat2str(Mlad), MoRef, MoStd, aTop);

sc = A.buildScenarios(cfg);
opts = struct('cacheDir', P.cache, 'verbose', false);
lab = @(t) mat2str(t);
F1 = struct('floe', names, 'Nstar', NaN, 'ladder', [], 'p', [], 'unc', NaN, 'cands', [], 'dN', NaN, 'dMo', NaN, ...
    'dW', NaN, 'sigma', []);

for i = 1:numel(names)
    s = sc(strcmp({sc.name}, names{i}));
    if numel(s) ~= 1, error('diagnoseTruncationF1:floe', 'A1 scenario %s not found exactly once.', names{i}); end
    scn = A.buildScn(s, cfg, opts);
    out = evaluateScenario(scn, cfg);                                   % the ruler (cached A1 fields)
    Ns = angularOrderFor(aTop, scn.p(3));
    feat = @(tr) A.features(s, scn, out.omega, cfg.emm.nodeSpacing, [s.frfOptions, {'Truncation', tr}], opts);
    dist = @(b) A.projectBias(out, b);
    kTop = dispersionRoots(aTop, 0, 0, 1);                              % [k0; k1]; k0 = i p0
    fprintf('\n---- %s:  R = %.4f, p0 R at the band top = %.2f, N* = %d;  [50 10 10]: sigma_lnbeta %.4f%%, d_W %.3f, %s\n', ...
        names{i}, scn.p(3), imag(kTop(1)) * scn.p(3), Ns, 100 * out.sigma(1), out.dW, out.class);

    % ---- reference ladder
    nL = numel(Mlad); fL = cell(1, nL); tL = zeros(1, nL);
    for k = 1:nL
        tc = tic; fL{k} = feat([Mlad(k) 10 Ns MoRef]); tL(k) = toc(tc);
        fprintf('  ladder %-22s %7.0f s\n', lab([Mlad(k) 10 Ns MoRef]), tL(k));
        if exist('OCTAVE_VERSION', 'builtin'), fflush(stdout); else, drawnow; end
    end
    sStep = zeros(1, nL - 1);
    for k = 1:nL - 1, sStep(k) = norm(fL{k + 1} - fL{k}) / norm(fL{end}); end
    pL = log2(sStep(1:end - 1) ./ sStep(2:end));
    ext = cell(1, numel(pL));
    for j = 1:numel(pL), ext{j} = fL{j + 2} + (fL{j + 2} - fL{j + 1}) / (2^pL(j) - 1); end
    fInf = ext{end};
    if numel(ext) >= 2, qU = dist(ext{end} - ext{end - 1}); unc = qU.d; else, unc = NaN; end
    fprintf('  steps |df|/|f|: %s   rates p: %s   reference uncertainty d = %.3f\n', ...
        mat2str(sStep, 3), mat2str(pL, 3), unc);

    % ---- candidates against f_inf
    cand = {[50 10 10], [50 10 Ns], [50 10 Ns MoStd], [200 10 Ns MoStd], [400 10 Ns MoStd], [800 10 Ns MoRef]};
    if QUICK, cand = {[50 10 10], [50 10 Ns], [50 10 Ns MoStd], [100 10 Ns MoStd], [200 10 Ns MoRef]}; end
    C = struct('tr', cand, 'd', NaN, 'betaPct', NaN, 'RPct', NaN, 'dTot', NaN, 'relF', NaN, 'seconds', NaN);
    fprintf('  %-22s %9s %11s %11s %9s %10s\n', 'candidate', 'd_trunc', 'dbeta/beta', 'dR/R', 'd_tot', '|df|/|f|');
    for c = 1:numel(cand)
        tr = cand{c}; tc = tic;
        if isequal(tr, [50 10 10])
            f = out.f0;
        else
            k = find(cellfun(@(t) isequal(t, tr), arrayfun(@(m) [m 10 Ns MoRef], Mlad, 'UniformOutput', false)), 1);
            if ~isempty(k), f = fL{k}; else, f = feat(tr); end
        end
        b = f(:) - fInf(:);
        q = dist(b); qt = dist(out.bW(:) + b);
        C(c).d = q.d; C(c).betaPct = 100 * (exp(q.dtheta(1)) - 1); C(c).RPct = 100 * (exp(q.dtheta(2)) - 1);
        C(c).dTot = qt.d; C(c).relF = norm(b) / norm(fInf); C(c).seconds = toc(tc);
        fprintf('  %-22s %9.3f %+10.4f%% %+10.5f%% %9.3f %10.2e\n', lab(tr), C(c).d, C(c).betaPct, C(c).RPct, ...
            C(c).dTot, C(c).relF);
    end

    % ---- N and Mo checks at M = 400 (or the largest quick level)
    Mc = Mlad(end - 1);
    fBase = feat([Mc 10 Ns MoStd]);
    qN = dist(feat([Mc 10 Ns + 5 MoStd]) - fBase);
    qMo = dist(fL{nL - 1} - fBase);
    fprintf('  N check  %s -> %s: d = %.4f\n', lab([Mc 10 Ns MoStd]), lab([Mc 10 Ns + 5 MoStd]), qN.d);
    fprintf('  Mo check %s -> %s: d = %.4f\n', lab([Mc 10 Ns MoStd]), lab([Mc 10 Ns MoRef]), qMo.d);

    F1(i).Nstar = Ns; F1(i).ladder = struct('M', Mlad, 'steps', sStep, 'seconds', tL); F1(i).p = pL; F1(i).unc = unc;
    F1(i).cands = rmfield(C, 'seconds'); F1(i).dN = qN.d; F1(i).dMo = qMo.d; F1(i).dW = out.dW; F1(i).sigma = out.sigma;
    save(fullfile(dDir, sprintf('truncationF1_%s.mat', stamp)), 'F1', 'cfg', '-v7');   % partial result, every floe
end

%% ---- summary ---------------------------------------------------------------------------------------------------------
fprintf('\nSUMMARY  d_trunc against the reference (reference uncertainty in brackets); G1 target d <= 0.1\n');
fprintf('  %-26s', 'candidate'); fprintf('%26s', names{:}); fprintf('\n');
for c = 1:numel(F1(1).cands)
    fprintf('  %-26s', lab(F1(1).cands(c).tr));
    for i = 1:numel(F1)
        if c <= numel(F1(i).cands), fprintf('%16.3f (%6.3f) ', F1(i).cands(c).d, F1(i).unc); end
    end
    fprintf('\n');
end
fprintf('  (for floes other than p0 the N entries are each floe''s own N*: see the per-floe tables)\n');
fprintf('\n  saved %s\n', fullfile(dDir, sprintf('truncationF1_%s.mat', stamp)));
diary('off');
catch err
    diary('off');
    rethrow(err);
end

function s = ternaryStr(c, a, b)
if c, s = a; else, s = b; end
end
