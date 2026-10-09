%% testM1MonteCarloSpotChecks.m
% Pre-flight checks for study M1. All must pass before runM1MonteCarloSpotChecks is run. Uses the cached
% fields of D1, A1, R2 and R1 (seconds), one direct EMM evaluation for check 11 (about 30 s), and an
% analytic field for the Monte Carlo, checkpoint and solver mechanics (no EMM).
%
%   1  cfg.M1: four cases, pre-registered acceptance values unchanged, the clarification logged
%   2  g3EstimatorCore is G3's code: all seven function bodies identical to runTransmissibilityTwinInversion.m
%   3  scenario identities (no results files): D1 cfg.M1.D1 layout and reference on the D1 grid; A1 =
%      a1Tools.buildScn of cfg.M1.A1scenario; R2 = buildR2Seas's cfg.M1.R2sea (normalised Hs, not 0.05 m);
%      R1 = E1 default + R1 sea (names read from cfg since 2026-10-09)
%   4  every analytic target reproduces its source study (sigma, d_W, class) and the source identity
%      (D1 optimum, R2 selection, A1 node grid) matches the pre-registration
%   5  estimator path (case 4, one record): G3's segment Welch reproduces extractTransmissibility; the
%      linearised estimate and q are finite; no accepted bin rejected
%   6  records are deterministic in the seed (same seed bitwise equal, other seed different); the seed
%      schedule never repeats across cases and records, and wave and noise seeds never collide
%   7  linear algebra on a real target (case 4): z = f0 gives dtheta = 0 exactly; z = f0 + J delta gives
%      dtheta = delta; q = delta' F delta; G3's infoMatrix on the accepted bins equals the engine's F
%   8  checkpoint and resume (analytic field, 4 records, checkpoint every 2): an interrupted and resumed
%      batch is bitwise identical to an uninterrupted one; a checkpoint with another key is ignored
%   9  the borderline rule on known numbers, and summarise on the batch of check 8
%  10  nonlinear solver (analytic field): the model at the truth equals f0; from beta x 0.8 and x 1.2 the G3
%      solver converges within 0.2 sigma of the linearised estimate on the same record
%  11  the nonlinear model on the REAL EMM equals the engine's f0 for the A1 case (same node grid, nu,
%      truncation and band as the analytic target)
%  12  the provenance guard: uncommitted .m files are found in `git status --porcelain` output, other files
%      are not
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
P = setupStudyPaths();
cfg = defineStudyScenarios();
T = m1Tools(); G = g3EstimatorCore();
nPass = 0; nFail = 0;

%% ---- 1  configuration ------------------------------------------------------------------------------------------
M = cfg.M1;
ok = numel(M.caseNames) == 4 && numel(M.rules) == 4 && isequal(M.accept.empPred, [0.85 1.20]) && ...
    isequal(M.accept.coverage95, [0.88 0.99]) && M.accept.nonlinVsLinSigma == 0.2 && ...
    isequal(M.nonlinearBetaOffsets, [-0.20 0.20]) && M.nRecords == 100 && M.nRecordsExtended == 200 && ...
    any(arrayfun(@(c) ~isempty(strfind(c.field, 'M1.caseNames')), cfg.meta.changes));
[nPass, nFail] = check(ok, sprintf('1  cfg.M1: %s; acceptance [%.2f %.2f], [%.2f %.2f], %.1f sigma unchanged; clarification logged', ...
    strjoin(M.caseNames, ', '), M.accept.empPred, M.accept.coverage95, M.accept.nonlinVsLinSigma), nPass, nFail);

%% ---- 2  G3 code reused verbatim ------------------------------------------------------------------------------------
g3File = fullfile(fileparts(P.dr), 'Output-Only Twin', 'runTransmissibilityTwinInversion.m');
src = fileread(g3File); cpy = fileread(which('g3EstimatorCore'));
same = cellfun(@(n) strcmp(fnText(src, n), fnText(cpy, n)) && ~isempty(fnText(src, n)), G.names);
msg = sprintf('2  g3EstimatorCore: %d of %d functions identical to %s', nnz(same), numel(same), G.source);
if ~all(same), msg = [msg ': DIFFER ' strjoin(G.names(~same), ', ')]; end
[nPass, nFail] = check(all(same), msg, nPass, nFail);

%% ---- 3  scenario identities --------------------------------------------------------------------------------------
cs = resolveM1Scenarios(cfg, P, false);
geom = d1LayoutTools().buildGeometry(cfg, cfg.p0.vec);
s1 = cs(1).scn;
ok1 = isequal(s1.layout, sort(cfg.M1.D1.layout)) && s1.layout(s1.ref) == cfg.M1.D1.ref && isequal(s1.sensors, geom.sensorsND) && ...
    s1.Snn == cfg.noise.lsm6dsv16x && isequal(s1.p, cfg.p0.vec);
A = a1Tools(); sc = A.buildScenarios(cfg); ia = strcmp({sc.name}, cfg.M1.A1scenario);
ok2 = isequal(cs(2).scn, A.buildScn(sc(ia), cfg, struct('cacheDir', P.cache)));
[seas, mkScn] = buildR2Seas(cfg); ir = strcmp({seas.name}, cfg.M1.R2sea);
s3 = cs(3).scn;
ok3 = isequal(s3, mkScn(seas(ir), cfg.welch.L, P.cache)) && nnz(ia) == 1 && nnz(ir) == 1 && s3.sea.omegaP == seas(ir).omegaP && s3.sea.gammaJ == seas(ir).gammaJ && ...
    abs(s3.sea.Hs - seas(ir).Hs) == 0 && abs(s3.sea.Hs - cfg.sea.Hs) > 1e-6 && s3.L == 2048;
s4 = cs(4).scn;
ok4 = s4.Snn == cfg.noise.lsm6dsv16x && isequal(s4.sea, cfg.R1.sea) && isequal(cfg.R1.sea, cfg.sea) && ...
    isempty(setdiff(fieldnames(s4), {'Snn', 'sea', 'cacheDir'}));
[nPass, nFail] = check(ok1 && ok2 && ok3 && ok4, sprintf(['3  identities: D1 %s, A1 %s, R2 %s (Hs %.5f m), R1 %s ' ...
    '(%s)'], tf(ok1), tf(ok2), tf(ok3), s3.sea.Hs, tf(ok4), cs(3).note), nPass, nFail);

%% ---- 4  analytic targets reproduce the source studies ---------------------------------------------------------------
fprintf('  ... analytic targets of the four cases (cached fields)\n');
cs = resolveM1Scenarios(cfg, P);
ctx = cell(1, 4); r4 = cell(1, 4); line = {};
for q = 1:4
    ctx{q} = T.context(cs(q), cfg);
    r4{q} = T.compareSource(ctx{q}.tgt, cs(q).source, cfg);
    line{end + 1} = sprintf('%s %.1e/%.1e', cs(q).study, r4{q}.relSigma, r4{q}.dDW); %#ok<SAGROW>
end
ok = all(cellfun(@(r) r.ok, r4));
[nPass, nFail] = check(ok, sprintf('4  targets vs source studies (rel. sigma / d_W): %s; classes %s', strjoin(line, ', '), ...
    strjoin(cellfun(@(c) c.tgt.class, ctx, 'UniformOutput', false), ', ')), nPass, nFail);

%% ---- 5  estimator path, one record ---------------------------------------------------------------------------------
c4 = ctx{4}; t4 = c4.tgt;
Y = synthesiseTwinRecords(c4.twin, t4.N, t4.dt, 'Seed', T.seed(cfg, 4, 1), c4.specArgs{:}, 'NoiseStd', c4.noiseStd);
Y = Y(t4.layout, :);
Pw = G.welchSegments(Y, t4.dt, t4.L, t4.overlap, t4.kBins, t4.ref, c4.others);
z = G.zFromSegments(Pw, 1:Pw.K, t4.Snn);
ex = extractTransmissibility(Y, t4.dt, 'Reference', t4.ref, 'SegmentLength', t4.L, 'Overlap', t4.overlap, ...
    'Band', c4.scn.bandEst, 'NoisePSD', t4.Snn);
[hit, at] = ismember(t4.kBins, ex.binIndex);
zx = reshape(stackTransmissibility(ex.T(:, at(hit))), t4.pFeat, []);
e5 = T.estimateRecord(c4, T.seed(cfg, 4, 1));
ok = all(hit) && max(abs(z(:) - zx(:))) < 1e-10 * max(abs(zx(:))) && isequal(e5.z, z) && all(isfinite(e5.dtheta)) && ...
    isfinite(e5.q) && e5.nBad == 0;
[nPass, nFail] = check(ok, sprintf(['5  estimator: segment Welch = extractTransmissibility on %d bins (%d segments); ' ...
    'record 1: dtheta [%+.4f%% %+.4f%%] (%.2f, %.2f sigma), q %.2f, rejected bins %d'], numel(t4.kBins), Pw.K, ...
    100 * e5.dtheta, e5.dtheta ./ t4.sigma, e5.q, e5.nBad), nPass, nFail);

%% ---- 6  seeds -------------------------------------------------------------------------------------------------------
eA = T.estimateRecord(c4, T.seed(cfg, 4, 1)); eB = T.estimateRecord(c4, T.seed(cfg, 4, 2));
seeds = arrayfun(@(q) T.seed(cfg, q, 1:M.nRecordsExtended), 1:4, 'UniformOutput', false); seeds = [seeds{:}];
ok = isequal(eA.z, e5.z) && ~isequal(eA.z, eB.z) && numel(unique(seeds)) == numel(seeds) && ...
    isempty(intersect(seeds, seeds + 1e6)) && max(seeds) < min(seeds) + 1e6;
[nPass, nFail] = check(ok, sprintf(['6  seeds: same seed bitwise equal, next seed differs; %d record seeds (%d to %d) ' ...
    'unique, no wave/noise collision'], numel(seeds), min(seeds), max(seeds)), nPass, nFail);

%% ---- 7  linear algebra on a real target -----------------------------------------------------------------------------
f0 = reshape(t4.f0, t4.pFeat, []);
e0 = T.linearisedFromZ(t4, G, f0);
delta = [0.004; -0.0007];
zl = f0 + reshape(t4.J * delta, t4.pFeat, []);
eD = T.linearisedFromZ(t4, G, zl);
Fg = G.infoMatrix(t4.J(t4.rows, :), t4.blocks(:, :, t4.validMask));
ok = all(e0.dtheta == 0) && max(abs(eD.dtheta - delta)) < 1e-9 * max(abs(delta)) && ...
    abs(eD.q - delta.' * t4.F * delta) < 1e-9 * eD.q && norm(Fg - t4.F) < 1e-10 * norm(t4.F);
[nPass, nFail] = check(ok, sprintf(['7  linear algebra: z = f0 -> dtheta [%g %g]; z = f0 + J delta -> error %.1e; ' ...
    'q %.4f = delta'' F delta; G3 infoMatrix vs engine F %.1e'], e0.dtheta, max(abs(eD.dtheta - delta)), eD.q, ...
    norm(Fg - t4.F) / norm(t4.F)), nPass, nFail);

%% ---- 8  checkpoint and resume (analytic field) ----------------------------------------------------------------------
frf = analyticFRF(cfg.p0.vec);
cA = struct('index', 9, 'name', 'analytic_test', 'scn', struct('Snn', cfg.noise.lsm6dsv16x, 'frf', frf));
cxA = T.context(cA, cfg);
ckf = [tempname() '.mat'];
cfgS = cfg; cfgS.M1.nRecords = 4;
s0 = T.newState(cxA, cfgS);
full = T.mcBatch(cxA, s0, 4, struct('cfg', cfgS));
part = T.mcBatch(cxA, s0, 2, struct('cfg', cfgS, 'file', ckf, 'every', 2));
[res, okL] = T.loadCheckpoint(ckf, T.checkpointKey(cxA, cfgS));
res = T.mcBatch(cxA, res, 4, struct('cfg', cfgS, 'file', ckf, 'every', 2));
keyX = T.checkpointKey(cxA, cfgS); keyX.seedBase = keyX.seedBase + 1;
[~, okX] = T.loadCheckpoint(ckf, keyX);
cfgX = cfgS; cfgX.M1.solver.tolSigma = 0.02;                 % a changed solver setting must also invalidate it
[~, okS] = T.loadCheckpoint(ckf, T.checkpointKey(cxA, cfgX));
cfgX = cfgS; cfgX.M1.accept.coverage95 = [0.90 0.99];       % and a changed acceptance band
[~, okA] = T.loadCheckpoint(ckf, T.checkpointKey(cxA, cfgX));
okX = okX || okS || okA;
ok = okL && part.nDone == 2 && res.nDone == 4 && isequal(res.dtheta(:, 1:4), full.dtheta(:, 1:4)) && ...
    isequal(res.q(1:4), full.q(1:4)) && isequal(res.seeds(1:4), full.seeds(1:4)) && ~okX;   % (unfilled columns are NaN)
[nPass, nFail] = check(ok, sprintf(['8  checkpoint: 2 + resume 2 records bitwise equal to 4 uninterrupted (analytic field, ' ...
    'sigma [%.3f%% %.3f%%]); stale key (seed base, solver setting, acceptance band) ignored'], 100 * cxA.tgt.sigma), nPass, nFail);
delete(ckf);

%% ---- 9  borderline rule and summary -------------------------------------------------------------------------------------
B = @(r, c, n) T.isBorderline(r, c, n, cfg);
cases9 = {[1.00 1.00], 0.95, 100, false; [0.87 1.00], 0.95, 100, true; [1.00 1.00], 0.97, 100, true; ...
    [1.00 1.00], 0.96, 100, false; [0.70 1.00], 0.95, 100, false; [1.00 1.18], 0.95, 100, true; ...
    [1.00 1.00], 0.99, 200, true};
got = cellfun(@(r, c, n) B(r, c, n), cases9(:, 1), cases9(:, 2), cases9(:, 3));
sm = T.summarise(full, cxA.tgt, cfgS);
ok = isequal(got(:), [cases9{:, 4}].') && sm.n == 4 && all(isfinite(sm.empPred)) && isnan(sm.passNonlinear) && ~sm.pass;
[nPass, nFail] = check(ok, sprintf(['9  borderline rule on %d known cases (emp/pred 0.87 and 1.18, coverage 0.97 and ' ...
    '0.99: extend; 0.70: clear fail, no extension); summarise runs (n %d, no nonlinear yet -> not passed)'], ...
    size(cases9, 1), sm.n), nPass, nFail);

%% ---- 10  nonlinear solver (analytic field) ------------------------------------------------------------------------------
fun = T.makeModel(cxA);
fT = fun(cxA.tgt.thetaTrue);
nlA = T.nonlinear(cxA, full, -0.20, fun, cfgS); nlB = T.nonlinear(cxA, full, 0.20, fun, cfgS);
ok = max(abs(fT - cxA.tgt.f0)) < 1e-12 * max(abs(cxA.tgt.f0)) && nlA.pass && nlB.pass;
[nPass, nFail] = check(ok, sprintf(['10 solver (analytic): model(theta_true) = f0; beta x 0.8: %s in %d it, NL - lin ' ...
    '[%+.4f %+.4f] sigma; x 1.2: %s in %d it, [%+.4f %+.4f] sigma'], tf(nlA.converged), nlA.it, nlA.gapSigma, ...
    tf(nlB.converged), nlB.it, nlB.gapSigma), nPass, nFail);

%% ---- 11  nonlinear model on the real EMM (A1 case) ------------------------------------------------------------------------
fprintf('  ... one direct EMM evaluation of the nonlinear model (A1 case, about 30 s)\n');
tic; fA = T.makeModel(ctx{2}); f11 = fA(ctx{2}.tgt.thetaTrue); s11 = toc;
rel = max(abs(f11 - ctx{2}.tgt.f0)) / max(abs(ctx{2}.tgt.f0));
[nPass, nFail] = check(rel < 1e-10, sprintf(['11 nonlinear model on the EMM at the A1 truth = engine f0: max rel. ' ...
    'difference %.1e (frfOptions {%s}, node grid %.2g, %.0f s)'], rel, optStr(ctx{2}.scn.frfOptions), ...
    ctx{2}.scn.nodeSpacing, s11), nPass, nFail);

%% ---- 12  provenance guard -------------------------------------------------------------------------------------------------
[c1, f1] = T.cleanTree(sprintf([' M Inverse Mapping/a.m\n?? b.m\n M Results/M1/x.csv\nR  old.m -> new dir/new.m\n' ...
    '?? "space name.m"\n']));
[c2, f2] = T.cleanTree(sprintf(' M x.csv\n?? 1__MY_SUMMARIES/R1.pdf\n'));
[c3] = T.cleanTree('');
ok = ~c1 && isequal(f1, {'Inverse Mapping/a.m', 'b.m', 'new dir/new.m', 'space name.m'}) && c2 && isempty(f2) && c3;
[nPass, nFail] = check(ok, sprintf('12 provenance guard: %d .m files found, csv/pdf ignored, empty status clean', numel(f1)), ...
    nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testM1MonteCarloSpotChecks:failed', '%d check(s) failed: do not run M1.', nFail); end

%% ================================================================================================
function [nPass, nFail] = check(ok, msg, nPass, nFail)
if ok
    fprintf('PASS  %s\n', msg); nPass = nPass + 1;
else
    fprintf('FAIL  %s\n', msg); nFail = nFail + 1;
end
end

function s = tf(b)
if b, s = 'ok'; else, s = 'MISMATCH'; end
end

function s = optStr(c)
s = strjoin(cellfun(@(x) num2str(x), c, 'UniformOutput', false), ' ');
end

function t = fnText(src, name)
% the text of one function, from its "function" line up to the next "function" line (or the end)
t = '';
src = strrep(src, sprintf('\r\n'), sprintf('\n'));
k = regexp(src, ['(?m)^function [^\n=]*=\s*' name '\('], 'once');
if isempty(k), return; end
rest = src(k:end);
nx = regexp(rest(2:end), '(?m)^function ', 'once');
if ~isempty(nx), rest = rest(1:nx); end
t = regexprep(rest, '\s+$', '');
end

function frf = analyticFRF(p0)
% smooth, non-degenerate test field: acceleration FRF at points [r theta] depending on beta and R
frf = @(w, pp) fieldAt(w, pp, p0);
end

function H = fieldAt(w, pp, p0)
persistent sens
if isempty(sens)
    cfg = defineStudyScenarios(); sens = cfg.layout.level2C.frac .* [cfg.p0.R 1];
end
w = w(:).'; b = (pp(1) / p0(1))^0.35; x = sens(:, 1) / pp(3);
H = -(ones(size(sens, 1), 1) * w.^2) .* (1 + 0.6 * (x.^2 .* cos(sens(:, 2))) * ones(1, numel(w)) .* ...
    exp(1i * (0.9 * b * x) * w) + 0.25 * b * (x * w));
end