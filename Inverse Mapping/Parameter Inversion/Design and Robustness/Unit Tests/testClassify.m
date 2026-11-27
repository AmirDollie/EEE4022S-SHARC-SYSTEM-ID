%% testClassify.m
% Checks of Engine/identifiabilityClass (paper plan F10, 2026-11-27). No EMM solves; under a minute.
%   1  thesis rule = the thesis classifier, copied verbatim below, on 20 000 random and edge-case points
%   2  thesis rule = the copies still inside d1LayoutTools (classCode) and d3Tools (classOf)
%   3  defaults: cfg.classRule is 'thesis'; an unknown rule is refused
%   4  accuracy rule on hand-computed cases (e_beta, ladder, calibration flag, kappa gate, singular F)
%   5  the two rules disagree exactly where they should (tiny sigma, harmless bias)
%   6  regression: the newest saved R1 run, every one of its 41 sweep points reclassified from its raw record
%   7  v1Tools.classify (used by V1 and V2) = identifiabilityClass; it refuses the accuracy rule without the bias
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
P = setupStudyPaths();
cfg = defineStudyScenarios();
nPass = 0; nFail = 0;
names = cfg.class.names;

%% 1  thesis rule = the verbatim thesis classifier
rng(20261127);
n = 20000;
pts = randomPoints(n);
mism = 0;
for i = 1:n
    a = identifiabilityClass(pts(i), cfg, 'thesis');
    b = thesisReference(pts(i), cfg);
    mism = mism + ~strcmp(a, b);
end
cnt = arrayfun(@(k) sum(arrayfun(@(p) strcmp(thesisReference(p, cfg), names{k}), pts)), 1:3);
[nPass, nFail] = check(mism == 0, sprintf(['1  thesis rule = thesis classifier on %d points (%d I / %d M / %d N, ' ...
    'incl. thresholds hit exactly, singular, Inf and NaN): %d mismatches'], n, cnt, mism), nPass, nFail);

%% 2  thesis rule = the copies still in d1LayoutTools and d3Tools
T1 = d1LayoutTools(); T3 = d3Tools();
m1 = 0; m3 = 0;
for i = 1:n
    c = find(strcmp(identifiabilityClass(pts(i), cfg, 'thesis'), names));
    q = struct('sigma', pts(i).sigma, 'singular', pts(i).singular, 'dW', pts(i).dTot, 'kappa', pts(i).condF);
    m1 = m1 + (T1.classCode(q, cfg) ~= c);
    m3 = m3 + (T3.classOf(pts(i).sigma, pts(i).dTot, pts(i).condF, pts(i).singular, cfg) ~= c);
end
[nPass, nFail] = check(m1 == 0 && m3 == 0, sprintf('2  = d1LayoutTools.classCode (%d mismatches) and d3Tools.classOf (%d)', ...
    m1, m3), nPass, nFail);

%% 3  defaults
raw0 = struct('sigma', [0.0005; 0.0001], 'singular', false, 'condF', 80, 'dTot', 0.4, 'dthetaTot', [0.0002; 0]);
[~, ~, i0] = identifiabilityClass(raw0, cfg);
cfgNo = rmfield(cfg, 'classRule'); [~, ~, i1] = identifiabilityClass(raw0, cfgNo);
try
    identifiabilityClass(raw0, cfg, 'nonsense'); refused = false;
catch err
    refused = strcmp(err.identifier, 'identifiabilityClass:rule');
end
[nPass, nFail] = check(strcmp(cfg.classRule, 'thesis') && strcmp(i0.rule, 'thesis') && strcmp(i1.rule, 'thesis') && refused, ...
    '3  cfg.classRule = ''thesis''; no classRule -> thesis; an unknown rule is refused', nPass, nFail);

%% 4  accuracy rule, hand cases (sigma and bias in ln units)
A = cfg.classAccuracy;
cases = { ...  % sigma_lnb   sigma_lnR  bias_lnb  bias_lnR  condF   singular  dTot   -> code  ladder  calibrated
    'p0 (R1, LSM6DSV16X)',      0.000500, 0.000084, 0.000200, 0.00002, 83,   false, 0.436,  1, 0.01, true;
    'e just under 5%',          0.030,    0.001,    0.0399,   0,       100,  false, 1.4,    1, 0.05, false;
    'e = 12% (Marginal)',       0.012,    0.001,    0.1194,   0,       100,  false, 10,     2, 0.20, false;
    'e just over 20%',          0.12,     0.01,     0.1601,   0,       100,  false, 1.4,    3, NaN,  false;
    'V1c -20 mm (+61.8%)',      0.0005,   0.0001,   log(1.618), 0.096, 83,   false, 1164,   3, NaN,  false;
    'small e, kappa > 1e8',     0.001,    0.001,    0,        0,       2e8,  false, 0.1,    3, 0.01, true;
    'singular F',               Inf,      Inf,      NaN,      NaN,     Inf,  true,  NaN,    3, NaN,  false};
ok4 = true; msg = {};
for i = 1:size(cases, 1)
    c = cases(i, :);
    raw = struct('sigma', [c{2}; c{3}], 'singular', c{7}, 'condF', c{6}, 'dTot', c{8}, 'dthetaTot', [c{4}; c{5}]);
    [nm, code, inf4] = identifiabilityClass(raw, cfg, 'accuracy');
    lad = inf4.ladder;
    okI = code == c{9} && strcmp(nm, names{c{9}}) && (isequaln(lad, c{10}) || abs(lad - c{10}) < 1e-12) && ...
        inf4.calibrated == c{11};
    if ~c{7}, okI = okI && abs(inf4.eBeta - hypot(c{2}, c{4})) < 1e-15; end
    ok4 = ok4 && okI;
    if ~okI, msg{end + 1} = sprintf('%s: got %s, ladder %g, calibrated %d', c{1}, nm, lad, inf4.calibrated); end %#ok<SAGROW>
end
[nPass, nFail] = check(ok4, sprintf('4  accuracy rule on %d hand cases (thresholds %g/%g, ladder %s)%s', size(cases, 1), ...
    A.identifiableMaxE, A.marginalMaxE, mat2str(A.ladder), sprintf(' | %s', msg{:})), nPass, nFail);

%% 5  where the rules disagree: V1c-like 0.44 mm radial error at p0 (|dbeta/beta| ~ 1%, d_sys ~ 25)
raw5 = struct('sigma', [0.0005; 0.000084], 'singular', false, 'condF', 83, 'dTot', 25, 'dthetaTot', [0.00995; 0.002]);
[nT, ~, ~] = identifiabilityClass(raw5, cfg, 'thesis');
[nA, ~, iA] = identifiabilityClass(raw5, cfg, 'accuracy');
[nPass, nFail] = check(strcmp(nT, names{3}) && strcmp(nA, names{1}) && ~iA.calibrated && abs(iA.ladder - 0.01) < 1e-12, ...
    sprintf(['5  1%% bias in beta at sigma 0.05%%: thesis %s (d_tot 25), accuracy %s with e_beta %.2f%% ' ...
    '(ladder %g%%), flagged not calibrated'], nT, nA, 100 * iA.eBeta, 100 * iA.ladder), nPass, nFail);

%% 6  regression on the newest saved R1 run
L = dir(fullfile(P.results, 'R1', 'R1_20*.mat'));
if isempty(L)
    fprintf('  SKIP  6  no saved R1 run in %s\n', fullfile(P.results, 'R1'));
else
    [~, k] = max([L.datenum]); S = load(fullfile(L(k).folder, L(k).name));
    G = S.results.grid;
    if ~isfield(G, 'raw') || isempty(G.raw)
        fprintf('  SKIP  6  %s has no raw records (saved before 2026-11-26)\n', L(k).name);
    else
        bad6 = 0;
        for j = 1:numel(G.raw)
            bad6 = bad6 + ~strcmp(identifiabilityClass(G.raw{j}, cfg, 'thesis'), G.raw{j}.class);
        end
        acc = cellfun(@(r) identifiabilityClass(r, cfg, 'accuracy'), G.raw, 'UniformOutput', false);
        [nPass, nFail] = check(bad6 == 0, sprintf(['6  %s: %d of %d saved points reclassified to their saved class ' ...
            '(accuracy rule, for information: %d I / %d M / %d N)'], L(k).name, numel(G.raw) - bad6, numel(G.raw), ...
            sum(strcmp(acc, names{1})), sum(strcmp(acc, names{2})), sum(strcmp(acc, names{3}))), nPass, nFail);
    end
end

%% 7  v1Tools.classify (V1, V2) goes through identifiabilityClass
V = v1Tools();
out = struct('sigma', [0.0005; 0.000084], 'singular', false, 'condF', 83, 'sMax', 0.0005);
same = true;
for d = [0.1 0.5 0.7 1 3]
    same = same && strcmp(V.classify(out, d, cfg), identifiabilityClass(setfield(out, 'dTot', d), cfg, 'thesis')); %#ok<SFLD>
end
cfgA = cfg; cfgA.classRule = 'accuracy';
try
    V.classify(out, 0.4, cfgA); refused7 = false;
catch err
    refused7 = strcmp(err.identifier, 'v1Tools:classifyBias');
end
okA = strcmp(V.classify(out, 25, cfgA, [0.00995; 0.002]), names{1});
[nPass, nFail] = check(same && refused7 && okA, ['7  v1Tools.classify = identifiabilityClass (thesis); under ' ...
    '''accuracy'' it needs dtheta_tot and uses it'], nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0, error('testClassify:failed', '%d check(s) failed.', nFail); end

%% ---- helpers ---------------------------------------------------------------------------------------------
function c = thesisReference(m, cfg)
% VERBATIM copy of the classifier in evaluateScenario before 2026-11-27 (and v1Tools, runR3, runR4)
I = cfg.class.identifiable; M = cfg.class.marginal;
if m.singular || ~all(isfinite(m.sigma))
    c = cfg.class.names{3};
    return
end
d = m.dTot;
if m.sMax <= I.maxSigma && d <= I.maxBiasDistance && m.condF <= I.maxCondF
    c = cfg.class.names{1};
elseif m.sMax <= M.maxSigma && d <= M.maxBiasDistance && m.condF <= M.maxCondF
    c = cfg.class.names{2};
else
    c = cfg.class.names{3};
end
end

function pts = randomPoints(n)
% log-uniform sigma, d and kappa across all three classes, with 10% of each coordinate set EXACTLY on a threshold
% and a few singular / Inf / NaN points
edgeS = [0.01 0.05]; edgeD = [0.5 1]; edgeK = 1e8;
pts = repmat(struct('sigma', [0; 0], 'sMax', 0, 'singular', false, 'condF', 0, 'dTot', 0, 'dthetaTot', [0; 0]), n, 1);
for i = 1:n
    s = 10 .^ (-4 + 3.5 * rand(2, 1));
    if rand < 0.1, s(randi(2)) = edgeS(randi(2)); end
    d = 10 ^ (-2 + 3 * rand); if rand < 0.1, d = edgeD(randi(2)); end
    k = 10 ^ (1 + 8 * rand); if rand < 0.05, k = edgeK; end
    sg = false;
    u = rand;
    if u < 0.01, sg = true; s = [Inf; Inf]; d = NaN; k = Inf;
    elseif u < 0.02, s(randi(2)) = Inf;
    elseif u < 0.03, d = NaN;
    end
    pts(i).sigma = s; pts(i).sMax = max(s); pts(i).singular = sg; pts(i).condF = k; pts(i).dTot = d;
    pts(i).dthetaTot = d * s .* sign(randn(2, 1)) / sqrt(2);
end
end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end
