function [name, code, info] = identifiabilityClass(raw, cfg, rule)
%IDENTIFIABILITYCLASS  Class of one evaluated point from its raw engine numbers, under a named rule.
%   [name, code, info] = identifiabilityClass(raw, cfg)          rule = cfg.classRule ('thesis' by default)
%   [name, code, info] = identifiabilityClass(raw, cfg, rule)    rule = 'thesis' | 'accuracy'
%
%   RAW is an evaluateScenario output, a rawOf record, or any struct with the same fields:
%     sigma [2x1] (ln beta, ln R), singular, condF       always required
%     dTot                                               thesis rule (falls back to dW if absent)
%     dthetaTot [2x1]                                    accuracy rule (falls back to dthetaW (+ dthetaSys))
%     sMax                                               optional (recomputed as max(sigma) if absent)
%
%   RULE 'thesis' (the pre-registered rule, cfg.class; bit-for-bit the classifier that was duplicated in
%   evaluateScenario, v1Tools, d1LayoutTools, d3Tools, runR3, runR4 before 2026-11-27):
%     singular F or a non-finite sigma                            -> Not identifiable
%     s_max <= 1%, d_tot <= 0.5, kappa <= 1e8                     -> Identifiable
%     s_max <= 5%, d_tot <= 1,   kappa <= 1e8                     -> Marginal
%     otherwise                                                   -> Not identifiable
%
%   RULE 'accuracy' (paper plan F10; cfg.classAccuracy). The thesis rule measures bias in Fisher units
%   (d = bias / sigma), so when sigma is tiny (0.05% at p0) a physically negligible bias (9 um of sensor
%   position, 0.008 mm of draught) fails it. The accuracy rule asks the question a user of the method asks
%   instead, "how far from the truth is beta likely to be?", and keeps the Fisher distance as a separate flag:
%     e_beta = sqrt(sigma_lnbeta^2 + bias_lnbeta^2)   (ln units ~ relative error; bias = dthetaTot(1))
%     singular F, non-finite sigma, or kappa > maxCondF           -> Not identifiable
%     e_beta <= 5%                                                -> Identifiable
%     e_beta <= 20%                                               -> Marginal
%     otherwise                                                   -> Not identifiable
%     info.calibrated = d_tot <= 0.5: the reported sigma is an honest error bar (bias small against it).
%     info.ladder = the smallest of [1 5 10 20]% that e_beta meets (NaN if none): finer reporting levels.
%   The class names (and codes 1, 2, 3) are the thesis ones, so every table, plot and code path that uses
%   cfg.class.names works unchanged under either rule; info.rule says which rule produced them.
%
%   INFO: rule, eBeta, eR (same construction for ln R; reported, not classified), biasBeta, biasR,
%         calibrated, ladder, dTot, sMax
%
% Added 2026-11-27 (paper plan F10, steps 1-3). Lives in Design and Robustness/Engine/.

if nargin < 3 || isempty(rule)
    if isfield(cfg, 'classRule') && ~isempty(cfg.classRule), rule = cfg.classRule; else, rule = 'thesis'; end
end
rule = lower(char(rule));

sigma = raw.sigma(:);
singular = isfield(raw, 'singular') && ~isempty(raw.singular) && logical(raw.singular);
condF = fieldOr(raw, 'condF', NaN);
sMax = fieldOr(raw, 'sMax', NaN);
if ~isfinite(sMax), sMax = max(sigma(1:min(2, end))); end
dTot = fieldOr(raw, 'dTot', NaN);
if ~isfield(raw, 'dTot') || isempty(raw.dTot)
    dTot = fieldOr(raw, 'dW', NaN);                  % Welch-only records (D1, D3a): d_tot = d_W
end
dth = totalBias(raw);

bad = singular || ~all(isfinite(sigma));
info = struct('rule', rule, 'eBeta', NaN, 'eR', NaN, 'biasBeta', dth(1), 'biasR', dth(2), ...
    'calibrated', ~bad && dTot <= cfg.class.identifiable.maxBiasDistance, 'ladder', NaN, 'dTot', dTot, 'sMax', sMax);

switch rule
    case 'thesis'
        I = cfg.class.identifiable; M = cfg.class.marginal;
        if bad
            code = 3;
        elseif sMax <= I.maxSigma && dTot <= I.maxBiasDistance && condF <= I.maxCondF
            code = 1;
        elseif sMax <= M.maxSigma && dTot <= M.maxBiasDistance && condF <= M.maxCondF
            code = 2;
        else
            code = 3;
        end
    case 'accuracy'
        A = cfg.classAccuracy;
        info.calibrated = ~bad && dTot <= A.calibrationMaxD;
        if ~bad
            info.eBeta = sqrt(sigma(1)^2 + dth(1)^2);
            info.eR = sqrt(sigma(2)^2 + dth(2)^2);
            k = find(info.eBeta <= A.ladder, 1);
            if ~isempty(k), info.ladder = A.ladder(k); end
        end
        if bad || ~(condF <= A.maxCondF) || ~isfinite(info.eBeta)
            code = 3;
        elseif info.eBeta <= A.identifiableMaxE
            code = 1;
        elseif info.eBeta <= A.marginalMaxE
            code = 2;
        else
            code = 3;
        end
    otherwise
        error('identifiabilityClass:rule', 'Unknown class rule ''%s'' (use ''thesis'' or ''accuracy'').', rule);
end
name = cfg.class.names{code};
end

function dth = totalBias(raw)
% dthetaTot if present; else dthetaW plus dthetaSys (whichever exist); else zeros (no bias known)
dth = [0; 0];
if isfield(raw, 'dthetaTot') && numel(raw.dthetaTot) >= 2 && all(isfinite(raw.dthetaTot(1:2)))
    dth = raw.dthetaTot(1:2); dth = dth(:); return
end
if isfield(raw, 'dthetaW') && numel(raw.dthetaW) >= 2 && all(isfinite(raw.dthetaW(1:2)))
    dth = raw.dthetaW(1:2); dth = dth(:);
end
if isfield(raw, 'dthetaSys') && numel(raw.dthetaSys) >= 2 && all(isfinite(raw.dthetaSys(1:2)))
    dth = dth + reshape(raw.dthetaSys(1:2), [], 1);
end
end

function v = fieldOr(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
