function r = rawOf(out)
%RAWOF  Compact record of the engine numbers a class is computed from, for saving with study results.
%   r = rawOf(out)   with out from evaluateScenario (or any struct carrying the same fields).
%
%   Added 2026-11-26 (paper plan F10). Studies store r at EVERY point they evaluate (every sweep-grid
%   point, layout, sea, error size), not only the class or the boundary they report, so that:
%     - the thesis classes can be recomputed exactly (regression test of classify(raw, rule));
%     - any other class rule (accuracy class + calibration flag, tolerance ladder) is post-processing;
%     - boundaries and tolerances can be re-derived under a new rule by interpolation.
%   Fields (missing ones in OUT are stored as NaN/empty, so r always has the same shape):
%     sigma [2x1] (ln beta, ln R), corr, condF, singular, F [2x2], sMax,
%     dthetaW [2x1], dW, dthetaSys [2x1], dSys, dthetaTot [2x1], dTot, class, nUsedBins
    r = struct('sigma', col(get(out, 'sigma', [NaN; NaN])), 'corr', get(out, 'corr', NaN), ...
        'condF', get(out, 'condF', NaN), 'singular', get(out, 'singular', false), ...
        'F', get(out, 'F', NaN(2)), 'sMax', get(out, 'sMax', NaN), ...
        'dthetaW', col(get(out, 'dthetaW', [NaN; NaN])), 'dW', get(out, 'dW', NaN), ...
        'dthetaSys', col(get(out, 'dthetaSys', [NaN; NaN])), 'dSys', get(out, 'dSys', NaN), ...
        'dthetaTot', col(get(out, 'dthetaTot', [NaN; NaN])), 'dTot', get(out, 'dTot', NaN), ...
        'class', get(out, 'class', ''), 'nUsedBins', get(out, 'nUsedBins', NaN));
end

function v = get(s, name, default)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function v = col(v)
    v = v(:);
end
