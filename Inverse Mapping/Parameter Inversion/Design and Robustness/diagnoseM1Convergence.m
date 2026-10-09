%% diagnoseM1Convergence.m
% M1 post-run DIAGNOSTIC (not a registered criterion; the M1 verdict is unchanged): why did a nonlinear
% inversion end "NOT converged" while matching the linearised estimate?
%
% WHY. g3EstimatorCore.gaussNewton (G3, verbatim) sets converged = true only when an ACCEPTED step is below
% tolSigma (0.01 sigma) in every component. It also stops when, with a freshly computed Jacobian, none of the
% step lengths [1 0.5 0.25 0.125] lowers the objective ("at the minimum", see its comment), but that exit
% leaves converged = false. This script shows, from the stored iteration histories only (no EMM solves),
% which exit each inversion took: the size of every accepted step in sigma units, the final objective, and
% the distance between the two starts' end points.
%
% Reads the newest Results/M1/M1_*.mat. Prints only; writes nothing.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
cfg = defineStudyScenarios();
Lst = dir(fullfile(P.results, 'M1', 'M1_20*.mat'));
if isempty(Lst), error('M1diag:none', 'No M1 results found.'); end
[~, k] = max([Lst.datenum]); R = load(fullfile(Lst(k).folder, Lst(k).name)); res = R.results;
fprintf('\n==== M1 diagnostic: nonlinear termination (%s) ====\n', Lst(k).name);
fprintf('tolSigma %.3g (converged = an accepted step below this in every component)\n', cfg.M1.solver.tolSigma);

for q = 1:numel(res.names)
    st = res.state{q}; sig = res.targets{q}.sigma(:);
    if isempty(st.nl), continue; end
    fprintf('\n%s (sigma [%.4f%% %.4f%%])\n', res.names{q}, 100 * sig);
    for j = 1:numel(st.nl)
        n = st.nl(j); h = n.history;                       % rows: [ln beta; ln R; objective], one column per accepted iterate
        steps = diff(h(1:2, :), 1, 2) ./ sig;              % accepted steps in sigma units
        fprintf('  beta x %.2f: %s, %d it, %d accepted steps\n', 1 + n.betaOffset, ...
            ternary(n.converged, 'converged', 'NOT converged'), n.it, size(steps, 2));
        for s = 1:size(steps, 2)
            fprintf('    step %d: [%+.4g %+.4g] sigma, objective %.6f\n', s, steps(:, s), h(3, s + 1));
        end
        last = max(abs(steps(:, end)));
        if n.converged
            fprintf('    exit: last accepted step %.3g sigma < %.3g (converged)\n', last, cfg.M1.solver.tolSigma);
        else
            fprintf(['    exit: last accepted step %.3g sigma >= %.3g, then no descent with a fresh Jacobian ' ...
                '("at the minimum" branch, flagged not converged)\n'], last, cfg.M1.solver.tolSigma);
        end
        fprintf('    end point - linearised estimate: [%+.2e %+.2e] sigma\n', n.gapSigma);
    end
    if numel(st.nl) >= 2
        d = (st.nl(1).theta - st.nl(2).theta) ./ sig;
        fprintf('  end points of the two starts differ by [%+.2e %+.2e] sigma; objectives %.8f and %.8f\n', d, ...
            st.nl(1).chi2, st.nl(2).chi2);
    end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end