function G = g3EstimatorCore()
%G3ESTIMATORCORE  The G3 estimator and solver operations, VERBATIM, so M1 can call them.
%   G = g3EstimatorCore();   then e.g.  P = G.welchSegments(Y, dt, L, ov, kBins, ref, others)
%
%   WHY THIS FILE EXISTS. Gate G3 (runTransmissibilityTwinInversion.m, Output-Only Twin/) validated the
%   finite-record estimator, the linearised estimate and the safeguarded Gauss-Newton solver, but keeps
%   them as LOCAL functions of a script, which nothing else can call. M1 must reuse that code unchanged.
%   The seven functions below are copied character for character from that script (extracted by a
%   program, not retyped); testM1MonteCarloSpotChecks compares every function body against the G3
%   source and fails on any difference. G3 itself is not modified.
%
%   FUNCTIONS (signatures and conventions exactly as in G3)
%   P = G.welchSegments(Y, dt, L, ov, kBins, ref, others)   per-segment scaled cross-spectra
%   z = G.zFromSegments(P, segIdx, Snn)                     corrected features (p x nF); bins failing
%                                                           den > 0.1 S_rr are NaN
%   F = G.infoMatrix(J, blocks)                             sum_k J_k' Sigma_k^-1 J_k
%   g = G.gradTerm(J, blocks, r)                            J' Sigma^-1 r (block diagonal Sigma)
%   v = G.objective(blocks, r)                              r' Sigma^-1 r
%   J = G.jacobianCD(fun, theta, h)                         central differences
%   g = G.gaussNewton(fun, z, blocks, theta, J, sigma, h, maxIt, stepSigma, tolSigma)
%                                                           safeguarded modified Gauss-Newton
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
G.welchSegments = @welchSegments;
G.zFromSegments = @zFromSegments;
G.infoMatrix = @infoMatrix;
G.gradTerm = @gradTerm;
G.objective = @objective;
G.jacobianCD = @jacobianCD;
G.gaussNewton = @gaussNewton;
G.source = 'Output-Only Twin/runTransmissibilityTwinInversion.m';
G.names = {'welchSegments', 'zFromSegments', 'infoMatrix', 'gradTerm', 'objective', 'jacobianCD', 'gaussNewton'};
end

%% ================================================================================================
%  Verbatim from runTransmissibilityTwinInversion.m (do not edit: the M1 tester compares them)
%% ================================================================================================
function P = welchSegments(Y, dt, L, ov, kBins, ref, others)
% Per-segment scaled cross-spectra, identical conventions to extractTransmissibility.
n = size(Y, 2);
w = 0.5 * (1 - cos(2 * pi * (0:L - 1).' / L));
U = sum(w.^2);
hop = max(1, round(L * (1 - ov)));
starts = 1:hop:(n - L + 1);
K = numel(starts);
scale = 2 * dt / (U * 2 * pi);
nF = numel(kBins);
P.Sjr = complex(zeros(numel(others), nF, K)); P.Srr = zeros(1, nF, K); P.K = K;
for s = 1:K
    seg = Y(:, starts(s):starts(s) + L - 1);
    seg = seg - mean(seg, 2);
    X = fft(seg .* w.', [], 2);
    X = X(:, kBins + 1);
    P.Sjr(:, :, s) = scale * X(others, :) .* conj(X(ref, :));
    P.Srr(1, :, s) = scale * abs(X(ref, :)).^2;
end
end

function z = zFromSegments(P, segIdx, Snn)
% Corrected transmissibility features (p x nF) from a subset of segments.
Sjr = mean(P.Sjr(:, :, segIdx), 3); Srr = mean(P.Srr(:, :, segIdx), 3);
den = Srr - Snn;
T = Sjr ./ den;
T(:, ~(den > 0.1 * Srr)) = complex(NaN, NaN);
z = reshape(stackTransmissibility(T), 2 * size(T, 1), []);
end

function F = infoMatrix(J, blocks)
[p, ~, nF] = size(blocks);
F = zeros(size(J, 2));
for k = 1:nF
    Jk = J((k - 1) * p + (1:p), :);
    F = F + Jk.' * (blocks(:, :, k) \ Jk);
end
end

function g = gradTerm(J, blocks, r)
% J' Sigma^-1 r for blkdiag Sigma; r is p x nF or a stacked vector.
[p, ~, nF] = size(blocks);
r = reshape(r, p, nF);
g = zeros(size(J, 2), 1);
for k = 1:nF
    Jk = J((k - 1) * p + (1:p), :);
    g = g + Jk.' * (blocks(:, :, k) \ r(:, k));
end
end

function v = objective(blocks, r)
[p, ~, nF] = size(blocks);
r = reshape(r, p, nF);
v = 0;
for k = 1:nF
    v = v + r(:, k).' * (blocks(:, :, k) \ r(:, k));
end
end

function J = jacobianCD(fun, theta, h)
f1 = fun(theta);
J = zeros(numel(f1), numel(theta));
for k = 1:numel(theta)
    e = zeros(size(theta)); e(k) = h;
    J(:, k) = (fun(theta + e) - fun(theta - e)) / (2 * h);
end
end

function g = gaussNewton(fun, z, blocks, theta, J, sigma, h, maxIt, stepSigma, tolSigma)
% Safeguarded modified Gauss-Newton (J reused, recomputed when needed).
evals = 0;
f = fun(theta); evals = evals + 1;
r = z(:) - f; obj = objective(blocks, r);
g.jRecomputes = 0; g.converged = false; g.history = [theta; obj];
sinceJ = 0; jFresh = true;                                  % J(theta0) was computed at the start point
for it = 1:maxIt
    F = infoMatrix(J, blocks);
    step = F \ gradTerm(J, blocks, r);
    accepted = false;
    for a = [1 0.5 0.25 0.125]
        thNew = theta + a * step;
        fNew = fun(thNew); evals = evals + 1;
        rNew = z(:) - fNew; objNew = objective(blocks, rNew);
        if objNew < obj, accepted = true; break; end
    end
    if ~accepted
        if jFresh
            break                                            % J is current and still no descent: at the minimum
        end
        J = jacobianCD(fun, theta, h); evals = evals + 4; g.jRecomputes = g.jRecomputes + 1; sinceJ = 0; jFresh = true;
        continue
    end
    theta = thNew; r = rNew; obj = objNew; sinceJ = sinceJ + 1; jFresh = false;
    g.history(:, end + 1) = [theta; obj];
    if all(abs(a * step) ./ sigma < tolSigma)
        g.converged = true; break
    end
    if any(abs(a * step) ./ sigma > stepSigma) || sinceJ >= 3
        J = jacobianCD(fun, theta, h); evals = evals + 4; g.jRecomputes = g.jRecomputes + 1; sinceJ = 0; jFresh = true;
    end
end
g.theta = theta; g.obj = obj; g.it = it; g.evals = evals;
end