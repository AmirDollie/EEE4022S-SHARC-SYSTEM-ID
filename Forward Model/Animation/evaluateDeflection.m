function eta = evaluateDeflection(data, r, theta)
%EVALUATEDEFLECTION Cheap per-point deflection evaluation from
%precomputed data (Step 7 only). R and THETA may be same-sized arrays.
% F1 (2026-11-25): Ihat_n(kappa r) now from besselIRatio (overflow-safe). Nothing else changed.
% (Production file lives in Forward Model/Animation/.)
    addpath('..');

    kappaFull = data.kappa;
    dphiFull = data.dphiFull;
    AfullMat = data.AfullMat;
    N = data.N;
    R = data.R;
    alpha = data.alpha;

    etaSum = zeros(size(r));
    for nIdx = 1:(2*N+1)
        n = nIdx - N - 1;
        innerSum = zeros(size(r));
        for idx = 1:length(kappaFull)
            km = kappaFull(idx);
            IhatVal = besselIRatio(n, km*r, km*R);
            innerSum = innerSum + dphiFull(idx)*AfullMat(idx,nIdx)*IhatVal;
        end
        etaSum = etaSum + innerSum .* exp(1i*n*theta);
    end

    eta = etaSum / (1i*alpha);
end
