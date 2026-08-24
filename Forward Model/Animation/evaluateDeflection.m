function eta = evaluateDeflection(data, r, theta)
%EVALUATEDEFLECTION Cheap per-point deflection evaluation from
%precomputed data (Step 7 only). R and THETA may be same-sized arrays.
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
            IhatVal = besseli(n, km*r) ./ besseli(n, km*R);
            innerSum = innerSum + dphiFull(idx)*AfullMat(idx,nIdx)*IhatVal;
        end
        etaSum = etaSum + innerSum .* exp(1i*n*theta);
    end

    eta = etaSum / (1i*alpha);
end