function k = openWaterRoots(alpha, Mo)
%OPENWATERROOTS  Open-water dispersion roots [k0; k1; ...; kMo] for large Mo.   (F1, 2026-11-25)
%   Same output as dispersionRoots(alpha, 0, 0, Mo), built for Mo in the tens of thousands.
%   k0 = i p0 comes from dispersionRoots. For m >= 1 the root of f(k) = k sin k + alpha cos k lies in
%   ((m - 1/2) pi, m pi), where f changes sign. fzero is used for m <= 50; beyond that a vectorised
%   Newton iteration from k ~ m pi - atan(alpha / (m pi)) converges in a few steps, and any root that
%   lands outside its bracket is redone with fzero.
%   Verified against dispersionRoots to 3e-16 relative at Mo = 500 (diag_F1_Mo_ladder.m, Part 0).
    k0 = dispersionRoots(alpha, 0, 0, 1);
    k0 = k0(1);
    f  = @(x) x .* sin(x) + alpha * cos(x);
    fp = @(x) (1 - alpha) * sin(x) + x .* cos(x);
    km = zeros(Mo, 1);
    mSmall = min(50, Mo);
    for i = 1:mSmall
        km(i) = fzero(f, [(i - 0.5) * pi, i * pi]);
    end
    if Mo > mSmall
        mm = (mSmall + 1:Mo).';
        x = mm * pi - atan(alpha ./ (mm * pi));
        for it = 1:8
            x = x - f(x) ./ fp(x);
        end
        bad = ~(x > (mm - 0.5) * pi & x < mm * pi) | ~isfinite(x);
        for i = find(bad).'
            x(i) = fzero(f, [(mm(i) - 0.5) * pi, mm(i) * pi]);
        end
        km(mSmall + 1:end) = x;
    end
    k = [k0; km];
end
