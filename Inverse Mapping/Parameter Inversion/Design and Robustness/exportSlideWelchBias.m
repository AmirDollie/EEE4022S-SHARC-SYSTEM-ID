%EXPORTSLIDEWELCHBIAS  Data for talk slide 4 (slide_welch_bias.pdf): pointwise EMM transmissibility
%   against the Welch-expected transmissibility (ratio of Hann-kernel-smoothed spectra) at the
%   production L = 2048 and a short L, plus the normalised spectral kernels.
%
%   Same case as slide 3: p0, Level 2C [s2 s6 s14 s26], reference s26, j = s14, G3 JONSWAP sea,
%   N = 65536, dt = 0.1 s. Uses the frozen engine only (evaluateScenario, synthesiseTwinRecords,
%   welchExpectedTransmissibility), so the smoothing is applied to numerator and denominator
%   spectra separately on the twin's own fine grid (full synthesis band with its taper), not to T.
%
%   Place in Inverse Mapping/Parameter Inversion/Design and Robustness/ and run. No new EMM solves
%   if the p0 field is cached. Writes, next to this script:
%     slide_welch_bias.csv    omega, ReT, ImT, ReW_<Llong>, ImW_<Llong>, ReW_<Lshort>, ImW_<Lshort>
%     slide_welch_kernel.csv  dw, K_<Llong>, K_<Lshort>   (each normalised to peak 1)
%     slide_welch_bins.csv    omega bins actually retained at each L (L, omega)
%   and prints d_W at both L for the talk notes.

P = setupStudyPaths();
cfg = defineStudyScenarios();
Llong = cfg.welch.L;            % 2048, production
Lshort = 256;                   % try 512 if 256 looks too dramatic
jName = 14;                     % sensor j shown (s14)

%% ---- field and twin info from the engine (p0, Level 2C, ref s26) ----
outL = evaluateScenario(struct('L', Llong, 'cacheDir', P.cache), cfg);
outS = evaluateScenario(struct('L', Lshort, 'cacheDir', P.cache), cfg);
scn = outL.scn; fld = outL.field;
lay = cfg.layout.level2C.idx;                       % [2 6 14 26]
assert(isequal(scn.layout, 1:4) && lay(scn.ref) == 26, 'expected Level 2C with reference s26');
others = setdiff(1:numel(scn.layout), scn.ref);
jIdx = find(lay(others) == jName);                  % position of s14 among the non-reference rows

sea = scn.sea;
[~, ~, info] = synthesiseTwinRecords(fld.twin, scn.N, scn.dt, 'Seed', 1, 'Spectrum', 'jonswap', ...
    'Hs', sea.Hs, 'PeakFrequency', sea.omegaP, 'PeakEnhancement', sea.gammaJ);
infoL = info; infoL.H = info.H(scn.layout, :);

%% ---- curves on a dense grid inside the retained band ----
w = linspace(cfg.band.basin(1), cfg.band.basin(2), 400);
H = fld.H(w); H = H(scn.layout, :);
T = H(others(jIdx), :) ./ H(scn.ref, :);                                    % pointwise EMM
TWl = welchExpectedTransmissibility(w, infoL, Llong, scn.ref, others, 0);  % corrected, no Snn left
TWs = welchExpectedTransmissibility(w, infoL, Lshort, scn.ref, others, 0);
TWl = TWl(jIdx, :); TWs = TWs(jIdx, :);

here = fileparts(mfilename('fullpath'));
M = [w(:), real(T(:)), imag(T(:)), real(TWl(:)), imag(TWl(:)), real(TWs(:)), imag(TWs(:))];
fid = fopen(fullfile(here, 'slide_welch_bias.csv'), 'w');
fprintf(fid, 'omega,ReT,ImT,ReW_%d,ImW_%d,ReW_%d,ImW_%d\n', Llong, Llong, Lshort, Lshort);
fprintf(fid, '%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', M.'); fclose(fid);

%% ---- spectral kernels |W_L(dw dt)|^2, normalised to peak 1 ----
dw = linspace(-1.2, 1.2, 1201);
Kl = abs(hannDTFT(dw * scn.dt, Llong)).^2; Kl = Kl / max(Kl);
Ks = abs(hannDTFT(dw * scn.dt, Lshort)).^2; Ks = Ks / max(Ks);
fid = fopen(fullfile(here, 'slide_welch_kernel.csv'), 'w');
fprintf(fid, 'dw,K_%d,K_%d\n', Llong, Lshort);
fprintf(fid, '%.10g,%.10g,%.10g\n', [dw; Kl; Ks]); fclose(fid);

%% ---- retained bins (for markers) ----
fid = fopen(fullfile(here, 'slide_welch_bins.csv'), 'w'); fprintf(fid, 'L,omega\n');
fprintf(fid, '%d,%.10g\n', [Llong * ones(1, numel(outL.omega)); outL.omega(:).']);
fprintf(fid, '%d,%.10g\n', [Lshort * ones(1, numel(outS.omega)); outS.omega(:).']); fclose(fid);

relL = abs(TWl - T) ./ abs(T); relS = abs(TWs - T) ./ abs(T);
fprintf('L = %4d: d_W = %.3f, max |Tbar - T|/|T| on s14 = %.2e, max |dphi| = %.3f deg\n', Llong, outL.dW, ...
    max(relL), max(abs(angle(TWl ./ T))) * 180 / pi);
fprintf('L = %4d: d_W = %.3f, max |Tbar - T|/|T| on s14 = %.2e, max |dphi| = %.3f deg\n', Lshort, outS.dW, ...
    max(relS), max(abs(angle(TWs ./ T))) * 180 / pi);
fprintf('sigma(ln beta) at L = %d: %.4g, at L = %d: %.4g\n', Llong, outL.sigma(1), Lshort, outS.sigma(1));

function W = hannDTFT(nu, L)      % as in welchExpectedTransmissibility.m (periodic Hann)
W = 0.5 * dirichletSum(nu, L) - 0.25 * dirichletSum(nu - 2 * pi / L, L) - 0.25 * dirichletSum(nu + 2 * pi / L, L);
end
function D = dirichletSum(x, L)
den = 1 - exp(-1i * x);
D = (1 - exp(-1i * x * L)) ./ den;
small = abs(den) < 1e-9; D(small) = L;
end