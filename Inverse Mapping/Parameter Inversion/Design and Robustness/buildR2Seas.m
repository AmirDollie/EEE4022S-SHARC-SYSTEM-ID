function [seas, mkScn] = buildR2Seas(cfg)
%BUILDR2SEAS  The 16 pre-registered R2 sea states, normalised to a fixed retained-band variance.
%   [seas, mkScn] = buildR2Seas(cfg)
%
%   Shared by runR2SpectralShape and testR2SpectralShape, so the tester checks the construction the
%   study actually uses (the only reason this is a separate file).
%
%   SEAS. 15 JONSWAP seas, omega_p in cfg.R2.omegaP x gamma_J in cfg.R2.gammaJ, plus the bimodal sea
%   cfg.R2.bimodal (two JONSWAP components). The reference sea is cfg.R2.referenceSea (= cfg.sea, the
%   G3/E1 sea: Hs 0.05 m, omega_p 5, gamma_J 3.3).
%
%   NORMALISATION (cfg.R2.normalisation = 'fixed_inband_elevation_variance'). Every sea has the same
%   incident elevation variance inside the retained estimator band cfg.R2.normalisationBand as the
%   reference sea, so R2 changes the SHAPE of the spectrum and not the in-band signal level (R1 covers
%   the level). The JONSWAP of synthesiseTwinRecords is normalised over the full spectrum to Hs, so its
%   in-band variance is proportional to Hs^2: each unimodal sea is passed as the struct form
%   (Hs, omegaP, gammaJ) with Hs_q = Hs_ref sqrt(V_ref / V_q), V = in-band variance at Hs = 1. The
%   reference sea therefore has Hs exactly cfg.sea.Hs and takes the identical code path as the E1
%   default. The bimodal sea is passed as a spectrum handle: each component is scaled so that its
%   in-band variance is (energy fraction) x (reference in-band variance).
%   The spectra are taken from synthesiseTwinRecords itself (its info.spectrum), never re-implemented.
%
%   seas(q): name, type ('jonswap' | 'bimodal'), omegaP, gammaJ (bimodal: per component), Hs (jonswap),
%     sea (what evaluateScenario receives: struct or handle), S (handle, elevation PSD per rad/s),
%     Vinband (retained-band variance, m^2), scale (Hs_q / Hs_ref; bimodal: NaN), isReference,
%     componentVinband (bimodal: per component)
%   mkScn(seaStruct, L, cacheDir): the evaluateScenario input. Only the sea and L differ from the E1
%     default (p0, Level 2C, reference s26, N, dt, band, node grid, truncation all default);
%     noise cfg.R2.noiseCase.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

band = cfg.R2.normalisationBand;
ref = cfg.R2.referenceSea;
unitS = @(wp, gJ) jonswapUnitHs(wp, gJ);
Vof = @(S) integral(S, band(1), band(2), 'RelTol', 1e-12, 'AbsTol', 0);

Sref1 = unitS(ref.omegaP, ref.gammaJ);
VrefUnit = Vof(Sref1);
Vref = VrefUnit * ref.Hs^2;                         % the target in-band variance (m^2)

seas = struct('name', {}, 'type', {}, 'omegaP', {}, 'gammaJ', {}, 'Hs', {}, 'sea', {}, 'S', {}, ...
    'Vinband', {}, 'scale', {}, 'isReference', {}, 'componentVinband', {});
for gJ = cfg.R2.gammaJ
    for wp = cfg.R2.omegaP
        S1 = unitS(wp, gJ);
        isRef = wp == ref.omegaP && gJ == ref.gammaJ;
        if isRef
            Hs = ref.Hs;                            % exactly: the E1 default sea
        else
            Hs = ref.Hs * sqrt(VrefUnit / Vof(S1));
        end
        sea = struct('Hs', Hs, 'omegaP', wp, 'gammaJ', gJ);
        S = @(w) Hs^2 * S1(w);
        seas(end + 1) = struct('name', sprintf('wp%.4g_gJ%.2g', wp, gJ), 'type', 'jonswap', 'omegaP', wp, ...
            'gammaJ', gJ, 'Hs', Hs, 'sea', sea, 'S', S, 'Vinband', Vof(S), 'scale', Hs / ref.Hs, ...
            'isReference', isRef, 'componentVinband', []); %#ok<AGROW>
    end
end
B = cfg.R2.bimodal;
a = zeros(1, 2); Sc = cell(1, 2);
for c = 1:2
    Sc{c} = unitS(B.omegaP(c), B.gammaJ(c));
    a(c) = B.energyFraction(c) * Vref / Vof(Sc{c});
end
Sb = @(w) a(1) * Sc{1}(w) + a(2) * Sc{2}(w);
seas(end + 1) = struct('name', sprintf('bimodal_wp%.4g_%.4g', B.omegaP), 'type', 'bimodal', 'omegaP', B.omegaP, ...
    'gammaJ', B.gammaJ, 'Hs', NaN, 'sea', Sb, 'S', Sb, 'Vinband', Vof(Sb), 'scale', NaN, 'isReference', false, ...
    'componentVinband', [Vof(@(w) a(1) * Sc{1}(w)), Vof(@(w) a(2) * Sc{2}(w))]);

Snn = cfg.noise.(cfg.R2.noiseCase);
mkScn = @(s, L, cacheDir) struct('Snn', Snn, 'sea', s.sea, 'L', L, 'cacheDir', cacheDir);
end

%% ================================================================================================
function S = jonswapUnitHs(wp, gJ)
% the engine's own JONSWAP at Hs = 1 (synthesiseTwinRecords info.spectrum), from a tiny analytic twin:
% no EMM solve, and the record length does not affect the spectrum handle
persistent twin
if isempty(twin)
    spec = struct('p', [1e-4 1e-3 0.4], 'sensors', [0.1 0; 0.2 0], 'nodeSpacing', 0.5, 'band', [2.9784 8.5012], ...
        'checkPoints', 0, 'frf', @(w) ones(2, numel(w)));
    twin = synthesiseTwinRecords(spec);
end
[~, ~, info] = synthesiseTwinRecords(twin, 256, 0.1, 'Spectrum', 'jonswap', 'Hs', 1, 'PeakFrequency', wp, ...
    'PeakEnhancement', gJ, 'Seed', 1);
S = info.spectrum;
end