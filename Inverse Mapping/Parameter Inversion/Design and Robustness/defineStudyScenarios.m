function cfg = defineStudyScenarios()
%DEFINESTUDYSCENARIOS  Pre-registered configuration of the final experiment programme.
%   cfg = defineStudyScenarios()
%
%   WHAT WILL BE TESTED, fixed before any sweep is run (SHARC_Experiment_Programme.tex, rev. 3).
%   This file is data: scenario lists, sweep grids, admissibility rules and class thresholds.
%   It runs no EMM solve, Fisher calculation, optimisation or plot. The only arithmetic is the
%   deterministic conversion of finite-depth physical floes to the EMM's scaling (same formulas
%   as screenFloeRegime.m), so every study reads identical (beta, gamma, R) values.
%   Deep-water (MIZ) cases are stored physically; their EMM scaling needs the dispersion
%   relation and is produced by screenFloeRegime.m in the study script.
%
%   Division of labour:
%     defineStudyScenarios.m  what is tested (this file)
%     evaluateScenario.m      how one scenario is evaluated
%     run<Study>.m            loops, optimisation, outputs
%
%   CHANGE CONTROL: after the pre-registration date below, any change is recorded in
%   cfg.meta.changes with a date and reason; never edit a value silently.
%
%   EMM SCALING (from the Forward Model, confirmed in computeSensorFRF.m and screenFloeRegime.m)
%     length scale = water depth H;  alpha = H omega^2 / g
%     beta  = D / (rhoW g H^4) = (l_f / H)^4,   D = E h^3 / (12 (1 - nu^2))
%     gamma = draught / H = rhoI h / (rhoW H)
%     R     = Rphys / H;   sensor radius r = rphys / H
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

%% ---- Metadata ---------------------------------------------------------------------------------
cfg.meta.created = '2026-09-30';
cfg.meta.version = '1.0';
cfg.meta.status = 'PRE-REGISTERED';
cfg.meta.programme = 'SHARC_Experiment_Programme.tex rev. 3';
cfg.meta.changes = struct('date', {}, 'field', {}, 'reason', {});
cfg.meta.changes(end + 1) = struct('date', '2026-10-01', 'field', 'E1.tolD (added)', 'reason', ...
    ['AMENDMENT, not pre-registered. The pre-registration fixed the d_W targets (0.26, 0.91, 3.68) but no ' ...
    'acceptance tolerance. After the first analytic comparison a 10% tolerance was introduced; it reflects the ' ...
    'sampling uncertainty of the 200-realisation Monte Carlo covariance and its inverse. The observed analytic ' ...
    'values are 6-7% below the Monte Carlo-derived targets.']);
cfg.meta.changes(end + 1) = struct('date', '2026-10-01', 'field', 'R4.noiseCases, headlineNoise, toleranceSearch, toleranceN (added)', ...
    'reason', ['CLARIFICATION, no value changed. These choices were made in runR4GammaTolerance before it was ' ...
    'run (all noise cases, headline = cfg.noise.default, tolerance search eps in [-0.30, 0.30] on 1201 points). ' ...
    'They are moved here so the configuration alone records what R4 tested.']);

%% ---- Physics and the basin reference p0 -------------------------------------------------------
cfg.physics.g = 9.81;
cfg.physics.Hbasin = 1.88;          % m, basin depth (EMM default WaterDepth)
cfg.physics.rhoWbasin = 1000;       % kg/m^3, fresh-water basin (runFloeRegimeScreen)
cfg.physics.rhoWsea = 1025;         % kg/m^3, MIZ
cfg.physics.nuEMM = 0.30;           % Poisson ratio used for p0

cfg.p0 = struct('beta', 4.6985e-5, 'gamma', 1.4548e-3, 'R', 0.3830);
cfg.p0.vec = [cfg.p0.beta, cfg.p0.gamma, cfg.p0.R];
cfg.p0.Rphys = cfg.p0.R * cfg.physics.Hbasin;                                            % 0.720 m
cfg.p0.D = cfg.p0.beta * cfg.physics.rhoWbasin * cfg.physics.g * cfg.physics.Hbasin^4;  % 5.76 N m
cfg.p0.m = cfg.p0.gamma * cfg.physics.Hbasin * cfg.physics.rhoWbasin;                   % 2.74 kg/m^2

%% ---- Frozen estimator and forward-model settings ----------------------------------------------
cfg.acq.dt = 0.1;
cfg.acq.Nref = 65536;               % 109 min
cfg.welch = struct('L', 2048, 'overlap', 0.5, 'window', 'hann-periodic', 'binStep', 2, 'NeffRef', 59.7);
cfg.band.basin = [3.178 8.301];     % rad/s, retained band (every 2nd L = 2048 bin, 84 bins)
cfg.band.estimate = [3 8.5];        % rad/s, estimator band before bin selection
cfg.emm.truncation = [50 10 10];    % [M P N]
cfg.emm.truncationCheck = [70 15 15];
cfg.emm.alphaValidated = [1.7 13.85];
cfg.emm.nodeSpacing = 0.2;          % rad/s, solver model grid (validated at p0, A, B only)
cfg.emm.nodeSpacingCheck = 0.1;     % rad/s, record grid; A1 checks 0.2 against 0.1
cfg.emm.maxInterpBiasDistance = 0.1;
cfg.emm.jacobianStep = 1e-3;
cfg.sea = struct('Hs', 0.05, 'omegaP', 5.0, 'gammaJ', 3.3);   % G3 JONSWAP

%% ---- Sensor layouts --------------------------------------------------------------------------------
% Level 2C candidate grid (level2C_finalCheck.m): centre + r = [0.3 0.5 0.7 0.9] R x 8 angles.
% Stored as FRACTIONS of R and angles; converted to physical coordinates ONCE per scenario at the
% true R and then held fixed (never moved with a trial R during differentiation).
rFrac = [0.3 0.5 0.7 0.9];
th = (0:7) * pi / 4;
[TT, RR] = meshgrid(th, rFrac);
RR = RR.'; TT = TT.';
cfg.layout.gridFrac = [0 0; RR(:), TT(:)];         % 33 x 2 [r/R, theta], index order as level2C
cfg.layout.level2C.idx = [2 6 14 26];             % s2 s6 s14 s26
cfg.layout.level2C.frac = [0.3 0; 0.3 pi; 0.5 pi; 0.9 0];
cfg.layout.level2C.ref = 4;                       % s26 (position within the layout)

%% ---- Classification (exhaustive, applied in order) -------------------------------------------------
% s_max = max(sigma_lnbeta, sigma_lnR). Bias: the evaluator adds the bias VECTORS,
%   dtheta_tot = dtheta_W + dtheta_sys,  d_tot = sqrt(dtheta_tot' F dtheta_tot),
% so systematic errors can add or partly cancel; only the thresholds live here.
cfg.class.names = {'Identifiable', 'Marginal', 'Not identifiable'};
cfg.class.identifiable = struct('maxSigma', 0.01, 'maxBiasDistance', 0.5, 'maxCondF', 1e8);
cfg.class.marginal = struct('maxSigma', 0.05, 'maxBiasDistance', 1.0, 'maxCondF', 1e8);
cfg.class.referenceRecordN = cfg.acq.Nref;
cfg.class.toleranceLevels = struct('d', [0.5 1], 'relBeta', 0.01);   % first-crossing rule

%% ---- Noise cases (one-sided, per rad/s, vertical, basin band) --------------------------------------
names = {'datasheet', 'fused3', 'lsm6dsv16x', 'g3', 'icm42688p', 'lsm6ds3trc'};
Snn = [6.5e-8, 5.0e-7, 8.0e-7, 1.416e-6, 1.52e-6, 1.24e-5];
labels = {'Datasheet grade (65 ug/rtHz)', 'Three on-board IMUs fused', 'LSM6DSV16X', ...
    'G3 twin assumption (10% stress)', 'ICM-42688-P', 'LSM6DS3TR-C'};
cfg.noise.names = names;
cfg.noise.labels = labels;
cfg.noise.Snn = Snn;
cfg.noise.sigma = sqrt(pi * Snn / cfg.acq.dt);    % per-sample std at dt = 0.1 s
for i = 1:numel(names), cfg.noise.(names{i}) = Snn(i); end
cfg.noise.default = 'lsm6dsv16x';

%% ---- E1 validation gate targets (from U1 / L2048 / G3 result files) ----------------------------------
cfg.E1.sigmaF = [0.00462 0.00088];                % setting F at p0
cfg.E1.dPred = struct('A512', 3.68, 'B1024', 0.91, 'F2048', 0.26);
cfg.E1.sigmaTruthA = [0.00498 0.00095];
cfg.E1.sigmaTruthB = [0.00431 0.00082];
cfg.E1.tolSigma = 0.05;
cfg.E1.tolSigmaCovMC = 0.20;                      % per-element, MC standard error ~10%
cfg.E1.tolD = 0.10;                               % AMENDMENT 2026-10-01, see cfg.meta.changes
cfg.E1.truthA = [1.2 * cfg.p0.beta, cfg.p0.gamma, 0.97 * cfg.p0.R];
cfg.E1.truthB = [0.8 * cfg.p0.beta, cfg.p0.gamma, 1.03 * cfg.p0.R];

%% ---- R4: gamma uncertainty ---------------------------------------------------------------------------
cfg.R4.epsGamma = [-0.20 -0.10 -0.05 -0.02 -0.01 -0.005 0.005 0.01 0.02 0.05 0.10 0.20];
cfg.R4.nonlinearCheck = [-0.20 -0.05 0.05 0.20];
cfg.R4.noiseCases = names;                         % CLARIFICATION 2026-10-01, see cfg.meta.changes
cfg.R4.headlineNoise = 'lsm6dsv16x';               % the noise case the headline numbers and panel (a) use
cfg.R4.toleranceSearch = [-0.30 0.30];             % eps range searched for the first-crossing tolerances
cfg.R4.toleranceN = 1201;                          % points on that grid (0.05% spacing)

%% ---- D2: record duration -----------------------------------------------------------------------------
cfg.D2.durationMin = [10 20 30 60 109 180 240];
cfg.D2.L = [512 1024 2048 4096];
cfg.D2.minSegments = 8;                           % K = floor(2N/L) - 1 at 50% overlap
cfg.D2.noiseCases = names;

%% ---- R3: noise-PSD mis-specification -----------------------------------------------------------------
cfg.R3.epsNoise = [-1.0 -0.5 -0.3 -0.1 0.1 0.3 0.5];   % -1 = no correction
cfg.R3.colouredCase = 'lsm6ds3trc_lowfreq_shape';       % measured spectral shape, basin level
cfg.R3.noiseCases = names;

%% ---- D1: sensor count and layout -----------------------------------------------------------------------
cfg.D1.sensorCounts = [2 3 4];
cfg.D1.sensorCountsStretch = [5 6];
cfg.D1.maxRadiusFrac = 0.90;                      % edge clearance as Level 2C
cfg.D1.buoyDiameter = 0.10;                       % m, PROVISIONAL (SHARC buoy not final)
cfg.D1.minSpacing = cfg.D1.buoyDiameter;          % m, centre-to-centre: two buoys cannot overlap
cfg.D1.edgeClearance = cfg.D1.buoyDiameter / 2;   % m, buoy must sit wholly on the floe
cfg.D1.minSpacingSensitivity = [0.05 0.15];       % m, reported alongside, not a second optimum
cfg.D1.optimiseReference = true;
cfg.D1.criteria = struct('primary', 'A', 'secondary', {{'D', 'E'}});
cfg.D1.noiseCases = {'lsm6dsv16x', 'icm42688p', 'g3'};
cfg.D1.reportUnconstrained = true;
cfg.D1.compareGramian = true;                     % Level 2C lambda_min optimum per n_s
cfg.D1.worthItFraction = 0.20;                    % extra sensor must cut sigma_lnbeta by >= 20%
cfg.D1.floe = 'p0';

%% ---- A1: floe-property sweep ---------------------------------------------------------------------------
% (a) parameter-space ring and corners (basin, p0 gamma, nu = 0.3)
cfg.A1.ringBetaRatio = [0.25 0.5 2 4];
cfg.A1.ringRRatio = [0.5 0.75 1.5 2];
cfg.A1.cornerRatios = [0.5 0.5; 0.5 2; 2 0.5; 2 2];         % [beta/beta0, R/R0]
param = struct('name', {}, 'group', {}, 'p', {});
for b = cfg.A1.ringBetaRatio
    param(end + 1) = struct('name', sprintf('ring_beta_x%g', b), 'group', 'ring', ...
        'p', [b * cfg.p0.beta, cfg.p0.gamma, cfg.p0.R]); %#ok<AGROW>
end
for r = cfg.A1.ringRRatio
    param(end + 1) = struct('name', sprintf('ring_R_x%g', r), 'group', 'ring', ...
        'p', [cfg.p0.beta, cfg.p0.gamma, r * cfg.p0.R]); %#ok<AGROW>
end
for c = 1:size(cfg.A1.cornerRatios, 1)
    cr = cfg.A1.cornerRatios(c, :);
    param(end + 1) = struct('name', sprintf('corner_beta_x%g_R_x%g', cr), 'group', 'corner', ...
        'p', [cr(1) * cfg.p0.beta, cfg.p0.gamma, cr(2) * cfg.p0.R]); %#ok<AGROW>
end
cfg.A1.param = param;

% (b) physical classes. Material values are nominal; flagged ones are provisional.
PP = struct('E', 1.5e9, 'nu', 0.40, 'rho', 910);   % polypropylene: SG 0.91 (builder), E typical 1.3-1.8 GPa
H = cfg.physics.Hbasin; rhoW = cfg.physics.rhoWbasin; g = cfg.physics.g;
phys = struct('name', {}, 'class', {}, 'deepWater', {}, 'Rphys', {}, 'h', {}, 'E', {}, 'nu', {}, ...
    'rhoI', {}, 'H', {}, 'rhoW', {}, 'provisional', {}, 'note', {});
phys(end + 1) = mk('basin_p0', 'basin', false, cfg.p0.Rphys, NaN, NaN, cfg.physics.nuEMM, NaN, H, rhoW, false, ...
    'EMM reference; enters via D and m, not h and E');
phys(end + 1) = mk('pp3mm_R072', 'tank', false, 0.72, 0.003, PP.E, PP.nu, PP.rho, H, rhoW, false, ...
    'p0-radius PP analogue: draught matched to p0, D 4.02 vs 5.76 N m (beta 0.70 beta0)');
phys(end + 1) = mk('pp5mm_D40', 'tank', false, 0.20, 0.005, PP.E, PP.nu, PP.rho, H, rhoW, false, ...
    'one 5 mm PP sheet, 40 cm diameter (the model floe''s core sheet)');
phys(end + 1) = mk('pp10mm_D40', 'tank', false, 0.20, 0.010, PP.E, PP.nu, PP.rho, H, rhoW, false, ...
    'two bonded 5 mm PP sheets, 40 cm diameter');
phys(end + 1) = mk('pp5mm_D100', 'tank', false, 0.50, 0.005, PP.E, PP.nu, PP.rho, H, rhoW, false, ...
    '5 mm PP sheet, 1 m diameter');
hModel = 0.04;                                     % volume-equivalent flat plate of a 6 cm lens (2/3 of max)
for Eeff = [0.5e9 1.5e9 3.0e9]
    phys(end + 1) = mk(sprintf('modelFloe_E%.1fGPa', Eeff / 1e9), 'tank', false, 0.22, hModel, Eeff, PP.nu, ...
        PP.rho, H, rhoW, true, ['SHARC model floe: 44 cm, 6 cm max (elliptical), SG 0.91, E unknown; ' ...
        'h and R provisional until CAD and weighing (Fri 2 Oct); thin-plate validity marginal']); %#ok<AGROW>
end
MIZ = struct('E', 6e9, 'nu', 0.3, 'rho', 917);
for Rm = [10 25 50]
    phys(end + 1) = mk(sprintf('miz_R%d_h1', Rm), 'miz', true, Rm, 1.0, MIZ.E, MIZ.nu, MIZ.rho, Inf, ...
        cfg.physics.rhoWsea, false, 'screening class B2; regime/domain screening unless checks pass'); %#ok<AGROW>
end
phys(end + 1) = mk('miz_R25_h05_E3', 'miz', true, 25, 0.5, 3e9, 0.3, MIZ.rho, Inf, cfg.physics.rhoWsea, false, ...
    'screening class B2; regime/domain screening unless checks pass');
for i = 1:numel(phys)
    t = toEMM(phys(i), cfg, g);
    if i == 1, physOut = t; else, physOut(i) = t; end %#ok<AGROW>
end
cfg.A1.physical = physOut;

% Band rule. The floe-relative rule (same kR range as p0) was dropped at pre-registration: for the
% 20-22 cm floes it needs omega up to ~15 rad/s, alpha ~ 40, far outside the validated EMM range.
% Finite-depth (basin/tank) cases therefore use the validated basin band itself (alpha in [1.94, 13.2]).
% A1 then answers "which floes are identifiable under the fixed, validated basin excitation band?",
% an engineering question. It does NOT isolate floe properties: changing R changes kR as well.
cfg.A1.bandRule.finiteDepth = 'basin band, same alpha window as p0';
cfg.A1.bandRule.finiteDepthBand = cfg.band.basin;
cfg.A1.bandRule.mizPeriods = [6 12.5];            % s; effective depth set at 12.5 s (screenFloeRegime)
cfg.A1.requireTruncationCheck = true;
cfg.A1.requireNodeGridCheck = true;
cfg.A1.layoutFrac = cfg.layout.level2C.frac;      % sensors placed at Level 2C fractions of each true R
cfg.A1.noiseCase = 'lsm6dsv16x';

%% ---- R1: SNR axis ---------------------------------------------------------------------------------------
cfg.R1.Snn = logspace(-8, -4, 41);
cfg.R1.sea = cfg.sea;
cfg.R1.marked = names;

%% ---- R2: spectral shape (normalised to fixed incident variance inside the retained band) ------------
cfg.R2.omegaP = [3.5 4.25 5.0 6.0 7.5];
cfg.R2.gammaJ = [1.0 3.3 7.0];
cfg.R2.normalisation = 'fixed_inband_elevation_variance';   % isolates shape from SNR (R1 covers SNR)
cfg.R2.referenceSea = cfg.sea;                               % sets the in-band variance
cfg.R2.bimodal = struct('omegaP', [3.5 6.5], 'gammaJ', [7.0 3.3], 'energyFraction', [0.4 0.6]);
cfg.R2.noiseCase = 'lsm6dsv16x';
cfg.R2.checkL = [1024 2048 4096];

%% ---- V1: forward-model mismatch -------------------------------------------------------------------------
cfg.V1.truncationTruth = cfg.emm.truncationCheck;
cfg.V1.truncationInverse = cfg.emm.truncation;
cfg.V1.depthErrorFrac = [-0.05 -0.02 0.02 0.05];
cfg.V1.sensorRadialError = [-0.02 -0.01 0.01 0.02];       % m (divide by H for the EMM)

%% ---- V2: directional forcing ----------------------------------------------------------------------------
cfg.V2.angleDeg = [15 30 60 90];
cfg.V2.secondEnergyFraction = [0.05 0.10 0.25 0.50];
cfg.V2.coherenceDiagnostic = true;

%% ---- M1: selection RULES (scenarios resolved after D1, A1, R2) ------------------------------------------
cfg.M1.nRecords = 100;
cfg.M1.nRecordsExtended = 200;
cfg.M1.rules = {'D1 optimum at n_s = 3', 'A1 case nearest the Identifiable/Marginal boundary', ...
    'R2 narrowest or bimodal sea', 'p0 at LSM6DSV16X noise'};
cfg.M1.nonlinearBetaOffsets = [-0.20 0.20];
cfg.M1.accept = struct('empPred', [0.85 1.20], 'coverage95', [0.88 0.99], 'nonlinVsLinSigma', 0.2);

%% ---- H1: calibration and co-located transmissibility (hardware) ---------------------------------------
cfg.H1.gLocal = 9.796;                            % m/s^2, Cape Town
cfg.H1.holdTrim = 10;                             % s trimmed from each end of a static hold
cfg.H1.accept = struct('relGain', 0.0042, 'delay', 0.37e-3, 'refNoiseRatio', 1e-3);
end

%% ================================================================================================
function s = mk(name, cls, deep, Rphys, h, E, nu, rhoI, H, rhoW, prov, note)
s = struct('name', name, 'class', cls, 'deepWater', deep, 'Rphys', Rphys, 'h', h, 'E', E, 'nu', nu, ...
    'rhoI', rhoI, 'H', H, 'rhoW', rhoW, 'provisional', prov, 'note', note);
end

function s = toEMM(s, cfg, g)
% Physical -> EMM scaling for finite depth (screenFloeRegime.m formulas). Deep water: NaN here,
% computed with screenFloeRegime in the study script (needs the dispersion relation).
if strcmp(s.name, 'basin_p0')
    s.D = cfg.p0.D; s.m = cfg.p0.m;
else
    s.D = s.E * s.h^3 / (12 * (1 - s.nu^2));
    s.m = s.rhoI * s.h;
end
s.draught = s.m / s.rhoW;
s.lf = (s.D / (s.rhoW * g))^(1/4);
s.RoverLf = s.Rphys / s.lf;
if s.deepWater
    s.p = [NaN NaN NaN];
else
    s.p = [s.D / (s.rhoW * g * s.H^4), s.draught / s.H, s.Rphys / s.H];
end
end