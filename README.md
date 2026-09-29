# SHARC Hydroelastic Forward Model, SSI, Inverse Mapping and IMU Characterisation

MATLAB research code for the SHARC buoy inversion pipeline. The repository now contains
these connected stages:

1. a hydroelastic Forward Model for a single floating circular elastic plate, based on
   Montiel's (2012) eigenfunction matching method (EMM);
2. an output-only Stochastic Subspace Identification (SSI) pipeline for recovering
   frequencies, damping and spatial mode information from multi-sensor time series;
3. an Inverse Mapping stage that extracts modal-response features from the EMM and
   recovers the parameter vector `(beta, gamma, R)` through bounded nonlinear inversion;
4. an **output-only transmissibility inversion** (the "SSI bridge" and its successor): sensor
   FRFs, an information analysis that motivated the pivot from SSI poles to
   transmissibilities `T_j = H_j/H_r`, a validated Welch transmissibility estimator, and a
   closed synthetic inverse loop for `theta = [ln beta, ln R]` with `gamma` fixed
   (see *SSI bridge and the pivot to transmissibilities* and *Output-only twin* below);
5. **Mission Physics**: dimensional screening of tank and Antarctic MIZ floe-wave regimes
   before any EMM run; and
6. **IMU Characterisation**: parsing, timing audit, noise PSD/coherence and Allan deviation
   for the real SHARC board, validated on synthetic records with known answers.

Stages 4-6 were added from 25 September 2026 onwards; their sections follow the original
Inverse Mapping section, which is kept unchanged as the record of the EMM-feature route.

The Forward Model remains the physical foundation of the repository. Given plate
parameters and an incident wave frequency it returns the complex deflection `eta(r,theta)`,
and the JONSWAP extension superposes those single-frequency responses into synthetic
time series. SSI and the inverse-mapping code are built around that model without
modifying its core equations.

Equation numbers in the Forward Model section refer to Montiel's 2012 Otago thesis,
Chapter 2 and Appendices A.1/B.1, unless stated otherwise.

## Repository layout

```
test_repo/
  Forward Model/
    (core EMM pipeline functions, described below)
    Animation/
      precomputeDeflectionData.m
      evaluateDeflection.m
      plotDeflectionSurface.m
      plotDeflectionTimeSeries.m
      runDeflectionAnimation.m
    Unit-Tests/
      (one tester per core function)
    Validation Experiments/
      singleFrequencyModalProjectionTest.m
    JONSWAP/
      (spectral extension functions)
      Unit-Tests-JONSWAP/
        (one tester per JONSWAP function)

  SSI/
    (SSI-DATA, modal tracking and Gramian sensor-placement tools)
    Unit Tests SSI/
      (unit/integration tests for SSI and Gramian utilities)
    Validation Experiments/
      (Forward-Model-to-SSI validation ladder, Levels 1-3B)

  Inverse Mapping/
    (modal feature extraction, Jacobian analysis and anchor prediction)
    Inverse Mapping Unit Tests/
      (feature/Jacobian/tracker tests)
    Parameter Inversion/
      inverseResidual.m
      solveInverse.m
      probeValidRegion.m
      Tests/
      EMM Twin/
        runProbeValidRegionEMM.m
        runEMMTwinInversion.m
        runFeaturePerturbationStudy.m
      Output-Only Twin/            (renamed from SSI Twin/ on 28 Sep)
        runSensorFRFReconnaissance.m, runSensorFRFSensitivity.m,
        runTransmissibilityInformation.m, runOutputOnlyInformationBound.m,
        runTransmissibilityEstimatorTest.m, checkSegmentSweepCorrelation.m,
        runTransmissibilityCovariance.m, runL2048BiasVarianceCheck.m,
        runTransmissibilityTwinInversion.m, runMeasurementTolerance.m
      Results/
        (saved probe, twin, perturbation, information, estimator, covariance,
         inversion and tolerance outputs; twinCache/ holds cached EMM node solves)
    computeSensorFRF.m, featureInformation.m, transmissibilityFeatureVector.m
      (new top-level Inverse Mapping functions, tests in Inverse Mapping Unit Tests/)

  Transmissibility/
    extractTransmissibility.m, stackTransmissibility.m, unstackTransmissibility.m,
    synthesiseTwinRecords.m
    Unit Tests Transmissibility/

  Mission Physics/
    screenFloeRegime.m, waveNumber.m, runFloeRegimeScreen.m
    Unit Tests/
    Results/

  IMU Characterisation/
    Processing/
      loadImuRecord.m, checkImuTimestamps.m, synthesiseImuNoise.m, reportImuRecord.m,
      runShellRateSweep.m, alignImuSensors.m, computeNoisePSD.m, allanDeviation.m,
      runStaticNoiseCharacterisation.m
      Unit Tests/
    Data/
      Raw/                         (untouched PuTTY logs from the SHARC board)
    Results/
```

`SSI/` also gained `stochasticSensorSynthesis.m` (with its tester) and
`Validation Experiments/level4A_broadbandKnownOscillator.m`.

## Quickstart

```matlab
addpath('Forward Model');

alpha = 10; beta = 1e-2; gamma = 0.1; R = 1; nu = 0.3;
M = 20; P = 10; N = 15;

eta = deflection(0.5, pi/3, alpha, beta, gamma, R, nu, M, P, N);
```

or, for a full animated plot:

```matlab
cd 'Forward Model/Animation'
runDeflectionAnimation
```

## Non-dimensional parameters

Lengths are non-dimensionalised by the water depth H, time by sqrt(H/g). All
functions take non-dimensional inputs:

| Symbol | Meaning | Definition |
|---|---|---|
| `alpha` | non-dim frequency | H·ω²/g |
| `beta`  | non-dim flexural rigidity | D/(ρ₀gH⁴), D = Eh³/12(1-ν²) |
| `gamma` | non-dim draught | d/H, d = (ρ/ρ₀)h (Archimedes) |
| `R`     | non-dim plate radius | R_physical/H |
| `nu`    | Poisson's ratio | fixed at 0.3 throughout |
| `M`     | vertical mode truncation | tower modes 0..M |
| `P`     | Gegenbauer mode truncation | modes 0..P |
| `N`     | angular mode truncation | n = -N..N |

**β is extremely sensitive to the assumed water depth H (∝H⁻⁴).** There is no
single "typical" value - it must be computed from the actual physical scenario.
See `Unit-Tests` discussion notes for worked examples spanning realistic sea-ice
and lab-experiment parameter ranges.

## Truncation guidance

- `N ≈ 10` is adequate at `α=10`; `N ≈ 20` needed by `α=20` (roughly `N` should
  scale with `κR`, the acoustic size of the disc).
- `M` requirements depend strongly on `β`: a softer plate (smaller β) excites more
  spatial/flexural detail and needs more tower modes to converge. Always verify
  convergence (e.g. compare `M` and `2M`) before trusting a new parameter regime -
  do not assume a truncation that worked for one `(α,β,γ,R)` transfers to another.
- `P < M` is required for `D0g`/`Dg` to remain well-conditioned.

---

## Core pipeline functions (`Forward Model/`)

Functions are listed in pipeline order. Each depends only on functions above it.

### Step 0 - Roots and vertical basis

**`dispersionFunction.m`** - `[D, Dp] = dispersionFunction(xi, alpha, beta, gamma)`
The single entire function whose roots are the vertical wavenumbers for *both*
the open-water and disc-covered dispersion relations (`beta=0, gamma=0` gives
open water). Returns the function value `D` and, optionally, its derivative `Dp`
(closed form, verified against finite differences). Used internally by
`dispersionRoots.m` and `bennettsCoefficients.m` - not normally called directly.

**`dispersionRoots.m`** - `xi = dispersionRoots(alpha, beta, gamma, M)`
Solves for the vertical wavenumbers. Pass `beta=0, gamma=0` for open water,
returning `[k0; k1; ...; kM]` (M+1 roots, k0 purely imaginary/travelling, rest
real/evanescent). Pass the actual `beta, gamma` for the disc-covered case,
returning `[kappa_{-2}; kappa_{-1}; kappa0; kappa1; ...; kappaM]` (M+3 roots,
including the complex-conjugate damped pair). This is usually the *first* call
in any script using this model.

**`verticalEigenfunction.m`** - `phi = verticalEigenfunction(z, xi, gamma)`
The normalised vertical mode shape, `cos(xi(z+1))/cos(xi(1-gamma))`. One formula
serves both open water (`gamma=0`) and disc-covered regions.

**`edgeDerivative.m`** - `dphi = edgeDerivative(xi, gamma)`
`φ'_m(-γ)`, the vertical derivative of the mode shape at the plate's submerged
edge - used throughout the matching system and in the final deflection sum.

### Step 1 - Bennetts (damped mode) reduction

**`bennettsCoefficients.m`** - `V = bennettsCoefficients(alpha, beta, gamma, kappa)`
Builds the `(M+1)×2` matrix `V` expressing the two damped-mode eigenfunctions
as a linear combination of the ordinary tower modes (Bennetts et al. 2007
reduction). Column 1 reconstructs `φ₋₂`, column 2 reconstructs `φ₋₁` (matches
`kappa`'s slot ordering, not the "natural" j-index ordering - see in-file
comments). Carries a confirmed sign correction relative to the thesis as
literally printed (traced to Bennetts 2007 eq. 4.4).

### Step 2 - Edge map

**`edgeConditionRows.m`** - `[row1, row2] = edgeConditionRows(kappaVal, R, n, nu)`
The two bracketed expressions from the free-edge conditions (eq. A.3a/A.3b),
evaluated at a single vertical wavenumber. Shared building block for `edgeMap.m`
and the end-to-end residual check in `scatteringMatricesTester.m`.

**`edgeMap.m`** - `Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa)`
Uses the plate's free-edge conditions (zero bending moment, zero shear) to
express the two damped-mode amplitudes in terms of the tower amplitudes:
`Ae_n = Ln * A_n`. Reduces the per-n unknown count from M+3 to M+1. Returns a
`2×(M+1)` matrix. One value of `n` at a time.

### Step 3 - Radial matrices

**`radialMatrix.m`** - `En = radialMatrix(n, r, R, kappa, V, Ln)`
The corrected radial matrix in the disc-covered region, `E_n(r) = 𝓔_n(r) +
V·𝓔̃_n(r)·Ln` (eq. A.5) - folds the damped-mode contribution into the tower via
`V` and `Ln`. Satisfies the exact identity `E_n(R) = I + V·Ln`.

**`scaledBesselDerivative.m`** - `val = scaledBesselDerivative(besselType, n, K, r, R)`
`d/dr[X̂_n(Kr)]` for the scaled Bessel function `X̂_n(Kr)=X_n(Kr)/X_n(KR)`, for
`X=I` (bounded) or `X=K` (radiating). Shared building block for the two
functions below.

**`openWaterRadialDerivatives.m`** - `[EIprime, E0prime] = openWaterRadialDerivatives(n, R, k)`
`d/dr` of the open-water diagonal radial matrices `E^(I)_n`, `E^(0)_n`, evaluated
at `r=R`. Diagonal - open water has no damped-mode mixing.

**`radialMatrixDerivative.m`** - `Enprime = radialMatrixDerivative(n, R, kappa, V, Ln)`
`d/dr[E_n(r)]` at `r=R` for the disc-covered region - the derivative counterpart
of `radialMatrix.m`, same `V`/`Ln` structure.

### Step 4 - Matching system

**`innerProductMatrix.m`** - `Dmat = innerProductMatrix(kappaTower, gamma)`
Inner-product matrix of tower modes, `[D]_ij = ∫φ_iφ_j dz`. Pass `gamma=0` with
the open-water tower for `D0` (comes out exactly diagonal - orthogonality);
pass the actual `gamma` with the disc-covered tower for `D` (not diagonal).
Closed form independently re-derived and confirmed (corrects a transcription
error in the printed thesis).

**`gegenbauerCrossMatrix.m`** - `Dg = gegenbauerCrossMatrix(rootVector, gammaDraught, P, normScale)`
Cross inner-product between a vertical mode basis and the Gegenbauer singular
basis (used to resolve the velocity singularity at the plate's submerged
corner). Pass the open-water tower with `normScale=1` for `D0g`; the
disc-covered tower with `normScale=1-gamma` for `Dg`. Returns `(M+1)×(P+1)`.

**`gegenbauerBasisFunction.m`** - `val = gegenbauerBasisFunction(p, gammaDraught, y)`
Validation-only helper evaluating the actual weighted Gegenbauer basis function
(eq. 8), used only by `gegenbauerCrossMatrixTester.m` for an independent
quadrature check. Requires the Symbolic Math Toolbox (`gegenbauerC`). Not
called anywhere in the main pipeline.

**`matchingMatrices.m`** - `[MIn, M0, Mmat] = matchingMatrices(EIprime, E0prime, Enprime, D0, D, D0g, Dg)`
Solves the velocity-continuity matching equations (B.1a/B.1b) for the amplitude
unknowns in terms of the incident wave and the interface unknown `U^(n)`.
Includes `rcond` checks that warn if any required inversion is close to
singular.

**`scatteringMatrices.m`** - `[MU, Sn, Sn0] = scatteringMatrices(D0g, Dg, En_R, MIn, M0, Mmat)`
Solves the pressure-continuity condition (B.1c) for the interface unknown
`U^(n)`, then assembles the actual scattering matrices: `An = Sn*AIn`,
`A0n = Sn0*AIn`. This is the payoff of the entire matching system.

### Steps 5-7 - Forcing and deflection

**`incidentAmplitude.m`** - `AIn = incidentAmplitude(alpha, n, k0, R, M)`
The incident plane-wave amplitude vector (eq. 2.22) - a unit-amplitude wave
travelling in the `-x` direction. Only the travelling (m=0) tower mode is
forced; every other tower mode is zero.

**`deflection.m`** - `eta = deflection(r, theta, alpha, beta, gamma, R, nu, M, P, N)`
The reference (simple, unoptimised) implementation of the full pipeline:
builds every piece above from scratch and returns the plate deflection at a
single point `(r,theta)`. Good for one-off checks; expensive if called
repeatedly for many points (see `Animation/` for the fast path).

---

## Animation (`Forward Model/Animation/`)

`deflection.m` rebuilds the entire pipeline on every call - fine for a single
point, wasteful for a plot needing thousands of points. These four functions
split the expensive, `(r,θ)`-independent part (Steps 0-6) from the cheap part
that actually depends on position (Step 7), so a grid or animation only pays
the expensive cost once.

**`precomputeDeflectionData.m`** - `data = precomputeDeflectionData(alpha, beta, gamma, R, nu, M, P, N)`
Runs Steps 0-6 once, for every angular order `n=-N..N`, and packages the
result into a struct for repeated fast evaluation.

**`evaluateDeflection.m`** - `eta = evaluateDeflection(data, r, theta)`
Cheap per-point (or per-grid, `r`/`theta` may be same-sized arrays) evaluation
from precomputed data. Regression-tested against `deflection.m` - exact
agreement.

**`plotDeflectionSurface.m`** - `plotDeflectionSurface(data, nr, ntheta, nFrames)`
Animated 3D surface of the deflection over one full non-dimensional wave
period (`tau` in `[0, 2*pi)`). Uses a diverging blue-white-red colormap
centred at `eta=0`, appropriate for a signed quantity.

**`plotDeflectionTimeSeries.m`** - `plotDeflectionTimeSeries(data, points, nTau)`
1D time traces of `eta` vs. non-dimensional phase at a list of chosen
`[r, theta]` points.

**`runDeflectionAnimation.m`**
Driver script - set parameters, precompute once, render both plots. Start
here for a first look at any new parameter set.

---

## Unit-Tests (`Forward Model/Unit-Tests/`)

One tester per core function, following the same pattern throughout: check
against a closed-form identity where one exists, cross-check against
independent numerical methods (finite differences, quadrature, contour
integration) where it doesn't, and check convergence (does the answer stop
changing as truncation increases) wherever a series or expansion is involved.
Every tester begins with `addpath('..')` (and, where needed, `addpath('../Animation')`)
to resolve the core functions from this subfolder.

Run any individual tester directly by name from within `Unit-Tests/`, e.g.:

```matlab
dispersionTester
```

Key testers and what they confirm:

- **`dispersionTester.m`** - roots satisfy the dispersion relation to ~1e-10;
  validated across 7 `(alpha,beta,gamma)` cases including a hard basin-of-
  attraction stress test for the damped-pair Newton solver.
- **`edgeDerivativeTest.m`** - matches finite-difference derivative, real and
  complex arguments.
- **`bennettsCoefficientsTester.m`** - reconstructs `phi_{-j}(z)` from the tower
  at both the edge and interior points; confirms convergence rate and the
  sign correction.
- **`edgeMapTester.m`** - correct shape/finiteness across multiple `n` and
  parameter sets, including `n=0`.
- **`radialMatrixTester.m`** - confirms the exact identity `En(R)=I+V*Ln` to
  machine precision.
- **`radialDerivativeTester.m`** - all three derivative functions checked
  against finite differences (~1e-9 to 1e-10).
- **`innerProductMatrixTester.m`**, **`gegenbauerCrossMatrixTester.m`** -
  closed forms checked against numerical quadrature; `D0` confirmed exactly
  diagonal (orthogonality), `D` confirmed non-diagonal.
- **`matchingMatricesTester.m`** - shape/finiteness, `rcond` guards exercised.
- **`scatteringMatricesTester.m`** - **the highest-value single test in the
  repository.** Forward-propagates a real incident wave through the entire
  matching system, then plugs the recovered amplitudes back into the
  *original* free-edge conditions (A.3a/A.3b) via `edgeConditionRows.m`.
  Residuals came back at machine precision (~1e-16), confirming the entire
  chain (Steps 0-4) is jointly self-consistent.
- **`deflectionTester.m`** - confirms `eta(0,theta)` is exactly
  theta-independent (a hard structural consequence of `I_n(0)=0` for `n≠0`,
  not an approximation) and that the deflection converges in `N` by
  `N≈10` at `alpha=10`, matching the stated truncation guidance.
- **`deflectionRefactorTester.m`** - confirms the fast `Animation/` path and
  the reference `deflection.m` agree exactly.

## JONSWAP extension (Forward Model/JONSWAP/)

Extends the model from a single monochromatic wave to a realistic irregular
sea state described by a JONSWAP spectrum. Because the underlying model is
linear, the response to a spectrum is exactly the superposition of the
responses to many discrete monochromatic components: no change was made to
any core pipeline function above. Every file in this folder is orchestration
around the existing, unmodified pipeline.

**`waveSpectrum.m`** - third-party function (Thor I. Fossen, MSS toolbox),
used as-is. Computes several standard wave spectra; only case 7 (JONSWAP via
Hs, peak frequency w0, peakedness gamma) is used here. Verified independently
that it correctly recovers its own input Hs from the zeroth moment of the
returned spectrum. Only cases 1-7 were checked as dependency-free; case 8
(Torsethaugen) calls a separate function, torsetSpectrum.m, not included and
not needed for this project.

**`discretizeSpectrum.m`** - `[omega, a, epsilon] = discretizeSpectrum(Hs, w0, gammaJONSWAP, wMin, wMax, numBins)`
Converts a continuous JONSWAP spectrum into a finite set of discrete
frequency components via the standard random-phase/amplitude method:
amplitude a_i = sqrt(2*S(omega_i)*domega), so that each component's
time-averaged variance (a_i^2/2) matches the variance the continuous
spectrum assigns to that frequency band; phases epsilon_i drawn independently
and uniformly on [0, 2*pi). Validated: recovered Hs converges cleanly as
numBins increases, and matches an independent continuous-integral
calculation of the same spectrum to 4 decimal places.

**`precomputeSpectralData.m`** - `specData = precomputeSpectralData(Hs, w0, gammaJONSWAP, wMin, wMax, numBins, H, beta, gamma, R, nu, M, P, N)`
Discretizes the spectrum, converts each bin's angular frequency into its own
non-dim alpha (alpha_i = H*omega_i^2/g), then calls the existing
precomputeDeflectionData.m once per bin, unmodified. Returns a struct with
the per-bin frequencies, amplitudes, phases, alphas, and precomputed
amplitude data. Validated: exact-match regression against calling
precomputeDeflectionData.m directly at one bin's alpha; correct shapes and
finiteness confirmed. Timing: roughly 1s/bin at M=50,P=10,N=10 (measured, not
a guarantee at other truncation levels).

**`evaluateSpectralDeflection.m`** - `zeta = evaluateSpectralDeflection(specData, r, theta, t)`
The cheap evaluation step: sums each bin's contribution,
a_i * Re{eta_i(r,theta) * exp(i*(omega_i*t + epsilon_i))}, using the
existing evaluateDeflection.m per bin. Returns the actual physical deflection
in metres (see note on units below), at real physical time t in seconds (not
the non-dim tau used in the single-frequency Animation/ tools). Warns if r
exceeds the plate radius. Validated: exact single-bin collapse onto
deflection.m called directly; exact multi-bin sum regression against a
hand-computed reference; exact theta-independence at r=0 carried through
the whole time series (a hard structural consequence of I_n(0)=0 for n!=0,
holding per-bin and therefore for the sum); r>R warning confirmed to fire.

**Units note:** eta (from the core pipeline) is non-dimensional, the
response to a unit non-dim incident amplitude. zeta (from this extension) is
in real physical metres, since each bin's amplitude a_i is already a
physical quantity in metres (from waveSpectrum's own units) and the two
factors of H that would otherwise appear (non-dimensionalising a_i, then
redimensionalising the result) cancel exactly, given the model's linearity.

**`plotSpectralDeflectionSurface.m`** - `plotSpectralDeflectionSurface(specData, nr, ntheta, nFrames, tMax, exaggeration, saveFilename)`
Animated 3D surface of the JONSWAP-forced deflection, in the same x/y style
as the single-frequency plotDeflectionSurface.m (non-dim r,theta axes).
Colour always reflects the true physical deflection in mm; the plotted
height is deliberately exaggerated (factor stated in the title) since real
deflections are typically too small, relative to the plate radius, to see
under true-scale axis proportions. EXAGGERATION and SAVEFILENAME are
optional: omit or pass [] for EXAGGERATION to auto-choose it; provide
SAVEFILENAME to save an MP4 instead of playing live, with no change to the
underlying computation either way.

**`runSpectralDeflectionAnimation.m`** - driver script. Set the spectrum and
floe parameters, call precomputeSpectralData.m once, then
plotSpectralDeflectionSurface.m. Start here for a first look at a new
JONSWAP scenario.

### Fixed truncation, chosen deliberately

M, P, N are held fixed across every frequency bin in a spectrum, rather than
adapted per bin, for simplicity. The current setting, M=50,P=10,N=10, was
confirmed (via a doubling check against M=70,P=15,N=15, diff of order 1e-6)
adequate across the wave-basin-scale test range used so far, alpha in
roughly [1.7, 13.9]. This was tested only at wave-basin scale (H=1.88m,
R=0.383) deliberately, to validate the new spectral machinery against
already-trusted single-frequency results before attempting a realistic
ocean-scale scenario. Realistic SCALE-deployment depths and periods push
alpha far beyond this range (into the tens or hundreds, since alpha is
structurally the same quantity as the deep-water parameter kH) and would
need this truncation check re-run before trusting results at that scale.

### Unit-Tests-JONSWAP (Forward Model/JONSWAP/Unit-Tests-JONSWAP/)

Same validation philosophy as the core Unit-Tests folder. Files anchor their
own paths using `fileparts(mfilename('fullpath'))` rather than plain
relative addpath calls, since this folder sits two levels below Forward
Model rather than one.

- **`discretizeSpectrumTester.m`** - variance/Hs convergence as numBins
  increases; cross-check against independent continuous integration; phase
  range check.
- **`precomputeSpectralDataTester.m`** - shape/finiteness checks; exact-match
  regression against a direct precomputeDeflectionData.m call; timing
  measurement.
- **`evaluateSpectralDeflectionTester.m`** (two parts) - single-bin collapse
  check against deflection.m; incident-wave-only statistical check (5
  trials, confirming Hs recovery is unbiased across repeated random-phase
  draws, not just a single lucky or unlucky run); genuine multi-bin sum
  regression against a hand-computed reference; theta-independence at r=0
  carried through the time domain; r>R warning firing check.

## SSI extension (SSI/)

Phase 2 of the SHARC pipeline: given a time series (real IMU data, or a
synthetic signal from the Forward Model above), reconstructs a discrete-time
state-space model and extracts candidate natural frequencies, damping ratios,
and mode shapes via Stochastic Subspace Identification (SSI-DATA). This is an
output-only method - no knowledge of the forcing is required or used, which
matters because the real wave forcing on a deployed buoy is never directly
measured. Nothing in this folder modifies any Forward Model file; the two
connect only through `syntheticSensorData.m`.

**Core algorithm, in pipeline order:**

**`buildHankelMatrix.m`** - `[Yp_ref, Yf, H] = buildHankelMatrix(y, i, refIdx)`
Stacks a multi-channel time series into the "past" and "future" block Hankel
matrices SSI is built on. `refIdx` optionally restricts which channels serve
as references (trims the past block only; every channel still gets its shape
reconstructed via the future block). Validated: exact match against a
hand-computed small case, reference-channel restriction confirmed, input
validation confirmed to fire correctly.

**`hankelProjection.m`** - `[P, Q1, R] = hankelProjection(Yp_ref, Yf)`
Projects the future block onto the past block's row space via QR
factorization (numerically safer than the direct `pinv`-based formula).
Validated: agrees with the direct formula to ~1e-13, QR factorization
structural checks and projection idempotency confirmed to ~1e-14.

**`extractStateSpace.m`** - `[A, C, singularValues] = extractStateSpace(P, l, n)`
SVD of the projection, truncated to model order `n`, recovering the state
matrix `A` (via the shift-invariance of the extended observability matrix)
and output matrix `C`. Validated: exact recovery of a known noise-free
2-state oscillator; a genuine twin test against real Forward Model output
recovered the true forcing frequency to 0.0000% error.

**`modalParameters.m`** - `[omega, zeta, Phi, lambda_c, isOscillatory] = modalParameters(A, C, dt)`
Converts `A`'s eigenvalues to natural frequency, damping ratio, and mode
shape (`Phi = C*Psi`). Validated against a known system (exact recovery) and
an independent algebraic check bypassing SSI entirely. **Carries a fixed
bug**: a real, negative discrete eigenvalue produces a spurious `+i*pi/dt`
term after `log()` (MATLAB's complex branch cut), indistinguishable from a
genuine oscillatory pole with no conjugate partner unless checked for
explicitly - `isOscillatory` tests `imag(lambda)` (before the branch cut),
not `imag(lambda_c)` (after it), specifically to avoid this. Caught by a
hard invariant check in `sweepModelOrdersTester.m` (see below), not by any
plausibility check.

**`computeMAC.m`** - `MACval = computeMAC(Phi_a, Phi_b)`
Modal Assurance Criterion, `MAC = |phi_a^H phi_b|^2 / [(phi_a^H phi_a)(phi_b^H phi_b)]`,
scale-and-phase-invariant similarity between two (possibly complex) mode
shapes. Validated: self-MAC exactly 1, orthogonal-MAC exactly 0, exact
complex-scale invariance confirmed, and a squared-numerator bug found in a
third-party implementation during this project explicitly confirmed absent
here. **Important limitation, proven not just stated**: with a single sensor
channel, MAC between any two nonzero scalars is always exactly 1 - it
provides zero discrimination until at least two channels are available.

**`sweepModelOrders.m`** - `results = sweepModelOrders(P, l, dt, nMax)`
Runs `extractStateSpace` + `modalParameters` at every order 1 to `nMax`,
keeping one representative per complex-conjugate pole pair (a real
oscillatory mode always appears as a pair; keeping both would double-count
every mode). Validated on a genuine multi-mode system against the hard
invariant that a real system of order `n` has at most `floor(n/2)` conjugate
pairs - this check caught the `modalParameters.m` branch-cut bug directly.

**`classifyPoleStability.m`** - `match = classifyPoleStability(prevResults, currResults, freqTol, dampTol, macTol)`
Matches poles at order `n` against order `n-1` via a hard frequency gate,
then damping and MAC as softer criteria, assigned via greedy one-to-one
matching (ranked by MAC) so no pole is double-claimed. Assigns a
hierarchical class: 0 (no match), 1 (frequency only), 2 (+damping), 3
(+MAC, "fully stable"). Validated: full class hierarchy exercised with exact
expected values, hard frequency gate confirmed, best-MAC-wins assignment
confirmed (not just non-duplication).

**`buildModalBranches.m`** - `branches = buildModalBranches(results, matches)`
Threads per-transition matches into full branches: one pole tracked across
consecutive orders. Deliberately does not bridge gaps - a pole failing to
match at one order ends its branch; reappearance starts a new one, with no
memory of the earlier branch. A branch continues through class 1/2 matches,
not only class 3 (branch membership = "is there a correspondence"; the
`class` field records the *quality* of that correspondence separately).
Validated: hand-constructed known branch structure (exact match), confirmed
against the noise-free 3-mode integration test that all 3 true modes form
single unbroken branches spanning their full valid order range.

**`computeBranchPersistence.m`** - `persistence = computeBranchPersistence(branches, windowSize, minStableInWindow)`
Sliding-window persistence: a branch is persistent if *at least one* window
of `windowSize` consecutive transitions contains at least `minStableInWindow`
class-3 links - deliberately not requiring an unbroken run throughout,
since that was shown empirically to under-report genuine, repeatedly-stable
modes under noise. `isPersistent` means "convincing evidence exists
somewhere in this branch's history," not "stable for its entire lifetime."
Validated on hand-constructed cases (isolated non-class-3 transitions
tolerated, sparse hits correctly rejected, too-short branches handled
without error).

**`syntheticSensorData.m`** - `y = syntheticSensorData(specData, sensorLocations, tVec)`
The bridge to the Forward Model: calls the existing, unmodified
`evaluateSpectralDeflection.m` once per sensor location and stacks the
results as rows, producing exactly the `l x n` matrix `buildHankelMatrix.m`
expects. Pure orchestration, no new physics. Validated: exact regression
against direct calls (single- and multi-sensor), distinct sensor locations
confirmed to give distinct signals, clean handoff into `buildHankelMatrix.m`
confirmed.

### Gramian sensor-placement utilities

The SSI folder now also contains a finite-horizon observability-based sensor-placement
branch, developed for Level 2C. These functions do **not** estimate a model from data;
they construct a known two-tone state-space model directly from the Forward Model so
candidate sensor layouts can be compared without circularly running SSI first.

**`buildKnownOscillatorA.m`** - constructs the exact real `4x4` block-diagonal rotation
matrix for two undamped frequencies, with state ordering
`[cos(omega1*t); sin(omega1*t); cos(omega2*t); sin(omega2*t)]`.

**`buildCandidateCrows.m`** - maps a candidate sensor's two complex, normalised Forward
Model responses into the matching real `1x4` output row `C_j`. Normalisation is kept
outside this function deliberately so Level 2C can separate spatial geometry from
amplitude imbalance.

**`finiteObservabilityGramian.m`** - finite-horizon discrete-time observability Gramian
`W_o = sum_{k=0}^{N-1} (A^k)' C' C A^k`. The finite-horizon form is required because
the known oscillator is undamped and therefore does not admit the usual convergent
infinite-horizon Lyapunov Gramian.

**`gramianMetrics.m`** - computes `logDet`, `trace`, `lambdaMin`, condition number and
Gramian eigenvalues, with explicit numerical-rank handling.

**`greedySensorSelection.m`** - generic forward greedy selection under any scalar
Gramian objective. It evaluates each candidate **together with the already selected
sensors**, rather than ranking sensors independently.

**`exhaustiveSensorSelection.m`** - exhaustive fixed-cardinality sensor search used to
check whether the greedy Level 2C selections are actually globally optimal on the
candidate grid.

### Two findings worth knowing before using this code

**`freqTol` has no single safe default.** With one sensor channel, MAC
provides no discrimination (see above), so widening `freqTol` enough to
tolerate realistic noise can wrongly merge two genuinely distinct, closely
spaced modes - demonstrated directly on a hand-built two-mode system. A
second channel with genuinely different modal weighting between the two
modes resolved this cleanly at the identical tolerance. The appropriate
value depends on sensor count and configuration, not a fixed constant.

**Branch fragmentation under noise is a real, diagnosed behaviour, not a
bug.** `buildModalBranches.m`'s no-gap-bridging design means a genuine mode
can fragment into several short branches under noise if its frequency
estimate drifts more than `freqTol` allows between some consecutive orders -
confirmed via a controlled `freqTol` sweep on the same dataset, which
reconnected the fragments into single stable branches once the gate was
widened appropriately, with no change to any other function. Deliberately
not "fixed" by adding gap-bridging to `buildModalBranches.m`, since the
fragmentation was fully explained by frequency tolerance, not a deficiency
in branch reconstruction itself.

## Unit Tests SSI (`SSI/Unit Tests SSI/`)

Same validation philosophy as the core Unit-Tests folder: one tester per
function, checking exact identities where they exist, hand-computed cases
where they don't, and - specifically for this folder - hard structural
invariants (e.g. the `floor(n/2)` conjugate-pair bound) rather than only
plausibility checks, since that discipline is what caught the
`modalParameters.m` branch-cut bug.

## Validation Experiments (`SSI/Validation Experiments/`)

Experiments connecting the SSI algorithm to genuine Forward Model physics,
answering "does SSI work" rather than "is this one function correct" - not
unit tests, but deliberate scientific validation scripts, run in sequence:

- **`level1MonochromaticTest.m`** - a single known incident frequency, 4
  sensors (centre, mid-radius, two near-edge points). **Important
  interpretive caveat, stated explicitly in the file**: `evaluateSpectralDeflection`
  produces a forced steady-state response, not a free-vibration record, so
  recovering the true frequency confirms SSI can extract it from real
  spatial physics - it does NOT confirm that frequency is a natural
  hydroelastic mode. Result: frequency recovered to ~1e-15 relative error at
  every model order tested, damping at the numerical noise floor (~1e-14),
  exactly as predicted for a forced tone with no damping. A negative control
  (checking the same result against a deliberately wrong target frequency)
  correctly failed as it should.
- **`level2MultiFrequencyTest.m`** - two known, deliberately close
  frequencies (2.70, 2.80 rad/s), built by hand (bypassing JONSWAP's
  automatic binning) for exact control over inputs. Computes the
  Forward-Model-only cross-frequency MAC *before* involving SSI at all,
  since this sets the ceiling on what any method could possibly achieve
  spatially. Result: both frequencies recovered exactly; SSI's cross-MAC
  matched the Forward-Model-only cross-MAC to 4 decimal places. The
  intended success criterion is `MAC_SSI,cross ≈ MAC_FM,cross` (SSI
  faithfully reproducing whatever the true physical similarity is),
  not an arbitrary fixed threshold.
- **`level2A_frequencyRegimeStudy.m`** - fixed sensors, varying frequency
  pair. Found that spatial distinguishability between two forcing
  frequencies is strongly frequency-dependent and non-monotonic (not simply
  a function of frequency separation) - across 5 pairs tested, SSI matched
  the true Forward-Model-only cross-MAC to 4 decimal places every time.
- **`level2B_sensorPlacementStudy.m`** - fixed frequencies (2.70, 6.00 rad/s,
  chosen from 2A as genuinely distinguishable), varying sensor placement
  across 4 configurations. Found that placement matters enormously and not
  in the naively expected direction: a radially-spread configuration
  performed far worse (MAC 0.86) than every edge-clustered configuration
  (MAC 0.13-0.27). Confirmed exact SSI-to-Forward-Model agreement at every
  common model order, not just the final one. **Scope note**: specific to
  this one frequency pair and these four configurations, not yet a general
  sensor-placement principle.

- **`level2C_sensorPlacementStudy.m`** - finite-horizon observability-Gramian
  sensor-placement study for a four-sensor limit, using the already well-separated
  `2.70/6.00 rad/s` pair so geometry is isolated from the frequency-proximity and
  amplitude-imbalance questions. Compares trace, log-determinant, minimum-eigenvalue
  and condition-number objectives, checks horizon sensitivity, then validates distinct
  finalist layouts through the full SSI pipeline.
- **`level2C_exhaustiveRecheck.m`** / **`level2C_finalCheck.m`** - exhaustive checks of
  the Level 2C candidate grid, added to separate greedy-algorithm behaviour from the
  underlying Gramian objective and to test horizon dependence at the actual 80 s SSI
  record length. These scripts are the current reference for the Gramian placement
  conclusions rather than the initial greedy run alone.
- **`level3A_frequencyProximityStudy.m`** - frequency-resolution study: holds sensor
  placement and amplitudes fixed and moves the second forced tone toward `2.70 rad/s`,
  using a coarse-then-refined sweep. The Fourier-bin scale `2*pi/T` is printed only as
  a reference scale, **not** assumed to be SSI's resolution limit.
- **`level3B_amplitudeImbalanceStudy.m`** - amplitude-imbalance study at the
  spatially and temporally easy `2.70/6.00 rad/s` pair, so disappearance of a component
  can be attributed to amplitude rather than an already difficult frequency pair. It
  tracks frequency error and MAC diagnostics before branch loss in both reciprocal
  directions.

### Standing caveat across every level above

All of these validate forced-response fidelity - does SSI faithfully
recover what a forced spectral input actually produces - not natural mode
recovery. None of them, however cleanly they pass, demonstrate that SSI can
find the floe's own natural hydroelastic modes from a free-vibration
response; that requires a separate experiment with a known free-decay
signal, not yet built.

### Still open

- A missing mode in a noisy 3-mode hand-built test case is still not fully explained
  (candidate causes include modal weighting, noise realisation, record length and
  numerical conditioning); it still needs systematic testing across noise seeds.
- The Forward-Model-to-SSI validation ladder through the current Level 3 scripts is
  still deterministic. Measured IMU noise has not yet been characterised and injected
  into this pipeline.
- The central bridge question is now explicit: the five EMM inverse observables are
  **forced-response extrema/near-zeros of dry-mode projection coefficients**, whereas
  SSI returns poles, damping and wet-system mode shapes. Those quantities must be
  connected or replaced by an SSI-compatible feature vector before the EMM inverse can
  be driven by SSI output.

---

## Inverse Mapping (`Inverse Mapping/`)

Phase 3 of the current SHARC workflow: map a compact set of Forward-Model response
features back to the EMM parameter vector

```text
p = [beta, gamma, R]^T.
```

The present inverse is intentionally local and synthetic. The reference/test parameter
set is

```text
beta0  = 4.6985e-5
gamma0 = 1.4548e-3
R0     = 0.3830
```

and is a **method-development point, not a claim about a representative real floe**.
The current five inverse observables are frequencies of three maxima and two near-zeros
of the EMM response after projection onto dry circular-plate modes. They are useful EMM
observables, but they are not automatically SSI poles or wet natural frequencies.

### Dry-mode projection and modal feature extraction

The modal basis itself lives in `Forward Model/`, because it is also useful outside the
inverse:

- **`circularPlateMode.m`** - dry free-edge circular-plate mode shapes, including rigid
  heave (`n=0,j=0`) and rigid tilt (`n=1,j=0`) plus flexural families.
- **`projectEMMOntoPlateModes.m`** - projects a complex EMM deflection field onto those
  dry modes. This gives modal coefficients `A_{n,j}`; they are coordinates of the wet
  **forced response in a dry basis**, not automatically wet natural modes.
- **`Forward Model/Validation Experiments/singleFrequencyModalProjectionTest.m`** -
  bridge test for the modal projection at one frequency.

The feature-extraction files in `Inverse Mapping/` are:

**`modalRAOReconnaissance.m`** - sweeps physical frequency and projects the EMM response
onto a selected set of dry modes. The current candidate features were found from the
responses of `A20`, `A01`, `A11`, `A21` and `A00`.

**`modalFeatureRefinement.m`** - refines the candidate feature locations on a local fine
frequency grid and saves the validated reference anchors in
`modalFeatureRefinementResults.mat`.

**`quadFit3.m`** / **`fitFeatureFrom3Points.m`** - three-point quadratic interpolation
of `|A|^2`, returning a continuous sub-grid vertex and flags for a degenerate fit,
extrapolation outside the bracket, or wrong curvature for the expected feature kind.

**`loadFeatureAnchors.m`** - loads the frozen three-frequency brackets and enforces the
fixed feature order

```text
[A20, A01, A11, A21, A00].
```

**`modalFeatureVector.m`** - the actual five-feature forward map used by the inverse.
For each trial `(beta,gamma,R)` it re-solves the EMM at the three anchors for each
feature, projects onto the relevant dry mode, fits the quadratic and returns the five
feature frequencies. `Strict=true` turns any invalid feature fit into an error so a bad
point cannot silently contaminate a Jacobian or nonlinear solve.

### Local Jacobian and identifiability

**`centralDiffJacobian.m`** - generic central finite-difference Jacobian helper.

**`computeSensitivityJacobian.m`** - evaluates the five-feature Jacobian with relative
parameter steps. Three stored runs (`eps = 1e-3, 1e-4, 1e-5`) are kept in the repository;
`1e-4` is the working step used downstream.

**`runSensitivityJacobian.m`** - driver for the production EMM Jacobian calculation.

**`analyzeIdentifiability.m`** - scales the Jacobian into approximately logarithmic
coordinates,

```text
J_s = diag(1./f0) * J * diag(p0),
```

and analyses its SVD. At the reference point the five-feature map is locally full rank,
with singular values approximately

```text
[1.195, 0.1124, 0.00408]
```

and condition number about `293`. The weakest right-singular direction is almost pure
`gamma`, so `gamma` is locally much less constrained by these feature locations than
`beta` or `R`. This is a local, scaled, noise-free result; it is not a global uniqueness
claim.

### Predicted-anchor tracking

The original feature extractor used the same three frequencies forever. That is smooth
for tiny finite-difference steps, but the first nonlinear validity probe showed it is a
poor tracker: the fixed brackets are already left by about a 1% change in `R`.

**`predictFeatureAnchors.m`** - moves each three-point bracket smoothly using the saved
reference scaled Jacobian:

```text
f_pred = f0 .* exp(J_s * log(p./p0))
shift  = f_pred - f0
```

and translates all three anchors for each feature rigidly by that shift. The prediction
only decides **where to look**; the feature returned by `modalFeatureVector.m` is still
the quadratic vertex fitted to fresh EMM evaluations at the moved anchors.

**`anchorPredictionModel.mat`** - saved full-precision `(p0,f0,J_s)` model used by the
predicted-anchor mode.

**`modalFeatureVector.m`** now supports

```matlab
'AnchorMode','fixed'      % default; original Jacobian/identifiability map
'AnchorMode','predicted'  % nonlinear inversion tracker
```

while preserving fixed mode unchanged. `predictFeatureAnchorsTester.m` checks identity
at the reference point, rigid translation, smoothness, option handling, real Forward
Model wiring and finite-difference behaviour.

A predicted-anchor validity probe kept `beta` and `gamma` valid at every tested axis
point to `+-20%`, kept `R` valid to `+-10%` (first tested failure at `+-20%`), and passed
all eight corners of the demonstrated box

```text
[beta/beta0, gamma/gamma0, R/R0]
in [0.84, 0.84, 0.92] .. [1.16, 1.16, 1.08].
```

That box is an **outer demonstrated region at the tested points/corners**, not proof that
every interior point is valid. One corner is tight (`A00` bracket margin about `0.031`),
so the nonlinear twin uses smaller bounds inside it.

### Generic parameter-inversion infrastructure (`Inverse Mapping/Parameter Inversion/`)

These functions are forward-model-agnostic: the EMM enters only through a function
handle, so the same solver machinery can later be reused for an SSI-compatible feature
map.

**`inverseResidual.m`** - converts scaled variables `x = p./pRef` back to physical
parameters, calls a supplied deterministic forward map, and returns relative feature
residuals `(f-fObs)./fObs`. Supports an optional feature mask.

**`solveInverse.m`** - bounded `lsqnonlin` wrapper using trust-region-reflective,
explicit **central differences at step `1e-4`** in scaled coordinates, iterate logging,
strict error classification and result/history packaging. Expected forward-map validity
errors are recorded as `leftValidRegion`; unexpected errors are re-thrown.

**`probeValidRegion.m`** - maps the contiguous axis neighbourhood around a reference
point, forms a conservative candidate box and checks all `2^n` corners. It reports only
validity at the tested points; it does not claim a mathematically continuous valid
region between them.

**`Parameter Inversion/Tests/`** - dedicated tests for `inverseResidual`, `solveInverse`
and `probeValidRegion`, using analytic toy problems and explicit failure-path checks.

### EMM twin inversion (`Parameter Inversion/EMM Twin/`)

**`runProbeValidRegionEMM.m`** - EMM-specific driver for fixed- or predicted-anchor
validity probes. Saved reference and predicted runs are in `Parameter Inversion/Results/`.

**`runEMMTwinInversion.m`** - first nonlinear inverse twin test. A hidden parameter truth
is passed through the EMM/predicted-anchor map once to create exact synthetic features,
then `solveInverse` attempts to recover all three parameters without receiving the
truth. The current truth and scaled bounds are

```text
x_true = [1.08, 0.93, 1.04]
lb     = [0.90, 0.85, 0.95]
ub     = [1.10, 1.15, 1.05].
```

Three distinct starting points all converged to the same hidden truth with residuals at
or near the numerical floor. The across-start spread in recovered scaled parameters was
approximately

```text
[6.1e-7, 1.65e-5, 4.5e-8]   % beta, gamma, R
```

with the largest spread, as expected, in the weak `gamma` direction. This is deliberately
a same-model/noise-free **inverse-crime integration test**; it validates the implemented
chain, not experimental robustness.

**`runFeaturePerturbationStudy.m`** - controlled feature-error sensitivity around the
solved twin. The baseline assumes independent equal relative errors on the five feature
frequencies and compares Monte-Carlo propagation with the local linear inverse. At the
solved truth the log-parameter Jacobian has approximately

```text
sigma = [1.161, 0.1108, 0.003677],  kappa ~= 316
```

and the RMS error amplification per unit relative feature standard deviation is
approximately

```text
beta: 12.35     gamma: 271.8     R: 1.12.
```

The script also perturbs the observation along the weakest left-singular direction and
runs selected full nonlinear inversions. At perturbation amplitudes `1e-4` and `3e-4`,
the nonlinear `gamma` response followed the local linear prediction to within about 4%,
confirming the local sensitivity model over the practically relevant range tested here.
This is still an idealised **feature-error study**, not a measured IMU/SSI noise model.

### Which features carry the inverse information

The full-precision solution Jacobian has also been checked under all three-, four- and
five-feature row subsets. Every subset of three or more tested features remains
numerically rank three, but practical conditioning varies strongly.

The key local finding is that the two near-zero features `A21` and `A00` jointly carry
most of the information that stabilises `gamma`. The maxima-only set
`[A20,A01,A11]` remains full rank but increases the predicted `gamma` error amplification
from about `272*sigma_f` to about `967*sigma_f`. By contrast, the best three-feature
subsets containing both near-zeros remain near `280-322*sigma_f`. `A20` is the least
consequential single deletion at this reference truth.

This matters directly for the next stage: **SSI does not automatically return these five
EMM feature frequencies.** An SSI-compatible inverse feature vector must be defined and
its own Jacobian/rank/conditioning assessed rather than assuming forced-response maxima
or near-zeros are equivalent to SSI poles.

### Inverse Mapping tests and saved results

`Inverse Mapping/Inverse Mapping Unit Tests/` currently contains tests for the modal
reconnaissance, feature vector, sensitivity Jacobian and predicted-anchor tracker. The
`Parameter Inversion/Tests/` folder covers the generic inverse solver infrastructure.

`Inverse Mapping/Parameter Inversion/Results/` stores the current reference/predicted
valid-region probes, the completed noise-free twin runs, and the feature-perturbation
`.mat`, `.txt`, `.png` and `.fig` outputs. These are reproducibility artefacts for the
current synthetic reference problem; they should not be mistaken for real-floe results.

### Current inverse-mapping status / next bridge

Established locally at the reference/test problem:

- smooth five-feature extraction and a stable local Jacobian;
- full local rank but a strongly `gamma`-dominated weak direction;
- Jacobian-predicted anchors suitable for the tested nonlinear neighbourhood;
- bounded three-parameter recovery from several starts on exact synthetic data;
- a quantified feature-error-to-parameter-error map, checked nonlinearly; and
- a feature-subset analysis showing that `A21`/`A00` information is important for
  keeping `gamma` practically constrained.

Not yet established:

- global identifiability away from this synthetic reference neighbourhood;
- realistic uncertainty from the IMU + SSI processing chain;
- equivalence between the present EMM forced-response features and anything SSI returns;
- inversion of a physically characterised real floe.

The immediate next task is therefore to define an **SSI-compatible feature vector**
`f_SSI(beta,gamma,R)` that can be generated from the Forward Model and estimated from SSI
output, then compute its Jacobian and compare its rank/conditioning with the present EMM
feature map.

> **Update (25-27 Sep):** this next task was carried out. The SSI-compatible feature vector
> turned out to be the sensor **transmissibilities**, not SSI poles; see the next two sections.

---

## SSI bridge and the pivot to transmissibilities (25 Sep)

The question after the EMM-feature inverse was: what can a **measured, output-only** record
actually deliver? These files answer it in order, each feeding the next.

### Broadband synthesis and Level 4A (`SSI/`)

**`stochasticSensorSynthesis.m`** - `[y, t, info] = stochasticSensorSynthesis(omega, H, dt, ...)`
Multi-channel stochastic time series from a sampled frequency response by frequency-domain
synthesis. One complex input coefficient per frequency is **shared by all sensors**, which
preserves the spatial coherence of a single incident field (independent per-sensor draws
would destroy the mode-shape information). `'gaussian'` (circular complex Gaussian input,
the stationary white-input assumption SSI is built on) or `'randomPhase'` (the older
forced-tone style, kept for comparison). Returns the expected second-order statistics.
Tested in `SSI/Unit Tests SSI/stochasticSensorSynthesisTester.m`.

**`Validation Experiments/level4A_broadbandKnownOscillator.m`** - does broadband synthesis
let SSI recover the **damped** poles of a known analytical modal system (no EMM solves)?
Six routes on the same system and record length, including a time-domain white-noise
control and a band-passed control. Result: the time-domain control passes from
`N >= 1024` (single mode) / `N >= 2048` (two close modes), but the band-limited routes,
including the synthesis route under test, do not pass consistently at any record length;
damping in particular is biased low by band limitation. This removed SSI damping as a
dependable observable for band-limited wave forcing.

### Sensor FRFs and the information analysis (`Inverse Mapping/`, `Output-Only Twin/`)

**`computeSensorFRF.m`** - `[H, info] = computeSensorFRF(omega, sensors, p, ...)`
EMM frequency response (displacement or acceleration) at fixed sensor locations, one EMM
solve per frequency, flagging frequencies outside the validated `alpha` range
`[1.7, 13.85]`. The time convention matches `evaluateSpectralDeflection` and
`stochasticSensorSynthesis`.

**`runSensorFRFReconnaissance.m`** - acceleration FRFs of the four Level 2C sensors
(`[0.3R,0; 0.3R,pi; 0.5R,pi; 0.9R,0]`) over 3-8.5 rad/s at `p0`. **No isolated resonance
peaks** exist in the sensor responses (only broad features, e.g. a 6 dB non-rigid bump near
5.5 rad/s); the non-rigid fraction grows from 0.04 at 3 rad/s to 0.67 at 8.5 rad/s. A local
pole fit has nothing local to fit.

**`runSensorFRFSensitivity.m`** - log-parameter derivatives of the complex sensor FRFs, an
upper bound on what any feature derived from them can carry. Weak direction almost pure
`gamma` (condition number about 150; joint std per unit relative FRF noise density:
`beta 3.21, gamma 39.3, R 0.406`). Decision: **fix `gamma`** (from measured thickness) and
invert `theta = [ln beta, ln R]`.

**`featureInformation.m`** - `out = featureInformation(Jf, Sf, Jr, Sr)`: fraction of a
reference observable's Fisher information that a candidate feature set keeps, via the
generalised eigenproblem `F_f v = lambda F_r v` (`0 <= lambda <= 1` when the noise is
propagated consistently; a warning fires otherwise).

**`runTransmissibilityInformation.m`** - transmissibilities `T_j = H_j/H_r` cancel the
unknown incident spectrum under the single-incident-field model. Complex `T` keeps
`lambda = [0.83, 0.097]` of the complex-FRF information on `[ln beta, ln R]` (std ratios
x1.72 and x2.98); magnitude alone keeps far less; the choice of reference sensor does not
matter to first order (checked to 4e-16).

**`runOutputOnlyInformationBound.m`** - what **any** output-only method (SSI poles included)
can extract, given what is known about the wave spectrum. With an arbitrary unknown
spectrum, transmissibilities carry all of it (variant E equals T-only to 2.5e-16); a known
spectrum would help (`lambda` up to `[0.99, 0.43]`), a JONSWAP-parameterised one modestly.
**Conclusion: transmissibilities are the output-only observable**, which is why the
inversion below uses them rather than SSI poles.

---

## Transmissibility estimator (`Transmissibility/`) and model side

**`extractTransmissibility.m`** - `est = extractTransmissibility(Y, dt, ...)`
Output-only Welch estimate `T_hat_j = S_jr / S_rr` (H1-type ratio, Hann, overlap, one-sided
PSDs **per rad/s**, white noise of std `sigma` giving `sigma^2 dt/pi`). Optional
reference-noise correction `S_rr - S_nn`; bins where the corrected denominator collapses
are flagged invalid and never used silently. Returns the frequency grid on which the model
must be evaluated.

**`stackTransmissibility.m` / `unstackTransmissibility.m`** - frequency-major real stacking,
per frequency `[Re T_1..T_m; Im T_1..T_m]`.

**`Inverse Mapping/transmissibilityFeatureVector.m`** - model side:
`theta -> computeSensorFRF -> T -> stacked f_T` on the estimator's grid (`ThetaIdx` default
`[1 3]`, `gamma` fixed at `p0`; reference sensor s26, index 4).

**`synthesiseTwinRecords.m`** - cached EMM twin records: EMM node solves on a frequency
grid (cached by a key including points, reference, spacing, band and grid), interpolated
to the synthesis grid, then JONSWAP-forced multi-sensor records with optional sensor noise.

Testers: `Transmissibility/Unit Tests Transmissibility/` and
`Inverse Mapping/Inverse Mapping Unit Tests/` (`computeSensorFRFTester`,
`featureInformationTester`, `transmissibilityFeatureVectorTester`).

---

## Output-only twin (`Inverse Mapping/Parameter Inversion/Output-Only Twin/`)

The synthetic validation ladder for the transmissibility inversion at `p0`, sensors as
above, reference s26, `theta = [ln beta, ln R]`, depth 1.88 m, validated band
`alpha in [1.7, 13.85]` (about 2.98-8.50 rad/s). Full write-ups are the Saturday (G2) and
Sunday (G3) reports in the project documents.

**G2 - estimator (`runTransmissibilityEstimatorTest.m`, `checkSegmentSweepCorrelation.m`)**
Monte Carlo against the exact twin: the apparent bias is fully explained by the **Welch
spectral-window (resolution) bias**, which the model can predict
(`T_pred = sum K G H_j H_r^* / sum K G |H_r|^2`); variance scales as `1/nEff`; the
Bendat-Piersol variance matches Monte Carlo (ratio 1.06); the reference-noise correction
removes the noise part of the bias. Neighbouring bins are correlated (lag 1), so
block-diagonal covariances overstate information unless bins are decimated.

**U1 - covariance (`runTransmissibilityCovariance.m`, `runL2048BiasVarianceCheck.m`)**
Finite-record covariance and bias projected into `theta`. The initial `L = 1024` setting
left a parameter bias `d_bias = 0.96` sigma (flagged); `L = 2048` cut it to about 0.26 with no
precision cost. **Frozen production setting F:** `N = 65536`, `dt = 0.1 s`, `L = 2048`,
50% overlap, every second bin in `[3.178, 8.301]` rad/s (84 frequencies, 504 features),
reference-noise correction on, block-diagonal `Sigma_k`. One 109-minute record gives
`sigma(ln beta) ~ 0.005`, `sigma(ln R) ~ 0.0009` (about 0.5% and 0.09%).

**U2 / G3 - nonlinear inversion (`runTransmissibilityTwinInversion.m`)**
Truth points A `(1.2 beta0, 0.97 R0)` and B `(0.8 beta0, 1.03 R0)`, node-cached EMM model
(0.2 rad/s spacing, interpolation error 1.7e-5), safeguarded Gauss-Newton from `p0`.
All gates pass: nonlinear solves converge in 4-5 iterations and agree with the linearised
solution to 0.003 sigma; oracle 95% coverage 0.915-0.92; a record-derived **pooled block
jackknife** covariance (31 blocks, pooled over +-2 bins, Hartlap-type correction) reaches
0.92-0.93, i.e. oracle level from one record; the mean offset equals the predicted Welch
bias (about +0.1% in `beta`). Same-model twin, so an **inverse-crime** result: it validates
the chain, not model error.

**Measurement tolerances (`runMeasurementTolerance.m`)**
Perturbs one sensor at a time (timing offset, clock-rate mismatch synced at start or at
record centre, relative scale error, both signs) and projects the exact feature change
through the U2 `J`, `Sigma`, `F`. Analytic model verified against time-warped twin records
(offset and scale exact, rate within 0.4% in `d`). Tolerance = largest error in the
contiguous safe interval from zero, reported as a percentage parameter bias (independent
of overall noise scaling, unlike `d_bias`). For **1% parameter bias** (stricter of A/B):

| Relative error between sensors | Tolerance |
|---|---|
| timing offset `dt0` | 0.37 ms |
| rate mismatch, clocks aligned at start | 0.11 ppm |
| rate mismatch, mean offset removed | 8.3 ppm |
| relative gain `s_j - s_r` | 0.42% |

`beta` sets every tolerance (at 1% `beta` bias, `R` is biased 0.08-0.29%). These are
**basin-twin requirements**, to be recomputed for the tank or MIZ Jacobian and the measured
noise. Consequence: relative timing must be known to about 0.3 ms, and clock offset and
drift must be **estimated and removed** (a shared clock or sync anchors), not just bounded.

---

## Mission Physics (`Mission Physics/`)

Screening only: no EMM runs, no identifiability verdicts.

**`waveNumber.m`** - open-water wavenumber from `omega^2 = g k tanh(kH)` (Newton from
Eckart's approximation; `H = Inf` gives deep water).

**`screenFloeRegime.m`** - `out = screenFloeRegime(floe, T, 'Depth', H, ...)`. The floe is
either physical `(R, h, E, nu, rhoIce)` or `(R, D, m)` (how the basin `p0` enters:
`D = beta rho g H^4`, `m = gamma H rho`). Returns `D, l_f, R/l_f`, and per period `k, lambda,
kH, kR, k l_f`, the ice-covered wavenumber `kappa` from the same dispersion relation as
`dispersionFunction.m`, `kappa/k`, the infinite-plate stiffness sensitivity
`d ln kappa / d ln D`, and the case in EMM scaling (`alpha, beta, gamma, R/H`; for deep
water an effective depth with `min(k,kappa) H = pi` at the longest period) with an
in-validated-range flag.

**`runFloeRegimeScreen.m`** - basin `p0`, an illustrative tank floe, and MIZ floes
(`h = 1 m, E = 6 GPa`, `R = 10/25/50 m`; `h = 0.5 m, E = 3 GPa, R = 25 m`) over 6-15 s, with a
sampling-independent log-range overlap table against `p0`. Findings: 1 m MIZ ice overlaps
`p0` substantially in `k l_f` (0.27-1.71 vs 0.15-1.15), a similar infinite-plate regime;
finite-floe similarity still depends on `kR` and `R/l_f` (`p0` has `R/l_f = 4.6`, needing
`R ~ 70 m` in 1 m ice). Stiffness sensitivity collapses toward long swell (below 0.02
beyond about 12.7 s for 1 m ice, about 9 s for 0.5 m ice). The 6 s end maps to `alpha ~ 19.6`,
beyond the validated 13.85; cutting the MIZ band near 12.5 s would keep it inside.

**`Unit Tests/testScreenFloeRegime.m`** - 23 checks, including dispersion residuals, deep
and shallow limits, scaling laws, `d ln kappa/d ln D` against finite differences, the `p0`
round trip, and **parity with `Forward Model/dispersionRoots.m`** (`kappa H` equal to the
travelling root to about 1e-15).

---

## IMU Characterisation (`IMU Characterisation/`)

The measurement chain on the real SHARC board: a Zephyr shell streaming three co-located
IMUs (LSM6DS3TR-C, LSM6DSV16X, ICM-42688-P; +-2 g, +-500 dps) over a USB virtual COM port,
logged with PuTTY ("All session output", local echo and line editing on **Auto**) via
`imu stream on <rate>`. Raw logs go to `Data/Raw/` untouched; results go to `Results/`.

### Conventions used by every function here

- **White noise**: one-sided amplitude density `N` (units/sqrt(Hz), the datasheet form):
  sample variance `N^2 fs/2`, one-sided PSD `N^2`, Allan deviation `N/sqrt(2 tau)`.
  The IEEE 952 coefficient is `N_IEEE = N/sqrt(2) = sigma_A(1 s)`.
- **PSD per rad/s** (the transmissibility and twin convention): `S_w = S_f/(2 pi)`; white
  noise of std `sigma` at `dt` gives `sigma^2 dt/pi`.
- **Flicker**: one-sided PSD `B^2/(2 pi f)`; Allan plateau `sqrt(ln2/pi) B = 0.470 B`
  (confirmed empirically, 0.467-0.475 over four records), **not** the IEEE factor 0.664.
- **Random walk**: increments `K sqrt(dt)`; Allan deviation `K sqrt(tau/3)`.
- Shell timestamps are **logging-path times** (MCU uptime when the log line was created,
  1 ms resolution) until the firmware confirms otherwise; there is no sample counter.

### Functions (`Processing/`)

**`loadImuRecord.m`** - parses a PuTTY log into per-sensor `t`, `acc` (m/s^2), `gyr`
(rad/s): strips ANSI codes and prompts, keeps only strictly matching lines, counts
malformed lines, Zephyr `messages dropped` notices (total and during streaming), and
backward time jumps (reboots). Vectorised parsing; rows kept in **arrival order** so a
reboot never interleaves two sessions.

**`checkImuTimestamps.m`** - timing audit: unwraps timestamps and counters (whole counter
wraps resolved from the timestamps), finds duplicates and gaps, fits `t_n = a + b n`
(actual rate and ppm error in device-clock units), timing residual, apparent resolution
(NaN when undetermined), windowed and quadratic drift, and host-vs-device clock rate.

**`synthesiseImuNoise.m`** - known-answer records with three clocks (IMU sampling, device
timestamp, host arrival), white + flicker + random-walk noise per channel, gravity on `z`,
and injected drops and duplicates.

**`reportImuRecord.m`** - one-command first look at one log: commanded rate (from the log
command or `_<rate>Hz` in the name), received rate `(N-1)/span` and estimated tick rate,
median / 99th-percentile / max interval, gap-like intervals (not proven losses), channel
statistics, pairwise correlation of the co-located sensors, figures and a text summary.

**`runShellRateSweep.m`** - compares logs at several commanded rates. Separate verdicts:
**loss-free** (no drop notices during streaming, malformed lines at most 0.5% of attempted
IMU lines, 99th-percentile interval at most 1.5 periods, all three sensors present, row
balance at least 0.98) and **on-rate** (received/commanded at least 0.95).

**`alignImuSensors.m`** - matches the three sensors tick by tick (they are read about 1 ms
apart in each loop tick); a tick missing in any sensor is dropped from all.

**`computeNoisePSD.m`** - one-sided Welch PSD/ASD (per Hz and per rad/s), cross-spectra
and coherence, equivalent degrees of freedom with the Welch overlap correction, 95% PSD
intervals, the expected independent-noise coherence and the single-bin 95% coherence
threshold `1 - 0.05^(1/(K_eff-1))`, band summaries for the basin (0.47-1.35 Hz) and MIZ
(0.08-0.17 Hz) bands, and a refusal to run on non-uniform timing (jitter over 5% or a gap
over 1.5 periods).

**`allanDeviation.m`** - overlapping Allan deviation on a log grid, approximate 95%
intervals, and noise terms fitted only where the local slope matches (`-1/2`, `0`, `+1/2`);
a term with no such region is returned as NaN. Averaging times with fewer than 4
equivalent degrees of freedom are shown but not fitted.

**`runStaticNoiseCharacterisation.m`** - driver for a long static record (newest
`*long*.log` in `Data/Raw/`): load, keep the longest reboot-free session, timing audit,
alignment, **longest clean contiguous block** (split where an interval exceeds 1.5x the
median or time does not increase; missing ticks are never deleted and joined), then PSD,
coherence and Allan deviation for all axes of all three IMUs, bias drift, and a summary
table with datasheet ratios (sensor x axis table: accelerometers 90 / 60 / 65-65-70
ug/sqrt(Hz); ICM-42688-P gyro 2.8 mdps/sqrt(Hz)).

### Tests (`Processing/Unit Tests/`)

| Tester | Checks | What it establishes |
|---|---|---|
| `testImuTimestamps.m` | 17 | rate error (37.0004 ppm exact), jitter, drift (also with very asymmetric gaps), gaps and duplicates at the right samples, 16-bit counter aliasing across a 70,000-sample gap, timestamp wrap, host clock, resolution including the undetermined case, white-noise scaling |
| `testNoisePSD.m` | 22 | absolute PSD level, per-rad/s level equal to `sigma^2 dt/pi`, bin scatter and CI coverage, Parseval, axis and Nyquist handling, coherence (shared, independent, 95% threshold exceeded in about 5% of bins), flicker and random-walk slopes, band rms, timing refusal, alignment |
| `testAllanDeviation.m` | 19 | exact agreement with the brute-force definition, white (`N` and `N_IEEE`), random walk, flicker plateau 0.470 B, a combined three-term record over five seeds, grid, linearity, gap refusal |
| `testStaticNoiseCharacterisation.m` | 9 | end-to-end on a synthetic Zephyr-format log with an outage and a reboot: correct session and block selection, white levels recovered on all 9 accelerometer channels via PSD (within 2.1%) and Allan (within 2.9%), coherence at chance level |

### Results so far (29 Sep, diagnostic shell path)

- **Shell rate sweep** (two passes): **10 Hz is loss-free and reproducible** (received
  9.754 Hz per sensor, 0 drop notices, 99th-percentile interval 104 ms); 20 Hz is
  borderline (drop notices, stalls to about 160 ms); 50 Hz saturates (thousands of drops,
  received 23.7 / 8.6 / 15.8 Hz, the LSM6DSV16X lines consistently most under-represented).
  This characterises the diagnostic logging path, **not** the IMU output data rate.
- The loop runs about 2.5% below the command with a timing residual of about 2-3 ms,
  consistent with a sleep-based loop rather than a hardware timer.
- **Quiet-bench noise** (74 s at 10 Hz, z-axis): LSM6DSV16X about 0.16-0.25 mg/sqrt(Hz),
  ICM-42688-P about 0.26-0.33, LSM6DS3TR-C about 0.73-0.87; roughly 2-10x the datasheet
  densities (aliasing from an unknown internal ODR/filter is the leading hypothesis).
  Coherence between sensors is at chance level on a quiet bench (an earlier log taken
  while typing on the same desk showed 0.8-0.9 correlation: bench motion). PSD and Allan
  white levels agree within about 5-8%.
- Static `|a|` of 10.18 / 10.53 / 9.93 m/s^2 (LSM6DS3TR-C / LSM6DSV16X / ICM-42688-P),
  repeatable to about 0.003 m/s^2: persistent scale errors far above the 0.42% relative-gain
  tolerance, so calibration is required. Gyro values are limited by the 3-decimal print
  step in either direction and are not yet quotable.
- The ICM-42688-P supports an external clock input, a 2 kB FIFO with 18/19-bit data and
  programmable filters: the strongest candidate for the array and a possible route to
  shared-clock synchronisation between boards.

---

## Current status and next steps (29 Sep 2026)

Established: the output-only transmissibility inversion is validated in a synthetic twin
at the basin reference, with numerical requirements on the measurement chain; MIZ regimes
are screened; the IMU processing chain is complete and validated on synthetic data, and the
diagnostic shell path is characterised.

Next:

1. long static 10 Hz record, analysed with `runStaticNoiseCharacterisation`;
2. firmware questions (internal ODR and filter, polling vs data-ready/FIFO, meaning of the
   timestamps, sample counter, ICM external clock routing, shared timing anchor between
   boards, board count and arrival date, temperature output);
3. multi-orientation gravity calibration (at least 12 orientations, local
   `g ~ 9.796 m/s^2`);
4. clock-map estimation and resampling across boards (synthetic until a second board is
   available), then comparison of measured relative timing and gain against the
   tolerances above;
5. Stage 3: EMM sensitivity and Fisher information at the selected tank/MIZ cases with the
   measured sensor noise.
