# Forward Model - Single Floating Elastic Disc (EMM)

Phase 1 of the SHARC buoy inversion pipeline: a MATLAB implementation of Montiel's
(2012) eigenfunction matching method (EMM) for the hydroelastic response of a single
floating circular elastic plate (sea ice floe) with finite draught, forced by a
plane incident wave. Given the plate's physical parameters and an incident wave
frequency, the model returns the plate's deflection η(r,θ).

Equation numbers throughout refer to Montiel's 2012 Otago thesis, Chapter 2 and
Appendices A.1/B.1, unless stated otherwise.

## Repository layout

```
test_repo/
  Forward Model/
    (core pipeline functions, described below)
    Animation/
      precomputeDeflectionData.m
      evaluateDeflection.m
      plotDeflectionSurface.m
      plotDeflectionTimeSeries.m
      runDeflectionAnimation.m
    Unit-Tests/
      (one tester per core function, described below)
    JONSWAP/
      (spectral extension functions, described below)
      Unit-Tests-JONSWAP/
        (one tester per JONSWAP function, described below)
```

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
