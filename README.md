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
  SSI/
    (Stochastic Subspace Identification functions, described below)
    Unit-Tests-SSI/
      (one tester per SSI function, described below)
    Validation-Experiments/
      (Forward-Model-to-SSI validation ladder scripts, described below)
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

## Unit-Tests-SSI (SSI/Unit-Tests-SSI/)

Same validation philosophy as the core Unit-Tests folder: one tester per
function, checking exact identities where they exist, hand-computed cases
where they don't, and - specifically for this folder - hard structural
invariants (e.g. the `floor(n/2)` conjugate-pair bound) rather than only
plausibility checks, since that discipline is what caught the
`modalParameters.m` branch-cut bug.

## Validation-Experiments (SSI/Validation-Experiments/)

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

### Standing caveat across every level above

All of these validate forced-response fidelity - does SSI faithfully
recover what a forced spectral input actually produces - not natural mode
recovery. None of them, however cleanly they pass, demonstrate that SSI can
find the floe's own natural hydroelastic modes from a free-vibration
response; that requires a separate experiment with a known free-decay
signal, not yet built.

### Still open

- A missing mode in a noisy 3-mode hand-built test case, not yet explained
  (candidate causes: modal weighting, noise realisation, record length,
  numerical conditioning) - needs testing across multiple noise seeds.
- No noise has been added to any Forward-Model-based validation level
  (1 through 2B) - all Forward Model validation so far is deterministic and
  noise-free.
