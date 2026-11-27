# F10: thesis classes vs accuracy classes

*Generated from `reclassifyStudies` run `Results/F10/reclassify_20261127_105913` (27 Nov 2026), on the study reruns at the converged truncation `[800 10 25 20000]` (code commit 50db4e7). No EMM solve: every number below is recomputed from the raw records the studies saved.*

## 1. The two rules

Both rules use the same three class names (Identifiable / Marginal / Not identifiable). They are implemented once, in `Engine/identifiabilityClass.m`. `cfg.classRule` selects the rule and defaults to `'thesis'`.

| | Thesis rule (pre-registered) | Accuracy rule (F10) |
|---|---|---|
| Error measure | s_max = max(σ_lnβ, σ_lnR); bias as the Fisher distance d_tot = √(δθᵀ F δθ) | e_β = √(σ²_lnβ + bias²_lnβ) (ln units ≈ relative error in β) |
| Identifiable | s_max ≤ 1 %, d_tot ≤ 0.5, κ ≤ 1e8 | e_β ≤ 5 %, κ ≤ 1e8 |
| Marginal | s_max ≤ 5 %, d_tot ≤ 1, κ ≤ 1e8 | e_β ≤ 20 %, κ ≤ 1e8 |
| Bias in Fisher units | sets the class | reported as a **calibration flag**: d_tot ≤ 0.5, meaning σ is an honest error bar |
| Extra reporting | — | ladder: the smallest of 1 / 5 / 10 / 20 % that e_β meets |

**Why a second rule.** The thesis rule measures bias in units of σ. At p0, σ_lnβ = 0.05 %, so a bias that is physically negligible counts as many σ and fails the rule. Examples are 9 µm of sensor position or 0.008 mm of draught. The accuracy rule answers the question a user actually asks: how far from the true β is the estimate likely to be? The Fisher distance stays available to say whether the stated error bar can be trusted.

**Thesis regression.** Every saved point was reclassified under the thesis rule and compared with its saved class. There were 0 mismatches across A1, R1, R3, R4, V1, V1a and V2, and the D3a/D3b counts were reproduced exactly. `testClassify` also checks the classifier against the old code on 20 000 points.

## 2. Class counts per study

| Study | Points | Thesis I / M / N | Accuracy I / M / N | Changed | Calibrated (d_tot ≤ 0.5) |
|---|---:|---:|---:|---:|---:|
| A1 floe sweep | 24 | 16 / 2 / 6 | 19 / 3 / 2 | 6 | 22 |
| R1 noise sweep | 41 | 24 / 7 / 10 | 41 / 0 / 0 | 17 | 24 |
| R3 noise-floor error | 42 | 31 / 4 / 7 | 42 / 0 / 0 | 11 | 31 |
| R4 draught (γ) error | 72 | 4 / 6 / 62 | 72 / 0 / 0 | 68 | 4 |
| V1 depth and radial error (all evaluated points) | 66 | 17 / 16 / 33 | 62 / 0 / 4 | 45 | 17 |
| V1a truncation, 4 floes | 4 | 2 / 1 / 1 | 4 / 0 / 0 | 2 | 2 |
| V2 directional forcing (all evaluated points) | 147 | 23 / 39 / 85 | 133 / 7 / 7 | 117 | 23 |

The V1 and V2 counts include the bisection steps, which crowd around the thesis crossings, so read them as "points evaluated", not as a designed sample.

## 3. Tolerances (headline noise, LSM6DSV16X, p0)

| What must be known | Thesis d ≤ 0.5 | Thesis \|Δβ/β\| ≤ 1 % | Accuracy e_β ≤ 1 % | e_β ≤ 5 % | e_β ≤ 20 % |
|---|---:|---:|---:|---:|---:|
| Draught γ (R4), relative | 0.31 % (±0.008 mm) | 12.1 % | 12.3 % (±0.34 mm) | > 30 % (search range) | > 30 % |
| Sensor radial position (V1c) | ±8.9 µm | −0.443 / +0.447 mm | −0.434 / +0.453 mm | ≈ ±2.2 mm\* | ≈ −8.6 / +9.2 mm\* |
| Water depth (V1b) | not crossed in ±5 % | not crossed | not crossed | not crossed | not crossed |
| Noise floor (R3) | not crossed (LSM6DS3TR-C: 65 %) | not crossed | not crossed | not crossed | not crossed |

\* Log-interpolated between the evaluated points at 1 mm and 10 mm. Treat these as estimates; a short V1c grid refinement would pin them down.

R4 for all six noise cases: the e_β ≤ 1 % tolerance on γ is 11.6–12.3 %, against 0.12–0.85 % under the thesis rule.

## 4. Directional forcing (V2): largest second-direction energy fraction f₂

| φ | Thesis d_sys ≤ 0.5 | Thesis \|Δβ\| ≤ 1 % | e_β ≤ 1 % | e_β ≤ 5 % | e_β ≤ 20 % | Coherence detects at |
|---:|---:|---:|---:|---:|---:|---:|
| 15° | 8.2e-4 | 0.033 | 0.033 | 0.165 | > 0.5 | 1.2e-3 |
| 30° | 2.2e-4 | 0.0087 | 0.0085 | 0.043 | 0.168 | 8.2e-5 |
| 60° | 9.6e-5 | 0.0047 | 0.0046 | 0.023 | 0.090 | 9.2e-6 |
| 90° | 5.7e-5 | 0.0059 | 0.0059 | 0.029 | 0.109 | < 5e-6 |

Under the thesis rule, the coherence test gives no early warning at φ = 15°, because the bias fails first. Under the accuracy rule, coherence detects the second direction **before** β degrades at every angle: at 15° it detects at 0.0012 and β reaches 1 % at 0.033.

## 5. Sensor layout (D3, LSM6DSV16X)

| | n_s = 2 | n_s = 3 | n_s = 4 |
|---|---:|---:|---:|
| D3a (heading unknown): layouts Identifiable at every heading, thesis | 28 / 528 | 1 760 / 5 456 | 14 844 / 40 920 |
| D3a, accuracy | 464 / 528 | 5 456 / 5 456 | 40 920 / 40 920 |
| D3a, best worst-heading e_β | 0.35 % | 0.10 % | 0.07 % |
| D3b (directional seas): feasible layouts, thesis | 0 | 0 | 0 |
| D3b, accuracy Identifiable / Marginal | 0 / 0 | 0 / 0 | 0 / 20 |
| D3b, best worst-case e_β | 347 % | 32 % | 18.4 % |

The ICM-42688-P and G3 noise cases behave the same way: 464 / 528 layouts at n_s = 2, and every layout at n_s = 3 and 4. D3b has no κ stored per layout; the D3a worst κ is 7.9e3.

## 6. A1 floes whose class changes

| Floe | Thesis | Accuracy | e_β | Calibrated | Note |
|---|---|---|---:|:---:|---|
| ring_R_x2 | Marginal | Identifiable | 0.044 % | no | d_W = 0.67: accurate, but σ understates the error |
| corner_β×0.5_R×2 | Not identifiable | Identifiable | 0.049 % | no | d_W = 1.07 |
| pp5mm_D40 | Marginal | Identifiable | 1.7 % | yes | stiff and small (R/l_f = 0.96) |
| pp10mm_D40 | Not identifiable | Marginal | 12.8 % | yes | R/l_f = 0.57 |
| modelFloe_E3.0GPa | Not identifiable | (Marginal) | (7.0 %) | yes | **σ not converged with truncation: keep Not identifiable** |
| miz_R10_h1 | Not identifiable | (Marginal) | (12.0 %) | yes | **σ not converged: keep Not identifiable** |

modelFloe E0.5 and E1.5 GPa stay Not identifiable under both rules (e_β 35 % and 26 %).

## 7. Findings for the paper

1. **Positioning, draught and model truncation are not the limiting errors.** Under the thesis rule they looked fatal because σ is so small. Their real effect on β is below 1 % for ±0.44 mm of sensor position, 12 % of draught, and any converged truncation.
2. **Fewer sensors suffice when only accuracy is required.** With an unknown heading, 464 of 528 two-sensor layouts keep β within 5 % at every heading, and the best keeps it within 0.35 %. The thesis rule needed n_s = 4.
3. **Directional seas remain the real limitation.** No layout up to n_s = 4 keeps β within 5 % in every directional sea tested; the best four-sensor layout reaches 18 %. The gain from the new rule is that coherence now flags the problem before it matters, at every angle.
4. **Stiff, small floes are a property of the physics, not of the rule.** They stay poorly identifiable (R/l_f ≲ 0.6). pp5mm_D40 at R/l_f = 0.96 is now Identifiable to 1.7 %.
5. **The calibration flag carries the Fisher-distance information.** "Accurate but not calibrated" (for example ring_R_x2) means σ understates the error, and the reported uncertainty should be inflated to e_β.

## 8. Open items

- Thresholds 5 % / 20 %, and whether the accuracy class should also gate on e_R: to agree with R. Verrinder.
- V1c 5 % and 20 % tolerances: interpolated over 1–10 mm. A finer grid would take about 10 minutes of V1-style evaluations.
- The A1 model-check floes need σ converged (M = 1600 or Richardson) before any accuracy class is quoted for them.
- D1, D2 and R2 are not reclassified. D1 ranks by A-optimality, not by class. D2 and R2 did not save bias vectors and would need a rerun with `rawOf`, which is fast.
