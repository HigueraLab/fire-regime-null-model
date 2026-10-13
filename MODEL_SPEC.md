# fireSim 2026: model specification

This specification defines the null model of a fire regime used in FORS 333. The R version (`R/firesim.R`) and the web app both implement it. The Python mirror in `dev/` was used to calibrate defaults and produce benchmark values. Version 1.0 (2026-10-09). Developed with AI assistance from Claude (Opus 5.5, model `claude-opus-5-5`; Anthropic) and reviewed by P. Higuera.

## Notation
FRP = fire rotation period; mFRI = mean fire return interval (point scale); mCFRI = mean composite fire return interval; FP/FN = false positive/negative. Defaults (app and R): realistic mode, p = 0.5, mean fire size 2%, largest possible fire 50%, 2000 yr, giving FRP = 100 yr (the boundary between fire-scar and lake-sediment territory) and record ÷ FRP = 20 (expected mFRI estimate 95% of FRP). The benchmarks and validation in sections 7–8 use p = 0.25, mean size 2% (FRP = 200 yr).

## 1. Purpose
The model generates fire histories in which all variability in fire timing, location, and size is random: there are no trends, no spatial gradients, and no feedbacks. Students use it to (a) see how much fire-history statistics vary by chance alone, (b) compare area-based (FRP) and point-based (mFRI, mCFRI) metrics, (c) quantify false positives and false negatives in composite records, and (d) ask whether a real change in the fire regime can be detected against that background variability.

## 2. Landscape
- An *n* × *n* grid of cells (default *n* = 100; study area *A* = 10,000 cells; units are arbitrary).
- The landscape is homogeneous: every cell can burn, and nothing about a cell's history affects its probability of burning (no fuel feedback).
- Edges wrap around (a torus). A fire that crosses one edge continues on the opposite edge. This mimics a study area embedded in a larger, similar landscape and removes edge effects. On the plotted map, a single fire can therefore appear as pieces at opposite edges.

## 3. Fire occurrence (WHEN)
- In each year, one fire occurs with probability *p* (Bernoulli; default *p* = 0.5). There is at most one fire per year.
- The ignition point is a uniformly random cell (WHERE).

## 4. Fire size and shape (HOW BIG)
Both modes are parameterized by the **target mean fire size**, *s̄*, given as a fraction of the study area (default 0.02 = 2%). That makes the expected FRP identical across modes:

> **Expected FRP = 1 / (*p* × *s̄*)**; defaults give 1 / (0.5 × 0.02) = **100 yr**; the validation benchmarks use 1 / (0.25 × 0.02) = 200 yr.

**Simple mode (earlier classroom version; kept for comparison).** Each fire is a square centered on the ignition cell. Its side length is ceil(Exp(mean *m*)) cells, capped at *n*. The value of *m* is solved numerically so that E[min(side, *n*)²] = *s̄A*. (For reference, the 2020 classroom settings burned, on average, *s̄* ≈ 1.5% of the study area per fire for FS = 8 and ≈ 4.9% for FS = 15. These are the realized footprints, (L+1)²; the old code's own `L^2` accounting gave 1.36% and 4.6%.)

**Realistic mode (default; recommended for students).** Fire sizes follow a truncated power law (Pareto), with pdf ∝ *a*^(−β) on [1, *a*_max] cells. The default *a*_max is 0.5*A*. The exponent β is solved numerically to match *s̄*. Defaults give β ≈ 1.28 for *s̄* = 2%; β ≈ 0.99 for 6%; β ≈ 1.36 for 1.36%.
- Shape: an ellipse with a random orientation and a length-to-breadth ratio drawn from U(1, 3). The ratio is reduced for the largest fires so the long axis never exceeds the grid width, which keeps a fire from wrapping onto itself.
- Area is exact: the fire burns the round(*a*) cells closest to the ignition point in elliptical distance (ties broken at random). The realized area therefore equals the drawn area (the rounding bias is about 0.01%).
- Same mean, heavier tail: a few large fires account for most of the area burned, so FRP varies far more between runs than in simple mode.

## 5. Regime change (optional)
Starting in year *t*_change, the model switches to *p*_after and/or *s̄*_after. Metrics can be computed for any time window (for example, before vs. after the change).

## 6. Fire-history metrics (window from *start* to *end*, length *T* = end − start + 1)
| Metric | Definition | Scale |
|---|---|---|
| **FRP** | *T* / (Σ fire areas in the window / *A*). Areas are realized burned cells, so reburns count again. | Area |
| **mFRI** (mean fire return interval, pooled) | The mean of all intervals between successive fires at the same cell, with both fires inside the window, pooled across all cells. This is the point-based estimate a fire-scar study reports. In theory, mFRI = FRP. | Point |
| **mCFRI** (mean composite fire return interval) | The mean interval between composite fire years recorded by *k* sample points. The default is 4 points in a 2 × 2 grid spaced 4 cells apart around the landscape center. A year counts as a composite fire year if ≥ *min_rec* points record it. | Composite |
| **Sampled area** | The square spanned by the fire-scarred trees, extended beyond the outer trees by half the tree spacing (⌊spacing/2⌋ cells) on every side, so each tree sits at the center of the area it represents: ((*k* − 1)·spacing + 2⌊spacing/2⌋ + 1)² cells. Default: 9 trees (3 × 3), 4 cells apart, so 13 × 13 cells. | Composite |
| **Composite fire year** | A year in which at least *min_rec* trees were scarred (default 4 of 9). Applies to fire-scar records, i.e., low- to mixed-severity regimes. | Composite |
| **Truth** | The fire regime: expected FRP = expected mFRI at every point. For one run, the realized value is the FRP of the sampled area = *T* / (mean fires per point in the sampled area). | Area |
| **False positive (FP)** | Count of composite fire years in which most (≥ 50%) of the sampled area did not burn (includes false scars). Shortens the mCFRI. | Composite vs. area |
| **False negative (FN)** | Count of years in which most of the sampled area burned but no composite fire was registered (too few trees scarred, or burned trees did not scar). Lengthens the mCFRI. | Composite vs. area |
| Area-weighted FP/FN (advanced) | fp_area = Σ over composite years of the unburned share; fn_area = Σ over missed fire years of the burned share. Exactly: composite fire years = area-equivalent fires + fp_area − fn_area. | Composite vs. area |
| **Lake-sediment record** | The center cell records fires that burn it, with probability *lake_p_record* (default 1); false positives (charcoal peaks attributed to fire that was not from a fire) occur with annual probability *lake_p_false* (default 0). Lake mFRI = mean interval between recorded events. Cells have no units, so lake size is relative to cell and fire size: when fires typically span many cells, a one-cell lake is realistic. | Point |

Two optional recorder-error parameters make the FP/FN module more realistic:
- *p*_record: the probability that a burned tree records the fire (a scar forms). Default 1. The lake has its own *lake_p_record*.
- *p*_false: the annual probability of a tree false positive (a scar attributed to fire that was not from a fire). Default 0. The lake has its own *lake_p_false* (a charcoal peak attributed to fire).

Note: *T* / mean(*n_i*) across all cells is algebraically identical to FRP. This is why point-based and area-based estimates converge in a homogeneous landscape.

## 7. Analytical expectations used for validation
1. The mean realized fire area equals *s̄A* in both modes (exact in expectation; in realistic mode E[round(*a*)] = 199.99 cells for a target of 200). Sample means vary because fire areas are highly variable; in realistic mode the SD is about 617 cells.
2. The expected number of fires is *pT*.
3. Each cell burns in a given year with probability *q* = *p* *s̄*, independently across years, so point fire intervals are geometric with mean 1/*q* = FRP.
4. **The expected pooled mFRI estimate from a finite window** is less than the true mFRI (= FRP), because intervals longer than the window cannot be observed (truncation):
   E[mFRI_obs] ≈ Σ_L L (T − L) f(L) / Σ_L (T − L) f(L), where f(L) = q(1 − q)^(L−1).
   For FRP = 200, this gives 196 yr (*T* = 10,000; record ÷ FRP = 50), 178 yr (2000), 152 yr (1000), 110 yr (500), 68 yr (250), and 31 yr (100). As a ratio of the FRP it depends only on record ÷ FRP: 0.48 (2×), 0.76 (5×), 0.89 (10×), 0.95 (20×), 0.98 (50×).
   *T* / mean(fires per cell), averaged over all cells including unburned ones, equals the FRP exactly; it uses the censored time before the first and after the last fire.
5. The FRP estimator is a ratio, so for short windows its mean exceeds the expected FRP (Jensen's inequality), and its median does too in realistic mode.

## 8. Benchmarks from the Python mirror (p = 0.25, mean size 2%, 1000-yr runs)
| | Simple (200 reps) | Realistic (200 reps) | Expected |
|---|---|---|---|
| Fires per run | 248 ± 15 | 250 ± 14 | 250 |
| FRP (mean ± SD) | 203 ± 28 | 211 ± 43 | 200 |
| mFRI | 154 ± 17 | 158 ± 24 | 152 (theory 200) |
| mCFRI (4 pts, 4 apart) | 139 ± 72 | 151 ± 97 | none (benchmark only) |

How the 95% range of FRP changes with record length (100 reps):

| Record length | Simple | Realistic |
|---|---|---|
| 100 yr | 84–617 yr | 90–1315 yr |
| 250 yr | 123–346 yr | 113–527 yr |
| 500 yr | 130–297 yr | 124–387 yr |
| 1000 yr | 147–256 yr | 140–295 yr |

Regime change (*p* doubled from 0.25 to 0.5 at year 501): FRP was 207 ± 44 yr before the change and 102 ± 15 yr after.

### Fire-scar composite benchmarks (JS app: p = 0.3, mean size 2.22% [FRP 150], 4000 yr, 9 trees 4 apart, 100 runs)
At the current defaults (realistic, p = 0.5, mean size 2%, 2000 yr, FRP 100, 9 trees 4 apart, 13 × 13-cell sampled area, minimum 4; 150 runs): mCFRI 101 yr, FP 3.8, FN 3.9 whole-area fires; minimum 2 / 3 / 4 / 5 / 6 gives mCFRI 69 / 83 / 101 / 111 / 128 yr. Balanced by design. (The table below used the earlier 9 × 9 area, before the half-spacing extension; mCFRI values are unaffected by the area definition, FP/FN are.)

| | Simple | Realistic | Truth |
|---|---|---|---|
| mCFRI, minimum 1 / 2 / 3 / 4 / 5 / 7 / 9 trees | 75 / 90 / 113 / **155** / 186 / 333 / 333 yr | 81 / 103 / 121 / **146** / 159 / 215 / 287 yr | 150 yr |
| FP / FN per point at minimum 4 (default) | 4.6 / 5.5 | 3.7 / 3.7 | none |
| Composite fire years at minimum 4 | 25.7 | 27.1 | 26.7 fires per point |
| FRP of the sampled area | 155 yr | 152 yr | 150 yr |
| Lake mFRI | 152 yr | 148 yr | 150 yr |

(Squares aligned with the grid burn trees in rectangular blocks, so in simple mode some minimums give identical results, e.g., 5 = 6 and 7 = 8 = 9.) FP/FN were first defined against a single unsampled center point, then by whole-area burns. Both versions were replaced on 2026-10-12 with the per-point definitions above. The Python mirror in `dev/` still uses the original center-point definition.

## 9. Changes from fireSim.m (2012/2020)
- mCFRI is now computed from fire **years**. The old code differenced fire *index numbers*, which made it about 3.3× too short (the mean number of years between landscape fires).
- Fire area is the number of burned cells. The old code used `L^2` for a footprint of `(L+1)^2` cells, which overstated FRP by 7–14%.
- mFRI is now interval-based. The old code averaged *T* / *n_i* over burned cells, which is biased high (mean of reciprocals; unburned cells dropped): 278 yr vs. a true 216 yr in a test of the 2020 settings. That metric has been removed.
- Occurrence is Bernoulli(*p*), as the slides describe. The old code used rounded Weibull intervals.
- Fires are centered on the ignition cell. In the old code the "center" was actually a corner; with wrap-around edges this has no statistical effect.

## 10. Known simplifications (state these to students)
- There is no fuel or vegetation feedback (real fire probability depends on time since the last fire), no spatial heterogeneity, no climate forcing, and at most one fire per year.
- Square and elliptical shapes stand in for irregular fire perimeters. Real fires leave unburned islands; these do not.
- The default power-law exponents (about 1.0–1.4) are set by matching the mean fire size, not fitted to data. Published frequency–area exponents for western U.S. regions are roughly in this range, but check specific values before quoting them to students.
