# =============================================================================
# validate_firesim.R -- checks firesim.R against analytical expectations and
# against benchmarks from the independent Python mirror (dev/firesim_mirror.py).
#
# Run from this folder:   source("validate_firesim.R")
# Runtime: about 2-5 min on a laptop. Every check prints PASS or FAIL.
# Tolerances are set at roughly 3-4 standard errors, so an occasional
# borderline FAIL can occur by chance; re-run with a different seed before
# concluding there is a bug.
#
# AI-use acknowledgement: Developed with AI assistance from Claude (Opus 5.5,
# model claude-opus-5-5; Anthropic), October 2026, and reviewed by P. Higuera.
# =============================================================================

source("firesim.R")
seed0 <- 20261009
n_pass <- 0; n_fail <- 0
check <- function(label, ok, detail = "") {
  if (isTRUE(ok)) n_pass <<- n_pass + 1 else n_fail <<- n_fail + 1
  cat(sprintf("[%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, detail))
}
rel_err <- function(obs, exp) abs(obs - exp) / exp
t_start <- Sys.time()

# ---- 1. Size calibration: realized mean fire area = target -----------------------
set.seed(seed0)
for (mode in c("simple", "realistic")) {
  par <- fs_size_params(mode, 0.02, 100)
  a <- replicate(20000, length(fs_footprint(par, fs_draw_size(par),
                                            sample.int(100, 1) - 1,
                                            sample.int(100, 1) - 1, 100)))
  se <- sd(a) / sqrt(length(a))
  tol <- 4 * se
  check(sprintf("Mean fire area, %s mode", mode),
        abs(mean(a) - 200) < tol,
        sprintf("(%.1f cells; target 200; tol +/- %.1f)", mean(a), tol))
  if (mode == "realistic")
    cat(sprintf("       beta = %.3f (Python mirror: 1.276); max area = %d cells\n",
                par$beta, max(a)))
  if (mode == "simple")
    cat(sprintf("       m = %.3f (Python mirror: 9.747)\n", par$m))
}

# ---- 2. Footprint geometry: no duplicated cells, area exact -----------------------
set.seed(seed0 + 1)
par <- fs_size_params("realistic", 0.02, 100)
ok <- TRUE
for (i in 1:500) {
  s <- fs_draw_size(par)
  ids <- fs_footprint(par, s, 0, 99, 100)          # corner: tests wrap-around
  if (anyDuplicated(ids) || length(ids) != max(1, round(s)) ||
      any(ids < 1 | ids > 10000)) ok <- FALSE
}
check("Realistic footprints: unique cells, exact area, valid ids", ok)
par <- fs_size_params("simple", 0.02, 100)
ids <- fs_footprint(par, 7, 99, 0, 100)
check("Simple footprint: 7 x 7 square wraps correctly",
      length(unique(ids)) == 49 && all(ids >= 1 & ids <= 10000))

# ---- 3. Long runs: number of fires and FRP converge on expectations ---------------
for (mode in c("simple", "realistic")) {
  rr <- run_replicates(10, sim_args = list(n_years = 10000, p = 0.25, mean_size_frac = 0.02, mode = mode),
                       seed = seed0 + 2, verbose = FALSE)
  nf <- mean(rr$n_fires); frp <- mean(rr$frp)
  check(sprintf("Fires per 10,000 yr, %s", mode),
        abs(nf - 2500) < 4 * sqrt(10000 * 0.25 * 0.75 / 10),
        sprintf("(%.0f; expected 2500)", nf))
  check(sprintf("Long-run FRP, %s", mode), rel_err(frp, 200) < 0.06,
        sprintf("(%.1f yr; expected 200)", frp))
  check(sprintf("Long-run mFRI ~ FRP, %s", mode),
        rel_err(mean(rr$mfri), expected_mfri_obs(200, 10000)) < 0.05,
        sprintf("(%.1f yr; expected %.1f)", mean(rr$mfri),
                expected_mfri_obs(200, 10000)))
}

# ---- 4. 1000-yr runs: means and spread vs. Python benchmarks ----------------------
# Python benchmarks (200 reps). The SD tolerance is wider in realistic mode
# because heavy-tailed fire sizes make the SD itself noisy (bootstrap: ~20%
# false-FAIL rate at a 25% tolerance, <2% at 45%).
bench <- list(simple    = c(frp = 202.9, frp_sd = 27.6, mfri = 154.2, mcfri = 139.1),
              realistic = c(frp = 210.9, frp_sd = 42.6, mfri = 158.3, mcfri = 151.4))
sd_tol <- c(simple = 0.25, realistic = 0.45)
for (mode in c("simple", "realistic")) {
  rr <- run_replicates(100, sim_args = list(n_years = 1000, p = 0.25, mean_size_frac = 0.02, mode = mode),
                       metric_args = list(sample_pts = make_sample_points(100, 2, 4), min_rec = 1),  # Python benchmark layout
                       seed = seed0 + 3, verbose = FALSE)
  b <- unname(bench[[mode]]); names(b) <- names(bench[[mode]])
  b <- as.list(b)
  check(sprintf("FRP mean (1000 yr), %s", mode),
        abs(mean(rr$frp) - b$frp) < 4 * b$frp_sd / sqrt(100),
        sprintf("(%.1f; Python %.1f)", mean(rr$frp), b$frp))
  check(sprintf("FRP SD (1000 yr), %s", mode),
        rel_err(sd(rr$frp), b$frp_sd) < sd_tol[[mode]],
        sprintf("(%.1f; Python %.1f)", sd(rr$frp), b$frp_sd))
  # tolerance: 4% model approximation (spatial correlation, rounding) + 4 SE
  e_mfri <- expected_mfri_obs(200, 1000)
  check(sprintf("mFRI vs. truncation formula, %s", mode),
        abs(mean(rr$mfri) - e_mfri) < 0.04 * e_mfri + 4 * sd(rr$mfri) / 10,
        sprintf("(%.1f; formula %.1f; Python %.1f)", mean(rr$mfri),
                expected_mfri_obs(200, 1000), b$mfri))
  check(sprintf("mCFRI < mFRI on average, %s", mode),
        mean(rr$mcfri, na.rm = TRUE) < mean(rr$mfri),
        sprintf("(mCFRI %.1f; Python %.1f)", mean(rr$mcfri, na.rm = TRUE), b$mcfri))
}

# ---- 5. Logic checks for composite / FP-FN ----------------------------------------
sim <- simulate_fire_regime(1000, p = 0.25, mean_size_frac = 0.02, seed = seed0 + 4)
ctr <- matrix(c(50, 50), ncol = 2)
m1 <- fire_history_metrics(sim, sample_pts = ctr, spacing = 1, min_rec = 1)   # one tree = whole sampled area
ty <- fs_cell_years(sim$cell, sim$year, 50, 50, 100)
check("One tree (tree = sampled area): FP = FN = 0", m1$FP == 0 && m1$FN == 0)
check("One tree: mCFRI equals that tree's mFRI",
      (length(ty) < 2 && is.na(m1$mcfri)) || isTRUE(all.equal(m1$mcfri, mean(diff(ty)))))
m5 <- fire_history_metrics(sim, min_rec = 2, p_record = 0.7, p_false = 0.01, seed = 2)
check("Accounting identity: composite fires = area-equivalent fires + fp_area - fn_area",
      isTRUE(all.equal(m5$n_composite_fires, m5$area_equiv_fires + m5$fp_area - m5$fn_area)))
m2 <- fire_history_metrics(sim, min_rec = 1)
m3 <- fire_history_metrics(sim, min_rec = 4)
check("Stricter filter (min_rec 4 vs 1): fewer composite fires, FP not higher, FN not lower",
      m3$n_composite_fires <= m2$n_composite_fires && m3$FP <= m2$FP &&
        m3$FN >= m2$FN)
m4 <- fire_history_metrics(sim, min_rec = 1, p_false = 0.02, seed = 1)
check("Non-fire scars (p_false > 0) add false positives", m4$FP >= m2$FP)
check("FRP identity: T / mean(fires per cell) == FRP from fire_history_metrics()",
      isTRUE(all.equal(fire_history_metrics(sim)$frp,
                       1000 / mean(tabulate(sim$cell, 1e4)))))
pts <- make_sample_points(100, k_side = 2, spacing = 3)
check("make_sample_points spacing (k_side 2, spacing 3)",
      length(unique(pts[, 1])) == 2 && diff(sort(unique(pts[, 1]))) == 3)
sim0 <- simulate_fire_regime(50, p = 0, seed = seed0 + 6)   # no fires at all
m0 <- fire_history_metrics(sim0)
check("Zero-fire run: no error, FRP = Inf, n_fires = 0",
      is.infinite(m0$frp) && m0$n_fires == 0)

# ---- 6. Regime change -------------------------------------------------------------
rr <- run_replicates(30, sim_args = list(n_years = 1000, p = 0.25, mean_size_frac = 0.02, mode = "simple", t_change = 501, p_after = 0.5),
                     windows = list(c(1, 500), c(501, 1000)),
                     seed = seed0 + 5, verbose = FALSE)
f1 <- mean(rr$frp[rr$start == 1]); f2 <- mean(rr$frp[rr$start == 501])
check("Regime change: FRP before ~200, after ~100",
      rel_err(f1, 200) < 0.12 && rel_err(f2, 100) < 0.12,
      sprintf("(%.0f -> %.0f yr; Python 207 -> 102)", f1, f2))

cat(sprintf("\n%d passed, %d failed. Runtime %.1f min.\n", n_pass, n_fail,
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
