# =============================================================================
# firesim.R -- A null model of a fire regime (FORS 333, Fire Ecology)
#
# P. Higuera, University of Montana. Original MATLAB version (fireSim.m) 2012;
# R translation and revision 2026 (see MODEL_SPEC.md for the full specification
# and a list of changes from the MATLAB code).
#
# AI-use acknowledgement: The R translation, revisions, validation code, and
# documentation were developed with AI assistance from Claude (Opus 5.5,
# model claude-opus-5-5; Anthropic), October 2026, and reviewed by P. Higuera.
#
# Base R only; no packages required. Load with:  source("firesim.R")
#
# Main functions
#   simulate_fire_regime()   simulate a fire history on an n x n landscape
#   fire_history_metrics()   FRP, mFRI, mCFRI, false positives/negatives
#   run_replicates()         repeat simulation + metrics, return a data frame
#   plot_fire_history()      4-panel summary figure (as in the 2020 lecture)
#   expected_frp(), expected_mfri_obs()   analytical expectations
#
# Coordinates are 0-based (0..n-1) in x and y; cell id = x * n + y + 1 so that
# ids index R vectors directly (1..n^2).
# =============================================================================


# ---- 1. Fire-size distributions ------------------------------------------------

# Expected burned area (cells) in SIMPLE mode: squares with side ceil(Exp(mean m)),
# capped at n. The side length is geometric: P(side = k) = q (1 - q)^(k - 1),
# with q = 1 - exp(-1/m).
fs_simple_mean_cells <- function(m, n) {
  q <- 1 - exp(-1 / m)
  k <- 1:20000
  sum(pmin(k, n)^2 * q * (1 - q)^(k - 1))
}

# Mean of a truncated power law (Pareto) with pdf proportional to a^(-beta) on
# [amin, amax].
fs_pareto_mean <- function(beta, amin, amax) {
  if (abs(beta - 1) < 1e-9) return((amax - amin) / log(amax / amin))
  if (abs(beta - 2) < 1e-9) return(log(amax / amin) / (1 / amin - 1 / amax))
  num <- (amax^(2 - beta) - amin^(2 - beta)) / (2 - beta)
  den <- (amax^(1 - beta) - amin^(1 - beta)) / (1 - beta)
  num / den
}

# Simple bisection: find x in [lo, hi] with f(x) = target.
# increasing = TRUE if f increases with x.
fs_bisect <- function(f, lo, hi, target, increasing = TRUE, iter = 200) {
  for (i in seq_len(iter)) {
    mid <- (lo + hi) / 2
    v <- f(mid)
    if ((v < target) == increasing) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

# Solve for the size-distribution parameter that gives the target mean fire size
# (mean_size_frac = mean fire size as a fraction of the study area).
fs_size_params <- function(mode, mean_size_frac, n, amax_frac = 0.5) {
  A <- n^2
  target <- mean_size_frac * A
  if (mode == "simple") {
    if (mean_size_frac > 0.8) stop("mean_size_frac must be <= 0.8 in simple mode")
    m <- fs_bisect(function(m) fs_simple_mean_cells(m, n), 0.05, 500, target,
                   increasing = TRUE)
    return(list(mode = "simple", m = m))
  }
  if (mode == "realistic") {
    amin <- 1
    amax <- amax_frac * A
    if (amax_frac <= 0 || amax_frac > 0.9)
      stop("amax_frac must be in (0, 0.9] so fire areas stay exact")
    if (target <= 2 || target >= 0.95 * (amin + amax) / 2)
      stop("mean_size_frac out of range for realistic mode (need > 2 cells and ",
           "well below amax/2)")
    beta <- fs_bisect(function(b) fs_pareto_mean(b, amin, amax), 0, 4, target,
                      increasing = FALSE)
    return(list(mode = "realistic", beta = beta, amin = amin, amax = amax))
  }
  stop("mode must be 'simple' or 'realistic'")
}

# Draw one fire size: side length (simple) or area in cells (realistic).
fs_draw_size <- function(par) {
  if (par$mode == "simple") return(ceiling(rexp(1, rate = 1 / par$m)))
  b <- par$beta; a0 <- par$amin; a1 <- par$amax
  u <- runif(1)
  if (abs(b - 1) < 1e-9) return(a0 * (a1 / a0)^u)
  e <- 1 - b
  (a0^e + u * (a1^e - a0^e))^(1 / e)          # inverse-CDF sampling
}


# ---- 2. Fire footprint -----------------------------------------------------------

# Returns the ids of cells burned by one fire centered on (cx, cy).
fs_footprint <- function(par, size, cx, cy, n) {
  if (par$mode == "simple") {
    s <- max(1, min(size, n))                    # guard: 0:(s-1) misbehaves if s = 0
    off <- (0:(s - 1)) - (s - 1) %/% 2
    g <- expand.grid(dx = off, dy = off)
  } else {
    k <- max(1, round(size))
    lb <- runif(1, 1, 3)                         # length-to-breadth ratio
    amaj_max <- (n - 3) / 2                      # keep long axis < grid width
    lb <- max(1, min(lb, pi * amaj_max^2 / size))
    a <- sqrt(size * lb / pi)                    # semi-major axis
    b <- sqrt(size / (pi * lb))                  # semi-minor axis
    th <- runif(1, 0, pi)                        # orientation
    h <- min(ceiling(a) + 1, (n - 1) %/% 2)
    off <- -h:h
    g <- expand.grid(dx = off, dy = off)
    u <-  g$dx * cos(th) + g$dy * sin(th)
    v <- -g$dx * sin(th) + g$dy * cos(th)
    r2 <- (u / a)^2 + (v / b)^2                  # elliptical distance
    k <- min(k, length(r2))
    keep <- order(r2, runif(length(r2)))[1:k]    # k closest cells, random ties
    g <- g[keep, , drop = FALSE]
  }
  x <- (cx + g$dx) %% n                          # wrap-around edges (torus)
  y <- (cy + g$dy) %% n
  x * n + y + 1
}


# ---- 3. Simulation ---------------------------------------------------------------

#' Simulate a fire history under a null model of a fire regime.
#'
#' @param n_years  length of simulation (yr)
#' @param p        annual probability of a fire somewhere in the study area
#' @param mean_size_frac  mean fire size, as a fraction of the study area
#' @param mode     "simple" (squares, exponential side length) or
#'                 "realistic" (ellipses, power-law sizes)
#' @param n        landscape is n x n cells
#' @param t_change optional year when the regime changes (NULL = no change)
#' @param p_after, mean_size_frac_after  regime parameters from t_change onward
#' @param amax_frac largest possible fire (realistic mode), fraction of area
#' @param seed     optional random seed for reproducibility
#' @return a list of class "firesim"
simulate_fire_regime <- function(n_years = 2000, p = 0.5, mean_size_frac = 0.02,
                                 mode = c("realistic", "simple"), n = 100,
                                 t_change = NULL, p_after = NULL,
                                 mean_size_frac_after = NULL, amax_frac = 0.5,
                                 seed = NULL) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)
  stopifnot(p >= 0, p <= 1, mean_size_frac > 0, n >= 10, n_years >= 1)
  if (is.null(t_change) && (!is.null(p_after) || !is.null(mean_size_frac_after)))
    warning("p_after / mean_size_frac_after ignored because t_change is NULL")

  par1 <- fs_size_params(mode, mean_size_frac, n, amax_frac)
  par2 <- par1
  if (is.null(p_after)) p_after <- p
  if (is.null(mean_size_frac_after)) mean_size_frac_after <- mean_size_frac
  stopifnot(p_after >= 0, p_after <= 1, mean_size_frac_after > 0)
  if (!is.null(t_change))
    par2 <- fs_size_params(mode, mean_size_frac_after, n, amax_frac)

  cells <- vector("list", n_years)
  fx <- fy <- rep(NA_real_, n_years)            # ignition coordinates
  for (yr in seq_len(n_years)) {
    after <- !is.null(t_change) && yr >= t_change
    pp  <- if (after) p_after else p
    par <- if (after) par2 else par1
    if (runif(1) < pp) {
      size <- fs_draw_size(par)
      cx <- sample.int(n, 1) - 1
      cy <- sample.int(n, 1) - 1
      cells[[yr]] <- fs_footprint(par, size, cx, cy, n)
      fx[yr] <- cx
      fy[yr] <- cy
    }
  }
  burned <- !vapply(cells, is.null, logical(1))
  fyears <- which(burned)
  areas <- lengths(cells[burned])
  fires_df <- data.frame(
    year = fyears,
    period = if (is.null(t_change)) rep(1L, length(fyears)) else
      ifelse(fyears >= t_change, 2L, 1L),
    cx = fx[burned], cy = fy[burned],
    area_cells = areas, area_frac = areas / n^2)
  # as.integer() keeps these as integer(0), not NULL, if no fire occurred at all
  cell_vec <- as.integer(unlist(cells[burned], use.names = FALSE))
  year_vec <- as.integer(rep(fyears, times = areas))

  structure(list(
    n = n, n_years = n_years, mode = mode,
    p = p, mean_size_frac = mean_size_frac,
    t_change = t_change, p_after = p_after,
    mean_size_frac_after = mean_size_frac_after,
    size_par = list(before = par1, after = par2),
    fires = fires_df,          # one row per fire
    cell = cell_vec,           # burned cell ids (long format) ...
    year = year_vec            # ... and the year each was burned
  ), class = "firesim")
}

print.firesim <- function(x, ...) {
  cat("Null-model fire history (", x$mode, " mode)\n", sep = "")
  cat("  ", x$n, " x ", x$n, " cells; ", x$n_years, " yr; ", nrow(x$fires),
      " fires\n", sep = "")
  cat("  p = ", x$p, "; mean fire size = ", 100 * x$mean_size_frac,
      "% of study area; expected FRP = ",
      round(expected_frp(x$p, x$mean_size_frac)), " yr\n", sep = "")
  if (!is.null(x$t_change))
    cat("  Regime change in yr ", x$t_change, ": p = ", x$p_after,
        "; mean size = ", 100 * x$mean_size_frac_after,
        "%; expected FRP = ",
        round(expected_frp(x$p_after, x$mean_size_frac_after)), " yr\n", sep = "")
  invisible(x)
}


# ---- 4. Sample points ------------------------------------------------------------

#' A k_side x k_side grid of sample points ("trees"), centered on the landscape
#' center, spaced `spacing` cells apart. Default: 9 trees (3 x 3), 4 cells apart.
make_sample_points <- function(n = 100, k_side = 3, spacing = 4) {
  ctr <- n %/% 2
  off <- (0:(k_side - 1)) * spacing - ((k_side - 1) * spacing) %/% 2
  g <- expand.grid(x = ctr + off, y = ctr + off)
  g$x <- g$x %% n
  g$y <- g$y %% n
  as.matrix(g)
}

# Fire years at one cell (x, y; 0-based), within the supplied long-format vectors
fs_cell_years <- function(cell, year, x, y, n) sort(year[cell == x * n + y + 1])

#' Cell ids of the sampled area: the square spanned by the trees, extended beyond
#' the outer trees by half the tree spacing (floor(spacing / 2) cells) on every side,
#' so each tree sits at the center of the area it represents.
sampled_area_cells <- function(sample_pts, n, spacing = 4) {
  e <- spacing %/% 2
  g <- expand.grid(x = (min(sample_pts[, 1]) - e):(max(sample_pts[, 1]) + e),
                   y = (min(sample_pts[, 2]) - e):(max(sample_pts[, 2]) + e))
  sort(unique((g$x %% n) * n + (g$y %% n) + 1))
}


# ---- 5. Fire-history metrics -----------------------------------------------------

#' Fire-history statistics (FRP; mean fire return interval, mFRI; mean composite
#' fire return interval, mCFRI; false positives/negatives) for one time window of a simulated fire history.
#'
#' @param sim       output of simulate_fire_regime()
#' @param start,end first and last year of the window (record length T_len = end - start + 1)
#' @param sample_pts matrix (x, y; 0-based) of sample points for the composite record
#' @param spacing   tree spacing (cells); the sampled area extends spacing/2 beyond the outer trees
#' @param min_rec   a composite fire year requires >= min_rec points recording it
#' @param p_record  probability a burned sample point records the fire (scar forms)
#' @param p_false   annual probability a tree records a false positive (a scar attributed to fire that was not from a fire)
#' @param lake_p_record probability the lake records a fire that burns the lake (center) cell
#' @param lake_p_false  annual probability of a lake false positive (a charcoal peak attributed to fire that was not from a fire)
#' @param seed      optional seed (only matters if p_record < 1 or p_false > 0)
#' @return one-row data.frame (attribute "intervals" holds all pooled intervals)
#'
#' False positives and negatives (fire-scar composite vs. the sampled area: the square
#' spanned by the trees, extended by half the tree spacing). Truth is the fire regime: expected FRP = expected
#' mFRI at every point. Both are averages over all points in the sampled area:
#'   FP = count of composite fire years in which most (>= 50%) of the sampled area did
#'        not burn (includes false scars).
#'   FN = count of years in which most of the sampled area burned but no composite fire
#'        was registered (too few trees scarred, or, if p_record < 1, trees did not scar).
#'   Area-weighted versions (fp_area, fn_area; whole-area fires) satisfy exactly:
#'   n_composite_fires = area_equiv_fires + fp_area - fn_area.
fire_history_metrics <- function(sim, start = 1, end = sim$n_years,
                                 sample_pts = make_sample_points(sim$n),
                                 spacing = 4, min_rec = 4, p_record = 1, p_false = 0,
                                 lake_p_record = 1, lake_p_false = 0,
                                 seed = NULL) {
  if (!is.null(seed)) {
    # use a local seed without disturbing the caller's random-number stream
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  n <- sim$n; A <- n^2
  stopifnot(start >= 1, end <= sim$n_years, end >= start)
  T_len <- end - start + 1                        # record length (yr)
  sample_pts <- matrix(sample_pts, ncol = 2)      # allows c(x, y) for one point

  # Area-based: fire rotation period
  f <- sim$fires
  in_w <- f$year >= start & f$year <= end
  area <- sum(f$area_cells[in_w])
  frp <- if (area > 0) T_len / (area / A) else Inf

  # Point-based: intervals between successive fires at the same cell
  w  <- sim$year >= start & sim$year <= end
  cc <- sim$cell[w]
  yy <- sim$year[w]
  n_burns <- tabulate(cc, nbins = A)              # fires per cell
  o  <- order(cc, yy)
  same <- diff(cc[o]) == 0
  ints <- diff(yy[o])[same]
  mfri <- if (length(ints) > 0) mean(ints) else NA_real_

  # Composite record from sample points
  rec <- list()
  for (i in seq_len(nrow(sample_pts))) {
    yrs <- fs_cell_years(cc, yy, sample_pts[i, 1], sample_pts[i, 2], n)
    if (p_record < 1) yrs <- yrs[runif(length(yrs)) < p_record]  # missed scars
    fake <- if (p_false > 0) (start:end)[runif(T_len) < p_false] else integer(0)
    rec[[i]] <- union(yrs, fake)
  }
  all_rec <- unlist(rec, use.names = FALSE)
  if (length(all_rec) > 0) {
    tab <- table(all_rec)
    comp <- sort(as.integer(names(tab)[tab >= min_rec]))
  } else comp <- integer(0)
  mcfri <- if (length(comp) > 1) mean(diff(comp)) else NA_real_

  # Sampled-area truth, false positives / negatives
  area_ids <- sampled_area_cells(sample_pts, n, spacing)
  in_area <- cc %in% area_ids
  f_area <- if (any(in_area)) tapply(in_area[in_area], yy[in_area], length) / length(area_ids) else numeric(0)
  area_years <- as.integer(names(f_area))
  f_comp <- ifelse(comp %in% area_years, f_area[as.character(comp)], 0)
  fn_years <- setdiff(area_years, comp)
  lake_true <- fs_cell_years(cc, yy, n %/% 2, n %/% 2, n)   # fires at the lake (center cell)
  lake_kept <- if (lake_p_record < 1) lake_true[runif(length(lake_true)) < lake_p_record] else lake_true
  lake_fake <- if (lake_p_false > 0) setdiff((start:end)[runif(T_len) < lake_p_false], lake_true) else integer(0)
  lake_years <- sort(c(lake_kept, lake_fake))                 # lake-sediment record
  FP <- sum(f_comp < 0.5)        # composite fire registered, but most of the sampled area did not burn
  FN <- sum(f_area[as.character(fn_years)] >= 0.5)   # most of the sampled area burned, no composite fire
  fp_area <- sum(1 - f_comp)                      # area-weighted versions (whole-area fires)
  fn_area <- sum(f_area[as.character(fn_years)])
  area_equiv <- sum(f_area)
  tree_ints <- unlist(lapply(seq_len(nrow(sample_pts)), function(i)
    diff(fs_cell_years(cc, yy, sample_pts[i, 1], sample_pts[i, 2], n))))

  out <- data.frame(
    start = start, end = end, T_len = T_len,
    n_fires = sum(in_w),
    frp = frp,
    mfri = mfri,
    n_intervals = length(ints),
    frac_cells_ge2_fires = mean(n_burns >= 2),
    frac_cells_unburned = mean(n_burns == 0),
    mcfri = mcfri,
    n_composite_fires = length(comp),
    tree_mfri = if (length(tree_ints) > 0) mean(tree_ints) else NA_real_,
    area_cells = length(area_ids),
    area_equiv_fires = area_equiv,
    frp_sampled_area = if (area_equiv > 0) T_len / area_equiv else Inf,
    FP = FP, FN = FN, fp_area = fp_area, fn_area = fn_area,
    n_lake_true_fires = length(lake_true),
    n_lake_fires = length(lake_years),
    lake_FP = length(lake_fake), lake_FN = length(lake_true) - length(lake_kept),
    lake_mfri = if (length(lake_years) > 1) mean(diff(lake_years)) else NA_real_
  )
  attr(out, "intervals") <- ints
  attr(out, "composite_years") <- comp
  attr(out, "area_burned_share") <- f_area     # named by year
  attr(out, "fp_years") <- comp[f_comp < 1 - 1e-12]      # composite years that missed some points
  attr(out, "fn_years") <- sort(fn_years)
  out
}


# ---- 6. Analytical expectations --------------------------------------------------

# Expected FRP (yr) for annual fire probability p and mean fire size (fraction)
expected_frp <- function(p, mean_size_frac) 1 / (p * mean_size_frac)

# Expected pooled mFRI estimate observed in a window of T_len yr, when point
# intervals are geometric with mean frp. Long intervals are truncated by the
# window, so this is < frp unless T_len >> frp.
expected_mfri_obs <- function(frp, T_len) {
  if (T_len < 2) return(NA_real_)
  q <- 1 / frp
  L <- 1:(T_len - 1)
  fL <- q * (1 - q)^(L - 1)
  w <- (T_len - L) * fL
  sum(L * w) / sum(w)
}


# Expected FRP for a window [start, end] of a simulation, allowing for a regime
# change inside the window: T_len / (sum over years of p_y * mean_size_y).
expected_frp_window <- function(sim, start = 1, end = sim$n_years) {
  yrs <- start:end
  after <- if (is.null(sim$t_change)) rep(FALSE, length(yrs)) else yrs >= sim$t_change
  burn <- ifelse(after, sim$p_after * sim$mean_size_frac_after,
                 sim$p * sim$mean_size_frac)
  length(yrs) / sum(burn)
}


# ---- 7. Replicates ---------------------------------------------------------------

#' Run many simulations and compute metrics for one or more windows.
#'
#' @param n_reps  number of replicate fire histories
#' @param windows list of c(start, end) windows; default = whole simulation
#' @param seed    seed for the whole batch (reproducible). Do NOT put a seed in
#'                sim_args or metric_args: that would make every replicate identical.
#' @param verbose print progress every 10 replicates
#' @param sim_args list of arguments passed to simulate_fire_regime()
#' @param metric_args list of extra arguments passed to fire_history_metrics()
#' @return data.frame, one row per replicate x window
run_replicates <- function(n_reps = 50, sim_args = list(), windows = NULL,
                           metric_args = list(), seed = 1, verbose = TRUE) {
  if (!is.null(sim_args$seed) || !is.null(metric_args$seed))
    stop("Remove 'seed' from sim_args/metric_args; use run_replicates(seed = ...)")
  if (!is.null(windows) && !is.list(windows))
    stop("windows must be a list, e.g. list(c(1, 500), c(501, 1000))")
  set.seed(seed)
  res <- vector("list", n_reps)
  for (r in seq_len(n_reps)) {
    sim <- do.call(simulate_fire_regime, sim_args)
    win <- if (is.null(windows)) list(c(1, sim$n_years)) else windows
    rows <- lapply(win, function(wd)
      do.call(fire_history_metrics,
              c(list(sim = sim, start = wd[1], end = wd[2]), metric_args)))
    res[[r]] <- cbind(rep = r, do.call(rbind, rows))
    if (verbose && r %% 10 == 0) message("replicate ", r, " of ", n_reps)
  }
  do.call(rbind, res)
}


# ---- 8. Plotting -----------------------------------------------------------------

#' Four-panel summary of one fire history (layout follows the 2020 lecture):
#' map of number of fires per cell, annual area burned, distribution of fires
#' per cell, and cumulative proportion of the study area burned.
plot_fire_history <- function(sim, start = 1, end = sim$n_years,
                              sample_pts = make_sample_points(sim$n), spacing = 4) {
  n <- sim$n; A <- n^2
  m <- fire_history_metrics(sim, start, end, sample_pts, spacing = spacing)
  w <- sim$year >= start & sim$year <= end
  nb <- tabulate(sim$cell[w], nbins = A)
  # cell id = x*n + y + 1, so filling a matrix by row gives rows = x, cols = y
  M <- matrix(nb, nrow = n, ncol = n, byrow = TRUE)
  f <- sim$fires[sim$fires$year >= start & sim$fires$year <= end, ]

  op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), family = "sans")
  on.exit(par(op))

  # (1) Map
  image(0:(n - 1), 0:(n - 1), M, col = gray.colors(max(nb) + 1, 1, 0),
        breaks = seq(-0.5, max(nb) + 0.5, 1),
        xlab = "x distance", ylab = "y distance", asp = 1,
        main = sprintf("Fires per cell; mCFRI = %s yr",
                       ifelse(is.na(m$mcfri), "NA", round(m$mcfri))))
  points(sample_pts[, 1], sample_pts[, 2], pch = 15, col = "red", cex = 0.8)
  e <- spacing %/% 2                                                 # sampled area
  rect(min(sample_pts[, 1]) - e - 0.5, min(sample_pts[, 2]) - e - 0.5,
       max(sample_pts[, 1]) + e + 0.5, max(sample_pts[, 2]) + e + 0.5, lty = 2, lwd = 1.5)
  points(n %/% 2, n %/% 2, pch = 4, col = "darkgreen", cex = 1.4, lwd = 2)   # lake

  # (2) Annual area burned
  aab <- numeric(end - start + 1)
  aab[f$year - start + 1] <- f$area_frac
  plot(start:end, aab, type = "h", xlab = "Year of simulation",
       ylab = "Proportion of study area",
       main = sprintf("Annual area burned (%d fires)", nrow(f)))

  # (3) Fires per cell
  hist(nb, breaks = seq(-0.5, max(nb) + 0.5, 1), col = "gray70",
       xlab = "Number of fires", ylab = "Number of cells",
       main = sprintf("mFRI = %s yr",
                      ifelse(is.na(m$mfri), "NA", round(m$mfri))))

  # (4) Cumulative area burned
  plot(start:end, cumsum(aab), type = "s", xlab = "Year of simulation",
       ylab = "Cumulative proportion burned",
       main = sprintf("FRP = %s yr (expected %s)", round(m$frp),
                      round(expected_frp_window(sim, start, end))))
  abline(h = 0:ceiling(sum(aab)), col = "gray85", lty = 3)
  invisible(m)
}
