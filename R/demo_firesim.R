# =============================================================================
# demo_firesim.R -- Quick demonstration of the null model of a fire regime
#
# Open this file in RStudio and press Source. It loads firesim.R, simulates one
# 2000-year fire history with the default settings, draws two four-panel
# figures (the whole record, then only the last 100 years), and prints the
# fire-history statistics for each. Use the arrows in the Plots pane to move
# between the two figures. If R reports "figure margins too large," enlarge
# the Plots pane and run the file again.
#
# P. Higuera, University of Montana. Developed with AI assistance from Claude
# (Opus 5.5, model claude-opus-5-5; Anthropic), October 2026, and reviewed by
# P. Higuera. Base R only; no packages required.
# =============================================================================

# Load the model from the folder that holds this file (falls back to the
# working directory when lines are run one at a time).
demo_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (is.null(demo_dir) || !file.exists(file.path(demo_dir, "firesim.R"))) demo_dir <- "."
if (!file.exists(file.path(demo_dir, "firesim.R")))
  stop("Cannot find firesim.R. Set the working directory to the R/ folder ",
       "(Session > Set Working Directory > To Source File Location) and try again.")
source(file.path(demo_dir, "firesim.R"))

# ---- Settings (change these and press Source again) ----------------------------
seed       <- 1      # same seed + same settings = same fire history
n_years    <- 2000   # length of the simulated record (yr)
p          <- 0.5    # chance of one fire somewhere on the landscape each year
mean_size  <- 0.02   # mean fire size, as a fraction of the study area
short_len  <- 100    # length of the short record shown in figure 2 (yr)

# ---- Simulate one fire history ---------------------------------------------------
sim <- simulate_fire_regime(n_years = n_years, p = p, mean_size_frac = mean_size,
                            seed = seed)
frp_expected <- expected_frp(p, mean_size)

report <- function(m, label) {
  cat(sprintf("\n%s (years %d-%d, %d yr):\n", label, m$start, m$end, m$T_len))
  cat(sprintf("  Fires:                                  %d\n", m$n_fires))
  cat(sprintf("  Fire rotation period (FRP):             %.0f yr\n", m$frp))
  cat(sprintf("  Mean fire return interval (mFRI):       %.0f yr  (expected for this record length: %.0f yr)\n",
              m$mfri, expected_mfri_obs(frp_expected, m$T_len)))
  cat(sprintf("  Composite interval (mCFRI, 9 trees, >= 4 scarred): %.0f yr\n", m$mcfri))
  cat(sprintf("  False positives / false negatives:      %d / %d\n", m$FP, m$FN))
}

cat(sprintf("Null model: p = %.2f, mean fire size = %.1f%% of the study area, seed = %d\n",
            p, 100 * mean_size, seed))
cat(sprintf("Expected FRP = mFRI = 1 / (p x mean fire size) = %.0f yr\n", frp_expected))

# Figure 1: the whole record
m_all <- plot_fire_history(sim)
report(m_all, "Whole record")

# Figure 2: only the most recent short_len years
m_short <- plot_fire_history(sim, start = n_years - short_len + 1, end = n_years)
report(m_short, sprintf("Last %d years", short_len))

cat("\nQuestions: How close is each estimate to the expected value?",
    "\nWhy does the mFRI fall below the FRP in the short record?",
    "\nTry another seed, or set short_len to 500, and press Source again.\n")
