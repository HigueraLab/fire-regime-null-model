# Null Model of a Fire Regime

An interactive exercise for fire ecology courses (e.g., FORS 333 and FORS 540 at the University of Montana; FOR 326 at the University of Idaho). Students simulate a fire regime in which all variability is due to chance alone, then use the results to compare area-based and point-based fire-history metrics, to see how false positives and false negatives in fire-scar and lake-sediment records affect accuracy and interpretation, and to experience how hard it is to detect a true change in a fire regime against expected random variability.

**Open the exercise:** <https://higueralab.github.io/fire-regime-null-model/>

## What the exercise covers

The exercise has an introduction and three tabs:

- **About the model.** The three rules of the null model (when, where, and how big fires are) and its simplifying assumptions.

1.  **Simulate a fire regime.** Run a 2,000-year fire history on a uniform 100 × 100-cell landscape. Watch fires accumulate, then compare the fire rotation period (FRP) with the mean fire return interval (mFRI). Restrict the analysis window to the last 100 to 1,000 years to see how short records make the FRP erratic and bias the mFRI low. Optionally, change the fire regime partway through the run.
2.  **Fire history from fire-scarred trees and sediment-charcoal records.** Sample the simulated history with a grid of fire-scarred trees and a lake at the center of the landscape. Choose how many trees must record a fire for it to count in the composite record, and see how that choice trades false positives for false negatives and shifts the mean composite fire return interval (mCFRI). Add recorder errors: trees that burn without scarring, and charcoal peaks that do not come from fire.
3.  **Summarize statistics from many runs.** Repeat the simulation 5 to 1,000 times to see the distribution of each metric, test how record length affects the FRP, and ask how often a doubling of fire activity can be detected before vs. after a regime change. Results can be downloaded as a CSV file.

**Terms.** The fire rotation period (FRP) is the time required to burn an area equal in size to the study area. The mean fire return interval (mFRI) is the average time between fires at a point. The mean composite fire return interval (mCFRI) is the average time between fires recorded by a group of sample points, such as fire-scarred trees. A false positive is a recorded fire that did not burn most of the sampled area; a false negative is a fire that burned most of the sampled area but was not recorded.

**Seeds.** Each run uses a random seed. With the same seed and settings, everyone gets the same fire history, so class results and answer keys can be checked.

## The model

The web app runs entirely in `index.html`; no data are loaded from the web. The full rules, metric definitions, and analytical expectations are in [`MODEL_SPEC.md`](MODEL_SPEC.md), which the web app and the R version both implement.

- **When:** each year, one fire occurs somewhere on the landscape with probability *p* (default 0.5).
- **Where:** the ignition point is a random cell. The landscape is uniform, and fires that cross an edge continue on the opposite edge.
- **How big:** fire sizes follow a truncated power law (default mean 2% of the study area; largest possible fire 50%), with elliptical shapes. A simpler square-fire mode is available for comparison.
- **Expected values:** the chance that a given cell burns in a year is *p* × mean fire size, so FRP = mFRI = 1 / (*p* × mean fire size) = 100 years at the defaults.
- **Simplifications:** there are no fuel or vegetation feedbacks, no spatial variation in fuels or topography, no climate forcing, and at most one fire per year. Records are perfect unless the recorder-error settings are changed.

## Using the exercise

- **Online:** open the link above in any current browser, on a laptop or a phone.
- **Offline:** download `index.html` and open it in a browser. It needs no internet connection or installation.

## R version

The `R/` folder holds a reference implementation of the same model for graduate and advanced students. It uses base R only; no packages are required.

| File | Purpose |
|------------------------------------|------------------------------------|
| `R/demo_firesim.R` | **Start here.** Open in RStudio and press Source: it simulates one fire history, draws two figures (the whole 2,000-year record and the last 100 years), and prints the statistics for each. |
| `R/firesim.R` | The model functions. Sourcing this file only loads the functions; it produces no output on its own. Simulates fire histories and calculates FRP, mFRI, mCFRI, false positives and negatives, and the lake-sediment record. |
| `R/validate_firesim.R` | Checks the R code against analytical expectations and benchmarks from an independent Python implementation (27 tests; about 3 minutes). Prints PASS or FAIL for each test. |
| `R/firesim_walkthrough.qmd` | A Quarto walkthrough with five exercises: one fire history; FRP vs. mFRI in theory and practice; record length; fire-scar composites, false positives, and false negatives; and detecting a change in the fire regime. |

To check the code, set the working directory to `R/` and run:

``` r
source("validate_firesim.R")    # expect all 27 tests to PASS
```

To use the functions directly:

``` r
source("firesim.R")
sim <- simulate_fire_regime(seed = 1)    # defaults: realistic mode, 2000 yr, p = 0.5, mean fire size 2% (FRP = 100 yr)
m <- plot_fire_history(sim)              # four-panel figure; returns the statistics
m[, c("frp", "mfri", "mcfri", "FP", "FN")]
```

If R reports "figure margins too large," enlarge the RStudio Plots pane and run the plot again.

To render the walkthrough to a self-contained HTML file, run `quarto render firesim_walkthrough.qmd` from the `R/` folder (requires [Quarto](https://quarto.org)).

The R code and the web app use different random-number generators, so the same seed gives different fire histories in each. Their statistics agree across many runs.

## License

- **Code** (HTML, CSS, JavaScript, and R): [MIT License](LICENSE).
- **Text and figures:** [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

## How to cite

Higuera, P. E. (2026). *Null Model of a Fire Regime* [Interactive teaching tool]. FORS 333: Fire Ecology, University of Montana. <https://github.com/HigueraLab/fire-regime-null-model>

## Development note

The exercise is based on a MATLAB null model of a fire regime written by P. Higuera (2012) for FOR 326 at the University of Idaho. The 2026 web and R versions were developed by P. Higuera with AI assistance (Claude, Anthropic), and reviewed by the author.

## Contact

Philip Higuera, Department of Ecosystem and Conservation Sciences, W.A. Franke College of Forestry and Conservation, University of Montana. [philip.higuera\@umontana.edu](mailto:philip.higuera@umontana.edu){.email}
