# fireSim 2026

A null model of a fire regime for Fire Ecology (FORS 333, Univ. of Montana). This
replaces `fireSim.m` (2012/2020), developed in Matlab (FOR 330, Univ. of Idaho).

| File | Purpose |
|---|---|
| `MODEL_SPEC.md` | Model rules, metric definitions, analytical expectations, benchmarks, and changes from the MATLAB code. The R code and the web app both implement this spec. |
| `R/demo_firesim.R` | Start here: press Source to simulate one fire history, draw two figures, and print the statistics. |
| `R/firesim.R` | Reference implementation (base R, no packages). Sourcing it only defines functions. |
| `R/validate_firesim.R` | Checks the R code against analytical expectations and the Python benchmarks. Prints PASS/FAIL. |
| `R/firesim_walkthrough.qmd` | Graduate/advanced walkthrough with 5 exercises (renders to a self-contained HTML file). |
| `dev/` | Independent Python mirror used to calibrate defaults and produce benchmarks. Not for students. |

## Run it
In R or RStudio, set the working directory to the `R/` folder and run:

```r
source("validate_firesim.R")    # about 2-5 min; expect all PASS
```

Then, in a terminal or the RStudio Terminal tab:

```
quarto render firesim_walkthrough.qmd
```

A quick start, without the walkthrough:

```r
source("firesim.R")
sim <- simulate_fire_regime(seed = 1)      # defaults: realistic, 2000 yr, p = 0.5, mean size 2% (FRP 100)
plot_fire_history(sim)
```

## Web app
`web/index.html` is a single self-contained file with no internet or server requirement. Double-click it to open it in any browser.

The source lives in `dev/web_src/` (`firesim_core.js`, `app.js`, `app_template.html`). After editing, rebuild with `python build.py`, and check the model with `node validate_core.js`.

To publish on GitHub Pages:
1. Create a public repository (e.g., `fors333-firesim`) and upload `web/index.html` to its root.
2. Go to Settings → Pages and set Source to "Deploy from a branch," branch `main`, folder `/ (root)`.
3. The app is live at `https://<username>.github.io/fors333-firesim/` within a minute or two. Post that link in Canvas.

## Acknowledgement of AI use
The R translation, model revisions, validation code, Python mirror, web app, and documentation were developed with AI assistance from Claude (Opus 5.5, model `claude-opus-5-5`; Anthropic) in October 2026, and reviewed by P. Higuera. The original model, `fireSim.m`, was written by P. Higuera (2012).

## Defaults
- p = 0.5 fires/yr, a mean fire size of 2% of the study area, and a largest possible fire of 50%, so the expected FRP is 1 / (0.5 × 0.02) = 100 yr: the boundary between fire-scar (FRP < ~100 yr) and lake-sediment (FRP > ~100 yr) territory, in mixed-severity regimes.
- The record is 2000 yr, 20 times the FRP, so the FRP and mFRI estimates sit close to theory (mFRI = FRP; the expected mFRI estimate is 95% of the FRP). Shorter analysis windows show the bias and variability that come with short records.
- Notation: FRP = fire rotation period; mFRI = mean fire return interval; mCFRI = mean composite fire return interval.
- Realistic mode (power-law sizes, ellipses) is the default and recommended for students; simple mode (squares) is kept for comparison.
- Fire-scar module (tab 2): 9 fire-scarred trees, 4 cells apart; the sampled area extends half the tree spacing beyond the outer trees (13 × 13 cells); a composite fire year requires at least 4 scarred trees. The default minimum (4) gives an mCFRI close to the true FRP. False positive: composite fire registered, but most of the sample area did not burn. False negative: most of the sample area burned, but no composite fire registered. Applies to low- to mixed-severity regimes. The center cell is a lake-sediment record.
