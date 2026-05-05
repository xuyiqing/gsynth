# gsynth 1.5.0

* **New: placebo test**. Arguments `placeboTest = TRUE` and `placebo.period` (e.g. `c(-2, 0)`) hold out pre-treatment periods as a placebo region. Introduced in `vignette("01-gsynth", package = "gsynth")` §Placebo test, applied to IFE-EM and matrix completion in `vignette("02-ife-mc", package = "gsynth")`. (`gsynth()` does not expose `carryoverTest`: the GSC framework does not allow treatment reversals, which fect's carryover test requires. Use `fect::fect()` directly for carryover diagnostics on reversal-friendly panels.)
* **New: `ci.method` argument on `gsynth()`**. `"normal"` (default; Wald: point ± z × SE) or `"basic"` (reflected-pivot bootstrap interval, Davison and Hinkley 1997). The point estimate and SE are unaffected by the choice; only the CI columns differ. New in fect v2.4.2.
* **New: loading-overlap plot type** for the GSC estimator. `plot(fit, type = "loading.overlap")` shows treated and control loadings in the first two factor dimensions, with the convex hull of control loadings shaded; treated points outside the hull indicate the GSC counterfactual is extrapolating. See `vignette("01-gsynth", package = "gsynth")` §Loading-overlap diagnostic.
* **New: modern plot recipe** (from **fect** 2.3.x). New `plot.gsynth()` arguments `legacy.style = FALSE`, `highlight = NULL`, `highlight.fill = FALSE`. The pre-2.3 visual recipe is reproducible with `legacy.style = TRUE`.
* **New: per-role weight arguments** `W.est` (outcome-model fit) and `W.agg` (aggregation), plus `weight = ...` continues to apply to both roles. From **fect** 2.3.1.
* **New: rolling-window CV is now the default**. New arguments `cv.method = "rolling"` (default; pass `"block"` for the pre-v1.5 design), `cv.prop = 0.1`, `cv.buffer = 1`. The default `k` is flipped from 5 to 20. From **fect** 2.3.0.
* **Breaking change**: `inference = "parametric"` combined with `estimator = "ife"` or `estimator = "mc"` now errors instead of silently falling back to nonparametric bootstrap. This aligns gsynth's wrapper with **fect** v2.2.0's hard gate. Migrate to `inference = "bootstrap"` or `inference = "jackknife"` for those estimators.
* **Soft-deprecated**: `estimator = "ife"`, `estimator = "mc"`, `EM = TRUE`, and `effect()`. Each emits a one-time-per-session console message pointing to the recommended replacement. Removal is targeted for v2.0.0. **For IFE-EM, matrix completion, or post-hoc estimands (cumulative ATT, APTT, log-ATT), use [`fect::fect()`](https://yiqingxu.org/packages/fect/) and `fect::estimand()` directly.** See `vignette("02-ife-mc", package = "gsynth")` for migration recipes.
* **Dependency**: `fect (>= 2.4.2)`.

# gsynth 1.4.0

* **Breaking change**: The `estimator` parameter now maps directly to estimation methods. `estimator = "gsynth"` (new default) uses the generalized synthetic control method (Xu 2017); `estimator = "ife"` uses the IFE-EM algorithm (Gobillon & Magnac 2016); `estimator = "mc"` uses matrix completion (Athey et al. 2021). The legacy `EM = TRUE` parameter is equivalent to `estimator = "ife"` when `estimator` is not explicitly specified.
* Fixed loadings plot warnings from `GGally::ggpairs` (upstream fix in **fect**).
* Replaced pkgdown vignette with Quarto book tutorial.
* Updated DESCRIPTION: added BugReports, Encoding, GitHub URL.

# gsynth 1.3.1

* Updated DESCRIPTION and documentation.

# gsynth 1.3.0

We merged all functionalities of **gsynth** into the package **fect**. **gsynth** is now a wrapper of **fect**. Please check [fect User Manual](https://yiqingxu.org/packages/fect/) for updates. We maintain
**gsynth** solely for backward compatibility.

# gsynth 1.1.7

Using normal approximation instead of the percentile method to obtain confidence intervals based on bootstrapped standard errors. 

# gsynth 1.1.4

1. Import function *felm* from **lfe** to fit two-way fixed effects model as the 
starting value for estimation of interactive fixed effects model with unbalanced 
panel data.
2. Add cluster bootstrap option for uncertainty estimates.
3. Add jackknife uncertainty estimates.
3. Add a new function *cumuEff* for calculation of sub-group and cumulative 
treatment effects.
 
# gsynth 1.0.9

1. Function `panelView()` is removed from **gsynth** and becomes an independent package [**panelview**](https://yiqingxu.org/packages/panelview/). 
2. Implement the matrix completion method.
2. Fix bugs.
3. Change the color scheme.

# gsynth 1.0.8

1. Add a function `panelView()` to visualize raw data and data structure before estimation.
2. Fix bugs.

# gsynth 1.0.7

1. Add "implied weights" of control units for each treated unit to the output of the main function (`wgt.implied`).
2. Add a plot to visualize missing data and treatment status (`type = "missing"`).
3. Accommodate unbalanced panels.
