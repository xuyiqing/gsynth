# gsynth 1.5.0

## Changes that affect results

* **One `r` with `CV = TRUE` searches `r` to 5 again**, as documented and as in gsynth 1.2.1. In gsynth 1.3.1 to 1.4.0 a single `r` was fitted without a search, so the default call (`r = 0`) fitted no factors. Default estimates change: on `turnout` with two-way fixed effects and `seed = 1`, the ATT moves from 0.83 (r = 0) to 3.40 (r = 1). Pass `CV = FALSE` to fix `r` (#26, #34, #102).
* **`weight` weights only the averaged effects** (`att.avg`, `att` and their uncertainty estimates), as documented and as in gsynth 1.2.1; it is the same as `W.agg`. In gsynth 1.3.x and 1.4.0 it was silently ignored. Use `W.est` to weight the model fit (#101).
* **Several fect 2.4.6 fixes also change gsynth results**: for example, `wgt.implied` now uses the formula of Xu (2017), `criterion = "pc"` picks the `r` with the lowest PC, and covariates that the fixed effects absorb or that are collinear are dropped with a warning (coefficient `NA`). See fect's NEWS.

## Bug fixes

* **`effect()` intervals**: after a fit with the default (parametric) inference, `effect()` gave confidence intervals centered near zero. With fect 2.4.6 it reads the inference from the fit (`fit$vartype`) (#37, #60, #75).
* **The stored call is the call as typed**: gsynth no longer adds a `vartype` argument to it, so `print()` shows the call as typed and `update()` works when `inference` was given.
* **`plot(type = "raw")` and `plot(type = "missing")`** work after a `Y =`/`D =` call, with a formula stored in a variable, or with `as.formula()`. They use `id` and return the plot instead of drawing it, so `plot(fit, type = "raw") + ggplot2::theme_bw()` draws once (#26, #58).
* **`inference = "parametric"` with `estimator = "ife"` or `"mc"`, or with `EM = TRUE`,** stops with a message in gsynth's argument names instead of fect's (#47).
* `print()` returns the fit invisibly.
* Deprecation and startup messages point to https://yiqingxu.org/packages/gsynth/02-ife-mc.html. The package ships no vignettes, so the old `vignette()` pointers found nothing.

## Documentation

* `?gsynth`: the Value section lists the elements a fit actually holds and says that `eff` and `Y.ct` are T x N matrices that include the control units (#15, #92, #102). `validX` is described correctly (#80). The arguments explain `r` and `CV` (#34), `weight`, `W.est` and `W.agg`, `parallel` and `cores` (#9), and the rules for covariates (#33).
* `?plot.gsynth`, `?print.gsynth` and `?effect` gain a Value section. `?plot.gsynth` corrects `nfactors` (loadings plot only) and `id` (ignored by the gap plot).

## New features and other changes

* **New: placebo test**. Arguments `placeboTest = TRUE` and `placebo.period` (e.g. `c(-2, 0)`) hold out pre-treatment periods as a placebo region. Introduced in the [Gsynth chapter](https://yiqingxu.org/packages/gsynth/01-gsynth.html) of the online manual (section "Placebo test") and applied to IFE-EM and matrix completion in the [IFE-EM & MC chapter](https://yiqingxu.org/packages/gsynth/02-ife-mc.html). (`gsynth()` does not expose `carryoverTest`: the GSC framework does not allow treatment reversals, which fect's carryover test requires. Use `fect::fect()` directly for carryover diagnostics on reversal-friendly panels.)
* **New: `ci.method` argument on `gsynth()`**. `"normal"` (default; Wald: point ± z × SE) or `"basic"` (reflected-pivot bootstrap interval, Davison and Hinkley 1997). The point estimate and SE are unaffected by the choice; only the CI columns differ. New in fect v2.4.2.
* **New: loading-overlap plot type** for the GSC estimator. `plot(fit, type = "loading.overlap")` shows treated and control loadings in the first two factor dimensions, with the convex hull of control loadings shaded; treated points outside the hull indicate the GSC counterfactual is extrapolating. See the [Gsynth chapter](https://yiqingxu.org/packages/gsynth/01-gsynth.html) of the online manual (section "Loading-overlap diagnostic").
* **New: modern plot recipe** (from **fect** 2.3.x). New `plot.gsynth()` arguments `legacy.style = FALSE`, `highlight = NULL`, `highlight.fill = FALSE`. The pre-2.3 visual recipe is reproducible with `legacy.style = TRUE`.
* **New: per-role weight arguments** `W.est` (model fit) and `W.agg` (averaging), from **fect** 2.3.1. `weight` is the same as `W.agg` (see above).
* **New: rolling-window CV is now the default**. New arguments `cv.method = "rolling"` (default; pass `"block"` for the pre-v1.5 design), `cv.prop = 0.1`, `cv.buffer = 1`. The default `k` is flipped from 5 to 20. From **fect** 2.3.0.
* **New: `cv.rule`, `cv.nobs` and `cv.donut` arguments**, passed to **fect**. `cv.rule` picks the number of factors from the cross-validation errors: `"1se"` (default), `"min"` (lowest error) or `"1pct"`. Before, `gsynth()` had no such argument, so `gsynth(..., cv.rule = "min")` stopped with an "unused argument" error. For `estimator = "gsynth"`, **fect** applies the rule from version 2.4.6 on; earlier versions always used `"1se"` there (xuyiqing/fect#146). `cv.nobs` and `cv.donut` set the held-out runs of block cross-validation (`cv.method = "block"`).
* **Breaking change**: `inference = "parametric"` combined with `estimator = "ife"` or `estimator = "mc"` now errors instead of silently falling back to nonparametric bootstrap. This aligns gsynth's wrapper with **fect** v2.2.0's hard gate. Migrate to `inference = "bootstrap"` or `inference = "jackknife"` for those estimators.
* **Soft-deprecated**: `estimator = "ife"`, `estimator = "mc"`, `EM = TRUE`, and `effect()`. Each emits a one-time-per-session console message pointing to the recommended replacement. Removal is targeted for v2.0.0. **For IFE-EM, matrix completion, or post-hoc estimands (cumulative ATT, APTT, log-ATT), use [`fect::fect()`](https://yiqingxu.org/packages/fect/) and `fect::estimand()` directly.** See the [IFE-EM & MC chapter](https://yiqingxu.org/packages/gsynth/02-ife-mc.html) of the online manual for migration recipes.
* **Dependency**: `fect (>= 2.4.6)` and `panelView (>= 1.3.1)`.
* **Install fect with the 2.4.6 fixes**: gsynth 1.5.0 needs the correctness fixes of fect 2.4.6. fect 2.4.6 on CRAN will include them. Until then, install fect from GitHub `dev` (`devtools::install_github("xuyiqing/fect", ref = "dev")`) after the fect fix PR is merged there. A fect `dev` build from before that merge also reports version 2.4.6 but lacks the fixes, so reinstall it.

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
