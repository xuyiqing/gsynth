# gsynth 1.5.0

Development version, not yet on CRAN.

## Changes that affect results

* `gsynth()` with `CV = TRUE` (the default) and one number for `r` now chooses the number of factors by cross-validation from `r` to 5, as documented and as in gsynth 1.2.1. From 1.3.1 to 1.4.0 it fitted `r` factors (0 by default) without a search, so default calls now give different results. Pass `CV = FALSE` to keep a fixed `r` (#104; also #26, #34, #102).
* `weight` now weights only the averaged effects, as documented and as in gsynth 1.2.1. Use `W.est` to weight the model fit. Estimates change only when some control observations are missing (#101).
* With **fect** 2.4.7 the implied weights `wgt.implied` change: they now use the formula from Xu (2017), so the weighted control outcomes rebuild the factor part of each treated unit's counterfactual, and their rows and columns are named by unit id (#17, #82).
* `effect()` returns the cumulative effect as the running sum of the per-period effects, as documented (from **fect** 2.4.7). Before, it was k times the average effect over periods 1 to k, which differs when the number of treated units changes over time (#75).
* With **fect** 2.4.7, `estimator = "ife"` with `time.component.from = "nevertreated"`, `CV = FALSE` and `se = TRUE` gives the never-treated estimate and standard errors, the same as `estimator = "gsynth"`. Before, asking for standard errors silently switched to the not-yet-treated model.
* With **fect** 2.4.7, `effect()` builds the intervals and p-values of parametric fits (the default) and jackknife fits from normal critical values, as `fit$est.att` does. Before, it used t critical values, so its intervals are now slightly narrower: by 0.6% for a parametric fit at the default `nboots = 200` (2.5% at 50, 0.1% at 1000). The standard errors of parametric fits do not change.
* With **fect** 2.4.7, `effect()` on jackknife fits (`inference = "jackknife"`) computes its standard errors with the same jackknife formula as `fit$est.att`. Before, they were too large by the factor sqrt(N/(N - 1)) for N units (1% with 50 units). For the cumulative effect, **fect**'s `att.cumu()` gave a standard error about sqrt(N - 1) times too small; it now agrees with `effect()`.
* With **fect** 2.4.7, `effect()` on fits with `weight` (or `W.agg`) uses those weights in the estimates and standard errors, so its per-period effects equal `fit$att`. Before, it used equal weights.
* With **fect** 2.4.7, `ci.method = "basic"` with parametric inference (the default) gives valid p-values for the effects, which agree with the basic intervals. Before, they were near 1 whatever the estimate (on `simdata`, 0.97 for an ATT of 5.54 with S.E. 0.25).

## Bug fixes

* `effect()` gave confidence intervals centered near zero after a fit with the default (parametric) inference. With **fect** 2.4.7 its intervals come from the inference stored on the fit, `fit$vartype` (#37, #60, #75).
* The stored call no longer gets an extra `vartype` argument, so `print()` shows the call as typed and `update()` works when `inference` was given.
* `plot(type = "raw")` and `plot(type = "missing")` now work after a `Y =`/`D =` call or with a formula stored in a variable, use `id`, and return the plot instead of drawing it (#58, #26).
* `inference = "parametric"` with `estimator = "ife"`/`"mc"` or `EM = TRUE` now stops with a message in gsynth's own argument names (#47).
* `inference = "jackknife"` with `ci.method = "basic"` now stops with a message in gsynth's own argument names; only `ci.method = "normal"` works with the jackknife. Before, **fect** stopped the call with a message about its own argument `vartype`.
* `print()` returns the fit invisibly.
* `fect::fect_mspe()` now refits a gsynth fit with `gsynth()`, so it scores the fit's own model (with **fect** 2.4.7). Before, it refitted **fect**'s default fixed-effects model, and it stopped on a fit made with a gsynth-only argument such as `inference`.

## New features and other changes

* **New: placebo test**. Arguments `placeboTest = TRUE` and `placebo.period` (e.g. `c(-2, 0)`) hold out pre-treatment periods as a placebo region. Introduced in the [Gsynth chapter](https://yiqingxu.org/packages/gsynth/01-gsynth.html) (section "Placebo test") and applied to IFE-EM and matrix completion in the [IFE-EM and MC chapter](https://yiqingxu.org/packages/gsynth/02-ife-mc.html). (`gsynth()` does not expose `carryoverTest`: the GSC framework does not allow treatment reversals, which fect's carryover test requires. Use `fect::fect()` directly for carryover diagnostics on reversal-friendly panels.)
* **New: `ci.method` argument on `gsynth()`**. `"normal"` (default; Wald: point ± z × SE) or `"basic"` (reflected-pivot bootstrap interval, Davison and Hinkley 1997). The point estimate and SE are unaffected by the choice; the intervals and p-values differ. `"basic"` works with parametric and nonparametric inference, not with the jackknife (see `?gsynth`). New in fect v2.4.2.
* **New: loading-overlap plot type** for the GSC estimator. `plot(fit, type = "loading.overlap")` shows treated and control loadings in the first two factor dimensions, with the convex hull of control loadings shaded; treated points outside the hull indicate the GSC counterfactual is extrapolating. See the [Gsynth chapter](https://yiqingxu.org/packages/gsynth/01-gsynth.html), section "Loading-overlap diagnostic".
* **New: modern plot recipe** (from **fect** 2.3.x). New `plot.gsynth()` arguments `legacy.style = FALSE`, `highlight = NULL`, `highlight.fill = FALSE`. The pre-2.3 visual recipe is reproducible with `legacy.style = TRUE`.
* **New: per-role weight arguments** `W.est` (model fit) and `W.agg` (averaging), from **fect** 2.3.1. `weight` is the same as `W.agg`.
* **New: rolling-window CV is now the default**. New arguments `cv.method = "rolling"` (default; pass `"block"` for the pre-v1.5 design), `cv.prop = 0.1`, `cv.buffer = 1`. The default `k` is flipped from 5 to 20. From **fect** 2.3.0.
* **New: `cv.rule`, `cv.nobs` and `cv.donut` arguments**, passed to **fect**. `cv.rule` picks the number of factors from the cross-validation errors: `"1se"` (default), `"min"` (lowest error) or `"1pct"`. Before, `gsynth()` had no such argument, so `gsynth(..., cv.rule = "min")` stopped with an "unused argument" error. For `estimator = "gsynth"`, **fect** applies the rule from version 2.4.7 on; earlier versions always used `"1se"` there (xuyiqing/fect#146). `cv.nobs` and `cv.donut` set the held-out runs of block cross-validation (`cv.method = "block"`).
* **Breaking change**: `inference = "parametric"` combined with `estimator = "ife"` or `estimator = "mc"` now errors instead of silently falling back to nonparametric bootstrap. This aligns gsynth's wrapper with **fect** v2.2.0's hard gate. Migrate to `inference = "nonparametric"` (alias `"bootstrap"`) or `inference = "jackknife"` for those estimators.
* **Soft-deprecated**: `estimator = "ife"`, `estimator = "mc"`, `EM = TRUE`, and `effect()`. Each emits a one-time-per-session console message pointing to the recommended replacement. Removal is targeted for v2.0.0. **For IFE-EM, matrix completion, or post-hoc estimands (cumulative ATT, APTT, log-ATT), use [`fect::fect()`](https://yiqingxu.org/packages/fect/) and `fect::estimand()` directly.** See the [IFE-EM and MC chapter](https://yiqingxu.org/packages/gsynth/02-ife-mc.html) for migration recipes.
* **Dependency**: `fect (>= 2.4.7)`. Several **fect** 2.4.7 fixes also change gsynth results; see fect's NEWS.

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
