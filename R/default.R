## Synthetic Control for Multiple Treated Units
## (Causal Inference with Interactive Fixed Effects Models)
## Authors: Yiqing Xu, Stanford University; Licheng Liu; Ziyi Liu; Shiyun Hu

## MAIN FUNCTION
## gsynth.default()

## DEPENDENT FUNCTIONS
## fect()

## METHODS
## print.gsynth()
## plot.gsynth()

#####################################################################
## Soft-deprecation message infrastructure
#####################################################################

## Each deprecation key fires once per session.
.gsynth_deprecation_seen <- local({
    seen <- character(0)
    function(key) {
        if (key %in% seen) return(invisible(FALSE))
        seen <<- c(seen, key)
        invisible(TRUE)
    }
})

.deprecation_message <- function(key, msg) {
    if (.gsynth_deprecation_seen(key)) {
        packageStartupMessage(msg)
    }
    invisible(NULL)
}

#####################################################################
## A Shell Function
#####################################################################

## default function

gsynth <- function(formula = NULL, data, # a data frame (long-form)
                           Y, # outcome
                           D, # treatment
                           X = NULL, # time-varying covariates
                           na.rm = FALSE, # remove missing values
                           index, # c(unit, time) indicators
                           weight = NULL,
                           force = "unit", # fixed effects demeaning
                           cl = NULL,
                           r = 0, # number of factors
                           lambda = NULL, ## mc method: regularization parameter
                           nlambda = 10, ## mc method: regularization parameter
                           CV = TRUE, # cross-validation
                           criterion = "mspe", # mspe or pc
                           k = 20, # cross-validation folds (was 5 in v1.4.0)
                           cv.method = "rolling", # NEW: "rolling" or "block"
                           cv.prop = 0.1, # NEW: per-fold unit-sampling fraction
                           cv.buffer = 1, # NEW: past-side buffer for rolling CV
                           EM = FALSE, # legacy; use estimator = "ife" instead
                           estimator = "gsynth", # gsynth/ife/mc method
                           time.component.from = NULL, # NEW: NULL lets fect dispatch by estimator
                           se = FALSE, # report uncertainties
                           nboots = 200, # number of bootstraps
                           inference = NULL, # type of inference
                           parallel = TRUE, # parallel computing
                           cores = NULL, # number of cores
                           tol = 0.001, # tolerance level
                           seed = NULL, # set seed
                           min.T0 = 5,
                           alpha = 0.05,
                           normalize = FALSE,
                           W.est = NULL, # NEW: per-role weight (estimation)
                           W.agg = NULL, # NEW: per-role weight (aggregation)
                           ci.method = "normal", # NEW (v1.5.0): "normal" (Wald) or "basic" (reflected pivot) for est.* slots
                           placeboTest = FALSE, # NEW (v1.5.0): hold out pre-treatment periods as placebo
                           placebo.period = NULL # NEW (v1.5.0): event-time periods to use for placebo (e.g., c(-2, 0))
                           ) {

    ##-------------------------------##
    ## Method routing (with soft-deprecation for non-GSC paths)
    ##-------------------------------##

    ## estimator = "gsynth" -> Xu (2017): factors from control group only
    ## estimator = "ife"    -> Gobillon & Magnac (2016): EM using treated pre-treatment data
    ## estimator = "mc"     -> Athey et al. (2021): matrix completion
    ## EM = TRUE (legacy)   -> same as estimator = "ife"
    ##
    ## In gsynth v1.5.0, estimator %in% c("ife", "mc") and EM = TRUE are
    ## SOFT-DEPRECATED. Each fires a one-time-per-session message pointing
    ## the user to fect::fect(). Removal is targeted for v2.0.0.

    method <- estimator
    if (isTRUE(EM) && missing(estimator)) {
        method <- "ife" # legacy backward compatibility
        .deprecation_message(
            key = "EM_TRUE",
            msg = paste0(
                "gsynth(EM = TRUE): soft-deprecated in v1.5.0, ",
                "removal in v2.0.0. ",
                "Use fect::fect(method = \"ife\", ",
                "time.component.from = \"notyettreated\", ...) directly. ",
                "See vignette(\"02-ife-mc\", package = \"gsynth\")."
            )
        )
    }

    if (method == "ife") {
        .deprecation_message(
            key = "estimator_ife",
            msg = paste0(
                "gsynth(estimator = \"ife\"): soft-deprecated in v1.5.0, ",
                "removal in v2.0.0. ",
                "For IFE-EM, use fect::fect(method = \"ife\", ",
                "time.component.from = \"notyettreated\", ...) directly. ",
                "See vignette(\"02-ife-mc\", package = \"gsynth\") ",
                "for migration recipes."
            )
        )
    }
    if (method == "mc") {
        .deprecation_message(
            key = "estimator_mc",
            msg = paste0(
                "gsynth(estimator = \"mc\"): soft-deprecated in v1.5.0, ",
                "removal in v2.0.0. ",
                "For matrix completion, use fect::fect(method = \"mc\", ",
                "time.component.from = \"notyettreated\", ...) directly. ",
                "See vignette(\"02-ife-mc\", package = \"gsynth\") ",
                "for migration recipes."
            )
        )
    }

    ##-------------------------------##
    ## Inference routing
    ##-------------------------------##

    ## v1.5.0 behavior: do NOT silently coerce inference = "parametric"
    ## to "bootstrap" for IFE-EM / MC. fect's hard gate (v2.2.0+) errors
    ## on this combination, and the wrapper passes the error through.
    ##
    ## Public commitment: gsynth-note (Xu 2026), Appendix A.5 ("Hard gate
    ## on IFE-EM + parametric") states that this combination must error
    ## rather than auto-fallback, on the grounds that silent coercion
    ## would mask the inferential issue from analysts who deliberately
    ## chose IFE-EM.

    if (is.null(inference)) {
        if (method %in% c("ife", "mc")) {
            inference <- "bootstrap"
        } else { # gsynth
            inference <- "parametric"
        }
    }
    if (inference == "nonparametric") {
        inference <- "bootstrap"
    }

    ##-------------------------------##
    ## time.component.from: route by estimator if not explicitly set
    ##-------------------------------##

    ## "gsynth" -> nevertreated (factor estimation on controls only;
    ##   the GSC objective).
    ## "ife" / "mc" -> notyettreated (factor estimation includes treated
    ##   pre-treatment cells via EM imputation).
    if (is.null(time.component.from)) {
        time.component.from <- if (method == "gsynth") {
            "nevertreated"
        } else {
            "notyettreated"
        }
    }

    ##-------------------------------##
    ## Pass-through to fect::fect()
    ##-------------------------------##

    output <- fect::fect(
        formula = formula, data = data, method = method,
        Y = Y, D = D, X = X,
        na.rm = na.rm, index = index, cl = cl,
        force = force, r = r, lambda = lambda, nlambda = nlambda,
        CV = CV, criterion = criterion, k = k,
        cv.method = cv.method, cv.prop = cv.prop, cv.buffer = cv.buffer,
        time.component.from = time.component.from,
        se = se, nboots = nboots, vartype = inference,
        ci.method = ci.method,
        placeboTest = placeboTest, placebo.period = placebo.period,
        parallel = parallel, cores = cores, tol = tol, seed = seed,
        min.T0 = min.T0, alpha = alpha, normalize = normalize,
        W = weight, W.est = W.est, W.agg = W.agg,
        keep.sims = TRUE
    )

    ##-------------------------------##
    ## Storage
    ##-------------------------------##

    output$call <- match.call()
    output$call$vartype <- output$call$inference # Name compatible with fect
    output$data <- data # Save original long-form data, to utilize panelView
    class(output) <- "gsynth"
    return(output)

} ## Program GSynth ends
