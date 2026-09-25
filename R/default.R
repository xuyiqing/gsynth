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
                           cv.nobs = 3, # block CV: consecutive periods per held-out run
                           cv.donut = 1, # block CV: unscored periods at each end of a run
                           cv.rule = "1se", # rule for picking r from the CV errors
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
                "See https://yiqingxu.org/packages/gsynth/02-ife-mc.html ",
                "for migration recipes."
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
                "See https://yiqingxu.org/packages/gsynth/02-ife-mc.html ",
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
                "See https://yiqingxu.org/packages/gsynth/02-ife-mc.html ",
                "for migration recipes."
            )
        )
    }

    ##-------------------------------##
    ## Inference routing
    ##-------------------------------##

    ## v1.5.0 behavior: do NOT silently coerce inference = "parametric"
    ## to "bootstrap" for IFE-EM / MC. fect's hard gate (v2.2.0+) errors
    ## on this combination; gsynth stops first, with the same condition
    ## and a message in gsynth's argument names (see below).
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
    ## Hard gate: parametric inference with IFE-EM / MC
    ##-------------------------------##

    ## Same condition as fect's gate (so the same calls stop), checked
    ## before any fitting and worded in gsynth's arguments. fect allows
    ## estimator = "ife" with time.component.from = "nevertreated", and
    ## estimator = "gsynth" (always never-treated controls).
    if (isTRUE(as.logical(se)) && identical(inference, "parametric") &&
        (method == "mc" ||
         (method == "ife" && !identical(time.component.from, "nevertreated")))) {
        stop("inference = \"parametric\" is not available for estimator = \"ife\" ",
             "or \"mc\", or for EM = TRUE. ",
             "Use inference = \"nonparametric\" or \"jackknife\".", call. = FALSE)
    }

    ##-------------------------------##
    ## Weights: `weight` weights the averaged effects only
    ##-------------------------------##

    ## As documented and as in gsynth 1.2.1, `weight` is the aggregation
    ## weight (fect's W.agg). It does not enter the model fit; W.est does.
    if (!is.null(weight)) {
        if (!is.null(W.agg) && !identical(weight, W.agg)) {
            stop("`weight` and `W.agg` name different columns. `weight` sets ",
                 "the weights for averaging treatment effects (the same role ",
                 "as `W.agg`): give one of them, or the same column in both. ",
                 "Use `W.est` for weights in the model fit.", call. = FALSE)
        }
        W.agg <- weight
    }

    ##-------------------------------##
    ## Number of factors: one r with CV searches r to 5
    ##-------------------------------##

    ## ?gsynth documents that CV = TRUE chooses the number of factors from
    ## r to 5 (as gsynth <= 1.2.1 did). fect searches r..r for one r, so
    ## the default call (r = 0) would fit no factors without a search.
    ## Widen one r below 5 to c(r, 5). A range, CV = FALSE, r >= 5 and
    ## estimator = "mc" (which has no r) pass through unchanged.
    cv.r.end <- 5 # documented upper end of the search for one r
    if (method %in% c("gsynth", "ife") && isTRUE(as.logical(CV)) &&
        is.numeric(r) && length(r) == 1L && !is.na(r) &&
        r >= 0 && r < cv.r.end) {
        r <- c(r, cv.r.end)
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
        cv.nobs = cv.nobs, cv.donut = cv.donut, cv.rule = cv.rule,
        time.component.from = time.component.from,
        se = se, nboots = nboots, vartype = inference,
        ci.method = ci.method,
        placeboTest = placeboTest, placebo.period = placebo.period,
        parallel = parallel, cores = cores, tol = tol, seed = seed,
        min.T0 = min.T0, alpha = alpha, normalize = normalize,
        W.est = W.est, W.agg = W.agg,
        keep.sims = TRUE
    )

    ##-------------------------------##
    ## Storage
    ##-------------------------------##

    ## The call is stored as typed (print() and update() use it). The
    ## inference actually used is recorded in output$vartype, which fect
    ## sets for every se = TRUE fit and its readers (effect(), plots) use.
    output$call <- match.call()
    if (isTRUE(as.logical(se))) {
        vt <- output$vartype
        if (!(is.character(vt) && length(vt) == 1L)) output$vartype <- inference
    }
    output$data <- data # Save original long-form data, to utilize panelView
    class(output) <- "gsynth"
    return(output)

} ## Program GSynth ends
