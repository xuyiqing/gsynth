## Synthetic Control for Multiple Treated Units
## (Causal Inference with Interactive Fixed Effects Models)
## Version
## Authors: Yiqing Xu, University of California, San Diego; Licheng Liu (Thu)
## Last Modified by Shiyun Hu
## Date: 2025.2.11

## MAIN FUNCTION
## gsynth.default()

## DEPENDENT FUNCTIONS
## fect()

## METHODS
## print.gsynth()
## plot.gsynth()

#####################################################################
## A Shell Function
#####################################################################

## default function

gsynth <- function(formula = NULL,data, # a data frame (long-form)
                           Y, # outcome
                           D, # treatment
                           X = NULL, # time-varying covariates
                           na.rm = FALSE, # remove missing values
                           index, # c(unit, time) indicators
                           weight = NULL,
                           force = "unit", # fixed effects demeaning
                           # cl = NULL,
                           r = 0, # nubmer of factors
                           lambda = NULL, ## mc method: regularization parameter
                           nlambda = 10, ## mc method: regularization parameter
                           CV = TRUE, # cross-validation
                           criterion = "mspe", # mspe or pc
                           k = 5, # cross-validation times
                           EM = FALSE, # EM algorithm
                           estimator = "ife", # ife/mc method
                           se = FALSE, # report uncertainties
                           nboots = 200, # number of bootstraps
                           inference = "parametric", # type of inference
                           # cov.ar = 1,
                           parallel = TRUE, # parallel computing
                           cores = NULL, # number of cores
                           tol = 0.001, # tolerance level
                           seed = NULL, # set seed
                           min.T0 = 5,
                           alpha = 0.05,
                           normalize = FALSE
                           ) {
    method <- "gsynth"
    if (EM == TRUE) {
        method <- "ife" # Gobillon & Magnac (2016)
    }
    if (estimator == "mc") {
        method <- "mc" # Athey et al. (2021)
    }
    if (inference == "nonparametric") {
        inference <- "bootstrap"
        warning("Using bootstrap for nonparametric inference.")
    }
    output <- fect::fect(formula = formula, data = data, method = method, Y = Y, D = D, X = X,
        na.rm = na.rm, index = index,
        force = force, r = r, lambda = lambda, nlambda = nlambda, CV = CV,
        criterion = criterion, k = k, se = se, nboots = nboots, vartype = inference,
        parallel = parallel, cores = cores, tol = tol, seed = seed, min.T0 = min.T0,
        alpha = alpha, normalize = normalize, need_cumu = TRUE)

    ##-------------------------------##
    ## storage
    ##-------------------------------##

    output$call = match.call()
    output$call$vartype <- output$call$inference # Name is compatible with `fect`
    output$data <- data # Save origninal long-form data, to utilize panelView
    class(output) <- "gsynth"
    return(output)

} ## Program GSynth ends
