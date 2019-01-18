#######################################################
## METHODS
#######################################################

##########
## Print
##########
## a gsynth object
print.gsynth <- function(x,  
                         ...) {
    
    cat("Call:\n")
    print(x$call, digits = 4)
    
    if (is.null(x$est.avg) == TRUE) { # no uncertainties
        cat("\nAverage Treatment Effect on the Treated:\n")
        print(x$att.avg, digits = 4)
        cat("\n   ~ by Period (including Pre-treatment Periods):\n")
        print(x$att, digits = 4)
        if (is.null(x$X) == FALSE) {
            cat("\nCoefficients for the Covariates:\n")
            print(x$beta, digits = 4)
        }
        cat("\nUncertainty estimates not available.\n")
    } else {
        cat("\nAverage Treatment Effect on the Treated:\n")
        print(x$est.avg, digits = 4)
        cat("\n   ~ by Period (including Pre-treatment Periods):\n")
        print(x$est.att, digits = 4)
        if (is.null(x$X) == FALSE) {
            cat("\nCoefficients for the Covariates:\n")
            print(x$est.beta, digits = 4)
        }
    }
}