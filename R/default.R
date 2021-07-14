## Synthetic Control for Multiple Treated Units
## (Causal Inference with Interactive Fixed Effects Models)
## Version 1.0.9
## Authors: Yiqing Xu, University of California, San Diego; Licheng Liu (Thu)
## Date: 2019.1.18

## MAIN FUNCTION
## gsynth.formula()
## gsynth.default()

## DEPENDENT FUNCTIONS
## synth.core()
## synth.em()
## synth.em.cv()
## synth.mc()
## synth.boot()

## METHODS
## print.gsynth()
## plot.gsynth()


#####################################################################
## A Shell Function
#####################################################################

## generic function

gsynth <- function(formula = NULL, data, # a data frame (long-form)
                   Y, # outcome
                   D, # treatment
                   X = NULL, # time-varying covariates
                   na.rm = FALSE, # remove missing values
                   index, # c(unit, time) indicators
                   weight = NULL, # weighting
                   force = "unit", # fixed effects demeaning
                   cl = NULL,
                   r = 0, # nubmer of factors
                   lambda = NULL, # mc method: regularization parameter
                   nlambda = 10, ## mc method: regularization parameter
                   CV = TRUE, # cross-validation
                   criterion = "mspe", # mspe or pc
                   k = 5, # cross-validation times
                   EM = FALSE, # EM algorithm
                   estimator = "ife", # ife/mc method
                   se = FALSE, # report uncertainties
                   nboots = 200, # number of bootstraps
                   inference = "nonparametric", # type of inference
                   cov.ar = 1,
                   parallel = TRUE, # parallel computing
                   cores = NULL, # number of cores
                   tol = 0.001, # tolerance level
                   seed = NULL, # set seed
                   min.T0 = 5,
                   alpha = 0.05,
                   normalize = FALSE
                   ) {
    UseMethod("gsynth")
}

## formula method

gsynth.formula <- function(formula = NULL,data, # a data frame (long-form)
                           Y, # outcome
                           D, # treatment
                           X = NULL, # time-varying covariates
                           na.rm = FALSE, # remove missing values
                           index, # c(unit, time) indicators
                           weight = NULL,
                           force = "unit", # fixed effects demeaning
                           cl = NULL,
                           r = 0, # nubmer of factors
                           lambda = NULL, # mc method: regularization parameter
                           nlambda = 10, ## mc method: regularization parameter
                           CV = TRUE, # cross-validation
                           criterion = "mspe", # mspe or pc
                           k = 5, # cross-validation times
                           EM = FALSE, # EM algorithm
                           estimator = "ife", # ife/mc method
                           se = FALSE, # report uncertainties
                           nboots = 200, # number of bootstraps
                           inference = "nonparametric", # type of inference
                           cov.ar = 1,
                           parallel = TRUE, # parallel computing
                           cores = NULL, # number of cores
                           tol = 0.001, # tolerance level
                           seed = NULL, # set seed
                           min.T0 = 5,
                           alpha = 0.05,
                           normalize = FALSE
                           ) {
    ## parsing
    varnames <- all.vars(formula)
    Yname <- varnames[1]
    Dname <- varnames[2]
    if (length(varnames) > 2) {
        Xname <- varnames[3:length(varnames)]
    } else {
        Xname <- NULL
    }

    namesData <- colnames(data)
    for (i in 1:length(varnames)) {
        if(!varnames[i] %in% namesData) {
            stop(paste("variable \"", varnames[i],"\" is not in the data set.", sep = ""))
        }
    }

    ## run the model
    out <- gsynth.default(formula = NULL, data = data, Y = Yname,
                          D = Dname, X = Xname,
                          na.rm, index, weight, force, cl, r, lambda, nlambda,
                          CV, criterion, k, EM, estimator, se, nboots,
                          inference, cov.ar,
                          parallel, cores, tol, seed, min.T0, alpha,
                          normalize)

    out$call <- match.call()
    out$formula <- formula
    ## print(out)
    return(out)

}

## default function

gsynth.default <- function(formula = NULL,data, # a data frame (long-form)
                           Y, # outcome
                           D, # treatment
                           X = NULL, # time-varying covariates
                           na.rm = FALSE, # remove missing values
                           index, # c(unit, time) indicators
                           weight = NULL,
                           force = "unit", # fixed effects demeaning
                           cl = NULL,
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
                           inference = "nonparametric", # type of inference
                           cov.ar = 1,
                           parallel = TRUE, # parallel computing
                           cores = NULL, # number of cores
                           tol = 0.001, # tolerance level
                           seed = NULL, # set seed
                           min.T0 = 5,
                           alpha = 0.05,
                           normalize = FALSE
                           ) {

    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------##
    ## library(ggplot2)

    if (is.data.frame(data) == FALSE || length(class(data)) > 1) {
         data <- as.data.frame(data)
         ## warning("Not a data frame.")
    }
    ## index
    if (length(index) != 2 | sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
    }

    unique_label <- unique(paste(data[,index[1]],"_",data[,index[2]],sep=""))
    if (length(unique_label)!= dim(data)[1]) {
        stop("Some records may be duplicated or wrongly marked in the data set. Check the index.")
    }

    ## force
    if (force == "none") { # force = 0 "none": no additive fixed effects imposed
        force <- 0
    } else if (force == "unit") { # force = 1 "unit": unit fixed-effect (default)
        force <- 1
    } else if (force == "time") { # force = 2 "time": time fixed-effect
        force <- 2
    } else if (force == "two-way") { # force = 3 "two-way": two-way fixed-effect
        force <- 3
    }
    if (!force %in% c(0, 1, 2, 3)) {
        stop("\"force\" option misspecified; choose from c(\"none\", \"unit\", \"time\", \"two-way\").")
    }

    ## estimator
    if (!estimator %in% c("ife","mc")) {
        stop("\"estimator\" must be either \"ife\" or \"mc\".")
    }
    if (estimator == "mc") {
        MC <- TRUE
        ## inference <- "nonparametric"
    } else {
        MC <- FALSE
    }

    ## r
    if ( MC == FALSE & r[1] < 0) {
        stop("\"r\" option misspecified. The number of factors must be non-negative.")
    }

    ## lambda
    if ( MC == TRUE & !is.null(lambda)) {
        if (sum(lambda < 0) > 0) {
            stop("\"lambda\" option misspecified. It must be non-negative.")
        }
    }

    ## CV
    if (CV == TRUE) {
        if (MC == FALSE) {
            if (length(r) == 2 & r[1] > r[2]) {
                stop("\"r\" option misspecified.")
            }
        } else {
            if (nlambda <= 0) {
                stop("\"nlambda\" option misspecified.")
            }
        }
    } else {
        if (MC == TRUE && is.null(lambda)) {
            stop("The value of \"lambda\" should be specified.")
        }
    }

    ## criterion
    if (!criterion %in% c("mspe", "pc")) {
        stop("\"criterion\" option misspecified; choose from c(\"mspe\", \"pc\").\n")
    }

    if (length(r) == 1) {
        if (r>=5) {
            r.end <- r
        } else {
            r.end <- 5
        }
    } else {
        r.end <- r[2]; r <- r[1]
    }

    ## EM
    if (is.logical(EM) == FALSE & !EM%in%c(0, 1)) {
        stop("EM is not a logical flag.")
    }

    ## se
    if (is.logical(se) == FALSE & !se%in%c(0, 1)) {
        stop("se is not a logical flag.")
    }

    ## normalize
    if (is.logical(normalize) == FALSE & !normalize%in%c(0, 1)) {
        stop("normalize is not a logical flag.")
    }

    ## inference
    if (inference == "para") {
        inference <- "parametric"
    }
    if (inference == "nonpara") {
        inference <- "nonparametric"
    }
    if (!inference %in% c("parametric", "nonparametric", "jackknife")) {
        stop("\"inference\" option misspecified; choose from c(\"parametric\", \"nonparametric\", \"jackknife\").")
    }

    if (inference == "parametric" && !is.null(cl)) {
        cl <- NULL
        cat("\nFor clustered bootsrap, please use the nonparametric procedure.\n")
    }

    # Do not attempt to use the parametric bootstrap if the number of treated units is fewer than 40
    n_treated = length(unique(data[data[,D] ==  1, index[1]]))
    if (inference == "nonparametric" && n_treated < 40 && se) {
        stop("Nonparametric bootstrap is inappropriate when there are fewer than 40 treated units")
    }

    ## nboots
    if (se == TRUE & nboots <= 0) {
        stop("\"nboots\" option misspecified. Try, for example, nboots = 200.")
    }

    ## parallel & cores
    if (parallel == TRUE) {
        if (is.null(cores) == FALSE) {
            if (cores <= 0) {
                stop("\"cores\" option misspecified. Try, for example, cores = 2.")
            }
        }
    }

    ## tol
    if (tol <= 0) {
        stop("\"tol\" option misspecified. Try using the default option.")
    }

    ## cl
    if (alpha <= 0 || alpha >= 1) {
        stop("\"alpha\" should be in the range of 0 and 1. Try, for example, alpha = 0.05.")
    }

    ## mc inference
    if (MC == TRUE && se == 1) {
        if (inference == "parametric") {
            stop("For matrix completion method, please use nonparametric bootstrap.")
        }
    }

    ## seed
    if (is.null(seed) == FALSE) {
        if (is.numeric(seed) == FALSE) {
            stop("seed should be a number.")
        }
    }

    ## remove missing values
    if (is.logical(na.rm) == FALSE & !na.rm%in%c(0, 1)) {
        stop("na.rm is not a logical flag.")
    }

    ## data <- data[,c(index, Y, D, X, weight)] ## some variables may not be used
    ## if (na.rm == TRUE) {
      ## data <- data[,c(index, Y, D, X, weight)] ## some variables may not be used
    ##     data <- na.omit(data)
    ## }

    ## select variable that are to be used
    if (!is.null(cl)) {
        if (cl %in% index) {
            data <- data[,c(index, Y, D, X, weight)]
        } else {
            data <- data[,c(index, Y, D, X, cl, weight)]
        }
    } else {
        data <- data[,c(index, Y, D, X, weight)] ## some variables may not be used
    }

    if (na.rm == TRUE) {
      ## data <- data[,c(index, Y, D, X, weight)] ## some variables may not be used
        data <- na.omit(data)
    }

    ##-------------------------------##
    ## Parsing raw data
    ##-------------------------------##

    ##store variable names
    data.old <- data
    Yname <- Y
    Dname <- D
    Xname <- X
    Wname <- weight
    clname <- cl

    if (!is.null(clname)) {
        if (!clname %in% index) {
            data[, clname] <- as.numeric(as.factor(data[, clname]))
        }
    }

    ## normalize
    norm.para <- NULL
    if (normalize == TRUE) {
        sd.Y <- sd(as.matrix(data[,Yname]))
        data[,c(Yname, Xname)] <- data[,c(Yname, Xname)]/sd.Y
        norm.para <- sd.Y ## normalized parameter
    }

    ## check index and treatment indicator
    if (!(class(data[, Dname]) %in% c("numeric", "integer"))) {
        ## data[, Dname] <- as.numeric(as.character(data[, Dname]))
        stop("Treatment indicator should be a numberic value.")
    }

    if (class(data[, index[1]]) == "factor") {
        data[, index[1]] <- as.character(data[, index[1]])
    }

    if (class(data[, index[2]]) == "factor") {
        data[, index[2]] <- as.character(data[, index[2]])
    }

    id <- index[1]
    time <- index[2]
    TT <- length(unique(data[,time]))
    N <- length(unique(data[,id]))
    p <- length(Xname)
    id.series <- unique(sort(data[,id])) ## unit id
    time.uni <- unique(sort(data[,time])) ## period

    ## sort data
    data <- data[order(data[,id], data[,time]), ]

    ## check missingness
    if (sum(is.na(data[, Yname])) > 0) {
        stop(paste("Missing values in variable \"", Yname,"\".", sep = ""))
    }
    if (sum(is.na(data[, Dname])) > 0) {
        stop(paste("Missing values in variable \"", Dname,"\".", sep = ""))
    }

    if (!(1%in%data[, Dname] & 0%in%data[,Dname] & length(unique(data[,Dname])) == 2)) {
        stop(paste("Error values in variable \"", Dname,"\".", sep = ""))
    }

    if (p > 0) {
        for (i in 1:p) {
            if (sum(is.na(data[, Xname[i]])) > 0) {
                stop(paste("Missing values in variable \"", Xname[i],"\".", sep = ""))
            }
            if (sum(tapply(data[, Xname[i]], data[, id], var), na.rm = TRUE) == 0) {
                stop(paste("Variable \"",Xname[i], "\" is unit-invariant. Try to remove it.", sep = ""))
            }
            if (sum(tapply(data[, Xname[i]], data[, time], var), na.rm = TRUE) == 0) {
                stop(paste("Variable \"",Xname[i], "\" is time-invariant. Try to remove it.", sep = ""))
            }
        }
    }

    if (sum(is.na(data[, id])) > 0) {
        stop(paste("Missing values in variable \"", id,"\".", sep = ""))
    }
    if (sum(is.na(data[, time])) > 0) {
        stop(paste("Missing values in variable \"", time,"\".", sep = ""))
    }

    ## check balanced panel and fill unbalanced panel
    if (dim(data)[1] < TT*N) {
        data[,time] <- as.numeric(as.factor(data[,time]))
        ## ob <- "time_ob_ls"

        ## while (ob %in% colnames(data)) {
        ##     ob <- paste(ob, ob, sep = "_")
        ## }

        ## data[, ob] <- data[, time]
        ## for (i in 1:N) {
        ##     data[data[,id] == id.series[i], ob] <- data[data[,id] == id.series[i],time] + (i - 1) * TT
        ## }

        ob.indicator <- data[,time]
        id.indicator <- table(data[, id])
        sub.start <- 1
        for (i in 1:(N - 1)) {
            sub.start <- sub.start + id.indicator[i]
            sub.end <- sub.start + id.indicator[i+1] - 1
            ob.indicator[sub.start:sub.end] <- ob.indicator[sub.start:sub.end] + i * TT
        }

        variable <- c(Yname, Dname, Xname, Wname, clname)

        data_I <- matrix(0, N * TT, 1)
        data_I[ob.indicator, 1] <- 1
        data_ub <- as.matrix(data[, variable])
        data <- data_ub_adj(data_I, data_ub)
        colnames(data) <- variable
    }

    ## index matrix that indicates if data is observed
    I <- matrix(1, TT, N)
    Y.ind <- matrix(data[, Yname], TT, N)
    I[is.nan(Y.ind)] <- 0
    I.old <- I

    if (0%in%I) {
        data[is.nan(data)] <- 0
    }

    ##treatment indicator
    D <- matrix(data[, Dname], TT, N)

    ## once treated, always treated
    ## careful unbalanced case: calculate treated units
    D <- apply(D, 2, function(vec){cumsum(vec)})
    D <- ifelse(D > 0, 1, 0)

    ## weighting variable
    if (is.null(Wname)) {
        W <- NULL
    } else {
        W <- matrix(data[, Wname], TT, N)
    }

    ##outcome variable
    Y <- matrix(data[, Yname], TT, N)

    ## time-varying covariates
    X <- array(0, dim = c(TT, N, p))
    ## xp <- rep(0, p) ## label invariant x
    ## x.pos <- 0

    if (p > 0) {
        ## x.pos <- 1:p
        for (i in 1:p) {
            X[,,i] <- matrix(data[, Xname[i]], TT, N)
            ## if (force %in% c(1,3)) {
            ##     if (!0%in%I) {
            ##         tot.var.unit <- sum(apply(X[, , i], 2, var))
            ##     } else {
            ##         Xi <- X[,,i]
            ##         Xi[which(I == 0)] <- NA
            ##         tot.var.unit <- sum(apply(Xi, 2, var, na.rm = TRUE))
            ##     }
            ##     if(!is.na(tot.var.unit)) {
            ##         if (tot.var.unit == 0) {
                        ## time invariant covar can be removed
            ##             xp[i] <- 1
            ##             cat(paste("Variable \"", Xname[i],"\" is time-invariant.\n", sep = ""))
            ##         }
            ##     }
            ## }
            ## if (force %in% c(2, 3)) {
            ##     if (!0%in%I) {
            ##         tot.var.time <- sum(apply(X[, , i], 1, var))
            ##     } else {
            ##         Xi <- X[,,i]
            ##         Xi[which(I == 0)] <- NA
            ##         tot.var.time <- sum(apply(Xi, 1, var, na.rm = TRUE))
            ##     }
            ##     if (!is.na(tot.var.time)) {
            ##         if (tot.var.time == 0) {
                        ## can be removed in inter_fe
            ##             xp[i] <- 1
            ##             cat(paste("Variable \"", Xname[i],"\" has no cross-sectional variation.\n", sep = ""))
            ##         }
            ##     }
            ## }
        }
    }

    if (!is.null(clname)) {
        if (clname %in% index) {
            cl <- 1:N
        } else {
            cl <- matrix(data[, clname], TT, N)
            v.cl <- c()
            for (i in 1:N) {
                if (sum(is.na(cl[,i])) > 0) {
                    v.cl <- c(v.cl, na.omit(cl[,i])[1])
                } else {
                    v.cl <- c(v.cl, cl[1, i])
                }
            }
            cl <- v.cl
        }
    } else {
        cl <- NULL
    }

    tr <- D[TT,] == 1     # cross-sectional: treated unit
    ## ----------------------------------------------------------- ##
    II <- I
    II[which(D==1)] <- 0 ## regard treated values as missing

    ## 1. remove units that have too few observations
    T0 <- apply(II, 2, sum)
    T0.min <- min(T0)

    if (sum(T0[which(tr == 1)] >= min.T0) == 0) {
        stop ("All treated units have been removed.\n")
    }
    ## T0.min : minimum T0,  min.T0: manually set
    ## rm.tr.id: relative location of treated units (within all treated units)
    ## that will be removed
    if (T0.min < min.T0) {
        cat("Some treated units has too few pre-treatment periods. \nThey will be automatically removed.\n")
    }

    rm.id <- sort(which(T0 < min.T0))
    ## rm.id <- which(T0 < min.T0) ## removed id
    ## rem.id <- which(T0 >= min.T0) ## remaining id
    ## rem.id <- setdiff(1:N, rm.id)

    if (length(rm.id) == N) {
        stop("All units have been removed.\n")
    }

    if (length(rm.id) > 0) {
        X.old <- X
        if (p > 0) {
            X <- array(0,dim = c(TT, (N - length(rm.id)), p))
            for (i in 1:p) {
                subX <- X.old[, , i]
                X[, , i] <- as.matrix(subX[, -rm.id])
            }
        } else {
            X <- array(0,dim = c(TT, (N - length(rm.id)), 0))
        }

        # N <- N - length(rm.id)
        Y <- as.matrix(Y[,-rm.id])
        D <- as.matrix(D[,-rm.id])
        I <- as.matrix(I[,-rm.id]) ## after removing
        II <- as.matrix(II[,-rm.id])
        T0 <- T0[-rm.id]
        ## tr <- tr[-rm.id]
        if (!is.null(cl)) {
            cl <- cl[-rm.id]
        }
    }

    ## 2. check if some periods when all units are missing
    I.use <- apply(II, 1, sum)
    if (0%in%I.use) {
        for (i in 1:TT) {
            if (I.use[i] == 0) {
                cat("There are not any observations under control at ",time.uni[i],", drop that period.\n")
            }
        }
        TT <- TT - sum(I.use == 0)
        time.uni <- time.uni[-which(I.use == 0)]

        I <- I[-which(I.use == 0),] ## remove that period
        II <- II[-which(I.use == 0),] ## remove that period
        D <- D[-which(I.use == 0),] ## remove that period
        Y <- Y[-which(I.use == 0),] ## remove that period

        X.old <- X
        if (p > 0) {
            X <- array(0,dim = c(TT, (N - length(rm.id)), p))
            for (i in 1:p) {
                subX <- X.old[, , i]
                X[, , i] <- as.matrix(subX[-which(I.use == 0),])
            }
        } else {
            X <- array(0,dim = c(TT, (N - length(rm.id)), 0))
        }
    }

    ## 3. check candidate factor numbers for cross-validation
    if (MC == FALSE) { ## factor model
        T0.min.2 <- min(T0)
        if ( (length(r) == 1) & (!CV) ) {
            con1 <- (T0.min.2 < r.end) & (force%in%c(0,2))
            con2 <- (T0.min.2 <= r.end) & (force%in%c(1,3))
            if (con1) {
                T0.min.e <- r.end
            }
            if (con2) {
                T0.min.e <- r.end + 1
            }
            if (con1 | con2) {
                stop("Some treated units has too few pre-treatment periods. Please set a larger value for min.T0 to remove them. Equal or greater than ",T0.min.e," is recommended.\n")
            }
        }
        if (CV) {
            con1 <- (T0.min.2 < r.end + 1) & (force%in%c(0,2))
            con2 <- (T0.min.2 < r.end + 2) & (force%in%c(1,3))
            if (con1) {
                T0.min.e <- r.end + 1
            }
            if (con2) {
                T0.min.e <- r.end + 2
            }
            if (con1 | con2) {
                stop("Some treated units has too few pre-treatment periods. Please set a larger value for min.T0 to remove them. Equal or greater than ",T0.min.e," is recommended. Or you can set a smaller range of factor numbers.\n")
            }
        }
    }

    ## if (is.null(X) == TRUE) {
    ##     X <- array(0, dim = c(1, 1, 0))
    ##     p <- 0
    ## } else {
    ##     p <- dim(X)[3]
    ## }
    ## for AR1, burn the first period
    AR1 <- FALSE
    ## if (AR1 == TRUE) {
    ##     Y.first <- Y[1,]
    ##     Y.lag <- Y[1:(T-1),]
    ##     Y <- Y[2:T,]
    ##     D <- D[2:T,]
    ##     if (p == 0) {
    ##         X <- array(NA, dim=c((T-1),N,1))
    ##         X[,,1] <- Y.lag
    ##     } else {
    ##         X.first <- X[1,,]
    ##         X.sav <- X[2:T,,]
    ##         X <- array(NA,dim=c((T-1),N,(p+1)))
    ##         X[,,1] <- Y.lag
    ##         X[,,2:(p+1)] <- X.sav
    ##     }
    ##     T <- T-1
    ## }

    ##-------------------------------##
    ## Register clusters
    ##-------------------------------##

    if (se == TRUE & parallel==TRUE) {

        if (is.null(cores) == TRUE) {
            cores <- detectCores()
        }

        para.clusters <- future::makeClusterPSOCK(cores)
        registerDoParallel(para.clusters)
        if (is.null(seed) == FALSE) {
            registerDoRNG(seed)
        }


        cat("Parallel computing ...\n")
    }

    ##-------------------------------##
    ## run main program
    ##-------------------------------##

    if (se == FALSE) {
        if (MC == FALSE) {
            if (EM == FALSE) { # the algorithm suggested in the paper
                out <- synth.core(Y = Y, X = X, D = D, I = I, W = W,
                                  r = r, r.end = r.end, force = force,
                                  CV = CV, criterion = criterion, tol = tol,
                                  AR1 = AR1, norm.para = norm.para)

            } else { # EM algorithm
                if (CV == FALSE) {
                    out <- synth.em(Y = Y, X = X, D = D, I = I, W = W,
                                    r = r, force = force,
                                    tol = tol, AR1 = AR1, norm.para = norm.para)


                } else { # cross-validation
                    out <- synth.em.cv(Y = Y,X = X, D = D, I = I, W = W,
                                       r = r, r.end = r.end,
                                       criterion = criterion, force = force,
                                       tol = tol, AR1 = AR1, norm.para = norm.para)

                }
            }
        } else {
            out <- synth.mc(Y = Y, X = X, D = D, I = I, W = W, lambda = lambda,
                            nlambda = nlambda, force = force, CV = CV, k = k,
                            tol = tol, AR1 = AR1, norm.para = norm.para)
        }
    } else  {
        if (is.null(seed) == FALSE) {
            set.seed(seed)
        }
        out <- synth.boot(Y = Y, X = X, D = D, cl = cl, I=I, W = W, EM = EM,
                          r = r, r.end = r.end, lambda = lambda,
                          nlambda = nlambda, force = force,
                          CV = CV, criterion = criterion, k = k,
                          tol = tol, MC = MC,
                          nboots = nboots, inference = inference,
                          cov.ar = cov.ar,
                          parallel = parallel, cores = cores,
                          AR1 = AR1, norm.para = norm.para,
                          alpha = alpha)

    }

    if (se == TRUE & parallel == TRUE) {
        stopCluster(para.clusters)
        ##closeAllConnections()
    }

    if ( (out$validX == 0) & (p!=0) ) {
        warning("Multi-colinearity among covariates. Try removing some of them.\r")
    }


    ##-------------------------------##
    ## storage
    ##-------------------------------##

    iname.old <- iname <- unique(sort(data.old[, id]))
    ## tname.old <- tname <- unique(sort(data.old[,time]))
    if (!0%in%I.use) {
        tname.old <- tname <- unique(sort(data.old[, time]))
    } else {
        tname.old <- tname <- unique(sort(data.old[, time]))[which(I.use != 0)]
    }

    if (length(rm.id) > 0) {
        ## tr.remove.id <- iname[rm.tr.id.s]
        iname <- iname[-rm.id]
    }

    ## remaining control and treated units
    id.tr <- which(tr == 1)
    id.co <- which(tr == 0)
    if (length(rm.id) > 0) {
        id.tr <- setdiff(id.tr, rm.id)
        id.co <- setdiff(id.co, rm.id)
    }

    obs.missing <- matrix(1, TT, N) ## control group:1

    tr.pre <- out$pre
    tr.post <- out$post

    tr.pre[which(tr.pre == 1)] <- 2 ## pre 2
    tr.post[which(tr.post == 1)] <- 3 ## post 3
    obs.missing[, id.tr] <- tr.pre + tr.post

    if (length(rm.id) > 0) {
        obs.missing[which(I.old == 0)] <- 0 ## I: after removing I.old: total
        obs.missing[, rm.id] <- 4 ## removed 4
    } else {
        obs.missing[which(I == 0)] <- 0 ## missing 0 ## I: total
    }
    ## obs.missing[which(obs.missing==1)] <- "control"
    ## obs.missing[which(obs.missing==2)] <- "pre"
    ## obs.missing[which(obs.missing==3)] <- "post"
    ## obs.missing[which(obs.missing==4)] <- "removed"
    ## obs.missing[which(obs.missing==0)] <- "missing"

    if (!is.null(norm.para)) {
        Y <-  Y * norm.para[1]
    }

    colnames(obs.missing) <- unique(sort(data.old[, id]))
    colnames(Y) <- iname
    if (!is.null(out$res.co)) {
        colnames(out$res.co) <- iname[which(out$tr == 0)]
        rownames(out$res.co) <- tname
    } else {
        colnames(out$res) <- iname
        rownames(out$res) <- tname
    }
    colnames(out$Y.co) <- iname[which(out$tr == 0)]
    colnames(out$Y.ct) <- colnames(out$Y.tr) <- colnames(out$I.tr) <- colnames(out$D.tr) <- colnames(out$post) <- colnames(out$pre) <- iname[which(out$tr == 1)]
    rownames(out$Y.ct) <- rownames(out$Y.tr) <- rownames(out$I.tr) <- rownames(out$D.tr) <- rownames(out$Y.co) <- rownames(out$post) <- rownames(out$pre) <- rownames(out$Y.bar) <- rownames(Y) <- rownames(obs.missing) <- tname


    if (AR1 == TRUE) {
        tname <- tname[-1]
    }
    Xname.tmp <- Xname
    if (AR1 == TRUE) {
        Xname.tmp <- c(paste(Yname, "_lag", sep=""), Xname)
    }

    ##  ----------- add label ---------- ##
    ## beta
    rownames(out$beta) <- Xname.tmp
    if (se == TRUE) {
        rownames(out$est.beta) <- Xname.tmp
        rownames(out$beta.boot) <- Xname.tmp
        if (class(out$eff.boot) == "array") {
            dimnames(out$eff.boot)[[2]] <- iname[which(out$tr == 1)]
        }

    }
    ## eff
    colnames(out$eff) <- iname[which(out$tr == 1)]
    rownames(out$eff) <- tname
    if (!is.null(out$eff.cnt)) {
        colnames(out$eff.cnt) <- iname[which(out$tr == 1)]
    }
    if (out$sameT0 == TRUE) {
        names(out$att) <- tname
    }

    ## individual eff
    if (!is.null(out$est.ind)) {
        dimnames(out$est.ind)[[3]] <- iname[which(out$tr == 1)]
    }

    ## cross validation
    if (!is.null(out$CV.out)) {
        rownames(out$CV.out) <- rep("", dim(out$CV.out)[1])
    }
    ## ife
    if (!is.null(out$factor)) {
        rownames(out$factor) <- tname
        rownames(out$lambda.tr) <- iname[which(out$tr == 1)]
        rownames(out$lambda.co) <- iname[which(out$tr == 0)]
        colnames(out$lambda.tr) <- colnames(out$lambda.co) <- colnames(out$factor) <- sapply(1:dim(out$factor)[2], function(i){paste("r", i, sep = "")})
    }

    ## add fe
    if (!is.null(out$xi)) {
        rownames(out$xi) <- tname
        colnames(out$xi) <- ""
    }

    if (!is.null(out$alpha.tr)) {
        ## if (class(out$alpha.tr) != "matrix") {
        ##     out$alpha.tr <- as.matrix(out$alpha.tr)
        ## }
        rownames(out$alpha.tr) <- iname[which(out$tr == 1)]
        colnames(out$alpha.tr) <- ""
    }

    if (!is.null(out$alpha.co)) {
        rownames(out$alpha.co) <- iname[which(out$tr == 0)]
        colnames(out$alpha.co) <- ""
    }

    ## group
    names(out$tr) <- iname

    if (MC == FALSE) {
        if (out$r.cv > 0) {
            colnames(out$wgt.implied) <- iname[which(out$tr == 1)]
            rownames(out$wgt.implied) <- iname[which(out$tr == 0)]
        }
    }

    output <- c(list(Y.dat = Y,
                     Y = Yname,
                     D = Dname,
                     X = Xname,
                     W = Wname,
                     index = index,
                     id = iname,
                     time = tname,
                     obs.missing = obs.missing,
                     id.tr = iname[which(out$tr == 1)],
                     id.co = iname[which(out$tr == 0)]),
                     out)

    if (!is.null(Wname)) {
        output <- c(output, list(W = Wname))
    }

    if (length(rm.id) > 0) {
        removed.id <- iname.old[rm.id]
        output <- c(output,list(removed.id = removed.id))
        cat("list of removed units:", removed.id)
        cat("\n\n")
    }
    output <- c(output, list(call = match.call()))
    class(output) <- "gsynth"
    return(output)

} ## Program GSynth ends
