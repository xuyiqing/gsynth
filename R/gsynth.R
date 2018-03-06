## Synthetic Control for Multiple Treated Units
## (Causal Inference with Interactive Fixed Effects Models)
## Version 1.03
## Author: Yiqing Xu, University of California, San Diego
## Date: 2018.1.20

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

## preView()

#####################################################################
## A Shell Function
#####################################################################

## generic function

gsynth <- function(formula = NULL,data, # a data frame (long-form)
                   Y, # outcome
                   D, # treatment 
                   X = NULL, # time-varying covariates
                   na.rm = FALSE, # remove missing values
                   index, # c(unit, time) indicators
                   weight = NULL, # weighting
                   force = "unit", # fixed effects demeaning
                   r = 0, # nubmer of factors
                   lambda = NULL, # mc method: regularization parameter
                   nlambda = 10, ## mc method: regularization parameter
                   CV = TRUE, # cross-validation
                   EM = FALSE, # EM algorithm
                   MC = FALSE, # mc method
                   se = FALSE, # report uncertainties
                   nboots = 200, # number of bootstraps
                   inference = "nonparametric", # type of inference
                   cov.ar = 1,
                   parallel = FALSE, # parallel computing
                   cores = NULL, # number of cores
                   tol = 0.001, # tolerance level
                   seed = NULL, # set seed
                   min.T0 = 5,
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
                           r = 0, # nubmer of factors
                           lambda = NULL, # mc method: regularization parameter
                           nlambda = 10, ## mc method: regularization parameter
                           CV = TRUE, # cross-validation
                           EM = FALSE, # EM algorithm
                           MC = FALSE, # mc method 
                           se = FALSE, # report uncertainties
                           nboots = 200, # number of bootstraps
                           inference = "nonparametric", # type of inference
                           cov.ar = 1,
                           parallel = FALSE, # parallel computing
                           cores = NULL, # number of cores
                           tol = 0.001, # tolerance level
                           seed = NULL, # set seed
                           min.T0 = 5,
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
    ## run the model
    out <- gsynth.default(formula = NULL, data = data, Y = Yname,
                          D = Dname, X = Xname,
                          na.rm, index, weight, force, r, lambda, nlambda, 
                          CV, EM, MC, se, nboots, 
                          inference, cov.ar, 
                          parallel, cores, tol, seed, min.T0, 
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
                           r = 0, # nubmer of factors
                           lambda = NULL, ## mc method: regularization parameter
                           nlambda = 10, ## mc method: regularization parameter
                           CV = TRUE, # cross-validation
                           EM = FALSE, # EM algorithm 
                           MC = FALSE, # mc method
                           se = FALSE, # report uncertainties
                           nboots = 200, # number of bootstraps
                           inference = "nonparametric", # type of inference
                           cov.ar = 1,
                           parallel = FALSE, # parallel computing
                           cores = NULL, # number of cores
                           tol = 0.001, # tolerance level
                           seed = NULL, # set seed
                           min.T0 = 5,
                           normalize = FALSE
                           ) {  
    
    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------## 
    ## library(ggplot2)

    if (is.data.frame(data) == FALSE) {
        data <- as.data.frame(data)
        warning("Not a data frame.")
    }
    ## index
    if (length(index) != 2 | sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
    }
    
    unique_label <- unique(paste(data[,index[1]],"_",data[,index[2]],sep=""))
    if (length(unique_label)!= dim(data)[1]) {
        stop("Some records may be replicated or wrongly marked in the data set.")
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

    ## MC
    if (is.logical(MC) == FALSE & !MC%in%c(0, 1)) {
        stop("MC is not a logical flag.")
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
    if (!inference %in% c("parametric", "nonparametric")) {
        stop("\"inference\" option misspecified; choose from c(\"parametric\", \"nonparametric\").")
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

    ## mc inference
    if (MC == TRUE && se == 1) {
        if (inference=="parametric") {
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

    if (na.rm == TRUE) {
    	data <- data[,c(index, Y, D, X, weight)] ## some variables may not be used
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

    ## normalize
    norm.para <- NULL
    if (normalize == TRUE) {
        sd.Y <- sd(as.matrix(data[,Yname]))
        data[,c(Yname, Xname)] <- data[,c(Yname, Xname)]/sd.Y
        norm.para <- sd.Y ## normalized parameter
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
        }
    }

    if (sum(is.na(data[, id])) > 0) {
        stop(paste("Missing values in variable \"", id,"\".", sep = ""))
    }
    if (sum(is.na(data[, time])) > 0) {
        stop(paste("Missing values in variable \"", time,"\".", sep = ""))
    } 

    ## check balanced panel and fill unbalanced panel
    if (var(table(data[,id])) + var(table(data[, time])) > 0 | TT == N) {
        
        data[,time]<-as.numeric(as.factor(data[,time]))
        ob <- "time_ob_ls"
        
        while(ob%in%colnames(data)){
            ob <- paste(ob,ob,sep="")
        }

        data[,ob]<-data[,time]
        for (i in 1:N) {
            data[data[,id] == id.series[i], ob] <- data[data[,id] == id.series[i],time] + (i - 1) * TT  
        }

        variable <- c(Yname, Dname, Xname, Wname)

        data_I <- matrix(0, N * TT, 1)
        data_I[c(data[,ob]), 1] <- 1
        data_ub <- as.matrix(data[, variable])
        data <- data_ub_adj(data_I, data_ub)
        colnames(data) <- variable
    }

    ## index matrix that indicates if data is observed 
    I <- matrix(1, TT, N)
    Y.ind <- matrix(data[, Yname], TT, N)
    I[is.nan(Y.ind)] <- 0

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
    tr <- D[TT,] == 1     # cross-sectional: treated unit
    
    I.tr.use <- apply(as.matrix(I[, which(!tr)]), 1, sum) ## check if at some periods all control units are missing
    if (0%in%I.tr.use) {
        for (i in 1:TT) {
            if (I.tr.use[i] == 0) {
                cat("There are not any observations in control group at ",time.uni[i],", drop observations in treated group units at that period.\n")
            }
        }
        TT <- TT - sum(I.tr.use == 0)
        time.uni <- time.uni[-which(I.tr.use == 0)]
        
        I <- I[-which(I.tr.use == 0),] ## remove that period
        D <- D[-which(I.tr.use == 0),] ## remove that period
        Y <- Y[-which(I.tr.use == 0),] ## remove that period

        ## Y[which(I==0)] <- 0
        data <- data[-which(rep(I.tr.use, N) == 0),] ## remove that period
    }

    pre <- as.matrix(D[,which(tr == 1)] == 0 & I[,which(tr == 1)] == 1) # a matrix indicating before treatment
    T0 <- apply(pre, 2, sum) 
    T0.min <- min(T0)
    id.tr <- which(tr == 1)
    id.co <- which(tr == 0)
    
    if (MC == FALSE) { ## factor model
        if ( (length(r) == 1) & (!CV) ) {
 
            if (sum(T0 >= min.T0) == 0) {
                stop ("All treated units have been removed.\n")
            }       
 
            T0.min.2 <- min(T0[which(T0 >= min.T0)])

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
        
            if (sum(T0 >= min.T0) == 0) {
                stop ("All treated units have been removed.\n")
            }

            T0.min.2 <- min(T0[which(T0 >= min.T0)])

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

    ## T0.min : minimum T0,  min.T0: manually set
    ## rm.tr.id: relative location of treated units (within all treated units) 
    ## that will be removed 
    if (T0.min < min.T0) {
        cat("Some treated units has too few pre-treatment periods. \nThey will be automatically removed.\n")
    }

    rm.tr.id <- rep(0,dim(pre)[2])
    for (i in 1:dim(pre)[2]) {
        if(T0[i] < min.T0) {
            rm.tr.id[i] <- 1
        }
    }

    if (sum(rm.tr.id) == dim(pre)[2]) {
        stop("All treated units have been removed.\n")
    }

    if (1 %in% rm.tr.id) {
        tr.id <- which(tr == 1) ## tr.id: location of all treated units
        rm.tr.pos <- which(rm.tr.id == 1)
        rm.tr.id.s <- tr.id[rm.tr.pos] ## location of treated units that will be removed
        id.tr <- id.tr[-rm.tr.pos] ## remaining treated units
    }
    
    ## time-varying covariates
    X <- array(0, dim = c(TT, N, p))

    if (p > 0) {
        for (i in 1:p) {
            X[,,i] <- matrix(data[, Xname[i]], TT, N)
            if (force %in% c(1,3)) {
                if (!0%in%I) {
                    tot.var.unit <- sum(apply(X[, , i], 2, var))
                } else {
                    Xi <- X[,,i]
                    Xi[which(I == 0)] <- NA
                    tot.var.unit <- sum(apply(Xi, 2, var, na.rm = TRUE))
                }
                if(!is.na(tot.var.unit)) {
                    if (tot.var.unit == 0) {
                        ## time invariant covar can be removed
                        cat(paste("Variable \"", Xname[i],"\" is time-invariant.\n", sep = ""))   
                    }
                }
            }
            if (force %in% c(2, 3)) {
                if (!0%in%I) {
                    tot.var.time <- sum(apply(X[, , i], 1, var))
                } else {
                    Xi <- X[,,i]
                    Xi[which(I == 0)] <- NA
                    tot.var.time <- sum(apply(Xi, 1, var, na.rm = TRUE))
                } 
                if (!is.na(tot.var.time)) {
                    if (tot.var.time == 0) {
                        ## can be removed in inter_fe
                        cat(paste("Variable \"", Xname[i],"\" has no cross-sectional variation.\n", sep = ""))
                    }
                }
            } 
        } 
    }

    if (1 %in% rm.tr.id) {

        X.old <- X
        if (p > 0) {
            X <- array(0,dim = c(TT, (N - length(rm.tr.id.s)), p))
            for (i in 1:p) {
                subX <- X.old[, , i]
                X[, , i] <- as.matrix(subX[, -rm.tr.id.s])
            }
        } else {
            X <- array(0,dim = c(TT, (N - length(rm.tr.id.s)), 0))
        }

        Y <- as.matrix(Y[,-rm.tr.id.s])
        D <- as.matrix(D[,-rm.tr.id.s])
        I.old <- I ## total I
        I <- as.matrix(I[,-rm.tr.id.s]) ## after removing
    }    

    if (is.null(dim(X)[3]) == TRUE) {
        p <- 0
    } else {
        p <- dim(X)[3]
    }
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
        para.clusters <- makeCluster(cores)
        registerDoParallel(para.clusters)
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
                                  CV = CV, tol = tol, 
                                  AR1 = AR1, norm.para = norm.para)

            } else { # EM algorithm
                if (CV == FALSE) { 
                    out <- synth.em(Y = Y, X = X, D = D, I = I, W = W,
                                    r = r, force = force,
                                    tol = tol, AR1 = AR1, norm.para = norm.para)

                
                } else { # cross-validation
                    out <- synth.em.cv(Y = Y,X = X, D = D, I = I, W = W,
                                       r = r, r.end = r.end, force = force,
                                       tol = tol, AR1 = AR1, norm.para = norm.para)

                } 
            }
        } else {
            out <- synth.mc(Y = Y, X = X, D = D, I = I, W = W, lambda = lambda,
                            nlambda = nlambda, force = force, CV = CV,
                            tol = tol, AR1 = AR1, norm.para = norm.para)
        } 
    } else  {
        if (is.null(seed) == FALSE) {
            set.seed(seed)
        }
        out <- synth.boot(Y = Y, X = X, D = D, I=I, W = W, EM = EM,
                          r = r, r.end = r.end, lambda = lambda,
                          nlambda = nlambda, force = force,
                          CV = CV, tol = tol, MC = MC,
                          nboots = nboots, inference = inference,
                          cov.ar = cov.ar,
                          parallel = parallel, cores = cores,           
                          AR1 = AR1, norm.para = norm.para)

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
    
    iname.old <- iname <- unique(sort(data.old[,id]))
    ## tname.old <- tname <- unique(sort(data.old[,time]))
    if (!0%in%I.tr.use) {
        tname.old <- tname <- unique(sort(data.old[,time]))
    } else {
        tname.old <- tname <- unique(sort(data.old[,time]))[which(I.tr.use != 0)]
    }

    if (1%in%rm.tr.id) {
        tr.remove.id <- iname[rm.tr.id.s]
        iname <- iname[-rm.tr.id.s]
    }

    obs.missing <-matrix(1, TT, N) ## control group:1
    
    tr.pre <- out$pre
    tr.post <- out$post

    tr.pre[which(tr.pre == 1)] <- 2 ## pre 2
    tr.post[which(tr.post == 1)] <- 3 ## post 3
    obs.missing[,id.tr] <- tr.pre + tr.post

    if (1%in%rm.tr.id) {
        obs.missing[which(I.old == 0)] <- 0 ## I: after removing I.old: total
        obs.missing[,rm.tr.id.s] <- 4 ## removed 4 
    } else {
        obs.missing[which(I == 0)] <- 0 ## missing 0 ## I: total    
    }
    ## obs.missing[which(obs.missing==1)] <- "control"
    ## obs.missing[which(obs.missing==2)] <- "pre"
    ## obs.missing[which(obs.missing==3)] <- "post"
    ## obs.missing[which(obs.missing==4)] <- "removed"
    ## obs.missing[which(obs.missing==0)] <- "missing"

    colnames(obs.missing) <- unique(sort(data.old[,id]))
    
    rownames(obs.missing) <- tname
    
    if (AR1 == TRUE) {
        tname <- tname[-1]
    } 
    Xname.tmp <- Xname
    if (AR1 == TRUE) {
        Xname.tmp <- c(paste(Yname, "_lag", sep=""), Xname)
    }
    rownames(out$beta) <- Xname.tmp
    if (se == TRUE) {
        rownames(out$est.beta) <- Xname.tmp
    } 
    colnames(out$eff) <- iname[which(out$tr == 1)]
    rownames(out$eff) <- tname

    if (MC == FALSE) {
        if (out$r.cv>0) {
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
                
    if (1 %in% rm.tr.id) {
        output <- c(output,list(tr.remove.id = tr.remove.id))
        cat("list of removed treated units:",tr.remove.id)
        cat("\n\n")
    }
    output <- c(output, list(call = match.call()))
    class(output) <- "gsynth"
    return(output)
    
} ## Program GSynth ends 


###################################################################
## Core Function
###################################################################
synth.core<-function(Y, # Outcome variable, (T*N) matrix
                     X, # Explanatory variables:  (T*N*p) array
                     D, #  Indicator for treated unit (tr==1) 
                     I,
                     W = NULL,
                     r=0, # initial number of factors considered if CV==1
                     r.end,
                     force,
                     CV = 1, # cross-validation
                     tol, # tolerance level
                     AR1 = 0,
                     beta0 = NULL, # starting value 
                     norm.para = NULL,
                     boot = 0) { # bootstrap: needn't to calculate weight  
    
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  
    na.pos <- NULL
    ## unit id and time
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(X) == FALSE) {p <- dim(X)[3]} else {p <- 0}

     
    ## treatement indicator
    tr <- D[TT,] == 1  ## cross-sectional: treated unit
    co <- D[TT,] == 0
    I.tr <- as.matrix(I[,tr]) ## maybe only 1 treated unit
    I.co <- I[,co]

    if (!0%in%I.tr) {
        ## a (TT*Ntr) matrix, time dimension: before treatment
        pre <- as.matrix(D[,which(tr == 1)] == 0)
        post <- as.matrix(D[,which(tr == 1)] == 1)   
    } else {
        pre <- as.matrix(D[,which(tr == 1)] == 0 & I[,which(tr == 1)] == 1)
        post <- as.matrix(D[,which(tr == 1)] == 1 & I[,which(tr == 1)] == 1)
    }

    D.tr <- as.matrix(D[,which(tr == 1)])
    T0.ub <- apply(D.tr == 0, 2, sum) 
    T0.ub.min <- min(T0.ub) ## unbalanced data

    if (!is.null(W)) {
        W.tr <- as.matrix(W[,which(tr == 1)])
    }

    Ntr <- sum(tr)
    Nco <- N - Ntr
    ## careful: only valid for balanced panel
    T0 <- apply(pre, 2, sum) 
    T0.min <- min(T0)
    sameT0 <- length(unique(T0)) == 1 ## treatment kicks in at the same time 
                                      ## unbalanced case needs more conditions
    if (!0%in%I.tr) {
        DID <- sameT0
    } else {
        DID <- length(unique(T0.ub)) == 1
    }
    
    id <- 1:N
    time <- 1:TT
    id.tr <- which(tr == 1) ## treated id
    id.co <- which(tr == 0)
    
    pre.v <- as.vector(pre)  ## vectorized "pre-treatment" indicator
    id.tr.pre.v <- rep(id, each = TT)[which(pre.v == 1)]  ## vectorized pre-treatment grouping variable for the treated
    time.pre <- split(rep(time, Ntr)[which(pre.v == 1)], id.tr.pre.v) ## a list of pre-treatment periods

    ## parsing data
    Y.tr <- as.matrix(Y[,id.tr])
    Y.co <- as.matrix(Y[,id.co])

    if (p == 0) {
        X.tr <- array(0, dim = c(TT, Ntr, 0))
        X.co <- array(0, dim = c(TT, Nco, 0)) 
    } else {
        X.tr <- array(NA, dim=c(TT, Ntr, p))
        X.co <- array(NA, dim=c(TT, Nco, p))
        for (j in 1:p) {
            X.tr[, , j] <- X[, id.tr, j]
            X.co[, , j] <- X[, id.co, j]
        }
    } 

    if (is.null(beta0) == TRUE ) {
        beta0 <- matrix(0, p, 1)
    }
    
    ##-------------------------------##
    ## Main Algorithm
    ##-------------------------------##

    validX <- 1 ## no multi-colinearity
    
    if (CV == FALSE) { ## case: CV==0
        
        ## inter.fe on the control group
        if(!0%in%I.co){
            ## if (force!=0) {
                est.co.best <- inter_fe(Y.co, X.co, r, force = force, beta0 = beta0, tol)
            ## } else {
                ## est.co.best<-inter_fe(Y.co, abind(I.co,X.co,along=3), r, force=0, beta0 = beta0)  
            ## }    
        } else {
            ## if (force!=0) {
              est.co.best <- inter_fe_ub(Y.co, X.co, I.co, r, force = force, tol)
            ## } else {
            ##   est.co.best<-inter_fe_ub(Y.co, abind(I.co,X.co,along=3), I.co, r, force=0, beta0 = beta0)  
            ## }
        }
                
        if (p > 0) {
            na.pos <- is.nan(est.co.best$beta)
            
        } 
        r.cv <- r.max <- r
        
    }  else if (CV == TRUE) { 
        
        ##-------------------------------##
        ## Cross-validation of r
        ##-------------------------------##
        
        ## starting r    
        if ((r > (T0.min-1) & force%in%c(0,2)) | (r > (T0.min-2) & force%in%c(1,3))) {
            cat("Warning: r is too big compared with T0; reset to 0.\n")
            r <- 0
        }
        
        ## initial values
        cat("Cross-validating ...","\r")
        
        ## store all MSPE
        if (force%in%c(0, 2)) {
            r.max <- max(min((T0.min-1), r.end), 0)
        } else {
            r.max <- max(min((T0.min-2), r.end), 0)
        }
        if (r.max == 0) {
            r.cv <- 0
            cat("Cross validation cannot be performed since available pre-treatment records of treated units are too few. So set r.cv = 0.\n ")
            if (!0%in%I.co) {
                ## if (force!=0) {
                    est.co.best <- inter_fe(Y.co, X.co, 0, force = force, beta0 = beta0, tol) 
                ## } else {
                ##     est.co.best<-inter_fe(Y.co, abind(I.co, X.co, along=3), 0, force=0, beta0 = beta0) 
                ## }  
            } else {
                ## if (force!=0) {
                    est.co.best <- inter_fe_ub(Y.co, X.co, I.co, 0, force = force, tol)
                ## } else {
                ##     est.co.best<-inter_fe_ub(Y.co, abind(I.co, X.co, along=3), I.co, 0, force=0, beta0 = beta0)
                ## }
            }

        } else {
            CV.out <- matrix(NA, (r.max - r + 1), 4)
            colnames(CV.out) <- c("r", "sigma2", "IC", "MSPE")
            CV.out[,"r"] <- c(r:r.max)
            CV.out[,"MSPE"] <- 1e20
        
            for (i in 1:dim(CV.out)[1]) { ## cross-validation loop starts 
  
                ## inter FE based on control, before & after 
                r <- CV.out[i, "r"]
                if (!0%in%I.co) {
                    ## if (force!=0) {
                        est.co <- inter_fe(Y = Y.co, X = X.co, r,
                                           force = force, beta0 = beta0, tol)
                    ## } else {
                    ##     est.co<-inter_fe(Y = Y.co, X = abind(I.co, X.co, along=3), 
                    ##                      r, force = 0, beta0 = beta0)                        
                    ## }
                } else {
                    ## if (force!=0) {
                        est.co <- inter_fe_ub(Y = Y.co, X = X.co, I = I.co, r,
                                              force = force, tol)
                    ## } else {
                    ##     est.co<-inter_fe_ub(Y = Y.co, X = abind(I.co, X.co, along=3), 
                    ##                         I = I.co, r, force = 0, beta0 = beta0)                        
                    ## }
                }       
   
                if (p > 0) {
                    na.pos <- is.nan(est.co$beta)
                    ## if (est.co$validX == 0) {
                    ##     beta <- matrix(0, p, 1) 
                    ## }
                    beta <- est.co$beta
                    beta[is.nan(est.co$beta)] <- 0 ## time invariant covar
                } 
            
                if (is.null(norm.para)) {
                    sigma2 <- est.co$sigma2
                    IC <- est.co$IC
                } else {
                    sigma2 <- est.co$sigma2 * (norm.para[1]^2)
                    IC <- est.co$IC - log(est.co$sigma2) + log(sigma2)
                }
               
                if (r!=0) {
                    ## factor T*r
                    F.hat <- as.matrix(est.co$factor)
                    ## the last column is for alpha_i
                    if (force%in%c(1, 3)) {F.hat <- cbind(F.hat, rep(1,TT))}  
                }      
                
                ## take out the effect of X (nothing to do with CV)
                U.tr <- Y.tr
               
                if (p > 0) {
                    for (j in 1:p) {
                        U.tr <- U.tr - X.tr[, , j] * beta[j]
                    }
                }
            
                ## take out grant mean and time fixed effects (nothing to do with CV)
                U.tr <- U.tr - matrix(est.co$mu, TT, Ntr) ## grand mean
                if (force%in%c(2, 3)) {
                    U.tr <- U.tr - matrix(est.co$xi, TT, Ntr, byrow = FALSE)
                } ## time fixed effects

                ## reset unbalanced data
                if (0%in%I.tr) {
                    U.tr[which(I.tr == 0)] <- 0
                }
            
                ## save for the moment       
                U.sav <- U.tr
            
                ## leave-one-out cross-validation
                sum.e2 <- num.y <- 0
                for (lv in unique(unlist(time.pre))) { ## leave one out using the pre-treatment period
   
                    U.tr <- U.sav
                    ## take out lv from U.tr.pre
                    if ( max(T0) == T0.min & (!0%in%I.tr) ) {
                        U.lv <- as.matrix(U.tr[setdiff(c(1:T0.min), lv), ]) ## setdiff : x
                    } else {
                        U.tr.pre.v <- as.vector(U.tr)[which(pre.v == 1)]    ## pre-treatment residual in a vector
                        U.tr.pre <- split(U.tr.pre.v, id.tr.pre.v) ##  a list of pretreatment residuals
                        if (!0%in%I.tr) {
                            U.lv<-lapply(U.tr.pre, function(vec){return(vec[-lv])}) ## a list    
                        } else {
                            ## U.tr.pre.sav <- U.tr.pre
                            for (i.tr in 1:Ntr) {
                                U.tmp <- U.tr.pre[[i.tr]]
                                U.tr.pre[[i.tr]] <- U.tmp[!time.pre[[i.tr]] == lv]
                            }                        
                            U.lv <- U.tr.pre
                            ## U.tr.pre <- U.tr.pre.sav    
                        }        
                    }

                    if (r == 0) {            
                        if (force%in%c(1,3)) { ## take out unit fixed effect
                            if ( max(T0) == T0.min & (!0%in%I.tr) ) {
                                alpha.tr.lv <- colMeans(U.lv)
                                U.tr <- U.tr - matrix(alpha.tr.lv, TT, Ntr, byrow=TRUE)
                            } else {
                                alpha.tr.lv <- sapply(U.lv,mean)
                                U.tr <- U.tr - matrix(alpha.tr.lv, TT, Ntr, byrow=TRUE)
                                ######
                                ## if (0%in%I.tr) {
                                ##     U.tr[which(I.tr==0)] <- 0
                                ## }
                                ###### not necessary
                            }
                        } 
                        e <- U.tr[which(time == lv),] ## that period
                    } else {  ## case: r>0
                        ## take out the effect of factors
                        F.lv <- as.matrix(F.hat[which(time != lv), ])
                        if ( max(T0) ==T0.min & (!0%in%I.tr) ) {
                            F.lv.pre <- F.hat[setdiff(c(1:T0.min), lv), ]
                            lambda.lv <- try(
                                solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% U.lv,
                                silent = TRUE)
                            if('try-error' %in% class(lambda.lv)) {
                                break
                            }
                            ## lamda.lv r*N matrix
                        } else {
                            if (!0%in%I.tr) {
                                lambda.lv <- try(as.matrix(sapply(U.lv, function(vec){
                                    F.lv.pre <- as.matrix(F.lv[1:length(vec),])
                                    l.lv.tr <- solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% vec 
                                    return(l.lv.tr) ## a vector of each individual lambdas
                                })), silent = TRUE)
                                if('try-error' %in% class(lambda.lv)) {
                                    break
                                } else {
                                    if( (r == 1) & (force%in%c(0,2)) ){
                                        lambda.lv <- t(lambda.lv)
                                    }    
                                }
                            } else { ## unbalanced data r*N
                                if (force%in%c(1,3)) {
                                    lambda.lv <- matrix(NA, (r+1), Ntr)
                                } else {
                                    lambda.lv <- matrix(NA, r, Ntr)
                                }
                                test <- try(
                                    for (i.tr in 1:Ntr) {
                                        ## F.lv.pre <- as.matrix(F.lv[unlist(time.pre[i.tr]),])
                                        F.lv.pre <- as.matrix(F.hat[setdiff(time.pre[[i.tr]], lv),])
                                        lambda.lv[,i.tr] <- solve(t(F.lv.pre) %*% F.lv.pre) %*% t(F.lv.pre) %*% as.matrix(U.lv[[i.tr]])
                                    }, silent = TRUE)
                                if('try-error' %in% class(test)) {
                                    break
                                }
                            } 
                        }
                        ##if (r>1|(r==1&force%in%c(1,3))) {
                            lambda.lv <- t(lambda.lv) ## N*r
                        ##}
                        ## error term (the left-out period)
                        e <- U.tr[which(time == lv),] - c(F.hat[which(time == lv),] %*% t(lambda.lv)) 
                    }
                    if (sameT0 == FALSE | 0%in%I.tr) { # those who are actually not treated
                        e <- e[which(pre[which(time == lv),] == TRUE)]    
                    }
                    ## sum up
                    sum.e2 <- sum.e2 + t(e) %*% e
                    num.y <- num.y + length(e)
                
                } ## end of leave-one-out
            
                MSPE <- ifelse(num.y == 0, Inf, sum.e2/num.y)
                if (!is.null(norm.para)) {
                    MSPE <- MSPE * (norm.para[1]^2)
                }

                if ((min(CV.out[,"MSPE"]) - MSPE) > tol * min(CV.out[,"MSPE"])) {
                    ## at least 5% improvement for MPSE
                    est.co.best <- est.co  ## interFE result with the best r
                    r.cv <- r
                } else {
                    if (r == r.cv + 1) cat("*")
                } 
                CV.out[i, 2:4] <- c(sigma2, IC, MSPE)
                cat("\n r = ",r, "; sigma2 = ",
                    sprintf("%.5f",sigma2), "; IC = ",
                    sprintf("%.5f",IC), "; MSPE = ",
                    sprintf("%.5f",MSPE), sep="")
            
            } ## end of while: search for r_star over
           
        
            if (r > (T0.min-1)) {cat(" (r hits maximum)")}
            cat("\n\n r* = ",r.cv, sep="")
            cat("\n\n") 
        
            MSPE.best <- min(CV.out[,"MSPE"])
        }
        
    } ## End of Cross-Validation

    validX <- est.co.best$validX
    
    ##-------------------------------##
    ## ATT and Counterfactuals 
    ##-------------------------------##
    
    ## variance of the error term
    if (is.null(norm.para)) {
        sigma2 <- est.co.best$sigma2   
        IC <- est.co.best$IC
    } else {
        sigma2 <- est.co.best$sigma2 * (norm.para[1]^2)
        IC <- est.co.best$IC - log(est.co.best$sigma2) + log(sigma2)       
    }
 
    ## ## take out the effect of X
    U.tr <- Y.tr
    res.co <- est.co.best$residuals


    if (p>0) {
        beta <- est.co.best$beta
        if (est.co.best$validX == 0) {
            beta <- matrix(0, p, 1) 
        } else {
            beta <- est.co.best$beta
            beta[is.nan(est.co.best$beta)] <- 0
        }
        ## if (!is.null(norm.para)) {
        ##     est.co.best$beta <- est.co.best$beta/norm.para[2]    
        ## }
        for (j in 1:p) {U.tr<-U.tr-X.tr[, , j] * beta[j]}
    } else {
        beta <- NA
    }
    
    mu <- est.co.best$mu 
    U.tr <- U.tr - matrix(mu, TT, Ntr) ## grand mean
    Y.fe.bar <- rep(mu, TT)
    
    if (force%in%c(2,3)) {
        xi <- est.co.best$xi ## a (TT*1) matrix
        U.tr <- U.tr - matrix(c(xi), TT, Ntr, byrow = FALSE) ## will be adjusted at last
        Y.fe.bar <- Y.fe.bar + xi
    }
    if ( max(T0) == T0.min & (!0%in%I.tr) ) {
        U.tr.pre <- as.matrix(U.tr[1:T0.min,])
    } else {
        ## not necessary to reset utr for ub data for pre.v doesn't include them
        U.tr.pre.v <- as.vector(U.tr)[which(pre.v == 1)] # pre-treatment residual in a vector
        U.tr.pre <- split(U.tr.pre.v, id.tr.pre.v) ##  a list of pretreatment residuals
    }
     
    
    ## the error structure
    if (r.cv == 0) {
        if (force%in%c(1,3)) { ## take out unit fixed effect
            if ((max(T0) == T0.min) & (!0%in%I.tr)) {
                alpha.tr <- colMeans(U.tr.pre)
                U.tr <- U.tr - matrix(alpha.tr, TT, Ntr, byrow = TRUE)
            } else {
                alpha.tr <- sapply(U.tr.pre, mean)
                U.tr <- U.tr - matrix(alpha.tr, TT, Ntr, byrow = TRUE)
            }
        }     
        eff <- U.tr  ## and that's it!

    } else { ## r.cv>0
        ## Factors
        F.hat <- as.matrix(est.co.best$factor)
        if (force %in% c(1, 3)) {F.hat <- cbind(F.hat, rep(1,TT))}
                                    # the last column is for alpha_i
         
        ## Lambda_tr (Ntr*r) or (Ntr*(r+1))
        if ( max(T0) == T0.min & (!0%in%I.tr)) {
            F.hat.pre <- F.hat[1:T0.min,]
            lambda.tr <- try(solve(t(F.hat.pre) %*% F.hat.pre) %*% t(F.hat.pre) %*% U.tr.pre,
                           silent = TRUE)
            if('try-error' %in% class(lambda.tr)) {
                return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1)))
                ## stop("Error occurs. Please set a smaller value of factor number.")
            }
        } else {
            if (!0%in%I.tr) {
                lambda.tr <- try(as.matrix(sapply(U.tr.pre, function(vec) {
                    F.hat.pre <- as.matrix(F.hat[1:length(vec),])
                    l.tr <- solve(t(F.hat.pre)%*%F.hat.pre)%*%t(F.hat.pre)%*%vec
                    return(l.tr) ## a vector of each individual lambdas
                })), silent =TRUE)
                if('try-error' %in% class(lambda.tr)) {
                    return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1)))
                    ## stop("Error occurs. Please set a smaller value of factor number.")
                }
                if ( (r.cv == 1) & (force%in%c(0,2)) ) {
                    lambda.tr <- t(lambda.tr)    
                }
            } else {
                if (force%in%c(1,3)) {
                    lambda.tr <- matrix(NA, (r.cv+1), Ntr)
                } else {
                    lambda.tr <- matrix(NA, r.cv, Ntr)
                }
                test <- try(
                    for (i.tr in 1:Ntr) {
                        F.hat.pre <- as.matrix(F.hat[time.pre[[i.tr]],])
                        lambda.tr[,i.tr] <- solve(t(F.hat.pre)%*%F.hat.pre)%*%t(F.hat.pre)%*%as.matrix(U.tr.pre[[i.tr]])
                    }, silent =TRUE
                )
                if('try-error' %in% class(lambda.tr)) {
                    return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1)))
                    ## stop("Error occurs. Please set a smaller value of factor number.")
                }
            }
        }
        ## if ((r.cv>1) | (r.cv==1 & force%in%c(1,3))) {
            lambda.tr <- t(lambda.tr) ## Ntr * r
        ## }

        ## predicting the treatment effect
        eff <- U.tr - F.hat%*%t(lambda.tr) 

        ## for storage
        if (force%in%c(1,3)) {
            alpha.tr <- lambda.tr[, (r.cv+1), drop = FALSE]
            lambda.tr <- lambda.tr[, 1:r.cv, drop = FALSE] 
        }

        if (boot == 0) {
            inv.tr <- try(
                ginv(t(as.matrix(lambda.tr))), silent = TRUE
            )
            if (!'try-error' %in% class(inv.tr)) {
                wgt.implied <- t(inv.tr%*%t(as.matrix(est.co.best$lambda)))
            }
        }

    } ## end of r!=0 case

    if (0%in%I.tr) {
        eff[which(I.tr == 0)] <- 0 ## adjust    
    } ## missing data will be adjusted to NA finally
   
    ##-------------------------------##
    ## Summarize
    ##-------------------------------##  
    
    ## counterfactuals and averages
    Y.ct <- as.matrix(Y.tr - eff)
    
    if (is.null(W)) {
        if (!0%in%I) {
            Y.tr.bar <- rowMeans(Y.tr)
            Y.ct.bar <- rowMeans(Y.ct)
            Y.co.bar <- rowMeans(Y.co)  
        } else {
            Y.tr.bar <- rowSums(Y.tr)/rowSums(I.tr)
            Y.ct.bar <- rowSums(Y.ct)/rowSums(I.tr)
            Y.co.bar <- rowSums(Y.co)/rowSums(I.co)
        }
    } else {
        Y.tr.bar <- rowSums(Y.tr * W.tr)/rowSums(W.tr)
        Y.ct.bar <- rowSums(Y.ct * W.tr)/rowSums(W.tr)
        Y.co.bar <- rowSums(Y.co)/rowSums(I.co)
    }

    ##Y.tr and Y.ct
    Y.bar <- cbind(Y.tr.bar, Y.ct.bar, Y.co.bar)
    colnames(Y.bar) <- c("Y.tr.bar", "Y.ct.bar", "Y.co.bar")
    
    ## ATT and average outcomes
    if (DID == TRUE) { ## diff-in-diffs: same timing
        if (is.null(W)) {
            if (!0%in%I.tr) {
                att <- rowMeans(eff)
            } else {
                att <- rowSums(eff)/rowSums(I.tr)
            }
        } else {
            att <- rowSums(eff * W.tr)/rowSums(W.tr)
        }
    } else { ## diff timing, centered the att
        if (!0%in%I.tr) {
            eff.cnt <- Y.tr.center <- matrix(NA, TT, Ntr)
            if (!is.null(W)) {
                W.tr.center <- matrix(NA, TT, Ntr)
            }
            for (j in 1:Ntr) {
                eff.cnt[1:(TT+T0.min-T0[j]), j] <- eff[(T0[j]-T0.min+1):TT, j]  
                Y.tr.center[1:(TT+T0.min-T0[j]), j] <- Y.tr[(T0[j]-T0.min+1):TT, j]
                if (!is.null(W)) {
                    W.tr.center[1:(TT+T0.min-T0[j]), j] <- W.tr[(T0[j]-T0.min+1):TT, j]   
                }
            }
            if (is.null(W)) {
                att <- apply(eff.cnt, 1, mean, na.rm = TRUE)
                Y.tr.cnt <- apply(Y.tr.center, 1, mean, na.rm = TRUE)
                Y.ct.cnt <- Y.tr.cnt - att
            } else {
                eff.cnt[which(is.na(eff.cnt))] <- 0
                Y.tr.center[which(is.na(Y.tr.center))] <- 0
                W.tr.center[which(is.na(W.tr.center))] <- 0
                att <- rowSums(eff.cnt * W.tr.center)/rowSums(W.tr.center)
                Y.tr.center <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
                Y.ct.cnt <- Y.tr.cnt - att
            }
        } else {
            T0.ub <- apply(as.matrix(D[,which(tr == 1)] == 0), 2, sum) 
            T0.ub.min <- min(T0.ub)
            eff.cnt <- Y.tr.center <- matrix(NA, TT, Ntr)
            eff[which(I.tr == 0)] <- NA
            Y.tr[which(I.tr == 0)] <- NA
            if (!is.null(W)) {
                W.tr.center <- matrix(NA, TT, Ntr)
                W.tr[which(I.tr == 0)] <- NA
            }
            for (j in 1:Ntr) {
                eff.cnt[1:(TT+T0.ub.min-T0.ub[j]), j] <- eff[(T0.ub[j]-T0.ub.min+1):TT, j]  
                Y.tr.center[1:(TT+T0.ub.min-T0.ub[j]),j] <- Y.tr[(T0.ub[j]-T0.ub.min+1):TT, j]
                if (!is.null(W)) {
                    W.tr.center[1:(TT+T0.ub.min-T0.ub[j]), j] <- W.tr[(T0.ub[j]-T0.ub.min+1):TT, j]   
                }
            }
            if (is.null(W)) {
                att <- apply(eff.cnt, 1, mean, na.rm = TRUE)
                Y.tr.cnt <- apply(Y.tr.center, 1, mean, na.rm = TRUE)
                Y.ct.cnt <- Y.tr.cnt - att
            } else {
                eff.cnt[which(is.na(eff.cnt))] <- 0
                Y.tr.center[which(is.na(Y.tr.center))] <- 0
                W.tr.center[which(is.na(W.tr.center))] <- 0
                att <- rowSums(eff.cnt * W.tr.center)/rowSums(W.tr.center)
                Y.tr.center <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
                Y.ct.cnt <- Y.tr.cnt - att
            }
        }
    }
    eff[which(is.na(eff))] <- 0 ## to calulate att
    if (is.null(W)) {
        att.avg <- sum(eff * post)/sum(post)
    } else {
        att.avg <- sum(eff * post * W.tr)/sum(post * W.tr)
    }

    ## AR1: calculate accumulative effect
    if (AR1 == TRUE) {
        rho <- est.co.best$beta[1]
        if (length(beta) > 1) {
            beta <- beta[-1]
        } 
        eff.tmp <- eff * D[, id.tr]
        eff.acc <- matrix(0, TT, Ntr)
        for (t in (T0.min + 1):TT) {
            for (i in 0:(t-T0.min-1)) {
                eff.acc[t,] <- eff.acc[t,] + eff.tmp[t-i,] * (rho^i)
            }
        }      
    } 
    
    ## final adjust unbalanced output
    if (0%in%I) {
        eff[which(I.tr == 0)] <- NA
        Y.ct[which(I.tr == 0)] <- NA
        Y.tr[which(I.tr == 0)] <- NA
        res.co[which(I.co == 0)] <- NA
        Y.co[which(I.co == 0)] <- NA
    }
    ## adjust beta: invariant covar
    if (p > 0) {
        if( sum(na.pos) > 0 ) {
            beta[na.pos] <- NA
        }
    }

    ## final adjustment
    if (!is.null(norm.para)) {
        mu <- mu * norm.para[1]
        ## if (p>0) {
        ##     beta <- beta*norm.para[1]/norm.para[2:length(norm.para)]
        ## }
        if (r.cv > 0) {
            est.co.best$lambda <- est.co.best$lambda * norm.para[1]
            lambda.tr <- lambda.tr * norm.para[1]
        }
        if (force%in%c(1, 3)) {
            est.co.best$alpha <- est.co.best$alpha * norm.para[1]
            alpha.tr <- alpha.tr * norm.para[1]
        }
        if (force%in%c(2,3)) {
            xi <- xi * norm.para[1]
        }
        res.co <- res.co * norm.para[1] 
        Y.tr <- Y.tr * norm.para[1] 
        Y.ct <- Y.ct * norm.para[1]
        Y.co <- Y.co * norm.para[1]
        eff <- eff * norm.para[1]
        Y.bar <- Y.bar * norm.para[1]
        att <- att * norm.para[1]
        att.avg <- att.avg * norm.para[1]
        if ( sameT0 == FALSE | 0%in%I.tr ) {
            eff.cnt <- eff.cnt * norm.para[1]
            Y.tr.cnt <- Y.tr.cnt * norm.para[1]
            Y.ct.cnt <- Y.ct.cnt * norm.para[1]
        }
    }

    T0<-apply(as.matrix(D[,which(tr == 1)] == 0), 2, sum) ## for plot

    
    ##-------------------------------##
    ## Storage 
    ##-------------------------------##  

    ##control group residuals
    out<-list(
        ## main results
        D.tr = D.tr,
        I.tr = I.tr,
        Y.tr = Y.tr,
        Y.ct = Y.ct,
        Y.co = Y.co, 
        eff = eff,
        Y.bar = Y.bar,
        att = att,
        att.avg = att.avg,
        ## supporting
        force = force,
        DID = DID,
        T = TT,
        N = N,
        p = p,
        Ntr = Ntr,
        Nco = Nco,
        T0 = T0,
        tr = tr,
        pre = pre,
        post = post,
        r.cv = r.cv, 
        IC = IC,
        beta = beta,
        est.co = est.co.best,
        mu = mu,
        validX = validX
    )

    out <- c(out,list(sigma2 = sigma2, res.co=res.co))
    

    if ( DID == FALSE ) {
        out<-c(out,list(eff.cnt = eff.cnt,
                        Y.tr.cnt = Y.tr.cnt,
                        Y.ct.cnt = Y.ct.cnt))
    }
    if (CV == 1 & r.max != 0) {
        out<-c(out, list(MSPE = MSPE.best,
                         CV.out = CV.out))
    }
    if (r.cv>0) {
        out<-c(out,list(
                       niter = est.co.best$niter,
                       factor = as.matrix(est.co.best$factor),
                       lambda.co = as.matrix(est.co.best$lambda),
                       lambda.tr = as.matrix(lambda.tr) ## Ntr*r
                   ))
        if (boot == 0) {
            if (!'try-error' %in% class(inv.tr)) {
                out <- c(out, list(wgt.implied = wgt.implied))
            }
        } 
    } 
    if (force==1) {
        out<-c(out, list(
                         alpha.tr = alpha.tr,
                         alpha.co = est.co.best$alpha))
    } else if (force == 2) {
        out<-c(out,list(xi = xi))
    } else if (force == 3) {
        out<-c(out,list(
                        alpha.tr = alpha.tr,
                        alpha.co = est.co.best$alpha,
                        xi = xi))
    }
    if (AR1 == TRUE) {
        out<-c(out,list(rho = rho,
                        eff.acc = eff.acc))
    }

    return(out)
} ## Core functions ends

###################################################################
## EM Algorithm
###################################################################

synth.em<-function(Y, # Outcome variable, (T*N) matrix
                   X, # Explanatory variables:  (T*N*p) array
                   D, # indicator for treated unit (tr==1) 
                   I,
                   W = NULL,
                   r = 0, # number of factors
                   force, # specifying fixed effects
                   tol = 1e-5,
                   AR1 = 0,
                   norm.para,
                   boot = 0
                   ){

    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  
    
    ## unit id and time
    TT <-dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(X) == FALSE) {p <- dim(X)[3]} else {p <- 0}
     
    ## treatement indicator
    tr <- D[TT,] == 1  ## cross-sectional: treated unit
    co <- D[TT,] == 0
    I.tr <- as.matrix(I[, tr])
    I.co <- I[, co]

    if (!0 %in% I.tr) {
        ## a (TT*Ntr) matrix, time dimension: before treatment
        pre <- as.matrix(D[,which(tr == 1)] == 0)
        post <- as.matrix(D[,which(tr == 1)] == 1)   
    } else {
        pre <- as.matrix(D[,which(tr == 1)] == 0 & I[,which(tr==1)] == 1)
        post <- as.matrix(D[,which(tr==1)] == 1 & I[,which(tr==1)] == 1)
    }
    
    Ntr <- sum(tr)
    Nco <- N - Ntr
    ## careful: only valid for balanced panel
    T0 <- apply(pre, 2, sum) 
    T0.min <- min(T0)
    sameT0 <- length(unique(T0)) == 1 ## treatment kicks in at the same time 

    D.tr <- D[,which(tr == 1)]
    T0.ub <- apply(as.matrix(D[,which(tr == 1)] == 0), 2, sum) 
    T0.ub.min <- min(T0.ub) ## unbalanced data

    if (!0%in%I.tr) {
        DID <- sameT0
    } else {
        DID <- length(unique(T0.ub)) == 1
    }

    if (!is.null(W)) {
        W.tr <- as.matrix(W[,which(tr == 1)])
    }

    
    id <- 1:N
    time <- 1:TT
    id.tr <- which(tr == 1) ## treated id
    id.co <- which(tr == 0)
    
    pre.v <- as.vector(pre)  ## vectorized "pre-treatment" indicator
    id.tr.pre.v <- rep(id, each = TT)[which(pre.v == 1)]  ## vectorized pre-treatment grouping variable for the treated
    time.pre <- split(rep(time, Ntr)[which(pre.v == 1)], id.tr.pre.v) ## a list of pre-treatment periods

    ## parsing data
    Y.tr <- as.matrix(Y[,tr])
    Y.co <- Y[,!tr]
    
    
    ##-------------------------------##
    ## Main Algorithm
    ##-------------------------------##

    init<-synth.core(Y = Y, X = X, D = D, I = I, W = W,
                     r = r, force = force,
                     CV = 0, tol = tol, AR1 = AR1, 
                     norm.para = NULL, boot = boot)
    
    ## throw out error: may occur during bootstrap
    if(length(init) == 2 || length(init) == 3) {
        return(init)
    }
    
    eff0 <- init$eff
    eff0[is.na(eff0)] <- 0
    Y.ct <- init$Y.ct
    Y.ct[is.na(Y.ct)] <- 0

    if (p > 0) {
        beta0 <- init$beta
        if (NA%in%beta0) {
            beta0 <- as.matrix(beta0[which(!is.na(beta0))])
        }
    } else {
        beta0 <- matrix(0, 0, 1)
    }
    
    diff <- 100
    trace.diff <- c()
    niter <- 0

    while (niter <= 500 & diff > tol) {

        ## E step
        Y.e <- Y  # T*N
        Y.e.tr.tmp <- Y.e[,id.tr]
        Y.e.tr.tmp[which(post == 1)] <- Y.ct[which(post == 1)]
        Y.e[,id.tr] <- Y.e.tr.tmp 

        ## M step
        if (!0%in%I) {
            ## if (force!=0) {
                est<-inter_fe(Y.e, X, r, force=force, beta0 = beta0, tol)
            ## } else {
            ##     est<-inter_fe(Y.e, abind(I,X,along=3), r, force=0, beta0 = beta0)
            ## }
        } else {
            ## if (force!=0) {
                est<-inter_fe_ub(Y.e, X, I, r, force=force, tol)
            ## } else {
            ##     est<-inter_fe_ub(Y.e, abind(I,X,along=3), I, r, force=0, beta0 = beta0)
            ## }
        }
        Y.ct <- as.matrix(Y.e[,id.tr] - est$residuals[,id.tr]) # T * Ntr

        eff <- as.matrix(Y.tr - Y.ct)  # T * Ntr
        diff <- norm(eff0-eff, type="F")

        eff0 <- eff

        trace.diff <- c(trace.diff,diff)
        niter <- niter + 1  
    }
    

    ## variance of the error term
    if (is.null(norm.para)) {
        sigma2<-est$sigma2   
        IC<-est$IC
    } else {
        sigma2<-est$sigma2*(norm.para[1]^2)
        IC <- est$IC-log(est$sigma2) + log(sigma2)       
    }


    ##-------------------------------##
    ## Summarize
    ##-------------------------------##  
    
    ## counterfactuals and averages
    if (is.null(W)) {
        if (!0%in%I) {
            Y.tr.bar <- rowMeans(Y.tr)
            Y.ct.bar <- rowMeans(Y.ct)
            Y.co.bar <- rowMeans(Y.co)  
        } else {
            Y.tr.bar <- rowSums(Y.tr)/rowSums(I.tr)
            Y.ct.bar <- rowSums(Y.ct)/rowSums(I.tr)
            Y.co.bar <- rowSums(Y.co)/rowSums(I.co)
        }
    } else {
        Y.tr.bar <- rowSums(Y.tr * W.tr)/rowSums(W.tr)
        Y.ct.bar <- rowSums(Y.ct * W.tr)/rowSums(W.tr)
        Y.co.bar <- rowSums(Y.co)/rowSums(I.co)
    }

    ##Y.tr and Y.ct
    Y.bar <- cbind(Y.tr.bar,Y.ct.bar,Y.co.bar)
    colnames(Y.bar) <- c("Y.tr.bar","Y.ct.bar","Y.co.bar")

    ## ATT and average outcomes
    if (DID == TRUE) { ## diff-in-diffs: same timing
        if (is.null(W)) {
            if (!0%in%I.tr) {
                att <- rowMeans(eff)
            } else {
                att <- rowSums(eff)/rowSums(I.tr)
            }
        } else {
            att <- rowSums(eff * W.tr)/rowSums(W.tr)
        }
    } else { ## diff timing, centered the att
        if (!0%in%I.tr) {
            eff.cnt <- Y.tr.center <- matrix(NA, TT, Ntr)
            if (!is.null(W)) {
                W.tr.center <- matrix(NA, TT, Ntr)
            }
            for (j in 1:Ntr) {
                eff.cnt[1:(TT+T0.min-T0[j]), j] <- eff[(T0[j]-T0.min+1):TT, j]  
                Y.tr.center[1:(TT+T0.min-T0[j]), j] <- Y.tr[(T0[j]-T0.min+1):TT, j]
                if (!is.null(W)) {
                    W.tr.center[1:(TT+T0.min-T0[j]), j] <- W.tr[(T0[j]-T0.min+1):TT, j]   
                }
            }
            if (is.null(W)) {
                att <- apply(eff.cnt, 1, mean, na.rm = TRUE)
                Y.tr.cnt <- apply(Y.tr.center, 1, mean, na.rm = TRUE)
                Y.ct.cnt <- Y.tr.cnt - att
            } else {
                eff.cnt[which(is.na(eff.cnt))] <- 0
                Y.tr.center[which(is.na(Y.tr.center))] <- 0
                W.tr.center[which(is.na(W.tr.center))] <- 0
                att <- rowSums(eff.cnt * W.tr.center)/rowSums(W.tr.center)
                Y.tr.center <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
                Y.ct.cnt <- Y.tr.cnt - att
            }
        } else {
            T0.ub <- apply(as.matrix(D[,which(tr == 1)] == 0), 2, sum) 
            T0.ub.min <- min(T0.ub)
            eff.cnt <- Y.tr.center <- matrix(NA, TT, Ntr)
            eff[which(I.tr == 0)] <- NA
            Y.tr[which(I.tr == 0)] <- NA
            if (!is.null(W)) {
                W.tr.center <- matrix(NA, TT, Ntr)
                W.tr[which(I.tr == 0)] <- NA
            }
            for (j in 1:Ntr) {
                eff.cnt[1:(TT+T0.ub.min-T0.ub[j]), j] <- eff[(T0.ub[j]-T0.ub.min+1):TT, j]  
                Y.tr.center[1:(TT+T0.ub.min-T0.ub[j]),j] <- Y.tr[(T0.ub[j]-T0.ub.min+1):TT, j]
                if (!is.null(W)) {
                    W.tr.center[1:(TT+T0.ub.min-T0.ub[j]), j] <- W.tr[(T0.ub[j]-T0.ub.min+1):TT, j]   
                }
            }
            if (is.null(W)) {
                att <- apply(eff.cnt, 1, mean, na.rm = TRUE)
                Y.tr.cnt <- apply(Y.tr.center, 1, mean, na.rm = TRUE)
                Y.ct.cnt <- Y.tr.cnt - att
            } else {
                eff.cnt[which(is.na(eff.cnt))] <- 0
                Y.tr.center[which(is.na(Y.tr.center))] <- 0
                W.tr.center[which(is.na(W.tr.center))] <- 0
                att <- rowSums(eff.cnt * W.tr.center)/rowSums(W.tr.center)
                Y.tr.center <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
                Y.ct.cnt <- Y.tr.cnt - att
            }
        }
    }
    eff[which(is.na(eff))] <- 0 ## to calulate att
    if (is.null(W)) {
        att.avg <- sum(eff * post)/sum(post)
    } else {
        att.avg <- sum(eff * post * W.tr)/sum(post * W.tr)
    }

    ## fixed effects
    mu<-est$mu
    
    if (force%in%c(1,3)) {
        alpha <- est$alpha
        alpha.tr <- alpha[id.tr]
        alpha.co <- alpha[id.co]
    }
    if (force%in%c(2,3)) {
        xi<-est$xi 
    }
    ## factors
    if (r > 0) {
        lambda <- est$lambda
        lambda.tr <- lambda[id.tr,,drop=FALSE]
        lambda.co <- lambda[id.co,]
        if (boot == 0) {
            inv.tr <- try(
                ginv(t(as.matrix(lambda.tr))), silent = TRUE
            )
            if (!'try-error' %in% class(inv.tr)) {
                wgt.implied <- t(inv.tr%*%t(as.matrix(est.co.best$lambda)))
            }
        }
    }
    ## AR1: calculate accumulative effect
    if (AR1 == TRUE) {
        rho<-est$beta[1]
        if (length(beta)>1) {
            beta<-beta[-1]
        }  
        eff.tmp<-eff*D[,id.tr]
        eff.acc<-matrix(0,TT,Ntr)
        for (t in (T0.min+1):TT) {
            for (i in 0:(t-T0.min-1)) {
                eff.acc[t,]<-eff.acc[t,]+eff.tmp[t-i,]*(rho^i)
            }
        }      
    } 

    if (p > 0) {
        beta <- est$beta
        beta[is.nan(beta)] <- NA
    } else {
        beta <- NA
    }
    
    res.co <- est$residuals[,id.co]
    
    
    if (0%in%I) {
        eff[which(I.tr==0)] <- NA 
        Y.ct[which(I.tr==0)] <- NA
        Y.tr[which(I.tr==0)] <- NA
        res.co[which(I.co==0)] <- NA
        Y.co[which(I.co==0)] <- NA
    }

        ## final adjustment
    if (!is.null(norm.para)) {
        mu <- mu*norm.para[1]
        ## if (p>0) {
        ##     beta <- beta*norm.para[1]/norm.para[2:length(norm.para)]
        ## }
        if (r>0) {
            lambda.tr <- lambda.tr*norm.para[1]
            lambda.co <- lambda.co*norm.para[1]
        }
        if (force%in%c(1,3)) {
            alpha.tr <- alpha.tr*norm.para[1]
            alpha.co <- alpha.co*norm.para[1]
        }
        if (force%in%c(2,3)) {
            xi <- xi*norm.para[1]
        }
        res.co <- res.co*norm.para[1] 
        Y.tr <- Y.tr*norm.para[1] 
        Y.ct <- Y.ct*norm.para[1]
        Y.co <- Y.co*norm.para[1]
        eff <- eff*norm.para[1]
        Y.bar <- Y.bar*norm.para[1]
        att <- att*norm.para[1]
        att.avg <- att.avg*norm.para[1]
        if ( sameT0 == FALSE | 0%in%I.tr ) {
            eff.cnt <- eff.cnt*norm.para[1]
            Y.tr.cnt <- Y.tr.cnt*norm.para[1]
            Y.ct.cnt <- Y.ct.cnt*norm.para[1]
        }
    }

    T0<-apply(as.matrix(D[,which(tr==1)]==0),2,sum) 

    
    ##-------------------------------##
    ## Storage 
    ##-------------------------------##  
    
    out<-list(
        ## main results
        D.tr = D.tr,
        I.tr = I.tr,
        Y.tr=Y.tr,
        Y.ct=Y.ct,
        Y.co=Y.co, 
        eff=eff,
        Y.bar = Y.bar,
        att=att,
        att.avg=att.avg,
        ## supporting
        force=force,
        DID=DID,
        T=TT,
        N=N,
        p=p,
        Ntr=Ntr,
        Nco=Nco,
        T0=T0,
        tr=tr,
        pre=pre,
        post = post,
        r.cv=r,
       ## res.co=res.co,  ##control group residuals 
        beta = beta,
        niter = niter,
        IC = IC,
        mu = mu,
        validX = est$validX
    )

    out <- c(out,list(sigma2 = sigma2, res.co = res.co))
    
    
    if (DID==FALSE) {
        out<-c(out,list(eff.cnt=eff.cnt, ##
                        Y.tr.cnt=Y.tr.cnt, ##
                        Y.ct.cnt=Y.ct.cnt)) ##
    }
    
    if (r > 0) {
        out<-c(out,list(factor=as.matrix(est$factor),
                        lambda.co=as.matrix(lambda.co),
                        lambda.tr=as.matrix(lambda.tr)
                        )) 
        if (boot == 0) {
            if (!'try-error' %in% class(inv.tr)) {
                out <- c(out, list(wgt.implied = wgt.implied))
            }
        }
    } 
    if (force == 1) {
        out<-c(out,list(
                        alpha.tr=alpha.tr,
                        alpha.co=alpha.co))
    } else if (force == 2) {
        out<-c(out,list(xi = xi))
    } else if (force == 3) {
        out<-c(out,list(
                        alpha.tr = alpha.tr,
                        alpha.co = alpha.co,
                        xi = xi))
    }
    if (AR1 == TRUE) {
        out<-c(out,list(rho = rho,
                        eff.acc = eff.acc))
    }
    return(out)
    
    
} ## EM function ends


###################################################################
## EM Algorithm
###################################################################

synth.em.cv<-function(Y,  # Outcome variable, (T*N) matrix
                      X, # Explanatory variables:  (T*N*p) array
                      D, # indicator for treated unit (tr==1) 
                      I,
                      W = NULL,
                      r = 0, # number of factors: starting point
                      r.end = 5, # end point
                      force, # specifying fixed effects
                      tol=1e-5,
                      AR1 = 0,
                      norm.para){
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  
    
    
    ## unit id and time
    TT <-dim(Y)[1]
    N<-dim(Y)[2]
    if (is.null(X)==FALSE) {p<-dim(X)[3]} else {p<-0}
     
    ## treatement indicator
    tr<-D[TT,]==1  ## cross-sectional: treated unit
    co<-D[TT,]==0
    I.tr<-as.matrix(I[,tr])
    I.co<-I[,co]
    D.tr<-D[,tr]

    if (!0%in%I.tr) {
        ## a (TT*Ntr) matrix, time dimension: before treatment
        pre <- as.matrix(D[,which(tr==1)]==0)   
    } else {
        pre <- as.matrix(D[,which(tr==1)]==0&I[,which(tr==1)]==1)
    }
    
    Ntr<-sum(tr)
    Nco<-N-Ntr
    ## careful: only valid for balanced panel
    T0<-apply(pre,2,sum) 
    T0.min<-min(T0)
    sameT0<-length(unique(T0))==1 ## treatment kicks in at the same time 
    
    id<-1:N
    time<-1:TT
    id.tr<-which(tr==1) ## treated id
    id.co<-which(tr==0)
    
    pre.v<-as.vector(pre)  ## vectorized "pre-treatment" indicator
    id.tr.pre.v<-rep(id,each=TT)[which(pre.v==1)]  ## vectorized pre-treatment grouping variable for the treated
    time.pre<-split(rep(time,Ntr)[which(pre.v==1)],id.tr.pre.v) ## a list of pre-treatment periods

    ## parsing data
    Y.tr<-as.matrix(Y[,id.tr])
    Y.co<-Y[,id.co]

    ##-------------------------------##
    ## Cross-validation of r
    ##-------------------------------##
    
    ## starting r    
    if (r > (T0.min-1)) {
        cat("Warning: r is too big compared with T0; reset to 0.\n")
        r <- 0
    }

    ## store all MSPE 
    if (force%in%c(0,2)) {
        r.max<-max(min((T0.min-1),r.end),0)
    } else {
        r.max<-max(min((T0.min-2),r.end),0)
    }
    if (r.max==0) {
        stop("Cross validation cannot be performed since available pre-treatment records of treated units are too few. r.cv = 0.\n ")
    } else {
        CV.out<-matrix(NA,(r.max-r+1),4)
        colnames(CV.out)<-c("r","sigma2","IC","MSPE")
        CV.out[,"r"]<-c(r:r.max)
        CV.out[,"MSPE"]<-1e20
        cat("Cross-validating ...","\r")
        for (i in 1:dim(CV.out)[1]) { ## cross-validation loop starts
      
            r <- CV.out[i,"r"]
            est<-synth.em(Y = Y,X = X, D = D, I = I, W = W, r = r, force = force,
                          tol = tol, AR1 = AR1, norm.para = norm.para, boot = 0)
            sigma2<-est$sigma2
            IC<-est$IC
        
            ## leave-one-out cross-validation
            sum.e2<-num.y<-0
            for (lv in unique(unlist(time.pre))){ ## leave one out using the pre-treatment period

                D.cv <- D
                D.cv[which(time == lv), id.tr] <- 1 # set the left-out period to treated
                ###########################
                ## if (0%in%I.tr) {
                ##     D.cv[which(I==0)] <- 0
                ## } not necessary! synth.em can detect pre and post period for ub data
                out <- synth.em(Y = Y, X = X, D = D.cv, I = I, W = W, r = r, force = force,
                                tol = tol, AR1 = AR1, norm.para = norm.para, boot = 0)

                e <- out$eff[which(time == lv),]
              
                if (sameT0 == FALSE|0%in%I.tr) { # those who are actually not treated
                    e<-e[which(pre[which(time==lv),]==TRUE)]    
                }
                ## sum up
                sum.e2<-sum.e2+t(e)%*%e
                num.y<-num.y + length(e) 
            } ## end of leave-one-out
        
            MSPE<-sum.e2/num.y
            if ((min(CV.out[,"MSPE"]) - MSPE) > tol*min(CV.out[,"MSPE"])) {
                ## at least 5% improvement for MPSE
                est.best<-est  ## interFE result with the best r
                r.cv<-r
            } else {
                if (r==r.cv+1) cat("*")
            }
            CV.out[i,2:4]<-c(sigma2,IC,MSPE)
            cat("\n r = ",r,"; sigma2 = ",
                sprintf("%.5f",sigma2),"; IC = ",
                sprintf("%.5f",IC),"; MSPE = ",
                sprintf("%.5f",MSPE),sep="") 
        } ## end of while: search for r_star over
     
        if (r>(T0.min-1)) {
            cat(" (r hits maximum)")
        }
        cat("\n\n r* = ", r.cv, sep="") 
        cat("\n\n") 
        MSPE.best <- min(CV.out[,"MSPE"])
    }
    

    ##-------------------------------##
    ## Storage 
    ##-------------------------------##  
    
    out<-c(est.best, list(MSPE = MSPE.best, CV.out = CV.out)) 
    return(out) 
    
} ## EM cross validation ends


###################################################################
## Matrix Completion Function
###################################################################
synth.mc<-function(Y, # Outcome variable, (T*N) matrix
                   X, # Explanatory variables:  (T*N*p) array
                   D, #  Indicator for treated unit (tr==1) 
                   I,
                   W = NULL,
                   lambda = NULL,
                   nlambda = 10,
                   force,
                   CV = 1,
                   hasF = 1,
                   tol, # tolerance level
                   AR1 = 0,
                   norm.para = NULL){  
    
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  
    na.pos <- NULL
    ## CV <- is.null(lambda)
    ## unit id and time
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(X) == FALSE) {p <- dim(X)[3]} else {p <- 0}
     
    ## treatement indicator
    tr <- D[TT,] == 1  ## cross-sectional: treated unit
    co <- D[TT,] == 0
    I.tr <- as.matrix(I[,tr]) ## maybe only 1 treated unit
    I.co <- I[,co]

    ## replicate data
    YY <- Y
    II <- I
    ## treat post-treatment period treated units as missing
    YY[which(D==1)] <- 0
    II[which(D==1)] <- 0

    if (!0%in%I.tr) {
        ## a (TT*Ntr) matrix, time dimension: before treatment
        pre <- as.matrix(D[,which(tr == 1)] == 0)
        post <- as.matrix(D[,which(tr == 1)] == 1)   
    } else {
        pre <- as.matrix(D[,which(tr == 1)] == 0 & I[,which(tr == 1)] == 1)
        post <- as.matrix(D[,which(tr == 1)] == 1 & I[,which(tr == 1)] == 1)
    }

    D.tr <- as.matrix(D[,which(tr == 1)])
    T0.ub <- apply(D.tr == 0, 2, sum) 
    T0.ub.min <- min(T0.ub) ## unbalanced data

    if (!is.null(W)) {
        W.tr <- as.matrix(W[,which(tr == 1)])
    }

    Ntr <- sum(tr)
    Nco <- N - Ntr
    ## careful: only valid for balanced panel
    T0 <- apply(pre, 2, sum) 
    T0.min <- min(T0)
    sameT0 <- length(unique(T0)) == 1 ## treatment kicks in at the same time 
                                      ## unbalanced case needs more conditions
    if (!0%in%I.tr) {
        DID <- sameT0
    } else {
        DID <- length(unique(T0.ub)) == 1
    }
    
    id <- 1:N
    time <- 1:TT
    id.tr <- which(tr == 1) ## treated id
    id.co <- which(tr == 0)

    ## parsing data
    Y.tr <- as.matrix(Y[,id.tr])
    Y.co <- as.matrix(Y[,id.co])
    
    ##-------------------------------##
    ## Main Algorithm
    ##-------------------------------##

    validX <- 1 ## no multi-colinearity
    
    if (CV == FALSE) { ## case: CV==0 or no factor  
        ## matrix completion
        est.best <- inter_fe_mc(YY, X, II, hasF, lambda[1], force, tol) 

        if (p > 0) {
            na.pos <- is.nan(est.best$beta)
        } 

        lambda.cv <- lambda[1]    
    } else { 
        
        ##-------------------------------##
        ## Cross-validation of lambda
        ##-------------------------------##
        
        ## initial values
        cat("Cross-validating ...","\r")

        tot.id <- which(c(II)==1) ## observed control data
        cv.count <- ceiling((sum(II)*sum(II))/(N*TT))

        if (is.null(lambda) || length(lambda) == 1) {
            ## create the hyper-parameter sequence
            ## lambda.max <- log10(max(svd(Y)$d)*2/(N*TT-sum(II)))
            lambda.max <- log10(max(svd(Y)$d))
            lambda <- rep(NA, nlambda)
            lambda.by <- 3/(nlambda - 2)
            for (i in 1:(nlambda - 1)) {
                lambda[i] <- 10^(lambda.max - (i - 1) * lambda.by)
            }
            lambda[nlambda] <- 0
        }
        
        ## store all MSPE
        CV.out <- matrix(NA, length(lambda), 3)
        colnames(CV.out) <- c("lambda", "sigma2", "MSPE")
        CV.out[,"lambda"] <- c(lambda)
        CV.out[,"MSPE"] <- 1e20
        for (i in 1:length(lambda)) {    
            k <- 5
            SSE <- 0
            for (ii in 1:k) {
                YY.cv <- YY
                II.cv <- II
                repeat{
                    cv.id <- sample(tot.id, as.integer(sum(II) - cv.count), replace = FALSE)
                    II.cv[cv.id] <- 0
                    con1 <- sum(apply(II.cv, 1, sum) > 0) == TT
                    con2 <- sum(apply(II.cv, 2, sum) > 0) == N
                    if (con1 & con2) {
                        break
                    }
                }
                YY.cv[cv.id] <- 0
                est.cv.fit <- inter_fe_mc(YY.cv, X, II.cv, 1, lambda[i], force, tol)$fit
                SSE <- SSE + sum((YY[cv.id]-est.cv.fit[cv.id])^2)
            }
            MSPE <- SSE/(k*(sum(II) - cv.count))

            est.cv <- inter_fe_mc(YY, X, II, 1, lambda[i], force, tol) ## overall
            sigma2 <- est.cv$sigma2 

            if(!is.null(norm.para)){
                MSPE <- MSPE*(norm.para[1]^2)
                sigma2 <- sigma2*(norm.para[1]^2)
            }

            if ((min(CV.out[,"MSPE"]) - MSPE) > tol*min(CV.out[,"MSPE"])) {
                ## at least 5% improvement for MPSE
                est.best <- est.cv  
                lambda.cv <- lambda[i]
            } else {
                if (i > 1) {
                    if (lambda.cv == lambda[i-1]) cat("*")
                }
            }
            CV.out[i, "MSPE"] <- MSPE
            CV.out[i, "sigma2"] <- sigma2 

            cat("\n lambda = ",
            sprintf("%.5f",lambda[i]),"; sigma2 = ",
            sprintf("%.5f",sigma2),"; MSPE = ",
            sprintf("%.5f",MSPE),sep="")

        } 
        cat("\n\n lambda* = ",lambda.cv, sep="")
        cat("\n\n")
        MSPE.best <- min(CV.out[,"MSPE"])
    } ## End of Cross-Validation

    validX <- est.best$validX
    validF <- est.best$validF
    
    ##-------------------------------##
    ## ATT and Counterfactuals 
    ##-------------------------------##
    
    ## variance of the error term
    if (is.null(norm.para)) {
        sigma2 <- est.best$sigma2   
    } else {
        sigma2 <- est.best$sigma2 * (norm.para[1]^2)       
    }
 
    ## effect
    Y.fit <- est.best$fit
    Y.ct <- as.matrix(Y.fit[,tr])
    if (0%in%I.tr) {
        Y.ct[which(I.tr==0)] <- 0 ## adjust    
    }
    eff <- Y.tr - Y.ct
    res <- est.best$residuals


    if (p>0) {
        beta <- est.best$beta
        if (est.best$validX == 0) {
            beta <- matrix(0, p, 1) 
        } else {
            beta <- est.best$beta
            beta[is.nan(est.best$beta)] <- 0
        }
    } else {
        beta <- NA
    }
    
    mu <- est.best$mu 
    Y.fe.bar <- rep(mu, TT)
    
    if (force%in%c(2,3)) {
        xi <- est.best$xi ## a (TT*1) matrix
        Y.fe.bar <- Y.fe.bar + xi
    }

    if (0%in%I.tr) {
        eff[which(I.tr == 0)] <- 0 ## adjust    
    } ## missing data will be adjusted to NA finally
   
    ##-------------------------------##
    ## Summarize
    ##-------------------------------##  
    
    ## counterfactuals and averages
    
    if (is.null(W)) {
        if (!0%in%I) {
            Y.tr.bar <- rowMeans(Y.tr)
            Y.ct.bar <- rowMeans(Y.ct)
            Y.co.bar <- rowMeans(Y.co)  
        } else {
            Y.tr.bar <- rowSums(Y.tr)/rowSums(I.tr)
            Y.ct.bar <- rowSums(Y.ct)/rowSums(I.tr)
            Y.co.bar <- rowSums(Y.co)/rowSums(I.co)
        }
    } else {
        Y.tr.bar <- rowSums(Y.tr * W.tr)/rowSums(W.tr)
        Y.ct.bar <- rowSums(Y.ct * W.tr)/rowSums(W.tr)
        Y.co.bar <- rowSums(Y.co)/rowSums(I.co)
    }

    ##Y.tr and Y.ct
    Y.bar <- cbind(Y.tr.bar, Y.ct.bar, Y.co.bar)
    colnames(Y.bar) <- c("Y.tr.bar", "Y.ct.bar", "Y.co.bar")
    
    ## ATT and average outcomes
    if (DID == TRUE) { ## diff-in-diffs: same timing
        if (is.null(W)) {
            if (!0%in%I.tr) {
                att <- rowMeans(eff)
            } else {
                att <- rowSums(eff)/rowSums(I.tr)
            }
        } else {
            att <- rowSums(eff * W.tr)/rowSums(W.tr)
        }
    } else { ## diff timing, centered the att
        if (!0%in%I.tr) {
            eff.cnt <- Y.tr.center <- matrix(NA, TT, Ntr)
            if (!is.null(W)) {
                W.tr.center <- matrix(NA, TT, Ntr)
            }
            for (j in 1:Ntr) {
                eff.cnt[1:(TT+T0.min-T0[j]), j] <- eff[(T0[j]-T0.min+1):TT, j]  
                Y.tr.center[1:(TT+T0.min-T0[j]), j] <- Y.tr[(T0[j]-T0.min+1):TT, j]
                if (!is.null(W)) {
                    W.tr.center[1:(TT+T0.min-T0[j]), j] <- W.tr[(T0[j]-T0.min+1):TT, j]   
                }
            }
            if (is.null(W)) {
                att <- apply(eff.cnt, 1, mean, na.rm = TRUE)
                Y.tr.cnt <- apply(Y.tr.center, 1, mean, na.rm = TRUE)
                Y.ct.cnt <- Y.tr.cnt - att
            } else {
                eff.cnt[which(is.na(eff.cnt))] <- 0
                Y.tr.center[which(is.na(Y.tr.center))] <- 0
                W.tr.center[which(is.na(W.tr.center))] <- 0
                att <- rowSums(eff.cnt * W.tr.center)/rowSums(W.tr.center)
                Y.tr.center <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
                Y.ct.cnt <- Y.tr.cnt - att
            }
        } else {
            T0.ub <- apply(as.matrix(D[,which(tr == 1)] == 0), 2, sum) 
            T0.ub.min <- min(T0.ub)
            eff.cnt <- Y.tr.center <- matrix(NA, TT, Ntr)
            eff[which(I.tr == 0)] <- NA
            Y.tr[which(I.tr == 0)] <- NA
            if (!is.null(W)) {
                W.tr.center <- matrix(NA, TT, Ntr)
                W.tr[which(I.tr == 0)] <- NA
            }
            for (j in 1:Ntr) {
                eff.cnt[1:(TT+T0.ub.min-T0.ub[j]), j] <- eff[(T0.ub[j]-T0.ub.min+1):TT, j]  
                Y.tr.center[1:(TT+T0.ub.min-T0.ub[j]),j] <- Y.tr[(T0.ub[j]-T0.ub.min+1):TT, j]
                if (!is.null(W)) {
                    W.tr.center[1:(TT+T0.ub.min-T0.ub[j]), j] <- W.tr[(T0.ub[j]-T0.ub.min+1):TT, j]   
                }
            }
            if (is.null(W)) {
                att <- apply(eff.cnt, 1, mean, na.rm = TRUE)
                Y.tr.cnt <- apply(Y.tr.center, 1, mean, na.rm = TRUE)
                Y.ct.cnt <- Y.tr.cnt - att
            } else {
                eff.cnt[which(is.na(eff.cnt))] <- 0
                Y.tr.center[which(is.na(Y.tr.center))] <- 0
                W.tr.center[which(is.na(W.tr.center))] <- 0
                att <- rowSums(eff.cnt * W.tr.center)/rowSums(W.tr.center)
                Y.tr.center <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
                Y.ct.cnt <- Y.tr.cnt - att
            }
        }
    }
    eff[which(is.na(eff))] <- 0 ## to calulate att
    if (is.null(W)) {
        att.avg <- sum(eff * post)/sum(post)
    } else {
        att.avg <- sum(eff * post * W.tr)/sum(post * W.tr)
    }

    ## AR1: calculate accumulative effect
    if (AR1 == TRUE) {
        rho <- est.best$beta[1]
        if (length(beta) > 1) {
            beta <- beta[-1]
        } 
        eff.tmp <- eff * D[, id.tr]
        eff.acc <- matrix(0, TT, Ntr)
        for (t in (T0.min + 1):TT) {
            for (i in 0:(t-T0.min-1)) {
                eff.acc[t,] <- eff.acc[t,] + eff.tmp[t-i,] * (rho^i)
            }
        }      
    } 
    
    ## final adjust unbalanced output
    if (0%in%I) {
        eff[which(I.tr == 0)] <- NA
        Y.ct[which(I.tr == 0)] <- NA
        Y.tr[which(I.tr == 0)] <- NA
        res[which(II == 0)] <- NA
        Y.co[which(I.co == 0)] <- NA
    }
    ## adjust beta: invariant covar
    if (p > 0) {
        if( sum(na.pos) > 0 ) {
            beta[na.pos] <- NA
        }
    }

    ## final adjustment
    if (!is.null(norm.para)) {
        mu <- mu * norm.para[1]
        ## if (p>0) {
        ##     beta <- beta*norm.para[1]/norm.para[2:length(norm.para)]
        ## }
        if (force%in%c(1, 3)) {
            est.best$alpha <- est.best$alpha * norm.para[1]
        }
        if (force%in%c(2,3)) {
            xi <- xi * norm.para[1]
        }
        res <- res * norm.para[1] 
        Y.tr <- Y.tr * norm.para[1] 
        Y.ct <- Y.ct * norm.para[1]
        Y.co <- Y.co * norm.para[1]
        eff <- eff * norm.para[1]
        Y.bar <- Y.bar * norm.para[1]
        att <- att * norm.para[1]
        att.avg <- att.avg * norm.para[1]
        if ( sameT0 == FALSE | 0%in%I.tr ) {
            eff.cnt <- eff.cnt * norm.para[1]
            Y.tr.cnt <- Y.tr.cnt * norm.para[1]
            Y.ct.cnt <- Y.ct.cnt * norm.para[1]
        }
    }

    T0 <- apply(as.matrix(D[,which(tr == 1)] == 0), 2, sum) ## for plot

    
    ##-------------------------------##
    ## Storage 
    ##-------------------------------##  

    ##control group residuals
    out<-list(
        ## main results
        D.tr = D.tr,
        I.tr = I.tr,
        Y.tr = Y.tr,
        Y.ct = Y.ct,
        Y.co = Y.co, 
        eff = eff,
        Y.bar = Y.bar,
        att = att,
        att.avg = att.avg,
        ## supporting
        force = force,
        DID = DID,
        T = TT,
        N = N,
        p = p,
        Ntr = Ntr,
        Nco = Nco,
        T0 = T0,
        tr = tr,
        pre = pre,
        post = post,
        lambda.cv = lambda.cv, 
        beta = beta,
        est = est.best,
        mu = mu,
        validX = validX,
        validF = validF,
        niter = est.best$niter
    )

    out <- c(out,list(sigma2 = sigma2, res=res))
    

    if ( DID == FALSE ) {
        out<-c(out,list(eff.cnt = eff.cnt,
                        Y.tr.cnt = Y.tr.cnt,
                        Y.ct.cnt = Y.ct.cnt))
    }
    if (CV) {
        out<-c(out, list(MSPE = MSPE.best,
                         CV.out = CV.out))
    } 
    if (force==1) {
        out<-c(out, list(alpha= est.best$alpha))
    } else if (force == 2) {
        out<-c(out,list(xi = xi))
    } else if (force == 3) {
        out<-c(out,list(alpha= est.best$alpha,
                        xi = xi))
    }
    if (AR1 == TRUE) {
        out<-c(out,list(rho = rho,
                        eff.acc = eff.acc))
    }

    return(out)
} ## Core functions ends


###############################################
## Inference 
###############################################

synth.boot<-function(Y,
                     X,
                     D, ## input
                     I,
                     W = NULL, 
                     EM, ## EM algorithm
                     r=0, r.end,
                     lambda = NULL,
                     nlambda = 10,
                     MC = FALSE,
                     force,
                     CV, ## cross validation
                     nboots,
                     tol,
                     inference, ## c("parametric","nonparametric")
                     cov.ar=1, 
                     AR1 = FALSE,
                     norm.para,
                     parallel = TRUE,
                     cores = NULL){
    
    
    na.pos <- NULL
    TT<-dim(Y)[1]
    N<-dim(Y)[2]
    if (is.null(X)==FALSE) {
        p<-dim(X)[3]
    } else {
        p<-0
    }

    ## treatement indicator
    tr<-D[TT,]==1  ## cross-sectional: treated unit
    co <- D[TT,] == 0 
    I.tr <- as.matrix(I[,tr])
    I.co <- I[,co]
    D.tr <- D[,tr]
    pre <- as.matrix(D[,which(tr==1)]==0&I[,which(tr==1)]!=0)
    post <- as.matrix(D[,which(tr==1)]!=0&I[,which(tr==1)]!=0)
                                         
    Ntr<-sum(tr)
    Nco<-N-Ntr
    T0<-apply(pre,2,sum) 
    T0.min<-min(T0)
    sameT0<-length(unique(T0))==1 ## treatment kicks in at the same time

    ## to calculate treated numbers
    T0.ub<-apply(as.matrix(D[,which(tr==1)]==0),2,sum) 
    T0.ub.min<-min(T0.ub)

    id<-1:N
    time<-1:TT
    id.tr<-which(tr==1) ## treated id
    id.co<-which(tr==0)

    ## vectorized "pre-treatment" indicator
    pre.v<-as.vector(pre)
    ## vectorized pre-treatment grouping variable for the treated
    id.tr.pre.v<-rep(id,each=TT)[which(pre.v==1)]
    ## create a list of pre-treatment periods
    time.pre<-split(rep(time,Ntr)[which(pre.v==1)],id.tr.pre.v) 
    
    ## estimation
    if (MC == FALSE) {
        if (EM == FALSE) {
            out<-synth.core(Y = Y, X = X, D = D, I=I, W=W, r = r, r.end = r.end, 
                            force = force, CV = CV, tol=tol,
                            AR1 = AR1, norm.para= norm.para, boot = 0)
        } else { # the case with EM
            if (CV == FALSE) {
                out<-synth.em(Y = Y,X = X, D = D, I=I, W=W, r = r, force = force,
                              tol = tol, AR1 = AR1, norm.para = norm.para, boot = 0)
            } else {
                out<-synth.em.cv(Y = Y, X = X, D = D, I=I, W=W, r = r, r.end = r.end,
                                 force = force, tol=tol,
                                 AR1 = AR1, norm.para = norm.para)
            }
        }
        ## for parametric bootstarp: some control group units may not be suitble
        if (inference == "parametric") {
            co.pre <- apply(as.matrix(I.co[1:T0.ub.min,]),2,sum)
            co.post <- apply(as.matrix(I.co[(max(T0.ub)+1):TT,]),2,sum)
            if (force%in%c(1,3)) {
                valid.co <- id.co[(co.pre >= out$r.cv+1)&(co.post >= 1)]
            } else {
                valid.co <- id.co[(co.pre >= out$r.cv)&(co.post >= 1)]
            }
        }
    } else {
        out<-synth.mc(Y = Y, X = X, D = D, I=I, W=W, 
                      lambda = lambda, nlambda = nlambda, 
                      force = force, tol=tol,
                      AR1 = AR1, norm.para= norm.para)
    }


    ## output
    validX <- out$validX
    eff<-out$eff
    att<-out$att
    att.avg<-out$att.avg
    DID <- out$DID

    if (p > 0) {
        beta <- out$beta
        if (NA%in%beta) {
            if (sum(is.na(beta))<p) {
                beta.it <- as.matrix(beta[which(!is.na(beta))])
            } else {
                beta.it <- matrix(0,0,1)
            }
        } else {
            beta.it <- beta
        }
    } else {
        beta.it <- beta <-matrix(0,0,1)
    }
    
    if (MC == FALSE) {
        error.co<-out$res.co ## error terms (T*Nco) contains NA unbalanced data
    } else {
        error <- out$res
    }
    

    Y.tr.bar=out$Y.bar[,1]
    Y.ct.bar=out$Y.bar[,2]
    Y.co.bar=out$Y.bar[,3]
    
    ## bootstrapped estimates
    eff.boot<-array(0,dim=c(TT,Ntr,nboots))  ## to store results
    att.boot<-matrix(0,TT,nboots)
    att.avg.boot<-matrix(0,nboots,1)
    if (p>0) {
        beta.boot<-matrix(0,p,nboots)
    }
    
    if (inference=="nonparametric") { ## nonparametric bootstrap

        cat("\rBootstrapping ...\n")
        if (MC == FALSE) {
            if (EM == FALSE) {
                one.nonpara <- function(){

                    repeat {
                        fake.co <- sample(id.co,Nco, replace=TRUE)
                        if (sum(apply(as.matrix(I[,fake.co]),1,sum)>=1)==TT) {
                            break
                        }
                    }

                    boot.id<-c(sample(id.tr,Ntr,replace=TRUE), fake.co)
                
                    X.boot<-X[,boot.id,,drop=FALSE]
                    W.boot <- NULL
                    if (!is.null(W)) {
                        W.boot <- W[,boot.id]
                    }
                    boot<-synth.core(Y[,boot.id], X.boot, D[,boot.id], I=I[,boot.id],
                                     W = W.boot, force = force, r = out$r.cv, CV=0,
                                     tol = tol, AR1 = AR1,
                                     beta0 = beta.it, norm.para = norm.para, boot = 1)
                    return(boot)
                
                } 
            } else { # the case of EM
                one.nonpara <- function(){
                
                    repeat {
                        fake.co <- sample(id.co,Nco, replace=TRUE)
                        if (sum(apply(as.matrix(I[,fake.co]),1,sum)>=1)==TT) {
                            break
                        }
                    }

                    boot.id<-c(sample(id.tr,Ntr,replace=TRUE), fake.co)
                
                    X.boot<-X[,boot.id,,drop=FALSE]
                    W.boot <- NULL
                    if (!is.null(W)) {
                        W.boot <- W[,boot.id]
                    }
                    boot<-synth.em(Y = Y[,boot.id], X = X.boot, D = D[,boot.id], I=I[,boot.id],
                                   W = W.boot, force = force, r = out$r.cv,
                                   tol = tol, AR1 = AR1, norm.para = norm.para, boot = 1)
                    return(boot)
                
                } 
            }
        } else { ## mc
            one.nonpara <- function(){
                repeat {
                    fake.co <- sample(id.co,Nco, replace=TRUE)
                    boot.id<-c(sample(id.tr,Ntr,replace=TRUE), fake.co)
                    con1 <- sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT
                    con2 <- sum(apply(as.matrix(I[,boot.id]),2,sum)>=1)==N
                    if (con1 & con2) {
                        break
                    }
                }
                
                X.boot<-X[,boot.id,,drop=FALSE]
                W.boot <- NULL
                if (!is.null(W)) {
                    W.boot <- W[,boot.id]
                }
                boot<-synth.mc(Y[,boot.id], X.boot, D[,boot.id], I=I[,boot.id],
                               W = W.boot, force = force, 
                               lambda = out$lambda.cv, hasF = out$validF, 
                               CV = 0, tol = tol, AR1 = AR1, norm.para = norm.para)
                return(boot)
            
            } 
        }
        ## computing
        if (parallel == TRUE) { 
            boot.out <- foreach(j=1:nboots, 
                                .inorder = FALSE,
                                .export = c("synth.core","synth.em","synth.mc"),
                                .packages = c("gsynth")
                                ) %dopar% {
                                    return(one.nonpara())
                                }

            for (j in 1:nboots) { 
                att.boot[,j]<-boot.out[[j]]$att
                att.avg.boot[j,]<-boot.out[[j]]$att.avg  
                if (p>0) {
                    beta.boot[,j]<-boot.out[[j]]$beta
                } 
            } 
        } else {
            for (j in 1:nboots) { 
                boot <- one.nonpara() 
                att.boot[,j]<-boot$att
                att.avg.boot[j,]<-boot$att.avg
                if (p>0) {
                    beta.boot[,j]<-boot$beta
                }
                ## report progress
                if (j%%100==0)  {
                    cat(".")   
                }  
            }  
        } 
        ## end of bootstrapping
        cat("\r")
        
    } else if (inference=="parametric") { ## end of non-parametric
        
        ## library("mvtnorm") ## generate multivariate normal distribution residual
        if (EM == FALSE) { # the case without EM
            ## y fixed
            if (is.null(norm.para)) {
                error.co<-out$res.co ## contains NA unbalanced data 
                Y.fixed<-Y
                Y.fixed[,id.tr]<-as.matrix(out$Y.ct)
                Y.fixed[,id.co]<-Y.fixed[,id.co]-error.co
            } else {
                error.co<-out$res.co/norm.para[1]
                Y.fixed<-Y
                Y.fixed[,id.tr]<-as.matrix(out$Y.ct/norm.para[1])
                Y.fixed[,id.co]<-Y.fixed[,id.co]-error.co
            }
            
            Y.fixed[which(I==0)] <- 0

            draw.error <- function() {
                ## draw 1 prediction error at a time      
                repeat {
                    fake.tr<-sample(id.co,1,replace=FALSE)
                    if (fake.tr%in%valid.co) {
                        break
                    }
                }
            
                id.co.rest<-id.co[which(!id.co%in%fake.tr)]
                ## resample control, to smooth CV prediction error
                repeat {
                    id.co.pseudo <- sample(id.co.rest, Nco, replace=TRUE)
                    if (sum(apply(as.matrix(I[,id.co.pseudo]),1,sum)>=1)==TT) {
                        break
                    }
                }
                      
                id.pseudo<-c(rep(fake.tr,Ntr),id.co.pseudo)  ## Ntr + ...
                I.id.pseudo<-I[,id.pseudo] 
                
                ## obtain the prediction eror
                D.pseudo<-D[,c(id.tr,id.co.pseudo)]  ## fake.tr + control left
                Y.pseudo<-Y[,id.pseudo]
                X.pseudo<-X[,id.pseudo,,drop=FALSE]
                W.pseudo <- NULL
                if (!is.null(W)) {
                    W.boot <- W[,id.pseudo]
                }
                ## output
                synth.out <- synth.core(Y = Y.pseudo, X = X.pseudo, D = D.pseudo,
                                        I = I.id.pseudo, W = W.pseudo,
                                        force = force, r = out$r.cv, CV = 0,
                                        tol = tol, AR1 = AR1, beta0 = beta.it,
                                        norm.para = norm.para, boot = 1)
                if (is.null(norm.para)) {
                    output <- synth.out$eff
                } else {
                    output <- synth.out$eff/norm.para[1]
                }
                           
                return(as.matrix(output)) ## TT * Ntr
                
            }

            cat("\rSimulating errors ...")
            if (parallel == TRUE) {
                error.tr <- foreach(j = 1:nboots,
                                    .combine = function(...) abind(...,along=3),
                                    .multicombine=TRUE,
                                    .export = c("synth.core"),
                                    .packages = c("gsynth"),
                                    .inorder = FALSE)  %dopar% {
                                        return(draw.error())
                                    } 
            } else {
                error.tr<-array(NA,dim=c(TT,Ntr,nboots))
                for (j in 1:nboots) {
                    error.tr[,,j] <- draw.error()
                    if (j%%100==0) {
                        cat(".")
                    }
                }
            }            
            
            if (0%in%I) {
                ## calculate vcov of ep_tr
                error.tr.adj <- array(NA,dim=c(TT,nboots,Ntr))
                for(i in 1:Ntr){
                    error.tr.adj[,,i] <- error.tr[,i,]
                }
                vcov_tr<-array(NA,dim=c(TT,TT,Ntr))
                for(i in 1:Ntr){
                    vcov_tr[,,i]<-res.vcov(res=error.tr.adj[,,i],
                                           cov.ar=cov.ar)
                    vcov_tr[,,i][is.na(vcov_tr[,,i])|is.nan(vcov_tr[,,i])] <- 0
                }
                
                ## calculate vcov of e_co
                vcov_co <- res.vcov(res=error.co,cov.ar=cov.ar)
                vcov_co[is.na(vcov_co)|is.nan(vcov_co)] <- 0
            }

            one.boot <- function(){
                ## boostrap ID
                repeat {
                    fake.co <- sample(id.co,Nco, replace=TRUE)
                    if (sum(apply(as.matrix(I[,fake.co]),1,sum)>=1)==TT) {
                        break
                    }
                }
                id.boot<-c(id.tr, fake.co)
                
                ## get the error for the treated and control
                error.tr.boot<-matrix(NA,TT,Ntr)
                if (0%in%I) {
                    
                    for (w in 1:Ntr) {
                        error.tr.boot[,w]<-t(rmvnorm(n=1,rep(0,TT),vcov_tr[,,w],method="svd"))
                    }
                    
                    error.tr.boot[which(I.tr==0)] <- 0
                    
                    error.co.boot <- 
                        t(rmvnorm(n=Nco,rep(0,TT),vcov_co,method="svd"))

                    error.co.boot[which(as.matrix(I[,fake.co])==0)] <- 0
                    
                } else {
                    for (w in 1:Ntr) {
                        error.tr.boot[,w]<-error.tr[,w,sample(1:nboots,1,replace=TRUE)]
                    }
                    error.co.boot<-error.co[,sample(1:Nco,Nco,replace=TRUE)] 
                    
                }

                Y.boot<-Y.fixed[,id.boot]
                Y.boot[,1:Ntr]<- as.matrix(Y.boot[,1:Ntr] + error.tr.boot)
                ## new treated: conterfactual+effect+ (same) new error
                Y.boot[,(Ntr+1):length(id.boot)]<-
                Y.boot[,(Ntr+1):length(id.boot)] + error.co.boot
                
                X.boot<-X[,id.boot,,drop=FALSE] 
                D.boot<-D[,id.boot] 
                I.boot<-I[,id.boot]
                W.boot <- NULL
                if (!is.null(W)) {
                    W.boot <- W[,id.boot]
                }
                
                ## re-estimate the model 
                boot <- synth.core(Y.boot, X.boot, D.boot, I=I.boot,
                                   W = W.boot, force = force, r = out$r.cv,
                                   CV = 0, tol = tol, AR1 = AR1,
                                   beta0 = beta.it, norm.para = norm.para, boot = 1)

                b.out <- list(eff = boot$eff + out$eff,
                              att = boot$att + out$att,
                              att.avg = boot$att.avg + out$att.avg)
                if (p>0) {
                    b.out <- c(b.out, list(beta = boot$beta))
                }
                return(b.out)
            }     
        } else { # the case of EM
            ## y fixed
            if (is.null(norm.para)) {
                error.co<-out$res.co ## contains NA unbalanced data 
                Y.fixed<-Y
                Y.fixed[,id.tr]<-as.matrix(out$Y.ct)
                Y.fixed[,id.co]<-Y.fixed[,id.co]-error.co
            } else {
                error.co<-out$res.co/norm.para[1]
                Y.fixed<-Y
                Y.fixed[,id.tr]<-as.matrix(out$Y.ct/norm.para[1])
                Y.fixed[,id.co]<-Y.fixed[,id.co]-error.co
            }
        
            Y.fixed[which(I==0)] <- 0

            if (0%in%I) {
                vcov_co <- res.vcov(res=error.co,cov.ar=cov.ar)
                vcov_co[is.na(vcov_co)|is.nan(vcov_co)] <- 0
            }
            
            one.boot <- function() {

                ## sample errors
                error.id <- sample(1:Nco, N, replace = TRUE)
                
                ## produce the new outcome data
                if (0%in%I) {
                    error.boot <- 
                        t(rmvnorm(n=N,rep(0,TT),vcov_co,method="svd"))
                    Y.boot <- Y.fixed + error.boot   
                } else {
                    Y.boot<-Y.fixed + error.co[,error.id]
                }
                                
                ## re-estimate the model
                boot<-synth.em(Y.boot, X, D, I=I, W=W, force=force, r=out$r.cv,
                               tol=tol, AR1 = AR1, norm.para = norm.para, boot = 1)

                b.out <- list(eff = boot$eff + out$eff,
                              att = boot$att + out$att,
                              att.avg = boot$att.avg + out$att.avg)
                if (p>0) {
                    b.out <- c(b.out, list(beta = boot$beta))
                } 
                return(b.out) 
            }
             
        } # the end of the EM case

        ## computing
        cat("\rBootstrapping ...\n")
        if (parallel == TRUE) { 
            boot.out <- foreach(k=1:nboots,
                                .inorder = FALSE,
                                .export = c("synth.core","synth.em"),
                                .packages = c("gsynth")
                                ) %dopar% {
                                    return(one.boot())
                                }
            for (j in 1:nboots) {
                eff.boot[,,j]<-boot.out[[j]]$eff
                att.boot[,j]<-boot.out[[j]]$att
                att.avg.boot[j,]<-boot.out[[j]]$att.avg
                if (p>0) {
                    beta.boot[,j]<-boot.out[[j]]$beta
                }
            }
        } else {
            for (j in 1:nboots) {
                boot.out <- one.boot()
                eff.boot[,,j]<-boot.out$eff
                att.boot[,j]<-boot.out$att
                att.avg.boot[j,]<-boot.out$att.avg
                if (p>0) {
                    beta.boot[,j]<-boot.out$beta
                }
                if (j%%100==0) {
                    cat(".")
                }
            }
        } 
        cat("\r")
    }
     
   
    ####################################
    ## Variance and CIs
    ####################################

    ## function to get two-sided p-values
    get.pvalue <- function(vec) {
        if (NaN%in%vec|NA%in%vec) {
            nan.pos <- is.nan(vec)
            na.pos <- is.na(vec)
            pos <- c(which(nan.pos),which(na.pos))
            vec.a <- vec[-pos]
            a <- sum(vec.a >= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2
            b <- sum(vec.a <= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2  
        } else {
            a <- sum(vec >= 0)/length(vec) * 2
            b <- sum(vec <= 0)/length(vec) * 2  
        }
        return(min(as.numeric(min(a, b)),1))
    }
    
    ## ATT estimates
    CI.att <- t(apply(att.boot, 1, function(vec) quantile(vec,c(0.025,0.975), na.rm=TRUE)))
    se.att <- apply(att.boot, 1, function(vec) sd(vec, na.rm=TRUE))
    pvalue.att <- apply(att.boot, 1, get.pvalue)

    if (DID == TRUE) {
        ntreated <- apply(post, 1, sum)
    } else {
        if (!0%in%I.tr) {
            rawcount <- apply(1-pre, 1, sum)
            ntreated <- c(rep(0, T0.min), rev(rawcount[(T0.min + 1): TT]))
        } else {
            Itr.sub <- I.tr[(T0.ub.min+1):TT,]
            Itr.count <- matrix(0,(TT-T0.ub.min),Ntr)
            for(i in 1:Ntr){
                Itr.count[1:(TT-T0.ub[i]),i] <- Itr.sub[(T0.ub[i]+1-T0.ub.min):(TT-T0.ub.min),i]
            }
            rawcount <- apply(Itr.count, 1, sum)
            ## ntreated <- c(rep(0, T0.ub.min), rev(rawcount[(T0.ub.min + 1): TT]))
            ntreated <- c(rep(0, T0.ub.min), rawcount)
        }
    }
    est.att <- cbind(att, se.att, CI.att, pvalue.att, ntreated)
    colnames(est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                           "p.value", "n.Treated")
    if (DID == TRUE) {
        rownames(est.att) <- time
    } else {
        if (!0%in%I.tr) {
            rownames(est.att) <- c(1:TT) - min(T0)
        } else {
            rownames(est.att) <- c(1:TT) - min(T0.ub)
        }
    }
    
    ## average (over time) ATT
    CI.avg <- quantile(att.avg.boot, c(0.025,0.975), na.rm=TRUE)
    se.avg <- sd(att.avg.boot, na.rm=TRUE)
    pvalue.avg <- get.pvalue(att.avg.boot)
    est.avg <- t(as.matrix(c(att.avg, se.avg, CI.avg, pvalue.avg)))
    colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")
    
    ## individual effects
    if (inference == "parametric") {
        CI.ind <- apply(eff.boot,c(1,2),function(vec)
            quantile(vec,c(0.025,0.975), na.rm=TRUE)) ## 2*T*Ntr
        est.ind <- array(NA,dim=c(TT, 5, Ntr)) ## eff, se, CI.lower, CI.upper
        est.ind[,1,] <- eff
        est.ind[,2,] <- apply(eff.boot,c(1,2),sd)
        est.ind[,3,] <- CI.ind[1,,]
        est.ind[,4,] <- CI.ind[2,,]
        est.ind[,5,] <- apply(eff.boot,c(1,2),get.pvalue)
    }

    
    ## regression coefficents
    if (p>0) {
        CI.beta<-t(apply(beta.boot, 1, function(vec)
            quantile(vec,c(0.025, 0.975), na.rm=TRUE)))
        se.beta<-apply(beta.boot, 1, function(vec)sd(vec,na.rm=TRUE))
        pvalue.beta <- apply(beta.boot, 1, get.pvalue)
        beta[na.pos] <- NA
        est.beta<-cbind(beta, se.beta, CI.beta, pvalue.beta)
        colnames(est.beta)<-c("beta", "S.E.", "CI.lower", "CI.upper", "p.value")
    }
  
    ##storage
    result<-list(inference = inference,
                 est.att = est.att,
                 est.avg = est.avg,
                 att.boot = att.boot
                 )
    if (p>0) {
        result <- c(result,list(beta.boot = beta.boot))
    }
    
    if (inference == "parametric") {
        result<-c(result,list(est.ind = est.ind,
                              eff.boot = eff.boot))
    }
    if (p>0) {
        result<-c(result,list(est.beta = est.beta))
    } 

    return(c(out,result))

    
} ## end of synth.boot()

 
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


##########
## Plot
##########
#x a gsynth object
# type of the plot; axes limits; axes labels; 
# show raw data in "counterfactual" mode # ("none","band","all")
# main: whether to show the title;
# nfactors: whose loadings to be plotted 
# id: individual plot
plot.gsynth <- function(x,  
                        type = "gap", 
                        xlim = NULL, 
                        ylim = NULL,
                        xlab = NULL, 
                        ylab = NULL,
                        legendOff = FALSE,
                        raw = "none", 
                        main = NULL,
                        nfactors = NULL, 
                        id = NULL,
                        axis.adjust = FALSE,
                        ...){


    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------##  

    outcome <- NULL
    ATT <- NULL
    CI.lower <- NULL
    CI.upper <- NULL
    co5 <- NULL
    co95 <- NULL
    tr5 <- NULL
    tr95 <- NULL
    group <- NULL
    L1 <- NULL
    out <- NULL

    if (class(x)!="gsynth") {
        stop("Not a \"gsynth\" object.")
    }
    if (!type %in% c("gap","counterfactual","factors","missing","loadings","raw")) {
        stop("\"type\" option misspecified.")
    }
    if (is.null(x$factor) & type == "factors") {
        stop("No factors to be plotted.")
    }
    if (is.null(x$lambda.tr) & type == "factors") {
        stop("No loadings to be plotted.")
    }

    if (is.null(xlim)==FALSE) {
        if (is.numeric(xlim)==FALSE) {
            stop("Some element in \"xlim\" is not numeric.")
        } else {
            if (length(xlim)!=2) {
                stop("xlim must be of length 2.")
            }
        }
    }
    if (is.null(ylim)==FALSE) {
        ## if (type!="missing") {
            if (is.numeric(ylim)==FALSE) {
                stop("Some element in \"ylim\" is not numeric.")
            } else {
                if (length(ylim)!=2) {
                    stop("ylim must be of length 2.")
                }
            }
        ## } else {
        ##     m.l <- length(ylim)
        ##     for (i in 1:m.l) {
        ##         if (!ylim[m.l]%in%x$id) {
        ##             stop("Some specified units are not in the data.")
        ##         }
        ##     }
        ## }
    }

    if (is.null(xlab)==FALSE) {
        if (is.character(xlab) == FALSE) {
            stop("\"xlab\" is not a string.")
        } else {
            xlab <- xlab[1]
        }   
    }
    if (is.null(ylab)==FALSE) {
        if (is.character(ylab) == FALSE) {
            stop("\"ylab\" is not a string.")
        } else {
            ylab <- ylab[1]
        }   
    }
    if (is.logical(legendOff) == FALSE & is.numeric(legendOff)==FALSE) {
        stop("\"legendOff\" is not a logical flag.")
    }
    if (type == "counterfactual") {
        if (! raw %in% c("none","band","all")) {
            cat("\"raw\" option misspecifed. Reset to \"none\".")
            raw <- "none" 
        }
        if (is.null(id)==FALSE) {
            if (length(id)>1) {
               stop("More than 1 element in \"id\".") 
            }
        } 
    }
    if (is.null(main)==FALSE) {
        if (is.character(main) == FALSE) {
            stop("\"main\" is not a string.")
        } else {
            main <- main[1]
        }   
    }
    if (is.null(nfactors)==FALSE) {
        if (is.numeric(nfactors)==FALSE) {
            stop("\"nfactors\" is not a positive integer.")
        } else {
            nfactors <- nfactors[1]
            if (nfactors%%1!=0 | nfactors<=0) {
                stop("\"nfactors\" is not a positive integer.")
            }  
        } 
    }

    if (axis.adjust==TRUE) {
        angle <- 45
        x.v <- 1
        x.h <- 1
    } else {
        angle <- 0
        x.v <- 0
        if (type=="missing") {
            x.h <- 0.5
        } else {
            x.h <- 0
        }
    }
    
    ##-------------------------------##
    ## Plotting
    ##-------------------------------##  

    I.tr <- x$I.tr
    D.tr <- x$D.tr
    Y.tr <- x$Y.tr
    Y.co <- x$Y.co
    Y.ct <- x$Y.ct
    tb <- x$est.att
    Yb <- x$Y.bar[,1:2] ## treated average and counterfactual average
    tr <- x$tr
    pre <- x$pre
    post <- x$post
    # I.tr <- x$I.tr
    TT <- x$T
    T0 <- x$T0 ## notice
    p <- x$p
    ## m <- x$m
    Ntr <- x$Ntr
    Nco <- x$Nco
    N <- x$N 
    force <- x$force
    F.hat <- x$factor
    L.tr <- x$lambda.tr

    ## time.label <- x$time
    ## T.b <- 1:TT
    if (!is.null(L.tr)) {
        r <- dim(L.tr)[2]
    } else {
        r <- 0
    }
    
    if (type!="missing") {
        if (is.null(id)==TRUE) {
            id <- x$id.tr
        }
    } else {
        if (is.null(id)==TRUE) {
            id <- colnames(x$obs.missing)
        }
        m.l <- length(id)
            for (i in 1:m.l) {
                if (!id[i]%in%colnames(x$obs.missing)) {
                    stop("Some specified units are not in the data.")
                }
        }
    }

    ## parameters
    line.width <- c(1.2,0.5)
  
    ## type of plots
    if (type == "raw"| type == "counterfactual" | 
        type == "factors" |  length(id) == 1 | type =="missing" | 
        type=="loadings") {
        time <- x$time
        if (!is.numeric(time[1])) {
            time <- 1:TT
        }
        
        if (type!="missing") {
            if (length(id) == 1) {
                time.bf <- time[T0[which(id == x$id.tr)]]
            } else {
                time.bf <- time[unique(T0)]
            }
        }

        ## periods to show
        if (length(xlim) != 0) {
            ## if(is.numeric(time[1])){
                show <- which(time>=xlim[1]& time<=xlim[2])
            ## } else {
            ##     xlim[1] <- which(x$time>=xlim[1])[1]
            ##     xlim[2] <- which(x$time<=xlim[2])[length(which(x$time<=xlim[2]))]
            ##     show <- which(time>=xlim[1]& time<=xlim[2])

            ## }
        } else {
            show <- 1:length(time)
        }     
    }

    if (type == "gap")  { ## variable treatment timing
        time <- c(1:TT) - min(T0)
        time.bf <- 0 ## before treatment

        if (length(xlim) != 0) {
            show <- which(time>=xlim[1]& time<=xlim[2])     
        } else {
            show <- 1:length(time)    
        }
    }

    nT <- length(show)
    time.label <- x$time[show]

    ## if (axis.adjust==FALSE) {
    ##     n.period <- length(show)
    ## } else {
    ##     n.period <- length(show) ## min(length(show),20)
    ## }

    ## if (axis.adjust==TRUE) {
    ##     n.period <- n.period - 1
    ##     T.n <- (nT-1)%/%n.period
    ##     T.res <- (nT-1)%%n.period
    ##     T.b <- seq(from=1,to=T.n*n.period+1,by=T.n)
    ##     if (T.res!=0) {
    ##         T.j <- 1
    ##         for(i in (n.period-T.res+2):(n.period+1)) {
    ##             T.b[i] <- T.b[i] + T.j
    ##             T.j <- T.j + 1
    ##         }
    ##     }
        ## T.b <- show[T.b]
    ## } else {
        T.b <- 1:length(show)
    ## }
 

    ## legend on/off
    if (legendOff == TRUE) {
        legend.pos <- "none"
    } else {
        legend.pos <- "bottom"
    }

    ############  START  ###############
    
    if (type == "raw") {
        ## axes labels
        if (is.null(xlab)==TRUE) {
            xlab <- x$index[2]
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab)==TRUE) {
            ylab <- x$Yname
        } else if (ylab == "") {
            ylab <- NULL
        }
            
        pst <- D.tr
        for (i in 1:Ntr){
            pst[T0[i],i] <- 1 ## paint the period right before treatment
        }
        time.pst <- c(pst[show,] * time[show])
        time.pst <- time.pst[which(c(pst[show,])==1)]
        Y.tr.pst <- c(Y.tr[show,])[which(pst[show,]==1)]
        id.tr.pst <- matrix(rep(1:Ntr,each=TT),TT,Ntr,byrow=FALSE)[show,]
        id.tr.pst <- c(id.tr.pst)[which(pst[show,]==1)]

        data <- cbind.data.frame("time" = c(rep(time[show], N), time.pst),
                                 "outcome" = c(c(Y.tr[show,]),
                                               c(Y.co[show,]),
                                               Y.tr.pst),
                                 "type" = c(rep("tr",(Ntr*nT)),
                                            rep("co",(Nco*nT)),
                                            rep("tr.pst",length(Y.tr.pst))),
                                 "id" = c(rep(1:N,each = nT), id.tr.pst*(-1)))
        
        ## theme
        p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
            theme(legend.position = legend.pos,
                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.h),
                  plot.title = element_text(size=20,
                                            hjust = 0.5,
                                            face="bold",
                                            margin = margin(10, 0, 10, 0)))

        
        
        if (x$DID==TRUE) {
            p <- p + geom_vline(xintercept=time.bf,colour="white",size = 2) +
                annotate("rect", xmin= time.bf, xmax= Inf,
                         ymin=-Inf, ymax=Inf, alpha = .3) 
        }
        
        ## main
        p <- p + geom_line(aes(time, outcome,
                               colour = type,
                               size = type,
                               linetype = type,
                               group = id))

        ## legend
        set.limits = c("tr","tr.pst","co")
        set.labels = c("Treated (Pre)",
                       "Treated (Post)",
                       "Controls")
        set.colors = c("#FC8D6280","red","#99999950")
        set.linetypes = c("solid","solid","solid")
        set.linewidth = c(0.5, 0.5, 0.5)
        
        p <- p + scale_colour_manual(limits = set.limits,
                                     labels = set.labels,
                                     values =set.colors) +
            scale_linetype_manual(limits = set.limits,
                                  labels = set.labels,
                                  values = set.linetypes) +
            scale_size_manual(limits = set.limits,
                              labels = set.labels,
                              values = set.linewidth) +
            guides(linetype = guide_legend(title=NULL, ncol=3),
                   colour = guide_legend(title=NULL, ncol=3),
                   size = guide_legend(title=NULL, ncol=3)) 
        if (!is.numeric(time.label)) {
            p <- p + 
            scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
        }

        ## title
        if (is.null(main) == TRUE) {
            p <- p + ggtitle("Raw Data")
        } else if (main!="") {
            p <- p + ggtitle(main)
        }

        ## ylim
        if (is.null(ylim) == FALSE) {
            p <- p + coord_cartesian(ylim = ylim)
        }
        suppressWarnings(print(p))
        
    } else if (type == "gap") { 
        
        if (length(id) == 1 & !(id[1] %in% x$id.tr)) { ## error
            cat(paste(id,"not in the treatment group"))
        } else { ## no error

            ## axes labels
            if (is.null(xlab) == TRUE) {
                if (x$DID == TRUE) {
                    xlab <- x$index[2]
                } else {
                    xlab <- paste("Time relative to Treatment")
                }
            } else if (xlab == "") {
                xlab <- NULL
            }
            if (is.null(ylab) == TRUE) {
                ylab <- "Coefficient"
            } else if (ylab == "") {
                ylab <- NULL
            }
            
            ## title
            if (length(id) == 1) { ## id specified
                maintext <- paste(x$index[1],"=",id) 
            }  else {
                maintext <- "Estimated Average Treatment Effect"
            } 
            
            ## contruct data for plotting
            if (is.null(x$est.att)==TRUE) { 
                cat("Uncertainty estimates not available.\n")
                if (length(id) == 1) { ## id specified
                    data <- cbind.data.frame(time, x$eff)[show,]
                    colnames(data) <- c("time","ATT")
                } else {
                    data <- cbind.data.frame(time, ATT = x$att)[show,] 
                } 
            } else {
                if (length(id) == 1) { ## id specified
                    id <- which(x$id.tr == id)
                    tb <- x$est.ind[,,id]
                    time.bf <- time[T0[id]] 
                    time <- time - time.bf
                    time.bf <- 0
                    if (!is.null(tb)) {
                        colnames(tb) <- c("ATT", "S.E.", "CI.lower", "CI.upper","p.value")
                    } else {
                        tb <- as.matrix(x$eff[,id])
                        colnames(tb) <- "ATT" 
                    } 
                }
                data <- cbind.data.frame(time, tb)[show,]
            }
             
            ## plotting
            p <- ggplot(data) +
                geom_vline(xintercept = time.bf, colour="white",size = 2) +
                geom_hline(yintercept = 0, colour="white",size = 2) +
                ## annotate("rect", xmin= time.bf, xmax= Inf,
                ##          ymin=-Inf, ymax=Inf, alpha = .1,
                ##          fill = "yellow") +
                xlab(xlab) +  ylab(ylab) +
                theme(legend.position = legend.pos,
                      plot.title = element_text(size=20,
                                                hjust = 0.5,
                                                face="bold",
                                                margin = margin(10, 0, 10, 0)))
           
            
            ## point estimates
            p <- p + geom_line(aes(time, ATT), size = 1.2)
             
            ## confidence intervals
            if (is.null(x$est.att)==FALSE) {
                if(!(is.null(x$est.ind)&length(id) == 1)) {
                    p <- p + geom_ribbon(aes(x = time, ymin=CI.lower, ymax=CI.upper),alpha=0.2)
                }
            }
            
            ## title
            if (is.null(main) == TRUE) {
                p <- p + ggtitle(maintext)
            } else if (main!=""){
                p <- p + ggtitle(main)
            }

            ## ylim
            if (is.null(ylim) == FALSE) {
                p <- p + coord_cartesian(ylim = ylim)
            }
            suppressWarnings(print(p))
        }  ## end of "gap" (in case of no id error)
       
        
    } else if (type=="counterfactual") { 

        if (length(id) == 1|length(x$id.tr) == 1|x$DID==TRUE) { 
            if (length(id)==1 & !(id[1]%in%x$id.tr)) { ## error
            
                cat(paste(id,"not in the treatment group"))
            
            } else { ## one treated unit case

                ## axes labels
                if (is.null(xlab)==TRUE) {
                    xlab <- x$index[2]
                } else if (xlab == "") {
                    xlab <- NULL
                }
                if (is.null(ylab)==TRUE) {
                    ylab <- x$Yname
                } else if (ylab == "") {
                    ylab <- NULL
                }
             
                if (length(id) == 1 | length(x$id.tr) == 1) { ## one treated unit
  
                    if (is.null(id) == TRUE) {
                        id <- x$id.tr
                    }
                    maintext <- paste("Treated and Counterfactual (",id,")",sep="") 
                    tr.info <- Y.tr[,which(id==x$id.tr)]
                    ct.info <- Y.ct[,which(id==x$id.tr)] 
                    if (raw == "none") { 
                        data <- cbind.data.frame("time" = rep(time[show],2),
                                                 "outcome" = c(tr.info[show],
                                                               ct.info[show]),
                                                 "type" = c(rep("tr",nT),
                                                            rep("ct",nT)))
                        ## theme
                        p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour="white",size = 2) +
                            annotate("rect", xmin= time.bf, xmax= Inf,
                                     ymin=-Inf, ymax=Inf, alpha = .3) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0))) 
                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type)) 
                        ## legend
                        set.limits = c("tr","ct")
                        set.labels = c("Treated", "Estimated Y(0)")
                        set.colors = c("red","steelblue")
                        set.linetypes = c("solid","longdash")
                        set.linewidth = rep(line.width[1],2)
                        p <- p + scale_colour_manual(limits = set.limits,
                                                     labels = set.labels,
                                                     values =set.colors) +
                            scale_linetype_manual(limits = set.limits,
                                                  labels = set.labels,
                                                  values = set.linetypes) +
                            scale_size_manual(limits = set.limits,
                                              labels = set.labels,
                                              values = set.linewidth) +
                            guides(linetype = guide_legend(title=NULL, ncol=2),
                                   colour = guide_legend(title=NULL, ncol=2),
                                   size = guide_legend(title=NULL, ncol=2))
                        if (!is.numeric(time.label)) {
                            p <- p + 
                            scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                        } 

                    
                    } else if  (raw == "band") {

                        Y.co.90 <- t(apply(Y.co, 1, quantile, prob=c(0.05,0.95), na.rm = TRUE)) 
                        data <- cbind.data.frame("time" = rep(time[show],2),
                                                 "outcome" = c(tr.info[show],
                                                               ct.info[show]),
                                                 "type" = c(rep("tr",nT),
                                                            rep("ct",nT)))

                        data.band <- cbind.data.frame(time, Y.co.90)[show,]
                        colnames(data.band) <- c("time","co5","co95")

                    
                        ## theme 
                        p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour="white",size = 2) +
                            annotate("rect", xmin= time.bf, xmax= Inf,
                                     ymin=-Inf, ymax=Inf, alpha = .3) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0)))

                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type))

                        ## band
                        p <- p + geom_ribbon(data = data.band,
                                        aes(ymin = co5, ymax = co95, x=time),
                                        alpha = 0.15, fill = "steelblue")

                        set.limits = c("tr","co.band","ct")
                        set.labels = c("Treated", "Controls (5-95% Quantiles)",
                                       "Estimated Y(0)")
                        set.colors = c("red","#4682B480","steelblue")
                        set.linetypes = c("solid","solid","longdash")
                        set.linewidth = c(line.width[1],4,line.width[1])

                        p <- p + scale_colour_manual(limits = set.limits,
                                                     labels = set.labels,
                                                     values =set.colors) +
                            scale_linetype_manual(limits = set.limits,
                                                  labels = set.labels,
                                                  values = set.linetypes) +
                            scale_size_manual(limits = set.limits,
                                              labels = set.labels,
                                              values = set.linewidth) +
                            guides(linetype = guide_legend(title=NULL, ncol=3),
                                   colour = guide_legend(title=NULL, ncol=3),
                                   size = guide_legend(title=NULL, ncol=3)) 

                        if (!is.numeric(time.label)) {
                            p <- p + 
                            scale_x_continuous(expand = c(0, 0), breaks = T.b, labels = time.label[T.b])
                        }
                    
                    } else if (raw == "all") { ## plot all the raw data
                    
                        data <- cbind.data.frame("time" = rep(time[show],(2 + Nco)),
                                                 "outcome" = c(tr.info[show],
                                                               ct.info[show],
                                                               c(Y.co[show,])),
                                                 "type" = c(rep("tr",nT),
                                                            rep("ct",nT),
                                                            rep("raw.co",(Nco * nT))),
                                                 "id" = c(rep("tr",nT),
                                                          rep("ct",nT),
                                                          rep(c(x$id.co), each = nT)))
                    
                        ## theme
                        p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour="white",size = 2) +
                            annotate("rect", xmin= time.bf, xmax= Inf,
                                     ymin=-Inf, ymax=Inf, alpha = .3) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0)))
                    
                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type,
                                               group = id))

                        ## legend
                        set.limits = c("tr","raw.co","ct")
                        set.labels = c("Treated","Controls","Estimated Y(0)")
                        set.colors = c("red","#99999950","steelblue")
                        set.linetypes = c("solid","solid","longdash")
                        set.linewidth = c(line.width[1],line.width[2],line.width[1])
                    
                        p <- p + scale_colour_manual(limits = set.limits,
                                                     labels = set.labels,
                                                     values =set.colors) +
                            scale_linetype_manual(limits = set.limits,
                                                  labels = set.labels,
                                                  values = set.linetypes) +
                            scale_size_manual(limits = set.limits,
                                              labels = set.labels,
                                              values = set.linewidth) +
                            guides(linetype = guide_legend(title=NULL, ncol=3),
                                   colour = guide_legend(title=NULL, ncol=3),
                                   size = guide_legend(title=NULL, ncol=3)) 

                        if (!is.numeric(time.label)) {
                            p <- p + 
                            scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                        }                        
                     
                    } 
                
                } else { # begin multiple treated unit case
                    maintext <- "Treated and Counterfactual Averages"
                    if (raw == "none") {
                        data <- cbind.data.frame("time" = rep(time[show],2),
                                                 "outcome" = c(Yb[show,1],
                                                               Yb[show,2]),
                                                 "type" = c(rep("tr",nT),
                                                            rep("co",nT))) 
                        ## theme 
                        p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour="white",size = 2) +
                            annotate("rect", xmin= time.bf, xmax= Inf,
                                     ymin=-Inf, ymax=Inf, alpha = .3) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0)))
                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type))

                        ## legend
                        set.limits = c("tr","co")
                        set.labels = c("Treated Average",
                                       "Estimated Y(0) Average")
                        set.colors = c("red","steelblue")
                        set.linetypes = c("solid","longdash")
                        set.linewidth = rep(line.width[1],2)
                        p <- p + scale_colour_manual(limits = set.limits,
                                                     labels = set.labels,
                                                     values =set.colors) +
                            scale_linetype_manual(limits = set.limits,
                                                  labels = set.labels,
                                                  values = set.linetypes) +
                            scale_size_manual(limits = set.limits,
                                              labels = set.labels,
                                              values = set.linewidth) +
                            guides(linetype = guide_legend(title=NULL, ncol=2),
                                   colour = guide_legend(title=NULL, ncol=2),
                                   size = guide_legend(title=NULL, ncol=2)) 

                        if (!is.numeric(time.label)) {
                            p <- p + 
                            scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                        }
                    
                    } else if  (raw == "band") {
                    
                        Y.tr.90 <- t(apply(Y.tr, 1, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                        Y.co.90 <- t(apply(Y.co, 1, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                    
                        data <- cbind.data.frame("time" = rep(time[show],2),
                                                 "outcome" = c(Yb[show,1],
                                                               Yb[show,2]),
                                                 "type" = c(rep("tr",nT),
                                                            rep("co",nT)))

                        data.band <- cbind.data.frame(time, Y.tr.90, Y.co.90)[show,]
                        colnames(data.band) <- c("time","tr5","tr95","co5","co95")
                    
                        ## theme 
                        p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour="white",size = 2) +
                            annotate("rect", xmin= time.bf, xmax= Inf,
                                     ymin=-Inf, ymax=Inf, alpha = .3) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0)))
                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type))
                        ## band
                        p <- p + geom_ribbon(data = data.band,
                                             aes(ymin = tr5, ymax = tr95, x=time),
                                             alpha = 0.15, fill = "red") +
                            geom_ribbon(data = data.band,
                                        aes(ymin = co5, ymax = co95, x=time),
                                        alpha = 0.15, fill = "steelblue")

                        set.limits = c("tr","co","tr.band","co.band")
                        set.labels = c("Treated Average",
                                       "Estimated Y(0) Average",
                                       "Treated 5-95% Quantiles",
                                       "Controls 5-95% Quantiles")
                        set.colors = c("red","steelblue","#FF000030","#4682B480")
                        set.linetypes = c("solid","longdash","solid","solid")
                        set.linewidth = c(rep(line.width[1],2),4,4)

                        p <- p + scale_colour_manual(limits = set.limits,
                                                     labels = set.labels,
                                                     values =set.colors) +
                            scale_linetype_manual(limits = set.limits,
                                                  labels = set.labels,
                                                  values = set.linetypes) +
                            scale_size_manual(limits = set.limits,
                                              labels = set.labels,
                                              values = set.linewidth) +
                            guides(linetype = guide_legend(title=NULL, ncol=2),
                                   colour = guide_legend(title=NULL, ncol=2),
                                   size = guide_legend(title=NULL, ncol=2)) 

                        if (!is.numeric(time.label)) {
                            p <- p + 
                            scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                        }
                    
                    } else if (raw == "all") { ## plot all the raw data
                    
                        data <- cbind.data.frame("time" = rep(time[show],(2 + N)),
                                                 "outcome" = c(Yb[show,1],
                                                               Yb[show,2],
                                                               c(Y.tr[show,]),
                                                               c(Y.co[show,])),
                                                 "type" = c(rep("tr",nT),
                                                            rep("co",nT),
                                                            rep("raw.tr",(Ntr * nT)),
                                                            rep("raw.co",(Nco * nT))),
                                                 "id" = c(rep("tr",nT),
                                                          rep("co",nT),
                                                          rep(c(x$id.tr,x$id.co),
                                                              each = nT))) 
                        ## theme
                        p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour="white",size = 2) +
                            annotate("rect", xmin= time.bf, xmax= Inf,
                                     ymin=-Inf, ymax=Inf, alpha = .3) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0))) 
                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type,
                                               group = id))
                        ## legend
                        set.limits = c("tr","co","raw.tr","raw.co")
                        set.labels = c("Treated Average",
                                       "Estimated Y(0) Average",
                                       "Treated Raw Data",
                                       "Controls Raw Data")
                        set.colors = c("red","steelblue","#FC8D6280","#99999950")
                        set.linetypes = c("solid","longdash","solid","solid")
                        set.linewidth = rep(line.width,each=2)
                    
                        p <- p + scale_colour_manual(limits = set.limits,
                                                     labels = set.labels,
                                                     values =set.colors) +
                            scale_linetype_manual(limits = set.limits,
                                                  labels = set.labels,
                                                  values = set.linetypes) +
                            scale_size_manual(limits = set.limits,
                                              labels = set.labels,
                                              values = set.linewidth) +
                            guides(linetype = guide_legend(title=NULL, ncol=2),
                                   colour = guide_legend(title=NULL, ncol=2),
                                   size = guide_legend(title=NULL, ncol=2)) 

                        if (!is.numeric(time.label)) {
                            p <- p + 
                            scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                        }
                    }

                } # end multiple treated unit case

                ## title
                if (is.null(main) == TRUE) {
                    p <- p + ggtitle(maintext)
                } else if (main!="") {
                    p <- p + ggtitle(main)
                }
            
                ## ylim
                if (is.null(ylim) == FALSE) {
                    p <- p + coord_cartesian(ylim = ylim)
                }
                suppressWarnings(print(p))
            }
        } else {
            maintext <- "Treated and Counterfactual Averages"

            ## axes labels
            if (is.null(xlab)==TRUE) {
                xlab <- paste("Time relative to Treatment")
            } else if (xlab == "") {
                xlab <- NULL
            }
            if (is.null(ylab)==TRUE) {
                ylab <- x$Yname
            } else if (ylab == "") {
                ylab <- NULL
            }
            
            xx <- ct.adjsut(x$Y.tr, x$Y.ct, x$T0)

            time <- xx$timeline
            Yb <- xx$Yb
            Y.tr.aug <- xx$Y.tr.aug
            ## Y.ct.aug <- xx$Y.ct.aug
            time.bf <- 0 ## before treatment

            if (!is.null(xlim)) {
                show <- which(time>=xlim[1]& time<=xlim[2])
            } else {
                show <- 1:length(time)
            }
            nT <- length(show)

            if (raw == "none") {
                data <- cbind.data.frame("time" = rep(time[show],2),
                                         "outcome" = c(Yb[show,1],
                                                       Yb[show,2]),
                                         "type" = c(rep("tr",nT),
                                                    rep("co",nT))) 
                ## theme 
                p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                    geom_vline(xintercept=time.bf,colour="white",size = 2) +
                    annotate("rect", xmin= time.bf, xmax= Inf,
                                ymin=-Inf, ymax=Inf, alpha = .3) +
                    theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0)))
                ## main
                p <- p + geom_line(aes(time, outcome,
                                       colour = type,
                                       size = type,
                                       linetype = type))

                ## legend
                set.limits = c("tr","co")
                set.labels = c("Treated Average",
                               "Estimated Y(0) Average")
                set.colors = c("red","steelblue")
                set.linetypes = c("solid","longdash")
                set.linewidth = rep(line.width[1],2)
                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                    scale_linetype_manual(limits = set.limits,
                                          labels = set.labels,
                                          values = set.linetypes) +
                    scale_size_manual(limits = set.limits,
                                      labels = set.labels,
                                      values = set.linewidth) +
                    guides(linetype = guide_legend(title=NULL, ncol=2),
                            colour = guide_legend(title=NULL, ncol=2),
                            size = guide_legend(title=NULL, ncol=2)) 
                    
            } else if  (raw == "band") {
                    
                Y.tr.90 <- t(apply(Y.tr.aug, 1, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                ## Y.co.90 <- t(apply(Y.co, 1, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                    
                data <- cbind.data.frame("time" = rep(time[show],2),
                                         "outcome" = c(Yb[show,1],
                                                       Yb[show,2]),
                                         "type" = c(rep("tr",nT),
                                                    rep("co",nT)))

                data.band <- cbind.data.frame(time, Y.tr.90)[show,]
                colnames(data.band) <- c("time","tr5","tr95")
                    
                ## theme 
                p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                    geom_vline(xintercept=time.bf,colour="white",size = 2) +
                    annotate("rect", xmin= time.bf, xmax= Inf,
                             ymin=-Inf, ymax=Inf, alpha = .3) +
                    theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0)))
                ## main
                p <- p + geom_line(aes(time, outcome,
                                       colour = type,
                                       size = type,
                                       linetype = type))
                ## band
                p <- p + geom_ribbon(data = data.band,
                                     aes(ymin = tr5, ymax = tr95, x=time),
                                         alpha = 0.15, fill = "red")

                set.limits = c("tr","co","tr.band")
                set.labels = c("Treated Average",
                               "Estimated Y(0) Average",
                                "Treated 5-95% Quantiles")
                set.colors = c("red","steelblue","#FF000030")
                set.linetypes = c("solid","longdash","solid")
                set.linewidth = c(rep(line.width[1],2),4)

                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                    scale_linetype_manual(limits = set.limits,
                                          labels = set.labels,
                                          values = set.linetypes) +
                    scale_size_manual(limits = set.limits,
                                      labels = set.labels,
                                      values = set.linewidth) +
                    guides(linetype = guide_legend(title=NULL, ncol=2),
                           colour = guide_legend(title=NULL, ncol=2),
                           size = guide_legend(title=NULL, ncol=2)) 
                    
            } else if (raw == "all") { ## plot all the raw data
                    
                data <- cbind.data.frame("time" = rep(time[show],(2 + Ntr)),
                                         "outcome" = c(Yb[show,1],
                                                       Yb[show,2],
                                                       c(Y.tr.aug[show,])),
                                         "type" = c(rep("tr",nT),
                                                    rep("co",nT),
                                                    rep("raw.tr",(Ntr * nT))),
                                          "id" = c(rep("tr",nT),
                                                  rep("co",nT),
                                                  rep(c(x$id.tr),
                                                      each = nT))) 
                ## theme
                p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                    geom_vline(xintercept=time.bf,colour="white",size = 2) +
                    annotate("rect", xmin= time.bf, xmax= Inf,
                             ymin=-Inf, ymax=Inf, alpha = .3) +
                    theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0))) 
                ## main
                p <- p + geom_line(aes(time, outcome,
                                       colour = type,
                                       size = type,
                                       linetype = type,
                                       group = id))
                ## legend
                set.limits = c("tr","co","raw.tr")
                set.labels = c("Treated Average",
                               "Estimated Y(0) Average",
                               "Treated Raw Data")
                set.colors = c("red","steelblue","#FC8D6280")
                set.linetypes = c("solid","longdash","solid")
                set.linewidth = c(rep(line.width[1],2),line.width[2])
                    
                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                    scale_linetype_manual(limits = set.limits,
                                          labels = set.labels,
                                          values = set.linetypes) +
                    scale_size_manual(limits = set.limits,
                                      labels = set.labels,
                                      values = set.linewidth) +
                    guides(linetype = guide_legend(title=NULL, ncol=2),
                           colour = guide_legend(title=NULL, ncol=2),
                           size = guide_legend(title=NULL, ncol=2)) 
            }

            ## title
            if (is.null(main) == TRUE) {
                p <- p + ggtitle(maintext)
            } else if (main!="") {
                p <- p + ggtitle(main)
            }

            ## ylim
            if (is.null(ylim) == FALSE) {
                p <- p + coord_cartesian(ylim = ylim)
            }
            suppressWarnings(print(p))
        }

    } else if (type=="factors") {
        
        if (x$r.cv==0) {
            cat("No factors included in the model.\n")
        } else {
            ## axes labels
            if (is.null(xlab)==TRUE) {
                xlab <- x$index[2]
            } else if (xlab == "") {
                xlab <- NULL
            }
            if (is.null(ylab)==TRUE) {
                ylab <- "Estimate"
            } else if (ylab == "") {
                ylab <- NULL
            }
            ## title
            if (is.null(main) == TRUE) {
                main <- "Latent Factors"
            } else if (main=="") {
                main <- NULL
            }
            ## prapre data
            L.co<-x$lambda.co
            norm<-sqrt(diag(t(L.co)%*%L.co)/(x$N-x$Ntr))
            data <- cbind.data.frame("time" = rep(time[show],r),
                                     "factor" = c(F.hat[show,])*rep(norm,each=nT),
                                     "group" = as.factor(c(rep(1:r,each=nT))))
            ## theme
            p <- ggplot(data) + xlab(xlab) +  ylab(ylab) + ggtitle(main) +
                geom_hline(yintercept=0,colour="white",size = 2) +
                theme(legend.position = legend.pos,
                      axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                      plot.title = element_text(size=20,
                                                hjust = 0.5,
                                                face="bold",
                                                margin = margin(10, 0, 10, 0)))  
            ## main plot
            p <- p + geom_line(aes(time, factor,
                                   colour = group,
                                   group = group), size = 1.2)
            ## legend
            p <- p + guides(colour = guide_legend(title="Factor(s)", ncol=4)) 

            if (!is.numeric(time.label)) {
                p <- p + 
                    scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
            }
           
            ## ylim
            if (is.null(ylim) == FALSE) {
                p <- p + coord_cartesian(ylim = ylim)
            }
            suppressWarnings(print(p))
        }
        
    } else if (type=="loadings") {

        
        if (x$r.cv==0) {
            cat("No factors are included in the model.\n") 
        } else {
            ## number of loadings to be plotted
            if (is.null(nfactors)==TRUE) {
                nfactors<-min(x$r.cv,4) 
            } else if (nfactors>x$r.cv) {
                cat("Too many factors specified. ")
                nfactors<-min(x$r.cv,4) 
            }
            if (nfactors == 1) {
                cat("Loadings for the first factor are shown...\n")
            } else if (nfactors < x$r.cv) {
                cat(paste("Loadings for the first",nfactors,"factors are shown...\n"))
            }
            

            ## title
            if (is.null(main) == TRUE) {
                main <- "Factor Loadings"
            } else if (main=="") {
                main <- NULL
            }
            
            ## prepare data
            L.hat <- rbind(x$lambda.tr, x$lambda.co)
            Lname <- Llabel <- c()
            for (i in 1:r) {
                Lname<-c(Lname,paste("L",i,sep=""))
                Llabel<-c(Llabel,paste("Factor",i))
            }
            colnames(L.hat) <- Lname
            rownames(L.hat) <- c()
            data <- cbind.data.frame(L.hat,
                          "id"=c(x$id.tr, x$id.co),
                          "group"=as.factor(c(rep("Treated",Ntr),
                                              rep("Control",Nco))))

            if (nfactors == 1) {
                p <- ggplot(data, aes(x=group, y=L1, fill = group)) +
                    geom_boxplot(alpha = 0.7) +
                    coord_flip() + guides(fill=FALSE) +
                    xlab("") + ylab("Factor Loading")
            } else {
                
                if (x$Ntr < 5) {
                    my_dens <- function(data, mapping, ...) {
                        ggplot(data = data, mapping = mapping) +
                            geom_density(..., fill = "gray", alpha = 0.7, color = "gray50")
                    }
                    p <- ggpairs(data, mapping = aes(color = group),
                                 columns = 1:nfactors,
                                 columnLabels = Llabel[1:nfactors],
                                 diag = list(continuous = my_dens),
                                 title = main)
                } else {
                    my_dens <- function(data, mapping, ...) {
                        ggplot(data = data, mapping = mapping) +
                            geom_density(..., alpha = 0.7, color = NA)
                    }
                    p <- ggpairs(data, mapping = aes(color = group, fill = group),
                                 columns = 1:nfactors,
                                 columnLabels = Llabel[1:nfactors],
                                 diag = list(continuous = my_dens),
                                 title = main) +
                        theme(plot.title = element_text(hjust = 0.5))
                }
            } 
            suppressWarnings(print(p))
        }
           
    } else if (type=="missing") {
        
        if (is.null(xlab)==TRUE) {
            xlab <- x$index[2]
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab)==TRUE) {
            ylab <- x$index[1]
        } else if (ylab == "") {
            ylab <- NULL
        }
        if (is.null(main)==TRUE) {
            main <- "Treatment Status"
        } else if (main == "") {
            main <- NULL
        }

        m <- x$obs.missing
        if (!is.null(id)) {
            m <- as.matrix(m[show,which(colnames(m)%in%id)])
        } else {
            m <- as.matrix(m[show,])
            ## ylim <- colnames(m)
        }

        all <- unique(c(m))
        col <- col2 <- breaks <- label <- NULL
        if (0%in%all) {
            col <- c(col,"#FFFFFF")
            col2 <- c(col2, "0"=NA)
            breaks <- c(breaks,0)
            label <- c(label,"Missing")
        }
        if (1%in%all) {
            col <- c(col,"#B0C4DE")
            col2 <- c(col2, "1"=NA)
            breaks <- c(breaks,1)
            label <- c(label,"Controls")
        }
        if (2%in%all) {
            col <- c(col,"#4671D5")
            col2 <- c(col2, "2"=NA)
            breaks <- c(breaks,2)
            label <- c(label,"Treated (Pre)")
        }
        if (3%in%all) {
            col <- c(col,"#06266F")
            col2 <- c(col2, "3"=NA)
            breaks <- c(breaks,3)
            label <- c(label,"Treated (Post)")
        }
        if (4%in%all) {
            col <- c(col,"#A9A9A9")
            col2 <- c(col2, "4"="red")
            breaks <- c(breaks,4)
            label <- c(label,"Treated (Removed)")
        }


        T <- dim(m)[1]
        N <- dim(m)[2]
        units <- rep(rev(1:N), each = T)
        period <- rep(1:T, N)
        res <- c(m)
        data <- cbind.data.frame(units=units, period=period, res=res)
        data[,"res"] <- as.factor(data[,"res"])

        N.b <- 1:N
        
        p <- ggplot(data, aes(x = period, y = units,
                              fill = res), position = "identity") 
        p <- p + geom_tile(colour="gray90", size=0.1, stat="identity") 
  
        p <- p +
            labs(x = xlab, y = ylab, 
                ## fill = "Value", 
                title=main) +
            theme_bw() + 
            scale_fill_manual(NA, breaks = breaks, values = col, labels=label)

        if(4%in%all) {
            p <- p + geom_point(aes(colour=res),size=0.5)
            p <- p + scale_color_manual(NA, breaks=breaks,
                                        values=col2, labels=label)
        }

        p <- p +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(fill=NA,color="gray90", size=0.5, linetype="solid"),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_text(color="black", size=14),
              axis.title=element_text(size=12),
              axis.text.x = element_text(size = 8, angle = angle, hjust=x.h, vjust=x.v),
              axis.text.y = element_text(size = 8),
              plot.background = element_rect(fill = "grey90"),
              legend.background = element_rect(fill = "grey90"),
              legend.position = legend.pos,
              legend.title=element_blank(),
              plot.title = element_text(size=20,
                                        hjust = 0.5,
                                        face="bold",
                                        margin = margin(10, 0, 10, 0))) +
        scale_x_continuous(expand = c(0, 0), breaks = T.b, labels = time.label[T.b]) +
        scale_y_continuous(expand = c(0, 0), breaks = N.b, labels = rev(sort(id)))
        
        if(length(all)>=4) {
            p <- p + guides(fill=guide_legend(nrow=2,byrow=TRUE))
        }
        suppressWarnings(print(p))
    }    
}

## counterfactual adjust
ct.adjsut <- function (Y.tr,
                       Y.ct, 
                       T0) {
    T <- dim(Y.tr)[1]
    N <- dim(Y.tr)[2]
    ## T.end <- T - min(T0)
    ## T.start <-
    T.m <- matrix(rep(1:T,N),T,N) - matrix(rep(T0,each=T),T,N)
    timeline <- min(T.m):max(T.m)
    Y.tr.aug <- matrix(NA,length(timeline),N)
    Y.ct.aug <- matrix(NA,length(timeline),N)
    for(i in 1:N) {
        Y.tr.aug[which(timeline%in%T.m[,i]),i] <- Y.tr[,i]
        Y.ct.aug[which(timeline%in%T.m[,i]),i] <- Y.ct[,i]
    }
    Y.tr.bar <- apply(Y.tr.aug, 1, mean, na.rm=TRUE)
    Y.ct.bar <- apply(Y.ct.aug, 1, mean, na.rm=TRUE)
    Yb <- cbind(Y.tr.bar,Y.ct.bar)
    return(list(timeline=timeline,
                Y.tr.aug=Y.tr.aug,
                Y.ct.aug=Y.ct.aug,
                Yb=Yb))
 
}

###################################
## parametric bootstrap for ub data
###################################
res.vcov <- function(res, ## TT*Nboots
                     cov.ar = 1) {
    T <- dim(res)[1]
    I <- is.na(res)
    count <- matrix(NA,T,T)

    res[is.na(res)] <- 0
    vcov <- res%*%t(res)


    for (i in 1:T) {
        for (j in 1:T) {
            if (i > j) {
                count[i, j] <- count[j, i]
            } else {
                if ((j-i) <= cov.ar) {
                  II <- I[i,] + I[j,]
                  count[i, j] <- min( 1/sum(II==0), 1)
                } else {
                  count[i, j] <- 0
                }
            }
        }
    }
    vcov <- vcov*count
    return(vcov)
}




