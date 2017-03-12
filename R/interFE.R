# Interactive Fixed Effect Model
# Version 1.02
# Yiqing Xu (MIT), 2014.12.8

## generic function
interFE <- function(formula=NULL,
                    data, # a data frame
                    Y, # outcome variable
                    X, # covariates
                    index, # id and time indicators
                    r = 0, # number of factors
                    force = "none", # additived fixed effects
                    se = TRUE, # standard error
                    nboots = 500, # number of bootstrap runs
                    seed = NULL) {
    UseMethod("interFE")
}

## formula method
interFE.formula <- function(formula=NULL, data, # a data frame
                            Y, # outcome variable
                            X, # covariates
                            index, # id and time indicators
                            r = 0, # number of factors
                            force = "none", # additived fixed effects
                            se = TRUE, # standard error
                            nboots = 500, # number of bootstrap runs
                            seed = NULL) {
    ## parsing
    varnames <- all.vars(formula)
    Yname <- varnames[1]
    Xname <- varnames[2:length(varnames)]
    
    ## run the model
    out <- interFE.default(formula=NULL, data = data, Y = Yname, X = Xname, 
                           index, # id and time indicators
                           r, # number of factors
                           force, # additived fixed effects
                           se, # standard error
                           nboots, # number of bootstrap runs
                           seed)
    out$call <- match.call()
    out$formula <- formula
    print(out)
    return(out)
}


print.interFE <- function(x, # a gsynth object
                         ...) {
    cat("Call:\n")
    print(x$call, digits = 4)
    cat("\nEstimated Coefficients:\n")
    print(x$est.table, digits = 4) 
}


###################################
# panel interactive fixed effects
###################################

interFE.default <- function(formula=NULL, data, # a data frame
                            Y, # outcome variable
                            X, # covariates
                            index, # id and time indicators
                            r = 0, # number of factors
                            force = "none", # additived fixed effects
                            se = TRUE, # standard error
                            nboots = 500, # number of bootstrap runs
                            seed = NULL
                            ){ 
    
    ##-------------------------------#
    ## Parameters
    ##-------------------------------#  

    ## index
    if (length(index) != 2 | sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
    }
    if (force == "none") { # no additive fixed effects imposed
        force <- 0
    } else if (force == "unit") { # unit fixed-effect
        force <- 1
    } else if (force == "time") { # time fixed-effect
        force <- 2
    } else if (force == "two-way") { # two-way fixed-effect 
        force <- 3
    }
    if (!force %in% c(0, 1, 2, 3)) {
        stop("\"force\" option misspecified; choose from c(\"none\", \"unit\", \"time\", \"two-way\").")
    } 
  
    ##-------------------------------#
    ## Parsing raw data
    ##-------------------------------#

    ## store variable names
    yname <- Y
    xname <- X
    id <- index[1]
    time <- index[2] 

    ## dimensions
    T <- length(unique(data[,time]))
    N <- length(unique(data[,id]))
    p<-length(xname) 
    
    
    ## parse data
    Y <- matrix(data[,yname],T,N)
    
    ## check time-varying covariates
    if (p==0) {
        X<-c()
    } else {
        X<-array(0,dim=c(T, N, p))
        for (i in 1: p) {
            X[,,i] <- matrix(data[, xname[i]], T, N)
            tot.var.unit <- sum(apply(X[, , i], 2, var))
            if (tot.var.unit == 0) {
                stop(paste("Variable \"", xname[i],"\" is time-invariant.", sep = ""))   
            }
            if (force %in% c(2, 3)) {
                tot.var.time <- sum(apply(X[, , i], 1, var))
                if (tot.var.time == 0) {
                    stop(paste("Variable \"", xname[i],"\" has no cross-sectional variation.", sep = ""))
                }
            } 
        } 
    } 
  
    ##-------------------------------#
    ## Estimation
    ##-------------------------------# 

    ## estimates
    out<-inter_fe(Y = Y, X = X, r = r, beta0 = as.matrix(rep(0,p)),
                  force = force)
    beta<-as.matrix(out$beta)
    mu <- out$mu
    

    ##-------------------------------#
    ## Standard Errors
    ##-------------------------------#

    ## function to get two-sided p-values
    get.pvalue <- function(vec){
        a <- sum(vec >= 0)/nboots * 2
        b <- sum(vec <= 0)/nboots * 2
        return(as.numeric(min(a, b)))
    }
    
    if (se == TRUE) {
        if (is.null(seed) == FALSE) {
            set.seed(seed)
        }
        ## to store results
        est.boot <- matrix(NA,nboots,(p+1))
        cat("Bootstraping")
        for (i in 1:nboots) {
            smp<-sample(1:N, N , replace=TRUE)
            Y.boot<-Y[,smp]
            X.boot<-X[,smp,,drop=FALSE]
            inter.out <- inter_fe(Y=Y.boot, X=X.boot, r=r,
                                  force=force, beta0 = beta)
            est.boot[i,]<- c(c(inter.out$beta), inter.out$mu)
            if (i%%100==0) {cat(".")}
        }
        cat("\r")
        ## T*2: lower,upper
        CI<-t(apply(est.boot,2,function(vec)
            quantile(vec,c(0.025,0.975))))
        SE<-apply(est.boot,2,sd)
        pvalue <- apply(est.boot, 2, get.pvalue)
         
        ## estimate table
        est.table<-cbind(c(beta,mu), SE, CI, pvalue)
        colnames(est.table) <- c("Coef","S.E.","CI.lower","CI.upper", "p.value")
    } else {
        est.table <- as.matrix(c(beta,mu))
    }
    rownames(est.table) <- c(xname,"_const")
    
    ##-------------------------------#
    ## Storage
    ##-------------------------------# 
   
    out<-c(out, list(dat.Y = Y,
                     dat.X = X,
                     Y = yname,
                     X = xname,
                     index = c(id,time)))
    if (se == TRUE) {
        out <- c(out,list(est.table = est.table,
                          est.boot = est.boot # bootstrapped coef.
                          ))
    } else {
        out <- c(out, list(est.table = est.table))
    }
    class(out) <- "interFE"
    return(out)

}




