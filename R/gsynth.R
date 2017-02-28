## Synthetic Control for Multiple Treated Units
## (Causal Inference with Interactive Fixed Effects Models)
## Version 1.03
## Author: Yiqing Xu, University of California, San Diego
## Date: 2016.7.28

## MAIN FUNCTION
## gsynth.formula()
## gsynth.default()

## DEPENDENT FUNCTIONS
## synth.core()
## synth.em()
## synth.em.cv()
## synth.boot()

## METHODS
## print.gsynth()
## plot.gsynth()

#####################################################################
## A Shell Function
#####################################################################

## generic function
gsynth <- function(x, ...) {
    UseMethod("gsynth")
}

## formula method
gsynth.formula <- function(formula, data, ...) {
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
    out <- gsynth.default(data = data, Y = Yname, D = Dname, X = Xname, ...)
    out$call <- match.call()
    out$formula <- formula
    print(out)
    return(out)
}

## default function
gsynth.default <- function(data, # a data frame (long-form)
                           Y, # outcome
                           D, # treatment 
                           X = NULL, # time-varying covariates
                           na.rm = FALSE, # remove missing values
                           index, # c(unit, time) indicators
                           force = "unit", # fixed effects demeaning
                           r = 0, # nubmer of factors
                           CV = TRUE, # cross-validation
                           EM = FALSE, # EM algorithm 
                           se = FALSE, # report uncertainties
                           nboots = 200, # number of bootstraps
                           inference = "parametric", # type of inference
                           cluster = NULL, #  clustering variable for block bootstrap
                           parallel = TRUE, # parallel computing
                           cores = NULL, # number of cores
                           tol = 0.001, # tolerance level
                           seed = NULL # set seed
                           ){  
    
    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------##  

    if (is.data.frame(data) == FALSE) {
        stop("Not a data frame.")
    }
    ## index
    if (length(index) != 2 | sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
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
    if (!force %in% c(0, 1,2,3)) {
        stop("\"force\" option misspecified; choose from c(\"none\", \"unit\", \"time\", \"two-way\").")
    } 
    ## r
    if (r[1] < 0) {
        stop("\"r\" option misspecified. The number of factors must be non-negative.")
    }
    ## CV
    if (CV == TRUE) {
        if (length(r) == 2 & r[1] > r[2]) {
            stop("\"r\" option misspecified.")
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
    if (is.logical(EM) == FALSE & is.numeric(EM)==FALSE) {
        stop("EM is not a logical flag.")
    }
    ## se
    if (is.logical(se) == FALSE & is.numeric(se)==FALSE) {
        stop("se is not a logical flag.")
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
    if (se == TRUE & nboots <= 0) {
        stop("\"nboots\" option misspecified. Try, for example, nboots = 200.")
    }
    ## cluster
    if (is.null(cluster) == FALSE) {
        if (!cluster %in% colnames(data)) {
            stop(paste("\"cluster\" option misspecified: varible", cluster, "not found."))
        }    
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
    ## seed
    if (is.null(seed)==FALSE) {
        if (is.numeric(seed)==FALSE) {
            stop("seed should be a number.")
        }
    }
    ## remove missing values
    if (is.logical(na.rm) == FALSE & is.numeric(na.rm)==FALSE) {
        stop("na.rm is not a logical flag.")
    } 
    if (na.rm == TRUE) {
        data <- na.omit(data[,c(Y, D, X, Z, FE)])
    } 
    ##-------------------------------##
    ## Parsing raw data
    ##-------------------------------##  

    ##store variable names
    Yname <- Y
    Dname <- D
    Xname <- X
    
    id <- index[1];
    time <- index[2];
    T <- length(unique(data[,time]))
    N <- length(unique(data[,id]))
    p <- length(Xname)
    
    ## check balanced panel
    if (var(table(data[,id])) + var(table(data[, time])) > 0) {
        stop("The panel is not balanced.")
    }
    
    ## check missingness
    if (sum(is.na(data[, Yname])) > 0) {
        stop(paste("Missing values in variable \"", Yname,"\".", sep = ""))
    }
    if (sum(is.na(data[, Dname])) > 0) {
        stop(paste("Missing values in variable \"", Dname,"\".", sep = ""))
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
    
    
    ##treatment indicator
    D<-matrix(data[,D],T,N)  
    
    ##outcome variable
    Y<-matrix(data[,Y],T,N)
    tr<-D[T,]==1     # cross-sectional: treated unit
    pre<-as.matrix(D[,which(tr==1)]==0) # a matrix indicating before treatment
    T0<-apply(pre,2,sum) 
    T0.min<-min(T0)
    id.co<-which(tr==0)
    
    ## time-varying covariates
    X <- array(0,dim=c(T,N,p))
    if (p > 0) {
        for (i in 1:p) {
            X[,,i] <- matrix(data[, Xname[i]], T, N)
            tot.var.unit <- sum(apply(X[, , i], 2, var))
            if (tot.var.unit == 0) {
                stop(paste("Variable \"", Xname[i],"\" is time-invariant.", sep = ""))   
            }
            if (force %in% c(2, 3)) {
                tot.var.time <- sum(apply(X[, , i], 1, var))
                if (tot.var.time == 0) {
                    stop(paste("Variable \"", Xname[i],"\" has no cross-sectional variation.", sep = ""))
                }
            } 
        } 
    }
    if (is.null(dim(X)[3])==TRUE) {
        p<-0
    } else {
        p<-dim(X)[3]
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

    ## for block bootstrap
    if (is.null(cluster) == TRUE) {
        cl.id <- NULL    
    } else {    
        if (cluster == id) {
            cl.id <- NULL  
        } else {
            cl.id <- matrix(data[,cluster],T,N)[1,]
            cl.id <- as.numeric(as.factor(cl.id))
        }
    }

    ##-------------------------------##
    ## Register clusters
    ##-------------------------------##
    
    if (se == TRUE & parallel==TRUE) {
        library(foreach)
        library(doParallel)
        library(abind)
        if (is.null(cores)==TRUE) {
            cores <- detectCores()
        }
        para.clusters <- makeCluster(cores)
        registerDoParallel(para.clusters)
        cat("Parallel computing ... ")
    }
    
    ##-------------------------------##
    ## run main program
    ##-------------------------------##

    if (se == FALSE) {
        if (EM == FALSE) { # the algorithm suggested in the paper 
            out<-synth.core(Y = Y, X = X, D = D,
                            r = r, r.end = r.end, force = force,
                            CV = CV, tol = tol, AR1 = AR1) 
        } else { # EM algorithm
            if (CV == FALSE) { 
                out<-synth.em(Y = Y, X = X, D = D,
                              r = r, force = force,
                              tol = tol, AR1 = AR1)
                
            } else { # cross-validation
                out<-synth.em.cv(Y = Y,X = X, D = D,
                                 r = r, r.end = r.end, force = force,
                                 tol = tol, AR1 = AR1) 
            } 
        } 
    } else  {
        if (is.null(seed) == FALSE) {
            set.seed(seed)
        }
        out<-synth.boot(Y = Y, X = X, D = D, EM = EM,
                        r = r, r.end = r.end, force = force,
                        CV = CV, tol = tol,
                        nboots = nboots, inference = inference,
                        parallel = parallel, cores = cores,
                        cl.id = cl.id, 
                        AR1 = AR1) 
    } 

    if (se == TRUE & parallel == TRUE) {
        stopCluster(para.clusters)
        ##closeAllConnections()
    }
    
    
    ##-------------------------------##
    ## storage
    ##-------------------------------## 
    
    iname<-unique(data[,id])
    tname<-unique(data[,time])
    if (AR1 == TRUE) {
        tname <- tname[-1]
    } 
    Xname.tmp <- Xname
    if (AR1 == TRUE) {
        Xname.tmp<-c(paste(Yname,"_lag",sep=""),Xname)
    }
    rownames(out$beta)<-Xname.tmp
    if (se == TRUE) {
        rownames(out$est.beta)<-Xname.tmp
    } 
    colnames(out$eff) <- iname[which(out$tr==1)]
    rownames(out$eff) <- tname
   
    output <- c(list(Y.dat = Y,
                     Y = Yname,
                     D = Dname,
                     X = Xname,
                     index = index,
                     id = iname,
                     time = tname,
                     id.tr = iname[which(out$tr==1)],
                     id.co = iname[which(out$tr==0)]),
                out,
                list(call = match.call()))
    class(output) <- "gsynth"
    return(output)
    
} ## Program GSynth ends 


###################################################################
## Core Function
###################################################################

synth.core<-function(Y, # Outcome variable, (T*N) matrix
                     X, # Explanatory variables:  (T*N*p) array
                     D, #  Indicator for treated unit (tr==1) 
                     r=0, # initial number of factors considered if CV==1
                     r.end,
                     force,
                     CV = 1, # cross-validation
                     tol, # tolerance level
                     AR1 = 0,
                     beta0 = NULL # starting value 
                     ){  
    
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  
    
    ## unit id and time
    T<-dim(Y)[1]
    N<-dim(Y)[2]
    if (is.null(X)==FALSE) {p<-dim(X)[3]} else {p<-0}
     
    ## treatement indicator
    tr<-D[T,]==1  ## cross-sectional: treated unit
    pre<-as.matrix(D[,which(tr==1)]==0) ## a (T*Ntr) matrix, time dimension: before treatment
    
    Ntr<-sum(tr)
    Nco<-N-Ntr
    T0<-apply(pre,2,sum) 
    T0.min<-min(T0)
    sameT0<-length(unique(T0))==1 ## treatment kicks in at the same time
    
    id<-1:N
    time<-1:T
    id.tr<-which(tr==1) ## treated id
    id.co<-which(tr==0)
    
    pre.v<-as.vector(pre)  ## vectorized "pre-treatment" indicator
    id.tr.pre.v<-rep(id,each=T)[which(pre.v==1)]  ## vectorized pre-treatment grouping variable for the treated
    time.pre<-split(rep(time,Ntr)[which(pre.v==1)],id.tr.pre.v) ## a list of pre-treatment periods

    ## parsing data
    Y.tr<-as.matrix(Y[,tr])
    Y.co<-Y[,!tr]
    if (p==0) {
        X.tr<-array(0,dim=c(T,Ntr,0))
        X.co<-array(0,dim=c(T,Nco,0)) 
    } else {
        X.tr<-array(NA,dim=c(T,Ntr,p))
        X.co<-array(NA,dim=c(T,Nco,p))
        for (j in 1:p) {
            X.tr[,,j]<-X[,tr,j]
            X.co[,,j]<-X[,!tr,j]
        }
    } 

    if (is.null(beta0) == TRUE ) {
        beta0 <- matrix(0, p, 1)
    }
    
    ##-------------------------------##
    ## Main Algorithm
    ##-------------------------------##
    
    if (CV == FALSE) { ## case: CV==0
        
        ## inter.fe on the control group
        est.co.best<-tryCatch(inter_fe(Y.co, X.co, r, force=force, beta0 = beta0),
                              error = function(e) {
                                  stop("Error: Multi-collinearity among covariates. Try dropping (some of) them.")
                              })
        r.cv<-r
        
    }  else if (CV == TRUE) { 
        
        ##-------------------------------##
        ## Cross-validation of r
        ##-------------------------------##
        
        ## starting r    
        if (r>(T0.min-1)) {
            cat("Warning: r is too big compared with T0; reset to 0.\n")
            r<-0
        }
        
        ## initial values
        cat("Cross-validating ...","\r")
        
        ## store all MSPE
        r.max<-min((T0.min-1),r.end)
        CV.out<-matrix(NA,(r.max-r+1),4)
        colnames(CV.out)<-c("r","sigma2","IC","MSPE")
        CV.out[,"r"]<-c(r:r.max)
        CV.out[,"MSPE"]<-1e7
        
        
        for (i in 1:dim(CV.out)[1]) { ## cross-validation loop starts 
            
            ## inter FE based on control, before & after 
            r<-CV.out[i,"r"]
            est.co<-tryCatch(inter_fe(Y = Y.co, X = X.co, r, force = force, beta0 = beta0),
                             error = function(e){ stop("Error: Multi-collinearity among covariates. Try dropping (some of) them.")})
            sigma2<-est.co$sigma2
            IC<-est.co$IC
            if (r!=0) {
                F.hat<-as.matrix(est.co$factor)
                if (force%in%c(1,3)) {F.hat<-cbind(F.hat,rep(1,T))} ## the last column is for alpha_i
            }      
            
            ## take out the effect of X (nothing to do with CV)
            U.tr<-Y.tr
            if (p>0) {for (j in 1:p) {
                          U.tr<-U.tr-X.tr[,,j]*est.co$beta[j]}
            }
            
            ## take out grant mean and time fixed effects (nothing to do with CV)
            if (force%in%c(1,2,3)) {
                U.tr<-U.tr-matrix(est.co$mu,T,Ntr)
            } ## grand mean
            if (force%in%c(2,3)) {
                U.tr<-U.tr-matrix(est.co$xi,T,Ntr,byrow=FALSE)
            } ## time fixed effects
            
            ## save for the moment       
            U.sav<-U.tr
            
            ## leave-one-out cross-validation
            sum.e2<-num.y<-0
            for (lv in unique(unlist(time.pre))){ ## leave one out using the pre-treatment period

                U.tr<-U.sav
                ## take out lv from U.tr.pre
                if (max(T0)==T0.min) {
                    U.lv<-as.matrix(U.tr[setdiff(c(1:T0.min),lv),])
                } else {
                    U.tr.pre.v<-as.vector(U.tr)[which(pre.v==1)]    ## pre-treatment residual in a vector
                    U.tr.pre<-split(U.tr.pre.v, id.tr.pre.v) ##  a list of pretreatment residuals
                    U.lv<-lapply(U.tr.pre,function(vec){return(vec[-lv])}) ## a list      
                }

                if (r==0) {            
                    if (force%in%c(1,3)) { ## take out unit fixed effect
                        if (max(T0)==T0.min) {
                            alpha.tr.lv<-colMeans(U.lv)
                            U.tr<-U.tr-matrix(alpha.tr.lv,T,Ntr,byrow=TRUE)
                        } else {
                            alpha.tr.lv<-sapply(U.lv,mean)
                            U.tr<-U.tr-matrix(alpha.tr.lv,T,Ntr,byrow=TRUE)
                        }
                    } 
                    e<-U.tr[which(time==lv),] ## that period
                } else {  ## case: r>0
                    ## take out the effect of factors
                    F.lv<-as.matrix(F.hat[which(time!=lv),])
                    if (max(T0)==T0.min) {
                        F.lv.pre<-F.hat[setdiff(c(1:T0.min),lv),]
                        lambda.lv<-solve(t(F.lv.pre)%*%F.lv.pre)%*%t(F.lv.pre)%*%U.lv
                    } else {
                        lambda.lv<-as.matrix(sapply(U.lv,function(vec){
                            F.lv.pre<-as.matrix(F.lv[1:length(vec),])
                            l.lv.tr<-solve(t(F.lv.pre)%*%F.lv.pre)%*%t(F.lv.pre)%*%vec 
                            return(l.lv.tr) ## a vector of each individual lambdas
                        }))
                    }
                    if (r>1|(r==1&force%in%c(1,3))) {
                        lambda.lv<-t(lambda.lv)
                    }
                    ## error term (the left-out period)
                    e<-U.tr[which(time==lv),] - c(F.hat[which(time==lv),]%*%t(lambda.lv)) 
                }
                if (sameT0 == FALSE) { # those who are actually not treated
                    e<-e[which(pre[which(time==lv),]==TRUE)]    
                }
                ## sum up
                sum.e2 <- sum.e2+t(e)%*%e
                num.y <- num.y+length(e)
                
            } ## end of leave-one-out
            
            MSPE<-sum.e2/num.y
            if ((min(CV.out[,"MSPE"]) - MSPE) > tol*min(CV.out[,"MSPE"])) {
                ## at least 5% improvement for MPSE
                est.co.best<-est.co  ## interFE result with the best r
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
        
        if (r>(T0.min-1)) {cat(" (r hits maximum)")}
        cat("\n\n r* = ",r.cv, sep="") 
        
        MSPE.best<-min(CV.out[,"MSPE"])  
        
    } ## End of Cross-Validation  
    
    
    ##-------------------------------##
    ## ATT and Counterfactuals 
    ##-------------------------------##
    
    ## variance of the error term
    sigma2<-est.co.best$sigma2
    IC<-est.co.best$IC

    
    ## ## take out the effect of X
    U.tr<-Y.tr
    if (p>0) {
        beta<-est.co.best$beta
        for (j in 1:p) {U.tr<-U.tr-X.tr[,,j]*beta[j]}
    }
   if (force%in%c(1,2,3)) {
        mu<-est.co.best$mu 
        U.tr<-U.tr-matrix(mu,T,Ntr) ## grand mean
        Y.fe.bar<-rep(mu,T)
    }
    if (force%in%c(2,3)) {
        xi<-est.co.best$xi ## a (T*1) matrix
        U.tr<-U.tr-matrix(c(xi),T,Ntr,byrow=FALSE)
        Y.fe.bar<-Y.fe.bar+xi
    }
    if (max(T0)==T0.min) {
        U.tr.pre<-as.matrix(U.tr[1:T0.min,])
    } else {
        U.tr.pre.v<-as.vector(U.tr)[which(pre.v==1)] # pre-treatment residual in a vector
        U.tr.pre<-split(U.tr.pre.v, id.tr.pre.v) ##  a list of pretreatment residuals
    }
     
    
    ## the error structure
    if (r.cv==0) {
        if (force%in%c(1,3)) { ## take out unit fixed effect
            if (max(T0)==T0.min) {
                alpha.tr<-colMeans(U.tr.pre)
                U.tr<-U.tr-matrix(alpha.tr,T,Ntr,byrow=TRUE)
            } else {
                alpha.tr<-sapply(U.tr.pre,mean)
                U.tr<-U.tr-matrix(alpha.tr,T,Ntr,byrow=TRUE)
            }
        }     
        eff<-U.tr  ## and that's it!
        
    } else { ## r.cv>0
        
        ## Factors
        F.hat<-as.matrix(est.co.best$factor)
        if (force%in%c(1,3)) {F.hat<-cbind(F.hat,rep(1,T))}
                                        # the last column is for alpha_i
         
        ## Lambda_tr (Ntr*r) or (Ntr*(r+1))
        if (max(T0)==T0.min) {
            F.hat.pre<-F.hat[1:T0.min,]
            lambda.tr<-solve(t(F.hat.pre)%*%F.hat.pre)%*%t(F.hat.pre)%*%U.tr.pre
        } else {
            lambda.tr<-as.matrix(sapply(U.tr.pre,function(vec){
                F.hat.pre<-as.matrix(F.hat[1:length(vec),])
                l.tr<-solve(t(F.hat.pre)%*%F.hat.pre)%*%t(F.hat.pre)%*%vec
                return(l.tr) ## a vector of each individual lambdas
            }))
        }
        if ((r.cv>1) | (r.cv==1 & force%in%c(1,3))) {
            lambda.tr<-t(lambda.tr) ## Ntr * r
        }

        ## predicting the treatment effect
        eff<-U.tr-F.hat%*%t(lambda.tr)   

        ## for storage
        if (force%in%c(1,3)) {
            alpha.tr<-lambda.tr[,(r.cv+1), drop = FALSE]
            lambda.tr<-lambda.tr[,1:r.cv, drop = FALSE] 
        }

    } ## end of r!=0 case


    ## AR1: calculate accumulative effect
    if (AR1 == TRUE) {
        rho<-est.co.best$beta[1]
        if (length(beta)>1) {
            beta<-beta[-1]
        } 
        eff.tmp<-eff*D[,id.tr]
        eff.acc<-matrix(0,T,Ntr)
        for (t in (T0.min+1):T) {
            for (i in 0:(t-T0.min-1)) {
                eff.acc[t,] <- eff.acc[t,]+eff.tmp[t-i,]*(rho^i)
            }
        }      
    }  

   
    ##-------------------------------##
    ## Summarize
    ##-------------------------------##  
    
    ## counterfactuals and averages
    Y.ct <- as.matrix(Y.tr-eff)
    Y.tr.bar <- rowMeans(Y.tr)
    Y.ct.bar <- rowMeans(Y.ct)
    Y.co.bar <- rowMeans(Y.co)

    ##Y.tr and Y.ct
    Y.bar <- cbind(Y.tr.bar,Y.ct.bar,Y.co.bar)
    colnames(Y.bar) <- c("Y.tr.bar","Y.ct.bar","Y.co.bar")
    
    ## ATT and average outcomes
    if (sameT0 == TRUE) { ## diff-in-diffs: same timing
        att <- rowMeans(eff)
    }  else { ## diff timing, centered the att
        eff.cnt <- Y.tr.center<-matrix(NA,T,Ntr)
        for (j in 1:Ntr) {
            eff.cnt[1:(T+T0.min-T0[j]), j] <- eff[(T0[j]-T0.min+1):T,j]  
            Y.tr.center[1:(T+T0.min-T0[j]),j] <- Y.tr[(T0[j]-T0.min+1):T,j]
        }
        att <- apply(eff.cnt, 1, mean, na.rm=TRUE)
        Y.tr.cnt <- apply(Y.tr.center, 1, mean, na.rm=TRUE)
        Y.ct.cnt <- Y.tr.cnt-att
    }
    att.avg<-sum(eff * (1 - pre))/sum(1 - pre)

    
    ##-------------------------------##
    ## Storage 
    ##-------------------------------##  

    ##control group residuals
    out<-list(
        ## main results
        Y.tr = Y.tr,
        Y.ct = Y.ct,
        Y.co = Y.co, 
        eff = eff,
        Y.bar = Y.bar,
        att = att,
        att.avg = att.avg,
        ## supporting
        force = force,
        DID = sameT0,
        T = T,
        N = N,
        p = p,
        Ntr = Ntr,
        Nco = Nco,
        T0 = T0,
        tr = tr,
        pre = pre,
        r.cv = r.cv,
        res.co = est.co.best$res,  
        sigma2 = sigma2,
        IC = IC,
        beta = beta,
        est.co = est.co.best
    )

    if (sameT0 == FALSE) {
        out<-c(out,list(eff.cnt = eff.cnt,
                        Y.tr.cnt = Y.tr.cnt,
                        Y.ct.cnt = Y.ct.cnt))
    }
    if (CV == 1) {
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
    } 
    if (force==1) {
        out<-c(out, list(mu = mu,
                         alpha.tr = alpha.tr,
                         alpha.co = est.co.best$alpha))
    } else if (force == 2) {
        out<-c(out,list(mu = mu,xi = xi))
    } else if (force == 3) {
        out<-c(out,list(mu = mu,
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
                   r = 0, # number of factors
                   force, # specifying fixed effects
                   tol=1e-5,
                   AR1 = 0
                   ){

    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  

    ## unit id and time
    T<-dim(Y)[1]
    N<-dim(Y)[2]
    if (is.null(X)==FALSE) {
        p<-dim(X)[3]
    } else {
        p<-0
    }

    ## treatement indicator
    tr<-D[T,]==1  ## cross-sectional: treated unit
    pre<-as.matrix(D[,which(tr==1)]==0) ## a (T*Ntr) matrix, time dimension: before treatment
    
    Ntr<-sum(tr)
    Nco<-N-Ntr
    T0<-apply(pre,2,sum) 
    T0.min<-min(T0)
    sameT0 <- length(unique(T0))==1 ## treatment kicks in at the same time
    
    id<-1:N
    time<-1:T
    id.tr<-which(tr==1) ## treated id
    id.co<-which(tr==0)
    
    pre.v<-as.vector(pre)  ## vectorized "pre-treatment" indicator
    id.tr.pre.v<-rep(id,each=T)[which(pre.v==1)]  ## vectorized pre-treatment grouping variable for the treated
    time.pre<-split(rep(time,Ntr)[which(pre.v==1)],id.tr.pre.v) ## a list of pre-treatment periods

    ## parsing data
    Y.tr<-as.matrix(Y[,tr])
    Y.co<-Y[,!tr]
    
    
    ##-------------------------------##
    ## Main Algorithm
    ##-------------------------------##

    init<-synth.core(Y = Y, X = X, D = D,
                     r = r, force = force,
                     CV = 0, tol = tol, AR1 = AR1)

    
    eff0 <- init$eff
    Y.ct <- init$Y.ct
    if (p > 0) {
        beta0 <- init$beta
    } else {
        beta0 <- matrix(0, 0, 0)
    }
    
    diff <- 100
    trace.diff <- c()
    niter <- 0

    while (niter <= 500 & diff > tol) {

        ## E step
        Y.e <- Y  # T*N
        Y.e[which(D==1)] <- Y.ct[which(pre==0)]  

        ## M step
        est<-inter_fe(Y.e, X, r, force=force, beta0 = beta0)
        Y.ct <- as.matrix(Y.e[,id.tr] - est$residuals[,id.tr]) # T * Ntr
        eff <- as.matrix(Y.tr - Y.ct)  # T * Ntr

        ## difference
        diff <- norm(eff0-eff, type="F")
        eff0 <- eff
        trace.diff <- c(trace.diff,diff)
        niter <- niter + 1
        
    }
    ## variance of the error term
    sigma2<-est$sigma2
    IC<-est$IC

    ##-------------------------------##
    ## Summarize
    ##-------------------------------##  
    
    ## counterfactuals and averages
    Y.tr.bar<-rowMeans(Y.tr)
    Y.ct.bar<-rowMeans(Y.ct)
    Y.co.bar<-rowMeans(Y.co)

    ##Y.tr and Y.ct
    Y.bar <- cbind(Y.tr.bar,Y.ct.bar,Y.co.bar)
    colnames(Y.bar) <- c("Y.tr.bar","Y.ct.bar","Y.co.bar")

    
    ## ATT and average outcomes
    if (sameT0==TRUE) { ## diff-in-diffs: same timing
        att<-rowMeans(eff)
        eff.cnt <- eff
        Y.tr.cnt <- Y.tr.bar
        Y.ct.cnt <- Y.ct.bar
    }  else { ## diff timing, centered the att
        eff.cnt<-Y.tr.center<-matrix(NA,T,Ntr)
        for (j in 1:Ntr) {
            eff.cnt[1:(T+T0.min-T0[j]),j]<-eff[(T0[j]-T0.min+1):T,j]  
            Y.tr.center[1:(T+T0.min-T0[j]),j]<-Y.tr[(T0[j]-T0.min+1):T,j]  
        }
        att<-apply(eff.cnt,1,mean,na.rm=TRUE)
        Y.tr.cnt<-apply(Y.tr.center,1,mean,na.rm=TRUE)
        Y.ct.cnt<-Y.tr.cnt-att
    }
    att.avg<-sum(eff*(1-pre))/sum(1-pre)

    ## fixed effects
    if (force%in%c(1,2,3)) {
        mu<-est$mu
    }
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
    }
    ## AR1: calculate accumulative effect
    if (AR1 == TRUE) {
        rho<-est$beta[1]
        if (length(beta)>1) {
            beta<-beta[-1]
        } 
        eff.tmp<-eff*D[,id.tr]
        eff.acc<-matrix(0,T,Ntr)
        for (t in (T0.min+1):T) {
            for (i in 0:(t-T0.min-1)) {
                eff.acc[t,]<-eff.acc[t,]+eff.tmp[t-i,]*(rho^i)
            }
        }      
    }  
    
    ##-------------------------------##
    ## Storage 
    ##-------------------------------##  
    
    out<-list(
        ## main results
        Y.tr=Y.tr,
        Y.ct=Y.ct,
        Y.co=Y.co, 
        eff=eff,
        Y.bar = Y.bar,
        att=att,
        att.avg=att.avg,
        ## supporting
        force=force,
        DID=sameT0,
        T=T,
        N=N,
        p=p,
        Ntr=Ntr,
        Nco=Nco,
        T0=T0,
        tr=tr,
        pre=pre,
        r.cv=r,
        res.co=est$residuals[,id.co],  ##control group residuals 
        beta = est$beta,
        niter = niter,
        sigma2 = sigma2,
        IC = IC
    )
    
    out<-c(out,list(eff.cnt=eff.cnt, ##
                    Y.tr.cnt=Y.tr.cnt, ##
                    Y.ct.cnt=Y.ct.cnt)) ##
    
    if (r > 0) {
        out<-c(out,list(factor=as.matrix(est$factor),
                        lambda.co=as.matrix(lambda.co),
                        lambda.tr=as.matrix(lambda.tr)
                        )) 
    } 
    if (force == 1) {
        out<-c(out,list(mu=mu,
                        alpha.tr=alpha.tr,
                        alpha.co=alpha.co))
    } else if (force == 2) {
        out<-c(out,list(mu = mu, xi = xi))
    } else if (force == 3) {
        out<-c(out,list(mu = mu,
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

synth.em.cv<-function(Y, # Outcome variable, (T*N) matrix
                      X, # Explanatory variables:  (T*N*p) array
                      D, # indicator for treated unit (tr==1) 
                      r = 0, # number of factors: starting point
                      r.end = 5, # end point
                      force, # specifying fixed effects
                      tol=1e-5, AR1 = 0){
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  

    ## unit id and time
    T<-dim(Y)[1]
    N<-dim(Y)[2]
    if (is.null(X)==FALSE) {
        p<-dim(X)[3]
    } else {
        p<-0
    }

    ## treatement indicator
    tr<-D[T,]==1  ## cross-sectional: treated unit
    pre<-as.matrix(D[,which(tr==1)]==0) ## a (T*Ntr) matrix, time dimension: before treatment
    
    Ntr<-sum(tr)
    Nco<-N-Ntr
    T0<-apply(pre,2,sum) 
    T0.min <- min(T0)
    sameT0 <- length(unique(T0))==1 ## treatment kicks in at the same time
    
    id<-1:N
    time<-1:T
    id.tr<-which(tr==1) ## treated id
    id.co<-which(tr==0)
    
    pre.v<-as.vector(pre)  ## vectorized "pre-treatment" indicator
    id.tr.pre.v<-rep(id,each=T)[which(pre.v==1)]  ## vectorized pre-treatment grouping variable for the treated
    time.pre<-split(rep(time,Ntr)[which(pre.v==1)],id.tr.pre.v) ## a list of pre-treatment periods

    ## parsing data
    Y.tr<-as.matrix(Y[,tr])
    Y.co<-Y[,!tr]

    ##-------------------------------##
    ## Cross-validation of r
    ##-------------------------------##
    
    ## starting r    
    if (r > (T0.min-1)) {
        cat("Warning: r is too big compared with T0; reset to 0.\n")
        r <- 0
    }
    
    ## store all MSPE
    r.max<-min((T0.min-1),r.end)
    CV.out<-matrix(NA,(r.max-r+1),4)
    colnames(CV.out)<-c("r","sigma2","IC","MSPE")
    CV.out[,"r"]<-c(r:r.max)
    CV.out[,"MSPE"]<-1e7
    cat("Cross-validating ...","\r")
    
    for (i in 1:dim(CV.out)[1]) { ## cross-validation loop starts
      
        r <- CV.out[i,"r"]
        est<-synth.em(Y = Y,X = X, D = D, r = r, force = force,
                      tol = tol, AR1 = AR1)
        sigma2<-est$sigma2
        IC<-est$IC
        
        ## leave-one-out cross-validation
        sum.e2<-num.y<-0
        for (lv in unique(unlist(time.pre))){ ## leave one out using the pre-treatment period

            D.cv <- D
            D.cv[which(time == lv), id.tr] <- 1 # set the left-out period to treated
            out<-synth.em(Y = Y,X = X, D = D.cv, r = r, force = force,
                          tol = tol, AR1 = AR1)

            e <- out$eff[which(time == lv),]
            if (sameT0 == FALSE) { # those who are actually not treated
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
    MSPE.best <- min(CV.out[,"MSPE"])  

    
    ##-------------------------------##
    ## Storage 
    ##-------------------------------##  
    
    out<-c(est.best, list(MSPE = MSPE.best, CV.out = CV.out)) 
    return(out) 
    
} ## EM cross validation ends



###############################################
## Inference 
###############################################

synth.boot<-function(Y,
                     X,
                     D, # input
                     EM, # EM algorithm
                     r=0, r.end,
                     force,
                     CV, # cross validation
                     nboots,
                     tol,
                     inference, # c("parametric","nonparametric")
                     AR1 = FALSE,
                     cl.id = NULL,  # a vector of cluster id,
                     parallel = TRUE,
                     cores = NULL){
    
    
    T<-dim(Y)[1]
    N<-dim(Y)[2]
    if (is.null(X)==FALSE) {
        p<-dim(X)[3]
    } else {
        p<-0
    }

    ## treatement indicator
    tr<-D[T,]==1    ## cross-sectional: indicating the treated units
    pre<-as.matrix(D[,which(tr==1)]==0)
                                         
    Ntr<-sum(tr)
    Nco<-N-Ntr
    T0<-apply(pre,2,sum) 
    T0.min<-min(T0)
    sameT0<-length(unique(T0))==1 ## treatment kicks in at the same time
    
    id<-1:N
    time<-1:T
    id.tr<-which(tr==1) ## treated id
    id.co<-which(tr==0)

    ## vectorized "pre-treatment" indicator
    pre.v<-as.vector(pre)
    ## vectorized pre-treatment grouping variable for the treated
    id.tr.pre.v<-rep(id,each=T)[which(pre.v==1)]
    ## create a list of pre-treatment periods
    time.pre<-split(rep(time,Ntr)[which(pre.v==1)],id.tr.pre.v) 
    
    ## estimation
    if (EM == FALSE) {
        out<-synth.core(Y = Y,X = X, D = D, r = r, r.end = r.end,
                        force = force, CV = CV, tol=tol, AR1 = AR1)
    } else { # the case with EM
        if (CV == FALSE) {
            out<-synth.em(Y = Y,X = X, D = D, r = r, force = force,
                          tol = tol, AR1 = AR1)
        } else {
            out<-synth.em.cv(Y = Y,X = X, D = D, r = r, r.end = r.end,
                             force = force, tol=tol, AR1 = AR1)
        }
    }
    ## output 
    eff<-out$eff
    att<-out$att
    att.avg<-out$att.avg
    if (p > 0) {
        beta<-out$beta
    } else {
        beta<-matrix(0,0,1)
    }
    error.co<-out$res.co ## error terms (T*Nco)
   
    Y.tr.bar=out$Y.bar[,1]
    Y.ct.bar=out$Y.bar[,2]
    Y.co.bar=out$Y.bar[,3]
    
    ## cluster structure
    if (is.null(cl.id)==FALSE) { 
        cl.tr <- unique(cl.id[id.tr])
        cl.co <- unique(cl.id[id.co])
        cl.all <- unique(cl.id)
        Ntr.cl <- length(cl.tr)  ## number of treated clusters
        Nco.cl <- length(cl.co)  ## number of control clusters
        N.cl <- length(cl.all)
        cl.list <- split(id,cl.id)
    }
    
    ## bootstrapped estimates
    eff.boot<-array(NA,dim=c(T,Ntr,nboots))  ## to store results
    att.boot<-matrix(NA,T,nboots)
    att.avg.boot<-matrix(NA,nboots,1)
    if (p>0) {
        beta.boot<-matrix(NA,p,nboots)
    }
    
    if (inference=="nonparametric") { ## nonparametric bootstrap

        cat("\rBootstrapping ...")
        if (EM == FALSE) {
            one.nonpara <- function(){
                if (is.null(cl.id)==FALSE) { ## cluster structure,
                    boot.cl<-c(sample(cl.tr,Ntr.cl,replace=TRUE),
                               sample(cl.co, Nco.cl, replace=TRUE))
                    ## bootstrapped clusters
                    boot.id<-unlist(cl.list[boot.cl])
                } else {
                    boot.id<-c(sample(id.tr,Ntr,replace=TRUE),
                               sample(id.co,Nco, replace=TRUE))
                }
                X.boot<-X[,boot.id,,drop=FALSE]
                boot<-synth.core(Y[,boot.id], X.boot,D[,boot.id],
                                 force = force, r = out$r.cv, CV=0,
                                 tol = tol, AR1 = AR1,
                                 beta0 = beta)
                return(boot)
            } 
        } else { # the case of EM
            one.nonpara <- function(){
                if (is.null(cl.id)==FALSE) { ## cluster structure,
                    boot.cl<-c(sample(cl.tr,Ntr.cl,replace=TRUE),
                               sample(cl.co, Nco.cl, replace=TRUE))
                    ## bootstrapped clusters
                    boot.id<-unlist(cl.list[boot.cl])
                } else {
                    boot.id<-c(sample(id.tr,Ntr,replace=TRUE),
                               sample(id.co,Nco, replace=TRUE))
                }
                X.boot<-X[,boot.id,,drop=FALSE]
                boot<-synth.em(Y[,boot.id], X.boot,D[,boot.id],
                               force = force, r = out$r.cv, 
                               tol = tol, AR1 = AR1,
                               beta0 = beta)
                return(boot)
            } 
        }
        ## computing
        if (parallel == TRUE) { 
            boot.out <- foreach(j=1:nboots, 
                                .inorder = FALSE,
                                .export = c("synth.core","synth.em"),
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

        if (EM == FALSE) { # the case without EM
            
            ## y fixed
            error.co<-out$res.co
            Y.fixed<-Y
            Y.fixed[,id.tr]<-as.matrix(out$Y.ct)
            Y.fixed[,id.co]<-Y.fixed[,id.co]-error.co

            draw.error <- function(){
                ## draw 1 prediction error at a time      
                if (is.null(cl.id)==TRUE) { 
                    fake.tr<-sample(id.co,1,replace=FALSE)
                    id.co.rest<-id.co[which(!id.co%in%fake.tr)]
                    ## resample control, to smooth CV prediction error
                    id.co.pseudo<-sample(id.co.rest,Nco,replace=TRUE)
                    
                } else { # clustered 
                    ## sample a unit
                    fake.tr<-sample(id.co,1,replace=FALSE)
                    cl.rest<-setdiff(cl.co, which(sapply(cl.list, function(vec)
                        fake.tr%in%vec)==TRUE))  
                    ## resample the rest of the clusters
                    cl.pseudo<-sample(cl.rest,Nco.cl,replace=TRUE)
                    ## resample control, to smooth CV prediction error
                    id.co.pseudo<-unlist(cl.list[cl.pseudo]) 
                }     
                id.pseudo<-c(rep(fake.tr,Ntr),id.co.pseudo)  ## Ntr + ...
                
                ## obtain the prediction eror
                D.pseudo<-D[,c(id.tr,id.co.pseudo)]  ## fake.tr + control left
                Y.pseudo<-Y[,id.pseudo]
                X.pseudo<-X[,id.pseudo,,drop=FALSE]

                ## output
                tryit <- tryCatch(synth.core(Y = Y.pseudo,
                                       X = X.pseudo,
                                       D = D.pseudo,
                                       force = force,
                                       r = out$r.cv,
                                       CV = 0,
                                       tol = tol,
                                       AR1 = AR1,
                                       beta0 = beta),
                                  error = function(e) {
                                      stop("Error: Some covariates do not have enough variation after resampling. Try dropping them.")
                                  }) 
                return(as.matrix(tryit$eff)) ## T * Ntr
                
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
                error.tr<-array(NA,dim=c(T,Ntr,nboots))
                for (j in 1:nboots) {
                    error.tr[,,j] <- draw.error()
                    if (j%%100==0) {
                        cat(".")
                    }
                }
            }
 
            
            one.boot <- function(){
                ## boostrap ID
                if (is.null(cl.id)==TRUE) {
                    id.boot<-c(id.tr,sample(id.co,Nco,replace=TRUE))
                } else { # cluster structure
                    boot.cl<-sample(cl.co,Nco.cl,replace=TRUE)
                    ## bootstrapped clusters for controls
                    id.boot.co<-unlist(cl.list[boot.cl])
                    id.boot<-c(id.tr,id.boot.co)
                }
                ## get the error for the treated and control
                error.tr.boot<-matrix(NA,T,Ntr)
                for (w in 1:Ntr) {
                    error.tr.boot[,w]<-error.tr[,w,sample(1:nboots,1,replace=TRUE)]
                }
                error.co.boot<-error.co[,sample(1:Nco,(length(id.boot)-Ntr),
                                                replace=TRUE)] 
                
                Y.boot<-Y.fixed[,id.boot]
                Y.boot[,1:Ntr]<- as.matrix(Y.boot[,1:Ntr] + error.tr.boot)
                ## new treated: conterfactual+effect+ (same) new error
                Y.boot[,(Ntr+1):length(id.boot)]<-
                    Y.boot[,(Ntr+1):length(id.boot)] + error.co.boot
                X.boot<-X[,id.boot,,drop=FALSE] 
                D.boot<-D[,id.boot] 
                
                ## re-estimate the model 
                boot <- tryCatch(synth.core(Y.boot,
                                            X.boot,
                                            D.boot,
                                            force = force,
                                            r = out$r.cv,
                                            CV = 0,
                                            tol = tol,
                                            AR1 = AR1,
                                            beta0 = beta),
                                 error = function(e) {
                                     stop("Error: Some covariates do not have enough variation after resampling. Try dropping them.")
                                 })
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
            error.co<-out$res.co
            Y.fixed <- Y
            Y.fixed[,id.tr]<-as.matrix(out$Y.ct)
            Y.fixed[,id.co]<-Y.fixed[,id.co] - out$res.co

            one.boot <- function(){

                ## sample errors
                if (is.null(cl.id)==TRUE) {
                    error.id <- sample(1:Nco, N, replace = TRUE)
                } else { # cluster structure
                    boot.cl <- sample(cl.co, Nco.cl*2, replace=TRUE)
                    ## bootstrapped clusters for controls
                    error.id <- unlist(cl.list[boot.cl])[1:N] 
                }
                
                ## produce the new outcome data
                Y.boot<-Y.fixed + error.co[,error.id]
                
                ## re-estimate the model
                boot<-synth.em(Y.boot, X, D, force=force, r=out$r.cv,
                               tol=tol, AR1 = AR1)

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
        cat("\rBootstrapping ...")
        if (parallel == TRUE) { 
            boot.out <- foreach(k=1:nboots,
                                .inorder = FALSE,
                                .export = c("synth.core","synth.em"),
                                .packages = c("gsynth")) %dopar% {
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
        
    }  ## end of parametric bootstrap
   
    ####################################
    ## Variance and CIs
    ####################################

    ## function to get two-sided p-values
    get.pvalue <- function(vec){
        a <- sum(vec >= 0)/nboots * 2
        b <- sum(vec <= 0)/nboots * 2
        return(as.numeric(min(a, b)))
    }
    
    ## ATT estimates
    CI.att <- t(apply(att.boot, 1, function(vec) quantile(vec,c(0.025,0.975))))
    se.att <- apply(att.boot, 1, sd)
    pvalue.att <- apply(att.boot, 1, get.pvalue)
    if (sameT0 == TRUE) {
        ntreated <- apply((1 - pre), 1, sum)
    } else {
        rawcount <- apply((1 - pre), 1, sum)
        ntreated <- c(rep(0, T0.min), rev(rawcount[(T0.min + 1): T]))
    }
    est.att <- cbind(att, se.att, CI.att, pvalue.att, ntreated)
    colnames(est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                           "p.value", "n.Treated")
    if (sameT0 == TRUE) {
        rownames(est.att) <- time
    } else {
        rownames(est.att) <- c(1:T) - min(T0)
    }
    
    ## average (over time) ATT
    CI.avg <- quantile(att.avg.boot, c(0.025,0.975))
    se.avg <- sd(att.avg.boot)
    pvalue.avg <- get.pvalue(att.avg.boot)
    est.avg <- t(as.matrix(c(att.avg, se.avg, CI.avg, pvalue.avg)))
    colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")
    
    ## individual effects
    if (inference == "parametric") {
        CI.ind <- apply(eff.boot,c(1,2),function(vec)
            quantile(vec,c(0.025,0.975))) ## 2*T*Ntr
        est.ind <- array(NA,dim=c(T, 5, Ntr)) ## eff, se, CI.lower, CI.upper
        est.ind[,1,] <- eff
        est.ind[,2,] <- apply(eff.boot,c(1,2),sd)
        est.ind[,3,] <- CI.ind[1,,]
        est.ind[,4,] <- CI.ind[2,,]
        est.ind[,5,] <- apply(eff.boot,c(1,2),get.pvalue)
    }

    
    ## regression coefficents
    if (p>0) {
        CI.beta<-t(apply(beta.boot, 1, function(vec)
            quantile(vec,c(0.025, 0.975))))
        se.beta<-apply(beta.boot, 1, sd)
        pvalue.beta <- apply(beta.boot, 1, get.pvalue)
        est.beta<-cbind(beta, se.beta, CI.beta, pvalue.beta)
        colnames(est.beta)<-c("beta", "S.E.", "CI.lower", "CI.upper", "p.value")
    }
    
  
    ##storage
    result<-list(inference = inference,
                 est.att = est.att,
                 est.avg = est.avg
                 )
    
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

print.gsynth <- function(out, # a gsynth object
                         ...) {
    cat("Call:\n")
    print(out$call, digits = 4)
    
    if (is.null(out$est.avg) == TRUE) { # no uncertainties
        cat("\nAverage Treatment Effect on the Treated:\n")
        print(out$att.avg, digits = 4)
        cat("\n   ~ by Period (including Pre-treatment Periods):\n")
        print(out$att, digits = 4)
        if (is.null(out$X) == FALSE) {
            cat("\nCoefficients for the Covariates:\n")
            print(out$beta, digits = 4)
        }
        cat("\nUncertainty estimates not available.\n")
    } else {
        cat("\nAverage Treatment Effect on the Treated:\n")
        print(out$est.avg, digits = 4)
        cat("\n   ~ by Period (including Pre-treatment Periods):\n")
        print(out$est.att, digits = 4)
        if (is.null(out$X) == FALSE) {
            cat("\nCoefficients for the Covariates:\n")
            print(out$est.beta, digits = 4)
        }
    }
    
}


##########
## Plot
##########

plot.gsynth <- function(out, # a gsynth object
                        type = "gap", # type of the plot
                        xlim = NULL, ylim = NULL, # axes limits
                        xlab = NULL, ylab = NULL, # axes labels
                        legendOff = FALSE,
                        raw = "band", # show raw data in "counterfactual" mode
                                        # ("none","band","all")
                        main = NULL, # whether to show the title
                        nfactors = NULL, # whose loadings to be plotted
                        id = NULL # individual plot
                        ){


    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------##  

    if (class(out)!="gsynth") {
        stop("Not a \"gsynth\" object.")
    }
    if (!type %in% c("gap","counterfactual","factors","loadings","raw")) {
        stop("\"type\" option misspecified.")
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
        if (is.numeric(ylim)==FALSE) {
            stop("Some element in \"ylim\" is not numeric.")
        } else {
            if (length(ylim)!=2) {
                stop("ylim must be of length 2.")
            }
        }
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
            stop("\"raw\" option misspecifed.") 
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
    
    ##-------------------------------##
    ## Plotting
    ##-------------------------------##  

    require(ggplot2)
    Y.tr <- out$Y.tr
    Y.co <- out$Y.co
    Y.ct <- out$Y.ct
    tb <- out$est.att
    Yb <- out$Y.bar[,1:2] ## treated average and counterfactual average
    tr <- out$tr
    pre <- out$pre
    T <- out$T
    T0 <- out$T0
    p <- out$p
    m <- out$m
    Ntr <- out$Ntr
    Nco <- out$Nco
    N <- out$N 
    force <- out$force
    F.hat <- out$factor
    L.tr <- out$lambda.tr
    if (!is.null(L.tr)) {
        r <- dim(L.tr)[2]
    } else {
        r <- 0
    }
    if (is.null(id)==TRUE) {
        id <- out$id.tr
    }

    ## parameters
    line.width <- c(1.2,0.5)
  
    ## type of plots
    if (type == "raw"| type == "counterfactual" | type == "factors" |
        out$DID == TRUE | length(id) == 1) {
        time <- out$time
        if (length(id) == 1) {
            time.bf <- time[T0[which(id == out$id.tr)]]
        } else {
            time.bf <- time[unique(T0)]
        }
       
    } else if (type == "gap")  { ## variable treatment timing
        time <- c(1:T) - min(T0)
        time.bf <- 0 ## before treatment
    }

    ## periods to show
    if (length(xlim) != 0) {
        show <- which(time>=xlim[1]& time<=xlim[2])
    } else {
        show <- 1:length(time)
    }
    nT <- length(show) 

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
            xlab <- out$index[2]
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab)==TRUE) {
            ylab <- out$Yname
        } else if (ylab == "") {
            ylab <- NULL
        }
        pst <- (1 - out$pre)
        for (i in 1:Ntr){
            pst[T0[i],i] <- 1 ## paint the period right before treatment
        }
        time.pst <- c(pst[show,] * time[show])
        time.pst <- time.pst[which(c(pst[show,])==1)]
        Y.tr.pst <- c(Y.tr[show,])[which(pst[show,]==1)]
        id.tr.pst <- matrix(rep(1:Ntr,each=T),T,Ntr,byrow=FALSE)[show,]
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
                  plot.title = element_text(size=20,
                                            hjust = 0.5,
                                            face="bold",
                                            margin = margin(10, 0, 10, 0)))

        
        
        if (out$DID==TRUE) {
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
        print(p)
        
    } else if (type == "gap") { 
        
        if (length(id) == 1 & !(id[1] %in% out$id.tr)) { ## error
            cat(paste(id,"not in the treatment group"))
        } else { ## no error

            ## axes labels
            if (is.null(xlab) == TRUE) {
                if (out$DID == TRUE) {
                    xlab <- out$index[2]
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
                maintext <- paste(out$index[1],"=",id) 
            }  else {
                maintext <- "Estimated Average Treatment Effect"
            } 
            
            ## contruct data for plotting
            if (is.null(out$est.att)==TRUE) { 
                cat("Uncertainty estimates not available.\n")
                if (length(id) == 1) { ## id specified
                    data <- cbind.data.frame(time, out$eff)[show,]
                    colnames(data) <- c("time","ATT")
                } else {
                    data <- cbind.data.frame(time, ATT = out$att)[show,] 
                } 
            } else {
                if (length(id) == 1) { ## id specified
                    id <- which(out$id.tr == id)
                    tb <- out$est.ind[,,id]
                    time.bf <- time[T0[id]] 
                    colnames(tb) <- c("ATT", "S.E.", "CI.lower", "CI.upper","p.value") 
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
            if (is.null(out$est.att)==FALSE) {
                p <- p + geom_ribbon(aes(x = time, ymin=CI.lower, ymax=CI.upper),alpha=0.2)
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
            print(p)
        }  ## end of "gap" (in case of no id error)
       
        
    } else if (type=="counterfactual") {

         ## axes labels
        if (is.null(xlab)==TRUE) {
            xlab <- out$index[2]
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab)==TRUE) {
            ylab <- out$Yname
        } else if (ylab == "") {
            ylab <- NULL
        } 

        if (length(id)==1 & !(id[1]%in%out$id.tr)) { ## error
            
            cat(paste(id,"not in the treatment group"))
            
        } else { ## one treated unit case

             
            if (length(id) == 1 | length(out$id.tr) == 1) { ## one treated unit

                if (is.null(id) == TRUE) {
                    id <- out$id.tr
                }
                maintext <- paste("Treated and Counterfactual (",id,")",sep="") 
                tr.info <- Y.tr[,which(id==out$id.tr)]
                ct.info <- Y.ct[,which(id==out$id.tr)] 
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
                    
                } else if  (raw == "band") {

                    Y.co.90 <- t(apply(Y.co, 1, quantile, prob=c(0.05,0.95))) 
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
                                                      rep(c(out$id.co), each = nT)))
                    
                    ## theme
                    p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                        geom_vline(xintercept=time.bf,colour="white",size = 2) +
                        annotate("rect", xmin= time.bf, xmax= Inf,
                                 ymin=-Inf, ymax=Inf, alpha = .3) +
                        theme(legend.position = legend.pos,
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
                    set.labels = c("Treated Avaerge",
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
                    
                    Y.tr.90 <- t(apply(Y.tr, 1, quantile, prob=c(0.05,0.95)))
                    Y.co.90 <- t(apply(Y.co, 1, quantile, prob=c(0.05,0.95)))
                    
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
                    set.labels = c("Treated Avaerge",
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
                                                      rep(c(out$id.tr,out$id.co),
                                                          each = nT))) 
                    ## theme
                    p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                        geom_vline(xintercept=time.bf,colour="white",size = 2) +
                        annotate("rect", xmin= time.bf, xmax= Inf,
                                 ymin=-Inf, ymax=Inf, alpha = .3) +
                        theme(legend.position = legend.pos,
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
                    set.labels = c("Treated Avaerge",
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
        }
        print(p)

    } else if (type=="factors") {
        
        if (out$r.cv==0) {
            cat("No factors included in the model.\n")
        } else {
            ## axes labels
            if (is.null(xlab)==TRUE) {
                xlab <- out$index[2]
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
            L.co<-out$lambda.co
            norm<-sqrt(diag(t(L.co)%*%L.co)/(out$N-out$Ntr))
            data <- cbind.data.frame("time" = rep(time[show],r),
                                     "factor" = c(F.hat[show,])*rep(norm,each=nT),
                                     "group" = as.factor(c(rep(1:r,each=nT))))
            ## theme
            p <- ggplot(data) + xlab(xlab) +  ylab(ylab) + ggtitle(main) +
                geom_hline(yintercept=0,colour="white",size = 2) +
                theme(legend.position = legend.pos,
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
           
            ## ylim
            if (is.null(ylim) == FALSE) {
                p <- p + coord_cartesian(ylim = ylim)
            }
            print(p)
        }
        
    } else if (type=="loadings") {

        
        if (out$r.cv==0) {
            cat("No factors are included in the model.\n") 
        } else {
            ## number of loadings to be plotted
            if (is.null(nfactors)==TRUE) {
                nfactors<-min(out$r.cv,4) 
            } else if (nfactors>out$r.cv) {
                cat("Too many factors specified. ")
                nfactors<-min(out$r.cv,4) 
            }
            if (nfactors == 1) {
                cat("Loadings for the first factor are shown...\n")
            } else if (nfactors < out$r.cv) {
                cat(paste("Loadings for the first",nfactors,"factors are shown...\n"))
            }
            

            ## title
            if (is.null(main) == TRUE) {
                main <- "Factor Loadings"
            } else if (main=="") {
                main <- NULL
            }
            
            ## prepare data
            L.hat <- rbind(out$lambda.tr, out$lambda.co)
            Lname <- Llabel <- c()
            for (i in 1:r) {
                Lname<-c(Lname,paste("L",i,sep=""))
                Llabel<-c(Llabel,paste("Factor",i))
            }
            colnames(L.hat) <- Lname
            rownames(L.hat) <- c()
            data <- cbind.data.frame(L.hat,
                          "id"=c(out$id.tr, out$id.co),
                          "group"=as.factor(c(rep("Treated",Ntr),
                                              rep("Control",Nco))))

            if (nfactors == 1) {
                p <- ggplot(data, aes(x=group, y=L1, fill = group)) +
                    geom_boxplot(alpha = 0.7) +
                    coord_flip() + guides(fill=FALSE) +
                    xlab("") + ylab("Factor Loading")
            } else {
                require(GGally)
                if (out$Ntr < 5) {
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
            print(p) 
        }
           
    }## fig: loadings

   
    
}
