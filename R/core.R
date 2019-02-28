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

    Y0.co <- NULL
    if (0 %in% I.co) { ## initial fit 
        ## Y0.co <- as.matrix(Y0[, id.co])
        data.ini <- matrix(NA, Nco*TT, (p+3))
        data.ini[,1] <- c(Y.co)
        data.ini[,2] <- rep(1:Nco, each = TT)
        data.ini[,3] <- rep(1:TT, Nco)
        if (p > 0) {
          for (i in 1:p) {
              data.ini[, (3 + i)] <- c(X.co[, , i])
          }
        }
        
        initialOut <- try(initialFit(data.ini, force, which(c(I.co) == 1)), silent = TRUE)
        if('try-error' %in% class(initialOut)) {
            return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1), eff = matrix(NA, TT, Ntr)))
            ## stop("Error occurs. Please set a smaller value of factor number.")
        }

        Y0.co <- initialOut$Y0
        if (p > 0) {
            beta0 <- initialOut$beta0
        }
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
              est.co.best <- inter_fe_ub(Y.co, Y0.co, X.co, I.co, beta0, r, force = force, tol)
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
                est.co.best <- inter_fe(Y.co, X.co, 0, force = force, beta0 = beta0, tol) 
            } else {
                est.co.best <- inter_fe_ub(Y.co, Y0.co, X.co, I.co, beta0, 0, force = force, tol)
            }

        } else {
            r.old <- r ## save the minimal number of factors 
            
            cat("Cross-validating ...","\r")
            CV.out <- matrix(NA, (r.max - r.old + 1), 5)
            colnames(CV.out) <- c("r", "sigma2", "IC", "PC", "MSPE")
            CV.out[,"r"] <- c(r.old:r.max)
            CV.out[,"MSPE"] <- CV.out[,"PC"] <- 1e20
            r.pc <- est.co.pc.best <- NULL
        
            for (i in 1:dim(CV.out)[1]) { ## cross-validation loop starts 
  
                ## inter FE based on control, before & after 
                r <- CV.out[i, "r"]
                if (!0%in%I.co) {
                    est.co <- inter_fe(Y = Y.co, X = X.co, r,
                                       force = force, beta0 = beta0, tol)
                } else {
                    est.co <- inter_fe_ub(Y = Y.co, Y0 = Y0.co, X = X.co, I = I.co, 
                                          beta0, r, force = force, tol)
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
                    PC <- est.co$PC
                } else {
                    sigma2 <- est.co$sigma2 * (norm.para[1]^2)
                    IC <- est.co$IC - log(est.co$sigma2) + log(sigma2)
                    PC <- est.co$PC * (norm.para[1]^2)
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

                if (PC < min(CV.out[,"PC"])) {
                    r.pc <- r
                    est.co.pc.best <- est.co
                } 
                CV.out[i, 2:5] <- c(sigma2, IC, PC, MSPE)
                cat("\n r = ",r, "; sigma2 = ",
                    sprintf("%.5f",sigma2), "; IC = ",
                    sprintf("%.5f",IC), "; PC = ",
                    sprintf("%.5f",PC), "; MSPE = ",
                    sprintf("%.5f",MSPE), sep="")
            
            } ## end of while: search for r_star over

            MSPE.best <- min(CV.out[,"MSPE"])

            ## compare 
            if (r.cv > r.pc) {
                cat("\n\n Factor number selected via cross validation may be larger than the true number. Using the PC criterion.\n\n ")
                r.cv <- r.pc 
                est.co.best <- est.co.pc.best
            }
        
            if (r > (T0.min-1)) {cat(" (r hits maximum)")}
            cat("\n\n r* = ",r.cv, sep="")
            cat("\n\n") 
        
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
        PC <- est.co.best$PC
    } else {
        sigma2 <- est.co.best$sigma2 * (norm.para[1]^2)
        IC <- est.co.best$IC - log(est.co.best$sigma2) + log(sigma2) 
        PC <- est.co.best$PC * (norm.para[1]^2)      
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
                alpha.tr <- as.matrix(colMeans(U.tr.pre))
                U.tr <- U.tr - matrix(alpha.tr, TT, Ntr, byrow = TRUE)
            } else {
                alpha.tr <- as.matrix(sapply(U.tr.pre, mean))
                U.tr <- U.tr - matrix(alpha.tr, TT, Ntr, byrow = TRUE)
            }
        }     
        eff <- U.tr  ## and that's it!
    ## r.cv>0
    } else {  
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
                if('try-error' %in% class(test)) {
                    return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1), eff = matrix(NA, TT, Ntr)))
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
            alpha.tr <- as.matrix(lambda.tr[, (r.cv+1), drop = FALSE])
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
    eff.cnt <- NULL
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
                Y.tr.cnt <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
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
                Y.tr.cnt <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
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
        est.co.best$residuals[which(I.co == 0)] <- NA
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
        est.co.best$mu <- est.co.best$mu * norm.para[1]

        ## sigma2 IC PC have been adjusted before
        est.co.best$sigma2 <- est.co.best$sigma2 * (norm.para[1]^2)
        est.co.best$IC <- est.co.best$IC - log(est.co.best$sigma2) + log(sigma2) 
        est.co.best$PC <- est.co.best$PC * (norm.para[1]^2)  
        
        if (r.cv > 0) {
            est.co.best$lambda <- est.co.best$lambda * norm.para[1]
            lambda.tr <- lambda.tr * norm.para[1]
        }
        if (force%in%c(1, 3)) {
            est.co.best$alpha <- est.co.best$alpha * norm.para[1]
            alpha.tr <- alpha.tr * norm.para[1]
        }
        if (force%in%c(2,3)) {
            est.co.best$xi <- est.co.best$xi * norm.para[1]
            xi <- xi * norm.para[1]
        }
        res.co <- res.co * norm.para[1] 
        est.co.best$residuals <- est.co.best$residuals * norm.para[1]
        est.co.best$fit <- est.co.best$fit * norm.para[1]
        
        Y.tr <- Y.tr * norm.para[1] 
        Y.ct <- Y.ct * norm.para[1]
        Y.co <- Y.co * norm.para[1]
        eff <- eff * norm.para[1]
        Y.bar <- Y.bar * norm.para[1]
        att <- att * norm.para[1]
        att.avg <- att.avg * norm.para[1]
        if ( !is.null(eff.cnt) ) {
            eff.cnt <- eff.cnt * norm.para[1]
            Y.tr.cnt <- Y.tr.cnt * norm.para[1]
            Y.ct.cnt <- Y.ct.cnt * norm.para[1]
        }
    }

    T0<-apply(as.matrix(D[,which(tr == 1)] == 0), 2, sum) ## for plot

    
    ##-------------------------------##
    ## Storage 
    ##-------------------------------##  

    names(att) <- c(1:TT) - min(T0.ub)

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
        sameT0 = DID,
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
        PC = PC,
        beta = beta,
        est.co = est.co.best,
        mu = mu,
        validX = validX
    )

    out <- c(out,list(sigma2 = sigma2, res.co=res.co))
    

    if ( DID == FALSE ) {
        names(Y.ct.cnt) <- names(Y.tr.cnt) <- rownames(eff.cnt) <- c(1:TT) - min(T0.ub)
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
                   beta0 = NULL,
                   norm.para,
                   boot = 0
                   ) {
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

    D.tr <- as.matrix(D[,which(tr == 1)])
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

    if (is.null(beta0) == TRUE ) {
        beta0 <- matrix(0, p, 1)
    }

    II <- I == 1 & D == 0 ## indicator

    
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

    #init <- synth.core(Y = Y, X = X, D = D, I = I, Y0 = Y0, W = W,
    #                   r = r, force = force,
    #                   CV = 0, tol = tol, AR1 = AR1, beta0 = beta0, 
    #                   norm.para = NULL, boot = boot)
    
    ## throw out error: may occur during bootstrap
    #if(length(init) == 2 || length(init) == 3) {
    #    return(init)
    #}
    
    #eff0 <- init$eff
    #eff0[is.na(eff0)] <- 0
    #Y.ct <- init$Y.ct
    #Y.ct[is.na(Y.ct)] <- 0

    Y0 <- NULL
    ## initial fit 
    data.ini <- matrix(NA, N*TT, (p+3))
    data.ini[,1] <- c(Y)
    data.ini[,2] <- rep(1:N, each = TT)
    data.ini[,3] <- rep(1:TT, N)
    if (p > 0) {
      for (i in 1:p) {
          data.ini[, (3 + i)] <- c(X[, , i])
      }
    }

    initialOut <- try(initialFit(data.ini, force, which(c(II) == 1)), silent = TRUE)
    if('try-error' %in% class(initialOut)) {
        return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1), eff = matrix(NA, TT, Ntr)))
        ## stop("Error occurs. Please set a smaller value of factor number.")
    }

    Y0 <- initialOut$Y0
    beta0 <- initialOut$beta0
    
    #diff <- 100
    #trace.diff <- c()
    #niter <- 0

    #while (niter <= 500 & diff > tol) {

        ## E step
        #if (!0 %in% I) {
        #    Y.e <- Y  # T*N
        #    Y.e.tr.tmp <- Y.e[,id.tr]
        #    Y.e.tr.tmp[which(post == 1)] <- Y.ct[which(post == 1)]
        #    Y.e[,id.tr] <- Y.e.tr.tmp 
        #} else {
        #    Y.e <- Y
        #    if (niter == 0) {
        #        Y.e[which(I==0)] <- Y0[which(I==0)]
        #        Y.e.tr.tmp <- Y.e[,id.tr]
        #        Y.e.tr.tmp[which(post == 1)] <- Y.ct[which(post == 1)]
        #        Y.e[,id.tr] <- Y.e.tr.tmp
        #    } else {
        #        Y.e[which(I==0)] <- est$fit[which(I==0)]
        #    }
        #}
        
    ## E step
    YY <- Y
    ## YY[which(II == 0)] <- 0 
        

        ## M step
        #if (!0%in%I) {
            ## if (force!=0) {
        #        est <- inter_fe(Y.e, X, r, force=force, beta0 = beta0, tol)
            ## } else {
            ##     est<-inter_fe(Y.e, abind(I,X,along=3), r, force=0, beta0 = beta0)
            ## }
        #} else {
            ## if (force!=0) {
    est <- inter_fe_ub(YY, Y0, X, II, beta0, r, force=force, tol)
            ## } else {
            ##     est<-inter_fe_ub(Y.e, abind(I,X,along=3), I, r, force=0, beta0 = beta0)
            ## }
        #}
        ## Y.ct <- as.matrix(Y.e[,id.tr] - est$residuals[,id.tr]) # T * Ntr
    Y.ct <- as.matrix(est$fit[,id.tr])

    eff <- as.matrix(Y.tr - Y.ct)  # T * Ntr
    eff[which(I.tr==0)] <- 0
    #    diff <- norm(eff0-eff, type="F")

    #    eff0 <- eff

    #    trace.diff <- c(trace.diff,diff)
    #    niter <- niter + 1  
    #}

    ## PC criterion
    

    ## variance of the error term
    if (is.null(norm.para)) {
        sigma2<-est$sigma2   
        IC<-est$IC
        PC <- est$PC
    } else {
        sigma2<-est$sigma2*(norm.para[1]^2)
        IC <- est$IC-log(est$sigma2) + log(sigma2)  
        PC <- est$PC*(norm.para[1]^2)     
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
    eff.cnt <- NULL
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
                Y.tr.cnt <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
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
                Y.tr.cnt <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
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
        alpha.tr <- as.matrix(alpha[id.tr])
        alpha.co <- as.matrix(alpha[id.co])
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
                wgt.implied <- t(inv.tr%*%t(lambda.co))
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
    
    res <- est$residuals
    res.co <- as.matrix(est$residuals[,id.co])
    
    
    if (0%in%I) {
        eff[which(I.tr==0)] <- NA 
        Y.ct[which(I.tr==0)] <- NA
        Y.tr[which(I.tr==0)] <- NA
        res[which(II == 0)] <- NA
        Y.co[which(I.co==0)] <- NA
        res.co[which(I.co==0)] <- NA
        est$residuals[which(II == 0)] <- NA

    }

        ## final adjustment
    if (!is.null(norm.para)) {
        mu <- mu*norm.para[1]
        est$mu <- est$mu * norm.para[1]

        ## sigma2 IC PC have been adjusted before
        est$sigma2 <- est$sigma2 * (norm.para[1]^2)
        est$IC <- est$IC - log(est$sigma2) + log(sigma2) 
        est$PC <- est$PC * (norm.para[1]^2)


        if (r>0) {
            est$lambda <- est$lambda * norm.para[1]
            lambda.tr <- lambda.tr*norm.para[1]
            lambda.co <- lambda.co*norm.para[1]
        }
        if (force%in%c(1,3)) {
            est$alpha <- est$alpha * norm.para[1]
            alpha.tr <- alpha.tr*norm.para[1]
            alpha.co <- alpha.co*norm.para[1]
        }
        if (force%in%c(2,3)) {
            est$xi <- est$xi * norm.para[1]
            xi <- xi*norm.para[1]
        }

        res <- res * norm.para[1] 
        res.co <- res.co * norm.para[1] 
        est$residuals <- est$residuals * norm.para[1]
        est$fit <- est$fit * norm.para[1]


        Y.tr <- Y.tr*norm.para[1] 
        Y.ct <- Y.ct*norm.para[1]
        Y.co <- Y.co*norm.para[1]
        eff <- eff*norm.para[1]
        Y.bar <- Y.bar*norm.para[1]
        att <- att*norm.para[1]
        att.avg <- att.avg*norm.para[1]
        if ( !is.null(eff.cnt) ) {
            eff.cnt <- eff.cnt*norm.para[1]
            Y.tr.cnt <- Y.tr.cnt*norm.para[1]
            Y.ct.cnt <- Y.ct.cnt*norm.para[1]
        }
    }

    T0<-apply(as.matrix(D[,which(tr==1)]==0),2,sum) 

    
    ##-------------------------------##
    ## Storage 
    ##-------------------------------##  

    names(att) <- c(1:TT) - min(T0.ub)
    
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
        sameT0=DID,
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
        ## niter = niter,
        IC = IC,
        PC = PC,
        mu = mu,
        validX = est$validX
    )

    out <- c(out,list(sigma2 = sigma2, res = res, res.co = res.co))
    
    
    if (DID==FALSE) {
        names(Y.ct.cnt) <- names(Y.tr.cnt) <- rownames(eff.cnt) <- c(1:TT) - min(T0.ub)
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
                      beta0 = NULL,
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

    if (is.null(beta0) == TRUE ) {
        beta0 <- matrix(0, p, 1)
    }

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

        r.old <- r ## save the minimal number of factors 
            
        cat("Cross-validating ...","\r")
        CV.out<-matrix(NA,(r.max-r.old+1),5)
        colnames(CV.out)<-c("r","sigma2","IC","PC","MSPE")
        CV.out[,"r"]<-c(r.old:r.max)
        CV.out[,"MSPE"]<-CV.out[,"PC"]<-1e20

        for (i in 1:dim(CV.out)[1]) { ## cross-validation loop starts
      
            r <- CV.out[i,"r"]
            est<-synth.em(Y = Y, X = X, D = D, I = I, W = W, r = r, force = force,
                          tol = tol, AR1 = AR1, beta0 = beta0, norm.para = norm.para, boot = 0)
            sigma2<-est$sigma2
            IC<-est$IC
            PC <- est$PC
        
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
                                tol = tol, AR1 = AR1, beta0 = beta0, norm.para = norm.para, boot = 0)

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

            if (PC < min(CV.out[,"PC"])) {
                r.pc <- r
                est.pc.best <- est
            }
            CV.out[i,2:5]<-c(sigma2,IC,PC,MSPE)
            cat("\n r = ",r,"; sigma2 = ",
                sprintf("%.5f",sigma2),"; IC = ",
                sprintf("%.5f",IC),"; PC = ",
                sprintf("%.5f",PC),"; MSPE = ",
                sprintf("%.5f",MSPE),sep="") 
        } ## end of while: search for r_star over

        MSPE.best <- min(CV.out[,"MSPE"])
        PC.best <- min(CV.out[,"PC"])

        ## compare 
        if (r.cv > r.pc) {
            cat("\n\n Factor number selected via cross validation may be larger than the true number. Using the PC criterion.\n\n ")
            r.cv <- r.pc 
            est.best <- est.pc.best
        }

        if (r>(T0.min-1)) {
            cat(" (r hits maximum)")
        }
        cat("\n\n r* = ", r.cv, sep="") 
        cat("\n\n") 
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
                   k = 5, 
                   hasF = 1,
                   tol, # tolerance level
                   AR1 = 0,
                   beta0 = NULL,
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

    oci <- which(c(II) == 1)

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

    if (is.null(beta0) == TRUE ) {
        beta0 <- matrix(0, p, 1)
    }
    
    id <- 1:N
    time <- 1:TT
    id.tr <- which(tr == 1) ## treated id
    id.co <- which(tr == 0)

    ## parsing data
    Y.tr <- as.matrix(Y[,id.tr])
    Y.co <- as.matrix(Y[,id.co])


    Y0 <- NULL
    ## initial fit 
    data.ini <- matrix(NA, N*TT, (p+3))
    data.ini[,1] <- c(Y)
    data.ini[,2] <- rep(1:N, each = TT)
    data.ini[,3] <- rep(1:TT, N)
    if (p > 0) {
      for (i in 1:p) {
          data.ini[, (3 + i)] <- c(X[, , i])
      }
    }

    initialOut <- try(initialFit(data.ini, force, oci), silent = TRUE)
    if('try-error' %in% class(initialOut)) {
        return(list(att = rep(NA, TT), att.avg = NA, beta = matrix(NA, p, 1), eff = matrix(NA, TT, Ntr)))
        ## stop("Error occurs. Please set a smaller value of factor number.")
    }

    Y0 <- initialOut$Y0
    beta0 <- initialOut$beta0
    
    ##-------------------------------##
    ## Main Algorithm
    ##-------------------------------##

    validX <- 1 ## no multi-colinearity
    
    if (CV == FALSE) { ## case: CV==0 or no factor  
        ## matrix completion
        est.best <- inter_fe_mc(YY, Y0, X, II, beta0, hasF, lambda[1], force, tol) 

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

        ## tot.id <- which(c(II)==1) ## observed control data
        cv.count <- ceiling((sum(II)*sum(II))/sum(I))

        if (is.null(lambda) || length(lambda) == 1) {
            ## create the hyper-parameter sequence
            ## lambda.max <- log10(max(svd(Y)$d)*2/(N*TT-sum(II)))
            ## Y.l <- YY - Y0
            ## Y.l[which(II == 0)] <- 0

            lambda.max <- log10(max(svd(YY)$d))
            lambda <- rep(NA, nlambda)
            lambda.by <- 3/(nlambda - 2)
            for (i in 1:(nlambda - 1)) {
                lambda[i] <- 10^(lambda.max - (i - 1) * lambda.by)
            }
            lambda[nlambda] <- 0
        }
        
        ## store all MSPE
        CV.out <- matrix(NA, length(lambda), 2)
        colnames(CV.out) <- c("lambda", "MSPE")
        CV.out[,"lambda"] <- c(lambda)
        CV.out[,"MSPE"] <- 1e20

        ociCV <- matrix(NA, cv.count, k) ## store indicator
        rmCV <- matrix(NA, (length(oci) - cv.count), k) ## removed indicator
        Y0CV <- array(NA, dim = c(TT, N, k)) ## store initial Y0
        if (p > 0) {
            beta0CV <- array(NA, dim = c(p, 1, k)) 
        } else {
            beta0CV <- array(0, dim = c(1, 0, k)) ## store initial beta0
        }

        for (i in 1:k) {
            cv.n <- 0
            repeat{
                cv.n <- cv.n + 1
                cv.id <- sample(oci, as.integer(sum(II) - cv.count), replace = FALSE)
                II.cv <- II
                II.cv[cv.id] <- 0
                con1 <- sum(apply(II.cv, 1, sum) > 0) == TT
                con2 <- sum(apply(II.cv, 2, sum) > 0) == N
                if (con1 & con2) {
                    break
                }
                if (cv.n > 100) {
                    stop("Some units have too few pre-treatment observations. Try to remove them.")
                }
            }
            rmCV[,i] <- cv.id
            ocicv <- setdiff(oci, cv.id)
            ociCV[,i] <- ocicv

            initialOutCv <- initialFit(data = data.ini, force = force, oci = ocicv)
            Y0CV[,,i] <- initialOutCv$Y0
                
            if (p > 0) {
                beta0cv <- initialOutCv$beta0
                beta0CV[,,i] <- beta0cv
            }
        }


        for (i in 1:length(lambda)) {    
            ## k <- 5
            SSE <- 0
            for (ii in 1:k) {
                II.cv <- II
                II.cv[rmCV[,ii]] <- 0
                YY.cv <- YY
                YY.cv[rmCV[,ii]] <- 0
                est.cv.fit <- inter_fe_mc(YY.cv, as.matrix(Y0CV[,,ii]), X, II.cv, as.matrix(beta0CV[,,ii]), 1, lambda[i], force, tol)$fit
                SSE <- SSE + sum((YY[rmCV[,ii]]-est.cv.fit[rmCV[,ii]])^2)
            }
            MSPE <- SSE/(k*(sum(II) - cv.count))

            est.cv <- inter_fe_mc(YY, Y0, X, II, beta0, 1, lambda[i], force, tol) ## overall

            if(!is.null(norm.para)){
                MSPE <- MSPE*(norm.para[1]^2)
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

            cat("\n lambda = ",
            sprintf("%.5f",lambda[i]),"; MSPE = ",
            sprintf("%.5f",MSPE), sep="")

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
 
    ## effect
    Y.fit <- est.best$fit
    Y.ct <- as.matrix(Y.fit[,tr])
    if (0%in%I.tr) {
        Y.ct[which(I.tr==0)] <- 0 ## adjust    
    }
    eff <- Y.tr - Y.ct
    res <- est.best$residuals
    res.co <- as.matrix(res[,id.co])


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
    eff.cnt <- NULL
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
                Y.tr.cnt <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
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
                Y.tr.cnt <- rowSums(Y.tr.center * W.tr.center)/rowSums(W.tr.center)
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
        res.co[which(I.co==0)] <- NA
        est.best$residuals[which(II == 0)] <- NA
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
        est.best$mu <- est.best$mu * norm.para[1]

        if (force%in%c(1, 3)) {
            est.best$alpha <- est.best$alpha * norm.para[1]
        }
        if (force%in%c(2,3)) {
            est$xi <- est$xi * norm.para[1]
            xi <- xi * norm.para[1]
        }

        res <- res * norm.para[1]
        res.co <- res.co * norm.para[1]
        est.best$residuals <- est.best$residuals * norm.para[1]

        Y.tr <- Y.tr * norm.para[1] 
        Y.ct <- Y.ct * norm.para[1]
        Y.co <- Y.co * norm.para[1]
        eff <- eff * norm.para[1]
        Y.bar <- Y.bar * norm.para[1]
        att <- att * norm.para[1]
        att.avg <- att.avg * norm.para[1]
        if ( !is.null(eff.cnt) ) {
            eff.cnt <- eff.cnt * norm.para[1]
            Y.tr.cnt <- Y.tr.cnt * norm.para[1]
            Y.ct.cnt <- Y.ct.cnt * norm.para[1]
        }
    }

    T0 <- apply(as.matrix(D[,which(tr == 1)] == 0), 2, sum) ## for plot

    
    ##-------------------------------##
    ## Storage 
    ##-------------------------------##
    names(att) <- c(1:TT) - min(T0.ub)  

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
        sameT0 = DID,
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

    out <- c(out,list(res = res, res.co = res.co))
    

    if ( DID == FALSE ) {
        names(Y.ct.cnt) <- names(Y.tr.cnt) <- rownames(eff.cnt) <- c(1:TT) - min(T0.ub)
        out<-c(out,list(eff.cnt = eff.cnt,
                        Y.tr.cnt = Y.tr.cnt,
                        Y.ct.cnt = Y.ct.cnt))
    }
    if (CV) {
        out<-c(out, list(MSPE = MSPE.best,
                         CV.out = CV.out))
    } 
    if (force==1) {
        out<-c(out, list(alpha.tr = as.matrix(est.best$alpha[id.tr,]), alpha.co = as.matrix(est.best$alpha[id.co,])))
    } else if (force == 2) {
        out<-c(out,list(xi = xi))
    } else if (force == 3) {
        out<-c(out,list(alpha.tr = as.matrix(est.best$alpha[id.tr,]), alpha.co = as.matrix(est.best$alpha[id.co,]),
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
                     k,
                     nboots,
                     tol,
                     inference, ## c("parametric","nonparametric")
                     cov.ar=1, 
                     AR1 = FALSE,
                     beta0 = NULL,
                     norm.para,
                     parallel = TRUE,
                     conf.lvl = 0.95,
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
                            AR1 = AR1, beta0 = beta0, norm.para= norm.para, boot = 0)
        } else { # the case with EM
            if (CV == FALSE) {
                out<-synth.em(Y = Y,X = X, D = D, I=I, W=W, r = r, force = force,
                              tol = tol, AR1 = AR1, beta0 = beta0, norm.para = norm.para, boot = 0)
            } else {
                out<-synth.em.cv(Y = Y, X = X, D = D, I=I, W=W, r = r, r.end = r.end,
                                 force = force, tol=tol,
                                 AR1 = AR1, beta0 = beta0, norm.para = norm.para)
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
                      force = force, tol=tol, CV = CV, k = k, 
                      AR1 = AR1, beta0 = beta0, norm.para= norm.para)
    }


    ## output
    validX <- out$validX
    eff<-out$eff
    att<-out$att
    att.avg<-out$att.avg
    DID <- out$sameT0

    if (p > 0) {
        beta <- out$beta
        if (!0 %in% I.co) {
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
          beta.it <- beta0
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
                    boot <- try(synth.core(Y[,boot.id], X.boot, D[,boot.id], I=I[,boot.id],
                                       W = W.boot, force = force, r = out$r.cv, CV=0,
                                       tol = tol, AR1 = AR1,
                                       beta0 = beta.it, norm.para = norm.para, boot = 1), silent = TRUE)
                    if ('try-error' %in% class(boot)) {
                        boot0 <- list(att.avg = NA, 
                                      beta = NA,
                                      att = NA)
                        return(boot0)
                    } else {
                        return(boot)
                    }
                
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
                    boot <- try(synth.em(Y = Y[,boot.id], X = X.boot, D = D[,boot.id], I=I[,boot.id],
                                   W = W.boot, force = force, r = out$r.cv,
                                   tol = tol, AR1 = AR1, beta0 = beta.it, norm.para = norm.para, boot = 1), silent = TRUE)
                    if ('try-error' %in% class(boot)) {
                        boot0 <- list(att.avg = NA, 
                                      beta = NA,
                                      att = NA)
                        return(boot0)
                    } else {
                        return(boot)
                    }
                
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
                boot <- try(synth.mc(Y[,boot.id], X.boot, D[,boot.id], I=I[,boot.id],
                               W = W.boot, force = force, 
                               lambda = out$lambda.cv, hasF = out$validF, 
                               CV = 0, tol = tol, AR1 = AR1, beta0 = beta.it, norm.para = norm.para), silent = TRUE)
                if ('try-error' %in% class(boot)) {
                    boot0 <- list(att.avg = NA, 
                                  beta = NA,
                                  att = NA)
                    return(boot0)
                } else {
                    return(boot)
                }
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
                boot <- try(synth.core(Y.boot, X.boot, D.boot, I=I.boot, 
                                   W = W.boot, force = force, r = out$r.cv,
                                   CV = 0, tol = tol, AR1 = AR1,
                                   beta0 = beta.it, norm.para = norm.para, boot = 1), silent = TRUE)

                if ('try-error' %in% class(boot)) {
                    boot0 <- list(eff = NA,
                                  att.avg = NA, 
                                  beta = NA,
                                  att = NA)
                    return(boot0)
                } else {
                    b.out <- list(eff = boot$eff + out$eff,
                                  att = boot$att + out$att,
                                  att.avg = boot$att.avg + out$att.avg)
                    if (p>0) {
                        b.out <- c(b.out, list(beta = boot$beta))
                    }
                    return(b.out)
                }
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
                boot <- try(synth.em(Y.boot, X, D, I=I, W=W, force=force, r=out$r.cv,
                               tol=tol, AR1 = AR1, beta0 = beta.it, norm.para = norm.para, boot = 1), silent = TRUE)



                if ('try-error' %in% class(boot)) {
                    boot0 <- list(eff = NA,
                                  att.avg = NA, 
                                  beta = NA,
                                  att = NA)
                    return(boot0)
                } else {
                    b.out <- list(eff = boot$eff + out$eff,
                              att = boot$att + out$att,
                              att.avg = boot$att.avg + out$att.avg)
                    if (p>0) {
                        b.out <- c(b.out, list(beta = boot$beta))
                    }
                    return(b.out)
                }
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
    conf.lvl.lb <- (1 - conf.lvl)/2
    conf.lvl.ub <- conf.lvl.lb + conf.lvl

    CI.att <- t(apply(att.boot, 1, function(vec) 
        quantile(vec,c(conf.lvl.lb, conf.lvl.ub), na.rm=TRUE)))
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
    CI.avg <- quantile(att.avg.boot, c(conf.lvl.lb, conf.lvl.ub), na.rm=TRUE)
    se.avg <- sd(att.avg.boot, na.rm=TRUE)
    pvalue.avg <- get.pvalue(att.avg.boot)
    est.avg <- t(as.matrix(c(att.avg, se.avg, CI.avg, pvalue.avg)))
    colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")
    rownames(est.avg) <- ""
    
    ## individual effects
    if (inference == "parametric") {
        CI.ind <- apply(eff.boot,c(1,2),function(vec)
            quantile(vec,c(conf.lvl.lb, conf.lvl.ub), na.rm=TRUE)) ## 2*T*Ntr
        est.ind <- array(NA,dim=c(TT, 5, Ntr)) ## eff, se, CI.lower, CI.upper
        est.ind[,1,] <- eff
        est.ind[,2,] <- apply(eff.boot,c(1,2),sd)
        est.ind[,3,] <- CI.ind[1,,]
        est.ind[,4,] <- CI.ind[2,,]
        est.ind[,5,] <- apply(eff.boot,c(1,2),get.pvalue)

        dimnames(est.ind)[[1]] <- rownames(est.att)
        dimnames(est.ind)[[2]] <- c("EFF", "S.E.", "CI.lower", "CI.upper", "p.value")
    }

    colboot <- sapply(1:nboots, function(i){paste("boot",i,sep="")})

    
    ## regression coefficents
    if (p>0) {
        CI.beta<-t(apply(beta.boot, 1, function(vec)
            quantile(vec,c(conf.lvl.lb, conf.lvl.ub), na.rm=TRUE)))
        se.beta<-apply(beta.boot, 1, function(vec)sd(vec,na.rm=TRUE))
        pvalue.beta <- apply(beta.boot, 1, get.pvalue)
        ## beta[na.pos] <- NA
        est.beta<-cbind(beta, se.beta, CI.beta, pvalue.beta)
        colnames(est.beta)<-c("beta", "S.E.", "CI.lower", "CI.upper", "p.value")
        colnames(beta.boot) <- colboot
    }

    rownames(att.boot) <- rownames(est.att)
    colnames(att.boot) <- colboot

    dimnames(eff.boot)[[1]] <- rownames(est.att)
    dimnames(eff.boot)[[3]] <- colboot


  
    ##storage
    result<-list(inference = inference,
                 est.att = est.att,
                 est.avg = est.avg,
                 att.boot = att.boot,
                 eff.boot = eff.boot
                 )
    if (p>0) {
        result <- c(result,list(beta.boot = beta.boot))
    }
    
    if (inference == "parametric") {
        result<-c(result,list(est.ind = est.ind))
    }
    if (p>0) {
        result<-c(result,list(est.beta = est.beta))
    } 

    return(c(out,result))

    
} ## end of synth.boot()


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

###################################
## plm for initial values
###################################
initialFit <- function(data,
                       force, 
                       oci) {

    p <- dim(data)[2] - 3

    data2 <- as.data.frame(data)
    colnames.data2 <- c("y","id","time")

    data2[,2] <- as.factor(data2[,2])
    data2[,3] <- as.factor(data2[,3])

    x.name <- NULL
    if (p > 0) {
        for (i in 1:p) {
            x.name <- c(x.name, paste("x", i, sep = ""))
        }
        colnames.data2 <- c(colnames.data2, x.name)
    }

    colnames(data2) <- colnames.data2

    if (p > 0) {
        if (force == 1) {
            Fit.formula <- paste0("y ~", paste0(x.name, collapse = "+"), "+ id")
        } 
        else if (force == 2) {
            Fit.formula <- paste0("y ~", paste0(x.name, collapse = "+"), "+ time")
        }
        else if (force == 3) {
            Fit.formula <- paste0("y ~", paste0(x.name, collapse = "+"), "+ id + time")
        }
        else {
            Fit.formula <- paste0("y ~", paste0(x.name, collapse = "+"))
        }
    } else {
        if (force == 1) {
            Fit.formula <- paste0("y ~ id")
        } 
        else if (force == 2) {
            Fit.formula <- paste0("y ~ time")
        }
        else if (force == 3) {
            Fit.formula <- paste0("y ~ id + time")
        }
        else {
            Fit.formula <- paste0("y ~ 1")
        }
    }
    N <- length(unique(data[,2]))
    T <- length(unique(data[,3]))
    lfit <- lm(as.formula(Fit.formula), data = data2[oci,])
    Y0 <- matrix(predict(lfit, data2), T, N)
    if (p > 0) {
        beta0 <- as.matrix(lfit$coefficients[2:(p+1)])
        if (sum(is.na(beta0)) > 0) {
            beta0[which(is.na(beta0))] <- 0
        }
    } else {
        beta0 <- matrix(0, 1, 1)
    }
    result <- list(Y0 = Y0, beta0 = beta0)
    return(result)
}
