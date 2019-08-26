## an auxiliary function for estimaing cumulative treatment effect based on gsynth 

## 1. cumulative effect for a specified event window
## 2. averaging treatment effect in a sub-group


cumuEff <- function(x,                ## a gsynth object 
	                cumu = TRUE,      ## whether to calculate cumulative effect
	                id = NULL,        ## units to be averaged on 
	                period = NULL) {  ## event window 

	eff <- x$eff 
	D <- x$D.tr
	I <- x$I.tr 
	inference <- x$inference ## bootstrap type 
	tr.all <- colnames(eff) ## all treated units 

	TT <- dim(eff)[1]
	N <- dim(eff)[2]

	if (sum(is.na(c(D))) > 0) {
		D[which(is.na(D))] <- 0
	}
	D.old <- D
	D <- apply(D, 2, function(vec){cumsum(vec)})
	minT0 <- min(apply(D, 2, function(vec){sum(vec == 0)}))
	D <- D.old

	period.raw <- c(0, TT - minT0) 
	if (!is.null(period)) {
		if (period[2] > period.raw[2]) {
		    stop(paste("Ending period should not be greater than ", period.raw[2], sep = ""))
	    }
	} else {
		period <- period.raw
	}

	id.pos <- c()
	if (!is.null(id)) {
		## if (is.null(x$est.ind)) {
		## 	stop("No estimation results at unit level.")
		## } else {
			for (i in 1:length(id)) {
				if (!id[i] %in% tr.all) {
					stop(paste(id[i], " is not in the treated units.\n", sep = ""))
				} else {
					id.pos <- c(id.pos, which(tr.all == id[i]))
				}
			}
		## }
	}

	catt.boot <- NULL

	## obtain effect 
    if (is.null(id)) {
    	catt <- getEffect(D, I, eff, cumu, period)
    	if (is.null(x$est.avg)) {
    		cat("No uncertainty estimates.")
    	} else {
    		if (sum(c(D.boot[,,1])) == 0) {
    			cat("Cannot get uncertainty estimates.")
    		} else {
    			nboots <- length(x$att.avg.boot)
	    		D.boot <- x$Dtr.boot
	    		I.boot <- x$Itr.boot
	    		eff.boot <- x$eff.boot

	    		catt.boot <- matrix(NA, period[2] - period[1] + 1, nboots)

	    		if (class(D.boot) == "array") {
	    			for (i in 1:nboots) {
		    			catt.boot[, i] <- getEffect(D.boot[,,i], I.boot[,,i], 
		    				                        eff.boot[,,i], cumu, period)
		    		}
	    		} else {
	    			for (i in 1:nboots) {
		    			catt.boot[, i] <- getEffect(D.boot[[i]], I.boot[[i]], 
		    				                        eff.boot[[i]], cumu, period)
		    		}
	    		}

    		}    		
    		
    	}
    } else {
    	catt <- getEffect(as.matrix(D[, id.pos]), as.matrix(I[, id.pos]), 
    		              as.matrix(eff[, id.pos]), cumu, period)

    	if (is.null(x$est.avg)) {
    		cat("No uncertainty estimates.")
    	} else {
    		if (is.null(x$est.ind)) {
				stop("No uncertainty estimates at unit level.")
			} else {
				nboots <- length(x$att.avg.boot)
	    		D.boot <- x$Dtr.boot
	    		I.boot <- x$Itr.boot
	    		eff.boot <- x$eff.boot

	    		catt.boot <- matrix(NA, period[2] - period[1] + 1, nboots)

	    		for (i in 1:nboots) {
	    			subD <- D.boot[,,i]
	    			subI <- I.boot[,,i]
	    			subeff <- eff.boot[,,i]
	    			
	    			catt.boot[, i] <- getEffect(as.matrix(subD[, id.pos]), as.matrix(subI[, id.pos]), 
	    				                        as.matrix(subeff[, id.pos]), cumu, period)
	    		}

			}
    		
    	}
    }

    ## non-parametric test
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

    if (!is.null(catt.boot)) {
    	CI.att <- t(apply(catt.boot, 1, function(vec) 
	        quantile(vec,c(0.025, 0.975), na.rm=TRUE)))
	    se.att <- apply(catt.boot, 1, function(vec) sd(vec, na.rm=TRUE))
	    pvalue.att <- apply(catt.boot, 1, get.pvalue)

	    est.catt <- cbind(catt, se.att, CI.att, pvalue.att)
	    colnames(est.catt) <- c("CATT", "S.E.", "CI.lower", "CI.upper", "p.value")
	    rownames(est.catt) <- period[1]:period[2]
    }

    out <- list(catt = catt)
    if (!is.null(catt.boot)) {
    	## out <- c(out, list(est.catt = est.catt, catt.boot = catt.boot))
    	out <- c(out, list(est.catt = est.catt))
    }

    return(out)

}


getEffect <- function(D,     
	                  I, 
	                  eff, 
	                  cumu, 
	                  period) {

	aeff <- rep(NA, period[2] - period[1] + 1)

	if (sum(is.na(D)) == length(c(D))) {
		return(aeff)
	}
	
	## D.old <- D
	if (sum(is.na(c(D))) > 0) {
		D[which(is.na(D))] <- 0
	}
	D <- apply(D, 2, function(vec){cumsum(vec)})

	TT <- dim(D)[1]
	N <- dim(D)[2]
	
	for (i in 1:N) {
		subd <- D[, i]
		t0 <- sum(subd == 0) 
		D[, i] <- 1:TT - t0
	}

	if (sum(c(I) == 0) > 0) {
		D[which(I == 0)] <- NA
	}

	vd <- c(D)
	veff <- c(eff)
	if (sum(is.na(vd)) > 0) {
		vd.rm <- which(is.na(vd))
	    vd <- vd[-vd.rm]
	    veff <- veff[-vd.rm]
	}

	uniT <- unique(vd)

	ts <- period[1]
	te <- min(period[2], max(vd))
	effT <- ts:te

	

	if (cumu == TRUE) {
		if (sum(!(effT %in% uniT)) == 0) {
			pos <- c()
			for (i in 1:length(effT)) {
				pos <- c(pos, which(vd == effT[i]))
				aeff[i] <- mean(veff[pos]) * i
			}
		}
	} else {
		## only average treatment effect
		ave <- as.numeric(tapply(veff, vd, mean))
		effT2 <- period[1]:period[2]
		for (i in 1:length(effT2)) {
			if (effT2[i] %in% uniT ) {
				aeff[i] <- ave[which(uniT == effT2[i])]
			}
		}
	}

	return(aeff)

}














