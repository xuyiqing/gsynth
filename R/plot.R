#######################################################
## METHODS
#######################################################

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
                        theme.bw = FALSE,
                        shade.post = NULL,
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

    scaleFUN <- function(x) sprintf("%.f", x) ## integer value at x axis

    if (class(x)!="gsynth") {
        stop("Not a \"gsynth\" object.")
    }
    if (!type %in% c("gap","counterfactual","ct","factors","missing","loadings","raw")) {
        stop("\"type\" option misspecified.")        
    }
    if (type == "ct") {
      type <- "counterfactual"
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

    ## color of axes
    if (theme.bw == TRUE) {
      line.color <- "#AAAAAA70"
    } else {
      line.color <- "white"
    }

    ## shade in the post-treatment period
    if (is.null(shade.post) == TRUE) {
      if (type %in% c("raw","counterfactual")) {
        shade.post <- TRUE
      }
      if (type %in% c("gap","factors")) {
        shade.post <- FALSE
      }    
    } else {
      if (!class(shade.post) %in% c("logical","numeric")) {
        stop("Wrong type for option \"shade.post\"")
      }
    }
    
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
        p <- ggplot(data) 
        if (theme.bw == TRUE) {
          p <- p + theme_bw()
        }
        ## labels and legend
        p <- p + xlab(xlab) +  ylab(ylab) +
            theme(legend.position = legend.pos,
                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.h),
                  plot.title = element_text(size=20,
                                            hjust = 0.5,
                                            face="bold",
                                            margin = margin(10, 0, 10, 0)))        
        
        if (x$sameT0==TRUE) {
            p <- p + geom_vline(xintercept=time.bf,colour=line.color,size = 2) 
            if (shade.post == TRUE) {
              p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
            }  
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
        if (theme.bw == FALSE) {
          set.colors = c("#FC8D6280","red","#99999950")
        } else {
          set.colors = c("#4671D565","#06266F","#5E5E5E50")
        }
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
        } else {
            p <- p + scale_x_continuous(labels=scaleFUN)
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
            stop(paste(id,"not in the treatment group"))
        } else { ## no error

            ## axes labels
            if (is.null(xlab) == TRUE) {
                if (x$sameT0 == TRUE) {
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
            if (length(id) == 1 && !is.null(x$est.ind)) { ## id specified
                maintext <- paste(x$index[1],"=",id) 
            }  else {
                maintext <- "Estimated ATT"
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
                if (length(id) == 1 && !is.null(x$est.ind)) { ## id specified
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
             
            ### plotting
            p <- ggplot(data) 
            ## black/white theme
            if (theme.bw == TRUE) {
              p <- p + theme_bw()
            }
            p <- p + geom_vline(xintercept = time.bf, colour=line.color,size = 2) +
                geom_hline(yintercept = 0, colour=line.color,size = 2) +
                xlab(xlab) +  ylab(ylab) +
                theme(legend.position = legend.pos,
                      plot.title = element_text(size=20,
                                                hjust = 0.5,
                                                face="bold",
                                                margin = margin(10, 0, 10, 0)))
            if (shade.post == TRUE) {
              p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
            }  
           
            
            ## point estimates
            p <- p + geom_line(aes(time, ATT), size = 1.2)
             
            ## confidence intervals
            if (is.null(x$est.att)==FALSE || !(is.null(x$est.ind)&length(id) == 1)) {
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
            
            suppressWarnings(print(p))
        }  ## end of "gap" (in case of no id error)
       
        
    } else if (type=="counterfactual") { 

        if (length(id) == 1|length(x$id.tr) == 1|x$sameT0==TRUE) { 
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
                        p <- ggplot(data) 
                        if (theme.bw == TRUE) {
                          p <- p + theme_bw()
                        }
                        p <- p + xlab(xlab) +  ylab(ylab) +
                        geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0))) 
                        if (shade.post == TRUE) {
                          p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
                        }        

                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type)) 
                        ## legend
                        set.limits = c("tr","ct")
                        set.labels = c("Treated", "Estimated Y(0)")
                        set.colors = c("black","steelblue")
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
                        } else {
                            p <- p + scale_x_continuous(labels=scaleFUN)
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
                        p <- ggplot(data) 
                        if (theme.bw == TRUE) {
                          p <- p + theme_bw()
                        }
                        p <- p + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0)))
                        if (shade.post == TRUE) {
                          p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
                        }      

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
                        set.colors = c("black","#4682B480","steelblue")
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
                                scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                        } else {
                            p <- p + scale_x_continuous(labels=scaleFUN)
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
                        p <- ggplot(data) 
                        if (theme.bw == TRUE) {
                          p <- p + theme_bw()
                        }
                        p <- p + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0)))
                        if (shade.post == TRUE) {
                          p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
                        }  
                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type,
                                               group = id))

                        ## legend
                        set.limits = c("tr","raw.co","ct")
                        set.labels = c("Treated","Controls","Estimated Y(0)")
                        set.colors = c("black","#4682B420","steelblue")
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
                        } else {
                            p <- p + scale_x_continuous(labels=scaleFUN)
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
                        p <- ggplot(data) 
                        if (theme.bw == TRUE) {
                          p <- p + theme_bw()
                        }
                        p <- p + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0)))
                        if (shade.post == TRUE) {
                          p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
                        }      
                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type))

                        ## legend
                        set.limits = c("tr","co")
                        set.labels = c("Treated Average",
                                       "Estimated Y(0) Average")
                        set.colors = c("black","steelblue")
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
                        } else {
                            p <- p + scale_x_continuous(labels=scaleFUN)
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
                        p <- ggplot(data) 
                        if (theme.bw == TRUE) {
                          p <- p + theme_bw()
                        }
                        p <- p + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0)))
                        if (shade.post == TRUE) {
                          p <- p + annotate("rect", xmin= time.bf, xmax= Inf, ymin=-Inf, ymax=Inf, alpha = .3) 
                        }      
                        ## main
                        p <- p + geom_line(aes(time, outcome,
                                               colour = type,
                                               size = type,
                                               linetype = type))
                        ## band
                        p <- p + geom_ribbon(data = data.band,
                                             aes(ymin = tr5, ymax = tr95, x=time),
                                             alpha = 0.15, fill = "black") +
                            geom_ribbon(data = data.band,
                                        aes(ymin = co5, ymax = co95, x=time),
                                        alpha = 0.15, fill = "steelblue")

                        set.limits = c("tr","co","tr.band","co.band")
                        set.labels = c("Treated Average",
                                       "Estimated Y(0) Average",
                                       "Treated 5-95% Quantiles",
                                       "Controls 5-95% Quantiles")
                        set.colors = c("black","steelblue","#77777750","#4682B480")
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
                        } else {
                            p <- p + scale_x_continuous(labels=scaleFUN)
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
                        p <- ggplot(data) 
                        if (theme.bw == TRUE) {
                          p <- p + theme_bw()
                        }
                        p <- p + xlab(xlab) +  ylab(ylab) +
                            geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                            theme(legend.position = legend.pos,
                                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                                  plot.title = element_text(size=20,
                                                            hjust = 0.5,
                                                            face="bold",
                                                            margin = margin(10, 0, 10, 0))) 
                        if (shade.post == TRUE) {
                          p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
                        }      
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
                        set.colors = c("black","steelblue","#77777750","#4682B420")
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
                        } else {
                            p <- p + scale_x_continuous(labels=scaleFUN)
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
                p <- ggplot(data) 
                if (theme.bw == TRUE) {
                  p <- p + theme_bw()
                }
                p <- p  + xlab(xlab) +  ylab(ylab) +
                    geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                    theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0)))
                if (shade.post == TRUE) {
                  p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
                }      
                ## main
                p <- p + geom_line(aes(time, outcome,
                                       colour = type,
                                       size = type,
                                       linetype = type))

                ## legend
                set.limits = c("tr","co")
                set.labels = c("Treated Average",
                               "Estimated Y(0) Average")
                set.colors = c("black","steelblue")
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
                } else {
                    p <- p + scale_x_continuous(labels=scaleFUN)
                }
                    
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
                p <- ggplot(data) 
                if (theme.bw == TRUE) {
                  p <- p + theme_bw()
                }
                p <- p  + xlab(xlab) +  ylab(ylab) +
                    geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                    theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0)))
                if (shade.post == TRUE) {
                  p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
                }      
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
                set.colors = c("black","steelblue","#77777750")
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

                if (!is.numeric(time.label)) {
                    p <- p + 
                        scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                } else {
                    p <- p + scale_x_continuous(labels=scaleFUN)
                } 
                    
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
                p <- ggplot(data) 
                if (theme.bw == TRUE) {
                  p <- p + theme_bw()
                }
                p <- p + xlab(xlab) +  ylab(ylab) +
                    geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                    theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0))) 
                if (shade.post == TRUE) {
                  p <- p + annotate("rect", xmin= time.bf, xmax= Inf,ymin=-Inf, ymax=Inf, alpha = .3) 
                }      
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
                set.colors = c("black","steelblue","#77777750")
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

                if (!is.numeric(time.label)) {
                    p <- p + 
                        scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                } else {
                    p <- p + scale_x_continuous(labels=scaleFUN)
                }
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
            ## prepare data
            L.co<-x$lambda.co
            norm<-sqrt(diag(t(L.co)%*%L.co)/(x$N-x$Ntr))
            data <- cbind.data.frame("time" = rep(time[show],r),
                                     "factor" = c(F.hat[show,])*rep(norm,each=nT),
                                     "group" = as.factor(c(rep(1:r,each=nT))))
            ## theme
            p <- ggplot(data) 
            if (theme.bw == TRUE) {
              p <- p + theme_bw()
            }
            p <- p + xlab(xlab) +  ylab(ylab) + ggtitle(main) +
                geom_hline(yintercept=0,colour=line.color,size = 2) +
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


            brew.colors <- c("black","steelblue","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9")
            set.colors = brew.colors[1:r]
            p <- p + scale_colour_manual(values =set.colors) 

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
