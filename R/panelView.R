#################################################################
## preview of data treatment status, missing values and raw data
#################################################################

panelView <- function(formula,
                      data, # a data frame (long-form)
                      index, # c(unit, time) indicators
                      na.rm = FALSE, # remove missing values
                      type = "missing",
                      by.group = FALSE, # (color pre-treatment treated differently)
                      xlim = NULL, 
                      ylim = NULL,
                      xlab = NULL, 
                      ylab = NULL,
                      legendOff = FALSE,
                      legend.labs = NULL,
                      main = NULL,
                      id = NULL,
                      show.id = NULL,
                      axis.adjust = FALSE,
                      axis.lab = "both",
                      axis.lab.gap = c(0, 0)
                    ) {  
    ## ------------------------------------- ##
    ##          part 1: parsing data
    ## ------------------------------------- ##
    varnames <- all.vars(formula)
    Y <- varnames[1]
    D <- varnames[2]
    if (length(varnames) > 2) {
        X <- varnames[3:length(varnames)]
    } else {
        X <- NULL
    }

    p <- length(X)

    ##store variable names
    id.old <- id
    data.old <- data
    Yname <- Y
    Dname <- D
    Xname <- X
    
    id <- index[1]
    time <- index[2]

    TT <- length(unique(data[,time]))
    N <- length(unique(data[,id]))
    id.series <- unique(sort(data[,id])) ## unit id
    time.uni <- unique(sort(data[,time])) ## period
    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------## 
    ## library(ggplot2)

    if (!type %in% c("missing", "raw")) {
        stop("\"type\" option misspecified.")
    }

    if (is.data.frame(data) == FALSE) {
        data <- as.data.frame(data)
        warning("Not a data frame.")
    }
    ## index
    if (length(index) != 2 | sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
    }
    
    unique_label <- unique(paste(data[,index[1]],"_",data[,index[2]],sep=""))
    if (length(unique_label)!=dim(data)[1]) {
        stop("Some records may be replicated or wrongly marked in the data set.")
    }

    ## remove missing values
    if (is.logical(na.rm) == FALSE & !na.rm%in%c(0, 1)) {
        stop("na.rm is not a logical flag.")
    } 

    ## axis.lab
    if (!axis.lab %in% c("both", "unit", "time", "off")) {
        stop("\"axis.lab\" option misspecified. Try, for example, axis.lab = c(\"both\", \"unit\", \"time\", \"off\").") 
    }
    ## time labels gap
    if (sum(axis.lab.gap < 0) > 0) {
        stop("\"gap\" should be equal to or greater than 0.\n")
    }

    ## legend labels
    if (is.null(legend.labs)==FALSE) {
        legend.labs <- as.character(legend.labs)
    }

    ## missing plot : y-axis
    if (!is.null(id.old)) {
        if (!is.null(show.id)) {
           warning("Using \"id\" option.\n") 
        }
    } else {
        if (!is.null(show.id)) {
            if (length(show.id) > N ) {
                stop("Length of \"show.id\" should not be larger than total number of units. \n")
            }
            if (!class(show.id) %in% c("numeric", "integer")) {
                stop("\"show.id\" option misspecified. Try, for example, show.id = 1:100. \n")
            }
            if (sum(show.id > N) > 0) {
                stop("Some specified units are not in the data.\n")
            }
            if (length(unique(show.id)) != length(show.id)) {
                stop("Repeated values in \"show.id\" option.")
            }
            id.old <- id.series[show.id] 
        }
    }


    ##-------------------------------##
    ## Parsing raw data
    ##-------------------------------##  

    ## sort data
    data <- data[order(data[,id], data[,time]), ]
    id.all <- time.all <- count <- coordin <- data.x <- x.na <- NULL
    
    ## if (p > 0) { ## covariates may contain missing values
        ## location
    ##     id.all <- as.matrix(data[, id])
    ##     time.all <- as.matrix(data[, time])

    ##     count <- dim(data)[1]
    ##     coordin <- matrix(NA, count, 2)
    ##     for (i in 1:count) {
    ## 	    coordin[i, 2] <- which(id.series == id.all[i])
    ##	    coordin[i, 1] <- which(time.uni == time.all[i]) ## obs.missing: T*N
    ##     }
    ##     data.x <- as.matrix(data[, X])
    ##     x.na <- apply(is.na(data.x), 1, sum)
    ##     x.na <- ifelse(x.na > 0, 1, 0)
    ## }

    if (na.rm == TRUE) {
        data <- data[,c(index, Y, D, X)] ## covariates are excluded
        data <- na.omit(data)
    } 

    ## check missingness
    if (sum(is.na(data[, Yname])) > 0) {
        stop(paste("Missing values in variable \"", Yname,"\".", sep = ""))
    }
    if (sum(is.na(data[, Dname])) > 0) {
        stop(paste("Missing values in variable \"", Dname,"\".", sep = ""))
    }

    if (!(1%in%data[, Dname]&0%in%data[,Dname]&length(unique(data[,Dname]))==2)) {
        stop(paste("Error values in variable \"", Dname,"\".", sep = ""))
    }

    if (sum(is.na(data[, id])) > 0) {
        stop(paste("Missing values in variable \"", id,"\".", sep = ""))
    }
    if (sum(is.na(data[, time])) > 0) {
        stop(paste("Missing values in variable \"", time,"\".", sep = ""))
    } 

    ## check balanced panel and fill unbalanced panel
    if (var(table(data[,id])) + var(table(data[, time])) > 0 | TT == N) {
        
        data[,time] <- as.numeric(as.factor(data[,time]))
        ob <- "time_ob_ls"
        
        while(ob%in%colnames(data)){
            ob <- paste(ob,ob,sep="")
        }

        data[,ob] <- data[,time]
        for (i in 1:N) {
            data[data[,id]==id.series[i],ob] <- data[data[,id]==id.series[i],time]+(i-1)*TT  
        }

        variable <- c(Yname, Dname)
       
        data_I <- matrix(0,N*TT,1)
        data_I[c(data[,ob]),1] <- 1
        data_ub <- as.matrix(data[,variable])
        data <- data_ub_adj(data_I,data_ub)
        colnames(data) <- variable
    }

    ## index matrix that indicates if data is observed 
    I <- matrix(1, TT, N)
    Y.ind <- matrix(data[,Yname], TT, N)
    I[is.nan(Y.ind)] <- 0

    if (0%in%I) {
        data[is.nan(data)] <- 0
    }
    
    ##treatment indicator
    D <- matrix(data[,Dname],TT,N)
    D.old <- D ## store D 

    ## once treated, always treated
    D <- apply(D, 2, function(vec){cumsum(vec)})
    D <- ifelse(D > 0, 1, 0)

    ## DID timing
    tr.old <- D[TT,]==1
    D.tr.old <- as.matrix(D[,which(tr.old==1)])
    T0 <- apply(D.tr.old==0,2,sum) 
    DID <- length(unique(T0))==1

    ## check DID mode
    if (sum(abs(D.old[which(I==1)] - D[which(I==1)])) == 0) {
        by.group <- by.group
    } else { ## FE mode
        if (by.group == FALSE & type == "raw") {
            warning("Treatment has reversals. Show by groups.\n")
        }
        by.group <- TRUE
    }

    ##-------------------------------##
    ## storage
    ##-------------------------------##    
    iname <- unique(sort(data.old[,id]))
    tname <- unique(sort(data.old[,time]))

    ##outcome variable
    Y <- matrix(data[,Yname],TT,N)
    Y[which(I==0)] <- NA

    ## plot a part of outcomes: raw plot
    if (!is.null(id.old)) {
        if (length(unique(id.old)) != length(id.old)) {
            stop("Repeated values in \"id\" option.")
        }

        id.pos <- rep(NA, length(id.old))
        for (i in 1:length(id.old)) {
            if (id.old[i]%in%iname) {
                id.pos[i] <- which(iname == id.old[i])
            } else {
                stop("Some specified units are not in the data.")
            }
        }
        N <- length(id.old)
        Y <- as.matrix(Y[, id.pos])
        I <- as.matrix(I[, id.pos])
        D <- as.matrix(D[, id.pos])
        D.old <- as.matrix(D.old[, id.pos])
        iname <- id.old

    }

    if (by.group == FALSE) {
        tr <- D[TT,]==1     # cross-sectional: treated unit
        pre <- as.matrix(D[,which(tr==1)]==0&I[,which(tr==1)]==1) # a matrix indicating before treatment
        post <- as.matrix(D[,which(tr==1)]==1&I[,which(tr==1)]==1)
        id.tr <- which(tr==1)
        id.co <- which(tr==0)

        D.tr <- as.matrix(D[,which(tr==1)])
        I.tr <- as.matrix(I[,which(tr==1)])
        Y.tr <- as.matrix(Y[,which(tr==1)])
        Y.co <- as.matrix(Y[,which(tr==0)])

        Ntr <- sum(tr)
        Nco <- N - Ntr

        obs.missing <- matrix(3, TT, N) ## control group:1
        pre[which(pre==1)] <- 1 ## pre 2
        post[which(post==1)] <- 2 ## post 3
        obs.missing[,id.tr] <- pre + post
        obs.missing[which(I==0)] <- 4 ## missing

        id.tr <- iname[which(tr==1)] ## re-define id
        id.co <- iname[which(tr==0)]
    } else {
        unit.type <- rep(NA, N) ## 1 for control; 2 for treated; 3 for reversal
        for (i in 1:N) {
            di <- D.old[, i]
            ii <- I[, i]
            if (length(unique(di[which(ii==1)])) == 1) { ## treated or control
                if (0 %in% unique(di[which(ii==1)])) {
                    unit.type[i] <- 1 ## control
                } else {
                    unit.type[i] <- 2 ## treated
                }
            } else {
                unit.type[i] <- 3 ## reversal
            }
        }
        
        obs.missing <- matrix(2, TT, N) ## not under treatment
        obs.missing[which(D.old==1)] <- 1 ## under treatment
        obs.missing[which(I==0)] <- 4 ## missing
    }
        
    ## obs.missing[which(obs.missing==1)] <- "control"
    ## obs.missing[which(obs.missing==2)] <- "pre"
    ## obs.missing[which(obs.missing==3)] <- "post"
    ## obs.missing[which(obs.missing==4)] <- "removed"
    ## obs.missing[which(obs.missing==0)] <- "missing"
    
    ## if (p > 0 & 1 %in% x.na ) {
    ##     for (i in 1:count) {
    ##     	if (x.na[i] == 1) {
    ##    		obs.missing[coordin[i, 1], coordin[i, 2]] <- -1 ## missing -1(covar missing)
    ##    	}
    ##    }
    ## }

    colnames(obs.missing) <- iname
    rownames(obs.missing) <- tname

    time <- tname 
    id <- id.old ## recover parameter from the function

    ## ------------------------------------- ##
    ##          part 2: plot
    ## ------------------------------------- ##
    outcome <- NULL ## global variable
    labels1 <- labels2 <- labels3 <- NULL
    
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

    if (is.null(main)==FALSE) {
        if (is.character(main) == FALSE) {
            stop("\"main\" is not a string.")
        } else {
            main <- main[1]
        }   
    }

     if (axis.adjust == TRUE) {
        angle <- 45
        x.v <- 1
        x.h <- 1
    } else {
        angle <- 0
        x.v <- 0
        if (type == "missing") {
            x.h <- 0.5
        } else {
            x.h <- 0
        }
    }
    
    if (type == "missing") {
        if (is.null(id) == TRUE) {
            ## if (is.null(show.id) == TRUE) {
                id <- colnames(obs.missing)
            ## } else {
            ##     id <- colnames(obs.missing)[show.id]
            ## }
        }
        m.l <- length(id)
        for (i in 1:m.l) {
            if (!id[i]%in%colnames(obs.missing)) {
                stop("Some specified units are not in the data.")
            }
        }
    }
  
    ## type of plots
    if (!is.numeric(time[1])) {
        time <- 1:TT
    }

    ## periods to show
    if (length(xlim) != 0) {
        show <- which(time>=xlim[1] & time<=xlim[2])
    } else {
        show <- 1:length(time)
    }     

    nT <- length(show)
    time.label <- tname[show]
    T.b <- 1:length(show)
 
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
            xlab <- index[2]
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab)==TRUE) {
            ylab <- Yname
        } else if (ylab == "") {
            ylab <- NULL
        }

        if (by.group == FALSE) { ## DID-type plot
            ## time-line
            if (Ntr >= 1) { ## with treatment group
                if (length(id) == 1) {
                    time.bf <- time[T0[which(id == id.tr)]]
                } else {
                    time.bf <- time[unique(T0)]
                }
                pst <- D.tr
                for (i in 1:Ntr) {
                    pst[T0[i], i] <- 1 ## paint the period right before treatment
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
            } else {
                data <- cbind.data.frame("time" = c(rep(time[show], N)),
                                         "outcome" = c(c(Y.tr[show,]),
                                                       c(Y.co[show,])),
                                         "type" = c(rep("tr",(Ntr*nT)),
                                                    rep("co",(Nco*nT))),
                                         "id" = c(rep(1:N,each = nT)))
            }
        
            ## theme
            p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                theme(legend.position = legend.pos,
                      axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.h),
                      plot.title = element_text(size=20,
                                                hjust = 0.5,
                                                face="bold",
                                                margin = margin(10, 0, 10, 0)))

        
        
            if (DID == TRUE && Ntr >= 1) {
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
            set.colors = c("#FC8D6280","red","#99999950")
            set.linetypes = c("solid","solid","solid")
            set.linewidth = c(0.5, 0.5, 0.5)
            if (!is.null(legend.labs)) {
               if (length(legend.labs) != 3) {
                   warning("Wrong number of labels in the legend. Using default.\n")
                   set.labels = c("Treated (Pre)","Treated (Post)","Controls")  
               } else {
                   set.labels <- legend.labs
               }
            } else {
                set.labels = c("Treated (Pre)","Treated (Post)","Controls") 
            }            
            
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
            ## end of raw plot
        
        } else { ## FE-type plot
            
            ## sequence: always control, always treated, reversal
            ## for (i in 1:(TT-1)) {
            ##     for (j in 1:N) {
            ##         if (D.old[i+1, j] == 1) {
            ##             D.plot[i, j] <- 1 ## link the last period
            ##         }
            ##     }
            ## }

            if (is.null(ylim) == TRUE) {
                ylim <- c(min(c(Y[show,]), na.rm = TRUE), max(c(Y[show,]), na.rm = TRUE))
            }

            if (is.null(main) == TRUE) {
                main <- "Raw Data"
            }

            if (!is.null(legend.labs)) {
               if (length(legend.labs) != 2) {
                   warning("Wrong number of labels in the legend. Using default.\n")
                   set.labels = c("Control", "Treatment") 
               } else {
                   set.labels <- legend.labs
               }
            } else {
                set.labels = c("Control", "Treatment") 
            } 

            

            if (1 %in% unit.type) {
                co.pos <- which(unit.type == 1)
                Nco <- length(co.pos)
                data1 <- cbind.data.frame("time" = c(rep(time[show], Nco)),
                                          "outcome" = c(Y[show, co.pos]),
                                          "type" = c(rep("co", (Nco*nT))),
                                          "id" = c(rep(1:Nco, each = nT)))
                limits1 <- c("co", "tr")
                #limits1 <- "co"
                colors1 <- c("#99999950", "#FC8D6280")
                #colors1 <- "#99999950"
                ## main1 <- main
                main1 <- "Always Under Control"
            }

            if (2 %in% unit.type) {
                tr.pos <- which(unit.type == 2)
                Ntr <- length(tr.pos)
                data2 <- cbind.data.frame("time" = c(rep(time[show], Ntr)),
                                          "outcome" = c(Y[show, tr.pos]),
                                          "type" = c(rep("tr",(Ntr*nT))),
                                          "id" = c(rep(1:Ntr,each = nT)))
                limits2 <- c("co", "tr")
                #limits2 <- "tr"
                colors2 <- c("#99999950", "#FC8D6280") 
                #colors2 <- "#FC8D6280" 
                ## main2 <- ifelse(1%in%unit.type, "", main)
                main2 <- "Always Under Treatment"
            }

            if (3 %in% unit.type) {
                rv.pos <- which(unit.type == 3)
                Nrv <- length(rv.pos)

                D.plot <- D.old
                #for (i in 1:(TT-1)) {
                #    for (j in 1:N) {
                #        if (D.old[i+1, j] == 1) {
                #            D.plot[i, j] <- 1 ## link the last period
                #        }
                #    }
                #}

                D.plot[which(D.plot == 0)] <- NA
                D.plot[which(I == 0)] <- NA

                D.rv <- as.matrix(D.plot[, rv.pos])
                Y.rv <- as.matrix(Y[, rv.pos])

                Y.trt <- Y.rv * D.rv
                Y.trt.show <- as.matrix(Y.trt[show,])
                time.trt.show <- time[show]
                ut.time <- ut.id <- NULL
                for (i in 1:Nrv) {
                    if (sum(is.na(Y.trt.show[,i])) != nT) {
                        ut.id <- c(ut.id, rep(i, nT - sum(is.na(Y.trt.show[,i]))))
                        ut.time <- c(ut.time, time.trt.show[which(!is.na(Y.trt.show[,i]))])
                    }
                }


                data3 <- cbind.data.frame("time" = c(rep(time[show], Nrv), ut.time),
                                          "outcome" = c(c(Y[show, rv.pos]),
                                                      c(Y.trt.show[which(!is.na(Y.trt.show))])),
                                          "type" = c(rep("nut",(Nrv*nT)),
                                                   rep("ut",length(ut.id))),
                                          "id" = c(rep(1:Nrv,each = nT), ut.id))
                limits3 <- c("nut", "ut")
                colors3 <- c("#99999950", "#FC8D6280")
                ##if (1%in%unit.type||2%in%unit.type) {
                ##    main3 <- ""   
                ##} else {
                ##    main3 <- main 
                ##}
                main3 <- "Treatment Status Changed"
            }
            
            subplot <- function (data, limits, labels, colors, main) {
                ## theme
                p <- ggplot(data) + xlab(xlab) +  ylab(ylab) +
                    theme(axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.h),
                          plot.title = element_text(size=10,
                                                    hjust = 0.5,
                                                    margin = margin(10, 0, 10, 0)))

                ## main
                p <- p + geom_line(aes(time, outcome,
                                       colour = type,
                                       size = type,
                                       linetype = type,
                                       group = id))

                ## legend
                set.limits = limits
                set.colors = colors
                set.linetypes = rep("solid", length(limits))
                set.linewidth = rep(0.5, length(limits))

                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                    scale_linetype_manual(limits = set.limits,
                                          labels = set.labels,
                                          values = set.linetypes) +
                    scale_size_manual(limits = set.limits,
                                      labels = set.labels,
                                      values = set.linewidth) +
                    guides(linetype = guide_legend(title=NULL, nrow=1),
                           colour = guide_legend(title=NULL, nrow=1),
                           size = guide_legend(title=NULL, nrow=1)) 
        
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
                return(p)
            }

            if (length(unique(unit.type))==1) {
                if (1%in%unit.type) {
                    p1 <- subplot(data1, limits1, labels1, colors1, main1)
                    if (legend.pos != "none") {
                        suppressWarnings(g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs)
                        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                        suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                         legend, nrow = 2, heights = c (1, 1/5)),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))
                    } else {
                        suppressWarnings(grid.arrange(p1 + theme(legend.position="none"),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))
                    }   
                }
                else if (2%in%unit.type) {
                    p2 <- subplot(data2, limits2, labels2, colors2, main2)
                    if (legend.pos != "none") {
                        suppressWarnings(g <- ggplotGrob(p2 + theme(legend.position="bottom"))$grobs)
                        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                        suppressWarnings(grid.arrange(arrangeGrob(p2 + theme(legend.position="none"),
                                         legend, nrow = 2, heights = c (1, 1/5)),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2)))) 
                    } else {
                        suppressWarnings(grid.arrange(p2 + theme(legend.position="none"),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))
                    }   
                }
                else if (3%in%unit.type) {
                    p3 <- subplot(data3, limits3, labels3, colors3, main3)
                    if (legend.pos != "none") {
                        suppressWarnings(g <- ggplotGrob(p3 + theme(legend.position="bottom"))$grobs)
                        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                        suppressWarnings(grid.arrange(arrangeGrob(p3 + theme(legend.position="none"),
                                         legend, nrow = 2, heights = c (1, 1/5)),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2)))) 
                    } else {
                        suppressWarnings(grid.arrange(p3 + theme(legend.position="none"),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))
                    }  
                }
            }
            else if (length(unique(unit.type))==2) {
                if (!1%in%unit.type) {
                    p2 <- subplot(data2, limits2, labels2, colors2, main2)
                    p3 <- subplot(data3, limits3, labels3, colors3, main3)
                    if (legend.pos != "none") {
                        suppressWarnings(g <- ggplotGrob(p2 + theme(legend.position="bottom"))$grobs)
                        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                        suppressWarnings(grid.arrange(arrangeGrob(p2 + theme(legend.position="none"), p3 + theme(legend.position="none"),
                                         legend, nrow = 3, heights = c (1, 1, 1/5)),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))  
                    } else {
                        suppressWarnings(grid.arrange(p2 + theme(legend.position="none"),
                                         p3 + theme(legend.position="none"),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))
                    }  
                }
                else if (!2%in%unit.type) {
                    p1 <- subplot(data1, limits1, labels1, colors1, main1)
                    p3 <- subplot(data3, limits3, labels3, colors3, main3)
                    if (legend.pos != "none") {
                        suppressWarnings(g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs)
                        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                        suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p3 + theme(legend.position="none"),
                                         legend, nrow = 3, heights = c (1, 1, 1/5)),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))  
                    } else {
                        suppressWarnings(grid.arrange(p1 + theme(legend.position="none"),
                                         p3 + theme(legend.position="none"),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))
                    }  
                }
                else if (!3%in%unit.type) {
                    p1 <- subplot(data1, limits1, labels1, colors1, main1)
                    p2 <- subplot(data2, limits2, labels2, colors2, main2)
                    if (legend.pos != "none") {
                        suppressWarnings(g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs)
                        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                        suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"),
                                         legend, nrow = 3, heights = c (1, 1, 1/5)),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))  
                    } else {
                        suppressWarnings(grid.arrange(p1 + theme(legend.position="none"),
                                         p2 + theme(legend.position="none"),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))
                    }   
                }
            }
            else {
                p1 <- subplot(data1, limits1, labels1, colors1, main1)
                p2 <- subplot(data2, limits2, labels2, colors2, main2)
                p3 <- subplot(data3, limits3, labels3, colors3, main3)
                if (legend.pos != "none") {
                    suppressWarnings(g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs)
                    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
                    suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"),
                                     p3 + theme(legend.position="none"), legend, nrow = 4, heights = c (1, 1, 1, 1/5)),
                                     top = textGrob(main, gp = gpar(fontsize=20,font=2))))
                } else {
                    suppressWarnings(grid.arrange(p1 + theme(legend.position="none"),
                                         p2 + theme(legend.position="none"),
                                         p3 + theme(legend.position="none"),
                                         top = textGrob(main, gp = gpar(fontsize=20,font=2))))
                }
            }
            ## end of raw plot
        }
        
    } else if (type=="missing") {
        
        if (is.null(xlab)==TRUE) {
            xlab <- index[2]
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab)==TRUE) {
            ylab <- index[1]
        } else if (ylab == "") {
            ylab <- NULL
        }
        if (is.null(main)==TRUE) {
            main <- "Treatment Status"
        } else if (main == "") {
            main <- NULL
        }

        m <- obs.missing
        
        if (!is.null(id)) {
            m <- as.matrix(m[show,which(colnames(m)%in%id)])
        } else {
            ## if (!is.null(show.id)) {
            ##     m <- as.matrix(m[show, c(show.id)])
            ## } else {
                m <- as.matrix(m[show,])
            ## }
        }

        all <- unique(c(m))
        col <- breaks <- label <- NULL
        if (3 %in% all) { ## DID case
            if (1%in%all) {
                col <- c(col,"#4671D5")
                breaks <- c(breaks,1)
                label <- c(label,"Treated (Pre)")
            }
            if (2%in%all) {
                col <- c(col,"#06266F")
                breaks <- c(breaks,2)
                label <- c(label,"Treated (Post)")
            }
            if (3%in%all) {
                col <- c(col,"#B0C4DE")
                breaks <- c(breaks,3)
                label <- c(label,"Controls")
            }
        } else {
            if (1%in%all) {
                col <- c(col,"#06266F")
                breaks <- c(breaks,1)
                label <- c(label,"Under Treatment")
            }
            if (2%in%all) {
                col <- c(col,"#B0C4DE")
                breaks <- c(breaks,2)
                label <- c(label,"Under Control")
            }
        }
        if (4%in%all) {
            col <- c(col,"#FFFFFF")
            breaks <- c(breaks,4)
            label <- c(label,"Missing")
        }
        if (!is.null(legend.labs)) {
            if (length(legend.labs) != length(all)) {
                warning("Wrong number of labels in the legend. Using default.\n")
            } else {
                label <- legend.labs
            }
        } 
        
        N <- dim(m)[2]
        units <- rep(rev(1:N), each = TT)
        period <- rep(1:TT, N)
        res <- c(m)
        data <- cbind.data.frame(units=units, period=period, res=res)
        data[,"res"] <- as.factor(data[,"res"])

        ## check if N >= 200
        if (dim(m)[2] >= 200) {
            if (axis.lab == "both") {
                axis.lab <- "time"
            }
            else if (axis.lab == "unit") {
                axis.lab <- "off"
            }
        }

        ## labels
        N.b <- 1:N
        if (axis.lab == "both") {
            if (length(axis.lab.gap)==2) {
                x.gap <- axis.lab.gap[1]
                y.gap <- axis.lab.gap[2] 
            } else {
                x.gap <- y.gap <- axis.lab.gap[1]
            }
        } else {
            x.gap <- y.gap <- axis.lab.gap[1]
        }

        if (x.gap != 0) {
            T.b <- seq(from = 1, to = length(show), by = (x.gap + 1))
        }
        if (y.gap != 0) {
            N.b <- seq(from = N, to = 1, by = -(y.gap + 1))
        }
        id <- rev(id)
        
        p <- ggplot(data, aes(x = period, y = units,
                              fill = res), position = "identity") 
        p <- p + geom_tile(colour="gray90", size=0.1, stat="identity") 
  
        p <- p +
            labs(x = xlab, y = ylab,  
                title=main) +
            theme_bw() + 
            scale_fill_manual(NA, breaks = breaks, values = col, labels=label)

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
                                        margin = margin(10, 0, 10, 0)))
        if (axis.lab == "both") {
            p <- p + scale_x_continuous(expand = c(0, 0), breaks = T.b, labels = time.label[T.b]) +
            scale_y_continuous(expand = c(0, 0), breaks = N.b, labels = id[N.b])
        }
        else if (axis.lab == "unit") {
            p <- p + scale_x_continuous(expand = c(0, 0), breaks = T.b, labels = NULL) +
            scale_y_continuous(expand = c(0, 0), breaks = N.b, labels = id[N.b])            
        }
        else if (axis.lab == "time") {
            p <- p + scale_x_continuous(expand = c(0, 0), breaks = T.b, labels = time.label[T.b]) +
            scale_y_continuous(expand = c(0, 0), breaks = N.b, labels = NULL)
        }
        else if (axis.lab == "off") {
            p <- p + scale_x_continuous(expand = c(0, 0), breaks = 1:length(show), labels = NULL) +
            scale_y_continuous(expand = c(0, 0), breaks = 1:N, labels = NULL)
        }
        
        if(length(all)>=4) {
            p <- p + guides(fill=guide_legend(nrow=2,byrow=TRUE))
        }
        suppressWarnings(print(p))
        ## end of missing plot
    }    
    ## if (type == "missing") {
    ##     out <- list(plot = p, obs.missing = obs.missing)
    ## } else {
    ##     out <- list(plot = p)
    ## }
    ## return(out)   
}