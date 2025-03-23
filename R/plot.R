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
plot.gsynth <- function(
    x,
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
    theme.bw = TRUE,
    shade.post = FALSE,
    ...){

  if (type %in% c("raw","missing")){
    # Implementing missing data visualization
    # Extract data from the gsynth object
    data <- x$data
    index <- x$index
    formula <- out$call$formula
    if (type=="missing"){

      # Set default main title if not provided
      if(is.null(main)) {
        main <- "Treatment Status and Missing Data"
      }

      # Call panelview with pre.post=TRUE to show treatment status
      panelView::panelview(
                    data = data,
                    formula,
                    index = index,
                    pre.post = TRUE,
                    main = main,
                    xlab = xlab,
                    ylab = ylab,
                    xlim = xlim,
                    ylim = ylim,
                    axis.adjust = axis.adjust)
    } else {

      # Set default main title if not provided
      if(is.null(main)) {
        main <- "Raw Data"
      }

      # Call panelview with type="outcome" to show outcome values
      panelView::panelview(
                    data = data,
                    formula,
                    index = index,
                    type = "outcome",
                    main = main,
                    xlab = xlab,
                    ylab = ylab,
                    xlim = xlim,
                    ylim = ylim,
                    legendOff = legendOff,
                    axis.adjust = axis.adjust)
    }
  } else {
    class(x) <- "fect"
    p <- fect:::plot.fect(x=x,
      type=type,
      xlim=xlim,
      ylim=ylim,
      xlab = xlab,
      ylab = ylab,
      legendOff=legendOff,
      raw=raw,
      main=main,
      nfactors = nfactors,
      id=id,
      axis.adjust = axis.adjust,
      theme.bw=theme.bw,
      shade.post = shade.post,
        ...)
    return(p)
  }
}


