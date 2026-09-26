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
    legacy.style = FALSE,
    highlight = NULL,
    highlight.fill = FALSE,
    ...){

  if (type %in% c("raw","missing")){
    # Raw data / missing-data plots are drawn by panelView from the data
    # and the variable names stored in the fit (not from x$call, which may
    # hold Y =/D = strings, a formula in a variable, or as.formula(...)).
    # [[ ]] matches names exactly (the fit also has "data.long", and a
    # second "Y"/"D" element holding the T x N matrices).
    if (is.null(x[["data"]]) || !is.character(x[["Y"]]) ||
        !is.character(x[["D"]]) || is.null(x[["index"]])) {
      stop("This fit does not store the data needed for type = \"raw\" or \"missing\".",
           call. = FALSE)
    }
    pv.args <- list(data = x[["data"]], Y = x[["Y"]], D = x[["D"]],
                    X = x[["X"]], index = x[["index"]],
                    xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
                    axis.adjust = axis.adjust, id = id)
    if (type == "missing") {
      pv.args$main <- if (is.null(main)) "Treatment Status and Missing Data" else main
      pv.args$pre.post <- TRUE # show treatment status
    } else {
      pv.args$main <- if (is.null(main)) "Raw Data" else main
      pv.args$type <- "outcome"
      pv.args$legendOff <- legendOff
    }
    # panelview() prints its plot before returning it. Send that print to
    # a throwaway device, then restore the user's devices, so the plot is
    # drawn once, when the returned object is printed (#58).
    old.dev <- grDevices::dev.cur()
    grDevices::pdf(file = NULL)
    tmp.dev <- grDevices::dev.cur()
    on.exit({
      if (tmp.dev %in% grDevices::dev.list()) grDevices::dev.off(tmp.dev)
      if (old.dev > 1L && old.dev %in% grDevices::dev.list()) grDevices::dev.set(old.dev)
    }, add = TRUE)
    p <- do.call(panelView::panelview, pv.args)
    return(p)
  } else {
    class(x) <- "fect"
    p <- fect::plot.fect(x=x,
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
      legacy.style = legacy.style,
      highlight = highlight,
      highlight.fill = highlight.fill,
        ...)
    return(p)
  }
}


