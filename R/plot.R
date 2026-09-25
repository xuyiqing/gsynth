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
    # Raw data and treatment status/missing data, drawn by panelView.
    # Build the panelView call from what the fit stores (the data and the
    # variable names), not from x$call, which may hold a formula stored in
    # a variable or no formula at all (Y = / D = calls). Exact matching:
    # x$data would fall back to x$data.long if `data` were missing; `Y` and
    # `D` return the first element of that name (the column name).
    pv.data <- x[["data"]]
    pv.Y <- x[["Y"]]
    pv.D <- x[["D"]]
    pv.index <- x[["index"]]
    if (is.null(pv.data) || !is.character(pv.Y) || !is.character(pv.D) ||
        is.null(pv.index)) {
      stop("This fit does not store the data needed for type = \"raw\" or ",
           "\"missing\".", call. = FALSE)
    }
    pv.args <- list(data = pv.data, Y = pv.Y, D = pv.D, X = x[["X"]],
                    index = pv.index,
                    xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
                    axis.adjust = axis.adjust, id = id)
    if (type == "missing") {
      # treatment status, pre/post, and missing cells
      pv.args$main <- if (is.null(main)) {
        "Treatment Status and Missing Data"
      } else {
        main
      }
      pv.args$pre.post <- TRUE
    } else {
      # outcome paths
      pv.args$main <- if (is.null(main)) "Raw Data" else main
      pv.args$type <- "outcome"
      pv.args$legendOff <- legendOff
    }
    # panelview() prints its plot before returning it. Send that print to a
    # throwaway device and return the plot, as for the other types: it is
    # drawn when printed (e.g. auto-printed at the console), once.
    old.dev <- grDevices::dev.cur()
    grDevices::pdf(file = NULL)
    tmp.dev <- grDevices::dev.cur()
    on.exit({
      if (tmp.dev %in% grDevices::dev.list()) grDevices::dev.off(tmp.dev)
      if (old.dev > 1L && old.dev %in% grDevices::dev.list()) {
        grDevices::dev.set(old.dev)
      }
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


