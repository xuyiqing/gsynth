# Calculate cumulative average treatment effect
# A wrapper function of fect:::cumuEff
#' @import fect
#' @param x A `fect` object with `method="gsynth"`.
#' @param cumu Boolean, whether to calculate cumulative ATT.
#' @param period c(start, end)
#' @param id ID of the units of interest
#' @param plot Whether to plot the cumulative effects
#' @return Cumulative effects
#' @export
effect <- function(x, cumu=TRUE, period=NULL, id=NULL, plot=FALSE){
  out <- fect:::effect(x=x, cumu=cumu, period=period, id=id, plot=plot)
  return(out)
}
