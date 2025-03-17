# Calculate cumulative average treatment effect
# A wrapper function of fect:::cumuEff
#' @import fect
#' @param x A `fect` object with `method="gsynth"`.
#' @param cumu Boolean, whether to calculate cumulative ATT.
#' @param period c(start, end)
#' @param id ID of the units of interest
#' @return Cumulative effects
#' @export
cumuEff <- function(x, cumu=TRUE, period, id=NULL){
  library(fect)
  out <- fect:::cumuEff(x=x, cumu=cumu, period=period, id=id)
  return(out)
}
