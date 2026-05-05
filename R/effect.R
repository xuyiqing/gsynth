# Calculate cumulative average treatment effect or group effects
# A wrapper function of fect:::cumuEff
#
# Soft-deprecated in v1.5.0. For new code, class-munge the gsynth fit
# to "fect" and call fect::estimand(fit, "att.cumu", ...). Removal not
# before gsynth v2.0.0.
#' @import fect
#' @param x A `fect` object with `method="gsynth"`.
#' @param cumu Boolean, whether to calculate cumulative ATT.
#' @param period c(start, end)
#' @param id ID of the units of interest
#' @param plot Whether to plot the cumulative effects
#' @return Cumulative effects
#' @export
effect <- function(x, cumu=TRUE, period=NULL, id=NULL, plot=FALSE){
  .deprecation_message(
    key = "effect",
    msg = paste0(
      "effect(): soft-deprecated in v1.5.0. ",
      "For new code, class-munge the gsynth fit to \"fect\" ",
      "and call fect::estimand(fit, \"att.cumu\", ...) directly. ",
      "Removal not before gsynth v2.0.0."
    )
  )
  out <- fect::effect(x=x, cumu=cumu, period=period, id=id, plot=plot)
  return(out)
}
