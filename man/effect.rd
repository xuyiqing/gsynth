\name{effect}
\alias{effect}
\title{Cumulative or Sub-group Treatment Effects}
\description{Calculates cumulative or sub-group treatment effects.

  \strong{Soft-deprecated in v1.5.0}: for cumulative or sub-group
  estimands, use \code{\link[fect]{estimand}} from \pkg{fect} (v2.4.0+).
  After class-munging the gsynth fit to \code{"fect"}, call
  \code{fect::estimand(fit, "att.cumu", ...)}. Each call to
  \code{effect()} emits a one-time-per-session message. Removal not
  before gsynth v2.0.0.
}
\usage{effect(x, cumu = TRUE, period = NULL, id = NULL, plot = FALSE)}
\arguments{
  \item{x}{a \code{\link{gsynth}} object fitted with \code{se = TRUE}.
    \code{effect()} stops if \code{x} has no bootstrap or jackknife
    results.}
  \item{cumu}{a logical flag. \code{TRUE} (default) gives cumulative
    effects: for each period, the average effect over the periods from the
    start of \code{period} up to that period, times the number of those
    periods. \code{FALSE} gives the average effect of each period.}
  \item{id}{a vector of treated-unit ids specifying a sub-group of treated
    units whose treatment effects are averaged. \code{NULL} (default) uses
    all treated units.}
  \item{period}{a two-element numeric vector: the first and the last
    event time (1 = the first treated period) over which treatment effects
    are accumulated. If left blank, all post-treatment periods are used.}
  \item{plot}{a logical flag indicating whether to plot the cumulative
    effects.}
}
\value{
  An object of class \code{"fect"}: the fit \code{x} with two added
  elements.
  \item{effect.est.avg}{the (cumulative) average treatment effect on the
    treated for each period in \code{period}.}
  \item{effect.est.att}{a matrix with one row per period and columns
    \code{ATT}, \code{S.E.}, \code{CI.lower}, \code{CI.upper} and
    \code{p.value}. The intervals follow the inference method stored in
    \code{x$vartype} (\pkg{fect} 2.4.6): the estimate plus or minus a t
    quantile times the standard error for parametric and jackknife fits,
    and quantiles of the bootstrap draws for bootstrap fits.}
  With \code{plot = TRUE} the plot is also drawn.
}
\author{
  Yiqing Xu <yiqingxu@stanford.edu>, Stanford University
}
\references{
  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25,
  Iss. 1, January 2017, pp. 57-76.
}
\seealso{
  \code{\link{gsynth}}, \code{\link[fect]{estimand}}.
}
