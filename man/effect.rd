\name{effect}
\alias{effect}
\title{Cumulative or Sub-group Treatment Effects}
\description{Calculates cumulative or sub-group treatment effects.

  \strong{Soft-deprecated in v1.5.0}: for cumulative or sub-group
  estimands, use \code{\link[fect]{estimand}} from \pkg{fect} (v2.4.0+)
  on the gsynth fit directly:
  \code{fect::estimand(fit, "att.cumu", ...)}. Each call to
  \code{effect()} emits a one-time-per-session message. Removal not
  before gsynth v2.0.0.
}
\usage{effect(x, cumu = TRUE, period = NULL, id = NULL, plot = FALSE)}
\arguments{
  \item{x}{a \code{\link{gsynth}} object fitted with \code{se = TRUE}.}
  \item{cumu}{a logical flag. If \code{TRUE} (default), the effect at
    each period is the cumulative effect: the running sum of the
    per-period ATTs, from the first period in \code{period} up to that
    period. If \code{FALSE}, the per-period ATTs.}
  \item{id}{a vector of ids specifying a sub-group of treated units that treatment
  effects are to be averaged on. \code{NULL} (default) uses all treated
  units.}
  \item{period}{a two-element numeric vector specifying the range of term during which treatment effects are to be accumulated. If left blank, atts at all post-treatment
  periods will be calculated.}
  \item{plot}{a logical flag indicating whether to plot the cumulative effects.}
}
\value{
  An object of class \code{"fect"}: the fit \code{x} with the two
  elements below added, so \code{print()} and \code{plot()} on it use
  \pkg{fect}'s methods. With \code{plot = TRUE} the plot is also drawn.
  \code{effect()} stops if \code{x} has no bootstrap results (fit it with
  \code{se = TRUE}).
  \item{effect.est.avg}{a numeric vector: the (cumulative) ATT for each
    period in \code{period}.}
  \item{effect.est.att}{a matrix with one row per period in
    \code{period} and columns \code{ATT}, \code{S.E.}, \code{CI.lower},
    \code{CI.upper} and \code{p.value}.}
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


