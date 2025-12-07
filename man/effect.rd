\name{effect}
\alias{effect}
\title{Cumulative or Sub-group Treatment Effects}
\description{Calculates cumulative or sub-group treatment effects}
\usage{effect(x, cumu = TRUE, period = NULL, id = NULL, plot = FALSE)}
\arguments{
  \item{x}{a \code{\link{gsynth}} object.}
  \item{cumu}{a logical flag indicating whether to calculate cumulative effects or not.}
  \item{id}{a string vector speicfying a sub-group of treated units that treatment
  effects are to be averaged on. }
  \item{period}{a two-element numeric vector specifying the range of term during which treatment effects are to be accumulated. If left blank, atts at all post-treatment
  periods will be calculated.}
  \item{plot}{a logical flag indicating whether to plot the cumulative effects.}
}
\value{
  \item{catt}{esimated (cumulative) atts.}
  \item{est.catt}{uncertainty estimates for \code{catt}.}
}

\author{
  Yiqing Xu <yiqingxu@stanfprd.edu>, Stanford University
}
\references{
  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25,
  Iss. 1, January 2017, pp. 57-76.
}
\seealso{
  \code{\link{gsynth}}
}


