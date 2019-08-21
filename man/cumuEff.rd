\name{cumuEff}
\alias{cumuEff}
\title{Calculate Cumulative or sub-gr Treatment Effects}
\description{Calculate Cumulative or sub-gr Treatment Effects}
\usage{cumuEff(x, cumu = TRUE, id = NULL, period = NULL)} 
\arguments{
  \item{x}{a \code{\link{gsynth}} object.}
  \item{cumu}{a logical flag indicating whether to calculate cumulative effects or not.}
  \item{id}{a string vector speicfying a sub-group of treated units that treatment 
  effects are to be averaged on. }
  \item{period}{a two-element numeric vector specifying the range of term during which treatment effects are to be accumulated. If left blank, atts at all post-treatment 
  periods will be calculated.}
}
\value{
  \item{catt}{esimated (cumulative) atts.}
  \item{est.catt}{uncertainty estimates for \code{catt}.}
}

\author{
  Yiqing Xu and Licheng Liu 
}
\references{  
  Jushan Bai. 2009. "Panel Data Models with Interactive Fixed
  Effects." Econometrica 77:1229--1279.

  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25, 
  Iss. 1, January 2017, pp. 57-76. Available at: \url{https://doi.org/10.1017/pan.2016.2}.
}
\seealso{
  \code{\link{gsynth}}
}


