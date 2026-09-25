\name{print.gsynth}
\alias{print.gsynth}
\title{Print Results}
\description{Print results of the generalized synthetic control method.}
\usage{\method{print}{gsynth}(x, \dots)}
\arguments{
  \item{x}{a \code{\link{gsynth}} object.}
  \item{\dots}{not used.}
}
\value{
  The object \code{x}, invisibly. The function is called for its side
  effect: it prints the call as typed, the average treatment effect on the
  treated, the effects by period and, with covariates, their coefficients,
  with uncertainty estimates when the fit has them.
}
\author{
  Yiqing Xu <yiqingxu@stanford.edu>, Stanford University
}
\references{Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25,
  Iss. 1, January 2017, pp. 57-76.

  For more details, see \url{https://yiqingxu.org/packages/gsynth/}.
}
\seealso{
  \code{\link{gsynth}} and \code{\link{plot.gsynth}}
}
