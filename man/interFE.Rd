\name{interFE}
\alias{interFE}
\title{Interactive Fixed Effects Models}
\description{Estimating interactive fixed effect models.}
\usage{interFE(formula = NULL, data, Y, X, index, r = 0, force = "none",
         se = TRUE, nboots = 500, seed = NULL, normalize = FALSE)
}
\arguments{
  \item{formula}{an object of class "formula": a symbolic description of the model to be fitted. }
  \item{data}{a data frame (must be with a dichotomous treatment but balanced
    is not reqiored).}
  \item{Y}{outcome.}
  \item{X}{time-varying covariates.}
  \item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2.}
  \item{r}{an integer specifying the number of factors.}
  \item{force}{a string indicating whether unit or time fixed effects will be
    imposed. Must be one of the following,
    "none", "unit", "time", or "two-way". The default is "unit".}
  \item{se}{a logical flag indicating whether uncertainty estimates will
    be produced via bootstrapping.}
  \item{nboots}{an integer specifying the number of bootstrap
    runs. Ignored if \code{se = FALSE}.}
  \item{seed}{an integer that sets the seed in random number
    generation. Ignored if  \code{se = FALSE} and \code{r} is specified.}
  \item{normalize}{a logic flag indicating whether to scale outcome and 
    covariates. Useful for accelerating computing speed when magnitude of data is large.The default is \code{normalize=FALSE}.}
}
\details{
  \code{interFE} estimates interactive fixed effect models proposed by
  Bai (2009).
}
\value{
  \item{beta}{estimated coefficients.}
  \item{mu}{estimated grand mean.}
  \item{factor}{estimated factors.}
  \item{lambda}{estimated factor loadings.}
  \item{VNT}{a diagonal matrix that consists of the r eigenvalues.}
  \item{niter}{the number of iteration before convergence.}
  \item{alpha}{estimated unit fixed effect (if \code{force} is "unit" or "two-way").}
  \item{xi}{estimated time fixed effect (if \code{force} is "time" or "two-way").}
  \item{residuals}{residuals of the estimated interactive fixed effect model.}
  \item{sigma2}{mean squared error of the residuals.}
  \item{IC}{the information criterion.}
  \item{ValidX}{a bool value specifying whether there are valid covariates.}
  \item{dat.Y}{a matrix storing data of the outcome variable.}
  \item{dat.X}{an array storing data of the independent variables.}
  \item{Y}{name of the outcome variable.}
  \item{X}{name of the time-varying control variables.}
  \item{index}{name of the unit and time indicators.}
  \item{est.table}{a table of the estimation results.}
  \item{est.boot}{a matrix storing results from bootstraps.}
}
\author{
  Yiqing Xu <yiqingxu@ucsd.edu>

  Licheng Liu <liulch.16@sem.tsinghua.edu.cn>
}
\references{Jushan Bai. 2009. "Panel Data Models with Interactive Fixed
  Effects." Econometrica 77:1229--1279.}
\seealso{
  \code{\link{print.interFE}} and \code{\link{gsynth}}
}
\examples{
library(gsynth)
data(gsynth)
d <- simdata[-(1:150),] # remove the treated units
out <- interFE(Y ~ X1 + X2, data = d, index=c("id","time"),
               r = 2, force = "two-way", nboots = 50)
}


