\name{gsynth}
\alias{gsynth}
\title{Generalized Synthetic Control Method}
\description{Implements the generalized synthetic control method based
on interactive fixed effect models.}
\usage{gsynth(formula, data, na.rm = FALSE, index, force = "unit", r = 0, 
       CV = TRUE, EM = FALSE, se = FALSE, nboots = 200, 
       inference = "parametric", cluster = NULL,
       parallel = TRUE, cores = NULL, tol = 0.001, seed = NULL)  
}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of
  the model to be fitted.}
\item{data}{a data frame (must be a balanced panel with a dichotomous
  treatment).}
\item{na.rm}{a logical flag indicating whether to list-wise delete
  missing data. The algorithm will report an error if missing data
  exist.}
\item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2.}
\item{force}{a string indicating whether unit or time fixed effects will be
    imposed. Must be one of the following,
    "none", "unit", "time", or "two-way". The default is "unit".}
\item{r}{an integer specifying the number of factors. If \code{CV = TRUE}, 
  the cross validation procedure will select the optimal number of
  factors from \code{r} to 5.}
\item{CV}{a logical flag indicating whether cross-validation will be
  performed to select the optimal number of factors. If \code{r} is not
  specified, the procedure will search through \code{r = 0} to \code{5}.}
\item{EM}{a logical flag indicating whether an Expectation Maximization
  algorithm will be used (Gobillon and Magnac 2016).}
\item{se}{a logical flag indicating whether uncertainty estimates will
  be produced.}
\item{nboots}{an integer specifying the number of bootstrap
  runs. Ignored if \code{se = FALSE}.}
\item{inference}{a string specifying which type of inferential method
  will be used, either "parametric" or "non-parametric". "parametric" is
  recommended when the number of treated units is small.}
\item{cluster}{a string indicating the clustering variable to be blocked
  on during bootstrapping.}
\item{parallel}{a logical flag indicating whether parallel computing
  will be used in bootstrapping and/or cross-validation. Ignored if
  \code{se = FALSE}.}
\item{cores}{an integer indicating the number of cores to be used in
  parallel computing. If not specified, the algorithm will use the
  maximum number of logical cores of your computer (warning: this
  could prevent you from multi-tasking on your computer).}
\item{tol}{a positive number indicating the tolerance level.}
\item{seed}{an integer that sets the seed in random number
  generation. Ignored if  \code{se = FALSE} and \code{r} is specified.}
}
\details{
  \code{gsynth} implements the generalized synthetic control method. It
  imputes counterfactuals for each treated unit using control group
  information based on a linear interactive fixed effects model that
  incorporates unit-specific intercepts interacted with time-varying
  coefficients. It generalizes the synthetic control method to the case
  of multiple treated units and variable treatment periods, and improves
  efficiency and interpretability. It allows the treatment to be
  correlated with unobserved unit and time heterogeneities under
  reasonable modelling assumptions. With a built-in cross-validation
  procedure, it avoids specification searches and thus is easy to
  implement. Data must be in form of a balanced panel with a dichotomous
  treatment. 
}
\value{
  \item{Y.dat}{a matrix storing data of the outcome variable.}
  \item{Y}{name of the outcome variable.}
  \item{D}{name of the treatment variable.}
  \item{X}{name of the time-varying control variables.}
  \item{index}{name of the unit and time indicators.}
  \item{id}{a vector of unit IDs.}
  \item{time}{a vector of time periods.}
  \item{id.tr}{a vector of IDs for the treatment units.}
  \item{id.co}{a vector of IDs for the control units.}
  \item{Y.tr}{data of the treated unit outcome.}
  \item{Y.ct}{predicted counterfactuals for the treated units.}
  \item{Y.co}{data of the control unit outcome.}
  \item{eff}{difference between actual outcome and predicted Y(0).}
  \item{Y.bar}{average values of Y.tr, Y.ct, and Y.co over time.}
  \item{att}{average treatment effect on the treated over time (it is averaged based on
    the timing of the treatment if it is different for each unit).}
  \item{att.avg}{average treatment effect on the treated.}
  \item{force}{user specified \code{force} option.}
  \item{DID}{TRUE if the timing of the treatment is the same.}
  \item{T}{the number of time periods.}
  \item{N}{the total number of units.}
  \item{p}{the number of time-varying observables.}
  \item{Ntr}{the number of treated units.}
  \item{Nco}{the number of control units.}
  \item{T0}{a vector that stores the timing of the treatment.}
  \item{tr}{a vector indicating treatment status for each unit.}
  \item{pre}{a matrix indicating the pre-treatment/non-treatment status.}
  \item{r.cv}{the number of factors included in the model -- either
    supplied by users or automatically chosen via cross-validation.}
  \item{res.co}{residuals of the control group units.}
  \item{beta}{coefficients of time-varying observables from the
    interactive fixed effect model.}
  \item{sigma2}{the mean squared error of interactive fixed effect model.}
  \item{IC}{the information criterion.}
  \item{est.co}{result of the interactive fixed effect model based on
    the control group data. An \code{\link{interFE}} object.}
  \item{eff.cnt}{difference between actual outcome and predicted Y(0);
    rearranged based on the timing of the treatment.}
  \item{Y.tr.cnt}{data of the treated unit outcome, rearranged based on the timing of the treatment.}
  \item{Y.ct.cnt}{data of the predicted Y(0), rearranged based on the timing of the treatment.}
  \item{MSPE}{mean squared prediction error of the cross-validated model.}
  \item{CV.out}{result of the cross-validation procedure.}
  \item{niter}{the number of iterations in the estimation of the
    interactive fixed effect model.}
  \item{factor}{estimated time-varying factors.}
  \item{lambda.co}{estimated loadings for the control group.}
  \item{lambda.tr}{estimated loadings for the treatment group.}
  \item{mu}{estimated ground mean. }
  \item{xi}{estimated time fixed effects. }
  \item{alpha.tr}{estimated unit fixed effects for the treated units.}
  \item{alpha.co}{estimated unit fixed effects for the control units.} 
}
\author{
  Yiqing Xu <yiqingxu@ucsd.edu>
}
\references{
  Laurent Gobillon and Thierry Magnac, 2016.
  "Regional Policy Evaluation: Interactive Fixed Effects and Synthetic Controls."
  The Review of Economics and Statistics, July 2016, Vol. 98, No. 3,
  pp. 535--551.

  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis,
  forthcoming. Available at: \url{http://dx.doi.org/10.1017/pan.2016.2}.
  
  For more details, see \url{http://yiqingxu.org/software/gsynth/gsynth_examples.html}. 
}
\seealso{
  \code{\link{plot.gsynth}} and \code{\link{print.gsynth}}
}
\examples{
library(gsynth)
data(gsynth)
out <- gsynth(Y ~ D + X1 + X2, data = simdata, 
              index = c("id","time"), force = "two-way",
              CV = TRUE, r = c(0, 5), se = FALSE)  
}


