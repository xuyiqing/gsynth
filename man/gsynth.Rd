\name{gsynth}
\alias{gsynth}
\title{Generalized Synthetic Control Method}
\description{Implements the generalized synthetic control method based
on interactive fixed effect models.}
\usage{gsynth(formula = NULL, data, Y, D, X = NULL, na.rm = FALSE,
       index, weight = NULL, force = "unit", cl = NULL, r = 0,
       lambda = NULL, nlambda = 10, CV = TRUE, criterion = "mspe",
       k = 20, cv.method = "rolling", cv.prop = 0.1, cv.buffer = 1,
       EM = FALSE, estimator = "gsynth", time.component.from = NULL,
       se = FALSE, nboots = 200,
       inference = NULL, parallel = TRUE,
       cores = NULL, tol = 0.001, seed = NULL, min.T0 = 5,
       alpha = 0.05, normalize = FALSE,
       W.est = NULL, W.agg = NULL,
       ci.method = "normal",
       placeboTest = FALSE, placebo.period = NULL)
}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of
  the model to be fitted.}
\item{data}{a data frame (must be with a dichotomous treatment but balanced
  is not required).}
\item{Y}{outcome.}
\item{D}{treatment.}
\item{X}{time-varying covariates.}
\item{na.rm}{a logical flag indicating whether to list-wise delete
  missing data. The algorithm will report an error if missing data
  exist.}
\item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2.}
\item{weight}{a string specifying the weighting variable (if any) to
    estimate the weighted average treatment effect. Default is
    \code{weight = NULL}. Applies to both estimation and aggregation
    roles. For per-role control, use \code{W.est} and \code{W.agg}.}
\item{force}{a string indicating whether unit or time fixed effects will be
    imposed. Must be one of the following,
    "none", "unit", "time", or "two-way". The default is "unit".}
\item{cl}{a string indicator the cluster variable. The default value is
  \code{NULL}. If \code{cl = NULL}, bootstrap will be blocked at unit level
  (only for non-parametric bootstrap).}
\item{r}{an integer specifying the number of factors. If \code{CV = TRUE},
  the cross validation procedure will select the optimal number of
  factors from \code{r} to 5.}
\item{lambda}{a single or sequence of positive numbers specifying the
  hyper-parameter sequence for matrix completion method. If \code{lambda}
  is a sequence and \code{CV = 1}, cross-validation will be performed.}
\item{nlambda}{an integer specifying the length of hyper-parameter sequence
  for matrix completion method. Default is \code{nlambda = 10}.}
\item{CV}{a logical flag indicating whether cross-validation will be
  performed to select the optimal number of factors or hyper-parameter
  in matrix completion algorithm. If \code{r} is not specified, the procedure
   will search through \code{r = 0} to \code{5}.}
\item{criterion}{a string specifying the criteria used for determining the number
  of factors. Choose from \code{c("mspe", "pc")}. "mspe" stands for the mean squared
  prediction error obtained through the loocv procedure, and "pc" stands for a kind
  of information criterion. If \code{criterion = "pc"}, the number of factors that
  minimize "pc" will be selected. Default is \code{criterion = "mspe"}.}
\item{k}{a positive integer specifying the number of cross-validation folds.
  Default is \code{k = 20} (changed from \code{k = 5} in v1.5.0 to align with
  \code{fect} v2.3.0+ defaults).}
\item{cv.method}{a string specifying the cross-validation masking strategy.
  Either \code{"rolling"} (default; rolling-window CV introduced in
  \code{fect} v2.3.0) or \code{"block"} (the pre-v1.5.0 design, retained
  for replication). New in v1.5.0.}
\item{cv.prop}{numeric. Fraction of eligible units sampled per fold under
  rolling CV (proportion of observations masked per round under block CV).
  Default \code{0.1}. New in v1.5.0.}
\item{cv.buffer}{integer. Past-side buffer (in time periods) for
  \code{cv.method = "rolling"}. Default \code{1}. New in v1.5.0.}
\item{EM}{a logical flag indicating whether an Expectation Maximization
  algorithm will be used (Gobillon and Magnac 2016). Soft-deprecated in
  v1.5.0; equivalent to \code{estimator = "ife"}. Removal targeted for
  v2.0.0. For IFE-EM, use \code{fect::fect(method = "ife", ...)} directly.}
\item{estimator}{a string that controls the estimation method:
  \code{"gsynth"} (generalized synthetic control, default), \code{"ife"}
  (interactive fixed effects with EM algorithm), or \code{"mc"} (matrix
  completion method). \code{"ife"} and \code{"mc"} are soft-deprecated in
  v1.5.0; removal targeted for v2.0.0. For IFE-EM or matrix completion,
  use \code{fect::fect(method = "ife", ...)} or
  \code{fect::fect(method = "mc", ...)} directly.}
\item{time.component.from}{a string specifying the sample for factor
  estimation: \code{"nevertreated"} (control units only; the GSC objective)
  or \code{"notyettreated"} (includes treated-unit pre-treatment cells).
  Default \code{NULL} routes by \code{estimator}: \code{"gsynth"} maps to
  \code{"nevertreated"}; \code{"ife"} and \code{"mc"} map to
  \code{"notyettreated"}. New in v1.5.0.}
\item{se}{a logical flag indicating whether uncertainty estimates will
  be produced.}
\item{nboots}{an integer specifying the number of bootstrap
  runs. Ignored if \code{se = FALSE}.}
\item{inference}{a string specifying which type of inferential method
  will be used: \code{"parametric"}, \code{"nonparametric"}, or
  \code{"jackknife"}. \code{"parametric"} is recommended when the number of
  treated units is small and is only valid for \code{estimator = "gsynth"}.
  \strong{Breaking change in v1.5.0:} \code{inference = "parametric"} with
  \code{estimator = "ife"} or \code{"mc"} now errors instead of silently
  falling back to nonparametric bootstrap, aligning with the hard gate in
  \code{fect} v2.2.0+ and the public commitment in Xu (2026), Appendix A.5.}
\item{parallel}{a logical flag indicating whether parallel computing
  will be used in bootstrapping and/or cross-validation. Ignored if
  \code{se = FALSE}.}
\item{cores}{an integer indicating the number of cores to be used in
  parallel computing. If not specified, the algorithm will use the
  maximum number of logical cores of your computer (warning: this
  could prevent you from multi-tasking on your computer).}
\item{tol}{a positive number indicating the tolerance level.}
\item{seed}{an integer that sets the seed in random number
  generation. Ignored if \code{se = FALSE} and \code{r} is specified.}
\item{min.T0}{an integer specifying the minimum value of pre-treatment
  periods. Treated units with pre-treatment periods less than that will
  be removed automatically. This item is important for unbalanced panels.
  If users want to perform cross validation procedure to select the optimal
  number of factors from \code{(r.min, r.max)}, they should set \code{min.T0}
  larger than \code{(r.max+1)} if no individual fixed effects or \code{(r.max+2)}
  otherwise. If there are too few pre-treatment periods among all treated units,
  a smaller value of \code{r.max} is recommended.}
\item{alpha}{a positive number in the range of 0 and 1 specifying significant
 levels for uncertainty estimates. The default value is \code{alpha = 0.05}.}
\item{normalize}{a logic flag indicating whether to scale outcome and
  covariates. Useful for accelerating computing speed when magnitude of
  data is large. The default is \code{normalize=FALSE}.}
\item{W.est}{a string specifying the weighting column for the
  outcome-model fit (per-role weight). Overrides \code{weight} for that
  role. From \code{fect} v2.3.1; new in gsynth v1.5.0.}
\item{W.agg}{a string specifying the weighting column for the ATT
  aggregation step (per-role weight). Overrides \code{weight} for that
  role. From \code{fect} v2.3.1; new in gsynth v1.5.0.}
\item{ci.method}{a string controlling the confidence-interval method
  used for the \code{est.*} slots. \code{"normal"} (default) uses Wald
  intervals (\eqn{\hat{\theta} \pm z_{1-\alpha/2}\,\widehat{\mathrm{SE}}});
  \code{"basic"} uses the reflected-pivot bootstrap interval (Davison and
  Hinkley 1997, Sec. 5.2.1) when \code{inference = "bootstrap"} or
  \code{"jackknife"}. For percentile, BC, and BCa intervals on
  alternative estimands (cumulative ATT, APTT, log-ATT), call
  \code{\link[fect]{estimand}} on the fit object after class-munging
  to \code{"fect"}. New in v1.5.0 (forwards \code{fect}'s
  \code{ci.method}, added in fect v2.4.2).}
\item{placeboTest}{a logical flag. If \code{TRUE}, a contiguous block of
  pre-treatment periods specified by \code{placebo.period} is held out
  from estimation and used as a placebo region. The fitted ATT is
  reported on the placebo periods, allowing the analyst to check
  whether the estimator produces non-zero "effects" in periods where
  no effect should exist. New in v1.5.0; forwards \code{fect}'s
  \code{placeboTest}.}
\item{placebo.period}{an integer vector specifying which event-time
  periods to use as the placebo block. For example,
  \code{placebo.period = c(-2, 0)} treats event-times \eqn{-2} through
  \eqn{0} as placebo. Default \code{NULL} lets \code{fect} pick a
  sensible default (typically the two periods immediately preceding
  treatment). Ignored unless \code{placeboTest = TRUE}.

  \pkg{gsynth} does not expose \code{fect}'s \code{carryoverTest}.
  The carryover test requires panels with treatment reversals, which
  the GSC estimator does not allow. For carryover diagnostics on
  reversal-friendly panels, use \code{\link[fect]{fect}} directly.}
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
  reasonable modeling assumptions. With a built-in cross-validation
  procedure, it avoids specification searches and thus is easy to
  implement. Data must be with a dichotomous treatment.

  \strong{Note (v1.5.0).} The \code{estimator = "ife"} (Gobillon-Magnac
  IFE-EM) and \code{estimator = "mc"} (matrix completion) paths are
  soft-deprecated in this release; each call emits a one-time-per-session
  message. For IFE-EM or matrix completion, use \code{fect::fect()}
  directly with \code{method = "ife"} or \code{method = "mc"} and
  \code{time.component.from = "notyettreated"}. Removal is targeted for
  gsynth v2.0.0. See \code{vignette("02-ife-mc", package = "gsynth")}.
}
\value{
  \item{Y.dat}{a matrix storing data of the outcome variable.}
  \item{Y}{name of the outcome variable.}
  \item{D}{name of the treatment variable.}
  \item{X}{name of the time-varying control variables.}
  \item{index}{name of the unit and time indicators.}
  \item{id}{a vector of unit IDs.}
  \item{time}{a vector of time periods.}
  \item{obs.missing}{a matrix storing status of each unit at each time point.
    \code{0} for missing, \code{1} for control group units, \code{2} for treat group
    units at pre-treatment period, \code{3} for treat group units at post-treatment
    period, and \code{4} for removed treated group units. Useful for unbalanced panel data.}
  \item{id.tr}{a vector of IDs for the treatment units.}
  \item{id.co}{a vector of IDs for the control units.}
  \item{removed.id}{a vector of IDs for units that are removed.}
  \item{D.tr}{a matrix of treatment indicator for the treated unit outcome.}
  \item{I.tr}{a matrix of observation indicator for the treated unit outcome.}
  \item{Y.tr}{data of the treated unit outcome.}
  \item{Y.ct}{predicted counterfactuals for the treated units.}
  \item{Y.co}{data of the control unit outcome.}
  \item{eff}{difference between actual outcome and predicted Y(0).}
  \item{Y.bar}{average values of Y.tr, Y.ct, and Y.co over time.}
  \item{att}{average treatment effect on the treated over time (it is averaged
    based on the timing of the treatment if it is different for each unit).}
  \item{att.avg}{average treatment effect on the treated.}
  \item{force}{user specified \code{force} option.}
  \item{sameT0}{TRUE if the timing of the treatment is the same.}
  \item{T}{the number of time periods.}
  \item{N}{the total number of units.}
  \item{p}{the number of time-varying observables.}
  \item{Ntr}{the number of treated units.}
  \item{Nco}{the number of control units.}
  \item{T0}{a vector that stores the timing of the treatment for balanced panel
    data.}
  \item{tr}{a vector indicating treatment status for each unit.}
  \item{pre}{a matrix indicating the pre-treatment/non-treatment status.}
  \item{post}{a matrix indicating the post-treatment status.}
  \item{r.cv}{the number of factors included in the model -- either
    supplied by users or automatically chosen via cross-validation.}
  \item{lambda.cv}{the optimal hyper-parameter in matrix completion method
  chosen via cross-validation.}
  \item{res.co}{residuals of the control group units.}
  \item{beta}{coefficients of time-varying observables from the
    interactive fixed effect model.}
  \item{sigma2}{the mean squared error of interactive fixed effect model.}
  \item{IC}{the information criterion.}
  \item{PC}{the proposed criterion for determining factor numbers.}
  \item{est.co}{result of the interactive fixed effect model based on
    the control group data.}
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
  \item{wgt.implied}{estimated weights of each of the control group unit for each of the treatment group unit.}
  \item{mu}{estimated ground mean. }
  \item{xi}{estimated time fixed effects.}
  \item{alpha.tr}{estimated unit fixed effects for the treated units.}
  \item{alpha.co}{estimated unit fixed effects for the control units.}
  \item{validX}{a logic value indicating if multicollinearity exists.}
  \item{inference}{a string indicating bootstrap procedure.}
  \item{est.att}{inference for \code{att}.}
  \item{est.att.avg}{inference for \code{att.avg}.}
  \item{est.beta}{inference for \code{beta}.}
  \item{est.ind}{inference for \code{att} of each treated unit.}
  \item{att.avg.boot}{bootstrap results for \code{att.avg}.}
  \item{att.boot}{bootstrap results for \code{att}.}
  \item{beta.boot}{bootstrap results for \code{beta}.}
}
\author{
  Yiqing Xu <yiqingxu@stanford.edu>, Stanford University
}
\references{
  Laurent Gobillon and Thierry Magnac, 2016.
  "Regional Policy Evaluation: Interactive Fixed Effects and Synthetic Controls."
  The Review of Economics and Statistics, July 2016, Vol. 98, No. 3,
  pp. 535--551.

  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25,
  Iss. 1, January 2017, pp. 57-76.

  Athey S, Bayati M, Doudchenko N, et al. Matrix completion methods for causal panel data models[J]. arXiv preprint arXiv:1710.10251, 2017.

  Yiqing Xu. 2026. "A Response to Green and Aronow."

  For more details, see \url{https://yiqingxu.org/packages/gsynth/}.

  For more details about the matrix completion method, see \url{https://github.com/susanathey/MCPanel}.
}
\seealso{
  \code{\link{plot.gsynth}}, \code{\link{print.gsynth}},
  \code{\link[fect]{fect}}.
}
\examples{
library(gsynth)
data(gsynth)
out <- gsynth(Y ~ D + X1 + X2, data = simdata, parallel = FALSE,
              index = c("id","time"), force = "two-way",
              CV = TRUE, r = c(0, 5), se = FALSE)
print(out)
}
