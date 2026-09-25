\name{gsynth}
\alias{gsynth}
\title{Generalized Synthetic Control Method}
\description{Implements the generalized synthetic control method based
on interactive fixed effect models.}
\usage{gsynth(formula = NULL, data, Y, D, X = NULL, na.rm = FALSE,
       index, weight = NULL, force = "unit", cl = NULL, r = 0,
       lambda = NULL, nlambda = 10, CV = TRUE, criterion = "mspe",
       k = 20, cv.method = "rolling", cv.prop = 0.1, cv.buffer = 1,
       cv.nobs = 3, cv.donut = 1, cv.rule = "1se",
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
\item{formula}{an object of class \code{"formula"}: \code{Y ~ D + X1 + X2},
  with the outcome on the left, then the treatment, then the covariates,
  all as plain column names of \code{data}. Transformations such as
  \code{log(Y)}, \code{factor()} terms and interactions are not allowed
  (\pkg{fect} 2.4.6 stops with a message): create the column first, for
  example \code{data$logY <- log(data$Y)}, and make dummy columns with
  \code{model.matrix()}. Intercept terms (\code{+ 0}, \code{- 1}) are
  accepted and ignored, because the model has its own fixed effects.}
\item{data}{a data frame (must be with a dichotomous treatment but balanced
  is not required).}
\item{Y}{name of the outcome column. Used when \code{formula} is not
  given.}
\item{D}{name of the treatment column (1 = treated, 0 = not treated). Used
  when \code{formula} is not given.}
\item{X}{names of the covariate columns. Used only when \code{formula} is
  not given: giving \code{X} together with a formula stops, so put the
  covariates in one place. Covariates, whether given here or in
  \code{formula}, must be numeric: logical columns count as 0/1, and factor
  or character columns stop (make numeric dummy columns first). A covariate
  that the fixed effects absorb (for example, one that is constant over
  time within units under unit fixed effects) or that is an exact linear
  combination of other covariates is dropped with a warning, and its
  coefficient is reported as \code{NA} (\pkg{fect} 2.4.6).}
\item{na.rm}{a logical flag indicating whether to list-wise delete
  missing data. The algorithm will report an error if missing data
  exist.}
\item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2. A factor time index is used
    in the order of its levels; a character time index must contain
    numbers (\pkg{fect} 2.4.6).}
\item{weight}{name of a numeric, non-negative weight column used to
    \emph{average} the treatment effects: \code{att}, \code{att.avg} and
    their uncertainty estimates. This is the same role as \code{W.agg}, and
    the same as in gsynth 1.2.1. It does not weight the model fit; use
    \code{W.est} for that. If \code{weight} and \code{W.agg} name different
    columns, gsynth stops. Rows whose weight is exactly 0 are dropped from
    the data. Default \code{weight = NULL} (no weights).}
\item{force}{a string indicating whether unit or time fixed effects will be
    imposed. Must be one of the following,
    "none", "unit", "time", or "two-way". The default is "unit".}
\item{cl}{name of the cluster column for the nonparametric bootstrap
  (\code{inference = "nonparametric"}), which then resamples whole
  clusters of units instead of single units. It must be a single column
  that is constant within each unit, with at least two clusters
  (\pkg{fect} 2.4.6). Parametric inference and the jackknife ignore it,
  with a warning. Default \code{NULL}: resample single units.}
\item{r}{the number of latent factors \emph{on top of} the fixed effects in
  \code{force}: with \code{force = "two-way"}, \code{r = 0} is the two-way
  fixed-effects model. Give one number or a range \code{c(r.min, r.max)}.
  With \code{CV = TRUE} and one number, cross-validation chooses from
  \code{r} to 5 (only \code{r} itself when \code{r >= 5}); with a range, it
  chooses from \code{r.min} to \code{r.max}. \pkg{fect} lowers the top of
  the search when there are too few pre-treatment periods or control
  units. With \code{CV = FALSE} the model uses \code{r} (the smaller value
  of a range). So the default call (\code{r = 0}, \code{CV = TRUE})
  searches 0 to 5. Not used by \code{estimator = "mc"}.}
\item{lambda}{a single or sequence of positive numbers specifying the
  hyper-parameter sequence for matrix completion method. If \code{lambda}
  is a sequence and \code{CV = 1}, cross-validation will be performed.}
\item{nlambda}{an integer specifying the length of hyper-parameter sequence
  for matrix completion method. Default is \code{nlambda = 10}.}
\item{CV}{a logical flag: whether to choose the number of factors (for
  \code{estimator = "mc"}, the penalty \code{lambda}) by cross-validation.
  Default \code{TRUE}. See \code{r} for the range searched. \code{criterion}
  and \code{cv.rule} set how the number is picked; by default it is the
  smallest \code{r} whose error is within one standard error of the
  lowest.}
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
\item{cv.nobs}{integer. Under \code{cv.method = "block"}, the number of
  consecutive periods held out in each run. Default \code{3}. New in v1.5.0.}
\item{cv.donut}{integer. Under \code{cv.method = "block"}, the number of
  periods at each end of a held-out run that are not scored; \code{2 * cv.donut}
  must be smaller than \code{cv.nobs}. Default \code{1}. New in v1.5.0.}
\item{cv.rule}{a string: the rule for choosing the number of factors from
  the cross-validation errors. \code{"1se"} (default) picks the smallest
  \code{r} whose error is within one fold standard error of the lowest;
  \code{"min"} picks the lowest error; \code{"1pct"} picks the smallest
  \code{r} within 1\% of the lowest. See \code{fect::fect}. For
  \code{estimator = "gsynth"} the rule takes effect with \pkg{fect} 2.4.6
  or later. New in v1.5.0.}
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
  runs. Ignored if \code{se = FALSE}. The jackknife ignores it and uses one
  run per unit.}
\item{inference}{the inference method: \code{"parametric"},
  \code{"nonparametric"} (alias \code{"bootstrap"}: resamples units, or
  the clusters given by \code{cl}), or \code{"jackknife"} (leaves out one
  unit at a time). \code{NULL} (default) uses \code{"parametric"} for
  \code{estimator = "gsynth"} and \code{"bootstrap"} for \code{"ife"} and
  \code{"mc"}. \code{"parametric"} is recommended when the number of
  treated units is small. The method used is stored in the output as
  \code{vartype}.
  \strong{Breaking change in v1.5.0:} with \code{se = TRUE},
  \code{inference = "parametric"} is not available for
  \code{estimator = "ife"} or \code{"mc"} (or \code{EM = TRUE}): gsynth
  stops with a message instead of silently falling back to the
  nonparametric bootstrap, aligning with the hard gate in \code{fect}
  v2.2.0+ and the public commitment in Xu (2026), Appendix A.5. The one
  exception is \code{estimator = "ife"} with
  \code{time.component.from = "nevertreated"}, which fits the same model
  as \code{estimator = "gsynth"}.}
\item{parallel}{\code{TRUE} (default) runs the bootstrap in parallel R
  worker processes, and cross-validation too on large panels.
  \code{FALSE} runs everything in the current R session. \code{"boot"} or
  \code{"cv"} runs only that step in parallel.}
\item{cores}{the number of parallel R worker processes (not CPU cores).
  \code{NULL} (default) uses the number of available cores minus 2, at
  most 8. With a multithreaded BLAS (OpenBLAS, MKL, Accelerate), each
  worker may also start several threads; set, for example,
  \code{OPENBLAS_NUM_THREADS=1} before starting R to avoid oversubscribing
  the machine.}
\item{tol}{a positive number indicating the tolerance level.}
\item{seed}{an integer seed for the random numbers used by
  cross-validation and the bootstrap, so that repeating a call gives the
  same results. \code{NULL} (default) uses the current state of R's random
  number generator (see \code{\link{set.seed}}).}
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
\item{W.est}{name of the weight column for the model fit. If \code{W.est}
  and \code{W.agg} (or \code{weight}) are both given, they must name the
  same column (a \pkg{fect} limitation). From \code{fect} v2.3.1; new in
  gsynth v1.5.0.}
\item{W.agg}{name of the weight column for averaging the treatment effects
  (\code{att}, \code{att.avg} and their uncertainty estimates). The same as
  \code{weight}. From \code{fect} v2.3.1; new in gsynth v1.5.0.}
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
  gsynth v2.0.0. See
  \url{https://yiqingxu.org/packages/gsynth/02-ife-mc.html} for migration
  recipes.
}
\value{
  An object of class \code{"gsynth"}: a list with the elements below (the
  output of \code{fect::fect()}, plus \code{data} and \code{call}). T is
  the number of periods and N the number of units used. T x N matrices
  have periods in rows (in the order of \code{rawtime}) and units in
  columns (in the order of \code{id}). Some elements exist only for some
  fits, as noted.

  \emph{Data and identifiers}
  \item{Y.dat, D.dat, I.dat}{T x N matrices of the outcome, the treatment
    (1 = treated) and the observation indicator (1 = observed, 0 =
    missing). Missing cells hold 0 in \code{Y.dat} and \code{D.dat}.}
  \item{Y, D}{names of the outcome and treatment columns. The list also
    holds a second element named \code{Y} and a second named \code{D}:
    copies of \code{Y.dat} and \code{D.dat}. \code{x$Y} and
    \code{x[["Y"]]} return the name.}
  \item{X}{names of the covariates given (\code{NULL} if none). A covariate
    that \pkg{fect} dropped keeps its name here. For \code{"ife"} and
    \code{"mc"} fits the list also holds a second \code{X}: the
    T x N x p array of covariate values.}
  \item{W}{name of the weight column (\code{NULL} if none). With weights,
    the list also holds a second \code{W}: the T x N matrix of weights.}
  \item{W.est.col, W.agg.col}{the weight column used in the model fit and
    the one used in averaging (\code{NULL} if none).}
  \item{W.in.fit, W.in.agg}{\code{TRUE} if weights enter the model fit and
    the averages. With \code{weight}, \code{W.in.fit} is \code{FALSE} and
    \code{W.in.agg} is \code{TRUE}.}
  \item{index}{names of the unit and time columns.}
  \item{id}{ids of the units used (the columns of the T x N matrices).}
  \item{rawtime}{the time periods (the rows of the T x N matrices).}
  \item{data}{the data frame passed to \code{gsynth()}. The
    \code{"raw"} and \code{"missing"} plots use it.}
  \item{data.long}{the long-form data used by \pkg{fect} (unit, time,
    outcome, treatment).}
  \item{call}{the matched call, as typed. The inference method used is in
    \code{vartype}.}
  \item{formula}{the model formula; \code{NULL} when \code{Y} and \code{D}
    were given as names.}

  \emph{Sample structure}
  \item{T, N, Ntr, Nco, p}{the numbers of periods, units, treated units,
    control units, and covariates in the model (after any were dropped).}
  \item{tr, co}{column positions of the treated and the control units in
    the T x N matrices. The treated units' ids are \code{id[tr]}.}
  \item{T.on}{T x N matrix of event time: 1 in the first treated period,
    0 in the last period before treatment, negative numbers before that;
    \code{NA} for never-treated units.}
  \item{obs.missing}{T x N matrix with the status of every unit in the
    data, including removed ones (columns named by unit id, rows by
    period): 1 = treated (under treatment); 2 = not under treatment
    (control units, and treated units before treatment); 3 = missing;
    4 = unit removed (a treated unit with fewer than \code{min.T0}
    pre-treatment periods); 5 = held out as a placebo period.}
  \item{sample}{T x N logical matrix, laid out like \code{obs.missing}:
    \code{TRUE} for the cells in the estimation sample (codes 1, 2 and
    5).}
  \item{unit.type}{for each unit used: 1 = never treated, 2 = treated in
    every observed period, 3 = treatment starts during the panel.}
  \item{remove.id}{ids of the treated units dropped because they have
    fewer than \code{min.T0} pre-treatment periods. Present only when
    units are removed (with \pkg{fect} 2.4.6, whenever units are
    removed). It replaces \code{removed.id}.}
  \item{I, II}{T x N matrices: the observed cells (the same as
    \code{I.dat}); the cells used to fit the model (observed, not under
    treatment, and not held out as placebo cells).}
  \item{hasRevs}{1 if the treatment ever switches off; always 0 for a fit
    that gsynth accepts.}

  \emph{Model}
  \item{method}{\code{"gsynth"}, \code{"ife"} or \code{"mc"}.}
  \item{time.component.from}{\code{"nevertreated"} (factors estimated from
    the control units only) or \code{"notyettreated"} (also from the
    treated units' pre-treatment cells).}
  \item{em}{\code{TRUE} if \pkg{fect} filled cells by EM while fitting
    (\code{"ife"} and \code{"mc"} fits).}
  \item{force}{the fixed effects, coded 0 = none, 1 = unit, 2 = time,
    3 = two-way.}
  \item{r.cv}{the number of factors used, chosen by cross-validation or
    given (\code{"gsynth"} and \code{"ife"} fits).}
  \item{CV.out}{the cross-validation table of an \code{estimator =
    "gsynth"} fit with \code{CV = TRUE}: one row per candidate \code{r},
    with \code{sigma2}, \code{IC}, \code{PC}, \code{MSPE} and other scores.
    It replaces the old \code{MSPE} element. For \code{"ife"} the table is
    \code{CV.out.ife}; for \code{"mc"} it is \code{CV.out.mc}, with one row
    per candidate \code{lambda}.}
  \item{lambda.cv, lambda.seq, lambda.norm, eigen.all}{\code{"mc"} fits:
    the chosen penalty; the candidate penalties; the chosen penalty
    relative to the largest singular value; the singular values.}
  \item{beta}{covariate coefficients: one row per covariate requested,
    named by covariate. \code{NA} for a covariate dropped because the
    fixed effects absorb it or because it is collinear with others
    (\pkg{fect} 2.4.6). A single \code{NA} when there are no
    covariates.}
  \item{mu, xi, alpha.tr, alpha.co}{the grand mean; the time fixed effects
    (T x 1, only with time fixed effects); the unit fixed effects of the
    treated and of the control units (only with unit fixed effects).}
  \item{alpha}{unit fixed effects of all units (\code{"ife"} and
    \code{"mc"} fits).}
  \item{factor}{T x r matrix of the estimated factors (when
    \code{r.cv > 0}).}
  \item{lambda.co, lambda.tr}{factor loadings of the control and of the
    treated units (when \code{r.cv > 0}).}
  \item{lambda}{factor loadings of all units (\code{"ife"} fits).}
  \item{wgt.implied}{Nco x Ntr matrix of implied weights, for
    \code{estimator = "gsynth"} fits with \code{r.cv > 0}. Rows are named
    by control-unit id and columns by treated-unit id (\pkg{fect} 2.4.6).
    Column j holds the weights of the control units whose loadings rebuild
    treated unit j's loadings (\code{t(lambda.co) \%*\% wgt.implied}
    equals \code{t(lambda.tr)}), so the weighted control outcomes rebuild
    the factor part of unit j's counterfactual. The weights can be
    negative and need not sum to one. \code{NULL} when \code{r.cv = 0}
    and for \code{"ife"} and \code{"mc"} fits. The values differ from
    those of earlier gsynth versions, which used another formula.}
  \item{loading.bound}{\code{"none"}: \pkg{fect}'s bounded-loading option,
    which gsynth does not expose (\code{"gsynth"} fits).}
  \item{sigma2, sigma2.fect}{the error variance of the fitted model; the
    error variance of the model without factors (used for the
    equivalence-test bound).}
  \item{IC, PC}{the information criterion and the "PC" criterion of the
    fitted model (\code{"gsynth"} and \code{"ife"} fits).}
  \item{rmse}{root mean squared error of the fit to the treated units'
    pre-treatment outcomes.}
  \item{niter}{the number of iterations used.}
  \item{validX}{1 if at least one covariate is left in the model after
    \pkg{fect} drops the covariates that the fixed effects absorb or that
    are collinear; 0 if none is left or none was given. It is not a
    collinearity flag.}
  \item{validF}{1 if the model has at least one factor, else 0.}
  \item{est}{the raw output of the model fit. For \code{"gsynth"} fits:
    \code{mu}, \code{beta}, \code{factor}, \code{lambda}, \code{VNT},
    \code{niter}, \code{alpha}, \code{xi}, \code{residuals},
    \code{sigma2}, \code{IC}, \code{PC} and \code{validX} of the fit on
    the control units.}

  \emph{Counterfactuals and effects}
  \item{Y.ct}{T x N matrix of the estimated Y(0) for \emph{all} units. At
    the treated units' treated cells (and at placebo cells held out by
    \code{placeboTest}) it is the imputed counterfactual; at the cells used
    to fit the model (control units, and treated units before treatment)
    it is the model's fitted value. \code{NA} at missing cells.}
  \item{eff}{T x N matrix \code{Y.dat - Y.ct} for \emph{all} units
    (columns named by unit id). At treated cells (\code{D.dat == 1}) it is
    the unit-period treatment effect; at the cells used to fit the model it
    is the model residual (see \code{res}). \code{NA} at missing cells. Only the treated cells enter \code{att}, \code{att.avg}
    and the \code{est.*} elements.}
  \item{eff.tr}{the treated units' columns of \code{eff} (\code{"gsynth"}
    fits).}
  \item{Y.tr.cnt}{T x Ntr matrix of the treated units' counterfactual
    Y(0), the same as \code{Y.ct[, tr]}, in calendar time (\code{"gsynth"}
    fits). This is not the event-time average that gsynth 1.2.1 and
    earlier stored under this name.}
  \item{Y.ct.cnt}{T x Nco matrix of the control units' fitted values, the
    same as \code{Y.ct[, co]}, in calendar time (\code{"gsynth"} fits).}
  \item{res, res.full}{\code{res}: the model residuals at cells not under
    treatment. For \code{"gsynth"} fits it is a T x Nco matrix for the
    control units (0 at missing cells); for \code{"ife"} and \code{"mc"}
    fits a T x N matrix (\code{NA} at treated cells). It replaces
    \code{res.co}. \code{res.full} (\code{"gsynth"} fits): T x N matrix of
    residuals at the cells used to fit the model, \code{NA} elsewhere.}
  \item{att, time, count}{the ATT by event time; the event times (relative
    to the treatment onset, 1 = first treated period; calendar time is in
    \code{rawtime}); the number of treated observations at each event
    time.}
  \item{att.avg}{the average treatment effect on the treated over all
    treated cells: the mean of \code{eff} over cells with
    \code{D.dat == 1}, weighted by \code{weight} (or \code{W.agg}) if
    given. \code{att} is weighted the same way.}
  \item{att.avg.unit}{the average of each treated unit's own average
    effect, so that every treated unit counts once (not weighted).}
  \item{eff.calendar, eff.calendar.fit, N.calendar}{the ATT by calendar
    period (columns \code{ATT-calendar} and \code{count}); the same with a
    loess-smoothed column (\code{Fitted}); the number of treated
    observations in each calendar period.}
  \item{eff.pre, eff.pre.equiv, pre.sd}{the treated units' pre-treatment
    residuals in long form (columns \code{eff}, \code{period},
    \code{unit}); the same from the model without factors; their standard
    deviation by event time.}
  \item{pre.periods, proportion, tost.threshold}{settings of the
    pre-trend tests: the event times used, the minimum share of treated
    units, and the equivalence bound.}
  \item{placeboTest, placebo.period}{the placebo settings used.}
  \item{att.placebo}{the ATT over the placebo periods (with
    \code{placeboTest = TRUE}).}

  \emph{Uncertainty estimates} (only with \code{se = TRUE})
  \item{vartype}{the inference method used: \code{"parametric"},
    \code{"bootstrap"} or \code{"jackknife"}. It replaces the old
    \code{inference} element.}
  \item{para.error}{for parametric inference, how the errors were drawn
    (\code{"empirical"}, \code{"ar"} or \code{"wild"}); \code{NULL}
    otherwise.}
  \item{est.avg}{one-row matrix for \code{att.avg} with columns
    \code{ATT.avg}, \code{S.E.}, \code{CI.lower}, \code{CI.upper} and
    \code{p.value}. It replaces \code{est.att.avg}.}
  \item{est.avg.unit}{the same for \code{att.avg.unit}.}
  \item{est.att}{\code{att} by event time (one row per event time) with
    columns \code{ATT}, \code{S.E.}, \code{CI.lower}, \code{CI.upper},
    \code{p.value} and \code{count}.}
  \item{est.att90, att.bound}{\code{est.att} with 90\% intervals; the 90\%
    interval alone (used by the equivalence test).}
  \item{est.beta}{\code{beta} with columns \code{Coef}, \code{S.E.},
    \code{CI.lower}, \code{CI.upper} and \code{p.value}.}
  \item{est.eff.calendar, est.eff.calendar.fit}{uncertainty estimates for
    \code{eff.calendar} and \code{eff.calendar.fit}.}
  \item{est.placebo}{uncertainty estimates for \code{att.placebo} (with
    \code{placeboTest = TRUE}).}
  \item{att.avg.boot, att.avg.unit.boot, att.boot, att.boot.original,
    att.count.boot, beta.boot, att.placebo.boot}{the bootstrap (or
    jackknife) draws of \code{att.avg}, \code{att.avg.unit}, \code{att}
    (one column per draw; \code{att.boot.original} also keeps failed
    draws), the counts, \code{beta} and \code{att.placebo}. Under
    parametric inference the draws are placebo draws centered at zero.}
  \item{att.vcov}{the covariance matrix of \code{att} across event
    times.}
  \item{eff.boot, D.boot, I.boot, colnames.boot}{for each draw, the
    unit-level effects, treatment and observation matrices (arrays with one
    slice per draw), and the positions in \code{id} of the units drawn
    (\code{colnames.boot[[b]]}, whose length is draw b's width). Under the
    cluster bootstrap the draws differ in width, and slices are padded with
    \code{NA} to the widest draw (\pkg{fect} 2.4.6). Used by
    \code{\link{effect}} and \code{fect::estimand()}.}
  \item{Y.avg}{the average treated and counterfactual outcomes by period
    with confidence bands (bootstrap and jackknife only; used by the
    counterfactual plot).}
  \item{test.out}{the pre-trend tests: \code{f.stat}, \code{f.p},
    \code{f.threshold}, \code{f.equiv.p}, \code{df1}, \code{df2},
    \code{N_bar}, \code{tost.equiv.p} and \code{tost.threshold}.}

  \emph{Internal} (used by \pkg{fect}; not needed for analysis)
  \item{Y.ct.full, T.on.balance, G, T.off, balance.period,
    obs.missing.balance, carryoverTest, carryover.period, carryover.rm,
    loo, binary, eff_debias, gamma, kappa, calendar.enp, rmCV, estCV}{
    internal elements of the \pkg{fect} output.}
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
  \code{\link{effect}}, \code{\link[fect]{fect}}.
}
\examples{
library(gsynth)
data(gsynth)
out <- gsynth(Y ~ D + X1 + X2, data = simdata, parallel = FALSE,
              index = c("id","time"), force = "two-way",
              CV = TRUE, r = c(0, 5), se = FALSE)
print(out)
}
