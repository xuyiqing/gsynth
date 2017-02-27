\name{plot.gsynth}
\alias{plot.gsynth}
\title{Plotting}
\description{Visualizes estimation results of the generalized synthetic control method.}
\usage{plot.gsynth(out, type = "gap", xlim = NULL, ylim = NULL,
            xlab = NULL, ylab = NULL, legendOff = FALSE, raw = "band",
            main = NULL, nfactors = NULL, id = NULL)
}
\arguments{
  \item{out}{a \code{\link{gsynth}} object.}
  \item{type}{a string that specifies the type of the plot. Must be one
    of the following: "gap" (plotting the average treatment effect on
    the treated; "raw" (plotting the raw data); "counterfactual"
    (plotting predicted Y(0)'s); "factors" (plotting estimated factors); and
    "loadings" (plotting the distribution of estimated factor
    loadings).}
  \item{xlim}{a two-element numeric vector specifying the range of x-axis.}
  \item{ylim}{a two-element numeric vector specifying the range of y-axis.}
  \item{xlab}{a string indicating the label of the x-axis.}
  \item{ylab}{a string indicating the label of the y-axis.}
  \item{legendOff}{a logical flag controlling whether to show the legend.}
  \item{raw}{a string indicating whether or how raw data for the outcome
    variable will be shown
    in the "counterfactual" plot. Ignored if \code{type} is not
    "counterfactual". Must be one of the following: "none" (not showing
    the raw data); "band" (showing the middle 90 percentiles of the raw
    data); and "all" (showing the raw data as they are).}
  \item{main}{a string that controls the title of the plot. If not
    supplied, no title will be shown.}
  \item{nfactors}{a positive integer that specifies the number of
    factors to be shown. The maximum number if 4. Ignored if \code{type}
  is not "factors"}
  \item{id}{a unit identifier of which the predicted counterfactual or
    the difference between actual and predicted counterfactual is to be shown. Ignored if \code{type} is neither
   "counterfactual" nor "gap".}
}
\details{
  \code{plot.gsynth} visualizes the raw data used by, or estimation
  results obtained from, the generalized synthetic control method.
}
\author{
  Yiqing Xu <yiqingxu@ucsd.edu>
}
\references{Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis,
  forthcoming. Available at: \url{http://dx.doi.org/10.1017/pan.2016.2}.

  See \url{http://yiqingxu.org/software/gsynth/gsynth_examples.html} for
  more detailed information.
}
\seealso{
  \code{\link{gsynth}} and \code{\link{print.gsynth}}
}



