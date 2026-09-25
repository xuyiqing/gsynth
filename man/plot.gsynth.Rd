\name{plot.gsynth}
\alias{plot.gsynth}
\title{Plotting}
\description{Visualizes estimation results of the generalized synthetic control method.}
\usage{\method{plot}{gsynth}(x, type = "gap", xlim = NULL, ylim = NULL,
            xlab = NULL, ylab = NULL, legendOff = FALSE, raw = "none",
            main = NULL, nfactors = NULL, id = NULL, axis.adjust = FALSE,
            theme.bw = TRUE, shade.post = FALSE,
            legacy.style = FALSE, highlight = NULL, highlight.fill = FALSE,
            \dots)
}
\arguments{
  \item{x}{a \code{\link{gsynth}} object.}
  \item{type}{a string that specifies the type of the plot. Must be one
    of the following: "gap" (the estimated average treatment effect on
    the treated by period); "counterfactual", or "ct" for short (the
    treated units' outcomes and their predicted Y(0)'s); "factors" (the
    estimated factors); "loadings" (the distribution of the estimated
    factor loadings); "loading.overlap" (a diagnostic scatter of treated-
    and control-unit loadings in the first two factor dimensions, with the
    convex hull of control loadings shaded; new in v1.5.0, useful for GSC
    overlap diagnostics); "raw" (the raw outcome data); "missing" (the
    treatment status of each unit in each period, and the missing data).
    "raw" and "missing" are drawn by \pkg{panelView} from the data stored
    in the fit; the other types are drawn by \pkg{fect}.}
  \item{xlim}{a two-element numeric vector specifying the range of x-axis. When
    class of time variable is string, must specify not original value but a
    counting number e.g. \code{xlim=c(1,30)}.}
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
  \item{main}{a string for the title. \code{NULL} (default) uses the plot
    type's default title where it has one (for example, "Raw Data" or
    "Treatment Status and Missing Data"); \code{main = ""} removes it.}
  \item{nfactors}{a positive integer: the number of factors whose loadings
    are shown when \code{type = "loadings"} (plus the unit fixed effects,
    when the model has them). The default, \code{NULL}, is
    \code{min(r.cv, 4)}; the same value is used when \code{nfactors} is
    larger than \code{r.cv}. Ignored for other types: the "factors" plot
    always shows all estimated factors.}
  \item{id}{unit id(s). For "raw" and "missing", the units to show. For
    "counterfactual", the treated unit(s) whose outcome and predicted Y(0)
    are shown (averaged over the units when several are given); an id that
    is not a treated unit stops with an error. Ignored for the other
    types: the "gap" plot always shows the average effect of all treated
    units.}
  \item{axis.adjust}{a logical flag indicating whether to adjust labels on x-axis.
    Useful when class of time variable is string and data magnitude is large.}
  \item{theme.bw}{a logical flag indicating whether to use a black/white theme.}
  \item{shade.post}{a logical flag controlling whether to shade the post-treatment periods.}
  \item{legacy.style}{a logical flag (default \code{FALSE}). When \code{TRUE},
    reproduces the pre-v1.5.0 visual recipe (legacy axis text sizes, colors,
    and plot elements). Useful for byte-identical rendering of figures
    embedded in papers under review. New in v1.5.0.}
  \item{highlight}{controls which test-period types receive colored-glyph
    highlighting. \code{NULL} (default) auto-detects from the fit;
    \code{TRUE}/\code{FALSE} enables/disables across all types; or a character
    subset of \code{c("placebo", "carryover", "carryover.rm")} to highlight
    selected types only. New in v1.5.0.}
  \item{highlight.fill}{a logical flag (default \code{FALSE}). When \code{TRUE},
    draws a lightened-tone background rectangle behind each highlighted
    period (slide / talk style). Default \code{FALSE} keeps the glyph-only
    minimal look. New in v1.5.0.}
  \item{\dots}{other arguments passed to \code{fect::plot.fect()} (not used
    for "raw" and "missing").}
}
\details{
  \code{plot.gsynth} visualizes the raw data used by, or estimation
  results obtained from, the generalized synthetic control method.
}
\value{
  A plot object: a \pkg{ggplot2} plot; for \code{type = "loadings"}, a
  \pkg{GGally} \code{ggmatrix}; for "raw" and "missing", the \pkg{ggplot2}
  plot built by \pkg{panelView}. Nothing is drawn until the object is
  printed; at the console it is printed automatically. Change it with
  \code{+}, for example \code{plot(out) + ggplot2::theme_bw()}.
}
\author{
  Yiqing Xu <yiqingxu@stanford.edu>, Stanford University
}
\references{Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25,
  Iss. 1, January 2017, pp. 57-76.

  See \url{https://yiqingxu.org/packages/gsynth/} for
  more detailed information.
}
\seealso{
  \code{\link{gsynth}} and \code{\link{print.gsynth}}
}
