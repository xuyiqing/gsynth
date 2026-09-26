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
    of the following: "gap" (plotting the average treatment effect on
    the treated); "raw" (plotting the raw data); "counterfactual", or "ct" for short,
    (plotting predicted Y(0)'s); "factors" (plotting estimated factors);
    "loadings" (plotting the distribution of estimated factor
    loadings); "loading.overlap" (a diagnostic scatter of treated- and
    control-unit loadings in the first two factor dimensions, with the
    convex hull of control loadings shaded; new in v1.5.0, useful for
    GSC overlap diagnostics); "missing" (plotting status of each unit at
    each time point). "raw" and "missing" are drawn by \pkg{panelView}
    from the data stored in the fit.}
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
    type's default title where it has one (for example "Raw Data" or
    "Treatment Status and Missing Data"); \code{""} removes it.}
  \item{nfactors}{a positive integer: the number of factors whose
    loadings are shown in the "loadings" plot. Default
    \code{min(r.cv, 4)}. The "factors" plot shows all factors.}
  \item{id}{unit id(s). For "raw" and "missing", the units to show. For
    "counterfactual", the treated unit(s) to show: with one unit, its
    outcome and counterfactual; with several, their averages. Ignored by
    the other types; the "gap" plot always shows the average effect over
    all treated units.}
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
  \item{\dots}{other argv.}
}
\details{
  \code{plot.gsynth} visualizes the raw data used by, or estimation
  results obtained from, the generalized synthetic control method.
}
\value{
  A \pkg{ggplot2} object (for \code{type = "loadings"}, a \pkg{GGally}
  \code{ggmatrix}; for "raw" and "missing", the plot built by
  \pkg{panelView}). Nothing is drawn until the object is printed; at the
  console it is printed automatically. Add layers or themes with
  \code{+}.
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



