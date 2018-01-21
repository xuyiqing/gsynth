\name{preView}
\alias{preView}
\title{pre-plotting}
\description{Visualizes outcomes or treatment status before running the generalized synthetic control method.}
\usage{preView(formula, data, index, na.rm = FALSE, type = "missing", 
               ByGroup = FALSE, xlim = NULL, ylim = NULL, xlab = NULL, 
               ylab = NULL, labels = NULL, legendOff = FALSE, main = NULL, 
               id = NULL, axis.adjust = FALSE, axis.lab = "both", 
               axis.lab.gap = c(0, 0))
}
\arguments{
  \item{formula}{an object of class "formula": a symbolic description of
  the model to be fitted.}
  \item{data}{a data frame (must be with a dichotomous treatment but balanced
  is not reqiored).}
  \item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2.}
  \item{na.rm}{a logical flag indicating whether to list-wise delete
  missing data. The algorithm will report an error if missing data
  exist.}
  \item{type}{a string that specifies the type of the plot. Must be either
    of the following two: "raw" (plotting the raw data); "missing" (plotting status of each unit at each time point).}
  \item{ByGroup}{a logic flag indicating whether the data set is a Diff-in-Diff 
    data.}
  \item{xlim}{a two-element numeric vector specifying the range of x-axis. When 
    class of time varible is string, must specify not original value but a counting number e.g. \code{xlim=c(1,30)}.}
  \item{ylim}{a two-element numeric vector specifying the range of y-axis.}
  \item{xlab}{a string indicating the label of the x-axis.}
  \item{ylab}{a string indicating the label of the y-axis.}
  \item{labels}{a vector indicating legend names for missing plot.}
  \item{legendOff}{a logical flag controlling whether to show the legend.}
  \item{main}{a string that controls the title of the plot. If not
    supplied, no title will be shown.}
  \item{id}{a unit identifier of which the predicted counterfactual or
    the difference between actual and predicted counterfactual is to be shown. 
    It can also be a vector specifying units to be plotted if \code{type=="missing"} when data magitude is large. Ignored if \code{type} is "raw".}
  \item{axis.adjust}{a logic flag indicating whether to adjust labels on x-axis. 
    Useful when class of time variable is string and data magnitude is large.}
  \item{axis.lab}{a string indicating whether labels on x-axis and y-axis will be shown. There are 4 types: "both": labels on both axes will be shown; "unit": only labels on y-axis; "time": only labels on x-axis; "none": no labels will be shown. Default is \code{axis.lab = "both"}.  }
  \item{axis.lab.gap}{a numeric vector setting gaps between labels on x or y-axis for "missing" plot. Default is \code{axis.lab.gap = c(0, 0)}, which means that all labels will be shown. Useful for data set with large N or T.}
}
\details{
  \code{preView} visualizes the raw data used by the generalized synthetic control method or the treatment status of the data set.
}
\value{
   \item{obs.missing}{a matrix stroing status of each unit at each time point. 
    For DID data, \code{1} for treat group units at pre-treatment period, \code{2} for treat group units at post-treatment period, \code{3} for control group units, and \code{4} for missing values. And for non-DID data, \code{1} for under treatment, \code{2} for not under treatment, and \code{4} for missing values. Useful for unbalanced panel data.}
  \item{plot}{a \code{ggplot2} object.}
}
\author{
  Yiqing Xu <yiqingxu@ucsd.edu>
  
  Licheng Liu <liulch.16@sem.tsinghua.edu.cn>
}
\references{Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25, 
  Iss. 1, January 2017, pp. 57-76. Available at: \url{https://doi.org/10.1017/pan.2016.2}.

  See \url{http://yiqingxu.org/software/gsynth/gsynth_examples.html} for
  more detailed information.
}
\examples{
library(gsynth)
data(gsynth)
out <- preView(Y ~ D + X1 + X2, data = simdata, 
              index = c("id","time"), type="raw") 
}




