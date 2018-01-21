\name{preView}
\alias{preView}
\title{Plotting data structure and raw outcome data before estimating a panel model}
\description{Visualizes outcomes or treatment status before running the generalized synthetic control method.}
\usage{preView(formula, data, index, na.rm = FALSE, type = "missing", 
               by.group = FALSE, xlim = NULL, ylim = NULL, xlab = NULL, 
               ylab = NULL, legendOff = FALSE, legend.labs = NULL, main = NULL, 
               id = NULL, axis.adjust = FALSE, axis.lab = "both", 
               axis.lab.gap = c(0, 0))
}
\arguments{
  \item{formula}{an object of class "formula": a symbolic description of
  the model to be fitted. The first variable on the right-hand-side is designated as the treatment indicator.}
  \item{data}{a data frame (must be with a dichotomous treatment);  the panel does not have to be balanced.}
  \item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2.}
  \item{na.rm}{a logical flag indicating whether to list-wise delete missing data. The algorithm will report an error if missing data exist.}
  \item{type}{a string that specifies the type of the plot. Must be either \code{"missing"} (default), which plots the treatment/missing-data status of each unit at each time point) or \code{"raw"}, which plots the raw outcome data}
  \item{by.group}{a logic flag indicating whether the data should be plotted in separate groups based on treatment status changes; only applies when \code{type=="raw"}.}
  \item{xlim}{a two-element numeric vector specifying the range of x-axis. When the class of time variable is string, must specify the range of strings to be shown, e.g. \code{xlim=c(1,30)}.}
  \item{ylim}{a two-element numeric vector specifying the range of y-axis.}
  \item{xlab}{a string indicating the label of the x-axis.}
  \item{ylab}{a string indicating the label of the y-axis.}
  \item{legendOff}{a logical flag controlling whether to show the legend.}
  \item{legend.labs}{a vector specifying the legend labels. Ignored when \code{legendOff=TRUE}.}
  \item{main}{a string that controls the title of the plot.}
  \item{id}{a vector specifying units to be shown in the \code{"missing"} plot. Useful when the number of units is very large.}
  \item{axis.adjust}{a logic flag indicating whether to adjust labels on the x-axis.  Useful when the class of time variable is string and there are many time periods.}
  \item{axis.lab}{a string indicating whether labels on the x- and y-axis will be shown. There are four options: \code{"both"} (default): labels on both axes will be shown; \code{"unit"}: only labels on y-axis will be shown; \code{"time"}: only labels on the x-axis will be shown; "none": no labels will be shown.}
  \item{axis.lab.gap}{a numeric vector setting the gaps between labels on the x- or y-axis for \code{"missing"} plot. Default is \code{axis.lab.gap = c(0, 0)}, which means that all labels will be shown. Useful for datasets with large N or T.}
}
\details{
  \code{preView} visualizes treatment status, missing values, and the raw outcome data of a  time-series cross-sectional dataset.
}
% \value{
%    \item{obs.missing}{a matrix showing status of each unit at each time point. 
%     For DID data, \code{1} for treat group units at pre-treatment period, \code{2} for treat group units at post-treatment period, \code{3} for control group units, and \code{4} for missing values. And for non-DID data, \code{1} for under treatment, \code{2} for not under treatment, and \code{4} for missing values. Useful for unbalanced panel data.}
%   \item{plot}{a \code{ggplot2} object.}
% }
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
preView(Y ~ D + X1 + X2 , data = simdata, index = c("id","time")) 
preView(Y ~ D + X1 + X2, data = simdata, index = c("id","time"), type = "raw") 
}




