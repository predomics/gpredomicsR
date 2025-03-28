% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visu.R
\name{plotAbundance}
\alias{plotAbundance}
\title{Plot Feature Abundance by Class}
\usage{
plotAbundance(
  features,
  X,
  y,
  topdown = TRUE,
  main = "",
  plot = TRUE,
  log_scale = FALSE,
  col.pt = c("deepskyblue4", "firebrick4"),
  col.bg = c("deepskyblue1", "firebrick1")
)
}
\arguments{
\item{features}{Character vector of feature names to plot.}

\item{X}{Data matrix or data frame (features in rows, samples in columns).}

\item{y}{Vector of class labels (factor for classification) or numeric values (for regression).}

\item{topdown}{Logical; if \code{TRUE}, features are displayed from top to bottom (default: \code{TRUE}).}

\item{main}{Plot title (default: \code{""}).}

\item{plot}{Logical; if \code{TRUE}, displays the plot, otherwise returns statistical results (default: \code{TRUE}).}

\item{log_scale}{Logical; if \code{TRUE}, plots abundance data on a log10 scale (default: \code{FALSE}).}

\item{col.pt}{Colors for data points (default: \code{c("deepskyblue4", "firebrick4")}).}

\item{col.bg}{Colors for boxplot backgrounds (default: \code{c("deepskyblue1", "firebrick1")}).}
}
\value{
If \code{plot = TRUE}, returns a ggplot object. Otherwise, returns a data frame of statistical test results.
}
\description{
Visualizes the abundance of selected features across different classes using \strong{boxplots}.
Supports both \strong{classification} (discrete classes) and \strong{regression} (continuous response).
}
\examples{
# Example: Classification
features <- c("feature1", "feature2")
X <- data.frame(feature1 = rnorm(100), feature2 = rnorm(100))
y <- sample(c(1, -1), 100, replace = TRUE)
plotAbundanceByClass(features, X, y, log_scale = TRUE)

}
