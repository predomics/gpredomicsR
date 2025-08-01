% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visu.R
\name{plotBarcode}
\alias{plotBarcode}
\title{Plot Barcode-style Heatmap of Feature Abundance}
\usage{
plotBarcode(
  X = NULL,
  y = NULL,
  main = "",
  ylabl = "",
  ylabr = "",
  fixed.scale = TRUE,
  data = NULL,
  select_features = NULL,
  select_samples = NULL
)
}
\arguments{
\item{X}{A numeric matrix or data frame (features x samples).}

\item{y}{A vector of sample class labels corresponding to columns in \code{X}.}

\item{main}{A character string for the plot title.}

\item{ylabl}{(Unused; reserved for future use).}

\item{ylabr}{(Unused; reserved for future use).}

\item{fixed.scale}{Logical; if \code{TRUE}, uses a fixed log-scale color mapping.
If \code{FALSE}, the color scale is computed dynamically from data range. Default is \code{TRUE}.}

\item{data}{Optional. A list containing \code{X}, \code{y}, and optionally \code{classes}.
If provided, overrides \code{X} and \code{y}.}

\item{select_features}{Optional vector of feature names to display (subset of rownames of \code{X}).}

\item{select_samples}{Optional vector of sample names to display (subset of colnames of \code{X}).}
}
\value{
A \code{ggplot2} object containing the barcode heatmap.
}
\description{
This function generates a barcode-style heatmap from a feature-by-sample matrix,
grouped by class labels. It is designed for visualizing abundance or intensity values,
and supports both fixed and dynamic log-scaled color mapping.
}
