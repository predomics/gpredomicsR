% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{populationToDataFrame}
\alias{populationToDataFrame}
\title{Convert a population of model objects into a dataframe for plotting}
\usage{
populationToDataFrame(
  pop,
  attributes = c("coeff", "indexes", "k", "auc", "epoch", "fit", "specificity",
    "sensitivity", "accuracy", "threshold", "language", "data_type", "data_type_minimum",
    "hash")
)
}
\arguments{
\item{pop}{A list of model objects (i.e., a population of models).}

\item{attributes}{A vector of attribute names to extract (default: model structure attributes).}
}
\value{
A data frame with attributes for each model.
}
\description{
Extracts attributes from each model in the population and creates a dataframe for further exploration.
}
