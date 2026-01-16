#' @section Data S3 methods:
#' These methods provide functionality for printing and extracting fields from Data objects.

#' Prints a Data object
#' @description This method prints a summary of a Data object.
#' @param data An object of class \code{Data} (Rust-backed external pointer).
#' @export
print.Data <- function(data) {
  cat(data$print(), "\n")
  invisible(data)
}

#' Extracts fields from a Data object
#' 
#' @description This method allows for the extraction of internal Rust-managed 
#' fields from a Data object into an R-compatible format.
#'
#' @param data An object of class \code{Data} (Rust-backed external pointer).
#' @param field A character string specifying the field name to extract.
#' 
#' @return A vector or list containing the requested field values.
#' @export
`[[.Data` <- function(data, field) {
  dataList <- data$get()
  if (field %in% names(dataList)) {
    return(dataList[[field]])
  }
  stop("Unknown parameter: ", field)
}

#' Gets the names of fields in a Data object
#'
#' @description This method returns the names of the fields in a Data object.
#'
#' @param data An object of class \code{Data} (Rust-backed external pointer).
#'
#' @return A character vector containing the names of the fields in the Data object.
#' @export
names.Data <- function(data) {
  return(names(data$get()))
}

#' Gets the dimensions of a Data object
#'
#' @description This method returns the dimensions of a Data object.
#'
#' @param data An object of class \code{Data} (Rust-backed external pointer).
#'
#' @return A numeric vector containing the dimensions of the Data object.
#' @export
dim.Data <- function(data) {
  res <- data$get()
  return(c(res$sample_len, res$feature_len))
}

#' @section Param S3 methods:
#' These methods provide functionality for printing, extracting and editing fields from Param objects.

#' Prints a Param object
#'
#' @description This method prints the contents of a Param object.
#'
#' @param param An object of class \code{Param} (Rust-backed external pointer).
#'
#' @export
print.Param <- function(param, ...) {
  cat(param$address(), "\n")
  print(param$get())
  invisible(param)
}

#' Sets a parameter in a Param object
#'
#' @description This method allows for setting parameters in a Param object
#' using the double bracket assignment operator.
#'
#' @param param An object of class \code{Param} (Rust-backed external pointer).
#' @param parameter A character string specifying the parameter name to set.
#' @param value The value to assign to the specified parameter.
#'
#' @return The updated Param object.
#' @export
`[[<-.Param` <- function(param, parameter, value) {
  if (is.character(value)) {
    param$set_string(parameter, value)
  } else if (is.numeric(value)) {
    param$set(parameter, value)
  } else if (is.logical(value)) {
    param$set_bool(parameter, value)
  }
  return(param)
}

#' Extracts a parameter from a Param object
#'
#' @description This method allows for the extraction of parameters from a Param object
#' using the double bracket operator.
#'
#' @param param An object of class \code{Param} (Rust-backed external pointer).
#' @param parameter A character string specifying the parameter name to extract.
#'
#' @return The value of the specified parameter.
#' @export
`[[.Param` <- function(param, parameter) {
    paramlist <- param$get()
    for (name in names(paramlist)) {
      if (parameter %in% names(paramlist[[name]])) {
        return(paramlist[[name]][[parameter]])
      }
    }
   stop("Unknown parameter: ", parameter)
}

#' @section Jury S3 Experiment methods:
#' These methods provide functionality for printing and extracting fields from Experiment objects.

#' Prints an Experiment object
#'
#' @description This method prints a summary of an Experiment object.
#'
#' @param exp An object of class \code{Experiment} (Rust-backed external pointer).
#'
#' @export
print.Experiment <- function(exp) {
  cat(exp$print(), "\n")
  invisible(exp)
}

#' Gets the number of generations in an Experiment
#'
#' @description This method returns the number of generations in an Experiment object.
#'
#' @param exp An object of class \code{Experiment} (Rust-backed external pointer).
#'
#' @return An integer representing the number of generations in the Experiment.
#' @export
length.Experiment <- function(exp) {
  return(exp$generation_number())
}

#' Extracts Generations from an Experiment object
#'
#' @description This method allows for indexing into an Experiment object to
#' retrieve specific Generations based on generation number and fold.
#'
#' @param exp An object of class \code{Experiment} (Rust-backed external pointer).
#' @param generation A numeric value specifying the generation number to extract.
#' @param fold A numeric value specifying the fold number (for cross-validation).
#' @param on_train A logical value indicating whether to access the population fitted on training set (default is
#' TRUE).
#'
#' @return A Generation object corresponding to the specified generation and fold.
#' @export
`[.Experiment` <- function(exp, generation, fold, on_train = TRUE) {
  if (!is.numeric(generation) || (length(generation) != 1)) {
  stop("Generation indexing only supports single numeric values.")
  }

  if (exp$get_n_folds() == 1) {
    return(exp$get_generation(generation))
  } else if (fold <= exp$get_n_folds()) {
    return(exp$get_fold_generation(fold, generation, on_train))
  } else {
    stop("Fold number exceeds the number of folds in this Experiment.")
  }
}

#' Plots Training History from an Experiment Object
#'
#' @description This method generates a plot of the training history from an Experiment object.
#'
#' @param x An object of class \code{Experiment} (Rust-backed external pointer).
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot object representing the training history.
#' @export
plot.Experiment <- function(x, ...) {
  plotHistory(x, x$get_data(TRUE), ...)
}

#' @section Individual S3 methods:
#' These methods provide functionality for printing and extracting fields from Individual objects.

#' Prints an Individual object
#' 
#' @description This method prints a summary of an Individual object.
#'
#' @param individual An object of class \code{Individual} (Rust-backed external pointer).
#'
#' @export
print.Individual <- function(individual) {
  individual$print()
  invisible(individual)
}

#' Extracts coefficients from an Individual object
#'
#' @description This method retrieves the coefficients of an Individual object.
#'
#' @param object An object of class \code{Individual} (Rust-backed external pointer).
#'
#' @return A named vector containing the coefficients of the Individual.
#' @export
coef.Individual <- function(object, ...) {
  res <- object$get()
  setNames(res$coefficients, res$features)
}

#' Extracts fields from an Individual object
#'
#' @description This method allows for the extraction of internal Rust-managed 
#' fields from an Individual object into an R-compatible format.
#' @param individual An object of class \code{Individual} (Rust-backed external pointer).
#' @param field A character string specifying the field name to extract.
#'
#' @return A vector or list containing the requested field values.
#' @export
`[[.Individual` <- function(individual, field) {
  individualList <- individual$get()
  if (field %in% names(individualList)) {
    return(individualList[[field]])
  }
  stop("Unknown field: ", field)
}

#' Predicts classes or scores for an Individual object
#' @param object An object of class \code{Individual}.
#' @param new_data A Data object for prediction.
#' @param type Character: "class" (default) or "score" (raw decision values).
#' @export
predict.Individual <- function(object, new_data, type = c("class", "score"), ...) {
  type <- match.arg(type)
  results <- object$predict(new_data)
  
  if (type == "score") {
    return(results$score)
  } else {
    return(results$class)
  }
}

#' @section Population S3 methods:
#' These methods provide functionality for printing and extracting fields from Population objects.

#' Prints a Population object
#'
#' @description This method prints a summary of a Population object.
#'
#' @param pop An object of class \code{Population} (Rust-backed external pointer).
#' 
#' @export
print.Population <- function(pop) {
  cat(pop$print(), "\n")
  invisible(pop)
}

#' Gets the Population size
#'
#' @description This method returns the size of a Population object.
#'
#' @param pop An object of class \code{Population} (Rust-backed external pointer).
#'
#' @return An integer representing the size of the Population.
#' @export
length.Population <- function(pop) {
  return(pop$get_size())
}

#' Extracts Individuals from a Population object
#'
#' @description This method allows for indexing into a Population object to
#' retrieve specific Individuals based on numeric indices or logical masks.
#'
#' @param pop An object of class \code{Population} (Rust-backed external pointer).
#' @param accessions A numeric vector of indices or a logical vector as a mask.
#'
#' @return A list of Individual objects corresponding to the specified indices or mask.
#' @export
`[.Population` <- function(pop, accessions) {
  if (is.numeric(accessions)) {
    pop$get_individual(accessions)
  } else if (is.logical(accessions) && length(accessions) == pop$get_size()) {
    pop$filter_by_mask(as.integer(accessions))
  } else {
    warning("Population indexing only supports numeric or logical vectors.")
  }
}

#' Extracts fields from a Population object
#'
#' @description This method allows for the extraction of specific fields from all
#' Individuals within a Population object.
#'
#' @param pop An object of class \code{Population} (Rust-backed external pointer).
#' @param field A character string specifying the field name to extract from each Individual.
#'
#' @return A vector containing the requested field values from all Individuals in the Population.
#' @export
`[[.Population` <- function(pop, field) {
  for (individual in pop$get()$individuals) {
    if (field %in% names(individual)) {
      return(sapply(pop$get()$individuals, function(ind) ind[[field]]))
    }
  }
}

#' @section Jury S3 methods:
#' These methods provide functionality for printing and extracting fields from Jury objects.

#' Prints a Jury object
#'
#' @description This method prints a summary of a Jury object.
#'
#' @param x An object of class \code{Jury} (Rust-backed external pointer).
#'
#' @export
print.Jury <- function(jury, ...) {
  cat(jury$print(), "\n")
  invisible(jury)
}

#' Extracts parameters from a Jury object
#'
#' @description This method allows for indexing into a Jury object to
#' retrieve specific Individuals based on numeric indices or logical masks.
#'
#' @param jury An object of class \code{Jury} (Rust-backed external pointer).
#' @param accessions A numeric vector of indices or a logical vector as a mask.
#'
#' @return A list of Individual objects corresponding to the specified indices or mask.
#' @export
`[.Jury` <- function(jury, accessions) {
  if (is.numeric(accessions)) {
    jury$get_population()$get_individual(accessions)
  } else if (is.logical(accessions) && length(accessions) == jury$get_population()$get_size()) {
    jury$get_population()$filter_by_mask(as.integer(accessions))
  } else {
    warning("Population indexing only supports numeric or logical vectors.")
  }
}

#' Predicts classes or scores for an Jury object
#' @param object An object of class \code{Jury}.
#' @param new_data A Data object for prediction.
#' @param type Character: "class" (default) or "score" (raw decision values).
#' @export
predict.Jury <- function(object, new_data, type = c("class", "score"), ...) {
  type <- match.arg(type)
  
  results <- tryCatch(
    object$predict(new_data),
    error = function(e) {
      stop("Failed to make predictions with Jury: ", conditionMessage(e), call. = FALSE)
    }
  )
  
  if (type == "score") {
    return(results$score)
  } else {
    return(results$class)
  }
}

#' Gets the Jury population size
#'
#' @description This method returns the size of the Jury population.
#'
#' @param jury An object of class \code{Jury} (Rust-backed external pointer).
#'
#' @return An integer representing the size of the Jury population.
#' @export
length.Jury <- function(jury) {
  return(jury$get_population()$get_size())
}

