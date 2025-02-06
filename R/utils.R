#' Evaluates whether an object is a model
#'
#' @description Checks whether an object conforms to the expected model structure.
#' @export
#' @param obj An object to test.
#' @return TRUE if the object is a model, FALSE otherwise.
isModel <- function(obj) {
  # Check if the object is NULL
  if (is.null(obj)) {
    return(FALSE)
  }
  
  # Check if the object is a list
  if (!is.list(obj)) {
    return(FALSE)
  }
  
  # Check if the object contains the expected attributes
  required_fields <- c("indexes", "features", "coeff", "auc", "hash")
  has_fields <- all(sapply(required_fields, function(field) !is.null(obj[[field]])))
  
  return(has_fields)
}


#' Evaluates whether an object is a population of models
#'
#' @description Checks whether an object is a list of valid model objects.
#' @param obj An object to test.
#' @return TRUE if the object is a population of models, FALSE otherwise.
#' @export
isPopulation <- function(obj) {
  # Check if the object exists and is a non-empty list
  if (is.null(obj) || !is.list(obj) || length(obj) == 0) {
    return(FALSE)
  }
  
  # Ensure all elements in the list are models
  return(all(sapply(obj, isModel)))
}


#' Evaluates whether an object is a model collection object
#'
#' @description Checks whether an object is a collection of model populations.
#' @param obj An object to test.
#' @return TRUE if the object is a model collection, FALSE otherwise.
#' @export
isModelCollection <- function(obj) {
  # Check if the object exists, is a list, and is non-empty
  if (is.null(obj) || !is.list(obj) || length(obj) == 0) {
    return(FALSE)
  }
  
  # Ensure all elements in the list are valid populations
  return(all(sapply(obj, isPopulation)))
}



#' Retrieve the Best Individual from a Population
#'
#' @description Sorts a population based on a given evaluation criterion and returns the best individual.
#' @param pop A list representing a population of models.
#' @param evalToFit A character string indicating the attribute used for sorting (default: "fit").
#' @return The best individual from the population or NULL if the population is empty or all values are NA.
#' @export
getTheBestIndividual <- function(pop, evalToFit = "fit") {
  # Validate input: check if the population is non-empty
  if (is.null(pop) || length(pop) == 0) {
    return(NULL)
  }
  
  # Sort the population by the specified attribute
  sorted_pop <- tryCatch(
    sortPopulation(pop, evalToOrder = evalToFit),
    error = function(e) {
      message("getTheBestIndividual: Unable to sort population. Returning NULL.")
      return(NULL)
    }
  )
  
  # Return the best individual (first in sorted population)
  if (!is.null(sorted_pop) && length(sorted_pop) > 0) {
    return(sorted_pop[[1]])
  } else {
    return(NULL)
  }
}


#' Get the best model from a classifier result
#'
#' @description Extracts a given attribute from a population of predomics objects.
#' @param element2get A character string specifying the attribute to extract.
#' @param toVec Logical; should the results be unlisted into a vector? (default: TRUE)
#' @param na.rm Logical; should NA values be removed? (default: TRUE)
#' @return A vector or list of extracted attributes, depending on `toVec`.
#' @export
populationGet_X <- function(element2get, toVec = TRUE, na.rm = TRUE) {
  
  # Custom function to extract the attribute from a population
  func <- function(pop) {
    
    # Return NA if the population is empty
    if (length(pop) == 0) {
      return(NA)
    }
    
    # Extract the desired attribute from each individual in the population
    res <- lapply(pop, function(indiv) {
      if (!is.list(indiv) || is.null(indiv[[element2get]])) {
        return(NA)
      } else {
        return(indiv[[element2get]])
      }
    })
    
    # Convert to a vector if requested
    if (toVec) {
      res <- unlist(res, use.names = FALSE)
    }
    
    # Remove NA values if requested
    if (na.rm) {
      res <- Filter(Negate(is.na), res)
    }
    
    return(res)
  }
  
  return(func)
}


#' Sort a Population Based on an Evaluation Criterion
#'
#' @description Sorts a population of models based on a specified attribute.
#' @param pop A list representing a population of models.
#' @param evalToOrder A character string indicating the attribute used for sorting (default: "fit").
#' @param decreasing Logical; should sorting be in decreasing order? (default: TRUE).
#' @return A sorted list of models based on the specified attribute.
#' @export
sortPopulation <- function(pop, evalToOrder = "fit", decreasing = TRUE) {
  # Validate input: check if the population is non-empty
  if (is.null(pop) || length(pop) == 0) {
    stop("sortPopulation: Population is empty or NULL.")
  }
  
  # Extract attribute values for sorting
  attr_values <- populationGet_X(element2get = evalToOrder, toVec = TRUE, na.rm = FALSE)(pop)
  
  # Handle cases where attr_values is NULL or all NA
  if (is.null(attr_values) || all(is.na(attr_values))) {
    stop(paste("sortPopulation: All values of", evalToOrder, "are NA or missing. Cannot sort."))
  }
  
  # Create an ordering index (ignoring NAs by setting them to the lowest rank)
  order_idx <- order(attr_values, decreasing = decreasing, na.last = TRUE)
  
  # Return sorted population
  return(pop[order_idx])
}




#' Convert a population of model objects into a dataframe for plotting
#'
#' @description Extracts attributes from each model in the population and creates a dataframe for further exploration.
#' @param pop A list of model objects (i.e., a population of models).
#' @param attributes A vector of attribute names to extract (default: model structure attributes).
#' @return A data frame with attributes for each model.
#' @export
populationToDataFrame <- function(
    pop, 
    attributes = c("coeff", "indexes", "k", "auc", "epoch", 
                   "fit", "specificity", "sensitivity", "accuracy", 
                   "threshold", "language", "data_type", "data_type_minimum", "hash")
) {
  
  # Validate input population
  if (!isPopulation(pop)) {
    stop("populationToDataFrame: Please provide a valid population object")
  }
  
  # Ensure requested attributes exist in models
  model_sample <- pop[[1]]
  available_attrs <- names(model_sample)
  missing_attrs <- attributes[!attributes %in% available_attrs]
  if (length(missing_attrs) > 0) {
    stop(paste("populationToDataFrame: unknown attributes:", paste(missing_attrs, collapse = ", ")))
  }
  
  # Fetch all attributes using populationGet_X
  extracted_data <- lapply(attributes, function(attr) {
    values <- populationGet_X(element2get = attr, toVec = FALSE, na.rm = FALSE)(pop)
    
    # Handle `coeff` and `indexes`: Convert to concatenated string if multiple values exist
    if (attr %in% c("coeff", "indexes")) {
      values <- lapply(values, function(x) {
        if (is.null(x)) return(NA_character_)
        return(paste(x, collapse = ", "))  # Convert to string
      })
    }
    
    # Ensure all attributes return a vector of length equal to population size
    if (length(values) != length(pop)) {
      warning(paste("populationToDataFrame: Attribute", attr, 
                    "returned", length(values), "values instead of", length(pop)))
      values <- rep(NA, length(pop))  # Pad with NAs to prevent misalignment
    }
    
    return(unlist(values))  # Convert list to vector
  })
  
  # Convert extracted list into a data frame
  df <- as.data.frame(extracted_data, stringsAsFactors = FALSE)
  colnames(df) <- attributes
  
  # Assign row names correctly
  rownames(df) <- paste("mod", seq_along(pop), sep = "_")
  
  return(df)
}


#' Validate dataset and response vector consistency
#'
#' @description Checks if X (dataset) and y (response vector) have coherent dimensions.
#' @param X A data matrix or dataframe.
#' @param y A response vector.
#' @return None; stops execution with an error if validation fails.
#' @export
check.X_y <- function(X, y) {
  
  # Check if X is a matrix or dataframe
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("check.X_y: X must be a matrix or a dataframe.")
  }
  
  # Validate that X and y have the same number of rows
  if (ncol(X) != length(y)) { 
    stop("check.X_y: The number of rows in X (", ncol(X), ") must match the length of y (", length(y), ").") 
  }
}


#' Transform the model object into a dense coefficient vector
#'
#' @description Converts a model object into a dense (long) format vector.
#' @param natts The number of attributes (features).
#' @param mod A model object.
#' @return A dense vector representing the model coefficients.
#' @export
modelToDenseVec <- function(natts, mod) {
  
  # Validate model
  if (!isModel(mod)) {
    stop("modelToDenseVec: The model object is not valid.")
  }
  
  # Validate number of attributes
  if (!is.numeric(natts) || length(natts) != 1 || natts <= 0) {
    stop("modelToDenseVec: `natts` must be a single positive integer representing the number of features.")
  }
  
  # Initialize dense coefficient vector
  dense_vector <- rep(0, natts)
  
  # Extract model indices and coefficients safely
  indices <- populationGet_X("indexes", toVec = TRUE, na.rm = FALSE)(list(mod))
  coefficients <- populationGet_X("coeff", toVec = TRUE, na.rm = FALSE)(list(mod))
  
  # Ensure indices and coefficients are correctly mapped
  if (length(indices) != length(coefficients)) {
    stop("modelToDenseVec: Mismatch between model indices and coefficients.")
  }
  
  # Assign coefficients to the correct positions in the dense vector
  dense_vector[indices] <- coefficients
  
  return(dense_vector)
}


#' Builds a list of dense vector coefficients from a list of models
#'
#' @param X A dataset (matrix or dataframe).
#' @param y A vector of labels.
#' @param list.models A list of model objects.
#' @return A list of dense vectors of coefficients.
#' @export
listOfModelsToListOfDenseVec <- function(X, y, list.models) {
  
  # Validate input population
  if (!isPopulation(list.models)) {
    stop("listOfModelsToListOfDenseVec: Please specify a valid population of model objects")
  }
  
  # Validate X and y
  check.X_y(X, y)
  
  # Convert models to dense vectors
  res <- lapply(list.models, function(model) {
    modelToDenseVec(nrow(X), model)
  })
  
  return(res)
}




#' Convert a list of models to a dense coefficient matrix
#'
#' @description Converts a list of models into a dense coefficient matrix, ensuring correct ordering and filtering.
#' @param X The dataset (matrix or dataframe).
#' @param y The class vector.
#' @param list.models A list of model objects.
#' @param rm.empty Logical; should models with all-zero coefficients be removed? (default: TRUE).
#' @param order.row Logical; should rows be ordered by occurrence? (default: TRUE).
#' @return A dense coefficient matrix with features as rows and models as columns.
#' @export
listOfModelsToDenseCoefMatrix <- function(X, y, list.models, rm.empty = TRUE, order.row = TRUE) {
  
  # Validate input population
  if (!isPopulation(list.models)) {
    stop("listOfModelsToDenseCoefMatrix: please provide a valid population of models")
  }
  
  # Validate X and y
  check.X_y(X, y)
  
  # Convert models to dense vectors
  pop.dense <- listOfModelsToListOfDenseVec(X = X, y = y, list.models = list.models)
  
  # Convert list to matrix
  pop.dense <- as.matrix(do.call(cbind, pop.dense))
  
  # Assign row and column names
  rownames(pop.dense) <- rownames(X)
  colnames(pop.dense) <- paste("model", seq_len(ncol(pop.dense)), sep = "_")
  
  # Handle NA row names
  if (any(is.na(rownames(pop.dense)))) {
    warning("listOfModelsToDenseCoefMatrix: Some features are NA in pop.dense ... omitting them")
    pop.dense <- pop.dense[!is.na(rownames(pop.dense)), ]
  }
  
  # Order rows by feature occurrence
  if (order.row) {
    pop.dense <- pop.dense[order(rowSums(pop.dense != 0, na.rm = TRUE), decreasing = TRUE), ]
    if (is.numeric(pop.dense)) { # Handle single model case
      pop.dense <- as.matrix(as.data.frame(pop.dense))
    }
  }
  
  # Remove empty rows (features with all zeros)
  if (rm.empty) {
    features.ord <- rownames(pop.dense)
    nonzero_rows <- rowSums(pop.dense != 0, na.rm = TRUE) != 0
    pop.dense <- pop.dense[nonzero_rows, ]
    
    if (is.vector(pop.dense)) { # Handle single model case
      pop.dense <- as.matrix(as.data.frame(pop.dense))
      rownames(pop.dense) <- features.ord[nonzero_rows]
      colnames(pop.dense) <- colnames(X)
    }
  }
  
  return(pop.dense)
}



#' Compute the Binomial Confidence Interval
#'
#' @description Computes a confidence interval for a binomial proportion using the normal approximation.
#' @param accuracy The observed accuracy or proportion.
#' @param n The total sample size.
#' @param p The p-value threshold (default: 0.05).
#' @return A named vector with lower bound ("inf"), accuracy, and upper bound ("sup").
#' @export
confInterBinomial <- function(accuracy, n, p = 0.05) {
  if (n <= 0) {
    stop("confInterBinomial: Sample size (n) must be greater than zero.")
  }
  
  # Compute the critical value for the normal approximation
  t_value <- -qnorm(p / 2)
  
  # Compute standard error
  std_error <- sqrt((accuracy * (1 - accuracy)) / n)
  
  # Compute confidence interval
  ci_range <- t_value * std_error
  lower_bound <- max(0, accuracy - ci_range)  # Ensure it doesn't go below 0
  upper_bound <- min(1, accuracy + ci_range)  # Ensure it doesn't go above 1
  
  res <- c(lower_bound, accuracy, upper_bound)
  names(res) <- c("inf", "accuracy", "sup")
  
  return(res)
}



#' Select the Top Significant Best Models from a Population
#' 
#' @description Selects the best part of a population where models are significantly **not different** from the best model.
#' @param pop A list of model objects (i.e., a population of models).
#' @param score The attribute of the model used for evaluation (default: "fit").
#' @param p The p-value threshold for statistical significance (default: 0.05).
#' @param k_penalty A penalty applied to the score based on the k-sparsity value (default: 0).
#' @param k_max A threshold for sparsity selection (default: 0, meaning no sparsity filtering).
#' @return A filtered sub-population of models or `NULL` if no models meet the criteria.
#' @export
selectBestPopulation <- function(pop, score = "fit", p = 0.05, k_penalty = 0, k_max = 0) {
  # Validate that pop is a population
  if (!isPopulation(pop)) {
    stop("selectBestPopulation: Please provide a valid population.")
  }
  
  # Sort population by score
  pop <- sortPopulation(pop, evalToOrder = score)
  
  # Extract scores and sparsity levels
  eval <- populationGet_X(element2get = score, toVec = TRUE, na.rm = FALSE)(pop)
  spar <- populationGet_X(element2get = "k", toVec = TRUE, na.rm = FALSE)(pop)
  
  # Apply penalty to evaluation score
  eval <- eval - (k_penalty * spar)
  
  # Compute confidence interval for the best model
  eval_ci <- confInterBinomial(accuracy = eval[1], n = length(pop), p = p)
  eval_threshold <- eval_ci["inf"]  # Lower bound of confidence interval
  
  # Select models above the significance threshold
  if (k_max == 0) {
    selected_indices <- eval > eval_threshold
  } else {
    selected_indices <- (eval > eval_threshold) & (spar <= k_max)
  }
  
  # Handle missing values
  selected_indices[is.na(selected_indices)] <- FALSE
  
  # Ensure at least one model is selected
  selected_pop <- pop[selected_indices]
  
  if (!isPopulation(selected_pop)) {
    warning("selectBestPopulation: No models were found after selection.")
    return(NULL)
  } else {
    return(selected_pop)
  }
}


