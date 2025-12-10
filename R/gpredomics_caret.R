#' Caret Integration for gpredomicsR
#'
#' @description
#' This module provides integration with the caret package, enabling
#' grid search and random search hyperparameter tuning for gpredomics models.
#'
#' @details
#' The integration allows users to leverage caret's train() function with
#' gpredomics algorithms (GA, BEAM, MCMC). Both grid search and random search
#' are supported for hyperparameter optimization.
#' 
#' There are two ways to use gpredomics with caret:
#' 
#' 1. Using train_gpredomics_caret() - Convenience wrapper that handles everything
#' 2. Using caret::train() with method = getModelInfo_gpredomics()
#'
#' @section Available Parameters:
#' The default tuning grid includes only commonly tuned parameters (population_size,
#' max_epochs, ga_kmin, ga_kmax, beam_kmin, beam_kmax, language, data_type). However, 
#' you can add **any parameter** supported by gpredomics to your custom grid!
#' 
#' Use `gpredomics_params_info()` to discover all 50+ available parameters organized by category:
#' - General: population_size, max_epochs, language, data_type, algo, fit
#' - GA: ga_kmin, ga_kmax, select_elite_pct, mutated_children_pct, forced_diversity_pct
#' - BEAM: beam_kmin, beam_kmax, best_models_criterion, max_nb_of_models, method
#' - Penalties: k_penalty, fr_penalty, bias_penalty, threshold_ci_penalty
#' - GA-specific: select_elite_pct, mutated_children_pct, forced_diversity_pct
#' - BEAM-specific: best_models_criterion, max_nb_of_models
#' - MCMC-specific: n_iter, n_burn, lambda
#' - Cross-validation: inner_folds, overfit_penalty, outer_folds
#' - Voting: fbm_ci_alpha, min_perf, min_diversity
#' - Data processing: feature_minimal_prevalence_pct, feature_selection_method
#' 
#' Example:
#' \code{
#' # Discover all parameters
#' all_params <- gpredomics_params_info()
#' 
#' # Create custom grid with any parameters
#' custom_grid <- expand.grid(
#'   population_size = 5000,
#'   max_epochs = 100,
#'   ga_kmin = 1,
#'   ga_kmax = 200,
#'   beam_kmin = 1,
#'   beam_kmax = 200,
#'   language = "bin",
#'   data_type = "raw",
#'   algo = "ga",                      # Choose algorithm (ga, beam, mcmc)
#'   fit = c("auc", "mcc"),            # Add fitness metric
#'   k_penalty = c(0, 0.0001),         # Add complexity penalty
#'   select_elite_pct = c(2, 5),       # Add GA parameter
#'   feature_selection_method = "wilcoxon"  # Add data parameter
#' )
#' }
#'
#' @section Rust Object Access Patterns:
#' This module correctly interfaces with Rust objects from lib.rs:
#' 
#' - Experiment: Use methods get_best_population(), get_population(generation)
#' - Population: Use method get_individual(index) with 1-based R indexing
#' - Population: Use method get_first_pct(pct) to get top percentage subset
#' - Individual: Use methods get(), get_metrics(), predict(data)
#' - Individual$predict() returns list with $class and $score components
#' - All internal attributes (threshold, auc, etc.) accessed via methods, not directly
#'
#' @section Prediction Strategy:
#' During cross-validation and final predictions, gpredomics uses the **best individual**
#' from the final population. This ensures:
#' 
#' - Consistent evaluation across all CV folds
#' - Direct comparison with other caret models
#' - Optimal performance based on the chosen fitness metric
#' 
#' The best individual is selected based on the `fit` parameter specified in your
#' gpredomics configuration (e.g., AUC, MCC, sensitivity, etc.).
#'
#' @importFrom pROC roc auc
#' @author gpredomicsR Team
#' @name caret-integration


#' Gpredomics Summary Function for Caret
#'
#' @description
#' Custom summary function that computes Gpredomics classification metrics:
#' Sensitivity, Specificity, Accuracy, MCC, NPV, PPV, F1-score, G-mean.
#' 
#' Use this with trainControl(summaryFunction = gpredomicsSummary) to get
#' comprehensive metrics during hyperparameter tuning.
#'
#' @param data Data.frame with observed classes (obs) and predicted classes (pred)
#' @param lev Character vector of class levels
#' @param model Model type (currently unused)
#'
#' @return Named vector of performance metrics
#' @export
#'
#' @examples
#' \dontrun{
#' ctrl <- trainControl(
#'   method = "cv",
#'   number = 5,
#'   summaryFunction = gpredomicsSummary
#' )
#' 
#' model <- train_gpredomics_caret(
#'   outcome ~ .,
#'   data = mydata,
#'   trControl = ctrl,
#'   metric = "MCC"  # Can optimize on any metric
#' )
#' }
gpredomicsSummary <- function(data, lev = NULL, model = NULL) {
  
  # Ensure required columns exist
  if (!all(c("obs", "pred") %in% names(data))) {
    stop("data must contain 'obs' and 'pred' columns")
  }
  
  # Compute AUC from scores if available, otherwise from predicted probabilities
  auc_val <- NA
  if (!is.null(lev) && length(lev) == 2) {
    # Check if we have probability columns (class probabilities from prob function)
    prob_col <- lev[2]  # Probability of positive class
    if (prob_col %in% names(data)) {
      # Use pROC to compute AUC from probabilities
      tryCatch({
        roc_obj <- pROC::roc(data$obs, data[[prob_col]], levels = lev, direction = "<", quiet = TRUE)
        auc_val <- as.numeric(pROC::auc(roc_obj))
      }, error = function(e) {
        auc_val <- NA
      })
    }
  }
  
  # Convert to binary 0/1 for confusion matrix calculations
  obs_binary <- as.integer(data$obs == lev[2])
  pred_binary <- as.integer(data$pred == lev[2])
  
  # Confusion matrix components
  tp <- sum(pred_binary == 1 & obs_binary == 1)  # True Positives
  tn <- sum(pred_binary == 0 & obs_binary == 0)  # True Negatives
  fp <- sum(pred_binary == 1 & obs_binary == 0)  # False Positives
  fn <- sum(pred_binary == 0 & obs_binary == 1)  # False Negatives
  
  # Base metrics
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else 0
  accuracy <- if ((tp + tn + fp + fn) > 0) (tp + tn) / (tp + tn + fp + fn) else 0
  
  # Additional metrics (matching Gpredomics)
  # Positive Predictive Value (Precision)
  ppv <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  
  # Negative Predictive Value
  npv <- if ((tn + fn) > 0) tn / (tn + fn) else 0
  
  # F1-score (harmonic mean of precision and recall)
  f1_score <- if ((ppv + sensitivity) > 0) {
    2 * (ppv * sensitivity) / (ppv + sensitivity)
  } else {
    0
  }
  
  # Matthews Correlation Coefficient
  mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (mcc_denom > 0) {
    (tp * tn - fp * fn) / mcc_denom
  } else {
    0
  }
  
  # G-mean (geometric mean of sensitivity and specificity)
  g_mean <- sqrt(sensitivity * specificity)
  
  # Return all metrics
  c(
    ROC = auc_val,
    Sens = sensitivity,
    Spec = specificity,
    Accuracy = accuracy,
    MCC = mcc,
    NPV = npv,
    PPV = ppv,
    F1 = f1_score,
    Gmean = g_mean
  )
}


#' Get Gpredomics Model Specification for Caret
#'
#' @description
#' Returns the model specification list that can be passed directly to caret::train().
#' This is the recommended way to use gpredomics with caret.
#'
#' @return List with model specification compatible with caret
#' @export
#'
#' @examples
#' \dontrun{
#' library(caret)
#' library(gpredomicsR)
#' 
#' # Get model specification
#' gp_model <- getModelInfo_gpredomics()
#' 
#' # Use directly with train()
#' ctrl <- trainControl(method = "cv", number = 5)
#' model <- train(Species ~ ., data = iris_binary, 
#'                method = gp_model, trControl = ctrl)
#' }
getModelInfo_gpredomics <- function() {
  list(
    library = "gpredomicsR",
    type = "Classification",
    parameters = data.frame(
      parameter = c("population_size", "max_epochs", "ga_kmin", "ga_kmax", "beam_kmin", "beam_kmax", "language", "data_type"),
      class = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "character"),
      label = c("Population Size", "Max Epochs", "GA Min Features", "GA Max Features", "BEAM Min Features", "BEAM Max Features", "Model Language", "Data Type"),
      stringsAsFactors = FALSE
    ),
    grid = gpredomics_grid,
    fit = gpredomics_train,
    predict = gpredomics_predict,
    prob = gpredomics_prob,
    sort = gpredomics_sort,
    levels = function(x) x$levels
  )
}


#' Extract Best Model from Experiment
#'
#' @description
#' Internal helper to extract the best performing individual from a gpredomics experiment.
#' This function always returns the single best individual from the final population.
#'
#' @param experiment A gpredomics Experiment object (Rust object from fit_on)
#' @return The best Individual object from the population
#' @keywords internal
.extract_best_model <- function(experiment) {
  if (is.null(experiment)) {
    stop("Experiment object is NULL")
  }
  
  # Get the best population from the experiment
  best_pop <- tryCatch({
    if (!is.null(experiment$get_best_population)) {
      experiment$get_best_population()
    } else if (!is.null(experiment$get_population)) {
      experiment$get_population(NULL)
    } else {
      stop("Experiment object does not have get_population or get_best_population methods")
    }
  }, error = function(e) {
    stop(paste("Cannot extract population from experiment:", e$message))
  })
  
  if (is.null(best_pop)) {
    stop("No population found in experiment")
  }
  
  # Extract the best individual (index 1, using 1-based R indexing)
  best_individual <- tryCatch({
    best_pop$get_individual(1)
  }, error = function(e) {
    stop(paste("Cannot extract best individual:", e$message))
  })
  
  if (is.null(best_individual)) {
    stop("No best individual found in population")
  }
  
  return(best_individual)
}

#' Predict Scores from Best Model in Experiment
#'
#' @description
#' Internal helper to get raw scores from the best model in an experiment.
#' Returns the full prediction object with $class and $score components.
#'
#' @param experiment A gpredomics Experiment object
#' @param newdata A data.frame or matrix with new samples (features in columns)
#' @return List with $class and $score components from predict()
#' @keywords internal
.predict_scores_from_experiment <- function(experiment, newdata) {
  # Get best individual using Rust methods
  best_model <- .extract_best_model(experiment)
  
  if (is.null(best_model)) {
    stop("Could not extract best model from experiment")
  }
  
  # Convert newdata to data.frame if matrix
  if (is.matrix(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  
  # Create a dummy y vector (required by as_gpredomics_data but not used for prediction)
  n_samples <- nrow(newdata)
  dummy_y <- setNames(rep(0L, n_samples), rownames(newdata))
  if (is.null(names(dummy_y))) {
    names(dummy_y) <- paste0("S", seq_len(n_samples))
  }
  
  # Create Data object for prediction
  pred_data <- as_gpredomics_data(
    df = newdata,
    y_vec = dummy_y,
    features_in_columns = TRUE,
    prior_weight = NULL,
    feature_penalty = NULL,
    feature_tags = NULL,
    sample_tags = NULL
  )
  
  # Get predictions from best individual (returns list with $class and $score)
  predictions <- best_model$predict(pred_data)
  return(predictions)
}

#' Predict from Best Model in Experiment
#'
#' @description
#' Internal helper to generate predictions from the best model in an experiment.
#' Uses the best individual from the final population for predictions.
#'
#' @param experiment A gpredomics Experiment object
#' @param newdata A data.frame or matrix with new samples (features in columns)
#' @return Predictions as a factor vector
#' @keywords internal
.predict_from_experiment <- function(experiment, newdata) {
  # Get predictions (class and score)
  predictions <- .predict_scores_from_experiment(experiment, newdata)
  
  # predictions$class contains 0/1 classifications based on threshold
  class_predictions <- as.vector(predictions$class)
  return(as.factor(class_predictions))
}


#' Train Gpredomics Model via Caret
#'
#' @description
#' Training function compatible with caret's train() interface.
#' This function is called internally by caret during hyperparameter tuning.
#'
#' @param x Feature matrix (samples in rows, features in columns)
#' @param y Response vector (binary: 0/1)
#' @param wts Sample weights (currently not used)
#' @param param Data.frame with a single row of hyperparameters
#' @param lev Class levels
#' @param last Logical, whether this is the final model
#' @param weights Sample weights (currently not used)
#' @param classProbs Logical, whether to compute class probabilities
#' @param ... Additional arguments
#'
#' @return A list containing the trained experiment and metadata
#' @export
gpredomics_train <- function(x, y, wts = NULL, param, lev, last, weights, classProbs, ...) {
  
  # Convert y to integer 0/1 (not numeric/double!)
  if (is.factor(y)) {
    # For factors, convert to integer codes then adjust to 0-based
    y_int <- as.integer(y) - 1L
  } else {
    # Ensure it's integer type, not numeric/double
    y_int <- as.integer(y)
  }
  
  # Ensure binary classification
  if (!all(y_int %in% c(0L, 1L))) {
    stop("gpredomics currently only supports binary classification (0/1 or factor with 2 levels)")
  }
  
  # Convert x to data.frame if matrix
  if (is.matrix(x)) {
    x <- as.data.frame(x)
  }
  
  # Ensure x has column names
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }
  
  # Create named y vector (required by as_gpredomics_data)
  if (is.null(names(y_int))) {
    if (!is.null(rownames(x))) {
      names(y_int) <- rownames(x)
    } else {
      names(y_int) <- paste0("S", seq_len(length(y_int)))
    }
  }
  
  # Create Data object using the correct function and parameters
  gp_data <- as_gpredomics_data(
    df = x,
    y_vec = y_int,  # Pass integer vector, not numeric
    features_in_columns = TRUE,  # caret provides samples in rows, features in columns
    prior_weight = NULL,
    feature_penalty = NULL,
    feature_tags = NULL,
    sample_tags = NULL
  )
  
  # Create default Param object
  gp_param <- Param$new()
  
  # IMPORTANT: Set seed from R's current RNG state if not specified in param
  # This ensures reproducibility across caret CV folds
  if ("seed" %in% names(param) && !is.na(param$seed)) {
    gp_param$set("seed", as.integer(param$seed))
  } else {
    # Use R's current seed state (set by set.seed() or caret internally)
    # Sample a random integer to create a seed based on R's RNG state
    current_seed <- sample.int(.Machine$integer.max, 1)
    gp_param$set("seed", current_seed)
  }
  
  # Set numeric hyperparameters from param grid
  # General parameters
  if ("data_type_epsilon" %in% names(param)) gp_param$set("data_type_epsilon", param$data_type_epsilon)
  if ("k_penalty" %in% names(param)) gp_param$set("k_penalty", param$k_penalty)
  if ("fr_penalty" %in% names(param)) gp_param$set("fr_penalty", param$fr_penalty)
  if ("bias_penalty" %in% names(param)) gp_param$set("bias_penalty", param$bias_penalty)
  if ("threshold_ci_penalty" %in% names(param)) gp_param$set("threshold_ci_penalty", param$threshold_ci_penalty)
  if ("threshold_ci_alpha" %in% names(param)) gp_param$set("threshold_ci_alpha", param$threshold_ci_alpha)
  if ("threshold_ci_n_bootstrap" %in% names(param)) gp_param$set("threshold_ci_n_bootstrap", param$threshold_ci_n_bootstrap)
  if ("threshold_ci_frac_bootstrap" %in% names(param)) gp_param$set("threshold_ci_frac_bootstrap", param$threshold_ci_frac_bootstrap)
  
  # Data parameters
  if ("feature_minimal_prevalence_pct" %in% names(param)) gp_param$set("feature_minimal_prevalence_pct", param$feature_minimal_prevalence_pct)
  if ("feature_maximal_adj_pvalue" %in% names(param)) gp_param$set("feature_maximal_adj_pvalue", param$feature_maximal_adj_pvalue)
  if ("feature_minimal_feature_value" %in% names(param)) gp_param$set("feature_minimal_feature_value", param$feature_minimal_feature_value)
  if ("feature_minimal_log_abs_bayes_factor" %in% names(param)) gp_param$set("feature_minimal_log_abs_bayes_factor", param$feature_minimal_log_abs_bayes_factor)
  if ("max_features_per_class" %in% names(param)) gp_param$set("max_features_per_class", param$max_features_per_class)
  
  # GA parameters
  if ("population_size" %in% names(param)) gp_param$set("population_size", param$population_size)
  if ("max_epochs" %in% names(param)) gp_param$set("max_epochs", param$max_epochs)
  if ("min_epochs" %in% names(param)) gp_param$set("min_epochs", param$min_epochs)
  if ("max_age_best_model" %in% names(param)) gp_param$set("max_age_best_model", param$max_age_best_model)
  if ("ga_kmin" %in% names(param)) gp_param$set("k_min", param$ga_kmin)
  if ("ga_kmax" %in% names(param)) gp_param$set("k_max", param$ga_kmax)
  # Legacy support for kmin/kmax (map to GA parameters)
  if ("kmin" %in% names(param) && !("ga_kmin" %in% names(param))) gp_param$set("k_min", param$kmin)
  if ("kmax" %in% names(param) && !("ga_kmax" %in% names(param))) gp_param$set("k_max", param$kmax)
  if ("select_elite_pct" %in% names(param)) gp_param$set("select_elite_pct", param$select_elite_pct)
  if ("select_niche_pct" %in% names(param)) gp_param$set("select_niche_pct", param$select_niche_pct)
  if ("select_random_pct" %in% names(param)) gp_param$set("select_random_pct", param$select_random_pct)
  if ("mutated_children_pct" %in% names(param)) gp_param$set("mutated_children_pct", param$mutated_children_pct)
  if ("mutated_features_pct" %in% names(param)) gp_param$set("mutated_features_pct", param$mutated_features_pct)
  if ("mutation_non_null_chance_pct" %in% names(param)) gp_param$set("mutation_non_null_chance_pct", param$mutation_non_null_chance_pct)
  if ("forced_diversity_pct" %in% names(param)) gp_param$set("forced_diversity_pct", param$forced_diversity_pct)
  if ("forced_diversity_epochs" %in% names(param)) gp_param$set("forced_diversity_epochs", param$forced_diversity_epochs)
  if ("random_sampling_pct" %in% names(param)) gp_param$set("random_sampling_pct", param$random_sampling_pct)
  if ("random_sampling_epochs" %in% names(param)) gp_param$set("random_sampling_epochs", param$random_sampling_epochs)
  
  # BEAM parameters
  if ("beam_kmin" %in% names(param)) gp_param$set("k_min", param$beam_kmin)
  if ("beam_kmax" %in% names(param)) gp_param$set("k_max", param$beam_kmax)
  if ("best_models_criterion" %in% names(param)) gp_param$set("best_models_criterion", param$best_models_criterion)
  if ("max_nb_of_models" %in% names(param)) gp_param$set("max_nb_of_models", param$max_nb_of_models)
  if ("method" %in% names(param)) gp_param$set_string("method", as.character(param$method))
  
  # CV parameters
  if ("overfit_penalty" %in% names(param)) gp_param$set("overfit_penalty", param$overfit_penalty)
  if ("inner_folds" %in% names(param)) gp_param$set("inner_folds", param$inner_folds)
  if ("resampling_inner_folds_epochs" %in% names(param)) gp_param$set("resampling_inner_folds_epochs", param$resampling_inner_folds_epochs)
  if ("outer_folds" %in% names(param)) gp_param$set("outer_folds", param$outer_folds)
  if ("cv_best_models_ci_alpha" %in% names(param)) gp_param$set("cv_best_models_ci_alpha", param$cv_best_models_ci_alpha)
  
  # Voting parameters
  if ("min_perf" %in% names(param)) gp_param$set("min_perf", param$min_perf)
  if ("min_diversity" %in% names(param)) gp_param$set("min_diversity", param$min_diversity)
  if ("fbm_ci_alpha" %in% names(param)) gp_param$set("fbm_ci_alpha", param$fbm_ci_alpha)
  if ("method_threshold" %in% names(param)) gp_param$set("method_threshold", param$method_threshold)
  if ("threshold_windows_pct" %in% names(param)) gp_param$set("threshold_windows_pct", param$threshold_windows_pct)
  
  # String parameters
  if ("algo" %in% names(param)) gp_param$set_string("algo", as.character(param$algo))
  if ("language" %in% names(param)) gp_param$set_string("language", as.character(param$language))
  if ("data_type" %in% names(param)) gp_param$set_string("data_type", as.character(param$data_type))
  if ("fit" %in% names(param)) gp_param$set_string("fit", as.character(param$fit))
  if ("feature_selection_method" %in% names(param)) gp_param$set_string("feature_selection_method", as.character(param$feature_selection_method))
  if ("stratify_by" %in% names(param)) gp_param$set_string("stratify_by", as.character(param$stratify_by))
  
  # ============================================================================
  # IMPORTANT: Cross-validation and data splitting strategy
  # ============================================================================
  # We disable Gpredomics internal CV (cv = FALSE) because:
  # 1. Caret handles train/test splitting via trainControl (method = "cv", "boot", etc.)
  # 2. For each CV fold, caret passes only the training subset to gpredomics_train()
  # 3. Gpredomics trains on this subset (no internal splitting)
  # 4. Caret evaluates performance on the held-out fold
  # 5. Metrics are averaged across all CV folds
  #
  # This avoids nested CV (which would be cv inside cv) and ensures consistent
  # evaluation across all hyperparameter combinations.
  # ============================================================================
  
  gp_param$set_bool("cv", FALSE)
  
  # Reduce verbosity during grid search
  if (!last) {
    gp_param$set("n_model_to_display", 0)
  }
  
  # Create running flag
  running_flag <- RunningFlag$new()
  
  # Train model
  experiment <- fit_on(
    data = gp_data,
    param = gp_param,
    running_flag = running_flag
  )
  
  # Get number of features (k) from best individual for reporting
  k <- tryCatch({
    best_pop <- if (!is.null(experiment$get_best_population)) {
      experiment$get_best_population()
    } else {
      experiment$get_population(NULL)
    }
    best_ind <- best_pop$get_individual(1)
    ind_info <- best_ind$get()
    length(ind_info$features)
  }, error = function(e) {
    NA_integer_
  })
  
  # Extract prediction mode settings from param if present
  # (These are passed via train_gpredomics_caret wrapper)
  # Return model object
  list(
    experiment = experiment,
    param = param,
    features = colnames(x),
    levels = lev,
    k = k
  )
}


#' Predict from Gpredomics Model via Caret
#'
#' @description
#' Prediction function compatible with caret's train() interface.
#' This function is called internally by caret for generating predictions.
#'
#' @param modelFit Model object returned by gpredomics_train()
#' @param newdata New data for prediction (data.frame or matrix)
#' @param preProc Pre-processing object (currently not used)
#' @param submodels List of submodels (currently not used)
#' @param ... Additional arguments (not used)
#'
#' @return Predictions (factor for classification with levels 0 and 1)
#' @export
gpredomics_predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL, ...) {
  
  # Ensure newdata has the same features as training
  if (!all(modelFit$features %in% colnames(newdata))) {
    missing <- setdiff(modelFit$features, colnames(newdata))
    stop(paste("Missing features in newdata:", paste(missing, collapse = ", ")))
  }
  
  # Reorder columns to match training
  newdata <- newdata[, modelFit$features, drop = FALSE]
  
  # Generate predictions using best individual
  preds <- .predict_from_experiment(
    modelFit$experiment, 
    newdata
  )
  
  # Convert to factor with correct levels
  factor(preds, levels = c(0, 1), labels = modelFit$levels)
}


#' Predict Probabilities from Gpredomics Model via Caret
#'
#' @description
#' Returns probability matrix for caret compatibility using gpredomics scores.
#' The scores from predict()$score are converted to pseudo-probabilities:
#' - Score represents the model's confidence (higher score = more confident in class 1)
#' - Scores are normalized to range 0 to 1 to act as probabilities
#'
#' @param modelFit Model object returned by gpredomics_train()
#' @param newdata New data for prediction (data.frame or matrix)
#' @param preProc Pre-processing object (currently not used)
#' @param submodels List of submodels (currently not used)
#' @param ... Additional arguments (not used)
#'
#' @return Matrix with probability for each class (columns = class levels)
#' @export
gpredomics_prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL, ...) {
  
  # Ensure newdata has the same features as training
  if (!all(modelFit$features %in% colnames(newdata))) {
    missing <- setdiff(modelFit$features, colnames(newdata))
    stop(paste("Missing features in newdata:", paste(missing, collapse = ", ")))
  }
  
  # Reorder columns to match training
  newdata <- newdata[, modelFit$features, drop = FALSE]
  
  # Get predictions with scores
  predictions <- .predict_scores_from_experiment(
    modelFit$experiment,
    newdata
  )
  
  # Extract scores (higher score = more confident in positive class)
  scores <- as.vector(predictions$score)
  
  # Convert scores to probabilities using sigmoid-like transformation
  # This maps scores to [0, 1] range where:
  # - score >> 0 → prob(class=1) ≈ 1
  # - score << 0 → prob(class=1) ≈ 0
  # - score ≈ 0 → prob(class=1) ≈ 0.5
  prob_class1 <- 1 / (1 + exp(-scores))
  prob_class0 <- 1 - prob_class1
  
  # Create probability matrix with columns for each class level
  # Caret expects columns named after class levels
  prob_matrix <- cbind(prob_class0, prob_class1)
  colnames(prob_matrix) <- modelFit$levels
  
  return(prob_matrix)
}


#' Generate Tuning Grid for Gpredomics
#'
#' @description
#' Grid generation function compatible with caret's train() interface.
#' Supports both grid search and random search.
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param len Number of parameter combinations (for random search)
#' @param search Type of search: "grid" or "random"
#'
#' @return Data.frame with hyperparameter combinations
#' @export
gpredomics_grid <- function(x, y, len = NULL, search = "grid") {
  
  if (search == "grid") {
    # Grid search: predefined combinations
    # Adjust complexity based on len
    if (is.null(len)) len <- 3
    
    if (len == 1) {
      # Default parameters aligned with param.yaml defaults
      grid <- expand.grid(
        population_size = 5000,
        max_epochs = 100,
        ga_kmin = 1,
        ga_kmax = 200,
        beam_kmin = 1,
        beam_kmax = 200,
        language = "bin",
        data_type = "raw",
        stringsAsFactors = FALSE
      )
    } else if (len == 2) {
      # Small grid
      grid <- expand.grid(
        population_size = c(1000, 5000),
        max_epochs = c(50, 100),
        ga_kmin = c(1, 5),
        ga_kmax = c(50, 200),
        beam_kmin = c(1, 5),
        beam_kmax = c(50, 200),
        language = c("bin", "ter"),
        data_type = c("raw", "log"),
        stringsAsFactors = FALSE
      )
    } else {
      # Full grid
      grid <- expand.grid(
        population_size = c(1000, 5000, 10000),
        max_epochs = c(50, 100, 150),
        ga_kmin = c(1, 5, 10),
        ga_kmax = c(50, 100, 200),
        beam_kmin = c(1, 5, 10),
        beam_kmax = c(50, 100, 200),
        language = c("bin", "ter", "ratio", "pow2"),
        data_type = c("raw", "log", "prev"),
        stringsAsFactors = FALSE
      )
    }
    
  } else {
    # Random search: random sampling
    if (is.null(len)) len <- 10
    
    grid <- data.frame(
      population_size = sample(c(1000, 2000, 5000, 10000), len, replace = TRUE),
      max_epochs = sample(50:150, len, replace = TRUE),
      ga_kmin = sample(1:10, len, replace = TRUE),
      ga_kmax = sample(50:200, len, replace = TRUE),
      beam_kmin = sample(1:10, len, replace = TRUE),
      beam_kmax = sample(50:200, len, replace = TRUE),
      language = sample(c("bin", "ter", "ratio", "pow2"), len, replace = TRUE),
      data_type = sample(c("raw", "log", "prev"), len, replace = TRUE),
      stringsAsFactors = FALSE
    )
    
    # Ensure kmax > kmin for both GA and BEAM
    grid$ga_kmax <- pmax(grid$ga_kmax, grid$ga_kmin + 10)
    grid$beam_kmax <- pmax(grid$beam_kmax, grid$beam_kmin + 10)
  }
  
  grid
}


#' Sort Tuning Results for Gpredomics
#'
#' @description
#' Function to sort tuning results. Required by caret.
#'
#' @param x Data.frame with tuning results
#' @param metric Performance metric name (optional, default uses first numeric column)
#' @param maximize Logical, whether to maximize the metric (default TRUE for AUC/Accuracy)
#'
#' @return Sorted data.frame
#' @export
gpredomics_sort <- function(x, metric = NULL, maximize = TRUE) {
  # If no metric specified, use first numeric column
  if (is.null(metric)) {
    metric <- names(x)[sapply(x, is.numeric)][1]
  }
  
  if (maximize) {
    x[order(-x[, metric]), ]
  } else {
    x[order(x[, metric]), ]
  }
}


#' Get Gpredomics Model Info for Caret
#'
#' @description
#' Returns model information for caret registration.
#'
#' @return List with model metadata
#' @export
gpredomics_info <- function() {
  list(
    library = "gpredomicsR",
    type = "Classification",
    parameters = data.frame(
      parameter = c("population_size", "max_epochs", "ga_kmin", "ga_kmax", "beam_kmin", "beam_kmax", "language", "data_type"),
      class = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "character"),
      label = c("Population Size", "Max Epochs", "GA Min Features", "GA Max Features", "BEAM Min Features", "BEAM Max Features", "Model Language", "Data Type"),
      stringsAsFactors = FALSE
    ),
    grid = gpredomics_grid,
    fit = gpredomics_train,
    predict = gpredomics_predict,
    prob = gpredomics_prob,
    sort = gpredomics_sort,
    levels = function(x) x$levels
  )
}


#' Register Gpredomics with Caret
#'
#' @description
#' Registers the gpredomics model with caret for use in train().
#' After registration, you can use method = "gpredomics" directly in train().
#' 
#' Note: This is optional. You can also use getModelInfo_gpredomics() or
#' train_gpredomics_caret() without manual registration.
#'
#' @return NULL (invisible). Prints confirmation message.
#' @export
#'
#' @examples
#' \dontrun{
#' library(caret)
#' library(gpredomicsR)
#' 
#' # Register the model
#' register_gpredomics_caret()
#' 
#' # Now you can use method = "gpredomics"
#' data(iris)
#' iris_binary <- iris[iris$Species != "setosa", ]
#' iris_binary$Species <- droplevels(iris_binary$Species)
#' 
#' ctrl <- trainControl(
#'   method = "cv",
#'   number = 5,
#'   summaryFunction = gpredomicsSummary
#' )
#' 
#' model <- train(
#'   Species ~ .,
#'   data = iris_binary,
#'   method = "gpredomics",  # Can use string after registration
#'   trControl = ctrl,
#'   tuneLength = 2
#' )
#' }
register_gpredomics_caret <- function() {
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required but not installed. Install it with: install.packages('caret')")
  }
  
  # Create model specification
  gpredomics_model <- getModelInfo_gpredomics()
  
  # Try to add to caret's internal model list
  tryCatch({
    # Access caret's model list environment
    modelList <- get("modelList", envir = asNamespace("caret"))
    modelList$gpredomics <- gpredomics_model
    assign("modelList", modelList, envir = asNamespace("caret"))
    
    message("✓ Gpredomics successfully registered with caret")
    message("  You can now use method = 'gpredomics' in caret::train()")
  }, error = function(e) {
    message("ℹ Note: Could not register with caret's internal list.")
    message("  Use method = getModelInfo_gpredomics() or train_gpredomics_caret() instead.")
  })
  
  invisible(NULL)
}


#' Get Gpredomics Parameters Information for Caret
#'
#' @description
#' Returns comprehensive information about all parameters that can be tuned
#' via caret's grid search. This helps users discover which parameters are
#' available beyond the default grid.
#'
#' @param category Character, filter by category: "all" (default), "ga", "beam", 
#'   "mcmc", "general", "data", "cv", "voting", or "penalties"
#'
#' @return Data.frame with parameter information including name, type, default value,
#'   description, and category
#' @export
#'
#' @examples
#' \dontrun{
#' # Get all available parameters
#' all_params <- gpredomics_params_info()
#' View(all_params)
#' 
#' # Get only GA-specific parameters
#' ga_params <- gpredomics_params_info("ga")
#' 
#' # Get parameters related to penalties
#' penalty_params <- gpredomics_params_info("penalties")
#' }
gpredomics_params_info <- function(category = "all") {
  
  # Define all available parameters with their defaults and descriptions
  params <- data.frame(
    parameter = character(),
    type = character(),
    default = character(),
    description = character(),
    category = character(),
    stringsAsFactors = FALSE
  )
  
  # General parameters
  general <- data.frame(
    parameter = c("population_size", "max_epochs", "min_epochs", "max_age_best_model",
                  "language", "data_type", "algo", "fit", "threshold_ci_n_bootstrap", "threshold_ci_alpha"),
    type = c("numeric", "numeric", "numeric", "numeric", 
             "character", "character", "character", "character", "numeric", "numeric"),
    default = c("5000", "100", "1", "100",
                "bin,ratio,pow2,ter", "raw,prev,log", "ga", "auc", "0", "0.05"),
    description = c(
      "Target number of models per generation",
      "Maximum number of generations before stopping",
      "Minimum number of generations to run",
      "Stop if best model reaches this age (after min_epochs)",
      "Model language(s): ter, bin, ratio, pow2 (comma-separated)",
      "Data transformation(s): raw, prev, log (comma-separated)",
      "Algorithm: ga, beam, or mcmc",
      "Fitness metric: auc, mcc, sensitivity, specificity, f1_score, g_mean",
      "Number of bootstrap samples for threshold CI (0 = disabled)",
      "Alpha for threshold confidence interval"
    ),
    category = rep("general", 10),
    stringsAsFactors = FALSE
  )
  
  # Penalty parameters
  penalties <- data.frame(
    parameter = c("k_penalty", "fr_penalty", "bias_penalty", "threshold_ci_penalty"),
    type = rep("numeric", 4),
    default = c("0.0001", "0", "0", "0.5"),
    description = c(
      "Penalty for model complexity (multiplied by k)",
      "False rate penalty (for sensitivity/specificity fit)",
      "Penalty for models with sens/spec < 0.5",
      "Penalty based on threshold rejection rate"
    ),
    category = rep("penalties", 4),
    stringsAsFactors = FALSE
  )
  
  # GA-specific parameters
  ga <- data.frame(
    parameter = c("ga_kmin", "ga_kmax",
                  "select_elite_pct", "select_niche_pct", "select_random_pct",
                  "mutated_children_pct", "mutated_features_pct", "mutation_non_null_chance_pct",
                  "forced_diversity_pct", "forced_diversity_epochs",
                  "random_sampling_pct", "random_sampling_epochs"),
    type = rep("numeric", 12),
    default = c("1", "200",
                "2", "20", "10", "80", "20", "20", "0", "10", "0", "1"),
    description = c(
      "Minimal number of features in GA models",
      "Maximum number of features in GA models (0 = no limit)",
      "% of best models retained (lower = more elitist)",
      "% of best models retained per language/data_type",
      "% of random models retained",
      "% of children submitted to mutation",
      "% of mutation per gene/feature",
      "% chance of meaningful mutation (adding variable)",
      "% threshold for diversity filtering (0 = disabled)",
      "Epoch gap between diversity filters",
      "% of samples used per generation (0 = all samples)",
      "Number of epochs with same random sample subset"
    ),
    category = rep("ga", 12),
    stringsAsFactors = FALSE
  )
  
  # BEAM parameters  
  beam <- data.frame(
    parameter = c("beam_kmin", "beam_kmax", "best_models_criterion", "max_nb_of_models", "method"),
    type = c(rep("numeric", 4), "character"),
    default = c("1", "200", "10", "20000", "LimitedExhaustive"),
    description = c(
      "Minimal number of features in BEAM models",
      "Maximum number of features in BEAM models (0 = no limit)",
      "If ≤1: alpha for FBM CI; if >1: % of best models to keep",
      "Max number of models at each epoch",
      "Search method: LimitedExhaustive or ParallelForward"
    ),
    category = rep("beam", 5),
    stringsAsFactors = FALSE
  )
  
  # Cross-validation parameters
  cv <- data.frame(
    parameter = c("inner_folds", "overfit_penalty", "resampling_inner_folds_epochs",
                  "outer_folds", "cv_best_models_ci_alpha"),
    type = rep("numeric", 5),
    default = c("5", "0", "0", "5", "0.05"),
    description = c(
      "Number of inner CV folds for overfit penalty",
      "Penalty for overfitting based on CV performance",
      "Resplit inner folds every X epochs (0 = no resplit)",
      "Number of outer CV folds",
      "Alpha for FBM CI on validation fold"
    ),
    category = rep("cv", 5),
    stringsAsFactors = FALSE
  )
  
  # Voting parameters
  voting <- data.frame(
    parameter = c("fbm_ci_alpha", "min_perf", "min_diversity", 
                  "method_threshold", "threshold_windows_pct"),
    type = rep("numeric", 5),
    default = c("0.05", "0.50", "10", "0.5", "5"),
    description = c(
      "Alpha for Family of Best Models CI",
      "Required min sensitivity AND specificity for judges",
      "Required diversity between judges",
      "Voting threshold (0.5 for majority, 1 for consensus)",
      "% window around threshold for rejection"
    ),
    category = rep("voting", 5),
    stringsAsFactors = FALSE
  )
  
  # Data processing parameters
  data <- data.frame(
    parameter = c("data_type_epsilon", "feature_minimal_prevalence_pct",
                  "feature_maximal_adj_pvalue", "feature_minimal_feature_value",
                  "feature_minimal_log_abs_bayes_factor", "max_features_per_class",
                  "feature_selection_method"),
    type = c(rep("numeric", 6), "character"),
    default = c("1e-5", "10", "1", "1e-4", "2", "0", "wilcoxon"),
    description = c(
      "Threshold for log/prevalence transformation",
      "Min prevalence % per class to retain feature",
      "Max corrected p-value to retain feature",
      "Min mean value to retain feature",
      "Min log absolute Bayes factor (Bayesian method)",
      "Max features per class (0 = all significant)",
      "Method: wilcoxon, student_t, bayesian_fisher"
    ),
    category = rep("data", 7),
    stringsAsFactors = FALSE
  )
  
  # Combine all parameters
  all_params <- rbind(general, penalties, ga, beam, cv, voting, data)
  
  # Filter by category if requested
  if (category != "all") {
    if (!category %in% unique(all_params$category)) {
      stop(paste("Unknown category:", category, 
                 "\nAvailable categories:", 
                 paste(unique(all_params$category), collapse = ", ")))
    }
    all_params <- all_params[all_params$category == category, ]
  }
  
  # Sort by category then parameter name
  all_params <- all_params[order(all_params$category, all_params$parameter), ]
  rownames(all_params) <- NULL
  
  return(all_params)
}


#' Train Gpredomics Model with Caret (Convenience Function)
#'
#' @description
#' Convenience wrapper around caret::train() for gpredomics models.
#' Automatically uses the best individual from the population for predictions.
#' 
#' **Probability predictions**: Gpredomics scores are converted to pseudo-probabilities
#' using a sigmoid transformation, enabling ROC-based metrics. The raw scores from 
#' predict()$score represent model confidence and are transformed to range 0 to 1.
#'
#' @param formula Formula specifying the model
#' @param data Data.frame containing the data
#' @param method Character, tuning method: "grid" or "random"
#' @param tuneLength Integer, number of parameter combinations to try
#' @param tuneGrid Optional data.frame with custom parameter grid
#' @param trControl trainControl object (defaults to 5-fold CV with gpredomicsSummary)
#' @param metric Performance metric to optimize (default: "ROC"). 
#'   Options: "ROC", "Accuracy", "MCC", "Sens", "Spec", "F1", "Gmean", "NPV", "PPV"
#' @param ... Additional arguments passed to caret::train()
#'
#' @return A train object from caret
#' @export
#'
#' @examples
#' \dontrun{
#' # Simple example
#' data(iris)
#' iris_binary <- iris[iris$Species != "setosa", ]
#' iris_binary$Species <- droplevels(iris_binary$Species)
#' 
#' model <- train_gpredomics_caret(
#'   Species ~ .,
#'   data = iris_binary,
#'   method = "grid",
#'   tuneLength = 2,
#'   metric = "Accuracy"
#' )
#' 
#' # With custom control
#' ctrl <- trainControl(
#'   method = "repeatedcv",
#'   number = 10,
#'   repeats = 3,
#'   search = "random",
#'   classProbs = TRUE,  # Required for ROC/AUC metrics
#'   summaryFunction = gpredomicsSummary
#' )
#' 
#' model <- train_gpredomics_caret(
#'   Species ~ .,
#'   data = iris_binary,
#'   trControl = ctrl,
#'   tuneLength = 20,
#'   metric = "MCC"  # Matthews Correlation Coefficient
#' )
#' }
train_gpredomics_caret <- function(formula, 
                                   data, 
                                   method = "grid",
                                   tuneLength = 3,
                                   tuneGrid = NULL,
                                   trControl = NULL,
                                   metric = "ROC",
                                   ...) {
  
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required. Install it with: install.packages('caret')")
  }
  
  # Default trainControl if not provided
  if (is.null(trControl)) {
    trControl <- caret::trainControl(
      method = "cv",
      number = 5,
      search = method,
      classProbs = TRUE,  # Enable class probabilities for ROC and AUC
      summaryFunction = gpredomicsSummary  # Use gpredomicsSummary for comprehensive metrics
    )
  } else {
    # If user provides trControl but forgot classProbs, enable it for ROC/AUC
    if (is.null(trControl$classProbs) || !trControl$classProbs) {
      warning("classProbs was not set to TRUE in trControl. Enabling it for probability-based metrics (ROC, AUC).")
      trControl$classProbs <- TRUE
    }
  }
  
  # Get the model specification
  gpredomics_spec <- getModelInfo_gpredomics()
  
  # Call caret::train with the model specification
  model <- caret::train(
    form = formula,
    data = data,
    method = gpredomics_spec,  # Pass the list directly
    trControl = trControl,
    tuneLength = tuneLength,
    tuneGrid = tuneGrid,
    metric = metric,
    ...
  )
  
  return(model)
}
