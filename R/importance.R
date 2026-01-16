#' Feature Importance Analysis via K-Fold Cross-Validation
#'
#' @description Performs feature importance analysis by training models on k-fold cross-validation,


# ============================================================================
# HELPER FUNCTIONS (Global scope to avoid closure serialization issues)
# ============================================================================

#' Extract features and coefficients from a population
#' @keywords internal
extract_population_features <- function(population) {
  pop_size <- population$get_size()
  
  if (pop_size == 0) {
    return(list())
  }
  
  all_features <- list()
  
  for (i in 1:pop_size) {
    ind <- population$get_individual(i)
    
    if (is.null(ind)) {
      warning(sprintf("Individual %d is NULL in population", i))
      next
    }
    
    # Get individual data
    ind_data <- ind$get()
    
    # Extract features and coefficients
    features <- ind_data$features
    coefficients <- ind_data$coefficients
    
    if (is.null(features) || is.null(coefficients) || length(features) == 0) {
      next
    }
    
    # Store feature-coefficient pairs
    for (j in seq_along(features)) {
      feature_name <- features[j]
      coef_value <- coefficients[j]
      
      if (!feature_name %in% names(all_features)) {
        all_features[[feature_name]] <- list(coefficients = c())
      }
      
      all_features[[feature_name]]$coefficients <- c(
        all_features[[feature_name]]$coefficients,
        coef_value
      )
    }
  }
  
  return(all_features)
}

#' Run a single fold with a specific seed
#' @keywords internal
run_single_fold <- function(fold, seed = NULL, default_param, fold_data_list, name_prefix, 
                             glog_level, best_scope, best_criterion, fit_on_valid, 
                             param_yaml_path = NULL, parallel = FALSE, X = NULL, y = NULL, 
                             fold_indices = NULL, features_in_columns = TRUE, k = NULL) {
  
  # Set log level (Rust logger) - needed in each parallel worker
  try(GLogger$level(level = glog_level), silent = TRUE)
  
  # In parallel, load Param from YAML and create data objects in each worker
  if (parallel) {
    param <- Param$load(param_yaml_path)
    
    # Create data objects in this worker
    train_idx <- which(fold_indices != fold)
    valid_idx <- which(fold_indices == fold)
    
    # Subset X based on data orientation
    if (features_in_columns) {
      X_train_fold <- X[train_idx, , drop = FALSE]
      X_valid_fold <- X[valid_idx, , drop = FALSE]
      y_train_fold <- y[train_idx]
      y_valid_fold <- y[valid_idx]
      
      if (!is.null(rownames(X_train_fold))) {
        names(y_train_fold) <- rownames(X_train_fold)
      }
      if (!is.null(rownames(X_valid_fold))) {
        names(y_valid_fold) <- rownames(X_valid_fold)
      }
    } else {
      X_train_fold <- X[, train_idx, drop = FALSE]
      X_valid_fold <- X[, valid_idx, drop = FALSE]
      y_train_fold <- y[train_idx]
      y_valid_fold <- y[valid_idx]
      
      if (!is.null(colnames(X_train_fold))) {
        names(y_train_fold) <- colnames(X_train_fold)
      }
      if (!is.null(colnames(X_valid_fold))) {
        names(y_valid_fold) <- colnames(X_valid_fold)
      }
    }
    
    fold_data <- list(
      train_gdata = as.gpredomics.data(X_train_fold, y_train_fold, features_in_columns),
      valid_gdata = as.gpredomics.data(X_valid_fold, y_valid_fold, features_in_columns)
    )
  } else {
    param <- default_param$clone()
    fold_data <- fold_data_list[[fold]]
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    param$set("seed", seed)
  }
  
  train_gdata <- fold_data$train_gdata
  valid_gdata <- fold_data$valid_gdata
  
  message(sprintf("Running fold %d/%d%s", fold, k, 
                  if (!is.null(seed)) sprintf(" (seed=%d)", seed) else ""))
  
  tryCatch({
    # Train model
    running_flag <- RunningFlag$new()
    exp_fold <- fit_on(train_gdata, param, running_flag)
    
    # Get final population
    final_pop <- exp_fold$get_best_population()
    
    if (is.null(final_pop) || final_pop$get_size() == 0) {
      warning(sprintf("Fold %d: No final population found", fold))
      return(list(success = FALSE, features = list()))
    }
    
    # If fit_on_valid=TRUE, re-evaluate and re-sort population on validation data
    if (fit_on_valid) {
      message(sprintf("  Re-evaluating population on validation data for fold %d", fold))
      final_pop <- final_pop$fit(valid_gdata, param)
    }
    
    # Select best models based on best_scope
    selected_pop <- NULL
    
    if (best_scope == "best_ind") {
      # Get only the best individual
      best_ind <- final_pop$get_individual(1)
      
      if (is.null(best_ind)) {
        warning(sprintf("Fold %d: Best individual is NULL", fold))
        return(list(success = FALSE, features = list()))
      }
      
      # Extract features directly from this individual
      ind_data <- best_ind$get()
      features <- ind_data$features
      coefficients <- ind_data$coefficients
      
      feature_list <- list()
      if (!is.null(features) && !is.null(coefficients) && length(features) > 0) {
        for (j in seq_along(features)) {
          feature_name <- features[j]
          coef_value <- coefficients[j]
          
          if (!feature_name %in% names(feature_list)) {
            feature_list[[feature_name]] <- list(coefficients = c())
          }
          
          feature_list[[feature_name]]$coefficients <- c(
            feature_list[[feature_name]]$coefficients,
            coef_value
          )
        }
      }
      
      selected_features <- feature_list
      
    } else if (best_scope == "fbm") {
      # Use Family of Best Models
      if (is.null(best_criterion)) best_criterion <- 0.05
      
      selected_pop <- final_pop$get_fbm(best_criterion)
      
      if (selected_pop$get_size() == 0) {
        warning(sprintf("Fold %d: FBM selection returned empty population", fold))
        return(list(success = FALSE, features = list()))
      }
      
      selected_features <- extract_population_features(selected_pop)
      
    } else if (best_scope == "topn") {
      # Use top N individuals
      if (is.null(best_criterion)) best_criterion <- 10
      n_top <- as.integer(best_criterion)
      
      pop_size <- final_pop$get_size()
      if (n_top > pop_size) {
        warning(sprintf("Fold %d: topn (%d) exceeds population size (%d), using all individuals", 
                        fold, n_top, pop_size))
        n_top <- pop_size
      }
      
      if (n_top == 0) {
        warning(sprintf("Fold %d: topn is 0", fold))
        return(list(success = FALSE, features = list()))
      }
      
      selected_pop <- final_pop$get_first_n(n_top)
      selected_features <- extract_population_features(selected_pop)
      
    } else if (best_scope == "pct") {
      # Use top X% of population
      if (is.null(best_criterion)) best_criterion <- 5
      pct <- best_criterion
      
      selected_pop <- final_pop$get_first_pct(pct)
      
      if (selected_pop$get_size() == 0) {
        warning(sprintf("Fold %d: Top pct selection returned empty population", fold))
        return(list(success = FALSE, features = list()))
      }
      
      selected_features <- extract_population_features(selected_pop)
    }
    
    # Cleanup
    if (exists("exp_fold")) rm(exp_fold)
    if (exists("running_flag")) rm(running_flag)
    if (exists("train_gdata")) rm(train_gdata)
    if (exists("valid_gdata")) rm(valid_gdata)
    if (exists("final_pop")) rm(final_pop)
    if (exists("selected_pop") && !is.null(selected_pop)) rm(selected_pop)
    if (exists("best_ind") && !is.null(best_ind)) rm(best_ind)
    
    gc(verbose = FALSE)
    
    return(list(success = TRUE, features = selected_features))
    
  }, error = function(e) {
    warning(sprintf("Fold %d failed: %s", fold, e$message))
    
    # Cleanup even on error
    if (exists("exp_fold")) rm(exp_fold)
    if (exists("running_flag")) rm(running_flag)
    if (exists("train_gdata")) rm(train_gdata)
    if (exists("valid_gdata")) rm(valid_gdata)
    gc(verbose = FALSE)
    
    return(list(success = FALSE, features = list()))
  })
}

#' Run a fold with multiple seeds
#' @keywords internal
run_fold_multiple_seeds <- function(fold, n_seeds, default_param, fold_data_list, name_prefix, 
                                     glog_level, best_scope, best_criterion, fit_on_valid, 
                                     param_yaml_path = NULL, parallel = FALSE, X = NULL, y = NULL, 
                                     fold_indices = NULL, features_in_columns = TRUE, k = NULL,
                                     n_cores = NULL, n_threads_per_model = 1) {
  
  if (n_seeds == 1) {
    # Single run without seed override
    return(run_single_fold(fold, seed = NULL, default_param, fold_data_list, name_prefix, 
                           glog_level, best_scope, best_criterion, fit_on_valid, 
                           param_yaml_path, FALSE, X, y, fold_indices, features_in_columns, k))
  }
  
  # Run multiple times with different seeds - parallelize if requested
  if (parallel) {
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package 'foreach' is required for parallel execution. Install it with install.packages('foreach')")
    }
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is required for parallel execution. Install it with install.packages('doParallel')")
    }
    
    if (is.null(n_cores)) {
      n_cores <- max(1, floor(parallel::detectCores() / 2))
    }
    
    message(sprintf("  Parallelizing %d seeds with %d cores for fold %d", n_seeds, n_cores, fold))
    
    cl <- parallel::makeCluster(n_cores, outfile = "")
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # CRITICAL FIX: Don't export fold_data_list or default_param (External Pointers can't be serialized)
    # Workers will load param from YAML and recreate data objects locally from X and y
    # Only export functions - foreach auto-detects variables from the calling environment
    seed_results_list <- foreach::foreach(
      seed_idx = 1:n_seeds,
      .packages = c("gpredomicsR"),
      .export = c("run_single_fold", "extract_population_features"),
      .noexport = c("fold_data_list", "default_param")
    ) %dopar% {
      
      # Configure threads
      if (!is.null(n_threads_per_model) && n_threads_per_model > 0) {
        options(gpredomics.threads.number = n_threads_per_model)
      }
      
      seed_value <- 41 + seed_idx  # 42, 43, 44, ...
      
      # CRITICAL: Pass NULL for both fold_data_list and default_param - worker recreates locally
      result <- run_single_fold(fold, seed = seed_value, 
                                default_param = NULL,  # Worker loads from YAML
                                fold_data_list = NULL,  # Worker recreates from X, y
                                name_prefix, glog_level, best_scope, best_criterion, fit_on_valid, 
                                param_yaml_path, TRUE, X, y, fold_indices, features_in_columns, k)
      
      return(result)
    }
    
    # Cluster will be stopped automatically by on.exit handler
    
  } else {
    # Sequential execution of seeds
    seed_results_list <- list()
    
    for (seed_idx in 1:n_seeds) {
      seed_value <- 41 + seed_idx  # 42, 43, 44, ...
      
      result <- run_single_fold(fold, seed = seed_value, default_param, fold_data_list, 
                                name_prefix, glog_level, best_scope, best_criterion, fit_on_valid, 
                                param_yaml_path, FALSE, X, y, fold_indices, features_in_columns, k)
      
      seed_results_list[[seed_idx]] <- result
    }
  }
  
  # Aggregate results from all seeds
  all_features <- list()
  success_count <- 0
  
  for (seed_idx in 1:length(seed_results_list)) {
    result <- seed_results_list[[seed_idx]]
    
    if (result$success) {
      success_count <- success_count + 1
      
      # Merge features from this seed
      for (feature_name in names(result$features)) {
        if (!feature_name %in% names(all_features)) {
          all_features[[feature_name]] <- list(
            coefficients = c(),
            seeds_seen = integer()  # NEW: Track which seeds saw this feature
          )
        }
        
        # Number of times the feature appears in THIS seed (e.g., 10 times if topn=10)
        n_occurrences <- length(result$features[[feature_name]]$coefficients)
        
        all_features[[feature_name]]$coefficients <- c(
          all_features[[feature_name]]$coefficients,
          result$features[[feature_name]]$coefficients
        )
        
        # NEW: Record the seed ID for each occurrence
        # If this is seed_idx 1, and it has 10 models, we add ten "1"s.
        all_features[[feature_name]]$seeds_seen <- c(
          all_features[[feature_name]]$seeds_seen,
          rep(seed_idx, n_occurrences)
        )
      }
    }
  }
  
  if (success_count == 0) {
    return(list(success = FALSE, features = list()))
  }
  
  return(list(success = TRUE, features = all_features, n_seeds_success = success_count))
}

# ============================================================================
# MAIN FUNCTION
# ============================================================================

#' Feature Importance Analysis via K-Fold Cross-Validation
#'
#' @description Performs feature importance analysis by training models on k-fold cross-validation,
#' selecting the best models on validation folds, and aggregating feature prevalence and coefficients
#' across all selected models. This provides insights into which features are consistently important
#' across different data splits.
#'
#' @param default_param A Rust Param object (created with Param$load()) that provides the
#'   configuration for model training.
#' @param X Feature matrix for training (samples in rows, features in columns). Required.
#' @param y Response vector (same length as nrow(X)) for training labels. Required.
#' @param k Number of folds for k-fold cross-validation. Must be >= 2. Required.
#' @param n_seeds Number of times to run each fold with different random seeds.
#'   If > 1, each fold is run multiple times with deterministic seeds (42 + i),
#'   providing more robust feature importance estimates. Default is 1 (single run per fold).
#' @param fit_on_valid Logical. If FALSE (default), selects best models based on training performance.
#'   If TRUE, re-evaluates the final population on validation data and re-sorts it before selecting 
#'   best models. This allows selecting models that generalize best to the validation fold rather than
#'   those that perform best on training data. Recommended: TRUE for better generalization.
#' @param features_in_columns Whether features are in columns (TRUE, default) or rows (FALSE).
#' @param name_prefix Optional prefix for naming each experiment run. Default is "importance".
#' @param glog_level Verbosity level for logging. Default is "warn".
#' @param best_scope Method to select best models: "best_ind" (best individual only, default),
#'   "fbm" (Family of Best Models), "pct" (top X% individuals), or "topn" (top N individuals).
#'   The "jury" option is not recommended for importance analysis as it aggregates predictions.
#' @param best_criterion Numeric value controlling the scope for "fbm", "pct", and "topn" best_scope options.
#'   For "fbm": alpha value for confidence interval (default 0.05). For "pct": percentage of top
#'   individuals to include (default 5). For "topn": number of top individuals to include (default 10).
#'   Ignored for "best_ind" scope.
#' @param parallel Whether to parallelize seed execution within each fold (requires doParallel). 
#'   Default is FALSE. Note: folds are always executed sequentially, but if parallel=TRUE and n_seeds>1,
#'   the different seeds for each fold will be executed in parallel. This is more memory-efficient 
#'   than parallelizing folds.
#' @param n_cores Number of cores to use if parallel = TRUE. Default is half of available cores.
#'   Make sure n_cores * n_threads_per_model <= total available cores to avoid oversubscription.
#' @param log_file Path to a log file where execution progress will be written. Default is "importance_log.txt".
#'   Set to NULL or "" to disable file logging.
#' @param n_threads_per_model Number of threads to use for each model (Rust level). Default is 1.
#'   In parallel mode, this controls how many threads each worker uses. Be careful with total thread count.
#'
#' @return A data frame with feature importance information, containing:
#'   \describe{
#'     \item{feature}{Feature name}
#'     \item{prevalence_folds}{Number of different folds in which the feature appears (cross-patient generalization)}
#'     \item{stability_seed_avg}{Average reproducibility across seeds within folds (0-1 scale). NORMALIZED metric
#'       based on unique seeds: measures the probability that a feature is discovered when running a seed in a fold.
#'       This is comparable across different best_scope settings (best_ind vs topn). High values (>0.8) indicate
#'       robust detection; low values suggest the feature is found only by chance in specific random initializations.
#'       Example: 0.5 means the feature was found by 50% of seeds when present in a fold.}
#'     \item{intra_seed_density}{Average number of models per seed containing this feature (when found). This metric
#'       measures saturation/intensity within a seed. ~1.0 for best_ind (one model per seed), ~10 for topn=10 if
#'       the feature saturates all selected models. Helps distinguish "consistently top-ranked" (high density) from
#'       "occasionally selected" (low density) features when using topn/pct/fbm selection methods.}
#'     \item{prevalence_pct}{Percentage of folds containing this feature (prevalence_folds / k * 100)}
#'     \item{sign_consistency}{Directional stability: score from 0 (chaotic, equal pos/neg) to 1 (always same sign).
#'       Computed as |prevalence_positive - prevalence_negative| / prevalence_total. Critical for biological interpretation:
#'       high values indicate consistent biological effect, low values suggest confounding or noise.}
#'     \item{mean_abs_coef}{Mean absolute coefficient: average magnitude when feature is present (coef_abs_sum / prevalence_total).
#'       Distinguishes "strong but rare drivers" from "frequent background noise".}
#'     \item{prevalence_total}{Total number of times the feature appears in selected models across all folds and seeds}
#'     \item{prevalence_positive}{Number of times the feature appears with a positive coefficient (>0)}
#'     \item{prevalence_negative}{Number of times the feature appears with a negative coefficient (<0)}
#'     \item{coef_abs_sum}{Sum of absolute values of all coefficients for this feature}
#'     \item{coef_total_sum}{Sum of all coefficients (positive + negative) for this feature}
#'     \item{coef_positive_sum}{Sum of all positive coefficients for this feature}
#'     \item{coef_negative_sum}{Sum of all negative coefficients for this feature}
#'   }
#'   The data frame is sorted by: (1) prevalence_folds (cross-patient generalization), (2) stability_seed_avg 
#'   (technical reproducibility - now normalized), (3) sign_consistency (directional stability), (4) coef_abs_sum (magnitude). 
#'   Features with high prevalence_folds, high stability_seed_avg AND high sign_consistency are the most reliable biomarkers.
#'   The normalized stability_seed_avg makes results comparable between best_ind and topn/pct/fbm selection methods.
#'
#' @details
#' The function:
#' 1. Creates stratified k-fold cross-validation splits
#' 2. For each fold:
#'    - Trains models on k-1 folds (repeated n_seeds times with different seeds)
#'    - If fit_on_valid=TRUE: Re-evaluates the final population on validation data and re-sorts by validation performance
#'    - If fit_on_valid=FALSE: Uses the population as-is (sorted by training performance)
#'    - Selects best models using the specified best_scope method (best_ind, fbm, pct, or topn)
#'    - Extracts features and coefficients from selected models
#' 3. Aggregates feature statistics across all folds and seeds
#' 4. Returns a sorted data frame with feature importance metrics
#'
#' When fit_on_valid=TRUE, the selection prioritizes models that generalize well to unseen data,
#' which typically results in more robust and stable feature importance estimates.
#'
#' @examples
#' \dontrun{
#' # Load default parameters
#' default_param <- Param$load("sample/param.yaml")
#'
#' # Run importance analysis with 5-fold cross-validation
#' importance_df <- compute_feature_importance(
#'   default_param = default_param,
#'   X = X_train,
#'   y = y_train,
#'   k = 5,
#'   fit_on_valid = TRUE  # Select models based on validation performance
#' )
#'
#' # View top features
#' head(importance_df)
#'
#' # Run with multiple seeds for more robust estimates
#' importance_robust <- compute_feature_importance(
#'   default_param = default_param,
#'   X = X_train,
#'   y = y_train,
#'   k = 5,
#'   n_seeds = 5,  # Run each fold 5 times
#'   fit_on_valid = TRUE
#' )
#'
#' # Use top 10% of models for importance analysis
#' importance_top_pct <- compute_feature_importance(
#'   default_param = default_param,
#'   X = X_train,
#'   y = y_train,
#'   k = 5,
#'   best_scope = "pct",
#'   best_criterion = 10,  # Top 10% of models
#'   fit_on_valid = TRUE
#' )
#' }
#'
#' @importFrom foreach %dopar%
#' @author Generated for gpredomicsR
#' @export
compute_feature_importance <- function(default_param,
                                       X,
                                       y,
                                       k,
                                       n_seeds = 1,
                                       fit_on_valid = FALSE,
                                       features_in_columns = TRUE,
                                       name_prefix = "importance",
                                       glog_level = "warn",
                                       best_scope = c("best_ind", "fbm", "pct", "topn"),
                                       best_criterion = NULL,
                                       parallel = FALSE,
                                       n_cores = NULL,
                                       log_file = "importance_log.txt",
                                       n_threads_per_model = 1) {
  
  # Start timing
  start_time <- Sys.time()
  
  # Initialize log file
  if (!is.null(log_file) && nchar(log_file) > 0) {
    cat(sprintf("=== Feature Importance Analysis Started: %s ===\n", Sys.time()), file = log_file)
    cat(sprintf("Log file: %s\n", log_file), file = log_file, append = TRUE)
  }
  
  best_scope <- match.arg(best_scope)
  
  # Validate n_seeds
  if (!is.numeric(n_seeds) || n_seeds < 1 || n_seeds != floor(n_seeds)) {
    stop("n_seeds must be a positive integer >= 1")
  }
  
  if (n_seeds > 1) {
    message(sprintf("Running each fold %d times with different seeds (42 to %d)", 
                    n_seeds, 41 + n_seeds))
  }
  
  # Set log level globally at the start (Rust logger)
  try(GLogger$level(level = glog_level), silent = TRUE)
  
  # Validate inputs
  if (!inherits(default_param, "Param")) {
    stop("default_param must be a Rust Param object (use Param$load() to create one)")
  }
  
  # Set default best_criterion based on best_scope
  if (is.null(best_criterion)) {
    best_criterion <- switch(best_scope,
                              fbm = 0.05,      # Default alpha for FBM
                              pct = 5,         # Default 5% for pct
                              topn = 10,       # Default 10 individuals for topn
                              best_ind = NA    # Not used for best_ind
    )
  }
  
  # Validate best_criterion
  if (best_scope %in% c("fbm", "pct", "topn") && (is.null(best_criterion) || is.na(best_criterion))) {
    stop(sprintf("best_criterion must be specified for best_scope = '%s'", best_scope))
  }
  
  # Validate X and y - required
  if (is.null(X) || is.null(y)) {
    stop("X and y are required. Provide training data as matrices/data.frames.")
  }
  
  if (features_in_columns && (nrow(X) != length(y))) {
    stop("X and y must have the same number of samples (nrow(X) must equal length(y))")
  } else if (!features_in_columns && (ncol(X) != length(y))) {
    stop("X and y must have the same number of samples (ncol(X) must equal length(y))")
  }
  
  # Validate k parameter - required and must be >= 2
  if (is.null(k) || !is.numeric(k) || k < 2) {
    stop("k must be a numeric value >= 2 for cross-validation")
  }
  
  n_samples <- if (features_in_columns) nrow(X) else ncol(X)
  if (k > n_samples) {
    stop(sprintf("k (%d) cannot be greater than the number of samples (%d)", k, n_samples))
  }
  
  # Create stratified k-fold cross-validation
  message(sprintf("Setting up %d-fold stratified cross-validation", k))
  
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required for stratified k-fold cross-validation. Install it with install.packages('caret')")
  }
  
  # Create stratified folds using caret
  fold_list <- caret::createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  
  # Convert to vector format (fold number for each sample)
  fold_indices <- integer(n_samples)
  for (i in 1:k) {
    fold_indices[fold_list[[i]]] <- i
  }
  
  message(sprintf("Fold indices created for %d folds", k))
  
  # Log configuration
  if (!is.null(log_file) && nchar(log_file) > 0) {
    cat(sprintf("Number of folds: %d\n", k), file = log_file, append = TRUE)
    cat(sprintf("Number of seeds per fold: %d\n", n_seeds), file = log_file, append = TRUE)
    cat(sprintf("Fit on validation: %s\n", fit_on_valid), file = log_file, append = TRUE)
    cat(sprintf("Best scope: %s\n", best_scope), file = log_file, append = TRUE)
    cat(sprintf("Parallel execution: %s\n", parallel), file = log_file, append = TRUE)
    if (parallel && !is.null(n_cores)) {
      cat(sprintf("Number of cores: %d\n", n_cores), file = log_file, append = TRUE)
    }
    cat(sprintf("\n--- Starting Execution ---\n\n"), file = log_file, append = TRUE)
  }
  
  # Create temporary YAML file for parameter passing (used in parallel mode)
  param_yaml_path <- NULL
  if (parallel && n_seeds > 1) {
    param_yaml_path <- default_param$save(tempfile(pattern = "gpred_param_", fileext = ".yaml"))
    
    # Note: On Windows, temp files may be locked by workers; unlink might fail silently (harmless)
    on.exit(if(!is.null(param_yaml_path) && file.exists(param_yaml_path)) {
      try(unlink(param_yaml_path), silent = TRUE)
    }, add = TRUE)
  }
  
  # Pre-create gpredomics.data objects for all folds
  # These will be used in sequential mode or passed to parallel workers
  message("Pre-creating data objects for all folds...")
  fold_data_list <- vector("list", k)
  
  for (fold in 1:k) {
    train_idx <- which(fold_indices != fold)
    valid_idx <- which(fold_indices == fold)
    
    # Subset X based on data orientation
    if (features_in_columns) {
      # Samples in rows, features in columns - subset rows
      X_train_fold <- X[train_idx, , drop = FALSE]
      X_valid_fold <- X[valid_idx, , drop = FALSE]
      y_train_fold <- y[train_idx]
      y_valid_fold <- y[valid_idx]
      
      # Align names using row names
      if (!is.null(rownames(X_train_fold))) {
        names(y_train_fold) <- rownames(X_train_fold)
      }
      if (!is.null(rownames(X_valid_fold))) {
        names(y_valid_fold) <- rownames(X_valid_fold)
      }
    } else {
      # Features in rows, samples in columns - subset columns
      X_train_fold <- X[, train_idx, drop = FALSE]
      X_valid_fold <- X[, valid_idx, drop = FALSE]
      y_train_fold <- y[train_idx]
      y_valid_fold <- y[valid_idx]
      
      # Align names using column names
      if (!is.null(colnames(X_train_fold))) {
        names(y_train_fold) <- colnames(X_train_fold)
      }
      if (!is.null(colnames(X_valid_fold))) {
        names(y_valid_fold) <- colnames(X_valid_fold)
      }
    }
    
    # Create and store data objects
    fold_data_list[[fold]] <- list(
      train_gdata = as.gpredomics.data(X_train_fold, y_train_fold, features_in_columns),
      valid_gdata = as.gpredomics.data(X_valid_fold, y_valid_fold, features_in_columns)
    )
  }
  
  message(sprintf("Data objects created for %d folds", k))
  
  # Helper functions are now defined at global scope to avoid closure capture
  # (extract_population_features, run_single_fold, run_fold_multiple_seeds)
  
  # Execute importance analysis - folds are sequential, seeds are parallelized
  if (!is.null(n_threads_per_model) && n_threads_per_model > 0 && !parallel) {
    options(gpredomics.threads.number = n_threads_per_model)
  }
  
  if (parallel && n_seeds == 1) {
    message("Note: parallel=TRUE but n_seeds=1. No parallelization will occur. Consider increasing n_seeds.")
  }
  
  results_list <- lapply(1:k, function(fold) {
    timestamp <- format(Sys.time(), "%H:%M:%S")
    
    start_msg <- sprintf("[%s] START Fold %d/%d", timestamp, fold, k)
    cat(paste0(start_msg, "\n"))
    if (!is.null(log_file) && nchar(log_file) > 0) {
      try(cat(paste0(start_msg, "\n"), file = log_file, append = TRUE), silent = TRUE)
    }
    
    # run_fold_multiple_seeds will parallelize seeds if parallel=TRUE and n_seeds>1
    res <- run_fold_multiple_seeds(fold, n_seeds, default_param, fold_data_list, name_prefix, 
                                    glog_level, best_scope, best_criterion, fit_on_valid, 
                                    param_yaml_path = param_yaml_path, parallel = parallel,
                                    X = X, y = y, fold_indices = fold_indices,
                                    features_in_columns = features_in_columns, k = k,
                                    n_cores = n_cores, n_threads_per_model = n_threads_per_model)
    
    end_timestamp <- format(Sys.time(), "%H:%M:%S")
    status <- if (res$success) "SUCCESS" else "FAILURE"
    end_msg <- sprintf("[%s] FINISH Fold %d/%d | Status: %s", end_timestamp, fold, k, status)
    
    cat(paste0(end_msg, "\n"))
    if (!is.null(log_file) && nchar(log_file) > 0) {
      try(cat(paste0(end_msg, "\n"), file = log_file, append = TRUE), silent = TRUE)
    }
    
    return(res)
  })
  
  # Aggregate features across all folds
  message("Aggregating feature importance across all folds...")
  
  all_features <- list()
  success_count <- 0
  
  for (i in seq_along(results_list)) {
    res <- results_list[[i]]  # Result from fold 'i' (may contain multiple seeds)
    
    if (res$success) {
      success_count <- success_count + 1
      
      # Merge features from this fold
      for (feature_name in names(res$features)) {
        if (!feature_name %in% names(all_features)) {
          all_features[[feature_name]] <- list(
            coefficients = c(),
            folds_seen_raw = integer(),  # Store ALL fold occurrences (e.g., 1,1,1,2,2)
            unique_seeds_count_by_fold = numeric()  # NEW: Ratio of unique seeds per fold
          )
        }
        
        # Coefficients from this fold (may be multiple if n_seeds > 1)
        new_coeffs <- res$features[[feature_name]]$coefficients
        count_in_this_fold <- length(new_coeffs)
        
        # Add coefficients from this fold
        all_features[[feature_name]]$coefficients <- c(
          all_features[[feature_name]]$coefficients,
          new_coeffs
        )
        
        # CRITICAL CORRECTION:
        # If the feature appears 5 times in fold 'i' (via 5 seeds), we add 'i' 5 times.
        # This allows us to calculate the density per fold later.
        all_features[[feature_name]]$folds_seen_raw <- c(
          all_features[[feature_name]]$folds_seen_raw,
          rep(i, count_in_this_fold)
        )
        
        # KEY CALCULATION: How many UNIQUE seeds saw the feature in this fold?
        # res$features[[feature_name]]$seeds_seen contains e.g., c(1, 1, 1, 2, 2, 5)
        # unique(...) -> c(1, 2, 5) -> length -> 3 unique seeds.
        if (!is.null(res$features[[feature_name]]$seeds_seen)) {
          n_unique_seeds <- length(unique(res$features[[feature_name]]$seeds_seen))
        } else {
          # Fallback for single seed case (no seeds_seen tracking)
          n_unique_seeds <- 1
        }
        
        # Store this number (e.g., 3) for this fold.
        # Divide by n_seeds (e.g., 10) to get a proper 0-1 ratio.
        ratio_for_this_fold <- n_unique_seeds / n_seeds
        
        all_features[[feature_name]]$unique_seeds_count_by_fold <- c(
          all_features[[feature_name]]$unique_seeds_count_by_fold,
          ratio_for_this_fold
        )
      }
    }
  }
  
  message(sprintf("Importance analysis completed: %d/%d folds successful", success_count, k))
  
  if (success_count == 0) {
    warning("All folds failed. Returning empty data frame.")
    return(data.frame(
      feature = character(),
      prevalence_folds = integer(),
      stability_seed_avg = numeric(),
      intra_seed_density = numeric(),
      prevalence_pct = numeric(),
      sign_consistency = numeric(),
      mean_abs_coef = numeric(),
      prevalence_total = integer(),
      prevalence_positive = integer(),
      prevalence_negative = integer(),
      coef_abs_sum = numeric(),
      coef_total_sum = numeric(),
      coef_positive_sum = numeric(),
      coef_negative_sum = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Build final importance data frame
  message("Building feature importance data frame...")
  
  importance_list <- list()
  
  for (feature_name in names(all_features)) {
    coeffs <- all_features[[feature_name]]$coefficients
    folds_raw <- all_features[[feature_name]]$folds_seen_raw
    
    # 1. Analyze stability per fold (How many seeds found the feature in each fold?)
    # table(folds_raw) gives us: Fold 1 -> 10 times, Fold 2 -> 3 times...
    counts_per_fold <- table(folds_raw)
    
    # Number of UNIQUE folds (this is your previous prevalence_folds)
    unique_folds_count <- length(counts_per_fold)
    
    # 2. Calculate Seed Stability (CORRECTED AND NORMALIZED)
    # We directly stored the ratios (0 to 1) calculated properly per fold
    # These ratios are based on UNIQUE seeds, making them comparable across best_scope
    ratios <- all_features[[feature_name]]$unique_seeds_count_by_fold
    
    # stability_seed_avg: Average of ratios across folds where the feature appears
    # If TopN=10 or BestInd=1, this number is now comparable!
    # 1.0 means "Every time we ran seeds in this fold, we found it (at least once)"
    stability_seed_avg <- mean(ratios)
    
    # 3. Standard metrics
    prevalence_total <- length(coeffs)
    prevalence_pct <- (unique_folds_count / k) * 100
    prevalence_positive <- sum(coeffs > 0)
    prevalence_negative <- sum(coeffs < 0)
    coef_positive_sum <- sum(coeffs[coeffs > 0])
    coef_negative_sum <- sum(coeffs[coeffs < 0])
    coef_total_sum <- sum(coeffs)
    coef_abs_sum <- sum(abs(coeffs))
    
    # Sign consistency (0 = random, 1 = perfect)
    sign_consistency <- abs(prevalence_positive - prevalence_negative) / prevalence_total
    
    # Mean absolute coefficient (impact when present)
    mean_abs_coef <- coef_abs_sum / prevalence_total
    
    # 4. NEW: Intra-seed density (intensity)
    # How many models per seed contain this feature on average?
    # Formula: prevalence_total / (prevalence_folds × stability_seed_avg × n_seeds)
    # This tells us: when a seed finds the feature, how saturated is it?
    # ~1.0 for best_ind, ~10 for topn=10 if feature is in all selected models
    denominator <- unique_folds_count * stability_seed_avg * n_seeds
    if (denominator > 0) {
      intra_seed_density <- prevalence_total / denominator
    } else {
      intra_seed_density <- NA
    }
    
    importance_list[[feature_name]] <- list(
      feature = feature_name,
      prevalence_folds = unique_folds_count,  # Generalization (Cross-Patient)
      stability_seed_avg = stability_seed_avg,  # Reproducibility (Cross-Run) - NOW NORMALIZED
      intra_seed_density = intra_seed_density,  # Intensity (models per seed)
      prevalence_pct = prevalence_pct,
      sign_consistency = sign_consistency,
      mean_abs_coef = mean_abs_coef,
      prevalence_total = prevalence_total,
      prevalence_positive = prevalence_positive,
      prevalence_negative = prevalence_negative,
      coef_abs_sum = coef_abs_sum,
      coef_total_sum = coef_total_sum,
      coef_positive_sum = coef_positive_sum,
      coef_negative_sum = coef_negative_sum
    )
  }
  
  # Convert to data frame
  importance_df <- do.call(rbind, lapply(importance_list, function(x) {
    data.frame(
      feature = x$feature,
      prevalence_folds = x$prevalence_folds,
      stability_seed_avg = x$stability_seed_avg,
      intra_seed_density = x$intra_seed_density,
      prevalence_pct = x$prevalence_pct,
      sign_consistency = x$sign_consistency,
      mean_abs_coef = x$mean_abs_coef,
      prevalence_total = x$prevalence_total,
      prevalence_positive = x$prevalence_positive,
      prevalence_negative = x$prevalence_negative,
      coef_abs_sum = x$coef_abs_sum,
      coef_total_sum = x$coef_total_sum,
      coef_positive_sum = x$coef_positive_sum,
      coef_negative_sum = x$coef_negative_sum,
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(importance_df) <- NULL
  
  # Reorganize columns for readability
  importance_df <- importance_df[, c(
    "feature", 
    "prevalence_folds",        # Generalization (Cross-Patient) - The King
    "stability_seed_avg",      # Reproducibility (Cross-Run) - The Duke (NOW NORMALIZED)
    "intra_seed_density",      # Intensity (models per seed when found)
    "prevalence_pct",          # Percentage of folds
    "sign_consistency",        # Reliability (Directional Stability) - The Queen
    "mean_abs_coef",           # Impact (Magnitude) - The Force
    "prevalence_total",        # Raw frequency
    "prevalence_positive", "prevalence_negative",
    "coef_abs_sum", "coef_total_sum", "coef_positive_sum", "coef_negative_sum"
  )]
  
  # INTELLIGENT SORT:
  # 1. First, features that generalize across most folds (Biological Robustness)
  # 2. Then, features that are most reproducible within those folds (Technical Stability)
  # 3. Then, directional consistency (same sign)
  # 4. Finally, total magnitude
  importance_df <- importance_df[order(importance_df$prevalence_folds, 
                                       importance_df$stability_seed_avg,
                                       importance_df$sign_consistency,
                                       importance_df$coef_abs_sum, 
                                       decreasing = TRUE), ]
  
  rownames(importance_df) <- NULL
  
  # Calculate execution time
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Log completion
  if (!is.null(log_file) && nchar(log_file) > 0) {
    cat(sprintf("\n=== Feature Importance Analysis Completed: %s ===\n", end_time), 
        file = log_file, append = TRUE)
    cat(sprintf("Total execution time: %.2f seconds (%.2f minutes)\n", 
                execution_time, execution_time / 60), file = log_file, append = TRUE)
    cat(sprintf("Successful folds: %d/%d\n", success_count, k), file = log_file, append = TRUE)
    cat(sprintf("Total features found: %d\n", nrow(importance_df)), file = log_file, append = TRUE)
  }
  
  message(sprintf("\nFeature importance analysis complete. Found %d features.", nrow(importance_df)))
  message(sprintf("Execution time: %.2f seconds (%.2f minutes)", execution_time, execution_time / 60))
  
  return(importance_df)
}


#' Repeated Feature Importance Analysis via K-Fold Cross-Validation
#'
#' @description Performs repeated feature importance analysis by running k-fold cross-validation
#' multiple times with different data splits (different seeds for fold creation). This approach,
#' known as Repeated K-Fold Cross-Validation, reduces bias from specific data splits and provides
#' more robust feature importance estimates, especially critical for small datasets where the random
#' assignment of samples to folds can significantly impact results.
#'
#' @param n_repeats Number of times to repeat the entire k-fold CV process with different data splits.
#'   Each repeat uses a different seed for fold creation. Default is 10. Higher values (20-30) provide
#'   more stable estimates but increase computation time. For small datasets (N < 200), 20+ repeats
#'   are recommended.
#' @param base_seed Starting seed for fold creation. Each repeat will use base_seed + i (where i is
#'   the repeat number). Default is 1000. This ensures reproducibility across runs.
#' @param ... All other arguments passed to \code{\link{compute_feature_importance}}. This includes:
#'   default_param, X, y, k, n_seeds, fit_on_valid, features_in_columns, name_prefix, glog_level,
#'   best_scope, best_criterion, parallel, n_cores, log_file, n_threads_per_model.
#'
#' @return A data frame with aggregated feature importance information across all repeats, containing:
#'   \describe{
#'     \item{feature}{Feature name}
#'     \item{prevalence_folds_mean}{Mean number of folds where feature appears across repeats}
#'     \item{prevalence_folds_sd}{Standard deviation of prevalence_folds across repeats}
#'     \item{stability_seed_avg_mean}{Mean seed stability across repeats}
#'     \item{stability_seed_avg_sd}{Standard deviation of seed stability across repeats}
#'     \item{prevalence_pct_mean}{Mean percentage of folds containing the feature}
#'     \item{sign_consistency_mean}{Mean sign consistency across repeats}
#'     \item{sign_consistency_sd}{Standard deviation of sign consistency}
#'     \item{mean_abs_coef_mean}{Mean of the mean absolute coefficients across repeats}
#'     \item{mean_abs_coef_sd}{Standard deviation of mean absolute coefficients}
#'     \item{prevalence_total_mean}{Mean total prevalence across repeats}
#'     \item{n_repeats_seen}{Number of repeats (out of n_repeats) where the feature appeared at least once}
#'     \item{repeat_stability}{Proportion of repeats where feature was detected (n_repeats_seen / n_repeats)}
#'     \item{coef_direction_majority}{Predominant direction: "positive", "negative", or "mixed" (if inconsistent)}
#'   }
#'   The data frame is sorted by: (1) repeat_stability (features appearing in most repeats), 
#'   (2) prevalence_folds_mean (cross-patient generalization), (3) stability_seed_avg_mean 
#'   (technical reproducibility), (4) sign_consistency_mean (directional stability).
#'   
#'   Features with high repeat_stability, high prevalence_folds_mean, and high sign_consistency_mean
#'   are the most robust biomarkers, having survived multiple independent data splits.
#'
#' @details
#' This function implements Repeated K-Fold Cross-Validation for feature importance:
#' 
#' 1. For each of n_repeats iterations:
#'    - Sets a different seed for fold creation (base_seed + i)
#'    - Calls compute_feature_importance() which creates different k-fold splits
#'    - Collects feature importance metrics from each repeat
#' 
#' 2. Aggregates results across all repeats:
#'    - Calculates mean and standard deviation for key metrics
#'    - Computes repeat_stability: proportion of repeats detecting each feature
#'    - Determines predominant coefficient direction across repeats
#' 
#' 3. Returns consolidated importance estimates
#' 
#' Why this is critical for small datasets (N < 200):
#' - Reduces "split bias": features shouldn't be important only because of lucky fold assignments
#' - Identifies truly robust biomarkers that generalize across different data partitions
#' - Provides uncertainty estimates (standard deviations) for all metrics
#' - Features surviving multiple repeats are more likely to replicate in independent studies
#' 
#' Computational cost: n_repeats × cost of single run. Use parallel=TRUE and consider running
#' overnight for thorough analysis with high n_repeats.
#'
#' @examples
#' \dontrun{
#' # Load default parameters
#' default_param <- Param$load("sample/param.yaml")
#'
#' # Run repeated importance analysis (20 repeats for robust estimates)
#' importance_repeated <- compute_feature_importance_repeated(
#'   n_repeats = 20,
#'   base_seed = 1000,
#'   default_param = default_param,
#'   X = X_train,
#'   y = y_train,
#'   k = 5,
#'   n_seeds = 5,
#'   fit_on_valid = TRUE,
#'   parallel = TRUE,
#'   n_cores = 8
#' )
#'
#' # View top robust features (appearing in most repeats)
#' head(importance_repeated)
#'
#' # Filter for highly stable features (present in >80% of repeats)
#' robust_features <- importance_repeated[importance_repeated$repeat_stability > 0.8, ]
#'
#' # Check features with consistent direction and high stability
#' gold_standard <- importance_repeated[
#'   importance_repeated$repeat_stability > 0.9 & 
#'   importance_repeated$sign_consistency_mean > 0.8,
#' ]
#' }
#'
#' @seealso \code{\link{compute_feature_importance}}
#' @author Generated for gpredomicsR
#' @export
compute_feature_importance_repeated <- function(n_repeats = 10,
                                                base_seed = 1000,
                                                ...) {
  
  # Start timing
  start_time <- Sys.time()
  
  # Validate inputs
  if (!is.numeric(n_repeats) || n_repeats < 1 || n_repeats != floor(n_repeats)) {
    stop("n_repeats must be a positive integer >= 1")
  }
  
  if (n_repeats == 1) {
    message("n_repeats = 1. Consider using compute_feature_importance() directly instead.")
  }
  
  message(sprintf("=== Starting Repeated Feature Importance Analysis ==="))
  message(sprintf("Number of repeats: %d", n_repeats))
  message(sprintf("Base seed: %d (using seeds %d to %d for fold creation)", 
                  base_seed, base_seed + 1, base_seed + n_repeats))
  
  # Store results from each repeat
  repeat_results <- list()
  
  # Run importance analysis n_repeats times with different fold splits
  for (i in 1:n_repeats) {
    current_seed <- base_seed + i
    
    message(sprintf("\n--- Repeat %d/%d (fold creation seed: %d) ---", i, n_repeats, current_seed))
    
    # Set seed for fold creation (this affects caret::createFolds)
    set.seed(current_seed)
    
    # Run single importance analysis
    tryCatch({
      result <- compute_feature_importance(...)
      
      # Add repeat identifier
      result$repeat_id <- i
      
      repeat_results[[i]] <- result
      
      message(sprintf("Repeat %d completed: %d features found", i, nrow(result)))
      
    }, error = function(e) {
      warning(sprintf("Repeat %d failed: %s", i, e$message))
      repeat_results[[i]] <- NULL
    })
  }
  
  # Filter out failed repeats
  repeat_results <- Filter(Negate(is.null), repeat_results)
  
  if (length(repeat_results) == 0) {
    stop("All repeats failed. Cannot aggregate results.")
  }
  
  n_successful <- length(repeat_results)
  message(sprintf("\n=== Aggregating results from %d successful repeats ===", n_successful))
  
  # Combine all results
  all_features_data <- do.call(rbind, repeat_results)
  
  # Get unique features across all repeats
  unique_features <- unique(all_features_data$feature)
  
  message(sprintf("Total unique features found: %d", length(unique_features)))
  
  # Aggregate metrics for each feature
  aggregated_list <- list()
  
  for (feat in unique_features) {
    feat_data <- all_features_data[all_features_data$feature == feat, ]
    
    n_repeats_seen <- nrow(feat_data)
    repeat_stability <- n_repeats_seen / n_successful
    
    # Calculate mean and sd for key metrics
    prevalence_folds_mean <- mean(feat_data$prevalence_folds, na.rm = TRUE)
    prevalence_folds_sd <- sd(feat_data$prevalence_folds, na.rm = TRUE)
    
    stability_seed_avg_mean <- mean(feat_data$stability_seed_avg, na.rm = TRUE)
    stability_seed_avg_sd <- sd(feat_data$stability_seed_avg, na.rm = TRUE)
    
    prevalence_pct_mean <- mean(feat_data$prevalence_pct, na.rm = TRUE)
    
    sign_consistency_mean <- mean(feat_data$sign_consistency, na.rm = TRUE)
    sign_consistency_sd <- sd(feat_data$sign_consistency, na.rm = TRUE)
    
    mean_abs_coef_mean <- mean(feat_data$mean_abs_coef, na.rm = TRUE)
    mean_abs_coef_sd <- sd(feat_data$mean_abs_coef, na.rm = TRUE)
    
    prevalence_total_mean <- mean(feat_data$prevalence_total, na.rm = TRUE)
    
    # Determine predominant direction
    total_pos <- sum(feat_data$prevalence_positive, na.rm = TRUE)
    total_neg <- sum(feat_data$prevalence_negative, na.rm = TRUE)
    
    if (total_pos > total_neg * 2) {
      coef_direction_majority <- "positive"
    } else if (total_neg > total_pos * 2) {
      coef_direction_majority <- "negative"
    } else {
      coef_direction_majority <- "mixed"
    }
    
    aggregated_list[[feat]] <- list(
      feature = feat,
      prevalence_folds_mean = prevalence_folds_mean,
      prevalence_folds_sd = prevalence_folds_sd,
      stability_seed_avg_mean = stability_seed_avg_mean,
      stability_seed_avg_sd = stability_seed_avg_sd,
      prevalence_pct_mean = prevalence_pct_mean,
      sign_consistency_mean = sign_consistency_mean,
      sign_consistency_sd = sign_consistency_sd,
      mean_abs_coef_mean = mean_abs_coef_mean,
      mean_abs_coef_sd = mean_abs_coef_sd,
      prevalence_total_mean = prevalence_total_mean,
      n_repeats_seen = n_repeats_seen,
      repeat_stability = repeat_stability,
      coef_direction_majority = coef_direction_majority
    )
  }
  
  # Convert to data frame
  aggregated_df <- do.call(rbind, lapply(aggregated_list, function(x) {
    data.frame(
      feature = x$feature,
      prevalence_folds_mean = x$prevalence_folds_mean,
      prevalence_folds_sd = x$prevalence_folds_sd,
      stability_seed_avg_mean = x$stability_seed_avg_mean,
      stability_seed_avg_sd = x$stability_seed_avg_sd,
      prevalence_pct_mean = x$prevalence_pct_mean,
      sign_consistency_mean = x$sign_consistency_mean,
      sign_consistency_sd = x$sign_consistency_sd,
      mean_abs_coef_mean = x$mean_abs_coef_mean,
      mean_abs_coef_sd = x$mean_abs_coef_sd,
      prevalence_total_mean = x$prevalence_total_mean,
      n_repeats_seen = x$n_repeats_seen,
      repeat_stability = x$repeat_stability,
      coef_direction_majority = x$coef_direction_majority,
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(aggregated_df) <- NULL
  
  # Replace NA values with 0 (can occur when feature appears in only one repeat)
  # This provides cleaner output while maintaining data integrity
  aggregated_df[is.na(aggregated_df)] <- 0
  
  # Sort by: (1) repeat stability, (2) fold prevalence, (3) seed stability, (4) sign consistency
  aggregated_df <- aggregated_df[order(aggregated_df$repeat_stability,
                                        aggregated_df$prevalence_folds_mean,
                                        aggregated_df$stability_seed_avg_mean,
                                        aggregated_df$sign_consistency_mean,
                                        decreasing = TRUE), ]
  
  rownames(aggregated_df) <- NULL
  
  # Calculate total execution time
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  message(sprintf("\n=== Repeated Feature Importance Analysis Completed ==="))
  message(sprintf("Total execution time: %.2f seconds (%.2f minutes)", total_time, total_time / 60))
  message(sprintf("Successful repeats: %d/%d", n_successful, n_repeats))
  message(sprintf("Total unique features: %d", nrow(aggregated_df)))
  message(sprintf("Features in >80%% of repeats: %d", sum(aggregated_df$repeat_stability > 0.8)))
  message(sprintf("Features in >50%% of repeats: %d", sum(aggregated_df$repeat_stability > 0.5)))
  
  return(aggregated_df)
}
