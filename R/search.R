#' Grid Search for Gpredomics Hyperparameter Tuning
#'
#' @description Performs a comprehensive grid search over all modifiable parameters in gpredomics.
#' Each parameter combination from the grid is merged with default parameters from a
#' Rust Param object, and gpredomics is executed for each unique configuration.
#'
#' @param grid A data frame created with expand.grid() containing all parameter combinations
#'   to test. Each column represents a parameter, and each row is a unique combination.
#'   Column names must match parameter names in the Param object.
#' @param default_param A Rust Param object (created with Param$load()) that provides the
#'   baseline configuration. Only parameters specified in `grid` will be modified.
#' @param name_prefix Optional prefix for naming each experiment run. Default is "grid_search".
#' @param glog_level Verbosity level for logging. Default is "warn".
#' @param metric Metric name to sort results by. Default is "auc". Only metrics from compute_metrics are allowed:
#'   auc, accuracy, sensitivity, specificity, rejection_rate, mcc, npv, ppv, f1_score, g_mean.
#'   The function automatically prepends "valid_" to the metric name (e.g., "auc" becomes "valid_auc").
#'   You can also specify "train_" prefix explicitly. Sorting is done in descending order (higher is better).
#' @param parallel Whether to run grid search in parallel (requires doParallel). Default is FALSE.
#' @param n_cores Number of cores to use if parallel = TRUE. Default is half of available cores.
#' @param best_scope Method to compute performance metrics: "jury" (default, uses voting/jury system),
#'   "best_ind" (best individual only), "fbm" (Family of Best Models),
#'   "pct" (average performance of top X% individuals), or "topn" (average performance of top N individuals).
#'   For "fbm", "pct", and "topn", the alpha/percentage/count value is controlled by best_criterion.
#' @param best_criterion Numeric value controlling the scope for "fbm", "pct", and "topn" best_scope options.
#'   For "fbm": alpha value for confidence interval (default 0.05). For "pct": percentage of top
#'   individuals to average (default 5). For "topn": number of top individuals to average (default 10).
#'   Ignored for "best_ind" and "jury" scopes.
#' @param aggregation Aggregation method for "fbm" and "pct" scopes: "mean" (default) or "median".
#' @param X Feature matrix for training (samples in rows, features in columns). Required.
#' @param y Response vector (same length as nrow(X)) for training labels. Required.
#' @param k Number of folds for k-fold cross-validation. If NULL or 1, no cross-validation is performed
#'   (models are validated on training data). Default is NULL.
#' @param cv_aggregation Aggregation method for cross-validation results: "mean" (default) or "median".
#' @param features_in_columns Whether features are in columns (TRUE, default) or rows (FALSE).
#'
#' @importFrom foreach %dopar%
#'
#' @return A list of class "gpredomics_grid_search" containing:
#'   \describe{
#'     \item{results}{Data frame with all parameter combinations and their metrics, sorted by the specified metric}
#'     \item{best}{Experiment object trained with best parameters on full training data (NULL if training failed)}
#'     \item{best_params}{Named list containing the best parameter combination and its metrics}
#'     \item{metric}{The metric column name used for sorting}
#'     \item{n_success}{Number of successful runs}
#'     \item{n_total}{Number of parameter combinations tested}
#'     \item{n_total_grid}{Total number of combinations in the grid (before random sampling)}
#'     \item{execution_time}{Total execution time in seconds}
#'     \item{random_search}{Logical indicating if random search was used}
#'   }
#'
#' @details
#' The function:
#' 1. Takes a Rust Param object as the baseline configuration
#' 2. Iterates over each row of the expand.grid data frame
#' 3. For each combination, modifies the default parameters accordingly
#' 4. Saves a temporary YAML file and runs gpredomics
#' 5. Collects metrics and results from each run
#' 6. Returns aggregated results with the best configuration identified
#'
#' @param get_train If TRUE, includes train_* columns in the result. Default is FALSE (only valid_* columns are returned).
#' @param random_search_n If > 0, performs random search by sampling this many parameter combinations
#'   from the grid instead of testing all combinations. This is often more efficient than exhaustive
#'   grid search, especially for large grids. Default is 0 (exhaustive grid search).
#' @param n_seeds Number of times to run each parameter combination with different random seeds.
#'   If > 1, each configuration is run multiple times with deterministic seeds (42 + i),
#'   and the performance metrics are averaged across seeds. This improves result robustness
#'   and reduces variance from random initialization. Default is 1 (single run per configuration).
#'   When n_seeds > 1, standard deviation columns (e.g., valid_auc_sd) are also included in results.
#' @param log_file Path to a log file where execution progress will be written. Default is "grid_search_log.txt".
#'   Set to NULL or "" to disable file logging. Logs include start/finish times for each configuration,
#'   worker IDs (in parallel mode), parameter values, and success/failure status. This provides real-time
#'   monitoring of long-running grid searches.
#' @param n_threads_per_model Number of threads to use for each model (Rust level). Default is 1.
#'   In parallel mode, this controls how many threads each worker uses. Be careful: if you use n_cores=10
#'   workers and n_threads_per_model=4, you'll spawn 40 threads total, which can cause context switching
#'   overhead. Best practice: keep this at 1 in parallel mode, or use formula: n_cores * n_threads_per_model <= total_cores.
#'   This sets the gpredomics.threads.number option for the Rust backend.
#'
#' @examples
#' \dontrun{
#' # Load default parameters
#' default_param <- Param$load("sample/param.yaml")
#'
#' # Create grid with expand.grid
#' my_grid <- expand.grid(
#'   population_size = c(1000, 5000, 10000),
#'   max_epochs = c(50, 100),
#'   language = c("bin", "ter", "ratio"),
#'   data_type = c("raw", "log"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Run grid search with 5-fold cross-validation
#' results <- grid_search(
#'   grid = my_grid,
#'   default_param = default_param,
#'   X = X_train,
#'   y = y_train,
#'   k = 5,
#'   metric = "auc"  # automatically uses "valid_auc" column
#' )
#'
#' # View results (sorted by metric)
#' print(head(results$results))
#'
#' # Best parameters and metrics are in results$best_params
#' print(results$best_params)
#'
#' # Access best metric value
#' print(results$best_params$valid_auc)
#'
#' # Use the trained best model for predictions
#' predictions <- results$best$get_best_population()$get_individual(1)$predict(test_data)
#'
#' # For more robust results, run each configuration multiple times with different seeds
#' results_robust <- grid_search(
#'   grid = my_grid,
#'   default_param = default_param,
#'   X = X_train,
#'   y = y_train,
#'   k = 5,
#'   metric = "auc",
#'   n_seeds = 5  # Run each config 5 times (seeds 42-46) and average the results
#' )
#'
#' # Results will include mean metrics and standard deviations (e.g., valid_auc_sd)
#' head(results_robust$results)
#' }
#'
#' @author Generated for gpredomicsR
#' @export
grid_search <- function(grid,
                        default_param,
                        X,
                        y,
                        k = NULL,
                        cv_aggregation = c("mean", "median"),
                        features_in_columns = TRUE,
                        name_prefix = "grid_search",
                        glog_level = "warn",
                        metric = "auc",
                        parallel = FALSE,
                        n_cores = NULL,
                        best_scope = c("best_ind", "fbm", "pct", "topn", "jury"),
                        best_criterion = 5,
                        aggregation = c("mean", "median"),
                        get_train = FALSE,
                        random_search_n = 0,
                        n_seeds = 1,
                        log_file = "grid_search_log.txt",
                        n_threads_per_model = 1) {
  
  # Start timing
  start_time <- Sys.time()
  
  # Initialize log file
  if (!is.null(log_file) && nchar(log_file) > 0) {
    cat(sprintf("=== Grid Search Started: %s ===\n", Sys.time()), file = log_file)
    cat(sprintf("Log file: %s\n", log_file), file = log_file, append = TRUE)
  }
  
  best_scope <- match.arg(best_scope)
  aggregation <- match.arg(aggregation)
  cv_aggregation <- match.arg(cv_aggregation)
  
  # Validate n_seeds
  if (!is.numeric(n_seeds) || n_seeds < 1 || n_seeds != floor(n_seeds)) {
    stop("n_seeds must be a positive integer >= 1")
  }
  
  if (n_seeds > 1) {
    message(sprintf("Running each parameter combination %d times with different seeds (42 to %d)", 
                    n_seeds, 41 + n_seeds))
  }
  
  # Set log level globally at the start (Rust logger)
  try(GLogger$level(level = glog_level), silent = TRUE)
  
  # Validate inputs
  if (!is.data.frame(grid)) {
    stop("grid must be a data frame created with expand.grid()")
  }
  
  if (nrow(grid) == 0) {
    stop("grid cannot be empty")
  }
  
  if (!inherits(default_param, "Param")) {
    stop("default_param must be a Rust Param object (use Param$load() to create one)")
  }
  
  # Set default best_criterion based on best_scope
  if (is.null(best_criterion)) {
    best_criterion <- switch(best_scope,
                              fbm = 0.05,      # Default alpha for FBM
                              pct = 5,         # Default 5% for pct
                              topn = 10,       # Default 10 individuals for topn
                              jury = NA,       # Not used for jury
                              best_ind = NA    # Not used for best_ind
    )
                        
  }
  
  # Validate best_criterion
  if (best_scope %in% c("fbm", "pct", "topn") && (is.null(best_criterion) || is.na(best_criterion))) {
    stop(sprintf("best_criterion must be specified for best_scope = '%s'", best_scope))
  }
  
  # Validate X and y - now required
  if (is.null(X) || is.null(y)) {
    stop("X and y are required for grid_search. Provide training data as matrices/data.frames.")
  }
  
  if (features_in_columns && (nrow(X) != length(y))) {
    stop("X and y must have the same number of samples (nrow(X) must equal length(y))")
  } else if (!features_in_columns && (ncol(X) != length(y))) {
    stop("X and y must have the same number of samples (ncol(X) must equal length(y))")
  }
  
  # Validate k parameter
  if (!is.null(k)) {
    n_samples_check <- if (features_in_columns) nrow(X) else ncol(X)
    if (!is.numeric(k) || k < 2) {
      stop("k must be a numeric value >= 2 for cross-validation")
    }
    if (k > n_samples_check) {
      stop(sprintf("k (%d) cannot be greater than the number of samples (%d)", k, n_samples_check))
    }
  }
  
  # Create fold indices if k-fold CV is requested
  fold_indices <- NULL
  fold_data_list <- NULL
  
  if (!is.null(k) && k >= 2) {
    message(sprintf("Setting up %d-fold stratified cross-validation", k))
    
    if (!requireNamespace("caret", quietly = TRUE)) {
      stop("Package 'caret' is required for stratified k-fold cross-validation. Install it with install.packages('caret')")
    }
    
    # Create stratified folds using caret
    fold_list <- caret::createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
    
    # Convert to vector format (fold number for each sample)
    n_samples <- if (features_in_columns) nrow(X) else ncol(X)
    fold_indices <- integer(n_samples)
    for (i in 1:k) {
      fold_indices[fold_list[[i]]] <- i
    }
    
    message(sprintf("Fold indices created for %d folds", k))
  } else {
    message("No cross-validation - models will be validated on training data")
  }
  
  # For non-parallel execution, pre-create data objects
  param_yaml_path <- NULL
  if (!parallel) {
    if (!is.null(k) && k >= 2) {
      # Pre-create gpredomics.data objects for all folds
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
    } else {
      # Pre-create single data object for no-CV case
      # Ensure names are aligned based on features_in_columns
      if (features_in_columns) {
        # Samples in rows, features in columns
        if (!is.null(rownames(X))) {
          names(y) <- rownames(X)
        }
      } else {
        # Features in rows, samples in columns
        if (!is.null(colnames(X))) {
          names(y) <- colnames(X)
        }
      }
      
      fold_data_list <- list(
        list(
          train_gdata = as.gpredomics.data(X, y, features_in_columns),
          valid_gdata = as.gpredomics.data(X, y, features_in_columns)
        )
      )
    }
  } else {
    # For parallel execution, create a temporary YAML file
    param_yaml_path <- default_param$save(tempfile(pattern = "gpred_param_", fileext = ".yaml"))
    
    on.exit(if(!is.null(param_yaml_path) && file.exists(param_yaml_path)) {
      unlink(param_yaml_path)
    }, add = TRUE)
  }
  
  # Use grid as-is (already expanded)
  grid_df <- grid
  n_combinations_total <- nrow(grid_df)
  
  # Random search: sample n rows if random_search_n > 0
  if (random_search_n > 0) {
    if (random_search_n >= n_combinations_total) {
      message(sprintf("random_search_n (%d) >= total combinations (%d), using full grid", 
                      random_search_n, n_combinations_total))
      n_combinations <- n_combinations_total
    } else {
      message(sprintf("Random search: sampling %d combinations from %d total", 
                      random_search_n, n_combinations_total))
      sampled_indices <- sample(1:n_combinations_total, random_search_n, replace = FALSE)
      grid_df <- grid_df[sampled_indices, , drop = FALSE]
      n_combinations <- random_search_n
    }
  } else {
    n_combinations <- n_combinations_total
  }
  
  message(sprintf("Exploration will test %d parameter combinations", n_combinations))
  
  # Log total configurations
  if (!is.null(log_file) && nchar(log_file) > 0) {
    cat(sprintf("Total configurations to test: %d\n", n_combinations), file = log_file, append = TRUE)
    cat(sprintf("Parallel execution: %s\n", parallel), file = log_file, append = TRUE)
    if (parallel && !is.null(n_cores)) {
      cat(sprintf("Number of cores: %d\n", n_cores), file = log_file, append = TRUE)
    }
    cat(sprintf("\n--- Starting Execution ---\n\n"), file = log_file, append = TRUE)
  }
  
  # Function to run a single experiment with given parameters and a specific seed
  run_single_experiment <- function(i, grid_row, default_param, fold_data_list, name_prefix, glog_level, best_scope, best_criterion, aggregation, cv_aggregation, param_yaml_path = NULL, parallel = FALSE, X = NULL, y = NULL, fold_indices = NULL, features_in_columns = TRUE, k = NULL, seed = NULL) {
    # Set log level (Rust logger) - needed in each parallel worker
    try(GLogger$level(level = glog_level), silent = TRUE)
    
    # In parallel, load Param from YAML in each worker
    if (parallel) {
      param <- Param$load(param_yaml_path)
      
      # Create data objects in this worker (can't serialize external pointers)
      if (!is.null(k) && k >= 2) {
        fold_data_list <- vector("list", k)
        
        for (fold in 1:k) {
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
          
          fold_data_list[[fold]] <- list(
            train_gdata = as.gpredomics.data(X_train_fold, y_train_fold, features_in_columns),
            valid_gdata = as.gpredomics.data(X_valid_fold, y_valid_fold, features_in_columns)
          )
        }
      } else {
        # No CV case
        if (features_in_columns) {
          if (!is.null(rownames(X))) {
            names(y) <- rownames(X)
          }
        } else {
          if (!is.null(colnames(X))) {
            names(y) <- colnames(X)
          }
        }
        
        fold_data_list <- list(
          list(
            train_gdata = as.gpredomics.data(X, y, features_in_columns),
            valid_gdata = as.gpredomics.data(X, y, features_in_columns)
          )
        )
      }
    } else {
      param <- default_param$clone()
    }
    
    # Set seed if provided
    if (!is.null(seed)) {
      param$set("seed", seed)
    }
    
    # Enable jury/voting if best_scope is jury
    if (best_scope == "jury") {
      param$set_bool("vote", TRUE)
    }
    # Update parameters based on grid row
    for (param_name in names(grid_row)) {
      param_value <- grid_row[[param_name]]
      if (is.factor(param_value)) {
        param_value <- as.character(param_value)
      }
      if (is.logical(param_value)) {
        param$set_bool(param_name, param_value)
      } else if (is.character(param_value)) {
        param$set_string(param_name, param_value)
      } else if (is.numeric(param_value)) {
        param$set(param_name, param_value)
      } else {
        warning(sprintf("Parameter '%s' has unsupported type, skipping", param_name))
      }
    }
    experiment_name <- sprintf("%s_%04d", name_prefix, i)
    message(sprintf("[%d/%d] Running experiment: %s", i, nrow(grid_df), experiment_name))
    tryCatch({
      if (length(fold_data_list) > 1) {
        k_folds <- length(fold_data_list)
        fold_metrics_list <- list()
        for (fold in 1:k_folds) {
          message(sprintf("  Fold %d/%d", fold, k_folds))
          fold_data <- fold_data_list[[fold]]
          train_gdata <- fold_data$train_gdata
          valid_gdata <- fold_data$valid_gdata
          running_flag <- RunningFlag$new()
          exp_fold <- fit_on(train_gdata, param, running_flag)
          fold_metrics <- extract_metrics(exp_fold, best_scope = best_scope, best_criterion = best_criterion, aggregation = aggregation, train_gdata = train_gdata, valid_data_obj = valid_gdata)
          fold_metrics_list[[fold]] <- fold_metrics
        }
        metrics <- aggregate_cv_metrics(fold_metrics_list, cv_aggregation)
      } else {
        fold_data <- fold_data_list[[1]]
        train_gdata <- fold_data$train_gdata
        valid_gdata <- fold_data$valid_gdata
        running_flag <- RunningFlag$new()
        exp <- fit_on(train_gdata, param, running_flag)
        metrics <- extract_metrics(exp, best_scope = best_scope, best_criterion = best_criterion, aggregation = aggregation, train_gdata = train_gdata, valid_data_obj = valid_gdata)
      }
      # Explicit memory cleanup to release Rust resources
      # Remove heavy objects created during the experiment
      if (exists("exp")) rm(exp)
      if (exists("exp_fold")) rm(exp_fold)
      if (exists("running_flag")) rm(running_flag)
      if (exists("train_gdata")) rm(train_gdata)
      if (exists("valid_gdata")) rm(valid_gdata)
      if (exists("fold_data")) rm(fold_data)
      
      # Force garbage collector to immediately release Rust memory
      gc(verbose = FALSE)
      
      result <- c(as.list(grid_row), metrics, list(experiment_name = experiment_name, success = TRUE, error = NA))
      return(list(params = grid_row, metrics = metrics, success = TRUE))
    }, error = function(e) {
      warning(sprintf("Experiment %d failed: %s", i, e$message))
      
      # Cleanup even on error to avoid memory leaks
      if (exists("exp")) rm(exp)
      if (exists("exp_fold")) rm(exp_fold)
      if (exists("running_flag")) rm(running_flag)
      if (exists("train_gdata")) rm(train_gdata)
      if (exists("valid_gdata")) rm(valid_gdata)
      if (exists("fold_data")) rm(fold_data)
      gc(verbose = FALSE)
      
      return(list(params = grid_row, metrics = list(), success = FALSE, error = e$message))
    })
  }
  
  # Function to run multiple seeds and aggregate results
  run_multiple_seeds <- function(i, grid_row, n_seeds, default_param, fold_data_list, name_prefix, glog_level, best_scope, best_criterion, aggregation, cv_aggregation, param_yaml_path = NULL, parallel = FALSE, X = NULL, y = NULL, fold_indices = NULL, features_in_columns = TRUE, k = NULL) {
    
    if (n_seeds == 1) {
      # No multiple seeds, run single experiment without seed override
      return(run_single_experiment(i, grid_row, default_param, fold_data_list, name_prefix, glog_level, best_scope, best_criterion, aggregation, cv_aggregation, param_yaml_path, parallel, X, y, fold_indices, features_in_columns, k, seed = NULL))
    }
    
    # Run multiple times with different seeds
    seed_results <- list()
    all_metrics_list <- list()
    success_count <- 0
    
    for (seed_idx in 1:n_seeds) {
      seed_value <- 41 + seed_idx  # 42, 43, 44, ...
      
      if (seed_idx == 1) {
        message(sprintf("[%d/%d] Running experiment with %d seeds (seeds %d to %d)", 
                        i, nrow(grid_df), n_seeds, 42, 41 + n_seeds))
      } else {
        message(sprintf("  Seed %d/%d (seed=%d)", seed_idx, n_seeds, seed_value))
      }
      
      result <- run_single_experiment(i, grid_row, default_param, fold_data_list, 
                                       sprintf("%s_seed%d", name_prefix, seed_value), 
                                       glog_level, best_scope, best_criterion, 
                                       aggregation, cv_aggregation, param_yaml_path, 
                                       parallel, X, y, fold_indices, features_in_columns, 
                                       k, seed = seed_value)
      
      seed_results[[seed_idx]] <- result
      
      if (result$success) {
        success_count <- success_count + 1
        all_metrics_list[[length(all_metrics_list) + 1]] <- result$metrics
      }
    }
    
    # If all seeds failed, return failure
    if (success_count == 0) {
      errors <- sapply(seed_results, function(r) if (!r$success) r$error else NA)
      error_msg <- paste(na.omit(errors), collapse = "; ")
      return(list(params = grid_row, metrics = list(), success = FALSE, 
                  error = sprintf("All %d seeds failed. Errors: %s", n_seeds, error_msg)))
    }
    
    # Aggregate metrics across successful seeds
    # Get all metric names from successful runs
    all_metric_names <- unique(unlist(lapply(all_metrics_list, names)))
    
    # Calculate mean and seed sd for each metric
    aggregated_metrics <- list()
    
    for (metric_name in all_metric_names) {
      values <- sapply(all_metrics_list, function(m) {
        val <- m[[metric_name]]
        if (is.null(val) || is.na(val)) NA else val
      })
      
      # Remove NAs for aggregation
      values <- values[!is.na(values)]
      
      if (length(values) > 0) {
        # 1. Mean value (Mean of Means OR Mean of SDs from folds)
        aggregated_metrics[[metric_name]] <- mean(values)
        
        # 2. Standard deviation BETWEEN SEEDS (only if n_seeds > 1)
        # Only create _seed_sd for performance metrics, not for existing sd/mad/size metrics
        if (length(values) > 1) {
          # Skip if this is already an sd, mad, size, or n_features metric
          if (!grepl("_sd$", metric_name) && !grepl("_mad$", metric_name) && 
              !grepl("_size$", metric_name) && metric_name != "n_features" &&
              metric_name != "n_generations" && metric_name != "n_folds") {
            # This is a performance metric (auc, accuracy, etc.) - track seed variance
            seed_sd_name <- paste0(metric_name, "_seed_sd")
            aggregated_metrics[[seed_sd_name]] <- sd(values)
          }
        }
      }
    }
    
    # Add metadata about seeds
    aggregated_metrics$n_seeds_success <- success_count
    aggregated_metrics$n_seeds_total <- n_seeds
    
    return(list(params = grid_row, metrics = aggregated_metrics, success = TRUE))
  }
  
  # Execute grid search (parallel or sequential)
  if (parallel) {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is required for parallel execution. Install it with install.packages('doParallel')")
    }
    if (is.null(n_cores)) {
      n_cores <- max(1, floor(parallel::detectCores() / 2))
    }
    message(sprintf("Running grid search in parallel with %d cores", n_cores))
    cl <- parallel::makeCluster(n_cores, outfile = "")
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    results_list <- foreach::foreach(
      i = 1:n_combinations,
      .packages = c("gpredomicsR"),
      .export = c("extract_metrics", "aggregate_population_metrics", "aggregate_cv_metrics", "param_yaml_path", "grid_df", "name_prefix", "glog_level", "best_scope", "best_criterion", "aggregation", "cv_aggregation", "X", "y", "fold_indices", "features_in_columns", "k", "n_seeds", "run_single_experiment", "run_multiple_seeds", "default_param", "log_file", "n_combinations", "n_threads_per_model")
    ) %dopar% {
      
      # --- 1. Setup output capture (sink) ---
      worker_id <- Sys.getpid()
      
      # Create unique temporary file for this worker
      worker_log_temp <- tempfile(pattern = paste0("worker_", worker_id, "_log_"))
      
      # Open connection to temporary file
      log_con <- file(worker_log_temp, open = "wt")
      
      # Redirect all output (stdout and stderr/warnings) to this file
      sink(log_con, type = "output")
      sink(log_con, type = "message")
      
      # --- 2. Begin execution (all output is now captured) ---
      
      # A. Configure threads
      if (!is.null(n_threads_per_model) && n_threads_per_model > 0) {
        options(gpredomics.threads.number = n_threads_per_model)
      }
      
      # B. Log start message
      current_params <- grid_df[i, , drop = FALSE]
      param_str <- paste(names(current_params), current_params, sep = "=", collapse = ", ")
      cat(sprintf("START Config %d/%d | Params: {%s}\n", i, n_combinations, param_str))
      
      # C. Handle beam search detection
      algo_name <- ""
      if ("algo" %in% names(grid_df)) {
        algo_name <- as.character(grid_df[i, "algo"])
      } else {
        tryCatch({
          param_temp <- Param$load(param_yaml_path)
          val <- param_temp[["algo"]]
          if (!is.null(val)) algo_name <- val
        }, error = function(e) {})
      }
      
      current_n_seeds <- n_seeds
      if (!is.null(algo_name) && nchar(algo_name) > 0 && grepl("beam", tolower(algo_name)) && n_seeds > 1) {
        cat("Beam search detected: forcing n_seeds=1\n")
        current_n_seeds <- 1
      }
      
      # D. Execute computation
      # Note: no log_file_path here since sink() captures everything
      res <- run_multiple_seeds(i, grid_df[i, , drop = FALSE], current_n_seeds, NULL, NULL, name_prefix, glog_level, best_scope, best_criterion, aggregation, cv_aggregation, param_yaml_path = param_yaml_path, parallel = TRUE, X = X, y = y, fold_indices = fold_indices, features_in_columns = features_in_columns, k = k)
      
      # E. Log completion message
      status <- if (res$success) "SUCCESS" else "FAILURE"
      cat(sprintf("FINISH Config %d/%d | Status: %s\n", i, n_combinations, status))
      
      # --- 3. Stop sink and flush output ---
      
      # Stop output redirection
      sink(type = "message")
      sink(type = "output")
      close(log_con)
      
      # Read temporary file and write to main log
      if (file.exists(worker_log_temp)) {
        # Read captured output
        captured_output <- readLines(worker_log_temp, warn = FALSE)
        unlink(worker_log_temp) # Delete temp file
        
        # Add [Worker ID] prefix to each line
        timestamp <- format(Sys.time(), "%H:%M:%S")
        prefix <- sprintf("[%s] [Worker %s] ", timestamp, worker_id)
        formatted_output <- paste0(prefix, captured_output)
        
        # Atomic write (entire block) to main log file
        if (!is.null(log_file) && nchar(log_file) > 0) {
          try({
            cat(paste(formatted_output, collapse = "\n"), "\n", file = log_file, append = TRUE)
            # Add visual separator
            cat(paste0(paste(rep("-", 80), collapse = ""), "\n"), file = log_file, append = TRUE)
          }, silent = TRUE)
        }
        
        # Write to console (for live monitoring if outfile="")
        cat(paste(formatted_output, collapse = "\n"), "\n")
      }
      
      return(res)
    }
  } else {
    # Sequential execution
    # Apply option globally before starting loop
    if (!is.null(n_threads_per_model) && n_threads_per_model > 0) {
      options(gpredomics.threads.number = n_threads_per_model)
    }
    
    results_list <- lapply(1:n_combinations, function(i) {
      # --- A. Prepare start message ---
      timestamp <- format(Sys.time(), "%H:%M:%S")
      
      # Format parameters for display
      current_params <- grid_df[i, , drop = FALSE]
      param_str <- paste(names(current_params), current_params, sep = "=", collapse = ", ")
      
      start_msg <- sprintf("[%s] START Config %d/%d | Params: {%s}", 
                           timestamp, i, n_combinations, param_str)
      
      # --- B. Write logs (console + file) ---
      cat(paste0(start_msg, "\n"))
      if (!is.null(log_file) && nchar(log_file) > 0) {
        try(cat(paste0(start_msg, "\n"), file = log_file, append = TRUE), silent = TRUE)
      }
      
      # --- C. Detect beam search ---
      # Detect if the algorithm is deterministic (Beam Search)
      is_deterministic <- FALSE
      algo_name <- ""
      
      # 1. Check if 'algo' parameter is in the grid
      if ("algo" %in% names(grid_df)) {
        algo_name <- as.character(grid_df[i, "algo"])
      } else {
        # 2. Otherwise, check in default parameters
        algo_name <- tryCatch({
          default_param[["algo"]]
        }, error = function(e) "")
      }
      
      # 3. Detect if it's Beam Search (deterministic)
      if (!is.null(algo_name) && nchar(algo_name) > 0 && grepl("beam", tolower(algo_name))) {
        is_deterministic <- TRUE
      }
      
      # 4. Adjust n_seeds dynamically for this iteration
      current_n_seeds <- n_seeds
      if (is_deterministic && n_seeds > 1) {
        current_n_seeds <- 1
        # Optional: log message for the user
        # message(sprintf("Configuration %d uses Beam Search (deterministic): forcing n_seeds=1", i))
      }
      
      # Execute the experiment
      res <- run_multiple_seeds(i, grid_df[i, , drop = FALSE], current_n_seeds, default_param, fold_data_list, name_prefix, glog_level, best_scope, best_criterion, aggregation, cv_aggregation, param_yaml_path = NULL, parallel = FALSE)
      
      # --- D. Write completion message ---
      end_timestamp <- format(Sys.time(), "%H:%M:%S")
      status <- if (res$success) "SUCCESS" else "FAILURE"
      end_msg <- sprintf("[%s] FINISH Config %d/%d | Status: %s", 
                         end_timestamp, i, n_combinations, status)
      
      cat(paste0(end_msg, "\n"))
      if (!is.null(log_file) && nchar(log_file) > 0) {
        try(cat(paste0(end_msg, "\n"), file = log_file, append = TRUE), silent = TRUE)
      }
      
      return(res)
    })
  }
  
  # Aggregate results
  message("Aggregating results...")
  
  all_metrics <- list()
  success_count <- 0
  
  # First, collect all unique column names from all runs
  all_cols <- unique(unlist(lapply(results_list, function(res) {
    if (res$success) {
      c(names(res$params), names(res$metrics), "success")
    } else {
      c(names(res$params), "success", "error")
    }
  })))
  
  for (i in 1:length(results_list)) {
    res <- results_list[[i]]
    
    # Create a list with all columns initialized to NA
    combined <- setNames(as.list(rep(NA, length(all_cols))), all_cols)
    
    if (res$success) {
      success_count <- success_count + 1
      # Fill in parameters and metrics
      params_list <- as.list(res$params)
      for (name in names(params_list)) {
        combined[[name]] <- params_list[[name]]
      }
      for (name in names(res$metrics)) {
        combined[[name]] <- res$metrics[[name]]
      }
      combined$success <- TRUE
    } else {
      # Fill in only parameters for failed runs
      params_list <- as.list(res$params)
      for (name in names(params_list)) {
        combined[[name]] <- params_list[[name]]
      }
      combined$success <- FALSE
      combined$error <- res$error
    }
    
    all_metrics[[i]] <- combined
  }
  
  # Convert to data frame - create manually to avoid rbind issues
  results_df <- data.frame(matrix(ncol = length(all_cols), nrow = length(all_metrics)))
  colnames(results_df) <- all_cols
  
  for (i in seq_along(all_metrics)) {
    for (col in all_cols) {
      val <- all_metrics[[i]][[col]]
      results_df[i, col] <- if (is.null(val)) NA else val
    }
  }
  
  message(sprintf("Grid search completed: %d/%d runs successful", success_count, n_combinations))

  # Remove train_* columns if get_train is FALSE
  if (!get_train) {
    train_cols <- grep("^train_", names(results_df), value = TRUE)
    results_df <- results_df[, setdiff(names(results_df), train_cols), drop = FALSE]
  }

  # Sort by metric
  # Only allow metrics from compute_metrics (with train_ or valid_ prefix)
  # Valid metrics: auc, accuracy, sensitivity, specificity, rejection_rate, mcc, npv, ppv, f1_score, g_mean
  valid_compute_metrics <- c("auc", "accuracy", "sensitivity", "specificity", 
                             "rejection_rate", "mcc", "npv", "ppv", 
                             "f1_score", "g_mean")
  
  # Check if metric starts with "train_" or "valid_" (user provided full name)
  if (grepl("^(train_|valid_)", metric)) {
    metric_col <- metric
    # Validate that the base metric is from compute_metrics
    base_metric <- sub("^(train_|valid_)", "", metric)
    if (!base_metric %in% valid_compute_metrics) {
      stop(sprintf("Invalid metric '%s'. Only metrics from compute_metrics are allowed: %s", 
                   metric, paste(valid_compute_metrics, collapse = ", ")))
    }
  } else if (metric %in% valid_compute_metrics) {
    # Add "valid_" prefix for standard metrics
    metric_col <- paste0("valid_", metric)
  } else {
    stop(sprintf("Invalid metric '%s'. Only metrics from compute_metrics are allowed: %s (with optional 'train_' or 'valid_' prefix)", 
                 metric, paste(valid_compute_metrics, collapse = ", ")))
  }
  
  if (metric_col %in% names(results_df)) {
    # All compute_metrics are sorted in descending order (higher is better)
    results_df <- results_df[order(results_df[[metric_col]], decreasing = TRUE, na.last = TRUE), ]
    rownames(results_df) <- NULL
    message(sprintf("Results sorted by '%s' (decreasing)", metric_col))
  } else {
    warning(sprintf("Metric '%s' (mapped to '%s') not found in results, returning unsorted", metric, metric_col))
  }

  # Train final model with best parameters on full training data
  best_experiment <- NULL
  best_params <- NULL
  
  if (success_count > 0 && nrow(results_df) > 0) {
    message("\nTraining final model with best parameters on full training data...")
    
    # Get best parameter combination (first row after sorting)
    best_row <- results_df[1, ]
    best_params <- as.list(best_row)
    
    # Get parameter names (exclude metrics and metadata)
    param_names <- names(best_row)[!grepl("^(train_|valid_|n_features|n_generations|n_folds|overfit|jury_size|fbm_size|top_size|topn_size|experiment_name|success|error|n_seeds_)", names(best_row)) & !grepl("_sd$", names(best_row))]
    
    # Create a new Param object with best parameters
    final_param <- default_param$clone()
    
    # Enable jury/voting if best_scope is jury
    if (best_scope == "jury") {
      final_param$set_bool("vote", TRUE)
    }
    
    # Apply best parameters
    for (param_name in param_names) {
      param_value <- best_row[[param_name]]
      if (is.factor(param_value)) {
        param_value <- as.character(param_value)
      }
      if (is.logical(param_value)) {
        final_param$set_bool(param_name, param_value)
      } else if (is.character(param_value)) {
        final_param$set_string(param_name, param_value)
      } else if (is.numeric(param_value)) {
        final_param$set(param_name, param_value)
      }
    }
    
    # Create full training data
    if (features_in_columns) {
      if (!is.null(rownames(X))) {
        names(y) <- rownames(X)
      }
    } else {
      if (!is.null(colnames(X))) {
        names(y) <- colnames(X)
      }
    }
    
    full_train_data <- as.gpredomics.data(X, y, features_in_columns)
    
    # Train final model
    tryCatch({
      message("Training final model on full training data...")
      running_flag <- RunningFlag$new()
      best_experiment <- fit_on(full_train_data, final_param, running_flag)
      message("Final model training completed successfully")
    }, error = function(e) {
      warning(sprintf("Failed to train final model: %s", e$message))
      best_experiment <- NULL
    })
  }

  # Calculate total execution time
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Log completion
  if (!is.null(log_file) && nchar(log_file) > 0) {
    cat(sprintf("\n=== Grid Search Completed: %s ===\n", end_time), file = log_file, append = TRUE)
    cat(sprintf("Total execution time: %.2f seconds (%.2f minutes)\n", execution_time, execution_time / 60), file = log_file, append = TRUE)
    cat(sprintf("Successful configurations: %d/%d\n", success_count, n_combinations), file = log_file, append = TRUE)
    cat(sprintf("Best metric (%s): %.4f\n", metric_col, best_params[[metric_col]]), file = log_file, append = TRUE)
  }

  # Return as a list with results and best experiment
  result <- list(
    results = results_df,
    best = best_experiment,
    best_params = best_params,
    metric = metric_col,
    n_success = success_count,
    n_total = n_combinations,
    n_total_grid = n_combinations_total,
    execution_time = execution_time,
    random_search = random_search_n > 0,
    n_seeds = n_seeds
  )
  
  class(result) <- c("gpredomics_grid_search", "list")
  
  return(result)
}

#' Print method for grid search results
#' @param x A gpredomics_grid_search object
#' @param ... Additional arguments (unused)
#' @export
print.gpredomics_grid_search <- function(x, ...) {
  cat("Grid Search Results\n")
  cat("===================\n\n")
  
  if (x$random_search) {
    cat(sprintf("Search type: Random search (%d/%d combinations tested)\n", 
                x$n_total, x$n_total_grid))
  } else {
    cat(sprintf("Search type: Exhaustive grid search (%d combinations)\n", x$n_total))
  }
  
  if (!is.null(x$n_seeds) && x$n_seeds > 1) {
    cat(sprintf("Seeds per configuration: %d (seeds 42 to %d)\n", x$n_seeds, 41 + x$n_seeds))
  }
  
  cat(sprintf("Total runs: %d successful, %d failed\n", 
              x$n_success, x$n_total - x$n_success))
  cat(sprintf("Execution time: %.2f seconds (%.2f min)\n", 
              x$execution_time, x$execution_time / 60))
  cat(sprintf("Sorted by: %s\n\n", x$metric))
  
  if (length(x$best_params) > 0) {
    cat("Best configuration:\n")
    cat("-------------------\n")
    
    # Print parameters (excluding metrics)
    param_names <- names(x$best_params)[!grepl("^(train_|valid_|n_features|n_generations|n_folds|overfit|jury_size|fbm_size|top_size|topn_size|experiment_name|success|error|n_seeds_)", names(x$best_params)) & !grepl("_sd$", names(x$best_params))]
    for (pname in param_names) {
      cat(sprintf("  %s: %s\n", pname, x$best_params[[pname]]))
    }
    
    cat("\nBest metrics:\n")
    cat("-------------\n")
    metric_names <- names(x$best_params)[grepl("^(valid_|train_|overfit|n_features|n_seeds_)", names(x$best_params)) | grepl("_sd$", names(x$best_params))]
    for (mname in metric_names) {
      val <- x$best_params[[mname]]
      if (is.numeric(val)) {
        cat(sprintf("  %s: %.4f\n", mname, val))
      } else {
        cat(sprintf("  %s: %s\n", mname, val))
      }
    }
    
    cat("\nFinal model:\n")
    cat("------------\n")
    if (!is.null(x$best)) {
      cat("  Status: Trained successfully on full training data\n")
      cat("  Access with: $best\n")
    } else {
      cat("  Status: Training failed\n")
    }
  } else {
    cat("No successful runs\n")
  }
  
  cat(sprintf("\nFull results table: %d rows × %d columns\n", 
              nrow(x$results), ncol(x$results)))
  cat("Access with: $results\n")
  
  invisible(x)
}


#' Extract Metrics from Gpredomics Experiment
#'
#' Internal function to extract key metrics from a gpredomics experiment object
#' using predictions from models according to the specified scope.
#'
#' @param exp An Experiment object returned by fit_on()
#' @param best_scope Method to compute metrics: "best_ind", "fbm", "pct", "topn", or "jury"
#' @param best_criterion Criterion value (alpha for fbm, percentage for pct, count for topn)
#' @param aggregation Aggregation method: "mean" or "median"
#' @param train_gdata Training data as Gpredomics Data object (for train metrics)
#' @param valid_data_obj Optional validation data as Gpredomics Data object (if NULL, no validation metrics computed)
#'
#' @return A named list of metrics with train_ and valid_ prefixes
#' @keywords internal
extract_metrics <- function(exp, best_scope = "jury", best_criterion = NULL, aggregation = "mean", train_gdata = NULL, valid_data_obj = NULL) {
  
  metrics <- list()
  
  # Initialize metadata with default values to avoid missing fields
  metrics$jury_size <- NA
  metrics$effective_scope <- best_scope # Default assumption: operation succeeds
  
  tryCatch({
    # exp is directly the Rust Experiment object from fit_on()
    exp
    
    if (is.null(exp)) {
      warning("No Rust experiment object found")
      return(metrics)
    }
    
    # Get final population
    final_pop <- exp$get_best_population()

    if (is.null(final_pop) || final_pop$get_size() == 0) {
      warning("No final population found")
      return(metrics)
    }
    
    # Select models based on best_scope
    if (best_scope == "best_ind") {
      # Use only the best individual
      best_ind <- final_pop$get_individual(1)

      if (is.null(best_ind)) {
        warning("Best individual is NULL")
        return(metrics)
      }
      
      # Compute metrics on train data
      if (!is.null(train_gdata)) {
        train_metrics <- best_ind$compute_metrics(train_gdata)
        metrics$train_auc <- train_metrics$auc
        metrics$train_accuracy <- train_metrics$accuracy
        metrics$train_sensitivity <- train_metrics$sensitivity
        metrics$train_specificity <- train_metrics$specificity
        metrics$train_rejection_rate <- train_metrics$rejection_rate
        metrics$train_mcc <- train_metrics$mcc
        metrics$train_npv <- train_metrics$npv
        metrics$train_ppv <- train_metrics$ppv
        metrics$train_f1_score <- train_metrics$f1_score
        metrics$train_g_mean <- train_metrics$g_mean
      }
      
      # Compute metrics on valid data
      if (!is.null(valid_data_obj)) {
        valid_metrics <- best_ind$compute_metrics(valid_data_obj)
        metrics$valid_auc <- valid_metrics$auc
        metrics$valid_accuracy <- valid_metrics$accuracy
        metrics$valid_sensitivity <- valid_metrics$sensitivity
        metrics$valid_specificity <- valid_metrics$specificity
        metrics$valid_rejection_rate <- valid_metrics$rejection_rate
        metrics$valid_mcc <- valid_metrics$mcc
        metrics$valid_npv <- valid_metrics$npv
        metrics$valid_ppv <- valid_metrics$ppv
        metrics$valid_f1_score <- valid_metrics$f1_score
        metrics$valid_g_mean <- valid_metrics$g_mean
      }
      
      # Add k from best individual
      metrics$n_features <- best_ind$get()$k
      
      # Cleanup intermediate Rust objects
      rm(best_ind, final_pop)
      
      # Metadata for best_ind
      metrics$jury_size <- NA # Not applicable for best_ind
      metrics$effective_scope <- "best_ind"
      
    } else if (best_scope == "fbm") {
      # Use Family of Best Models
      if (is.null(best_criterion)) best_criterion <- 0.05
      
      fbm_pop <- final_pop$get_fbm(best_criterion)
      fbm_size <- fbm_pop$get_size()
      
      if (fbm_size == 0) {
        warning("FBM selection returned empty population, falling back to best individual")
        fallback_metrics <- extract_metrics(exp, best_scope = "best_ind", best_criterion = NULL, aggregation = aggregation, train_gdata = train_gdata, valid_data_obj = valid_data_obj)
        
        # Keep fallback metrics but override metadata
        fallback_metrics$fbm_size <- 0
        fallback_metrics$effective_scope <- "best_ind" # Fallback from fbm
        
        return(fallback_metrics)
      }
      
      # Aggregate metrics from individual models in FBM
      metrics <- aggregate_population_metrics(fbm_pop, train_gdata, valid_data_obj, aggregation)
      metrics$fbm_size <- fbm_size
      metrics$effective_scope <- "fbm"
      
      # Cleanup intermediate Rust objects
      rm(fbm_pop, final_pop)
      
    } else if (best_scope == "jury") {
      # Use jury/voting system
      jury <- exp$get_jury()
      
      # Jury failure case
      if (is.null(jury) || jury$get_population()$get_size() == 0) {
        warning("Jury not available, falling back to best individual")
        
        # Recursive call to obtain best_ind metrics
        fallback_metrics <- extract_metrics(exp, best_scope = "best_ind", best_criterion = NULL, aggregation = aggregation, train_gdata = train_gdata, valid_data_obj = valid_data_obj)
        
        # Keep fallback metrics but override metadata
        fallback_metrics$jury_size <- 0 # 0 clearly indicates no jury available
        fallback_metrics$effective_scope <- "best_ind" # Report actual scope used
        
        return(fallback_metrics)
      }
      
      # Compute metrics on train data with jury
      if (!is.null(train_gdata)) {
        train_metrics <- jury$compute_metrics(train_gdata)
        metrics$train_auc <- train_metrics$auc
        metrics$train_accuracy <- train_metrics$accuracy
        metrics$train_sensitivity <- train_metrics$sensitivity
        metrics$train_specificity <- train_metrics$specificity
        metrics$train_rejection_rate <- train_metrics$rejection_rate
        metrics$train_mcc <- train_metrics$mcc
        metrics$train_npv <- train_metrics$npv
        metrics$train_ppv <- train_metrics$ppv
        metrics$train_f1_score <- train_metrics$f1_score
        metrics$train_g_mean <- train_metrics$g_mean
      }
      
      # Compute metrics on valid data with jury
      if (!is.null(valid_data_obj)) {
        valid_metrics <- jury$compute_metrics(valid_data_obj)
        metrics$valid_auc <- valid_metrics$auc
        metrics$valid_accuracy <- valid_metrics$accuracy
        metrics$valid_sensitivity <- valid_metrics$sensitivity
        metrics$valid_specificity <- valid_metrics$specificity
        metrics$valid_rejection_rate <- valid_metrics$rejection_rate
        metrics$valid_mcc <- valid_metrics$mcc
        metrics$valid_npv <- valid_metrics$npv
        metrics$valid_ppv <- valid_metrics$ppv
        metrics$valid_f1_score <- valid_metrics$f1_score
        metrics$valid_g_mean <- valid_metrics$g_mean
      }
      
      # Jury success case
      # Add jury information
      metrics$jury_size <- jury$get_population()$get_size()
      metrics$effective_scope <- "jury" # Record that jury was actually used
      
      # Average k from jury members
      k_values <- sapply(1:jury$get_population()$get_size(), function(i) jury$get_population()$get_individual(i)$get()$k)
      metrics$n_features <- mean(k_values)
      
      # Cleanup intermediate Rust objects
      rm(jury, final_pop)
      
    } else if (best_scope == "topn") {
      # Use top N individuals from population
      if (is.null(best_criterion)) best_criterion <- 10
      n_top <- as.integer(best_criterion)
      
      # Ensure n_top doesn't exceed population size
      pop_size <- final_pop$get_size()
      if (n_top > pop_size) {
        warning(sprintf("topn (%d) exceeds population size (%d), using all individuals", n_top, pop_size))
        n_top <- pop_size
      }
      
      if (n_top == 0) {
        warning("topn is 0, falling back to best individual")
        fallback_metrics <- extract_metrics(exp, best_scope = "best_ind", best_criterion = NULL, aggregation = aggregation, train_gdata = train_gdata, valid_data_obj = valid_data_obj)
        fallback_metrics$topn_size <- 0
        fallback_metrics$effective_scope <- "best_ind"
        return(fallback_metrics)
      }
      
      # Get top N individuals using get_first_n method
      top_pop <- final_pop$get_first_n(n_top)
      top_size <- top_pop$get_size()
      
      if (top_size == 0) {
        warning("topn selection returned empty population, falling back to best individual")
        fallback_metrics <- extract_metrics(exp, best_scope = "best_ind", best_criterion = NULL, aggregation = aggregation, train_gdata = train_gdata, valid_data_obj = valid_data_obj)
        fallback_metrics$topn_size <- 0
        fallback_metrics$effective_scope <- "best_ind"
        return(fallback_metrics)
      }
      
      # Aggregate metrics from individual models in top N
      metrics <- aggregate_population_metrics(top_pop, train_gdata, valid_data_obj, aggregation)
      metrics$topn_size <- top_size
      metrics$effective_scope <- "topn"
      
      # Cleanup intermediate Rust objects
      rm(top_pop, final_pop)
      
    } else if (best_scope == "pct") {
      # Use top X% of population
      if (is.null(best_criterion)) best_criterion <- 5
      pct <- best_criterion
      
      # Get top pct% using native method
      top_pop <- final_pop$get_first_pct(pct)
      top_size <- top_pop$get_size()
      
      if (top_size == 0) {
        warning("Top pct selection returned empty population, falling back to best individual")
        fallback_metrics <- extract_metrics(exp, best_scope = "best_ind", best_criterion = NULL, aggregation = aggregation, train_gdata = train_gdata, valid_data_obj = valid_data_obj)
        fallback_metrics$top_size <- 0
        fallback_metrics$effective_scope <- "best_ind"
        return(fallback_metrics)
      }
      
      # Aggregate metrics from individual models in top pct
      metrics <- aggregate_population_metrics(top_pop, train_gdata, valid_data_obj, aggregation)
      metrics$top_size <- top_size
      metrics$effective_scope <- "pct"
      
      # Cleanup intermediate Rust objects
      rm(top_pop, final_pop)
    }
    
    # Add number of generations if available
    if (!is.null(exp$generation_number())) {
      tryCatch({
        metrics$n_generations <- exp$generation_number()
      }, error = function(e) {})
    }
    
    # Final garbage collection to release all Rust memory
    gc(verbose = FALSE)
    
  }, error = function(e) {
    warning(sprintf("Error extracting metrics: %s", e$message))
    gc(verbose = FALSE)  # Cleanup even on error
  })
  
  return(metrics)
}


#' Aggregate Metrics from Population
#'
#' Internal function to compute and aggregate metrics from individual models in a population.
#'
#' @param population Population object
#' @param train_data_obj Training data as Gpredomics Data object (or NULL)
#' @param valid_data_obj validation data as Gpredomics Data object (or NULL)
#' @param aggregation Aggregation method: "mean" or "median"
#'
#' @return Named list of aggregated metrics with train_ and valid_ prefixes
#' @keywords internal
aggregate_population_metrics <- function(population, train_data_obj, valid_data_obj, aggregation = "mean") {
  
  metrics <- list()
  pop_size <- population$get_size()
  
  if (pop_size == 0) {
    return(metrics)
  }
  
  # Choose aggregation function
  agg_fn <- if (aggregation == "median") median else mean
  
  # Collect metrics from train data
  if (!is.null(train_data_obj)) {
    train_auc <- c()
    train_accuracy <- c()
    train_sensitivity <- c()
    train_specificity <- c()
    train_rejection_rate <- c()
    train_mcc <- c()
    train_npv <- c()
    train_ppv <- c()
    train_f1_score <- c()
    train_g_mean <- c()
    
    for (i in 1:pop_size) {
      ind <- population$get_individual(i)
      
      if (is.null(ind)) {
        warning(sprintf("Individual %d is NULL in population (train)", i))
        next
      }
      
      m <- ind$compute_metrics(train_data_obj)

      train_auc <- c(train_auc, m$auc)
      train_accuracy <- c(train_accuracy, m$accuracy)
      train_sensitivity <- c(train_sensitivity, m$sensitivity)
      train_specificity <- c(train_specificity, m$specificity)
      train_rejection_rate <- c(train_rejection_rate, m$rejection_rate)
      train_mcc <- c(train_mcc, m$mcc)
      train_npv <- c(train_npv, m$npv)
      train_ppv <- c(train_ppv, m$ppv)
      train_f1_score <- c(train_f1_score, m$f1_score)
      train_g_mean <- c(train_g_mean, m$g_mean)
    }

    metrics$train_auc <- agg_fn(train_auc[!is.na(train_auc)])
    metrics$train_accuracy <- agg_fn(train_accuracy[!is.na(train_accuracy)])
    metrics$train_sensitivity <- agg_fn(train_sensitivity[!is.na(train_sensitivity)])
    metrics$train_specificity <- agg_fn(train_specificity[!is.na(train_specificity)])
    metrics$train_rejection_rate <- agg_fn(train_rejection_rate[!is.na(train_rejection_rate)])
    metrics$train_mcc <- agg_fn(train_mcc[!is.na(train_mcc)])
    metrics$train_npv <- agg_fn(train_npv[!is.na(train_npv)])
    metrics$train_ppv <- agg_fn(train_ppv[!is.na(train_ppv)])
    metrics$train_f1_score <- agg_fn(train_f1_score[!is.na(train_f1_score)])
    metrics$train_g_mean <- agg_fn(train_g_mean[!is.na(train_g_mean)])

  }
  
  # Collect metrics from valid data
  if (!is.null(valid_data_obj)) {
    valid_auc <- c()
    valid_accuracy <- c()
    valid_sensitivity <- c()
    valid_specificity <- c()
    valid_rejection_rate <- c()
    valid_mcc <- c()
    valid_npv <- c()
    valid_ppv <- c()
    valid_f1_score <- c()
    valid_g_mean <- c()
    
    for (i in 1:pop_size) {
      ind <- population$get_individual(i)
      
      if (is.null(ind)) {
        warning(sprintf("Individual %d is NULL in population (valid)", i))
        next
      }
      
      m <- ind$compute_metrics(valid_data_obj)
      
      valid_auc <- c(valid_auc, m$auc)
      valid_accuracy <- c(valid_accuracy, m$accuracy)
      valid_sensitivity <- c(valid_sensitivity, m$sensitivity)
      valid_specificity <- c(valid_specificity, m$specificity)
      valid_rejection_rate <- c(valid_rejection_rate, m$rejection_rate)
      valid_mcc <- c(valid_mcc, m$mcc)
      valid_npv <- c(valid_npv, m$npv)
      valid_ppv <- c(valid_ppv, m$ppv)
      valid_f1_score <- c(valid_f1_score, m$f1_score)
      valid_g_mean <- c(valid_g_mean, m$g_mean)
    }
    
    metrics$valid_auc <- agg_fn(valid_auc[!is.na(valid_auc)])
    metrics$valid_accuracy <- agg_fn(valid_accuracy[!is.na(valid_accuracy)])
    metrics$valid_sensitivity <- agg_fn(valid_sensitivity[!is.na(valid_sensitivity)])
    metrics$valid_specificity <- agg_fn(valid_specificity[!is.na(valid_specificity)])
    metrics$valid_rejection_rate <- agg_fn(valid_rejection_rate[!is.na(valid_rejection_rate)])
    metrics$valid_mcc <- agg_fn(valid_mcc[!is.na(valid_mcc)])
    metrics$valid_npv <- agg_fn(valid_npv[!is.na(valid_npv)])
    metrics$valid_ppv <- agg_fn(valid_ppv[!is.na(valid_ppv)])
    metrics$valid_f1_score <- agg_fn(valid_f1_score[!is.na(valid_f1_score)])
    metrics$valid_g_mean <- agg_fn(valid_g_mean[!is.na(valid_g_mean)])
  }
  
  # Add average k
  k_values <- sapply(1:population$get_size(), function(i) { return(population$get_individual(i)$get()$k) })
  metrics$n_features <- agg_fn(k_values)

  return(metrics)
}


#' Aggregate Metrics Across CV Folds
#'
#' Internal function to aggregate metrics from multiple cross-validation folds.
#'
#' @param fold_metrics_list List of metric lists, one per fold
#' @param cv_aggregation Aggregation method: "mean" or "median"
#'
#' @return Named list of aggregated metrics with standard deviations
#' @keywords internal
aggregate_cv_metrics <- function(fold_metrics_list, cv_aggregation = "mean") {
  
  if (length(fold_metrics_list) == 0) {
    return(list())
  }
  
  # Choose aggregation function
  agg_fn <- if (cv_aggregation == "median") median else mean
  
  # Get all metric names
  all_metric_names <- unique(unlist(lapply(fold_metrics_list, names)))
  
  aggregated <- list()
  
  # Define metrics to aggregate
  metric_names <- c("train_auc", "train_accuracy", "train_sensitivity", "train_specificity", 
                   "train_rejection_rate", "train_mcc", "train_npv", "train_ppv", 
                   "train_f1_score", "train_g_mean",
                   "valid_auc", "valid_accuracy", "valid_sensitivity", "valid_specificity", 
                   "valid_rejection_rate", "valid_mcc", "valid_npv", "valid_ppv", 
                   "valid_f1_score", "valid_g_mean",
                   "n_features", "age", "jury_size", "fbm_size", "top_size", "topn_size", "n_generations")
  
  # Count number of juries actually used (to detect fallbacks)
  jury_success_count <- 0
  for (fold_metrics in fold_metrics_list) {
    if ("effective_scope" %in% names(fold_metrics) && fold_metrics$effective_scope == "jury") {
      jury_success_count <- jury_success_count + 1
    } else if ("jury_size" %in% names(fold_metrics) && !is.na(fold_metrics$jury_size) && fold_metrics$jury_size > 0) {
      # Fallback if effective_scope doesn't exist but jury_size > 0
      jury_success_count <- jury_success_count + 1
    }
  }
  aggregated$jury_success_count <- jury_success_count
  
  for (metric_name in metric_names) {
    # Collect values across folds
    values <- sapply(fold_metrics_list, function(m) {
      if (metric_name %in% names(m)) m[[metric_name]] else NA
    })
    
    # Remove NAs
    values <- values[!is.na(values)]
    
    if (length(values) > 0) {
      aggregated[[metric_name]] <- agg_fn(values)
      
      if (length(values) > 1) {
        if (cv_aggregation == "median") {
          mad_name <- paste0(metric_name, "_mad")
          aggregated[[mad_name]] <- mad(values, constant = 1.4826)
        } else {
          sd_name <- paste0(metric_name, "_sd")
          aggregated[[sd_name]] <- sd(values)
        }
      }
    }
  }
  
  # Add number of folds
  aggregated$n_folds <- length(fold_metrics_list)
  
  # Add jury success rate if applicable
  if (jury_success_count > 0) {
    aggregated$jury_success_rate <- jury_success_count / length(fold_metrics_list)
  }
  
  # Calculate overfit metric (difference between train and valid AUC)
  if ("train_auc" %in% names(aggregated) && "valid_auc" %in% names(aggregated)) {
    aggregated$overfit <- aggregated$train_auc - aggregated$valid_auc
  }
  
  return(aggregated)
  
}

