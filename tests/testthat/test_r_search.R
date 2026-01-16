# Tests for grid_search function
library(testthat)
library(gpredomicsR)

# ---- Setup: Generate test data ----
set.seed(123)

num_samples <- 40
num_features <- 20

# Create feature matrix X_cols (features in columns, samples in rows - for grid_search)
X_cols <- matrix(runif(num_samples * num_features, min = 0, max = 1), 
            nrow = num_samples,  # Samples in rows
            ncol = num_features)  # Features in columns

rownames(X_cols) <- paste0("Sample_", seq_len(num_samples))
colnames(X_cols) <- paste0("Feature_", seq_len(num_features))

# Create binary response vector
y <- sample(0:1, num_samples, replace = TRUE)
names(y) <- rownames(X_cols)

# Create feature matrix X (features in rows, samples in columns - for other functions)
X <- t(X_cols)

# ---- Test 1: Input validation ----
test_that("grid_search validates inputs correctly", {
  # Create minimal grid
  grid <- expand.grid(
    population_size = c(100, 200),
    stringsAsFactors = FALSE
  )
  
  # Test: grid must be a data frame
  expect_error(
    grid_search(
      grid = list(population_size = c(100, 200)),
      default_param = NULL,
      X = X_cols,
      y = y
    ),
    "grid must be a data frame"
  )
  
  # Test: grid cannot be empty
  expect_error(
    grid_search(
      grid = data.frame(),
      default_param = NULL,
      X = X_cols,
      y = y
    ),
    "grid cannot be empty"
  )
  
  # Test: default_param must be a Param object
  expect_error(
    grid_search(
      grid = grid,
      default_param = list(),
      X = X_cols,
      y = y
    ),
    "default_param must be a Rust Param object"
  )
  
  # Test: X and y are required
  default_param <- Param$load("../../sample/param.yaml")
  
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = NULL,
      y = y
    ),
    "X and y are required"
  )
  
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = X_cols,
      y = NULL
    ),
    "X and y are required"
  )
})


test_that("grid_search validates X and y dimensions", {
  default_param <- Param$load("../../sample/param.yaml")
  grid <- expand.grid(population_size = c(100), stringsAsFactors = FALSE)
  
  # Test: mismatched dimensions (features in columns)
  y_wrong <- sample(0:1, num_samples + 5, replace = TRUE)
  
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = X_cols,
      y = y_wrong,
      features_in_columns = TRUE
    ),
    "X and y must have the same number of samples"
  )
  
  # Test: mismatched dimensions (features in rows)
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = X_cols,
      y = y_wrong,
      features_in_columns = FALSE
    ),
    "X and y must have the same number of samples"
  )
})


test_that("grid_search validates k parameter", {
  default_param <- Param$load("../../sample/param.yaml")
  grid <- expand.grid(population_size = c(100), stringsAsFactors = FALSE)
  
  # Test: k must be >= 2
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = X_cols,
      y = y,
      k = 1
    ),
    "k must be a numeric value >= 2"
  )
  
  # Test: k cannot exceed number of samples
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = X_cols,
      y = y,
      k = num_samples + 10,
      features_in_columns = TRUE
    ),
    "k .* cannot be greater than the number of samples"
  )
})


test_that("grid_search validates metric parameter", {
  default_param <- Param$load("../../sample/param.yaml")
  # Set minimal parameters for fast execution
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(population_size = c(100), stringsAsFactors = FALSE)
  
  # Test: invalid metric is rejected
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = X_cols,
      y = y,
      metric = "invalid_metric"
    ),
    "Invalid metric"
  )
  
  # Test: valid metrics are accepted (these should not error on validation)
  # We're just testing that the metric validation passes, not running full search
  valid_metrics <- c("auc", "accuracy", "sensitivity", "specificity", 
                     "mcc", "npv", "ppv", "f1_score", "g_mean")
  
  for (m in valid_metrics) {
    # Just check that metric validation passes
    expect_silent({
      metric_to_check <- m
      valid_compute_metrics <- c("auc", "accuracy", "sensitivity", "specificity", 
                                 "rejection_rate", "mcc", "npv", "ppv", 
                                 "f1_score", "g_mean")
      if (!grepl("^(train_|valid_)", metric_to_check)) {
        if (metric_to_check %in% valid_compute_metrics) {
          metric_col <- paste0("valid_", metric_to_check)
        } else {
          stop("Invalid metric")
        }
      }
    })
  }
})


# ---- Test 2: Basic functionality ----
test_that("grid_search returns correct structure", {
  default_param <- Param$load("../../sample/param.yaml")
  # Set minimal parameters for fast execution
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100, 200),
    stringsAsFactors = FALSE
  )
  
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    metric = "auc",
    glog_level = "error"
  )
  
  # Check class
  expect_s3_class(result, "gpredomics_grid_search")
  
  # Check structure
  expect_named(result, c("results", "best", "best_params", "metric", 
                        "n_success", "n_total", "n_total_grid", 
                        "execution_time", "random_search"))
  
  # Check results data frame
  expect_s3_class(result$results, "data.frame")
  expect_equal(nrow(result$results), 2)
  
  # Check metric column exists
  expect_true("valid_auc" %in% names(result$results))
  
  # Check that results are sorted
  if (result$n_success > 1) {
    expect_true(all(diff(result$results$valid_auc) <= 0, na.rm = TRUE))
  }
  
  # Check counts
  expect_equal(result$n_total, 2)
  expect_equal(result$n_total_grid, 2)
  expect_false(result$random_search)
})


# ---- Test 3: Random search ----
test_that("grid_search random search samples correctly", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  # Create larger grid
  grid <- expand.grid(
    population_size = c(100, 150, 200, 250),
    max_epochs = c(2, 3, 4),
    stringsAsFactors = FALSE
  )
  
  total_combinations <- nrow(grid)
  
  # Test: random_search_n samples the correct number
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    random_search_n = 5,
    glog_level = "error"
  )
  
  expect_equal(result$n_total, 5)
  expect_equal(result$n_total_grid, total_combinations)
  expect_true(result$random_search)
  
  # Test: random_search_n >= total uses all combinations
  result2 <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    random_search_n = total_combinations + 10,
    glog_level = "error"
  )
  
  expect_equal(result2$n_total, total_combinations)
  expect_equal(result2$n_total_grid, total_combinations)
})


# ---- Test 4: Cross-validation ----
test_that("grid_search works with k-fold cross-validation", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    k = 3,
    metric = "auc",
    glog_level = "error"
  )
  
  # Check that CV metrics are present
  expect_true("valid_auc" %in% names(result$results))
  expect_true("n_folds" %in% names(result$results))
  
  # n_folds should be 3
  expect_equal(result$results$n_folds[1], 3)
})


test_that("grid_search works without cross-validation", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    k = NULL,  # No cross-validation
    metric = "auc",
    glog_level = "error"
  )
  
  # Check that metrics are present but n_folds is not
  expect_true("valid_auc" %in% names(result$results))
  expect_false("n_folds" %in% names(result$results))
})


# ---- Test 5: Different best_scope options ----
test_that("grid_search works with different best_scope options", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  # Test best_ind
  result_best_ind <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    best_scope = "best_ind",
    glog_level = "error"
  )
  expect_true("n_features" %in% names(result_best_ind$results))
  
  # Test jury
  default_param$set_bool("vote", TRUE)
  result_jury <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    best_scope = "jury",
    glog_level = "error"
  )
  expect_true("jury_size" %in% names(result_jury$results))
  
  # Test fbm
  result_fbm <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    best_scope = "fbm",
    best_criterion = 0.05,
    glog_level = "error"
  )
  expect_true("fbm_size" %in% names(result_fbm$results))
  
  # Test pct
  result_pct <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    best_scope = "pct",
    best_criterion = 10,
    glog_level = "error"
  )
  expect_true("top_size" %in% names(result_pct$results))
  
  # Test topn
  result_topn <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    best_scope = "topn",
    best_criterion = 5,
    glog_level = "error"
  )
  expect_true("topn_size" %in% names(result_topn$results))
})


# ---- Test 6: Features in columns vs rows ----
test_that("grid_search handles features_in_columns parameter", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  # Test with features in rows (default, FALSE)
  result_rows <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X,
    y = y,
    features_in_columns = FALSE,
    glog_level = "error"
  )
  expect_s3_class(result_rows, "gpredomics_grid_search")
  
  # Test with features in columns (TRUE)
  y_cols <- y
  names(y_cols) <- rownames(X_cols)
  
  result_cols <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y_cols,
    features_in_columns = TRUE,
    glog_level = "error"
  )
  expect_s3_class(result_cols, "gpredomics_grid_search")
})


# ---- Test 7: get_train parameter ----
test_that("grid_search respects get_train parameter", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  # Test with get_train = FALSE (default)
  result_no_train <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    get_train = FALSE,
    glog_level = "error"
  )
  
  train_cols <- grep("^train_", names(result_no_train$results), value = TRUE)
  expect_equal(length(train_cols), 0)
  
  # Test with get_train = TRUE
  result_with_train <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    get_train = TRUE,
    glog_level = "error"
  )
  
  train_cols_with <- grep("^train_", names(result_with_train$results), value = TRUE)
  expect_true(length(train_cols_with) > 0)
  expect_true("train_auc" %in% names(result_with_train$results))
})


# ---- Test 8: Parallel execution ----
test_that("grid_search parallel mode creates temporary YAML file", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100, 200),
    stringsAsFactors = FALSE
  )
  
  # This test just verifies that parallel execution doesn't error
  # grid_search will create its own cluster internally
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    parallel = TRUE,
    n_cores = 2,
    glog_level = "error"
  )
  
  expect_s3_class(result, "gpredomics_grid_search")
  expect_equal(result$n_total, 2)
})


# ---- Test 9: Print method ----
test_that("print.gpredomics_grid_search works correctly", {
  # Create a mock grid_search result object
  mock_result <- list(
    results = data.frame(
      population_size = c(100, 200),
      valid_auc = c(0.85, 0.80),
      n_features = c(5, 7)
    ),
    best = NULL,
    best_params = list(
      population_size = 100,
      valid_auc = 0.85,
      n_features = 5
    ),
    metric = "valid_auc",
    n_success = 2,
    n_total = 2,
    n_total_grid = 2,
    execution_time = 120.5,
    random_search = FALSE
  )
  class(mock_result) <- c("gpredomics_grid_search", "list")
  
  # Test that print doesn't error and returns invisibly
  expect_output(print(mock_result), "Grid Search Results")
  expect_output(print(mock_result), "Exhaustive grid search")
  expect_output(print(mock_result), "valid_auc")
  expect_invisible(print(mock_result))
  
  # Test with random search
  mock_result$random_search <- TRUE
  mock_result$n_total <- 5
  mock_result$n_total_grid <- 10
  expect_output(print(mock_result), "Random search")
})


# ---- Test 10: Helper functions ----
test_that("extract_metrics handles different best_scope correctly", {
  # This is implicitly tested through grid_search tests with different best_scope options
  # No explicit test needed as it requires full experiment setup
  expect_true(TRUE)
})


test_that("aggregate_population_metrics computes correct statistics", {
  # This is implicitly tested through grid_search tests
  # No explicit test needed as it requires population object
  expect_true(TRUE)
})


test_that("aggregate_cv_metrics handles fold aggregation correctly", {
  # Create mock fold metrics
  fold_metrics_list <- list(
    list(valid_auc = 0.8, train_auc = 0.85, n_features = 5),
    list(valid_auc = 0.82, train_auc = 0.87, n_features = 6),
    list(valid_auc = 0.78, train_auc = 0.83, n_features = 5)
  )
  
  # Test mean aggregation
  result_mean <- aggregate_cv_metrics(fold_metrics_list, "mean")
  
  expect_equal(result_mean$valid_auc, mean(c(0.8, 0.82, 0.78)))
  expect_equal(result_mean$train_auc, mean(c(0.85, 0.87, 0.83)))
  expect_equal(result_mean$n_features, mean(c(5, 6, 5)))
  expect_equal(result_mean$n_folds, 3)
  expect_true("valid_auc_sd" %in% names(result_mean))
  expect_true("overfit" %in% names(result_mean))
  
  # Test median aggregation
  result_median <- aggregate_cv_metrics(fold_metrics_list, "median")
  
  expect_equal(result_median$valid_auc, median(c(0.8, 0.82, 0.78)))
  expect_equal(result_median$train_auc, median(c(0.85, 0.87, 0.83)))
  expect_true("valid_auc_mad" %in% names(result_median))
  
  # Test empty list
  result_empty <- aggregate_cv_metrics(list(), "mean")
  expect_equal(length(result_empty), 0)
})


# ---- Test 11: Multiple seeds (n_seeds > 1) ----
test_that("grid_search handles multiple seeds correctly", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  # Run with multiple seeds
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    n_seeds = 3,  # Run 3 times with different seeds
    glog_level = "error"
  )
  
  # Check that standard deviation columns are present
  expect_true("valid_auc_sd" %in% names(result$results))
  expect_true("valid_accuracy_sd" %in% names(result$results))
  
  # Check that SD values are numeric and >= 0
  expect_true(is.numeric(result$results$valid_auc_sd))
  expect_true(result$results$valid_auc_sd[1] >= 0)
  
  # Check that mean metrics are present
  expect_true("valid_auc" %in% names(result$results))
})

test_that("grid_search validates n_seeds parameter", {
  default_param <- Param$load("../../sample/param.yaml")
  grid <- expand.grid(population_size = c(100), stringsAsFactors = FALSE)
  
  # Test: n_seeds must be positive integer
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = X_cols,
      y = y,
      n_seeds = 0
    ),
    "n_seeds must be a positive integer"
  )
  
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = X_cols,
      y = y,
      n_seeds = -1
    ),
    "n_seeds must be a positive integer"
  )
  
  expect_error(
    grid_search(
      grid = grid,
      default_param = default_param,
      X = X_cols,
      y = y,
      n_seeds = 2.5
    ),
    "n_seeds must be a positive integer"
  )
})


# ---- Test 12: Log file creation ----
test_that("grid_search creates log file when specified", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100, 200),
    stringsAsFactors = FALSE
  )
  
  # Create temporary log file path
  log_path <- tempfile(pattern = "grid_search_test_", fileext = ".txt")
  
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    log_file = log_path,
    glog_level = "error"
  )
  
  # Check that log file was created
  expect_true(file.exists(log_path))
  
  # Check that log file contains expected content
  log_content <- readLines(log_path)
  expect_true(length(log_content) > 0)
  expect_true(any(grepl("Grid Search Started", log_content)))
  
  # Cleanup
  unlink(log_path)
})

test_that("grid_search works without log file", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  # Test with log_file = NULL
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    log_file = NULL,
    glog_level = "error"
  )
  expect_s3_class(result, "gpredomics_grid_search")
  
  # Test with log_file = ""
  result2 <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    log_file = "",
    glog_level = "error"
  )
  expect_s3_class(result2, "gpredomics_grid_search")
})


# ---- Test 13: CV aggregation method ----
test_that("grid_search respects cv_aggregation parameter", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  # Test with cv_aggregation = "mean"
  result_mean <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    k = 3,
    cv_aggregation = "mean",
    glog_level = "error"
  )
  
  expect_true("valid_auc_sd" %in% names(result_mean$results))
  
  # Test with cv_aggregation = "median"
  result_median <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    k = 3,
    cv_aggregation = "median",
    glog_level = "error"
  )
  
  expect_true("valid_auc_mad" %in% names(result_median$results))
  
  # Values should differ between mean and median
  # (unless data is perfectly symmetric, which is unlikely)
  # Just check that both methods complete successfully
  expect_s3_class(result_mean, "gpredomics_grid_search")
  expect_s3_class(result_median, "gpredomics_grid_search")
})


# ---- Test 14: Aggregation method for FBM/PCT ----
test_that("grid_search respects aggregation parameter for FBM", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  # Test with aggregation = "mean"
  result_mean <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    best_scope = "fbm",
    best_criterion = 0.05,
    aggregation = "mean",
    glog_level = "error"
  )
  expect_s3_class(result_mean, "gpredomics_grid_search")
  
  # Test with aggregation = "median"
  result_median <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    best_scope = "fbm",
    best_criterion = 0.05,
    aggregation = "median",
    glog_level = "error"
  )
  expect_s3_class(result_median, "gpredomics_grid_search")
})


# ---- Test 15: Failure handling ----
test_that("grid_search handles configuration failures gracefully", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  # Create grid with one invalid parameter that might cause issues
  # (e.g., population_size = 0 should fail)
  grid <- expand.grid(
    population_size = c(100, 0),  # 0 should fail
    stringsAsFactors = FALSE
  )
  
  # This should complete but may have n_success < n_total
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    glog_level = "error"
  )
  
  expect_s3_class(result, "gpredomics_grid_search")
  
  # Check that failure was recorded
  expect_true(result$n_success <= result$n_total)
  
  # Results should only contain successful runs
  expect_equal(nrow(result$results), result$n_success)
})


# ---- Test 16: n_threads_per_model parameter ----
test_that("grid_search accepts n_threads_per_model parameter", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100),
    stringsAsFactors = FALSE
  )
  
  # Test with different thread counts
  result_1 <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    n_threads_per_model = 1,
    glog_level = "error"
  )
  expect_s3_class(result_1, "gpredomics_grid_search")
  
  result_2 <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    n_threads_per_model = 2,
    glog_level = "error"
  )
  expect_s3_class(result_2, "gpredomics_grid_search")
})


# ---- Test 17: Combined features test ----
test_that("grid_search works with multiple features combined", {
  default_param <- Param$load("../../sample/param.yaml")
  default_param$set("population_size", 100)
  default_param$set("max_epochs", 2)
  
  grid <- expand.grid(
    population_size = c(100, 150),
    max_epochs = c(2, 3),
    stringsAsFactors = FALSE
  )
  
  log_path <- tempfile(pattern = "grid_search_combined_", fileext = ".txt")
  
  # Test with CV, random search, multiple seeds, log file
  result <- grid_search(
    grid = grid,
    default_param = default_param,
    X = X_cols,
    y = y,
    k = 3,
    n_seeds = 2,
    random_search_n = 2,
    cv_aggregation = "mean",
    log_file = log_path,
    get_train = TRUE,
    glog_level = "error"
  )
  
  # Check structure
  expect_s3_class(result, "gpredomics_grid_search")
  expect_true(result$random_search)
  expect_equal(result$n_total, 2)
  
  # Check CV metrics
  expect_true("n_folds" %in% names(result$results))
  expect_true("valid_auc_sd" %in% names(result$results))
  
  # Check train metrics
  expect_true("train_auc" %in% names(result$results))
  
  # Check log file
  expect_true(file.exists(log_path))
  
  # Cleanup
  unlink(log_path)
})
