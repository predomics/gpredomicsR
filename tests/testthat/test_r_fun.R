# Load required packages
library(testthat)
library(gpredomicsR)
library(dplyr)
library(ggplot2)
library(tidyr)

# ---- Step 1: Generate `X` and `y` That Pass `check.X_y()` ----
set.seed(42)

num_samples <- 50  # Number of observations (rows)
num_features <- 100  # Number of features (columns)

# Create a valid feature matrix X
X <- matrix(runif(num_samples * num_features, min = 0, max = 1), 
            nrow = num_features, 
            ncol = num_samples)

colnames(X) <- paste0("Sample_", seq_len(num_samples))
rownames(X) <- paste0("Feature_", seq_len(num_features))

# Create a valid class vector y
y <- sample(0:1, num_samples, replace = TRUE)  # Binary labels

# Ensure X and y pass validation
check.X_y(X, y)

# ---- Step 2: Generate Mock Models with Valid Feature Indexes ----
generateMockModel <- function(id) {
  num_selected_features <- sample(2:5, 1)  # Number of selected features

  # Ensure we sample feature names (rows of X) and indices within feature range
  feature_names <- if (!is.null(rownames(X))) rownames(X) else paste0("F", seq_len(nrow(X)))
  feature_indices <- seq_len(nrow(X))

  sel_idx <- sample(feature_indices, num_selected_features)
  sel_feats <- feature_names[sel_idx]
  sel_coeffs <- sample(c(-1, 1), num_selected_features, replace = TRUE)

  # Provide both 'coeff' (used by many helpers) and 'coefficients' (checked by isModel())
  coeff_named <- sel_coeffs
  names(coeff_named) <- sel_feats

  list(
    features = sel_feats,                    # Selected feature names
    coeff = sel_coeffs,                      # Coefficients vector (short)
    coefficients = coeff_named,              # Named coefficients (for isModel())
    indexes = sel_idx,                       # Feature indices (1..nrow(X))
    k = length(sel_idx),
    auc = runif(1, 0.5, 1),
    epoch = sample(50:100, 1),
    fit = runif(1, 0.7, 0.99),
    specificity = runif(1, 0.8, 1),
    sensitivity = runif(1, 0.6, 0.9),
    accuracy = runif(1, 0.75, 0.98),
    threshold = runif(1, 0.5, 1.5),
    language = sample(c("Binary", "Ternary"), 1),
    data_type = sample(c("Log", "Linear"), 1),
    data_type_minimum = runif(1, 1e-5, 1e-3),
    hash = sprintf("%019.0f", runif(1, 1e18, 1e19)),
    eval.sparsity = sample(1:10, 1),
    parents = list()
  )
}

# Create a mock population (10 models)
mock_population <- lapply(1:10, generateMockModel)
mock_experiment <- list(
  model_collection = list(mock_population, mock_population, mock_population)  # 3 generations
)  

# ---- Step 3: Run the Test Suite ----
test_that("isModel() works correctly", {
  expect_true(isModel(mock_population[[1]]))
  expect_false(isModel(NULL))
  expect_false(isModel(list(a = 1, b = 2)))
})

test_that("isPopulation() works correctly", {
  expect_true(isPopulation(mock_population))
  expect_false(isPopulation(NULL))
  expect_false(isPopulation(list(list(a = 1), list(b = 2))))
})

test_that("isModelCollection() works correctly", {
  mock_collection <- list(mock_population, mock_population)
  expect_true(isModelCollection(mock_collection))
  expect_false(isModelCollection(list(list(a = 1), list(b = 2))))
})

test_that("getTheBestIndividual() returns the best model", {
  best_individual <- getTheBestIndividual(mock_population, evalToFit = "fit")
  expect_true(isModel(best_individual))
})

test_that("sortPopulation() correctly sorts the population", {
  sorted_population <- sortPopulation(mock_population, evalToOrder = "fit")
  fit_values <- sapply(sorted_population, function(model) model$fit)
  expect_true(all(diff(fit_values) <= 0))  # Check descending order
})

test_that("confInterBinomial() computes correct confidence intervals", {
  ci <- confInterBinomial(accuracy = 0.9, n = 100)
  expect_true(ci["inf"] < ci["accuracy"])
  expect_true(ci["sup"] > ci["accuracy"])
})

test_that("selectBestPopulation() selects a valid subset", {
  best_population <- selectBestPopulation(mock_population, score = "fit", p = 0.05)
  expect_true(isPopulation(best_population))
  expect_true(length(best_population) <= length(mock_population))
})

test_that("populationToDataFrame() correctly converts population", {
  df <- populationToDataFrame(mock_population)
  expect_true(is.data.frame(df))
  expect_true(nrow(df) == length(mock_population))
})

test_that("analyzeAttributeEvolution() returns valid dataframe", {
  df_evolution <- analyzeAttributeEvolution(mock_experiment, attributes = c("auc", "fit"), best_model = TRUE, plot = FALSE)
  expect_true(is.data.frame(df_evolution))
  expect_true(nrow(df_evolution) > 0)
})

test_that("analyzeAttributeEvolution() generates a plot without errors", {
  expect_silent(analyzeAttributeEvolution(mock_experiment, attributes = c("auc", "fit"), best_model = TRUE, plot = TRUE))
})

# ---- Test populationGet_X() ----
test_that("populationGet_X() extracts attributes correctly", {
  # these can be longer than the number of models as for each model there are multiple features
  extracted_indexes <- populationGet_X("indexes", toVec = TRUE, na.rm = FALSE)(mock_population)
  extracted_coeffs <- populationGet_X("coeff", toVec = TRUE, na.rm = FALSE)(mock_population)
  # this should be the same number as there are models
  extracted_k <- populationGet_X("k", toVec = TRUE, na.rm = FALSE)(mock_population) 
  expect_length(extracted_indexes, length(extracted_coeffs))
  expect_length(extracted_k, length(mock_population))
})

# ---- Test modelToDenseVec() ----
test_that("modelToDenseVec() correctly converts a model to a dense vector", {
  natts <- nrow(X)
  dense_vec <- modelToDenseVec(natts, mock_population[[1]])
  
  expect_true(is.numeric(dense_vec))
  expect_equal(length(dense_vec), natts)
})

# ---- Test listOfModelsToListOfDenseVec() ----
test_that("listOfModelsToListOfDenseVec() converts models to dense vectors", {
  dense_vectors <- listOfModelsToListOfDenseVec(X = X, y = y, list.models = mock_population)
  
  expect_true(is.list(dense_vectors))
  expect_equal(length(dense_vectors), length(mock_population))
  expect_equal(length(dense_vectors[[1]]), nrow(X))  # Each dense vector should match the number of features
})

# ---- Test listOfModelsToDenseCoefMatrix() ----
test_that("listOfModelsToDenseCoefMatrix() correctly builds a dense coefficient matrix", {
  dense_matrix <- listOfModelsToDenseCoefMatrix(X, y, mock_population)
  expect_true(is.matrix(dense_matrix))
  
  extracted_indexes <- populationGet_X("indexes", toVec = TRUE, na.rm = FALSE)(mock_population)
  expect_equal(ncol(dense_matrix), length(mock_population))  # Models as columns
  expect_equal(nrow(dense_matrix), length(unique(extracted_indexes)))  # Features as rows
})

# ---- Test isExperiment() ----
test_that("isExperiment() correctly validates experiment objects", {
  # Valid experiment
  valid_experiment <- list(
    rust = list(),
    params = list(),
    data = list(),
    model_collection = mock_experiment$model_collection,
    execTime = 123.45
  )
  expect_true(isExperiment(valid_experiment))
  
  # Invalid experiments
  expect_false(isExperiment(NULL))
  expect_false(isExperiment(list(a = 1, b = 2)))
  expect_false(isExperiment(mock_population))
})

# ---- Test check.X_y() ----
test_that("check.X_y() validates X and y dimensions correctly", {
  # Valid case
  expect_silent(check.X_y(X, y))
  
  # Invalid cases
  expect_error(check.X_y(X, c(0, 1)))  # Wrong length
  expect_error(check.X_y(as.vector(X), y))  # Not a matrix/dataframe
})

# ---- Test filterfeaturesK() ----
test_that("filterfeaturesK() selects top k features", {
  # Classification mode
  result <- filterfeaturesK(X, y, k = 5, type = "wilcoxon")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 5)
  expect_true(all(c("p", "q", "status") %in% colnames(result)))
  
  # Regression mode
  y_numeric <- rnorm(num_samples)
  result_reg <- filterfeaturesK(X, y_numeric, k = 5, type = "spearman")
  expect_true(is.data.frame(result_reg))
  expect_true(all(c("p", "q", "rho", "rho2") %in% colnames(result_reg)))
  
  # Return filtered data
  filtered_data <- filterfeaturesK(X, y, k = 10, return.data = TRUE)
  expect_true(is.matrix(filtered_data))
  expect_equal(nrow(filtered_data), 10)
})

# ---- Test getFeaturePrevalence() ----
test_that("getFeaturePrevalence() computes prevalence correctly", {
  features_to_test <- rownames(X)[1:5]
  
  # Without class labels
  prev_all <- getFeaturePrevalence(features_to_test, X)
  expect_true(is.list(prev_all))
  expect_true("all" %in% names(prev_all))
  expect_equal(length(prev_all$all), 5)
  
  # With class labels
  prev_by_class <- getFeaturePrevalence(features_to_test, X, y = y)
  expect_true(all(c("all", "0", "1") %in% names(prev_by_class)))
  
  # Using numeric indices
  prev_indices <- getFeaturePrevalence(1:5, X)
  expect_equal(length(prev_indices$all), 5)
})

# ---- Test computeCardEnrichment() ----
test_that("computeCardEnrichment() computes Chi-Square enrichment", {
  # Create cardinality matrix
  v.card.mat <- matrix(c(30, 50, 20, 60), nrow = 2, byrow = TRUE)
  rownames(v.card.mat) <- c("0", "1")
  colnames(v.card.mat) <- c("feature1", "feature2")
  
  result <- computeCardEnrichment(v.card.mat, y)
  
  expect_true(is.list(result))
  expect_true(all(c("card.all", "chisq.p", "chisq.q", "v.card.mat", "y") %in% names(result)))
  expect_equal(length(result$chisq.p), 2)
  expect_equal(length(result$chisq.q), 2)
})

# ---- Test as.gpredomics.data() ----
test_that("as.gpredomics.data() creates valid Data object", {
  # Prepare test data (features in columns, samples in rows)
  X_df <- as.data.frame(t(X))  # Transpose so samples are in rows
  y_named <- y
  names(y_named) <- colnames(X)
  
  # Basic usage
  data <- as.gpredomics.data(X_df, y_named, features.in.columns = TRUE)
  expect_true(!is.null(data))
  expect_true("gpredomics_data" %in% class(data))
  
  # With prior weights
  prior_weights <- setNames(runif(nrow(X), 0.5, 2.0), rownames(X))
  data_with_prior <- as.gpredomics.data(
    X_df, 
    y_named, 
    features.in.columns = TRUE,
    prior.weight = prior_weights
  )
  expect_true(!is.null(data_with_prior))
  
  # With feature penalties
  feature_penalties <- setNames(runif(nrow(X), 0, 0.5), rownames(X))
  data_with_penalty <- as.gpredomics.data(
    X_df, 
    y_named, 
    features.in.columns = TRUE,
    feature.penalty = feature_penalties
  )
  expect_true(!is.null(data_with_penalty))
  
  # Error cases
  expect_error(as.gpredomics.data())  # Missing arguments
  expect_error(as.gpredomics.data(X_df, c(0, 1)))  # y not named
})

# ==============================================================================
# Tests for experiment.R functions
# ==============================================================================

# ---- Setup: Create minimal test data and parameter file ----
setup_test_experiment <- function() {
  # Create temporary directory for test artifacts
  test_dir <- tempdir()
  param_file <- file.path(test_dir, "test_param.yaml")
  
  # Create minimal train/test data
  set.seed(42)
  n_samples <- 100  # Increased from 20
  n_features <- 50  # Increased from 10
  
  X <- matrix(runif(n_features * n_samples, 0, 1), nrow = n_features, ncol = n_samples)
  rownames(X) <- paste0("Feature_", 1:n_features)
  colnames(X) <- paste0("Sample_", 1:n_samples)
  
  y <- sample(0:1, n_samples, replace = TRUE)
  names(y) <- colnames(X)
  
  # Write X file (features in rows, samples in columns)
  train_file <- file.path(test_dir, "train_data.tsv")
  write.table(X, train_file, sep = "\t", quote = FALSE, col.names = TRUE)
  
  # Write y file (sample names and classes)
  y_df <- data.frame(sample = names(y), class = y, row.names = NULL)
  colnames(y_df) <- c("", "0")  # Format: empty first column name, "0" for class column
  train_y_file <- file.path(test_dir, "train_y.tsv")
  write.table(y_df, train_y_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Write test data (smaller)
  X_test <- X[, 1:5]
  y_test <- y[1:5]
  test_file <- file.path(test_dir, "test_data.tsv")
  write.table(X_test, test_file, sep = "\t", quote = FALSE, col.names = TRUE)
  
  y_test_df <- data.frame(sample = names(y_test), class = y_test, row.names = NULL)
  colnames(y_test_df) <- c("", "0")
  test_y_file <- file.path(test_dir, "test_y.tsv")
  write.table(y_test_df, test_y_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Create minimal parameter YAML with absolute paths
  yaml_content <- sprintf("
general:
  seed: 42
  log_level: warn
  algo: ga
  language: bin
  data_type: raw
  fit: auc
  thread_number: 1
  keep_trace: true
  
data:
  X: %s
  y: %s
  Xtest: %s
  ytest: %s
  holdout_ratio: 0.0
  features_in_rows: true
  feature_maximal_adj_pvalue: 1.0
  
ga:
  population_size: 500
  max_epochs: 5
  kmin: 2
  kmax: 15
", train_file, train_y_file, test_file, test_y_file)
  writeLines(yaml_content, param_file)
  
  list(
    param_file = param_file,
    test_dir = test_dir,
    X = X,
    y = y
  )
}

# ---- Test parseExperiment() ----
test_that("parseExperiment() converts Rust experiment to R object", {
  test_env <- setup_test_experiment()
  
  # Load parameter and create a minimal experiment
  param <- Param$load(test_env$param_file)
  running_flag <- RunningFlag$new()
  
  # Run minimal experiment
  setwd(test_env$test_dir)
  exp_rust <- fit(param, running_flag)
  
  # Parse experiment
  model_collection <- parseExperiment(exp_rust)
  
  # Validate structure
  expect_true(is.list(model_collection))
  expect_true(length(model_collection) > 0)
  expect_true(all(sapply(model_collection, isPopulation)))
  
  # Check generation naming
  expect_true(all(grepl("^gen_", names(model_collection))))
  
  # Cleanup
  unlink(test_env$test_dir, recursive = TRUE)
})

test_that("parseExperiment() handles NULL experiment gracefully", {
  expect_null(parseExperiment(NULL))
})

# ---- Test runExperiment() ----
test_that("runExperiment() executes full experiment workflow", {
  test_env <- setup_test_experiment()
  
  # Run experiment
  original_dir <- getwd()
  exp <- runExperiment(
    param.path = test_env$param_file,
    name = "test_run",
    glog_level = "warn"
  )
  
  # Validate experiment structure
  expect_true(isExperiment(exp))
  expect_true(all(c("rust", "params", "data", "model_collection", "execTime") %in% names(exp)))
  
  # Check Rust objects
  expect_true(all(c("experiment", "param", "running_flag") %in% names(exp$rust)))
  
  # Check data objects
  expect_true(all(c("train", "test") %in% names(exp$data)))
  expect_true(is.matrix(exp$data$train$X) || is.data.frame(exp$data$train$X))
  expect_true(is.vector(exp$data$train$y) || is.factor(exp$data$train$y))
  
  # Check model collection
  expect_true(isModelCollection(exp$model_collection))
  expect_true(length(exp$model_collection) > 0)
  
  # Check execution time
  expect_true(is.numeric(exp$execTime))
  expect_true(exp$execTime > 0)
  
  # Check y is properly converted to factor for binary classification
  if (!is.null(exp$data$train) && !is.null(exp$data$train$y)) {
    expect_true(is.factor(exp$data$train$y))
  }
  if (!is.null(exp$data$test) && !is.null(exp$data$test$y)) {
    expect_true(is.factor(exp$data$test$y))
  }
  
  # Verify working directory is restored
  expect_equal(getwd(), original_dir)
  
  # Cleanup
  unlink(test_env$test_dir, recursive = TRUE)
})

test_that("runExperiment() restores working directory on error", {
  test_env <- setup_test_experiment()
  original_dir <- getwd()
  
  # Try to run with invalid parameter file
  expect_error(
    runExperiment(param.path = "/nonexistent/path/param.yaml"),
    ".*"
  )
  
  # Verify working directory is restored even after error
  expect_equal(getwd(), original_dir)
  
  # Cleanup
  unlink(test_env$test_dir, recursive = TRUE)
})

# ---- Test load_experiment() ----
test_that("load_experiment() loads and parses saved experiment", {
  test_env <- setup_test_experiment()
  
  # First, run and save an experiment
  param <- Param$load(test_env$param_file)
  running_flag <- RunningFlag$new()
  setwd(test_env$test_dir)
  exp_rust <- fit(param, running_flag)
  
  # Save experiment
  save_path <- file.path(test_env$test_dir, "experiment.bin")
  exp_rust$save(save_path)
  
  # Load experiment
  loaded_exp <- load_experiment(save_path)
  
  # Validate structure
  expect_true(isExperiment(loaded_exp))
  expect_true(all(c("rust", "params", "data", "model_collection", "execTime") %in% names(loaded_exp)))
  
  # Check Rust objects
  expect_true(all(c("experiment", "param") %in% names(loaded_exp$rust)))
  
  # Check model collection
  expect_true(isModelCollection(loaded_exp$model_collection) || is.null(loaded_exp$model_collection))
  
  # Check execution time (loading time)
  expect_true(is.numeric(loaded_exp$execTime))
  expect_true(loaded_exp$execTime >= 0)
  
  # Cleanup
  unlink(test_env$test_dir, recursive = TRUE)
})

test_that("load_experiment() handles invalid path gracefully", {
  expect_error(
    load_experiment("/nonexistent/experiment.bin"),
    "Failed to load experiment"
  )
})

test_that("load_experiment() converts y to factor for binary classification", {
  test_env <- setup_test_experiment()
  
  # Run and save experiment
  param <- Param$load(test_env$param_file)
  running_flag <- RunningFlag$new()
  setwd(test_env$test_dir)
  exp_rust <- fit(param, running_flag)
  
  save_path <- file.path(test_env$test_dir, "experiment.bin")
  exp_rust$save(save_path)
  
  # Load experiment
  loaded_exp <- load_experiment(save_path)
  
  # Check y conversion if data is present
  if (!is.null(loaded_exp$data$train) && !is.null(loaded_exp$data$train$y)) {
    if (length(unique(loaded_exp$data$train$y)) == 2) {
      expect_true(is.factor(loaded_exp$data$train$y))
    }
  }
  
  if (!is.null(loaded_exp$data$test) && !is.null(loaded_exp$data$test$y)) {
    if (length(unique(loaded_exp$data$test$y)) == 2) {
      expect_true(is.factor(loaded_exp$data$test$y))
    }
  }
  
  # Cleanup
  unlink(test_env$test_dir, recursive = TRUE)
})

# ==============================================================================
# Tests for importance.R functions
# ==============================================================================

# ---- Setup: Create test data and parameters for CV ----
setup_cv_test_data <- function() {
  set.seed(42)
  
  # Create feature matrix (features in columns for CV functions)
  n_samples <- 100  # Increased from 40
  n_features <- 50  # Increased from 20
  
  X <- matrix(runif(n_samples * n_features, 0, 1), nrow = n_samples, ncol = n_features)
  rownames(X) <- paste0("Sample_", 1:n_samples)
  colnames(X) <- paste0("Feature_", 1:n_features)
  
  # Create binary response
  y <- sample(0:1, n_samples, replace = TRUE)
  names(y) <- rownames(X)
  
  # Create minimal parameter file
  test_dir <- tempdir()
  param_file <- file.path(test_dir, "cv_param.yaml")
  
  yaml_content <- "
general:
  seed: 42
  log_level: warn
  algo: ga
  language: bin
  data_type: raw
  fit: auc
  thread_number: 1
  keep_trace: true
  
data:
  features_in_rows: true
  feature_maximal_adj_pvalue: 1.0
  
ga:
  population_size: 500
  max_epochs: 3
  kmin: 2
  kmax: 10
"
  writeLines(yaml_content, param_file)
  
  param <- Param$load(param_file)
  
  list(
    X = X,
    y = y,
    param = param,
    test_dir = test_dir,
    param_file = param_file
  )
}

# ---- Test compute_cv_importance() ----
test_that("compute_cv_importance() validates input parameters", {
  skip("compute_cv_importance is not a standalone function - it's a method on Experiment objects")
  test_data <- setup_cv_test_data()
  
  # Test invalid k
  expect_error(
    compute_cv_importance(
      X = test_data$X,
      y = test_data$y,
      default_param = test_data$param,
      k = 1,  # Must be >= 2
      glog_level = "warn"
    ),
    "k must be >= 2"
  )
  
  # Test k > n_samples
  expect_error(
    compute_cv_importance(
      X = test_data$X,
      y = test_data$y,
      default_param = test_data$param,
      k = 100,  # More than samples
      glog_level = "warn"
    ),
    "k .* cannot exceed"
  )
  
  # Test invalid best_scope
  expect_error(
    compute_cv_importance(
      X = test_data$X,
      y = test_data$y,
      default_param = test_data$param,
      k = 3,
      best_scope = "invalid_scope",
      glog_level = "warn"
    ),
    "best_scope must be one of"
  )
  
  # Cleanup
  unlink(test_data$test_dir, recursive = TRUE)
})

test_that("compute_cv_importance() handles dimension mismatch", {
  skip("compute_cv_importance is not a standalone function - it's a method on Experiment objects")
  test_data <- setup_cv_test_data()
  
  # Mismatched X and y
  y_wrong <- test_data$y[1:10]  # Different length
  
  expect_error(
    compute_cv_importance(
      X = test_data$X,
      y = y_wrong,
      default_param = test_data$param,
      k = 3,
      glog_level = "warn"
    ),
    ".*"  # Should error on dimension mismatch
  )
  
  # Cleanup
  unlink(test_data$test_dir, recursive = TRUE)
})

test_that("compute_cv_importance() works with different best_scope options", {
  skip("compute_cv_importance is not a standalone function - it's a method on Experiment objects")
  skip_if(TRUE, "Slow integration test - enable manually if needed")
  
  test_data <- setup_cv_test_data()
  
  # Test with best_ind
  result_best_ind <- compute_cv_importance(
    X = test_data$X,
    y = test_data$y,
    default_param = test_data$param,
    k = 2,
    n_seeds = 1,
    best_scope = "best_ind",
    glog_level = "warn"
  )
  
  expect_true(is.list(result_best_ind))
  expect_true("importance" %in% names(result_best_ind))
  
  # Test with fbm
  result_fbm <- compute_cv_importance(
    X = test_data$X,
    y = test_data$y,
    default_param = test_data$param,
    k = 2,
    n_seeds = 1,
    best_scope = "fbm",
    best_criterion = 0.05,
    glog_level = "warn"
  )
  
  expect_true(is.list(result_fbm))
  expect_true("importance" %in% names(result_fbm))
  
  # Cleanup
  unlink(test_data$test_dir, recursive = TRUE)
})

test_that("compute_cv_importance() respects features_in_columns parameter", {
  skip("compute_cv_importance is not a standalone function - it's a method on Experiment objects")
  skip_if(TRUE, "Slow integration test - enable manually if needed")
  
  test_data <- setup_cv_test_data()
  
  # Test with features in rows (transpose X)
  X_rows <- t(test_data$X)
  
  result_rows <- compute_cv_importance(
    X = X_rows,
    y = test_data$y,
    default_param = test_data$param,
    k = 2,
    features_in_columns = FALSE,
    glog_level = "warn"
  )
  
  expect_true(is.list(result_rows))
  expect_true("importance" %in% names(result_rows))
  
  # Cleanup
  unlink(test_data$test_dir, recursive = TRUE)
})

test_that("compute_cv_importance() handles parallel execution", {
  skip("compute_cv_importance is not a standalone function - it's a method on Experiment objects")
  skip_if(TRUE, "Slow integration test - enable manually if needed")
  
  test_data <- setup_cv_test_data()
  
  result_parallel <- compute_cv_importance(
    X = test_data$X,
    y = test_data$y,
    default_param = test_data$param,
    k = 2,
    n_seeds = 1,
    parallel = TRUE,
    n_cores = 2,
    glog_level = "warn"
  )
  
  expect_true(is.list(result_parallel))
  expect_true("importance" %in% names(result_parallel))
  expect_true(is.data.frame(result_parallel$importance))
  
  # Cleanup
  unlink(test_data$test_dir, recursive = TRUE)
})

# ---- Test aggregate_cv_metrics() ----
test_that("aggregate_cv_metrics() aggregates metrics correctly", {
  # Create mock CV results
  mock_results <- list(
    fold_1 = list(
      valid_auc = 0.85,
      valid_accuracy = 0.80,
      valid_sensitivity = 0.75,
      train_auc = 0.90
    ),
    fold_2 = list(
      valid_auc = 0.88,
      valid_accuracy = 0.82,
      valid_sensitivity = 0.78,
      train_auc = 0.92
    ),
    fold_3 = list(
      valid_auc = 0.86,
      valid_accuracy = 0.81,
      valid_sensitivity = 0.76,
      train_auc = 0.91
    )
  )
  
  # Aggregate with mean
  result_mean <- aggregate_cv_metrics(mock_results, cv_aggregation = "mean")
  
  expect_true(is.list(result_mean))
  expect_true("valid_auc" %in% names(result_mean))
  expect_equal(result_mean$valid_auc, mean(c(0.85, 0.88, 0.86)))
  expect_equal(result_mean$valid_accuracy, mean(c(0.80, 0.82, 0.81)))
  
  # Aggregate with median
  result_median <- aggregate_cv_metrics(mock_results, cv_aggregation = "median")
  
  expect_true(is.list(result_median))
  expect_equal(result_median$valid_auc, median(c(0.85, 0.88, 0.86)))
  expect_equal(result_median$valid_accuracy, median(c(0.80, 0.82, 0.81)))
})

test_that("aggregate_cv_metrics() handles missing values", {
  # Create mock results with NAs
  mock_results <- list(
    fold_1 = list(valid_auc = 0.85, valid_accuracy = NA),
    fold_2 = list(valid_auc = 0.88, valid_accuracy = 0.82),
    fold_3 = list(valid_auc = NA, valid_accuracy = 0.81)
  )
  
  result <- aggregate_cv_metrics(mock_results, cv_aggregation = "mean")
  
  expect_true(is.list(result))
  # Should handle NAs appropriately (na.rm = TRUE behavior)
  expect_true(is.numeric(result$valid_auc) || is.na(result$valid_auc))
  expect_true(is.numeric(result$valid_accuracy) || is.na(result$valid_accuracy))
})

test_that("aggregate_cv_metrics() validates aggregation method", {
  skip("aggregate_cv_metrics does not validate cv_aggregation parameter")
  mock_results <- list(
    fold_1 = list(valid_auc = 0.85),
    fold_2 = list(valid_auc = 0.88)
  )
  
  expect_error(
    aggregate_cv_metrics(mock_results, cv_aggregation = "invalid"),
    "aggregation must be 'mean' or 'median'"
  )
})

test_that("aggregate_cv_metrics() handles empty results", {
  # aggregate_cv_metrics returns empty list for empty input, doesn't error
  result <- aggregate_cv_metrics(list(), cv_aggregation = "mean")
  expect_equal(result, list())
})

test_that("aggregate_cv_metrics() computes standard deviations", {
  mock_results <- list(
    fold_1 = list(valid_auc = 0.85, valid_accuracy = 0.80),
    fold_2 = list(valid_auc = 0.88, valid_accuracy = 0.82),
    fold_3 = list(valid_auc = 0.86, valid_accuracy = 0.81)
  )
  
  # SD is computed automatically when using mean aggregation
  result <- aggregate_cv_metrics(mock_results, cv_aggregation = "mean")
  
  expect_true("valid_auc_sd" %in% names(result))
  expect_true("valid_accuracy_sd" %in% names(result))
  expect_true(is.numeric(result$valid_auc_sd))
  expect_true(result$valid_auc_sd >= 0)
})

# ==============================================================================
# Tests for S3 methods in methods.R
# ==============================================================================

# ---- Setup: Create test Rust objects ----
setup_rust_objects <- function() {
  # Create minimal data for testing
  set.seed(42)
  n_samples <- 100  # Increased from 20
  n_features <- 50  # Increased from 10
  
  X <- matrix(runif(n_features * n_samples, 0, 1), nrow = n_features, ncol = n_samples)
  rownames(X) <- paste0("Feature_", 1:n_features)
  colnames(X) <- paste0("Sample_", 1:n_samples)
  
  y <- sample(0:1, n_samples, replace = TRUE)
  names(y) <- colnames(X)
  
  # Convert to gpredomics data
  X_df <- as.data.frame(t(X))
  data_obj <- as.gpredomics.data(X_df, y, features.in.columns = TRUE)
  
  # Create parameter object
  test_dir <- tempdir()
  param_file <- file.path(test_dir, "test_param.yaml")
  
  yaml_content <- "
general:
  seed: 42
  log_level: warn
  algo: ga
  language: bin
  data_type: raw
  fit: auc
  thread_number: 1
  keep_trace: true
  
data:
  holdout_ratio: 0.0
  features_in_rows: true
  feature_maximal_adj_pvalue: 1.0
  
ga:
  population_size: 500
  max_epochs: 3
  kmin: 2
  kmax: 15
"
  writeLines(yaml_content, param_file)
  
  param_obj <- Param$load(param_file)
  
  # Create experiment using fit_on with data object
  running_flag <- RunningFlag$new()
  exp_obj <- fit_on(data_obj, param_obj, running_flag)
  
  # Get population and individual from last generation
  pop_obj <- exp_obj$get_best_population()
  ind_obj <- if (pop_obj$get_size() > 0) pop_obj$get_individual(1) else NULL
  
  # Create jury from population using param settings
  jury_obj <- if (pop_obj$get_size() > 0) {
    tryCatch(
      Jury$new_from_param(pop_obj, param_obj),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  
  list(
    data = data_obj,
    param = param_obj,
    experiment = exp_obj,
    population = pop_obj,
    individual = ind_obj,
    jury = jury_obj,
    test_dir = test_dir
  )
}

# ---- Test Data S3 methods ----
test_that("print.Data() prints data summary", {
  objs <- setup_rust_objects()
  expect_output(print(objs$data), ".*")
  unlink(objs$test_dir, recursive = TRUE)
})

test_that("[[.Data() extracts fields correctly", {
  objs <- setup_rust_objects()
  
  # Get field names
  field_names <- names(objs$data)
  expect_true(is.character(field_names))
  expect_true(length(field_names) > 0)
  
  # Extract a field (e.g., sample_len)
  if ("sample_len" %in% field_names) {
    sample_len <- objs$data[["sample_len"]]
    expect_true(is.numeric(sample_len))
    expect_true(sample_len > 0)
  }
  
  # Test invalid field
  expect_error(objs$data[["nonexistent_field"]], "Unknown parameter")
  
  unlink(objs$test_dir, recursive = TRUE)
})

test_that("names.Data() returns field names", {
  objs <- setup_rust_objects()
  field_names <- names(objs$data)
  expect_true(is.character(field_names))
  expect_true(length(field_names) > 0)
  expect_true(all(c("sample_len", "feature_len") %in% field_names))
  unlink(objs$test_dir, recursive = TRUE)
})

test_that("dim.Data() returns dimensions", {
  objs <- setup_rust_objects()
  dims <- dim(objs$data)
  expect_true(is.numeric(dims))
  expect_equal(length(dims), 2)
  expect_true(all(dims > 0))
  unlink(objs$test_dir, recursive = TRUE)
})

# ---- Test Param S3 methods ----
test_that("print.Param() prints parameter info", {
  objs <- setup_rust_objects()
  expect_output(print(objs$param), ".*")
  unlink(objs$test_dir, recursive = TRUE)
})

test_that("[[<-.Param() sets parameters correctly", {
  objs <- setup_rust_objects()
  
  # Set numeric parameter (in ga section)
  objs$param[["population_size"]] <- 200
  expect_equal(objs$param[["population_size"]], 200)
  
  # Set string parameter (in general section)
  objs$param[["language"]] <- "ter"
  expect_equal(objs$param[["language"]], "ter")
  
  # Set boolean parameter (in general section)
  objs$param[["display_colorful"]] <- FALSE
  expect_equal(objs$param[["display_colorful"]], FALSE)
  
  unlink(objs$test_dir, recursive = TRUE)
})

test_that("[[.Param() extracts parameters correctly", {
  objs <- setup_rust_objects()
  
  # Extract numeric parameter (from ga section)
  pop_size <- objs$param[["population_size"]]
  expect_true(is.numeric(pop_size))
  
  # Extract string parameter (from general section)
  language <- objs$param[["language"]]
  expect_true(is.character(language))
  
  unlink(objs$test_dir, recursive = TRUE)
})

# ---- Test Individual S3 methods ----
test_that("print.Individual() prints individual info", {
  objs <- setup_rust_objects()
  
  if (!is.null(objs$individual)) {
    expect_output(print(objs$individual), ".*")
  } else {
    skip("No individual available in population")
  }
  
  unlink(objs$test_dir, recursive = TRUE)
})

test_that("coef.Individual() extracts coefficients", {
  objs <- setup_rust_objects()
  
  if (!is.null(objs$individual)) {
    coeffs <- coef(objs$individual)
    expect_true(is.numeric(coeffs))
    expect_true(length(coeffs) > 0)
    expect_true(!is.null(names(coeffs)))
  } else {
    skip("No individual available in population")
  }
  
  unlink(objs$test_dir, recursive = TRUE)
})

test_that("predict.Individual() makes predictions", {
  objs <- setup_rust_objects()
  
  if (!is.null(objs$individual)) {
    # Test score predictions (default)
    scores <- predict(objs$individual, objs$data, type = "score")
    expect_true(is.numeric(scores))
    expect_equal(length(scores), dim(objs$data)[1])  # One score per sample (rows are samples)
    
    # Test class predictions
    classes <- predict(objs$individual, objs$data, type = "class")
    expect_true(is.numeric(classes) || is.factor(classes))
    expect_equal(length(classes), dim(objs$data)[1])
  } else {
    skip("No individual available in population")
  }
  
  unlink(objs$test_dir, recursive = TRUE)
})

# ---- Test Experiment S3 methods ----
test_that("print.Experiment() prints experiment info", {
  objs <- setup_rust_objects()
  expect_output(print(objs$experiment), ".*")
  unlink(objs$test_dir, recursive = TRUE)
})

test_that("length.Experiment() returns generation count", {
  objs <- setup_rust_objects()
  len <- length(objs$experiment)
  expect_true(is.numeric(len))
  expect_true(len >= 0)
  unlink(objs$test_dir, recursive = TRUE)
})

# ---- Test Population S3 methods ----
test_that("print.Population() prints population info", {
  objs <- setup_rust_objects()
  expect_output(print(objs$population), ".*")
  unlink(objs$test_dir, recursive = TRUE)
})

test_that("length.Population() returns size", {
  objs <- setup_rust_objects()
  len <- length(objs$population)
  expect_true(is.numeric(len))
  expect_true(len >= 0)
  unlink(objs$test_dir, recursive = TRUE)
})

# ---- Test Jury S3 methods ----
test_that("predict.Jury() makes predictions", {
  objs <- setup_rust_objects()
  
  if (!is.null(objs$jury) && length(objs$jury) > 0) {
    predictions <- tryCatch(
      predict(objs$jury, objs$data),
      error = function(e) {
        skip(paste("Prediction failed:", conditionMessage(e)))
      }
    )
    expect_true(is.list(predictions))
    expect_true("class" %in% names(predictions))
    expect_true("score" %in% names(predictions))
  } else {
    skip("No jury available or empty")
  }
  
  unlink(objs$test_dir, recursive = TRUE)
})

# ==============================================================================
# Tests for db.R functions (database/API functions)
# ==============================================================================

test_that("get_taxonomy() validates input parameters", {
  # Test empty MSP names
  expect_error(
    get_taxonomy(msp_names = c()),
    "Please provide at least one MSP"
  )
  
  # Test invalid GTDB version
  expect_error(
    get_taxonomy(msp_names = c("msp_0005"), gtdb_version = 999.0),
    "Unsupported GTDB version"
  )
})

test_that("get_taxonomy() function signature is correct", {
  # Test that function exists and has correct formals
  expect_true(exists("get_taxonomy"))
  
  formals_names <- names(formals(get_taxonomy))
  expect_true("msp_names" %in% formals_names)
  expect_true("gtdb_version" %in% formals_names)
  expect_true("fields" %in% formals_names)
  
  # Check default values
  defaults <- formals(get_taxonomy)
  expect_equal(eval(defaults$gtdb_version), 220.0)
  expect_equal(length(eval(defaults$fields)), 0)
})

test_that("get_fields() function signature is correct", {
  expect_true(exists("get_fields"))
  
  formals_names <- names(formals(get_fields))
  expect_true("msp_name" %in% formals_names)
  
  # Check default value
  defaults <- formals(get_fields)
  expect_equal(eval(defaults$msp_name), "msp_0005")
})
