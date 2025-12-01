# Load required packages
library(testthat)
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
