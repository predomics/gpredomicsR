# Additional tests for visualization and utility functions
library(testthat)
library(gpredomicsR)
library(ggplot2)

# ---- Setup: Generate test data ----
set.seed(456)

num_samples <- 50
num_features <- 100

# Create feature matrix X
X <- matrix(runif(num_samples * num_features, min = 0, max = 1), 
            nrow = num_features, 
            ncol = num_samples)

colnames(X) <- paste0("Sample_", seq_len(num_samples))
rownames(X) <- paste0("Feature_", seq_len(num_features))

# Create binary response vector
y <- sample(0:1, num_samples, replace = TRUE)

# Generate mock models
generateMockModel <- function(id) {
  num_selected_features <- sample(2:5, 1)
  feature_names <- if (!is.null(rownames(X))) rownames(X) else paste0("F", seq_len(nrow(X)))
  feature_indices <- seq_len(nrow(X))
  
  sel_idx <- sample(feature_indices, num_selected_features)
  sel_feats <- feature_names[sel_idx]
  sel_coeffs <- sample(c(-1, 1), num_selected_features, replace = TRUE)
  
  coeff_named <- sel_coeffs
  names(coeff_named) <- sel_feats
  
  list(
    features = sel_feats,
    coeff = sel_coeffs,
    coefficients = coeff_named,
    indexes = sel_idx,
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

mock_population <- lapply(1:20, generateMockModel)
mock_experiment <- list(
  rust = list(),
  params = list(),
  data = list(),
  model_collection = list(
    gen_1 = mock_population,
    gen_2 = mock_population,
    gen_3 = mock_population
  ),
  execTime = 123.45
)

# ---- Test visualization functions ----
test_that("printModel() formats models correctly", {
  model <- mock_population[[1]]
  
  # Test short format
  output_short <- printModel(model, method = "short", score = "fit")
  expect_true(is.character(output_short))
  expect_true(nchar(output_short) > 0)
  
  # Test str format
  output_str <- printModel(model, method = "str", score = "fit")
  expect_true(is.character(output_str))
  
  # Test with invalid score
  expect_null(printModel(model, method = "short", score = "invalid_score"))
  
  # Test with invalid object
  expect_null(printModel(list(a = 1, b = 2), method = "short"))
})

test_that("printPopulation() displays population info", {
  # Test digested format
  expect_output(printPopulation(mock_population, method = "digested", score = "fit"))
  
  # Test short format
  expect_output(printPopulation(mock_population, method = "short", score = "fit"))
  
  # Test str format
  expect_output(printPopulation(mock_population, method = "str", score = "fit"))
  
  # Test with invalid object
  expect_null(printPopulation(list(list(a = 1)), method = "short"))
})

test_that("printModelCollection() displays collection info", {
  model_collection <- list(
    gen_1 = mock_population[1:10],
    gen_2 = mock_population[11:20]
  )
  
  # Test short format
  expect_output(printModelCollection(model_collection, method = "short"))
  
  # Test long format
  expect_output(printModelCollection(model_collection, method = "long"))
  
  # Test with invalid object
  expect_null(printModelCollection(list(list(a = 1)), method = "short"))
})

test_that("printExperiment() displays experiment info", {
  # Test with valid experiment
  expect_output(printExperiment(mock_experiment))
  
  # Test with invalid object
  expect_null(printExperiment(list(a = 1, b = 2)))
})

test_that("plotHistory() generates plot without errors", {
  # Create a simple history data frame
  history_df <- data.frame(
    epoch = 1:10,
    best_fit = seq(0.7, 0.95, length.out = 10),
    mean_fit = seq(0.65, 0.9, length.out = 10)
  )
  
  # This should not error
  expect_silent(p <- plotHistory(mock_experiment, attributes = c("fit", "auc"), plot = FALSE))
})

test_that("analyzePopulationFeatures() analyzes features correctly", {
  # This function analyzes feature usage in a population
  result <- analyzePopulationFeatures(mock_population, X = X, y = y)
  
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  expect_true("feature" %in% names(result))
  expect_true("n_models" %in% names(result))
})

test_that("plotAbundanceByClass() generates plot without errors", {
  # Create data for plotting
  features_to_plot <- rownames(X)[1:5]
  
  # Should not error
  expect_silent(
    p <- plotAbundanceByClass(
      features = features_to_plot,
      X = X,
      y = y,
      plot = FALSE
    )
  )
})

test_that("plotPrevalence() generates plot without errors", {
  features_to_plot <- rownames(X)[1:10]
  
  # Should not error
  expect_silent(
    p <- plotPrevalence(
      features = features_to_plot,
      X = X,
      y = y,
      plot = FALSE
    )
  )
})

test_that("plotBarcode() generates plot without errors", {
  # Mock individual with features
  individual <- mock_population[[1]]
  
  # Should not error (plot = FALSE to avoid displaying)
  expect_silent(
    p <- plotBarcode(
      individual = individual,
      X = X,
      y = y,
      plot = FALSE
    )
  )
})

test_that("plotFeatureModelCoeffs() generates plot without errors", {
  # Should not error
  expect_silent(
    p <- plotFeatureModelCoeffs(
      population = mock_population,
      X = X,
      y = y,
      plot = FALSE
    )
  )
})

test_that("plotClassHeatmap() generates plot without errors", {
  features_to_plot <- rownames(X)[1:20]
  
  # Should not error
  expect_silent(
    p <- plotClassHeatmap(
      features = features_to_plot,
      X = X,
      y = y,
      plot = FALSE
    )
  )
})

test_that("plotImportanceHeatmap() generates plot without errors", {
  # Should not error
  expect_silent(
    p <- plotImportanceHeatmap(
      population = mock_population,
      X = X,
      y = y,
      top_n = 10,
      plot = FALSE
    )
  )
})

test_that("plot_population_heatmap() generates plot without errors", {
  # Should not error
  expect_silent(
    p <- plot_population_heatmap(
      population = mock_population,
      X = X,
      y = y,
      plot = FALSE
    )
  )
})

test_that("plotIndividualWaterfall() generates plot without errors", {
  individual <- mock_population[[1]]
  
  # Should not error
  expect_silent(
    p <- plotIndividualWaterfall(
      individual = individual,
      X = X,
      y = y,
      plot = FALSE
    )
  )
})

# ---- Test additional utility functions ----
test_that("printy() prints objects appropriately", {
  # Test with different types
  expect_output(printy(mock_population[[1]], type = "model"))
  expect_output(printy(mock_population, type = "population"))
  expect_output(printy(mock_experiment, type = "experiment"))
})

test_that("get_fields() extracts parameter fields correctly", {
  param_list <- list(
    general = list(seed = 123, algo = "ga"),
    ga = list(max_epochs = 100, population_size = 1000)
  )
  
  # Extract a nested field
  result <- get_fields(param_list, "max_epochs")
  expect_equal(result, 100)
  
  # Extract from first level
  result2 <- get_fields(param_list, "seed")
  expect_equal(result2, 123)
})

test_that("get_taxonomy() extracts taxonomy information", {
  # Create mock feature names with taxonomy
  feature_names <- c(
    "k__Bacteria|p__Firmicutes|c__Bacilli",
    "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria",
    "k__Bacteria|p__Firmicutes|c__Clostridia"
  )
  
  # Extract at different levels
  result_phylum <- get_taxonomy(feature_names, level = "p")
  expect_true(is.character(result_phylum))
  expect_equal(length(result_phylum), 3)
  expect_true(any(grepl("Firmicutes", result_phylum)))
  
  result_class <- get_taxonomy(feature_names, level = "c")
  expect_true(is.character(result_class))
  expect_true(any(grepl("Bacilli", result_class)))
})

test_that("confInterBinomial() calculates confidence intervals", {
  # Test with typical values
  ci <- confInterBinomial(accuracy = 0.85, n = 100)
  
  expect_true(is.numeric(ci))
  expect_true("inf" %in% names(ci))
  expect_true("sup" %in% names(ci))
  expect_true("accuracy" %in% names(ci))
  expect_true(ci["inf"] < ci["accuracy"])
  expect_true(ci["sup"] > ci["accuracy"])
  
  # Test edge cases
  ci_low <- confInterBinomial(accuracy = 0.05, n = 50)
  expect_true(ci_low["inf"] >= 0)
  expect_true(ci_low["inf"] < ci_low["accuracy"])
  
  ci_high <- confInterBinomial(accuracy = 0.95, n = 50)
  expect_true(ci_high["sup"] <= 1)
  expect_true(ci_high["sup"] > ci_high["accuracy"])
})

test_that("selectBestPopulation() filters population by confidence", {
  # Select best models with 95% confidence
  best_pop <- selectBestPopulation(mock_population, score = "fit", p = 0.05)
  
  expect_true(isPopulation(best_pop))
  expect_true(length(best_pop) <= length(mock_population))
  expect_true(length(best_pop) > 0)
  
  # All selected models should have high fit scores
  fit_scores <- sapply(best_pop, function(m) m$fit)
  expect_true(all(fit_scores >= min(fit_scores)))
})

test_that("sortPopulation() sorts by different criteria", {
  # Sort by fit
  sorted_fit <- sortPopulation(mock_population, evalToOrder = "fit")
  fit_values <- sapply(sorted_fit, function(m) m$fit)
  expect_true(all(diff(fit_values) <= 1e-10))  # Descending order
  
  # Sort by auc
  sorted_auc <- sortPopulation(mock_population, evalToOrder = "auc")
  auc_values <- sapply(sorted_auc, function(m) m$auc)
  expect_true(all(diff(auc_values) <= 1e-10))  # Descending order
  
  # Sort by k (sparsity)
  sorted_k <- sortPopulation(mock_population, evalToOrder = "k")
  k_values <- sapply(sorted_k, function(m) m$k)
  expect_true(all(diff(k_values) >= 0))  # Ascending order for k
})

test_that("getTheBestIndividual() returns best model", {
  # Get best by fit
  best_fit <- getTheBestIndividual(mock_population, evalToFit = "fit")
  expect_true(isModel(best_fit))
  
  # Should be the model with highest fit
  all_fits <- sapply(mock_population, function(m) m$fit)
  expect_equal(best_fit$fit, max(all_fits))
  
  # Get best by auc
  best_auc <- getTheBestIndividual(mock_population, evalToFit = "auc")
  expect_true(isModel(best_auc))
  all_aucs <- sapply(mock_population, function(m) m$auc)
  expect_equal(best_auc$auc, max(all_aucs))
})

test_that("populationGet_X() extracts attributes correctly", {
  # Extract k values
  k_extractor <- populationGet_X("k", toVec = TRUE, na.rm = FALSE)
  k_values <- k_extractor(mock_population)
  
  expect_true(is.numeric(k_values))
  expect_equal(length(k_values), length(mock_population))
  
  # Extract features (should be longer due to multiple features per model)
  features_extractor <- populationGet_X("features", toVec = TRUE, na.rm = FALSE)
  features <- features_extractor(mock_population)
  
  expect_true(is.character(features))
  expect_true(length(features) > length(mock_population))
  
  # Extract coefficients
  coeff_extractor <- populationGet_X("coeff", toVec = TRUE, na.rm = FALSE)
  coeffs <- coeff_extractor(mock_population)
  
  expect_true(is.numeric(coeffs))
})

test_that("modelToDenseVec() converts model to dense vector", {
  model <- mock_population[[1]]
  natts <- nrow(X)
  
  dense_vec <- modelToDenseVec(natts, model)
  
  expect_true(is.numeric(dense_vec))
  expect_equal(length(dense_vec), natts)
  
  # Check that selected features have non-zero coefficients
  for (i in seq_along(model$indexes)) {
    idx <- model$indexes[i]
    expect_equal(dense_vec[idx], model$coeff[i])
  }
})

test_that("listOfModelsToListOfDenseVec() converts multiple models", {
  dense_list <- listOfModelsToListOfDenseVec(X = X, y = y, list.models = mock_population)
  
  expect_true(is.list(dense_list))
  expect_equal(length(dense_list), length(mock_population))
  
  # Each element should be a dense vector
  for (vec in dense_list) {
    expect_true(is.numeric(vec))
    expect_equal(length(vec), nrow(X))
  }
})

test_that("listOfModelsToDenseCoefMatrix() creates coefficient matrix", {
  dense_matrix <- listOfModelsToDenseCoefMatrix(X, y, mock_population)
  
  expect_true(is.matrix(dense_matrix))
  expect_equal(ncol(dense_matrix), length(mock_population))
  
  # Number of rows should be the number of unique features used
  all_indexes <- unlist(lapply(mock_population, function(m) m$indexes))
  expect_equal(nrow(dense_matrix), length(unique(all_indexes)))
})

test_that("populationToDataFrame() converts population correctly", {
  df <- populationToDataFrame(mock_population)
  
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), length(mock_population))
  
  # Check expected columns
  expect_true("fit" %in% names(df))
  expect_true("auc" %in% names(df))
  expect_true("k" %in% names(df))
  expect_true("sensitivity" %in% names(df))
  expect_true("specificity" %in% names(df))
  expect_true("accuracy" %in% names(df))
})

test_that("analyzeAttributeEvolution() tracks attribute changes", {
  df_evo <- analyzeAttributeEvolution(
    mock_experiment,
    attributes = c("fit", "auc", "k"),
    best_model = TRUE,
    plot = FALSE
  )
  
  expect_true(is.data.frame(df_evo))
  expect_true(nrow(df_evo) > 0)
  expect_true("generation" %in% names(df_evo))
  expect_true("fit" %in% names(df_evo) || "auc" %in% names(df_evo))
})

test_that("check.X_y() validates input data", {
  # Valid case
  expect_silent(check.X_y(X, y))
  
  # Invalid: mismatched dimensions
  y_wrong <- sample(0:1, num_samples + 10, replace = TRUE)
  expect_error(check.X_y(X, y_wrong))
  
  # Invalid: X not a matrix
  expect_error(check.X_y(as.vector(X), y))
})

test_that("filterfeaturesK() selects top features", {
  # Classification mode
  result <- filterfeaturesK(X, y, k = 10, type = "wilcoxon")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 10)
  expect_true(all(c("p", "q", "status") %in% colnames(result)))
  
  # Test return.data option
  filtered_X <- filterfeaturesK(X, y, k = 10, return.data = TRUE)
  expect_true(is.matrix(filtered_X))
  expect_equal(nrow(filtered_X), 10)
  
  # Regression mode
  y_numeric <- rnorm(num_samples)
  result_reg <- filterfeaturesK(X, y_numeric, k = 10, type = "spearman")
  expect_true(is.data.frame(result_reg))
  expect_true("rho" %in% colnames(result_reg))
})

test_that("getFeaturePrevalence() computes prevalence", {
  features_test <- rownames(X)[1:10]
  
  # Without class labels
  prev_all <- getFeaturePrevalence(features_test, X)
  expect_true(is.list(prev_all))
  expect_true("all" %in% names(prev_all))
  expect_equal(length(prev_all$all), 10)
  
  # With class labels
  prev_by_class <- getFeaturePrevalence(features_test, X, y = y)
  expect_true(all(c("all", "0", "1") %in% names(prev_by_class)))
  expect_equal(length(prev_by_class$all), 10)
})

test_that("computeCardEnrichment() performs chi-square test", {
  # Create cardinality matrix
  v.card.mat <- matrix(c(25, 15, 35, 25), nrow = 2, byrow = TRUE)
  rownames(v.card.mat) <- c("0", "1")
  colnames(v.card.mat) <- c("Feature_1", "Feature_2")
  
  result <- computeCardEnrichment(v.card.mat, y)
  
  expect_true(is.list(result))
  expect_true("chisq.p" %in% names(result))
  expect_true("chisq.q" %in% names(result))
  expect_true("card.all" %in% names(result))
  expect_equal(length(result$chisq.p), 2)
})
