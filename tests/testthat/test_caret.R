library(testthat)
library(gpredomicsR)
library(caret)

# Load caret interface
source("../../R/gpredomics_caret.R")

# ==============================================================================
# Minimal test data
# ==============================================================================

setup_test_data <- function() {
  set.seed(123)
  
  # Create simple but valid data
  n_samples <- 40
  n_features <- 20
  
  X <- matrix(
    runif(n_features * n_samples, 0, 10),
    nrow = n_features,
    ncol = n_samples
  )
  rownames(X) <- paste0("Feature_", 1:n_features)
  colnames(X) <- paste0("Sample_", 1:n_samples)
  
  y <- sample(0:1, n_samples, replace = TRUE)
  names(y) <- paste0("Sample_", 1:n_samples)  # Add names to y
  
  # Convert to Data for Gpredomics
  gp_data <- as.gpredomics.data(X, y, TRUE)
  
  # Create data.frame for caret
  df <- as.data.frame(t(X))
  df$Class <- factor(ifelse(y == 1, "Pos", "Neg"), levels = c("Neg", "Pos"))
  
  list(
    X = X,
    y = y,
    gp_data = gp_data,
    df = df
  )
}

# ==============================================================================
# Test 1: Direct comparison Gpredomics vs caret (fixed seed)
# ==============================================================================

test_that("caret produces same results as direct Gpredomics with fixed seed", {
  data <- setup_test_data()
  
  # Identical parameters
  seed_value <- 456
  n_gen <- 5
  pop_size <- 10
  k_val <- 5
  
  # --- Approach 1: Direct Gpredomics ---
  set.seed(seed_value)
  param_direct <- Param$new()
  param_direct$set("max_epochs", n_gen)
  param_direct$set("population_size", pop_size)
  param_direct$set("k_min", k_val)
  param_direct$set("k_max", k_val)
  param_direct$set("cv", FALSE)  # No internal CV
  param_direct$set("seed", seed_value)
  
  exp_direct <- Experiment$new(data$gp_data, param_direct)
  exp_direct$run()
  
  best_pop_direct <- exp_direct$get_best_population()
  best_ind_direct <- best_pop_direct$get_individual(1)
  
  # --- Approach 2: Via caret (no CV for direct comparison) ---
  set.seed(seed_value)
  
  # trainControl without resampling
  ctrl <- trainControl(
    method = "none",
    summaryFunction = gpredomicsSummary,
    classProbs = TRUE,
    savePredictions = FALSE
  )
  
  # Grid with same parameters
  tune_grid <- expand.grid(
    max_epochs = n_gen,
    population_size = pop_size,
    ga_kmin = k_val,
    ga_kmax = k_val,
    beam_kmin = k_val,
    beam_kmax = k_val
  )
  
  model_caret <- train_gpredomics_caret(
    Class ~ .,
    data = data$df,
    method = getModelInfo_gpredomics(prediction_mode = "best"),
    trControl = ctrl,
    tuneGrid = tune_grid,
    metric = "ROC",
    seed = seed_value  # Seed passed explicitly
  )
  
  # --- Comparisons ---
  
  # 1. Verify experiment is stored
  expect_true(!is.null(model_caret$finalModel$experiment))
  expect_true(isExperiment(model_caret$finalModel$experiment))
  
  # 2. Verify number of generations
  pop_direct <- exp_direct$get_population()
  pop_caret <- model_caret$finalModel$experiment$get_population()
  
  expect_equal(length(pop_direct$get_history()), length(pop_caret$get_history()))
  
  # 3. Verify best individual has same characteristics
  best_ind_caret <- pop_caret$get_individual(1)
  
  expect_equal(best_ind_direct$get_k(), best_ind_caret$get_k())
  
  # 4. Verify selected features (should be identical with same seed)
  features_direct <- best_ind_direct$get_features()
  features_caret <- best_ind_caret$get_features()
  
  expect_equal(sort(features_direct), sort(features_caret))
  
  # 5. Verify coefficients (should be identical)
  coeffs_direct <- best_ind_direct$get_coefficients()
  coeffs_caret <- best_ind_caret$get_coefficients()
  
  expect_equal(coeffs_direct[sort(names(coeffs_direct))], 
               coeffs_caret[sort(names(coeffs_caret))])
})

# ==============================================================================
# Test 2: Verification of correct parameter transmission
# ==============================================================================

test_that("All grid parameters are correctly transmitted to Gpredomics", {
  data <- setup_test_data()
  
  # Grid with multiple parameters
  tune_grid <- expand.grid(
    max_epochs = 3,
    population_size = 8,
    ga_kmin = 4,
    ga_kmax = 4,
    beam_kmin = 4,
    beam_kmax = 4,
    select_elite_pct = 0.25,
    mutated_children_pct = 0.15
  )
  
  ctrl <- trainControl(
    method = "none",
    summaryFunction = gpredomicsSummary,
    classProbs = TRUE
  )
  
  set.seed(789)
  model <- train_gpredomics_caret(
    Class ~ .,
    data = data$df,
    method = getModelInfo_gpredomics(prediction_mode = "best"),
    trControl = ctrl,
    tuneGrid = tune_grid,
    metric = "Accuracy"
  )
  
  # Get the Param used via experiment
  # Note: Cannot access Param directly, but can verify
  # results (k, population size, etc.)
  
  best_pop <- model$finalModel$experiment$get_best_population()
  best_ind <- best_pop$get_individual(1)
  
  # Verify k
  expect_equal(best_ind$get_k(), tune_grid$ga_kmin)
  
  # Verify population size (number of individuals)
  population <- model$finalModel$experiment$get_population()
  history <- population$get_history()
  
  # Last generation should have 'population_size' individuals
  last_gen <- history[[length(history)]]
  expect_equal(length(last_gen), tune_grid$population_size)
  
  # Verify number of generations
  expect_equal(length(history), tune_grid$max_epochs)
})

# ==============================================================================
# Test 3: Prediction modes (best, fbm, jury, percentage)
# ==============================================================================

test_that("All 4 prediction modes produce coherent results", {
  data <- setup_test_data()
  
  # Train a simple model
  tune_grid <- expand.grid(
    max_epochs = 5,
    population_size = 10,
    ga_kmin = 5,
    ga_kmax = 5,
    beam_kmin = 5,
    beam_kmax = 5
  )
  
  ctrl <- trainControl(
    method = "none",
    summaryFunction = gpredomicsSummary,
    classProbs = TRUE
  )
  
  set.seed(999)
  
  # Test each mode
  modes <- c("best", "fbm", "jury", "percentage")
  
  for (mode in modes) {
    model <- train_gpredomics_caret(
      Class ~ .,
      data = data$df,
      method = getModelInfo_gpredomics(
        prediction_mode = mode,
        fbm_alpha = 0.05,
        percentage = 20
      ),
      trControl = ctrl,
      tuneGrid = tune_grid,
      metric = "ROC"
    )
    
    # Verify mode is stored
    expect_equal(model$finalModel$gpredomics_settings$prediction_mode, mode)
    
    # Verify predictions work
    preds_class <- predict(model, newdata = data$df, type = "raw")
    preds_prob <- predict(model, newdata = data$df, type = "prob")
    
    expect_true(is.factor(preds_class))
    expect_equal(length(preds_class), nrow(data$df))
    
    expect_true(is.data.frame(preds_prob))
    expect_equal(nrow(preds_prob), nrow(data$df))
    expect_true(all(c("Neg", "Pos") %in% colnames(preds_prob)))
    
    # Verify probabilities sum to 1
    expect_true(all(abs(rowSums(preds_prob) - 1) < 1e-10))
  }
})

# ==============================================================================
# Test 4: Identical predictions with mode override
# ==============================================================================

test_that("Prediction mode override produces coherent results", {
  data <- setup_test_data()
  
  tune_grid <- expand.grid(
    max_epochs = 5,
    population_size = 12,
    ga_kmin = 6,
    ga_kmax = 6,
    beam_kmin = 6,
    beam_kmax = 6
  )
  
  ctrl <- trainControl(
    method = "none",
    summaryFunction = gpredomicsSummary,
    classProbs = TRUE
  )
  
  set.seed(111)
  
  # Train with "best" mode
  model_best <- train_gpredomics_caret(
    Class ~ .,
    data = data$df,
    method = getModelInfo_gpredomics(prediction_mode = "best"),
    trControl = ctrl,
    tuneGrid = tune_grid,
    metric = "ROC"
  )
  
  # Predict with "best" mode (default)
  preds_best_1 <- predict(model_best, newdata = data$df, type = "prob")
  
  # Predict with "best" mode (explicit)
  preds_best_2 <- predict(model_best, newdata = data$df, type = "prob", 
                          mode = "best")
  
  # Should be identical
  expect_equal(preds_best_1, preds_best_2)
  
  # Predict with "jury" mode (override)
  preds_jury <- predict(model_best, newdata = data$df, type = "prob", 
                        mode = "jury")
  
  # Should NOT be identical (unless by chance)
  expect_false(identical(preds_best_1, preds_jury))
  
  # But dimensions should be the same
  expect_equal(dim(preds_best_1), dim(preds_jury))
})

# ==============================================================================
# Test 5: Metrics computed correctly
# ==============================================================================

test_that("gpredomicsSummary computes all metrics correctly", {
  # Test data for metrics
  set.seed(222)
  n <- 100
  
  # Create predictions and observations
  obs <- factor(sample(c("Neg", "Pos"), n, replace = TRUE), 
                levels = c("Neg", "Pos"))
  
  # Simulated probabilities
  Pos_prob <- runif(n)
  Neg_prob <- 1 - Pos_prob
  
  data <- data.frame(
    obs = obs,
    pred = factor(ifelse(Pos_prob > 0.5, "Pos", "Neg"), levels = c("Neg", "Pos")),
    Neg = Neg_prob,
    Pos = Pos_prob
  )
  
  # Compute metrics
  metrics <- gpredomicsSummary(data, lev = c("Neg", "Pos"))
  
  # Verify all metrics are present
  expected_metrics <- c("ROC", "Sens", "Spec", "Accuracy", 
                       "MCC", "NPV", "PPV", "F1", "Gmean")
  
  expect_true(all(expected_metrics %in% names(metrics)))
  expect_equal(length(metrics), 9)
  
  # Verify metrics are in [0, 1] (except MCC in [-1, 1])
  expect_true(all(metrics[c("ROC", "Sens", "Spec", "Accuracy", 
                           "NPV", "PPV", "F1", "Gmean")] >= 0 & 
                 metrics[c("ROC", "Sens", "Spec", "Accuracy", 
                           "NPV", "PPV", "F1", "Gmean")] <= 1))
  
  expect_true(metrics["MCC"] >= -1 && metrics["MCC"] <= 1)
  
  # Verify manual accuracy calculation
  cm <- table(data$pred, data$obs)
  manual_acc <- sum(diag(cm)) / sum(cm)
  expect_equal(as.numeric(metrics["Accuracy"]), manual_acc, tolerance = 1e-10)
})

# ==============================================================================
# Test 6: Seed handling with CV
# ==============================================================================

test_that("Seed reproduces identical results with CV", {
  data <- setup_test_data()
  
  tune_grid <- expand.grid(
    max_epochs = 3,
    population_size = 8,
    ga_kmin = 4,
    ga_kmax = 4,
    beam_kmin = 4,
    beam_kmax = 4
  )
  
  # 3-fold CV
  ctrl <- trainControl(
    method = "cv",
    number = 3,
    summaryFunction = gpredomicsSummary,
    classProbs = TRUE,
    savePredictions = FALSE
  )
  
  seed_val <- 333
  
  # First run
  set.seed(seed_val)
  model_1 <- train_gpredomics_caret(
    Class ~ .,
    data = data$df,
    method = getModelInfo_gpredomics(prediction_mode = "best"),
    trControl = ctrl,
    tuneGrid = tune_grid,
    metric = "ROC",
    seed = seed_val
  )
  
  # Second run with same seed
  set.seed(seed_val)
  model_2 <- train_gpredomics_caret(
    Class ~ .,
    data = data$df,
    method = getModelInfo_gpredomics(prediction_mode = "best"),
    trControl = ctrl,
    tuneGrid = tune_grid,
    metric = "ROC",
    seed = seed_val
  )
  
  # CV metrics should be identical
  expect_equal(model_1$results$ROC, model_2$results$ROC, tolerance = 1e-10)
  expect_equal(model_1$results$Sens, model_2$results$Sens, tolerance = 1e-10)
  expect_equal(model_1$results$Spec, model_2$results$Spec, tolerance = 1e-10)
  
  # Predictions on new data should be identical
  preds_1 <- predict(model_1, newdata = data$df, type = "prob")
  preds_2 <- predict(model_2, newdata = data$df, type = "prob")
  
  expect_equal(preds_1, preds_2)
})

# ==============================================================================
# Test 7: Detailed population comparison
# ==============================================================================

test_that("caret population identical to direct Gpredomics population", {
  data <- setup_test_data()
  
  seed_value <- 555
  n_gen <- 4
  pop_size <- 10
  k_val <- 5
  
  # --- Direct Gpredomics ---
  set.seed(seed_value)
  param_direct <- Param$new()
  param_direct$set("max_epochs", n_gen)
  param_direct$set("population_size", pop_size)
  param_direct$set("k_min", k_val)
  param_direct$set("k_max", k_val)
  param_direct$set("cv", FALSE)
  param_direct$set("seed", seed_value)
  
  exp_direct <- Experiment$new(data$gp_data, param_direct)
  exp_direct$run()
  
  pop_direct <- exp_direct$get_population()
  history_direct <- pop_direct$get_history()
  
  # --- Via caret ---
  set.seed(seed_value)
  
  ctrl <- trainControl(
    method = "none",
    summaryFunction = gpredomicsSummary,
    classProbs = TRUE
  )
  
  tune_grid <- expand.grid(
    max_epochs = n_gen,
    population_size = pop_size,
    ga_kmin = k_val,
    ga_kmax = k_val,
    beam_kmin = k_val,
    beam_kmax = k_val
  )
  
  model_caret <- train_gpredomics_caret(
    Class ~ .,
    data = data$df,
    method = getModelInfo_gpredomics(prediction_mode = "best"),
    trControl = ctrl,
    tuneGrid = tune_grid,
    metric = "ROC",
    seed = seed_value
  )
  
  pop_caret <- model_caret$finalModel$experiment$get_population()
  history_caret <- pop_caret$get_history()
  
  # --- Detailed comparisons ---
  
  # 1. Same number of generations
  expect_equal(length(history_direct), length(history_caret))
  expect_equal(length(history_direct), n_gen)
  
  # 2. Same population size in each generation
  for (gen in seq_along(history_direct)) {
    expect_equal(length(history_direct[[gen]]), length(history_caret[[gen]]))
    expect_equal(length(history_direct[[gen]]), pop_size)
  }
  
  # 3. Compare all individuals from last generation
  last_gen_direct <- history_direct[[n_gen]]
  last_gen_caret <- history_caret[[n_gen]]
  
  for (i in seq_along(last_gen_direct)) {
    ind_direct <- pop_direct$get_individual(i)
    ind_caret <- pop_caret$get_individual(i)
    
    # Same k
    expect_equal(ind_direct$get_k(), ind_caret$get_k())
    
    # Same features (order may differ)
    expect_setequal(ind_direct$get_features(), ind_caret$get_features())
    
    # Same coefficients
    coeffs_direct <- ind_direct$get_coefficients()
    coeffs_caret <- ind_caret$get_coefficients()
    
    # Sort by names for comparison
    expect_equal(coeffs_direct[sort(names(coeffs_direct))],
                 coeffs_caret[sort(names(coeffs_caret))])
  }
})

# ==============================================================================
# Test 8: Identical predictions between approaches
# ==============================================================================

test_that("caret predictions identical to direct Gpredomics predictions", {
  data <- setup_test_data()
  
  seed_value <- 2468
  
  # Common configuration
  n_gen <- 5
  pop_size <- 10
  k_val <- 6
  
  # --- Direct Gpredomics ---
  set.seed(seed_value)
  param_direct <- Param$new()
  param_direct$set("max_epochs", n_gen)
  param_direct$set("population_size", pop_size)
  param_direct$set("k_min", k_val)
  param_direct$set("k_max", k_val)
  param_direct$set("cv", FALSE)
  param_direct$set("seed", seed_value)
  
  exp_direct <- Experiment$new(data$gp_data, param_direct)
  exp_direct$run()
  
  best_pop_direct <- exp_direct$get_best_population()
  best_ind_direct <- best_pop_direct$get_individual(1)
  
  # Predictions with best individual
  preds_direct <- best_ind_direct$predict(data$gp_data)
  
  # --- Via caret ---
  set.seed(seed_value)
  
  ctrl <- trainControl(
    method = "none",
    summaryFunction = gpredomicsSummary,
    classProbs = TRUE
  )
  
  tune_grid <- expand.grid(
    max_epochs = n_gen,
    population_size = pop_size,
    ga_kmin = k_val,
    ga_kmax = k_val,
    beam_kmin = k_val,
    beam_kmax = k_val
  )
  
  model_caret <- train_gpredomics_caret(
    Class ~ .,
    data = data$df,
    method = getModelInfo_gpredomics(prediction_mode = "best"),  # "best" mode to compare with best individual
    trControl = ctrl,
    tuneGrid = tune_grid,
    metric = "ROC",
    seed = seed_value
  )
  
  # Predictions via caret
  preds_caret_prob <- predict(model_caret, newdata = data$df, type = "prob")
  preds_caret_class <- predict(model_caret, newdata = data$df, type = "raw")
  
  # --- Comparisons ---
  
  # Rust predictions return list(class, score)
  # class: integer vector 0/1
  # score: probability vector for class 1
  
  expect_equal(length(preds_direct$class), nrow(data$df))
  expect_equal(length(preds_direct$score), nrow(data$df))
  
  # Compare predicted classes
  # caret: factor "Neg"/"Pos"
  # Gpredomics: integers 0/1
  caret_class_int <- ifelse(preds_caret_class == "Pos", 1, 0)
  expect_equal(preds_direct$class, caret_class_int)
  
  # Compare scores/probabilities
  # caret: "Pos" column of probs data.frame
  # Gpredomics: score vector
  expect_equal(preds_direct$score, preds_caret_prob$Pos, tolerance = 1e-10)
})

# ==============================================================================
# Test 9: Grid search with multiple parameters
# ==============================================================================

test_that("Grid search correctly explores all parameters", {
  data <- setup_test_data()
  
  # Grid 2x2 = 4 combinations
  tune_grid <- expand.grid(
    max_epochs = c(3, 5),
    population_size = c(8, 12),
    ga_kmin = 5,
    ga_kmax = 5,
    beam_kmin = 5,
    beam_kmax = 5
  )
  
  ctrl <- trainControl(
    method = "none",
    summaryFunction = gpredomicsSummary,
    classProbs = TRUE
  )
  
  set.seed(777)
  model <- train_gpredomics_caret(
    Class ~ .,
    data = data$df,
    method = getModelInfo_gpredomics(prediction_mode = "best"),
    trControl = ctrl,
    tuneGrid = tune_grid,
    metric = "ROC"
  )
  
  # Verify results contains 4 rows (4 combinations)
  expect_equal(nrow(model$results), 4)
  
  # Verify all combinations are present
  expect_setequal(model$results$max_epochs, c(3, 5))
  expect_setequal(model$results$population_size, c(8, 12))
  
  # Verify bestTune points to one of the combinations
  expect_true(model$bestTune$max_epochs %in% c(3, 5))
  expect_true(model$bestTune$population_size %in% c(8, 12))
  expect_equal(model$bestTune$ga_kmin, 5)
  expect_equal(model$bestTune$ga_kmax, 5)
})

# ==============================================================================
# Test 10: Verification of string parameters (formulas, operators)
# ==============================================================================

test_that("String parameters are correctly transmitted", {
  data <- setup_test_data()
  
  # Parameters with language and data_type
  tune_grid <- expand.grid(
    max_epochs = 3,
    population_size = 8,
    ga_kmin = 4,
    ga_kmax = 4,
    beam_kmin = 4,
    beam_kmax = 4,
    language = "ter",  # String parameter
    data_type = "log"  # String parameter
  )
  
  ctrl <- trainControl(
    method = "none",
    summaryFunction = gpredomicsSummary,
    classProbs = TRUE
  )
  
  set.seed(888)
  
  # Should work without error
  expect_error({
    model <- train_gpredomics_caret(
      Class ~ .,
      data = data$df,
      method = getModelInfo_gpredomics(prediction_mode = "best"),
      trControl = ctrl,
      tuneGrid = tune_grid,
      metric = "ROC"
    )
  }, NA)  # NA = no error expected
})