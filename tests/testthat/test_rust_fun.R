# ==============================================================================
# Unit tests for Rust structures in gpredomicsR
# ==============================================================================

library(testthat)
library(gpredomicsR)

# ==============================================================================
# 1. Tests for RunningFlag
# ==============================================================================

test_that("RunningFlag: creation and state management", {
  flag <- RunningFlag$new()
  
  # Test: initial state should be "running"
  expect_true(flag$is_running())
  
  # Test: stopping the flag
  flag$stop()
  expect_false(flag$is_running())
  
  # Test: resetting the flag
  flag$reset()
  expect_true(flag$is_running())
})

test_that("RunningFlag: multiple stops do not fail", {
  flag <- RunningFlag$new()
  flag$stop()
  expect_silent(flag$stop()) 
  expect_false(flag$is_running())
})

# ==============================================================================
# 2. Tests for Param
# ==============================================================================

test_that("Param: loading from YAML file", {
  skip_if_not(file.exists("../../sample/param.yaml"), 
              "param.yaml file not found")
  
  param <- Param$load("../../sample/param.yaml")
  
  # Test: Param object should be created
  expect_true(!is.null(param))
  
  # Test: get() should return an R object with parameters
  param_data <- param$get()
  expect_true(is.list(param_data))
})

test_that("Param: parameter modification", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  
  param <- Param$load("../../sample/param.yaml")
  
  # Test: set() for numeric value
  expect_silent(param$set("max_epochs", 50))
  
  # Test: set_string() for string value
  expect_silent(param$set_string("X", "test_path.tsv"))
  
  # Test: set_bool() for boolean value
  expect_silent(param$set_bool("gpu", FALSE))
})

test_that("Param: address() returns memory address", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  
  param <- Param$load("../../sample/param.yaml")
  addr <- param$address()
  
  expect_true(is.character(addr))
  expect_true(grepl("^0x", addr))  # Hexadecimal format
})

# ==============================================================================
# 3. Tests for Data
# ==============================================================================

test_that("Data: creation via as_gpredomics_data", {
  # Create synthetic data
  set.seed(123)
  X <- matrix(runif(100 * 20), nrow = 100, ncol = 20)
  rownames(X) <- paste0("Feature_", 1:100)
  colnames(X) <- paste0("Sample_", 1:20)
  
  y <- factor(sample(c("control", "case"), 20, replace = TRUE),
              levels = c("control", "case"),
              ordered = TRUE)
  names(y) <- colnames(X)
  
  # Test: Data creation
  data <- as_gpredomics_data(as.data.frame(X), y, FALSE, NULL, NULL, NULL, NULL)
  
  expect_true(!is.null(data))
  
  # Test: get() returns the data
  data_content <- data$get()
  expect_true(is.list(data_content))
  expect_true("X" %in% names(data_content))
  expect_true("y" %in% names(data_content))
  expect_true("classes" %in% names(data_content))
})

test_that("Data: handling transposed data", {
  set.seed(456)
  X <- matrix(runif(50 * 10), nrow = 50, ncol = 10)
  rownames(X) <- paste0("Feature_", 1:50)
  colnames(X) <- paste0("Sample_", 1:10)
  
  y <- factor(sample(0:1, 10, replace = TRUE))
  names(y) <- colnames(X)
  
  # Test with transposed data
  X_t <- as.data.frame(t(X))
  data_t <- as_gpredomics_data(X_t, y, TRUE, NULL, NULL, NULL, NULL)
  
  data_content <- data_t$get()
  expect_equal(nrow(data_content$X), 50)  # Should have 50 features
  expect_equal(ncol(data_content$X), 10)  # And 10 samples
})

test_that("Data: address() and print() work correctly", {
  set.seed(789)
  X <- matrix(runif(20 * 5), nrow = 20, ncol = 5)
  rownames(X) <- paste0("F", 1:20)
  colnames(X) <- paste0("S", 1:5)
  y <- factor(sample(0:1, 5, replace = TRUE))
  names(y) <- colnames(X)
  
  data <- as_gpredomics_data(as.data.frame(X), y, FALSE, NULL, NULL, NULL, NULL)
  
  # Test address
  addr <- data$address()
  expect_true(is.character(addr))
  expect_true(grepl("^0x", addr))
  
  # Test print
  print_out <- data$print()
  expect_true(is.character(print_out))
  expect_true(nchar(print_out) > 0)
})


# ---- Test data loader ----
test_that("as_gpredomics_data and Experiment$load_data build valid gpredomics_data objects", {
    # Launch a quick experiment to control that data is the same 
    try(GLogger$level(level = glog_level), silent = TRUE)
    running_flag <- RunningFlag$new()
    paramRust <- Param$load("../../sample/param.yaml")
    paramRust$set("max_epochs", 10)
    paramRust$set_string("X", "../../sample/Xtrain.tsv")
    paramRust$set_string("y", "../../sample/Ytrain.tsv")
    paramRust$set_string("Xtest", "../../sample/Xtest.tsv")
    paramRust$set_string("ytest", "../../sample/Ytest.tsv")

    # Run the genetic algorithm
    expRust <- fit(paramRust, running_flag)

    # Test as_gpredomics_data
    X <- read.table(paste(getwd(), "../../sample/Xtrain.tsv", sep = "/"), h=T, row.names=1)
    colnames(X) <- gsub("\\.", "-", colnames(X), , perl = TRUE)

    y_df <- read.table(paste(getwd(), "../../sample/Ytrain.tsv", sep = "/"), h=T)
    y <- factor(ifelse(y_df[,1]==0, "healthy", "cirrhosis"),
            levels = c("healthy", "cirrhosis"),
            ordered = TRUE)
    names(y) <- rownames(y_df)

    data <- as_gpredomics_data(X, y, FALSE, NULL, NULL, NULL, NULL)

    X_t <- as.data.frame(t(X))
    data_t <- as_gpredomics_data(X_t, y, TRUE, NULL, NULL, NULL, NULL)

    expect_equal(data$get()$X, expRust$train_data()$get()$X)
    expect_equal(data$get()$y, expRust$train_data()$get()$y)
    expect_equal(data$get()$classes, expRust$train_data()$get()$classes)
    expect_equal(data$get(), data_t$get())

    # Test load_data - should work with compatible data (same structure as training)
    # Note: load_data creates a NEW Data object, it doesn't modify the experiment
    x_path <- normalizePath("../../sample/Xtest.tsv", mustWork = TRUE)
    y_path <- normalizePath("../../sample/Ytest.tsv", mustWork = TRUE)
    
    # This should work - loading test data that's compatible with train data
    # features_in_columns = FALSE because data has features in rows (like training data)
    data_load <- expRust$load_data(x_path, y_path, TRUE)
    
    expect_true(!is.null(data_load))
    expect_true(!is.null(data_load$get()$X))
    expect_true(!is.null(data_load$get()$y))
    expect_equal(data_load$get()$classes, expRust$train_data()$get()$classes)
})

# ==============================================================================
# 4. Tests for Individual
# ==============================================================================

test_that("Individual: extraction and metrics from Experiment", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  # Create a small experiment
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 5)
  param$set("population_size", 1000)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: retrieving an individual
  ind <- exp$individual(generation = 1, order = 1)
  expect_true(!is.null(ind))
  
  # Test: get() returns individual information
  ind_data <- ind$get()
  expect_true(is.list(ind_data))
  expect_true("features" %in% names(ind_data))
  expect_true("auc" %in% names(ind_data))
  expect_true("fit" %in% names(ind_data))
  expect_true("sensitivity" %in% names(ind_data))
  expect_true("specificity" %in% names(ind_data))
})

test_that("Individual: get_metrics, compute_metrics", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  ind <- exp$individual(1, 1)
  test_data <- exp$test_data()
  
  # Test get_metrics (base metrics already computed)
  ind_metrics_base <- ind$get_metrics()
  expect_true(!is.null(ind_metrics_base))
  expect_true(is.numeric(ind_metrics_base$auc))
  expect_true(is.numeric(ind_metrics_base$fit))
  
  # Test compute_metrics (all metrics computed on test data)
  ind_metrics <- ind$compute_metrics(test_data)
  expect_true(!is.null(ind_metrics))
  expect_true(is.numeric(ind_metrics$auc))
  expect_true(is.numeric(ind_metrics$accuracy))
  expect_true(is.numeric(ind_metrics$mcc))
})

test_that("Individual: predictions", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  ind <- exp$individual(1, 1)
  train_data <- exp$train_data()
  
  # Test predict() - returns list with class and score
  pred_result <- ind$predict(train_data)
  expect_true(is.list(pred_result))
  expect_true("class" %in% names(pred_result))
  expect_true("score" %in% names(pred_result))
  expect_true(is.numeric(pred_result$class))
  expect_true(is.numeric(pred_result$score))
  expect_true(all(pred_result$class %in% c(0, 1)))
  expect_true(length(pred_result$class) == length(pred_result$score))
  
  # Verify prediction length matches number of samples
  td <- train_data$get()
  n_samples <- ncol(td$X)
  if (!is.null(n_samples)) {
    expect_equal(length(pred_result$class), n_samples)
    expect_equal(length(pred_result$score), n_samples)
  }
})

test_that("Individual: set_threshold modifies threshold", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  ind <- exp$individual(1, 1)
  
  # Get initial threshold
  initial_threshold <- ind$get()$threshold
  
  # Modify threshold
  new_threshold <- 0.75
  new_ind <- ind$set_threshold(new_threshold)
  
  # Verify threshold has changed
  updated_threshold <- new_ind$get()$threshold
  expect_equal(updated_threshold, new_threshold)
  expect_true(updated_threshold != initial_threshold)
})

test_that("Individual: to_string, address, print", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  ind <- exp$individual(1, 1)
  
  # Test to_string
  str_repr <- ind$to_string()
  expect_true(is.character(str_repr))
  expect_true(nchar(str_repr) > 0)
  
  # Test address
  addr <- ind$address()
  expect_true(is.character(addr))
  expect_true(grepl("^0x", addr))
  
  # Test print (produces output, so we expect output rather than silence)
  expect_output(ind$print())
})

# ==============================================================================
# 5. Tests for Experiment
# ==============================================================================

test_that("Experiment: creation and basic methods", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: access to training data
  train_data <- exp$train_data()
  expect_true(!is.null(train_data))
  
  # Test: generation number and population size
  gen_num <- exp$generation_number()
  expect_true(is.numeric(gen_num))
  expect_true(gen_num > 0)
  
  pop_size <- exp$population_size(1)
  expect_true(is.numeric(pop_size))
  expect_true(pop_size > 0)
})

test_that("Experiment: get_population", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: population of a specific generation (generation 0)
  gen0_pop <- exp$get_population(1)
  expect_true(!is.null(gen0_pop))
  
  # Test: population of the last generation
  last_gen <- exp$generation_number() - 1
  final_pop <- exp$get_population(last_gen)
  expect_true(!is.null(final_pop))
})

test_that("Experiment: get_param", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: retrieving parameters
  exp_param <- exp$get_param()
  expect_true(!is.null(exp_param))
  param_data <- exp_param$get()
  expect_true(is.list(param_data))
})

test_that("Experiment: address and print", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test address
  addr <- exp$address()
  expect_true(is.character(addr))
  expect_true(grepl("^0x", addr))
  
  # Test print
  print_out <- exp$print()
  expect_true(is.character(print_out))
  expect_true(nchar(print_out) > 0)
})

# ==============================================================================
# 6. Tests for Population
# ==============================================================================

test_that("Population: get returns individuals", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set("population_size", 1500)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test: get returns a list with individuals
  pop_data <- pop$get()
  expect_true(is.list(pop_data))
  expect_true("individuals" %in% names(pop_data))
  expect_true(length(pop_data$individuals) > 0)
})

test_that("Population: filters (auc, fit, sensitivity, specificity)", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set("population_size", 2000)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test filter_by_auc
  filtered_auc <- pop$filter_by_auc(0.5)
  expect_true(!is.null(filtered_auc))
  # Verify that individuals returned have auc >= 0.5 when available
  auc_inds <- filtered_auc$get()$individuals
  if (length(auc_inds) > 0) {
    extract_metric <- function(ind, name) {
      if (is.list(ind) && !is.null(ind[[name]])) return(ind[[name]])
      if (!is.null(ind$get)) {
        gv <- tryCatch(ind$get(), error = function(e) NULL)
        if (!is.null(gv) && !is.null(gv[[name]])) return(gv[[name]])
      }
      return(NA)
    }
    vals <- vapply(auc_inds, function(i) extract_metric(i, "auc"), FUN.VALUE = numeric(1))
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) expect_true(all(vals >= 0.5))
  }
  
  # Test filter_by_fit
  filtered_fit <- pop$filter_by_fit(0.6)
  expect_true(!is.null(filtered_fit))
  fit_inds <- filtered_fit$get()$individuals
  if (length(fit_inds) > 0) {
    vals <- vapply(fit_inds, function(i) extract_metric(i, "fit"), FUN.VALUE = numeric(1))
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) expect_true(all(vals >= 0.6))
  }
  
  # Test filter_by_sensitivity
  filtered_sens <- pop$filter_by_sensitivity(0.4)
  expect_true(!is.null(filtered_sens))
  sens_inds <- filtered_sens$get()$individuals
  if (length(sens_inds) > 0) {
    vals <- vapply(sens_inds, function(i) extract_metric(i, "sensitivity"), FUN.VALUE = numeric(1))
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) expect_true(all(vals >= 0.4))
  }
  
  # Test filter_by_specificity
  filtered_spec <- pop$filter_by_specificity(0.4)
  expect_true(!is.null(filtered_spec))
  spec_inds <- filtered_spec$get()$individuals
  if (length(spec_inds) > 0) {
    vals <- vapply(spec_inds, function(i) extract_metric(i, "specificity"), FUN.VALUE = numeric(1))
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) expect_true(all(vals >= 0.4))
  }
})

test_that("Population: filter_by_k", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")

  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test: filtering by number of features
  filtered_k <- pop$filter_by_k(min_k = 2, max_k = 10)
  expect_true(!is.null(filtered_k))
})

test_that("Population: get_fbm (Family of Best Models)", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 4)
  param$set("population_size", 2000)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test get_fbm
  fbm <- pop$get_fbm(alpha = 0.05)
  expect_true(!is.null(fbm))
  
  # FBM should be smaller or equal to original population
  fbm_data <- fbm$get()
  pop_data <- pop$get()
  expect_true(length(fbm_data$individuals) <= length(pop_data$individuals))
})

test_that("Population: predict_score_matrix and predict_class_matrix", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set("population_size", 1000)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  train_data <- exp$train_data()
  
  # Test predict_score_matrix
  scores_matrix <- pop$predict_score_matrix(train_data)
  expect_true(is.data.frame(scores_matrix))
  expect_true(nrow(scores_matrix) > 0)
  expect_true(ncol(scores_matrix) > 0)
  # Try to verify orientation: number of columns or rows should match number of individuals
  pop_count <- length(pop$get()$individuals)
  expect_true((ncol(scores_matrix) == pop_count) || (nrow(scores_matrix) == pop_count))
  
  # Test predict_class_matrix
  classes_matrix <- pop$predict_class_matrix(train_data)
  expect_true(is.data.frame(classes_matrix))
  # Values should be 0/1 or NA
  flat_vals <- unlist(classes_matrix)
  flat_vals <- flat_vals[!is.na(flat_vals)]
  if (length(flat_vals) > 0) expect_true(all(flat_vals %in% c(0, 1)))
  pop_count <- length(pop$get()$individuals)
  expect_true((ncol(classes_matrix) == pop_count) || (nrow(classes_matrix) == pop_count))
})

test_that("Importance: per-individual vs population matrix consistency", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))

  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set("population_size", 100)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")

  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  pop <- pop$get_fbm(alpha = 0.05) 
  train_data <- exp$train_data()

  # Use a small number of permutations for speed in tests
  n_perm <- 50
  seed_val <- 20251105
  used_only <- FALSE

  # Compute the full importance matrix: rows = features, cols = individuals
  imp_mat <- pop$compute_importance_matrix(train_data, n_perm = n_perm, used_only = used_only, seed = seed_val)
  expect_true(!is.null(imp_mat))
  # Ensure we have rownames and at least one column
  rn <- rownames(imp_mat)
  expect_true(!is.null(rn))
  expect_true(ncol(imp_mat) > 0)

  # Compute importance for the first individual using the same seed and params
  ind <- pop$get_individual(1)
  expect_true(!is.null(ind))
  ind_imp_raw <- ind$compute_importance(train_data, n_perm = n_perm, seed = seed_val, used_only = used_only)
  expect_true(!is.null(ind_imp_raw))

  # Align individual importance to matrix row order: fill NA for missing features
  mat_feat_names <- rn
  ind_vec <- rep(NA_real_, length(mat_feat_names))
  names(ind_vec) <- mat_feat_names
  if (!is.null(names(ind_imp_raw))) {
    common <- intersect(names(ind_imp_raw), mat_feat_names)
    ind_vec[common] <- as.numeric(ind_imp_raw[common])
  } else {
    # if unnamed, attempt to coerce by length (best effort)
    tmp <- as.numeric(ind_imp_raw)
    if (length(tmp) == length(ind_vec)) ind_vec[] <- tmp
  }

  # Extract the first individual's column from the matrix and compare
  mat_col <- as.numeric(imp_mat[, 1])
  names(mat_col) <- rownames(imp_mat)

  # Both should be numeric vectors of same length
  expect_equal(length(mat_col), length(ind_vec))
  # Compare values; allow a small numerical tolerance
  expect_equal(mat_col, ind_vec, tolerance = 1e-8)
})

test_that("Population: get_individual", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test: extracting an individual
  ind <- pop$get_individual(1)
  expect_true(!is.null(ind))
  ind_data <- ind$get()
  expect_true(is.list(ind_data))
  expect_true("features" %in% names(ind_data))
})

test_that("Population: address and print", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test address
  addr <- pop$address()
  expect_true(is.character(addr))
  expect_true(grepl("^0x", addr))
  
  # Test print - returns a string
  print_out <- pop$print()
  expect_true(is.character(print_out))
  expect_true(nchar(print_out) > 0)
})

# ==============================================================================
# 7. Tests for Jury
# ==============================================================================

test_that("Jury: creation from population and param", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set("population_size", 1500)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test: creating a Jury
  jury <- Jury$new_from_param(pop, param)
  expect_true(!is.null(jury))
})

test_that("Jury: get returns jury information", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  jury <- Jury$new_from_param(pop, param)
  
  # Test get
  jury_data <- jury$get()
  expect_true(is.list(jury_data))
  expect_true("experts" %in% names(jury_data))
  expect_true("voting_method" %in% names(jury_data))
  expect_true("auc" %in% names(jury_data))
})

test_that("Jury: predict", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  jury <- Jury$new_from_param(pop, param)
  train_data <- exp$train_data()
  
  # Test predictions - fit() instead of evaluate()
  jury$fit(train_data)
  predictions <- jury$predict(train_data)
  expect_true(is.list(predictions))
  expect_true("class" %in% names(predictions))
  expect_true("score" %in% names(predictions))
  # Verify lengths equal number of samples
  td <- train_data$get()
  n_samples <- ncol(td$X)
  if (!is.null(n_samples)) {
    expect_equal(length(predictions$class), n_samples)
    expect_equal(length(predictions$score), n_samples)
  }
})

test_that("Jury: get_population", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  jury <- Jury$new_from_param(pop, param)
  
  # Test extracting population
  jury_pop <- jury$get_population()
  expect_true(!is.null(jury_pop))
})

test_that("Jury: address and print", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  jury <- Jury$new_from_param(pop, param)
  
  # Test address
  addr <- jury$address()
  expect_true(is.character(addr))
  expect_true(grepl("^0x", addr))
  
  # Test print - returns a string
  print_out <- jury$print()
  expect_true(is.character(print_out))
  expect_true(nchar(print_out) > 0)
})

# ==============================================================================
# 8. Tests for GLogger
# ==============================================================================

test_that("GLogger: creation with new", {
  logger <- GLogger$new()
  expect_true(!is.null(logger))
})

test_that("GLogger: creation with specific level", {
  logger_info <- GLogger$level("info")
  expect_true(!is.null(logger_info))
  
  logger_debug <- GLogger$level("debug")
  expect_true(!is.null(logger_debug))
  
  logger_error <- GLogger$level("error")
  expect_true(!is.null(logger_error))
})

test_that("GLogger: set_level modifies level", {
  logger <- GLogger$new()
  
  # Test level changes
  expect_silent(logger$set_level("debug"))
  expect_silent(logger$set_level("warn"))
  expect_silent(logger$set_level("error"))
})

test_that("GLogger: creation from Param", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  
  param <- Param$load("../../sample/param.yaml")
  logger <- GLogger$get(param)
  
  expect_true(!is.null(logger))
})

# ==============================================================================
# 9. Integration tests
# ==============================================================================

test_that("Integration: complete simple workflow", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  # 1. Create a logger
  logger <- GLogger$level("info")
  
  # 2. Load parameters
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set("population_size", 1000)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  # 3. Create a flag
  flag <- RunningFlag$new()
  
  # 4. Run experiment
  exp <- fit(param, flag)
  
  # 5. Extract final population
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # 6. Filter by AUC
  filtered_pop <- pop$filter_by_auc(0.5)
  
  # 7. Create a jury
  jury <- Jury$new_from_param(filtered_pop, param)
  
  # 8. Make predictions - use fit() instead of evaluate()
  train_data <- exp$train_data()
  jury$fit(train_data)
  predictions <- jury$predict(train_data)
  
  # Verifications
  expect_true(!is.null(predictions))
  expect_true(is.list(predictions))
  expect_true(length(predictions$class) > 0)
})

test_that("Integration: FBM and feature importance", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 4)
  param$set("population_size", 2000)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Extract FBM
  fbm <- pop$get_fbm(alpha = 0.05)
  expect_true(!is.null(fbm))
  
  # Compute importance (if method is available)
  # Note: compute_importance may require specific parameters
  train_data <- exp$train_data()
  
  # Test with an individual
  ind <- fbm$get_individual(1)
  importance <- tryCatch({
    ind$compute_importance(train_data, n_perm = 100, seed = 42, used_only = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(importance)) {
    expect_true(is.numeric(importance) || is.data.frame(importance))
  }
})

# ======================================================================
# 10. Tests for serialization: JSON, MSPack (.mp) and BIN
# These tests attempt to save an Experiment to different formats and reload
# The backend may not support all formats; tests will skip if a format is
# not supported (caught as an error during save/load).
# ======================================================================

test_that("Experiment: serialization JSON", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))

  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")

  exp <- fit(param, running_flag)

  temp_file <- tempfile(fileext = ".json")

  # Save as JSON and assert success (do not skip silently)
  exp$save(temp_file)
  expect_true(file.exists(temp_file))
  expect_true(file.info(temp_file)$size > 0)

  # Load back and compare
  exp_loaded <- Experiment$load(temp_file)
  expect_equal(exp$generation_number(), exp_loaded$generation_number())
  expect_equal(exp$population_size(1), exp_loaded$population_size(1))

  unlink(temp_file)
})

test_that("Experiment: serialization MSPack (.mp)", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))

  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")

  exp <- fit(param, running_flag)

  temp_file <- tempfile(fileext = ".mp")

  # Save as MSPack and assert success (do not skip silently)
  exp$save(temp_file)
  expect_true(file.exists(temp_file))
  expect_true(file.info(temp_file)$size > 0)

  # Load back and compare
  exp_loaded <- Experiment$load(temp_file)
  expect_equal(exp$generation_number(), exp_loaded$generation_number())
  expect_equal(exp$population_size(1), exp_loaded$population_size(1))

  unlink(temp_file)
})

test_that("Experiment: serialization BIN (.bin)", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))

  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")

  exp <- fit(param, running_flag)

  temp_file <- tempfile(fileext = ".bin")

  # Save as BIN and assert success (do not skip silently)
  exp$save(temp_file)
  expect_true(file.exists(temp_file))
  expect_true(file.info(temp_file)$size > 0)

  # Load back and compare
  exp_loaded <- Experiment$load(temp_file)
  expect_equal(exp$generation_number(), exp_loaded$generation_number())
  expect_equal(exp$population_size(1), exp_loaded$population_size(1))

  unlink(temp_file)
})

# ======================================================================
# Test: Verify Population matrix predictions vs individual predictions
# ======================================================================
test_that("Population: predictions consistent between matrix and individuals", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))

  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")

  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  train_data <- exp$train_data()

  # Predictions for the entire population
  pop_predictions <- pop$predict_class_matrix(train_data)

  # Verify predictions for each individual
  for (i in seq_len(length(pop$get()$individuals))) {
    ind <- pop$get_individual(i)
    ind_predictions <- ind$predict(train_data)$class

    # Compare individual predictions with the corresponding column
    expect_equal(as.numeric(pop_predictions[, i]), ind_predictions)
  }
})


# ======================================================================
# Test: Mathematical verification of metrics for Jury and Individual
# ======================================================================
test_that("Mathematical verification of metrics for Jury and Individual", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))

  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")

  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  train_data <- exp$train_data()

  ind <- pop$get_individual(1)
  ind_predictions_result <- ind$predict(train_data)
  ind_metrics <- ind$compute_metrics(train_data)

  true_labels <- train_data$get()$y
  predicted_labels <- ind_predictions_result$class

  tp <- sum(true_labels == 1 & predicted_labels == 1)
  tn <- sum(true_labels == 0 & predicted_labels == 0)
  fp <- sum(true_labels == 0 & predicted_labels == 1)
  fn <- sum(true_labels == 1 & predicted_labels == 0)

  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)

  expect_equal(ind_metrics$accuracy, accuracy, tolerance = 1e-8)
  expect_equal(ind_metrics$sensitivity, sensitivity, tolerance = 1e-8)
  expect_equal(ind_metrics$specificity, specificity, tolerance = 1e-8)

  jury <- Jury$new_from_param(pop, param)
  jury$fit(train_data)
  jury_metrics <- jury$get()

  jury_predictions <- jury$predict(train_data)$class

  tp_jury <- sum(true_labels == 1 & jury_predictions == 1)
  tn_jury <- sum(true_labels == 0 & jury_predictions == 0)
  fp_jury <- sum(true_labels == 0 & jury_predictions == 1)
  fn_jury <- sum(true_labels == 1 & jury_predictions == 0)

  accuracy_jury <- (tp_jury + tn_jury) / (tp_jury + tn_jury + fp_jury + fn_jury)
  sensitivity_jury <- tp_jury / (tp_jury + fn_jury)
  specificity_jury <- tn_jury / (tn_jury + fp_jury)

  expect_equal(jury_metrics$accuracy, accuracy_jury, tolerance = 1e-8)
  expect_equal(jury_metrics$sensitivity, sensitivity_jury, tolerance = 1e-8)
  expect_equal(jury_metrics$specificity, specificity_jury, tolerance = 1e-8)
})

test_that("Population: creation from Individuals vector and data consistency", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))

  # Load parameters
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set("population_size", 1000)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")

  # Create a small experiment
  running_flag <- RunningFlag$new()
  exp <- fit(param, running_flag)

  # Retrieve individuals from the first generation
  individuals <- c(exp$get_population(1)$get_individual(1),
                   exp$get_population(1)$get_individual(1),
                   exp$get_population(1)$get_individual(2))

  # Create a new Population from these individuals
  new_population <- Population$from_individuals(individuals)

  # Test: Population object is created
  expect_true(!is.null(new_population))

  # Test: Data consistency in the new Population
  pop_data <- new_population$get()
  expect_true(is.list(pop_data))
  expect_true("individuals" %in% names(pop_data))
  expect_equal(length(pop_data$individuals), length(individuals))

  # Verify that the data of each individual matches
  for (i in seq_along(individuals)) {
    original_data <- individuals[[i]]$get()
    new_data <- pop_data$individuals[[i]]
    expect_equal(original_data, new_data)
  }
})

# ==============================================================================
# 11. Tests for Cross-Validation (CV) functionality
# ==============================================================================

test_that("CV: experiment with CV enabled creates fold structure", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set("population_size", 500)
  param$set_bool("cv", TRUE)
  param$set("outer_folds", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: CV folds are created
  n_folds <- exp$get_n_folds()
  expect_equal(n_folds, 3)
  
  # Test: Best population is available
  best_pop <- exp$get_best_population()
  expect_true(!is.null(best_pop))
  expect_true(length(best_pop$get()$individuals) > 0)
})

test_that("CV: get_fold_data returns correct train/validation data", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_bool("cv", TRUE)
  param$set("outer_folds", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: Get training data for fold 1
  fold0_train <- exp$get_fold_data(fold = 1, train = TRUE)
  expect_true(!is.null(fold0_train))
  expect_true(is.list(fold0_train$get()))
  
  # Test: Get validation data for fold 1
  fold0_valid <- exp$get_fold_data(fold = 1, train = FALSE)
  expect_true(!is.null(fold0_valid))
  expect_true(is.list(fold0_valid$get()))
  
  # Test: Train and validation data have different sizes
  train_samples <- ncol(fold0_train$get()$X)
  valid_samples <- ncol(fold0_valid$get()$X)
  expect_true(train_samples != valid_samples)
  
  # Test: Total samples should equal original data
  total_data <- exp$train_data()
  total_samples <- ncol(total_data$get()$X)
  expect_equal(train_samples + valid_samples, total_samples)
})

test_that("CV: get_fold_generation retrieves correct population", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_bool("cv", TRUE)
  param$set("outer_folds", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: Get population from fold 1, generation 1
  fold0_gen0 <- exp$get_fold_generation(fold = 1, generation = 1, train=TRUE)
  expect_true(!is.null(fold0_gen0))
  expect_true(length(fold0_gen0$get()$individuals) > 0)
  
  # Test: Get population from fold 2, last generation
  last_gen <- as.integer(exp$generation_number()[1])
  fold1_last <- exp$get_fold_generation(fold = 2, generation = last_gen, train=TRUE)
  expect_true(!is.null(fold1_last))
  
  # Test: Population sizes should be consistent
  expect_true(length(fold0_gen0$get()$individuals) > 0)
  expect_true(length(fold1_last$get()$individuals) > 0)
})

test_that("CV: multiple folds have independent populations", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_bool("cv", TRUE)
  param$set("outer_folds", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  
  # Get populations from different folds
  fold0_pop <- exp$get_fold_generation(fold = 1, generation = 1, train=TRUE)
  fold1_pop <- exp$get_fold_generation(fold = 2, generation = 1, train=TRUE)
  fold2_pop <- exp$get_fold_generation(fold = 3, generation = 1, train=TRUE)
  
  # Test: All populations exist
  expect_true(!is.null(fold0_pop))
  expect_true(!is.null(fold1_pop))
  expect_true(!is.null(fold2_pop))
  
  # Test: Populations have different individuals (by hash)
  fold0_hashes <- sapply(fold0_pop$get()$individuals, function(ind) ind$hash)
  fold1_hashes <- sapply(fold1_pop$get()$individuals, function(ind) ind$hash)
  
  # At least some individuals should be different across folds
  expect_false(all(fold0_hashes %in% fold1_hashes))
})

test_that("CV: compute_cv_importance_matrix returns valid importance matrix", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set("population_size", 500)
  param$set_bool("cv", TRUE)
  param$set("outer_folds", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: Compute CV importance matrix
  imp_matrix <- exp$compute_cv_importance_matrix(
    n_perm = 50,
    used_only = TRUE,
    seed = 42,
    aggregation = "mean",
    scaled = TRUE
  )
  
  # Test: Matrix structure
  expect_true(is.matrix(imp_matrix) || is.data.frame(imp_matrix))
  expect_true(nrow(imp_matrix) > 0)
  expect_equal(ncol(imp_matrix), 2)  # 2 folds
  
  # Test: Column names are fold identifiers
  expect_true(all(grepl("^Fold_", colnames(imp_matrix))))
  
  # Test: Row names are feature names
  expect_true(!is.null(rownames(imp_matrix)))
  expect_true(length(rownames(imp_matrix)) > 0)
  
  # Test: Values are numeric
  expect_true(is.numeric(imp_matrix[1,1]))
})

test_that("CV: error handling for invalid fold access", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_bool("cv", TRUE)
  param$set("outer_folds", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: Invalid fold index should error
  expect_error(exp$get_fold_data(fold = 10, train = TRUE))
  expect_error(exp$get_fold_generation(fold = 10, generation = 0))
  
  # Test: Invalid generation should error
  expect_error(exp$get_fold_generation(fold = 0, generation = 999))
})

test_that("CV: experiment without CV fails fold operations", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_bool("cv", FALSE)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  
  # Test: get_n_folds should return 1
  n_folds <- exp$get_n_folds()
  expect_equal(n_folds, 1)
  
  # Test: Fold operations should fail
  expect_error(exp$get_fold_data(fold = 0, train = TRUE))
})

# ==============================================================================
# 12. Tests for Pruning (Feature Importance-based)
# ==============================================================================

test_that("Pruning: Individual prune_by_threshold reduces features", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Get an individual with multiple features
  ind <- pop$get_individual(1)
  original_k <- ind$get()$k
  
  # Skip test if individual has only 1 feature
  skip_if(original_k <= 1, "Individual has only 1 feature, cannot prune")
  
  # Test: Prune by threshold
  pruned_ind <- ind$prune_by_threshold(
    threshold = 0.01,
    n_perm = 50,
    seed = 42,
    min_k = 1
  )
  
  pruned_k <- pruned_ind$get()$k
  
  # Test: Pruned individual should have fewer or equal features
  expect_true(pruned_k <= original_k)
  
  # Test: At least min_k features should remain
  expect_true(pruned_k >= 1)
  
  # Test: Metrics are recomputed
  expect_true(is.numeric(pruned_ind$get()$auc))
  expect_true(is.numeric(pruned_ind$get()$fit))
})

test_that("Pruning: Individual prune_by_quantile reduces features", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  ind <- pop$get_individual(1)
  original_k <- ind$get()$k
  
  skip_if(original_k <= 2, "Individual has too few features, cannot prune")
  
  # Test: Prune by quantile (25th percentile)
  pruned_ind <- ind$prune_by_quantile(
    quantile = 0.25,
    eps = 0.0,
    n_perm = 50,
    seed = 42,
    min_k = 1
  )
  
  pruned_k <- pruned_ind$get()$k
  
  # Test: Pruned individual should have fewer or equal features
  expect_true(pruned_k <= original_k)
  expect_true(pruned_k >= 1)
  
  # Test: Metrics are valid
  expect_true(!is.na(pruned_ind$get()$auc))
})

test_that("Pruning: min_k parameter enforces minimum features", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  ind <- pop$get_individual(1)
  original_k <- ind$get()$k
  
  skip_if(original_k <= 3, "Individual has too few features")
  
  # Test: Prune with high threshold but enforce min_k = 3
  pruned_ind <- ind$prune_by_threshold(
    threshold = 1000.0,  # Very high threshold
    n_perm = 50,
    seed = 42,
    min_k = 3
  )
  
  pruned_k <- pruned_ind$get()$k
  
  # Test: At least min_k features should remain
  expect_true(pruned_k >= 3)
})

test_that("Pruning: Population prune_by_threshold prunes all individuals", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set("population_size", 200)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Take a subset of population
  pop_subset <- pop$filter_by_auc(0.5)
  n_individuals <- length(pop_subset$get()$individuals)
  
  skip_if(n_individuals == 0, "No individuals passed AUC filter")
  
  # Test: Prune entire population
  pruned_pop <- pop_subset$prune_by_threshold(
    threshold = 0.01,
    n_perm = 30,
    seed = 42,
    min_k = 1
  )
  
  # Test: Population still has same number of individuals
  expect_equal(length(pruned_pop$get()$individuals), n_individuals)
  
  # Test: All individuals have valid metrics
  for (i in seq_len(n_individuals)) {
    ind <- pruned_pop$get_individual(i)
    expect_true(is.numeric(ind$get()$auc))
    expect_true(is.numeric(ind$get()$fit))
    expect_true(ind$get()$k >= 1)
  }
})

test_that("Pruning: Population prune_by_quantile prunes all individuals", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set("population_size", 200)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Take a subset
  pop_subset <- pop$filter_by_auc(0.5)
  n_individuals <- length(pop_subset$get()$individuals)
  
  skip_if(n_individuals == 0, "No individuals passed filter")
  
  # Test: Prune by quantile
  pruned_pop <- pop_subset$prune_by_quantile(
    quantile = 0.30,
    eps = 0.05,
    n_perm = 30,
    seed = 42,
    min_k = 1
  )
  
  # Test: Population integrity
  expect_equal(length(pruned_pop$get()$individuals), n_individuals)
  
  # Test: All individuals have valid metrics
  aucs <- sapply(pruned_pop$get()$individuals, function(ind) ind$auc)
  expect_true(all(!is.na(aucs)))
  expect_true(all(aucs >= 0))
  expect_true(all(aucs <= 1))
})

test_that("Pruning: pruned features are deterministic with same seed", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  ind <- pop$get_individual(1)
  
  skip_if(ind$get()$k <= 1, "Individual has only 1 feature")
  
  # Test: Prune twice with same seed
  pruned1 <- ind$prune_by_threshold(
    threshold = 0.01,
    n_perm = 50,
    seed = 123,
    min_k = 1
  )
  
  pruned2 <- ind$prune_by_threshold(
    threshold = 0.01,
    n_perm = 50,
    seed = 123,
    min_k = 1
  )
  
  # Test: Results should be identical
  expect_equal(pruned1$get()$k, pruned2$get()$k)
  expect_equal(pruned1$get()$features, pruned2$get()$features)
})

test_that("Pruning: different seeds produce potentially different results", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  ind <- pop$get_individual(1)
  
  skip_if(ind$get()$k <= 2, "Individual has too few features")
  
  # Prune with different seeds
  pruned1 <- ind$prune_by_threshold(
    threshold = 0.01,
    n_perm = 50,
    seed = 111,
    min_k = 1
  )
  
  pruned2 <- ind$prune_by_threshold(
    threshold = 0.01,
    n_perm = 50,
    seed = 222,
    min_k = 1
  )
  
  # Test: Results may differ (not guaranteed, but likely)
  # At minimum, they should both be valid
  expect_true(is.numeric(pruned1$get()$auc))
  expect_true(is.numeric(pruned2$get()$auc))
  expect_true(pruned1$get()$k >= 1)
  expect_true(pruned2$get()$k >= 1)
})

test_that("Pruning: importance calculation consistency", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  ind <- pop$get_individual(1)
  train_data <- exp$train_data()
  
  # Test: Compute importance manually
  imp <- ind$compute_importance(
    data = train_data,
    n_perm = 50,
    seed = 42,
    used_only = TRUE
  )
  
  # Test: Importance values are numeric
  expect_true(is.numeric(imp))
  expect_true(all(!is.na(imp)))
  
  # Test: Feature names match
  expect_true(!is.null(names(imp)))
  expect_true(all(names(imp) %in% train_data$get()$features))
})

test_that("Pruning: edge case with single feature individual", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set("kmin", 1)
  param$set("kmax", 1)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Find an individual with k=1
  single_feature_ind <- NULL
  for (i in seq_len(min(10, length(pop$get()$individuals)))) {
    ind <- pop$get_individual(i)
    if (ind$get()$k == 1) {
      single_feature_ind <- ind
      break
    }
  }
  
  skip_if(is.null(single_feature_ind), "No single-feature individual found")
  
  # Test: Pruning should preserve the single feature
  pruned <- single_feature_ind$prune_by_threshold(
    threshold = 0.01,
    n_perm = 30,
    seed = 42,
    min_k = 1
  )
  
  expect_equal(pruned$get()$k, 1)
  expect_equal(pruned$get()$features, single_feature_ind$get()$features)
})

# ==============================================================================
# 13. Additional tests for missing coverage
# ==============================================================================

# ------------------------------------------------------------------------------
# Tests for Param
# ------------------------------------------------------------------------------

test_that("Param: creation with new()", {
  param <- Param$new()
  
  # Test: Param object should be created
  expect_true(!is.null(param))
  
  # Test: get() should return a list
  param_data <- param$get()
  expect_true(is.list(param_data))
})

# ------------------------------------------------------------------------------
# Tests for Data
# ------------------------------------------------------------------------------

test_that("Data: creation with new()", {
  data <- Data$new()
  
  # Test: Data object should be created
  expect_true(!is.null(data))
  
  # Test: get() should return a list
  data_content <- data$get()
  expect_true(is.list(data_content))
})

test_that("Data: train_test_split with stratify_by", {
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  # Load data
  X <- read.table("../../sample/Xtrain.tsv", h=TRUE, row.names=1)
  colnames(X) <- gsub("\\.", "-", colnames(X), perl = TRUE)
  
  y_df <- read.table("../../sample/Ytrain.tsv", h=TRUE)
  y <- factor(ifelse(y_df[,1]==0, "healthy", "cirrhosis"),
              levels = c("healthy", "cirrhosis"),
              ordered = TRUE)
  names(y) <- rownames(y_df)
  
  # Create sample annotations with stratification variable
  sample_tags <- data.frame(
    sample = names(y),
    group = sample(c("A", "B"), length(y), replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  data <- as_gpredomics_data(X, y, FALSE, 
                            prior_weight = NULL, 
                            feature_penalty = NULL, 
                            feature_tags = NULL, 
                            sample_tags = sample_tags)
  
  # Test: Split with stratification
  split_result <- data$train_test_split(test_ratio = 0.3, stratify_by = "group", seed = 42)
  
  expect_true(is.list(split_result))
  expect_true("train" %in% names(split_result))
  expect_true("test" %in% names(split_result))
  
  train_size <- ncol(split_result$train$get()$X)
  test_size <- ncol(split_result$test$get()$X)
  
  expect_true(test_size > 0)
  expect_true(train_size > 0)
  expect_equal(train_size + test_size, ncol(X))
})

test_that("Data: train_test_split with different ratios", {
  set.seed(123)
  X <- matrix(runif(100 * 20), nrow = 100, ncol = 20)
  rownames(X) <- paste0("Feature_", 1:100)
  colnames(X) <- paste0("Sample_", 1:20)
  
  y <- factor(sample(c("control", "case"), 20, replace = TRUE),
              levels = c("control", "case"))
  names(y) <- colnames(X)
  
  data <- as_gpredomics_data(as.data.frame(X), y, FALSE, NULL, NULL, NULL, NULL)
  
  # Test with ratio 0.2
  split1 <- data$train_test_split(test_ratio = 0.2, stratify_by = NULL, seed = 42)
  test_size1 <- ncol(split1$test$get()$X)
  expect_equal(test_size1, 4)  # 20 * 0.2 = 4
  
  # Test with ratio 0.5
  split2 <- data$train_test_split(test_ratio = 0.5, stratify_by = NULL, seed = 42)
  test_size2 <- ncol(split2$test$get()$X)
  # Allow for rounding: 0.5 * 20 = 10, but may be 10 or 11
  expect_true(test_size2 >= 9 && test_size2 <= 11)
})

test_that("Data: train_test_split error handling", {
  set.seed(456)
  X <- matrix(runif(50 * 10), nrow = 50, ncol = 10)
  rownames(X) <- paste0("Feature_", 1:50)
  colnames(X) <- paste0("Sample_", 1:10)
  
  y <- factor(sample(0:1, 10, replace = TRUE))
  names(y) <- colnames(X)
  
  data <- as_gpredomics_data(as.data.frame(X), y, FALSE, NULL, NULL, NULL, NULL)
  
  # Test: Invalid ratio (too low)
  expect_error(data$train_test_split(test_ratio = 0.0, stratify_by = NULL))
  
  # Test: Invalid ratio (too high)
  expect_error(data$train_test_split(test_ratio = 1.0, stratify_by = NULL))
  
  # Test: Invalid ratio (negative)
  expect_error(data$train_test_split(test_ratio = -0.1, stratify_by = NULL))
})

# ------------------------------------------------------------------------------
# Tests for Individual
# ------------------------------------------------------------------------------

test_that("Individual: evaluate method", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  ind <- pop$get_individual(1)
  
  # Test: evaluate returns numeric scores
  scores <- ind$evaluate()
  expect_true(is.numeric(scores))
  expect_true(length(scores) > 0)
})

# ------------------------------------------------------------------------------
# Tests for Experiment
# ------------------------------------------------------------------------------

test_that("Experiment: get_history with different scopes", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 4)
  param$set("population_size", 1000)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  train_data <- exp$train_data()
  
  # Test: History with scope "best"
  hist_best <- exp$get_history(train_data, "best")
  expect_true(is.data.frame(hist_best))
  expect_true("Generation" %in% names(hist_best))
  expect_true("AUC" %in% names(hist_best))
  expect_equal(nrow(hist_best), 4)
  
  # Test: History with scope "fbm"
  hist_fbm <- exp$get_history(train_data, "fbm")
  expect_true(is.data.frame(hist_fbm))
  expect_equal(nrow(hist_fbm), 4)
  
  # Test: History with scope "top5"
  hist_top5 <- exp$get_history(train_data, "top5")
  expect_true(is.data.frame(hist_top5))
  expect_equal(nrow(hist_top5), 4)
  
  # Test: History with scope "all"
  hist_all <- exp$get_history(train_data, "all")
  expect_true(is.data.frame(hist_all))
  expect_equal(nrow(hist_all), 4)
})

test_that("Experiment: get_history with invalid scope", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  train_data <- exp$train_data()
  
  # Test: Invalid scope should error
  expect_error(exp$get_history(train_data, "invalid_scope"))
})

# ------------------------------------------------------------------------------
# Tests for Population
# ------------------------------------------------------------------------------

test_that("Population: display_feature_prevalence", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  train_data <- exp$train_data()
  
  # Test: Should not error (output may not be captured in all environments)
  expect_no_error(pop$display_feature_prevalence(train_data, 10))
})

test_that("Population: filter_by_diversity", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 3)
  param$set("population_size", 1500)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test: Filter by diversity (no niche)
  filtered_pop <- pop$filter_by_diversity(min_diversity_pct = 50.0, by_niche = FALSE)
  expect_true(!is.null(filtered_pop))
  expect_true(length(filtered_pop$get()$individuals) > 0)
  
  # Test: Filtered population should be smaller or equal
  expect_true(length(filtered_pop$get()$individuals) <= length(pop$get()$individuals))
  
  # Test: Filter by diversity (with niche)
  filtered_niche <- pop$filter_by_diversity(min_diversity_pct = 30.0, by_niche = TRUE)
  expect_true(!is.null(filtered_niche))
  expect_true(length(filtered_niche$get()$individuals) > 0)
})

test_that("Population: filter_by_mask", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set("population_size", 500)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  pop_size <- length(pop$get()$individuals)
  
  # Create a mask: keep first half (must be integers)
  mask <- as.integer(c(rep(1, ceiling(pop_size/2)), rep(0, floor(pop_size/2))))
  
  # Test: Filter with mask
  filtered <- pop$filter_by_mask(mask)
  expect_true(!is.null(filtered))
  expect_equal(length(filtered$get()$individuals), ceiling(pop_size/2))
  
  # Test: Error with wrong mask length
  wrong_mask <- c(1, 0, 1)
  expect_error(pop$filter_by_mask(wrong_mask))
})

test_that("Population: get_first_pct explicit test", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set("population_size", 1000)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test: Get first 10%
  first_10 <- pop$get_first_pct(10.0)
  expect_true(!is.null(first_10))
  
  pop_size <- length(pop$get()$individuals)
  first_size <- length(first_10$get()$individuals)
  
  expect_true(first_size <= pop_size)
  expect_true(first_size >= 0.09 * pop_size)  # Allow some rounding
  expect_true(first_size <= 0.11 * pop_size)
})

test_that("Population: fit method recomputes metrics", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  test_data <- exp$test_data()
  
  # Store original AUC
  original_aucs <- sapply(pop$get()$individuals, function(ind) ind$auc)
  
  # Test: Refit on test data
  refitted_pop <- pop$fit(test_data, param)
  
  # Test: Population has been refitted
  expect_true(!is.null(refitted_pop))
  expect_equal(length(refitted_pop$get()$individuals), length(pop$get()$individuals))
  
  # Metrics should be recomputed (may differ on test data)
  new_aucs <- sapply(refitted_pop$get()$individuals, function(ind) ind$auc)
  expect_true(is.numeric(new_aucs))
  expect_true(all(!is.na(new_aucs)))
})

# ------------------------------------------------------------------------------
# Tests for Jury
# ------------------------------------------------------------------------------

test_that("Jury: from_population creation", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set("population_size", 500)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  
  # Test: Create jury from population
  jury <- Jury$from_population(pop, threshold = 0.5, window = 5.0)
  
  expect_true(!is.null(jury))
  
  # Test: Jury has experts
  jury_data <- jury$get()
  expect_true(length(jury_data$experts) > 0)
  expect_equal(jury_data$voting_threshold, 0.5)
  expect_equal(jury_data$threshold_window, 5.0)
})

test_that("Jury: get_metrics and compute_metrics", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  jury <- Jury$new_from_param(pop, param)
  
  train_data <- exp$train_data()
  jury$fit(train_data)
  
  # Test: get_metrics (base metrics only)
  base_metrics <- jury$get_metrics()
  expect_true(is.list(base_metrics))
  expect_true("auc" %in% names(base_metrics))
  expect_true("accuracy" %in% names(base_metrics))
  expect_true("sensitivity" %in% names(base_metrics))
  expect_true("specificity" %in% names(base_metrics))
  expect_true("rejection_rate" %in% names(base_metrics))
  
  # Test: compute_metrics (all metrics)
  test_data <- exp$test_data()
  all_metrics <- jury$compute_metrics(test_data)
  expect_true(is.list(all_metrics))
  expect_true("mcc" %in% names(all_metrics))
  expect_true("npv" %in% names(all_metrics))
  expect_true("ppv" %in% names(all_metrics))
  expect_true("f1_score" %in% names(all_metrics))
  expect_true("g_mean" %in% names(all_metrics))
})

test_that("Jury: print_self_report and print_report", {
  skip_if_not(file.exists("../../sample/param.yaml"))
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/Xtest.tsv"))
  skip_if_not(file.exists("../../sample/Ytest.tsv"))
  
  running_flag <- RunningFlag$new()
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set_string("X", "../../sample/Xtrain.tsv")
  param$set_string("y", "../../sample/Ytrain.tsv")
  param$set_string("Xtest", "../../sample/Xtest.tsv")
  param$set_string("ytest", "../../sample/Ytest.tsv")
  
  exp <- fit(param, running_flag)
  pop <- exp$get_population(exp$generation_number() - 1)
  jury <- Jury$new_from_param(pop, param)
  
  train_data <- exp$train_data()
  jury$fit(train_data)
  
  # Test: print_self_report outputs to console
  expect_output(jury$print_self_report())
  
  # Test: print_report with test data
  test_data <- exp$test_data()
  expect_output(jury$print_report(test_data))
})

# ------------------------------------------------------------------------------
# Tests for fit_on function
# ------------------------------------------------------------------------------

test_that("fit_on: run algorithm with Data object", {
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/param.yaml"))
  
  # Load data
  X <- read.table("../../sample/Xtrain.tsv", h=TRUE, row.names=1)
  colnames(X) <- gsub("\\.", "-", colnames(X), perl = TRUE)
  
  y_df <- read.table("../../sample/Ytrain.tsv", h=TRUE)
  y <- factor(ifelse(y_df[,1]==0, "healthy", "cirrhosis"),
              levels = c("healthy", "cirrhosis"),
              ordered = TRUE)
  names(y) <- rownames(y_df)
  
  data <- as_gpredomics_data(X, y, FALSE, NULL, NULL, NULL, NULL)
  
  # Load parameters
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set("population_size", 500)
  
  running_flag <- RunningFlag$new()
  
  # Test: fit_on with Data object
  exp <- fit_on(data, param, running_flag)
  
  expect_true(!is.null(exp))
  expect_true(exp$generation_number() > 0)
  expect_true(exp$population_size(1) > 0)
  
  # Test: Final population exists
  final_pop <- exp$get_best_population()
  expect_true(!is.null(final_pop))
  expect_true(length(final_pop$get()$individuals) > 0)
})

test_that("fit_on: comparison with fit on same data", {
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/param.yaml"))
  
  # Load data
  X <- read.table("../../sample/Xtrain.tsv", h=TRUE, row.names=1)
  colnames(X) <- gsub("\\.", "-", colnames(X), perl = TRUE)
  
  y_df <- read.table("../../sample/Ytrain.tsv", h=TRUE)
  y <- factor(ifelse(y_df[,1]==0, "healthy", "cirrhosis"),
              levels = c("healthy", "cirrhosis"),
              ordered = TRUE)
  names(y) <- rownames(y_df)
  
  data <- as_gpredomics_data(X, y, FALSE, NULL, NULL, NULL, NULL)
  
  # Same parameters for both
  param1 <- Param$load("../../sample/param.yaml")
  param1$set("max_epochs", 2)
  param1$set("population_size", 300)
  param1$set("seed", 12345)
  param1$set_string("X", "../../sample/Xtrain.tsv")
  param1$set_string("y", "../../sample/Ytrain.tsv")
  
  param2 <- Param$load("../../sample/param.yaml")
  param2$set("max_epochs", 2)
  param2$set("population_size", 300)
  param2$set("seed", 12345)
  
  running_flag1 <- RunningFlag$new()
  running_flag2 <- RunningFlag$new()
  
  # Run with fit
  exp1 <- fit(param1, running_flag1)
  
  # Run with fit_on
  exp2 <- fit_on(data, param2, running_flag2)
  
  # Both should produce valid experiments
  expect_true(!is.null(exp1))
  expect_true(!is.null(exp2))
  expect_equal(exp1$generation_number(), exp2$generation_number())
})

test_that("fit_on: with voting enabled", {
  skip_if_not(file.exists("../../sample/Xtrain.tsv"))
  skip_if_not(file.exists("../../sample/Ytrain.tsv"))
  skip_if_not(file.exists("../../sample/param.yaml"))
  
  # Load data
  X <- read.table("../../sample/Xtrain.tsv", h=TRUE, row.names=1)
  colnames(X) <- gsub("\\.", "-", colnames(X), perl = TRUE)
  
  y_df <- read.table("../../sample/Ytrain.tsv", h=TRUE)
  y <- factor(ifelse(y_df[,1]==0, "healthy", "cirrhosis"),
              levels = c("healthy", "cirrhosis"),
              ordered = TRUE)
  names(y) <- rownames(y_df)
  
  data <- as_gpredomics_data(X, y, FALSE, NULL, NULL, NULL, NULL)
  
  # Enable voting
  param <- Param$load("../../sample/param.yaml")
  param$set("max_epochs", 2)
  param$set("population_size", 500)
  param$set_bool("vote", TRUE)
  
  running_flag <- RunningFlag$new()
  
  # Test: fit_on with voting
  exp <- fit_on(data, param, running_flag)
  
  expect_true(!is.null(exp))
  
  # Test: Jury should be created
  jury <- exp$get_jury()
  expect_true(!is.null(jury))
  expect_true(length(jury$get()$experts) > 0)
})
