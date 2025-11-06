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
  data <- as_gpredomics_data(as.data.frame(X), y, FALSE)
  
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
  data_t <- as_gpredomics_data(X_t, y, TRUE)
  
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
  
  data <- as_gpredomics_data(as.data.frame(X), y, FALSE)
  
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

    data <- as_gpredomics_data(X, y, FALSE)

    X_t <- as.data.frame(t(X))
    data_t <- as_gpredomics_data(X_t, y, TRUE)

    expect_equal(data$get()$X, expRust$train_data()$get()$X)
    expect_equal(data$get()$y, expRust$train_data()$get()$y)
    expect_equal(data$get()$classes, expRust$train_data()$get()$classes)
    expect_equal(data$get(), data_t$get())

    # Test load_data
    data_load <- expRust$load_data("../../sample/Xtrain.tsv", "../../sample/Ytrain.tsv") 
    expect_equal(expRust$train_data()$get()$X, data_load$get()$X)
    expect_equal(expRust$train_data()$get()$y, data_load$get()$y)
    expect_equal(expRust$train_data()$get()$classes, data_load$get()$classes)
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
  ind <- exp$individual(generation = 0, order = 0)
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

test_that("Individual: compute_auc, compute_metrics, compute_all", {
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
  ind <- exp$individual(0, 0)
  test_data <- exp$test_data()
  
  # Test compute_auc
  ind_auc <- ind$compute_auc(test_data)
  expect_true(!is.null(ind_auc))
  
  # Test compute_metrics
  ind_metrics <- ind$compute_metrics(test_data)
  expect_true(!is.null(ind_metrics))
  
  # Test compute_all
  ind_all <- ind$compute_all(test_data)
  expect_true(!is.null(ind_all))
  ind_all_data <- ind_all$get()
  expect_true(is.numeric(ind_all_data$auc))
  expect_true(is.numeric(ind_all_data$accuracy))
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
  ind <- exp$individual(0, 0)
  train_data <- exp$train_data()
  
  # Test predict()
  predictions <- ind$predict(train_data)
  expect_true(is.numeric(predictions))
  expect_true(all(predictions %in% c(0, 1)))
  # Verify prediction length matches number of samples
  td <- train_data$get()
  n_samples <- ncol(td$X)
  if (!is.null(n_samples)) expect_equal(length(predictions), n_samples)
  
  # Test predict_class_and_score()
  pred_result <- ind$predict_class_and_score(train_data)
  expect_true(is.list(pred_result))
  expect_true("class" %in% names(pred_result))
  expect_true("score" %in% names(pred_result))
  expect_true(length(pred_result$class) == length(pred_result$score))
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
  ind <- exp$individual(0, 0)
  
  # Get initial threshold
  initial_threshold <- ind$get()$threshold
  
  # Modify threshold
  new_threshold <- 0.75
  ind$set_threshold(new_threshold)
  
  # Verify threshold has changed
  updated_threshold <- ind$get()$threshold
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
  ind <- exp$individual(0, 0)
  
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
  
  pop_size <- exp$population_size(0)
  expect_true(is.numeric(pop_size))
  expect_true(pop_size > 0)
})

test_that("Experiment: get_generation returns individuals", {
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
  
  # Test: retrieving a generation
  generation <- exp$get_generation(0)
  expect_true(is.list(generation))
  expect_true(length(generation) > 0)
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
  gen0_pop <- exp$get_population(0)
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
  fbm <- pop$get_fbm(alpha = 0.05, min_pct_fallback = 5.0)
  expect_true(!is.null(fbm))
  
  # FBM should be smaller or equal to original population
  fbm_data <- fbm$get()
  pop_data <- pop$get()
  expect_true(length(fbm_data$individuals) <= length(pop_data$individuals))
})

test_that("Population: predict_scores_matrix and predict_classes_matrix", {
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
  
  # Test predict_scores_matrix
  scores_matrix <- pop$predict_scores_matrix(train_data)
  expect_true(is.data.frame(scores_matrix))
  expect_true(nrow(scores_matrix) > 0)
  expect_true(ncol(scores_matrix) > 0)
  # Try to verify orientation: number of columns or rows should match number of individuals
  pop_count <- length(pop$get()$individuals)
  expect_true((ncol(scores_matrix) == pop_count) || (nrow(scores_matrix) == pop_count))
  
  # Test predict_classes_matrix
  classes_matrix <- pop$predict_classes_matrix(train_data)
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
  pop <- pop$get_fbm(alpha = 0.05, min_pct_fallback = 5.0) 
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
  ind <- pop$get_individual(0)
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
  ind <- pop$get_individual(0)
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
  
  # Test print
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

test_that("Jury: predict_class_and_score", {
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
  
  # Test predictions
  jury$evaluate(train_data)
  predictions <- jury$predict_class_and_score(train_data)
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
  
  # Test print
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
  
  # 8. Make predictions
  train_data <- exp$train_data()
  jury$evaluate(train_data)
  predictions <- jury$predict_class_and_score(train_data)
  
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
  fbm <- pop$get_fbm(alpha = 0.05, min_pct_fallback = 5.0)
  expect_true(!is.null(fbm))
  
  # Compute importance (if method is available)
  # Note: compute_importance may require specific parameters
  train_data <- exp$train_data()
  
  # Test with an individual
  ind <- fbm$get_individual(0)
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
  expect_equal(exp$population_size(0), exp_loaded$population_size(0))

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
  expect_equal(exp$population_size(0), exp_loaded$population_size(0))

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
  expect_equal(exp$population_size(0), exp_loaded$population_size(0))

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
  pop_predictions <- pop$predict_classes_matrix(train_data)

  # Verify predictions for each individual
  for (i in seq_len(length(pop$get()$individuals))) {
    ind <- pop$get_individual(i - 1)
    ind_predictions <- ind$predict(train_data)

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

  ind <- pop$get_individual(0)
  ind_predictions <- ind$predict(train_data)
  ind_metrics <- ind$compute_all(train_data)

  true_labels <- train_data$get()$y
  predicted_labels <- ind_predictions

  tp <- sum(true_labels == 1 & predicted_labels == 1)
  tn <- sum(true_labels == 0 & predicted_labels == 0)
  fp <- sum(true_labels == 0 & predicted_labels == 1)
  fn <- sum(true_labels == 1 & predicted_labels == 0)

  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)

  expect_equal(ind_metrics$get()$accuracy, accuracy, tolerance = 1e-8)
  expect_equal(ind_metrics$get()$sensitivity, sensitivity, tolerance = 1e-8)
  expect_equal(ind_metrics$get()$specificity, specificity, tolerance = 1e-8)

  jury <- Jury$new_from_param(pop, param)
  jury$evaluate(train_data)
  jury_metrics <- jury$get()

  jury_predictions <- jury$predict_class_and_score(train_data)$class

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
  individuals <- c(exp$get_population(0)$get_individual(0),
                   exp$get_population(0)$get_individual(1),
                   exp$get_population(0)$get_individual(2))

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
