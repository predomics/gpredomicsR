# Caret wrapper for Gpredomics - THESE FUNCTIONS REMAIN EXPERIMENTAL
# This script demonstrates different ways to use gpredomics with caret
# for hyperparameter tuning, from simple to advanced approaches.

rextendr::document()

library(readr)
library(tibble)
library(dplyr)
library(caret)
library(gpredomicsR)

options(gpredomics.threads.number = 8)

#---------------------------------------------
# Discovering Available Parameters
#---------------------------------------------

# Get information about ALL available parameters
# (not just the ones in the default grid!)
all_params <- gpredomics_params_info()
print(all_params)

# Get parameters by category
ga_params <- gpredomics_params_info("ga")
penalty_params <- gpredomics_params_info("penalties")
data_params <- gpredomics_params_info("data")

# The default grid only includes a few common parameters
# But you can add ANY of these parameters to your custom grid

#---------------------------------------------
# Loading Qin2014 Data
#---------------------------------------------

# Load training data (samples in columns, features in rows)
df_X_train <- read_table("sample/Xtrain.tsv")
df_X_train <- df_X_train %>% column_to_rownames("msp_name")

# Transpose to get samples in rows (required by caret)
X_train <- as.data.frame(t(df_X_train))

# Load training labels
df_y_train <- read_table("sample/Ytrain.tsv", col_names = FALSE, skip = 1)
y_train <- factor(
  ifelse(df_y_train[, 2] == 0, "healthy", "cirrhosis"),
  levels = c("healthy", "cirrhosis")
)

# Create training data frame
train_data <- X_train
train_data$Class <- y_train

# Load test data
df_X_test <- read_table("sample/Xtest.tsv")
df_X_test <- df_X_test %>% column_to_rownames("msp_name")
X_test <- as.data.frame(t(df_X_test))

# Load test labels
df_y_test <- read_table("sample/Ytest.tsv", col_names = FALSE, skip = 1)
y_test <- factor(
  ifelse(df_y_test[, 2] == 0, "healthy", "cirrhosis"),
  levels = c("healthy", "cirrhosis")
)

# Create test data frame
test_data <- X_test
test_data$Class <- y_test

#---------------------------------------------
# Example 1: Simple Training (No Cross-Validation)
#---------------------------------------------

# Train with default parameters, no CV
# Uses method="none" in trainControl to skip cross-validation
ctrl_none <- trainControl(
  method = "none",
  classProbs = TRUE,
  summaryFunction = gpredomicsSummary
)

# Define single parameter set (aligned with param.yaml defaults)
params <- data.frame(
  population_size = 5000,
  max_epochs = 100,
  ga_kmin = 1,
  ga_kmax = 200,
  beam_kmin = 1,
  beam_kmax = 200,
  language = "bin",
  data_type = "raw"
)

model_no_cv <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_none,
  tuneGrid = params,
  metric = "ROC"
)

# Make predictions
preds_no_cv <- predict(model_no_cv, test_data)

# Evaluate
confusionMatrix(preds_no_cv, test_data$Class)

#---------------------------------------------
# Example 2: Training with Cross-Validation
#---------------------------------------------

# Train with 5-fold cross-validation
# This evaluates performance more robustly
ctrl_cv <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = gpredomicsSummary
)

# Single parameter set (aligned with param.yaml defaults)
params_cv <- data.frame(
  population_size = 5000,
  max_epochs = 100,
  ga_kmin = 1,
  ga_kmax = 200,
  beam_kmin = 1,
  beam_kmax = 200,
  language = "bin",
  data_type = "raw"
)

model_cv <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_cv,
  tuneGrid = params_cv,
  metric = "ROC"
)

# Check CV results
model_cv$results

# Make predictions on test set
preds_cv <- predict(model_cv, test_data)
confusionMatrix(preds_cv, test_data$Class)

#---------------------------------------------
# Example 3: Random Search for Hyperparameters
#---------------------------------------------

# Use random search to explore hyperparameter space
# Faster than exhaustive grid search for large parameter spaces
ctrl_random <- trainControl(
  method = "cv",
  number = 3,
  search = "random",
  classProbs = TRUE,
  summaryFunction = gpredomicsSummary
)

model_random <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_random,
  tuneLength = 10,  # Test 10 random combinations
  metric = "ROC"
)

# View all tested parameter combinations
model_random$results

# Best parameters found
model_random$bestTune

# Make predictions
preds_random <- predict(model_random, test_data)
confusionMatrix(preds_random, test_data$Class)

#---------------------------------------------
# Example 4: Grid Search with Custom Grid
#---------------------------------------------

# Define custom grid for systematic exploration
# Test different population sizes and epoch counts
ctrl_grid <- trainControl(
  method = "cv",
  number = 3,
  classProbs = TRUE,
  summaryFunction = gpredomicsSummary
)

custom_grid <- expand.grid(
  population_size = c(1000, 5000),
  max_epochs = c(50, 100),
  ga_kmin = c(1, 5),
  ga_kmax = c(50, 200),
  beam_kmin = c(1, 5),
  beam_kmax = c(50, 200),
  language = c("bin", "ter"),
  data_type = c("raw", "log")
)

# This will test 2 × 2 × 2 × 2 × 2 × 2 = 64 combinations
model_grid <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_grid,
  tuneGrid = custom_grid,
  metric = "ROC"
)

# View results sorted by ROC
model_grid$results %>%
  arrange(desc(ROC)) %>%
  head(10)

# Plot results
plot(model_grid, metric = "ROC")

# Best configuration
model_grid$bestTune

# Final predictions
preds_grid <- predict(model_grid, test_data)

# Evaluation
confusionMatrix(preds_grid, test_data$Class)

#---------------------------------------------
# Example 5: Advanced - Compare Multiple Models
#---------------------------------------------

# Use resamples() to compare different configurations
# Train multiple models for comparison
ctrl_compare <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = gpredomicsSummary,
  savePredictions = "final",
  index = createFolds(train_data$Class, k = 5)  # Use same folds for fair comparison
)

# Small model (fast, reduced size for testing)
grid_small <- data.frame(
  population_size = 1000,
  max_epochs = 50,
  ga_kmin = 1,
  ga_kmax = 50,
  beam_kmin = 1,
  beam_kmax = 50,
  language = "bin",
  data_type = "raw"
)

model_small <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_compare,
  tuneGrid = grid_small,
  metric = "ROC"
)

# Medium model (default parameters)
grid_medium <- data.frame(
  population_size = 5000,
  max_epochs = 100,
  ga_kmin = 1,
  ga_kmax = 200,
  beam_kmin = 1,
  beam_kmax = 200,
  language = "bin",
  data_type = "raw"
)

model_medium <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_compare,
  tuneGrid = grid_medium,
  metric = "ROC"
)

# Large model (more epochs and population)
grid_large <- data.frame(
  population_size = 10000,
  max_epochs = 150,
  ga_kmin = 1,
  ga_kmax = 200,
  beam_kmin = 1,
  beam_kmax = 200,
  language = "bin",
  data_type = "raw"
)

model_large <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_compare,
  tuneGrid = grid_large,
  metric = "ROC"
)

# Compare models using resamples
comparison <- resamples(list(
  Small = model_small,
  Medium = model_medium,
  Large = model_large
))

# Summary statistics
summary(comparison)

# Visual comparison
bwplot(comparison, metric = "ROC")
dotplot(comparison, metric = "ROC")

# Statistical tests
diff_results <- diff(comparison)
summary(diff_results)

#---------------------------------------------
# Example 6: Testing Different Model Languages and Data Types
#---------------------------------------------

# This example explores different model languages and data transformations
# Language options: "ter" (ternary), "bin" (binary), "ratio", "pow2" (power of 2)
# Data type options: "log" (log-transform), "raw" (no transform), "prev" (prevalence)

ctrl_lang <- trainControl(
  method = "cv",
  number = 3,
  classProbs = TRUE,
  summaryFunction = gpredomicsSummary
)

# Test all combinations of language and data_type
grid_lang_data <- expand.grid(
  population_size = 5000,
  max_epochs = 100,
  ga_kmin = 1,
  ga_kmax = 200,
  beam_kmin = 1,
  beam_kmax = 200,
  language = c("bin", "ter", "ratio", "pow2"),
  data_type = c("raw", "log", "prev"),
  stringsAsFactors = FALSE
)

# This creates 4 × 3 = 12 combinations
model_lang_data <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_lang,
  tuneGrid = grid_lang_data,
  metric = "ROC"
)

# View results by language and data_type
model_lang_data$results %>%
  arrange(desc(ROC)) %>%
  select(language, data_type, ROC, Sens, Spec, MCC) %>%
  head(12)

# Best combination
model_lang_data$bestTune

# Plot performance by language
ggplot(model_lang_data$results, aes(x = language, y = ROC, fill = data_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "Performance by Language and Data Type",
    x = "Model Language",
    y = "ROC AUC",
    fill = "Data Type"
  )

# Final predictions with best configuration
preds_lang <- predict(model_lang_data, test_data)
confusionMatrix(preds_lang, test_data$Class)

#---------------------------------------------
# Example 7: Advanced - Custom Parameters (fit, algo, penalties, etc.)
#---------------------------------------------

# The grid in caret.R only shows DEFAULT tunable parameters
# You can add ANY parameter supported by gpredomics to your custom grid!

# Example: Test different fitness metrics
ctrl_fit <- trainControl(
  method = "cv",
  number = 3,
  classProbs = TRUE,
  summaryFunction = gpredomicsSummary
)

grid_fit_metrics <- expand.grid(
  population_size = 5000,
  max_epochs = 100,
  ga_kmin = 1,
  ga_kmax = 200,
  beam_kmin = 1,
  beam_kmax = 200,
  language = "bin",
  data_type = "raw",
  fit = c("auc", "mcc", "f1_score", "g_mean"),  # Test different fitness metrics
  stringsAsFactors = FALSE
)

model_fit_comparison <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_fit,
  tuneGrid = grid_fit_metrics,
  metric = "MCC"
)

# View results by fitness metric
model_fit_comparison$results %>%
  arrange(desc(MCC)) %>%
  select(fit, MCC, Sens, Spec, F1) %>%
  head()

# Example: Add penalties to the grid
grid_penalties <- expand.grid(
  population_size = 5000,
  max_epochs = 100,
  ga_kmin = 1,
  ga_kmax = 200,
  beam_kmin = 1,
  beam_kmax = 200,
  language = "bin",
  data_type = "raw",
  k_penalty = c(0, 0.0001, 0.001),  # Penalize model complexity
  fr_penalty = c(0, 0.05),          # False rate penalty
  stringsAsFactors = FALSE
)

model_penalties <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_fit,
  tuneGrid = grid_penalties,
  metric = "MCC"
)

# View impact of penalties
model_penalties$results %>%
  arrange(desc(MCC)) %>%
  select(k_penalty, fr_penalty, MCC, Sens, Spec) %>%
  head(10)

# Example: Test different algorithms (GA, BEAM, MCMC)
grid_algorithms <- expand.grid(
  population_size = 5000,
  max_epochs = 100,
  ga_kmin = 1,
  ga_kmax = 200,
  beam_kmin = 1,
  beam_kmax = 200,
  language = "bin",
  data_type = "raw",
  algo = c("ga", "beam"),  # Compare GA vs BEAM
  stringsAsFactors = FALSE
)

model_algo_comparison <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_fit,
  tuneGrid = grid_algorithms,
  metric = "MCC"
)

# Compare algorithms
model_algo_comparison$results %>%
  select(algo, MCC, Sens, Spec, F1) %>%
  arrange(desc(MCC))

# You can combine multiple parameters in one grid!
grid_comprehensive <- expand.grid(
  population_size = c(5000, 10000),
  max_epochs = c(100, 150),
  ga_kmin = 1,
  ga_kmax = 200,
  beam_kmin = 1,
  beam_kmax = 200,
  language = c("bin", "ter"),
  data_type = c("raw", "log"),
  fit = c("auc", "mcc"),
  k_penalty = c(0, 0.0001),
  algo = "ga",
  stringsAsFactors = FALSE
)

# This creates 2 × 2 × 2 × 2 × 2 × 2 = 64 combinations
model_comprehensive <- train(
  Class ~ .,
  data = train_data,
  method = getModelInfo_gpredomics(),
  trControl = ctrl_fit,
  tuneGrid = grid_comprehensive,
  metric = "MCC"
)

# Analyze best configurations
model_comprehensive$results %>%
  arrange(desc(MCC)) %>%
  select(algo, fit, language, data_type, k_penalty, MCC, Sens, Spec) %>%
  head(10)

# Final model with best configuration
model_comprehensive$bestTune
