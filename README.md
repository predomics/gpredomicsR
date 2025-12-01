# Preliminary

[gpredomicsR](https://github.com/predomics/gpredomicsR) is an R package that wraps 
the [gpredomics](https://github.com/predomics/gpredomics) rust library, which 
implements accelerated heuristics to find sparse interpretable models described in 
[Prifti et al, Gigascience, 2020](https://academic.oup.com/gigascience/article/9/3/giaa010/5801229). 
The main algorithm is a genetic algorithm which applies numerous operators to evolve `individuals`, 
which are classification models. The fit of the models is evaluated on a penalized AUC.

# Installation

## From source

Open an R session at the package root directory:

```R
# Build and install in one step
devtools::install()

# Or build then install separately
devtools::build(".")
install.packages("../gpredomicsR_0.1.7.tar.gz", repos = NULL, type = "source")
```

## System Requirements

- **Rust**: Cargo (Rust's package manager) and rustc must be installed
- **Linux**: Vulkan libraries and drivers for GPU support (see [gpredomics](https://github.com/predomics/gpredomics) README for details)
- **macOS (Apple Silicon)**: Native GPU support via WGPU (Metal backend)
- **Windows**: GPU support may require additional setup

# Using

## Creating the environment

After loading the package, you need to:
1. Load parameters from a YAML file
2. Create a `RunningFlag` object for algorithm control
3. (Optional) Set up logging

```R
library(gpredomicsR)

# Load parameters
param <- Param$load("sample/param.yaml")

# Create running flag (allows to stop the algorithm)
running_flag <- RunningFlag$new()

# Optional: Configure logging
glog <- GLogger$level("info")  # Options: "trace", "debug", "info", "warn", "error"

# Access parameters as R list (optional)
rparam <- param$get()
```

**Note**: The `param` object is a Rust pointer wrapper. Use `param$get()` to retrieve parameters as an R list.

## Running the algorithm

The main function is `fit()`, which runs the genetic algorithm with the provided parameters.

```R
# Run the genetic algorithm
exp <- fit(param, running_flag)

# Basic experiment information
n_generations <- exp$generation_number()
print(paste("Algorithm ran for", n_generations, "generations"))

# Access specific individuals
# Note: generation and order indices are 1-based in R
best_individual <- exp$individual(generation = n_generations, order = 1)

# Get the final population
final_pop <- exp$get_population(n_generations)

# Access training and test data
train_data <- exp$train_data()
test_data <- exp$test_data()

# Convert experiment to R-friendly format
exp_results <- parseExperiment(exp)
```

**Key Methods**:
- `exp$generation_number()`: Number of generations run
- `exp$individual(generation, order)`: Get a specific individual (1-based indexing)
- `exp$get_population(generation)`: Get entire population for a generation
- `exp$train_data()` / `exp$test_data()`: Access data objects

## Analyzing the results

### Evolution analysis

Track how model attributes evolve across generations:

```R
# Analyze evolution of best models (default)
perf_df <- analyzeAttributeEvolution(
  exp, 
  attributes = c("auc", "fit", "k", "epoch"),
  best_model = TRUE,
  plot = TRUE
)

# View the data
head(perf_df)
```

### Population analysis

Analyze feature usage and importance across a population:

```R
# Get final population
final_pop <- exp$get_population(exp$generation_number())

# Filter population by performance
filtered_pop <- final_pop$filter_by_auc(0.7)

# Get Family of Best Models (FBM)
fbm <- final_pop$get_fbm(alpha = 0.05)

# Compute feature importance matrix
importance_matrix <- fbm$compute_importance_matrix(
  data = exp$train_data(),
  n_perm = 100,
  used_only = TRUE,
  seed = 42
)

# Predictions
predictions <- fbm$predict_class_matrix(exp$test_data())
scores <- fbm$predict_score_matrix(exp$test_data())
```

### Individual analysis

Examine specific models in detail:

```R
# Get best individual
best_ind <- exp$individual(generation = exp$generation_number(), order = 1)

# View individual details
best_ind$print()

# Get metrics
metrics <- best_ind$get_metrics()
print(metrics)

# Make predictions
preds <- best_ind$predict(exp$test_data())

# Feature importance for this individual
imp <- best_ind$compute_importance(
  data = exp$train_data(),
  n_perm = 100,
  seed = 42
)

# Prune low-importance features
pruned_ind <- best_ind$prune_by_threshold(
  threshold = 0.01,
  n_perm = 50,
  min_k = 1
)
```

### Jury/Ensemble methods

Create an ensemble of top models:

```R
# Create a jury from filtered population
jury <- Jury$new_from_param(fbm, param)

# Fit jury on data
jury$fit(exp$train_data())

# Make predictions
jury_preds <- jury$predict(exp$test_data())

# Get jury metrics
jury_metrics <- jury$get()
```

## Advanced Features

### Cross-Validation (CV)

Enable cross-validation in your parameter file or programmatically:

```R
# Enable CV in parameters
param$set_bool("cv", TRUE)
param$set("outer_folds", 5)

# Run experiment with CV
exp <- fit(param, running_flag)

# Access CV results
n_folds <- exp$get_n_folds()
best_pop <- exp$get_best_population()

# Get fold-specific data and populations
fold_train <- exp$get_fold_data(fold = 1, train = TRUE)
fold_valid <- exp$get_fold_data(fold = 1, train = FALSE)
fold_pop <- exp$get_fold_generation(fold = 1, generation = 10, train = TRUE)

# Compute CV importance matrix
cv_importance <- exp$compute_cv_importance_matrix(
  n_perm = 50,
  aggregation = "mean",
  scaled = TRUE
)
```

### GPU Acceleration

Enable GPU for faster computation:

```R
# In param.yaml: gpu: true
# Or programmatically:
param$set_bool("gpu", TRUE)
```

**Platform-specific requirements**:
- **Apple Silicon (M1/M2/M3)**: Native support via Metal (WGPU backend)
- **Linux**: Requires Vulkan drivers and libraries
- **Windows**: May require additional DirectX setup

### Thread Control

Control parallelization for computationally intensive operations:

```R
# Set number of threads
options(gpredomics.threads.number = 8)

# Methods that benefit from multi-threading:
# - Population$compute_importance_matrix()
# - Population$fit_on()
# - Individual$compute_importance()
```

If not set, the package defaults to available CPU cores or the value in `param.yaml`.

### Logging

Configure logging verbosity:

```R
# Set log level directly
glog <- GLogger$level("info")

# Or create then modify
glog <- GLogger$new()
glog$set_level("debug")

# Available levels: "trace", "debug", "info", "warn", "error"
```

**Note**: Logger can only be created once per session.

### Serialization

Save and load experiments:

```R
# Save experiment
exp$save("results/my_experiment.json")   # JSON format
exp$save("results/my_experiment.mp")     # MessagePack format
exp$save("results/my_experiment.bin")    # Binary format

# Load experiment
loaded_exp <- Experiment$load("results/my_experiment.json")
```

## Example Workflow

Complete example using sample data:

```R
library(gpredomicsR)

# 1. Setup
param <- Param$load("sample/param.yaml")
param$set("max_epochs", 50)
param$set("population_size", 1000)
param$set_string("X", "sample/Xtrain.tsv")
param$set_string("y", "sample/Ytrain.tsv")
param$set_string("Xtest", "sample/Xtest.tsv")
param$set_string("ytest", "sample/Ytest.tsv")

running_flag <- RunningFlag$new()
glog <- GLogger$level("info")

# 2. Run algorithm
exp <- fit(param, running_flag)

# 3. Analyze results
final_pop <- exp$get_population(exp$generation_number())
fbm <- final_pop$get_fbm(alpha = 0.05)

# 4. Feature importance
importance <- fbm$compute_importance_matrix(
  data = exp$train_data(),
  n_perm = 100
)

# 5. Create ensemble
jury <- Jury$new_from_param(fbm, param)
jury$fit(exp$train_data())
predictions <- jury$predict(exp$test_data())

# 6. Save results
exp$save("results/experiment.json")
```

## Documentation

For detailed documentation of all methods and classes:

```R
# R help system
?fit
?Individual
?Population
?Experiment
?Jury
?Data
?Param

# List all exported functions
help(package = "gpredomicsR")
```
