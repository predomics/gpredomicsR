\name{caret-integration}
\alias{caret-integration}
\alias{gpredomics_train}
\alias{gpredomics_predict}
\alias{gpredomics_prob}
\alias{gpredomics_grid}
\alias{gpredomics_sort}
\alias{gpredomics_info}
\alias{register_gpredomics_caret}
\alias{train_gpredomics_caret}

\title{Caret Integration for gpredomicsR}

\description{
Functions to integrate gpredomics with the caret package for hyperparameter
tuning using grid search and random search.
}

\usage{
register_gpredomics_caret()

train_gpredomics_caret(
  formula, 
  data, 
  method = "grid",
  tuneLength = 3,
  tuneGrid = NULL,
  trControl = NULL,
  metric = "ROC",
  ...
)

gpredomics_train(x, y, wts = NULL, param, lev, last, weights, classProbs, ...)
gpredomics_predict(modelFit, newdata, preProc = NULL, submodels = NULL)
gpredomics_prob(modelFit, newdata, preProc = NULL, submodels = NULL)
gpredomics_grid(x, y, len = NULL, search = "grid")
gpredomics_sort(x, metric, maximize)
gpredomics_info()
}

\arguments{
  \item{formula}{A formula specifying the model (e.g., \code{y ~ .})}
  \item{data}{A data.frame containing the training data}
  \item{method}{Character string: "grid" for grid search or "random" for random search}
  \item{tuneLength}{Integer specifying the number of parameter combinations to try}
  \item{tuneGrid}{Optional data.frame with custom hyperparameter grid}
  \item{trControl}{A \code{trainControl} object from caret (defaults to 5-fold CV)}
  \item{metric}{Performance metric to optimize (default: "ROC")}
  \item{...}{Additional arguments passed to \code{caret::train()}}
  \item{x}{Feature matrix for training}
  \item{y}{Response vector (binary: 0/1 or factor)}
  \item{wts}{Sample weights (currently not used)}
  \item{param}{Data.frame row with hyperparameters}
  \item{lev}{Class levels}
  \item{last}{Logical indicating if this is the final model}
  \item{weights}{Sample weights (currently not used)}
  \item{classProbs}{Logical indicating if class probabilities should be computed}
  \item{modelFit}{Fitted model object}
  \item{newdata}{New data for prediction}
  \item{preProc}{Pre-processing object (currently not used)}
  \item{submodels}{List of submodels (currently not used)}
  \item{len}{Number of parameter combinations for random search}
  \item{search}{Type of search: "grid" or "random"}
  \item{maximize}{Logical indicating if metric should be maximized}
}

\details{
The caret integration provides two main workflows:

\strong{1. Automatic registration and training:}
Use \code{train_gpredomics_caret()} which automatically handles registration
and calls \code{caret::train()} with appropriate settings.

\strong{2. Manual registration:}
Call \code{register_gpredomics_caret()} once, then use \code{caret::train()}
directly with \code{method = "gpredomics"}.

\strong{Hyperparameters available for tuning:}
\itemize{
  \item \code{population_size}: Size of the genetic algorithm population
  \item \code{max_epochs}: Maximum number of generations
  \item \code{kmin}: Minimum number of features in models
  \item \code{kmax}: Maximum number of features in models
}

Additional parameters can be added to the grid including:
\code{min_epochs}, \code{select_elite_pct}, \code{mutated_children_pct},
\code{k_penalty}, \code{fr_penalty}, \code{algo}, \code{language}, \code{fit}.

\strong{Grid Search:}
Evaluates all combinations in a predefined grid. Use \code{search = "grid"}
in \code{trainControl()}.

\strong{Random Search:}
Randomly samples from the parameter space. Use \code{search = "random"}
in \code{trainControl()}. Often more efficient for large parameter spaces.
}

\value{
\code{register_gpredomics_caret()} returns NULL invisibly and prints a confirmation message.

\code{train_gpredomics_caret()} returns a \code{train} object from caret containing:
\itemize{
  \item Best hyperparameter values
  \item Cross-validation results
  \item Final fitted model
  \item Resampling performance metrics
}

The other functions are internal to the caret interface and return objects
compatible with caret's expectations.
}

\examples{
\dontrun{
library(gpredomicsR)
library(caret)

# Load example data
data(iris)
iris_binary <- iris[iris$Species != "setosa", ]
iris_binary$Species <- droplevels(iris_binary$Species)

# Method 1: Using convenience wrapper
model <- train_gpredomics_caret(
  Species ~ .,
  data = iris_binary,
  method = "grid",
  tuneLength = 2,
  metric = "ROC"
)

print(model)
plot(model)

# Method 2: Manual registration with custom grid
register_gpredomics_caret()

custom_grid <- expand.grid(
  population_size = c(50, 100, 200),
  max_epochs = c(50, 100),
  kmin = c(3, 5),
  kmax = c(10, 15)
)

ctrl <- trainControl(
  method = "cv",
  number = 5,
  search = "grid",
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

model <- train(
  Species ~ .,
  data = iris_binary,
  method = "gpredomics",
  trControl = ctrl,
  tuneGrid = custom_grid,
  metric = "ROC"
)

# Method 3: Random search
ctrl_random <- trainControl(
  method = "cv",
  number = 10,
  search = "random",
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

model_random <- train(
  Species ~ .,
  data = iris_binary,
  method = "gpredomics",
  trControl = ctrl_random,
  tuneLength = 20,  # Try 20 random combinations
  metric = "ROC"
)

# Compare results
plot(model_random)
print(model_random$bestTune)

# Make predictions
predictions <- predict(model, newdata = iris_binary)
probabilities <- predict(model, newdata = iris_binary, type = "prob")
}
}

\seealso{
\code{\link{fit}}, \code{\link{fit_on}}, \code{\link{Param}}, 
\code{\link[caret]{train}}, \code{\link[caret]{trainControl}}
}

\author{gpredomicsR Development Team}
