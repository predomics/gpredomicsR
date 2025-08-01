#' Run a Predomics Experiment
#'
#' This function executes an experiment using the `gpredomicsR` package by loading
#' experiment parameters, setting up a logging mechanism, and running a genetic
#' algorithm. It returns a structured list containing the experiment results.
#'
#' @param param.path A string specifying the path to the parameter file (YAML format).
#'   Default is `"sample/param.yaml"`.
#' @param name an optional string to name the experiment.
#' @param glog_level the verbose level for glog (default: "debug").
#'
#' @return A list containing:
#' \item{rust}{A list with Rust objects: the experiment, parameters, running flag, and logger.}
#' \item{params}{The parameter settings loaded from the YAML file.}
#' \item{data}{A list with training and testing datasets.}
#' \item{model_collection}{An R representation of the experiment's model evolution.}
#' \item{execTime}{The execution time of the experiment in minutes.}
#'
#' @details
#' The function follows these steps:
#' - Saves the current working directory and switches to the experiment directory.
#' - Initializes or retrieves a global logger to prevent duplicate loggers.
#' - Loads parameters from the YAML file.
#' - Runs the genetic algorithm to generate and evaluate models.
#' - Constructs and returns an experiment object containing key components.
#'
#' @examples
#' \dontrun{
#' exp <- runExperiment(param.path = "sample/param.yaml")
#' print(exp)
#' }
#'
#' @author Edi Prifti (IRD)
#' @export
runExperiment <- function(param.path = "sample/param.yaml", name = "", glog_level = "debug") {
  startingTime <- Sys.time()
  
  # Extract the directory path from param.path
  dir.path <- dirname(param.path)
  
  # Save the current working directory
  original.dir <- getwd()
  
  # Ensure working directory is restored even if an error occurs
  on.exit(setwd(original.dir), add = TRUE)
  
  # Try to use an existing global logger, otherwise create a new one
  if (!exists("global_glog", envir = .GlobalEnv)) {
    message("Initializing new GLogger instance...")
    glog_attempt <- tryCatch(GLogger$level(level = glog_level), error = function(e) NULL)
    
    if (!is.null(glog_attempt)) {
      assign("global_glog", glog_attempt, envir = .GlobalEnv)
    } else {
      warning("Logger initialization failed. Proceeding without logging.")
    }
  }
  
  # Get the global logger (even if initialization failed, this avoids re-initialization)
  glog <- get("global_glog", envir = .GlobalEnv, inherits = FALSE)
  
  running_flag <- RunningFlag$new()
  
  # Load parameters
  paramRust <- Param$load(param.path)
  
  # Change to the new directory
  setwd(dir.path)
  
  # Run the genetic algorithm
  expRust <- fit(paramRust, running_flag)
  
  # Create experiment structure
  experiment <- list(
    rust = list(
      experiment = expRust,
      param = paramRust,
      running_flag = running_flag,
      glog = glog
    ),
    params = paramRust$get(),
    data = list(
      train = expRust$get_data_robj(train = TRUE),
      test = expRust$get_data_robj(train = FALSE)
    ),
    model_collection = parseExperiment(expRust),  # Convert Rust pointer to R object
    execTime = as.numeric(Sys.time() - startingTime, units = "mins")
  )
  
  # fix the class
  # for the TRAIN dataset
  # Ensure it's a vector
  if (!is.vector(experiment$data$train$y)) {
    experiment$data$train$y <- as.vector(experiment$data$train$y)  # Force conversion
  }
  # Convert binary numeric to factor
  unique_vals <- unique(experiment$data$train$y)
  if (length(table(experiment$data$train$y)) == 2) {
    experiment$data$train$y <- as.factor(experiment$data$train$y)
  }
  
  # same thing for the TEST dataset
  # Ensure it's a vector
  if (!is.vector(experiment$data$test$y)) {
    experiment$data$test$y <- as.vector(experiment$data$test$y)  # Force conversion
  }
  # Convert binary numeric to factor
  unique_vals <- unique(experiment$data$test$y)
  if (length(table(experiment$data$test$y)) == 2) {
    experiment$data$test$y <- as.factor(experiment$data$test$y)
  }
  
  return(experiment)
}


load_experiment <- function(path){
  startingTime <- Sys.time()
  
  expRust <- Experiment$load(path)
  
  paramRust <- expRust$get_param()
  
  experiment <- list(
    rust = list(
      experiment = expRust,
      param = paramRust#,
      # running_flag = running_flag,
      # glog = glog
    ),
    params = paramRust$get(),
    data = list(
      train = expRust$get_data_robj(train = TRUE),
      test = expRust$get_data_robj(train = FALSE)
    ),
    model_collection = parseExperiment(expRust),  # Convert Rust pointer to R object
    execTime = as.numeric(Sys.time() - startingTime, units = "mins")
  )
}




#' @title parseExperiment
#' @description Parse the experiment object to extract the individuals
#' @param exp The experiment object
#' @return A list of generations, each generation containing a list of individuals
#' @export
parseExperiment <- function(exp)
{
  # check validity of the experiment
  if (is.null(exp))
  {
    print("Experiment is null")
    return(NULL)
  }
  
  nb_generations <- exp$generation_number()
  generations <- list()
  for (i in 1:nb_generations)
  {
    generations[[paste0("gen_",i)]] <- exp$get_generation(i-1) # the index in rust is 0-based
  }
  return(generations)
}


#' Analyze evolution of model attributes over generations
#'
#' @description Extracts and visualizes how selected model attributes change over generations in an experiment.
#' @param exp A list representing the experiment data, containing populations over generations.
#' @param attributes A vector of attribute names to analyze (default: c("auc", "fit", "k", "epoch")).
#' @param best_model Logical; should the analyses focus on the best model? (default: TRUE).
#' @param plot Logical; should a plot be generated? (default: TRUE).
#' @return A dataframe containing the evolution of attributes over generations, or a plot if `plot = TRUE`.
#' @export
analyzeAttributeEvolution <- function(exp, attributes = c("auc", "fit", "k", "epoch"), best_model = TRUE, plot = TRUE) {
  require(dplyr)
  require(ggplot2)
  require(tidyr)
  
  # Validate input experiment
  if (!is.list(exp) || length(exp) == 0) {
    stop("analyzeEvolutionBestModel: Experiment must be a non-empty list.")
  }
  
  if (!is.list(exp$model_collection) || length(exp$model_collection) == 0) {
    stop("analyzeEvolutionBestModel: Experiment should contain a model collection (list of populations).")
  }
  
  # Ensure all elements of exp are valid populations
  if (!all(sapply(exp$model_collection, isPopulation))) {
    stop("analyzeEvolutionBestModel: Some elements in the experiment are not valid populations.")
  }
  
  # Extract generations
  generations <- seq_along(exp$model_collection)  # Assume experiment list is ordered by generation
  
  # Function to extract attribute data from populations
  getAttributeData <- function(attribute) {
    attribute_data <- lapply(seq_along(exp$model_collection), function(gen_idx) {
      pop <- exp$model_collection[[gen_idx]]
      
      # Select best model if needed
      if (best_model) {
        pop <- list(getTheBestIndividual(pop, evalToFit = "fit"))  # Extract best model
      }
      
      data <- populationGet_X(element2get = attribute, toVec = TRUE, na.rm = FALSE)(pop)
      
      if (is.null(data) || length(data) == 0) {
        return(data.frame())
      }
      
      # Convert to data frame
      df <- data.frame(
        generation = gen_idx,
        individual = seq_along(data),
        attribute = attribute,
        value = data
      )
      
      return(df)
    })
    
    do.call(rbind, attribute_data)
  }
  
  # Extract data for all attributes
  list_data <- lapply(attributes, getAttributeData)
  
  # Combine into a single dataframe
  perf.df.long <- do.call(rbind, list_data)
  
  # Handle case where no data is available
  if (nrow(perf.df.long) == 0 || all(is.na(perf.df.long$value))) {
    message("No data available for analysis.")
    return(NULL)
  }
  
  # Convert generation to ordered factor for proper plotting
  perf.df.long <- perf.df.long %>%
    mutate(generation = factor(generation, levels = sort(unique(generation), decreasing = FALSE)))
  
  if (plot) {
    p <- ggplot(perf.df.long, aes(x = generation, y = value, group = individual, color = ifelse(best_model, attribute, individual))) +
      geom_line(alpha = 0.6) +
      theme_minimal() +
      facet_wrap(~ attribute, scales = "free_y") +
      labs(title = "Evolution of Model Attributes Across Generations",
           x = "Generation",
           y = "Attribute Value",
           color = ifelse(best_model, "Attribute", "Individual")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    if (!best_model) {
      p <- p + theme(legend.position = "none")  # Hide legend if not focusing on best model
    }
    
    print(p)  # Display the plot
  } else {
    return(perf.df.long)
  }
}
