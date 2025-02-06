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
  
  # Ensure all elements of exp are valid populations
  if (!all(sapply(exp, isPopulation))) {
    stop("analyzeEvolutionBestModel: Some elements in the experiment are not valid populations.")
  }
  
  # Extract generations
  generations <- seq_along(exp)  # Assume experiment list is ordered by generation
  
  # Function to extract attribute data from populations
  getAttributeData <- function(attribute) {
    attribute_data <- lapply(seq_along(exp), function(gen_idx) {
      pop <- exp[[gen_idx]]
      
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
