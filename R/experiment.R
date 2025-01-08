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



#' Retrieve Specific Attribute for Each Individual in a Population
#'
#' This function iterates through a list representing a population and extracts a specified attribute from each individual. If the attribute is not present, `NA` is assigned for that individual.
#'
#' @param pop A list of lists, where each inner list represents an individual and should contain named attributes.
#' @param attribute A character string specifying the name of the attribute to extract from each individual. The default is "fit".
#'
#' @return A vector containing the values of the specified attribute for each individual in the population. If an individual does not have the attribute, the corresponding entry in the vector will be `NA`.
#'
#' @examples
#' # Assuming a list of individuals with attributes
#' individuals <- list(
#'   list(fit = 0.95, auc = 0.90),
#'   list(fit = 0.89),
#'   list(auc = 0.92)  # This individual does not have 'fit'
#' )
#' # Extract 'fit' attribute
#' fits <- getPopulationAttribute(individuals, attribute = "fit")
#' print(fits)  # Prints: c(0.95, 0.89, NA)
#'
#' @export
getPopulationAttribute <- function(pop, attribute = "fit") {
  # Initialize a vector to store the results
  results <- vector("list", length(pop))
  
  # Loop over each individual in the population
  for (i in seq_along(pop)) {
    individual <- pop[[i]]
    
    # Check if the individual is a list and has the required attribute
    if (is.list(individual) && attribute %in% names(individual)) {
      results[[i]] <- individual[[attribute]]
    } else {
      results[[i]] <- NA  # Assign NA if the attribute is not present
    }
  }
  
  # Simplify the list to a vector, handling types automatically
  do.call(c, results)
}


# Document the function parameters and expected input
#' @param exp A list containing generations and associated data
#' @param attributes Character vector specifying which attributes to analyze
#' @param plot Logical indicating if results should be plotted
#' @return A data frame of combined attributes or NULL if plot is TRUE
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @examples
#' @export
#' exp <- list(generations = list(generation1, generation2, ...))
#' analyzeEvolution(exp, attributes = c("auc", "fit", "k"), plot = TRUE)
analyzeEvolution <- function(exp, attributes = c("auc", "fit", "k", "n"), plot = TRUE) {
  require(dplyr)
  require(ggplot2)
  require(tidyr)
  
  # Function to fetch and bind data for a given attribute
  getAttributeData <- function(attribute, generations) {
    attribute_data <- lapply(generations, getPopulationAttribute, attribute = attribute)
    do.call(rbind, attribute_data)
  }
  
  # Extract generations from experiment data
  generations <- parseExperiment(exp)
  
  # Using Map to apply function over each attribute
  list_data <- setNames(Map(getAttributeData, attributes, MoreArgs = list(generations = generations)), attributes)
  
  # Create the performance dataframe
  perf.df <- data.frame(
    auc = list_data$auc[,1],
    fit = list_data$fit[,1],
    k = list_data$k[,1],
    age = list_data$n[,1],
    generation = seq_along(generations)
  )
  
  if (plot) {
    # Transform perf.df to long format
    perf.df.long <- perf.df %>%
      tidyr::pivot_longer(
        cols = -generation,
        names_to = "measure",
        values_to = "value"
      )
    
    # Plotting the data
    ggplot(perf.df.long, aes(x = generation, y = value, group = measure, color = measure)) +
      geom_line() +
      theme_minimal() +
      facet_grid(measure ~ ., scales = "free_y") +
      labs(title = "Best individual markers across generations",
           x = "Generations",
           y = "Attributes") + 
      theme(legend.position = "none")
    
  } else {
    return(perf.df)
  }
}


#' Analyze Evolution Across All Models
#'
#' This function processes experimental data across generations, extracting specified
#' attributes and optionally plotting their evolution over time using facets. Each facet represents
#' an attribute across all generations.
#'
#' @param exp A list containing generations and associated data.
#' @param attributes Character vector specifying which attributes to analyze.
#' @param plot Logical indicating if results should be plotted.
#' @return A data frame of combined attributes or NULL if plot is TRUE.
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @examples
#' exp <- list(generations = list(generation1, generation2, ...))
#' analyzeEvolutionAllModels(exp, attributes = c("auc", "fit", "k", "n"), plot = TRUE)
#' @export
analyzeEvolutionAllModels <- function(exp, attributes = c("auc", "fit", "k", "n"), plot = TRUE) {
  require(dplyr, quietly = TRUE)
  require(ggplot2, quietly = TRUE)
  require(tidyr, quietly = TRUE)
  
  data.list <- list()
  for (i in seq_along(attributes))
  {
    data <- sapply(generations, getPopulationAttribute, attribute = attributes[i])
    data <- t(data)
    
    data.long <- as.data.frame(data) %>%
      rownames_to_column(var = "generation") %>%
      pivot_longer(cols = -generation, names_to = "individual", values_to = "fit") %>%
      mutate(generation = factor(fit.df.long$generation,
                                 levels = unique(fit.df.long$generation[order(as.numeric(gsub("gen_", "", fit.df.long$generation)))])))
    data.long$attribute <- attributes[i]
    
    data.list[[attributes[i]]] <- data.long
  }
  
  data <- do.call(rbind, data.list)
  rownames(data) <- NULL
  
  
  if (plot) {
    # Plotting the data with facets for each attribute
    p <- ggplot(data, aes(x = generation, y = fit, group = individual, color = individual)) +
      geom_line(alpha = 0.5) +  # Set alpha for line transparency
      theme_minimal() +
      facet_grid(attribute ~ ., scales = "free_y") +  # Facet by attribute
      labs(title = "Evolution of Model Attributes Across Generations",
           x = "Generation",
           y = "Attribute Value",
           color = "Attribute") +
      theme(legend.position = "none",  # Remove legend if not needed
            axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
    
    print(p)  # Display the plot
  } else {
    return(data)
  }
}
