#' Print Model Information
#'
#' This function prints information about a given model in different formats.
#' Supports short, long, and structured outputs.
#'
#' @param mod A model object.
#' @param method Print format: "short" (default), "long", "str".
#' @param score Score attribute to display (default: "fit").
#' @return A formatted string summarizing the model.
#' @export
printModel <- function(mod, method = "short", score = "fit") {
  
  if (!isModel(mod)) {
    print("printModel: please provide a valid model object.")
    return(NULL)
  }
  
  if (!score %in% names(mod)) {
    print("printModel: provided score not found in model attributes.")
    return(NULL)
  }
  
  switch(method,
         
         # ---- Short Format ----
         short = {
           if (all(!is.null(mod$coefficients), length(mod$coefficients) > 0)) {
             ind.pos <- mod$coeff > 0
             ind.neg <- mod$coeff < 0
             
             term.pos <- if (any(ind.pos)) paste0("+", sort(mod$features[ind.pos]), collapse = " ") else "0"
             term.neg <- if (any(ind.neg)) paste0("+", sort(mod$features[ind.neg]), collapse = " ") else "0"
             
             res <- paste0("(", term.pos, ")", " - ", "(", term.neg, ")", " ≥ ", format(mod$threshold, scientific = TRUE), " then ", "1")
             
             mod.fit <- mod[[score]]
             fit_info <- if (!is.na(mod.fit)) paste("|(F=", signif(mod.fit, 4), ")", sep="") else ""
             
             res <- paste0(res, " {", "|K=", mod$k, "|La=", mod$language, "|Data=", mod$data_type, fit_info, "}", sep="")
             
           } else {
             res <- "No coefficients available."
           }
         },
         
         
         # ---- Structure Output ----
         str = {
           res <- capture.output(str(mod))
           res <- paste(res, collapse = "\n")
         },
         
         {
           warning("Invalid method! Choose from: short, long, str")
           return(NULL)
         }
  )
  
  return(res)
}


#' Print Information about a Population of Models
#'
#' Prints detailed information about a population of models using different display formats.
#'
#' @param obj A population of models (list of model objects).
#' @param method Print format: "digested", "short", "long", "str".
#' @param score Score attribute to display (default: "fit").
#' @param indent String indentation for structured printing.
#' @return None. Prints the information directly.
#' @export
printPopulation <- function(obj, method = "short", score = "fit", indent = "") {
  
  # Check if the object is a valid population
  if (!isPopulation(obj)) {
    warning("printPopulation: the object provided is not a valid population of models.")
    return(NULL)
  }
  
  switch(method,
         
         # ---- Digested Summary ----
         digested = {
           sparsity_values <- populationGet_X("k")(obj)
           unique_k_values <- unique(sparsity_values)
           
           # Define attribute name
           attribute_name <- if (length(unique_k_values) == 1) {
             paste0("k_", unique_k_values)
           } else {
             "k_mixed"
           }
           
           # Convert population to dataframe
           pop_df <- populationToDataFrame(obj)
           
           # Extract summary info (first 5 models max)
           summary_info <- paste(
             paste(pop_df$language, signif(pop_df$fit, 2), pop_df$k, sep = "_")[1:min(5, length(obj))],
             collapse = "; "
           )
           
           # Print population summary
           cat(paste(indent, attribute_name, ": ", length(obj), " models ... ", summary_info, "\n", sep=""))
         },
         
         # ---- Short / Long / Structured Model Output ----
         short =, 
         str ={
           for (i in seq_along(obj)) {
             mod <- obj[[i]]
             cat(paste0(i, ": ", printModel(mod, method = method, score = score), "\n"))
           }
         },
         
         # ---- Invalid Method Warning ----
         {
           warning("printPopulation: please provide a valid method (digested/short/long/str).")
         }
  )
}



#' Print Information about a Model Collection
#'
#' Prints detailed information about a collection of model populations.
#'
#' @param obj A model collection object (list of populations).
#' @param indent Indentation string for hierarchical formatting (default: "\t--- ").
#' @param method Print format: "short" (summary) or "long" (detailed, per population).
#' @return None. Prints the information directly.
#' @export
printModelCollection <- function(obj, indent = "\t--- ", method = "short") {
  
  # Check if the object is a valid model collection
  if (!isModelCollection(obj)) {
    warning("printModelCollection: the object provided is not a valid model collection.")
    return(NULL)
  }
  
  switch(method,
         
         # ---- Short Summary ----
         short = {
           summary_text <- paste(names(obj), unlist(lapply(obj, length)), sep = ": ")
           cat(paste(summary_text, collapse = " | "), "\n")
         },
         
         # ---- Detailed Population Info ----
         long = {
           for (i in seq_along(obj)) {
             cat(paste0(indent, "Population ", names(obj)[i], " (", length(obj[[i]]), " models):\n"))
             printPopulation(obj = obj[[i]], method = "digested", indent = paste(indent, "   "))
           }
         },
         
         # ---- Invalid Method Warning ----
         {
           warning("printModelCollection: please provide a valid method (short/long).")
         }
  )
}


#' Print Information about an Experiment Object
#'
#' This function prints details about an experiment object, including execution time,
#' parameters, dataset, and model collection.
#'
#' @param obj An experiment object created by `runExperiment()`.
#' @param indent A string for indentation to structure the output.
#'
#' @return None. The function prints experiment details to the console.
#' @export
printExperiment <- function(obj, indent = "\t--- ") {
  
  if (!isExperiment(obj)) {
    warning("printExperiment: The provided object is not a valid experiment.")
    return(NULL)
  }
  
  cat("\n========== Experiment Summary ==========\n")
  
  # Print execution time
  cat(paste(indent, "Execution Time (mins):", signif(obj$execTime, 2), "\n"))
  
  # Print Parameters with Proper Nesting
  printNestedList <- function(lst, indent) {
    for (name in names(lst)) {
      value <- lst[[name]]
      if (is.list(value)) {
        cat(paste0(indent, name, ":\n"))
        printNestedList(value, paste0(indent, "   "))  # Recursively print sublists
      } else {
        cat(paste0(indent, name, ": ", value, "\n"))
      }
    }
  }
  
  cat("\n========== Parameters ==========\n")
  printNestedList(obj$params, indent)
  
  # Print Rust-based Experiment Info
  cat("\n========== Rust Objects ==========\n")
  cat(paste(indent, "Experiment Pointer:", if (!is.null(obj$rust$experiment)) "Available" else "NULL", "\n"))
  cat(paste(indent, "Running Flag:", if (!is.null(obj$rust$running_flag)) "Active" else "NULL", "\n"))
  
  # Print Data Overview
  cat("\n========== Dataset ==========\n")
  if (!is.null(obj$data$train)) {
    cat(paste(indent, "Training Samples:", nrow(obj$data$train), "\n"))
    cat(paste(indent, "Training Features:", ncol(obj$data$train), "\n"))
  } else {
    cat(paste(indent, "Training Data: NULL\n"))
  }
  
  if (!is.null(obj$data$test)) {
    cat(paste(indent, "Test Samples:", nrow(obj$data$test), "\n"))
    cat(paste(indent, "Test Features:", ncol(obj$data$test), "\n"))
  } else {
    cat(paste(indent, "Test Data: NULL\n"))
  }
  
  # Print Model Collection
  cat("\n========== Model Collection ==========\n")
  if (isModelCollection(obj$model_collection)) {
    printModelCollection(obj$model_collection, method = "short")
  } else {
    cat(paste(indent, "No valid models found.\n"))
  }
}


#' Print Summary of Predomics Object
#'
#' This function prints a summary of a given object, identifying its type
#' (model, population, classifier, experiment, or model collection) and calling
#' the appropriate print function to display relevant information about the object.
#'
#' @param obj An object that can be of type model, population, classifier,
#'   experiment, or model collection.
#'
#' @return None. The function prints the summary of the object directly to the console.
#'
#' @details The function checks the type of the provided object using `isModel`,
#' `isPopulation`, `isClf`, `isExperiment`, and `isModelCollection` functions.
#' Based on the object type, it prints a summary:
#' - **Model**: Calls `printModel` with a detailed description of the model.
#' - **Population**: Calls `printPopulation`, showing a summary of the population of models.
#' - **Model Collection**: Calls `printModelCollection` to show a summary of a collection of models.
# #' - **Experiment**: Calls `printExperiment` to display experiment details.
# #' - **Classifier**: Calls `printClassifier` for classifier details.
#'
#' If the object type is not recognized, an error message is printed.
#'
#' @examples
#' \dontrun{
#' # Assuming 'model', 'population', 'classifier', 'experiment', and 'model_collection' are valid objects
#' printy(model)
#' printy(population)
#' printy(model_collection)
#' }
#'
#' @author Edi Prifti (IRD)
#' @export
printy <- function(obj) {
  
  if (isModel(obj)) {
    cat("Summary of Model object:\n")
    printModel(mod = obj, method = "short")
    
  } else if (isPopulation(obj)) {
    cat(sprintf("Summary of a population with %d models:\n", length(obj)))
    
    # Print details of up to 5 models
    printPopulation(obj = obj[1:min(5, length(obj))], method = "short")
    
    # Indicate that more models exist
    if (length(obj) > 5) cat("... and more models in the population\n")
    
  } else if (isModelCollection(obj)) {
    cat(sprintf("Summary of a Model Collection with %d populations:\n", length(obj)))
    printModelCollection(obj = obj, method = "short")
    
  } else if (isExperiment(obj)) {
    cat("Summary of Experiment object:\n")
    printExperiment(obj)
  } 
  else {
    stop("printy: The provided object is not a valid Predomics object.")
  }
}



#' Plot Feature Coefficients Across Models
#'
#' This function visualizes the coefficients of features across different models
#' in a heatmap. The colors represent positive, negative, and zero coefficients.
#' It supports sorting features, rotating axis labels, and using log scaling.
#'
#' @param feat.model.coeffs A numeric matrix or data frame where **rows** represent
#'   features, **columns** represent models, and values represent the feature coefficients.
#' @param topdown Logical; if `TRUE`, features are sorted from **top to bottom** (default: `TRUE`).
#' @param main A string for the **title** of the plot (default: `""`).
#' @param col A vector of **colors** for the heatmap (default: `c("deepskyblue1", "white", "firebrick1")`).
#' @param vertical.label Logical; if `TRUE`, feature labels are **rotated vertically** (default: `TRUE`).
#' @param log.scale Logical; if `TRUE`, **log transformation** is applied to the coefficients (default: `FALSE`).
#'
#' @return A `ggplot2` object displaying the heatmap of feature coefficients.
#'
#' @examples
#' # Example usage
#' features <- c("feature1", "feature2", "feature3")
#' models <- c("model1", "model2", "model3")
#' coeffs <- matrix(runif(9, -1, 1), nrow = 3, ncol = 3)
#' rownames(coeffs) <- features
#' colnames(coeffs) <- models
#'
#' # Plot feature model coefficients
#' plotFeatureModelCoeffs(coeffs, main = "Feature Coefficients Heatmap", log.scale = TRUE)
#'
#' @author Edi Prifti (IRD)
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plotFeatureModelCoeffs <- function(feat.model.coeffs, topdown = TRUE, main = "", 
                                   col = c("deepskyblue1", "white", "firebrick1"), 
                                   vertical.label = TRUE, log.scale = FALSE) {
  
  # Ensure input is a matrix
  if (!is.matrix(feat.model.coeffs)) {
    feat.model.coeffs <- as.matrix(feat.model.coeffs)
  }
  
  # Handle log scale transformation (avoid log(0))
  if (log.scale) {
    feat.model.coeffs <- sign(feat.model.coeffs) * log1p(abs(feat.model.coeffs))
  }
  
  # Sort features in top-down order if required
  if (topdown) {
    feat.model.coeffs <- feat.model.coeffs[nrow(feat.model.coeffs):1, , drop = FALSE]
  }
  
  # Convert matrix to long format for ggplot
  data.m <- reshape2::melt(feat.model.coeffs)
  colnames(data.m) <- c("feature", "model", "value")
  
  # Adjust color mapping based on available coefficient values
  unique_vals <- unique(data.m$value)
  if (length(unique_vals) < 3) {
    col <- col[c(-1, 0, 1) %in% unique_vals]
  }
  
  # Create heatmap plot
  p <- ggplot(data.m, aes(x = model, y = feature, fill = value)) + 
    geom_tile(color = "darkgray") + 
    scale_fill_gradientn(colors = col) +
    theme_bw() +
    ggtitle(main) +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 9),
      axis.text.x = if (vertical.label) element_text(angle = 90, hjust = 1) else element_text(angle = 0)
    )
  
  return(p)
}



#' Analyze Features in a Population of Models
#'
#' This function analyzes features in a population of models, allowing for the
#' visualization and examination of feature importance, prevalence, and model
#' coefficients. It can generate a variety of plots to understand the
#' distribution and importance of features in the given population.
#'
#' @param pop A population of models, typically obtained from
#'   `modelCollectionToPopulation` or similar functions.
#' @param X The data matrix containing features (rows represent features,
#'   columns represent samples).
#' @param y The response variable (class labels or continuous values depending
#'   on the model).
#' @param res_clf The classifier used for the analysis, typically a result from
#'   a classification experiment.
#' @param makeplot Logical. If `TRUE`, the function generates plots and saves
#'   them as a PDF. If `FALSE`, it returns the analysis results without
#'   plotting.
#' @param name A string representing the name of the analysis or output (used
#'   for saving files).
#' @param ord.feat A string indicating the ordering method for features. Options
#'   are:
#'   - "prevalence": Order by the prevalence of features across models.
#'   - "importance": Order by feature importance based on cross-validation.
#'   - "hierarchical": Order by hierarchical clustering of the feature-to-model coefficient matrix.
#' @param make.network Logical. If `TRUE`, generates a network of feature
#'   co-occurrence across the population of models.
#' @param network.layout A string indicating the layout of the network. Default
#'   is "circular". Other options may include "fr" for Fruchterman-Reingold
#'   layout.
#' @param network.alpha A numeric value controlling the alpha transparency of
#'   the network plot.
#' @param verbose Logical. If `TRUE`, prints additional information during
#'   execution.
#' @param pdf.dims A vector of two numbers specifying the width and height of
#'   the PDF output (in inches).
#' @param filter.perc A numeric value between 0 and 1 specifying the minimum
#'   prevalence of a feature to be included in the analysis.
#' @param k_penalty A penalty value for model selection in the population
#'   filtering.
#' @param k_max The maximum number of models to include in the final population
#'   after filtering.
#'
#' @return If `makeplot = TRUE`, returns a PDF with visualizations of feature
#'   importance, prevalence, and model coefficients. If `makeplot = FALSE`,
#'   returns a list of the analysis results including the normalized scores and
#'   feature importance.
#'
#' @details The function performs a variety of analyses on a population of
#' models:
#' - It filters models based on feature prevalence.
#' - It orders features by various metrics such as prevalence, importance, or hierarchical clustering.
#' - It generates plots of feature prevalence, model coefficients, and other characteristics.
#' - If requested, it also generates a network of feature co-occurrence across the models.
#'
#' @examples
#' \dontrun{
#' # Assuming 'pop' is a valid population of models, 'X' is the feature matrix, and 'y' is the response variable
#' analyzePopulationFeatures(pop = pop, X = X, y = y, res_clf = res_clf, makeplot = TRUE, name = "population_analysis")
#' }
#'
#' @author Edi Prifti (IRD)
#' @export
analyzePopulationFeatures <- function(pop, X, y, res_clf, makeplot = TRUE, name = "", ord.feat = "importance", 
                                      make.network = TRUE, network.layout = "circular", network.alpha = 0.0001, 
                                      verbose = TRUE, pdf.dims = c(width = 25, height = 20), filter.perc = 0.05, 
                                      k_penalty = 0.75/100, 
                                      k_max = 0)
{
  pop <- selectBestPopulation(pop, p = 0.05, k_penalty = k_penalty, k_max = k_max)
  if(verbose) print(paste("There are",length(pop), "models in this population"))
  
  if(length(pop) == 1)
  {
    print("analyzePopulationFeatures: only one model after filtering. Plot can not be built... returing empty handed.")
    return(NULL)
  }
  
  pop.df <- populationToDataFrame(pop)
  pop.noz <- listOfModelsToDenseCoefMatrix(clf = res_clf$classifier, X = X, y = y, list.models = pop)
  if(verbose) print(paste("Pop noz object is created with", nrow(pop.noz), "features and", ncol(pop.noz), "models"))
  
  # clean not conform bin/bininter models
  tocheck <- pop.df$language == "bin" | pop.df$language == "bininter"
  todelete <- apply(pop.noz[,tocheck] < 0, 2, any)
  
  if(any(todelete))
  {
    if(verbose) print(paste(sum(todelete)," bin/bininter models contain negative coefficients ... deleting them"))
  }
  
  if(length(todelete) != 0)
  {
    # clean population and recompute things
    pop <- pop[!todelete]  
  }
  
  # make the feature annots
  fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = res_clf$classifier)
  
  # filter features that are very rare in the models
  pop.noz <- filterFeaturesByPrevalence(X = fa$pop.noz, perc.prevalence = filter.perc * 100)
  
  if(!is.null(dim(pop.noz)))
  {
    if(nrow(pop.noz)==0)
    {
      pop.noz <- filterFeaturesByPrevalence(fa$pop.noz, perc.prevalence = 0)
    }  
  }
  if(verbose) print(paste("After filtering pop noz object is created with", nrow(pop.noz), "features and", ncol(pop.noz), "models"))
  
  if(verbose) print(paste("Ordering features by", ord.feat))
  switch(ord.feat,
         prevalence=
         {
           # Here we order based on the prevalence of features in the models
           prev <- getFeaturePrevalence(features = rownames(pop.noz), X = X, y = y, prop=TRUE)
           ind.feat <- order(prev$all, prev$`1`, prev$`-1`, decreasing = TRUE)
           features <- rownames(pop.noz)[ind.feat]
         },
         importance=
         {
           # Here we order based on a the crossval feature importance
           if(isExperiment(res_clf))
           {
             lr <- list(res_clf)
             names(lr) <- paste(res_clf$classifier$learner, res_clf$classifier$params$language, sep=".")
             
             feat.import <- mergeMeltImportanceCV(list.results = lr, 
                                                  filter.cv.prev = 0, 
                                                  min.kfold.nb = FALSE, 
                                                  learner.grep.pattern = "*", 
                                                  nb.top.features = NULL,
                                                  feature.selection = rownames(pop.noz),
                                                  scaled.importance = FALSE,
                                                  make.plot = TRUE)
           }
           
           if(is.null(feat.import))
           {
             print("analyzeImportanceFeatures: no feature importance data found... returning empty handed.")
             return(NULL)
           }
           
           # the most important features along with the order
           features <- rev(levels(feat.import$summary$feature))
           
         },
         hierarchical=
         {
           # Here we order based on a hierarchial clustering on the feature to model coefficient matrix
           hc.feat <- hclust(d = dist((pop.noz), method = "manhattan"), method = "ward.D")
           ind.feat <- hc.feat$order
           features <- rownames(pop.noz)[ind.feat]
         },
         {
           warning('This method does not exist !')
         }
  )
  
  if(ord.feat == "importance")
  {
    g6 <- feat.import$g
    
    g7 <- plotPrevalence(features = features, X, y)
    g7 <- g7 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    #plot abundance
    g8 <- plotAbundanceByClass(features = features, X, y)
    g8 <- g8 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    #tmp <- pop.noz; colnames(tmp) <- gsub("metal_","",colnames(pop.dense.noz))
    g9 <- plotFeatureModelCoeffs(feat.model.coeffs = pop.noz[features, ], 
                                 vertical.label = FALSE, col = c("deepskyblue1", "white", "firebrick1"))
    g9 <- g9 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    
    if(makeplot)
    {
      if(verbose) print(paste("Making plots in a dedicated pdf"))
      
      pdf(file=paste("population features",name,".pdf", sep=""), w=pdf.dims[1], h=pdf.dims[2])
      
      grid.arrange(g6, g9, g8, g7, ncol=4, widths = c(2,1,1,1))
      
      if(make.network)
      {
        if(verbose) print(paste("Making the network of co-occurance of features in the population of models"))
        #fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf)
        
        try(makeFeatureModelPrevalenceNetworkCooccur(pop.noz = pop.noz, 
                                                     feature.annot = fa$feature.df[rownames(pop.noz),], 
                                                     alpha = network.alpha, 
                                                     verbose = verbose, 
                                                     layout = network.layout),
            silent = TRUE)
      }
      dev.off()
      
    }else
    {
      return(grid.arrange(g6, g9, g8, g7, ncol=4, widths = c(2,1,1,1)))
    }
  }else
  {
    g7 <- plotPrevalence(features = features, X, y)
    g7 <- g7 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    #plot abundance
    g8 <- plotAbundanceByClass(features = features, X, y)
    g8 <- g8 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    #tmp <- pop.noz; colnames(tmp) <- gsub("metal_","",colnames(pop.dense.noz))
    g9 <- plotFeatureModelCoeffs(feat.model.coeffs = pop.noz[features, ], 
                                 vertical.label = FALSE, col = c("deepskyblue1", "white", "firebrick1"))
    
    if(makeplot)
    {
      if(verbose) print(paste("Making plots in a dedicated pdf"))
      
      pdf(file=paste("population features",name,".pdf", sep=""), width = pdf.dims[1], height = pdf.dims[2])
      
      grid.arrange(g9, g8, g7, ncol=3, widths = c(2,1,1))
      
      if(make.network)
      {
        if(verbose) print(paste("Making the network of co-occurance of features in the population of models"))
        #fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf)
        
        try(makeFeatureModelPrevalenceNetworkCooccur(pop.noz = pop.noz, 
                                                     feature.annot = fa$feature.df[rownames(pop.noz),], 
                                                     alpha = network.alpha, 
                                                     verbose = verbose, 
                                                     layout = network.layout),
            silent = TRUE)
      }
      dev.off()
      
    }else
    {
      return(grid.arrange(g9, g8, g7, ncol=3, widths = c(2,1,1)))
    }
  }
}



#' Prints as text the detail on a given experiment along with summarized results (if computed)
#'
#' @description This function will use the miic package to compute the co-occurance of features in a population of models
#' @param pop.noz: a data.frame of in features in the rows and models in the columns. 
#' This table contains the feature coefficients in the models and is obtained by makeFeatureAnnot()
#' @param feature.annot: a data frame with annotation on features obtained by makeFeatureAnnot()
#' @param cor.th: a threshold abtained on the partial correlation value
#' @param verbose: print out information during run
#' @param layout: the network layout by default is circular (layout_in_circle) and will be a weighted Fruchterman-Reingold otherwise
#' @return plots a graph
#' @export
makeFeatureModelPrevalenceNetworkMiic <- function(pop.noz, 
                                                  feature.annot, 
                                                  cor.th = 0.3, 
                                                  verbose = TRUE, 
                                                  layout = "circlular")
{
  require(igraph)
  require(miic)
  
  if(verbose) print("Creating the miic object")
  mobj <- miic(inputData = as.data.frame(t(pop.noz)))
  #miic.plot(mobj, igraphLayout = igraph::layout.fruchterman.reingold)
  
  #-----
  # load the edge information for spectral3off2 network
  edges <- mobj$retained.edges.summary
  if(verbose) print(paste(nrow(edges),"edges are found"))
  #dim(edges) # 350 edges
  colnames(edges)[1:2] <- c("from","to")
  rownames(edges) <- paste(edges$from, edges$to, sep=" => ")
  
  # clean those edges that have no link
  edges <- edges[abs(edges$partial_correlation) > cor.th,]
  if(verbose) print(paste(nrow(edges),"edges are kept after filtering by absolute correlation threshold", cor.th))
  
  
  # BUILD NETWORK
  #-------------------------------------------------------------------------------------------
  # ANNOTATION of the edges
  allnodes <- unique(c(edges$from,edges$to))
  edges.annot <- feature.annot[allnodes,] 
  edges.annot <- data.frame(rownames(edges.annot),edges.annot);
  colnames(edges.annot)[match("name",colnames(edges.annot))] <- "name_long"; 
  colnames(edges.annot)[1] <- "name"
  
  # create the igraph object
  gD <- graph.data.frame(d = edges, directed = TRUE, 
                         vertices = edges.annot)
  
  if(verbose) print(paste("The network is built"))
  
  # Calculate degree for all nodes
  degAll <- igraph::degree(gD, v = V(gD), mode = "all")
  gD <- set.vertex.attribute(gD, "degree", index = V(gD)$name, value = degAll)
  # Calculate betweenness for all nodes
  betAll <- igraph::betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
  betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll)); rm(betAll)
  # Add new node/edge attributes based on the calculated node properties/similarities
  gD <- set.vertex.attribute(gD, "betweenness", index = V(gD)$name, value = betAll.norm)
  
  # Calculate edge properties and add to the network
  E(gD)$color <- c("#DC143C","#A6A6A6")[as.factor(factor(sign(E(gD)$infOrt), levels=c('-1','1')))]
  E(gD)$infOrt[E(gD)$infOrt==  1] <- 0; 
  E(gD)$infOrt[E(gD)$infOrt== -2] <- 1
  
  #Calculate Dice similarities between all pairs of nodes
  dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
  # The following function will transform a square matrix to an edge driven one and add values to each edge
  F1 <- function(x) {data.frame(dice = dsAll[which(V(gD)$name == as.character(x$from)), which(V(gD)$name == as.character(x$to))])}
  # library(plyr) => took this out to force changing to tidyr
  edges.ext <- ddply(edges, .variables=c("from", "to"), function(x) data.frame(F1(x))); dim(edges.ext)
  
  gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
  E(gD)[as.character(edges.ext$from) %--% as.character(edges.ext$to)]$similarity <- as.numeric(edges.ext$dice)
  
  # Check the attributes
  summary(gD)
  
  if(layout == "circular")
  {
    l <- layout_in_circle(gD, order = V(gD))  
  }else
  {
    l <- layout_with_fr(gD, weights=E(gD)$weight)
  }
  
  plot(gD,
       vertex.label = V(gD)$name_long,
       vertex.color = c("deepskyblue1", "firebrick1")[factor(V(gD)$wilcox.class)],
       vertex.size=V(gD)$mod.prev/max(V(gD)$mod.prev)*10 + 3,
       edge.arrow.size=.2,
       asp=TRUE,
       rescale=TRUE,
       layout=l,
       edge.arrow.mode = E(gD)$infOrt,
       vertex.label.cex = 0.7,
       vertex.label.dist=0,
       edge.width = log(E(gD)$log_confidence)
  )
  
  #return(gD)
}


#' Plot Feature Abundance by Class
#'
#' Visualizes the abundance of selected features across different classes using **boxplots**. 
#' Supports both **classification** (discrete classes) and **regression** (continuous response).
#'
#' @param features Character vector of feature names to plot.
#' @param X Data matrix or data frame (features in rows, samples in columns).
#' @param y Vector of class labels (factor for classification) or numeric values (for regression).
#' @param topdown Logical; if `TRUE`, features are displayed from top to bottom (default: `TRUE`).
#' @param main Plot title (default: `""`).
#' @param plot Logical; if `TRUE`, displays the plot, otherwise returns statistical results (default: `TRUE`).
#' @param log_scale Logical; if `TRUE`, plots abundance data on a log10 scale (default: `FALSE`).
#' @param col.pt Colors for data points (default: `c("deepskyblue4", "firebrick4")`).
#' @param col.bg Colors for boxplot backgrounds (default: `c("deepskyblue1", "firebrick1")`).
#'
#' @return If `plot = TRUE`, returns a ggplot object. Otherwise, returns a data frame of statistical test results.
#' 
#' @examples
#' # Example: Classification
#' features <- c("feature1", "feature2")
#' X <- data.frame(feature1 = rnorm(100), feature2 = rnorm(100))
#' y <- sample(c(1, -1), 100, replace = TRUE)
#' plotAbundanceByClass(features, X, y, log_scale = TRUE)
#'
#' @import ggplot2
#' @import tidyr
#' @export
plotAbundance <- function(features, X, y, topdown = TRUE, 
                                 main = "", plot = TRUE, log_scale = FALSE,
                                 col.pt = c("deepskyblue4", "firebrick4"), 
                                 col.bg = c("deepskyblue1", "firebrick1")) {
  
  check.X_y(X, y)
  
  # Convert raw class if needed
  if (is.raw(y)) y <- as.integer(y)
  if (!is.vector(y)) y <- as.vector(y)
  
  # Ensure proper trait type
  if (is.numeric(y) && length(unique(y)) > 2) {
    mode <- "regression"
  } else {
    y <- as.factor(y)
    mode <- "classification"
  }
  
  if (any(is.na(match(features, rownames(X))))) {
    stop("plotAbundanceByClass: Some features are not found in the dataset.")
  }
  
  X <- as.matrix(X)
  dat <- X[features, , drop = FALSE]
  
  if (mode == "classification") {
    lev <- unique(y)
    datl1 <- X[features, y == lev[1], drop = FALSE]
    datl2 <- X[features, y == lev[2], drop = FALSE]
    
    dat.test <- filterfeaturesK(dat, y, k = nrow(dat), sort = FALSE)
    
    if (!plot) return(dat.test)
    
    qvals <- ifelse(dat.test$q < 0.05, "*", "")
    
    datl1_long <- as.data.frame(datl1) %>%
      tibble::rownames_to_column("feature") %>%
      pivot_longer(-feature, names_to = "observation", values_to = "abundance") %>%
      mutate(class = lev[1])
    
    datl2_long <- as.data.frame(datl2) %>%
      tibble::rownames_to_column("feature") %>%
      pivot_longer(-feature, names_to = "observation", values_to = "abundance") %>%
      mutate(class = lev[2])
    
    dat.reshape <- rbind(datl1_long, datl2_long)
    
  } else {  # Regression
    dat.test <- filterfeaturesK(dat, y, k = nrow(dat), sort = FALSE)
    
    if (!plot) return(dat.test)
    
    qvals <- ifelse(dat.test$q < 0.05, "*", "")
    
    dat.reshape <- as.data.frame(dat) %>%
      tibble::rownames_to_column("feature") %>%
      pivot_longer(-feature, names_to = "observation", values_to = "abundance") %>%
      mutate(class = "all")
  }
  
  dat.reshape$feature <- factor(dat.reshape$feature, levels = if (topdown) rev(features) else features)
  
  p <- ggplot(dat.reshape, aes(x = feature, y = abundance, fill = class, color = class)) +
    geom_boxplot() + coord_flip() +
    theme_bw() +
    scale_color_manual(values = col.pt) +
    scale_fill_manual(values = col.bg) +
    theme(legend.position = "none") +
    ggtitle(main)
  
  if (log_scale) {
    p <- p + scale_y_log10() + ylab("abundance (log10)")
  }
  
  p <- p + annotate("text", y = max(dat.reshape$abundance, na.rm = TRUE) * 1.1, 
                    x = 1:length(qvals), label = qvals, color = "gray", size = 7)
  
  return(p)
}


#' Plot Feature Coefficients Across Models
#'
#' This function visualizes the coefficients of features across different models
#' in a heatmap. It supports sorting features, rotating labels, and displaying 
#' coefficients in **log scale** when needed.
#'
#' @param feat.model.coeffs A numeric **matrix** where rows are **features**, 
#'   columns are **models**, and values represent the coefficients.
#' @param topdown Logical; if `TRUE`, features are sorted from **top to bottom** (default: `TRUE`).
#' @param main A string for the **title** of the plot (default: `""`).
#' @param col A vector of **colors** for the heatmap (default: `c("deepskyblue1", "white", "firebrick1")`).
#' @param vertical.label Logical; if `TRUE`, feature labels are **rotated vertically** (default: `TRUE`).
#' @param log.scale Logical; if `TRUE`, **log transformation** is applied to coefficients (default: `FALSE`).
#'
#' @return A `ggplot2` object displaying the heatmap of feature coefficients.
#'
#' @examples
#' # Example usage
#' features <- c("feature1", "feature2", "feature3")
#' models <- c("model1", "model2", "model3")
#' coeffs <- matrix(runif(9, -1, 1), nrow = 3, ncol = 3)
#' rownames(coeffs) <- features
#' colnames(coeffs) <- models
#'
#' # Plot feature model coefficients
#' plotFeatureModelCoeffs(coeffs, main = "Feature Coefficients Heatmap", log.scale = TRUE)
#'
#' @author Edi Prifti (IRD)
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plotFeatureModelCoeffs <- function(feat.model.coeffs, topdown = TRUE, main = "", 
                                   col = c("deepskyblue1", "white", "firebrick1"), 
                                   vertical.label = TRUE, log.scale = FALSE) {
  
  # Ensure input is a matrix
  if (!is.matrix(feat.model.coeffs)) {
    feat.model.coeffs <- as.matrix(feat.model.coeffs)
  }
  
  # Handle log scale transformation (avoid log(0))
  if (log.scale) {
    feat.model.coeffs <- sign(feat.model.coeffs) * log1p(abs(feat.model.coeffs))
  }
  
  # Sort features in top-down order if required
  if (topdown) {
    feat.model.coeffs <- feat.model.coeffs[nrow(feat.model.coeffs):1, , drop = FALSE]
  }
  
  # Convert matrix to long format for ggplot
  data.m <- reshape2::melt(feat.model.coeffs)
  colnames(data.m) <- c("feature", "model", "value")
  
  # Adjust color mapping based on available coefficient values
  unique_vals <- unique(data.m$value)
  if (length(unique_vals) < 3) {
    col <- col[c(-1, 0, 1) %in% unique_vals]
  }
  
  # Create heatmap plot
  p <- ggplot(data.m, aes(x = model, y = feature, fill = value)) + 
    geom_tile(color = "darkgray") + 
    theme_bw() +
    ggtitle(main) +
    scale_fill_gradientn(colors = col) +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 9),
      axis.text.x = if (vertical.label) element_text(angle = 90, hjust = 1) else element_text(angle = 0)
    )
  
  # Apply log scale on x-axis if enabled
  if (log.scale) {
    p <- p + scale_fill_viridis_c(option = "C", trans = "log10")  # Adjust fill colors for log
  }
  
  return(p)
}


#' Plot Feature Prevalence and Enrichment
#'
#' This function visualizes the prevalence of features across different groups
#' in a dataset and computes feature enrichment, providing optional statistical 
#' tests for enrichment. If enrichment data is available, significance markers are added.
#'
#' @param features A character vector of feature names to be plotted.
#' @param X A data matrix or data frame where **rows = features** and **columns = samples**.
#' @param y A vector of class labels (e.g., `1` and `-1` for binary classification) 
#'   corresponding to the columns in `X`.
#' @param topdown Logical; if `TRUE`, features are displayed in descending order (default: `TRUE`).
#' @param main A string for the title of the plot.
#' @param plot Logical; if `TRUE`, the function **displays the plot**. If `FALSE`, 
#'   it **returns the enrichment statistics** instead of plotting.
#' @param col.pt Colors for points in the plot (default: `c("deepskyblue4", "firebrick4")`).
#' @param col.bg Colors for bars in the plot (default: `c("deepskyblue1", "firebrick1")`).
#' @param zero.value The value treated as zero in prevalence calculations (default: `0`).
#'
#' @return If `plot = TRUE`, returns a `ggplot2` object displaying feature prevalence.  
#' If `plot = FALSE`, returns the **enrichment results** (statistical test output).
#'
#' @examples
#' features <- c("feature1", "feature2", "feature3")
#' X <- matrix(sample(0:1, 300, replace = TRUE), nrow = 3)
#' rownames(X) <- features
#' y <- sample(c(1, -1), 100, replace = TRUE)
#'
#' # Plot feature prevalence
#' plotPrevalence(features, X, y, main = "Feature Prevalence Plot")
#'
#' # Get enrichment statistics without plotting
#' plotPrevalence(features, X, y, plot = FALSE)
#'
#' @author Edi Prifti (IRD)
#' @export
plotPrevalence <- function(features, X, y, topdown = TRUE, main = "", plot = TRUE, 
                           col.pt = c("deepskyblue4", "firebrick4"), 
                           col.bg = c("deepskyblue1", "firebrick1"),
                           zero.value = 0) {
  
  # Compute feature prevalence
  v.prop <- getFeaturePrevalence(features = features, X = X, y = y, prop = TRUE, zero.value = zero.value)
  v.card <- getFeaturePrevalence(features = features, X = X, y = y, prop = FALSE, zero.value = zero.value)
  
  # Compute enrichment
  prev.enrichment <- computeCardEnrichment(v.card.mat = do.call(rbind, v.card), y = y)
  
  # Internal function to reshape data
  meltScoreList <- function(v, topdown = TRUE) {
    if (!is.list(v)) stop("meltScoreList: Input should be a list of vectors.")
    
    melted_data <- do.call(rbind, lapply(seq_along(v), function(i) {
      data.frame(feature = names(v[[i]]), 
                 prevalence = as.numeric(v[[i]]), 
                 group = names(v)[i], 
                 stringsAsFactors = FALSE)
    }))
    
    # Arrange feature order
    feature_order <- if (topdown) rev(names(v$all)) else names(v$all)
    melted_data$feature <- factor(melted_data$feature, levels = feature_order)
    
    # Invert prevalence for class "-1"
    melted_data$prevalence[melted_data$group == "-1"] <- 
      -melted_data$prevalence[melted_data$group == "-1"]
    
    return(melted_data)
  }
  
  # Prepare data for plotting
  v.prop.melt <- meltScoreList(v = v.prop, topdown = topdown)
  v.prop.melt$prevalence <- v.prop.melt$prevalence * 100  # Convert to percentage
  
  if (!plot) return(prev.enrichment)
  
  # Generate the plot
  p <- ggplot(v.prop.melt, aes(x = feature, y = prevalence, fill = group)) + 
    geom_bar(data = subset(v.prop.melt, group == "all"), stat = "identity") + 
    coord_flip() + 
    geom_point(data = subset(v.prop.melt, group %in% c("0", "1")), 
               aes(x = feature, y = prevalence, color = group, shape = group)) + 
    scale_color_manual(values = c("all" = "gray90", "0" = col.pt[1], "1" = col.pt[2])) +
    scale_fill_manual(values = c("all" = "gray90", "0" = col.bg[1], "1" = col.bg[2])) +
    scale_shape_manual(values = c(25, 24)) + 
    theme_bw() + 
    theme(legend.position = "none", axis.text = element_text(size = 9)) + 
    ggtitle(main)
  
  # Add significance markers if available
  if (!is.null(prev.enrichment$chisq.q)) {
    qvals <- ifelse(prev.enrichment$chisq.q < 0.05, "*", "")
    p <- p + annotate("text", y = rep(101, length(qvals)), 
                      x = seq_along(qvals) - 0.3, label = qvals, color = "gray", size = 7)
  }
  
  return(p)
}



#' Plot Barcode-style Heatmap of Feature Abundance
#'
#' This function generates a barcode-style heatmap from a feature-by-sample matrix,
#' grouped by class labels. It is designed for visualizing abundance or intensity values,
#' and supports both fixed and dynamic log-scaled color mapping.
#'
#' @param X A numeric matrix or data frame (features x samples).
#' @param y A vector of sample class labels corresponding to columns in `X`.
#' @param main A character string for the plot title.
#' @param ylabl (Unused; reserved for future use).
#' @param ylabr (Unused; reserved for future use).
#' @param fixed.scale Logical; if `TRUE`, uses a fixed log-scale color mapping.
#'   If `FALSE`, the color scale is computed dynamically from data range. Default is `TRUE`.
#' @param data Optional. A list containing `X`, `y`, and optionally `classes`.
#'   If provided, overrides `X` and `y`.
#' @param select_features Optional vector of feature names to display (subset of rownames of `X`).
#' @param select_samples Optional vector of sample names to display (subset of colnames of `X`).
#'
#' @return A `ggplot2` object containing the barcode heatmap.
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate filter
#' @importFrom scales log_trans label_number squish
#' @export
plotBarcode <- function(X = NULL, y = NULL, main = "", ylabl = "", ylabr = "",
                        fixed.scale = TRUE, data = NULL,
                        select_features = NULL, select_samples = NULL) {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(scales)
  
  # --- Load from data object if provided ---
  if (!is.null(data)) {
    if (!all(c("X", "y") %in% names(data))) {
      stop("Data object must contain at least 'X' and 'y'")
    }
    X <- data$X
    y <- data$y
    
    if (is.raw(y)) y <- as.integer(y)
    if (!is.vector(y)) y <- as.vector(y)
    y <- as.integer(y)
    y <- factor(y)
    
    if (!is.null(data$classes)) {
      class_labels <- as.character(data$classes)
      if (length(class_labels) < length(levels(y))) {
        stop("Not enough class labels provided for the unique values in y.")
      }
      levels(y) <- class_labels
    }
  } else {
    if (is.raw(y)) y <- as.integer(y)
    if (!is.vector(y)) y <- as.vector(y)
    y <- as.integer(y)
    y <- factor(y)
  }
  
  if (is.null(X) || is.null(y)) stop("Both X and y must be provided.")
  
  # --- Subset features and samples if requested ---
  if (!is.null(select_features)) {
    X <- X[rownames(X) %in% select_features, , drop = FALSE]
  }
  if (!is.null(select_samples)) {
    X <- X[, colnames(X) %in% select_samples, drop = FALSE]
    y <- y[colnames(X) %in% select_samples]
  }
  
  # --- Sanity check ---
  if (ncol(X) != length(y)) {
    stop("Length of y must match number of columns in X after selection.")
  }
  
  if (is.null(colnames(X))) colnames(X) <- paste0("Sample", seq_len(ncol(X)))
  if (is.null(rownames(X))) rownames(X) <- paste0("Feature", seq_len(nrow(X)))
  
  feature_levels <- rownames(X)
  sample_levels <- colnames(X)
  
  # --- Reshape data to long format ---
  df_long <- as.data.frame(X) %>%
    mutate(Feature = rownames(X)) %>%
    pivot_longer(cols = -Feature, names_to = "Sample", values_to = "value")
  
  df_long$Group <- y[match(df_long$Sample, colnames(X))]
  
  df_long <- df_long %>%
    filter(!is.na(Group)) %>%
    mutate(
      Feature = factor(Feature, levels = rev(feature_levels)),
      Sample = factor(Sample, levels = sample_levels),
      Group = factor(Group, levels = levels(y))
    )
  
  # Replace 0 with NA for white tiles
  df_long$value[df_long$value == 0] <- NA
  
  # --- Define color scale ---
  # --- Color scale ---
  if (fixed.scale) {
    breaks <- c(1e-07, 4e-07, 1.6e-06, 6.4e-06,
                2.56e-05, 0.0001024, 0.0004096, 0.0016384)
    
    fill_scale <- scale_fill_gradientn(
      colors = c("white", "deepskyblue", "blue", "green3", "yellow",
                 "orange", "red", "orangered2", "darkred"),
      trans = log_trans(base = 4),
      breaks = breaks,
      labels = scales::trans_format("log", math_format(10^.x)),
      limits = range(breaks),
      oob = squish,
      na.value = "white",
      guide = guide_colorbar(title = expression("Value (log"[4]*")"))
    )
  } else {
    non_zero_vals <- df_long$value[!is.na(df_long$value)]
    if (length(non_zero_vals) == 0) {
      stop("All values are zero; cannot compute dynamic log scale.")
    }
    
    min_val <- min(non_zero_vals)
    max_val <- max(non_zero_vals)
    log_min <- floor(log(min_val, base = 4))
    log_max <- ceiling(log(max_val, base = 4))
    breaks <- 4 ^ seq(log_min, log_max)
    
    fill_scale <- scale_fill_gradientn(
      colors = colorRampPalette(c("white", "blue", "green", "yellow", "red", "darkred"))(length(breaks) - 1),
      trans = log_trans(base = 4),
      breaks = breaks,
      labels = label_number(format = "e", accuracy = 1),
      limits = c(min_val, max_val),
      oob = squish,
      na.value = "white",
      guide = guide_colorbar(title = "Value (log₄)")
    )
  }
  
  # --- Build the plot ---
  p <- ggplot(df_long, aes(x = Sample, y = Feature, fill = value)) +
    geom_tile() +
    facet_wrap(~Group, scales = "free_x") +
    fill_scale +
    labs(title = main, x = "", y = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  return(p)
}

#' Plot a Heatmap of a Model Population with Signed Color Gradient
#'
#' This function extracts a model population from a `gpredomics` experiment object,
#' optionally selects the top-performing models, computes the dense coefficient matrix,
#' and visualizes it as a clustered heatmap. The color scale is centered on zero,
#' with negative values in red, positive values in green, and zero in white.
#'
#' @param exp A `gpredomics` experiment object containing the model collection and training data.
#' @param pop_index An integer specifying which population index to use from `exp$model_collection`
#'   (default is the last population).
#' @param focus_fbm Logical. If TRUE (default), selects the top models using `selectBestPopulation()`.
#'   If FALSE, uses the full population without filtering.
#' @param score A character string indicating which score to use for model selection (default is `"fit"`).
#' @param p A numeric value between 0 and 1 indicating the proportion of top models to retain (default is `0.05`).
#' @param clustering_distance Distance metric used for hierarchical clustering (`"euclidean"` by default).
#' @param clustering_method Method used for clustering (`"complete"` by default).
#' @param scale Scaling method passed to `pheatmap` (`"none"` by default).
#' @param ... Additional arguments passed to `pheatmap()` (e.g., `annotation_col`, `fontsize`, etc.).
#'
#' @return Invisibly returns the output of `pheatmap()`. The heatmap is displayed as a side effect.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Plot heatmap using filtered best models
#'   plot_population_heatmap(exp)
#'
#'   # Plot heatmap using all models from generation 5
#'   plot_population_heatmap(exp, pop_index = 5, focus_fbm = FALSE)
#' }
plot_population_heatmap <- function(exp,
                                    pop_index = length(exp$model_collection),
                                    focus_fbm = TRUE,
                                    score = "fit",
                                    p = 0.05,
                                    clustering_distance = "euclidean",
                                    clustering_method = "complete",
                                    scale = "none",
                                    ...) {
  library(pheatmap)
  
  # Step 1: Extract population
  if (pop_index < 1 || pop_index > length(exp$model_collection)) {
    stop("Invalid pop_index: must be between 1 and ", length(exp$model_collection))
  }
  
  pop <- exp$model_collection[[pop_index]]
  
  # Step 2: Optionally select best models
  list.models <- if (focus_fbm) {
    selectBestPopulation(pop, score = score, p = p)
  } else {
    pop
  }
  
  # Step 3: Compute dense coefficient matrix
  dense_matrix <- listOfModelsToDenseCoefMatrix(
    X = exp$data$train$X,
    y = exp$data$train$y,
    list.models = list.models
  )
  
  # Step 4: Check value range
  vmin <- min(dense_matrix, na.rm = TRUE)
  vmax <- max(dense_matrix, na.rm = TRUE)
  
  if (vmin == vmax) {
    warning("Matrix has constant values — slight jitter added to avoid error.")
    vmin <- vmin - 1e-6
    vmax <- vmax + 1e-6
  }
  
  # Step 5: Compute color breaks and palette
  n_colors <- max(20, min(100, nrow(dense_matrix) * ncol(dense_matrix)))
  neg_breaks <- seq(vmin, 0, length.out = ceiling(n_colors / 2))
  pos_breaks <- seq(0, vmax, length.out = floor(n_colors / 2) + 1)
  breaks <- unique(c(neg_breaks, pos_breaks[-1]))
  
  if (length(breaks) < 2) {
    stop("Breaks could not be computed correctly. Matrix may have insufficient variation.")
  }
  
  neg_colors <- colorRampPalette(c("darkred", "white"))(sum(breaks < 0))
  pos_colors <- colorRampPalette(c("white", "darkgreen"))(sum(breaks > 0))
  custom_colors <- c(neg_colors, pos_colors)
  
  # Step 6: Plot
  res <- pheatmap(dense_matrix,
                  color = custom_colors,
                  breaks = breaks,
                  clustering_distance_rows = clustering_distance,
                  clustering_distance_cols = clustering_distance,
                  clustering_method = clustering_method,
                  scale = scale,
                  ...)
  
  invisible(res)
}

#' Plot a precomputed genealogy
#'
#' @param td A list with components $nodes, $edges, $metadata as returned by Individual$get_genealogy_vis().
#' @param layout Character: "tree" (default), "kk", "fr", "dendrogram", "circle".
#' @param node_size_base Numeric base size for nodes (default 2).
#' @param text_repel Logical, use repel labels (default TRUE).
#' @return A ggraph plot object.
#' @export
plot_genealogy <- function(td, layout = "tree", node_size_base = 2, text_repel = TRUE) {
  if (is.null(td) || is.null(td$nodes) || is.null(td$edges)) {
    stop("`td` must be a list with $nodes and $edges as returned by get_genealogy_vis().")
  }
  nodes <- td$nodes
  edges <- td$edges

  if (!"name" %in% names(nodes)) {
    if ("id" %in% names(nodes)) {
      names(nodes)[names(nodes) == "id"] <- "name"
    } else {
      stop("Nodes need a 'name' or 'id' column.")
    }
  }
  nodes$name <- as.character(nodes$name)
  if ("label" %in% names(nodes)) nodes$label <- as.character(nodes$label)

  if (!all(c("from","to") %in% names(edges))) {
    stop("Edges need 'from' and 'to' columns.")
  }
  edges$from <- as.character(edges$from)
  edges$to   <- as.character(edges$to)

  used_ids <- unique(c(edges$from, edges$to))
  missing  <- setdiff(used_ids, nodes$name)
  if (length(missing) > 0) {
    add <- data.frame(
      name       = missing,
      label      = missing,
      auc        = NA_real_,
      k          = NA_integer_,
      generation = NA_integer_,
      stringsAsFactors = FALSE
    )
    for (cn in setdiff(names(add), names(nodes))) nodes[[cn]] <- NA
    for (cn in setdiff(names(nodes), names(add))) add[[cn]] <- NA
    nodes <- rbind(nodes[names(add)], add)
  }

  g <- igraph::graph_from_data_frame(
    d = edges[, c("from","to", setdiff(names(edges), c("from","to"))), drop = FALSE],
    vertices = nodes,
    directed = TRUE
  )

  lay <- match.arg(layout, c("tree","kk","fr","dendrogram","circle"))
  if ("k" %in% names(nodes)) {
    nodes$k[is.na(nodes$k)] <- 0L
  }

  p <- ggraph::ggraph(g, layout = lay) +
    ggraph::geom_edge_link(
      arrow   = grid::arrow(length = grid::unit(2, "mm")),
      end_cap = ggraph::circle(3, "mm"),
      alpha   = 0.5, colour = "grey55"
    ) +
    ggraph::geom_node_point(ggplot2::aes(color = .data$auc, size = .data$k)) +
    ggplot2::scale_color_viridis_c(name = "AUC", option = "C", na.value = "grey80") +
    ggplot2::scale_size_continuous(name = "k", range = c(node_size_base, node_size_base * 3)) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::labs(
      title = paste("Genealogy", if (!is.null(td$metadata$root_id)) paste("of", td$metadata$root_id) else ""),
      x = NULL, y = NULL
    )
    
  if ("label" %in% names(nodes)) {
    if (text_repel) {
      p <- p + ggraph::geom_node_text(
        ggplot2::aes(label = .data$label),
        repel = TRUE, size = 3, bg.color = "white", segment.color = "grey60"
      )
    } else {
      p <- p + ggraph::geom_node_text(ggplot2::aes(label = .data$label), size = 3, vjust = -1)
    }
  }

  p <- p + ggplot2::guides(
    color = ggplot2::guide_colorbar(barheight = 6, barwidth = 0.4),
    size  = ggplot2::guide_legend(override.aes = list(color = "grey30"))
  )

  return(p)
}