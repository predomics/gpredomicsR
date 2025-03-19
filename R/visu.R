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
           if (all(!is.null(mod$coeff), length(mod$coeff) > 0)) {
             ind.pos <- mod$coeff > 0
             ind.neg <- mod$coeff < 0
             
             term.pos <- if (any(ind.pos)) paste0("+", sort(mod$indexes[ind.pos]), collapse = " ") else "0"
             term.neg <- if (any(ind.neg)) paste0("+", sort(mod$indexes[ind.neg]), collapse = " ") else "0"
             
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
plotAbundanceByClass <- function(features, X, y, topdown = TRUE, 
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
    p <- p + scale_y_log10()
  }
  
  p <- p + annotate("text", y = max(dat.reshape$abundance, na.rm = TRUE) * 1.1, 
                    x = 1:length(qvals), label = qvals, color = "gray", size = 7)
  
  return(p)
}
