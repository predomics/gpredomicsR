#' Get the models from a classifier result for each k-sparsity
#'
#' @description Get the N best models from a classifier result for each k-sparsity.
#' @param obj: the classifier result output from the function fit. This can also be a ModelCollection or Population object
#' @param significance: if TRUE, (default:FALSE) a statistical test will be applied to find the lowest threshold that will delimit the window
#' of the best models. If FALSE, the models will be selected according to the rest of the criteria.
#' @param by.k.sparsity: if TRUE (default:TRUE), the filtering will be performed for each sparsity level
#' @param k.penalty: (default:0), it will penalize the models with large sparsity if different, when by.k.sparsity is set to TRUE
#' @param n.best: the number of best models to be returned for each sparsity if by.k.sparsity is set to TRUE or for the whole population 
#' otherwise (default:5).
#' @param nbest: the number of best models we wish to get from the population, per each sparsity or not. If there are less best models then this
#' number, less will be returned
#' @param single.best: if TRUE, this will return the best model of all (default:FALSE) and the n.best will be set to 1. 
#' @param single.best.cv: if single.best is TRUE, we could chose the best model based on data from cross validation (default:TRUE) and in this 
#' case obj should be an experiment or from empirical results not in CV.
#' @param single.best.k: if single.best is TRUE, we could chose the best model of a given sparsity that is specified by a number here. 
#' If this value is specified (default:NULL), then this will de-actvate single.best.cv.
#' @param max.min.prevalence: if TRUE (default:FALSE), the best models will be selected based on their performance but also on the prevalence of 
#' the features that compose it.
#' @param X: the dataset to be learned (default:NULL). This is neeeded when max.min.prevalence is set to TRUE.
#' @param verbose: provide more information about execution (default = FALSE)
#' @param evalToOrder: which attribute of the model object should we use to order the models and select them (default:fit_)
#' @param return.population: if set to TRUE (default:FALSE), the result will be send as a population of models
#' @param unique.control: if set to TRUE (default:TRUZ), we correct the population so that no dupplication of models takes place
#' @return a list of model objects or a model when it is a single one or a model collection
#' @export
getNBestModels <- function(obj, 
                           significance = FALSE, 
                           by.k.sparsity = TRUE,
                           k.penalty = 0,
                           n.best = 5,
                           single.best = FALSE,
                           single.best.cv = TRUE,
                           single.best.k = NULL,
                           max.min.prevalence = FALSE,
                           X = NULL,
                           verbose = FALSE, 
                           evalToOrder = "fit_",
                           return.population = FALSE,
                           unique.control = TRUE
)
{
  
  # sanity check
  if(!isExperiment(obj) & !isModelCollection(obj) & !isPopulation(obj))
  {
    warning("getNBestModels: please provide a valid experiment, modelCollection or population object ... returning NULL")
    return(NULL)
  }
  
  # if an experiment
  if(isExperiment(obj = obj)) 
  {
    if(verbose) print(paste0("getNBestModels: the object is an experiment"))
    mc              <- obj$classifier$models
  }
  
  # convert to a model collection if it is not
  if(isPopulation(obj = obj)) # if a population
  {
    if(verbose) print(paste0("getNBestModels: the object is an population of models"))
    mc              <- listOfModels2ModelCollection(pop = obj)
  }
  # if a modelCollection
  if(isModelCollection(obj = obj)) 
  {
    if(verbose) print(paste0("getNBestModels: the object is a model collection"))
    mc              <- obj
  }
  
  if(!isModelCollection(mc))
  {
    warning("getNBestModels: the object is not a valid model collection. Returning empty handed.")
    return(NULL)
  }
  
  if(single.best) 
  {
    n.best          <- 1
    if(verbose) print(paste0("getNBestModels: single best"))
  }
  
  # set up switch variables that are no needed to parameterize but that we could in the future
  penalty_by_kfold <- FALSE
  
  if(by.k.sparsity | single.best)
  {
    ####################################################################
    # # if by.k.sparsity
    ####################################################################
    
    if(verbose) print(paste0("getNBestModels: by k sparsity"))
    res <- list()
    pop.valids <- c()
    # for each k_sparsity
    for(i in 1:length(mc))
    {
      if(unique.control)
      {
        pop           <- unique(mc[[i]])
      }else
      {
        pop           <- mc[[i]]  
      }
      
      if(verbose) print(paste0("getNBestModels: the population of sparsity ", i," contains ", length(pop), " models"))
      
      if(significance)
      {
        # restrict the population to the confidence interval
        pop         <- selectBestPopulation(pop = pop, score = evalToOrder, p = 0.05, k_penalty = k.penalty)
        if(verbose) print(paste0("getNBestModels: after significance selection with k.penalty ", k.penalty," it contains ", length(pop), " models"))
      }
      
      # if we wish to select best models with max min prevalence
      if(max.min.prevalence)
      {
        if(!is.null(X))
        {
          eval      <- populationGet_X(element2get = evalToOrder, toVec = TRUE, na.rm = FALSE)(pop)
          best.eval <- max(eval, na.rm = T)
          # epsilon is used to be able to select a window of best models
          epsilon   <- sqrt(best.eval * (1 - best.eval) / ncol(X))
          pop       <- pop[eval >= (best.eval - epsilon) & !is.na(eval)]
          pop       <- getMaxMinPrevalenceModel(pop = pop, X = X, selected = 0)
          if(verbose) print(paste0("getNBestModels: after max.min.prevalence it contains ", length(pop), " models"))
        }
      }
      
      pop           <- sortPopulation(pop = pop, evalToOrder = evalToOrder)
      pop           <- pop[1:min(n.best, length(pop))]
      if(verbose) print(paste0("getNBestModels: the final population contains ", length(pop), " models"))
      res[[i]]      <- pop
      # mark valididity
      pop.valids <- c(pop.valids, isPopulation(pop))
      
    } # end for loop
    names(pop.valids) <- names(mc)
    names(res)      <- names(mc)[pop.valids]
    
    mc <- mc[pop.valids]
    
    if(!isModelCollection(mc))
    {
      warning("digestModelCollection: after treating the mc object no result is available. Returning NULL")
      return(NULL)
    }
    
    if(single.best)
    {
      single.best.cv <- FALSE
      # set best model type switch
      if(isExperiment(obj) & !myAssertNotNullNorNa(single.best.k))
      {
        if(!is.null(obj$crossVal))
        {
          single.best.cv <- TRUE  
        }
      }
      
      k_catalogue <- paste0("k_",obj$classifier$params$sparsity)
      spar        <- populationGet_X("eval.sparsity", toVec = TRUE, na.rm = FALSE)(modelCollectionToPopulation(res))
      
      # Here we are left with two options
      if(single.best.cv)
      {
        # get the best model from crossval information
        if(verbose) print(paste0("getNBestModels: single.best.cv mode ... returning the best model"))
        if(obj$classifier$params$objective == "auc" & !(evalToOrder == "accuracy_" | evalToOrder == "auc_"))
        {
          evalToOrder <- "accuracy_"
        }
        if(obj$classifier$params$objective == "cor" & !(evalToOrder == "cor_"))
        {
          evalToOrder <- "cor_"
        }
        
        # Abreviations for the results
        key <- data.frame(real=c("auc_","accuracy_","cor_"), abbrv=c("auc","acc","cor")); rownames(key) <- key$real
        emp.name <- paste("empirical", as.character(key[evalToOrder,]$abbrv), sep=".")
        gen.name <- paste("generalization", as.character(key[evalToOrder,]$abbrv), sep=".")
        # for each classifier
        emp.data <- obj$crossVal$scores[[emp.name]][k_catalogue, ]
        gen.data <- obj$crossVal$scores[[gen.name]][k_catalogue, ]
        # plot for debug
        # par(mfrow=c(2,1)); image(as.matrix(t(emp.data))); image(as.matrix(t(gen.data)))
        
        if(!is.null(emp.data) & !is.null(gen.data))
        {
          # compute the penalty data
          emp.data.penalty          <- emp.data
          k                         <- as.numeric(gsub("k_","",k_catalogue))
          
          
          # if we want to compute the penalty by k_fold
          if(penalty_by_kfold)
          {
            for(j in 1:nrow(emp.data))
            {
              emp.data.penalty[j,]  <- emp.data[j,] - penalty[i] * k[j]
            }
            # select the k_sparse for each k_fold
            ind <- apply(emp.data.penalty, 2, which.max)
            k_sparse <- rownames(emp.data.penalty)[ind]
            
            best_empirical <- diag(as.matrix(emp.data[ind,]))
            best_generalization <- diag(as.matrix(gen.data[ind,]))
            
          }else # otherwise we compute a meaned penalty
          {
            mean_score <- rowMeans(emp.data.penalty, na.rm = TRUE)
            mean_score_penalty  <- mean_score - k.penalty * k
            #plot(mean_score, mean_score_penalty, ylim=c(0,1),xlim=c(0.5,1))
            
            # make sure to be in the space of available sparsity models in the emperical models
            ind <- which.max(mean_score_penalty[names(mc)])
            k_sparse <- rep(rownames(emp.data.penalty[names(mc),])[ind], ncol(emp.data))
            best_empirical <- as.numeric(emp.data[names(mc),][ind,])
            best_generalization <- as.numeric(gen.data[names(mc),][ind,])
            
            # if no values are found in empirical
            if(length(ind) == 0)
            {
              best_empirical <- logical(0)
              best_generalization <- logical(0)
            }
            
            # => TEST if(all(is.na(best_empirical))) best_empirical <- logical(0)
            # => TEST if(all(is.na(best_generalization))) best_generalization <- logical(0)
            # plot(best_empirical, best_generalization, ylim=c(0.5,1),xlim=c(0.5,1))
            # boxplot(list(best_empirical,best_generalization), ylim=ylim)
          }
          res.k.cv <- data.frame(best_empirical, best_generalization, k_sparse)
        }
        else
        {
          res.k.cv <- data.frame(best_empirical = NA, best_generalization = NA)[-1,]
        }
        best.k <- as.numeric(gsub("k_","",as.character(unique(res.k.cv$k_sparse))))
        
      }else
      {
        # get the best model from empirical information
        if(verbose) print(paste0("getNBestModels: single.best mode ... returning the best model"))
        
        eval <- NULL
        for(i in 1:length(res))
        {
          eval.i <- populationGet_X(evalToOrder)(res[[i]])
          if(length(eval.i) > 0)
          {
            eval <- c(eval, eval.i)
          }else
          {
            eval <- c(eval, NA)
          }
        }
        # eval <- unlist(lapply(res, function(x){populationGet_X(evalToOrder)(x)[[1]]}))
        # apply the k_penalty
        eval <- eval - (k.penalty * spar)
        
        best.k <- as.numeric(gsub("k_","",names(eval[which.max(eval)])))
      }
      
      # select the best model for a given k
      if(!myAssertNotNullNorNa(single.best.k))
      {
        # set from above if it does not exist
        single.best.k <- best.k
      }else
      {
        if(single.best.k == 0)
        {
          # this is a special case, and means that there will not be a selection but the maximum number of variables will be taken into account
          #k   <- as.numeric(gsub("k_","",k_catalogue))
          single.best.k <- max(spar)
        }
      }
      
      k_spar <- paste0("k_", single.best.k)
      
      # otherwise when single.best.k is specified this will be the first choice
      if(length(single.best.k) == 1 & is.numeric(single.best.k))
      {
        if(k_spar %in% names(mc)) # check if we have it in the model collection
        {
          if(verbose) print(paste0("getNBestModels: single.best.k mode with k_spar ", k_spar, " returning the best model"))
          pop     <- mc[[k_spar]]
          eval    <- populationGet_X(element2get = evalToOrder, toVec = TRUE, na.rm = FALSE)(pop)
          mod     <- pop[[which.max(eval)]]
          if(return.population)
          {
            return(list(mod))
          }else
          {
            return(mod)  
          }
        }else # not found
        {
          if(verbose) print(paste0("getNBestModels: single.best.k mode with k_spar ", k_spar, " not found in the results"))
          if(return.population)
          {
            return(list(NA))
          }else
          {
            return(NA)  
          }
        }
      }else
      {
        print(paste0("getNBestModels: single.best.k mode but ",k_spar, " is not found in the model collection. Executing the default settings."))
      }
    } # end of single.best.k 
    
    if(return.population)
    {
      if(verbose) print(paste0("getNBestModels: returning a population of models"))
      # Transform the model collection onto a population
      return(modelCollectionToPopulation(mod.collection = res))
    }else
    {
      if(verbose) print(paste0("getNBestModels: returning a model collection"))
      # a model collection
      return(res)
    }
    
  }else 
  {
    ####################################################################
    # # else not by.k.sparsity
    ####################################################################
    if(verbose) print(paste0("getNBestModels: regardless of k sparsity"))
    # first convert the pop
    if(unique.control)
    {
      pop             <- unique(modelCollectionToPopulation(mc))
    }else
    {
      pop             <- modelCollectionToPopulation(mc)  
    }
    
    if(verbose) print(paste0("getNBestModels: the population with all sparsities contains ", length(pop), " models"))
    
    if(significance)
    {
      # restrict the population to the confidence interval
      pop           <- selectBestPopulation(pop = pop, score = evalToOrder, p = 0.05, k_penalty = k.penalty)
      if(verbose) print(paste0("getNBestModels: after significance selection with k.penalty ", k.penalty," it contains ", length(pop), " models"))
    }
    
    # if we wish to select best models with max min prevalence
    if(max.min.prevalence)
    {
      if(!is.null(X))
      {
        eval        <- populationGet_X(element2get = evalToOrder, toVec = TRUE, na.rm = FALSE)(pop)
        k           <- populationGet_X(element2get = "eval.sparsity", toVec = TRUE, na.rm = TRUE)(pop)
        eval.penalty<- eval - (k*k.penalty)
        best.eval   <- max(eval.penalty)
        epsilon     <- sqrt(best.eval * (1 - best.eval) / ncol(X))
        pop         <- pop[eval.penalty >= (best.eval - epsilon)]
        pop         <- getMaxMinPrevalenceModel(pop = pop, X = X, selected = 0)
        if(verbose) print(paste0("getNBestModels: after max.min.prevalence it contains ", length(pop), " models"))
      }
    } # end max.min.prevalence
    
  } # end by.k.sparsity ifelse
  
  
  if(return.population)
  {
    if(verbose) print(paste0("getNBestModels: returning a population of models"))
    return(pop)
  }else
  {
    if(verbose) print(paste0("getNBestModels: returning a model collection"))
    return(listOfModels2ModelCollection(pop))
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

