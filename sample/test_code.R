
# recompile
rextendr::document()


library(gpredomicsR)

# # other libraries
library(ggplot2)
library(tidyverse)


#---------------------------------------------
# Create a new experiment
#---------------------------------------------
glog <- GLogger$new()
running_flag <- RunningFlag$new()

setwd("sample")
param <- Param$load("param.yaml")

exp <- ga(param, running_flag)
# exp$generation_number()
# get a simple version of the individual
# (ind <- exp$get_individual(99,0))
# get a full individual
# (ind <- exp$get_individual_full(generation = 99,order = 0, verbose = FALSE))

# get all the individuals as they have evolved
generations <- parseExperiment(exp)

# get a specific population, for instance the latest
pop <- exp$get_generation(99)
# pop <- generations[[length(generations)]]


#---------------------------------------------
# Analyze the experiment
#---------------------------------------------
perf.df <- analyzeEvolution(exp, attributes = c("auc", "fit", "k", "n"), plot = FALSE)
analyzeEvolution(exp, attributes = c("auc", "fit", "k", "n"), plot = TRUE)
analyzeEvolutionAllModels(exp, attributes = c("auc", "fit", "k", "n"), plot = TRUE)
  




pop.final <- pop$get_all_individuals(99)


# glog$set_log_file("test.log")
# glog$set_log_level("INFO")

individual <- list()
individual$names  <- c("f1","f2") # vecteur
individual$coefs <- c() # vecteur
individual$auc <- c() # valeur
