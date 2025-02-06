
#---------------------------------------------
# Recompile and load libraries
#---------------------------------------------
rextendr::document()


library(gpredomicsR)

# # other libraries
library(ggplot2)
library(tidyverse)
library(GGally)
library(momr)

#---------------------------------------------
# Create a new experiment
#---------------------------------------------
glog <- GLogger$new()
running_flag <- RunningFlag$new()

setwd("sample")
param <- Param$load("param.yaml")
# Get the parameters in an r object
rparam <- param$get()

system.time(expRust <- ga(param, running_flag))
# exp$generation_number()
# get a simple version of the individual
# (ind <- exp$get_individual(99,0))
# get a full individual
# (ind <- exp$get_individual_full(generation = 99,order = 0, verbose = FALSE))


#---------------------------------------------
# Data extraction
#---------------------------------------------
train_data = expRust$get_data_robj(train = TRUE)
test_data = expRust$get_data_robj(train = FALSE)
# check the data barcodes
momr::plotBarcode(train_data$X/1e3)
momr::plotBarcode(test_data$X/1e3)
# load the annotation data
annot <- read_tsv("hs_10_4_1990_MSP.20240626.gtdb_r220_annot_long.tsv")


#---------------------------------------------
# Data extraction
#---------------------------------------------
train_data = expRust$get_data_robj(train = TRUE)
test_data = expRust$get_data_robj(train = FALSE)
momr::plotBarcode(train_data$X/1000)
momr::plotBarcode(test_data$X/1000)


#---------------------------------------------
# Analyze the experiment
#---------------------------------------------
# get all the individuals as they have evolved, this is an R object. exp on the other hand is a rust pointer.
exp <- parseExperiment(expRust)
# get a specific population, for instance the latest
pop <- expRust$get_generation(0) 
# or get it like this 
pop <- exp$gen_1
pop.final <- expRust$get_generation(expRust$generation_number()-1)

# check the number of individuals per generation
barplot(unlist(lapply(exp,length)), las = 2, names.arg = 1:length(exp), col = "lightblue", 
        main = "Number of individuals per generation", xlab = "Generation", ylab = "Number of individuals")

table(populationGet_X("language")(pop))
# Binary Ternary 
# 499     500 
table(populationGet_X("data_type")(pop))
# Log Raw 
# 500 499


# the best models
perf.df <- analyzeAttributeEvolution(exp, attributes = c("auc", "fit", "k", "epoch","specificity","sensitivity","accuracy"), plot = FALSE)
# perf.df <- analyzeAttributeEvolution(exp, attributes = c("auc", "fit", "k", "epoch"), plot = FALSE)
analyzeAttributeEvolution(exp, attributes = c("fit"), plot = TRUE, best_model = FALSE)
# 
# # all the models
# analyzeAttributeEvolution(exp, attributes = c("auc", "fit", "k", "epoch","specificity","sensitivity","accuracy"), plot = TRUE)
# 


#---------------------------------------------
# The family of the best models (FBM)
#---------------------------------------------
pop.fbm <- selectBestPopulation(pop.final, score = "fit", p = 0.05)
df <- populationToDataFrame(pop.fbm)

df %>% 
  ggplot(aes(x = as.character(k), y = fit, fill = language)) +
  geom_boxplot() +
  facet_grid(~data_type, scales = "free_y") +
  labs(title = "Boxplot of the best models", x = "Language", y = "Fit") +
  theme_bw()


plot(table(df$k))
plot(table(df$data_type, df$language))

GGally::ggpairs(df, columns = c("k", "fit", "data_type", "sensitivity", "specificity", "accuracy"), 
                aes(color = language)) + 
  scale_fill_manual(values = c("orange", "firebrick")) + 
  scale_color_manual(values = c("orange", "firebrick")) +
  theme_bw()


populationGet_X("indexes")(pop.fbm)

# get the feature model prevalence matrix
dense_matrix <- listOfModelsToDenseCoefMatrix(X = train_data$X, y = train_data$y, pop.fbm)
print(dense_matrix)





#---------------------------------------------
# Overall operations
#---------------------------------------------



# glog$set_log_file("test.log")
# glog$set_log_level("INFO")
