
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

# Run experiment
exp <- runExperiment(param.path = "sample/param.yaml")
# saveRDS(exp, paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_exp.rds"))

#---------------------------------------------
# Data extraction
#---------------------------------------------
train_data = exp$rust$experiment$get_data_robj(train = TRUE)
test_data = exp$rust$experiment$get_data_robj(train = FALSE)
# check the data barcodes
plotBarcode(train_data$X/1e3)
plotBarcode(test_data$X/1e3)
# load the annotation data
annot <- read_tsv("sample/hs_10_4_1990_MSP.20240626.gtdb_r220_annot_long.tsv")


#---------------------------------------------
# Analyze the experiment
#---------------------------------------------
pop.final <- exp$model_collection[[length(exp$model_collection)]] # the final population

# check models
mod <- pop.final[[1]]
printy(mod)
modrust <- exp$rust$experiment$individual(exp$rust$experiment$generation_number()-1,0)$get()
printy(modrust)




# TODO make a function that analyses the overall experiment evolution
# check the number of individuals per generation
barplot(unlist(lapply(exp$model_collection,length)), las = 2, names.arg = 1:length(exp$model_collection), col = "lightblue", 
        main = "Number of individuals per generation", xlab = "Generation", ylab = "Number of individuals")


table(populationGet_X("language")(pop.final))
# Binary Ternary 
# 390     514 
table(populationGet_X("data_type")(pop.final))
# Log Prevalence        Raw 
# 396        243        265 


# the best models
perf.df <- analyzeAttributeEvolution(exp, attributes = c("auc", "fit", "k", "epoch","specificity","sensitivity","accuracy"), plot = FALSE)
# perf.df <- analyzeAttributeEvolution(exp, attributes = c("auc", "fit", "k", "epoch"), plot = FALSE)
analyzeAttributeEvolution(exp$model_collection, attributes = c("fit"), plot = TRUE, best_model = TRUE)
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

ggpairs(df, columns = c("k", "fit", "data_type", "accuracy"), 
                aes(color = language)) + 
  scale_fill_manual(values = c("orange", "firebrick")) + 
  scale_color_manual(values = c("orange", "firebrick")) +
  theme_bw()


indexes <- unique(populationGet_X("indexes")(pop.fbm))

# get the feature model prevalence matrix
dense_matrix <- listOfModelsToDenseCoefMatrix(X = train_data$X, y = train_data$y, pop.fbm)
print(dense_matrix)

mod <- pop.fbm[[1]]
printy(mod)
printy(pop.fbm)
printy(exp)


plotAbundanceByClass(features = rownames(dense_matrix), X = train_data$X, y = train_data$y, log_scale = FALSE)



#---------------------------------------------
# Overall operations
#---------------------------------------------



# glog$set_log_file("test.log")
# glog$set_log_level("INFO")
