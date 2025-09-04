
#---------------------------------------------
# Recompile and load libraries
#---------------------------------------------
rextendr::document()


library(gpredomicsR)

# # other libraries
library(tidyverse)
library(GGally)
# library(momr)

#---------------------------------------------
# Create a new experiment
#---------------------------------------------

# Run experiment
system.time(exp <- runExperiment(param.path = "sample/param.yaml", name = "test_experiment", glog_level = "info"))
# saveRDS(exp, paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_exp.rds"))

# Save the experiment
exp$rust$experiment$save("sample/test_experiment.mp")

exp2 <- load_experiment("sample/test_experiment.mp")

tmp <- exp$rust$experiment$get_data(train = TRUE)$get()
Experiment$get_data()

#---------------------------------------------
# Data extraction
#---------------------------------------------
train_data = exp$rust$experiment$get_data_robj(train = TRUE)
test_data = exp$rust$experiment$get_data_robj(train = FALSE)

# # check the data barcodes
# plotBarcode(X = train_data$X, y = as.numeric(train_data$y), fixed.scale = FALSE)
# plotBarcode(X = train_data$X[1:30,1:10], y = train_data$y[1:10], fixed.scale = FALSE)

plotBarcode(data = train_data, fixed.scale = FALSE)
plotBarcode(data = test_data, fixed.scale = FALSE)
plotBarcode(data = train_data, select_features = c("msp_0937",  "msp_0938",  "msp_0939", "msp_0069", "msp_0005"), fixed.scale = FALSE)
# plotBarcode(data = test_data, fixed.scale = TRUE)

#---------------------------------------------
# Jury extraction
#---------------------------------------------
jury <- exp$rust$experiment$get_jury()
rust_train <- exp$rust$experiment$get_data(train = TRUE)
rust_test <- exp$rust$experiment$get_data(train = FALSE)

# Evaluate on Train Data
jury$evaluate(rust_train)

# Get R object (data is simply used to obtain the names of the features and does not influence the jury object)
jury$get(rust_train) 
jury$get(rust_test)

# Get Jury metrics on a specific dataset
jury$compute_new_metrics(rust_train)
jury$compute_new_metrics(rust_test)

# Get Jury class and scores on a specific dataset
jury$evaluate_class_and_score(rust_train)
jury$evaluate_class_and_score(rust_test)

# Print Jury
jury$display_train(rust_train, exp$rust$param)
jury$display_train_and_test(rust_train, rust_test, exp$rust$param)

# load the annotation data
annot <- read_tsv("sample/hs_10_4_1990_MSP.20240626.gtdb_r220_annot_long.tsv")


#---------------------------------------------
# Analyze the experiment
#---------------------------------------------
# extract the final evolved population
pop.final <- exp$model_collection[[length(exp$model_collection)]] # the final population

# check models
mod <- pop.final[[1]]
printy(mod)
# # extract directly the same model from the rust object
# modrust <- exp$rust$experiment$individual(exp$rust$experiment$generation_number()-1,0)$get()
# printy(modrust)


# TODO make a function that analyses the overall experiment evolution
# check the number of individuals per generation
barplot(unlist(lapply(exp$model_collection,length)), las = 2, names.arg = 1:length(exp$model_collection), col = "lightblue", 
        main = "Number of individuals per generation", xlab = "Generation", ylab = "Number of individuals")


table(populationGet_X("language")(pop.final))
# Binary    Pow2   Ratio Ternary 
# 1133    1172    1374    1183
table(populationGet_X("data_type")(pop.final))
# Log Prevalence        Raw 
# 1563       1787       1512 

plot(table(populationGet_X("data_type")(pop.final),populationGet_X("language")(pop.final)))
     

#---------------------------------------------
# Recompute all the statistics for all the populations
#---------------------------------------------


#---------------------------------------------
# Add the FBM selection to analyzeAttributeEvolution
#---------------------------------------------


# the best models
perf.df <- analyzeAttributeEvolution(exp, attributes = c("auc", "fit", "k", "epoch","specificity","sensitivity","accuracy"), plot = FALSE)
# perf.df <- analyzeAttributeEvolution(exp, attributes = c("auc", "fit", "k", "epoch"), plot = FALSE)
analyzeAttributeEvolution(exp, attributes = c("fit"), plot = TRUE, best_model = TRUE)
 
# all the models
analyzeAttributeEvolution(exp, attributes = c("auc", "fit", "k", "epoch","specificity","sensitivity","accuracy"), plot = TRUE)
 


#---------------------------------------------
# The family of the best models (FBM)
#---------------------------------------------
# plot the FMB coefficients heatmap
# plot_population_heatmap(exp, pop_index = 99, focus_fbm = FALSE)
# plot_population_heatmap(exp, pop_index = 99, cluster_rows = FALSE)


pop.fbm <- selectBestPopulation(pop.final, score = "fit", p = 0.05)
df <- populationToDataFrame(pop.fbm)
# get the feature model prevalence matrix
dense_matrix <- listOfModelsToDenseCoefMatrix(X = exp$data$train$X, y = exp$data$train$y, pop.fbm)
rownames(dense_matrix)



# plot the abundance of the selected features
plotAbundance(features = rownames(dense_matrix), X = exp$data$train$X, y = exp$data$train$y, log_scale = TRUE)
plotAbundance(features = rownames(dense_matrix), X = exp$data$train$X, y = exp$data$train$y, log_scale = FALSE)
# plot the prevalence of the selected features
plotPrevalence(features = rownames(dense_matrix), X = exp$data$train$X, y = exp$data$train$y)
# plot the heatmap of the coefficients of the FBM
plot_population_heatmap(exp, pop_index = 99, cluster_rows = FALSE)
# plot the barcode of the selected features
plotBarcode(data = train_data, select_features = rownames(dense_matrix), fixed.scale = FALSE)


# check annotation
View(annot[annot$msp_name %in% rownames(dense_matrix),])


get_taxonomy(msp_names = rownames(dense_matrix), gtdb_version = 220.0, fields = c("phylum","class","order","family","genus","species"))

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
  scale_fill_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") +
  theme_bw()


indexes <- unique(populationGet_X("indexes")(pop.fbm))


#---------------------------------------------
# The best models voting
#---------------------------------------------








#---------------------------------------------
# Focus on one specific model
#---------------------------------------------
mod <- pop.fbm[[1]]
printy(mod)
printy(pop.fbm)
printy(exp)

gmod <- exp$rust$experiment$individual(exp$rust$experiment$generation_number()-1,0)$get()

printy(mod)

computeClass <- function(X, y, mod){
  # compute the class of the model
  # X: the data
  # y: the labels
  # model: the model
  # return: the class of the model
  y_pred <- predict(mod, X)
  y_pred <- ifelse(y_pred > 0.5, 1, 0)
  return(y_pred)
}


#---------------------------------------------
# Overall operations
#---------------------------------------------



# glog$set_log_file("test.log")
# glog$set_log_level("INFO")
