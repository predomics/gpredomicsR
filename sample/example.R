#---------------------------------------------
# Recompile and load libraries
#---------------------------------------------
rextendr::document()

library(readr)
library(tibble)
library(dplyr)
options(gpredomics.threads.number = 8)

#---------------------------------------------
# Loading data
#---------------------------------------------

# X should be a data frame
df_X_train <- read_table("sample/Xtrain.tsv")
df_X_train <- df_X_train %>% column_to_rownames("msp_name")

# y an ordered factor with names=samples
df_y_train <- read_table("sample/Ytrain.tsv", col_names = FALSE, skip=1)
y_train <- factor(ifelse(df_y_train[,2]==0, "healthy", "cirrhosis"), levels=c("healthy", "cirrhosis"), ordered = T)
names(y_train) <- df_y_train[,1][[1]]

gpred_train <- as.gpredomics.data(
  X = df_X_train, 
  y = y_train, 
  features.in.columns = FALSE
)

# Same with test set
df_X_test <-  read_table("sample/Xtest.tsv")
df_X_test <- df_X_test %>% column_to_rownames("msp_name")
df_y_test <- read_table("sample/Ytest.tsv", col_names = FALSE, skip=1)
y_test <- factor(ifelse(df_y_test[,2]==0, "healthy", "cirrhosis"), levels=c("healthy", "cirrhosis"), ordered = T)
names(y_test) <- df_y_test[,1][[1]]
gpred_test <- as.gpredomics.data(
  X = df_X_test, 
  y = y_test, 
  features.in.columns = FALSE
)

#---------------------------------------------
# Launching the experiment
#---------------------------------------------

# Loading param and creating the Running flag
param <- Param$load("sample/param.yaml")
running <- RunningFlag$new()

# Launch the experiment
exp <- fit_on(gpred_train, param, running)

#---------------------------------------------
# Get experiment results
#---------------------------------------------

# Get the best population (last generation)
pop <- exp$get_best_population()

# Print the number of generation
exp$generation_number()

# Get a specific generation
exp$get_generation(20)

#---------------------------------------------
# Dealing with individuals
#---------------------------------------------

# Get the best model
best <- pop$get_individual(1)

# Get the best model metrics on training set
best$get_metrics()

# Get the best model metrics on other set
best$compute_metrics(gpred_test)

# Make prediction with the best model
best$predict(gpred_train)$class
best$predict(gpred_test)$class

# Compute the feature importance (MDA)
imp <- best$compute_importance(gpred_train, n_perm=1000, seed=4815162342, used_only=TRUE)

# Prune Individual based on a MDA threshold
pruned_best <- best$prune_by_threshold(0, n_perm=1000, seed=4815162342, min_k = 10)

# Better AUC, worst MCC
best$compute_metrics(gpred_test)$auc
pruned_best$compute_metrics(gpred_test)$auc

best$compute_metrics(gpred_test)$mcc
pruned_best$compute_metrics(gpred_test)$mcc

# Prune Individual based on a MDA quantile
pruned_best_2 <- best$prune_by_quantile(quantile=0.05, eps=0, n_perm=1000, seed=4815162342, min_k = 10)

# Create a Population from a vector/list of Individuals
ex_pop <- Population$from_individuals(c(best, pruned_best))

# Fit Individual on new data / on new param
fit_mcc <- param
fit_mcc$set_string("fit", "mcc")
best$fit(gpred_train, param)

plotBarcode(data = gpred_train$get(), select_features = best$get()$features, fixed.scale = TRUE)

# Study the Individual genealogy
gen <- best$get_genealogy(exp, max_depth = 4)
plot_genealogy(gen, node_vars = list(color="language"))

# Explain a specific prediction with waterfall plot
# Show feature contributions for first sample
plotIndividualWaterfall(best, gpred_train, sample = 1)

# Show contributions for a specific sample with custom title
plotIndividualWaterfall(best, gpred_train, sample = 5, 
                        main = "Feature contributions for sample 5")

# Show only top 10 features by absolute contribution
plotIndividualWaterfall(best, gpred_test, sample = 3, top_n = 10)

#---------------------------------------------
# Dealing with population
#---------------------------------------------

# Filter population by metrics
pop$filter_by_auc(0.999)
pop$filter_by_specificity(0.999)
pop$filter_by_sensitivity(0.999)

# Filter according to a 1/0 mask
# For example, to remove ratio:
mask <- as.integer(sapply(pop$get()$individuals, function(x) { x["language"]})!="Ratio")
pop$filter_by_mask(mask)

# Or to remove models including the feature msp_0005:
mask <- as.integer(sapply(pop$get()$individuals, function(x) { !"msp_0005" %in% x$features }))
pop$filter_by_mask(mask)

# Filter population by diversity 
# by_niche represents a filter by common model type:
# filter linear (bin/ter/pow2) & filter ratio apart
pop$filter_by_diversity(10, by_niche=FALSE) 

# Get the Family of Best Models
fbm <- pop$get_fbm(0.05)

# Fit FBM on new param or data 
fbm$fit(gpred_train, param)

#---------------------------------------------
# Dealing with Jury
#---------------------------------------------
jury <- Jury$from_population(fbm, 0.0, 0.0)

# Get Jury metrics on a specific dataset
jury$compute_metrics(gpred_train)
jury$compute_metrics(gpred_test)

# Get Jury class and scores on a specific dataset
jury$predict(gpred_train)$class
jury$predict(gpred_test)$class

# Print Jury report (Gpredomics binary style)
jury$print_self_report()
jury$print_report(gpred_test)

# Non-gpredomics method can also be used to compute meta-model from a Population
# To do that, extract each Individual scores or predictions 
jury_pop <- jury$get_population()

jury_pop$predict_score_matrix(gpred_train)
jury_pop$predict_class_matrix(gpred_train)

pimp <- jury_pop$get_first_pct(10)$compute_importance_matrix(gpred_train, 100, T, 4815162342)
plotImportanceHeatmap(pimp)

#---------------------------------------------
# Dealing with Cross-validation 
#---------------------------------------------

cv_param <- Param$load("sample/param.yaml")
cv_param$set_bool("cv", TRUE)
cv_param$set("outer_folds", 3)

exp_cv <- fit_on(gpred_train, cv_param, running)

# Get the best population (fold merged FBMs)
exp_cv$get_best_population()

# Get the fold numbers
exp_cv$get_n_folds()

# Get a fold data
exp_cv$get_fold_data(1, train=TRUE) # train
exp_cv$get_fold_data(1, train=FALSE) # valid

# Get gen 100 of the first fold with metrics computed on train (k-1)
exp_cv$get_fold_generation(1, 100, train=TRUE)  

# Get CV importance
exp_cv$compute_cv_importance(n_perm = 100, aggregation = "median", scaled = TRUE, seed = 4815162342, compact = FALSE)

#---------------------------------------------
# Get R object
#---------------------------------------------

param$get()
gpred_train$get()
best$get()
pop$get()
jury$get()
#exp$get() #also is possible but huge

#---------------------------------------------
# Deal with features annotations 
#---------------------------------------------

msps <- rownames(df_X_train)

weights <-rep(1, length(msps))

names(weights) <- msps
weights[1:1980] <- 1

penalties <-rep(0, length(msps))
names(penalties) <- msps
penalties["msp_0007"] <- 100000000

tags <- data.frame(feature=msps, name=paste(msps, "_name", sep=""))
samples <- data.frame(sample=colnames(df_X_train), test=rep(c("A", "B"), ncol(df_X_train)/2), test2=rep(c("C", "D"), ncol(df_X_train)/2))

gpred_train <- as.gpredomics.data(
  X = df_X_train, 
  y = y_train, 
  features.in.columns = FALSE,
  prior.weight = weights,
  feature.penalty = penalties,
  feature.tags = tags,
  sample.tags = samples
)
                      
# Loading param and creating the Running flag
param <- Param$load("sample/param.yaml")
running <- RunningFlag$new()

# Modifying a parameter
param$set("max_epochs", 100)
param$set("user_penalties_weight", 1)

# Launch the experiment
exp <- fit_on(gpred_train, param, running)



# Control the application of the penalty
pop <- exp$get_best_population()$get_first_pct(10)
msp_seven_count <- 0
for (ind in pop$get()$individuals) {
  if ("msp_0007" %in% ind$features) {
    msp_seven_count <- msp_seven_count+1
  }
}

# Severely penalized feature is not represented in our 10% best models
print(msp_seven_count)
