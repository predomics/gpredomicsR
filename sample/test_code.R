library(gpredomicsR)

glog <- GLogger$new()
running_flag <- RunningFlag$new()

setwd("sample")

param <- Param$get("param.yaml")
pop <- ga(param, running_flag)
pop$generation_number()
ind <- pop$get_individual(99,0)
ind <- pop$get_individual_full(99,0)

pop.final <- pop$get_all_individuals(99)


# glog$set_log_file("test.log")
# glog$set_log_level("INFO")

individual <- list()
individual$names  <- c("f1","f2") # vecteur
individual$coefs <- c() # vecteur
individual$auc <- c() # valeur
