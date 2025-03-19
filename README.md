# Preliminary

[gpredomicsR](https://github.com/predomics/gpredomicsR) is an R package that wraps 
the [gpredomics](https://github.com/predomics/gpredomics) rust library, which 
implements accelerated heuristics to find sparse interpretable models described in 
[Prifti et al, Gigascience, 2020](https://academic.oup.com/gigascience/article/9/3/giaa010/5801229). 
The main algorithm is a genetic algorithm which applies numerous operators to evolve `individuals`, 
which are classification models. The fit of the models is evaluated on a penalized AUC.

# Compiling

Open an R session at root and:
```R
devtools::build(".")
install.packages("../gpredomicsR_0.0.0.9000.tar.gz", repos = NULL, type = "source")
```

# Using

## Creating the environment
After loading the package, we need to load the parameters and create a `RunningFlag` object. 
`param` is a rust pointer and is need to call rust internal functions. 
However, it is possible to get it as an R object with `param$get()`.

```R
library(gpredomicsR)
running_flag <- RunningFlag$new()
param <- Param$load("param.yaml")
rparam <- param$get()
```

## Running the algorithm
The main heuristic based on genetic algorithm `ga()` is called with the `param` and `running_flag` objects.
This will return a `Experiment` object which will store the results of the algorithm.
The results in R are a list generations, each being a `Population` object, which contains the individuals of each generation.
Several functions allow to get specific information, such as the number of generations, the best individual of a generation, etc.

```R
# Create an experiment and run the genetic algorithm
exp <- fit(param, running_flag)
exp$generation_number()
exp$individual(99,0)
exp.results <- parseExperiment(exp)
```

## Analyzing the results
```R
# get the resulting attributes from all the evolve modesl
perf.df <- analyzeEvolution(exp, attributes = c("auc", "fit", "k", "n"), plot = FALSE)
analyzeEvolution(exp, attributes = c("auc", "fit", "k", "n"), plot = TRUE)
analyzeEvolutionAllModels(exp, attributes = c("auc", "fit", "k", "n"), plot = TRUE)
  
```

NB: you can additionally setup a GLogger object to add some display of info during ga() call
```R
glog <- GLogger$new()
```
NB2: the glog object can only be created once.
NB3: once created the verbosity level can be adjusted with : `glog$set_level("warning")` for instance.

*NEW*

- you may directly set a logger at a correct level `glog$set_level("trace")`
- you may use `ga/ga2/ga_no_overfit/ga2_no_overfit` sub algorithms of the `ga` family
- you may get all individuals at once in a a generation with `pop$get_all_individuals(98)` (which will give you a Dataframe of all the individuals) 

## GPU usage

You should have GPU set to true (directly in param.yaml or with `param$set_bool("gpu",TRUE)`). You may need additional libraries (Apple Silicons are natively supported by WGPU which gpredomics is using, Linux needs Vulkan libraries, with drivers, see gpredomics README.md for details).