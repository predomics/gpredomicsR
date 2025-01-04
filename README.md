# Compiling

Open an R session at root and:
```R
devtools::build(".")
install.packages("../gpredomicsR_0.0.0.9000.tar.gz", repos = NULL, type = "source")
```

# Using

```R
library(gpredomicsR)
running_flag <- RunningFlag$new()
param <- Param$get("param.yaml")
pop <- ga(param, running_flag)
pop$generation_number()
pop$get_individual(99,0)
```

NB: you can additionally setup a GLogger object to add some display of info during ga() call
```R
glog <- GLogger$new()
```
NB2: the glog object can only be created once.
NB3: once created the verbosity level can be adjusted with : `glog$set_level("warning")` for instance.

*NEW*

- you may directly set a logger at a correct level `glog$level("trace")`
- you may use `ga/ga2/ga_no_overfit/ga2_no_overfit` sub algorithms of the `ga` family
- you may get all individuals at once in a a generation with `pop$get_all_individuals(98)` (which will give you a Dataframe of all the individuals) 