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