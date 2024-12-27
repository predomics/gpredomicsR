# gpredomicsR

## Install

- Compile with:
```
R CMD build .
```

- Install with:
```
R CMD INSTALL gpredomicsR_0.1.0.tar.gz
```
(you'll need Rcpp package if you do not have it, it can be installed with `install.packages("Rcpp")` in an R session)

## Usage

```R
# Load the Rust package
library(gpredomicsR)

# Create a new RunningFlag instance
running_flag <- RunningFlag$new()

# Start the Rust function in a separate R thread
result <- parallel::mcparallel({
    run_ga_no_overfit("param.yaml", running_flag$get_arc())
})

# Check if the process is running
running_flag$is_running()

# Stop the process
running_flag$stop()

# Wait for the thread to complete
parallel::mccollect(result)
```

```R
library(gpredomicsR)

# 1. RunningFlag
rf <- new_runningflag()
rf    # prints RunningFlag: TRUE
stop_running(rf)
is_running(rf)   # FALSE
reset_running(rf)
rf    # RunningFlag: TRUE

# 2. Param
p1 <- new_param()
p2 <- get_param("path/to/param_file.toml")
set_feature_minimal_prevalence_pct(p2, 0.05)
p2   # prints "A gpredomics Param object"

# 3. GA run => Population
pop <- ga(p2, rf)
pop                # "A Population object with X generations."
get_individuals(pop, 0, 0)  # named int vector for the 0th individual's features
population_size(pop, 0)   
```