#' @useDynLib gpredomicsR, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL

# ============================
# RunningFlag S3 Wrapper
# ============================

#' Create a new RunningFlag.
#'
#' A RunningFlag tracks whether the GA is supposed to keep running.
#' @export
new_runningflag <- function() {
  rf <- .Call("gpredomicsR_RunningFlag__new")
  class(rf) <- "RunningFlag"
  rf
}

#' Stop the RunningFlag.
#'
#' Sets the internal flag to false.
#' @param rf A RunningFlag object
#' @export
stop_running <- function(rf) {
  .Call("gpredomicsR_RunningFlag__stop", rf)
  invisible(rf)
}

#' Check whether the RunningFlag is currently true/false.
#'
#' @param rf A RunningFlag object
#' @return `TRUE` if still running, `FALSE` otherwise
#' @export
is_running <- function(rf) {
  .Call("gpredomicsR_RunningFlag__is_running", rf)
}

#' Reset the RunningFlag to true.
#'
#' @param rf A RunningFlag object
#' @export
reset_running <- function(rf) {
  .Call("gpredomicsR_RunningFlag__reset", rf)
  invisible(rf)
}

#' Print method for RunningFlag.
#'
#' Shows whether it's currently running.
#'
#' @param x A RunningFlag object
#' @param ... unused
#' @export
print.RunningFlag <- function(x, ...) {
  cat("RunningFlag:", is_running(x), "\n")
  invisible(x)
}

# ============================
# Param S3 Wrapper
# ============================
# (You have an #[extendr] struct Param in Rust, but the module block didn't list `impl Param;`
#  so make sure that appears in your extendr_module! if you want these .Call() bindings.)

#' Create a new Param
#'
#' Creates an empty or default gpredomics::Param object
#' @export
new_param <- function() {
  p <- .Call("gpredomicsR_Param__new")
  class(p) <- "Param"
  p
}

#' Load a Param from a file
#'
#' Calls Param::get(file_path) in Rust, which may panic if the file is unsuitable.
#' @param file_path Path to a parameter file
#' @export
get_param <- function(file_path) {
  p <- .Call("gpredomicsR_Param__get", file_path)
  class(p) <- "Param"
  p
}

#' Set feature minimal prevalence percentage
#'
#' @param param A Param object
#' @param pct A numeric value
#' @export
set_feature_minimal_prevalence_pct <- function(param, pct) {
  .Call("gpredomicsR_Param__set_feature_minimal_prevalence_pct", param, pct)
  invisible(param)
}

#' Print method for Param
#'
#' (Param doesn't expose a custom `to_string` in Rust, so we just label it.)
#'
#' @param x A Param object
#' @param ... unused
#' @export
print.Param <- function(x, ...) {
  cat("A gpredomics Param object\n")
  invisible(x)
}

# ============================
# Population S3 Wrapper
# ============================

#' Run the GA
#'
#' This calls the Rust ga(param, running_flag) function, returning a Population.
#'
#' @param param A Param object
#' @param running_flag A RunningFlag object
#' @return A Population S3 object
#' @export
ga <- function(param, running_flag) {
  pop <- .Call("gpredomicsR_ga", param, running_flag)
  class(pop) <- "Population"
  pop
}

#' Get individuals from a Population
#'
#' @param pop A Population object
#' @param generation Which generation index (0-based)
#' @param order Which individual index (0-based)
#' @return A named integer vector indicating which features are present (and what value)
#' @export
get_individuals <- function(pop, generation, order) {
  .Call("gpredomicsR_Population__get_individuals", pop, generation, order)
}

#' Number of generations
#'
#' @param pop A Population object
#' @return integer number of generations
#' @export
generation_number <- function(pop) {
  .Call("gpredomicsR_Population__generation_number", pop)
}

#' Population size
#'
#' @param pop A Population object
#' @param generation Which generation index
#' @return integer population size for that generation
#' @export
population_size <- function(pop, generation) {
  .Call("gpredomicsR_Population__population_size", pop, generation)
}

#' Print method for Population
#'
#' @param x A Population object
#' @param ... ignored
#' @export
print.Population <- function(x, ...) {
  cat("A Population object with", generation_number(x), "generations.\n")
  invisible(x)
}
