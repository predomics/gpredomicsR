#' Create a new Individual object
#' @export
create_individual <- function() {
  ptr <- .Call("rust_create_individual")
  new_individual(ptr)
}

#' Evaluate an Individual
#' @param individual An Individual object
#' @param data A Data object (assuming Data object is defined similarly)
#' @export
evaluate_individual <- function(individual, data) {
  scores <- .Call("rust_evaluate_individual", individual$ptr, data$ptr)
  scores # Assuming it returns a vector of scores
}

#' Compute AUC for an Individual
#' @param individual An Individual object
#' @param data A Data object
#' @export
compute_auc_individual <- function(individual, data) {
  auc <- .Call("rust_compute_auc_individual", individual$ptr, data$ptr)
  auc
}