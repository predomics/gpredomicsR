# .onLoad <- function(libname, pkgname) {
#   try(GLogger$new(), silent = TRUE)
# }

.onLoad <- function(libname, pkgname) {
  try({
    .gp_logger_set(GLogger$new())
  }, silent = TRUE)
  registerS3method("print", "Data", print.Data)
  registerS3method("print", "Individual", print.Individual)
  registerS3method("print", "Population", print.Population)
  registerS3method("print", "Jury", print.Jury)
  registerS3method("print", "Experiment", print.Experiment)
  registerS3method("print", "Param", print.Param)
  
  op <- options()
  op.gpredomics <- list(
    gpredomics.threads.number = 4
  )
  toset <- !(names(op.gpredomics) %in% names(op))
  if (any(toset)) options(op.gpredomics[toset])
  invisible()
}

#' @export
print.Data <- function(x, ...) {
  cat(x$print(), "\n")
  invisible(x)
}

#' @export
print.Individual <- function(x, ...) {
  x$print()
  invisible(x)
}

#' @export
print.Population <- function(x, ...) {
  cat(x$print(), "\n")
  invisible(x)
}

#' @export
print.Jury <- function(x, ...) {
  cat(x$print(), "\n")
  invisible(x)
}

#' @export
print.Experiment <- function(x, ...) {
  cat(x$print(), "\n")
  invisible(x)
}

#' @export
print.Param <- function(x, ...) {
  cat(x$address(), "\n")
  print(x$get())
  invisible(x)
}