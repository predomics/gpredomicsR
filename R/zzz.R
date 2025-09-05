# .onLoad <- function(libname, pkgname) {
#   try(GLogger$new(), silent = TRUE)
# }

.onLoad <- function(libname, pkgname) {
  try({
    .gp_logger_set(GLogger$new())
  }, silent = TRUE)
}