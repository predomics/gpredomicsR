.gpredomics_env <- new.env(parent = emptyenv())

# ------------------------------------------------------------
# Log level management (R-side only)
# Levels: none < error < warn < info < debug
# ------------------------------------------------------------
.gpredomics_env$log_level <- "info"
.gpredomics_env$level_rank <- c(none = 0L, error = 1L, warn = 2L, info = 3L, debug = 4L)

#' Set the R-side log level
#' @param level one of c("none","error","warn","info","debug")
#' @keywords internal
gp_set_log_level <- function(level = c("info", "warn", "error", "debug", "none")) {
  level <- match.arg(level)
  .gpredomics_env$log_level <- level
  invisible(level)
}

#' Get the current R-side log level
#' @keywords internal
gp_get_log_level <- function() {
  lvl <- .gpredomics_env$log_level
  if (is.null(lvl) || !nzchar(lvl)) lvl <- "info"
  lvl
}

.gp_level_allowed <- function(level) {
  ranks <- .gpredomics_env$level_rank
  cur <- gp_get_log_level()
  (ranks[[level]] %||% 0L) <= (ranks[[cur]] %||% 3L)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------
# Optional storage for a GLogger instance (not required)
# ------------------------------------------------------------
.gp_logger_set <- function(gl = NULL) {
  assign("glogger", gl, envir = .gpredomics_env)
  invisible(gl)
}

.gp_logger_get <- function() {
  if (exists("glogger", envir = .gpredomics_env, inherits = FALSE)) {
    get("glogger", envir = .gpredomics_env, inherits = FALSE)
  } else {
    NULL
  }
}