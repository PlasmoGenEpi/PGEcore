#' Check that a suggested package is installed
#'
#' Used by wrappers that depend on optional R packages listed in Suggests.
#' Does not install the package.
#'
#' @param pkg Name of the package.
#' @param reason Optional short description of why it is needed.
#' @return Invisibly returns `TRUE` if the package is available.
#' @keywords internal
check_suggested_pkg <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    why <- if (is.null(reason)) {
      sprintf("Package '%s' is required for this functionality.", pkg)
    } else {
      sprintf("Package '%s' is required for %s.", pkg, reason)
    }
    stop(
      why,
      " Install it separately (for example via Conda or install.packages).",
      " PGEcore does not install specialised dependencies automatically.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Check that an external executable is on PATH
#'
#' External tools are not R package dependencies and are not installed by
#' PGEcore. Workflows must provide them in the process environment.
#'
#' @param tool Executable name as found by [Sys.which()].
#' @return Invisibly returns the absolute path to the executable.
#' @keywords internal
check_external_tool <- function(tool) {
  path <- Sys.which(tool)
  if (identical(unname(path), "") || identical(path, "")) {
    stop(
      sprintf(
        paste0(
          "Could not find executable '%s' on PATH. ",
          "PGEcore does not install external tools; install '%s' ",
          "separately and ensure it is available on PATH."
        ),
        tool,
        tool
      ),
      call. = FALSE
    )
  }
  invisible(path)
}
