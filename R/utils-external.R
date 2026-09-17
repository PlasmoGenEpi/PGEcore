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

#' Check several suggested packages
#'
#' @param pkgs Character vector of package names.
#' @param reason Optional short description of why they are needed.
#' @return Invisibly returns `TRUE` if all packages are available.
#' @keywords internal
check_suggested_pkgs <- function(pkgs, reason = NULL) {
  for (pkg in pkgs) {
    check_suggested_pkg(pkg, reason)
  }
  invisible(TRUE)
}

#' Require variantstring 1.x
#'
#' @param reason Passed to [check_suggested_pkg()].
#' @return Invisibly returns `TRUE` if a 1.x version is installed.
#' @keywords internal
check_variantstring_v1 <- function(reason = NULL) {
  check_suggested_pkg("variantstring", reason)
  ver <- as.character(utils::packageVersion("variantstring"))
  if (utils::compareVersion(ver, "1.0.0") < 0 ||
      utils::compareVersion(ver, "2.0.0") >= 0) {
    stop(
      "This functionality requires variantstring version 1.x.x, but version ",
      ver,
      " is installed.",
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
