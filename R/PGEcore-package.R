#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data
#' @importFrom optparse OptionParser make_option parse_args
## usethis namespace: end
NULL

# NSE symbols used in tidyr/dplyr pipelines (e.g. run_coiaf)
utils::globalVariables(c(
  "specimen_name",
  "snp_name",
  "seq_base",
  "wsmaf",
  "plmaf",
  "reads"
))
