#' @keywords internal
#' @useDynLib PGEcore, .registration = TRUE
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data :=
#' @importFrom stats median
#' @importFrom optparse OptionParser make_option parse_args
## usethis namespace: end
NULL

# NSE symbols used in tidyr/dplyr/validate pipelines
utils::globalVariables(c(
  "specimen_name",
  "snp_name",
  "seq_base",
  "wsmaf",
  "plmaf",
  "reads",
  "sample_id",
  "target_name",
  "allele",
  "chrom",
  "#chrom",
  "pos",
  "start",
  "end",
  "strand",
  "ref_base",
  "ref_seq",
  "seq",
  "gene",
  "gene_id",
  "aa_position",
  "aa_locus",
  "ref_aa",
  "aa",
  "target_value",
  "stats",
  "fws",
  "group_id",
  "unique_targets",
  "population",
  "target_count",
  "frequency",
  "freq",
  "sample_total",
  "host_id",
  "coi_estimate",
  "sequence",
  "haplo.set",
  "post.p",
  "Var1",
  "Var2",
  "i",
  "pair",
  "variant",
  "wsaf",
  "is_biallelic",
  "name",
  "length",
  "output$pred.haplo.set"
))
