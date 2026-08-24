test_that("IDM_wrapper requires exactly one input table and slaf_output", {
  expect_error(
    IDM_wrapper(slaf_output = "out.tsv"),
    "One and only one|Rmpfr|openxlsx"
  )
  expect_error(
    IDM_wrapper(
      allele_table_input = "a.tsv",
      aa_calls_input = "b.tsv",
      slaf_output = "out.tsv"
    ),
    "One and only one|Rmpfr|openxlsx"
  )
  expect_error(
    IDM_wrapper(allele_table_input = "a.tsv", slaf_output = NULL),
    "slaf_output"
  )
  expect_error(
    IDM_wrapper(
      allele_table_input = "a.tsv",
      slaf_output = "out.tsv",
      model = "NOPE"
    ),
    "model|Rmpfr|openxlsx"
  )
})

test_that("prepare_input_4_allele_table validates columns", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("specimen_name\ttarget_name\nS1\tL1", tmp)
  expect_error(
    prepare_input_4_allele_table(tmp),
    "validation|seq|reads"
  )
})

test_that("IDM_wrapper integration skips without Rmpfr", {
  skip_if_not_installed("Rmpfr")
  skip_if_not_installed("openxlsx")
  path <- system.file("extdata", "example_allele_table.tsv", package = "PGEcore")
  skip_if_not(nzchar(path) && file.exists(path))
  dat <- prepare_input_4_allele_table(path)
  expect_true(all(c("specimen_name", "locus", "variants") %in% names(dat)))
  out <- tempfile(fileext = ".tsv")
  on.exit(unlink(out), add = TRUE)
  res <- IDM_wrapper(
    allele_table_input = path,
    slaf_output = out,
    model = "OM"
  )
  expect_true(file.exists(out))
  expect_true(all(c("target_name", "seq", "freq") %in% names(readr::read_tsv(out, show_col_types = FALSE))))
  expect_true("variant" %in% names(res) || "freq" %in% names(res))
})

test_that("MultiLociBiallelicModel_wrapper validates required arguments", {
  expect_error(
    MultiLociBiallelicModel_wrapper(
      aa_calls = NULL,
      loci_group_table = "g.tsv",
      mlaf_output = "o.tsv"
    ),
    "required|variantstring"
  )
})

test_that("create_MultiLociBiallelicModel_input validates loci groups", {
  aa <- tempfile(fileext = ".tsv")
  lg <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(aa, lg)), add = TRUE)
  readr::write_tsv(
    tibble::tibble(
      specimen_name = "S1",
      gene_id = "g1",
      aa_position = 1L,
      ref_aa = "A",
      aa = "A"
    ),
    aa
  )
  readr::write_tsv(
    tibble::tibble(group_id = NA_character_, gene_id = "g1", aa_position = 1L),
    lg
  )
  expect_error(
    create_MultiLociBiallelicModel_input(aa, lg),
    "validation|loci_group"
  )
})

test_that("make_stave groups mutations by gene", {
  variant <- "pfdhfr_1_150:51:N;pfdhfr_1_150:59:C;pfdhps_400_550:437:A"
  expect_equal(
    make_stave(variant),
    "pfdhfr_1_150:51_59:N_C;pfdhps_400_550:437:A"
  )
})

test_that("MultiLociBiallelicModel_wrapper integration skips without variantstring", {
  skip_if_not_installed("variantstring")
  ver <- as.character(utils::packageVersion("variantstring"))
  skip_if(
    utils::compareVersion(ver, "1.0.0") < 0 ||
      utils::compareVersion(ver, "2.0.0") >= 0
  )

  aa <- tempfile(fileext = ".tsv")
  lg <- tempfile(fileext = ".tsv")
  out <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(aa, lg, out)), add = TRUE)
  readr::write_tsv(
    tibble::tibble(
      specimen_name = c("S1", "S1", "S2", "S2"),
      gene_id = "g1",
      aa_position = c(1L, 2L, 1L, 2L),
      ref_aa = "A",
      aa = c("A", "A", "T", "A")
    ),
    aa
  )
  readr::write_tsv(
    tibble::tibble(
      group_id = "grp",
      gene_id = "g1",
      aa_position = c(1L, 2L)
    ),
    lg
  )
  res <- MultiLociBiallelicModel_wrapper(
    aa_calls = aa,
    loci_group_table = lg,
    mlaf_output = out
  )
  expect_true(file.exists(out))
  expect_true(all(c("group_id", "variant", "freq") %in% names(res)))
})
