test_that("calculate_prevalence counts specimens per allele", {
  allele_table <- tibble::tibble(
    specimen_name = c("s1", "s1", "s2"),
    target_name = c("L1", "L1", "L1"),
    variant = c("A", "B", "A")
  )
  prev <- calculate_prevalence(allele_table)
  expect_equal(prev$prev[prev$variant == "A"], 1)
  expect_equal(prev$prev[prev$variant == "B"], 0.5)
  expect_equal(prev$sample_total[1], 2)
})

test_that("genotype prevalence counts each specimen once as its allele set", {
  allele_table <- tibble::tibble(
    specimen_name = c("s1", "s2", "s2", "s3"),
    target_name = "L1",
    variant = c("A", "A", "B", "B")
  )
  allele_prev <- calculate_prevalence(allele_table)
  expect_equal(allele_prev$prev[allele_prev$variant == "A"], 2 / 3)
  expect_equal(allele_prev$prev[allele_prev$variant == "B"], 2 / 3)

  geno <- calculate_genotype_prevalence(allele_table)
  expect_equal(sort(geno$variant), c("A", "A/B", "B"))
  expect_equal(geno$prev[geno$variant == "A"], 1 / 3)
  expect_equal(geno$prev[geno$variant == "A/B"], 1 / 3)
  expect_equal(geno$prev[geno$variant == "B"], 1 / 3)
  expect_equal(sum(geno$prev), 1)
})

test_that("observed genotypes sort alleles and ignore duplicate rows", {
  allele_table <- tibble::tibble(
    specimen_name = c("s1", "s1", "s1"),
    target_name = "L1",
    variant = c("B", "A", "A")
  )
  collapsed <- collapse_to_observed_genotype(allele_table)
  expect_equal(collapsed$variant, "A/B")
})

test_that("estimate_allele_prevalence_naive validates exclusive inputs", {
  expect_error(
    estimate_allele_prevalence_naive(),
    "One and only one"
  )
  expect_error(
    estimate_allele_prevalence_naive(
      aa_calls = "missing.tsv",
      allele_table = "also_missing.tsv"
    ),
    "One and only one"
  )

  aa_path <- system.file(
    "extdata",
    "example_aa_calls.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(aa_path) && file.exists(aa_path))
  expect_error(
    estimate_allele_prevalence_naive(aa_calls = aa_path, method = "nope"),
    "not a valid method"
  )
})

test_that("estimate_allele_prevalence_naive works on packaged example data", {
  aa_path <- system.file(
    "extdata",
    "example_aa_calls.tsv",
    package = "PGEcore"
  )
  mh_path <- system.file(
    "extdata",
    "example_allele_table.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(aa_path) && file.exists(aa_path))
  skip_if_not(nzchar(mh_path) && file.exists(mh_path))

  aa_out <- estimate_allele_prevalence_naive(aa_calls = aa_path)
  expect_true(
    all(c("variant", "prev", "sample_count", "sample_total") %in% names(aa_out))
  )
  expect_false("gene_id" %in% names(aa_out))

  aa_geno <- estimate_allele_prevalence_naive(
    aa_calls = aa_path,
    method = "genotype"
  )
  expect_true(
    all(c("variant", "prev", "sample_count", "sample_total") %in% names(aa_geno))
  )
  expect_true(any(grepl("/", aa_geno$variant)))

  mh_out <- estimate_allele_prevalence_naive(allele_table = mh_path)
  expect_true(
    all(
      c("target_name", "seq", "prev", "sample_count", "sample_total") %in%
        names(mh_out)
    )
  )

  tmp_out <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_out), add = TRUE)
  written <- estimate_allele_prevalence_naive(
    aa_calls = aa_path,
    output = tmp_out
  )
  expect_true(file.exists(tmp_out))
  expect_equal(nrow(written), nrow(aa_out))
})
