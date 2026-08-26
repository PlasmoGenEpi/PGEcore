test_that("filter_biallelic_calls keeps loci with at most two alleles", {
  aa <- tibble::tibble(
    gene_id = c("g1", "g1", "g1", "g2", "g2"),
    aa_position = c(1, 1, 1, 2, 2),
    ref_aa = c("A", "A", "A", "V", "V"),
    aa = c("A", "T", "G", "V", "I")
  )
  out <- filter_biallelic_calls(aa)
  expect_equal(unique(out$biallelic$gene_id), "g2")
  expect_equal(unique(out$nonbiallelic$gene_id), "g1")
  expect_true(all(out$nonbiallelic$allele_calls == 3))
  expect_true(all(out$biallelic$allele_calls == 2))
})

test_that("filter_biallelic_calls validates columns and overwrite", {
  expect_error(
    filter_biallelic_calls(data.frame(gene_id = "g")),
    "Missing required columns"
  )

  aa <- tibble::tibble(
    gene_id = "g1",
    aa_position = 1,
    ref_aa = "A",
    aa = "A"
  )
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  readr::write_tsv(aa, tmp)
  expect_error(
    filter_biallelic_calls(aa, output = tmp, overwrite = FALSE),
    "already exists"
  )
})

test_that("filter_biallelic_calls works on packaged example", {
  path <- system.file(
    "extdata",
    "example_aa_calls.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(path) && file.exists(path))
  out <- filter_biallelic_calls(path)
  expect_true(all(c("biallelic", "nonbiallelic") %in% names(out)))
  expect_true(nrow(out$biallelic) + nrow(out$nonbiallelic) >= 1)
})
