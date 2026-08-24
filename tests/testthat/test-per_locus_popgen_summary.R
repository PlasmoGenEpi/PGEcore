test_that("msa_method_binary maps methods to executables", {
  expect_equal(msa_method_binary("Muscle"), "muscle")
  expect_equal(msa_method_binary("ClustalW"), "clustalw")
  expect_equal(msa_method_binary("ClustalOmega"), "clustalo")
  expect_error(msa_method_binary("nope"), "msa_method")
})

test_that("calculate_popgen_stats returns zeros for identical sequences", {
  skip_if_not_installed("ape")
  skip_if_not_installed("msa")
  skip_if_not_installed("pegas")
  stats <- calculate_popgen_stats(c("ATGC", "ATGC"), msa_method = "Muscle")
  expect_equal(stats$Nucleotide_Diversity, 0)
  expect_equal(stats$Segregating_Sites, 0)
  expect_equal(stats$Tajima_D, 0)
  expect_null(stats$Tajima_D_pval_normal)
})

test_that("per_locus_popgen_summary summarises identical alleles without MSA", {
  skip_if_not_installed("ape")
  skip_if_not_installed("msa")
  skip_if_not_installed("pegas")

  tmp_in <- tempfile(fileext = ".tsv")
  tmp_out <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)
  readr::write_tsv(
    tibble::tibble(
      specimen_name = c("s1", "s2"),
      target_name = c("L1", "L1"),
      seq = c("ATGC", "ATGC"),
      reads = c(1, 2)
    ),
    tmp_in
  )
  res <- per_locus_popgen_summary(tmp_in, out = tmp_out)
  expect_true(file.exists(tmp_out))
  expect_equal(res$target_name, "L1")
  expect_equal(res$nucleotide_diversity, 0)
  expect_equal(res$segregating_sites, 0)
  expect_equal(res$tajima_d, 0)
})

test_that("per_locus_popgen_summary rejects unknown msa_method", {
  expect_error(
    per_locus_popgen_summary(tempfile(), msa_method = "nope"),
    "msa_method"
  )
})
