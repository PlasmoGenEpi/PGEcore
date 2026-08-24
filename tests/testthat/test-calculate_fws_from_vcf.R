test_that("derive_gds_path replaces vcf suffixes", {
  expect_equal(derive_gds_path("x.vcf"), "x.gds")
  expect_equal(derive_gds_path("x.vcf.gz"), "x.gds")
  expect_equal(derive_gds_path("x.VCF.GZ"), "x.gds")
  expect_equal(derive_gds_path("x.bcf"), "x.bcf.gds")
  expect_equal(derive_gds_path("x.vcf", "custom.gds"), "custom.gds")
})

test_that("gds_needs_conversion respects force and mtimes", {
  vcf <- tempfile(fileext = ".vcf")
  gds <- tempfile(fileext = ".gds")
  on.exit(unlink(c(vcf, gds)), add = TRUE)
  writeLines("##fileformat=VCFv4.2", vcf)
  expect_true(gds_needs_conversion(vcf, gds, force = FALSE))
  expect_true(gds_needs_conversion(vcf, gds, force = TRUE))
  writeLines("gds", gds)
  Sys.setFileTime(gds, file.mtime(vcf) + 10)
  expect_false(gds_needs_conversion(vcf, gds, force = FALSE))
  expect_true(gds_needs_conversion(vcf, gds, force = TRUE))
})

test_that("calculate_fws_from_vcf requires SeqArray and moimix", {
  err <- tryCatch(
    calculate_fws_from_vcf(input = tempfile()),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("SeqArray|moimix|not found", err))
})

test_that("calculate_fws_from_vcf errors when the VCF is missing", {
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("moimix")
  missing <- file.path(tempdir(), "no-such-file.vcf")
  expect_error(
    calculate_fws_from_vcf(input = missing),
    "not found"
  )
})
