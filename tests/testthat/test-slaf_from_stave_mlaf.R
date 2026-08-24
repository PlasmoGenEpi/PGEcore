test_that("slaf_from_stave_mlaf checks suggested package and columns", {
  err <- tryCatch(
    slaf_from_stave_mlaf(data.frame(x = 1)),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("variantstring|Missing required columns", err))
})

test_that("slaf_from_stave_mlaf expands packaged MLAF when variantstring is present", {
  skip_if_not_installed("variantstring")

  path <- system.file("extdata", "example_mlaf.tsv", package = "PGEcore")
  skip_if_not(nzchar(path) && file.exists(path))

  out <- slaf_from_stave_mlaf(path)
  expect_true(all(c("variant", "freq") %in% names(out)))
  expect_true(nrow(out) >= 1)
  expect_true(all(grepl(":", out$variant)))
})
