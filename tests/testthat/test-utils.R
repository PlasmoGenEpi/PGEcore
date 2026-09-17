test_that("check_optparse_required_args detects missing flags", {
  expect_error(
    check_optparse_required_args(list(output = "x"), c("snp_calls", "output")),
    "--snp_calls"
  )
})

test_that("validate_required_columns detects missing columns and empty data", {
  expect_error(
    validate_required_columns(data.frame(a = 1), c("a", "b"), "demo"),
    "Missing required columns"
  )
  expect_error(
    validate_required_columns(data.frame(a = integer()), "a", "demo"),
    "empty"
  )
})

test_that("set_decompose returns expected partitions", {
  out <- set_decompose(c("a", "b"), c("b", "c"))
  expect_equal(out$only_in_vector_a, "a")
  expect_equal(out$only_in_vector_b, "c")
  expect_equal(out$shared, "b")
})

test_that("convert_single_locus_table_to_stave builds variant ids", {
  df <- data.frame(
    gene_id = "PF3D7_0417200.1",
    aa_position = 51,
    aa = "I",
    prev = 0.5
  )
  out <- convert_single_locus_table_to_stave(df, additional_columns = "prev")
  expect_equal(out$variant, "PF3D7_0417200.1:51:I")
  expect_equal(out$prev, 0.5)
})

test_that("stop_if_output_exists respects overwrite", {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines("x", tmp)
  expect_error(stop_if_output_exists(tmp, overwrite = FALSE), "already exists")
  expect_true(stop_if_output_exists(tmp, overwrite = TRUE))
  expect_true(stop_if_output_exists(NULL, overwrite = FALSE))
})

test_that("ensure_output_directory creates and optionally replaces a directory", {
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  ensure_output_directory(tmp, overwrite_dir = FALSE)
  expect_true(dir.exists(tmp))
  expect_error(ensure_output_directory(tmp, overwrite_dir = FALSE), "already exist")
  ensure_output_directory(tmp, overwrite_dir = TRUE)
  expect_true(dir.exists(tmp))
})

test_that("stop_on_validate_fails reports failing rules", {
  skip_if_not_installed("validate")
  df <- data.frame(a = c(1, NA_real_))
  rules <- validate::validator(!is.na(a))
  expect_error(stop_on_validate_fails(df, rules, "demo"), "demo")
})