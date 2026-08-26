test_that("summarize_allele_table computes counts and singlets", {
  locus_data <- data.frame(
    sample_id = c("s1", "s1", "s2", "s3"),
    target_name = c("L1", "L1", "L1", "L2"),
    allele = c("a", "a", "b", "c"),
    stringsAsFactors = FALSE
  )
  out <- summarize_allele_table(locus_data)
  l1 <- out[out$target_name == "L1", ]
  expect_equal(l1$total_allele_count, 3)
  expect_equal(l1$unique_allele_count, 2)
  expect_equal(l1$allele_singlets, 1L)

  l2 <- out[out$target_name == "L2", ]
  expect_equal(l2$total_allele_count, 1)
  expect_equal(l2$unique_allele_count, 1)
  expect_equal(l2$allele_singlets, 1L)
})

test_that("create_locus_data validates required character columns", {
  tmp_ok <- tempfile(fileext = ".tsv")
  tmp_bad <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_ok, tmp_bad)), add = TRUE)

  write.table(
    data.frame(
      specimen_name = "s1",
      target_name = "L1",
      seq = "ACGT",
      stringsAsFactors = FALSE
    ),
    tmp_ok,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  captured <- capture.output({
    locus <- create_locus_data(tmp_ok)
  })
  expect_equal(names(locus), c("sample_id", "target_name", "allele"))
  expect_true(any(grepl("Reading input data", captured)))

  write.table(
    data.frame(
      specimen_name = NA_character_,
      target_name = "L1",
      seq = "ACGT",
      stringsAsFactors = FALSE
    ),
    tmp_bad,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
  expect_error(
    capture.output(create_locus_data(tmp_bad)),
    "validation checks"
  )
})

test_that("allele_per_locus_summary works on packaged example data", {
  path <- system.file(
    "extdata",
    "example_allele_table.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(path) && file.exists(path))

  tmp_out <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_out), add = TRUE)

  out <- capture.output({
    result <- allele_per_locus_summary(path, output = tmp_out)
  })
  expect_true(file.exists(tmp_out))
  expect_true(
    all(
      c(
        "target_name",
        "total_allele_count",
        "unique_allele_count",
        "allele_singlets"
      ) %in% names(result)
    )
  )
  expect_true(nrow(result) >= 1)
})
