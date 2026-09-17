test_that("slaf_from_mhaps_freqs aggregates and renormalises frequencies", {
  mhaps <- tibble::tibble(
    target_name = c("t1", "t1", "t2", "t2"),
    seq = c("AAA", "TTT", "AAA", "GGG"),
    freq = c(0.6, 0.4, 0.5, 0.5),
    sample_total = c(10, 10, 20, 20)
  )
  loci <- tibble::tibble(
    target_name = c("t1", "t1", "t2", "t2"),
    gene_id = "g1",
    aa_position = 10L,
    seq = c("AAA", "TTT", "AAA", "GGG"),
    aa = c("K", "F", "K", "G")
  )

  out <- slaf_from_mhaps_freqs(mhaps, loci)
  expect_true(all(c("slaf", "per_target_slaf") %in% names(out)))
  # Renormalisation is within (gene_id, aa_position, sample_total) groups.
  grouped_sum <- out$slaf |>
    dplyr::group_by(.data$sample_total) |>
    dplyr::summarise(s = sum(.data$freq), .groups = "drop")
  expect_equal(grouped_sum$s, rep(1, nrow(grouped_sum)), tolerance = 1e-8)
  expect_true(all(grepl("^g1:10:", out$slaf$variant)))
  expect_equal(max(out$slaf$sample_total), 20)
})

test_that("slaf_from_mhaps_freqs validates required columns", {
  expect_error(
    slaf_from_mhaps_freqs(data.frame(a = 1), data.frame(b = 1)),
    "Missing required columns"
  )
})

test_that("slaf_from_mhaps_freqs writes optional outputs", {
  mhaps <- tibble::tibble(
    target_name = "t1",
    seq = "AAA",
    freq = 1,
    sample_total = 5
  )
  loci <- tibble::tibble(
    target_name = "t1",
    gene_id = "g1",
    aa_position = 1L,
    seq = "AAA",
    aa = "K"
  )
  tmp_slaf <- tempfile(fileext = ".tsv")
  tmp_per <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_slaf, tmp_per)), add = TRUE)
  out <- slaf_from_mhaps_freqs(
    mhaps,
    loci,
    slaf_output = tmp_slaf,
    per_target_slaf_output = tmp_per
  )
  expect_true(file.exists(tmp_slaf))
  expect_true(file.exists(tmp_per))
  expect_equal(nrow(out$slaf), 1)
})
