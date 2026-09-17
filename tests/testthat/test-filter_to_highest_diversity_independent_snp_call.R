test_that("filter_to_highest_diversity_independent_snp_call keeps distant SNPs", {
  # Allele mixes chosen so he(100) > he(20000) > he(150); 150 is within
  # mindist of 100 and should be dropped.
  snp <- tibble::tibble(
    specimen_name = c(
      rep("s1", 5), rep("s2", 5), rep("s3", 5), rep("s4", 5)
    ),
    target_name = "t1",
    chrom = "chr1",
    pos = c(
      100, 100, 150, 20000, 20000,
      100, 100, 150, 20000, 20000,
      100, 100, 150, 20000, 20000,
      100, 100, 150, 20000, 20000
    ),
    snp_name = paste0("snp-", pos),
    ref_base = "A",
    seq_base = c(
      "A", "T", "A", "A", "T",
      "A", "T", "A", "A", "G",
      "A", "T", "A", "A", "T",
      "A", "T", "A", "A", "T"
    ),
    reads = 10,
    is_biallelic = TRUE
  )

  out <- filter_to_highest_diversity_independent_snp_call(
    snp,
    mindist_between_snps = 10000
  )
  kept_pos <- sort(unique(out$pos))
  expect_true(100 %in% kept_pos)
  expect_true(20000 %in% kept_pos)
  expect_false(150 %in% kept_pos)
})

test_that("filter_to_highest_diversity_independent_snp_call only_informative", {
  snp <- tibble::tibble(
    specimen_name = c("s1", "s2", "s1", "s2"),
    target_name = "t1",
    chrom = "chr1",
    pos = c(100L, 100L, 50000L, 50000L),
    snp_name = c("a", "a", "b", "b"),
    ref_base = "A",
    seq_base = c("A", "T", "A", "A"),
    reads = 10,
    is_biallelic = TRUE
  )
  out <- filter_to_highest_diversity_independent_snp_call(
    snp,
    mindist_between_snps = 1000,
    only_informative = TRUE
  )
  expect_equal(unique(out$pos), 100L)
})

test_that("filter_to_highest_diversity_independent_snp_call validates input", {
  expect_error(
    filter_to_highest_diversity_independent_snp_call(data.frame(x = 1)),
    "Missing required columns"
  )
})
