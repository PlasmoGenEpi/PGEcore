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

# 10 monoclonal specimens at five SNPs 100 kb apart plus one SNP ("stut") 100 bp
# from "good1"; specimens s1-s6 also carry a 5% stutter read at "stut", and
# s7-s10 are true mixtures with a minor allele at every SNP.
make_artifact_snp_table <- function() {
  snps <- tibble::tibble(
    snp_name = c("good1", "stut", "good2", "good3", "good4", "good5"),
    pos = c(100000, 100100, 200000, 300000, 400000, 500000)
  )
  specs <- paste0("s", 1:10)
  # "stut" is monomorphic apart from stutter; "good1" is moderately diverse
  major <- function(i, j) {
    if (snps$snp_name[j] == "stut") return("A")
    if (snps$snp_name[j] == "good1") return(if (i <= 8) "A" else "T")
    if ((i + j) %% 2 == 0) "A" else "T"
  }
  rows <- list()
  for (i in seq_along(specs)) {
    for (j in seq_len(nrow(snps))) {
      maj <- major(i, j)
      alt <- if (maj == "A") "T" else "A"
      stutter <- snps$snp_name[j] == "stut" && i <= 6
      mixed <- i >= 7
      rows[[length(rows) + 1]] <- tibble::tibble(
        specimen_name = specs[i], snp_name = snps$snp_name[j], pos = snps$pos[j],
        seq_base = if (stutter || mixed) c(maj, alt) else maj,
        reads = if (stutter) c(95, 5) else if (mixed) c(70, 30) else 100
      )
    }
  }
  dplyr::bind_rows(rows) |>
    dplyr::mutate(target_name = "t1", chrom = "chr1", ref_base = "A", is_biallelic = TRUE)
}

test_that("flag_artifact_snps flags a minor allele recurring in monoclonal specimens", {
  flags <- PGEcore:::flag_artifact_snps(make_artifact_snp_table())
  expect_setequal(flags$snp_name[flags$artifact], "stut")
  # true mixtures carry minor alleles everywhere and are not eligible
  expect_equal(flags$eligible[flags$snp_name == "good1"], 6)
})

test_that("drop_artifact_snps removes the artifact before ranking", {
  snp <- make_artifact_snp_table()
  kept_default <- unique(filter_to_highest_diversity_independent_snp_call(snp)$snp_name)
  kept_drop <- unique(filter_to_highest_diversity_independent_snp_call(snp, drop_artifact_snps = TRUE)$snp_name)
  # stutter inflates heterozygosity, so by default it wins its window over good1
  expect_true("stut" %in% kept_default)
  expect_false("good1" %in% kept_default)
  expect_false("stut" %in% kept_drop)
  expect_true("good1" %in% kept_drop)
})

test_that("exclude_snp_names drops the named SNPs", {
  out <- filter_to_highest_diversity_independent_snp_call(
    make_artifact_snp_table(),
    exclude_snp_names = "stut,good5"
  )
  expect_false(any(c("stut", "good5") %in% out$snp_name))
  expect_true("good1" %in% out$snp_name)
})

test_that("dropped_snps_output lists every dropped SNP, and is empty when none are", {
  out <- tempfile()
  suppressMessages(filter_to_highest_diversity_independent_snp_call(
    make_artifact_snp_table(),
    drop_artifact_snps = TRUE,
    exclude_snp_names = "good5",
    dropped_snps_output = out
  ))
  expect_setequal(readLines(out), c("good5", "stut"))
  expect_setequal(PGEcore:::parse_name_list_arg(out), c("good5", "stut"))

  none <- tempfile()
  filter_to_highest_diversity_independent_snp_call(make_artifact_snp_table(), dropped_snps_output = none)
  expect_equal(PGEcore:::parse_name_list_arg(none), character(0))
})
