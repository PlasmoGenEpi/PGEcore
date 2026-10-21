# prep_input_prop() must return allele-1 and allele-2 matrices on the same grid:
# McCOIL_prop reads both as flat vectors sized from a1 alone.
make_calls <- function() {
  # locus B is seen with its second allele before locus A is, so the two
  # pivot_wider subsets order their columns differently
  data.frame(
    specimen_name = c("s1", "s1", "s2", "s2", "s1", "s2"),
    snp_name      = c("A", "B", "A", "B", "A", "B"),
    ref_base      = c("A", "C", "A", "C", "A", "C"),
    seq_base      = c("A", "C", "T", "G", "T", "C"),
    reads         = c(10L, 20L, 30L, 40L, 5L, 7L)
  )
}

test_that("allele matrices share row and column order", {
  inp <- prep_input_prop(make_calls())
  expect_identical(rownames(inp$a1), rownames(inp$a2))
  expect_identical(colnames(inp$a1), colnames(inp$a2))
  expect_identical(dim(inp$a1), dim(inp$a2))
})

test_that("counts stay attached to their own locus, alternate allele first", {
  inp <- prep_input_prop(make_calls())
  # locus A: ref "A", alt "T", so T counts are allele 1
  expect_equal(inp$a1["s1", "A"], 5)
  expect_equal(inp$a2["s1", "A"], 10)
  expect_equal(inp$a1["s2", "A"], 30)
  expect_equal(inp$a2["s2", "A"], 0)
})

test_that("both models index the non-reference allele", {
  df <- make_calls()
  alleles <- mccoil_alt_alleles(df)
  expect_equal(alleles$alt_base[alleles$snp_name == "A"], "T")
  expect_equal(alleles$alt_base[alleles$snp_name == "B"], "G")
  expect_equal(alleles$ref_base, c("A", "C"))
  # categorical: homozygous non-reference scores 1
  cat_in <- prep_input_categorical(df)
  expect_equal(cat_in["s2", "A"], 1)
  expect_equal(cat_in["s1", "A"], 0.5)
  expect_equal(cat_in["s1", "B"], 0)
})

test_that("filter_biallelic drops multiallelic loci for both models", {
  df <- rbind(make_calls(), data.frame(
    specimen_name = "s2", snp_name = "A", ref_base = "A", seq_base = "C", reads = 5L
  ))
  expect_warning(df <- filter_biallelic(df), "1 multiallelic")
  expect_identical(colnames(prep_input_categorical(df)), "B")
  expect_identical(colnames(prep_input_prop(df)$a1), "B")
})

test_that("filter_biallelic drops loci whose reference base is absent", {
  df <- make_calls()
  df$ref_base[df$snp_name == "B"] <- "T"
  expect_warning(df <- filter_biallelic(df), "reference base was never observed: B")
  expect_identical(colnames(prep_input_categorical(df)), "A")
  expect_identical(colnames(prep_input_prop(df)$a1), "A")
})

test_that("ref_base must be consistent within a SNP", {
  df <- make_calls()
  df$ref_base[1] <- "G"
  expect_error(filter_biallelic(df), "ref_base differs")
})

test_that("SLAF rows name the non-reference allele as snp_name:seq_base", {
  draws <- array(
    c(1, 1, 0.2, 0.3, 0.6, 0.7),
    dim = c(2, 1, 3),
    dimnames = list(iteration = NULL, chain = NULL,
                    variable = c("coi[s1]", "freq[A]", "freq[B]"))
  )
  out <- format_mccoil_output(
    list(draws = draws, coi_names = "s1", freq_names = c("A", "B")),
    mccoil_alt_alleles(make_calls())
  )
  expect_equal(out$slaf$variant, c("A:T", "B:G"))
  expect_equal(out$slaf$ref_base, c("A", "C"))
  expect_equal(out$slaf$seq_base, c("T", "G"))
  expect_equal(out$slaf$freq, c(0.25, 0.65))
})

test_that("a locus missing an allele contributes zero rather than shifting columns", {
  d <- make_calls()
  d <- d[!(d$specimen_name == "s1" & d$snp_name == "A" & d$seq_base == "T"), ]
  inp <- prep_input_prop(d)
  expect_identical(colnames(inp$a1), colnames(inp$a2))
  # s1 now carries only the reference base "A", which is allele 2
  expect_equal(inp$a1["s1", "A"], 0)
  expect_equal(inp$a2["s1", "A"], 10)
})

test_that("a specimen-locus pair with no reads is marked missing, not zero", {
  d <- data.frame(
    specimen_name = c("s1", "s1", "s2", "s2"),
    snp_name      = c("A", "A", "B", "B"),
    ref_base      = c("A", "A", "C", "C"),
    seq_base      = c("A", "T", "C", "G"),
    reads         = c(10L, 5L, 7L, 3L)
  )
  inp <- prep_input_prop(d)
  # s1 has no reads at locus B and s2 none at locus A: both alleles go negative
  # so McCOIL_prop's own missing-data guard skips them
  expect_lt(inp$a1["s1", "B"], 0)
  expect_lt(inp$a2["s1", "B"], 0)
  expect_lt(inp$a1["s2", "A"], 0)
  expect_lt(inp$a2["s2", "A"], 0)
  # covered pairs keep their real counts, including a legitimate zero
  expect_equal(inp$a1["s1", "A"], 5)
  expect_equal(inp$a2["s1", "A"], 10)
})
