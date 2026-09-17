# prep_input_prop() must return allele-1 and allele-2 matrices on the same grid:
# McCOIL_prop reads both as flat vectors sized from a1 alone.
make_calls <- function() {
  # locus B is seen with its second allele before locus A is, so the two
  # pivot_wider subsets order their columns differently
  data.frame(
    specimen_name = c("s1", "s1", "s2", "s2", "s1", "s2"),
    snp_name      = c("A", "B", "A", "B", "A", "B"),
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

test_that("counts stay attached to their own locus", {
  inp <- prep_input_prop(make_calls())
  # s1 locus A: allele "A" 10 reads (idx 1), allele "T" 5 reads (idx 2)
  expect_equal(inp$a1["s1", "A"], 10)
  expect_equal(inp$a2["s1", "A"], 5)
})

test_that("a locus missing an allele contributes zero rather than shifting columns", {
  d <- make_calls()
  d <- d[!(d$specimen_name == "s1" & d$snp_name == "A" & d$seq_base == "T"), ]
  inp <- prep_input_prop(d)
  expect_identical(colnames(inp$a1), colnames(inp$a2))
  expect_equal(inp$a2["s1", "A"], 0)
})

test_that("a specimen-locus pair with no reads is marked missing, not zero", {
  d <- data.frame(
    specimen_name = c("s1", "s1", "s2", "s2"),
    snp_name      = c("A", "A", "B", "B"),
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
  expect_equal(inp$a1["s1", "A"], 10)
  expect_equal(inp$a2["s1", "A"], 5)
})
