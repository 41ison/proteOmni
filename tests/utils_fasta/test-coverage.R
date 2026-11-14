test_that("peptide_starts finds every non-overlapping occurrence", {
  expect_equal(peptide_starts("PEP", "XPEPXXPEP"), c(2L, 7L))
  expect_equal(peptide_starts("AA", "AAAA"), c(1L, 3L))
  # Case-insensitive, and the peptide is matched literally, not as a regex.
  expect_equal(peptide_starts("pep", "xPEPx"), 2L)
  expect_equal(peptide_starts("A.C", "ABCA.C"), 4L)
})

test_that("peptide_starts returns integer(0) when there is nothing to match", {
  expect_identical(peptide_starts("ZZ", "XPEPX"), integer(0))
  expect_identical(peptide_starts("", "XPEPX"), integer(0))
  expect_identical(peptide_starts(NA, "XPEPX"), integer(0))
  expect_identical(peptide_starts(c("PEP", "X"), "XPEPX"), integer(0))
})

test_that("protein_coverage accumulates depth over all occurrences", {
  cov <- protein_coverage("ABCDEFABC", c("ABC", "XYZ", "DE"))
  expect_equal(cov$depth, c(1, 1, 1, 1, 1, 0, 1, 1, 1))
  expect_equal(cov$mask, cov$depth > 0)
  expect_equal(cov$n_unmapped, 1L)
})

test_that("protein_coverage weights peptides by counts, NA counting as 1", {
  cov <- protein_coverage(
    "ABCDEFABC",
    c("ABC", "XYZ", "DE", "BCD"),
    c(2, 5, NA, 3)
  )
  expect_equal(cov$depth, c(2, 5, 5, 4, 1, 0, 2, 2, 2))
  expect_equal(sum(cov$mask), 8L)
  expect_equal(cov$n_unmapped, 1L)
})

test_that("protein_coverage with no peptides covers nothing", {
  cov <- protein_coverage("ABCDEF", character(0))
  expect_equal(cov$depth, numeric(6))
  expect_false(any(cov$mask))
  expect_equal(cov$n_unmapped, 0L)
})
