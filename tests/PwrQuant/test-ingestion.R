test_that("fixture is read with the expected shape and sample layout", {
  expect_equal(dim(pq_raw), c(6746L, 9L))
  expect_equal(
    colnames(pq_raw),
    paste0(rep(c("m0", "m1", "m2"), each = 3), "_r", 1:3)
  )
  expect_equal(as.vector(table(pq_groups)), c(3L, 3L, 3L))
  expect_true(is.numeric(pq_raw))
})

test_that("protein IDs are unquoted, unique UniProt accessions", {
  ids <- rownames(pq_raw)
  expect_false(any(grepl('"', ids, fixed = TRUE)))
  expect_false(anyDuplicated(ids) > 0)
  expect_true(all(is_uniprot_accession(ids)))
  expect_true(all(c("P00761", "P02769") %in% ids))
})

test_that("missingness profile matches the training file", {
  expect_equal(sum(rowSums(!is.na(pq_raw)) == 0), 976L)
  expect_equal(sum(stats::complete.cases(pq_raw)), 3608L)
  # MaxLFQ reports non-detects as NA, never 0
  expect_false(any(pq_raw == 0, na.rm = TRUE))
  expect_true(all(pq_raw > 0, na.rm = TRUE))
})

test_that("log2(x + 1) transform stays on a plausible log2 scale", {
  rng <- range(pq_log2, na.rm = TRUE)
  expect_gt(rng[1], 15)
  expect_lt(rng[2], 40)
  expect_identical(is.na(pq_log2), is.na(pq_raw))
})

test_that("groupwise valid-value filter keeps the expected rows", {
  expect_equal(nrow(pq_log2_complete), 3608L)
  expect_true(all(stats::complete.cases(pq_log2_complete)))
  for (g in unique(pq_groups)) {
    n_valid <- rowSums(!is.na(pq_log2_67[, pq_groups == g]))
    expect_true(all(n_valid >= 2))
  }
  expect_gt(nrow(pq_log2_67), nrow(pq_log2_complete))
})
