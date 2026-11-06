# A 400-protein slice keeps the KNN / BPCA / missForest runs fast.
imp_input <- pq_log2_67[seq_len(400), ]
observed <- !is.na(imp_input)

test_that("complete data is returned untouched", {
  expect_identical(
    groupwise_imputation(pq_log2_complete[1:50, ], pq_groups, verbose = FALSE),
    pq_log2_complete[1:50, ]
  )
})

test_that("MinProb fills every gap, keeps observed values, and draws low", {
  set.seed(7310)
  out <- groupwise_imputation(
    imp_input,
    pq_groups,
    method = "minprob",
    verbose = FALSE
  )
  expect_false(anyNA(out))
  expect_identical(out[observed], imp_input[observed])

  # Draws come from mean - 1.8 SD of the protein's own group values, so they
  # should sit below the observed group mean almost always.
  below <- c()
  for (g in unique(pq_groups)) {
    cols <- pq_groups == g
    grp_mean <- rowMeans(imp_input[, cols], na.rm = TRUE)
    miss <- which(!observed[, cols], arr.ind = TRUE)
    below <- c(below, out[, cols][miss] < grp_mean[miss[, "row"]])
  }
  expect_gt(mean(below), 0.95)
})

test_that("MinProb is reproducible under a seed", {
  set.seed(1904)
  a <- groupwise_imputation(
    imp_input,
    pq_groups,
    method = "minprob",
    verbose = FALSE
  )
  set.seed(1904)
  b <- groupwise_imputation(
    imp_input,
    pq_groups,
    method = "minprob",
    verbose = FALSE
  )
  expect_identical(a, b)
})

test_that("MinProb leaves rows that are all-NA within a group as NA", {
  x <- pq_log2[c("A0A0U1RRL7", "P02769"), ] # first is missing in all of m0
  set.seed(2266)
  out <- suppressWarnings(
    groupwise_imputation(x, pq_groups, method = "minprob", verbose = FALSE)
  )
  expect_true(all(is.na(out["A0A0U1RRL7", pq_groups == "m0"])))
  expect_false(anyNA(out["A0A0U1RRL7", pq_groups == "m1"]))
})

test_that("KNN imputation fills gaps and keeps observed values", {
  skip_if_not_installed("impute")
  out <- groupwise_imputation(
    imp_input,
    pq_groups,
    method = "knn",
    verbose = FALSE
  )
  expect_false(anyNA(out))
  expect_equal(out[observed], imp_input[observed])
  expect_equal(dimnames(out), dimnames(imp_input))
})

test_that("BPCA imputation fills gaps and keeps observed values", {
  skip_if_not_installed("pcaMethods")
  set.seed(3390)
  out <- groupwise_imputation(
    imp_input,
    pq_groups,
    method = "bpca",
    verbose = FALSE
  )
  expect_false(anyNA(out))
  expect_equal(out[observed], imp_input[observed], tolerance = 1e-8)
  expect_equal(dimnames(out), dimnames(imp_input))
})

test_that("missForest imputation fills gaps and keeps observed values", {
  skip_if_not_installed("missForest")
  set.seed(8127)
  sub <- imp_input[1:60, ]
  out <- groupwise_imputation(
    sub,
    pq_groups,
    method = "missforest",
    verbose = FALSE
  )
  expect_false(anyNA(out))
  expect_equal(out[!is.na(sub)], sub[!is.na(sub)])
})

test_that("unknown method is rejected", {
  expect_error(
    groupwise_imputation(
      imp_input,
      pq_groups,
      method = "zero",
      verbose = FALSE
    ),
    "Unknown imputation method"
  )
})
