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

  # Sample-wise: draws come from N(mean - 1.8 SD, 0.3 SD) of each sample's own
  # observed values, so they sit well below that sample's mean.
  for (j in seq_len(ncol(imp_input))) {
    col_obs <- imp_input[observed[, j], j]
    draws <- out[!observed[, j], j]
    if (!length(draws)) {
      next
    }
    expect_true(all(draws < mean(col_obs)))
    expect_lt(
      abs(mean(draws) - (mean(col_obs) - 1.8 * sd(col_obs))),
      0.5 * sd(col_obs)
    )
  }
})

test_that("MinProb draws do not depend on the protein's own values", {
  # Shifting one protein's observed values must not move its imputed draws.
  row <- which(rowSums(!observed) > 0)[1]
  shifted <- imp_input
  shifted[row, ] <- shifted[row, ] + 5
  set.seed(5108)
  a <- groupwise_imputation(imp_input, pq_groups, "minprob", verbose = FALSE)
  set.seed(5108)
  b <- groupwise_imputation(shifted, pq_groups, "minprob", verbose = FALSE)
  miss <- !observed[row, ]
  expect_equal(a[row, miss], b[row, miss], tolerance = 0.05)
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

test_that("KNN and missForest also leave all-NA-within-group rows as NA", {
  skip_if_not_installed("impute")
  skip_if_not_installed("missForest")
  x <- rbind(pq_log2["A0A0U1RRL7", , drop = FALSE], imp_input[1:60, ])
  set.seed(6620)
  for (m in c("knn", "missforest")) {
    out <- suppressWarnings(
      groupwise_imputation(x, pq_groups, method = m, verbose = FALSE)
    )
    expect_true(all(is.na(out["A0A0U1RRL7", pq_groups == "m0"])))
    expect_false(anyNA(out[-1, ]))
  }
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

test_that("KNN k is chosen by cross-validation, reproducibly, without touching the RNG", {
  skip_if_not_installed("impute")
  set.seed(2207)
  before <- .Random.seed
  k1 <- select_knn_k(imp_input, pq_groups)
  expect_identical(.Random.seed, before)
  expect_identical(select_knn_k(imp_input, pq_groups), k1)
  expect_true(k1 %in% c(3L, 5L, 8L, 10L, 15L, 20L))

  out <- groupwise_imputation(imp_input, pq_groups, "knn", verbose = FALSE)
  expect_identical(attr(out, "knn_k"), k1)
  fixed <- groupwise_imputation(
    imp_input,
    pq_groups,
    "knn",
    verbose = FALSE,
    knn_k = 7
  )
  expect_identical(attr(fixed, "knn_k"), 7L)
})

test_that("KNN k is bounded by the number of proteins, not samples", {
  skip_if_not_installed("impute")
  expect_identical(select_knn_k(imp_input[1:4, ], pq_groups), 3L)
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
