skip_if_not_installed("edgeR")

test_that("TMM output keeps shape, dimnames and missingness", {
  out <- tmm_normalize_log2(pq_log2)
  expect_equal(dim(out), dim(pq_log2))
  expect_equal(dimnames(out), dimnames(pq_log2))
  expect_identical(is.na(out), is.na(pq_log2))
})

test_that("TMM applies a single log2 shift per sample", {
  out <- tmm_normalize_log2(pq_log2)
  shift <- pq_log2 - out
  per_col_sd <- apply(shift, 2, sd, na.rm = TRUE)
  expect_true(all(per_col_sd < 1e-10))

  info <- attr(out, "tmm")
  expect_equal(
    unname(apply(shift, 2, function(x) x[!is.na(x)][1])),
    log2(info$scale_factor)
  )
})

test_that("TMM factor table is well formed and centred at geometric mean 1", {
  info <- attr(tmm_normalize_log2(pq_log2), "tmm")
  expect_equal(nrow(info), ncol(pq_log2))
  expect_equal(info$sample, colnames(pq_log2))
  expect_equal(sum(info$is_reference), 1L)
  expect_equal(info$method[info$is_reference], "reference")
  expect_true(all(info$method[!info$is_reference] == "TMM"))
  expect_equal(mean(log2(info$scale_factor)), 0, tolerance = 1e-12)
  expect_true(all(info$n_shared >= 50))
})

test_that("TMM recovers a known per-sample loading offset", {
  shifted <- pq_log2
  shifted[, "m1_r2"] <- shifted[, "m1_r2"] + 1
  f0 <- attr(tmm_normalize_log2(pq_log2), "tmm")
  f1 <- attr(tmm_normalize_log2(shifted), "tmm")
  # Doubling one sample should double its factor relative to the others.
  # Precision weights change slightly with the intensity scale, hence the
  # tolerance.
  rel <- log2(f1$scale_factor / f0$scale_factor)
  expect_equal(
    rel[colnames(pq_log2) == "m1_r2"] - mean(rel[-5]),
    1,
    tolerance = 0.02
  )
})

test_that("TMM falls back to median-M / no scaling with little overlap", {
  set.seed(5521)
  small <- matrix(rnorm(30 * 3, mean = 25), nrow = 30)
  colnames(small) <- c("a", "b", "c")
  info <- attr(tmm_normalize_log2(small), "tmm")
  expect_true(all(grepl("^median-M", info$method[!info$is_reference])))

  tiny <- small[1:8, ]
  expect_warning(
    info2 <- attr(tmm_normalize_log2(tiny), "tmm"),
    "scaling factor set to 1"
  )
  expect_true(all(
    info2$method[!info2$is_reference] == "none (insufficient overlap)"
  ))
})

test_that("TMM input validation", {
  expect_error(
    tmm_normalize_log2(pq_log2[, 1, drop = FALSE]),
    "at least two samples"
  )
  bad <- pq_log2_complete[1:100, 1:3]
  bad[-1, 3] <- NA
  expect_error(tmm_normalize_log2(bad), "fewer than two detected")
})
