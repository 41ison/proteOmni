test_that("compute_cv_mtx returns one CV column per group, in % units", {
  cv <- compute_cv_mtx(pq_raw, pq_groups)
  expect_equal(colnames(cv), c("CV_m0", "CV_m1", "CV_m2"))
  expect_equal(rownames(cv), rownames(pq_raw))

  x <- pq_raw["P02769", pq_groups == "m2"]
  expect_equal(cv["P02769", "CV_m2"], 100 * sd(x) / mean(x))
  expect_true(all(cv >= 0, na.rm = TRUE))
})

test_that("compute_cv_mtx gives NA/NaN for rows with < 2 values in a group", {
  cv <- compute_cv_mtx(pq_raw, pq_groups)
  # P00366 is entirely missing
  expect_true(all(is.na(unlist(cv["P00366", ]))))
  n_m0 <- rowSums(!is.na(pq_raw[, pq_groups == "m0"]))
  expect_true(all(is.na(cv$CV_m0[n_m0 < 2])))
  expect_false(anyNA(cv$CV_m0[n_m0 >= 2]))
})

test_that("compiled compute_cv_mtx_fast matches the interpreted version", {
  sub <- pq_raw[1:300, ]
  expect_identical(
    compute_cv_mtx_fast(sub, pq_groups),
    compute_cv_mtx(sub, pq_groups)
  )
})
