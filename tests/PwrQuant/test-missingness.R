miss_model <- estimate_missingness_model(pq_log2)

test_that("missingness model is estimated from the training matrix", {
  expect_type(miss_model, "list")
  expect_named(miss_model, c("coefficients", "overall_rate", "abundance_range"))
  expect_equal(miss_model$overall_rate, mean(is.na(pq_log2)))
  # MNAR: low-abundance proteins are missing more often
  expect_lt(miss_model$coefficients[2], 0)
  expect_equal(
    miss_model$abundance_range,
    range(rowMeans(pq_log2, na.rm = TRUE), na.rm = TRUE)
  )
})

test_that("complete data give no missingness model", {
  expect_null(estimate_missingness_model(pq_log2_complete))
})

test_that("simulated missingness respects the per-group minimum", {
  set.seed(6143)
  complete <- pq_log2_complete[1:1000, ]
  holed <- apply_simulated_missingness(
    complete,
    miss_model,
    pq_groups,
    min_valid = 2L
  )
  expect_equal(dim(holed), dim(complete))
  expect_gt(mean(is.na(holed)), 0)
  expect_equal(holed[!is.na(holed)], complete[!is.na(holed)])
  for (g in unique(pq_groups)) {
    expect_true(all(rowSums(!is.na(holed[, pq_groups == g])) >= 2))
  }
})

test_that("drop probability extremes behave as expected", {
  complete <- pq_log2_complete[1:200, ]
  never <- list(coefficients = c(-50, 0))
  always <- list(coefficients = c(50, 0))
  expect_false(anyNA(apply_simulated_missingness(complete, never, pq_groups)))

  set.seed(4472)
  holed <- apply_simulated_missingness(
    complete,
    always,
    pq_groups,
    min_valid = 2L
  )
  for (g in unique(pq_groups)) {
    expect_true(all(rowSums(!is.na(holed[, pq_groups == g])) == 2))
  }
})
