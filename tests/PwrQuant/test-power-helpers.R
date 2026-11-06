test_that("bh_effective_alpha follows the BH boundary", {
  expect_true(is.na(bh_effective_alpha(numeric(0))))
  expect_true(is.na(bh_effective_alpha(c(NA, NaN))))

  # No rejections -> Bonferroni-like floor fdr / m
  expect_equal(bh_effective_alpha(rep(0.5, 20), fdr = 0.05), 0.05 / 20)

  p <- c(0.001, 0.002, 0.2, 0.5, 0.9)
  n_rej <- sum(p.adjust(p, "BH") <= 0.05)
  expect_equal(bh_effective_alpha(p, 0.05), 0.05 * n_rej / length(p))
  # Non-finite p-values are dropped before counting m
  expect_equal(bh_effective_alpha(c(p, NA), 0.05), bh_effective_alpha(p, 0.05))
})

test_that("solve_ncp achieves the requested power", {
  alpha <- 0.01
  df <- 12
  ncp <- solve_ncp(alpha, df, power = 0.8)
  crit <- qt(1 - alpha / 2, df)
  achieved <- pt(crit, df, ncp, lower.tail = FALSE) + pt(-crit, df, ncp)
  expect_equal(achieved, 0.8, tolerance = 1e-5)
})

test_that("solve_ncp approaches the normal approximation at large df", {
  expect_equal(
    solve_ncp(0.05, Inf, 0.8),
    qnorm(0.975) + qnorm(0.8),
    tolerance = 1e-3
  )
})

test_that("solve_ncp is vectorised, monotone in df, and handles edge cases", {
  v <- solve_ncp(c(0.05, 0.05, 0.01), c(4, 40, 4))
  expect_length(v, 3)
  expect_equal(v[1], solve_ncp(0.05, 4))
  expect_gt(v[1], v[2]) # fewer df -> larger ncp needed
  expect_gt(v[3], v[1]) # stricter alpha -> larger ncp needed

  expect_true(is.na(solve_ncp(0, 10)))
  expect_true(is.na(solve_ncp(1, 10)))
  expect_true(is.na(solve_ncp(NA_real_, 10)))
  expect_equal(solve_ncp(0.5, 10, power = 0.4), 0) # power below alpha
})
