# Simulations run on a 400-protein fit with few replicates to stay fast; they
# check structure and invariants, not the precision of the MDD estimate.
sim_fit <- fit_pwrquant(pq_log2_complete[seq_len(400), ], pq_groups)
small_grid <- c(0.5, 1, 2, 4)

run_sim <- function(...) {
  simulate_pipeline_mdd(
    contrast_fit = sim_fit$contrast_fit,
    design = sim_fit$design,
    contrast_matrix = sim_fit$contrast_matrix,
    lfc_grid = small_grid,
    n_sim = 2L,
    ...
  )
}

test_that("grid search returns a full power curve and per-contrast MDD", {
  res <- run_sim(seed = 2718)
  expect_named(res, c("curve", "mdd", "settings"))
  expect_equal(nrow(res$curve), length(small_grid) * length(pq_contrasts))
  expect_true(all(res$curve$Power >= 0 & res$curve$Power <= 1))
  expect_true(all(res$curve$FDP >= 0 & res$curve$FDP <= 1))
  expect_true(all(res$curve$Missing_Rate == 0))
  expect_setequal(res$mdd$comparison, pq_contrasts)
  ok <- !is.na(res$mdd$MDD_Log2FC)
  expect_true(all(res$mdd$MDD_Log2FC[ok] >= min(small_grid)))
  expect_true(all(res$mdd$MDD_Log2FC[ok] <= max(small_grid)))
  expect_equal(res$settings$n_proteins, 400L)
  expect_equal(res$settings$n_samples, 9L)
})

test_that("power increases with the true effect size", {
  res <- run_sim(seed = 2718)
  by_lfc <- tapply(res$curve$Power, res$curve$True_LFC, mean)
  expect_gt(by_lfc[["4"]], by_lfc[["0.5"]])
})

test_that("a fixed seed makes the simulation reproducible", {
  expect_identical(run_sim(seed = 9051)$curve, run_sim(seed = 9051)$curve)
})

test_that("bisection evaluates fewer points than the grid", {
  res <- run_sim(seed = 2718, search = "bisect", bisect_iter = 3L)
  per_contrast <- table(res$curve$comparison)
  expect_true(all(per_contrast >= 2 & per_contrast <= 5))
  expect_setequal(res$mdd$comparison, pq_contrasts)
})

test_that("progress callback fires once per replicate", {
  n <- 0L
  run_sim(seed = 2718, progress = function() n <<- n + 1L)
  expect_equal(n, 2L * length(small_grid) * length(pq_contrasts))
})

test_that("missingness arm uses the fitted MNAR model and imputation", {
  mm <- estimate_missingness_model(pq_log2)
  res <- run_sim(seed = 2718, miss_model = mm, impute_method = "minprob")
  expect_true(all(res$curve$Missing_Rate > 0))
  expect_true(res$settings$missingness)
  expect_equal(res$settings$target_missing, mm$overall_rate)

  res_ls <- run_sim(seed = 2718, miss_model = mm, impute_method = NULL)
  expect_true(all(res_ls$curve$Missing_Rate > 0))
})

test_that("replicate sweep returns one MDD per contrast and n", {
  sweep <- simulate_mdd_by_replicates(
    sim_fit$contrast_fit,
    sim_fit$design,
    sim_fit$contrast_matrix,
    n_grid = c(4, 3),
    lfc_grid = small_grid,
    n_sim = 2L,
    seed = 2718
  )
  expect_equal(nrow(sweep), 2L * length(pq_contrasts))
  expect_equal(sort(unique(sweep$n_per_group)), c(3, 4))
  expect_equal(sweep$n_samples, sweep$n_per_group * 3)
})
