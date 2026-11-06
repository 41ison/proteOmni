ps <- extract_power_stats(pq_fit$contrast_fit, fdr = 0.05, target_power = 0.80)
n_prot <- nrow(pq_fit$contrast_fit$coefficients)

test_that("pipeline fit on the training matrix is well formed", {
  expect_equal(colnames(pq_fit$design), c("m0", "m1", "m2"))
  expect_equal(colnames(pq_fit$contrast_matrix), pq_contrasts)
  expect_equal(n_prot, nrow(pq_log2_67))
  expect_true(all(is.finite(pq_fit$contrast_fit$s2.post)))
})

test_that("extract_power_stats returns one row per protein and contrast", {
  expect_equal(nrow(ps), n_prot * length(pq_contrasts))
  expect_setequal(unique(ps$comparison), pq_contrasts)
  expect_true(all(
    c(
      "Protein",
      "logFC",
      "P.Value",
      "adj.P.Val",
      "Sigma",
      "SE",
      "df_total",
      "alpha_eff",
      "Conditional_MDD_Log2FC",
      "Above_MDD",
      "comparison"
    ) %in%
      names(ps)
  ))
  expect_equal(
    ps$Protein[ps$comparison == "m1-m0"],
    rownames(pq_fit$contrast_fit$coefficients)
  )
})

test_that("SE reproduces limma's moderated t statistic", {
  expect_equal(ps$t, unname(ps$logFC / ps$SE), tolerance = 1e-8)
})

test_that("alpha_eff is the BH boundary of each contrast", {
  for (cf in pq_contrasts) {
    d <- ps[ps$comparison == cf, ]
    expect_equal(unique(d$alpha_eff), bh_effective_alpha(d$P.Value, 0.05))
  }
})

test_that("conditional MDD attains the target power and flags consistently", {
  expect_true(all(ps$Conditional_MDD_Log2FC > 0))
  d <- ps[1, ]
  ncp <- d$Conditional_MDD_Log2FC / d$SE
  crit <- qt(1 - d$alpha_eff / 2, d$df_total)
  power <- pt(crit, d$df_total, ncp, lower.tail = FALSE) +
    pt(-crit, d$df_total, ncp)
  expect_equal(power, 0.80, tolerance = 1e-4)
  expect_identical(ps$Above_MDD, abs(ps$logFC) >= ps$Conditional_MDD_Log2FC)
})

test_that("the training data yield differential proteins", {
  n_sig <- tapply(ps$adj.P.Val <= 0.05, ps$comparison, sum)
  expect_true(all(n_sig > 0))
})

test_that("compiled extract_power_stats_fast matches", {
  expect_identical(extract_power_stats_fast(pq_fit$contrast_fit), ps)
})

test_that("TMM-normalised matrix runs through the pipeline", {
  skip_if_not_installed("edgeR")
  norm <- tmm_normalize_log2(pq_log2_67)
  fit <- fit_pwrquant(norm, pq_groups)
  out <- extract_power_stats(fit$contrast_fit)
  expect_equal(nrow(out), nrow(norm) * length(pq_contrasts))
  expect_false(anyNA(out$logFC))
})
