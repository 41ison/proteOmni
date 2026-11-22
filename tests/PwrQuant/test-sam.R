# SAM (Tusher et al. 2001): statistic, permutation null, delta / FDR, and the
# layout of the results that the module's tabs consume.

sam_sim <- function(p = 1200, n = 5, n_diff = 120, shift = 2, seed = 7243) {
  set.seed(seed)
  x <- matrix(
    rnorm(p * 2 * n, mean = 20, sd = rep(runif(p, 0.2, 1.2), 2 * n)),
    nrow = p,
    dimnames = list(paste0("P", seq_len(p)), paste0("s", seq_len(2 * n)))
  )
  x[seq_len(n_diff), seq_len(n)] <- x[seq_len(n_diff), seq_len(n)] + shift
  list(
    x = x,
    groups = rep(c("A", "B"), each = n),
    truth = seq_len(p) <= n_diff
  )
}

sim <- sam_sim()
w_ab <- c(A = 1, B = -1)
set.seed(915)
sam_fit <- sam_contrast(sim$x, sim$groups, w_ab, n_perm = 300)

test_that("r(i) and s(i) are the pooled-variance t-test's numerator and SE", {
  tt <- apply(sim$x, 1, function(v) {
    t <- stats::t.test(v[1:5], v[6:10], var.equal = TRUE)
    c(unname(t$estimate[1] - t$estimate[2]), t$stderr)
  })
  expect_equal(sam_fit$table$logFC, unname(tt[1, ]))
  expect_equal(sam_fit$table$SE, unname(tt[2, ]))
  expect_equal(
    sam_fit$table$d,
    sam_fit$table$logFC / (sam_fit$table$SE + sam_fit$table$s0)
  )
})

test_that("s0 = 0 reproduces the ordinary t statistic", {
  set.seed(915)
  fit0 <- sam_contrast(sim$x, sim$groups, w_ab, n_perm = 50, s0 = 0)
  tstat <- apply(sim$x, 1, function(v) {
    unname(stats::t.test(v[1:5], v[6:10], var.equal = TRUE)$statistic)
  })
  expect_equal(fit0$table$d, unname(tstat))
})

test_that("s0 is estimated as a quantile of s, with the small-data fallback", {
  s <- sam_fit$table$SE
  expect_true(sam_fit$table$s0[1] >= 0 && sam_fit$table$s0[1] <= max(s))
  small <- sam_estimate_s0(sam_fit$table$logFC[1:100], s[1:100])
  expect_equal(small, unname(stats::quantile(s[1:100], 0.05)))
})

test_that("two-group contrasts with few splits are enumerated completely", {
  labs <- sam_labellings(rep(c("A", "B"), each = 3), n_perm = 1000)
  expect_true(attr(labs, "complete"))
  expect_equal(ncol(labs), choose(6, 3))
  expect_equal(nrow(unique(t(labs))), choose(6, 3))
  expect_true(all(colSums(labs == "A") == 3))

  set.seed(1)
  labs <- sam_labellings(rep(c("A", "B"), each = 6), n_perm = 100)
  expect_false(attr(labs, "complete"))
  expect_equal(ncol(labs), 100)
  expect_true(all(colSums(labs == "A") == 6))
})

test_that("blocked labellings keep each group's batch composition", {
  lab <- rep(c("A", "B"), each = 4)
  batch <- rep(c("b1", "b2"), times = 4) # 2 + 2 of each batch per group
  labs <- sam_labellings(lab, n_perm = 1000, block = batch)
  # choose(4, 2)^2 = 36 balanced labellings, as in Tusher et al. (2001)
  expect_true(attr(labs, "complete"))
  expect_equal(attr(labs, "n_distinct"), 36)
  expect_equal(ncol(labs), 36)
  expect_equal(nrow(unique(t(labs))), 36)
  expect_true(all(colSums(labs[batch == "b1", ] == "A") == 2))
  expect_true(all(colSums(labs[batch == "b2", ] == "A") == 2))
  expect_true(any(apply(labs, 2, identical, lab)))

  # Unbalanced blocks, random draws, three groups
  lab3 <- rep(c("A", "B", "C"), times = c(6, 6, 6))
  batch3 <- rep(c("b1", "b1", "b2"), times = 6)
  set.seed(31)
  labs3 <- sam_labellings(lab3, n_perm = 200, block = batch3)
  expect_false(attr(labs3, "complete"))
  for (b in unique(batch3)) {
    for (g in unique(lab3)) {
      expect_true(all(
        colSums(labs3[batch3 == b, , drop = FALSE] == g) ==
          sum(lab3[batch3 == b] == g)
      ))
    }
  }
})

test_that("sam_contrast permutes within blocks, or falls back when it cannot", {
  x <- sim$x[1:600, ]
  batch <- rep(c("b1", "b2", "b1", "b2", "b1"), 2)
  set.seed(4)
  fit <- sam_contrast(x, sim$groups, w_ab, n_perm = 100, block = batch)
  expect_true(fit$summary$blocked)
  expect_equal(fit$summary$n_distinct_labellings, choose(6, 3) * choose(4, 2))
  # Blocking changes the null, never the observed statistic.
  set.seed(4)
  free <- sam_contrast(x, sim$groups, w_ab, n_perm = 100)
  expect_false(free$summary$blocked)
  expect_equal(fit$table$d, free$table$d)
  expect_false(isTRUE(all.equal(fit$table$d_expected, free$table$d_expected)))

  # Batch confounded with group, and a single batch: unrestricted.
  set.seed(4)
  conf <- sam_contrast(x, sim$groups, w_ab, n_perm = 100, block = sim$groups)
  expect_false(conf$summary$blocked)
  expect_equal(conf$table$d_expected, free$table$d_expected)
  set.seed(4)
  one <- sam_contrast(x, sim$groups, w_ab, n_perm = 100, block = rep("b", 10))
  expect_false(one$summary$blocked)
})

test_that("a batch effect widens the unrestricted null but not the blocked one", {
  set.seed(6602)
  p <- 1000
  batch <- rep(c("b1", "b2"), times = 4)
  groups <- rep(c("A", "B"), each = 4)
  x <- matrix(rnorm(p * 8, 20, 0.3), p, dimnames = list(paste0("P", 1:p), NULL))
  x[, batch == "b2"] <- x[, batch == "b2"] + rnorm(p, 0, 1.5) # per-protein shift
  x[1:100, groups == "A"] <- x[1:100, groups == "A"] + 1.5
  blk <- sam_contrast(x, groups, w_ab, n_perm = 1000, block = batch, s0 = 0)
  free <- sam_contrast(x, groups, w_ab, n_perm = 1000, s0 = 0)
  # Unrestricted labellings mix batches unevenly, which widens the null.
  expect_gt(
    stats::sd(free$table$d_expected),
    1.2 * stats::sd(blk$table$d_expected)
  )
  expect_gte(sum(blk$table$sam_called), sum(free$table$sam_called))
})

test_that("the delta table and q-values are internally consistent", {
  dt <- sam_fit$delta_table
  tb <- sam_fit$table
  expect_false(is.unsorted(dt$delta))
  # A larger delta never calls more proteins.
  expect_true(all(diff(dt$n_called) <= 0))
  # Every protein called at a delta lies beyond that delta's cutoffs.
  row <- dt[which.min(abs(dt$n_called - 100)), ]
  called <- tb$sam_delta >= row$delta
  expect_equal(sum(called), row$n_called)
  expect_true(all(tb$d[called] >= row$cut_up | tb$d[called] <= row$cut_down))
  expect_false(any(tb$d[!called] >= row$cut_up | tb$d[!called] <= row$cut_down))
  # q-value = lowest FDR among the deltas at which the protein is called.
  expect_true(all(tb$q.value >= 0 & tb$q.value <= 1))
  expect_equal(tb$q.value[which.max(tb$sam_delta)], min(dt$FDR, na.rm = TRUE))
  expect_equal(tb$sam_called, tb$q.value <= 0.05)
  expect_equal(sam_fit$summary$n_called, sum(tb$sam_called))
  expect_true(sam_fit$summary$FDR <= 0.05)
  expect_true(all(tb$P.Value > 0 & tb$P.Value <= 1))
})

test_that("SAM finds the spiked proteins at roughly the nominal FDR", {
  called <- sam_fit$table$sam_called
  expect_gt(mean(called[sim$truth]), 0.5)
  expect_lt(mean(!sim$truth[called]), 0.15)
  # Under the global null almost nothing is called.
  null <- sam_sim(n_diff = 0, seed = 2861)
  set.seed(77)
  fit_null <- sam_contrast(null$x, null$groups, w_ab, n_perm = 200)
  expect_lt(sum(fit_null$table$sam_called), 10)
})

test_that("a manual delta overrides the target FDR", {
  set.seed(915)
  fit <- sam_contrast(sim$x, sim$groups, w_ab, n_perm = 300, delta = 1.5)
  expect_equal(fit$summary$delta, 1.5)
  expect_equal(fit$table$sam_called, fit$table$sam_delta >= 1.5)
  expect_equal(fit$summary$n_called, sum(fit$table$sam_called))
})

test_that("the Tusher estimator is the mean false calls without pi0", {
  set.seed(915)
  tu <- sam_contrast(
    sim$x,
    sim$groups,
    w_ab,
    n_perm = 300,
    fdr_method = "tusher"
  )
  # Same permutations, same statistic, different FDR estimator.
  expect_equal(tu$table$d, sam_fit$table$d)
  expect_equal(tu$delta_table$n_called, sam_fit$delta_table$n_called)
  expect_true(all(tu$delta_table$n_false >= 0))
  expect_false(isTRUE(all.equal(tu$delta_table$FDR, sam_fit$delta_table$FDR)))
})

test_that("missing values are tolerated and untestable proteins return NA", {
  x <- sim$x
  x[1, 1:5] <- NA # one group entirely missing
  x[2, c(1, 6)] <- NA # still two values per group
  set.seed(915)
  fit <- sam_contrast(x, sim$groups, w_ab, n_perm = 100)
  expect_true(is.na(fit$table$d[1]))
  expect_true(is.na(fit$table$q.value[1]))
  expect_false(fit$table$sam_called[1] %in% TRUE)
  t2 <- stats::t.test(x[2, 2:5], x[2, 7:10], var.equal = TRUE)
  expect_equal(fit$table$SE[2], t2$stderr)
  expect_equal(fit$summary$n_tested, nrow(x) - 1)
  expect_false(anyNA(fit$table$q.value[-1]))
})

test_that("contrasts that do not compare groups are rejected", {
  expect_error(
    sam_contrast(sim$x, sim$groups, c(A = 1, B = 0)),
    "contrast between groups"
  )
})

test_that("run_sam returns the long layout the module expects", {
  cm <- build_contrast_matrix(pq_contrasts, c("m0", "m1", "m2"))
  rs <- run_sam(pq_log2_67, pq_groups, cm, n_perm = 100, seed = 4127)
  expect_equal(nrow(rs$results), nrow(pq_log2_67) * length(pq_contrasts))
  expect_equal(unique(rs$results$comparison), pq_contrasts)
  expect_true(all(
    c(
      "Protein",
      "logFC",
      "AveExpr",
      "d",
      "d_expected",
      "SE",
      "s0",
      "P.Value",
      "q.value",
      "sam_called",
      "comparison"
    ) %in%
      names(rs$results)
  ))
  expect_equal(rs$summary$comparison, pq_contrasts)
  # 3 vs 3: the 20 distinct splits are enumerated.
  expect_true(all(rs$summary$complete_enumeration))
  expect_true(all(rs$summary$n_perm == 20))

  # Same fold change as limma on the same matrix.
  lf <- pq_fit$contrast_fit$coefficients
  first <- rs$results[rs$results$comparison == pq_contrasts[1], ]
  expect_equal(
    first$logFC[match(rownames(lf), first$Protein)],
    unname(lf[, 1])
  )
  # Seeded, hence reproducible.
  expect_equal(
    run_sam(pq_log2_67, pq_groups, cm, n_perm = 100, seed = 4127)$results,
    rs$results
  )
  expect_equal(stat_result_columns(rs$results)$fdr, "q.value")
  expect_equal(stat_result_columns(data.frame(adj.P.Val = 1))$fdr, "adj.P.Val")
})

test_that("the SAM plots build from a results table", {
  skip_if_not_installed("ggplot2")
  withr::local_package("ggplot2")
  cm <- build_contrast_matrix(pq_contrasts[1], c("m0", "m1", "m2"))
  rs <- run_sam(pq_log2_67, pq_groups, cm, n_perm = 50, seed = 1)
  lr <- rs$results
  lr$status <- factor(
    ifelse(lr$sam_called %in% TRUE, "Increased", "Not significant"),
    levels = c("Decreased", "Not significant", "Increased")
  )
  expect_s3_class(ggplot_build(sam_qq_plot(lr, rs$summary)), "ggplot_built")
  expect_s3_class(ggplot_build(sam_scatter_plot(lr)), "ggplot_built")
})
