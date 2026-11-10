# Condition names are free text. These tests pin down that contrasts are
# resolved without ever treating a name as an R expression or a regex.

# Regex metacharacters, operators, spaces, a name containing the "-" separator,
# a name that is a substring of another, and one that looks like an alias.
odd_conds <- c(
  "WT (ctrl)",
  "KO-1",
  "KO",
  "a.*b|c+",
  "[5 ug/mL]",
  "grp_1",
  "Protein"
)

test_that("syntactic names reproduce limma::makeContrasts", {
  conds <- c("m0", "m1", "m2")
  specs <- c("m1-m0", "m2-m0", "(m1+m2)/2-m0")
  expect_equal(
    unclass(build_contrast_matrix(specs, conds)),
    unclass(limma::makeContrasts(contrasts = specs, levels = conds))
  )
})

test_that("every pair of free-form names resolves to the right groups", {
  pairs <- contrast_pair_labels(odd_conds)
  expect_equal(nrow(pairs), length(odd_conds) * (length(odd_conds) - 1))

  expect_equal(pairs$label, paste(pairs$group_a, "vs", pairs$group_b))
  # Unlike "-", the separator keeps hyphenated names distinguishable.
  expect_false(
    anyDuplicated(contrast_pair_labels(c("a-b", "c", "a", "b-c"))$label) > 0
  )

  cm <- build_contrast_matrix(pairs$label, odd_conds)
  expect_equal(rownames(cm), odd_conds)
  expect_equal(colnames(cm), pairs$label)
  for (i in seq_len(nrow(pairs))) {
    expect_equal(cm[pairs$group_a[i], i], 1)
    expect_equal(cm[pairs$group_b[i], i], -1)
    expect_equal(sum(cm[, i] != 0), 2)
  }

  groups <- contrast_group_table(cm)
  expect_equal(groups$comparison, pairs$label)
  expect_equal(groups$group_a, pairs$group_a)
  expect_equal(groups$group_b, pairs$group_b)
})

test_that("typed expressions work with free-form names", {
  cm <- build_contrast_matrix(
    c("KO-1 - WT (ctrl)", "(KO-1+KO)/2-WT (ctrl)"),
    odd_conds
  )
  expect_equal(unname(cm[c("WT (ctrl)", "KO-1", "KO"), 1]), c(-1, 1, 0))
  expect_equal(unname(cm[c("WT (ctrl)", "KO-1", "KO"), 2]), c(-1, 0.5, 0.5))
  # Not a plain A - B, so no single group on either side.
  expect_true(is.na(contrast_group_table(cm)$group_a[2]))
})

test_that("unresolvable or unsafe contrasts are rejected, not evaluated", {
  expect_error(
    build_contrast_matrix("KO-Missing", odd_conds),
    "Could not interpret"
  )
  expect_error(build_contrast_matrix("KO-KO", odd_conds), "Could not interpret")
  marker <- tempfile()
  expect_error(
    build_contrast_matrix(sprintf("file.create('%s')", marker), odd_conds),
    "Could not interpret"
  )
  expect_false(file.exists(marker))
  expect_error(build_contrast_matrix("A-", c("A", "")), "cannot be empty")
})

test_that("the limma pipeline and the MDD simulation run on free-form names", {
  renamed <- setNames(c("WT (ctrl)", "KO-1", "a.*b|c+"), c("m0", "m1", "m2"))
  groups <- unname(renamed[pq_groups])
  contrasts <- c("KO-1-WT (ctrl)", "a.*b|c+-WT (ctrl)", "a.*b|c+-KO-1")

  mat <- pq_log2_complete[seq_len(400), ]
  fit <- fit_pwrquant(mat, groups, contrasts = contrasts)
  ref <- fit_pwrquant(mat, pq_groups)

  ps <- extract_power_stats(fit$contrast_fit)
  ps_ref <- extract_power_stats(ref$contrast_fit)
  expect_equal(unique(ps$comparison), contrasts)
  expect_equal(ps$logFC, ps_ref$logFC)
  expect_equal(ps$adj.P.Val, ps_ref$adj.P.Val)

  sweep <- simulate_mdd_by_replicates(
    fit$contrast_fit,
    fit$design,
    fit$contrast_matrix,
    n_grid = c(3, 4),
    lfc_grid = c(0.5, 4),
    n_sim = 1L,
    seed = 6113
  )
  expect_setequal(sweep$comparison, contrasts)
  expect_equal(nrow(sweep), 2L * length(contrasts))
})

test_that("CV table keeps free-form condition names intact", {
  groups <- rep(c("CV_x (a|b)", "y+z"), each = 3)
  cv <- compute_cv_mtx(2^pq_log2_complete[1:5, 1:6], groups)
  expect_equal(sub("^CV_", "", names(cv)), unique(groups))
})
