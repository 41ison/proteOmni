# Smoke test for the MDD machinery in app/modules/mod_PwrQuant.r using a real
# matrix. Loads the statistical helpers only, not the Shiny module.

library(limma)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(stringr)

# Evaluate every top-level expression in the module. Defining a function never
# runs its body, so this pulls in the helpers without starting Shiny.
helper_env <- new.env(parent = globalenv())
for (e in parse("app/modules/mod_PwrQuant.r")) {
  try(eval(e, envir = helper_env), silent = TRUE)
}
needed <- c(
  "bh_effective_alpha",
  "solve_ncp",
  "extract_power_stats",
  "estimate_missingness_model",
  "apply_simulated_missingness",
  "simulate_pipeline_mdd",
  "simulate_mdd_by_replicates",
  "groupwise_imputation"
)
missing_fns <- needed[!vapply(needed, exists, logical(1), envir = helper_env)]
stopifnot(length(missing_fns) == 0)
list2env(as.list(helper_env), envir = globalenv())

# ── Data (comma-separated despite the .tsv extension) ───────────────────────
raw <- read.csv("matrix_human.tsv", check.names = FALSE)
mat <- as.matrix(raw[, -1])
rownames(mat) <- raw[[1]]
group_labels <- sub("_[^_]+$", "", colnames(mat))

log2_mat_raw <- log2(mat)
keep <- apply(log2_mat_raw, 1, function(x) {
  all(tapply(!is.na(x), group_labels, sum) >= 3)
})
log2_mat_raw <- log2_mat_raw[keep, ]

cat(sprintf(
  "Groups: %s | proteins retained: %d of %d | missing: %.1f%%\n",
  paste(
    names(table(group_labels)),
    table(group_labels),
    sep = "x",
    collapse = ", "
  ),
  nrow(log2_mat_raw),
  nrow(mat),
  100 * mean(is.na(log2_mat_raw))
))

# MinProb draws from the left tail, so unlike KNN it is stochastic: seed it or
# the fitted priors (and everything downstream) shift between runs.
set.seed(4817)
log2_mat <- groupwise_imputation(
  log2_mat_raw,
  group_labels,
  method = "minprob",
  verbose = FALSE
)
log2_mat <- normalizeBetweenArrays(as.matrix(log2_mat), method = "cyclicloess")

design <- model.matrix(~ 0 + factor(group_labels))
colnames(design) <- levels(factor(group_labels))
contrast_matrix <- makeContrasts(
  contrasts = c("SW620-SW480", "NOS2-V0"),
  levels = design
)
fit <- eBayes(
  contrasts.fit(lmFit(log2_mat, design), contrast_matrix),
  trend = TRUE
)
cat(sprintf(
  "df.residual = %d | df.prior = %.2f | s2.prior length = %d\n",
  unique(fit$df.residual)[1],
  fit$df.prior[1],
  length(fit$s2.prior)
))

# ── 1. Conditional MDD (diagnostic) ─────────────────────────────────────────
ps <- extract_power_stats(fit, fdr = 0.05, target_power = 0.80)
cat("\n[1] Conditional MDD by contrast:\n")
print(as.data.frame(
  ps |>
    group_by(comparison) |>
    summarise(
      n = n(),
      n_sig_fdr = sum(adj.P.Val <= 0.05),
      alpha_eff = dplyr::first(alpha_eff),
      df_total = dplyr::first(df_total),
      median_MDD = median(Conditional_MDD_Log2FC),
      .groups = "drop"
    )
))

# ── 2. Prospective MDD, fine grid, complete data ────────────────────────────
t0 <- Sys.time()
sim_fine <- simulate_pipeline_mdd(
  fit,
  design,
  contrast_matrix,
  group_labels = group_labels,
  trend = TRUE,
  pi1 = 0.10,
  fdr_cutoff = 0.05,
  target_power = 0.80,
  lfc_grid = seq(0.1, 2.0, by = 0.1),
  n_sim = 20,
  seed = 4817
)
cat(sprintf(
  "\n[2] Grid search (%.1f s):\n",
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
))
print(as.data.frame(sim_fine$mdd))

# ── 3. Bisection must agree with the grid (the sweep relies on it) ──────────
t0 <- Sys.time()
sim_bis <- simulate_pipeline_mdd(
  fit,
  design,
  contrast_matrix,
  group_labels = group_labels,
  trend = TRUE,
  pi1 = 0.10,
  fdr_cutoff = 0.05,
  target_power = 0.80,
  lfc_grid = c(0.05, 3.0),
  search = "bisect",
  bisect_iter = 8,
  n_sim = 20,
  seed = 4817
)
cat(sprintf(
  "\n[3] Bisection (%.1f s):\n",
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
))
print(as.data.frame(sim_bis$mdd))
cat(sprintf(
  "    max |grid - bisect| = %.3f log2FC\n",
  max(abs(sim_fine$mdd$MDD_Log2FC - sim_bis$mdd$MDD_Log2FC))
))

# ── 4. Replicate sweep, complete data vs simulated missingness ──────────────
cm1 <- contrast_matrix[, "SW620-SW480", drop = FALSE]
miss_model <- estimate_missingness_model(log2_mat_raw)
cat(sprintf(
  "\n[4] Missingness model: overall NA rate %.3f\n",
  mean(is.na(log2_mat_raw))
))

n_grid <- c(3, 4, 5, 6, 8, 10, 12, 15)
sweep_args <- list(
  contrast_fit = fit,
  design = design,
  contrast_matrix = cm1,
  n_grid = n_grid,
  trend = TRUE,
  pi1 = 0.10,
  fdr_cutoff = 0.05,
  target_power = 0.80,
  lfc_grid = c(0.05, 3.0),
  search = "bisect",
  bisect_iter = 8,
  n_sim = 10,
  seed = 4817
)

t0 <- Sys.time()
sweep_complete <- do.call(simulate_mdd_by_replicates, sweep_args)
cat(sprintf(
  "    complete-data sweep: %.1f s\n",
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
))

t0 <- Sys.time()
sweep_missing <- do.call(
  simulate_mdd_by_replicates,
  c(
    sweep_args,
    list(
      miss_model = miss_model,
      impute_method = "minprob",
      min_valid_per_group = 3L
    )
  )
)
cat(sprintf(
  "    missingness sweep:   %.1f s\n",
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
))

sweep <- bind_rows(
  sweep_complete |> mutate(scenario = "Complete data"),
  sweep_missing |> mutate(scenario = "Missingness + MinProb imputation")
)
cat("\n[4] Replicate sweep:\n")
print(as.data.frame(
  sweep |>
    mutate(
      MDD_log2 = round(MDD_Log2FC, 3),
      fold = round(2^MDD_Log2FC, 3),
      realised_missing = round(realised_missing, 4)
    ) |>
    select(
      scenario,
      n_per_group,
      n_samples,
      MDD_log2,
      fold,
      realised_missing,
      note
    ) |>
    arrange(scenario, n_per_group)
))

# Cost of missing data at matched replicate counts.
cat("\n[4b] MDD penalty from missingness:\n")
print(as.data.frame(
  sweep |>
    select(scenario, n_per_group, MDD_Log2FC) |>
    tidyr::pivot_wider(names_from = scenario, values_from = MDD_Log2FC) |>
    mutate(
      penalty_log2 = round(
        `Missingness + MinProb imputation` - `Complete data`,
        3
      ),
      pct_worse = round(
        100 * (`Missingness + MinProb imputation` / `Complete data` - 1),
        1
      )
    )
))

# 1/sqrt(n) reference: MDD should scale close to, but slightly worse than, this.
ref <- sweep |> filter(scenario == "Complete data") |> arrange(n_per_group)
cat("\n[4c] Observed vs 1/sqrt(n) scaling (anchored at n = 6):\n")
anchor <- ref$MDD_Log2FC[ref$n_per_group == 6]
print(data.frame(
  n_per_group = ref$n_per_group,
  observed = round(ref$MDD_Log2FC, 3),
  sqrt_law = round(anchor * sqrt(6 / ref$n_per_group), 3)
))

# Power is monotone in n, so the MDD must be non-increasing.
for (s in unique(sweep$scenario)) {
  d <- sweep |> filter(scenario == s) |> arrange(n_per_group)
  stopifnot(all(diff(d$MDD_Log2FC) <= 1e-8))
}
cat("\nMonotonicity in n: OK\n")

# ── 5. Reproducibility ──────────────────────────────────────────────────────
sweep_rep <- do.call(simulate_mdd_by_replicates, sweep_args)
stopifnot(isTRUE(all.equal(sweep_complete$MDD_Log2FC, sweep_rep$MDD_Log2FC)))
cat("Seeded reproducibility: OK\n")
