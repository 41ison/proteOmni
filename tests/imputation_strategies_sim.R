# Realised FDR and power of each PwrQuant imputation strategy, by simulation.
#
# Source of the table in the "Imputation" section of mod_PwrQuantDocs.r.
# Pushes simulated matrices through the module's own simulate_pipeline_mdd():
# variances drawn from the eBayes priors of the real matrix, known fold changes
# spiked into a fraction of proteins, intensity-dependent dropout fitted to the
# real matrix (optionally scaled up), then imputation -> lmFit -> eBayes -> BH.
#
# Run from the repository root:
#   Rscript app/tests/imputation_strategies_sim.R
# Full settings take ~15-25 min (~9 min without missForest). For a quick look:
#   PQ_SIM_NSIM=3 PQ_SIM_SKIP_MISSFOREST=1 Rscript app/tests/imputation_strategies_sim.R

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tidyr)
})

# ── Settings ────────────────────────────────────────────────────────────────
n_sim <- as.integer(Sys.getenv("PQ_SIM_NSIM", "15"))
skip_missforest <- nzchar(Sys.getenv("PQ_SIM_SKIP_MISSFOREST"))
seed <- 6153L
pi1 <- 0.10
fdr_cutoff <- 0.05
lfc_grid <- c(1.0, 2.0)
contrast <- "m1-m0"
# Dropout levels: NA = the matrix's own model; numbers = mean dropout
# probability before the per-group cap (>= 2 of 3 observed per group), so the
# realised rate is lower (about 11% and 18%).
dropout_targets <- c(NA, 0.15, 0.30)

# ── Module helpers (defining functions only; nothing Shiny runs) ────────────
helper_env <- new.env(parent = globalenv())
for (e in parse("app/modules/mod_PwrQuant.r")) {
  try(eval(e, envir = helper_env), silent = TRUE)
}
stopifnot(all(vapply(
  c(
    "groupwise_imputation",
    "select_knn_k",
    "estimate_missingness_model",
    "simulate_pipeline_mdd"
  ),
  exists,
  logical(1),
  envir = helper_env
)))

# ── Data: same ingestion and filter as the PwrQuant test fixture ────────────
df <- data.table::fread(
  "app/training/mtx_maxLFQ.tsv",
  header = TRUE,
  quote = ""
)
df <- df[!is.na(df[[1]]), ]
mtx <- as.matrix(as.data.frame(df)[, -1, drop = FALSE])
mode(mtx) <- "numeric"
rownames(mtx) <- make.unique(gsub('^"|"$', "", as.character(df[[1]])))
groups <- sub("_r[0-9]+$", "", colnames(mtx))
log2_mat <- log2(mtx + 1)

# >= 60% valid per group, i.e. at least 2 of 3 replicates.
keep <- Reduce(
  `&`,
  lapply(unique(groups), function(g) {
    rowMeans(!is.na(log2_mat[, groups == g, drop = FALSE])) >= 0.6
  })
)
log2_mat <- log2_mat[keep, ]

# ── Priors and missingness model from the real matrix ───────────────────────
# Priors come from the least-squares fit with NAs kept, so they are not
# shaped by any of the imputation methods being compared.
design <- model.matrix(~ 0 + factor(groups))
colnames(design) <- levels(factor(groups))
cm <- makeContrasts(contrasts = contrast, levels = design)
fit <- eBayes(contrasts.fit(lmFit(log2_mat, design), cm), trend = TRUE)
miss_model <- helper_env$estimate_missingness_model(log2_mat)

# Shift the GLM intercept so the mean dropout probability hits `target`,
# keeping the fitted intensity dependence (slope).
scaled_model <- function(target) {
  if (is.na(target)) {
    return(miss_model)
  }
  b <- miss_model$coefficients
  b[1] <- uniroot(
    function(a) mean(plogis(a + b[2] * fit$Amean)) - target,
    c(-50, 50)
  )$root
  modifyList(miss_model, list(coefficients = b))
}

# KNN k chosen on the real matrix, as the app does.
knn_k <- helper_env$select_knn_k(log2_mat, groups)

arms <- data.frame(
  arm = c(
    "none (ls)",
    "minprob",
    sprintf("knn (CV, k=%d)", knn_k),
    "knn (k=5)",
    "knn (k=10)",
    "bpca",
    "missforest"
  ),
  method = c(NA, "minprob", "knn", "knn", "knn", "bpca", "missforest"),
  k = c(NA, NA, knn_k, 5, 10, NA, NA),
  stringsAsFactors = FALSE
)
if (skip_missforest) {
  arms <- arms[arms$method %in% c(NA, "minprob", "knn", "bpca"), ]
}

run_arm <- function(i, target) {
  a <- arms[i, ]
  s <- helper_env$simulate_pipeline_mdd(
    fit,
    design,
    cm,
    group_labels = groups,
    trend = TRUE,
    pi1 = pi1,
    fdr_cutoff = fdr_cutoff,
    lfc_grid = lfc_grid,
    n_sim = n_sim,
    seed = seed,
    miss_model = scaled_model(target),
    impute_method = if (is.na(a$method)) NULL else a$method,
    knn_k = if (is.na(a$k)) NULL else a$k
  )
  data.frame(
    arm = a$arm,
    target = target,
    lfc = s$curve$True_LFC,
    power = s$curve$Power,
    fdp = s$curve$FDP,
    miss = s$curve$Missing_Rate
  )
}

t0 <- Sys.time()
results <- bind_rows(lapply(dropout_targets, function(tg) {
  bind_rows(lapply(seq_len(nrow(arms)), run_arm, target = tg))
}))
message(sprintf(
  "Done in %.1f min (n_sim = %d).",
  as.numeric(difftime(Sys.time(), t0, units = "mins")),
  n_sim
))

# Realised FDR / power (%), one row per dropout level and true log2FC.
summary_table <- results |>
  group_by(target, lfc) |>
  mutate(missing_pct = round(100 * mean(miss), 1)) |>
  ungroup() |>
  mutate(cell = sprintf("%.1f / %.1f", 100 * fdp, 100 * power)) |>
  select(missing_pct, lfc, arm, cell) |>
  pivot_wider(names_from = arm, values_from = cell) |>
  arrange(lfc, missing_pct)

cat(sprintf(
  "\nmtx_maxLFQ: %d proteins, %s | contrast %s | pi1 = %.2f | BH %.0f%% | KNN k = %d\n",
  nrow(log2_mat),
  paste(names(table(groups)), table(groups), sep = "x", collapse = ", "),
  contrast,
  pi1,
  100 * fdr_cutoff,
  knn_k
))
cat("Cells: realised FDR / power (%)\n\n")
print(as.data.frame(summary_table), row.names = FALSE)
