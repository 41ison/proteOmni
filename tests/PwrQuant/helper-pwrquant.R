# Shared setup for the mod_PwrQuant unit tests.
#
# testthat sources every helper-*.R file before the test files, so everything
# defined here (module helpers + fixtures) is visible to every test.

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(limma)
})

# Walk up from the test folder to the app directory (the one containing
# modules/mod_PwrQuant.r), so the suite works wherever the folder is placed.
find_app_dir <- function(start) {
  dir <- normalizePath(start, mustWork = TRUE)
  repeat {
    for (cand in c(dir, file.path(dir, "app"))) {
      if (file.exists(file.path(cand, "modules", "mod_PwrQuant.r"))) {
        return(cand)
      }
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      stop("Could not locate modules/mod_PwrQuant.r above ", start)
    }
    dir <- parent
  }
}

pq_app_dir <- find_app_dir(testthat::test_path())
pq_module_path <- file.path(pq_app_dir, "modules", "mod_PwrQuant.r")
pq_fixture_path <- file.path(pq_app_dir, "training", "mtx_maxLFQ.tsv")

# ── Load the module ─────────────────────────────────────────────────────────
# Sourcing only *defines* functions (UI, server, helpers); nothing Shiny runs.
# Any top-level error should fail the suite loudly, so no try() here.
source(pq_module_path, local = environment())

# ── Fixture: MaxLFQ matrix, 3 conditions (m0, m1, m2) x 3 replicates ────────
# Mirrors `raw_matrix_full()` in PwrQuant_server(): fread with quote = "",
# strip the literal quotes that leaves on the IDs, make IDs unique, coerce to
# numeric. test-server.R checks this copy against the real reactive.
read_pwrquant_matrix <- function(path) {
  df <- data.table::fread(path, header = TRUE, quote = "")
  df <- df[!is.na(df[[1]]), ]
  ids <- make.unique(gsub('^"|"$', "", as.character(df[[1]])))
  mtx <- as.matrix(as.data.frame(df)[, -1, drop = FALSE])
  mode(mtx) <- "numeric"
  rownames(mtx) <- ids
  mtx
}

pq_raw <- read_pwrquant_matrix(pq_fixture_path)
pq_groups <- sub("_r[0-9]+$", "", colnames(pq_raw))
# Server default (already_log2 unticked): log2(x + 1)
pq_log2 <- log2(pq_raw + 1)

# Groupwise valid-value filter, same rule as step 1b of the limma pipeline.
filter_min_valid <- function(log2_mat, group_labels, min_pct) {
  keep <- rep(TRUE, nrow(log2_mat))
  for (g in unique(group_labels)) {
    cols <- which(group_labels == g)
    keep <- keep &
      rowMeans(!is.na(log2_mat[, cols, drop = FALSE])) >= min_pct / 100
  }
  log2_mat[keep, , drop = FALSE]
}

pq_log2_complete <- filter_min_valid(pq_log2, pq_groups, 100)
pq_log2_67 <- filter_min_valid(pq_log2, pq_groups, 60) # >= 2 of 3 per group

pq_contrasts <- c("m1-m0", "m2-m0", "m2-m1")

# Mirrors steps 5-6 of the limma pipeline (design ~ 0 + group, contrasts,
# lmFit, optional removal of unestimable rows, eBayes).
fit_pwrquant <- function(
  log2_mat,
  group_labels,
  contrasts = pq_contrasts,
  trend = TRUE,
  method = "ls"
) {
  design_groups <- factor(group_labels)
  design <- stats::model.matrix(~ 0 + design_groups)
  colnames(design) <- str_remove(colnames(design), "design_groups")
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrasts,
    levels = design
  )
  fit <- limma::lmFit(as.matrix(log2_mat), design, method = method)
  cfit <- limma::contrasts.fit(fit, contrast_matrix)
  if (trend) {
    cfit <- cfit[which(!is.na(cfit$sigma) & !is.na(cfit$Amean)), ]
  }
  list(
    contrast_fit = limma::eBayes(cfit, trend = trend),
    design = design,
    contrast_matrix = contrast_matrix
  )
}

# The full-data fit is reused by several test files; compute it once.
pq_fit <- fit_pwrquant(pq_log2_67, pq_groups)
