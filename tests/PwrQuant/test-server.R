# Exercises the real ingestion reactives of PwrQuant_server() via
# shiny::testServer, so the fixture reader in helper-pwrquant.R cannot drift
# from the app's own parsing. Packages are attached only for this file
# (withr::local_package detaches them afterwards).

server_pkgs <- c(
  "shiny",
  "shinyjs",
  "shinydashboard",
  "tidyverse",
  "DT",
  "data.table"
)

# The upload also renders the QC plots, which emit ggplot warnings about the
# 976 all-NA proteins. They are expected and not what these tests check, so
# each testServer() call is wrapped in suppressWarnings(). (A wrapper function
# would break testServer's quoting of its expression.)

local_server_packages <- function(env = parent.frame()) {
  for (p in server_pkgs) {
    skip_if_not_installed(p)
  }
  for (p in server_pkgs) {
    suppressPackageStartupMessages(withr::local_package(p, .local_envir = env))
  }
}

fixture_upload <- data.frame(
  name = basename(pq_fixture_path),
  size = file.size(pq_fixture_path),
  type = "text/tab-separated-values",
  datapath = pq_fixture_path
)

test_that("uploading mtx_maxLFQ.tsv yields the expected raw and log2 matrices", {
  local_server_packages()
  suppressWarnings(shiny::testServer(PwrQuant_server, {
    session$setInputs(matrix_file = fixture_upload, already_log2 = FALSE)
    expect_equal(raw_matrix_full(), pq_raw)
    expect_equal(raw_matrix(), pq_raw)
    expect_equal(log2_input(), pq_log2)
  }))
})

test_that("already_log2 = TRUE bypasses the transform", {
  local_server_packages()
  suppressWarnings(shiny::testServer(PwrQuant_server, {
    session$setInputs(matrix_file = fixture_upload, already_log2 = TRUE)
    expect_equal(log2_input(), pq_raw)
  }))
})

test_that("limma runs end to end with free-form condition and display names", {
  local_server_packages()
  conds <- c(m0 = "WT (ctrl)", m1 = "KO-1", m2 = "a.*b|c+")
  suppressWarnings(shiny::testServer(PwrQuant_server, {
    session$setInputs(matrix_file = fixture_upload, already_log2 = FALSE)
    meta <- metadata_base()
    meta$Condition <- unname(conds[sub("_r[0-9]+$", "", meta$Sample)])
    meta$Display_Name <- paste0(meta$Condition, " [", seq_len(nrow(meta)), "]")
    meta_all(meta)

    contrasts <- c("KO-1 vs WT (ctrl)", "a.*b|c+ vs KO-1")
    session$setInputs(
      contrast_pairs = contrasts,
      limma_method = "robust",
      imputation_method = "minprob",
      norm_method = "none",
      trend_param = "TRUE",
      min_valid_pct = 0,
      sig_fdr = 0.05,
      sig_lfc = 0,
      run_limma = 1
    )
    res <- limma_results_ev()
    lr <- res$limma_results

    expect_equal(colnames(res$design), sort(unname(conds)))
    expect_equal(colnames(res$contrast_matrix), contrasts)
    expect_equal(unique(lr$comparison), contrasts)
    per_contrast <- lr[match(contrasts, lr$comparison), ]
    expect_equal(per_contrast$group_a, c("KO-1", "a.*b|c+"))
    expect_equal(per_contrast$group_b, c("WT (ctrl)", "KO-1"))

    # Imputation-driven = a group of the contrast is entirely missing.
    log2_mat <- res$log2_mat
    all_na <- function(cond) {
      cols <- meta$Sample[meta$Condition == cond]
      rowSums(!is.na(log2_mat[, cols, drop = FALSE])) == 0
    }
    first <- lr[lr$comparison == contrasts[1], ]
    expected <- (all_na("KO-1") | all_na("WT (ctrl)"))[first$Protein]
    expect_equal(first$imputation_driven, unname(expected))
    expect_true(any(first$imputation_driven))
    expect_false(anyNA(lr$imputation_driven))

    expect_equal(unname(get_display_names()), meta$Display_Name)
  }))
})

test_that("every Start limma press logs parameters and metadata assignment", {
  local_server_packages()
  logged <- character(0)
  withCallingHandlers(
    suppressWarnings(shiny::testServer(PwrQuant_server, {
      session$setInputs(matrix_file = fixture_upload, already_log2 = FALSE)
      session$setInputs(exclude_samples = "m2_r3")
      meta <- metadata_base()
      meta$Condition <- sub("_r[0-9]+$", "", meta$Sample)
      meta$Batch <- ifelse(meta$Condition == "m0", "A", "B")
      meta_all(meta)

      # Nothing reads limma_results_ev() here: the press alone must log.
      session$setInputs(
        contrast_pairs = "m1 vs m0",
        limma_method = "ls",
        norm_method = "none",
        trend_param = "TRUE",
        min_valid_pct = 0,
        sig_fdr = 0.05,
        sig_lfc = 0,
        run_limma = 1
      )
      session$setInputs(norm_method = "quantile", run_limma = 2)

      # A press that cannot run is recorded too.
      session$setInputs(contrast_pairs = "nope vs m0", run_limma = 3)
    })),
    message = function(m) {
      logged <<- c(logged, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  runs <- grep("PwrQuant limma parameters", logged, value = TRUE)
  expect_length(runs, 2)
  expect_match(runs[1], "Normalization\\s+: none")
  expect_match(runs[2], "Normalization\\s+: quantile")
  expect_match(runs[2], "PwrQuant sample metadata assignment", fixed = TRUE)
  expect_match(runs[2], "m0_r1\\s+m0_r1\\s+m0\\s+A")
  expect_match(runs[2], "m1_r2\\s+m1_r2\\s+m1\\s+B")
  expect_match(runs[2], "Removed samples : m2_r3", fixed = TRUE)
  expect_no_match(runs[2], "m2_r3\\s+m2_r3")
  expect_length(grep("PwrQuant limma not run", logged), 1)
})

test_that("excluding samples drops columns and enforces two kept samples", {
  local_server_packages()
  suppressWarnings(shiny::testServer(PwrQuant_server, {
    session$setInputs(matrix_file = fixture_upload, already_log2 = FALSE)
    session$setInputs(exclude_samples = c("m2_r3", "m0_r1"))
    expect_equal(
      colnames(raw_matrix()),
      setdiff(colnames(pq_raw), c("m2_r3", "m0_r1"))
    )

    session$setInputs(exclude_samples = colnames(pq_raw)[-1])
    expect_error(raw_matrix(), "At least two samples")
  }))
})
