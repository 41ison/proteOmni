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
