# Run the mod_PwrQuant unit tests.
#
# From the repository root:
#   Rscript app/tests/PwrQuant/run_tests.R
# or, in an R session (repo root or app/ as working directory):
#   source("app/tests/PwrQuant/run_tests.R")
#   testthat::test_dir("app/tests/PwrQuant")
#
# Fixture: app/training/mtx_maxLFQ.tsv (read in place, not copied).

local({
  # 1. Sourced: use the script's own location
  this_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error = function(e) NULL
  )
  # 2. Rscript: use --file=
  if (is.null(this_dir)) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
    if (length(file_arg)) this_dir <- dirname(normalizePath(file_arg))
  }
  # 3. Run line by line: look relative to the working directory
  if (is.null(this_dir)) {
    cands <- c("app/tests/PwrQuant", "tests/PwrQuant", ".")
    hit <- cands[file.exists(file.path(cands, "helper-pwrquant.R"))]
    if (!length(hit)) {
      stop(
        "Cannot find the PwrQuant test folder from ", getwd(),
        ". Set the working directory to the repo root or app/."
      )
    }
    this_dir <- hit[1]
  }
  testthat::test_dir(this_dir, reporter = "progress", stop_on_failure = TRUE)
})
