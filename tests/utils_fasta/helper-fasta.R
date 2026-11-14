# Shared setup for the utils_fasta.r unit tests.
#
# testthat sources every helper-*.R file before the test files, so everything
# defined here (FASTA / iBAQ helpers + fixtures) is visible to every test.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Walk up from the test folder to the app directory (the one containing
# modules/utils_fasta.r), so the suite works wherever the folder is placed.
find_app_dir <- function(start) {
  dir <- normalizePath(start, mustWork = TRUE)
  repeat {
    for (cand in c(dir, file.path(dir, "app"))) {
      if (file.exists(file.path(cand, "modules", "utils_fasta.r"))) {
        return(cand)
      }
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      stop("Could not locate modules/utils_fasta.r above ", start)
    }
    dir <- parent
  }
}

fa_app_dir <- find_app_dir(testthat::test_path())

# Sourcing only *defines* functions; the Shiny module in the file is not run.
source(file.path(fa_app_dir, "modules", "utils_fasta.r"), local = environment())

# ── Fixture: four small proteins with hand-countable tryptic digests ────────
#   P1TEST   MK | AAAAAARPEPTIDEK | GGGGGGR   (R-P is not cleaved; wrapped)
#   Q2TEST   AAAAAAK | CCCCCCK | DDDDDDK      (lower case in the file)
#   plain_id MK                               (no peptide of length >= 6)
#   P3TEST   MR | GGGGGGR                     (shares GGGGGGR with P1TEST)
fa_lines <- c(
  ">sp|P1TEST|ALBU_TEST Test albumin OS=Testus GN=ALB PE=1 SV=1",
  "MKAAAAAARPEPT",
  "IDEKGGGGGGR",
  ">tr|Q2TEST|Q2TEST_TEST Uncharacterized protein OS=Testus GN=XYZ1",
  "aaaaaakcccccckddddddk",
  ">plain_id no gene here",
  "MK",
  ">sp|P3TEST|DUP_TEST Duplicate peptide OS=Testus GN=DUP PE=1 SV=1",
  "MRGGGGGGR"
)

fa_ids <- c("P1TEST", "Q2TEST", "plain_id", "P3TEST")
fa_seqs <- setNames(
  c(
    "MKAAAAAARPEPTIDEKGGGGGGR",
    "AAAAAAKCCCCCCKDDDDDDK",
    "MK",
    "MRGGGGGGR"
  ),
  fa_ids
)

# Writes `lines` to a FASTA file that is deleted when the calling test ends.
local_fasta <- function(lines = fa_lines, env = parent.frame()) {
  withr::local_tempfile(lines = lines, fileext = ".fasta", .local_envir = env)
}
