# Shared setup for the mod_Sage.r unit tests.
#
# testthat sources every helper-*.R file before the test files, so everything
# defined here (module helpers + fixtures) is visible to every test.

# Walk up from the test folder to the app directory (the one containing
# modules/mod_Sage.r), so the suite works wherever the folder is placed.
find_app_dir <- function(start) {
  dir <- normalizePath(start, mustWork = TRUE)
  repeat {
    for (cand in c(dir, file.path(dir, "app"))) {
      if (file.exists(file.path(cand, "modules", "mod_Sage.r"))) {
        return(cand)
      }
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      stop("Could not locate modules/mod_Sage.r above ", start)
    }
    dir <- parent
  }
}

sg_app_dir <- find_app_dir(testthat::test_path())

# Sourcing only *defines* functions (UI, server, helpers); nothing Shiny runs.
# mod_Sage.r attaches its own packages and relies on the FASTA / iBAQ helpers
# of utils_fasta.r, which the app sources first. Any top-level error should
# fail the suite loudly, so no try() here.
source(file.path(sg_app_dir, "modules", "utils_fasta.r"), local = environment())
source(file.path(sg_app_dir, "modules", "mod_Sage.r"), local = environment())

# ── Fixture: FASTA as fasta_sequences_by_header_id() returns it ─────────────
# Named by the first word of the header, which is what Sage reports in the
# `proteins` column. Theoretical tryptic peptides (6-30 aa, no missed
# cleavage) are given on the right.
sg_sequences <- c(
  "sp|P1TEST|ALBU_TEST" = "MKAAAAAARPEPTIDEKGGGGGGR", # 2  (MK is too short)
  "tr|Q2TEST|Q2TEST_TEST" = "AAAAAAKCCCCCCKDDDDDDK", # 3
  "tiny_id" = "MK", # 0
  "plain_id" = "DDDDDDK", # 1
  "contam_sp|C9TEST|TRYP_TEST" = "AAAAAAKCCCCCCK" # 2
)
sg_genes <- c("ALB", "XYZ1", NA, NA, "PRSS")

# ── Fixture: filtered LFQ peptides, as filter_sage_lfq() returns them ───────
# Rows are deliberately out of file order.
sg_lfq <- data.frame(
  filename = c("f2.mzML", "f2.mzML", rep("f1.mzML", 6)),
  proteins = c(
    "sp|P1TEST|ALBU_TEST",
    "sp|P1TEST|ALBU_TEST",
    "sp|P1TEST|ALBU_TEST",
    "tr|Q2TEST|Q2TEST_TEST",
    "tiny_id",
    "plain_id",
    "contam_sp|C9TEST|TRYP_TEST",
    "sp|NOPE|NOT_TEST"
  ),
  peptide = c("p1", "p2", "p1", "p3", "p4", "p5", "p6", "p7"),
  intensity = c(100, 300, 50, 90, 10, 8, 40, 5)
)
