# ── filter_diann_report ─────────────────────────────────────────────────────
# Row 1 passes everything; rows 2-6 each fail exactly one filter; row 7 has a
# missing Q.Value.
qc_report <- data.frame(
  id = 1:7,
  Q.Value = c(0.001, 0.02, 0.001, 0.001, 0.001, 0.001, NA),
  PG.Q.Value = c(0.001, 0.001, 0.02, 0.001, 0.001, 0.001, 0.001),
  Lib.PG.Q.Value = c(0.001, 0.001, 0.001, 0.02, 0.001, 0.001, 0.001),
  Proteotypic = c(1, 1, 1, 1, 0, 1, 1),
  PG.MaxLFQ.Quality = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.4, 0.9)
)

test_that("filter_diann_report applies every default filter", {
  expect_equal(filter_diann_report(qc_report)$id, 1L)
})

test_that("each filter of filter_diann_report can be switched off on its own", {
  expect_equal(
    filter_diann_report(qc_report, q_value_cutoff = NULL)$id,
    c(1L, 2L, 7L)
  )
  expect_equal(
    filter_diann_report(qc_report, pg_q_value_cutoff = NULL)$id,
    c(1L, 3L)
  )
  expect_equal(
    filter_diann_report(qc_report, lib_pg_q_value_cutoff = NULL)$id,
    c(1L, 4L)
  )
  expect_equal(
    filter_diann_report(qc_report, proteotypic_only = FALSE)$id,
    c(1L, 5L)
  )
  expect_equal(
    filter_diann_report(qc_report, lfq_quality_cutoff = NULL)$id,
    c(1L, 6L)
  )
  expect_equal(
    filter_diann_report(qc_report, lfq_quality_cutoff = 0)$id,
    c(1L, 6L)
  )
})

test_that("filter_diann_report cutoffs are inclusive", {
  edge <- qc_report[1, ]
  edge$Q.Value <- 0.01
  edge$PG.MaxLFQ.Quality <- 0.5
  expect_equal(nrow(filter_diann_report(edge)), 1L)
})

test_that("filter_diann_report skips filters whose column is absent", {
  bare <- qc_report["id"]
  expect_equal(filter_diann_report(bare), bare)
  expect_equal(nrow(filter_diann_report(qc_report[0, ])), 0L)
})

# ── compute_ibaq ────────────────────────────────────────────────────────────
# Theoretical peptides (trypsin, 6-30 aa, no missed cleavage):
#   P1TEST 2, Q2TEST 3, plain_id 0, P3TEST 1
# Rows are deliberately out of run order; row 7 fails the Q.Value filter.
ibaq_report <- data.frame(
  Run = c("B", "B", "A", "A", "A", "B", "A", "A"),
  Protein.Group = c(
    "P1TEST",
    "P1TEST",
    "P1TEST",
    "Q2TEST",
    "plain_id",
    "P1TEST;Q2TEST",
    "Q2TEST",
    "NOTINFASTA"
  ),
  Genes = c("ALB", "ALB", "ALB", "XYZ1", "", "ALB;XYZ1", "XYZ1", "NIF"),
  Precursor.Normalised = c(100, 300, 50, 90, 10, 70, 999, 5),
  Precursor.Quantity = c(10, 30, 6, 9, 1, 7, 99, 1),
  Q.Value = c(rep(0.001, 6), 0.5, 0.001)
)

test_that("compute_ibaq divides summed intensity by theoretical peptide count", {
  ibaq <- compute_ibaq(ibaq_report, fa_seqs, log2_transform = FALSE)
  expect_s3_class(ibaq, "tbl_df")
  expect_equal(
    names(ibaq),
    c("protein_id", "Genes", "n_theoretical_peptides", "A", "B")
  )
  expect_equal(ibaq$protein_id, c("P1TEST", "Q2TEST"))
  expect_equal(ibaq$Genes, c("ALB", "XYZ1"))
  expect_equal(ibaq$n_theoretical_peptides, c(2L, 3L))
  expect_equal(ibaq$A, c(50 / 2, 90 / 3))
  # Q2TEST was not quantified in run B.
  expect_equal(ibaq$B, c((100 + 300) / 2, NA))
})

test_that("compute_ibaq reports log2(iBAQ + 1) by default", {
  raw <- compute_ibaq(ibaq_report, fa_seqs, log2_transform = FALSE)
  logged <- compute_ibaq(ibaq_report, fa_seqs)
  expect_equal(logged$A, log2(raw$A + 1))
  expect_equal(logged$B, log2(raw$B + 1))
})

test_that("compute_ibaq lists protein groups absent from the FASTA", {
  ibaq <- compute_ibaq(ibaq_report, fa_seqs)
  expect_setequal(
    attr(ibaq, "unmatched_proteins"),
    c("P1TEST;Q2TEST", "NOTINFASTA")
  )
  # plain_id is in the FASTA but has no observable peptide: dropped, not
  # reported as unmatched, and never a division by zero.
  expect_false("plain_id" %in% ibaq$protein_id)
  expect_true(all(
    is.finite(unlist(ibaq[c("A", "B")])) | is.na(unlist(ibaq[c("A", "B")]))
  ))
})

test_that("compute_ibaq long output has one row per protein and sample", {
  long <- compute_ibaq(
    ibaq_report,
    fa_seqs,
    log2_transform = FALSE,
    output_long = TRUE
  )
  expect_equal(
    names(long),
    c("protein_id", "Genes", "n_theoretical_peptides", "sample", "iBAQ")
  )
  expect_equal(long$protein_id, c("P1TEST", "P1TEST", "Q2TEST"))
  expect_equal(long$sample, c("A", "B", "A"))
  expect_equal(long$iBAQ, c(25, 200, 30))
  expect_setequal(
    attr(long, "unmatched_proteins"),
    c("P1TEST;Q2TEST", "NOTINFASTA")
  )
})

test_that("compute_ibaq passes the digest settings through", {
  ibaq <- compute_ibaq(
    ibaq_report,
    fa_seqs,
    max_missed_cleavages = 1,
    log2_transform = FALSE
  )
  expect_equal(ibaq$n_theoretical_peptides, c(4L, 5L))
  expect_equal(ibaq$A, c(50 / 4, 90 / 5))

  # Peptides of 8-22 aa: P1TEST keeps only AAAAAARPEPTIDEK, Q2TEST has none.
  narrow <- compute_ibaq(
    ibaq_report,
    fa_seqs,
    min_peptide_len = 8,
    max_peptide_len = 22,
    log2_transform = FALSE
  )
  expect_equal(narrow$protein_id, "P1TEST")
  expect_equal(narrow$A, 50)
})

test_that("compute_ibaq honours intensity_col, id_col and the quality filters", {
  qty <- compute_ibaq(
    ibaq_report,
    fa_seqs,
    intensity_col = "Precursor.Quantity",
    log2_transform = FALSE
  )
  expect_equal(qty$A, c(6 / 2, 9 / 3))
  expect_equal(qty$B, c(40 / 2, NA))

  # With the Q.Value filter off, row 7 contributes to Q2TEST in run A.
  unfiltered <- compute_ibaq(
    ibaq_report,
    fa_seqs,
    q_value_cutoff = NULL,
    log2_transform = FALSE
  )
  expect_equal(unfiltered$A, c(25, (90 + 999) / 3))

  # Protein IDs taken from another column.
  renamed <- ibaq_report
  names(renamed)[names(renamed) == "Protein.Group"] <- "Accession"
  by_acc <- compute_ibaq(
    renamed,
    fa_seqs,
    id_col = "Accession",
    log2_transform = FALSE
  )
  expect_equal(by_acc$protein_id, c("P1TEST", "Q2TEST"))
  expect_equal(by_acc$A, c(25, 30))
})

# Regression test: an id_col that is also an annotation column was selected
# twice and failed with "Element `Genes` doesn't exist".
test_that("compute_ibaq accepts an annotation column as id_col", {
  report <- cbind(
    ibaq_report,
    Protein.Names = paste0(ibaq_report$Genes, "_TEST")
  )
  by_gene <- compute_ibaq(
    report,
    setNames(fa_seqs[c("P1TEST", "Q2TEST")], c("ALB", "XYZ1")),
    id_col = "Genes",
    log2_transform = FALSE
  )
  expect_equal(
    names(by_gene),
    c("protein_id", "Protein.Names", "n_theoretical_peptides", "A", "B")
  )
  expect_equal(by_gene$protein_id, c("ALB", "XYZ1"))
  expect_equal(by_gene$Protein.Names, c("ALB_TEST", "XYZ1_TEST"))
  expect_equal(by_gene$A, c(25, 30))
  expect_equal(by_gene$B, c(200, NA))

  by_name <- compute_ibaq(
    report,
    setNames(fa_seqs[c("P1TEST", "Q2TEST")], c("ALB_TEST", "XYZ1_TEST")),
    id_col = "Protein.Names",
    log2_transform = FALSE,
    output_long = TRUE
  )
  expect_equal(
    names(by_name),
    c("protein_id", "Genes", "n_theoretical_peptides", "sample", "iBAQ")
  )
  expect_equal(by_name$iBAQ, c(25, 200, 30))
})

test_that("compute_ibaq works from a FASTA file end to end", {
  expect_equal(
    compute_ibaq(ibaq_report, ibaq_fasta_sequences(local_fasta())),
    compute_ibaq(ibaq_report, fa_seqs)
  )
})

test_that("compute_ibaq fails informatively on unusable input", {
  expect_error(
    compute_ibaq(ibaq_report[c("Run", "Genes")], fa_seqs),
    "Missing required columns.*Protein.Group, Precursor.Normalised"
  )
  expect_error(
    compute_ibaq(ibaq_report, fa_seqs, q_value_cutoff = 0),
    "No rows passed the quality filters"
  )
})
