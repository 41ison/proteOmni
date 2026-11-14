# One row of a compute_ibaq_sage() result, looked up by accession and file so
# the tests do not depend on how the rows happen to be sorted.
ibaq_row <- function(ibaq, protein_id, sample) {
  ibaq[ibaq$protein_id == protein_id & ibaq$sample == sample, , drop = FALSE]
}

test_that("compute_ibaq_sage divides summed intensity by theoretical peptide count", {
  ibaq <- compute_ibaq_sage(sg_lfq, sg_sequences, genes = sg_genes)
  expect_equal(
    names(ibaq),
    c(
      "protein_id",
      "entry_name",
      "gene",
      "n_theoretical_peptides",
      "sample",
      "iBAQ"
    )
  )
  # P1TEST in two files, Q2TEST, plain_id and the contaminant in one each.
  expect_equal(nrow(ibaq), 5L)

  # Untransformed: sum of peptide intensities / theoretical peptides.
  expect_equal(ibaq_row(ibaq, "P1TEST", "f1.mzML")$iBAQ, 50 / 2)
  expect_equal(ibaq_row(ibaq, "P1TEST", "f2.mzML")$iBAQ, (100 + 300) / 2)
  expect_equal(ibaq_row(ibaq, "Q2TEST", "f1.mzML")$iBAQ, 90 / 3)
  expect_equal(ibaq_row(ibaq, "plain_id", "f1.mzML")$iBAQ, 8 / 1)
  expect_equal(ibaq$iBAQ[which(ibaq$entry_name == "TRYP_TEST")], 40 / 2)

  # A protein not quantified in a file has no row for it.
  expect_equal(nrow(ibaq_row(ibaq, "Q2TEST", "f2.mzML")), 0L)
})

test_that("compute_ibaq_sage annotates accession, entry name and gene", {
  ibaq <- compute_ibaq_sage(sg_lfq, sg_sequences, genes = sg_genes)

  p1 <- ibaq_row(ibaq, "P1TEST", "f1.mzML")
  expect_equal(p1$entry_name, "ALBU_TEST")
  expect_equal(p1$gene, "ALB")
  expect_equal(p1$n_theoretical_peptides, 2L)

  q2 <- ibaq_row(ibaq, "Q2TEST", "f1.mzML")
  expect_equal(q2$entry_name, "Q2TEST_TEST")
  expect_equal(q2$gene, "XYZ1")
  expect_equal(q2$n_theoretical_peptides, 3L)

  # Identifiers without the db|accession|name layout have no entry name.
  plain <- ibaq_row(ibaq, "plain_id", "f1.mzML")
  expect_true(is.na(plain$entry_name))
  expect_true(is.na(plain$gene))

  # Sorted by accession, then file. dplyr::arrange() sorts in the C locale;
  # method = "radix" does the same, whereas plain order() follows the session
  # locale and would place the lower-case contaminant ID differently.
  expect_equal(
    order(ibaq$protein_id, ibaq$sample, method = "radix"),
    seq_len(nrow(ibaq))
  )
})

test_that("compute_ibaq_sage leaves gene as NA when no gene names are given", {
  ibaq <- compute_ibaq_sage(sg_lfq, sg_sequences)
  expect_true(all(is.na(ibaq$gene)))
  expect_type(ibaq$gene, "character")
})

test_that("compute_ibaq_sage drops and reports proteins it cannot quantify", {
  ibaq <- compute_ibaq_sage(sg_lfq, sg_sequences, genes = sg_genes)
  expect_equal(attr(ibaq, "unmatched_proteins"), "sp|NOPE|NOT_TEST")
  expect_false("NOPE" %in% ibaq$protein_id)
  # tiny_id is in the FASTA but has no observable peptide: dropped, not
  # reported as unmatched, and never a division by zero.
  expect_false("tiny_id" %in% ibaq$protein_id)
  expect_true(all(is.finite(ibaq$iBAQ)))
  expect_true(all(ibaq$n_theoretical_peptides > 0))
})

test_that("compute_ibaq_sage falls back to the accession when headers differ", {
  # Same accession, different entry name than in the FASTA header.
  d <- data.frame(
    filename = "f1.mzML",
    proteins = c("tr|Q2TEST|RENAMED_TEST", "sp|ZZZZZZ|ALBU_TEST"),
    peptide = c("p1", "p2"),
    intensity = c(60, 7)
  )
  ibaq <- compute_ibaq_sage(d, sg_sequences, genes = sg_genes)
  expect_equal(ibaq$protein_id, "Q2TEST")
  expect_equal(ibaq$entry_name, "RENAMED_TEST")
  expect_equal(ibaq$gene, "XYZ1")
  expect_equal(ibaq$iBAQ, 60 / 3)
  # Only the accession is used for the fallback, not the entry name.
  expect_equal(attr(ibaq, "unmatched_proteins"), "sp|ZZZZZZ|ALBU_TEST")
})

test_that("compute_ibaq_sage keeps contaminants unless asked to exclude them", {
  kept <- compute_ibaq_sage(sg_lfq, sg_sequences)
  expect_true("TRYP_TEST" %in% kept$entry_name)
  expect_equal(attr(kept, "contaminant_proteins"), character())

  dropped <- compute_ibaq_sage(
    sg_lfq,
    sg_sequences,
    exclude_contaminants = TRUE
  )
  expect_false("TRYP_TEST" %in% dropped$entry_name)
  expect_equal(
    attr(dropped, "contaminant_proteins"),
    "contam_sp|C9TEST|TRYP_TEST"
  )
  # The other proteins are unaffected.
  expect_equal(nrow(dropped), nrow(kept) - 1L)
  expect_equal(ibaq_row(dropped, "P1TEST", "f2.mzML")$iBAQ, 200)
})

test_that("compute_ibaq_sage recognises the contaminant prefixes", {
  d <- data.frame(
    filename = "f1.mzML",
    proteins = c(
      "contam_sp|C9TEST|TRYP_TEST",
      "CON_plain_id",
      "Cont_plain_id",
      # "con" inside the identifier is not a prefix
      "plain_id",
      "sp|P1TEST|ALBU_TEST"
    ),
    peptide = paste0("p", 1:5),
    intensity = c(1, 2, 3, 8, 50)
  )
  ibaq <- compute_ibaq_sage(d, sg_sequences, exclude_contaminants = TRUE)
  expect_setequal(
    attr(ibaq, "contaminant_proteins"),
    c("contam_sp|C9TEST|TRYP_TEST", "CON_plain_id", "Cont_plain_id")
  )
  expect_setequal(ibaq$protein_id, c("plain_id", "P1TEST"))

  custom <- compute_ibaq_sage(
    d,
    sg_sequences,
    exclude_contaminants = TRUE,
    contaminant_regex = "ALBU"
  )
  expect_equal(attr(custom, "contaminant_proteins"), "sp|P1TEST|ALBU_TEST")
  expect_false("P1TEST" %in% custom$protein_id)
})

test_that("compute_ibaq_sage stops when no peptide is left", {
  only_contam <- sg_lfq[grepl("^contam_", sg_lfq$proteins), ]
  expect_error(
    compute_ibaq_sage(only_contam, sg_sequences, exclude_contaminants = TRUE),
    "No quantified peptides left"
  )
  expect_error(
    compute_ibaq_sage(sg_lfq[0, ], sg_sequences),
    "No quantified peptides left"
  )
})

test_that("compute_ibaq_sage passes the digest settings through", {
  # One missed cleavage: P1TEST 4, Q2TEST 5 theoretical peptides.
  mc1 <- compute_ibaq_sage(sg_lfq, sg_sequences, max_missed_cleavages = 1)
  expect_equal(ibaq_row(mc1, "P1TEST", "f1.mzML")$n_theoretical_peptides, 4L)
  expect_equal(ibaq_row(mc1, "P1TEST", "f1.mzML")$iBAQ, 50 / 4)
  expect_equal(ibaq_row(mc1, "Q2TEST", "f1.mzML")$iBAQ, 90 / 5)

  # Peptides of 8-22 aa: only AAAAAARPEPTIDEK of P1TEST qualifies.
  narrow <- compute_ibaq_sage(
    sg_lfq,
    sg_sequences,
    min_peptide_len = 8,
    max_peptide_len = 22
  )
  expect_equal(unique(narrow$protein_id), "P1TEST")
  expect_equal(narrow$n_theoretical_peptides, c(1L, 1L))

  # Without the proline rule P1TEST is also cut at R-P: 3 peptides.
  no_pro <- compute_ibaq_sage(
    sg_lfq,
    sg_sequences,
    protease_regex = protease_regex_from_rule("KR")
  )
  expect_equal(ibaq_row(no_pro, "P1TEST", "f2.mzML")$iBAQ, 400 / 3)
})

test_that("compute_ibaq_sage works on the output of filter_sage_lfq()", {
  lfq <- data.frame(
    peptide = c("p1", "p2", "p3", "p4", "p5", "p6"),
    proteins = c(
      "sp|P1TEST|ALBU_TEST",
      "sp|P1TEST|ALBU_TEST",
      # shared peptide: cannot be attributed without protein inference
      "sp|P1TEST|ALBU_TEST;tr|Q2TEST|Q2TEST_TEST",
      "tr|Q2TEST|Q2TEST_TEST",
      "tr|Q2TEST|Q2TEST_TEST",
      "rev_sp|P1TEST|ALBU_TEST"
    ),
    is_decoy = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    q_value = c(0.001, 0.001, 0.001, 0.001, 0.2, 0.001),
    filename = "f1.mzML",
    # 0 = not quantified
    intensity = c(50, 0, 1000, 90, 500, 700)
  )
  ibaq <- compute_ibaq_sage(
    filter_sage_lfq(lfq),
    sg_sequences,
    genes = sg_genes
  )
  expect_equal(ibaq$protein_id, c("P1TEST", "Q2TEST"))
  expect_equal(ibaq$iBAQ, c(50 / 2, 90 / 3))
  expect_equal(attr(ibaq, "unmatched_proteins"), character())
})
