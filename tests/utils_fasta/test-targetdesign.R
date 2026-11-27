# Targeted proteomics assay design (SRM / MRM / PRM) — utils_fasta.r & mod_TargetDesign.r

test_that("digest_proteome_fasta digests FASTA entries with any protease rule", {
  seqs <- c(
    P1 = "MKAAAVVVKPEPTIDER",
    P2 = "MRLPEIVARGGGGGGR"
  )
  dig0 <- digest_proteome_fasta(seqs, max_missed = 0L, min_len = 6L, max_len = 30L)
  expect_true(all(c("AAAVVVK", "PEPTIDER", "LPEIVAR", "GGGGGGR") %in% dig0$peptide))
  expect_false("MK" %in% dig0$peptide)

  dig1 <- digest_proteome_fasta(seqs, max_missed = 1L, min_len = 6L, max_len = 30L)
  expect_true("AAAVVVKPEPTIDER" %in% dig1$peptide)

  expect_equal(nrow(digest_proteome_fasta(character(0))), 0L)
})

test_that("classify_targeted_proteotypicity checks intra-sequence and FASTA specificity with I=L", {
  intra <- classify_targeted_proteotypicity(
    c("LPEIVAR", "PEPTIDER", "LPELVAR"),
    digest_df = NULL,
    equate_il = TRUE
  )
  expect_equal(
    intra$proteotypicity,
    c("Repeated in sequence", "Unique (in sequence)", "Repeated in sequence")
  )

  intra_strict <- classify_targeted_proteotypicity(
    c("LPEIVAR", "LPELVAR"),
    digest_df = NULL,
    equate_il = FALSE
  )
  expect_equal(
    intra_strict$proteotypicity,
    c("Unique (in sequence)", "Unique (in sequence)")
  )

  db <- data.frame(
    protein = c("P02647", "P02647", "P99999"),
    peptide = c("THLAPYSDELR", "LPEIVAR", "LPELVAR"),
    stringsAsFactors = FALSE
  )
  cls <- classify_targeted_proteotypicity(
    c("THLAPYSDELR", "LPEIVAR", "ABSENTPEPK"),
    digest_df = db,
    equate_il = TRUE
  )
  expect_equal(cls$proteotypicity, c("Proteotypic", "Shared", "Not in FASTA"))
  expect_equal(cls$n_proteins, c(1L, 2L, 0L))
  expect_match(cls$mapped_proteins[2], "P02647")
  expect_match(cls$mapped_proteins[2], "P99999")
})

test_that("select_targeted_transitions prioritizes proline-enhanced y-ions above mz_2", {
  # THLAPYSDELR: 11 aa, [M+2H]2+ ~ 651.83, internal Pro at pos 5 -> y7 starts with P
  ions <- select_targeted_transitions("THLAPYSDELR", max_transitions = 4L)
  rec <- ions[ions$recommended, ]
  expect_equal(nrow(rec), 4L)
  expect_true(all(rec$type == "y" & rec$number >= 3L))
  expect_true(all(rec$mz > attr(ions, "mz_2")))
  # y7 ("PYSDELR") starts with Proline and should be ranked #1
  expect_equal(rec$ion[rec$transition_rank == 1L], "y7")

  skip_if_not_installed("ggplot2")
  expect_no_error(ggplot2::ggplot_build(plot_targeted_spectrum(ions)))
})

test_that("evaluate_targeted_candidates applies PMC4830481 rules and ranks tiers", {
  # Construct a multi-peptide protein covering Tier 1, Tier 2, and Tier 3 rules:
  #   1. MK                         (too short, protein N-term) -> Avoid
  #   2. THLAPYSDELR                (11 aa, clean internal tryptic) -> Recommended
  #   3. AVALCDEFGHK                (11 aa, internal alkylated Cys) -> Conditional
  #   4. QGAVALDEFGR                (11 aa, N-term Q + QG deamidation) -> Avoid
  #   5. AVALMDEFGHK                (11 aa, internal Met) -> Conditional
  #   6. AVALDPFGHIK                (11 aa, DP acid cleavage + ragged C-end) -> Avoid
  #   7. KRAVALDEFGHK               (dibasic KR ragged N-end) -> Avoid
  #   8. AVALXDEFGHK                (non-canonical X) -> Avoid
  prot <- paste0(
    "MK",
    "THLAPYSDELR",
    "AVALCDEFGHK",
    "QGAVALDEFGR",
    "AVALMDEFGHK",
    "AVALDPFGHIK",
    "KR",
    "AVALDEFGHK",
    "AVALXDEFGHK"
  )
  dig <- run_seq_digest(
    prot,
    "trypsin",
    list(max_missed = 0L, min_len = 8L, max_len = 22L, carbamidomethyl = TRUE)
  )
  cands <- evaluate_targeted_candidates(dig)
  by_seq <- split(cands, cands$sequence)

  expect_equal(by_seq[["THLAPYSDELR"]]$tier, "Recommended")
  expect_equal(by_seq[["THLAPYSDELR"]]$rank, 1L)

  expect_equal(by_seq[["AVALCDEFGHK"]]$tier, "Conditional")
  expect_match(by_seq[["AVALCDEFGHK"]]$flags, "Cys")

  expect_equal(by_seq[["AVALMDEFGHK"]]$tier, "Conditional")
  expect_match(by_seq[["AVALMDEFGHK"]]$flags, "Met")

  expect_equal(by_seq[["QGAVALDEFGR"]]$tier, "Avoid")
  expect_match(by_seq[["QGAVALDEFGR"]]$flags, "pyroglutamate")

  expect_equal(by_seq[["AVALDPFGHIK"]]$tier, "Avoid")
  expect_match(by_seq[["AVALDPFGHIK"]]$flags, "DP motif")

  expect_equal(by_seq[["AVALDEFGHK"]]$tier, "Avoid")
  expect_match(by_seq[["AVALDEFGHK"]]$flags, "Ragged N-end")

  expect_equal(by_seq[["AVALXDEFGHK"]]$tier, "Avoid")
  expect_match(by_seq[["AVALXDEFGHK"]]$flags, "Non-canonical")

  # Skyline transition export includes Recommended & Conditional candidates
  tr <- build_targeted_transition_list(
    cands,
    protein_name = "TEST_PROT",
    carbamidomethyl = TRUE
  )
  expect_true(all(tr$Tier %in% c("Recommended", "Conditional")))
  expect_true("THLAPYSDELR" %in% tr$PeptideSequence)
  expect_true(any(grepl("C\\[\\+57\\.0215\\]", tr$ModifiedSequence)))
})

test_that("TargetDesign_server sanitizes sequence, warns on B/Z/X/U/O/*, and checks FASTA proteotypicity", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("plotly")
  skip_if_not_installed("DT")
  suppressPackageStartupMessages(library(shiny))
  source(
    file.path(fa_app_dir, "modules", "mod_TargetDesign.r"),
    local = environment()
  )

  tmp_fa <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp_fa), add = TRUE)
  writeLines(
    c(
      ">sp|P02647|APOA1_HUMAN",
      "MKTHLAPYSDELRSHAREDPEPTIDEK",
      ">sp|P02652|APOA2_HUMAN",
      "MKSHAREDPEPTIDEK"
    ),
    tmp_fa
  )

  shiny::testServer(TargetDesign_server, {
    session$setInputs(
      sequence = ">test\nmk thlapy sdelr 123 sharedpeptidek*",
      enzyme = "trypsin",
      max_missed = 0,
      opt_min_len = 8,
      opt_max_len = 22,
      mz_min = 350,
      mz_max = 1250,
      min_transitions = 3,
      max_transitions = 5,
      carbamidomethyl = TRUE,
      allow_alkylated_cys = TRUE,
      penalize_met = TRUE,
      penalize_trp = TRUE,
      oxidation = FALSE,
      equate_il = TRUE,
      tier_filter = "all",
      fragment_charge = "1"
    )
    info_html <- output$seq_info$html
    expect_match(info_html, "Length: 28 aa", fixed = TRUE)
    expect_match(info_html, "Non-canonical residues detected (*)", fixed = TRUE)

    # Upload background proteome FASTA and evaluate
    session$setInputs(
      proteome_fasta = data.frame(
        name = "bg.fasta",
        size = file.info(tmp_fa)$size,
        type = "text/plain",
        datapath = tmp_fa,
        stringsAsFactors = FALSE
      ),
      evaluate = 1L
    )
    res <- design_result()
    expect_true(res$has_fasta)
    expect_equal(res$n_fasta_proteins, 2L)
    by_seq <- split(res$candidates, res$candidates$sequence)
    expect_equal(by_seq[["THLAPYSDELR"]]$proteotypicity, "Proteotypic")
    expect_equal(by_seq[["THLAPYSDELR"]]$tier, "Recommended")
    expect_equal(by_seq[["SHAREDPEPTIDEK"]]$proteotypicity, "Shared")
    expect_equal(by_seq[["SHAREDPEPTIDEK"]]$tier, "Avoid")
  })
})
