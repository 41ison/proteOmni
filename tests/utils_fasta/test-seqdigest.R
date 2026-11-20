# In silico digestion of a single sequence (Digestion tab of mod_SeqProfile.r).

sd_seq <- "MKAAAAAARPEPTIDEKGGGGGGR"

test_that("monoisotopic masses match known values", {
  # PEPTIDE: 799.3600 Da monoisotopic
  expect_equal(seqprofile_mono_mass("PEPTIDE"), 799.35997, tolerance = 1e-7)
  expect_equal(seqprofile_mono_mass("G"), 57.021464 + 18.010565)
  expect_true(is.na(seqprofile_mono_mass("PEPXIDE")))
  skip_if_not_installed("Peptides")
  peps <- c("PEPTIDE", "ACDEFGHIKLMNPQRSTVWY", "MK")
  expect_equal(
    seqprofile_mono_mass(peps),
    as.numeric(Peptides::mw(peps, monoisotopic = TRUE)),
    tolerance = 1e-6
  )
})

test_that("trypsin cuts after K/R but not before P", {
  res <- run_seq_digest(sd_seq, "trypsin")
  expect_equal(res$peptides$sequence, c("MK", "AAAAAARPEPTIDEK", "GGGGGGR"))
  expect_equal(res$peptides$start, c(1L, 3L, 18L))
  expect_equal(res$peptides$end, c(2L, 17L, 24L))
  expect_equal(res$peptides$length, c(2L, 15L, 7L))
  expect_equal(res$stats$n_sites, 2L)
  expect_equal(paste(res$peptides$sequence, collapse = ""), sd_seq)
})

test_that("each preset enzyme applies its rule", {
  seqs <- function(x, enzyme) run_seq_digest(x, enzyme)$peptides$sequence
  expect_equal(seqs("AAEAAEPAA", "gluc"), c("AAE", "AAEPAA"))
  expect_equal(seqs("AAKAAKA", "lysn"), c("AA", "KAA", "KA"))
  expect_equal(seqs("AADAAD", "aspn"), c("AA", "DAA", "D"))
  expect_equal(seqs("SFPSWSYSLS", "chymotrypsin"), c("SFPSW", "SY", "SL", "S"))
  expect_equal(seqs("SSASGPSVS", "elastase"), c("SSA", "SGPSV", "S"))
  expect_setequal(
    unname(seqprofile_enzyme_choices()),
    c(names(SEQPROFILE_ENZYMES), "custom")
  )
})

test_that("custom rules cut on either side and honour blocking residues", {
  c_side <- run_seq_digest(
    "AADAAEPAA",
    "custom",
    list(residues = "d, e", terminus = "C", blocked_by = "P")
  )
  expect_equal(c_side$peptides$sequence, c("AAD", "AAEPAA"))
  expect_match(c_side$rule, "C-terminal to D/E, except when followed by P")

  n_side <- run_seq_digest(
    "AADAAEPAA",
    "custom",
    list(residues = "DE", terminus = "N")
  )
  expect_equal(n_side$peptides$sequence, c("AA", "DAA", "EPAA"))

  expect_error(
    run_seq_digest(sd_seq, "custom", list(residues = "", terminus = "C")),
    class = "seqprofile_bad_param"
  )
  expect_error(
    run_seq_digest(sd_seq, "custom", list(residues = "K1", terminus = "C")),
    class = "seqprofile_bad_param"
  )
})

test_that("missed cleavages add the joined neighbouring peptides", {
  res <- run_seq_digest(sd_seq, "trypsin", list(max_missed = 1))
  expect_equal(nrow(res$peptides), 5L)
  expect_equal(res$peptides$missed_cleavages, c(0L, 1L, 0L, 1L, 0L))
  expect_true("MKAAAAAARPEPTIDEK" %in% res$peptides$sequence)
  # More missed cleavages than sites: capped at the whole sequence
  all_mc <- run_seq_digest(sd_seq, "trypsin", list(max_missed = 5))
  expect_equal(nrow(all_mc$peptides), 6L)
  expect_true(sd_seq %in% all_mc$peptides$sequence)
})

test_that("coverage counts residues in peptides of the chosen length range", {
  res <- run_seq_digest(sd_seq, "trypsin")
  expect_equal(res$peptides$in_range, c(FALSE, TRUE, TRUE))
  expect_equal(res$covered, c(FALSE, FALSE, rep(TRUE, 22)))
  expect_equal(res$stats$coverage_pct, 100 * 22 / 24)

  full <- run_seq_digest(sd_seq, "trypsin", list(min_len = 1, max_len = 100))
  expect_equal(full$stats$coverage_pct, 100)
  none <- run_seq_digest(sd_seq, "trypsin", list(min_len = 20, max_len = 30))
  expect_equal(none$stats$coverage_pct, 0)
  expect_equal(none$stats$n_in_range, 0L)
})

test_that("m/z values and carbamidomethyl Cys follow the neutral mass", {
  plain <- run_seq_digest("PEPTIDEKACCK", "trypsin")$peptides
  expect_equal(plain$sequence, c("PEPTIDEK", "ACCK"))
  # PEPTIDEK: [M+H]+ 928.4622, [M+2H]2+ 464.7347
  expect_equal(plain$mz_1[1], 928.4622, tolerance = 1e-7)
  expect_equal(plain$mz_2[1], 464.7347, tolerance = 1e-7)
  expect_equal(plain$mz_3, (plain$mono_mass + 3 * 1.007276467) / 3)

  cam <- run_seq_digest(
    "PEPTIDEKACCK",
    "trypsin",
    list(carbamidomethyl = TRUE)
  )
  expect_true(cam$params$carbamidomethyl)
  expect_equal(cam$peptides$mono_mass - plain$mono_mass, c(0, 2 * 57.021464))
  expect_equal(cam$peptides$mz_2[2] - plain$mz_2[2], 57.021464)
  expect_false(run_seq_digest("ACK", "trypsin")$params$carbamidomethyl)
})

test_that("Met oxidation adds 15.9949 Da per Met, alone or with CAM", {
  plain <- run_seq_digest("MAMCK", "trypsin")$peptides$mono_mass
  ox <- run_seq_digest("MAMCK", "trypsin", list(oxidation = TRUE))
  expect_true(ox$params$oxidation)
  expect_equal(ox$peptides$mono_mass - plain, 2 * 15.994915)
  both <- run_seq_digest(
    "MAMCK",
    "trypsin",
    list(oxidation = TRUE, carbamidomethyl = TRUE)
  )
  expect_equal(both$peptides$mono_mass - plain, 2 * 15.994915 + 57.021464)
})

test_that("the m/z filter keeps peptides with a 2+ or 3+ ion in range", {
  # GGGGGGR: 2+ 259.1 and 3+ 173.1, both below 350
  # AAAAAARPEPTIDEK: 2+ 755.9, in range
  res <- run_seq_digest(sd_seq, "trypsin", list(mz_min = 350, mz_max = 1500))
  expect_equal(res$peptides$in_range, c(FALSE, TRUE, FALSE))
  expect_equal(res$stats$coverage_pct, 100 * 15 / 24)

  # Only the 3+ ion (504.3) of AAAAAARPEPTIDEK is inside 450-550
  z3 <- run_seq_digest(sd_seq, "trypsin", list(mz_min = 450, mz_max = 550))
  expect_equal(round(z3$peptides$mz_3[2], 1), 504.3)
  expect_equal(z3$peptides$in_range, c(FALSE, TRUE, FALSE))

  # No limits: length filter only
  off <- run_seq_digest(sd_seq, "trypsin", list(mz_min = NA, mz_max = NA))
  expect_equal(off$peptides$in_range, c(FALSE, TRUE, TRUE))
  expect_true(is.na(off$params$mz_min))

  # A peptide without a mass cannot pass the m/z filter
  x <- run_seq_digest("AAAXAAAK", "trypsin", list(mz_min = 100, mz_max = 2000))
  expect_false(x$peptides$in_range)

  expect_error(
    run_seq_digest(sd_seq, "trypsin", list(mz_min = 350, mz_max = NA)),
    class = "seqprofile_bad_param"
  )
  expect_error(
    run_seq_digest(sd_seq, "trypsin", list(mz_min = 1500, mz_max = 350)),
    class = "seqprofile_bad_param"
  )
})

test_that("b/y fragment ions match known values", {
  ions <- seqprofile_fragment_ions("PEPTIDE")
  expect_equal(nrow(ions), 12L)
  expect_equal(ions$ion, c(paste0("b", 1:6), paste0("y", 1:6)))
  expect_true(all(ions$intensity == 1))
  mz <- stats::setNames(ions$mz, ions$ion)
  expect_equal(mz[["b2"]], 227.1026, tolerance = 1e-6)
  expect_equal(mz[["y1"]], 148.0604, tolerance = 1e-6)
  expect_equal(mz[["y2"]], 263.0874, tolerance = 1e-6)
  expect_equal(ions$fragment[ions$ion == "b3"], "PEP")
  expect_equal(ions$fragment[ions$ion == "y3"], "IDE")
  # Complementary pairs: b(i) + y(n - i) = [M+2H]2+ * 2
  m <- seqprofile_mono_mass("PEPTIDE")
  expect_equal(mz[["b2"]] + mz[["y5"]], m + 2 * 1.007276467)

  z2 <- seqprofile_fragment_ions("PEPTIDE", charge = 2)
  expect_equal(z2$mz, (ions$mz + 1.007276467) / 2)
  expect_equal(z2$ion[1], "b1 2+")

  mod <- seqprofile_fragment_ions(
    "ACMK",
    carbamidomethyl = TRUE,
    oxidation = TRUE
  )
  raw <- seqprofile_fragment_ions("ACMK")
  shift <- stats::setNames(mod$mz - raw$mz, raw$ion)
  expect_equal(shift[["b1"]], 0)
  expect_equal(shift[["b2"]], 57.021464)
  expect_equal(shift[["b3"]], 57.021464 + 15.994915)
  expect_equal(shift[["y1"]], 0)
  expect_equal(shift[["y2"]], 15.994915)

  expect_error(
    seqprofile_fragment_ions("K"),
    class = "seqprofile_short_sequence"
  )
  expect_error(
    seqprofile_fragment_ions("PEPXIDE"),
    class = "seqprofile_bad_param"
  )
  skip_if_not_installed("ggplot2")
  expect_no_error(ggplot2::ggplot_build(plot_fragment_spectrum(
    ions,
    "PEPTIDE"
  )))
})

test_that("edge cases and bad input are handled", {
  uncut <- run_seq_digest("AAAAAAA", "trypsin")
  expect_equal(uncut$peptides$sequence, "AAAAAAA")
  expect_equal(uncut$stats$n_sites, 0L)
  expect_equal(run_seq_digest("K", "trypsin")$peptides$sequence, "K")
  expect_error(run_seq_digest("", "trypsin"), class = "seqprofile_empty")
  expect_error(run_seq_digest(sd_seq, "pepsin"), class = "seqprofile_bad_param")
  expect_error(
    run_seq_digest(sd_seq, "trypsin", list(max_missed = NA)),
    class = "seqprofile_bad_param"
  )
  expect_error(
    run_seq_digest(sd_seq, "trypsin", list(min_len = 10, max_len = 5)),
    class = "seqprofile_bad_param"
  )
})

# ── SeqProfile_server(): digestion tab ──────────────────────────────────────
test_that("the module digests only on Digest and builds the coverage map", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("plotly")
  skip_if_not_installed("DT")
  suppressPackageStartupMessages(library(shiny))
  source(
    file.path(fa_app_dir, "modules", "mod_SeqProfile.r"),
    local = environment()
  )

  shiny::testServer(SeqProfile_server, {
    session$setInputs(
      sequence = "",
      enzyme = "trypsin",
      max_missed = 0,
      min_len = 6,
      max_len = 30
    )
    expect_error(digest_result(), class = "shiny.silent.error")
    session$setInputs(digest = 1)
    expect_error(digest_result(), "Paste a sequence")

    session$setInputs(sequence = sd_seq, digest = 2)
    expect_equal(digest_result()$enzyme_name, "Trypsin")
    expect_equal(digest_table_data()$Residues, c(2L, 15L, 7L))
    expect_true(all(
      c("m/z (1+)", "m/z (2+)", "m/z (3+)") %in% names(digest_table_data())
    ))
    expect_match(output$digest_coverage$html, "91.7 %", fixed = TRUE)

    # No spectrum until a peptide row is clicked
    expect_error(spectrum_ions(), "Click a peptide")
    session$setInputs(digest_table_rows_selected = 2)
    expect_equal(attr(spectrum_ions(), "peptide"), "AAAAAARPEPTIDEK")
    expect_equal(nrow(spectrum_ions()), 28L)
    session$setInputs(fragment_charge = "2")
    expect_equal(spectrum_ions()$ion[1], "b1 2+")
    session$setInputs(digest_table_rows_selected = NULL)

    # m/z filter and Met oxidation take effect on the next Digest
    session$setInputs(
      mz_filter = TRUE,
      mz_min = 350,
      mz_max = 1500,
      oxidation = TRUE,
      digest = 3
    )
    expect_equal(digest_table_data()$Counted_For_Coverage, c("no", "yes", "no"))
    expect_true(digest_result()$params$oxidation)
    expect_match(output$digest_coverage$html, "m/z 350-1500", fixed = TRUE)
    session$setInputs(mz_filter = FALSE, oxidation = FALSE)

    # Changing the enzyme does not re-digest until Digest is clicked again
    session$setInputs(
      enzyme = "custom",
      custom_residues = "E",
      custom_terminus = "N",
      custom_blocked = ""
    )
    expect_equal(digest_result()$enzyme_name, "Trypsin")
    session$setInputs(digest = 4)
    expect_equal(digest_result()$peptides$sequence[1], "MKAAAAAARP")
  })

  map <- as.character(seqprofile_coverage_map(
    "AAKAA",
    c(TRUE, TRUE, TRUE, FALSE, FALSE)
  ))
  expect_match(map, ">AAK<", fixed = TRUE)
  expect_match(map, ">AA<", fixed = TRUE)
})
