# Sequence profiling helpers in utils_fasta.r and the SeqProfile module server.

sp_ins_a <- "GIVEQCCTSICSLYQLENYCN" # human insulin A chain
sp_ins_b <- "FVNQHLCGSHLVEALYLVCGERGFFYTPKT" # human insulin B chain

# ── clean_protein_sequence() ────────────────────────────────────────────────
test_that("raw and FASTA input clean to the same sequence", {
  raw <- clean_protein_sequence(sp_ins_a)
  fa <- clean_protein_sequence(c(
    ">sp|P01308|INS_HUMAN Insulin A chain",
    "GIVEQCCTSI",
    "CSLYQLENYCN"
  ))
  expect_equal(raw$sequence, sp_ins_a)
  expect_equal(fa$sequence, sp_ins_a)
  expect_equal(fa$header, "sp|P01308|INS_HUMAN Insulin A chain")
  expect_true(is.na(raw$header))
  expect_equal(raw$length, 21L)
  expect_true(raw$is_valid)
  expect_equal(raw$n_entries, 1L)
})

test_that("lower case, spaces, digits and line breaks are removed", {
  cs <- clean_protein_sequence("  1 giveqcctsi 11\r\ncslyqlenyc n\n")
  expect_equal(cs$sequence, sp_ins_a)
  expect_true(cs$is_valid)
})

test_that("only the first of several FASTA entries is kept", {
  cs <- clean_protein_sequence(">a\nGIVEQ\n>b\nFVNQH\n>c\nAAAA")
  expect_equal(cs$sequence, "GIVEQ")
  expect_equal(cs$header, "a")
  expect_equal(cs$n_entries, 3L)
})

test_that("non-canonical residues are reported, illegal characters invalidate", {
  cs <- clean_protein_sequence("GIVXEQB*")
  expect_equal(cs$invalid_chars, sort(c("X", "B", "*")))
  expect_true(cs$is_valid)

  bad <- clean_protein_sequence("GIV-EQ?")
  expect_false(bad$is_valid)
  expect_setequal(bad$illegal_chars, c("-", "?"))
})

test_that("empty input is invalid, not an error", {
  for (x in list(NULL, "", "   \n ", ">header only")) {
    cs <- clean_protein_sequence(x)
    expect_equal(cs$sequence, "")
    expect_equal(cs$length, 0L)
    expect_false(cs$is_valid)
  }
})

# ── seqprofile_windows() ────────────────────────────────────────────────────
test_that("windows are centred, with the documented even-window convention", {
  odd <- seqprofile_windows("ABCDEFG", 3L)
  expect_equal(odd$position, 2:6)
  expect_equal(odd$residue, c("B", "C", "D", "E", "F"))
  expect_equal(odd$window_seq, c("ABC", "BCD", "CDE", "DEF", "EFG"))

  even <- seqprofile_windows("ABCDEFG", 4L)
  expect_equal(even$position, 2:5)
  expect_equal(even$window_seq[1], "ABCD")
})

test_that("a sequence shorter than the window raises a classed, readable error", {
  expect_error(
    seqprofile_windows("GIVEQ", 11L),
    "shorter than the window",
    class = "seqprofile_short_sequence"
  )
})

# ── extract_profile_regions() ───────────────────────────────────────────────
test_that("contiguous positions crossing the threshold are grouped", {
  prof <- data.frame(
    position = 1:10,
    residue = LETTERS[1:10],
    score = c(0, 2, 3, 1, -1, NA, 4, 4.5, -2, 1)
  )
  reg <- extract_profile_regions(prof, threshold = 1)
  expect_equal(reg$region_id, 1:3)
  expect_equal(reg$start, c(2L, 7L, 10L))
  expect_equal(reg$end, c(4L, 8L, 10L))
  expect_equal(reg$length, c(3L, 2L, 1L))
  expect_equal(reg$sequence, c("BCD", "GH", "J"))
  expect_equal(reg$peak_pos, c(3L, 8L, 10L))
  expect_equal(reg$peak_score, c(3, 4.5, 1))
  expect_equal(reg$mean_score, c(2, 4.25, 1))
})

test_that("direction = 'below', NA gaps and empty results are handled", {
  prof <- data.frame(
    position = 1:6,
    residue = LETTERS[1:6],
    score = c(-1, -3, NA, -2, 5, 5)
  )
  reg <- extract_profile_regions(prof, threshold = -1, direction = "below")
  expect_equal(reg$start, c(1L, 4L))
  expect_equal(reg$end, c(2L, 4L))
  expect_equal(reg$peak_pos, c(2L, 4L))

  none <- extract_profile_regions(prof, threshold = 100)
  expect_equal(nrow(none), 0L)
  expect_named(
    none,
    c(
      "region_id",
      "start",
      "end",
      "length",
      "sequence",
      "peak_pos",
      "peak_score",
      "mean_score",
      "note"
    )
  )
  expect_equal(nrow(extract_profile_regions(prof, threshold = NA)), 0L)
})

# ── run_seq_profile() ───────────────────────────────────────────────────────
test_that("every registered algorithm returns the full data contract", {
  skip_if_not_installed("Peptides")
  for (id in names(seqprofile_registry())) {
    res <- run_seq_profile(sp_ins_b, id)
    expect_named(
      res,
      c(
        "algorithm_id",
        "algorithm_name",
        "sequence",
        "params",
        "profile",
        "regions",
        "global_stats",
        "interpretation",
        "y_label",
        "direction",
        "baseline",
        "masked_residues"
      )
    )
    expect_equal(res$algorithm_id, id)
    expect_named(
      res$profile,
      c("position", "residue", "window_seq", "score", "above_thr")
    )
    expect_type(res$profile$score, "double")
    expect_false(anyNA(res$profile$score))
    expect_equal(nrow(res$profile), nchar(sp_ins_b) - res$params$window + 1L)
    expect_true(all(
      c(
        "seq_length",
        "mean_score",
        "sd_score",
        "min_score",
        "max_score",
        "n_regions",
        "pct_above"
      ) %in%
        names(res$global_stats)
    ))
    expect_named(
      res$interpretation,
      c("summary", "reading", "thresholds", "caveats", "reference")
    )
    expect_equal(res$global_stats$n_regions, nrow(res$regions))
  }
})

test_that("Kyte-Doolittle windows agree with Peptides and the in-repo GRAVY()", {
  skip_if_not_installed("Peptides")
  res <- run_seq_profile(sp_ins_b, "hydropathy_kd", list(window = 9))
  w <- res$profile$window_seq
  expect_equal(
    res$profile$score,
    Peptides::hydrophobicity(w, scale = "KyteDoolittle")
  )
  # GRAVY() carries its own copy of the scale: an independent cross-check
  expect_equal(res$profile$score, unname(vapply(w, GRAVY, numeric(1))))
})

test_that("Hopp-Woods and Lehninger charge match direct Peptides calls", {
  skip_if_not_installed("Peptides")
  hw <- run_seq_profile(sp_ins_b, "hydrophilicity_hw")
  expect_equal(hw$params$window, 6L)
  expect_equal(
    hw$profile$score,
    Peptides::hydrophobicity(hw$profile$window_seq, scale = "HoppWoods")
  )

  ch <- run_seq_profile(sp_ins_b, "charge_lehninger", list(pH = 5))
  expect_equal(
    ch$profile$score,
    Peptides::charge(ch$profile$window_seq, pH = 5, pKscale = "Lehninger") / 9
  )
  # Charge density falls as pH rises
  ch_hi <- run_seq_profile(sp_ins_b, "charge_lehninger", list(pH = 10))
  expect_true(all(ch_hi$profile$score < ch$profile$score))
})

test_that("the hydrophobic moment uses the chosen angle and window", {
  skip_if_not_installed("Peptides")
  a <- run_seq_profile(sp_ins_b, "hmoment_eisenberg", list(angle = 100))
  b <- run_seq_profile(sp_ins_b, "hmoment_eisenberg", list(angle = 160))
  expect_equal(
    a$profile$score[1],
    Peptides::hmoment(substr(sp_ins_b, 1, 11), angle = 100, window = 11)
  )
  expect_false(isTRUE(all.equal(a$profile$score, b$profile$score)))
  expect_true(all(a$profile$score >= 0))
  # No preset cutoff: no regions until the user supplies one
  expect_equal(nrow(a$regions), 0L)
  expect_true(is.na(a$global_stats$pct_above))
  expect_gt(
    nrow(
      run_seq_profile(
        sp_ins_b,
        "hmoment_eisenberg",
        list(threshold = min(a$profile$score))
      )$regions
    ),
    0L
  )
})

test_that("windows with non-canonical residues are masked, not mis-scored", {
  skip_if_not_installed("Peptides")
  seq_x <- "GIVEQCCTSIXSLYQLENYCN*"
  expect_no_warning(
    res <- run_seq_profile(seq_x, "hydropathy_kd", list(window = 5))
  )
  has_bad <- grepl("[X*]", res$profile$window_seq)
  expect_true(any(has_bad))
  expect_true(all(is.na(res$profile$score[has_bad])))
  expect_false(anyNA(res$profile$score[!has_bad]))
  expect_equal(res$masked_residues, sort(c("X", "*")))
  expect_equal(res$global_stats$n_masked, sum(has_bad))
  # A masked window splits a region instead of bridging it
  expect_true(all(res$regions$end < 11 | res$regions$start > 11))
})

test_that("bad input gives classed errors rather than R indexing errors", {
  skip_if_not_installed("Peptides")
  expect_error(
    run_seq_profile("GIVEQ", "hmoment_eisenberg"),
    class = "seqprofile_short_sequence"
  )
  expect_error(run_seq_profile("", "hydropathy_kd"), class = "seqprofile_empty")
  expect_error(
    run_seq_profile(sp_ins_a, "no_such_algorithm"),
    class = "seqprofile_bad_param"
  )
  expect_error(
    run_seq_profile(sp_ins_a, "hydropathy_kd", list(window = 2.5)),
    class = "seqprofile_bad_param"
  )
  expect_error(
    run_seq_profile(sp_ins_a, "charge_lehninger", list(pH = 20)),
    class = "seqprofile_bad_param"
  )
})

test_that("global stats and descriptors are consistent with the profile", {
  skip_if_not_installed("Peptides")
  res <- run_seq_profile(sp_ins_b, "hydropathy_kd")
  gs <- res$global_stats
  expect_equal(gs$seq_length, 30L)
  expect_equal(gs$max_score, max(res$profile$score))
  expect_equal(gs$max_pos, res$profile$position[which.max(res$profile$score)])
  expect_equal(gs$pct_above, 100 * mean(res$profile$score >= 0))

  gd <- seqprofile_global_descriptors(sp_ins_b)
  expect_equal(gd$gravy, GRAVY(sp_ins_b))
  expect_equal(gd$mw, Peptides::mw(sp_ins_b))
  expect_true(all(is.na(unlist(seqprofile_global_descriptors("GIVXEQ")))))
})

test_that("the profile plot builds, with and without regions", {
  skip_if_not_installed("Peptides")
  skip_if_not_installed("ggplot2")
  p <- plot_seq_profile(run_seq_profile(sp_ins_b, "hydropathy_kd"))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  p2 <- plot_seq_profile(run_seq_profile(sp_ins_b, "hmoment_eisenberg"))
  expect_no_error(ggplot2::ggplot_build(p2))
  all_na <- plot_seq_profile(run_seq_profile("XXXXXXXXXX", "hydropathy_kd"))
  expect_s3_class(all_na, "ggplot")
})

test_that("seqprofile_choices() groups every algorithm by category", {
  ch <- seqprofile_choices()
  expect_identical(names(ch), SEQPROFILE_CATEGORIES)
  expect_setequal(unlist(ch, use.names = FALSE), names(seqprofile_registry()))
  ids <- lapply(ch, function(g) unlist(g, use.names = FALSE))
  expect_identical(
    ids[["Physicochemical"]],
    c(
      "charge_lehninger",
      "hydropathy_kd",
      "hydrophilicity_hw",
      "hydrophilicity_parker"
    )
  )
  expect_identical(
    ids[["Structural / Membrane"]],
    c("hmoment_eisenberg", "surface_emini", "betaturn_chou_fasman")
  )
  expect_identical(
    ids[["Immunogenicity / Epitopes"]],
    c("antigenicity_kolaskar", "immunogenicity_calis")
  )
})

# ── SeqProfile_server() ─────────────────────────────────────────────────────
test_that("the module computes only on Predict and validates empty input", {
  skip_if_not_installed("Peptides")
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
      algorithm = "hydropathy_kd",
      window = 9,
      threshold = 0,
      table_view = "regions"
    )
    # Nothing is computed before the first click
    expect_error(result(), class = "shiny.silent.error")

    session$setInputs(predict = 1)
    expect_error(result(), "Paste a sequence")

    session$setInputs(sequence = "GIVEQ", predict = 2)
    expect_error(result(), "shorter than the window")

    session$setInputs(sequence = sp_ins_b, predict = 3)
    expect_equal(result()$algorithm_id, "hydropathy_kd")
    expect_equal(names(table_data())[1], "Region_ID")
    session$setInputs(table_view = "profile")
    expect_equal(nrow(table_data()), 22L)

    # Changing inputs does not recompute until Predict is clicked again
    session$setInputs(algorithm = "hydrophilicity_hw", window = 6)
    expect_equal(result()$algorithm_id, "hydropathy_kd")
    session$setInputs(predict = 4)
    expect_equal(result()$algorithm_id, "hydrophilicity_hw")
    expect_equal(result()$params$window, 6L)
  })
})
