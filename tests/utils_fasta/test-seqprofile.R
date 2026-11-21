# Sequence profiling helpers in utils_fasta.r and the SeqProfile module server.

sp_ins_a <- "GIVEQCCTSICSLYQLENYCN" # human insulin A chain
sp_ins_b <- "FVNQHLCGSHLVEALYLVCGERGFFYTPKT" # human insulin B chain
# Human preproinsulin (UniProt P01308, 110 aa): long enough for the 51-residue
# FoldIndex window; residues 1-24 are the signal peptide.
sp_preproins <- paste0(
  "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGG",
  "PGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)

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
    res <- run_seq_profile(sp_preproins, id)
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
    expect_equal(
      nrow(res$profile),
      nchar(sp_preproins) - res$params$window + 1L
    )
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
    c(
      "flexibility_ks",
      "tm_von_heijne",
      "idr_foldindex",
      "hmoment_eisenberg",
      "surface_emini",
      "betaturn_chou_fasman"
    )
  )
  expect_identical(
    ids[["Immunogenicity / Epitopes"]],
    c("antigenicity_jw", "antigenicity_kolaskar", "immunogenicity_calis")
  )
})

# ── Karplus-Schulz backbone flexibility ─────────────────────────────────────
test_that("Karplus-Schulz weights are triangular and sum to 1", {
  expect_equal(
    ks_triangular_weights(7),
    c(0.25, 0.5, 0.75, 1, 0.75, 0.5, 0.25) / 4
  )
  for (w in 3:12) {
    wt <- ks_triangular_weights(w)
    expect_length(wt, w)
    expect_equal(sum(wt), 1)
    expect_equal(wt, rev(wt))
    expect_true(all(wt > 0))
  }
})

test_that("Karplus-Schulz scores are the weighted sum of the flexibility scale", {
  expect_named(SEQPROFILE_SCALE_KS_FLEX, SEQPROFILE_AA, ignore.order = TRUE)
  # A homopolymer scores its own scale value, whatever the window
  expect_equal(score_karplus_schulz("AAAAAAA"), 0.946)
  expect_equal(score_karplus_schulz("DDD"), 1.089)
  # G I V E Q C C with weights 1 2 3 4 3 2 1 (/ 16), worked by hand
  expect_equal(
    score_karplus_schulz("GIVEQCC"),
    (1.042 +
      2 * 0.892 +
      3 * 0.927 +
      4 * 1.036 +
      3 * 1.025 +
      2 * 0.878 +
      0.878) /
      16
  )
  expect_equal(score_karplus_schulz("giveqcc"), score_karplus_schulz("GIVEQCC"))
  expect_true(is.na(score_karplus_schulz("GIVXQCC")))
  expect_equal(score_karplus_schulz(character(0)), numeric(0))
})

test_that("the Karplus-Schulz profile shades flexible stretches above 1", {
  res <- run_seq_profile(sp_preproins, "flexibility_ks")
  expect_equal(res$params$window, 7L)
  expect_equal(res$params$threshold, 1)
  expect_equal(res$direction, "above")
  expect_equal(res$profile$score, score_karplus_schulz(res$profile$window_seq))
  expect_identical(res$profile$above_thr, res$profile$score >= 1)
  expect_gt(nrow(res$regions), 0L)
  expect_true(all(res$regions$peak_score >= 1))
  expect_setequal(res$regions$note, "Flexible loop / segment")
  # The scale spans 0.862 (Met) to 1.089 (Asp)
  expect_true(all(res$profile$score > 0.862 & res$profile$score < 1.089))

  expect_equal(
    run_seq_profile(sp_ins_a, "flexibility_ks", list(window = 3))$params$window,
    3L
  )
  expect_error(
    run_seq_profile(sp_ins_a, "flexibility_ks", list(window = 2)),
    "at least 3",
    class = "seqprofile_bad_param"
  )
})

# ── von Heijne TopPred transmembrane helices ────────────────────────────────
test_that("TopPred weights form a trapezoid of the window's length", {
  # 21 residues: 5-residue ramps around an 11-residue flat core
  expect_equal(
    toppred_trapezoid_weights(21L),
    c(1:5, rep(6, 11), 5:1) / 96
  )
  for (w in 11:25) {
    wt <- toppred_trapezoid_weights(w)
    expect_length(wt, w)
    expect_equal(sum(wt), 1)
    expect_equal(wt, rev(wt))
  }
})

test_that("TopPred scores are positive for hydrophobic windows", {
  skip_if_not_installed("Peptides")
  # Homopolymers score the GES value of the residue (Leu 2.8, Ala 1.6, Lys -8.8)
  expect_equal(
    score_von_heijne_tm(c(strrep("L", 21), strrep("A", 21), strrep("K", 21))),
    c(2.8, 1.6, -8.8)
  )
  expect_true(is.na(score_von_heijne_tm(paste0(strrep("L", 20), "X"))))
  expect_equal(score_von_heijne_tm(character(0)), numeric(0))

  # Mixed window against the weights and the scale, both taken independently
  win <- substr(sp_preproins, 1, 21)
  aa <- strsplit(win, "")[[1]]
  ges <- -vapply(
    aa,
    function(a) Peptides::hydrophobicity(a, scale = "Engelman"),
    numeric(1)
  )
  expect_equal(
    score_von_heijne_tm(win),
    sum(ges * c(1:5, rep(6, 11), 5:1)) / 96
  )
})

test_that("TopPred regions are annotated with certainty and orientation", {
  # One region centred on residue 26: helix 16-36, flanks 1-15 and 37-51
  region <- function(peak, start = 26L, end = 26L) {
    data.frame(start = start, end = end, peak_score = peak)
  }
  p <- list(window = 21L)
  n_flank <- "AAKKRAAAAAAAAAA" # K+R = 3
  c_flank <- "AAAAAAAAAAAAAAR" # K+R = 1
  core <- strrep("L", 21)

  expect_equal(
    annotate_von_heijne_regions(
      region(2.8),
      paste0(n_flank, core, c_flank),
      p
    ),
    "Certain TM | N-in/C-out (N-flank K+R=3, C-flank K+R=1)"
  )
  expect_equal(
    annotate_von_heijne_regions(
      region(0.7),
      paste0(c_flank, core, n_flank),
      p
    ),
    "Putative TM | N-out/C-in (N-flank K+R=1, C-flank K+R=3)"
  )
  # The certain cutoff of 1.0 is inclusive
  expect_match(
    annotate_von_heijne_regions(
      region(1),
      paste0(c_flank, core, c_flank),
      p
    ),
    "^Certain TM \\| Balanced flanks \\(N-flank K\\+R=1, C-flank K\\+R=1\\)$"
  )
  # A helix at the very start or end has an empty flank, not an error
  expect_match(
    annotate_von_heijne_regions(
      region(2, start = 11L, end = 11L),
      paste0(core, n_flank),
      p
    ),
    "N-out/C-in (N-flank K+R=0, C-flank K+R=3)",
    fixed = TRUE
  )
  # One note per region
  expect_length(
    annotate_von_heijne_regions(
      region(c(2, 0.7), start = c(20L, 30L), end = c(22L, 31L)),
      paste0(n_flank, core, c_flank),
      p
    ),
    2L
  )
  expect_equal(
    annotate_von_heijne_regions(region(1)[0, ], core, p),
    character(0)
  )
})

test_that("the TopPred profile finds a hydrophobic helix and its signal-peptide look-alike", {
  skip_if_not_installed("Peptides")
  tm_seq <- paste0("AAKKRAAAAAAAAAA", strrep("L", 21), "AAAAAAAAAAAAAAR")
  res <- run_seq_profile(tm_seq, "tm_von_heijne")
  expect_equal(res$params$window, 21L)
  expect_equal(res$params$threshold, 0.6)
  expect_equal(res$baseline, 1)
  expect_equal(nrow(res$regions), 1L)
  # The window made only of Leu is centred on residue 26
  expect_equal(res$regions$peak_pos, 26L)
  expect_equal(res$regions$peak_score, 2.8)
  expect_match(res$regions$note, "^Certain TM \\| ")
  expect_equal(res$global_stats$max_pos, 26L)

  # Documented caveat: the preproinsulin signal peptide (1-24) passes as a helix
  sig <- run_seq_profile(sp_preproins, "tm_von_heijne")
  expect_equal(nrow(sig$regions), 1L)
  expect_true(sig$regions$start >= 1L && sig$regions$end <= 24L)

  # A threshold nothing reaches leaves no region to annotate
  none <- run_seq_profile(tm_seq, "tm_von_heijne", list(threshold = 99))
  expect_equal(nrow(none$regions), 0L)
  expect_type(none$regions$note, "character")

  expect_error(
    run_seq_profile(tm_seq, "tm_von_heijne", list(window = 9)),
    "at least 11",
    class = "seqprofile_bad_param"
  )
  expect_error(
    run_seq_profile(sp_ins_a, "tm_von_heijne", list(window = 25)),
    class = "seqprofile_short_sequence"
  )
})

# ── Jameson-Wolf antigenic index ────────────────────────────────────────────
test_that("the Jameson-Wolf index reaches its documented extremes", {
  skip_if_not_installed("Peptides")
  # Asp: hydrophilic (2), surface (1), flexible (1), turn (2)
  #   0.30 * 2 + 0.15 + 0.15 + 0.40 * 2 = 1.7
  # Ile: hydrophobic (-2), buried (0), rigid (0), no turn (0)
  #   0.30 * -2 = -0.6
  expect_equal(
    score_jameson_wolf(c("DDDDDDD", "IIIIIII")),
    c(1.7, -0.6)
  )
  expect_equal(score_jameson_wolf(character(0)), numeric(0))
})

test_that("the Jameson-Wolf index combines its four discretised components", {
  skip_if_not_installed("Peptides")
  res <- run_seq_profile(sp_preproins, "antigenicity_jw")
  expect_equal(res$params$window, 7L)
  expect_equal(res$params$threshold, 0)
  w <- res$profile$window_seq

  hw <- Peptides::hydrophobicity(w, scale = "HoppWoods")
  h <- c(-2, -1, 1, 2)[findInterval(hw, c(-0.4, 0, 0.5), left.open = TRUE) + 1L]
  # A 7-residue window scores Emini on its first six residues
  s <- as.numeric(score_emini_surface(substr(w, 1, 6)) >= 1)
  f <- as.numeric(score_karplus_schulz(w) >= 1)
  ct <- score_chou_fasman_turn(w)
  t <- (ct >= 1) + (ct >= 1.2)
  expect_equal(res$profile$score, 0.30 * h + 0.15 * s + 0.15 * f + 0.40 * t)

  expect_true(all(res$profile$score >= -0.6 & res$profile$score <= 1.7))
  expect_identical(res$profile$above_thr, res$profile$score >= 0)
  expect_setequal(res$regions$note, "Predicted antigenic peak")

  # Emini is read from the six central residues of a longer window
  w8 <- run_seq_profile(sp_preproins, "antigenicity_jw", list(window = 8))
  win <- w8$profile$window_seq
  hw8 <- Peptides::hydrophobicity(win, scale = "HoppWoods")
  h8 <- c(-2, -1, 1, 2)[
    findInterval(hw8, c(-0.4, 0, 0.5), left.open = TRUE) + 1L
  ]
  ct8 <- score_chou_fasman_turn(win)
  expect_equal(
    w8$profile$score,
    0.30 *
      h8 +
      0.15 * (score_emini_surface(substr(win, 2, 7)) >= 1) +
      0.15 * (score_karplus_schulz(win) >= 1) +
      0.40 * ((ct8 >= 1) + (ct8 >= 1.2))
  )

  expect_error(
    run_seq_profile(sp_preproins, "antigenicity_jw", list(window = 5)),
    "at least 6",
    class = "seqprofile_bad_param"
  )
})

# ── FoldIndex intrinsic disorder ────────────────────────────────────────────
test_that("FoldIndex follows the Uversky charge-hydropathy formula", {
  skip_if_not_installed("Peptides")
  # Hydrophobic and uncharged folds; charged and hydrophilic does not,
  # whatever the sign of the charge
  fi <- score_foldindex(c(strrep("I", 15), strrep("E", 15), strrep("K", 15)))
  expect_gt(fi[1], 0)
  expect_lt(fi[2], 0)
  expect_lt(fi[3], 0)
  expect_equal(score_foldindex(character(0)), numeric(0))

  res <- run_seq_profile(sp_preproins, "idr_foldindex", list(window = 15))
  w <- res$profile$window_seq
  # GRAVY() carries its own copy of the Kyte-Doolittle scale
  kd <- unname(vapply(w, GRAVY, numeric(1)))
  charge <- Peptides::charge(w, pH = 7.4, pKscale = "Lehninger") / 15
  expect_equal(
    res$profile$score,
    2.785 * (kd + 4.5) / 9 - abs(charge) - 1.151
  )

  # The charge term depends on pH
  acid <- run_seq_profile(
    sp_preproins,
    "idr_foldindex",
    list(window = 15, pH = 3)
  )
  expect_equal(acid$params$pH, 3)
  expect_false(isTRUE(all.equal(acid$profile$score, res$profile$score)))
  expect_equal(
    acid$profile$score,
    score_foldindex(acid$profile$window_seq, pH = 3)
  )
})

test_that("FoldIndex regions are the windows at or below the threshold", {
  skip_if_not_installed("Peptides")
  res <- run_seq_profile(sp_preproins, "idr_foldindex", list(window = 15))
  expect_equal(res$direction, "below")
  score <- res$profile$score
  expect_identical(res$profile$above_thr, score <= 0)
  expect_equal(res$global_stats$pct_above, 100 * mean(score <= 0))

  expect_gt(nrow(res$regions), 0L)
  expect_lt(nrow(res$regions), nrow(res$profile))
  covered <- unlist(Map(seq, res$regions$start, res$regions$end))
  expect_setequal(covered, res$profile$position[score <= 0])
  # The peak of a disordered region is its lowest score
  for (i in seq_len(nrow(res$regions))) {
    in_region <- res$profile$position %in%
      res$regions$start[i]:res$regions$end[i]
    expect_equal(res$regions$peak_score[i], min(score[in_region]))
  }
  expect_setequal(
    res$regions$note,
    "Predicted intrinsically disordered region (IDR)"
  )

  # At the highest score every window is at or below the threshold
  all_in <- run_seq_profile(
    sp_preproins,
    "idr_foldindex",
    list(window = 15, threshold = max(score))
  )
  expect_equal(nrow(all_in$regions), 1L)
  expect_equal(all_in$regions$length, nrow(all_in$profile))
  expect_equal(all_in$global_stats$pct_above, 100)
})

test_that("FoldIndex defaults to the 51-residue window of Prilusky et al.", {
  skip_if_not_installed("Peptides")
  res <- run_seq_profile(sp_preproins, "idr_foldindex")
  expect_equal(
    res$params[c("window", "pH", "threshold")],
    list(
      window = 51L,
      pH = 7.4,
      threshold = 0
    )
  )
  expect_equal(nrow(res$profile), 60L)
  expect_equal(res$profile$position[1], 26L)

  # Shorter than the default window: a readable error, not an indexing one
  expect_error(
    run_seq_profile(sp_ins_b, "idr_foldindex"),
    "shorter than the window of 51",
    class = "seqprofile_short_sequence"
  )
  expect_error(
    run_seq_profile(sp_preproins, "idr_foldindex", list(window = 10)),
    "at least 11",
    class = "seqprofile_bad_param"
  )
  expect_error(
    run_seq_profile(sp_preproins, "idr_foldindex", list(pH = -1)),
    class = "seqprofile_bad_param"
  )
})

# ── The four profiles above in the plot, the consensus and the module ───────
test_that("the newer profiles mask non-canonical residues and plot", {
  skip_if_not_installed("Peptides")
  seq_x <- paste0(
    substr(sp_preproins, 1, 54),
    "X",
    substr(sp_preproins, 56, 110)
  )
  windows <- c(
    flexibility_ks = 7,
    tm_von_heijne = 21,
    antigenicity_jw = 7,
    idr_foldindex = 15
  )
  for (id in names(windows)) {
    expect_no_warning(
      res <- run_seq_profile(seq_x, id, list(window = windows[[id]]))
    )
    has_x <- grepl("X", res$profile$window_seq, fixed = TRUE)
    expect_equal(sum(has_x), windows[[id]])
    expect_true(all(is.na(res$profile$score[has_x])))
    expect_false(anyNA(res$profile$score[!has_x]))
    # A masked window splits a region instead of bridging it
    expect_true(all(res$regions$end < 55 | res$regions$start > 55))

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- plot_seq_profile(
        run_seq_profile(sp_preproins, id, list(window = windows[[id]]))
      )
      expect_s3_class(p, "ggplot")
      expect_no_error(ggplot2::ggplot_build(p))
    }
  }
})

test_that("a consensus counts FoldIndex windows below their threshold", {
  skip_if_not_installed("Peptides")
  fi <- run_seq_profile(sp_preproins, "idr_foldindex", list(window = 15))
  ks <- run_seq_profile(sp_preproins, "flexibility_ks")
  cons <- seqprofile_consensus_regions(list(fi, ks))
  expected <- intersect(
    fi$profile$position[fi$profile$score <= 0],
    ks$profile$position[ks$profile$score >= 1]
  )
  expect_gt(length(expected), 0)
  expect_setequal(unlist(Map(seq, cons$start, cons$end)), expected)
  expect_true(all(cons$mean_idr_foldindex <= 0))
  expect_true(all(cons$mean_flexibility_ks >= 1))
})

test_that("the module labels FoldIndex windows as below the threshold", {
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
      sequence = sp_preproins,
      algorithm = c("idr_foldindex", "tm_von_heijne"),
      window_idr_foldindex = 15,
      table_view = "regions",
      predict = 1
    )
    expect_named(results(), c("idr_foldindex", "tm_von_heijne"))
    expect_equal(result()$direction, "below")
    expect_match(
      output$summary_strip$html,
      "Windows below threshold",
      fixed = TRUE
    )
    expect_setequal(
      table_data()$Note,
      "Predicted intrinsically disordered region (IDR)"
    )
    expect_true(all(table_data()$Peak_Score <= 0))

    # The TopPred note (certainty and orientation) reaches the regions table
    session$setInputs(summary_algorithm = "tm_von_heijne")
    expect_match(
      output$summary_strip$html,
      "Windows above threshold",
      fixed = TRUE
    )
    expect_match(table_data()$Note, "^(Certain|Putative) TM \\| ")

    # The default FoldIndex window is longer than a short peptide
    session$setInputs(
      sequence = sp_ins_b,
      algorithm = "idr_foldindex",
      window_idr_foldindex = 51,
      predict = 2
    )
    expect_error(result(), "shorter than the window of 51")
  })
})

# ── In-file scales: Kolaskar-Tongaonkar, Parker, Chou-Fasman, Emini, Calis ──
# Mean of a scale over each window, written without the helpers under test
sp_window_mean <- function(window_seq, scale) {
  vapply(
    strsplit(window_seq, ""),
    function(aa) sum(scale[aa]) / length(aa),
    numeric(1),
    USE.NAMES = FALSE
  )
}

test_that("every in-file scale has one value per canonical residue", {
  scales <- list(
    SEQPROFILE_SCALE_KOLASKAR,
    SEQPROFILE_SCALE_EMINI,
    SEQPROFILE_SCALE_PARKER,
    SEQPROFILE_SCALE_CHOU_FASMAN_TURN,
    SEQPROFILE_SCALE_CALIS
  )
  for (scale in scales) {
    expect_length(scale, 20L)
    expect_setequal(names(scale), SEQPROFILE_AA)
    expect_false(anyNA(scale))
  }
  # Extremes quoted in the reading guides
  expect_equal(SEQPROFILE_SCALE_PARKER[["W"]], -10)
  expect_equal(SEQPROFILE_SCALE_PARKER[["D"]], 10)
  expect_equal(names(which.max(SEQPROFILE_SCALE_KOLASKAR)), "C")
  expect_equal(names(which.max(SEQPROFILE_SCALE_CALIS)), "W")
  # Emini values are fractions; the anchors P1, P2 and P9 of Calis weigh nothing
  expect_true(all(SEQPROFILE_SCALE_EMINI > 0 & SEQPROFILE_SCALE_EMINI < 1))
  expect_length(SEQPROFILE_WEIGHTS_CALIS_9MER, 9L)
  expect_equal(SEQPROFILE_WEIGHTS_CALIS_9MER[c(1, 2, 9)], c(0, 0, 0))
})

test_that("seqprofile_scale_mean() averages a scale over each window", {
  scale <- c(A = 1, C = 3, D = -2)
  expect_equal(
    seqprofile_scale_mean(c("AC", "ACD", "DDDD"), scale),
    c(2, 2 / 3, -2)
  )
  expect_equal(seqprofile_scale_mean("acd", scale), 2 / 3)
  # A residue missing from the scale gives NA, not a partial mean
  expect_true(is.na(seqprofile_scale_mean("ACX", scale)))
  expect_true(is.na(seqprofile_scale_mean("", scale)))
  expect_equal(seqprofile_scale_mean(character(0), scale), numeric(0))
})

test_that("Kolaskar, Parker and Chou-Fasman scores are window means of their scales", {
  # G I V E Q C C, summed by hand from each scale
  expect_equal(
    score_kolaskar_tongaonkar("GIVEQCC"),
    (0.874 + 1.152 + 1.383 + 0.851 + 1.015 + 1.412 + 1.412) / 7
  )
  expect_equal(
    score_parker_hydrophilicity("GIVEQCC"),
    (5.7 - 8.0 - 3.7 + 7.8 + 6.0 + 1.4 + 1.4) / 7
  )
  expect_equal(
    score_chou_fasman_turn("GIVEQCC"),
    (1.56 + 0.47 + 0.50 + 0.74 + 0.98 + 1.19 + 1.19) / 7
  )
  # A homopolymer scores the scale value of its residue
  expect_equal(score_kolaskar_tongaonkar("VVVVVVV"), 1.383)
  expect_equal(score_parker_hydrophilicity("DDDDDDD"), 10)
  expect_equal(score_chou_fasman_turn("NNNNNNN"), 1.56)

  fns <- list(
    score_kolaskar_tongaonkar,
    score_parker_hydrophilicity,
    score_chou_fasman_turn
  )
  for (fn in fns) {
    expect_true(is.na(fn("GIVXQCC")))
    expect_equal(fn(character(0)), numeric(0))
    expect_length(fn(c("GIVEQCC", "IVEQCCT")), 2L)
  }

  # Documented caveat of Kolaskar-Tongaonkar: a hydrophobic core scores above
  # 1 although Parker calls the same stretch hydrophobic
  expect_gt(score_kolaskar_tongaonkar("VLCVLCV"), 1)
  expect_lt(score_parker_hydrophilicity("VLCVLCV"), 0)
  # Turn formers (Gly, Asn, Pro, Asp, Ser) against helix / strand formers
  expect_gt(score_chou_fasman_turn("GNPDSGN"), 1)
  expect_lt(score_chou_fasman_turn("VILMFAE"), 1)
})

test_that("the Emini score is the normalised product of surface fractions", {
  # Sn = prod(d) * 0.37^(-w): Lys 0.97, Cys 0.26
  expect_equal(score_emini_surface("KKKKKK"), (0.97 / 0.37)^6)
  expect_equal(score_emini_surface("CCCCCC"), (0.26 / 0.37)^6)
  expect_equal(
    score_emini_surface("GIVEQC"),
    0.48 * 0.34 * 0.36 * 0.84 * 0.84 * 0.26 / 0.37^6
  )
  # The exponent follows the window length
  expect_equal(score_emini_surface("KKKK"), (0.97 / 0.37)^4)
  # One buried residue lowers the score of an otherwise exposed window
  expect_equal(
    score_emini_surface("KKKCKK") / score_emini_surface("KKKKKK"),
    0.26 / 0.97
  )
  expect_gt(score_emini_surface("KKKKKK"), 1)
  expect_lt(score_emini_surface("CCCCCC"), 1)
  expect_equal(score_emini_surface("kkkkkk"), score_emini_surface("KKKKKK"))
  expect_true(is.na(score_emini_surface("KKKXKK")))
  expect_equal(score_emini_surface(character(0)), numeric(0))
})

test_that("the four IEDB B-cell profiles use their documented defaults", {
  specs <- list(
    antigenicity_kolaskar = list(
      window = 7L,
      threshold = 1,
      scale = SEQPROFILE_SCALE_KOLASKAR
    ),
    hydrophilicity_parker = list(
      window = 7L,
      threshold = 0,
      scale = SEQPROFILE_SCALE_PARKER
    ),
    betaturn_chou_fasman = list(
      window = 7L,
      threshold = 1,
      scale = SEQPROFILE_SCALE_CHOU_FASMAN_TURN
    )
  )
  for (id in names(specs)) {
    spec <- specs[[id]]
    res <- run_seq_profile(sp_preproins, id)
    expect_equal(res$params$window, spec$window)
    expect_equal(res$params$threshold, spec$threshold)
    expect_equal(res$direction, "above")
    # A 7-residue window is reported at its central (4th) residue
    expect_equal(res$profile$position[1], 4L)
    expect_equal(
      res$profile$score,
      sp_window_mean(res$profile$window_seq, spec$scale)
    )
    expect_identical(res$profile$above_thr, res$profile$score >= spec$threshold)
    expect_gt(nrow(res$regions), 0L)
    expect_lt(nrow(res$regions), nrow(res$profile))
  }

  em <- run_seq_profile(sp_preproins, "surface_emini")
  expect_equal(em$params$window, 6L)
  expect_equal(em$params$threshold, 1)
  # Even window: the score sits at the residue left of the midpoint
  expect_equal(em$profile$position[1], 3L)
  expect_equal(em$profile$window_seq[1], substr(sp_preproins, 1, 6))
  expect_equal(
    em$profile$score,
    vapply(
      strsplit(em$profile$window_seq, ""),
      function(aa) prod(SEQPROFILE_SCALE_EMINI[aa]) / 0.37^6,
      numeric(1)
    )
  )
  expect_true(all(em$profile$score > 0))
  expect_identical(em$profile$above_thr, em$profile$score >= 1)

  # The window size is a free parameter of all four
  wide <- run_seq_profile(
    sp_preproins,
    "antigenicity_kolaskar",
    list(window = 11)
  )
  expect_equal(
    wide$profile$score,
    sp_window_mean(wide$profile$window_seq, SEQPROFILE_SCALE_KOLASKAR)
  )
  expect_equal(nrow(wide$profile), 100L)
})

# ── Calis MHC-I immunogenicity ──────────────────────────────────────────────
test_that("parse_calis_mask() reads 1-based positions and the C-terminus", {
  expect_equal(parse_calis_mask("1,2,C", 9L), c(1L, 2L, 9L))
  expect_equal(parse_calis_mask("1, 2;c", 9L), c(1L, 2L, 9L))
  expect_equal(parse_calis_mask("1,1,cterm", 10L), c(1L, 10L))
  expect_equal(parse_calis_mask("C-term", 11L), 11L)
  # Numeric input; positions outside the peptide are dropped
  expect_equal(parse_calis_mask(c(1, 12, NA), 9L), 1L)
  expect_equal(parse_calis_mask("0,3,10", 9L), 3L)
  for (empty in list(NULL, character(0), "", "  ")) {
    expect_equal(parse_calis_mask(empty, 9L), integer(0))
  }
  expect_error(
    parse_calis_mask("1,z", 9L),
    "Invalid mask token 'Z'",
    class = "seqprofile_bad_param"
  )
  expect_error(parse_calis_mask("1.5", 9L), class = "seqprofile_bad_param")
})

test_that("the Calis score is the position-weighted sum of the scale", {
  # AAAAWAAAA: Trp (0.719) at P5 (0.30), Ala (0.127) at P3, P4, P6, P7, P8
  expect_equal(
    score_calis_immunogenicity("AAAAWAAAA"),
    0.30 * 0.719 + (0.10 + 0.31 + 0.29 + 0.26 + 0.18) * 0.127
  )
  # GILGFVFTL (influenza M1 58-66): L3 G4 F5 V6 F7 T8 carry the weights
  gil <- 0.10 *
    -0.036 +
    0.31 * 0.110 +
    0.30 * 0.380 +
    0.29 * 0.134 +
    0.26 * 0.380 +
    0.18 * 0.126
  expect_equal(score_calis_immunogenicity("GILGFVFTL"), gil)
  expect_equal(score_calis_immunogenicity("gilgfvftl"), gil)

  # P1, P2 and P9 already weigh nothing, so the default mask changes no 9-mer
  expect_equal(score_calis_immunogenicity("GILGFVFTL", mask = ""), gil)
  expect_equal(score_calis_immunogenicity("GILGFVFTL", mask = NULL), gil)
  # Masking a weighted position removes its term
  expect_equal(
    score_calis_immunogenicity("GILGFVFTL", mask = "4,5"),
    gil - 0.31 * 0.110 - 0.30 * 0.380
  )
  expect_equal(score_calis_immunogenicity("GILGFVFTL", mask = "3,4,5,6,7,8"), 0)

  # Anchor residues do not count: only P3-P8 differ between these peptides
  expect_equal(
    score_calis_immunogenicity("WWAAAAAAW"),
    score_calis_immunogenicity("KKAAAAAAK")
  )

  expect_true(is.na(score_calis_immunogenicity("GILGXVFTL")))
  expect_equal(score_calis_immunogenicity(character(0)), numeric(0))
  expect_error(
    score_calis_immunogenicity("GILGFVFT"),
    "at least 9",
    class = "seqprofile_bad_param"
  )
  expect_error(
    score_calis_immunogenicity("GILGFVFTL", mask = "1,x"),
    class = "seqprofile_bad_param"
  )
})

test_that("longer Calis peptides repeat the central weight", {
  # 10-mer: weights 0 0 .10 .31 .30 .30 .29 .26 .18 0
  expect_equal(
    score_calis_immunogenicity("AAAAWAAAAA"),
    0.30 * 0.719 + (0.10 + 0.31 + 0.30 + 0.29 + 0.26 + 0.18) * 0.127
  )
  # 11-mer: two extra central positions, here P6 (Trp) and P7 (Ala)
  expect_equal(
    score_calis_immunogenicity("AAAAAWAAAAA"),
    0.30 * 0.719 + (0.10 + 0.31 + 0.30 + 0.30 + 0.29 + 0.26 + 0.18) * 0.127
  )
  # "C" follows the peptide length, so the last residue never counts
  expect_equal(
    score_calis_immunogenicity("AAAAWAAAAW"),
    score_calis_immunogenicity("AAAAWAAAAA")
  )
  # Without a mask the last residue still weighs nothing (P9 weight is 0)
  expect_equal(
    score_calis_immunogenicity("AAAAWAAAAW", mask = ""),
    score_calis_immunogenicity("AAAAWAAAAA", mask = "")
  )
})

test_that("the Calis profile slides a 9-mer and passes the mask on", {
  res <- run_seq_profile(sp_preproins, "immunogenicity_calis")
  expect_equal(
    res$params[c("window", "mask", "threshold")],
    list(
      window = 9L,
      mask = "1,2,C",
      threshold = 0
    )
  )
  expect_equal(nrow(res$profile), 102L)
  expect_equal(res$profile$position[1], 5L)
  w9 <- SEQPROFILE_WEIGHTS_CALIS_9MER
  expect_equal(
    res$profile$score,
    vapply(
      strsplit(res$profile$window_seq, ""),
      function(aa) sum(SEQPROFILE_SCALE_CALIS[aa] * w9),
      numeric(1)
    )
  )
  expect_identical(res$profile$above_thr, res$profile$score >= 0)
  expect_setequal(res$regions$note, "Predicted immunogenic MHC-I peptide")

  masked <- run_seq_profile(
    sp_preproins,
    "immunogenicity_calis",
    list(mask = "1,2,5,C")
  )
  expect_equal(masked$params$mask, "1,2,5,C")
  expect_equal(
    res$profile$score - masked$profile$score,
    0.30 * unname(SEQPROFILE_SCALE_CALIS[substr(res$profile$window_seq, 5, 5)])
  )

  ten <- run_seq_profile(
    sp_preproins,
    "immunogenicity_calis",
    list(window = 10)
  )
  expect_equal(
    ten$profile$score,
    score_calis_immunogenicity(ten$profile$window_seq)
  )

  expect_error(
    run_seq_profile(sp_preproins, "immunogenicity_calis", list(window = 8)),
    "at least 9",
    class = "seqprofile_bad_param"
  )
  expect_error(
    run_seq_profile(sp_preproins, "immunogenicity_calis", list(mask = "1,z")),
    "Invalid mask token",
    class = "seqprofile_bad_param"
  )
  expect_error(
    run_seq_profile("GIVEQCCT", "immunogenicity_calis"),
    class = "seqprofile_short_sequence"
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
      window_hydropathy_kd = 9,
      threshold_hydropathy_kd = 0,
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
    session$setInputs(
      algorithm = "hydrophilicity_hw",
      window_hydrophilicity_hw = 6
    )
    expect_equal(result()$algorithm_id, "hydropathy_kd")
    session$setInputs(predict = 4)
    expect_equal(result()$algorithm_id, "hydrophilicity_hw")
    expect_equal(result()$params$window, 6L)
  })
})

# ── UniProt retrieval (no network: only parsing and validation) ─────────────
test_that("accessions are validated before any request is made", {
  expect_equal(seqprofile_clean_accession(" q07730 "), "Q07730")
  expect_equal(seqprofile_clean_accession("A0A024R161"), "A0A024R161")
  expect_equal(seqprofile_clean_accession("P05067-2"), "P05067-2")
  for (bad in list("", NULL, "insulin", "Q07730/../x", "Q0773")) {
    expect_error(
      seqprofile_clean_accession(bad),
      class = "seqprofile_bad_accession"
    )
  }
  expect_error(fetch_uniprot_fasta("?"), class = "seqprofile_bad_accession")
})

test_that("parse_uniprot_fasta() extracts the protein name and sequence", {
  txt <- ">sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1\nMALWMRLLPL\r\nLALLALWGPD\n"
  up <- parse_uniprot_fasta(txt)
  expect_equal(up$protein_name, "Insulin")
  expect_equal(up$sequence, "MALWMRLLPLLALLALWGPD")
  expect_false(grepl("\r", up$fasta))
  expect_equal(clean_protein_sequence(up$fasta)$sequence, up$sequence)
  expect_error(parse_uniprot_fasta(""), class = "seqprofile_fetch_failed")
  expect_error(
    parse_uniprot_fasta("<html>error</html>"),
    class = "seqprofile_fetch_failed"
  )
})

# ── seqprofile_consensus_regions() ──────────────────────────────────────────
test_that("consensus regions are the positions passing every threshold", {
  skip_if_not_installed("Peptides")
  kd <- run_seq_profile(sp_ins_b, "hydropathy_kd", list(window = 5))
  kt <- run_seq_profile(sp_ins_b, "antigenicity_kolaskar", list(window = 7))
  cons <- seqprofile_consensus_regions(list(kd, kt))
  expect_named(
    cons,
    c(
      "region_id",
      "start",
      "end",
      "length",
      "sequence",
      "mean_hydropathy_kd",
      "mean_antigenicity_kolaskar"
    )
  )

  pos_kd <- kd$profile$position[kd$profile$above_thr]
  pos_kt <- kt$profile$position[kt$profile$above_thr]
  expected <- intersect(pos_kd, pos_kt)
  expect_gt(length(expected), 0)
  covered <- unlist(Map(seq, cons$start, cons$end))
  expect_setequal(covered, expected)
  expect_equal(cons$length, cons$end - cons$start + 1L)
  expect_equal(cons$sequence, substring(sp_ins_b, cons$start, cons$end))
  first <- kd$profile$position %in% cons$start[1]:cons$end[1]
  expect_equal(cons$mean_hydropathy_kd[1], mean(kd$profile$score[first]))

  # A minimum length drops the shorter stretches and renumbers the rest
  cut <- max(cons$length)
  long <- seqprofile_consensus_regions(list(kd, kt), min_length = cut)
  expect_equal(long$start, cons$start[cons$length >= cut])
  expect_equal(long$region_id, seq_len(nrow(long)))
  expect_equal(
    nrow(seqprofile_consensus_regions(list(kd, kt), min_length = cut + 1L)),
    0L
  )
  for (bad in list(0, 2.5, NA, c(1, 2))) {
    expect_error(
      seqprofile_consensus_regions(list(kd, kt), min_length = bad),
      "minimum consensus length"
    )
  }

  # Consensus intervals are drawn as a second set of bands on each track
  skip_if_not_installed("ggplot2")
  n_rect <- function(p) {
    sum(vapply(
      p$layers,
      function(l) inherits(l$geom, "GeomRect"),
      logical(1)
    ))
  }
  expect_equal(n_rect(plot_seq_profile(kd)), 1L)
  shaded <- plot_seq_profile(kd, consensus = cons)
  expect_equal(n_rect(shaded), 2L)
  expect_no_error(ggplot2::ggplot_build(shaded))
  expect_equal(
    n_rect(plot_seq_profile(kd, consensus = cons[0, ])),
    1L
  )

  # Thresholds no window reaches give an empty table with the same columns
  none <- seqprofile_consensus_regions(list(
    kd,
    run_seq_profile(sp_ins_b, "antigenicity_kolaskar", list(threshold = 99))
  ))
  expect_equal(nrow(none), 0L)
  expect_named(none, names(cons))

  expect_error(
    seqprofile_consensus_regions(list(kd)),
    class = "seqprofile_bad_param"
  )
  expect_error(
    seqprofile_consensus_regions(list(
      kd,
      run_seq_profile(sp_ins_b, "hmoment_eisenberg")
    )),
    "Eisenberg"
  )
  expect_error(
    seqprofile_consensus_regions(list(
      kd,
      run_seq_profile(sp_ins_a, "antigenicity_kolaskar")
    )),
    "same sequence"
  )
})

# ── Helical wheel ───────────────────────────────────────────────────────────
test_that("the helical wheel reproduces the Peptides hydrophobic moment", {
  skip_if_not_installed("Peptides")
  pep <- "FLGKIWKLLAK"
  wheel <- seqprofile_helical_wheel(pep, start = 5L)
  expect_equal(wheel$position, 5:15)
  expect_equal(wheel$residue, strsplit(pep, "")[[1]])
  expect_equal(wheel$angle_deg[1:5], c(0, 100, 200, 300, 40))
  # Residue 1 at the top, residue 2 at 100 degrees clockwise
  expect_equal(c(wheel$x[1], wheel$y[1]), c(0, 1))
  expect_gt(wheel$x[2], 0)
  expect_lt(wheel$y[2], 0)
  expect_equal(
    attr(wheel, "moment")$uH,
    Peptides::hmoment(pep, angle = 100, window = 11)
  )

  # Residue 19 shares the angle of residue 1 and moves to an outer ring
  long <- seqprofile_helical_wheel(strrep("AL", 10))
  expect_equal(long$angle_deg[19], 0)
  expect_gt(long$y[19], long$y[1])

  expect_error(
    seqprofile_helical_wheel("AX"),
    class = "seqprofile_bad_param"
  )
  expect_error(
    seqprofile_helical_wheel("AAXAA"),
    class = "seqprofile_bad_param"
  )

  skip_if_not_installed("ggplot2")
  p <- plot_helical_wheel(wheel, title = "t", subtitle = "s")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("a wheel segment spans the windows behind a region", {
  skip_if_not_installed("Peptides")
  res <- run_seq_profile(
    sp_ins_b,
    "hmoment_eisenberg",
    list(threshold = 0)
  )
  # Every window passes a threshold of 0: one region covering all centres
  expect_equal(nrow(res$regions), 1L)
  seg <- seqprofile_wheel_segment(res, 1L)
  expect_false(seg$truncated)
  expect_equal(c(seg$start, seg$end), c(1L, 30L))
  expect_equal(seg$sequence, sp_ins_b)

  # Longer than max_len: 18 residues around the peak, inside the region span
  cut <- seqprofile_wheel_segment(res, 1L, max_len = 20L)
  expect_true(cut$truncated)
  expect_equal(nchar(cut$sequence), 18L)
  expect_gte(cut$start, 1L)
  expect_lte(cut$end, 30L)
  expect_true(res$regions$peak_pos[1] %in% cut$start:cut$end)

  expect_error(
    seqprofile_wheel_segment(res, 2L),
    class = "seqprofile_bad_param"
  )
})

# ── SeqProfile_server(): several algorithms, consensus, wheel, UniProt ──────
test_that("the module stacks algorithms, calls consensus and draws a wheel", {
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
      sequence = sp_ins_b,
      algorithm = c("antigenicity_kolaskar", "hydropathy_kd"),
      window_hydropathy_kd = 5,
      table_view = "regions",
      predict = 1
    )
    expect_named(results(), c("antigenicity_kolaskar", "hydropathy_kd"))
    expect_equal(results()$hydropathy_kd$params$window, 5L)
    # Untouched parameters fall back to the algorithm's defaults
    expect_equal(results()$antigenicity_kolaskar$params$window, 7L)

    # Two stacked tracks render (y and y2 axes of the subplot)
    expect_match(output$profile_plot, "yaxis2", fixed = TRUE)

    # The Summary box follows the chosen algorithm, the first by default
    expect_equal(result()$algorithm_id, "antigenicity_kolaskar")
    session$setInputs(summary_algorithm = "hydropathy_kd")
    expect_equal(result()$algorithm_id, "hydropathy_kd")

    session$setInputs(table_view = "consensus")
    expect_equal(
      names(table_data())[6:7],
      c(
        "Mean: Antigenicity (Kolaskar-Tongaonkar)",
        "Mean: Hydropathy (Kyte-Doolittle)"
      )
    )
    expect_gt(nrow(table_data()), 0)
    expect_equal(nrow(consensus_bands()), nrow(table_data()))

    # The minimum length filters table and plot bands without a new Predict
    longest <- max(table_data()$Length)
    session$setInputs(consensus_min_len = longest)
    expect_true(all(table_data()$Length >= longest))
    expect_equal(nrow(consensus_bands()), nrow(table_data()))
    session$setInputs(consensus_min_len = NA)
    expect_error(table_data(), "minimum consensus length")
    expect_null(consensus_bands())
    expect_no_error(output$profile_plot)
    session$setInputs(consensus_min_len = 1)

    # Consensus needs two algorithms and a threshold for each
    session$setInputs(
      algorithm = c("hydropathy_kd", "hmoment_eisenberg"),
      predict = 2
    )
    expect_error(table_data(), "None is set for")
    session$setInputs(algorithm = "hmoment_eisenberg", predict = 3)
    expect_error(table_data(), "Select two or three algorithms")

    # Helical wheel: Eisenberg at 100 degrees, regions view, one clicked row
    session$setInputs(table_view = "regions")
    expect_true(wheel_available())
    expect_error(wheel_data(), "enter a threshold")
    session$setInputs(threshold_hmoment_eisenberg = 0, predict = 4)
    expect_error(wheel_data(), "Click a region")
    session$setInputs(summary_table_rows_selected = 1)
    expect_equal(wheel_data()$position, 1:30)
    expect_equal(attr(wheel_data(), "segment")$sequence, sp_ins_b)
    session$setInputs(angle_hmoment_eisenberg = "160", predict = 5)
    expect_false(wheel_available())
    session$setInputs(
      algorithm = "hydropathy_kd",
      summary_table_rows_selected = NULL,
      predict = 6
    )
    expect_false(wheel_available())
  })
})

test_that("Fetch from UniProt fills the sequence box", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("plotly")
  skip_if_not_installed("DT")
  suppressPackageStartupMessages(library(shiny))
  source(
    file.path(fa_app_dir, "modules", "mod_SeqProfile.r"),
    local = environment()
  )
  # Shadows the real function for the module sourced above: no network
  requested <- character(0)
  fetch_uniprot_fasta <- function(accession, ...) {
    requested <<- c(requested, accession)
    acc <- seqprofile_clean_accession(accession)
    out <- parse_uniprot_fasta(paste0(">sp|", acc, "|X_TEST Test OS=T\nGIVEQ"))
    out$accession <- acc
    out
  }
  sent <- list()
  updateTextAreaInput <- function(session, inputId, ..., value = NULL) {
    sent[[inputId]] <<- value
  }
  showNotification <- function(...) NULL

  shiny::testServer(SeqProfile_server, {
    session$setInputs(accession = "q07730", fetch_uniprot = 1)
    expect_equal(requested, "q07730")
    expect_equal(sent$sequence, ">sp|Q07730|X_TEST Test OS=T\nGIVEQ")

    # A malformed accession is reported and leaves the sequence alone
    sent <<- list()
    session$setInputs(accession = "not an accession", fetch_uniprot = 2)
    expect_length(sent, 0)
  })
})
