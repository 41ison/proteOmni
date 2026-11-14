test_that("count_theoretical_peptides counts tryptic peptides within the length window", {
  expect_equal(
    count_theoretical_peptides(fa_seqs),
    c(P1TEST = 2L, Q2TEST = 3L, plain_id = 0L, P3TEST = 1L)
  )
  # Missed cleavages add the sums of consecutive fragments.
  expect_equal(
    unname(count_theoretical_peptides(fa_seqs, max_missed = 1)),
    c(4L, 5L, 0L, 2L)
  )
  expect_equal(
    unname(count_theoretical_peptides(fa_seqs, max_missed = 2)),
    c(5L, 6L, 0L, 2L)
  )
  # P1TEST with 2 missed cleavages: lengths 15, 7, 17, 22, 24
  expect_equal(
    unname(count_theoretical_peptides(
      fa_seqs["P1TEST"],
      min_len = 8,
      max_len = 22,
      max_missed = 2
    )),
    3L
  )
})

test_that("count_theoretical_peptides returns NA for empty or missing sequences", {
  expect_identical(
    count_theoretical_peptides(c(a = "", b = NA, c = "AAAAAAK")),
    c(a = NA_integer_, b = NA_integer_, c = 1L)
  )
})

test_that("count_theoretical_peptides agrees with in_silico_digest", {
  for (mc in 0:2) {
    digest <- in_silico_digest(fa_seqs, max_missed = mc)
    expect_equal(
      count_theoretical_peptides(fa_seqs, max_len = 50, max_missed = mc),
      setNames(tabulate(match(digest$protein, fa_ids), length(fa_ids)), fa_ids)
    )
  }
})

test_that("build_theoretical_peptide_counts returns one row per protein", {
  theo <- build_theoretical_peptide_counts(fa_seqs, max_missed = 1)
  expect_s3_class(theo, "tbl_df")
  expect_equal(theo$protein_id, fa_ids)
  expect_equal(theo$n_theoretical_peptides, c(4L, 5L, 0L, 2L))
})

test_that("in_silico_digest cleaves after K/R but not before P", {
  digest <- in_silico_digest(fa_seqs, max_missed = 0)
  expect_setequal(names(digest), c("protein", "peptide"))
  expect_setequal(
    digest$peptide[digest$protein == "P1TEST"],
    c("AAAAAARPEPTIDEK", "GGGGGGR")
  )
  expect_setequal(
    digest$peptide[digest$protein == "Q2TEST"],
    c("AAAAAAK", "CCCCCCK", "DDDDDDK")
  )
  # "MK" is below the minimum peptide length of 6.
  expect_false("plain_id" %in% digest$protein)
})

test_that("in_silico_digest adds missed-cleavage peptides", {
  digest <- in_silico_digest(fa_seqs["P1TEST"], max_missed = 2)
  expect_setequal(
    digest$peptide,
    c(
      "AAAAAARPEPTIDEK",
      "GGGGGGR",
      "MKAAAAAARPEPTIDEK",
      "AAAAAARPEPTIDEKGGGGGGR",
      "MKAAAAAARPEPTIDEKGGGGGGR"
    )
  )
  expect_equal(anyDuplicated(digest), 0L)
})

test_that("protease_regex_from_rule builds a C-terminal cleavage rule", {
  expect_equal(protease_regex_from_rule("KR", "P"), "(?<=[KR])(?![P])")
  expect_equal(
    protease_regex_from_rule("k, r", "p", sense = "c"),
    "(?<=[KR])(?![P])"
  )
  expect_equal(protease_regex_from_rule("KR"), "(?<=[KR])")
  expect_equal(
    count_theoretical_peptides(
      fa_seqs,
      protease_regex = protease_regex_from_rule("KR", "P")
    ),
    count_theoretical_peptides(fa_seqs)
  )
  # Without the proline rule P1TEST is also cut at R-P:
  # MK | AAAAAAR | PEPTIDEK | GGGGGGR
  expect_equal(
    unname(count_theoretical_peptides(
      fa_seqs["P1TEST"],
      protease_regex = protease_regex_from_rule("KR")
    )),
    3L
  )
})

test_that("protease_regex_from_rule builds an N-terminal cleavage rule", {
  split_with <- function(x, rx) strsplit(x, rx, perl = TRUE)[[1]]

  lys_n <- protease_regex_from_rule("K", sense = "N")
  expect_equal(split_with("AAKDDKCC", lys_n), c("AA", "KDD", "KCC"))
  # A cut residue at position 1 must not produce an empty first fragment.
  expect_equal(split_with("KAAKDD", lys_n), c("KAA", "KDD"))
  expect_equal(
    unname(count_theoretical_peptides(
      "AAAAAAKDDDDDDKCCCCCC",
      protease_regex = lys_n
    )),
    3L
  )

  blocked <- protease_regex_from_rule("K", "P", sense = "n")
  expect_equal(split_with("AAKDDPKCC", blocked), c("AA", "KDDPKCC"))
})

test_that("protease_regex_from_rule returns NULL without a cut residue", {
  expect_null(protease_regex_from_rule(""))
  expect_null(protease_regex_from_rule(NULL))
  expect_null(protease_regex_from_rule(NA))
  expect_null(protease_regex_from_rule("-"))
})

test_that("classify_peptides labels proteotypic, shared and unmapped peptides", {
  digest <- in_silico_digest(fa_seqs, max_missed = 0)
  out <- classify_peptides(
    c("AAAAAAK", "GGGGGGR", "NOTINFASTA", "AAAAAAK"),
    digest
  )
  expect_equal(out$peptide, c("AAAAAAK", "GGGGGGR", "NOTINFASTA"))
  expect_equal(out$classification, c("Proteotypic", "Shared", "Unmapped"))
  expect_equal(out$n_proteins, c(1L, 2L, NA))
  expect_equal(out$mapped_proteins[1], "Q2TEST")
  expect_setequal(
    strsplit(out$mapped_proteins[2], ";", fixed = TRUE)[[1]],
    c("P1TEST", "P3TEST")
  )
})

abundance_input <- data.frame(
  filename = c(rep("f1", 6), "f2", "f2"),
  pep = c("a", "b", "c", "d", "s", "e", "a", "n"),
  inten = c(1, 2, 3, 4, 100, 7, 10, NA),
  class = c(rep("Proteotypic", 4), "Shared", rep("Proteotypic", 3)),
  prot = c(rep("P1", 4), "P1;P2", "P2", "P1", "P2")
)

test_that("compute_protein_abundance sums proteotypic peptides per file", {
  mat <- compute_protein_abundance(
    abundance_input,
    "pep",
    "inten",
    "class",
    "prot"
  )
  expect_equal(names(mat), c("prot", "f1", "f2"))
  expect_equal(mat$prot, c("P1", "P2"))
  # Shared peptide (100) excluded; P2 has no usable value in f2.
  expect_equal(mat$f1, c(10, 7))
  expect_equal(mat$f2, c(10, NA))
})

test_that("compute_protein_abundance top3 keeps the three most intense peptides", {
  mat <- compute_protein_abundance(
    abundance_input,
    "pep",
    "inten",
    "class",
    "prot",
    method = "top3"
  )
  expect_equal(mat$f1, c(9, 7))
  expect_equal(mat$f2, c(10, NA))
})

test_that("compute_protein_abundance handles missing inputs", {
  expect_error(
    compute_protein_abundance(abundance_input, "pep", "nope", "class", "prot"),
    "Missing columns: nope"
  )
  expect_null(compute_protein_abundance(
    abundance_input[abundance_input$class == "Shared", ],
    "pep",
    "inten",
    "class",
    "prot"
  ))
  no_file <- abundance_input[1:4, names(abundance_input) != "filename"]
  expect_equal(
    names(compute_protein_abundance(no_file, "pep", "inten", "class", "prot")),
    c("prot", "Unknown")
  )
})
