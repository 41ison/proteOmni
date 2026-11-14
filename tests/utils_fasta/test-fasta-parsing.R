test_that("read_fasta_custom joins wrapped lines and names entries by accession", {
  seqs <- read_fasta_custom(local_fasta())
  expect_equal(names(seqs), fa_ids)
  # Sequences are returned as written; normalisation happens downstream.
  expect_equal(
    unname(seqs[c("P1TEST", "plain_id")]),
    c(unname(fa_seqs[1]), "MK")
  )
  expect_equal(unname(seqs["Q2TEST"]), tolower(unname(fa_seqs["Q2TEST"])))
})

test_that("read_fasta_custom handles files without headers or sequences", {
  expect_equal(read_fasta_custom(local_fasta(c("MKAAA", "CCC"))), character(0))
  expect_equal(
    read_fasta_custom(local_fasta(c(">a", ">b", "MK"))),
    c(a = "", b = "MK")
  )
})

# Regression test: the UniProt check looked for "sp" or "tr" anywhere in the
# first word, so ">trypsin_pig Porcine trypsin" was named "Porcine".
test_that("read_fasta_custom only treats sp|/tr| database tags as UniProt-style", {
  seqs <- read_fasta_custom(local_fasta(c(
    ">trypsin_pig Porcine trypsin",
    "MKAA",
    ">HSP90 heat shock protein",
    "CC",
    ">contam_sp|P02769|ALBU_BOVIN Albumin",
    "DD",
    ">gi|12345|ref|NP_000001.1| some protein",
    "EE",
    ">transporter|X1 desc",
    "FF"
  )))
  expect_equal(
    names(seqs),
    c("trypsin_pig", "HSP90", "P02769", "gi", "transporter")
  )
  expect_equal(unname(seqs), c("MKAA", "CC", "DD", "EE", "FF"))
})

test_that("fasta_seq_normalize strips whitespace and upper-cases", {
  expect_equal(fasta_seq_normalize(" mk aa\n"), "MKAA")
  expect_equal(fasta_seq_normalize(c("mk", "cc")), "MK")
})

test_that("fasta_header_tokens splits only the identifier field", {
  expect_equal(
    fasta_header_tokens(c(
      "sp|P02769|ALBU_BOVIN Albumin OS=Bos|taurus",
      "plain x"
    )),
    list(c("sp", "P02769", "ALBU_BOVIN"), "plain")
  )
})

test_that("fasta_index_exact requires a whole-token match", {
  headers <- c("sp|P02310|H1_TEST first", "sp|P0231|H2_TEST second")
  # A substring search for "P0231" would hit the first header.
  expect_equal(fasta_index_exact("P0231", headers), 2L)
  expect_equal(fasta_index_exact("P02310", headers), 1L)
  expect_equal(fasta_index_exact("H1_TEST", headers), 1L)
  expect_identical(fasta_index_exact("P023", headers), NA_integer_)
})

test_that("fasta_index_exact falls back to a case-insensitive match", {
  expect_equal(
    fasta_index_exact("p0231", c("sp|P02310|A x", "sp|P0231|B x")),
    2L
  )
  # An exact-case hit wins over an earlier case-insensitive one.
  expect_equal(fasta_index_exact("abc", c("ABC x", "abc y")), 2L)
})

test_that("protein_ids_contain matches complete identifiers in a protein field", {
  ids <- c("sp|P1|A_HUMAN; sp|P2|B_HUMAN", "sp|P10|C_HUMAN", NA, "P1")
  expect_equal(protein_ids_contain(ids, "P1"), c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(
    protein_ids_contain(ids, "B_HUMAN"),
    c(TRUE, FALSE, FALSE, FALSE)
  )
  expect_equal(protein_ids_contain(character(0), "P1"), logical(0))
})

test_that("fasta_accessions takes the UniProt accession or the first word", {
  expect_equal(
    fasta_accessions(c(
      ">sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus",
      "tr|A0A024RBG1|X_HUMAN",
      "ENSP00000354587.3 pep chromosome:GRCh38",
      "  plain  "
    )),
    c("P02769", "A0A024RBG1", "ENSP00000354587.3", "plain")
  )
})

test_that("fasta_header_ids returns the first word of each header", {
  expect_equal(
    fasta_header_ids(c(
      ">sp|P04075|ALDOA_HUMAN Fructose-bisphosphate aldolase A",
      "contam_sp|P02769|ALBU_BOVIN Albumin",
      "plain"
    )),
    c("sp|P04075|ALDOA_HUMAN", "contam_sp|P02769|ALBU_BOVIN", "plain")
  )
})

test_that("fasta_header_genes reads the GN= field and returns NA without one", {
  expect_equal(
    fasta_header_genes(c(fa_lines[c(1, 4, 6)], "x GN=TP53")),
    c("ALB", "XYZ1", NA, "TP53")
  )
})

test_that("fasta_sequences_by_header_id normalises and drops duplicate IDs", {
  # A named character vector stands in for the AAStringSet used by the app.
  fasta <- c("sp|P1|A_T desc" = "mk aa", "sp|P1|A_T other" = "CCC", "x" = "dd")
  expect_equal(
    fasta_sequences_by_header_id(fasta),
    c("sp|P1|A_T" = "MKAA", x = "DD")
  )
})

test_that("match_fasta_sequences matches the header ID, then the accession", {
  sequences <- c("sp|P1|A_T" = "MK", "sp|P2|B_T" = "CC", plain = "DD")
  expect_equal(
    match_fasta_sequences(
      protein = c("sp|P2|B_T", "tr|P1|OTHER_T", "nope", "plain"),
      accession = c("P2", "P1", "ZZ", "plain"),
      sequences = sequences
    ),
    c(2L, 1L, NA, 3L)
  )
})

test_that("ibaq_fasta_sequences reads a FASTA path into accession -> sequence", {
  expect_equal(ibaq_fasta_sequences(local_fasta()), fa_seqs)
})

test_that("ibaq_fasta_sequences accepts an already-parsed FASTA", {
  expected <- c(P1 = "MK", x = "DD")
  # seqinr::read.fasta(as.string = TRUE) style: list named by full header ID;
  # a repeated accession keeps its first sequence.
  parsed <- list("sp|P1|A_T" = "mk", "sp|P1|A2_T" = "CC", x = "d d")
  expect_equal(ibaq_fasta_sequences(parsed), expected)
  expect_equal(ibaq_fasta_sequences(unlist(parsed)), expected)
})
