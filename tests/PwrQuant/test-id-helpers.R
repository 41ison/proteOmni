test_that("is_uniprot_accession recognises accession formats", {
  expect_true(all(is_uniprot_accession(c("P00761", "A0A024RBG1", "Q9H0H5-2"))))
  expect_false(any(is_uniprot_accession(c("BSA", "ALBU_HUMAN", "P0076", ""))))
  expect_true(all(is_uniprot_accession(rownames(pq_raw))))
})

test_that("extract_uniprot_lookup_id handles common ID styles", {
  expect_equal(
    extract_uniprot_lookup_id(c(
      "sp|P02769|ALBU_BOVIN",
      "P00761;P00762",
      "P00761,Q9H0H5",
      "tr|A0A0B4J2D5|X OS=Homo sapiens GN=GATD3B",
      "P02769"
    )),
    c("P02769", "P00761", "P00761", "A0A0B4J2D5", "P02769")
  )
})

test_that("detect_go_keytype picks UNIPROT for the training IDs", {
  expect_equal(detect_go_keytype(rownames(pq_raw)), "UNIPROT")
  expect_equal(detect_go_keytype(c("ALB", "TP53", "GAPDH")), "SYMBOL")
  expect_equal(detect_go_keytype(c(NA, "")), "SYMBOL")
})

test_that("parse_ratio_fraction converts GeneRatio strings", {
  expect_equal(parse_ratio_fraction(c("3/12", "1/4", "bad")), c(0.25, 0.25, NA))
})

# Regression test: the lookup used `[[`, which errored with "subscript out of
# bounds" for unknown taxa before the informative stop() was reached.
test_that("ensure_orgdb rejects unsupported taxa without installing", {
  expect_error(ensure_orgdb("999999999"), "not in the supported")
})

test_that("organism tables are consistent", {
  expect_false(anyDuplicated(pwrquant_orgdb_organisms$taxon) > 0)
  expect_true(all(
    pwrquant_orgdb_organisms$taxon %in% names(pwrquant_orgdb_lookup)
  ))
  expect_equal(unname(pwrquant_orgdb_choices), pwrquant_orgdb_organisms$taxon)
})
