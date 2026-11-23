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

test_that("string_taxon_for maps OrgDb strain taxa to STRING's taxa", {
  expect_identical(string_taxon_for("559292"), 4932L)
  expect_identical(string_taxon_for("83333"), 511145L)
  expect_identical(string_taxon_for("10090"), 10090L)
  expect_identical(string_taxon_for(""), 9606L)
})

test_that("resolve_go_ids leaves human symbols and accessions unchanged", {
  skip_if_not_installed("org.Hs.eg.db")
  hs <- org.Hs.eg.db::org.Hs.eg.db
  res <- resolve_go_ids(c("ALB", "TP53", "GAPDH"), hs)
  expect_equal(res$keytype, "SYMBOL")
  expect_equal(unname(res$map), c("ALB", "TP53", "GAPDH"))
  expect_null(res$labels)
  expect_equal(resolve_go_ids(c("P02768", "P04637"), hs)$keytype, "UNIPROT")
})

# Regression test: org.Sc.sgd.db has no SYMBOL keytype, and yeast matrices mix
# gene names with systematic ORF names.
test_that("resolve_go_ids maps yeast gene names and ORFs to ORF keys", {
  skip_if_not_installed("org.Sc.sgd.db")
  sc <- org.Sc.sgd.db::org.Sc.sgd.db
  res <- resolve_go_ids(c("ACT1", "YBL036C", "NOT_A_GENE"), sc)
  expect_equal(res$keytype, "ORF")
  expect_equal(unname(res$map), c("YFL039C", "YBL036C", NA))
  expect_equal(res$labels[["YFL039C"]], "ACT1")
  expect_equal(to_go_keys(res, c("ACT1", NA, "NOT_A_GENE")), "YFL039C")
  expect_error(resolve_go_ids(c("ALB", "TP53"), sc), "could be matched")
  # UniProt accessions are relabelled with gene names.
  up <- resolve_go_ids("P60010", sc)
  expect_equal(unname(up$map), "YFL039C")
  expect_equal(up$labels[["YFL039C"]], "ACT1")
})

# Regression test: org.EcK12.eg.db has no UNIPROT keytype, so accessions are
# resolved through gene identifiers fetched from UniProt (stubbed here).
test_that("resolve_go_ids falls back to UniProt gene IDs without UNIPROT", {
  skip_if_not_installed("org.EcK12.eg.db")
  ec <- org.EcK12.eg.db::org.EcK12.eg.db
  expect_false("UNIPROT" %in% AnnotationDbi::keytypes(ec))
  stub <- function(accessions) {
    data.frame(
      id = c("P0A6F5", "P0A6F5", "P0A7V8"),
      candidate = c("93777681", "948665", "rpsD"),
      stringsAsFactors = FALSE
    )
  }
  res <- resolve_go_ids(c("P0A6F5", "P0A7V8", "P99999"), ec, stub)
  expect_equal(res$keytype, "ENTREZID")
  expect_equal(unname(res$map), c("948665", "947793", NA))
  expect_equal(res$labels[["948665"]], "groL")

  offline <- function(accessions) stop("offline")
  expect_error(
    resolve_go_ids(c("P0A6F5", "P0A7V8"), ec, offline),
    "UniProt lookup failed"
  )
})

# Regression test: enrichGO() needs ENTREZID, which org.Pf.plasmo.db cannot
# serve, so enrichment falls back to enricher() on the OrgDb's GO annotations.
test_that("enrichment works for an OrgDb without usable ENTREZID", {
  skip_if_not_installed("org.Pf.plasmo.db")
  skip_if_not_installed("clusterProfiler")
  pf <- suppressMessages(org.Pf.plasmo.db::org.Pf.plasmo.db)
  expect_false(orgdb_serves_entrez(pf))

  sets <- orgdb_go_term2gene(pf, "ORF", "BP")
  expect_named(sets$term2gene, c("term", "gene"))
  expect_true(all(sets$ontology == "BP"))

  # Genes of one term must come back enriched for that term.
  sizes <- table(sets$term2gene$term)
  term <- names(sizes)[sizes >= 20 & sizes <= 60][1]
  genes <- sets$term2gene$gene[sets$term2gene$term == term]
  res <- suppressMessages(run_per_contrast_enrichment(
    list(a = genes),
    pf,
    "ORF",
    ont = "ALL",
    labels = setNames(genes, genes)
  ))
  expect_true(term %in% res$ID)
  expect_equal(res$ONTOLOGY[res$ID == term], "BP")
  expect_equal(names(res)[1:3], c("ONTOLOGY", "ID", "Description"))
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
