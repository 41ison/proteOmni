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

# Regression test: the FASTA-header regex excluded the hyphen, so an isoform
# header fell back to the entry name ("NAME") instead of the accession.
test_that("extract_uniprot_lookup_id keeps isoform accessions", {
  ids <- extract_uniprot_lookup_id(c(
    "sp|P12345-2|NAME_ORG",
    "P12345-2",
    "sp|Q9H0H5-12|X_HUMAN;sp|P00761|TRYP_PIG"
  ))
  expect_equal(ids, c("P12345-2", "P12345-2", "Q9H0H5-12"))
  expect_true(all(is_uniprot_accession(ids)))
  expect_equal(strip_uniprot_isoform(ids), c("P12345", "P12345", "Q9H0H5"))
  expect_equal(strip_uniprot_isoform("A0A024RBG1"), "A0A024RBG1")
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

# ── UniProt GO annotation source ────────────────────────────────────────────

# Stub for fetch_uniprot_go()'s `fetcher`: serves rows of `known` and records
# every request in `calls`.
make_uniprot_go_stub <- function(known, release = "2026_04") {
  calls <- list()
  fetcher <- function(accessions) {
    calls[[length(calls) + 1L]] <<- accessions
    list(
      table = known[known$accession %in% accessions, , drop = FALSE],
      release = release
    )
  }
  list(fetcher = fetcher, calls = function() calls)
}

uniprot_go_known <- data.frame(
  accession = c("P11111", "P22222", "P33333"),
  gene = c("GENA", "", "GENC"),
  organism = c("Org one", "Org one", "Org two"),
  go_p = c(
    "protein folding [GO:0006457]; brand-new process [GO:9999999]",
    "",
    "protein folding [GO:0006457]"
  ),
  go_c = c("cytosol [GO:0005829]", "", ""),
  go_f = c("", "", "ATP binding [GO:0005524]"),
  stringsAsFactors = FALSE
)

local_uniprot_go_cache <- function(env = parent.frame()) {
  dir <- tempfile("uniprot_go_cache_test_")
  withr::defer(unlink(dir, recursive = TRUE), envir = env)
  dir
}

test_that("fetch_uniprot_go chunks, caches and records the release", {
  cache_dir <- local_uniprot_go_cache()
  stub <- make_uniprot_go_stub(uniprot_go_known)
  ids <- c("P11111-2", "P22222", "P33333", "P44444", "P11111", "GAPDH", NA)

  res <- fetch_uniprot_go(
    ids,
    fetcher = stub$fetcher,
    chunk_size = 2L,
    cache_dir = cache_dir
  )
  # Isoforms collapse to the canonical accession; non-accessions are never
  # sent (UniProt rejects the whole request for one malformed accession).
  expect_equal(res$accession, c("P11111", "P22222", "P33333", "P44444"))
  expect_equal(stub$calls(), list(c("P11111", "P22222"), c("P33333", "P44444")))
  expect_equal(res$found, c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(res$gene, c("GENA", "", "GENC", ""))
  expect_true(all(res$release == "2026_04"))
  expect_true(all(res$fetched == Sys.Date()))
  expect_true(file.exists(file.path(cache_dir, "uniprot_go_cache.rds")))

  # Second run: everything, including the unknown accession, is cached.
  again <- fetch_uniprot_go(ids, fetcher = stub$fetcher, cache_dir = cache_dir)
  expect_equal(again, res)
  expect_length(stub$calls(), 2)

  # A new R process (empty in-memory mirror) reads the disk cache.
  rm(list = ls(.uniprot_go_cache), envir = .uniprot_go_cache)
  from_disk <- fetch_uniprot_go(
    ids,
    fetcher = function(accessions) stop("must not be called"),
    cache_dir = cache_dir
  )
  expect_equal(from_disk, res)

  # Only accessions not seen before are requested.
  fetch_uniprot_go(
    c("P11111", "P55555"),
    fetcher = stub$fetcher,
    cache_dir = cache_dir
  )
  expect_equal(stub$calls()[[3]], "P55555")

  # Expired entries are fetched again; `refresh` forces it.
  fetch_uniprot_go(
    "P11111",
    fetcher = stub$fetcher,
    cache_dir = cache_dir,
    ttl_days = -1
  )
  fetch_uniprot_go(
    "P11111",
    fetcher = stub$fetcher,
    cache_dir = cache_dir,
    refresh = TRUE
  )
  expect_length(stub$calls(), 5)
})

test_that("fetch_uniprot_go fails clearly offline, unless cached", {
  cache_dir <- local_uniprot_go_cache()
  offline <- function(accessions) stop("Could not resolve host")
  expect_error(
    fetch_uniprot_go("P11111", fetcher = offline, cache_dir = cache_dir),
    "Could not download GO annotations from UniProt \\(Could not resolve host\\)"
  )
  stub <- make_uniprot_go_stub(uniprot_go_known)
  fetch_uniprot_go("P11111", fetcher = stub$fetcher, cache_dir = cache_dir)
  expect_equal(
    fetch_uniprot_go("P11111", fetcher = offline, cache_dir = cache_dir)$gene,
    "GENA"
  )
  # Expired entries are better than nothing when UniProt is unreachable.
  expect_message(
    stale <- fetch_uniprot_go(
      "P11111",
      fetcher = offline,
      cache_dir = cache_dir,
      ttl_days = -1
    ),
    "using expired cache entries"
  )
  expect_equal(stale$gene, "GENA")
})

test_that("uniprot_go_term2gene propagates annotations to GO.db ancestors", {
  skip_if_not_installed("GO.db")
  cache_dir <- local_uniprot_go_cache()
  stub <- make_uniprot_go_stub(uniprot_go_known)
  entries <- fetch_uniprot_go(
    c(uniprot_go_known$accession, "P44444"),
    fetcher = stub$fetcher,
    cache_dir = cache_dir
  )
  sets <- uniprot_go_term2gene(entries, "ALL")
  expect_named(sets$term2gene, c("term", "gene"))
  expect_named(sets$term2name, c("term", "name"))
  genes_of <- function(term) {
    sort(sets$term2gene$gene[sets$term2gene$term == term])
  }

  # Direct annotation, and its root, which UniProt does not list.
  expect_equal(genes_of("GO:0006457"), c("P11111", "P33333"))
  expect_equal(genes_of("GO:0008150"), c("P11111", "P33333"))
  expect_equal(genes_of("GO:0005575"), "P11111")
  expect_equal(genes_of("GO:0003674"), "P33333")
  expect_false("all" %in% sets$term2gene$term)
  expect_false(anyNA(sets$term2gene$term))
  expect_false(anyDuplicated(sets$term2gene) > 0)

  expect_equal(unname(sets$ontology["GO:0008150"]), "BP")
  expect_equal(unname(sets$ontology["GO:0005829"]), "CC")
  expect_equal(unname(sets$ontology["GO:0005524"]), "MF")
  expect_setequal(names(sets$ontology), sets$term2name$term)

  name_of <- function(term) sets$term2name$name[sets$term2name$term == term]
  expect_equal(name_of("GO:0008150"), "biological_process")
  # A GO ID missing from the installed GO.db keeps UniProt's description.
  expect_equal(name_of("GO:9999999"), "brand-new process")
  expect_equal(genes_of("GO:9999999"), "P11111")
  expect_equal(sets$missing_from_go_db, "GO:9999999")

  # Gene name, or the accession when UniProt has none.
  expect_equal(
    sets$labels[c("P11111", "P22222", "P44444")],
    c(P11111 = "GENA", P22222 = "P22222", P44444 = "P44444")
  )

  mf <- uniprot_go_term2gene(entries, "MF")
  expect_true(all(mf$ontology == "MF"))
  expect_equal(unique(mf$term2gene$gene), "P33333")
  expect_error(
    uniprot_go_term2gene(entries[entries$accession == "P22222", ], "ALL"),
    "GO annotations in UniProt"
  )
})

test_that("prepare_uniprot_go builds the universe and rejects gene symbols", {
  skip_if_not_installed("GO.db")
  cache_dir <- local_uniprot_go_cache()
  stub <- make_uniprot_go_stub(uniprot_go_known)
  expect_message(
    up <- prepare_uniprot_go(
      c("sp|P11111-2|A_ORG", "P22222;P33333", "P33333", "CONTAM_X"),
      fetcher = stub$fetcher,
      cache_dir = cache_dir
    ),
    "UniProt release\\(s\\): 2026_04"
  )
  expect_equal(up$universe, c("P11111", "P22222", "P33333"))
  expect_equal(
    up$map,
    c(
      "P11111-2" = "P11111",
      P22222 = "P22222",
      P33333 = "P33333",
      CONTAM_X = NA
    )
  )
  expect_equal(to_go_keys(up, c("P11111-2", "CONTAM_X")), "P11111")

  expect_error(
    prepare_uniprot_go(
      c("ALB", "TP53", "P11111"),
      fetcher = stub$fetcher,
      cache_dir = cache_dir
    ),
    "needs UniProt accessions"
  )
})

# Regression test: with ont = "ALL" the enricher() branch pooled BP, CC and MF
# in one call, inflating the background and the BH adjustment. enrichGO tests
# each ontology separately.
test_that("run_per_contrast_enrichment tests pre-built go_sets per ontology", {
  skip_if_not_installed("clusterProfiler")
  genes <- sprintf("P%05d", 1:300)
  set_of <- function(term, idx) {
    data.frame(term = term, gene = genes[idx], stringsAsFactors = FALSE)
  }
  # BP annotates all 300 genes, CC 150, MF 96, so backgrounds differ.
  t2g <- rbind(
    set_of("GO:BP1", 1:20),
    set_of("GO:BP2", 21:300),
    set_of("GO:BP3", seq(1, 300, by = 7)),
    set_of("GO:CC1", 1:20),
    set_of("GO:CC2", 21:150),
    set_of("GO:MF1", 5:24),
    set_of("GO:MF2", 25:100)
  )
  terms <- unique(t2g$term)
  go_sets <- list(
    term2gene = t2g,
    term2name = data.frame(
      term = terms,
      name = paste("name of", terms),
      stringsAsFactors = FALSE
    ),
    ontology = setNames(substr(terms, 4, 5), terms),
    labels = setNames(sub("^P", "gene", genes), genes)
  )
  hits <- genes[1:15]
  run <- function(ont, sets = go_sets) {
    # qvalue() cannot be estimated from two or three p-values and warns.
    suppressWarnings(suppressMessages(run_per_contrast_enrichment(
      list(a = hits, empty = character(0)),
      universe = genes,
      ont = ont,
      pvalue_cutoff = 1,
      qvalue_cutoff = 1,
      go_sets = sets
    )))
  }
  all_res <- run("ALL")
  expect_equal(names(all_res)[1:3], c("ONTOLOGY", "ID", "Description"))
  expect_equal(all_res$ONTOLOGY, substr(all_res$ID, 4, 5))
  expect_setequal(unique(all_res$ONTOLOGY), c("BP", "CC", "MF"))
  expect_equal(levels(all_res$Group), "a")

  # Each ontology has its own annotated background ...
  bg <- tapply(sub(".*/", "", all_res$BgRatio), all_res$ONTOLOGY, unique)
  expect_equal(bg[["BP"]], "300")
  expect_equal(bg[["CC"]], "150")
  expect_equal(bg[["MF"]], "96")

  # ... and BH adjustment: identical to running that ontology alone.
  for (o in c("BP", "CC", "MF")) {
    single <- run(o)
    expect_false("ONTOLOGY" %in% names(single))
    pooled <- all_res[all_res$ONTOLOGY == o, ]
    expect_equal(pooled$ID, single$ID)
    expect_equal(pooled$pvalue, single$pvalue)
    expect_equal(pooled$p.adjust, single$p.adjust)
  }

  # geneID is relabelled with the gene names carried by go_sets.
  bp1 <- strsplit(all_res$geneID[all_res$ID == "GO:BP1"], "/", fixed = TRUE)[[
    1
  ]]
  expect_setequal(bp1, sprintf("gene%05d", 1:15))
  expect_equal(all_res$Description[all_res$ID == "GO:BP1"], "name of GO:BP1")

  # An ontology without gene sets is skipped rather than failing.
  no_mf <- go_sets
  no_mf$term2gene <- t2g[!grepl("MF", t2g$term), ]
  expect_setequal(unique(run("ALL", no_mf)$ONTOLOGY), c("BP", "CC"))
  expect_error(
    run_per_contrast_enrichment(list(a = hits)),
    "orgdb.*go_sets"
  )
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
