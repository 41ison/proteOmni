test_that("cluster_heatmap_rows clusters a z-scored training slice", {
  mat <- pq_log2_complete[1:60, ]
  mat_z <- t(scale(t(mat)))
  res <- cluster_heatmap_rows(mat_z, k = 4)
  expect_equal(res$k, 4L)
  expect_equal(names(res$cluster), rownames(mat_z))
  expect_equal(levels(res$cluster), as.character(1:4))
  # Clusters are numbered in dendrogram order: cluster 1 is drawn first
  expect_equal(as.character(res$cluster[res$hc$order][1]), "1")
  expect_equal(unique(as.integer(res$cluster[res$hc$order])), 1:4)
})

test_that("cluster_heatmap_rows clamps k and handles NA/small input", {
  mat_z <- t(scale(t(pq_log2_complete[1:5, ])))
  expect_equal(cluster_heatmap_rows(mat_z, k = 50)$k, 4L)
  expect_equal(cluster_heatmap_rows(mat_z, k = 1)$k, 2L)
  expect_equal(cluster_heatmap_rows(mat_z, k = "x")$k, 4L)
  mat_z[1, 2] <- NA
  expect_length(cluster_heatmap_rows(mat_z, k = 2)$cluster, 5L)
  expect_null(cluster_heatmap_rows(mat_z[1:2, ]))
})

test_that("get_palette_colors returns n valid colours for every palette", {
  all_pals <- unname(unlist(pwrquant_palette_choices))
  for (pal in all_pals) {
    for (n in c(1L, 3L, 15L)) {
      cols <- get_palette_colors(pal, n)
      expect_length(cols, n)
      expect_false(anyNA(cols))
      expect_silent(grDevices::col2rgb(cols))
    }
  }
  expect_length(get_palette_colors("not_a_palette", 3), 3)
})

test_that("string_halo_colors maps regulation status", {
  expect_equal(
    string_halo_colors(c("Increased", "Decreased", "Not significant", NA)),
    c("#d62728", "#1f77b4", "#7f7f7f", "#7f7f7f")
  )
})

test_that("format_unsupported_species_error builds a readable message", {
  msg <- format_unsupported_species_error(list(
    taxon_species = 12345,
    org_name = "Sporothrix brasiliensis",
    suggestions = data.frame(
      official_name_NCBI = "Sporothrix schenckii",
      taxon_id = 29908
    )
  ))
  expect_match(msg, "Sporothrix brasiliensis (taxon 12345)", fixed = TRUE)
  expect_match(msg, "Sporothrix schenckii (taxon 29908)", fixed = TRUE)

  msg2 <- format_unsupported_species_error(list(
    taxon_species = 1,
    org_name = NA_character_,
    suggestions = NULL
  ))
  expect_match(msg2, "taxon 1")
  expect_match(msg2, "No closely related organism")
})
