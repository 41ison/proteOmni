## Diagnostic harness (temporary): load mod_PwrQuant in a testServer and measure
## how the reactive graph behaves when an abundance matrix is uploaded.
suppressWarnings(suppressPackageStartupMessages({
  libs <- c(
    "shiny",
    "shinydashboard",
    "shinyjs",
    "tidyverse",
    "DT",
    "data.table",
    "impute",
    "ggpointdensity",
    "ggtext",
    "ggsci",
    "viridis",
    "RColorBrewer",
    "ggfortify",
    "gridExtra",
    "grid",
    "GGally",
    "limma",
    "patchwork",
    "ComplexHeatmap",
    "compiler",
    "httr",
    "jsonlite"
  )
  for (l in libs) {
    library(l, character.only = TRUE)
  }
}))

options(shiny.reactlog = TRUE)
source("app/modules/mod_PwrQuant.r")

# MaxLFQ matrix, 3 conditions (m0, m1, m2) x 3 replicates.
f <- "app/training/mtx_maxLFQ.tsv"
fi <- data.frame(
  name = basename(f),
  size = file.size(f),
  type = "text/tab-separated-values",
  datapath = f
)

t0 <- Sys.time()
testServer(PwrQuant_server, {
  cat("--- setting matrix_file ---\n")
  ta <- Sys.time()
  session$setInputs(matrix_file = fi)
  cat(
    "upload flush seconds:",
    round(as.numeric(difftime(Sys.time(), ta, units = "secs")), 1),
    "\n"
  )
  cat("dim kept:", paste(dim(raw_matrix()), collapse = " x "), "\n")

  samples <- colnames(raw_matrix_full())
  m <- meta_all()
  m$Condition <- sub("_r[0-9]+$", "", m$Sample)
  meta_all(m)
  cat("\nconditions:\n")
  print(condition_replicates())

  session$setInputs(exclude_samples = character(0))
  cat("no-op selection -> n excluded:", length(excluded_samples()), "\n")

  # Leave a single m2 replicate so the thin-condition path is exercised.
  drop <- samples[grepl("^m2_", samples)][-1]
  session$setInputs(exclude_samples = drop)
  cat(
    "removed",
    length(drop),
    "-> dim kept:",
    paste(dim(raw_matrix()), collapse = " x "),
    "\n"
  )
  cat("thin conditions:\n")
  print(thin_conditions())
  cat("emptied:", paste(emptied_conditions(), collapse = ", "), "\n")
  cat("downstream metadata rows:", nrow(meta_edit_df()), "\n")

  # The raw log comes from shiny; the reactlog package only renders it.
  entries <- shiny::reactlog()
  field <- function(name) {
    vapply(
      entries,
      function(e) if (is.null(e[[name]])) NA_character_ else e[[name]][1],
      character(1)
    )
  }
  actions <- field("action")
  react_ids <- field("reactId")
  labels <- field("label")

  # "enter" events carry only a reactId; the label lives on the "define" (or a
  # later "updateNodeLabel") entry, so keep the last label seen for each id.
  has_label <- actions %in% c("define", "updateNodeLabel") & !is.na(labels)
  id_labels <- tapply(labels[has_label], react_ids[has_label], function(x) {
    x[length(x)]
  })
  # Unnamed reactives are labelled with their deparsed code, so flatten and
  # truncate it; the reactId keeps reactives with the same opening apart.
  id_labels[] <- sprintf(
    "%s [%s]",
    strtrim(gsub("\\s+", " ", id_labels), 60),
    names(id_labels)
  )

  entered <- react_ids[actions == "enter"]
  entered_labels <- unname(id_labels[entered])
  entered_labels[is.na(entered_labels)] <- entered[is.na(entered_labels)]
  cat("\n--- top reactive executions (enter events) ---\n")
  print(head(sort(table(entered_labels), decreasing = TRUE), 25))
})
cat(
  "total seconds:",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
  "\n"
)
