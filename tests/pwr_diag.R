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

f <- "mtx_prot_abundance.tsv"
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
  m$Condition <- sub("_r[0-9]+_b[0-9]+$", "", m$Sample)
  meta_all(m)
  cat("\nconditions:\n")
  print(condition_replicates())

  session$setInputs(exclude_samples = character(0))
  cat("no-op selection -> n excluded:", length(excluded_samples()), "\n")

  drop <- samples[grepl("^NaCl", samples)][-1]
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

  entries <- reactlog::reactlog()
  actions <- vapply(
    entries,
    function(e) if (is.null(e$action)) NA_character_ else e$action,
    character(1)
  )
  labels <- vapply(
    entries,
    function(e) if (is.null(e$label)) NA_character_ else e$label,
    character(1)
  )
  cat("\n--- top reactive executions (enter events) ---\n")
  print(head(sort(table(labels[actions == "enter"]), decreasing = TRUE), 25))
})
cat(
  "total seconds:",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
  "\n"
)
