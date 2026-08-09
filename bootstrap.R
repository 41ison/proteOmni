## proteOmni — dependency bootstrap
##
## Runs with BASE R ONLY. Nothing here may reference an add-on package before
## it has been installed, so this file can execute on a completely fresh R
## installation (no shiny, no BiocManager, no pak).
##
## Sourced by run.R before Shiny is loaded, and defensively by proteOmni.r so
## that launching the app directly from an IDE still works.

if (!isTRUE(getOption("proteOmni.bootstrapped"))) {
  options(repos = c(CRAN = "https://cran.rstudio.com/"))

  # ── Writable library ───────────────────────────────────────────────────────
  # Rscript does not offer to create the personal library the way interactive R
  # does, so a fresh install has nowhere writable to put packages.
  .omni_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(.omni_lib)) {
    .omni_lib <- file.path(
      path.expand("~"),
      "R",
      paste0(R.version$platform, "-library"),
      paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")
    )
  }
  .omni_lib <- strsplit(.omni_lib, .Platform$path.sep)[[1]][1]
  if (!dir.exists(.omni_lib)) {
    dir.create(.omni_lib, recursive = TRUE, showWarnings = FALSE)
  }
  if (dir.exists(.omni_lib) && file.access(.omni_lib, 2) == 0) {
    .libPaths(c(.omni_lib, .libPaths()))
  }
  if (file.access(.libPaths()[1], 2) != 0) {
    stop(
      "No writable R library found. Set R_LIBS_USER to a writable directory ",
      "and restart proteOmni.",
      call. = FALSE
    )
  }

  # ── Helpers ────────────────────────────────────────────────────────────────
  .omni_has <- function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  }

  .omni_missing <- character(0)

  .omni_install <- function(pkgs, installer) {
    for (pkg in pkgs) {
      if (.omni_has(pkg)) next
      message("[proteOmni] Installing ", pkg, " ...")
      ok <- tryCatch(
        {
          installer(pkg)
          .omni_has(pkg)
        },
        error = function(e) {
          message("[proteOmni] Failed to install ", pkg, ": ", conditionMessage(e))
          FALSE
        }
      )
      if (!ok) .omni_missing <<- c(.omni_missing, pkg)
    }
  }

  .omni_cran <- function(pkg) install.packages(pkg)
  .omni_bioc <- function(pkg) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }

  # ── CRAN ───────────────────────────────────────────────────────────────────
  # `pak` is listed explicitly: it is needed to install the GitHub-only `diann`.
  .omni_install(
    c(
      "shiny", "shinydashboard", "shinyjs", "fresh",
      "BiocManager", "pak", "devtools",
      "tidyverse", "tidytext", "janitor", "data.table",
      "ggpointdensity", "ggtext", "ggrepel", "ggseqlogo", "GGally", "ggsci",
      "ggfortify", "plotly", "viridis", "RColorBrewer", "scales", "patchwork",
      "lsa", "vegan", "seqinr", "zip", "DT", "colourpicker", "R6",
      "gridExtra", "lavaan", "naniar", "missForest",
      "httr", "jsonlite"
    ),
    .omni_cran
  )

  # arrow needs these set at install time to get a full-featured libarrow.
  Sys.setenv(LIBARROW_MINIMAL = "false", ARROW_WITH_ZSTD = "ON")
  .omni_install("arrow", .omni_cran)

  # ── Bioconductor ───────────────────────────────────────────────────────────
  if (.omni_has("BiocManager")) {
    .omni_install(
      c(
        "limma", "Biostrings", "sva", "impute", "pcaMethods",
        "ComplexHeatmap", "STRINGdb", "clusterProfiler", "enrichplot",
        "AnnotationDbi", "mixOmics",
        # Baseline OrgDb (human); other organisms are installed on demand by
        # the PwrQuant enrichment module via ensure_orgdb().
        "org.Hs.eg.db"
      ),
      .omni_bioc
    )
  } else {
    .omni_missing <- c(.omni_missing, "BiocManager (Bioconductor deps skipped)")
  }

  # ── GitHub ─────────────────────────────────────────────────────────────────
  if (!.omni_has("diann")) {
    .omni_install("diann", function(pkg) {
      if (.omni_has("pak")) {
        pak::pak("https://github.com/vdemichev/diann-rpackage", ask = FALSE)
      } else {
        remotes::install_github("vdemichev/diann-rpackage", upgrade = "never")
      }
    })
  }

  # ── Report ─────────────────────────────────────────────────────────────────
  if (length(.omni_missing)) {
    message(
      "\n[proteOmni] The following packages could not be installed:\n  - ",
      paste(unique(.omni_missing), collapse = "\n  - "),
      "\nproteOmni may fail to start. See the messages above for the cause ",
      "(commonly a missing system library or no internet connection).\n"
    )
  } else {
    message("[proteOmni] All dependencies are available.")
  }

  rm(.omni_lib, .omni_has, .omni_install, .omni_cran, .omni_bioc, .omni_missing)
  options(proteOmni.bootstrapped = TRUE)
}
