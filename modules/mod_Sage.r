## ============================================================
## mod_Sage.r  —  Sage DDA/DIA results viewer
## ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyjs)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(data.table)
  library(arrow)
  library(ggpointdensity)
  library(ggtext)
  library(DT)
})

# ── Helper functions ──────────────────────────────────────────────────────────

theme_sage <- function(...) {
  theme_bw(...) +
    theme(
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.x = element_text(
        angle = 65,
        hjust = 1,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.title = element_text(size = 11, face = "bold", color = "black"),
      legend.position = "bottom",
      legend.title = element_text(
        size = 10,
        face = "bold",
        color = "black",
        hjust = 0.5
      ),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

sp_wrap_sage <- function(sid, ui_el) {
  div(
    class = "plot-wrap",
    tags$div(
      class = "spinner-overlay",
      id = sid,
      icon("spinner", class = "fa-spin")
    ),
    ui_el
  )
}

# ── Sage folder discovery ─────────────────────────────────────────────────────

#' Locate Sage output files under a user-supplied folder.
#'
#' Searches recursively for \code{results.sage.parquet} and
#' \code{lfq.parquet}. The \code{.tsv} equivalents are accepted as a fallback
#' when no parquet exists. When several matches are found, the one with the
#' shortest path (closest to the given folder) is used, and the parquet is
#' preferred over the tsv within the same directory.
#'
#' @param path Character. Folder to search.
#' @return Named list with elements \code{results} and \code{lfq}
#'   (full path or \code{NA}).
SAGE_FILE_TARGETS <- c(results = "results\\.sage", lfq = "lfq")

find_sage_files <- function(path) {
  targets <- SAGE_FILE_TARGETS
  path <- path.expand(trimws(path))
  if (!nzchar(path) || !dir.exists(path)) {
    return(setNames(
      as.list(rep(NA_character_, length(targets))),
      names(targets)
    ))
  }
  lapply(targets, function(f) {
    hits <- list.files(
      path,
      pattern = paste0("^", f, "\\.(parquet|tsv)$"),
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(hits) == 0) {
      return(NA_character_)
    }
    # Shallowest directory first, then prefer parquet within it
    hits <- hits[order(nchar(dirname(hits)))]
    hits <- hits[dirname(hits) == dirname(hits[1])]
    is_pq <- grepl("\\.parquet$", hits)
    if (any(is_pq)) hits[is_pq][1] else hits[1]
  })
}

#' Read a Sage table from parquet or tsv according to its extension.
read_sage_table <- function(path) {
  if (grepl("\\.parquet$", path)) {
    as.data.frame(arrow::read_parquet(path))
  } else {
    as.data.frame(data.table::fread(path, sep = "\t", header = TRUE))
  }
}

# ── Sage LFQ -> protein abundance matrix ──────────────────────────────────────

#' Columns of the Sage LFQ table that are not per-sample intensities.
SAGE_LFQ_ID_COLS <- c(
  "peptide",
  "stripped_peptide",
  "charge",
  "proteins",
  "is_decoy",
  "q_value",
  "score",
  "spectral_angle"
)

#' Protein roll-up methods available for the LFQ matrix.
SAGE_ROLLUP_METHODS <- c(
  "Sum of peptide intensities" = "sum",
  "Median of peptide intensities" = "median",
  "Mean of top-3 peptides" = "top3"
)

#' Bring the Sage LFQ table to long format (one row per peptide x file).
#'
#' \code{lfq.parquet} is already long (\code{filename}, \code{intensity});
#' \code{lfq.tsv} is wide with one column per raw file.
sage_lfq_long <- function(lfq) {
  if (all(c("filename", "intensity") %in% names(lfq))) {
    return(lfq)
  }
  sample_cols <- setdiff(names(lfq), SAGE_LFQ_ID_COLS)
  if (length(sample_cols) == 0) {
    stop("No per-sample intensity columns found in the Sage LFQ table.")
  }
  lfq |>
    tidyr::pivot_longer(
      dplyr::all_of(sample_cols),
      names_to = "filename",
      values_to = "intensity"
    )
}

#' Build a protein abundance matrix from the Sage LFQ table.
#'
#' Peptides are filtered (decoys, q-value, optionally shared peptides) and
#' rolled up per protein and raw file. Zero intensities are treated as
#' "not quantified" (NA), as in MaxQuant. The output has a \code{Protein ID}
#' column followed by one numeric column per raw file, which is the layout
#' expected by the PwrQuant module.
#'
#' @param lfq            data.frame read from lfq.parquet / lfq.tsv.
#' @param q_max          Maximum LFQ q-value for a peptide to be kept.
#' @param remove_shared  Drop peptides mapping to more than one protein.
#' @param method         One of \code{SAGE_ROLLUP_METHODS}.
#' @param min_peptides   Minimum distinct peptides per protein (across all
#'   files) for the protein to be reported.
#' @param log2_transform Apply log2 to the abundance columns.
#' @return A data.frame.
build_sage_protein_matrix <- function(
  lfq,
  q_max = 0.01,
  remove_shared = TRUE,
  method = "sum",
  min_peptides = 1L,
  log2_transform = FALSE
) {
  filter_sage_lfq(lfq, q_max = q_max, remove_shared = remove_shared) |>
    rollup_sage_proteins(
      method = method,
      min_peptides = min_peptides,
      log2_transform = log2_transform
    )
}

#' Filter the Sage LFQ table to quantifiable target peptides (long format).
#'
#' Removes decoys, peptides above \code{q_max}, optionally shared peptides,
#' and zero / missing intensities (Sage writes 0 for "not quantified").
filter_sage_lfq <- function(lfq, q_max = 0.01, remove_shared = TRUE) {
  d <- sage_lfq_long(lfq)

  if ("is_decoy" %in% names(d)) {
    d <- d |> dplyr::filter(!is_decoy)
  }
  if ("q_value" %in% names(d)) {
    d <- d |> dplyr::filter(is.na(q_value) | q_value <= q_max)
  }
  if (remove_shared) {
    d <- d |> dplyr::filter(!grepl(";", proteins, fixed = TRUE))
  }
  d <- d |>
    dplyr::mutate(intensity = as.numeric(intensity)) |>
    dplyr::filter(!is.na(intensity), intensity > 0, !is.na(proteins))

  if (nrow(d) == 0) {
    stop("No peptides left after filtering the Sage LFQ table.")
  }
  d
}

#' Roll filtered LFQ peptides up to a protein x sample matrix.
#'
#' @param d Output of \code{filter_sage_lfq()}.
#' @return data.frame with \code{Protein ID} followed by one column per file.
rollup_sage_proteins <- function(
  d,
  method = "sum",
  min_peptides = 1L,
  log2_transform = FALSE
) {
  method <- match.arg(method, unname(SAGE_ROLLUP_METHODS))
  rollup <- switch(
    method,
    sum = function(x) sum(x),
    median = function(x) stats::median(x),
    top3 = function(x) {
      mean(sort(x, decreasing = TRUE)[seq_len(min(3L, length(x)))])
    }
  )

  prot <- d |>
    dplyr::group_by(proteins, filename) |>
    dplyr::summarise(abundance = rollup(intensity), .groups = "drop")

  n_pep <- d |>
    dplyr::group_by(proteins) |>
    dplyr::summarise(n_peptides = dplyr::n_distinct(peptide), .groups = "drop")
  keep <- n_pep$proteins[n_pep$n_peptides >= min_peptides]

  out <- prot |>
    dplyr::filter(proteins %in% keep) |>
    tidyr::pivot_wider(
      names_from = filename,
      values_from = abundance
    ) |>
    dplyr::rename(`Protein ID` = proteins) |>
    dplyr::arrange(`Protein ID`) |>
    as.data.frame()

  if (log2_transform) {
    val_cols <- setdiff(names(out), "Protein ID")
    out[val_cols] <- lapply(out[val_cols], log2)
  }
  out
}

#' Distribution of distinct peptides per protein in the filtered LFQ table.
#'
#' Bars are coloured by whether the protein passes the \code{min_peptides}
#' threshold used for the abundance matrix.
plot_sage_peptides_per_protein <- function(
  d,
  min_peptides = 1L,
  color = "#1b9e77",
  cap = 20L
) {
  n_pep <- d |>
    dplyr::group_by(proteins) |>
    dplyr::summarise(n_peptides = dplyr::n_distinct(peptide), .groups = "drop")
  n_keep <- sum(n_pep$n_peptides >= min_peptides)

  n_pep |>
    dplyr::mutate(
      bin = factor(
        ifelse(n_peptides >= cap, paste0(cap, "+"), as.character(n_peptides)),
        levels = c(as.character(seq_len(cap - 1L)), paste0(cap, "+"))
      ),
      status = ifelse(
        n_peptides >= min_peptides,
        "In matrix",
        "Below min peptides"
      )
    ) |>
    ggplot(aes(x = bin, fill = status)) +
    geom_bar(color = "white", linewidth = 0.25) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(
      values = c("In matrix" = color, "Below min peptides" = "grey70"),
      name = NULL
    ) +
    labs(
      title = sprintf(
        "Peptides per protein (%s of %s proteins in matrix)",
        format(n_keep, big.mark = ","),
        format(nrow(n_pep), big.mark = ",")
      ),
      x = "Distinct peptides per protein",
      y = "Number of proteins"
    ) +
    theme_sage() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "top"
    )
}

#' Missingness summary of a protein abundance matrix.
#'
#' Top: percentage of proteins missing in each sample. Bottom: number of
#' samples in which each protein is quantified.
plot_sage_missingness <- function(pm, color = "#1b9e77") {
  val_cols <- setdiff(names(pm), "Protein ID")
  mat <- as.matrix(pm[val_cols])
  n_samp <- length(val_cols)

  per_sample <- data.frame(
    sample = factor(val_cols, levels = val_cols),
    pct_missing = colMeans(is.na(mat)) * 100
  )
  per_protein <- data.frame(n_quant = rowSums(!is.na(mat))) |>
    dplyr::count(n_quant, name = "n_proteins") |>
    dplyr::mutate(n_quant = factor(n_quant, levels = 0:n_samp))

  p1 <- ggplot(per_sample, aes(x = sample, y = pct_missing)) +
    geom_col(fill = color, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = sprintf("%.2f%%", pct_missing)),
      vjust = -0.3,
      size = 3,
      fontface = "bold"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = sprintf(
        "Overall missingness: %.2f%%",
        mean(is.na(mat)) * 100
      ),
      x = NULL,
      y = "Proteins missing (%)"
    ) +
    theme_sage()

  p2 <- ggplot(per_protein, aes(x = n_quant, y = n_proteins)) +
    geom_col(fill = color, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = n_proteins),
      vjust = -0.3,
      size = 3,
      fontface = "bold"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x = "Number of samples in which the protein is quantified",
      y = "Number of proteins"
    ) +
    theme_sage() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1.2, 1))
}

# ── iBAQ from the Sage LFQ table ──────────────────────────────────────────────
#   iBAQ = summed peptide intensity of a protein in a file /
#          number of theoretically observable peptides of that protein
# Same definition as in mod_QC4DIANN / mod_PSManalyst; the in-silico digest
# helpers live in modules/utils_fasta.r.

#' Digest settings of the search, from the results.json Sage writes next to
#' results.sage.parquet.
#'
#' @param results_path Path to results.sage.parquet / .tsv.
#' @return NULL when results.json is absent or unreadable, otherwise a list
#'   with \code{cut}, \code{nocut}, \code{sense} ("C"/"N"),
#'   \code{missed_cleavages}, \code{min_len}, \code{max_len} and \code{fasta}
#'   (the FASTA path used in the search); fields are NA when missing.
read_sage_digest_params <- function(results_path) {
  f <- file.path(dirname(results_path), "results.json")
  if (!file.exists(f)) {
    return(NULL)
  }
  cfg <- tryCatch(
    jsonlite::fromJSON(f, simplifyVector = FALSE),
    error = function(e) NULL
  )
  db <- cfg$database
  if (is.null(db)) {
    return(NULL)
  }
  enz <- db$enzyme
  chr <- function(x) if (is.null(x)) NA_character_ else as.character(x)[1]
  int <- function(x) if (is.null(x)) NA_integer_ else as.integer(x)[1]
  list(
    file = f,
    cut = chr(enz$cleave_at),
    nocut = chr(enz$restrict),
    # Sage cleaves C-terminally unless c_terminal is explicitly false
    sense = if (isFALSE(enz$c_terminal)) "N" else "C",
    missed_cleavages = int(enz$missed_cleavages),
    min_len = int(enz$min_len),
    max_len = int(enz$max_len),
    fasta = chr(db$fasta)
  )
}

#' Protein identifiers treated as contaminants (common FASTA prefixes).
SAGE_CONTAMINANT_REGEX <- "^(contam|cont|con)_"

#' iBAQ per protein and raw file, long format (untransformed).
#'
#' Sage does no protein inference, so only peptides mapping to a single
#' protein can be attributed: \code{d} must come from
#' \code{filter_sage_lfq(remove_shared = TRUE)}.
#'
#' @param d         Output of \code{filter_sage_lfq()}.
#' @param sequences Named sequences from \code{fasta_sequences_by_header_id()}.
#' @param genes     Optional gene names parallel to \code{sequences}.
#' @return tibble with protein_id, entry_name, gene, n_theoretical_peptides,
#'   sample, iBAQ. Proteins absent from the FASTA are dropped and listed in
#'   the "unmatched_proteins" attribute; excluded contaminants in
#'   "contaminant_proteins".
compute_ibaq_sage <- function(
  d,
  sequences,
  genes = NULL,
  exclude_contaminants = FALSE,
  contaminant_regex = SAGE_CONTAMINANT_REGEX,
  protease_regex = "(?<=[KR])(?!P)",
  min_peptide_len = 6,
  max_peptide_len = 30,
  max_missed_cleavages = 0
) {
  contaminants <- character()
  if (isTRUE(exclude_contaminants)) {
    is_contam <- grepl(contaminant_regex, d$proteins, ignore.case = TRUE)
    contaminants <- unique(d$proteins[is_contam])
    d <- d[!is_contam, , drop = FALSE]
  }
  if (nrow(d) == 0) {
    stop("No quantified peptides left to compute iBAQ.")
  }

  intensity <- d |>
    dplyr::group_by(sample = filename, protein = proteins) |>
    dplyr::summarise(sum_intensity = sum(intensity), .groups = "drop")

  protein <- unique(intensity$protein)
  proteins <- data.frame(
    protein = protein,
    protein_id = fasta_accessions(protein),
    entry_name = ifelse(
      grepl("^[^|]+\\|[^|]+\\|.+", protein),
      sub("^[^|]+\\|[^|]+\\|", "", protein),
      NA_character_
    )
  )
  idx <- match_fasta_sequences(protein, proteins$protein_id, sequences)
  found <- !is.na(idx)
  proteins$gene <- if (is.null(genes)) NA_character_ else genes[idx]
  proteins$n_theoretical_peptides <- NA_integer_
  proteins$n_theoretical_peptides[found] <- unname(count_theoretical_peptides(
    sequences[idx[found]],
    protease_regex = protease_regex,
    min_len = min_peptide_len,
    max_len = max_peptide_len,
    max_missed = max_missed_cleavages
  ))

  out <- intensity |>
    dplyr::inner_join(proteins[found, , drop = FALSE], by = "protein") |>
    # Proteins without any observable peptide would divide by zero
    dplyr::filter(
      !is.na(n_theoretical_peptides),
      n_theoretical_peptides > 0
    ) |>
    dplyr::mutate(iBAQ = sum_intensity / n_theoretical_peptides) |>
    dplyr::select(
      protein_id,
      entry_name,
      gene,
      n_theoretical_peptides,
      sample,
      iBAQ
    ) |>
    dplyr::arrange(protein_id, sample)
  attr(out, "unmatched_proteins") <- protein[!found]
  attr(out, "contaminant_proteins") <- contaminants
  out
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — SIDEBAR UI
# ═══════════════════════════════════════════════════════════════════════════════
Sage_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      style = "padding:12px 16px 4px;color:#ffffff;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      icon("leaf", lib = "font-awesome"),
      " Sage DDA/DIA"
    ),
    tags$div(
      style = "padding:0 8px;",
      textInput(
        ns("sage_folder"),
        "Path to Sage output folder",
        value = "",
        placeholder = "/path/to/sage_search"
      ),
      tags$p(
        style = "color:#adb5bd;font-size:11px;margin-top:-6px;",
        "results.sage.parquet (PSMs) and lfq.parquet (quantification) are ",
        "located automatically, searching subfolders. The .tsv versions ",
        "are used when no parquet is found."
      ),
      tags$div(
        style = "text-align:center;",
        actionButton(
          ns("load_files"),
          "Load Sage Files",
          class = "btn-primary",
          style = "width:80%;font-weight:bold;margin-bottom:6px;"
        )
      ),
      uiOutput(ns("file_status"))
    ),
    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),
    sliderInput(
      ns("lda_filter"),
      "Min Sage Discriminant Score (LDA)",
      min = -100,
      max = 100,
      value = -100,
      step = 0.5
    ),
    sliderInput(
      ns("qval_filter"),
      "Max Peptide q-value",
      min = 0,
      max = 0.05,
      value = 0.01,
      step = 0.005
    ),
    sliderInput(
      ns("rt_tolerance"),
      "ΔRT artifact threshold (min)",
      min = 0,
      max = 5,
      value = 0.5,
      step = 0.1
    ),
    checkboxInput(
      ns("filter_decoy"),
      "Filter out decoys for plots",
      value = TRUE
    ),

    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),
    colourpicker::colourInput(
      ns("color_target"),
      "Target colour",
      value = "#1b9e77"
    ),
    colourpicker::colourInput(
      ns("color_decoy"),
      "Decoy / Line colour",
      value = "#d95f02"
    ),
    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),
    tags$div(
      style = "padding:12px 16px 4px;color:#ffffff;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      "Protein Abundance Matrix (LFQ)"
    ),
    div(
      style = "padding:0 8px;",
      tags$p(
        style = "color:#adb5bd;font-size:11px;",
        "Built from lfq.parquet. Decoys are removed and zero intensities ",
        "are treated as missing. Output is compatible with PwrQuant."
      ),
      sliderInput(
        ns("pm_qval"),
        "Max LFQ peptide q-value",
        min = 0,
        max = 0.05,
        value = 0.01,
        step = 0.005
      ),
      selectInput(
        ns("pm_method"),
        "Protein roll-up",
        choices = SAGE_ROLLUP_METHODS,
        selected = "sum"
      ),
      numericInput(
        ns("pm_min_peptides"),
        "Min peptides per protein",
        value = 1,
        min = 1,
        max = 10,
        step = 1
      ),
      checkboxInput(
        ns("pm_remove_shared"),
        "Remove shared peptides (multi-protein)",
        value = TRUE
      ),
      checkboxInput(
        ns("pm_log2"),
        "log2-transform values",
        value = FALSE
      ),
      downloadButton(
        ns("download_protein_matrix"),
        "⬇ Protein Abundance Matrix (.tsv)",
        class = "dl-btn",
        style = "width:100%;text-align:left;"
      )
    )
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — BODY UI
# ═══════════════════════════════════════════════════════════════════════════════
Sage_body_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      id = ns("tabs"),
      type = "tabs",

      # ── Interactive Plot Viewer ──────────────────────────────────────────────
      tabPanel(
        "Interactive Plot Viewer",
        fluidRow(infoBoxOutput(ns("info_box"), width = 12)),
        fluidRow(
          box(
            title = "Dynamic Plot View",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            fluidRow(
              column(
                4,
                selectInput(
                  ns("plot_select"),
                  "Select Graphic",
                  choices = c(
                    "Number of PSMs" = "plot_psm_counts",
                    "Proteins & Peptides by File" = "plot_id_counts",
                    "Sage Discriminant Score (LDA)" = "plot_lda",
                    "Charge State Density" = "plot_charge",
                    "Peptide Length Density" = "plot_length",
                    "Missed Cleavages" = "plot_missed",
                    "GRAVY Index Distribution" = "plot_gravy",
                    "pI Distribution" = "plot_pi",
                    "RT vs Mass Error (Da)" = "plot_rt_error",
                    "Fragment Error (ppm)" = "plot_frag_error",
                    "RT vs Precursor Error (ppm)" = "plot_rt_precursor",
                    "Precursor Mass Error Density (ppm)" = "plot_precursor_error",
                    "Peptide vs Protein q-value" = "plot_pep_prot_qval",
                    "Peptide vs Spectrum q-value" = "plot_qvals",
                    "Peptide Yield vs. FDR" = "plot_fdr_curve"
                  )
                )
              ),
              column(
                4,
                tags$div(
                  style = "margin-top:25px;",
                  actionButton(
                    ns("run_plot"),
                    "Plot Selected Graphic",
                    icon = icon("chart-bar"),
                    class = "btn-primary",
                    style = "width:100%;"
                  )
                )
              ),
              column(
                4,
                tags$div(
                  style = "margin-top:25px;",
                  downloadButton(
                    ns("download_plot"),
                    "Download Plot (.png)",
                    style = "width:100%;"
                  )
                )
              )
            ),
            sp_wrap_sage(ns("spi_main"), uiOutput(ns("dynamic_plot_ui")))
          )
        )
      ),

      # ── Modification Diagnostic ─────────────────────────────────────────
      tabPanel(
        "Modification Diagnostic",
        fluidRow(
          box(
            title = "RT Shift Profile — Modified vs Unmodified Peptides",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(
                4,
                selectInput(
                  ns("mod_diag_sample"),
                  "Select Sample",
                  choices = c("All samples" = "__all__"),
                  selected = "__all__"
                )
              ),
              column(
                4,
                numericInput(
                  ns("mod_diag_top_n"),
                  "Top N Peptides",
                  value = 30,
                  min = 5,
                  max = 100,
                  step = 5
                )
              ),
              column(
                4,
                tags$div(
                  style = "margin-top:25px;display:flex;gap:6px;",
                  downloadButton(
                    ns("download_mod_diag_plot"),
                    "Plot (.png)",
                    style = "flex:1;"
                  ),
                  downloadButton(
                    ns("download_mod_diag_table"),
                    "Table (.tsv)",
                    style = "flex:1;"
                  )
                )
              )
            ),
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp_moddiag"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("mod_diag_plot_ui"))
            )
          )
        ),
        fluidRow(
          box(
            title = "Modification Diagnostic Table",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            DT::dataTableOutput(ns("mod_diag_table"))
          )
        )
      ),

      # ── Protein Matrix (LFQ) ────────────────────────────────────────────
      tabPanel(
        "Protein Matrix",
        fluidRow(
          box(
            title = "Protein Matrix Figure",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            fluidRow(
              column(
                4,
                selectInput(
                  ns("pm_plot_select"),
                  "Select Figure",
                  choices = c(
                    "Peptides per protein" = "peptides",
                    "Missingness summary" = "missing"
                  )
                )
              ),
              column(
                4,
                tags$div(
                  style = "margin-top:25px;",
                  actionButton(
                    ns("pm_build"),
                    "Build Plot",
                    icon = icon("chart-bar"),
                    class = "btn-primary",
                    style = "width:100%;"
                  )
                )
              ),
              column(
                4,
                tags$div(
                  style = "margin-top:25px;",
                  downloadButton(
                    ns("pm_download_plot"),
                    "Download Plot (.png)",
                    style = "width:100%;"
                  )
                )
              )
            ),
            plotOutput(ns("pm_plot"), height = "600px")
          )
        ),
        fluidRow(
          box(
            title = "Protein Abundance Matrix — Preview",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            tags$p(
              style = "color:#6c757d;font-size:12px;",
              "Settings are taken from the 'Protein Abundance Matrix (LFQ)' ",
              "section of the sidebar. The download contains the full matrix."
            ),
            DT::dataTableOutput(ns("pm_table"))
          )
        )
      ),

      # ── iBAQ (from the LFQ table) ───────────────────────────────────────
      tabPanel(
        "iBAQ",
        fluidRow(uiOutput(ns("ibaq_banner"))),
        fluidRow(box(
          title = "iBAQ Controls",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          fluidRow(
            column(
              4,
              fileInput(
                ns("ibaq_fasta"),
                "Protein FASTA (search database)",
                accept = c(".fasta", ".fa", ".fas")
              )
            ),
            column(
              2,
              numericInput(
                ns("ibaq_top_n"),
                "Top N proteins (bar chart)",
                value = 20,
                min = 1,
                max = 100,
                step = 1
              )
            ),
            column(
              2,
              numericInput(
                ns("ibaq_rank_labels"),
                "Labels per facet (ranked)",
                value = 5,
                min = 0,
                max = 50,
                step = 1
              )
            ),
            column(
              2,
              style = "margin-top:25px;",
              checkboxInput(
                ns("ibaq_no_contam"),
                "Exclude contaminants",
                value = FALSE
              )
            ),
            column(
              2,
              style = "margin-top:25px;",
              checkboxInput(
                ns("ibaq_log2"),
                "log₂(iBAQ + 1) in table / export",
                value = TRUE
              )
            )
          ),
          fluidRow(
            column(
              2,
              numericInput(
                ns("ibaq_min_len"),
                "Min peptide length",
                value = 6,
                min = 1,
                max = 100,
                step = 1
              )
            ),
            column(
              2,
              numericInput(
                ns("ibaq_max_len"),
                "Max peptide length",
                value = 30,
                min = 1,
                max = 100,
                step = 1
              )
            ),
            column(
              2,
              numericInput(
                ns("ibaq_missed_cleavages"),
                "Missed cleavages",
                value = 0,
                min = 0,
                max = 5,
                step = 1
              )
            ),
            column(
              2,
              textInput(ns("ibaq_cut"), "Cleaved residues", value = "KR")
            ),
            column(
              2,
              textInput(ns("ibaq_nocut"), "Blocked by residue", value = "P")
            ),
            column(
              2,
              selectInput(
                ns("ibaq_sense"),
                "Cleavage side",
                choices = c("C-terminal" = "C", "N-terminal" = "N")
              )
            )
          ),
          uiOutput(ns("ibaq_digest_source")),
          actionButton(
            ns("run_ibaq"),
            "Compute iBAQ",
            icon = icon("calculator"),
            class = "btn-primary",
            style = "margin-right:8px;margin-bottom:8px;"
          ),
          downloadButton(
            ns("download_ibaq_table"),
            tagList(icon("download"), " iBAQ Matrix (.tsv)"),
            class = "dl-btn",
            style = "margin-right:8px;margin-bottom:8px;"
          ),
          downloadButton(
            ns("download_all_ibaq_plots"),
            tagList(icon("file-zipper"), " Download All Plots (.zip)"),
            class = "dl-btn",
            style = "margin-bottom:8px;"
          ),
          uiOutput(ns("ibaq_summary"))
        )),
        fluidRow(box(
          title = "iBAQ Distribution per Sample",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          plotOutput(ns("ibaq_dist_plot"), height = "500px")
        )),
        fluidRow(box(
          title = "Dynamic Range — Ranked Protein Abundance",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          uiOutput(ns("ibaq_rank_plot_ui"))
        )),
        fluidRow(box(
          title = "Most Abundant Proteins — Relative iBAQ",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          uiOutput(ns("ibaq_top_plot_ui"))
        )),
        fluidRow(box(
          title = "iBAQ Table",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          DT::dataTableOutput(ns("ibaq_table"))
        ))
      )
    )
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — SERVER
# ═══════════════════════════════════════════════════════════════════════════════
Sage_server <- function(id, fasta_digest) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    spin_ids <- paste0("spi", sprintf("%02d", 1:14))
    show_all <- function() lapply(spin_ids, function(s) shinyjs::show(id = s))
    hide_sp <- function(s) shinyjs::hide(id = s)
    rh <- function(fn, sid) {
      on.exit(hide_sp(sid), add = TRUE)
      fn()
    }

    # ── Dynamic plot height helper ─────────────────────────────────────────────
    # ncol=3 for most faceted plots; precursor_error uses ncol=1
    sage_plot_h <- reactive({
      d <- filtered_data()
      n <- dplyr::n_distinct(d$filename)
      max(400L, ceiling(n / 3L) * 350L)
    })
    sage_plot_h1 <- reactive({
      d <- filtered_data()
      n <- dplyr::n_distinct(d$filename)
      max(400L, n * 350L) # ncol=1 for precursor_error
    })

    # ── Load data from a Sage output folder ───────────────────────────────────
    sage_data_rv <- reactiveVal(NULL)
    lfq_data_rv <- reactiveVal(NULL)
    sage_files_rv <- reactiveVal(NULL)
    status_log <- reactiveVal(character(0))

    log_step <- function(msg, level = "info", notify = TRUE) {
      icon_chr <- switch(
        level,
        ok = "\u2714",
        warn = "\u26A0",
        error = "\u2716",
        "\u2139"
      )
      status_log(c(status_log(), paste(icon_chr, msg)))
      if (notify && level %in% c("warn", "error")) {
        showNotification(
          msg,
          type = if (level == "error") "error" else "warning"
        )
      }
    }

    sage_data <- reactive({
      req(sage_data_rv())
      sage_data_rv()
    })
    lfq_data <- reactive({
      req(lfq_data_rv())
      lfq_data_rv()
    })

    observeEvent(input$load_files, {
      sage_data_rv(NULL)
      lfq_data_rv(NULL)
      status_log(character(0))

      folder <- trimws(input$sage_folder)
      if (!nzchar(folder)) {
        log_step("Please enter the path to the Sage output folder.", "error")
        return()
      }
      if (!dir.exists(path.expand(folder))) {
        log_step(paste0("Folder not found: ", folder), "error")
        return()
      }

      files <- find_sage_files(folder)
      sage_files_rv(files)
      if (is.na(files$results)) {
        log_step(
          "results.sage.parquet (or .tsv) not found under this folder.",
          "error"
        )
        return()
      }
      log_step(paste0("PSMs: ", basename(files$results)), "ok", notify = FALSE)
      if (is.na(files$lfq)) {
        log_step(
          "lfq.parquet not found; protein matrix export unavailable.",
          "warn"
        )
      } else {
        log_step(paste0("LFQ: ", basename(files$lfq)), "ok", notify = FALSE)
      }

      show_all()
      withProgress(message = "Loading Sage results", value = 0, {
        incProgress(0.1, detail = basename(files$results))
        res <- tryCatch(
          {
            df <- read_sage_table(files$results)
            if (
              !"stripped_peptide" %in% names(df) && "peptide" %in% names(df)
            ) {
              df$stripped_peptide <- str_replace_all(
                df$peptide,
                "\\[.*?\\]|\\.|\\-",
                ""
              )
            }
            if (!"filename" %in% names(df)) {
              df$filename <- "Unknown"
            }
            df
          },
          error = function(e) {
            log_step(paste("Error reading Sage results:", e$message), "error")
            NULL
          }
        )
        sage_data_rv(res)
        if (!is.null(res)) {
          log_step(
            sprintf(
              "%s PSMs across %d files loaded.",
              format(nrow(res), big.mark = ","),
              dplyr::n_distinct(res$filename)
            ),
            "ok",
            notify = FALSE
          )
        }

        if (!is.na(files$lfq)) {
          incProgress(0.6, detail = basename(files$lfq))
          lfq <- tryCatch(
            read_sage_table(files$lfq),
            error = function(e) {
              log_step(paste("Error reading Sage LFQ:", e$message), "error")
              NULL
            }
          )
          lfq_data_rv(lfq)
          if (!is.null(lfq)) {
            log_step(
              sprintf(
                "%s LFQ peptide rows loaded.",
                format(nrow(lfq), big.mark = ",")
              ),
              "ok",
              notify = FALSE
            )
          }
        }
        setProgress(1)
      })
    })

    output$file_status <- renderUI({
      lines <- status_log()
      if (length(lines) == 0) {
        return(tags$p(
          style = "color:#adb5bd;font-size:11px;",
          "No folder loaded yet."
        ))
      }
      tags$div(
        style = "color:#adb5bd;font-size:11px;line-height:1.4;",
        lapply(lines, function(l) tags$div(l))
      )
    })

    observe({
      d <- sage_data()
      req(d, "sage_discriminant_score" %in% names(d))
      sc <- d$sage_discriminant_score[is.finite(d$sage_discriminant_score)]
      if (length(sc)) {
        sc_lo <- floor(min(sc) * 10) / 10
        sc_hi <- ceiling(max(sc) * 10) / 10
        if (input$lda_filter == -100) {
          updateSliderInput(
            session,
            "lda_filter",
            min = sc_lo,
            max = sc_hi,
            value = sc_lo
          )
        } else {
          updateSliderInput(session, "lda_filter", min = sc_lo, max = sc_hi)
        }
      }
    })

    # ── Filtered data based on decoy checkbox and filters
    filtered_data <- reactive({
      req(sage_data())
      d <- sage_data()
      if ("peptide_q" %in% names(d)) {
        d <- d |>
          dplyr::filter(is.na(peptide_q) | peptide_q <= input$qval_filter)
      }

      if ("sage_discriminant_score" %in% names(d)) {
        d <- d |>
          dplyr::filter(
            is.na(sage_discriminant_score) |
              sage_discriminant_score >= input$lda_filter
          )
      }

      if (input$filter_decoy && "is_decoy" %in% names(d)) {
        d <- d |> dplyr::filter(is_decoy == FALSE)
      }

      if ("stripped_peptide" %in% names(d)) {
        safe_seqs <- str_remove_all(
          d$stripped_peptide,
          "[^ACDEFGHIKLMNPQRSTVWY]"
        )
        d$gravy <- sapply(safe_seqs, GRAVY)
        d$pI <- calculate_pI(safe_seqs)
        d$MW <- calculate_MW(safe_seqs)
      }
      d
    })

    fasta_digest <- reactive({
      req(input$fasta_file)
      showNotification(
        "Reading and digesting FASTA...",
        id = "fasta_notif_sage",
        duration = NULL
      )
      seqs <- read_fasta_custom(input$fasta_file$datapath)
      df <- in_silico_digest(seqs, max_missed = input$missed_cleavages)
      removeNotification("fasta_notif_sage")
      df
    })

    mapped_data <- reactive({
      d <- filtered_data()
      if (isTruthy(input$fasta_file)) {
        dig <- fasta_digest()
        classes <- classify_peptides(d$stripped_peptide, dig)
        d <- dplyr::left_join(
          d,
          classes,
          by = c("stripped_peptide" = "peptide")
        )
      } else {
        d$classification <- "Unmapped"
        d$mapped_proteins <- NA_character_
      }
      d
    })

    output$info_box <- renderInfoBox({
      req(sage_data())
      total <- nrow(sage_data())
      filtered <- nrow(filtered_data())
      pct <- if (total > 0) round(filtered / total * 100, 1) else 0
      infoBox(
        "Retained PSMs",
        paste0(filtered, " / ", total, " (", pct, "%)"),
        icon = icon("leaf", lib = "font-awesome"),
        color = "green"
      )
    })

    vals_decoy <- reactive({
      c("TRUE" = input$color_decoy, "FALSE" = input$color_target)
    })

    plot_psm_counts_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      if (!"is_decoy" %in% names(d)) {
        d$is_decoy <- FALSE
      }

      d |>
        dplyr::count(filename, is_decoy, name = "n_psm") |>
        ggplot(aes(x = filename, y = n_psm, fill = as.character(is_decoy))) +
        geom_col(position = position_dodge(preserve = "single")) +
        geom_text(
          aes(label = n_psm),
          position = position_dodge(width = 0.9, preserve = "single"),
          vjust = -0.3,
          size = 3,
          fontface = "bold"
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(x = "File", y = "Number of PSMs", fill = "Is decoy?") +
        scale_fill_manual(values = vals_decoy()) +
        theme_sage() +
        theme(legend.position = "top")
    })

    plot_id_counts_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      summ <- d |>
        dplyr::group_by(filename) |>
        dplyr::summarise(
          n_peptides = dplyr::n_distinct(stripped_peptide),
          n_proteins = if ("proteins" %in% names(d)) {
            dplyr::n_distinct(proteins)
          } else {
            0
          },
          .groups = "drop"
        )

      p <- ggplot(summ, aes(x = filename))
      if ("proteins" %in% names(d)) {
        p <- p +
          geom_col(
            aes(y = n_proteins),
            fill = input$color_target,
            position = "dodge"
          ) +
          geom_text(
            aes(y = n_proteins, label = n_proteins),
            vjust = 1.5,
            size = 3,
            fontface = "bold",
            color = "white"
          )
      }
      p +
        geom_line(
          aes(y = n_peptides, group = 1),
          color = input$color_decoy,
          linewidth = 1
        ) +
        geom_point(aes(y = n_peptides), color = input$color_decoy, size = 2) +
        geom_text(
          aes(y = n_peptides, label = n_peptides),
          vjust = -1,
          size = 3,
          fontface = "bold",
          color = input$color_decoy
        ) +
        labs(
          title = "Sage IDs Validation",
          x = "File",
          y = "Proteins (bar) / Peptides (line)"
        ) +
        theme_sage() +
        theme(axis.text.x = element_text(angle = 65, hjust = 1))
    })

    plot_lda_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "sage_discriminant_score" %in% names(d))
      if (!"is_decoy" %in% names(d)) {
        d$is_decoy <- FALSE
      }

      ggplot(
        d,
        aes(x = sage_discriminant_score, fill = as.character(is_decoy))
      ) +
        geom_density(alpha = 0.6, color = "white", linewidth = 0.25) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        labs(
          x = "Sage discriminant score (LDA)",
          y = "Density",
          fill = "Is Decoy?"
        ) +
        scale_fill_manual(values = vals_decoy()) +
        facet_wrap(~filename, ncol = 3) +
        theme_sage() +
        theme(legend.position = "top")
    })

    plot_annotated_density_faceted <- function(d, col_sym, title, xlab, color) {
      val <- d[[col_sym]]
      if (all(is.na(val))) {
        return(
          ggplot() +
            annotate("text", x = 0, y = 0, label = "Not enough data") +
            theme_void()
        )
      }
      m_df <- d |>
        dplyr::filter(!is.na(!!sym(col_sym))) |>
        dplyr::group_by(filename) |>
        dplyr::summarise(m = median(!!sym(col_sym)), .groups = "drop")

      ggplot(d, aes(x = !!sym(col_sym))) +
        geom_density(fill = color, color = "black", alpha = 0.6) +
        geom_vline(
          data = m_df,
          aes(xintercept = m),
          linetype = "dashed",
          color = "red",
          linewidth = 1
        ) +
        geom_text(
          data = m_df,
          aes(x = m, y = Inf, label = paste("Median:", round(m, 2))),
          vjust = 2,
          hjust = -0.1,
          color = "red",
          fontface = "bold"
        ) +
        labs(title = title, x = xlab, y = "Density") +
        facet_wrap(~filename) +
        theme_sage()
    }

    plot_charge_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "charge" %in% names(d))
      ggplot(d, aes(x = charge)) +
        geom_density(
          alpha = 0.6,
          fill = input$color_target,
          color = "black",
          linewidth = 0.25
        ) +
        labs(x = "Charge state", y = "Density") +
        scale_x_continuous(
          breaks = seq(1, max(6, max(d$charge, na.rm = TRUE)))
        ) +
        facet_wrap(~filename) +
        theme_sage()
    })

    plot_length_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      ggplot(d, aes(x = nchar(stripped_peptide))) +
        geom_density(
          alpha = 0.6,
          fill = input$color_target,
          color = "black",
          linewidth = 0.25
        ) +
        labs(x = "Peptide length (AA)", y = "Density") +
        facet_wrap(~filename) +
        theme_sage()
    })

    plot_missed_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "missed_cleavages" %in% names(d))
      ggplot(d, aes(x = factor(missed_cleavages))) +
        geom_bar(
          alpha = 0.8,
          fill = input$color_target,
          color = "black",
          linewidth = 0.25
        ) +
        labs(x = "Number of missed cleavages", y = "Count") +
        facet_wrap(~filename) +
        theme_sage() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    })

    plot_gravy_obj <- reactive({
      d <- mapped_data()
      req(nrow(d) > 0, "gravy" %in% names(d))
      plot_annotated_density_faceted(
        d,
        "gravy",
        "GRAVY Index Distribution",
        "GRAVY Index",
        input$color_target
      )
    })

    plot_pi_obj <- reactive({
      d <- mapped_data()
      req(nrow(d) > 0, "pI" %in% names(d))
      if (all(is.na(d$pI))) {
        return(
          ggplot() +
            annotate(
              "text",
              x = 0,
              y = 0,
              label = "pI not available (Peptides package missing?)"
            ) +
            theme_void()
        )
      }
      plot_2d_gel(
        d,
        pi_col = "pI",
        mw_col = "MW",
        title = "Virtual 2D Gel",
        facet_col = "filename"
      )
    })

    plot_rt_error_obj <- reactive({
      d <- filtered_data()
      req(
        nrow(d) > 0,
        "rt" %in% names(d),
        "expmass" %in% names(d),
        "calcmass" %in% names(d)
      )
      ggplot(d, aes(x = rt, y = expmass - calcmass)) +
        geom_density2d_filled(show.legend = FALSE) +
        labs(x = "Retention time (min)", y = "Mass error (Da)") +
        facet_wrap(~filename, ncol = 3) +
        theme_sage()
    })

    plot_frag_error_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "fragment_ppm" %in% names(d))
      ggplot(d, aes(x = fragment_ppm)) +
        geom_histogram(
          binwidth = 1,
          fill = input$color_target,
          alpha = 0.8,
          color = "black",
          linewidth = 0.25
        ) +
        labs(x = "Fragment error (ppm)", y = "Count") +
        facet_wrap(~filename, ncol = 3) +
        theme_sage()
    })

    plot_rt_precursor_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "rt" %in% names(d), "precursor_ppm" %in% names(d))
      ggplot(d, aes(x = rt, y = precursor_ppm)) +
        geom_density2d_filled(show.legend = FALSE) +
        labs(x = "Retention time (min)", y = "Precursor mass error (ppm)") +
        facet_wrap(~filename, ncol = 3) +
        theme_sage()
    })

    plot_precursor_error_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "expmass" %in% names(d), "calcmass" %in% names(d))
      d$pm_err <- (d$expmass - d$calcmass) / d$calcmass * 1e6
      lims <- unname(quantile(d$pm_err, probs = c(0.01, 0.99), na.rm = TRUE))
      if (any(is.na(lims))) {
        lims <- c(-50, 50)
      }

      ggplot(d, aes(x = pm_err)) +
        geom_density(
          fill = input$color_target,
          color = "black",
          alpha = 0.6,
          linewidth = 0.25
        ) +
        coord_cartesian(xlim = lims) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        labs(x = "Precursor mass error (ppm)", y = "Density") +
        facet_wrap(~filename, ncol = 1) +
        theme_sage()
    })

    plot_pep_prot_qval_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "peptide_q" %in% names(d), "protein_q" %in% names(d))
      if (!"is_decoy" %in% names(d)) {
        d$is_decoy <- FALSE
      }
      ggplot(
        d,
        aes(
          x = -log10(peptide_q + 1e-10),
          y = -log10(protein_q + 1e-10),
          color = as.character(is_decoy)
        )
      ) +
        geom_point(alpha = 0.4, size = 1.5) +
        labs(
          x = "Peptide-level q-value (-log10)",
          y = "Protein-level q-value (-log10)",
          color = "Is Decoy?"
        ) +
        scale_color_manual(values = vals_decoy()) +
        facet_wrap(~filename, ncol = 3) +
        theme_sage() +
        theme(legend.position = "top")
    })

    plot_qvals_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "peptide_q" %in% names(d), "spectrum_q" %in% names(d))
      if (!"is_decoy" %in% names(d)) {
        d$is_decoy <- FALSE
      }
      ggplot(
        d,
        aes(
          x = -log10(peptide_q + 1e-10),
          y = -log10(spectrum_q + 1e-10),
          color = as.character(is_decoy)
        )
      ) +
        geom_point(alpha = 0.4, size = 1.5) +
        labs(
          x = "Peptide-level q-value (-log10)",
          y = "Spectrum-level q-value (-log10)",
          color = "Is Decoy?"
        ) +
        scale_color_manual(values = vals_decoy()) +
        facet_wrap(~filename, ncol = 3) +
        theme_sage() +
        theme(legend.position = "top")
    })

    plot_fdr_curve_obj <- reactive({
      d <- sage_data()
      req(d, nrow(d) > 0, "peptide_q" %in% names(d))
      td <- d |> dplyr::filter(!is.na(peptide_q), !is.na(filename))
      if (input$filter_decoy && "is_decoy" %in% names(td)) {
        td <- td |> dplyr::filter(is_decoy == FALSE)
      }
      q_sorted <- td |>
        dplyr::group_by(filename) |>
        dplyr::arrange(peptide_q, .by_group = TRUE) |>
        dplyr::mutate(cumulative_peptides = dplyr::row_number()) |>
        dplyr::ungroup()

      ggplot(q_sorted, aes(x = peptide_q, y = cumulative_peptides)) +
        geom_line(color = input$color_target, linewidth = 1) +
        geom_vline(
          xintercept = 0.01,
          linetype = "dashed",
          color = "red",
          linewidth = 0.7
        ) +
        geom_vline(
          xintercept = 0.05,
          linetype = "dashed",
          color = "orange",
          linewidth = 0.7
        ) +
        annotate(
          "text",
          x = 0.01,
          y = Inf,
          label = "1% FDR",
          vjust = 2,
          hjust = -0.1,
          color = "red",
          fontface = "bold",
          size = 3.5
        ) +
        annotate(
          "text",
          x = 0.05,
          y = Inf,
          label = "5% FDR",
          vjust = 2,
          hjust = -0.1,
          color = "orange",
          fontface = "bold",
          size = 3.5
        ) +
        scale_x_continuous(labels = scales::label_scientific()) +
        labs(
          title = "Peptide Yield vs. FDR",
          x = "Peptide q-value",
          y = "Cumulative peptide count"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_sage()
    })

    # ── renderUI wrappers for dynamic height ──────────────────────────────────
    output$dynamic_plot_ui <- renderUI({
      req(filtered_data())
      h <- if (input$plot_select == "plot_precursor_error") {
        sage_plot_h1()
      } else {
        sage_plot_h()
      }
      plotOutput(ns("dynamic_plot_out"), height = paste0(h, "px"))
    })

    current_plot_obj <- eventReactive(input$run_plot, {
      req(input$plot_select)
      shinyjs::show(id = "spi_main")

      switch(
        input$plot_select,
        "plot_psm_counts" = plot_psm_counts_obj(),
        "plot_id_counts" = plot_id_counts_obj(),
        "plot_lda" = plot_lda_obj(),
        "plot_charge" = plot_charge_obj(),
        "plot_length" = plot_length_obj(),
        "plot_missed" = plot_missed_obj(),
        "plot_gravy" = plot_gravy_obj(),
        "plot_pi" = plot_pi_obj(),
        "plot_rt_error" = plot_rt_error_obj(),
        "plot_frag_error" = plot_frag_error_obj(),
        "plot_rt_precursor" = plot_rt_precursor_obj(),
        "plot_precursor_error" = plot_precursor_error_obj(),
        "plot_pep_prot_qval" = plot_pep_prot_qval_obj(),
        "plot_qvals" = plot_qvals_obj(),
        "plot_fdr_curve" = plot_fdr_curve_obj()
      )
    })

    output$dynamic_plot_out <- renderPlot({
      on.exit(hide_sp("spi_main"), add = TRUE)
      req(current_plot_obj())
      current_plot_obj()
    })

    # ── Download Plots
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("Sage_", input$plot_select, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(current_plot_obj())
        ggsave(
          file,
          current_plot_obj(),
          width = 11,
          height = 8,
          bg = "white",
          device = "png"
        )
      }
    )

    # ── Protein abundance matrix (PwrQuant-compatible) ────────────────────────
    pm_min_peptides <- reactive({
      v <- suppressWarnings(as.integer(input$pm_min_peptides))
      if (is.na(v)) 1L else max(1L, v)
    })

    # Sidebar inputs only exist while the Sage sidebar is rendered; fall back
    # to the UI defaults otherwise.
    or_default <- function(x, default) if (is.null(x)) default else x

    # NB: shiny::validate / shiny::need are namespaced on purpose -
    # jsonlite::validate() masks shiny's in proteOmni.
    lfq_filtered <- reactive({
      shiny::validate(shiny::need(
        !is.null(lfq_data_rv()),
        "Load a Sage output folder containing lfq.parquet."
      ))
      filter_sage_lfq(
        lfq_data(),
        q_max = or_default(input$pm_qval, 0.01),
        remove_shared = isTRUE(or_default(input$pm_remove_shared, TRUE))
      )
    })

    protein_matrix <- reactive({
      rollup_sage_proteins(
        lfq_filtered(),
        method = or_default(input$pm_method, "sum"),
        min_peptides = pm_min_peptides(),
        log2_transform = isTRUE(input$pm_log2)
      )
    })

    pm_plot_obj <- eventReactive(input$pm_build, {
      col <- or_default(input$color_target, "#1b9e77")
      if (identical(input$pm_plot_select, "missing")) {
        pm <- protein_matrix()
        shiny::validate(shiny::need(
          nrow(pm) > 0,
          "No proteins pass the current filters."
        ))
        plot_sage_missingness(pm, color = col)
      } else {
        plot_sage_peptides_per_protein(
          lfq_filtered(),
          min_peptides = pm_min_peptides(),
          color = col
        )
      }
    })

    output$pm_plot <- renderPlot({
      shiny::validate(shiny::need(
        input$pm_build > 0,
        "Select a figure and click 'Build Plot'."
      ))
      pm_plot_obj()
    })

    output$pm_download_plot <- downloadHandler(
      filename = function() {
        paste0(
          "Sage_protein_matrix_",
          input$pm_plot_select,
          "_",
          Sys.Date(),
          ".png"
        )
      },
      content = function(file) {
        req(input$pm_build > 0)
        ggsave(
          file,
          pm_plot_obj(),
          width = 11,
          height = if (identical(input$pm_plot_select, "missing")) 10 else 7,
          bg = "white",
          device = "png"
        )
      }
    )

    output$pm_table <- DT::renderDataTable({
      pm <- protein_matrix()
      val_cols <- setdiff(names(pm), "Protein ID")
      DT::datatable(
        head(pm, 500),
        rownames = FALSE,
        caption = sprintf(
          "%s proteins x %d samples (%s%s, %.1f%% missing) - first 500 rows shown",
          format(nrow(pm), big.mark = ","),
          length(val_cols),
          names(SAGE_ROLLUP_METHODS)[SAGE_ROLLUP_METHODS == input$pm_method],
          if (isTRUE(input$pm_log2)) ", log2" else "",
          mean(is.na(as.matrix(pm[val_cols]))) * 100
        ),
        options = list(dom = "frtip", pageLength = 20, scrollX = TRUE),
        class = "display compact"
      ) |>
        DT::formatRound(
          columns = val_cols,
          digits = if (isTRUE(input$pm_log2)) 3 else 0
        )
    })

    output$download_protein_matrix <- downloadHandler(
      filename = function() {
        tag <- paste0("sage_", input$pm_method)
        if (isTRUE(input$pm_log2)) {
          tag <- paste0(tag, "_log2")
        }
        paste0("protein_abundance_matrix_", tag, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        if (is.null(lfq_data_rv())) {
          showNotification(
            "lfq.parquet is not loaded; load a Sage folder containing it first.",
            type = "error"
          )
          req(FALSE)
        }
        pm <- withProgress(
          message = "Building protein matrix...",
          value = 0.5,
          protein_matrix()
        )
        data.table::fwrite(pm, file = file, sep = "\t", na = "NA")
        log_step(
          sprintf(
            "Protein matrix downloaded: %s proteins x %d samples (%s%s).",
            format(nrow(pm), big.mark = ","),
            ncol(pm) - 1,
            names(SAGE_ROLLUP_METHODS)[SAGE_ROLLUP_METHODS == input$pm_method],
            if (isTRUE(input$pm_log2)) ", log2" else ""
          ),
          "ok",
          notify = FALSE
        )
      }
    )

    # ════════════════════════════════════════════════════════════════════════
    # iBAQ  (from lfq.parquet; compute_ibaq_sage is defined above)
    # ════════════════════════════════════════════════════════════════════════
    # Digest settings of the search (results.json) become the defaults of the
    # iBAQ controls whenever a folder is loaded; the user can still edit them.
    sage_digest <- reactive({
      files <- sage_files_rv()
      if (is.null(files) || is.na(files$results)) {
        return(NULL)
      }
      read_sage_digest_params(files$results)
    })

    observeEvent(sage_digest(), {
      p <- sage_digest()
      if (!is.na(p$min_len)) {
        updateNumericInput(session, "ibaq_min_len", value = p$min_len)
      }
      if (!is.na(p$max_len)) {
        updateNumericInput(session, "ibaq_max_len", value = p$max_len)
      }
      if (!is.na(p$missed_cleavages)) {
        updateNumericInput(
          session,
          "ibaq_missed_cleavages",
          value = p$missed_cleavages
        )
      }
      # Leave the enzyme controls alone when the rule has no cut residue
      # (non-enzymatic or no-digestion search): there is no digest to mirror.
      if (!is.null(protease_regex_from_rule(p$cut))) {
        updateTextInput(session, "ibaq_cut", value = p$cut)
        updateTextInput(
          session,
          "ibaq_nocut",
          value = if (is.na(p$nocut)) "" else p$nocut
        )
        updateSelectInput(session, "ibaq_sense", selected = p$sense)
      }
    })

    output$ibaq_digest_source <- renderUI({
      req(sage_files_rv())
      p <- sage_digest()
      msg <- if (is.null(p)) {
        "No results.json found next to the Sage results: digest settings are the app defaults."
      } else {
        paste0(
          "Digest defaults read from results.json. Edit the fields above to override them.",
          if (is.null(protease_regex_from_rule(p$cut))) {
            " The search has no cleavage rule (non-enzymatic), so the enzyme fields were left unchanged."
          }
        )
      }
      helpText(msg)
    })

    # FASTA for the in-silico digest: the upload in the iBAQ controls, otherwise
    # the search database recorded in results.json when it is still on disk.
    ibaq_fasta_source <- reactive({
      if (!is.null(input$ibaq_fasta)) {
        return(list(path = input$ibaq_fasta$datapath, source = "upload"))
      }
      p <- sage_digest()
      if (!is.null(p) && !is.na(p$fasta) && file.exists(p$fasta)) {
        return(list(path = p$fasta, source = "results.json"))
      }
      NULL
    })

    ibaq_sequences <- reactive({
      src <- ibaq_fasta_source()
      req(src)
      fa <- Biostrings::readAAStringSet(src$path)
      keep <- !duplicated(fasta_header_ids(names(fa)))
      list(
        seqs = fasta_sequences_by_header_id(fa),
        genes = fasta_header_genes(names(fa))[keep]
      )
    })

    output$ibaq_banner <- renderUI({
      src <- ibaq_fasta_source()
      msg <- if (is.null(lfq_data_rv())) {
        "Load a Sage output folder containing lfq.parquet to enable iBAQ calculation."
      } else if (is.null(src)) {
        p <- sage_digest()
        paste0(
          "Upload the protein FASTA used in the search (iBAQ Controls) to enable iBAQ calculation.",
          if (!is.null(p) && !is.na(p$fasta)) {
            paste0(" Sage searched: ", p$fasta)
          }
        )
      } else if (src$source == "results.json") {
        paste0(
          "Using the search database recorded in results.json: ",
          src$path,
          ". Upload a FASTA to use a different one."
        )
      }
      if (!is.null(msg)) {
        div(
          class = "plot-placeholder",
          icon("circle-info", style = "margin-right:10px;"),
          msg
        )
      }
    })

    # Untransformed iBAQ in long format. Uses the sidebar's LFQ peptide
    # q-value; decoys, unquantified and shared peptides are always removed
    # (Sage does no protein inference, so shared peptides cannot be assigned).
    ibaq_long <- eventReactive(input$run_ibaq, {
      if (is.null(lfq_data_rv()) || is.null(ibaq_fasta_source())) {
        showNotification(
          "Load a Sage folder with lfq.parquet and provide a protein FASTA to compute iBAQ.",
          type = "warning"
        )
        req(FALSE)
      }
      num <- function(x, default) if (isTruthy(x)) as.integer(x) else default
      min_len <- max(1L, num(input$ibaq_min_len, 6L))
      max_len <- max(min_len, num(input$ibaq_max_len, 30L))
      tryCatch(
        withProgress(message = "Computing iBAQ…", value = 0.3, {
          protease_regex <- protease_regex_from_rule(
            input$ibaq_cut,
            input$ibaq_nocut,
            input$ibaq_sense
          )
          if (is.null(protease_regex)) {
            stop("enter at least one cleaved residue (e.g. KR for trypsin).")
          }
          fa <- ibaq_sequences()
          d <- filter_sage_lfq(
            lfq_data(),
            q_max = or_default(input$pm_qval, 0.01),
            remove_shared = TRUE
          )
          compute_ibaq_sage(
            d,
            fa$seqs,
            genes = fa$genes,
            exclude_contaminants = isTRUE(input$ibaq_no_contam),
            protease_regex = protease_regex,
            min_peptide_len = min_len,
            max_peptide_len = max_len,
            max_missed_cleavages = min(
              5L,
              max(0L, num(input$ibaq_missed_cleavages, 0L))
            )
          )
        }),
        error = function(e) {
          showNotification(
            paste("iBAQ error:", conditionMessage(e)),
            type = "error",
            duration = 10
          )
          NULL
        }
      )
    })

    ibaq_top_n <- reactive({
      n <- suppressWarnings(as.integer(input$ibaq_top_n))
      if (length(n) == 0 || is.na(n)) 20L else max(1L, min(100L, n))
    })

    # Labels per facet in the ranked-abundance plot; 0 turns labels off
    ibaq_rank_labels <- reactive({
      n <- suppressWarnings(as.integer(input$ibaq_rank_labels))
      if (length(n) == 0 || is.na(n)) 5L else max(0L, min(50L, n))
    })

    # Per-sample rank and relative iBAQ (riBAQ = share of the sample's total
    # iBAQ), labelled by gene, falling back to entry name, then protein ID.
    ibaq_per_sample <- reactive({
      d <- ibaq_long()
      req(d, nrow(d) > 0)
      labels <- d[!duplicated(d$protein_id), ]
      lab <- dplyr::coalesce(
        dplyr::na_if(labels$gene, ""),
        dplyr::na_if(labels$entry_name, ""),
        labels$protein_id
      )
      d$label <- make.unique(lab)[match(d$protein_id, labels$protein_id)]
      d |>
        dplyr::group_by(sample) |>
        dplyr::mutate(
          riBAQ = iBAQ / sum(iBAQ) * 100,
          rank = dplyr::row_number(dplyr::desc(iBAQ))
        ) |>
        dplyr::ungroup()
    })

    # Facet layout and pixel heights for the per-sample plots
    IBAQ_NCOL <- 3L
    ibaq_rank_h <- reactive({
      n <- dplyr::n_distinct(ibaq_per_sample()$sample)
      max(500L, ceiling(n / IBAQ_NCOL) * 300L)
    })
    ibaq_top_h <- reactive({
      n <- dplyr::n_distinct(ibaq_per_sample()$sample)
      max(500L, ceiling(n / IBAQ_NCOL) * (ibaq_top_n() * 14L + 70L))
    })

    output$ibaq_rank_plot_ui <- renderUI({
      plotOutput(ns("ibaq_rank_plot"), height = paste0(ibaq_rank_h(), "px"))
    })
    output$ibaq_top_plot_ui <- renderUI({
      plotOutput(ns("ibaq_top_plot"), height = paste0(ibaq_top_h(), "px"))
    })

    ibaq_wide <- reactive({
      d <- ibaq_long()
      req(d, nrow(d) > 0)
      if (isTRUE(input$ibaq_log2)) {
        d$iBAQ <- log2(d$iBAQ + 1)
      }
      d |>
        dplyr::arrange(sample) |>
        tidyr::pivot_wider(names_from = sample, values_from = iBAQ) |>
        dplyr::arrange(protein_id)
    })

    output$ibaq_summary <- renderUI({
      d <- ibaq_long()
      req(d)
      n_unmatched <- length(attr(d, "unmatched_proteins"))
      n_contam <- length(attr(d, "contaminant_proteins"))
      tags$p(
        style = "margin:6px 0 0 0;",
        icon("circle-info"),
        sprintf(
          " iBAQ computed for %s proteins across %s samples, from peptides unique to one protein.",
          scales::comma(dplyr::n_distinct(d$protein_id)),
          dplyr::n_distinct(d$sample)
        ),
        if (n_contam > 0) {
          sprintf(
            " %s contaminant protein(s) were excluded.",
            scales::comma(n_contam)
          )
        },
        if (n_unmatched > 0) {
          sprintf(
            " %s protein(s) were excluded because they are not in the FASTA (a FASTA that differs from the search database).",
            scales::comma(n_unmatched)
          )
        }
      )
    })

    ibaq_color <- reactive(or_default(input$color_target, "#1b9e77"))

    ibaq_dist_obj <- reactive({
      d <- ibaq_long()
      req(d, nrow(d) > 0)
      ggplot(d, aes(x = sample, y = log2(iBAQ + 1))) +
        geom_boxplot(fill = ibaq_color(), alpha = 0.75) +
        labs(x = NULL, y = "log₂(iBAQ + 1)") +
        theme_sage()
    })

    ibaq_rank_obj <- reactive({
      s <- ibaq_per_sample()
      ggplot(s, aes(x = rank, y = log2(iBAQ + 1))) +
        geom_point(color = ibaq_color(), alpha = 0.6, size = 0.8) +
        ggrepel::geom_text_repel(
          data = dplyr::filter(s, rank <= ibaq_rank_labels()),
          aes(label = label),
          size = 3,
          max.overlaps = Inf,
          min.segment.length = 0
        ) +
        labs(x = "Protein rank", y = "log₂(iBAQ + 1)") +
        facet_wrap(~sample, ncol = IBAQ_NCOL) +
        theme_sage() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    })

    ibaq_top_obj <- reactive({
      top <- ibaq_per_sample() |>
        dplyr::filter(rank <= ibaq_top_n())
      ggplot(
        top,
        aes(y = tidytext::reorder_within(label, riBAQ, sample), x = riBAQ)
      ) +
        geom_col(fill = ibaq_color(), color = "black", alpha = 0.85) +
        tidytext::scale_y_reordered() +
        labs(y = NULL, x = "Relative iBAQ (% of sample total)") +
        facet_wrap(~sample, ncol = IBAQ_NCOL, scales = "free_y") +
        theme_sage() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    })

    output$ibaq_dist_plot <- renderPlot({
      ibaq_dist_obj()
    })
    output$ibaq_rank_plot <- renderPlot({
      ibaq_rank_obj()
    })
    output$ibaq_top_plot <- renderPlot({
      ibaq_top_obj()
    })

    output$ibaq_table <- DT::renderDataTable({
      tbl <- ibaq_wide()
      DT::datatable(
        tbl,
        rownames = FALSE,
        filter = "top",
        options = list(pageLength = 20, scrollX = TRUE)
      ) |>
        DT::formatRound(columns = unique(ibaq_long()$sample), digits = 2)
    })

    output$download_ibaq_table <- downloadHandler(
      filename = function() {
        paste0(
          "Sage_",
          if (isTRUE(input$ibaq_log2)) "log2_" else "",
          "iBAQ_matrix_",
          Sys.Date(),
          ".tsv"
        )
      },
      content = function(file) {
        data.table::fwrite(ibaq_wide(), file = file, sep = "\t", na = "NA")
      }
    )

    output$download_all_ibaq_plots <- downloadHandler(
      filename = function() {
        paste0("Sage_iBAQ_plots_", Sys.Date(), ".zip")
      },
      content = function(file) {
        td <- tempdir()
        fps <- character()
        save_p <- function(nm, fn, w = 10, h = 8) {
          fp <- file.path(td, nm)
          tryCatch(
            {
              ggsave(
                fp,
                fn(),
                width = w,
                height = h,
                dpi = 300,
                bg = "white",
                limitsize = FALSE
              )
              fps <<- c(fps, fp)
            },
            error = function(e) message("Skip ", nm, ": ", conditionMessage(e))
          )
        }

        withProgress(message = "Saving plots...", value = 0, {
          save_p("iBAQ_distribution.png", function() ibaq_dist_obj(), w = 12)
          incProgress(1 / 3)
          # Same proportions as on screen (px / 80), capped for ggsave
          save_p(
            "iBAQ_ranked_abundance.png",
            function() ibaq_rank_obj(),
            w = 16,
            h = min(45, ibaq_rank_h() / 80)
          )
          incProgress(1 / 3)
          save_p(
            "iBAQ_top_proteins.png",
            function() ibaq_top_obj(),
            w = 16,
            h = min(45, ibaq_top_h() / 80)
          )
          incProgress(1 / 3)
        })

        req(length(fps) > 0)
        utils::zip(file, files = fps, flags = "-j")
      },
      contentType = "application/zip"
    )

    # ════════════════════════════════════════════════════════════════════════
    # MODIFICATION DIAGNOSTIC
    # ════════════════════════════════════════════════════════════════════════
    mod_diag_data <- reactive({
      d <- filtered_data()
      req(d, nrow(d) > 0)
      req(all(c("peptide", "stripped_peptide", "rt", "filename") %in% names(d)))

      d_mod <- d |>
        dplyr::filter(
          !is.na(peptide),
          !is.na(stripped_peptide),
          peptide != stripped_peptide,
          !is.na(rt)
        ) |>
        dplyr::select(stripped_peptide, peptide, rt, filename) |>
        dplyr::rename(
          peptide_seq = stripped_peptide,
          mod_seq = peptide,
          rt_mod = rt,
          sample_name = filename
        )

      d_unmod <- d |>
        dplyr::filter(!is.na(peptide), peptide == stripped_peptide) |>
        dplyr::group_by(stripped_peptide, filename) |>
        dplyr::summarise(
          rt_unmod = median(rt, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(peptide_seq = stripped_peptide, sample_name = filename)

      dplyr::inner_join(d_mod, d_unmod, by = c("peptide_seq", "sample_name")) |>
        dplyr::mutate(
          delta_rt = abs(rt_mod - rt_unmod),
          classification = dplyr::if_else(
            delta_rt <= input$rt_tolerance,
            "Suspected Artifact",
            "Sample-Derived"
          )
        )
    })

    # Populate the sample selector whenever new results are loaded
    observeEvent(sage_data_rv(), ignoreNULL = FALSE, {
      d <- sage_data_rv()
      samples <- if (is.null(d)) character(0) else sort(unique(d$filename))
      updateSelectInput(
        session,
        "mod_diag_sample",
        choices = c("All samples" = "__all__", setNames(samples, samples)),
        selected = "__all__"
      )
    })

    mod_diag_top_n <- reactive({
      v <- suppressWarnings(as.integer(input$mod_diag_top_n))
      if (length(v) == 0 || is.na(v)) 30L else max(1L, v)
    })

    # Paired data restricted to the selected sample
    mod_diag_selected <- reactive({
      paired <- mod_diag_data()
      sel <- input$mod_diag_sample
      if (!is.null(sel) && nzchar(sel) && sel != "__all__") {
        paired <- paired |> dplyr::filter(sample_name == sel)
      }
      paired
    })

    mod_diag_plot_obj <- reactive({
      paired <- mod_diag_selected()
      req(paired, nrow(paired) > 0)

      # Top N peptides by largest |dRT|, chosen within each sample
      top_pairs <- paired |>
        dplyr::group_by(sample_name, peptide_seq) |>
        dplyr::summarise(
          max_delta = max(delta_rt, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::group_by(sample_name) |>
        dplyr::slice_max(
          order_by = max_delta,
          n = mod_diag_top_n(),
          with_ties = FALSE
        ) |>
        dplyr::ungroup() |>
        dplyr::select(sample_name, peptide_seq)

      plot_data <- paired |>
        dplyr::semi_join(top_pairs, by = c("sample_name", "peptide_seq"))
      req(nrow(plot_data) > 0)
      n_facets <- dplyr::n_distinct(plot_data$sample_name)

      ggplot(plot_data, aes(y = reorder(peptide_seq, delta_rt))) +
        geom_segment(
          aes(
            x = rt_unmod,
            xend = rt_mod,
            yend = reorder(peptide_seq, delta_rt),
            color = classification
          ),
          linewidth = 0.7
        ) +
        geom_point(aes(x = rt_unmod), color = "grey50", size = 2) +
        geom_point(aes(x = rt_mod, color = classification), size = 2.5) +
        scale_color_manual(
          values = c(
            "Suspected Artifact" = "#e74c3c",
            "Sample-Derived" = "#2ecc71"
          ),
          name = "Classification"
        ) +
        labs(
          title = "RT shift profile \u2014 modified vs unmodified (Sage)",
          x = "Retention time (min)",
          y = "Peptide sequence (stripped)",
          caption = paste0(
            "\u0394RT threshold: ",
            input$rt_tolerance,
            " min | Artifact if |\u0394RT| \u2264 threshold"
          )
        ) +
        facet_wrap(
          ~sample_name,
          ncol = if (n_facets == 1) 1 else 2,
          scales = "free"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 7, face = "bold", color = "black"),
          axis.text.x = element_text(face = "bold", color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", face = "bold"),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "bottom",
          legend.text = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          plot.caption = element_text(size = 12, face = "bold")
        )
    })

    output$mod_diag_plot_ui <- renderUI({
      paired <- tryCatch(mod_diag_selected(), error = function(e) NULL)
      hide_sp("sp_moddiag")
      if (is.null(paired) || nrow(paired) == 0) {
        return(tags$p(
          style = "color:#adb5bd;text-align:center;padding:20px;",
          "No paired modified/unmodified peptides found. Load a Sage output folder and ensure results.sage.parquet contains 'peptide', 'stripped_peptide', 'rt', and 'filename' columns."
        ))
      }
      n_samples <- dplyr::n_distinct(paired$sample_name)
      n_rows <- if (n_samples == 1) 1L else ceiling(n_samples / 2)
      # ~16 px per peptide label keeps the y axis readable as Top N grows
      row_h <- max(400L, mod_diag_top_n() * 16L + 150L)
      dynamic_h <- min(n_rows * row_h, 6000L)
      plotOutput(ns("mod_diag_plot"), height = paste0(dynamic_h, "px"))
    })

    output$mod_diag_plot <- renderPlot({
      mod_diag_plot_obj()
    })

    # Table shown in the tab and exported by the TSV download
    mod_diag_table_data <- reactive({
      paired <- mod_diag_selected()
      req(paired, nrow(paired) > 0)
      paired |>
        dplyr::select(
          sample_name,
          peptide_seq,
          mod_seq,
          rt_unmod,
          rt_mod,
          delta_rt,
          classification
        ) |>
        dplyr::rename(peptide = peptide_seq, modified_peptide = mod_seq) |>
        dplyr::arrange(dplyr::desc(delta_rt))
    })

    mod_diag_file_tag <- function() {
      sel <- input$mod_diag_sample
      if (is.null(sel) || sel == "__all__") {
        "all_samples"
      } else {
        gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(sel))
      }
    }

    output$download_mod_diag_plot <- downloadHandler(
      filename = function() {
        paste0(
          "Sage_mod_diagnostic_",
          mod_diag_file_tag(),
          "_top",
          mod_diag_top_n(),
          "_",
          Sys.Date(),
          ".png"
        )
      },
      content = function(file) {
        p <- mod_diag_plot_obj()
        n_samples <- dplyr::n_distinct(p$data$sample_name)
        n_rows <- if (n_samples == 1) 1L else ceiling(n_samples / 2)
        # Same sizing rule as the on-screen plot, converted to inches
        row_in <- max(4, (mod_diag_top_n() * 16 + 150) / 100)
        ggsave(
          file,
          p,
          width = if (n_samples == 1) 9 else 14,
          height = min(n_rows * row_in, 49),
          bg = "white",
          device = "png",
          limitsize = FALSE
        )
      }
    )

    output$download_mod_diag_table <- downloadHandler(
      filename = function() {
        paste0(
          "Sage_mod_diagnostic_",
          mod_diag_file_tag(),
          "_",
          Sys.Date(),
          ".tsv"
        )
      },
      content = function(file) {
        data.table::fwrite(
          mod_diag_table_data(),
          file = file,
          sep = "\t",
          na = "NA"
        )
      }
    )

    output$mod_diag_table <- DT::renderDataTable({
      tbl <- mod_diag_table_data()
      DT::datatable(
        tbl,
        rownames = FALSE,
        filter = "top",
        options = list(pageLength = 20, scrollX = TRUE)
      ) |>
        DT::formatRound(
          columns = c("rt_unmod", "rt_mod", "delta_rt"),
          digits = 3
        )
    })
  })
}
