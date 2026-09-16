# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  proteOmni — MaxQuant MS/MS Spectrum + Evidence QC Module                  ║
# ║  File: mod_MaxQuantMSMS.r                                                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── 1. HELPER FUNCTIONS ─────────────────────────────────────────────────────
theme_mq <- function(...) {
  theme_bw(...) +
    theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text = element_text(face = "bold", color = "black", size = 12),
      axis.title = ggtext::element_markdown(
        size = 14,
        face = "bold",
        color = "black"
      ),
      legend.position = "bottom",
      legend.title = element_text(
        size = 12,
        face = "bold",
        color = "black",
        hjust = 0.5
      ),
      legend.title.position = "top",
      legend.key.height = unit(0.2, "line"),
      legend.key.width = unit(3, "line"),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold", size = 12),
      panel.border = element_rect(color = "black", fill = NA)
    )
}


# ── MaxQuant folder discovery ─────────────────────────────────────────────────

#' Locate the MaxQuant output files under a user-supplied folder.
#'
#' Accepts either the \code{combined} folder or the \code{combined/txt}
#' subfolder; searches recursively and prefers matches inside a \code{txt}
#' directory when more than one exists.
#'
#' @param path Character. Folder to search.
#' @return Named list with elements \code{msms}, \code{evidence},
#'   \code{summary}, \code{proteinGroups} (full path or \code{NA}).
find_mq_files <- function(path) {
  targets <- c(
    msms = "msms.txt",
    evidence = "evidence.txt",
    summary = "summary.txt",
    proteinGroups = "proteinGroups.txt"
  )
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
      pattern = paste0("^", f, "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(hits) == 0) {
      return(NA_character_)
    }
    in_txt <- grepl("/txt/", hits, fixed = TRUE)
    if (any(in_txt)) hits[in_txt][1] else hits[1]
  })
}


# ── summary.txt helpers ───────────────────────────────────────────────────────

#' Read the MaxQuant summary.txt file, dropping columns that are entirely empty.
#' @param path Character. Full path to summary.txt.
#' @return A data.table with one row per raw file (plus the MaxQuant Total row).
read_summary_file <- function(path) {
  dt <- data.table::fread(path)
  empty <- vapply(
    dt,
    function(x) all(is.na(x) | (is.character(x) & !nzchar(x))),
    logical(1)
  )
  dt[, .SD, .SDcols = !empty]
}


# ── proteinGroups.txt helpers ─────────────────────────────────────────────────

#' Quantity types that can be exported from proteinGroups.txt.
#' Names are display labels, values are the column prefixes.
MQ_QUANT_TYPES <- c(
  "Intensity" = "Intensity",
  "LFQ intensity" = "LFQ intensity",
  "iBAQ" = "iBAQ"
)

#' Read and filter proteinGroups.txt.
#'
#' Removes reverse hits, potential contaminants and groups only identified by
#' site, and adds a \code{Protein ID} column holding the first accession of
#' each protein group.
#'
#' @param path Character. Full path to proteinGroups.txt.
#' @return Filtered data.table with all original columns plus \code{Protein ID}.
read_protein_groups <- function(path) {
  pg <- data.table::fread(path)

  flagged <- function(col) {
    if (!col %in% names(pg)) {
      return(rep(FALSE, nrow(pg)))
    }
    x <- pg[[col]]
    !is.na(x) & as.character(x) == "+"
  }
  keep <- !(flagged("Reverse") |
    flagged("Potential contaminant") |
    flagged("Only identified by site"))
  pg <- pg[keep]
  pg[, `Protein ID` := sub(";.*$", "", `Protein IDs`)]
  pg
}


#' Per-sample column names for a quantity type.
#'
#' Sample-level columns are "<quant> <sample>"; the bare "<quant>" column is
#' the summed total and "iBAQ peptides" is a count, so both are excluded.
quant_sample_cols <- function(pg, quant) {
  cols <- grep(paste0("^", quant, " .+"), names(pg), value = TRUE)
  setdiff(cols, paste0(quant, " peptides"))
}


#' Quantity types present (with per-sample columns) in a proteinGroups table.
#' @return Named character vector, subset of \code{MQ_QUANT_TYPES}.
available_quant_types <- function(pg) {
  present <- vapply(
    MQ_QUANT_TYPES,
    function(q) length(quant_sample_cols(pg, q)) > 0,
    logical(1)
  )
  MQ_QUANT_TYPES[present]
}


#' Build a protein abundance matrix from a filtered proteinGroups table.
#'
#' @param pg         data.table from \code{read_protein_groups()}.
#' @param quant      Character. One of \code{MQ_QUANT_TYPES}.
#' @param zero_to_na Logical. Replace 0 (MaxQuant's "not quantified") with NA.
#' @param log2_transform Logical. Apply log2 to the abundance columns. Zeros
#'   are always converted to NA first, since log2(0) is -Inf.
#' @return A data.table with a \code{Protein ID} column followed by one
#'   abundance column per sample.
build_protein_matrix <- function(
  pg,
  quant = "Intensity",
  zero_to_na = FALSE,
  log2_transform = FALSE
) {
  cols <- quant_sample_cols(pg, quant)
  if (length(cols) == 0) {
    stop("No per-sample '", quant, "' columns found in proteinGroups.txt.")
  }
  out <- cbind(
    data.table::data.table(`Protein ID` = pg[["Protein ID"]]),
    pg[, .SD, .SDcols = cols]
  )
  data.table::setnames(out, cols, sub(paste0("^", quant, " "), "", cols))
  val_cols <- setdiff(names(out), "Protein ID")

  if (zero_to_na || log2_transform) {
    out[,
      (val_cols) := lapply(.SD, function(x) {
        x <- as.numeric(x)
        x[!is.na(x) & x == 0] <- NA_real_
        x
      }),
      .SDcols = val_cols
    ]
  }
  if (log2_transform) {
    out[, (val_cols) := lapply(.SD, log2), .SDcols = val_cols]
  }
  out
}


# ── msms.txt helpers ──────────────────────────────────────────────────────────

#' Read and pre-process the MaxQuant msms.txt file.
#' @param path Character. Full path to msms.txt.
#' @return A tibble ready for \code{tidy_msms()}.
read_msms_file <- function(path) {
  data.table::fread(
    path,
    select = c(
      "Raw file",
      "Charge",
      "m/z",
      "Retention time",
      "Sequence",
      "Gene Names",
      "Matches",
      "Intensities",
      "Masses",
      "Intensities2",
      "Masses2",
      "Number of matches"
    )
  )[, Charge := paste0(as.character(Charge), "+")]
}


#' Tidy a raw msms tibble for a single peptide sequence.
#' @param data        Tibble from \code{read_msms_file()}.
#' @param peptide_seq Character. Peptide sequence to subset.
#' @return Tidy tibble with one row per fragment ion.
tidy_msms <- function(data, peptide_seq) {
  # Ensure 'data' is a data.table
  data.table::setDT(data)
  data[
    Sequence == peptide_seq,
    .(
      `Raw file`,
      `Gene Names`,
      Sequence,
      `m/z`,
      `Retention time`,
      Charge,
      Match = unlist(strsplit(as.character(Matches), ";")),
      Intensity = as.numeric(unlist(strsplit(as.character(Intensities), ";"))),
      MZ = as.numeric(unlist(strsplit(as.character(Masses), ";")))
    )
  ][!is.na(MZ) & Intensity > 0]
}


#' Build the annotated MS/MS spectrum ggplot object.
#' @param tidy_data  Tibble from \code{tidy_msms()}.
#' @param label_size Numeric. Font size for ion annotations.
#' @return A ggplot object.
build_msms_spectrum <- function(tidy_data, label_size = 3) {
  if (nrow(tidy_data) == 0) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "Peptide not found in the uploaded file.",
          size = 6,
          colour = "grey50"
        ) +
        theme_void()
    )
  }

  gene_name <- dplyr::first(stats::na.omit(tidy_data$`Gene Names`))
  peptide <- dplyr::first(tidy_data$Sequence)

  ggplot(tidy_data, aes(x = MZ, y = Intensity)) +
    geom_segment(aes(xend = MZ, yend = 0), colour = "grey25", linewidth = 0.4) +
    geom_text(
      aes(label = Match, color = stringr::str_detect(Match, "^y")),
      vjust = -0.8,
      size = label_size,
      fontface = "bold",
      check_overlap = TRUE
    ) +
    scale_color_manual(
      values = c("TRUE" = "#d95f02", "FALSE" = "#1b9e77"),
      guide = "none"
    ) +
    facet_wrap(
      ~ `Raw file` + Charge + `Retention time`,
      ncol = 2,
      scales = "free",
      labeller = label_both
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15)),
      labels = scales::label_scientific()
    ) +
    scale_x_continuous(
      breaks = tidy_data$MZ,
      labels = scales::label_number(accuracy = 0.1)
    ) +
    labs(
      title = paste0("MS/MS Fragmentation: ", peptide, " (", gene_name, ")"),
      x = "*m/z*",
      y = "Intensity"
    ) +
    theme_mq() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(color = "black", face = "bold", size = 8),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 8,
        color = "black"
      ),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.25)
    )
}


# ── evidence.txt helpers ──────────────────────────────────────────────────────

#' Read the MaxQuant evidence.txt file.
#'
#' @param path Character. Full path to evidence.txt.
#' @return A tibble.
read_evidence_file <- function(path) {
  dt <- fread(
    path,
    select = c(
      "Raw file",
      "Sequence",
      "Length",
      "Modifications",
      "Missed cleavages",
      "Charge",
      "m/z",
      "Mass",
      "Mass error [ppm]",
      "Mass error [Da]",
      "Retention time",
      "Number of data points",
      "Type",
      "PEP",
      "Taxonomy names"
    )
  )

  dt[, `:=`(
    Charge = as.character(Charge),
    `Missed cleavages` = as.character(`Missed cleavages`)
  )]
}


# ── Individual evidence plot builders ─────────────────────────────────────────

plot_ev_mz_rt <- function(df) {
  df |>
    ggplot(aes(x = `Retention time`, y = `m/z`)) +
    ggpointdensity::geom_pointdensity(method = "kde2d", adjust = 3) +
    scale_color_viridis_c(option = "D", direction = -1) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Retention time (min)", y = "m/z", color = "Density") +
    theme_mq()
}

plot_ev_ndp_dist <- function(df) {
  df |>
    ggplot(aes(x = `Number of data points`)) +
    geom_density(
      fill = "#1b9e77",
      alpha = 0.8,
      color = "white",
      linewidth = 0.25
    ) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Number of data points", y = "Density") +
    theme_mq()
}

plot_ev_ndp_mz <- function(df) {
  df |>
    ggplot(aes(x = `m/z`, y = `Number of data points`)) +
    ggpointdensity::geom_pointdensity(method = "kde2d") +
    scale_color_viridis_c(option = "D", direction = -1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "m/z", y = "Number of data points", color = "Density") +
    theme_mq() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
}

plot_ev_length <- function(df) {
  df |>
    ggplot(aes(x = Length)) +
    geom_density(
      fill = "#1b9e77",
      alpha = 0.8,
      color = "white",
      linewidth = 0.25
    ) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Peptide length (number of amino acids)", y = "Density") +
    theme_mq()
}

plot_ev_modifications <- function(df) {
  df |>
    ggplot(aes(y = Modifications)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      position = position_stack(vjust = 0.95),
      size = 5,
      fontface = "bold"
    ) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(y = NULL, x = "Count") +
    theme_mq() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
}

plot_ev_missed_cleavages <- function(df) {
  df |>
    ggplot(aes(x = `Missed cleavages`)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      position = position_stack(vjust = 1),
      size = 5,
      fontface = "bold"
    ) +
    scale_x_discrete(limits = as.character(0:5)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Number of missed cleavages", y = "Count") +
    theme_mq()
}

plot_ev_id_type <- function(df) {
  df |>
    ggplot(aes(x = Type)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      position = position_stack(vjust = 0.95),
      size = 5,
      fontface = "bold"
    ) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Identification type", y = "Count") +
    theme_mq() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
}

plot_ev_charge_bar <- function(df) {
  df |>
    ggplot(aes(x = Charge)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      position = position_stack(vjust = 0.99),
      size = 5,
      fontface = "bold"
    ) +
    scale_x_discrete(limits = as.character(1:6)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Charge state", y = "Count") +
    theme_mq()
}

plot_ev_mz_dist <- function(df) {
  df |>
    ggplot(aes(x = `m/z`)) +
    geom_density(
      fill = "#1b9e77",
      alpha = 0.8,
      color = "white",
      linewidth = 0.25
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "*m/z*", y = "Density") +
    theme_mq() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
}

plot_ev_mass_dist <- function(df) {
  df |>
    ggplot(aes(x = Mass)) +
    geom_density(
      fill = "#1b9e77",
      alpha = 0.8,
      color = "white",
      linewidth = 0.25
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Mass (Da)", y = "Density") +
    theme_mq() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
}

plot_ev_mass_err_ppm <- function(df) {
  df |>
    ggplot(aes(x = `Mass error [ppm]`)) +
    geom_density(
      fill = "#1b9e77",
      alpha = 0.8,
      color = "white",
      linewidth = 0.25
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Mass error (ppm)", y = "Density") +
    theme_mq()
}

plot_ev_mass_err_da <- function(df) {
  df |>
    ggplot(aes(x = `Mass error [Da]`)) +
    geom_density(
      fill = "#1b9e77",
      alpha = 0.8,
      color = "white",
      linewidth = 0.25
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Mass error (Da)", y = "Density") +
    theme_mq() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
}

plot_ev_pep <- function(df) {
  df |>
    ggplot(aes(x = PEP)) +
    geom_density(
      fill = "#1b9e77",
      alpha = 0.8,
      color = "white",
      linewidth = 0.25
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Posterior error probability (PEP)", y = "Density") +
    theme_mq() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
}

plot_ev_taxonomy <- function(df) {
  df |>
    ggplot(aes(y = `Taxonomy names`)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      position = position_stack(vjust = 0.95),
      size = 5,
      fontface = "bold"
    ) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(y = NULL, x = "Peptide count") +
    theme_mq() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
}


# ── Shared facet-height helper ────────────────────────────────────────────────

#' Compute a sensible plot height (px) based on the number of raw files.
#' Assumes 3-column facet layout, ~350 px per row, minimum 400 px.
#'
#' @param df  A tibble with a \code{`Raw file`} column.
#' @param row_px  Height in pixels per facet row (default 350).
#' @return Integer pixel height.
facet_height_px <- function(df, row_px = 350L) {
  n_files <- dplyr::n_distinct(df$`Raw file`)
  n_rows <- ceiling(n_files / 3L)
  max(400L, n_rows * row_px)
}


# ── 2. SIDEBAR UI ───────────────────────────────────────────────────────────

MaxQuantMSMS_sidebar_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tags$div(
      id = ns("sidebar_content"),

      # ── Data Input ────────────────────────────────────────────────────────
      tags$div(class = "sidebar-section-label", "MaxQuant Output Folder"),

      tags$div(
        style = "padding:0 8px;",
        textInput(
          ns("mq_folder"),
          "Path to MaxQuant 'combined' folder",
          value = "",
          placeholder = "/path/to/combined"
        ),
        tags$p(
          style = "color:#adb5bd;font-size:11px;margin-top:-6px;",
          "msms.txt, evidence.txt, summary.txt and proteinGroups.txt are ",
          "located automatically in the txt/ subfolder."
        ),
        tags$div(
          style = "text-align:center;",
          actionButton(
            ns("load_files"),
            "Load MaxQuant Files",
            class = "btn-primary",
            style = "width:80%;font-weight:bold;margin-bottom:6px;"
          )
        ),
        uiOutput(ns("file_status"))
      ),

      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── MS/MS Parameters ──────────────────────────────────────────────────
      tags$div(class = "sidebar-section-label", "MS/MS Parameters"),

      selectizeInput(
        ns("peptide_seq"),
        "Select Peptide Sequence",
        choices = NULL,
        options = list(
          placeholder = "Load MaxQuant files first…",
          maxOptions = 5000,
          searchField = "value"
        )
      ),

      numericInput(
        ns("label_size"),
        "Ion Label Size",
        value = 3,
        min = 1,
        max = 8,
        step = 0.5
      ),

      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── Evidence QC Parameters ────────────────────────────────────────────
      tags$div(class = "sidebar-section-label", "Evidence QC Parameters"),

      selectInput(
        ns("ev_plot_select"),
        "Select QC Plot",
        choices = c(
          "m/z vs Retention Time" = "mz_rt",
          "Data Points Distribution" = "ndp_dist",
          "Data Points vs m/z" = "ndp_mz",
          "Peptide Length" = "length",
          "Modifications" = "modifications",
          "Missed Cleavages" = "missed_cleavages",
          "Identification Type" = "id_type",
          "Charge State Distribution" = "charge_bar",
          "m/z Distribution" = "mz_dist",
          "Mass Distribution" = "mass_dist",
          "Mass Error (ppm)" = "mass_err_ppm",
          "Mass Error (Da)" = "mass_err_da",
          "PEP Distribution" = "pep",
          "Taxonomy Names" = "taxonomy"
        ),
        selected = "charge_rt"
      ),

      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── Actions ───────────────────────────────────────────────────────────
      tags$div(
        style = "padding:0 8px;text-align:center;",
        actionButton(
          ns("run_msms"),
          "Plot MS/MS Spectrum",
          class = "btn-primary",
          style = "width:80%;font-weight:bold;margin-top:8px;margin-bottom:6px;"
        ),
        actionButton(
          ns("run_evidence"),
          "Plot Evidence QC",
          class = "btn-primary",
          style = "width:80%;font-weight:bold;margin-top:2px;margin-bottom:10px;"
        )
      ),

      tags$div(
        style = "padding:0 8px;",
        downloadButton(
          ns("download_msms_plot"),
          "\u2B07 MS/MS Plot (.pdf)",
          class = "dl-btn",
          style = "width:100%;text-align:left;margin-bottom:6px;"
        ),
        downloadButton(
          ns("download_msms_data"),
          "\u2B07 MS/MS Tidy Data (.tsv)",
          class = "dl-btn",
          style = "width:100%;text-align:left;margin-bottom:6px;"
        ),
        downloadButton(
          ns("download_ev_plot"),
          "\u2B07 Evidence Plot (.pdf)",
          class = "dl-btn",
          style = "width:100%;text-align:left;margin-bottom:6px;"
        ),
        downloadButton(
          ns("download_ev_data"),
          "\u2B07 Evidence Data (.tsv)",
          class = "dl-btn",
          style = "width:100%;text-align:left;margin-bottom:6px;"
        ),
        tags$hr(style = "border-color:#2d3741;margin:8px 0 4px 0;"),
        selectInput(
          ns("pg_quant"),
          "Protein abundance values",
          choices = MQ_QUANT_TYPES,
          selected = "Intensity"
        ),
        checkboxInput(
          ns("pg_zero_na"),
          "Replace zeros with NA",
          value = FALSE
        ),
        checkboxInput(
          ns("pg_log2"),
          "log2-transform values (zeros become NA)",
          value = FALSE
        ),
        downloadButton(
          ns("download_protein_matrix"),
          "\u2B07 Protein Abundance Matrix (.tsv)",
          class = "dl-btn",
          style = "width:100%;text-align:left;"
        )
      )
    ),
    tags$hr(style = "border-color:#2d3741;margin:4px 0;"),
    tags$div(class = "sidebar-section-label", "Modification Diagnostic"),
    tags$div(
      style = "padding:0 16px;",
      sliderInput(
        ns("rt_tolerance"),
        "ΔRT artifact threshold (min)",
        min = 0,
        max = 5,
        value = 0.5,
        step = 0.1
      )
    )
  )
}


# ── 3. BODY UI ──────────────────────────────────────────────────────────────

MaxQuantMSMS_body_ui <- function(id) {
  ns <- NS(id)

  tabsetPanel(
    id = ns("tabs"),
    type = "tabs",

    # ── Tab 1: Run Summary (summary.txt) ─────────────────────────────────
    tabPanel(
      title = tagList(icon("clipboard-list"), "Run Summary"),
      fluidRow(
        box(
          title = "Processing Status",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          uiOutput(ns("status_log_ui"))
        )
      ),
      fluidRow(
        box(
          title = "MaxQuant run summary (summary.txt)",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          DT::dataTableOutput(ns("run_summary_table"))
        )
      ),
      fluidRow(
        box(
          title = "Protein abundance matrix preview (proteinGroups.txt)",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          DT::dataTableOutput(ns("protein_matrix_table"))
        )
      )
    ),

    # ── Tab 2: MS/MS Spectrum ─────────────────────────────────────────────
    tabPanel(
      title = tagList(icon("chart-bar"), "MS/MS Spectrum"),
      fluidRow(
        box(
          title = "Annotated MS/MS Fragmentation Spectrum",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_spectrum"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("spectrum_ui"))
          )
        )
      )
    ),

    # ── Tab 2: Evidence QC ────────────────────────────────────────────────
    tabPanel(
      title = tagList(icon("microscope"), "Evidence QC"),
      fluidRow(
        box(
          title = uiOutput(ns("ev_plot_title")),
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_evidence"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("evidence_plot_ui"))
          )
        )
      )
    ),

    # ── Tab 3: MS/MS Tidy Data ────────────────────────────────────────────
    tabPanel(
      title = tagList(icon("table"), "MS/MS Data"),
      fluidRow(
        box(
          title = "Fragment ion data for selected peptide",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          DT::dataTableOutput(ns("tidy_table"))
        )
      )
    ),

    # ── Tab 4: Evidence Data Preview ─────────────────────────────────────
    tabPanel(
      title = tagList(icon("table"), "Evidence Data"),
      fluidRow(
        box(
          title = "Evidence table preview",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          DT::dataTableOutput(ns("evidence_table"))
        )
      )
    ),

    # ── Tab 5: File Summary ───────────────────────────────────────────────
    tabPanel(
      title = tagList(icon("list"), "MS/MS File Summary"),
      fluidRow(
        box(
          title = "Peptide sequences in msms.txt",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          DT::dataTableOutput(ns("summary_table"))
        )
      )
    ),

    # ── Tab 6: Modification Diagnostic ───────────────────────────────────────
    tabPanel(
      title = tagList(icon("flask"), "Modification Diagnostic"),
      fluidRow(
        box(
          title = "RT Shift Profile — Modified vs Unmodified Peptides",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
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
    )
  )
}


# ── 4. SERVER ────────────────────────────────────────────────────────────────

MaxQuantMSMS_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Spinner helpers ───────────────────────────────────────────────────
    hide_spinner <- function(sid) shinyjs::hide(id = sid)
    show_spinner <- function(sid) shinyjs::show(id = sid)

    # ════════════════════════════════════════════════════════════════════
    # 0. PROCESSING LOG + FILE LOADING
    # ════════════════════════════════════════════════════════════════════

    # Timestamped step log shown in the sidebar and the Run Summary tab so
    # the user can always tell whether the app is working or idle.
    status_log <- reactiveVal(
      data.frame(
        time = character(0),
        level = character(0),
        msg = character(0),
        stringsAsFactors = FALSE
      )
    )

    log_step <- function(
      msg,
      level = c("info", "ok", "warn", "error"),
      notify = TRUE
    ) {
      level <- match.arg(level)
      status_log(rbind(
        status_log(),
        data.frame(
          time = format(Sys.time(), "%H:%M:%S"),
          level = level,
          msg = msg,
          stringsAsFactors = FALSE
        )
      ))
      if (notify) {
        showNotification(
          msg,
          type = switch(
            level,
            info = "message",
            ok = "message",
            warn = "warning",
            error = "error"
          ),
          duration = if (level == "error") NULL else 4
        )
      }
    }

    mq_files <- reactiveVal(NULL)
    raw_msms_rv <- reactiveVal(NULL)
    raw_evidence_rv <- reactiveVal(NULL)
    raw_summary_rv <- reactiveVal(NULL)
    protein_groups_rv <- reactiveVal(NULL)

    # Read one file with progress + log; returns NULL on failure
    load_one <- function(label, path, reader, step, n_steps) {
      if (is.na(path)) {
        log_step(
          paste0(label, " not found — related tabs will stay empty."),
          "warn"
        )
        return(NULL)
      }
      size_mb <- round(file.info(path)$size / 1024^2, 1)
      log_step(
        sprintf("Reading %s (%s MB)…", label, size_mb),
        "info",
        notify = FALSE
      )
      setProgress(
        value = (step - 1) / n_steps,
        message = sprintf("Step %d/%d: reading %s", step, n_steps, label),
        detail = sprintf("%s MB — large files may take a minute", size_mb)
      )
      t0 <- Sys.time()
      res <- tryCatch(reader(path), error = function(e) {
        log_step(
          sprintf("Failed to read %s: %s", label, conditionMessage(e)),
          "error"
        )
        NULL
      })
      if (!is.null(res)) {
        secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
        log_step(
          sprintf(
            "%s loaded: %s rows in %s s",
            label,
            format(nrow(res), big.mark = ","),
            secs
          ),
          "ok"
        )
      }
      res
    }

    observeEvent(input$load_files, {
      shinyjs::disable("load_files")
      on.exit(shinyjs::enable("load_files"), add = TRUE)

      # Reset previous state
      raw_msms_rv(NULL)
      raw_evidence_rv(NULL)
      raw_summary_rv(NULL)
      protein_groups_rv(NULL)
      status_log(status_log()[0, ])

      folder <- trimws(input$mq_folder)
      if (!nzchar(folder)) {
        log_step(
          "Please enter the path to the MaxQuant 'combined' folder.",
          "error"
        )
        return()
      }
      if (!dir.exists(path.expand(folder))) {
        log_step(paste0("Folder not found: ", folder), "error")
        return()
      }

      log_step(
        paste0("Scanning ", folder, " for MaxQuant files…"),
        "info",
        notify = FALSE
      )
      files <- find_mq_files(folder)
      mq_files(files)

      found <- names(files)[!is.na(unlist(files))]
      missing <- setdiff(names(files), found)
      log_step(
        sprintf(
          "Found %d/4 files%s",
          length(found),
          if (length(missing)) {
            paste0(" (missing: ", paste(missing, collapse = ", "), ")")
          } else {
            ""
          }
        ),
        if (length(missing)) "warn" else "ok"
      )
      if (length(found) == 0) {
        return()
      }

      withProgress(message = "Loading MaxQuant files", value = 0, {
        # Small files first so the Run Summary tab populates quickly
        raw_summary_rv(load_one(
          "summary.txt",
          files$summary,
          read_summary_file,
          1,
          4
        ))
        protein_groups_rv(load_one(
          "proteinGroups.txt",
          files$proteinGroups,
          read_protein_groups,
          2,
          4
        ))
        # Offer only the quantity types actually present in this run
        if (!is.null(protein_groups_rv())) {
          avail <- available_quant_types(protein_groups_rv())
          if (length(avail) == 0) {
            log_step(
              "proteinGroups.txt has no per-sample Intensity / LFQ / iBAQ columns.",
              "warn"
            )
          } else {
            log_step(
              paste0(
                "Abundance values available: ",
                paste(names(avail), collapse = ", ")
              ),
              "info",
              notify = FALSE
            )
          }
          updateSelectInput(
            session,
            "pg_quant",
            choices = avail,
            selected = if ("Intensity" %in% avail) "Intensity" else avail[1]
          )
        }
        raw_evidence_rv(load_one(
          "evidence.txt",
          files$evidence,
          read_evidence_file,
          3,
          4
        ))
        raw_msms_rv(load_one("msms.txt", files$msms, read_msms_file, 4, 4))
        setProgress(1, message = "Loading complete", detail = "")
      })

      log_step("All available files processed. Ready for analysis.", "ok")
    })

    # Sidebar status block
    output$file_status <- renderUI({
      files <- mq_files()
      if (is.null(files)) {
        return(tags$p(
          style = "color:#adb5bd;font-size:11px;",
          icon("circle-info"),
          " No folder loaded yet."
        ))
      }
      loaded <- list(
        summary = raw_summary_rv(),
        proteinGroups = protein_groups_rv(),
        evidence = raw_evidence_rv(),
        msms = raw_msms_rv()
      )
      row <- function(nm) {
        if (is.na(files[[nm]])) {
          tags$li(
            style = "color:#e74c3c;",
            icon("xmark"),
            " ",
            nm,
            ".txt missing"
          )
        } else if (is.null(loaded[[nm]])) {
          tags$li(
            style = "color:#f39c12;",
            icon("spinner", class = "fa-spin"),
            " ",
            nm,
            ".txt found — loading…"
          )
        } else {
          tags$li(
            style = "color:#2ecc71;",
            icon("check"),
            " ",
            nm,
            ".txt (",
            format(nrow(loaded[[nm]]), big.mark = ","),
            " rows)"
          )
        }
      }
      tags$ul(
        style = "list-style:none;padding-left:4px;font-size:11px;margin-top:4px;",
        lapply(c("summary", "proteinGroups", "evidence", "msms"), row)
      )
    })

    # Full step log in the Run Summary tab
    output$status_log_ui <- renderUI({
      lg <- status_log()
      if (nrow(lg) == 0) {
        return(tags$p(
          style = "color:#adb5bd;",
          icon("circle-info"),
          " Enter the MaxQuant 'combined' folder path in the sidebar and click ",
          tags$b("Load MaxQuant Files"),
          ". Progress will be reported here."
        ))
      }
      colors <- c(
        info = "#3498db",
        ok = "#2ecc71",
        warn = "#f39c12",
        error = "#e74c3c"
      )
      icons <- c(
        info = "circle-info",
        ok = "check",
        warn = "triangle-exclamation",
        error = "xmark"
      )
      tags$ul(
        style = "list-style:none;padding-left:0;font-family:monospace;font-size:12px;margin:0;",
        lapply(rev(seq_len(nrow(lg))), function(i) {
          tags$li(
            style = paste0("color:", colors[lg$level[i]], ";"),
            tags$span(style = "color:#7f8c8d;", lg$time[i]),
            "  ",
            icon(icons[lg$level[i]]),
            " ",
            lg$msg[i]
          )
        })
      )
    })

    # ════════════════════════════════════════════════════════════════════
    # A. msms.txt REACTIVES
    # ════════════════════════════════════════════════════════════════════

    raw_msms <- reactive({
      req(raw_msms_rv())
    })

    # Populate peptide selector after upload
    observeEvent(raw_msms(), {
      peptides <- sort(unique(raw_msms()$Sequence))
      updateSelectizeInput(
        session,
        "peptide_seq",
        choices = peptides,
        selected = peptides[1],
        server = TRUE
      )
    })

    # Tidy MS/MS data — recalculate only on button click
    tidy_msms_data <- eventReactive(input$run_msms, {
      if (is.null(raw_msms_rv())) {
        log_step(
          "msms.txt is not loaded — load the MaxQuant folder first.",
          "warn"
        )
      }
      req(raw_msms(), nchar(input$peptide_seq) > 0)
      show_spinner("sp_spectrum")
      log_step(
        paste0("Building MS/MS spectrum for ", input$peptide_seq, "…"),
        "info",
        notify = FALSE
      )
      withProgress(message = "Tidying MS/MS data…", value = 0.5, {
        result <- tidy_msms(raw_msms(), input$peptide_seq)
        incProgress(0.5, detail = "Done.")
        result
      })
      log_step(
        sprintf("MS/MS spectrum ready (%d fragment ions).", nrow(result)),
        "ok"
      )
      result
    })

    # Number of facets for dynamic height
    n_msms_facets <- reactive({
      req(tidy_msms_data())
      tidy_msms_data() |>
        dplyr::distinct(`Raw file`, Charge, `Retention time`) |>
        nrow()
    })

    msms_plot_height_px <- reactive({
      n_rows <- ceiling(n_msms_facets() / 3)
      max(400L, n_rows * 350L)
    })

    # ════════════════════════════════════════════════════════════════════
    # B. evidence.txt REACTIVES
    # ════════════════════════════════════════════════════════════════════

    raw_evidence <- reactive({
      req(raw_evidence_rv())
    })

    # The currently selected evidence plot — recalculate on button click
    current_ev_plot <- eventReactive(input$run_evidence, {
      if (is.null(raw_evidence_rv())) {
        log_step(
          "evidence.txt is not loaded — load the MaxQuant folder first.",
          "warn"
        )
      }
      req(raw_evidence())
      show_spinner("sp_evidence")

      df <- raw_evidence()
      sel <- input$ev_plot_select
      log_step(
        paste0("Building evidence QC plot: ", sel, "…"),
        "info",
        notify = FALSE
      )

      withProgress(message = "Building evidence plot…", value = 0.3, {
        p <- switch(
          sel,
          mz_rt = plot_ev_mz_rt(df),
          ndp_dist = plot_ev_ndp_dist(df),
          ndp_mz = plot_ev_ndp_mz(df),
          length = plot_ev_length(df),
          modifications = plot_ev_modifications(df),
          missed_cleavages = plot_ev_missed_cleavages(df),
          id_type = plot_ev_id_type(df),
          charge_bar = plot_ev_charge_bar(df),
          mz_dist = plot_ev_mz_dist(df),
          mass_dist = plot_ev_mass_dist(df),
          mass_err_ppm = plot_ev_mass_err_ppm(df),
          mass_err_da = plot_ev_mass_err_da(df),
          pep = plot_ev_pep(df),
          taxonomy = plot_ev_taxonomy(df)
        )
        incProgress(0.7, detail = "Done.")
        p
      })
      log_step("Evidence QC plot built — rendering…", "ok", notify = FALSE)
      p
    })

    # Dynamic height for evidence plots
    ev_plot_height_px <- reactive({
      req(raw_evidence())
      facet_height_px(raw_evidence())
    })

    # ════════════════════════════════════════════════════════════════════
    # C0. OUTPUTS — Run Summary / Protein matrix
    # ════════════════════════════════════════════════════════════════════

    output$run_summary_table <- DT::renderDataTable({
      req(raw_summary_rv())
      DT::datatable(
        raw_summary_rv(),
        rownames = FALSE,
        extensions = "FixedColumns",
        options = list(
          dom = "frtip",
          pageLength = 25,
          scrollX = TRUE,
          fixedColumns = list(leftColumns = 1)
        ),
        class = "display compact"
      )
    })

    # Abundance matrix for the currently selected quantity type
    protein_matrix <- reactive({
      req(protein_groups_rv(), input$pg_quant)
      req(input$pg_quant %in% available_quant_types(protein_groups_rv()))
      build_protein_matrix(
        protein_groups_rv(),
        input$pg_quant,
        zero_to_na = isTRUE(input$pg_zero_na),
        log2_transform = isTRUE(input$pg_log2)
      )
    })

    # Human-readable description of the current matrix settings
    protein_matrix_label <- reactive({
      paste0(
        input$pg_quant,
        if (isTRUE(input$pg_log2)) {
          " (log2, zeros -> NA)"
        } else if (isTRUE(input$pg_zero_na)) {
          " (zeros -> NA)"
        } else {
          ""
        }
      )
    })

    # Keep the two checkboxes consistent: log2 implies zeros -> NA
    observeEvent(
      input$pg_log2,
      {
        if (isTRUE(input$pg_log2)) {
          updateCheckboxInput(session, "pg_zero_na", value = TRUE)
          shinyjs::disable("pg_zero_na")
        } else {
          shinyjs::enable("pg_zero_na")
        }
      },
      ignoreInit = TRUE
    )

    output$protein_matrix_table <- DT::renderDataTable({
      pm <- protein_matrix()
      DT::datatable(
        head(pm, 500),
        rownames = FALSE,
        caption = paste0(
          "Values: ",
          protein_matrix_label(),
          " — ",
          format(nrow(pm), big.mark = ","),
          " protein groups × ",
          ncol(pm) - 1,
          " samples (first 500 rows shown)"
        ),
        options = list(dom = "frtip", pageLength = 20, scrollX = TRUE),
        class = "display compact"
      ) |>
        DT::formatSignif(
          columns = setdiff(names(pm), "Protein ID"),
          digits = 4
        )
    })

    # ════════════════════════════════════════════════════════════════════
    # C. OUTPUTS — MS/MS
    # ════════════════════════════════════════════════════════════════════

    # Dynamic container so the plot height can scale with facets
    output$spectrum_ui <- renderUI({
      plotOutput(
        ns("msms_spectrum"),
        height = paste0(msms_plot_height_px(), "px")
      )
    })

    output$msms_spectrum <- renderPlot({
      on.exit(hide_spinner("sp_spectrum"), add = TRUE)
      req(tidy_msms_data())
      build_msms_spectrum(tidy_msms_data(), label_size = input$label_size)
    })

    output$tidy_table <- DT::renderDataTable({
      req(tidy_msms_data())
      DT::datatable(
        tidy_msms_data(),
        rownames = FALSE,
        options = list(dom = "frtip", pageLength = 20, scrollX = TRUE),
        class = "display compact"
      )
    })

    output$summary_table <- DT::renderDataTable({
      req(raw_msms())
      summary_df <- raw_msms() |>
        dplyr::count(Sequence, `Gene Names`, name = "n_spectra") |>
        dplyr::arrange(dplyr::desc(n_spectra))
      DT::datatable(
        summary_df,
        rownames = FALSE,
        options = list(dom = "frtip", pageLength = 25, scrollX = TRUE),
        class = "display compact"
      )
    })

    # ════════════════════════════════════════════════════════════════════
    # D. OUTPUTS — Evidence QC
    # ════════════════════════════════════════════════════════════════════

    # Dynamically update box title to reflect the chosen plot
    output$ev_plot_title <- renderUI({
      label_map <- c(
        mz_rt = "m/z vs Retention Time",
        ndp_dist = "Number of Data Points Distribution",
        ndp_mz = "Number of Data Points vs m/z",
        length = "Peptide Length Distribution",
        modifications = "Modification Distribution",
        missed_cleavages = "Missed Cleavage Distribution",
        id_type = "Identification Type",
        charge_bar = "Charge State Distribution",
        mz_dist = "m/z Distribution",
        mass_dist = "Mass Distribution",
        mass_err_ppm = "Mass Error Distribution (ppm)",
        mass_err_da = "Mass Error Distribution (Da)",
        pep = "PEP Distribution",
        taxonomy = "Taxonomy Names"
      )
      label_map[input$ev_plot_select]
    })

    # Dynamic container so the plot height can scale with raw files
    output$evidence_plot_ui <- renderUI({
      plotOutput(
        ns("evidence_plot"),
        height = paste0(ev_plot_height_px(), "px")
      )
    })

    output$evidence_plot <- renderPlot({
      on.exit(hide_spinner("sp_evidence"), add = TRUE)
      req(current_ev_plot())
      current_ev_plot()
    })

    output$evidence_table <- DT::renderDataTable({
      req(raw_evidence())
      DT::datatable(
        head(raw_evidence(), 500),
        rownames = FALSE,
        options = list(dom = "frtip", pageLength = 20, scrollX = TRUE),
        class = "display compact"
      )
    })

    # ════════════════════════════════════════════════════════════════════
    # E. DOWNLOAD HANDLERS
    # ════════════════════════════════════════════════════════════════════

    # MS/MS spectrum PDF
    output$download_msms_plot <- downloadHandler(
      filename = function() {
        paste0("msms_spectrum_", input$peptide_seq, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        p <- build_msms_spectrum(
          tidy_msms_data(),
          label_size = input$label_size
        )
        n_rows <- ceiling(n_msms_facets() / 3)
        ggplot2::ggsave(
          file,
          plot = p,
          device = "pdf",
          width = 16,
          height = max(5, n_rows * 4),
          units = "in"
        )
      }
    )

    # MS/MS tidy TSV
    output$download_msms_data <- downloadHandler(
      filename = function() {
        paste0("msms_tidy_", input$peptide_seq, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        data.table::fwrite(tidy_msms_data(), file = file, sep = ",", na = "NA")
      }
    )

    # Evidence plot PDF
    output$download_ev_plot <- downloadHandler(
      filename = function() {
        paste0("evidence_", input$ev_plot_select, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        p <- current_ev_plot()
        n_rows <- ceiling(
          dplyr::n_distinct(raw_evidence()$`Raw file`) / 3
        )
        ggplot2::ggsave(
          file,
          plot = p,
          device = "pdf",
          width = 16,
          height = max(5, n_rows * 4),
          units = "in"
        )
      }
    )

    # Evidence tidy TSV
    output$download_ev_data <- downloadHandler(
      filename = function() {
        paste0("evidence_data_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        data.table::fwrite(raw_evidence(), file = file, sep = ",", na = "NA")
      }
    )

    # Protein abundance matrix TSV
    output$download_protein_matrix <- downloadHandler(
      filename = function() {
        tag <- gsub("[^A-Za-z0-9]+", "_", input$pg_quant)
        if (isTRUE(input$pg_log2)) {
          tag <- paste0(tag, "_log2")
        } else if (isTRUE(input$pg_zero_na)) {
          tag <- paste0(tag, "_zeroNA")
        }
        paste0("protein_abundance_matrix_", tag, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        pm <- protein_matrix()
        data.table::fwrite(pm, file = file, sep = "\t", na = "NA")
        log_step(
          sprintf(
            "Protein abundance matrix (%s) downloaded: %s proteins x %d samples.",
            protein_matrix_label(),
            format(nrow(pm), big.mark = ","),
            ncol(pm) - 1
          ),
          "ok"
        )
      }
    )

    # ════════════════════════════════════════════════════════════════════════
    # MODIFICATION DIAGNOSTIC
    # ════════════════════════════════════════════════════════════════════════
    mod_diag_data <- reactive({
      d <- raw_evidence()
      req(d, nrow(d) > 0)
      req(all(
        c("Sequence", "Modifications", "Retention time", "Raw file") %in%
          names(d)
      ))

      d_mod <- d |>
        dplyr::filter(
          !is.na(Modifications),
          !Modifications %in% c("", "Unmodified")
        ) |>
        dplyr::select(
          `Sequence`,
          `Modifications`,
          `Retention time`,
          `Raw file`
        ) |>
        dplyr::rename(
          peptide = `Sequence`,
          mod_label = `Modifications`,
          rt_mod = `Retention time`,
          sample_name = `Raw file`
        )

      d_unmod <- d |>
        dplyr::filter(
          is.na(Modifications) | Modifications %in% c("", "Unmodified")
        ) |>
        dplyr::group_by(`Sequence`, `Raw file`) |>
        dplyr::summarise(
          rt_unmod = median(`Retention time`, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(peptide = `Sequence`, sample_name = `Raw file`)

      dplyr::inner_join(d_mod, d_unmod, by = c("peptide", "sample_name")) |>
        dplyr::mutate(
          delta_rt = abs(rt_mod - rt_unmod),
          classification = dplyr::if_else(
            delta_rt <= input$rt_tolerance,
            "Suspected Artifact",
            "Sample-Derived"
          )
        )
    })

    mod_diag_plot_obj <- reactive({
      paired <- mod_diag_data()
      req(paired, nrow(paired) > 0)

      top_peptides <- paired |>
        dplyr::group_by(sample_name, peptide) |>
        dplyr::summarise(
          max_delta = max(delta_rt, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::group_by(sample_name) |>
        dplyr::slice_max(order_by = max_delta, n = 30) |>
        dplyr::pull(peptide) |>
        unique()

      plot_data <- paired |> dplyr::filter(peptide %in% top_peptides)
      req(nrow(plot_data) > 0)

      ggplot(plot_data, aes(y = reorder(peptide, delta_rt))) +
        geom_segment(
          aes(
            x = rt_unmod,
            xend = rt_mod,
            yend = reorder(peptide, delta_rt),
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
          title = "RT shift profile \u2014 modified vs unmodified",
          x = "Retention time (min)",
          y = "Peptide sequence",
          caption = paste0(
            "\u0394RT threshold: ",
            input$rt_tolerance,
            " min | Artifact if |\u0394RT| \u2264 threshold"
          )
        ) +
        facet_wrap(~sample_name, ncol = 2, scales = "free") +
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
      paired <- tryCatch(mod_diag_data(), error = function(e) NULL)
      shinyjs::hide("sp_moddiag")
      if (is.null(paired) || nrow(paired) == 0) {
        return(tags$p(
          style = "color:#adb5bd;text-align:center;padding:20px;",
          "No paired modified/unmodified peptides found. Load the MaxQuant folder (evidence.txt) and ensure it contains 'Sequence', 'Modifications', and 'Retention time' columns."
        ))
      }
      n_samples <- dplyr::n_distinct(paired$sample_name)
      dynamic_h <- max(400L, min(n_samples * 500L, 2000L))
      plotOutput(ns("mod_diag_plot"), height = paste0(dynamic_h, "px"))
    })

    output$mod_diag_plot <- renderPlot({
      mod_diag_plot_obj()
    })

    output$mod_diag_table <- DT::renderDataTable({
      paired <- mod_diag_data()
      req(paired, nrow(paired) > 0)
      tbl <- paired |>
        dplyr::select(
          sample_name,
          peptide,
          mod_label,
          rt_unmod,
          rt_mod,
          delta_rt,
          classification
        ) |>
        dplyr::arrange(dplyr::desc(delta_rt))
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
