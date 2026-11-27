## ============================================================
## mod_EncyclopeDIA.r  —  EncyclopeDIA DIA results viewer
## ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyjs)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(data.table)
  library(readr)
  library(purrr)
  library(viridis)
  library(DT)
})

# ── Constants & Helper functions ──────────────────────────────────────────────

ENC_THEME_CHOICES <- c(
  "Black & White (theme_bw)" = "bw",
  "Classic (theme_classic)" = "classic",
  "Minimal (theme_minimal)" = "minimal",
  "Light (theme_light)" = "light",
  "Linedraw (theme_linedraw)" = "linedraw",
  "Grey (theme_grey)" = "grey"
)

ENC_FEATURE_COLS <- c(
  "id",
  "Label",
  "ScanNr",
  "primary",
  "evalue",
  "correlationToGaussian",
  "correlationToPrecursor",
  "HyperScore",
  "xCorrLib",
  "xCorrModel",
  "scribe",
  "numberOfMatchingPeaks",
  "isotopeDotProduct",
  "averageFragmentDeltaMasses",
  "averageParentDeltaMass",
  "deltaRT",
  "ms1MassError",
  "ms2MassError",
  "numMissedCleavage",
  "pepLength",
  "charge1",
  "charge2",
  "charge3",
  "charge4",
  "precursorMz",
  "precursorMass",
  "RTinMin",
  "sequence",
  "Proteins"
)

#' Normalize a user-entered folder or file path (quotes, whitespace, ~).
normalize_enc_path <- function(path) {
  if (is.null(path) || length(path) == 0 || is.na(path[1])) {
    return("")
  }
  p <- trimws(as.character(path[1]))
  p <- gsub('^["\']+|["\']+$', "", p)
  p <- trimws(p)
  if (!nzchar(p)) "" else path.expand(p)
}

#' Derive a clean sample stem shared by *.encyclopedia(2).txt and *.features.txt.
enc_sample_stem <- function(path) {
  b <- basename(path)
  b <- sub("\\.(encyclopedia2?|features)\\.txt$", "", b, ignore.case = TRUE)
  b <- sub("\\.(mzML|mzXML|raw|dia|wiff)(\\.gz)?$", "", b, ignore.case = TRUE)
  b
}

#' Locate EncyclopeDIA PSM (*.encyclopedia.txt / *.encyclopedia2.txt) and
#' feature (*.features.txt, including *.mzML.features.txt) tables.
find_encyclopedia_files <- function(path) {
  p <- normalize_enc_path(path)
  empty <- list(psm = character(0), features = character(0))
  if (!nzchar(p)) {
    return(empty)
  }

  if (file.exists(p) && !dir.exists(p)) {
    if (
      grepl("\\.encyclopedia2?\\.txt$", p, ignore.case = TRUE) &&
        !grepl("\\.decoy\\.txt$", p, ignore.case = TRUE)
    ) {
      sib_feat <- sub(
        "\\.encyclopedia2?\\.txt$",
        ".features.txt",
        p,
        ignore.case = TRUE
      )
      return(list(
        psm = p,
        features = if (file.exists(sib_feat)) sib_feat else character(0)
      ))
    }
    if (grepl("\\.features\\.txt$", p, ignore.case = TRUE)) {
      return(list(psm = character(0), features = p))
    }
    return(empty)
  }

  if (!dir.exists(p)) {
    return(empty)
  }

  psm_hits <- list.files(
    p,
    pattern = "\\.encyclopedia2?\\.txt$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  psm_hits <- psm_hits[!grepl("\\.decoy\\.txt$", psm_hits, ignore.case = TRUE)]

  # If both .encyclopedia.txt and .encyclopedia2.txt exist for the same stem in
  # the same folder, prefer .encyclopedia2.txt to avoid duplicate runs.
  if (length(psm_hits) > 1L) {
    is_v2 <- grepl("\\.encyclopedia2\\.txt$", psm_hits, ignore.case = TRUE)
    key <- paste0(dirname(psm_hits), "/", enc_sample_stem(psm_hits))
    ord <- order(!is_v2)
    psm_hits <- psm_hits[ord][!duplicated(key[ord])]
    psm_hits <- sort(psm_hits)
  }

  feat_hits <- list.files(
    p,
    pattern = "\\.features\\.txt$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  list(psm = sort(psm_hits), features = sort(feat_hits))
}

#' Vectorized Kyte-Doolittle GRAVY calculation (matches utils_fasta.r GRAVY()).
calculate_GRAVY_vec <- function(sequence) {
  seqs <- toupper(gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", as.character(sequence)))
  n <- length(seqs)
  if (n == 0L) {
    return(numeric(0))
  }
  hi <- c(
    A = 1.8,
    R = -4.5,
    N = -3.5,
    D = -3.5,
    C = 2.5,
    Q = -3.5,
    E = -3.5,
    G = -0.4,
    H = -3.2,
    I = 4.5,
    L = 3.8,
    K = -3.9,
    M = 1.9,
    F = 2.8,
    P = -1.6,
    S = -0.8,
    T = -0.7,
    W = -0.9,
    Y = -1.3,
    V = 4.2
  )
  len <- nchar(seqs)
  total <- numeric(n)
  for (aa in names(hi)) {
    k <- nchar(gsub(paste0("[^", aa, "]"), "", seqs))
    total <- total + k * hi[[aa]]
  }
  out <- ifelse(len > 0L & !is.na(sequence), total / len, NA_real_)
  out
}

#' Compute GRAVY, pI, and MW once per unique stripped peptide and map back.
compute_enc_peptide_props <- function(df) {
  if (is.null(df) || nrow(df) == 0L || !"stripped_peptide" %in% names(df)) {
    return(df)
  }
  uniq_raw <- unique(df$stripped_peptide)
  safe_uniq <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", toupper(as.character(uniq_raw)))
  gravy_u <- calculate_GRAVY_vec(safe_uniq)
  pi_u <- calculate_pI(safe_uniq)
  mw_u <- calculate_MW(safe_uniq)

  idx <- match(df$stripped_peptide, uniq_raw)
  df$gravy <- gravy_u[idx]
  df$pI <- pi_u[idx]
  df$MW <- mw_u[idx]
  df
}

#' Quietly read a TSV via fread while muffling Percolator's `pi_0=...` footer warning.
fread_enc_tsv <- function(path, ...) {
  withCallingHandlers(
    data.table::fread(path, sep = "\t", header = TRUE, ...),
    warning = function(w) {
      if (grepl("Discarded single-line footer: <<pi_0=", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Derive precursor charge from one-hot charge1..charge4 columns or PSMId.
derive_enc_charge <- function(df) {
  n <- nrow(df)
  chg <- rep(NA_integer_, n)
  for (z in 1:4) {
    col <- paste0("charge", z)
    if (col %in% names(df)) {
      hit <- !is.na(df[[col]]) & as.numeric(df[[col]]) > 0
      chg[is.na(chg) & hit] <- z
    }
  }
  if ("Charge" %in% names(df)) {
    miss <- is.na(chg)
    chg[miss] <- suppressWarnings(as.integer(df$Charge[miss]))
  } else if ("charge" %in% names(df) && !all(is.na(df$charge))) {
    miss <- is.na(chg)
    chg[miss] <- suppressWarnings(as.integer(df$charge[miss]))
  }
  if ("PSMId" %in% names(df)) {
    miss <- is.na(chg)
    if (any(miss)) {
      chg[miss] <- suppressWarnings(as.integer(
        stringr::str_extract(df$PSMId[miss], "(?<=\\+)\\d+$")
      ))
    }
  }
  chg
}

#' Parse one EncyclopeDIA target PSM file (*.encyclopedia.txt or *.encyclopedia2.txt).
read_encyclopedia_psm_file <- function(path, sample_name = NULL) {
  d <- as.data.frame(fread_enc_tsv(path))

  # Standardize column aliases across EncyclopeDIA v1 / v2 / Percolator exports
  if (!"PSMId" %in% names(d) && "id" %in% names(d)) {
    d$PSMId <- d$id
  }
  if (!"peptide" %in% names(d)) {
    if ("Peptide" %in% names(d)) {
      d$peptide <- d$Peptide
    } else if ("sequence" %in% names(d)) {
      d$peptide <- d$sequence
    }
  }
  if (!"proteinIds" %in% names(d)) {
    if ("Protein" %in% names(d)) {
      d$proteinIds <- d$Protein
    } else if ("Proteins" %in% names(d)) {
      d$proteinIds <- d$Proteins
    }
  }
  if (!"q-value" %in% names(d)) {
    if ("Qvalue" %in% names(d)) {
      d$`q-value` <- d$Qvalue
    } else if ("q_value" %in% names(d)) {
      d$`q-value` <- d$q_value
    }
  }
  if (!"score" %in% names(d) && "Score" %in% names(d)) {
    d$score <- d$Score
  }
  if (!"posterior_error_prob" %in% names(d)) {
    pep_alt <- intersect(c("PEP", "pep", "PosteriorErrorProbability"), names(d))
    if (length(pep_alt) > 0L) {
      d$posterior_error_prob <- d[[pep_alt[1]]]
    }
  }

  req_cols <- c("peptide", "proteinIds", "q-value")
  missing_cols <- setdiff(req_cols, names(d))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "%s is missing required column(s): %s",
        basename(path),
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Drop trailing Percolator `pi_0=null` row if fread kept it (e.g. 1-row files)
  if ("PSMId" %in% names(d)) {
    d <- d[!grepl("^pi_0=", as.character(d$PSMId)), , drop = FALSE]
  }
  d$`q-value` <- suppressWarnings(as.numeric(d$`q-value`))
  d <- d[!is.na(d$`q-value`) & !is.na(d$peptide) & nzchar(d$peptide), , drop = FALSE]

  fname <- if (!is.null(sample_name) && nzchar(sample_name)) {
    sample_name
  } else {
    enc_sample_stem(path)
  }

  d$filename <- fname
  d$stripped_peptide <- stringr::str_replace_all(
    as.character(d$peptide),
    "\\[[^\\]]*\\]|\\.|\\-",
    ""
  )
  d$pep_len <- nchar(d$stripped_peptide)
  d$n_mods <- stringr::str_count(as.character(d$peptide), "\\[[^\\]]+\\]")
  d$charge <- derive_enc_charge(d)
  d
}

#' Read an EncyclopeDIA *.features.txt table either to enrich known PSMIds or
#' as a standalone table with target-decoy competition q-values.
read_encyclopedia_feature_file <- function(
  path,
  target_ids = NULL,
  sample_name = NULL,
  max_qval = 0.25
) {
  hdr <- names(fread_enc_tsv(path, nrows = 0))
  wanted <- if (!is.null(target_ids)) {
    setdiff(ENC_FEATURE_COLS, c("sequence", "Proteins", "pepLength", "Label"))
  } else {
    ENC_FEATURE_COLS
  }
  avail <- intersect(wanted, hdr)
  if (!"id" %in% avail) {
    stop(
      sprintf("%s is missing the 'id' column expected in *.features.txt.", basename(path)),
      call. = FALSE
    )
  }

  dt <- fread_enc_tsv(path, select = avail)

  # Mode 1: Join feature metrics onto already-loaded PSMIds
  if (!is.null(target_ids)) {
    sub_dt <- dt[id %in% target_ids]
    rm(dt)
    if (nrow(sub_dt) == 0L) {
      return(NULL)
    }
    sub_dt <- unique(sub_dt, by = "id")
    sub_df <- as.data.frame(sub_dt)
    sub_df$charge_feat <- derive_enc_charge(sub_df)
    drop_cols <- intersect(
      c("sequence", "Proteins", "pepLength", "charge1", "charge2", "charge3", "charge4"),
      names(sub_df)
    )
    sub_df[, setdiff(names(sub_df), drop_cols), drop = FALSE]
  } else {
    # Mode 2: Standalone *.features.txt ingestion
    req_f <- c("sequence", "Proteins")
    miss_f <- setdiff(req_f, names(dt))
    if (length(miss_f) > 0L) {
      stop(
        sprintf(
          "%s is missing required feature column(s): %s",
          basename(path),
          paste(miss_f, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    score_col <- intersect(c("primary", "scribe", "HyperScore", "xCorrLib"), names(dt))
    if (length(score_col) > 0L) {
      sc <- score_col[1]
      data.table::setorderv(dt, cols = sc, order = -1L, na.last = TRUE)
    } else {
      sc <- NA_character_
    }

    if ("Label" %in% names(dt)) {
      dt <- unique(dt, by = c("sequence", "Label"))
      has_decoys <- any(dt$Label == -1L, na.rm = TRUE) && any(dt$Label == 1L, na.rm = TRUE)
      if (has_decoys && !is.na(sc)) {
        data.table::setorderv(dt, cols = sc, order = -1L, na.last = TRUE)
        cum_t <- cumsum(dt$Label == 1L)
        cum_d <- cumsum(dt$Label == -1L)
        fdr <- (cum_d + 1) / pmax(cum_t, 1L)
        dt[, `q-value` := rev(cummin(rev(fdr)))]
      } else {
        dt[, `q-value` := 0]
      }
      dt <- dt[is.na(Label) | Label >= 0L]
    } else {
      dt <- unique(dt, by = "sequence")
      dt[, `q-value` := 0]
    }

    if (nrow(dt) > 100000L && !is.null(max_qval) && is.finite(max_qval)) {
      dt <- dt[`q-value` <= max_qval]
    }

    d <- as.data.frame(dt)
    rm(dt)

    fname <- if (!is.null(sample_name) && nzchar(sample_name)) {
      sample_name
    } else {
      enc_sample_stem(path)
    }

    d$PSMId <- as.character(d$id)
    d$score <- if (!is.na(sc)) as.numeric(d[[sc]]) else NA_real_
    d$posterior_error_prob <- NA_real_
    d$peptide <- as.character(d$sequence)
    d$proteinIds <- as.character(d$Proteins)
    d$filename <- fname
    d$stripped_peptide <- stringr::str_replace_all(
      d$peptide,
      "\\[[^\\]]*\\]|\\.|\\-",
      ""
    )
    d$pep_len <- nchar(d$stripped_peptide)
    d$n_mods <- stringr::str_count(d$peptide, "\\[[^\\]]+\\]")
    d$charge <- derive_enc_charge(d)

    drop_cols <- intersect(c("charge1", "charge2", "charge3", "charge4"), names(d))
    d[, setdiff(names(d), drop_cols), drop = FALSE]
  }
}

#' Load an EncyclopeDIA directory (PSMs, features, or joined PSM + features).
load_encyclopedia_directory <- function(
  dir_path,
  mode = c("auto", "psm_only", "features_only"),
  progress_cb = function(value, detail) NULL
) {
  mode <- match.arg(mode)
  files <- find_encyclopedia_files(dir_path)

  has_psm <- length(files$psm) > 0L
  has_feat <- length(files$features) > 0L

  if (mode == "psm_only" && !has_psm) {
    stop(
      "No *.encyclopedia.txt or *.encyclopedia2.txt files found in the directory.",
      call. = FALSE
    )
  }
  if (mode == "features_only" && !has_feat) {
    stop("No *.features.txt files found in the directory.", call. = FALSE)
  }
  if (!has_psm && !has_feat) {
    stop(
      "No EncyclopeDIA result files (*.encyclopedia2.txt, *.encyclopedia.txt, or *.features.txt) found.",
      call. = FALSE
    )
  }

  effective_mode <- if (mode == "auto") {
    if (has_psm && has_feat) {
      "psm_plus_features"
    } else if (has_psm) {
      "psm_only"
    } else {
      "features_only"
    }
  } else {
    mode
  }

  if (effective_mode %in% c("psm_only", "psm_plus_features")) {
    psm_stems <- enc_sample_stem(files$psm)
    if (anyDuplicated(psm_stems)) {
      psm_stems <- paste0(basename(dirname(files$psm)), "/", psm_stems)
    }
    feat_map <- stats::setNames(files$features, enc_sample_stem(files$features))
    n_f <- length(files$psm)
    joined_any <- FALSE

    rows_list <- vector("list", n_f)
    for (i in seq_len(n_f)) {
      fp <- files$psm[i]
      st <- psm_stems[i]
      progress_cb(
        (i - 0.5) / (n_f + 0.5),
        paste0("Reading ", basename(fp), " (", i, "/", n_f, ")")
      )
      d <- read_encyclopedia_psm_file(fp, sample_name = st)

      if (effective_mode == "psm_plus_features" && "PSMId" %in% names(d)) {
        raw_stem <- enc_sample_stem(fp)
        ff <- unname(feat_map[raw_stem])
        if (length(ff) == 1L && !is.na(ff) && file.exists(ff)) {
          progress_cb(
            (i - 0.15) / (n_f + 0.5),
            paste0("Joining features: ", basename(ff))
          )
          fsub <- read_encyclopedia_feature_file(ff, target_ids = d$PSMId)
          if (!is.null(fsub) && nrow(fsub) > 0L) {
            d <- dplyr::left_join(d, fsub, by = c("PSMId" = "id"))
            if ("charge_feat" %in% names(d)) {
              miss_c <- is.na(d$charge)
              d$charge[miss_c] <- d$charge_feat[miss_c]
              d$charge_feat <- NULL
            }
            joined_any <- TRUE
          }
        }
      }
      rows_list[[i]] <- d
    }
    res <- dplyr::bind_rows(rows_list)
    rm(rows_list)
    n_files_loaded <- n_f
    has_features_loaded <- joined_any
  } else {
    feat_stems <- enc_sample_stem(files$features)
    if (anyDuplicated(feat_stems)) {
      feat_stems <- paste0(basename(dirname(files$features)), "/", feat_stems)
    }
    n_f <- length(files$features)
    rows_list <- vector("list", n_f)
    for (i in seq_len(n_f)) {
      ff <- files$features[i]
      st <- feat_stems[i]
      progress_cb(
        (i - 0.3) / (n_f + 0.5),
        paste0("Reading features: ", basename(ff), " (", i, "/", n_f, ")")
      )
      rows_list[[i]] <- read_encyclopedia_feature_file(
        ff,
        target_ids = NULL,
        sample_name = st
      )
    }
    res <- dplyr::bind_rows(rows_list)
    rm(rows_list)
    n_files_loaded <- n_f
    has_features_loaded <- TRUE
  }

  if ("filename" %in% names(res) && nrow(res) > 0L) {
    nat_levels <- stringr::str_sort(unique(as.character(res$filename)), numeric = TRUE)
    res$filename <- factor(as.character(res$filename), levels = nat_levels)
  }

  progress_cb(0.95, "Computing peptide physicochemical properties...")
  res <- compute_enc_peptide_props(res)
  progress_cb(1, "Done")

  list(
    data = res,
    mode_used = effective_mode,
    n_files = n_files_loaded,
    n_psm_files = length(files$psm),
    n_feature_files = length(files$features),
    has_features = has_features_loaded
  )
}

#' Build consistent ggplot2 styling from user-selected theme and font size.
theme_enc_custom <- function(theme_name = "bw", font_size = 13) {
  base_th <- get_da_theme(theme_name = theme_name, base_size = font_size)
  base_th +
    theme(
      text = element_text(size = font_size),
      plot.title = element_text(
        size = round(font_size * 1.1),
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = round(font_size * 0.9),
        hjust = 0.5
      ),
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold")
    )
}

sp_wrap_enc <- function(sid, ui_el) {
  div(
    class = "plot-wrap",
    shinyjs::hidden(
      tags$div(
        class = "spinner-overlay",
        id = sid,
        icon("spinner", class = "fa-spin")
      )
    ),
    ui_el
  )
}

enc_empty_plot <- function(msg) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = msg,
      size = 5,
      fontface = "bold",
      color = "#6c757d"
    ) +
    theme_void()
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — SIDEBAR UI
# ═══════════════════════════════════════════════════════════════════════════════
EncyclopeDIA_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      style = "padding:12px 16px 4px;color:#ffffff;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      icon("book-open", lib = "font-awesome"),
      " EncyclopeDIA DIA"
    ),
    tags$div(
      style = "padding: 0 16px; margin-bottom: 12px;",
      textInput(
        ns("encyclopedia_dir"),
        "Path to results directory:",
        placeholder = "/path/to/encyclopedia_search/"
      ),
      tags$p(
        style = "color:#adb5bd;font-size:11px;margin-top:-6px;",
        "Automatically locates *.encyclopedia2.txt, *.encyclopedia.txt, and ",
        "*.mzML.features.txt (searching subfolders)."
      ),
      selectInput(
        ns("input_mode"),
        "Input tables",
        choices = c(
          "Auto (PSMs + Feature metrics)" = "auto",
          "PSM results only (*.encyclopedia*.txt)" = "psm_only",
          "Feature tables only (*.features.txt)" = "features_only"
        ),
        selected = "auto"
      ),
      tags$div(
        style = "text-align:center;",
        actionButton(
          ns("load_dir_btn"),
          "Load Directory",
          icon = icon("folder-open"),
          class = "btn-primary",
          style = "width:85%;font-weight:bold;margin-bottom:6px;"
        )
      ),
      uiOutput(ns("file_status"))
    ),
    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),
    tags$div(
      style = "padding: 0 8px;",
      sliderInput(
        ns("qval_filter"),
        "Max q-value",
        min = 0,
        max = 0.05,
        value = 0.01,
        step = 0.005
      )
    ),
    Fasta_sidebar_ui(ns("fasta"))
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — BODY UI
# ═══════════════════════════════════════════════════════════════════════════════
EncyclopeDIA_body_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      infoBoxOutput(ns("info_box_files"), width = 6),
      infoBoxOutput(ns("info_box_psms"), width = 6)
    ),
    tabsetPanel(
      id = ns("enc_tabs"),
      type = "tabs",

      # ── Interactive Plot Viewer ──────────────────────────────────────────────
      tabPanel(
        "Interactive Plot Viewer",
        fluidRow(
          box(
            title = "Plot Controls",
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
                  choices = list(
                    "Identifications & FDR" = c(
                      "Protein and Peptide Identifications" = "plot_identifications",
                      "Score Distribution" = "plot_score",
                      "Posterior Error Probability (PEP)" = "plot_pep",
                      "q-value Distribution" = "plot_qval",
                      "Peptide Yield vs. FDR" = "plot_fdr_curve"
                    ),
                    "Peptide Properties" = c(
                      "Charge State" = "plot_charge",
                      "Modifications" = "plot_mods",
                      "Peptide Length" = "plot_length",
                      "GRAVY Index Distribution" = "plot_gravy",
                      "Virtual 2D Gel (pI vs MW)" = "plot_pi",
                      "Amino acid frequencies" = "plot_aa_freq"
                    ),
                    "Feature Metrics (*.features.txt)" = c(
                      "Precursor m/z vs Retention Time" = "plot_mz_rt",
                      "Retention Time Error (ΔRT)" = "plot_delta_rt",
                      "MS1 & MS2 Mass Error (ppm)" = "plot_mass_error",
                      "Missed Cleavages" = "plot_missed_cleavages",
                      "Feature Scores (Scribe / HyperScore / xCorrLib)" = "plot_feature_scores",
                      "Isotope Dot Product" = "plot_idotp"
                    )
                  ),
                  width = "100%"
                )
              ),
              column(
                2,
                colourpicker::colourInput(
                  ns("bar_color"),
                  "Primary / Bar colour",
                  value = "#1b9e77"
                )
              ),
              column(
                2,
                colourpicker::colourInput(
                  ns("line_color"),
                  "Secondary / Line colour",
                  value = "#d95f02"
                )
              ),
              column(
                2,
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
                2,
                tags$div(
                  style = "margin-top:25px;",
                  downloadButton(
                    ns("download_plot"),
                    tagList(icon("download"), " Download (.png)"),
                    class = "dl-btn",
                    style = "width:100%;margin:0;padding:7px 12px;"
                  )
                )
              )
            ),
            tags$hr(style = "margin:8px 0 12px;"),
            fluidRow(
              column(
                3,
                tags$div(
                  style = "margin-top:25px;",
                  checkboxInput(
                    ns("facet_enable"),
                    "Facet by file",
                    value = TRUE
                  )
                )
              ),
              column(
                3,
                numericInput(
                  ns("facet_ncol"),
                  "Facet columns",
                  value = 3,
                  min = 1,
                  max = 12,
                  step = 1,
                  width = "100%"
                )
              ),
              column(
                3,
                numericInput(
                  ns("font_size"),
                  "Font size",
                  value = 13,
                  min = 8,
                  max = 28,
                  step = 1,
                  width = "100%"
                )
              ),
              column(
                3,
                selectInput(
                  ns("plot_theme"),
                  "Theme",
                  choices = ENC_THEME_CHOICES,
                  selected = "bw",
                  width = "100%"
                )
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Dynamic Plot View",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_enc(ns("spi_main"), uiOutput(ns("dynamic_plot_ui")))
          )
        )
      ),

      # ── Data Table ────────────────────────────────────────────────────────
      tabPanel(
        "Data Table",
        fluidRow(
          box(
            title = "Filtered EncyclopeDIA Records",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            div(
              style = "margin:6px 0;",
              DT::dataTableOutput(ns("res_table"))
            )
          )
        )
      )
    )
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — SERVER
# ═══════════════════════════════════════════════════════════════════════════════
EncyclopeDIA_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    fasta_digest <- Fasta_server("fasta")

    hide_sp <- function(s) shinyjs::hide(id = s)

    enc_data <- reactiveVal(NULL)
    enc_meta <- reactiveVal(NULL)
    status_msg <- reactiveVal(NULL)
    active_plot_id <- reactiveVal("plot_identifications")

    enc_aesthetics <- reactive({
      ncol_val <- suppressWarnings(as.integer(input$facet_ncol %||% 3L))
      if (is.na(ncol_val) || ncol_val < 1L) {
        ncol_val <- 3L
      }
      font_size <- suppressWarnings(as.numeric(input$font_size %||% 13))
      if (is.na(font_size) || font_size < 6) {
        font_size <- 13
      }
      theme_name <- input$plot_theme %||% "bw"
      facet_on <- if (is.null(input$facet_enable)) TRUE else isTRUE(input$facet_enable)
      bar_col <- input$bar_color %||% "#1b9e77"
      line_col <- input$line_color %||% "#d95f02"

      list(
        facet = facet_on,
        ncol = ncol_val,
        font_size = font_size,
        label_size = max(2.5, round(font_size * 0.28, 2)),
        theme_name = theme_name,
        theme_obj = theme_enc_custom(theme_name = theme_name, font_size = font_size),
        bar_color = bar_col,
        line_color = line_col
      )
    })

    add_enc_facet <- function(p, aes_cfg, scales = "fixed") {
      if (isTRUE(aes_cfg$facet)) {
        p + facet_wrap(~filename, ncol = aes_cfg$ncol, scales = scales)
      } else {
        p
      }
    }

    observeEvent(input$load_dir_btn, {
      dir_path <- normalize_enc_path(input$encyclopedia_dir)
      if (!nzchar(dir_path) || (!dir.exists(dir_path) && !file.exists(dir_path))) {
        msg <- "Please enter a valid EncyclopeDIA directory or file path."
        status_msg(list(type = "error", text = msg))
        showNotification(msg, type = "error")
        return()
      }

      mode_sel <- input$input_mode %||% "auto"
      shinyjs::show(id = "spi_main")
      on.exit(hide_sp("spi_main"), add = TRUE)

      withProgress(
        message = "Loading EncyclopeDIA files...",
        value = 0.05,
        {
          loaded <- tryCatch(
            load_encyclopedia_directory(
              dir_path = dir_path,
              mode = mode_sel,
              progress_cb = function(val, detail) {
                setProgress(value = val, detail = detail)
              }
            ),
            error = function(e) {
              status_msg(list(type = "error", text = conditionMessage(e)))
              showNotification(
                paste("Error loading EncyclopeDIA files:", conditionMessage(e)),
                type = "error"
              )
              NULL
            }
          )

          if (!is.null(loaded) && !is.null(loaded$data) && nrow(loaded$data) > 0L) {
            enc_data(loaded$data)
            enc_meta(loaded[setdiff(names(loaded), "data")])
            active_plot_id(input$plot_select %||% "plot_identifications")

            mode_lbl <- switch(
              loaded$mode_used,
              psm_plus_features = sprintf(
                "%d PSM file(s) + feature metrics joined",
                loaded$n_files
              ),
              psm_only = sprintf("%d PSM file(s) loaded", loaded$n_files),
              features_only = sprintf("%d feature table(s) loaded", loaded$n_files),
              sprintf("%d file(s) loaded", loaded$n_files)
            )
            summary_txt <- sprintf(
              "%s (%s records)",
              mode_lbl,
              format(nrow(loaded$data), big.mark = ",")
            )
            status_msg(list(type = "ok", text = summary_txt))
            showNotification(
              paste("Successfully loaded:", summary_txt),
              type = "message"
            )
          }
        }
      )
    })

    observeEvent(input$run_plot, {
      req(input$plot_select)
      active_plot_id(input$plot_select)
    })

    output$file_status <- renderUI({
      st <- status_msg()
      if (is.null(st)) {
        return(NULL)
      }
      col <- if (identical(st$type, "ok")) "#2ecc71" else "#e74c3c"
      ic <- if (identical(st$type, "ok")) "check-circle" else "exclamation-circle"
      tags$div(
        style = sprintf(
          "padding:6px 8px;margin-top:4px;border-radius:4px;background:#1f2c38;color:%s;font-size:11px;line-height:1.35;",
          col
        ),
        icon(ic),
        " ",
        st$text
      )
    })

    filtered_data <- reactive({
      req(enc_data())
      q_max <- input$qval_filter %||% 0.01
      enc_data() |> dplyr::filter(`q-value` <= q_max)
    })

    mapped_data <- reactive({
      d <- filtered_data()
      dig <- fasta_digest()
      if (!is.null(dig)) {
        classes <- classify_peptides(d$stripped_peptide, dig)
        d <- dplyr::left_join(d, classes, by = c("stripped_peptide" = "peptide"))
      } else {
        d$classification <- "Unmapped"
        d$mapped_proteins <- NA_character_
      }
      d
    })

    output$info_box_files <- renderInfoBox({
      d <- enc_data()
      if (is.null(d)) {
        return(infoBox(
          "Files Loaded",
          "0",
          subtitle = "Load a directory from the sidebar",
          icon = icon("folder-open"),
          color = "light-blue"
        ))
      }
      n_files <- dplyr::n_distinct(d$filename)
      meta <- enc_meta()
      sub_txt <- if (!is.null(meta) && isTRUE(meta$has_features)) {
        if (identical(meta$mode_used, "features_only")) {
          "Feature tables (*.features.txt)"
        } else {
          "PSM tables + *.features.txt metrics"
        }
      } else {
        "PSM tables (*.encyclopedia*.txt)"
      }
      infoBox(
        "Files Loaded",
        n_files,
        subtitle = sub_txt,
        icon = icon("file-alt"),
        color = "light-blue"
      )
    })

    output$info_box_psms <- renderInfoBox({
      d <- enc_data()
      if (is.null(d)) {
        return(infoBox(
          "Filtered Records",
          "—",
          subtitle = paste("q-value <=", input$qval_filter %||% 0.01),
          icon = icon("filter"),
          color = "green"
        ))
      }
      total <- nrow(d)
      filtered <- nrow(filtered_data())
      pct <- if (total > 0L) round(filtered / total * 100, 1) else 0
      infoBox(
        "Filtered Records",
        paste0(
          format(filtered, big.mark = ","),
          " / ",
          format(total, big.mark = ","),
          " (",
          pct,
          "%)"
        ),
        subtitle = paste("q-value <=", input$qval_filter %||% 0.01),
        icon = icon("filter"),
        color = "green"
      )
    })

    empty_filtered_plot <- function() {
      all_d <- enc_data()
      min_q <- if (!is.null(all_d) && nrow(all_d) > 0L && "q-value" %in% names(all_d)) {
        suppressWarnings(min(all_d$`q-value`, na.rm = TRUE))
      } else {
        NA_real_
      }
      hint <- if (is.finite(min_q)) {
        sprintf(" (minimum q-value in loaded data: %.4g)", min_q)
      } else {
        ""
      }
      enc_empty_plot(sprintf(
        "No records pass Max q-value <= %s%s.",
        input$qval_filter %||% 0.01,
        hint
      ))
    }

    plot_identifications_obj <- reactive({
      req(enc_data())
      d <- filtered_data()
      if (nrow(d) == 0L) {
        return(empty_filtered_plot())
      }
      aes_cfg <- enc_aesthetics()

      summ <- d |>
        dplyr::group_by(filename, .drop = FALSE) |>
        dplyr::summarise(
          n_proteins = dplyr::n_distinct(proteinIds[!is.na(proteinIds) & nzchar(proteinIds)]),
          n_peptides = dplyr::n_distinct(stripped_peptide[!is.na(stripped_peptide) & nzchar(stripped_peptide)]),
          .groups = "drop"
        )

      ggplot(summ, aes(x = filename)) +
        geom_col(aes(y = n_proteins), fill = aes_cfg$bar_color, width = 0.7) +
        geom_text(
          aes(
            y = n_proteins,
            label = ifelse(n_proteins > 0L, as.character(n_proteins), "")
          ),
          vjust = 1.5,
          size = aes_cfg$label_size,
          fontface = "bold",
          color = "white"
        ) +
        geom_line(
          aes(y = n_peptides, group = 1),
          color = aes_cfg$line_color,
          linewidth = 1.2
        ) +
        geom_point(
          aes(y = n_peptides),
          color = aes_cfg$line_color,
          size = 3
        ) +
        geom_text(
          aes(y = n_peptides, label = n_peptides),
          vjust = -0.8,
          size = aes_cfg$label_size,
          fontface = "bold",
          color = aes_cfg$line_color
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(
          title = "EncyclopeDIA Identifications",
          subtitle = paste("q-value <=", input$qval_filter %||% 0.01),
          x = "File",
          y = "Number of proteins (bar) or peptides (line)"
        ) +
        aes_cfg$theme_obj +
        theme(
          axis.text.x = element_text(angle = 55, hjust = 1, color = "black")
        )
    })

    plot_annotated_density_enc <- function(
      d,
      col_sym,
      title,
      xlab,
      color,
      aes_cfg
    ) {
      if (!col_sym %in% names(d) || all(is.na(d[[col_sym]]))) {
        return(enc_empty_plot(paste0("Column '", col_sym, "' not available in loaded data.")))
      }
      d_valid <- d |> dplyr::filter(is.finite(.data[[col_sym]]))
      if (nrow(d_valid) < 2L) {
        return(enc_empty_plot("Not enough finite observations to plot density."))
      }

      if (isTRUE(aes_cfg$facet)) {
        m_df <- d_valid |>
          dplyr::group_by(filename) |>
          dplyr::summarise(m = stats::median(.data[[col_sym]], na.rm = TRUE), .groups = "drop")
      } else {
        m_df <- tibble::tibble(m = stats::median(d_valid[[col_sym]], na.rm = TRUE))
      }

      p <- ggplot(d_valid, aes(x = .data[[col_sym]])) +
        geom_density(fill = color, color = "black", alpha = 0.6) +
        geom_vline(
          data = m_df,
          aes(xintercept = m),
          linetype = "dashed",
          color = "red",
          linewidth = 0.9
        ) +
        geom_text(
          data = m_df,
          aes(x = m, y = Inf, label = paste("Median:", round(m, 3))),
          vjust = 2,
          hjust = -0.1,
          size = aes_cfg$label_size,
          color = "red",
          fontface = "bold"
        ) +
        labs(title = title, x = xlab, y = "Density") +
        aes_cfg$theme_obj

      add_enc_facet(p, aes_cfg)
    }

    plot_score_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      sc_col <- intersect(c("score", "Score", "scores", "primary"), names(d))
      if (length(sc_col) == 0L) {
        return(enc_empty_plot("Score column not found in loaded data."))
      }
      plot_annotated_density_enc(
        d,
        sc_col[1],
        "Score Distribution",
        sc_col[1],
        aes_cfg$bar_color,
        aes_cfg
      )
    })

    plot_pep_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      p_col <- intersect(
        c("posterior_error_prob", "PEP", "pep", "PosteriorErrorProbability"),
        names(d)
      )
      if (length(p_col) == 0L || all(is.na(d[[p_col[1]]]))) {
        return(enc_empty_plot(
          "Posterior Error Probability (PEP) is only present in Percolator PSM outputs (*.encyclopedia*.txt)."
        ))
      }
      plot_annotated_density_enc(
        d,
        p_col[1],
        "Posterior Error Probability (PEP)",
        p_col[1],
        aes_cfg$line_color,
        aes_cfg
      )
    })

    plot_qval_obj <- reactive({
      d <- enc_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      p <- ggplot(d, aes(x = `q-value`)) +
        geom_histogram(
          bins = 50,
          fill = aes_cfg$bar_color,
          color = "black",
          alpha = 0.8
        ) +
        geom_vline(
          xintercept = input$qval_filter %||% 0.01,
          color = "red",
          linetype = "dashed",
          linewidth = 0.9
        ) +
        labs(title = "q-value Distribution", x = "q-value", y = "Count") +
        aes_cfg$theme_obj
      add_enc_facet(p, aes_cfg)
    })

    plot_charge_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      if (!"charge" %in% names(d) || all(is.na(d$charge))) {
        return(enc_empty_plot("Charge state information is missing in loaded data."))
      }
      d_c <- d |> dplyr::filter(!is.na(charge))
      p <- ggplot(d_c, aes(x = factor(charge))) +
        geom_bar(
          alpha = 0.85,
          fill = aes_cfg$bar_color,
          color = "black",
          linewidth = 0.25
        ) +
        labs(
          title = "Precursor Charge State Distribution",
          x = "Charge state",
          y = "Count"
        ) +
        aes_cfg$theme_obj
      add_enc_facet(p, aes_cfg)
    })

    plot_length_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      plot_annotated_density_enc(
        d,
        "pep_len",
        "Peptide Length Distribution",
        "Peptide length (AA)",
        aes_cfg$bar_color,
        aes_cfg
      )
    })

    plot_mods_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      if (!"n_mods" %in% names(d) || all(is.na(d$n_mods))) {
        return(enc_empty_plot("Modification annotations are missing in loaded data."))
      }
      p <- ggplot(d, aes(x = factor(n_mods))) +
        geom_bar(
          alpha = 0.85,
          fill = aes_cfg$bar_color,
          color = "black",
          linewidth = 0.25
        ) +
        labs(
          title = "Modifications per Peptide",
          x = "Number of modifications per peptide",
          y = "Count"
        ) +
        aes_cfg$theme_obj
      add_enc_facet(p, aes_cfg)
    })

    plot_gravy_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      plot_annotated_density_enc(
        d,
        "gravy",
        "GRAVY Index Distribution",
        "GRAVY Index",
        aes_cfg$line_color,
        aes_cfg
      )
    })

    plot_pi_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "pI" %in% names(d))
      aes_cfg <- enc_aesthetics()
      p <- plot_2d_gel(
        d,
        pi_col = "pI",
        mw_col = "MW",
        title = "Virtual 2D Gel (pI vs MW)",
        facet_col = if (isTRUE(aes_cfg$facet)) "filename" else NULL,
        ncol = aes_cfg$ncol
      )
      p +
        aes_cfg$theme_obj +
        theme(legend.position = "bottom")
    })

    plot_aa_freq_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      base_aas <- c(
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
        "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
      )

      groups <- if (isTRUE(aes_cfg$facet)) unique(d$filename) else "__all__"
      obs_df <- lapply(groups, function(fname) {
        sub_d <- if (identical(fname, "__all__")) d else d[d$filename == fname, , drop = FALSE]
        pep_aa <- paste0(
          sub_d$stripped_peptide[!is.na(sub_d$stripped_peptide)],
          collapse = ""
        )
        if (nchar(pep_aa) == 0L) {
          return(NULL)
        }
        tc <- table(strsplit(pep_aa, "")[[1]])
        tc <- tc[names(tc) %in% base_aas]
        if (length(tc) == 0L) {
          return(NULL)
        }
        freq <- as.data.frame(
          tc / sum(tc, na.rm = TRUE) * 100,
          stringsAsFactors = FALSE
        )
        colnames(freq) <- c("AA", "Frequency")
        missing_aas <- setdiff(base_aas, freq$AA)
        if (length(missing_aas) > 0L) {
          freq <- rbind(
            freq,
            data.frame(
              AA = missing_aas,
              Frequency = 0,
              stringsAsFactors = FALSE
            )
          )
        }
        freq$filename <- fname
        freq
      })

      obs_freq <- do.call(rbind, obs_df)
      if (is.null(obs_freq) || nrow(obs_freq) == 0L) {
        return(enc_empty_plot("Not enough sequences to compute amino acid frequencies."))
      }
      obs_freq$AA <- factor(obs_freq$AA, levels = base_aas)

      p <- ggplot(obs_freq, aes(x = AA, y = Frequency)) +
        geom_col(fill = aes_cfg$bar_color, color = "black", alpha = 0.85) +
        labs(
          title = "Amino Acid Frequencies",
          x = "Amino acid",
          y = "Frequency (%)"
        ) +
        aes_cfg$theme_obj
      add_enc_facet(p, aes_cfg)
    })

    plot_fdr_curve_obj <- reactive({
      d <- enc_data()
      req(d, nrow(d) > 0, "q-value" %in% names(d))
      aes_cfg <- enc_aesthetics()

      q_sorted <- d |>
        dplyr::group_by(filename) |>
        dplyr::arrange(`q-value`, .by_group = TRUE) |>
        dplyr::mutate(cumulative_peptides = dplyr::row_number()) |>
        dplyr::ungroup()

      if (isTRUE(aes_cfg$facet)) {
        p <- ggplot(q_sorted, aes(x = `q-value`, y = cumulative_peptides)) +
          geom_line(color = aes_cfg$bar_color, linewidth = 1) +
          labs(
            title = "Peptide Yield vs. FDR",
            x = "q-value",
            y = "Cumulative peptide count"
          )
      } else {
        p <- ggplot(
          q_sorted,
          aes(x = `q-value`, y = cumulative_peptides, color = filename)
        ) +
          geom_line(linewidth = 1) +
          labs(
            title = "Peptide Yield vs. FDR",
            x = "q-value",
            y = "Cumulative peptide count",
            color = "File"
          )
      }

      p <- p +
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
          size = aes_cfg$label_size
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
          size = aes_cfg$label_size
        ) +
        scale_x_continuous(labels = scales::label_scientific()) +
        aes_cfg$theme_obj

      add_enc_facet(p, aes_cfg)
    })

    # ── Feature-table (*.mzML.features.txt) QC plots ─────────────────────────
    feature_missing_msg <- paste(
      "Feature metrics require *.features.txt.",
      "Select 'Auto' or 'Feature tables only' under Input tables and click Load Directory."
    )

    plot_mz_rt_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      if (!all(c("RTinMin", "precursorMz") %in% names(d))) {
        return(enc_empty_plot(feature_missing_msg))
      }
      dv <- d[is.finite(d$RTinMin) & is.finite(d$precursorMz), , drop = FALSE]
      if (nrow(dv) == 0L) {
        return(enc_empty_plot("No finite RT / precursor m/z values available."))
      }
      if (isTRUE(aes_cfg$facet)) {
        dv$.dens <- stats::ave(
          seq_len(nrow(dv)),
          dv$filename,
          FUN = function(i) point_density_2d(dv$RTinMin[i], dv$precursorMz[i])
        )
      } else {
        dv$.dens <- point_density_2d(dv$RTinMin, dv$precursorMz)
      }
      p <- ggplot(
        dv,
        aes(x = RTinMin, y = precursorMz, color = .data[[".dens"]])
      ) +
        geom_point(size = 0.8, alpha = 0.75) +
        scale_color_viridis_c(name = "Ion density", option = "D") +
        labs(
          title = "Precursor m/z vs Retention Time",
          x = "Retention time (min)",
          y = "Precursor m/z"
        ) +
        aes_cfg$theme_obj +
        theme(legend.position = "bottom")
      add_enc_facet(p, aes_cfg)
    })

    plot_delta_rt_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      if (!"deltaRT" %in% names(d) || all(is.na(d$deltaRT))) {
        return(enc_empty_plot(feature_missing_msg))
      }
      plot_annotated_density_enc(
        d,
        "deltaRT",
        "Retention Time Error (ΔRT)",
        "ΔRT (min)",
        aes_cfg$bar_color,
        aes_cfg
      )
    })

    plot_mass_error_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      err_cols <- intersect(c("ms1MassError", "ms2MassError"), names(d))
      if (length(err_cols) == 0L) {
        return(enc_empty_plot(feature_missing_msg))
      }
      long_err <- d |>
        dplyr::select(filename, dplyr::all_of(err_cols)) |>
        tidyr::pivot_longer(
          cols = dplyr::all_of(err_cols),
          names_to = "Level",
          values_to = "ppm_error"
        ) |>
        dplyr::filter(is.finite(ppm_error)) |>
        dplyr::mutate(
          Level = dplyr::recode(
            Level,
            ms1MassError = "MS1 Mass Error",
            ms2MassError = "MS2 Mass Error"
          )
        )
      if (nrow(long_err) == 0L) {
        return(enc_empty_plot("No finite mass error observations available."))
      }
      p <- ggplot(long_err, aes(x = ppm_error, fill = Level)) +
        geom_density(alpha = 0.55, color = "black", linewidth = 0.3) +
        scale_fill_manual(
          values = c(
            "MS1 Mass Error" = aes_cfg$bar_color,
            "MS2 Mass Error" = aes_cfg$line_color
          ),
          name = NULL
        ) +
        labs(
          title = "MS1 & MS2 Mass Error Distribution",
          x = "Absolute mass error (ppm)",
          y = "Density"
        ) +
        aes_cfg$theme_obj
      add_enc_facet(p, aes_cfg)
    })

    plot_missed_cleavages_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      if (!"numMissedCleavage" %in% names(d) || all(is.na(d$numMissedCleavage))) {
        return(enc_empty_plot(feature_missing_msg))
      }
      dm <- d |> dplyr::filter(!is.na(numMissedCleavage))
      p <- ggplot(dm, aes(x = factor(numMissedCleavage))) +
        geom_bar(
          fill = aes_cfg$bar_color,
          color = "black",
          alpha = 0.85,
          linewidth = 0.25
        ) +
        labs(
          title = "Missed Cleavages per Peptide",
          x = "Number of missed cleavages",
          y = "Count"
        ) +
        aes_cfg$theme_obj
      add_enc_facet(p, aes_cfg)
    })

    plot_feature_scores_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      sc_cols <- intersect(c("scribe", "HyperScore", "xCorrLib"), names(d))
      if (length(sc_cols) == 0L) {
        return(enc_empty_plot(feature_missing_msg))
      }
      long_sc <- d |>
        dplyr::select(filename, dplyr::all_of(sc_cols)) |>
        tidyr::pivot_longer(
          cols = dplyr::all_of(sc_cols),
          names_to = "Metric",
          values_to = "Value"
        ) |>
        dplyr::filter(is.finite(Value))
      if (nrow(long_sc) == 0L) {
        return(enc_empty_plot("No finite feature score values available."))
      }
      if (isTRUE(aes_cfg$facet)) {
        ggplot(long_sc, aes(x = filename, y = Value, fill = Metric)) +
          geom_boxplot(outlier.size = 0.4, alpha = 0.8) +
          facet_wrap(~Metric, scales = "free_y", ncol = min(3L, length(sc_cols))) +
          labs(
            title = "EncyclopeDIA Feature Scores by File",
            x = "File",
            y = "Score"
          ) +
          aes_cfg$theme_obj +
          theme(
            axis.text.x = element_text(angle = 55, hjust = 1, color = "black"),
            legend.position = "none"
          )
      } else {
        ggplot(long_sc, aes(x = Value, fill = Metric)) +
          geom_density(alpha = 0.65, color = "black") +
          facet_wrap(~Metric, scales = "free", ncol = min(3L, length(sc_cols))) +
          labs(
            title = "EncyclopeDIA Feature Score Distributions",
            x = "Score",
            y = "Density"
          ) +
          aes_cfg$theme_obj +
          theme(legend.position = "none")
      }
    })

    plot_idotp_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      aes_cfg <- enc_aesthetics()
      if (!"isotopeDotProduct" %in% names(d) || all(is.na(d$isotopeDotProduct))) {
        return(enc_empty_plot(feature_missing_msg))
      }
      plot_annotated_density_enc(
        d,
        "isotopeDotProduct",
        "Isotope Dot Product Distribution",
        "Isotope dot product",
        aes_cfg$bar_color,
        aes_cfg
      )
    })

    # ── Dynamic plot height & rendering ──────────────────────────────────────
    enc_plot_dims <- reactive({
      d <- enc_data()
      req(d)
      aes_cfg <- enc_aesthetics()
      pl <- active_plot_id()
      single_panel <- pl %in% c("plot_identifications", "plot_feature_scores") ||
        !isTRUE(aes_cfg$facet)
      n_files <- max(1L, dplyr::n_distinct(d$filename))
      row_h <- max(280L, as.integer(round(aes_cfg$font_size * 24)))
      if (single_panel) {
        list(px = max(440L, as.integer(round(row_h * 1.45))), w_in = 11, h_in = 7)
      } else {
        eff_ncol <- min(aes_cfg$ncol, n_files)
        eff_nrow <- ceiling(n_files / aes_cfg$ncol)
        list(
          px = max(420L, eff_nrow * row_h),
          w_in = max(8, eff_ncol * 3.8),
          h_in = max(6, eff_nrow * 3.2 + 1.2)
        )
      }
    })

    output$dynamic_plot_ui <- renderUI({
      if (is.null(enc_data())) {
        return(div(
          class = "plot-placeholder",
          "Load an EncyclopeDIA directory from the sidebar to view QC plots."
        ))
      }
      dims <- enc_plot_dims()
      plotOutput(ns("dynamic_plot_out"), height = paste0(dims$px, "px"))
    })

    current_plot_obj <- reactive({
      req(enc_data(), active_plot_id())
      pl <- active_plot_id()
      if (
        !pl %in% c("plot_qval", "plot_fdr_curve") &&
          nrow(filtered_data()) == 0L
      ) {
        return(empty_filtered_plot())
      }
      switch(
        pl,
        "plot_identifications" = plot_identifications_obj(),
        "plot_score" = plot_score_obj(),
        "plot_pep" = plot_pep_obj(),
        "plot_qval" = plot_qval_obj(),
        "plot_fdr_curve" = plot_fdr_curve_obj(),
        "plot_charge" = plot_charge_obj(),
        "plot_mods" = plot_mods_obj(),
        "plot_length" = plot_length_obj(),
        "plot_gravy" = plot_gravy_obj(),
        "plot_pi" = plot_pi_obj(),
        "plot_aa_freq" = plot_aa_freq_obj(),
        "plot_mz_rt" = plot_mz_rt_obj(),
        "plot_delta_rt" = plot_delta_rt_obj(),
        "plot_mass_error" = plot_mass_error_obj(),
        "plot_missed_cleavages" = plot_missed_cleavages_obj(),
        "plot_feature_scores" = plot_feature_scores_obj(),
        "plot_idotp" = plot_idotp_obj()
      )
    })

    output$dynamic_plot_out <- renderPlot({
      shinyjs::show(id = "spi_main")
      on.exit(hide_sp("spi_main"), add = TRUE)
      req(current_plot_obj())
      current_plot_obj()
    })

    output$res_table <- DT::renderDataTable({
      req(mapped_data())
      d <- mapped_data()
      num_cols <- names(d)[vapply(d, is.numeric, logical(1))]
      round_cols <- intersect(
        num_cols,
        c(
          "q-value",
          "score",
          "Score",
          "posterior_error_prob",
          "PEP",
          "gravy",
          "pI",
          "MW",
          "precursorMz",
          "precursorMass",
          "RTinMin",
          "deltaRT",
          "ms1MassError",
          "ms2MassError",
          "scribe",
          "HyperScore",
          "xCorrLib",
          "isotopeDotProduct"
        )
      )

      dt <- DT::datatable(
        d,
        extensions = "Buttons",
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ),
        class = "cell-border stripe",
        rownames = FALSE,
        filter = "top"
      )
      if (length(round_cols) > 0L) {
        dt <- dt |> DT::formatRound(columns = round_cols, digits = 4)
      }
      dt
    })

    # ── Plot Download Button ─────────────────────────────────────────────────
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("EncyclopeDIA_", active_plot_id(), "_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(current_plot_obj())
        dims <- enc_plot_dims()
        ggsave(
          file,
          current_plot_obj(),
          width = dims$w_in,
          height = dims$h_in,
          dpi = 300,
          bg = "white",
          device = "png"
        )
      }
    )
  })
}
