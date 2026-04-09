## ============================================================
## dash_deNovo.r  —  Casanovo de novo results viewer
##
## Architecture:
##   deNovo_sidebar_ui(id)  → sidebar controls
##   deNovo_body_ui(id)     → tab panels + plots
##   deNovo_ui(id)          → tagList(sidebar, body)  [for convenience]
##   deNovo_server(id)      → module server
##
## Standalone: run as shinyApp at the bottom.
## Embed in proteOmni: source this file, then call
##   deNovo_sidebar_ui("dnv")  in the sidebar
##   deNovo_body_ui("dnv")     in the body
##   deNovo_server("dnv")      in the server
## ============================================================

# ── Libraries ────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyjs)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggseqlogo)
  library(readr)
  library(tibble)
  library(ggpointdensity)
  library(viridis)
  library(colourpicker)
  library(DT)
})

# ── Helper functions ──────────────────────────────────────────────────────────

## Strip modification annotations from the sequence.
## Handles:
##   "[Carbamyl]-ELDDLEALVERR" → "ELDDLEALVERR"
##   "THM+15.995ELGGK"         → "THMELGGK"
##   "WGEAGAEYVVESTGVFTTM[Oxidation]EK" → "WGEAGAEYVVESTGVFTTMEK"
strip_sequence <- function(seq) {
  s <- as.character(seq)
  s <- sub("^\\[.*?\\]-?", "", s)
  s <- gsub("\\[.*?\\]", "", s)
  s <- gsub("[+-][0-9]+\\.?[0-9]*", "", s)
  s <- gsub("[^A-Z]", "", s)
  s
}

## Clean modification name from mzTab format:
## "19-Oxidation (M):UNIMOD:35" → "Oxidation (M)"
## "null" / NA                  → NA
parse_mod_name <- function(x) {
  if (is.na(x) || str_trim(x) %in% c("null", "NULL", "")) {
    return(NA_character_)
  }
  parts <- str_split(x, ",")[[1]]
  parts <- str_trim(parts)
  parts <- sub("^[0-9]+-", "", parts) # strip position prefix
  parts <- sub(":UNIMOD:.*$", "", parts) # strip UNIMOD suffix
  paste(parts, collapse = "; ")
}


## Parse a Casanovo .mztab file → list(metadata df, psm df)
read_mztab <- function(path) {
  lines <- readLines(path, warn = FALSE)
  # Strip carriage returns (Windows CRLF files)
  lines <- gsub("\r", "", lines)

  # Metadata
  mtd <- lines[startsWith(lines, "MTD")]
  meta <- do.call(rbind, strsplit(mtd, "\t")) |>
    as.data.frame(stringsAsFactors = FALSE)
  colnames(meta) <- c("prefix", "key", "value")
  meta <- meta[, c("key", "value")]

  # PSM header
  psh <- lines[startsWith(lines, "PSH")][1]
  cols <- strsplit(psh, "\t")[[1]]
  cols[1] <- "row_type"

  # PSM rows
  psm_lines <- lines[startsWith(lines, "PSM")]
  if (!length(psm_lines)) {
    return(list(metadata = meta, psm = NULL))
  }

  mat <- do.call(rbind, strsplit(psm_lines, "\t"))
  if (ncol(mat) < length(cols)) {
    mat <- cbind(mat, matrix("", nrow(mat), length(cols) - ncol(mat)))
  }
  psm <- as.data.frame(mat[, seq_len(length(cols))], stringsAsFactors = FALSE)
  colnames(psm) <- cols
  psm$row_type <- NULL

  num_cols <- c(
    "PSM_ID",
    "search_engine_score[1]",
    "retention_time",
    "charge",
    "exp_mass_to_charge",
    "calc_mass_to_charge"
  )
  for (col in intersect(num_cols, names(psm))) {
    psm[[col]] <- suppressWarnings(as.numeric(psm[[col]]))
  }

  psm <- psm |>
    mutate(across(
      where(is.character),
      ~ ifelse(. %in% c("null", "NULL", ""), NA, .)
    ))

  score_col <- which(names(psm) == "search_engine_score[1]")
  if (length(score_col) == 0) {
    score_col <- grep("search_engine_score", names(psm), fixed = FALSE)[1]
  }
  if (length(score_col) > 0 && !is.na(score_col)) {
    names(psm)[score_col] <- "score"
  } else {
    stop(paste(
      "Cannot find search_engine_score column. Available columns:",
      paste(names(psm), collapse = ", ")
    ))
  }

  psm <- psm |>
    mutate(
      stripped_sequence = sapply(sequence, strip_sequence, USE.NAMES = FALSE),
      peptide_length = nchar(stripped_sequence),
      mz_error = exp_mass_to_charge - calc_mass_to_charge,
      mod_name = sapply(modifications, parse_mod_name, USE.NAMES = FALSE),
      is_modified = !is.na(mod_name),
      aa_score_mean = sapply(
        `opt_ms_run[1]_aa_scores`,
        function(x) {
          if (is.na(x)) {
            return(NA_real_)
          }
          sc <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
          mean(sc, na.rm = TRUE)
        },
        USE.NAMES = FALSE
      ),
      gravy = sapply(stripped_sequence, GRAVY, USE.NAMES = FALSE),
      pI = sapply(stripped_sequence, calculate_pI, USE.NAMES = FALSE)
    )

  list(metadata = meta, psm = psm)
}

prep_seqlogo <- function(seqs, width) {
  seqs <- seqs[!is.na(seqs) & nchar(seqs) >= width]
  seqs <- substr(seqs, 1, width)
  seqs[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", seqs)]
}

theme_dn <- function(...) {
  theme_bw(...) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.title.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold"),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

## Spinner overlay wrapper
sp_wrap <- function(sid, ui_el) {
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

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — SIDEBAR UI
# ═══════════════════════════════════════════════════════════════════════════════
deNovo_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      style = "padding:12px 16px 4px;color:#ffffff;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      icon("fingerprint"),
      " Casanovo de Novo"
    ),
    tags$div(
      style = "padding: 0 16px; margin-bottom: 12px;",
      textInput(
        ns("mztab_dir"),
        "Path to results directory:",
        placeholder = "/path/to/mztab_files/"
      ),
      actionButton(
        ns("load_dir_btn"),
        "Load Directory",
        icon = icon("folder-open"),
        class = "btn-primary",
        style = "width:80%;"
      )
    ),
    sliderInput(
      ns("score_filter"),
      "Min Casanovo score",
      min = -1,
      max = 1,
      value = -1,
      step = 0.01
    ),
    sliderInput(
      ns("aa_score_filter"),
      "Min mean AA score",
      min = 0,
      max = 1,
      value = 0,
      step = 0.01
    ),
    colourpicker::colourInput(
      ns("plot_color"),
      "Plot colour",
      value = "#2980b9"
    ),
    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),
    tags$div(
      style = "padding:4px 16px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      "SeqLogo options"
    ),
    numericInput(
      ns("nterm_width"),
      "N-terminus width (AA)",
      value = 5,
      min = 2,
      max = 15
    ),
    numericInput(
      ns("cterm_width"),
      "C-terminus width (AA)",
      value = 5,
      min = 2,
      max = 15
    ),
    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),
    div(
      style = "padding:0 8px;",
      downloadButton(
        ns("download_plots"),
        "⬇ Download plots",
        class = "dl-btn",
        style = "width:100%;text-align:left;"
      )
    )
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — BODY UI
# ═══════════════════════════════════════════════════════════════════════════════
deNovo_body_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      id = ns("tabs"),
      type = "tabs",

      # ── Overview ─────────────────────────────────────────────────────────────
      tabPanel(
        "Overview",
        fluidRow(column(
          12,
          div(
            id = ns("meta_box"),
            style = "background:#1a2332;border:1px solid #2d3741;border-radius:6px;padding:16px;margin:12px 0;",
            h5(
              style = "color:#7fb3d5;margin:0 0 8px 0;font-weight:700;",
              icon("circle-info"),
              " Casanovo Run Summary"
            ),
            uiOutput(ns("metadata_ui"))
          )
        )),
        fluidRow(infoBoxOutput(ns("info_box"), width = 12)),
        fluidRow(
          box(
            title = "Score distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp01"), uiOutput(ns("plot_score_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Mean AA-score distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp02"), uiOutput(ns("plot_aa_score_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Peptide length distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp03"), uiOutput(ns("plot_length_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Charge state distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp04"), uiOutput(ns("plot_charge_ui")))
          )
        ),
        fluidRow(
          box(
            title = "m/z error (exp − calc)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp05"), uiOutput(ns("plot_mz_error_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Score vs mean AA score",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp06"), uiOutput(ns("plot_score_vs_aa_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Number of retained PSMs vs. score threshold",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp17"), uiOutput(ns("plot_elbow_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Retention time vs. m/z error",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp18"), uiOutput(ns("plot_rt_mz_error_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Retention time vs. PSM count",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp19"), uiOutput(ns("plot_rt_density_ui")))
          )
        )
      ),

      # ── Peptide Analysis ─────────────────────────────────────────────────────
      tabPanel(
        "Peptide Analysis",
        fluidRow(
          box(
            title = "Modifications summary",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp07"), uiOutput(ns("plot_mods_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Modified vs Unmodified PSMs",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp08"), uiOutput(ns("plot_mod_pie_ui")))
          )
        ),
        fluidRow(
          box(
            title = "GRAVY (Grand Average of Hydropathy)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp09"), uiOutput(ns("plot_gravy_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Isoelectric Point (pI) distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp15"), uiOutput(ns("plot_pi_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Median score by peptide length",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp10"), uiOutput(ns("plot_score_length_ui")))
          )
        ),
        fluidRow(
          box(
            title = "m/z error vs Casanovo score",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp11"), uiOutput(ns("plot_mz_score_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Amino acid frequencies (Identified peptides)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp16"), uiOutput(ns("plot_aa_freq_ui")))
          )
        )
      ),

      # ── Termini SeqLogos ──────────────────────────────────────────────────────
      tabPanel(
        "Termini SeqLogos",
        fluidRow(
          box(
            title = "N-terminus SeqLogo (first N AA from N-terminal)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp12"), plotOutput(ns("plot_nterm"), height = "280px"))
          ),
          box(
            title = "C-terminus SeqLogo (last N AA from C-terminal)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(ns("sp13"), plotOutput(ns("plot_cterm"), height = "280px"))
          )
        ),
        fluidRow(
          box(
            title = "N:C-terminus co-occurrence probability",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap(
              ns("sp14"),
              uiOutput(ns("plot_cooccurrence_ui"))
            )
          )
        )
      ),

      # ── PSM Table ────────────────────────────────────────────────────────────
      tabPanel(
        "PSM Table",
        fluidRow(column(
          12,
          div(
            style = "margin:12px 0;",
            DT::dataTableOutput(ns("psm_table"))
          )
        ))
      )
    )
  )
}

## Convenience wrapper (sidebar + body together) for minimal integrations
deNovo_ui <- function(id) {
  tagList(deNovo_sidebar_ui(id), deNovo_body_ui(id))
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — SERVER
# ═══════════════════════════════════════════════════════════════════════════════
deNovo_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    spin_ids <- paste0("sp", sprintf("%02d", 1:16))
    show_all <- function() lapply(spin_ids, function(s) shinyjs::show(id = s))
    hide_sp <- function(s) shinyjs::hide(id = s)
    rh <- function(fn, sid) {
      on.exit(hide_sp(sid), add = TRUE)
      fn()
    }

    # ── Load raw data ───────────────────────────────────────────────────────────
    raw_data_rv <- reactiveValues(metadata = NULL, psm = NULL)

    observeEvent(input$load_dir_btn, {
      dir_path <- trimws(input$mztab_dir)
      if (dir_path == "" || !dir.exists(dir_path)) {
        showNotification("Please enter a valid directory path.", type = "error")
        return()
      }
      show_all()
      withProgress(message = "Parsing mzTab files...", value = 0, {
        mz_files <- list.files(
          dir_path,
          pattern = "\\.mztab$",
          full.names = TRUE,
          recursive = TRUE
        )
        if (length(mz_files) == 0) {
          showNotification(
            "No .mztab files found in the directory.",
            type = "warning"
          )
          setProgress(1)
          hide_sp("spall")
          return()
        }

        all_meta <- psm_list <- list()
        n <- length(mz_files)

        for (i in seq_along(mz_files)) {
          fpath <- mz_files[i]
          fname <- str_extract(basename(fpath), "^[^\\.]+")
          setProgress(i / n, detail = paste("Reading", basename(fpath)))

          res <- tryCatch(read_mztab(fpath), error = function(e) NULL)
          if (!is.null(res)) {
            if (!is.null(res$psm) && nrow(res$psm) > 0) {
              res$psm$filename <- fname
              psm_list[[fname]] <- res$psm
            }
            if (!is.null(res$metadata)) {
              res$metadata$filename <- fname
              all_meta[[fname]] <- res$metadata
            }
          }
        }

        if (length(psm_list) > 0) {
          raw_data_rv$psm <- bind_rows(psm_list)
          raw_data_rv$metadata <- bind_rows(all_meta)
        } else {
          raw_data_rv$psm <- NULL
          raw_data_rv$metadata <- NULL
        }
        setProgress(1)
      })
    })

    raw_data <- reactive({
      req(raw_data_rv$psm)
      list(metadata = raw_data_rv$metadata, psm = raw_data_rv$psm)
    })

    # ── Update slider ranges once file is loaded ────────────────────────────────
    observeEvent(raw_data_rv$psm, {
      psm_dat <- raw_data_rv$psm
      sc <- psm_dat$score[is.finite(psm_dat$score)]
      if (length(sc)) {
        sc_lo <- floor(min(sc) * 100) / 100
        sc_hi <- ceiling(max(sc) * 100) / 100
        updateSliderInput(
          session,
          "score_filter",
          min = sc_lo,
          max = sc_hi,
          value = sc_lo
        )
      }
      aa <- psm_dat$aa_score_mean[is.finite(psm_dat$aa_score_mean)]
      if (length(aa)) {
        updateSliderInput(
          session,
          "aa_score_filter",
          max = ceiling(max(aa) * 100) / 100
        )
      }
    })

    # ── Filtered PSM data ───────────────────────────────────────────────────────
    data <- reactive({
      rd <- raw_data()
      req(rd, rd$psm)
      rd$psm |>
        filter(
          is.na(score) | score >= input$score_filter,
          is.na(aa_score_mean) | aa_score_mean >= input$aa_score_filter
        )
    })

    meta <- reactive({
      rd <- raw_data()
      req(rd)
      rd$metadata
    })

    # ── Metadata box UI ─────────────────────────────────────────────────────────
    output$metadata_ui <- renderUI({
      m <- meta()
      req(m)
      m <- m |> filter(filename == unique(m$filename)[1])
      key_map <- c(
        "mzTab-version" = "mzTab version",
        "mzTab-mode" = "Mode",
        "description" = "Description",
        "software[1]" = "Software",
        "psm_search_engine_score[1]" = "Score type",
        "fixed_mod[1]" = "Fixed mods",
        "variable_mod[1]" = "Variable mods"
      )
      badge_style <- "background:#243447;border:1px solid #2d3741;border-radius:4px;padding:4px 10px;font-size:11px;color:#e0e0e0;margin:3px;"
      badges <- lapply(names(key_map), function(k) {
        v <- m$value[m$key == k]
        if (!length(v) || is.na(v) || v == "") {
          return(NULL)
        }
        tags$span(
          style = badge_style,
          tags$b(style = "color:#7fb3d5;", key_map[k], ": "),
          v
        )
      })

      sk <- m$key[startsWith(m$key, "software[1]-setting")]
      sv <- m$value[startsWith(m$key, "software[1]-setting")]
      settings_el <- if (length(sk)) {
        tags$details(
          tags$summary(
            style = "color:#7fb3d5;cursor:pointer;font-size:11px;font-weight:700;margin-top:8px;",
            paste0("▶ Model settings (", length(sk), " parameters)")
          ),
          tags$table(
            style = "font-size:10px;color:#adb5bd;margin-top:6px;border-collapse:collapse;width:100%;",
            lapply(seq_along(sk), function(i) {
              tags$tr(
                style = if (i %% 2 == 0) "background:#1e2d3d;" else "",
                tags$td(
                  style = "padding:2px 8px;color:#7fb3d5;font-weight:bold;min-width:220px;",
                  sub("\\s*=.*$", "", sv[i])
                ),
                tags$td(
                  style = "padding:2px 8px;",
                  sub("^[^=]+=\\s*", "", sv[i])
                )
              )
            })
          )
        )
      } else {
        NULL
      }

      tagList(div(style = "display:flex;flex-wrap:wrap;", badges), settings_el)
    })

    # ── Info box ────────────────────────────────────────────────────────────────
    output$info_box <- renderInfoBox({
      d <- data()
      total <- nrow(raw_data()$psm)
      n <- nrow(d)
      pct <- if (total > 0) round(n / total * 100, 1) else 0
      infoBox(
        "PSM Filter",
        paste0(
          n,
          " / ",
          total,
          " PSMs  (",
          pct,
          "%)  |  score ≥ ",
          input$score_filter,
          "  |  mean AA score ≥ ",
          input$aa_score_filter
        ),
        icon = icon("filter"),
        color = "black"
      )
    })

    # ── PLOTS ───────────────────────────────────────────────────────────────────

    # ── Dynamic plot height ───────────────────────────────────────────────
    dn_plot_h <- reactive({
      d <- data()
      n <- dplyr::n_distinct(d$filename)
      max(400L, ceiling(n / 3L) * 350L)
    })

    make_dn_plot_ui <- function(plot_id) {
      renderUI({ req(data()); plotOutput(ns(plot_id), height = paste0(dn_plot_h(), "px")) })
    }

    output$plot_score_ui        <- make_dn_plot_ui("plot_score")
    output$plot_aa_score_ui     <- make_dn_plot_ui("plot_aa_score")
    output$plot_length_ui       <- make_dn_plot_ui("plot_length")
    output$plot_charge_ui       <- make_dn_plot_ui("plot_charge")
    output$plot_mz_error_ui     <- make_dn_plot_ui("plot_mz_error")
    output$plot_score_vs_aa_ui  <- make_dn_plot_ui("plot_score_vs_aa")
    output$plot_elbow_ui        <- make_dn_plot_ui("plot_elbow")
    output$plot_rt_mz_error_ui  <- make_dn_plot_ui("plot_rt_mz_error")
    output$plot_rt_density_ui   <- make_dn_plot_ui("plot_rt_density")

    output$plot_mods_ui         <- make_dn_plot_ui("plot_mods")
    output$plot_mod_pie_ui      <- make_dn_plot_ui("plot_mod_pie")
    output$plot_gravy_ui        <- make_dn_plot_ui("plot_gravy")
    output$plot_pi_ui           <- make_dn_plot_ui("plot_pi")
    output$plot_score_length_ui <- make_dn_plot_ui("plot_score_length")
    output$plot_mz_score_ui     <- make_dn_plot_ui("plot_mz_score")
    output$plot_aa_freq_ui      <- make_dn_plot_ui("plot_aa_freq")
    output$plot_cooccurrence_ui <- make_dn_plot_ui("plot_cooccurrence")

    output$plot_score <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          ggplot(d, aes(x = score)) +
            geom_histogram(
              bins = 80,
              fill = input$plot_color,
              color = "black",
              alpha = 0.85
            ) +
            geom_vline(
              xintercept = input$score_filter,
              color = "red",
              linetype = "dashed",
              linewidth = 0.8
            ) +
            labs(
              title = "Casanovo score distribution",
              x = "Casanovo Score",
              y = "Count",
              caption = "Red line = current filter threshold"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn()
        },
        "sp01"
      )
    })

    output$plot_aa_score <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          med_aa <- round(median(d$aa_score_mean, na.rm = TRUE), 3)
          q10 <- round(quantile(d$aa_score_mean, 0.10, na.rm = TRUE), 3)
          interp <- dplyr::case_when(
            med_aa >= 0.9 ~ "high-confidence (most residues well-predicted)",
            med_aa >= 0.7 ~ "moderate confidence (some residues uncertain)",
            TRUE ~ "low confidence (many ambiguous residue predictions)"
          )
          ggplot(d, aes(x = aa_score_mean)) +
            geom_histogram(
              bins = 60,
              fill = input$plot_color,
              color = "black",
              alpha = 0.85
            ) +
            geom_vline(
              xintercept = input$aa_score_filter,
              color = "red",
              linetype = "dashed",
              linewidth = 0.8
            ) +
            geom_vline(
              xintercept = med_aa,
              color = "orange",
              linetype = "dotted",
              linewidth = 0.8
            ) +
            annotate(
              "text",
              x = med_aa,
              y = Inf,
              vjust = 1.5,
              hjust = -0.1,
              label = paste0("median=", med_aa),
              color = "orange",
              size = 3.5,
              fontface = "bold"
            ) +
            labs(
              title = "Mean per-residue AA score",
              x = "Mean AA score (per PSM)",
              y = "Count",
              caption = paste0(
                "Each PSM receives one score per predicted amino acid (opt_ms_run[1]_aa_scores); ",
                "this plot shows the mean across residues in each peptide.\n",
                "Interpretation: scores near 1 = high model confidence; scores near 0 = uncertain prediction.\n",
                "Median = ",
                med_aa,
                " → ",
                interp,
                ".  ",
                "10th percentile = ",
                q10,
                " (filter threshold removes the lowest-scoring PSMs)."
              )
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn() +
            theme(
              plot.caption = element_text(
                size = 10,
                color = "gray40",
                hjust = 0
              )
            )
        },
        "sp02"
      )
    })

    output$plot_length <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          len_range <- range(d$peptide_length, na.rm = TRUE)
          ggplot(d, aes(x = peptide_length)) +
            geom_histogram(
              bins = len_range[2] - len_range[1] + 1,
              fill = input$plot_color,
              color = "black",
              alpha = 0.85
            ) +
            scale_x_continuous(breaks = seq(0, 100, 5)) +
            labs(
              title = "Peptide length distribution",
              x = "Length (AA)",
              y = "Count"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn()
        },
        "sp03"
      )
    })

    output$plot_charge <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          d |>
            dplyr::count(filename, charge) |>
            mutate(charge = factor(charge)) |>
            ggplot(aes(x = charge, y = n, fill = charge)) +
            geom_col(color = "black", show.legend = FALSE) +
            geom_text(
              aes(label = n),
              vjust = -0.4,
              fontface = "bold",
              size = 4.5
            ) +
            scale_fill_viridis_d(option = "D") +
            labs(
              title = "Charge state distribution",
              x = "Charge (z)",
              y = "Count"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn()
        },
        "sp04"
      )
    })

    output$plot_mz_error <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          ggplot(d, aes(x = mz_error)) +
            geom_histogram(
              bins = 80,
              fill = input$plot_color,
              color = "black",
              alpha = 0.85
            ) +
            geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
            labs(
              title = "m/z error distribution",
              x = "m/z error  (exp m/z − calc m/z)",
              y = "Count",
              caption = "Distribution centred at 0 indicates good mass accuracy"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn()
        },
        "sp05"
      )
    })

    output$plot_score_vs_aa <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          ggplot(d, aes(x = score, y = aa_score_mean)) +
            ggpointdensity::geom_pointdensity(size = 0.4, alpha = 0.6) +
            viridis::scale_color_viridis(option = "plasma", name = "Density") +
            geom_smooth(
              method = "lm",
              se = FALSE,
              color = "red",
              linetype = "dashed",
              linewidth = 0.8
            ) +
            labs(
              title = "Casanovo score vs mean AA score",
              x = "Casanovo score",
              y = "Mean per-residue AA score"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn() +
            theme(
              legend.position = "bottom",
              legend.key.width = unit(1.5, "cm"),
              legend.key.height = unit(0.25, "cm")
            )
        },
        "sp06"
      )
    })

    output$plot_elbow <- renderPlot({
      rh(
        function() {
          rd <- raw_data()
          req(rd, rd$psm)
          d_raw <- rd$psm |> dplyr::filter(!is.na(score))
          if (nrow(d_raw) < 10) {
            return(
              ggplot() +
                theme_void()
            )
          }

          df_list <- lapply(unique(d_raw$filename), function(fname) {
            sub_d <- d_raw[d_raw$filename == fname, ]
            if (nrow(sub_d) == 0) {
              return(NULL)
            }
            scores <- sort(sub_d$score, decreasing = TRUE)
            df_c <- data.frame(
              Score_Threshold = scores,
              PSMs_Retained = seq_along(scores),
              filename = fname
            )
            if (nrow(df_c) > 500) {
              idx <- round(seq(1, nrow(df_c), length.out = 500))
              df_c <- df_c[idx, ]
            }
            df_c
          })
          df_curve <- do.call(rbind, df_list)
          if (is.null(df_curve) || nrow(df_curve) == 0) {
            return(
              ggplot() +
                theme_void()
            )
          }

          ggplot(df_curve, aes(x = Score_Threshold, y = PSMs_Retained)) +
            geom_line(color = input$plot_color, linewidth = 1.2) +
            geom_vline(
              xintercept = input$score_filter,
              color = "red",
              linetype = "dashed",
              linewidth = 1
            ) +
            labs(
              title = "Number of retained PSMs vs. Casanovo score threshold",
              x = "Minimum Casanovo score threshold",
              y = "Number of retained PSMs",
              caption = "Red dashed line marks the current score filter threshold."
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn()
        },
        "sp17"
      )
    })

    output$plot_rt_mz_error <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          if (!"retention_time" %in% names(d) || all(is.na(d$retention_time))) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "Retention time not available",
                  color = "gray40"
                ) +
                theme_void()
            )
          }
          ggplot(d, aes(x = retention_time, y = mz_error)) +
            ggpointdensity::geom_pointdensity(size = 0.5, alpha = 0.5) +
            viridis::scale_color_viridis(option = "plasma", name = "Density") +
            geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
            geom_smooth(
              method = "gam",
              formula = y ~ s(x, bs = "cs"),
              se = TRUE,
              color = "white",
              linetype = "dashed",
              linewidth = 0.8
            ) +
            labs(
              title = "Retention time vs. m/z error",
              x = "Retention time",
              y = "m/z error (exp − calc)"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn() +
            theme(
              legend.position = "bottom",
              legend.key.width = unit(1.5, "cm")
            )
        },
        "sp18"
      )
    })

    output$plot_rt_density <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          if (!"retention_time" %in% names(d) || all(is.na(d$retention_time))) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "Retention time not available",
                  color = "gray40"
                ) +
                theme_void()
            )
          }
          ggplot(d, aes(x = retention_time)) +
            geom_histogram(
              bins = 80,
              fill = input$plot_color,
              color = "black",
              alpha = 0.85
            ) +
            labs(
              title = "PSM count along LC gradient",
              x = "Retention time",
              y = "Count"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn()
        },
        "sp19"
      )
    })

    output$plot_mods <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          mod_df <- d |>
            dplyr::filter(!is.na(mod_name)) |>
            tidyr::separate_rows(mod_name, sep = "; ") |>
            dplyr::count(filename, mod_name)
          if (!nrow(mod_df)) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "No modifications in filtered PSMs",
                  size = 5,
                  color = "gray50"
                ) +
                theme_void()
            )
          }

          ggplot(mod_df, aes(x = n, y = reorder(mod_name, n))) +
            geom_col(fill = input$plot_color, color = "black") +
            geom_text(
              aes(label = n),
              hjust = -0.1,
              fontface = "bold",
              size = 4
            ) +
            scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
            labs(title = "Modification summary", x = "Count", y = NULL) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn()
        },
        "sp07"
      )
    })

    output$plot_mod_pie <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          d |>
            dplyr::mutate(
              status = ifelse(is_modified, "Modified", "Unmodified")
            ) |>
            dplyr::count(filename, status) |>
            dplyr::group_by(filename) |>
            dplyr::mutate(
              pct = round(n / sum(n) * 100, 1),
              lbl = paste0(status, "\n", n, " (", pct, "%)")
            ) |>
            ggplot(aes(x = "", y = n, fill = status)) +
            geom_col(color = "black", width = 0.5) +
            geom_text(
              aes(label = lbl),
              position = position_stack(vjust = 0.5),
              fontface = "bold",
              size = 5,
              color = "white"
            ) +
            scale_fill_manual(
              values = c(Modified = "#e67e22", Unmodified = "#2980b9")
            ) +
            coord_polar("y") +
            labs(title = "Modified vs unmodified PSMs") +
            facet_wrap(~filename, ncol = 3) +
            theme_void() +
            theme(
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              legend.position = "none"
            )
        },
        "sp08"
      )
    })

    output$plot_gravy <- renderPlot({
      rh(
        function() {
          d <- data() |> filter(!is.na(gravy))
          if (nrow(d) == 0) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "GRAVY not available",
                  color = "gray40",
                  size = 5
                ) +
                theme_void()
            )
          }
          med_g <- round(median(d$gravy, na.rm = TRUE), 3)
          pct_pos <- round(mean(d$gravy > 0, na.rm = TRUE) * 100, 1)
          pct_neg <- round(mean(d$gravy < 0, na.rm = TRUE) * 100, 1)
          pct_amp <- round(mean(abs(d$gravy) <= 0.3, na.rm = TRUE) * 100, 1)
          pop_label <- dplyr::case_when(
            pct_pos > 60 ~ paste0(
              "predominantly HYDROPHOBIC (",
              pct_pos,
              "% GRAVY > 0)"
            ),
            pct_neg > 60 ~ paste0(
              "predominantly HYDROPHILIC (",
              pct_neg,
              "% GRAVY < 0)"
            ),
            pct_amp > 40 ~ paste0(
              "AMPHIPATHIC tendency (",
              pct_amp,
              "% within −0.3/+0.3)"
            ),
            TRUE ~ paste0(
              "mixed population (hydrophobic: ",
              pct_pos,
              "%, hydrophilic: ",
              pct_neg,
              "%)"
            )
          )
          ggplot(d, aes(x = gravy, fill = after_stat(x))) +
            geom_histogram(bins = 50, color = "black") +
            scale_fill_viridis_c(option = "C", name = "GRAVY") +
            geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
            geom_vline(
              xintercept = med_g,
              color = "orange",
              linetype = "dotted",
              linewidth = 0.9
            ) +
            annotate(
              "text",
              x = med_g,
              y = Inf,
              vjust = 1.5,
              hjust = ifelse(med_g >= 0, -0.1, 1.1),
              label = paste0("median=", med_g),
              color = "orange",
              size = 3.5,
              fontface = "bold"
            ) +
            labs(
              title = "GRAVY distribution",
              x = "GRAVY index",
              y = "Count",
              caption = paste0(
                "Positive = hydrophobic | Negative = hydrophilic | ~0 = amphipathic\n",
                "Population: ",
                pop_label,
                "."
              )
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn() +
            theme(
              legend.position = "bottom",
              legend.key.width = unit(2, "cm"),
              legend.key.height = unit(0.25, "cm"),
              plot.caption = element_text(
                size = 10,
                color = "gray30",
                hjust = 0
              )
            )
        },
        "sp09"
      )
    })

    output$plot_pi <- renderPlot({
      rh(
        function() {
          d <- data() |> filter(!is.na(pI))
          if (nrow(d) == 0) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "pI not available",
                  color = "gray40",
                  size = 5
                ) +
                theme_void()
            )
          }
          med_pi <- round(median(d$pI, na.rm = TRUE), 2)
          pct_acid <- round(mean(d$pI < 7, na.rm = TRUE) * 100, 1)
          pct_bas <- round(mean(d$pI > 7, na.rm = TRUE) * 100, 1)
          ggplot(d, aes(x = pI, fill = after_stat(x))) +
            geom_histogram(bins = 40, color = "black") +
            scale_fill_gradient2(
              low = "#e74c3c",
              mid = "#f0e6ff",
              high = "#2980b9",
              midpoint = 7,
              name = "pI"
            ) +
            geom_vline(
              xintercept = 7,
              color = "gray40",
              linetype = "dashed",
              linewidth = 0.7
            ) +
            geom_vline(
              xintercept = med_pi,
              color = "orange",
              linetype = "dotted",
              linewidth = 0.9
            ) +
            annotate(
              "text",
              x = med_pi,
              y = Inf,
              vjust = 1.5,
              hjust = ifelse(med_pi >= 7, -0.1, 1.1),
              label = paste0("median pI=", med_pi),
              color = "orange",
              size = 3.5,
              fontface = "bold"
            ) +
            labs(
              title = "Isoelectric point (pI) distribution",
              x = "Isoelectric point (pI)",
              y = "Count",
              caption = paste0(
                "Acidic (pI < 7): ",
                pct_acid,
                "%  |  Basic (pI > 7): ",
                pct_bas,
                "%\n",
                "Red = acidic end; Blue = basic end; dashed grey line = pI 7 (neutral)"
              )
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn() +
            theme(
              legend.position = "bottom",
              legend.key.width = unit(2, "cm"),
              legend.key.height = unit(0.25, "cm"),
              plot.caption = element_text(
                size = 10,
                color = "gray30",
                hjust = 0
              )
            )
        },
        "sp15"
      )
    })

    output$plot_score_length <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          d |>
            filter(!is.na(peptide_length)) |>
            group_by(filename, peptide_length) |>
            summarise(
              med = median(score, na.rm = TRUE),
              n = n(),
              .groups = "drop"
            ) |>
            ggplot(aes(x = peptide_length, y = med, size = n)) +
            geom_point(color = input$plot_color, alpha = 0.8) +
            geom_smooth(
              method = "loess",
              se = TRUE,
              color = "red",
              linetype = "dashed",
              linewidth = 0.8,
              show.legend = FALSE
            ) +
            scale_size_continuous(name = "# PSMs") +
            labs(
              title = "Median score by peptide length",
              x = "Peptide length (AA)",
              y = "Median Casanovo score"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn() +
            theme(legend.position = "bottom")
        },
        "sp10"
      )
    })

    output$plot_mz_score <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          ggplot(d, aes(x = score, y = mz_error)) +
            ggpointdensity::geom_pointdensity(size = 0.5, alpha = 0.5) +
            viridis::scale_color_viridis(option = "plasma", name = "Density") +
            geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
            labs(
              title = "m/z error vs Casanovo score",
              x = "Casanovo score",
              y = "m/z error (exp − calc)"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn() +
            theme(
              legend.position = "bottom",
              legend.key.width = unit(1.5, "cm"),
              legend.key.height = unit(0.25, "cm")
            )
        },
        "sp11"
      )
    })

    output$plot_nterm <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          w <- max(2, min(15, input$nterm_width))
          seq_lst <- split(d$stripped_sequence, d$filename)
          seq_lst <- lapply(seq_lst, function(s) prep_seqlogo(s, w))
          seq_lst <- seq_lst[sapply(seq_lst, length) >= 5]

          if (length(seq_lst) == 0) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "Not enough sequences for SeqLogo",
                  size = 5,
                  color = "gray50"
                ) +
                theme_dn()
            )
          }

          ggseqlogo::ggseqlogo(
            lapply(seq_lst, function(x) substr(x, 1, w)),
            method = "bits",
            seq_type = "AA"
          ) +
            labs(
              title = paste0("N-terminus SeqLogo (first ", w, " AA)"),
              x = "Position from N-terminus",
              y = "Bits"
            ) +
            theme_bw() +
            theme(
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              axis.text = element_text(
                size = 10,
                face = "bold",
                color = "black"
              ),
              axis.title = element_text(size = 11, face = "bold")
            )
        },
        "sp12"
      )
    })

    output$plot_cterm <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          w <- max(2, min(15, input$cterm_width))
          seq_lst <- split(d$stripped_sequence, d$filename)
          seq_lst <- lapply(seq_lst, function(seqs) {
            seqs <- seqs[!is.na(seqs) & nchar(seqs) >= w]
            cseqs <- substr(seqs, nchar(seqs) - w + 1, nchar(seqs))
            cseqs[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", cseqs)]
          })
          seq_lst <- seq_lst[sapply(seq_lst, length) >= 5]

          if (length(seq_lst) == 0) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "Not enough sequences for SeqLogo",
                  size = 5,
                  color = "gray50"
                ) +
                theme_void()
            )
          }

          ggseqlogo::ggseqlogo(seq_lst, method = "bits", seq_type = "AA") +
            labs(
              title = paste0("C-terminus SeqLogo (last ", w, " AA)"),
              x = paste0("Position (1 = AA at position -", w, " from C-term)"),
              y = "Bits"
            ) +
            theme_bw() +
            theme(
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              axis.text = element_text(
                size = 10,
                face = "bold",
                color = "black"
              ),
              axis.title = element_text(size = 11, face = "bold")
            )
        },
        "sp13"
      )
    })

    output$plot_cooccurrence <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          aas <- c(
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "V",
            "W",
            "Y"
          )
          td <- d |>
            filter(!is.na(stripped_sequence), nchar(stripped_sequence) >= 1) |>
            mutate(
              N_aa = substr(stripped_sequence, 1, 1),
              C_aa = substr(
                stripped_sequence,
                nchar(stripped_sequence),
                nchar(stripped_sequence)
              )
            ) |>
            filter(N_aa %in% aas, C_aa %in% aas)
          req(nrow(td) > 0)

          hd <- td |>
            dplyr::count(filename, N_aa, C_aa) |>
            group_by(filename) |>
            mutate(Probability = n / sum(n)) |>
            ungroup() |>
            tidyr::complete(
              filename,
              N_aa = aas,
              C_aa = aas,
              fill = list(Probability = 0)
            ) |>
            mutate(
              N_terminus = factor(N_aa, levels = aas),
              C_terminus = factor(C_aa, levels = aas)
            )

          ggplot(hd, aes(x = C_terminus, y = N_terminus, fill = Probability)) +
            geom_tile(color = "white", linewidth = 0.15) +
            scale_fill_gradient2(
              low = "white",
              mid = "#B8E6D9",
              high = "#1BB99A",
              midpoint = max(hd$Probability) / 2,
              name = "Probability"
            ) +
            labs(
              title = "N:C-terminus co-occurrence",
              x = "C-terminus AA",
              y = "N-terminus AA"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_minimal(base_size = 11) +
            theme(
              plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
              axis.text = element_text(
                size = 10,
                face = "bold",
                color = "black"
              ),
              axis.title = element_text(size = 12, face = "bold"),
              legend.title = element_text(size = 10, face = "bold"),
              legend.title.position = "top",
              legend.key.height = unit(0.25, "cm"),
              legend.key.width = unit(1, "cm"),
              legend.position = "bottom",
              panel.border = element_rect(color = "black", fill = NA)
            ) +
            coord_fixed()
        },
        "sp14"
      )
    })

    output$plot_aa_freq <- renderPlot({
      rh(
        function() {
          d <- data()
          req(nrow(d) > 0)
          base_aas <- c(
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "V",
            "W",
            "Y"
          )

          obs_df <- lapply(unique(d$filename), function(fname) {
            sub_d <- d[d$filename == fname, ]
            pep_aa <- paste0(
              sub_d$stripped_sequence[!is.na(sub_d$stripped_sequence)],
              collapse = ""
            )
            if (nchar(pep_aa) == 0) {
              return(NULL)
            }

            tc <- table(strsplit(pep_aa, "")[[1]])
            tc <- tc[names(tc) %in% base_aas]

            freq <- as.data.frame(
              tc / sum(tc, na.rm = TRUE) * 100,
              stringsAsFactors = FALSE
            )
            colnames(freq) <- c("AA", "Frequency")

            missing_aas <- setdiff(base_aas, freq$AA)
            if (length(missing_aas) > 0) {
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
          if (is.null(obs_freq) || nrow(obs_freq) == 0) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "Not enough sequences",
                  color = "gray50",
                  size = 5
                ) +
                theme_void()
            )
          }

          obs_freq$AA <- factor(obs_freq$AA, levels = base_aas)

          ggplot(obs_freq, aes(x = AA, y = Frequency)) +
            geom_col(fill = input$plot_color, color = "black", alpha = 0.85) +
            labs(
              title = "Amino acid frequencies (identified peptides)",
              x = "Amino acid",
              y = "Frequency (%)"
            ) +
            facet_wrap(~filename, ncol = 3) +
            theme_dn()
        },
        "sp16"
      )
    })

    # ── PSM Table ────────────────────────────────────────────────────────────────
    output$psm_table <- DT::renderDataTable({
      d <- data()
      req(nrow(d) > 0)

      d_show <- d |>
        # Ensure list columns do not break DataTables by converting them to string summaries if any exist
        dplyr::mutate(across(
          where(is.list),
          ~ sapply(., function(x) paste(head(x, 5), collapse = ","))
        )) |>
        dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
      DT::datatable(
        d_show,
        options = list(
          scrollX = TRUE,
          pageLength = 25,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ),
        rownames = FALSE,
        filter = "top",
        extensions = "Buttons",
        class = "cell-border stripe hover"
      )
    })

    # ── Download handler ─────────────────────────────────────────────────────────
    output$download_plots <- downloadHandler(
      filename = function() paste0("deNovo_plots_", Sys.Date(), ".zip"),
      content = function(file) {
        td <- tempdir()
        d <- data()
        fps <- character()
        save_p <- function(nm, fn, w = 10, h = 6) {
          fp <- file.path(td, nm)
          tryCatch(
            {
              ggsave(fp, fn(), width = w, height = h, dpi = 300, bg = "white")
              fps <<- c(fps, fp)
            },
            error = function(e) message("Skip: ", nm)
          )
        }
        withProgress("Saving plots...", value = 0, {
          save_p("01_score.png", function() {
            ggplot(d, aes(x = score)) +
              geom_histogram(
                bins = 80,
                fill = input$plot_color,
                color = "black"
              ) +
              labs(title = "Score distribution", x = "Score", y = "Count") +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("02_aa_score.png", function() {
            ggplot(d, aes(x = aa_score_mean)) +
              geom_histogram(
                bins = 60,
                fill = input$plot_color,
                color = "black"
              ) +
              labs(
                title = "Mean AA-score distribution",
                x = "Mean AA score",
                y = "Count"
              ) +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("03_length.png", function() {
            ggplot(d, aes(x = peptide_length)) +
              geom_histogram(
                bins = 30,
                fill = input$plot_color,
                color = "black"
              ) +
              labs(
                title = "Peptide length distribution",
                x = "Length (AA)",
                y = "Count"
              ) +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("04_charge.png", function() {
            d |>
              dplyr::count(filename, charge) |>
              mutate(charge = factor(charge)) |>
              ggplot(aes(x = charge, y = n, fill = charge)) +
              geom_col(color = "black", show.legend = FALSE) +
              scale_fill_viridis_d(option = "D") +
              labs(
                title = "Charge state distribution",
                x = "Charge",
                y = "Count"
              ) +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("05_mz_error.png", function() {
            ggplot(d, aes(x = mz_error)) +
              geom_histogram(
                bins = 80,
                fill = input$plot_color,
                color = "black"
              ) +
              geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
              labs(
                title = "m/z error distribution",
                x = "EXP − CALC m/z",
                y = "Count"
              ) +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("06_score_vs_aa.png", function() {
            ggplot(d, aes(x = score, y = aa_score_mean)) +
              ggpointdensity::geom_pointdensity(size = 0.4, alpha = 0.6) +
              viridis::scale_color_viridis(option = "plasma") +
              labs(
                title = "Score vs mean AA score",
                x = "Score",
                y = "Mean AA score"
              ) +
              facet_wrap(~filename, ncol = 3) +
              theme_dn() +
              theme(legend.position = "bottom")
          })
          incProgress(1 / 18)
          save_p("07_rt_mz_error.png", function() {
            if (
              !"retention_time" %in% names(d) || all(is.na(d$retention_time))
            ) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "No retention time"
                  ) +
                  theme_void()
              )
            }
            ggplot(d, aes(x = retention_time, y = mz_error)) +
              ggpointdensity::geom_pointdensity() +
              viridis::scale_color_viridis(option = "plasma") +
              geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
              labs(
                title = "Retention time vs. m/z error",
                x = "RT",
                y = "m/z error"
              ) +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("08_rt_density.png", function() {
            if (
              !"retention_time" %in% names(d) || all(is.na(d$retention_time))
            ) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "No retention time"
                  ) +
                  theme_void()
              )
            }
            ggplot(d, aes(x = retention_time)) +
              geom_histogram(
                bins = 80,
                fill = input$plot_color,
                color = "black"
              ) +
              labs(title = "PSM count along RT", x = "RT", y = "Count") +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("09_mods.png", function() {
            mod_df <- d |>
              filter(!is.na(mod_name)) |>
              separate_rows(mod_name, sep = "; ") |>
              dplyr::count(filename, mod_name)
            if (!nrow(mod_df)) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "No modifications"
                  ) +
                  theme_void()
              )
            }
            ggplot(mod_df, aes(x = n, y = reorder(mod_name, n))) +
              geom_col(fill = input$plot_color, color = "black") +
              geom_text(aes(label = n), hjust = -0.1) +
              labs(title = "Modification summary", x = "Count", y = NULL) +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("10_mod_pie.png", function() {
            d |>
              mutate(status = ifelse(is_modified, "Modified", "Unmodified")) |>
              dplyr::count(filename, status) |>
              ggplot(aes(x = "", y = n, fill = status)) +
              geom_col(color = "black") +
              coord_polar("y") +
              scale_fill_manual(
                values = c("Modified" = "#e67e22", "Unmodified" = "#2980b9")
              ) +
              labs(title = "Modified vs Unmodified") +
              facet_wrap(~filename, ncol = 3) +
              theme_void() +
              theme(plot.title = element_text(hjust = 0.5, face = "bold"))
          })
          incProgress(1 / 18)
          save_p("11_gravy.png", function() {
            td <- d |> filter(!is.na(gravy))
            if (nrow(td) == 0) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "GRAVY not available"
                  ) +
                  theme_dn()
              )
            }
            ggplot(td, aes(x = gravy, fill = after_stat(x))) +
              geom_histogram(bins = 50, color = "black") +
              scale_fill_viridis_c(option = "C") +
              geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
              labs(
                title = "GRAVY distribution",
                x = "GRAVY index",
                y = "Count"
              ) +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("12_pi.png", function() {
            td <- d |> filter(!is.na(pI))
            if (nrow(td) == 0) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "pI not available"
                  ) +
                  theme_dn()
              )
            }
            ggplot(td, aes(x = pI, fill = after_stat(x))) +
              geom_histogram(bins = 40, color = "black") +
              scale_fill_gradient2(
                low = "#e74c3c",
                mid = "#f0e6ff",
                high = "#2980b9",
                midpoint = 7
              ) +
              geom_vline(
                xintercept = 7,
                color = "gray40",
                linetype = "dashed"
              ) +
              labs(title = "Isoelectric point (pI)", x = "pI", y = "Count") +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("13_score_length.png", function() {
            td <- d |> filter(!is.na(peptide_length))
            if (nrow(td) == 0) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "Length not available"
                  ) +
                  theme_dn()
              )
            }
            td |>
              group_by(filename, peptide_length) |>
              summarise(
                med = median(score, na.rm = TRUE),
                n = n(),
                .groups = "drop"
              ) |>
              ggplot(aes(x = peptide_length, y = med, size = n)) +
              geom_point(color = input$plot_color) +
              geom_smooth(
                method = "loess",
                se = FALSE,
                color = "red",
                linetype = "dashed"
              ) +
              labs(
                title = "Median score by length",
                x = "Length",
                y = "Median score"
              ) +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("14_mz_score.png", function() {
            ggplot(d, aes(x = score, y = mz_error)) +
              ggpointdensity::geom_pointdensity() +
              viridis::scale_color_viridis(option = "plasma") +
              geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
              labs(title = "m/z error vs score", x = "Score", y = "m/z error") +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          save_p("15_aa_freq.png", function() {
            base_aas <- c(
              "A",
              "C",
              "D",
              "E",
              "F",
              "G",
              "H",
              "I",
              "K",
              "L",
              "M",
              "N",
              "P",
              "Q",
              "R",
              "S",
              "T",
              "V",
              "W",
              "Y"
            )

            obs_df <- lapply(unique(d$filename), function(fname) {
              sub_d <- d[d$filename == fname, ]
              pep_aa <- paste0(
                sub_d$stripped_sequence[!is.na(sub_d$stripped_sequence)],
                collapse = ""
              )
              if (nchar(pep_aa) == 0) {
                return(NULL)
              }

              tc <- table(strsplit(pep_aa, "")[[1]])
              tc <- tc[names(tc) %in% base_aas]

              freq <- as.data.frame(
                tc / sum(tc, na.rm = TRUE) * 100,
                stringsAsFactors = FALSE
              )
              colnames(freq) <- c("AA", "Frequency")

              missing_aas <- setdiff(base_aas, freq$AA)
              if (length(missing_aas) > 0) {
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
            if (is.null(obs_freq) || nrow(obs_freq) == 0) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "Not enough sequences",
                    color = "gray50",
                    size = 5
                  ) +
                  theme_void()
              )
            }

            obs_freq$AA <- factor(obs_freq$AA, levels = base_aas)

            ggplot(obs_freq, aes(x = AA, y = Frequency)) +
              geom_col(fill = input$plot_color, color = "black") +
              labs(title = "AA Frequencies", x = "AA", y = "Freq (%)") +
              facet_wrap(~filename, ncol = 3) +
              theme_dn()
          })
          incProgress(1 / 18)
          wn <- max(2, min(15, input$nterm_width))
          wc <- max(2, min(15, input$cterm_width))
          save_p("16_nterm.png", function() {
            df <- d |>
              filter(!is.na(stripped_sequence), nchar(stripped_sequence) >= wn)
            if (nrow(df) == 0) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "Not enough sequences"
                  ) +
                  theme_void()
              )
            }
            seqs <- sapply(
              df$stripped_sequence,
              function(x) {
                if (nchar(x) >= wn) {
                  return(substr(x, 1, wn))
                }
                paste0(substr(x, 1, nchar(x)), strrep("-", wn - nchar(x)))
              },
              USE.NAMES = FALSE
            )
            seq_lst <- split(seqs, df$filename)
            seq_lst <- seq_lst[sapply(seq_lst, length) >= 5]
            if (length(seq_lst) == 0) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "Not enough sequences"
                  ) +
                  theme_void()
              )
            }
            ggseqlogo::ggseqlogo(seq_lst, method = "bits", seq_type = "AA") +
              labs(
                title = paste0("N-term SeqLogo (", wn, " AA)"),
                x = "Pos",
                y = "Bits"
              ) +
              theme_bw()
          })
          incProgress(1 / 18)
          save_p("17_cterm.png", function() {
            df <- d |>
              filter(!is.na(stripped_sequence), nchar(stripped_sequence) >= wc)
            if (nrow(df) == 0) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "Not enough sequences"
                  ) +
                  theme_void()
              )
            }
            cseqs <- sapply(
              df$stripped_sequence,
              function(x) {
                substr(x, nchar(x) - wc + 1, nchar(x))
              },
              USE.NAMES = FALSE
            )
            valid <- grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", cseqs)
            seq_lst <- split(cseqs[valid], df$filename[valid])
            seq_lst <- seq_lst[sapply(seq_lst, length) >= 5]
            if (length(seq_lst) == 0) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "Not enough sequences"
                  ) +
                  theme_void()
              )
            }
            ggseqlogo::ggseqlogo(seq_lst, method = "bits", seq_type = "AA") +
              labs(
                title = paste0("C-term SeqLogo (", wc, " AA)"),
                x = "Pos",
                y = "Bits"
              ) +
              theme_bw()
          })
          incProgress(1 / 18)
          save_p("18_cooccurrence.png", function() {
            aas <- c(
              "A",
              "C",
              "D",
              "E",
              "F",
              "G",
              "H",
              "I",
              "K",
              "L",
              "M",
              "N",
              "P",
              "Q",
              "R",
              "S",
              "T",
              "V",
              "W",
              "Y"
            )
            td <- d |>
              filter(
                !is.na(stripped_sequence),
                nchar(stripped_sequence) >= 1
              ) |>
              mutate(
                N_aa = substr(stripped_sequence, 1, 1),
                C_aa = substr(
                  stripped_sequence,
                  nchar(stripped_sequence),
                  nchar(stripped_sequence)
                )
              ) |>
              filter(N_aa %in% aas, C_aa %in% aas)
            if (!nrow(td)) {
              return(
                ggplot() +
                  annotate(
                    "text",
                    x = 0.5,
                    y = 0.5,
                    label = "Not enough data"
                  ) +
                  theme_void()
              )
            }
            hd <- td |>
              dplyr::count(filename, N_aa, C_aa) |>
              group_by(filename) |>
              mutate(Probability = n / sum(n)) |>
              ungroup() |>
              tidyr::complete(
                filename,
                N_aa = aas,
                C_aa = aas,
                fill = list(Probability = 0)
              ) |>
              mutate(
                N_terminus = factor(N_aa, levels = aas),
                C_terminus = factor(C_aa, levels = aas)
              )

            ggplot(
              hd,
              aes(x = C_terminus, y = N_terminus, fill = Probability)
            ) +
              geom_tile(color = "white") +
              scale_fill_gradient2(
                low = "white",
                mid = "#B8E6D9",
                high = "#1BB99A",
                midpoint = max(hd$Probability) / 2
              ) +
              labs(title = "N:C co-occurrence", x = "C-term", y = "N-term") +
              facet_wrap(~filename, ncol = 3) +
              theme_minimal() +
              coord_fixed()
          })
          incProgress(1 / 18)
        })
        utils::zip(file, files = fps, flags = "-j")
      },
      contentType = "application/zip"
    )
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# STANDALONE APP  (skipped when sourced inside proteOmni with proteOmni_LOADED=TRUE)
# ═══════════════════════════════════════════════════════════════════════════════
if (!exists("proteOmni_LOADED")) {
  dark_css <- "
    body, .content-wrapper, .main-footer { background-color:#0a111c !important; color:#e0e0e0; }
    .skin-black .main-header .logo,
    .skin-black .main-header .navbar   { background-color:#0d1b2a !important; border-bottom:1px solid #1c2e42; }
    .skin-black .main-sidebar           { background-color:#111e2e !important; }
    .skin-black .sidebar-menu > li.active > a,
    .skin-black .sidebar-menu > li:hover > a { background:#182538 !important; }
    .box { background-color:#111e2e !important; border:1px solid #1c2e42 !important; color:#e0e0e0; }
    .box.box-primary { border-top-color:#2980b9 !important; }
    .box-header.with-border { border-bottom:1px solid #1c2e42; background:#182538 !important; }
    .box-title { color:#7fb3d5 !important; font-weight:700; }
    .nav-tabs-custom { background:#111e2e; border-color:#1c2e42; }
    .nav-tabs-custom>.nav-tabs>li.active>a { border-top:2px solid #2980b9; color:#7fb3d5; background:#182538; }
    .nav-tabs-custom>.nav-tabs>li>a { color:#adb5bd; background:#111e2e; }
    .plot-wrap { position:relative; min-height:250px; }
    .spinner-overlay { display:none; position:absolute; top:50%; left:50%;
                       transform:translate(-50%,-50%); color:#2980b9; font-size:2em; z-index:10; }
    .info-box, .info-box-content { background:#111e2e !important; color:#e0e0e0 !important; }
    .info-box .info-box-icon { background:#1f618d !important; }
    .form-control, input[type='number'] { background:#1a2332 !important; color:#e0e0e0 !important;
                                          border:1px solid #2d3741 !important; }
    .control-label, label { color:#adb5bd !important; font-weight:600; }
    .dl-btn { background-color:#1f618d !important; color:#fff !important;
              border:none !important; border-radius:4px !important; }
    .dl-btn:hover { background-color:#2980b9 !important; }
    .irs--shiny .irs-bar, .irs--shiny .irs-bar-edge { background:#2980b9 !important; border-color:#2980b9 !important; }
    .irs--shiny .irs-handle > i:first-child { background:#2980b9 !important; border-color:#2980b9 !important; }
    .irs--shiny .irs-from, .irs--shiny .irs-to,
    .irs--shiny .irs-single { background:#2980b9 !important; }
    table.dataTable { background:#111e2e !important; color:#e0e0e0 !important; }
    table.dataTable thead th,tr.odd,tr.even { background:#182538 !important; color:#7fb3d5 !important; }
    table.dataTable tbody tr.odd  { background:#111e2e !important; }
    table.dataTable tbody tr.even { background:#0f1b2b !important; }
    table.dataTable tbody tr:hover { background:#1c2e42 !important; }
    .dataTables_wrapper .dataTables_filter input,
    .dataTables_wrapper .dataTables_length select { background:#1a2332 !important; color:#e0e0e0 !important;
                                                     border:1px solid #2d3741 !important; }
    .dataTables_wrapper .dataTables_info,
    .dataTables_wrapper .dataTables_paginate .paginate_button { color:#adb5bd !important; }
    details summary::-webkit-details-marker { color:#7fb3d5; }
  "

  ui <- tagList(
    tags$head(tags$title("de Novo Viewer — Casanovo")),
    useShinyjs(),
    shinydashboard::dashboardPage(
      skin = "black",
      shinydashboard::dashboardHeader(
        title = tags$span(
          style = "color:#ffffff;",
          icon("fingerprint"),
          " de Novo Viewer — Casanovo"
        ),
        titleWidth = 310
      ),
      shinydashboard::dashboardSidebar(
        width = 290,
        tags$head(tags$style(HTML(dark_css))),
        deNovo_sidebar_ui("dnv")
      ),
      shinydashboard::dashboardBody(
        useShinyjs(),
        tags$head(tags$style(HTML(dark_css))),
        deNovo_body_ui("dnv")
      )
    )
  )

  server <- function(input, output, session) {
    deNovo_server("dnv")
  }

  shinyApp(ui = ui, server = server)
}
