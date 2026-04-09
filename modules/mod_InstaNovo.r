## ============================================================
## mod_InstaNovo.r  —  InstaNovo de novo results viewer
## ============================================================

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
  library(ggpointdensity)
  library(viridis)
  library(colourpicker)
  library(DT)
})

# ── Helper functions ──────────────────────────────────────────────────────────

strip_sequence_instanovo <- function(seq) {
  s <- as.character(seq)
  s <- gsub("[^A-Z]", "", s)
  s
}

read_instanovo_csv <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)

  req_cols <- c(
    "predictions",
    "log_probs",
    "precursor_charge",
    "delta_mass_ppm"
  )
  missing <- setdiff(req_cols, names(df))
  if (length(missing) > 0) {
    stop(paste(
      "Missing columns in InstaNovo CSV:",
      paste(missing, collapse = ", ")
    ))
  }

  df <- df |>
    mutate(
      score = as.numeric(log_probs),
      charge = as.numeric(precursor_charge),
      mz_error_ppm = as.numeric(delta_mass_ppm),
      stripped_sequence = sapply(
        predictions,
        strip_sequence_instanovo,
        USE.NAMES = FALSE
      ),
      peptide_length = nchar(stripped_sequence)
    )

  df <- df |> filter(!is.na(stripped_sequence), peptide_length > 0)

  list(psm = df)
}

prep_seqlogo_instanovo <- function(seqs, width) {
  seqs <- seqs[!is.na(seqs) & nchar(seqs) >= width]
  seqs <- substr(seqs, 1, width)
  seqs[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", seqs)]
}

theme_ins <- function(...) {
  theme_bw(...) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10, face = "bold", color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title.position = "top",
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold"),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

sp_wrap_ins <- function(sid, ui_el) {
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
InstaNovo_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      style = "padding:12px 16px 4px;color:#ffffff;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      icon("dna", lib = "font-awesome"),
      " InstaNovo de Novo"
    ),
    tags$div(
      style = "padding: 0 16px; margin-bottom: 12px;",
      textInput(
        ns("csv_dir"),
        "Path to results directory:",
        placeholder = "/path/to/csv_files/"
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
      "Min log_probs score",
      min = -100,
      max = 0,
      value = -100,
      step = 0.5
    ),
    colourpicker::colourInput(
      ns("plot_color"),
      "Plot colour",
      value = "#27ae60"
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
InstaNovo_body_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      id = ns("tabs"),
      type = "tabs",

      # ── Overview ─────────────────────────────────────────────────────────────
      tabPanel(
        "Overview",
        fluidRow(infoBoxOutput(ns("info_box"), width = 12)),
        fluidRow(
          box(
            title = "Score distribution (log_probs)",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi01"), uiOutput(ns("plot_score_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Peptide length distribution",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi02"), uiOutput(ns("plot_length_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Charge state distribution",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi03"), uiOutput(ns("plot_charge_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Mass error distribution (ppm)",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi04"), uiOutput(ns("plot_mz_error_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Number of retained PSMs vs. score threshold",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi05"), uiOutput(ns("plot_elbow_ui")))
          )
        )
      ),

      # ── Peptide Analysis ─────────────────────────────────────────────────────
      tabPanel(
        "Peptide Analysis",
        fluidRow(
          box(
            title = "Median score by peptide length",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi06"), uiOutput(ns("plot_score_length_ui")))
          )
        ),
        fluidRow(
          box(
            title = "ppm error vs InstaNovo score",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi07"), uiOutput(ns("plot_mz_score_ui")))
          )
        ),
        fluidRow(
          box(
            title = "GRAVY Index Distribution",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi11"), uiOutput(ns("plot_gravy_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Isoelectric Point (pI) Distribution",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi12"), uiOutput(ns("plot_pi_ui")))
          )
        ),
        fluidRow(
          box(
            title = "Amino acid frequencies (Identified peptides)",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(ns("spi08"), uiOutput(ns("plot_aa_freq_ui")))
          )
        )
      ),

      # ── Termini SeqLogos ──────────────────────────────────────────────────────
      tabPanel(
        "Termini SeqLogos",
        fluidRow(
          box(
            title = "N-terminus SeqLogo",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(
              ns("spi09"),
              plotOutput(ns("plot_nterm"), height = "280px")
            )
          ),
          box(
            title = "C-terminus SeqLogo",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            sp_wrap_ins(
              ns("spi10"),
              plotOutput(ns("plot_cterm"), height = "280px")
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

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — SERVER
# ═══════════════════════════════════════════════════════════════════════════════
InstaNovo_server <- function(id, fasta_digest) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    spin_ids <- paste0("spi", sprintf("%02d", 1:10))
    show_all <- function() lapply(spin_ids, function(s) shinyjs::show(id = s))
    hide_sp <- function(s) shinyjs::hide(id = s)
    rh <- function(fn, sid) {
      on.exit(hide_sp(sid), add = TRUE)
      fn()
    }

    # ── Load raw data ───────────────────────────────────────────────────────────
    raw_data_rv <- reactiveValues(psm = NULL)

    observeEvent(input$load_dir_btn, {
      dir_path <- trimws(input$csv_dir)
      if (dir_path == "" || !dir.exists(dir_path)) {
        showNotification("Please enter a valid directory path.", type = "error")
        return()
      }
      show_all()
      withProgress(message = "Parsing InstaNovo results...", value = 0, {
        csv_files <- list.files(
          dir_path,
          pattern = "\\.csv$",
          full.names = TRUE,
          recursive = TRUE
        )
        if (length(csv_files) == 0) {
          showNotification(
            "No .csv files found in the directory.",
            type = "warning"
          )
          setProgress(1)
          hide_sp("spi01")
          return()
        }

        psm_list <- list()
        n <- length(csv_files)

        for (i in seq_along(csv_files)) {
          fpath <- csv_files[i]
          fname <- str_extract(basename(fpath), "^[^\\.]+")
          setProgress(i / n, detail = paste("Reading", basename(fpath)))

          res <- tryCatch(read_instanovo_csv(fpath), error = function(e) NULL)
          if (!is.null(res) && !is.null(res$psm) && nrow(res$psm) > 0) {
            res$psm$filename <- fname
            psm_list[[fname]] <- res$psm
          }
        }

        if (length(psm_list) > 0) {
          raw_data_rv$psm <- bind_rows(psm_list)
        } else {
          raw_data_rv$psm <- NULL
        }
        setProgress(1)
      })
    })

    raw_data <- reactive({
      req(raw_data_rv$psm)
      list(psm = raw_data_rv$psm)
    })

    # ── Update slider limits ────────────────────────────────────────────────────
    observeEvent(raw_data_rv$psm, {
      psm_dat <- raw_data_rv$psm
      sc <- psm_dat$score[is.finite(psm_dat$score)]
      if (length(sc)) {
        sc_lo <- floor(min(sc) * 10) / 10
        sc_hi <- ceiling(max(sc) * 10) / 10
        if (sc_hi > 0) {
          sc_hi <- 0
        }
        updateSliderInput(
          session,
          "score_filter",
          min = sc_lo,
          max = sc_hi,
          value = sc_lo
        )
      }
    })

    # ── Filtered data
    data <- reactive({
      rd <- raw_data()
      req(rd, rd$psm)
      rd$psm |> filter(is.na(score) | score >= input$score_filter)
    })

    mapped_data <- reactive({
      d <- data()
      dig <- fasta_digest()
      if (!is.null(dig)) {
        classes <- classify_peptides(d$stripped_sequence, dig)
        d <- left_join(d, classes, by = c("stripped_sequence" = "peptide"))
      } else {
        d$classification <- "Unmapped"
        d$mapped_proteins <- NA_character_
      }

      safe_seqs <- str_remove_all(
        d$stripped_sequence,
        "[^ACDEFGHIKLMNPQRSTVWY]"
      )
      d$gravy <- sapply(safe_seqs, GRAVY)
      d$pI <- sapply(safe_seqs, calculate_pI)
      d
    })

    output$info_box <- renderInfoBox({
      d <- data()
      total <- nrow(raw_data()$psm)
      n <- nrow(d)
      pct <- if (total > 0) round(n / total * 100, 1) else 0
      infoBox(
        "InstaNovo Filter",
        paste0(
          n,
          " / ",
          total,
          " PSMs  (",
          pct,
          "%)  |  score ≥ ",
          input$score_filter
        ),
        icon = icon("filter"),
        color = "green"
      )
    })

    plot_annotated_density <- function(d, col_sym, title, xlab, color) {
      d <- d |> filter(!is.na(!!sym(col_sym)))
      if (nrow(d) == 0) {
        return(
          ggplot() +
            annotate("text", x = 0, y = 0, label = "Not enough data") +
            theme_void()
        )
      }
      med_df <- d |>
        group_by(filename) |>
        summarise(m = median(!!sym(col_sym), na.rm = TRUE), .groups = "drop")

      ggplot(d, aes(x = !!sym(col_sym))) +
        geom_density(fill = color, color = "black", alpha = 0.6) +
        geom_vline(
          data = med_df,
          aes(xintercept = m),
          linetype = "dashed",
          color = "red",
          linewidth = 1
        ) +
        geom_text(
          data = med_df,
          aes(x = m, y = Inf, label = paste("Median:", round(m, 2))),
          vjust = 2,
          hjust = -0.1,
          color = "red",
          fontface = "bold"
        ) +
        labs(title = title, x = xlab, y = "Density") +
        facet_wrap(~filename, ncol = 3) +
        theme_ins()
    }

    plot_score_obj <- reactive({
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
          title = "InstaNovo log_probs distribution",
          x = "InstaNovo Score (log_prob)",
          y = "Count"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_ins()
    })

    plot_length_obj <- reactive({
      d <- data()
      req(nrow(d) > 0)
      len_range <- range(d$peptide_length, na.rm = TRUE)
      ggplot(d, aes(x = peptide_length)) +
        geom_histogram(
          bins = max(1, len_range[2] - len_range[1] + 1),
          fill = input$plot_color,
          color = "black",
          alpha = 0.85
        ) +
        labs(
          title = "Peptide length distribution",
          x = "Length (AA)",
          y = "Count"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_ins()
    })

    plot_charge_obj <- reactive({
      d <- data()
      req(nrow(d) > 0)
      d |>
        dplyr::count(filename, charge) |>
        mutate(charge = factor(charge)) |>
        ggplot(aes(x = charge, y = n, fill = charge)) +
        geom_col(color = "black", show.legend = FALSE) +
        geom_text(aes(label = n), vjust = -0.4, fontface = "bold", size = 4.5) +
        scale_fill_viridis_d(option = "D") +
        labs(
          title = "Charge state distribution",
          x = "Charge (z)",
          y = "Count"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_ins()
    })

    plot_mz_error_obj <- reactive({
      d <- data()
      req(nrow(d) > 0)
      ggplot(d, aes(x = mz_error_ppm)) +
        geom_histogram(
          bins = 80,
          fill = input$plot_color,
          color = "black",
          alpha = 0.85
        ) +
        geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
        labs(
          title = "Mass error distribution (ppm)",
          x = "delta_mass_ppm",
          y = "Count"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_ins()
    })

    plot_elbow_obj <- reactive({
      rd <- raw_data()
      req(rd, rd$psm)
      td <- rd$psm |> dplyr::filter(!is.na(score), !is.na(filename))
      if (nrow(td) < 10) {
        return(
          ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = "Not enough scored PSMs") +
            theme_void()
        )
      }

      df_curve <- td |>
        group_by(filename) |>
        arrange(desc(score), .by_group = TRUE) |>
        mutate(PSMs_Retained = row_number()) |>
        ungroup() |>
        rename(Score_Threshold = score)

      df_curve <- df_curve |>
        group_by(filename) |>
        slice(seq(1, n(), length.out = min(n(), 500))) |>
        ungroup()

      ggplot(df_curve, aes(x = Score_Threshold, y = PSMs_Retained)) +
        geom_line(color = input$plot_color, linewidth = 1.2) +
        geom_vline(
          xintercept = input$score_filter,
          color = "red",
          linetype = "dashed",
          linewidth = 1
        ) +
        labs(
          title = "Number of retained PSMs vs. score threshold",
          x = "Minimum log_probs threshold",
          y = "Number of retained PSMs"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_ins()
    })

    plot_score_length_obj <- reactive({
      d <- data()
      req(nrow(d) > 0)
      d |>
        dplyr::filter(!is.na(peptide_length)) |>
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
          y = "Median InstaNovo score"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_ins() +
        theme(legend.position = "bottom")
    })

    plot_mz_score_obj <- reactive({
      d <- data()
      req(nrow(d) > 0)
      ggplot(d, aes(x = score, y = mz_error_ppm)) +
        ggpointdensity::geom_pointdensity(size = 0.5, alpha = 0.5) +
        viridis::scale_color_viridis(option = "plasma", name = "Density") +
        geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
        labs(
          title = "ppm error vs InstaNovo score",
          x = "InstaNovo score",
          y = "delta_mass_ppm"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_ins() +
        theme(legend.position = "bottom")
    })

    plot_aa_freq_obj <- reactive({
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
            annotate("text", x = 0.5, y = 0.5, label = "Not enough sequences") +
            theme_void()
        )
      }

      obs_freq$AA <- factor(obs_freq$AA, levels = base_aas)

      ggplot(obs_freq, aes(x = AA, y = Frequency)) +
        geom_col(fill = input$plot_color, color = "black", alpha = 0.85) +
        labs(
          title = "Amino acid frequencies",
          x = "Amino acid",
          y = "Frequency (%)"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_ins()
    })

    plot_gravy_obj <- reactive({
      d <- mapped_data()
      req(nrow(d) > 0, "gravy" %in% names(d))
      plot_annotated_density(
        d,
        "gravy",
        "GRAVY Index Distribution",
        "GRAVY Index",
        input$plot_color
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
              label = "pI not available (Peptides package missing)"
            ) +
            theme_void()
        )
      }
      plot_annotated_density(
        d,
        "pI",
        "Isoelectric Point (pI) Distribution",
        "pI",
        input$plot_color
      )
    })

    plot_nterm_obj <- reactive({
      d <- data()
      req(nrow(d) > 0)
      w <- max(2, min(15, input$nterm_width))
      df <- d |>
        filter(!is.na(stripped_sequence), nchar(stripped_sequence) >= w)
      if (nrow(df) == 0) {
        return(
          ggplot() +
            annotate(
              "text",
              x = 0.5,
              y = 0.5,
              label = "Not enough sequences",
              size = 5
            ) +
            theme_void()
        )
      }
      seqs <- sapply(
        df$stripped_sequence,
        function(x) {
          if (nchar(x) >= w) {
            return(substr(x, 1, w))
          }
          paste0(substr(x, 1, nchar(x)), strrep("-", w - nchar(x)))
        },
        USE.NAMES = FALSE
      )

      seq_lst <- split(seqs, df$filename)
      seq_lst <- seq_lst[sapply(seq_lst, length) >= 5]
      if (length(seq_lst) == 0) {
        return(
          ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = "Not enough sequences") +
            theme_void()
        )
      }
      ggseqlogo::ggseqlogo(seq_lst, method = "bits", seq_type = "AA") +
        labs(
          title = paste0("N-terminus SeqLogo (first ", w, " AA)"),
          x = "Position",
          y = "Bits"
        ) +
        theme_bw() +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    })

    plot_cterm_obj <- reactive({
      d <- data()
      req(nrow(d) > 0)
      w <- max(2, min(15, input$cterm_width))
      df <- d |>
        filter(!is.na(stripped_sequence), nchar(stripped_sequence) >= w)
      if (nrow(df) == 0) {
        return(
          ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = "Not enough sequences") +
            theme_void()
        )
      }

      cseqs <- sapply(
        df$stripped_sequence,
        function(x) {
          substr(x, nchar(x) - w + 1, nchar(x))
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
              label = "Not enough sequences",
              size = 5
            ) +
            theme_void()
        )
      }
      ggseqlogo::ggseqlogo(seq_lst, method = "bits", seq_type = "AA") +
        labs(
          title = paste0("C-terminus SeqLogo (last ", w, " AA)"),
          x = "Position",
          y = "Bits"
        ) +
        theme_bw() +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    })

    # ── Dynamic plot height ───────────────────────────────────────────────
    ins_plot_h <- reactive({
      d <- data()
      n <- dplyr::n_distinct(d$filename)
      max(400L, ceiling(n / 3L) * 350L)
    })

    make_ins_plot_ui <- function(plot_id) {
      renderUI({ req(data()); plotOutput(ns(plot_id), height = paste0(ins_plot_h(), "px")) })
    }

    output$plot_score_ui        <- make_ins_plot_ui("plot_score")
    output$plot_length_ui       <- make_ins_plot_ui("plot_length")
    output$plot_charge_ui       <- make_ins_plot_ui("plot_charge")
    output$plot_mz_error_ui     <- make_ins_plot_ui("plot_mz_error")
    output$plot_elbow_ui        <- make_ins_plot_ui("plot_elbow")
    output$plot_score_length_ui <- make_ins_plot_ui("plot_score_length")
    output$plot_mz_score_ui     <- make_ins_plot_ui("plot_mz_score")
    output$plot_aa_freq_ui      <- make_ins_plot_ui("plot_aa_freq")
    output$plot_gravy_ui        <- make_ins_plot_ui("plot_gravy")
    output$plot_pi_ui           <- make_ins_plot_ui("plot_pi")

    # ── Render UI plots
    output$plot_score <- renderPlot({
      rh(function() plot_score_obj(), "spi01")
    })
    output$plot_length <- renderPlot({
      rh(function() plot_length_obj(), "spi02")
    })
    output$plot_charge <- renderPlot({
      rh(function() plot_charge_obj(), "spi03")
    })
    output$plot_mz_error <- renderPlot({
      rh(function() plot_mz_error_obj(), "spi04")
    })
    output$plot_elbow <- renderPlot({
      rh(function() plot_elbow_obj(), "spi05")
    })
    output$plot_score_length <- renderPlot({
      rh(function() plot_score_length_obj(), "spi06")
    })
    output$plot_mz_score <- renderPlot({
      rh(function() plot_mz_score_obj(), "spi07")
    })
    output$plot_aa_freq <- renderPlot({
      rh(function() plot_aa_freq_obj(), "spi08")
    })
    output$plot_nterm <- renderPlot({
      rh(function() plot_nterm_obj(), "spi09")
    })
    output$plot_cterm <- renderPlot({
      rh(function() plot_cterm_obj(), "spi10")
    })
    output$plot_gravy <- renderPlot({
      rh(function() plot_gravy_obj(), "spi11")
    })
    output$plot_pi <- renderPlot({
      rh(function() plot_pi_obj(), "spi12")
    })

    # ── PSM Table
    output$psm_table <- DT::renderDataTable({
      d <- mapped_data()
      req(nrow(d) > 0)

      d_show <- d |>
        mutate(across(
          where(is.list),
          ~ sapply(., function(x) paste(head(x, 5), collapse = ","))
        )) |>
        mutate(across(where(is.numeric), ~ round(., 4)))

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

    # ── Download Plots
    output$download_plots <- downloadHandler(
      filename = function() {
        paste0("InstaNovo_Plots_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip")
      },
      content = function(file) {
        req(plot_score_obj())

        tmp_dir <- file.path(
          tempdir(),
          paste0("InstaNovo_Plots_", as.integer(Sys.time()))
        )
        dir.create(tmp_dir, showWarnings = FALSE)

        tryCatch(
          ggsave(
            file.path(tmp_dir, "01_Score_Dist.png"),
            plot_score_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "02_Length_Dist.png"),
            plot_length_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "03_Charge_Dist.png"),
            plot_charge_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "04_MZ_Error_Dist.png"),
            plot_mz_error_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "05_PSMs_vs_Threshold.png"),
            plot_elbow_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "06_Median_Score_Length.png"),
            plot_score_length_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "07_PPM_Error_vs_Score.png"),
            plot_mz_score_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "08_AA_Frequencies.png"),
            plot_aa_freq_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "09_N_Term_SeqLogo.png"),
            plot_nterm_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "10_C_Term_SeqLogo.png"),
            plot_cterm_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "11_GRAVY_Dist.png"),
            plot_gravy_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )
        tryCatch(
          ggsave(
            file.path(tmp_dir, "12_pI_Dist.png"),
            plot_pi_obj(),
            width = 11,
            height = 8,
            dpi = 300
          ),
          error = function(e) NULL
        )

        owd <- setwd(tmp_dir)
        on.exit(setwd(owd))
        utils::zip(zipfile = file, files = list.files(tmp_dir))
      }
    )
  })
}
