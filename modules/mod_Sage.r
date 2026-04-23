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
    fileInput(
      ns("sage_file"),
      "Choose results.sage.tsv or .parquet",
      accept = c(".tsv", ".parquet")
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
    ),
    actionButton(
      ns("run_plot"),
      "Plot Selected Graphic",
      icon = icon("chart-bar"),
      class = "btn-primary",
      style = "width:80%;margin-bottom:8px;"
    ),
    div(
      style = "padding:0 8px;",
      downloadButton(
        ns("download_plot"),
        "⬇ Download Plot (.png)",
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
                  choices = NULL
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

    # ── Load data
    sage_data <- reactive({
      req(input$sage_file)
      show_all()
      withProgress(message = "Parsing Sage results...", value = 0.5, {
        ext <- tools::file_ext(input$sage_file$name)
        res <- tryCatch(
          {
            df <- NULL
            if (ext == "parquet") {
              df <- arrow::read_parquet(input$sage_file$datapath)
            } else {
              df <- data.table::fread(
                input$sage_file$datapath,
                sep = "\t",
                header = TRUE
              )
            }

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
            showNotification(
              paste("Error reading Sage file:", e$message),
              type = "error"
            )
            NULL
          }
        )
        setProgress(1)
        res
      })
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
        d <- d |> filter(is.na(peptide_q) | peptide_q <= input$qval_filter)
      }

      if ("sage_discriminant_score" %in% names(d)) {
        d <- d |>
          filter(
            is.na(sage_discriminant_score) |
              sage_discriminant_score >= input$lda_filter
          )
      }

      if (input$filter_decoy && "is_decoy" %in% names(d)) {
        d <- d |> filter(is_decoy == FALSE)
      }

      if ("stripped_peptide" %in% names(d)) {
        safe_seqs <- str_remove_all(
          d$stripped_peptide,
          "[^ACDEFGHIKLMNPQRSTVWY]"
        )
        d$gravy <- sapply(safe_seqs, GRAVY)
        d$pI <- sapply(safe_seqs, calculate_pI)
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
        d <- left_join(d, classes, by = c("stripped_peptide" = "peptide"))
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
        group_by(filename) |>
        summarise(is_decoy = is_decoy, .groups = "drop") |>
        ggplot(aes(x = filename, fill = as.character(is_decoy))) +
        geom_bar(position = "dodge") +
        labs(x = "File", y = "Number of PSMs", fill = "Is decoy?") +
        scale_fill_manual(values = vals_decoy()) +
        theme_sage() +
        theme(legend.position = "top")
    })

    plot_id_counts_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      summ <- d |>
        group_by(filename) |>
        summarise(
          n_peptides = n_distinct(stripped_peptide),
          n_proteins = if ("proteins" %in% names(d)) {
            n_distinct(proteins)
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
        filter(!is.na(!!sym(col_sym))) |>
        group_by(filename) |>
        summarise(m = median(!!sym(col_sym)), .groups = "drop")

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
      plot_annotated_density_faceted(
        d,
        "pI",
        "Isoelectric Point (pI) Distribution",
        "pI",
        input$color_target
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

    mod_diag_plot_obj <- reactive({
      paired <- mod_diag_data()
      req(paired, nrow(paired) > 0)

      top_peptides <- paired |>
        dplyr::group_by(sample_name, peptide_seq) |>
        dplyr::summarise(
          max_delta = max(delta_rt, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::group_by(sample_name) |>
        dplyr::slice_max(order_by = max_delta, n = 30) |>
        dplyr::pull(peptide_seq) |>
        unique()

      plot_data <- paired |> dplyr::filter(peptide_seq %in% top_peptides)
      req(nrow(plot_data) > 0)

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
      hide_sp("sp_moddiag")
      if (is.null(paired) || nrow(paired) == 0) {
        return(tags$p(
          style = "color:#adb5bd;text-align:center;padding:20px;",
          "No paired modified/unmodified peptides found. Load a Sage results file and ensure it contains 'peptide', 'stripped_peptide', 'rt', and 'filename' columns."
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
          peptide_seq,
          mod_seq,
          rt_unmod,
          rt_mod,
          delta_rt,
          classification
        ) |>
        dplyr::rename(peptide = peptide_seq, modified_peptide = mod_seq) |>
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
