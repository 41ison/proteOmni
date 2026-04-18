## ============================================================
## mod_EncyclopeDIA.r  —  EncyclopeDIA DIA results viewer
## ============================================================

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
library(DT)

# ── Helper functions ──────────────────────────────────────────────────────────

theme_enc <- function(...) {
  theme_bw(...) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold"),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

sp_wrap_enc <- function(sid, ui_el) {
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
      actionButton(
        ns("load_dir_btn"),
        "Load Directory",
        icon = icon("folder-open"),
        class = "btn-primary",
        style = "width:80%;"
      )
    ),
    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),
    selectInput(
      ns("plot_select"),
      "Select Graphic",
      choices = c(
        "Protein and Peptide Identifications" = "plot_identifications",
        "Score Distribution" = "plot_score",
        "Posterior Error Probability (PEP)" = "plot_pep",
        "q-value Distribution" = "plot_qval",
        "Peptide Yield vs. FDR" = "plot_fdr_curve",
        "Charge State" = "plot_charge",
        "Modifications" = "plot_mods",
        "Peptide Length" = "plot_length",
        "GRAVY Index Distribution" = "plot_gravy",
        "pI Distribution" = "plot_pi",
        "Amino acid frequencies" = "plot_aa_freq"
      )
    ),
    actionButton(
      ns("run_plot"),
      "Plot Selected Graphic",
      icon = icon("chart-bar"),
      class = "btn-primary",
      style = "width:80%;margin-bottom:8px;"
    ),
    sliderInput(
      ns("qval_filter"),
      "Max q-value",
      min = 0,
      max = 0.05,
      value = 0.01,
      step = 0.005
    ),
    colourpicker::colourInput(
      ns("bar_color"),
      "Proteins Bar colour",
      value = "#1b9e77"
    ),
    colourpicker::colourInput(
      ns("line_color"),
      "Peptides Line colour",
      value = "#d95f02"
    ),
    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),
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
EncyclopeDIA_body_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      id = ns("enc_tabs"),
      # ── Interactive Plot Viewer ──────────────────────────────────────────────
      tabPanel(
        "Interactive Plot Viewer",
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
        fluidRow(column(
          12,
          div(
            style = "margin:12px 0;",
            DT::dataTableOutput(ns("res_table"))
          )
        ))
      )
    )
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE — SERVER
# ═══════════════════════════════════════════════════════════════════════════════
EncyclopeDIA_server <- function(id, fasta_digest) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    spin_ids <- paste0("spi", sprintf("%02d", 1:10))
    show_all <- function() lapply(spin_ids, function(s) shinyjs::show(id = s))
    hide_sp <- function(s) shinyjs::hide(id = s)
    rh <- function(fn, sid) {
      on.exit(hide_sp(sid), add = TRUE)
      fn()
    }

    enc_data <- reactiveVal(NULL)

    observeEvent(input$load_dir_btn, {
      dir_path <- trimws(input$encyclopedia_dir)
      if (dir_path == "" || !dir.exists(dir_path)) {
        showNotification("Please enter a valid directory path.", type = "error")
        return()
      }

      enc_files <- list.files(
        dir_path,
        pattern = "\\.encyclopedia2\\.txt$",
        recursive = TRUE,
        full.names = TRUE
      )
      if (length(enc_files) == 0) {
        showNotification(
          "No .encyclopedia2.txt files found in the directory.",
          type = "warning"
        )
        return()
      }

      show_all()
      withProgress(
        message = paste("Loading", length(enc_files), "EncyclopeDIA files..."),
        value = 0.5,
        {
          res <- tryCatch(
            {
              map_dfr(enc_files, function(f) {
                d <- data.table::fread(f, sep = "\t", header = TRUE)

                if (!"peptide" %in% names(d) && "Peptide" %in% names(d)) {
                  d$peptide <- d$Peptide
                }
                if (!"proteinIds" %in% names(d) && "Protein" %in% names(d)) {
                  d$proteinIds <- d$Protein
                }
                if (!"q-value" %in% names(d) && "Qvalue" %in% names(d)) {
                  d$`q-value` <- d$Qvalue
                }

                d <- d |>
                  dplyr::mutate(
                    filename = basename(f) |>
                      str_remove("\\.encyclopedia2\\.txt$"),
                    stripped_peptide = str_replace_all(
                      peptide,
                      "\\[.*?\\]|\\.|\\-",
                      ""
                    ),
                    pep_len = nchar(stripped_peptide)
                  )

                if ("PSMId" %in% names(d)) {
                  d$charge <- suppressWarnings(as.integer(str_extract(
                    d$PSMId,
                    "(?<=\\+)\\d+$"
                  )))
                  d$n_mods <- suppressWarnings(str_count(
                    d$PSMId,
                    "\\[\\+.*?\\]"
                  ))
                } else {
                  d$charge <- NA_integer_
                  d$n_mods <- NA_integer_
                }
                d
              })
            },
            error = function(e) {
              showNotification(
                paste("Error reading files:", e$message),
                type = "error"
              )
              NULL
            }
          )
          setProgress(1)
          if (!is.null(res)) {
            showNotification(
              paste("Successfully loaded", nrow(res), "records."),
              type = "message"
            )
            enc_data(res)
          }
        }
      )
    })

    filtered_data <- reactive({
      req(enc_data())
      d <- enc_data() |> dplyr::filter(`q-value` <= input$qval_filter)
      safe_seqs <- str_remove_all(d$stripped_peptide, "[^ACDEFGHIKLMNPQRSTVWY]")
      d$gravy <- sapply(safe_seqs, GRAVY)
      d$pI <- sapply(safe_seqs, calculate_pI)
      d
    })

    mapped_data <- reactive({
      d <- filtered_data()
      dig <- fasta_digest()
      if (!is.null(dig)) {
        classes <- classify_peptides(d$stripped_peptide, dig)
        d <- left_join(d, classes, by = c("stripped_peptide" = "peptide"))
      } else {
        d$classification <- "Unmapped"
        d$mapped_proteins <- NA_character_
      }
      d
    })

    output$info_box_files <- renderInfoBox({
      req(enc_data())
      n_files <- n_distinct(enc_data()$filename)
      infoBox(
        "Files Loaded",
        n_files,
        icon = icon("file-alt"),
        color = "light-blue"
      )
    })

    output$info_box_psms <- renderInfoBox({
      req(enc_data())
      total <- nrow(enc_data())
      filtered <- nrow(filtered_data())
      pct <- if (total > 0) round(filtered / total * 100, 1) else 0
      infoBox(
        "Filtered Records",
        paste0(filtered, " / ", total, " (", pct, "%)"),
        icon = icon("filter"),
        color = "green"
      )
    })

    plot_identifications_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)

      summ <- d |>
        group_by(filename) |>
        summarise(
          n_proteins = n_distinct(proteinIds),
          n_peptides = n_distinct(stripped_peptide),
          .groups = "drop"
        )

      ggplot(summ, aes(x = filename)) +
        geom_col(aes(y = n_proteins), fill = input$bar_color, width = 0.7) +
        geom_text(
          aes(y = n_proteins, label = n_proteins),
          vjust = 1.5,
          size = 4,
          fontface = "bold",
          color = "white"
        ) +
        geom_line(
          aes(y = n_peptides, group = 1),
          color = input$line_color,
          linewidth = 1.2
        ) +
        geom_point(
          aes(y = n_peptides),
          color = input$line_color,
          size = 3
        ) +
        geom_text(
          aes(y = n_peptides, label = n_peptides),
          vjust = -1,
          size = 4,
          fontface = "bold",
          color = input$line_color
        ) +
        labs(
          title = "EncyclopeDIA Identifications",
          subtitle = paste("q-value <=", input$qval_filter),
          x = "File",
          y = "Number of proteins (bar) or peptides (line)"
        ) +
        theme_enc() +
        theme(
          axis.text.x = element_text(angle = 65, hjust = 1, color = "black")
        )
    })

    plot_annotated_density_faceted_enc <- function(
      d,
      col_sym,
      title,
      xlab,
      color
    ) {
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
        facet_wrap(~filename, ncol = 3) +
        theme_enc()
    }

    plot_score_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      sc_col <- intersect(c("Score", "score", "scores"), names(d))
      if (length(sc_col) == 0) {
        return(
          ggplot() +
            annotate("text", x = 0, y = 0, label = "Score not in data") +
            theme_void()
        )
      }
      plot_annotated_density_faceted_enc(
        d,
        sc_col[1],
        "Score Distribution",
        sc_col[1],
        "#3498db"
      )
    })

    plot_pep_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      p_col <- intersect(
        c("PEP", "pep", "PosteriorErrorProbability", "posterior_error_prob"),
        names(d)
      )
      if (length(p_col) == 0) {
        return(
          ggplot() +
            annotate("text", x = 0, y = 0, label = "PEP not in data") +
            theme_void()
        )
      }
      plot_annotated_density_faceted_enc(
        d,
        p_col[1],
        "Posterior Error Probability (PEP)",
        p_col[1],
        "#e74c3c"
      )
    })

    plot_qval_obj <- reactive({
      d <- enc_data()
      req(nrow(d) > 0)
      ggplot(d, aes(x = `q-value`)) +
        geom_histogram(
          bins = 50,
          fill = "#3498db",
          color = "black",
          alpha = 0.8
        ) +
        geom_vline(
          xintercept = input$qval_filter,
          color = "red",
          linetype = "dashed",
          linewidth = 1
        ) +
        labs(title = "q-value Distribution", x = "q-value", y = "Count") +
        facet_wrap(~filename, ncol = 3) +
        theme_enc()
    })

    plot_charge_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "charge" %in% names(d))
      if (all(is.na(d$charge))) {
        return(
          ggplot() +
            annotate("text", x = 0, y = 0, label = "Charge missing") +
            theme_void()
        )
      }
      ggplot(d, aes(x = charge)) +
        geom_density(
          alpha = 0.6,
          fill = input$bar_color,
          color = "black",
          linewidth = 0.25
        ) +
        labs(x = "Charge state", y = "Density") +
        scale_x_continuous(
          breaks = seq(1, max(6, max(d$charge, na.rm = TRUE)))
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_enc()
    })

    plot_length_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0)
      ggplot(d, aes(x = pep_len)) +
        geom_density(
          alpha = 0.6,
          fill = input$bar_color,
          color = "black",
          linewidth = 0.25
        ) +
        labs(x = "Peptide length (AA)", y = "Density") +
        facet_wrap(~filename, ncol = 3) +
        theme_enc()
    })

    plot_mods_obj <- reactive({
      d <- filtered_data()
      req(nrow(d) > 0, "n_mods" %in% names(d))
      if (all(is.na(d$n_mods))) {
        return(
          ggplot() +
            annotate("text", x = 0, y = 0, label = "Mods missing") +
            theme_void()
        )
      }
      ggplot(d, aes(x = factor(n_mods))) +
        geom_bar(
          alpha = 0.8,
          fill = input$bar_color,
          color = "black",
          linewidth = 0.25
        ) +
        labs(x = "Number of modifications per peptide", y = "Count") +
        facet_wrap(~filename, ncol = 3) +
        theme_enc() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    })

    plot_gravy_obj <- reactive({
      d <- mapped_data()
      req(nrow(d) > 0, "gravy" %in% names(d))
      plot_annotated_density_faceted_enc(
        d,
        "gravy",
        "GRAVY Index Distribution",
        "GRAVY Index",
        input$line_color
      )
    })

    plot_pi_obj <- reactive({
      d <- mapped_data()
      req(nrow(d) > 0, "pI" %in% names(d))
      plot_annotated_density_faceted_enc(
        d,
        "pI",
        "Isoelectric Point (pI)",
        "pI",
        input$line_color
      )
    })

    plot_aa_freq_obj <- reactive({
      d <- mapped_data()
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
          sub_d$stripped_peptide[!is.na(sub_d$stripped_peptide)],
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
        geom_col(fill = input$bar_color, color = "black", alpha = 0.85) +
        labs(
          title = "Amino acid frequencies",
          x = "Amino acid",
          y = "Frequency (%)"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_enc()
    })

    plot_fdr_curve_obj <- reactive({
      d <- enc_data()
      req(d, nrow(d) > 0, "q-value" %in% names(d))
      q_sorted <- d |>
        dplyr::group_by(filename) |>
        dplyr::arrange(`q-value`, .by_group = TRUE) |>
        dplyr::mutate(cumulative_peptides = dplyr::row_number()) |>
        dplyr::ungroup()

      ggplot(q_sorted, aes(x = `q-value`, y = cumulative_peptides)) +
        geom_line(color = input$bar_color, linewidth = 1) +
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
          x = "q-value",
          y = "Cumulative peptide count"
        ) +
        facet_wrap(~filename, ncol = 3) +
        theme_enc()
    })

    # ── Dynamic plot height ──────────────────────────────────────────────────
    output$dynamic_plot_ui <- renderUI({
      req(filtered_data())
      d <- filtered_data()
      n <- dplyr::n_distinct(d$filename)
      h <- max(400L, ceiling(n / 3L) * 350L)
      plotOutput(ns("dynamic_plot_out"), height = paste0(h, "px"))
    })

    current_plot_obj <- eventReactive(input$run_plot, {
      req(input$plot_select)
      shinyjs::show(id = "spi_main")

      switch(
        input$plot_select,
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
        "plot_aa_freq" = plot_aa_freq_obj()
      )
    })

    output$dynamic_plot_out <- renderPlot({
      on.exit(hide_sp("spi_main"), add = TRUE)
      req(current_plot_obj())
      current_plot_obj()
    })

    output$res_table <- DT::renderDataTable({
      req(mapped_data())
      d <- mapped_data()
      num_cols <- names(d)[sapply(d, is.numeric)]
      round_cols <- intersect(
        num_cols,
        c("q-value", "PrecursorMz", "Score", "PEP", "gravy", "pI")
      )

      dt <- DT::datatable(
        d,
        options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip"),
        class = "cell-border stripe",
        rownames = FALSE,
        filter = "top"
      )
      if (length(round_cols) > 0) {
        dt <- dt |> formatRound(columns = round_cols, digits = 4)
      }
      dt
    })

    # ── Plot Download Button
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("EncyclopeDIA_", input$plot_select, "_", Sys.Date(), ".png")
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
  })
}
