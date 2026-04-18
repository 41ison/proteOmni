# ── Helper Functions ────────────────────────────────────────────────────────
compute_cv_mtx <- function(protein_matrix, group_labels) {
  unique_groups <- unique(group_labels)
  cv_results <- data.frame(row.names = rownames(protein_matrix))

  for (group in unique_groups) {
    group_samples <- which(group_labels == group)
    group_data <- protein_matrix[, group_samples, drop = FALSE]
    group_means <- rowMeans(group_data, na.rm = TRUE)
    group_sds <- apply(group_data, 1, sd, na.rm = TRUE)
    cv_values <- (group_sds / group_means) * 100
    cv_results[[paste0("CV_", group)]] <- cv_values
  }
  return(cv_results)
}

groupwise_imputation <- function(data, group_labels) {
  imputed_data <- data
  for (group in unique(group_labels)) {
    group_cols <- which(group_labels == group)
    group_data <- data[, group_cols, drop = FALSE]
    keep_rows <- !apply(group_data, 1, function(x) all(is.na(x)))
    group_data_filtered <- group_data[keep_rows, , drop = FALSE]
    if (nrow(group_data_filtered) > 0) {
      imputed_group <- missForest::missForest(
        as.data.frame(t(group_data_filtered)),
        verbose = FALSE
      )$ximp
      imputed_data[keep_rows, group_cols] <- t(imputed_group)
    }
  }
  return(imputed_data)
}

extract_power_stats <- function(contrast_fit, pwr_calc) {
  coef_names <- colnames(contrast_fit$coefficients)
  purrr::map(coef_names, \(coef) {
    limma::topTable(
      contrast_fit,
      coef = coef,
      number = Inf,
      sort.by = "none",
      adjust.method = "BH"
    ) |>
      tibble::rownames_to_column("Protein") |>
      dplyr::mutate(
        Sigma = contrast_fit$sigma,
        Min_Detectable_Log2FC = pwr_calc$d * Sigma,
        Is_reliable = abs(logFC) >= Min_Detectable_Log2FC,
        comparison = coef
      )
  }) |>
    dplyr::bind_rows()
}

# ── Valid Organism DBs ──────────────────────────────────────────────────────
org_db_choices <- c(
  "Human (org.Hs.eg.db)" = "org.Hs.eg.db",
  "Mouse (org.Mm.eg.db)" = "org.Mm.eg.db",
  "Rat (org.Rn.eg.db)" = "org.Rn.eg.db",
  "Zebrafish (org.Dr.eg.db)" = "org.Dr.eg.db",
  "C. elegans (org.Ce.eg.db)" = "org.Ce.eg.db",
  "Drosophila (org.Dm.eg.db)" = "org.Dm.eg.db",
  "Yeast (org.Sc.sgd.db)" = "org.Sc.sgd.db",
  "Arabidopsis (org.At.tair.db)" = "org.At.tair.db",
  "Bovine (org.Bt.eg.db)" = "org.Bt.eg.db",
  "Canine (org.Cf.eg.db)" = "org.Cf.eg.db",
  "Pig (org.Ss.eg.db)" = "org.Ss.eg.db",
  "Chicken (org.Gg.eg.db)" = "org.Gg.eg.db",
  "Macaque (org.Mmu.eg.db)" = "org.Mmu.eg.db"
)

# ── UI ───────────────────────────────────────────────────────────────────────
PwrQuant_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      id = ns("sidebar_content"),
      tags$div(
        style = "padding:12px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "Data Upload"
      ),
      fileInput(
        ns("matrix_file"),
        "Upload Abundance Matrix (.tsv/.csv)",
        accept = c(".tsv", ".csv", ".txt")
      ),
      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),
      tags$div(
        style = "padding:12px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "Limma Analysis"
      ),
      selectInput(
        ns("norm_method"),
        "Normalization Method",
        choices = c("scale", "quantile", "cyclicloess"),
        selected = "cyclicloess"
      ),
      selectInput(
        ns("limma_method"),
        "Limma Regression Method",
        choices = c("ls", "robust"),
        selected = "ls"
      ),
      selectInput(
        ns("trend_param"),
        "eBayes Trend",
        choices = c("TRUE", "FALSE"),
        selected = "TRUE"
      ),
      numericInput(
        ns("min_valid_pct"),
        "Min. valid values per group (%)",
        value = 0,
        min = 0,
        max = 100,
        step = 5
      ),
      uiOutput(ns("contrasts_ui")),
      tags$div(
        style = "padding:0 8px;text-align:center;",
        actionButton(
          ns("run_limma"),
          "Start limma",
          class = "btn-primary",
          style = "width:80%;font-weight:bold;margin-top:10px;margin-bottom:10px;"
        )
      ),
      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),
      tags$div(
        style = "padding:12px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "Overrepresentation Analysis"
      ),
      selectInput(
        ns("org_db"),
        "Organism Database (ORA)",
        choices = org_db_choices,
        selected = "org.Hs.eg.db"
      ),
      tags$div(
        style = "padding:0 8px;text-align:center;",
        actionButton(
          ns("run_ora"),
          "Start ORA",
          class = "btn-primary",
          style = "width:80%;font-weight:bold;margin-top:10px;margin-bottom:10px;"
        )
      ),
      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),
      tags$div(
        style = "padding:0 8px;",
        div(
          style = "margin-bottom:8px;",
          downloadButton(
            ns("download_limma"),
            "⬇ Download Limma Results",
            class = "dl-btn",
            style = "width:100%;text-align:left;"
          )
        ),
        div(
          style = "margin-bottom:8px;",
          downloadButton(
            ns("download_ora"),
            "⬇ Download ORA Results",
            class = "dl-btn",
            style = "width:100%;text-align:left;"
          )
        ),
        div(
          style = "margin-bottom:8px;",
          downloadButton(
            ns("download_all_plots"),
            "⬇ Download All Plots",
            class = "dl-btn",
            style = "width:100%;text-align:left;"
          )
        )
      )
    )
  )
}

PwrQuant_body_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    id = ns("tabs"),
    type = "tabs",
    tabPanel(
      "Metadata Mapping",
      fluidRow(
        box(
          title = "Sample metadata assignment",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          p(
            "Please assign each sample to a valid Condition and Batch. Double-click cells in the table to edit them. Ensure consistency in naming (e.g. 'Control' and 'Treatment')."
          ),
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_meta"),
              icon("spinner", class = "fa-spin")
            ),
            DT::dataTableOutput(ns("metadata_table"))
          )
        )
      )
    ),
    tabPanel(
      "Sparsity",
      fluidRow(
        box(
          title = "Data Sparsity Plot",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_sparse"),
              icon("spinner", class = "fa-spin")
            ),
            plotOutput(ns("sparsity_plot"), height = 600)
          )
        )
      )
    ),
    tabPanel(
      "Pre-processing QC",
      fluidRow(
        box(
          title = "Coefficient of Variation (CV) distributions",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_cv"),
              icon("spinner", class = "fa-spin")
            ),
            plotOutput(ns("cv_plot"), height = 500)
          )
        ),
        box(
          title = "Mean-Variance relationship",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          p(
            "Use this plot to decide whether to enable eBayes trend. If the line (the trend) is flat, trend = FALSE is fine. If the line has a distinct slope or curve, you need trend = TRUE to avoid biased results."
          ),
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_mv"),
              icon("spinner", class = "fa-spin")
            ),
            plotOutput(ns("mv_plot"), height = 500)
          )
        ),
        box(
          title = "Raw Abundance Distribution",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_raw"),
              icon("spinner", class = "fa-spin")
            ),
            plotOutput(ns("raw_boxplots"), height = 500)
          )
        ),
        box(
          title = "Normalized Abundance Distribution",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_norm"),
              icon("spinner", class = "fa-spin")
            ),
            plotOutput(ns("norm_boxplots"), height = 500)
          )
        )
      )
    ),
    tabPanel(
      "Differential Abundance",
      fluidRow(
        box(
          title = "Bland-Altman (MA) plots",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_dap"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("dap_plot_ui"))
          )
        ),
        box(
          title = "Volcano plots",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_volcano"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("volcano_plot_ui"))
          )
        ),
        box(
          title = "Protein labels",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          selectizeInput(
            ns("label_proteins"),
            "Select proteins to label in MA & Volcano plots",
            choices = NULL,
            multiple = TRUE,
            options = list(placeholder = "Type or select protein IDs...")
          )
        ),
        box(
          title = "Top differentially abundant proteins",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_barmir"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("barmir_plot_ui"))
          )
        ),
        box(
          title = "P-value distribution per contrast",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_pval_hist"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("pval_histogram_ui"))
          )
        )
      )
    ),
    tabPanel(
      "Correlation",
      fluidRow(
        box(
          title = "Advanced Contrast Correlation",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          selectInput(
            ns("adv_corr_x"),
            "X-axis contrast",
            choices = NULL
          ),
          selectInput(
            ns("adv_corr_y"),
            "Y-axis contrast",
            choices = NULL
          ),
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_adv_corr"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("adv_corr_plot_ui"))
          )
        )
      )
    ),
    tabPanel(
      "Power Statistics",
      fluidRow(
        box(
          title = "Statistical power and reliability mapping",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_pwr"),
              icon("spinner", class = "fa-spin")
            ),
            plotOutput(ns("limma_power_plot"), height = 600)
          )
        )
      )
    ),
    tabPanel(
      "Enrichment",
      fluidRow(
        box(
          title = "Overrepresentation Analysis (ORA) dotplots",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          tags$div(
            style = "padding:10px 14px;margin-bottom:8px;background:#2a3f54;border-radius:6px;color:#ddd;font-size:13px;",
            icon("info-circle", style = "color:#5bc0de;"),
            tags$b(" Protein identifier format: "),
            "For ORA to work, the row names in your abundance matrix should contain gene symbols. ",
            "Accepted formats: ",
            tags$ul(
              style = "margin-top:6px;margin-bottom:0;",
              tags$li(
                "UniProt FASTA headers with ",
                tags$code("GN=SYMBOL"),
                " tags (e.g. ",
                tags$code("sp|P12345|... GN=TP53 ..."),
                ")"
              ),
              tags$li(
                "Pipe-separated IDs (e.g. ",
                tags$code("sp|P12345|TP53_HUMAN"),
                ")"
              ),
              tags$li("Plain gene symbols (e.g. ", tags$code("TP53"), ")")
            )
          ),
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_ora"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("ora_plot_ui"))
          )
        )
      )
    ),
    tabPanel(
      "UpSet — Proteins by Condition",
      fluidRow(
        box(
          title = "UpSet plot — proteins detected per condition",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          p(
            "A protein is considered detected in a condition if it has at least one non-missing value across all samples belonging to that condition."
          ),
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_upset"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("upset_plot_ui"))
          )
        )
      )
    )
  )
}

# ── Server ───────────────────────────────────────────────────────────────────
PwrQuant_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    pq_plot_h <- reactive({
      lr <- limma_results_ev()$limma_results
      req(lr, "comparison" %in% names(lr))
      n <- length(unique(lr$comparison))
      max(400L, ceiling(n / 2) * 400L)
    })

    make_pq_plot_ui <- function(plot_id) {
      renderUI({
        req(limma_results_ev())
        plotOutput(ns(plot_id), height = paste0(pq_plot_h(), "px"))
      })
    }

    output$dap_plot_ui <- make_pq_plot_ui("dap_plot")
    output$volcano_plot_ui <- make_pq_plot_ui("volcano_plot")
    output$barmir_plot_ui <- make_pq_plot_ui("barmir_plot")
    output$pval_histogram_ui <- make_pq_plot_ui("pval_histogram")
    output$adv_corr_plot_ui <- make_pq_plot_ui("adv_corr_plot")

    ns <- session$ns

    # Placeholder for spinner IDs
    spin_ids <- c(
      "sp_meta",
      "sp_sparse",
      "sp_cv",
      "sp_mv",
      "sp_raw",
      "sp_norm",
      "sp_dap",
      "sp_volcano",
      "sp_barmir",
      "sp_corr",
      "sp_pwr",
      "sp_ora"
    )
    show_spinners <- function() {
      lapply(spin_ids, function(s) shinyjs::show(id = s))
    }
    hide_spinner <- function(sid) shinyjs::hide(id = sid)
    rh <- function(expr_fn, sid) {
      on.exit(hide_spinner(sid), add = TRUE)
      expr_fn()
    }

    # ── Matrix Ingestion ─────────────────────────────────────────────────────
    raw_matrix <- reactive({
      req(input$matrix_file)
      show_spinners()
      df <- data.table::fread(input$matrix_file$datapath, sep = "\t", header = TRUE)

      # Drop rows where the first column (representing IDs) is NA
      df <- df[!is.na(df[[1]]), ]

      # Robustly force the first column as rownames (resolving potential duplicates)
      df_clean <- df %>%
        dplyr::mutate(dplyr::across(1, ~ make.unique(as.character(.))))
      mtx <- df_clean %>%
        tibble::column_to_rownames(colnames(df_clean)[1]) %>%
        as.matrix()

      # Ensure all numeric
      mode(mtx) <- "numeric"
      return(mtx)
    })

    # ── Metadata Mapping ─────────────────────────────────────────────────────
    metadata_base <- reactive({
      req(raw_matrix())
      samples <- colnames(raw_matrix())

      # Construct default data.frame for DT with an editable Display_Name column
      df <- data.frame(
        Sample = samples,
        Display_Name = samples,
        Condition = "Control",
        Batch = "1",
        stringsAsFactors = FALSE
      )
      df
    })

    output$metadata_table <- DT::renderDataTable({
      rh(
        function() {
          req(metadata_base())
          DT::datatable(
            metadata_base(),
            editable = list(
              target = "cell",
              disable = list(columns = c(0)) # Disable editing of the Sample (original ID) column
            ),
            rownames = FALSE,
            options = list(dom = "t", paging = FALSE, ordering = FALSE),
            class = "display compact"
          )
        },
        "sp_meta"
      )
    })

    # Reactive value to capture the edited table
    meta_edit_df <- reactiveVal()

    observeEvent(metadata_base(), {
      meta_edit_df(metadata_base())
    })

    observeEvent(input$metadata_table_cell_edit, {
      info <- input$metadata_table_cell_edit
      df <- meta_edit_df()
      df[info$row, info$col + 1] <- info$value # +1 because rownames=FALSE but column index is 0-based in JS
      meta_edit_df(df)
    })

    # Helper: get display name lookup from metadata
    get_display_names <- reactive({
      req(meta_edit_df())
      meta <- meta_edit_df()
      setNames(meta$Display_Name, meta$Sample)
    })

    # ── Contrast UI ──────────────────────────────────────────────────────────
    output$contrasts_ui <- renderUI({
      req(meta_edit_df())
      conds <- unique(meta_edit_df()$Condition)
      if (length(conds) < 2) {
        return(tags$p(
          "Define at least 2 distinct conditions in the table above to configure contrasts.",
          style = "color:orange;"
        ))
      }

      tagList(
        selectizeInput(
          ns("contrast_pairs"),
          "Select Contrasts (e.g. Treat vs Control)",
          choices = NULL,
          multiple = TRUE,
          options = list(create = TRUE, placeholder = "Treat-Control")
        )
      )
    })

    observe({
      req(meta_edit_df())
      conds <- unique(meta_edit_df()$Condition)
      if (length(conds) >= 2) {
        pairs <- expand.grid(conds, conds, stringsAsFactors = FALSE)
        pairs <- pairs[pairs$Var1 != pairs$Var2, ]
        suggs <- paste0(pairs$Var1, "-", pairs$Var2)
        updateSelectizeInput(
          session,
          "contrast_pairs",
          choices = suggs,
          selected = suggs[1],
          server = TRUE
        )
      }
    })

    # ── Sparsity Plot ────────────────────────────────────────────────────────
    output$sparsity_plot <- renderPlot({
      rh(
        function() {
          req(raw_matrix(), meta_edit_df())
          pm <- as.data.frame(raw_matrix())
          disp_names <- get_display_names()
          colnames(pm) <- disp_names[colnames(pm)]
          naniar::vis_miss(pm) +
            labs(x = NULL, y = "Number of proteins") +
            theme(
              text = element_text(size = 20),
              axis.text.y = element_text(
                color = "black",
                vjust = 1,
                face = "bold"
              ),
              axis.text.x = element_text(
                angle = 90,
                color = "black",
                face = "bold"
              ),
              axis.title = element_text(color = "black", face = "bold"),
              line = element_blank(),
              panel.background = element_blank()
            )
        },
        "sp_sparse"
      )
    })

    # ── Live Preprocessing Metrics ───────────────────────────────────────────
    cv_results_live <- reactive({
      req(raw_matrix(), meta_edit_df())
      raw <- raw_matrix()
      meta <- meta_edit_df()
      group_labels <- meta$Condition[match(colnames(raw), meta$Sample)]

      log2_matrix <- log2(raw + 1)
      cv <- compute_cv_mtx(2^log2_matrix, group_labels) |>
        rownames_to_column("protein") |>
        pivot_longer(-protein, names_to = "condition", values_to = "CV") |>
        dplyr::mutate(condition = str_remove(condition, "CV_"))
      cv
    })

    # ── Live Mean-Variance Metrics (per protein, faceted by condition) ─────────
    mv_results_live <- reactive({
      req(raw_matrix(), meta_edit_df())
      raw <- raw_matrix()
      meta <- meta_edit_df()
      group_labels <- meta$Condition[match(colnames(raw), meta$Sample)]
      log2_matrix <- log2(raw + 1)

      # For each condition, compute per-protein mean and variance
      mv_list <- lapply(unique(group_labels), function(g) {
        g_cols <- which(group_labels == g)
        g_mat <- log2_matrix[, g_cols, drop = FALSE]
        data.frame(
          protein = rownames(log2_matrix),
          condition = g,
          mean_abundance = rowMeans(g_mat, na.rm = TRUE),
          variance = apply(g_mat, 1, var, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
      mv_df <- do.call(rbind, mv_list)
      mv_df[!is.na(mv_df$mean_abundance) & !is.na(mv_df$variance), ]
    })

    # ── Master Analytical Engine ─────────────────────────────────────────────
    limma_results_ev <- eventReactive(input$run_limma, {
      req(raw_matrix(), meta_edit_df(), input$contrast_pairs)

      withProgress(message = "Executing limma pipeline...", value = 0, {
        incProgress(0.1, detail = "Parsing matrix and metadata")
        # 1. Parse Data & Metadata
        raw <- raw_matrix()
        meta <- meta_edit_df()
        group_labels <- meta$Condition[match(colnames(raw), meta$Sample)]
        batch_labels <- meta$Batch[match(colnames(raw), meta$Sample)]

        # Prepare log2 transformation
        log2_matrix <- log2(raw + 1)

        # 1b. Groupwise missing value filter
        min_pct <- input$min_valid_pct
        if (!is.null(min_pct) && min_pct > 0) {
          keep_protein <- rep(TRUE, nrow(log2_matrix))
          for (g in unique(group_labels)) {
            g_cols <- which(group_labels == g)
            valid_prop <- rowMeans(!is.na(log2_matrix[, g_cols, drop = FALSE]))
            keep_protein <- keep_protein & (valid_prop >= min_pct / 100)
          }
          log2_matrix <- log2_matrix[keep_protein, , drop = FALSE]
        }

        # Identify totally missing (groupwise mask)
        rn <- rownames(log2_matrix)
        if (is.null(rn)) {
          rn <- as.character(seq_len(nrow(log2_matrix)))
        }
        group_na_mask <- data.frame(Protein = rn, stringsAsFactors = FALSE)
        for (g in unique(group_labels)) {
          g_cols <- which(group_labels == g)
          group_na_mask[[g]] <- apply(
            log2_matrix[, g_cols, drop = FALSE],
            1,
            function(x) all(is.na(x))
          )
        }

        # 2. Imputation Selection
        reg_method <- ifelse(input$limma_method == "robust", "robust", "ls")

        if (reg_method == "robust") {
          # add a coffee icon to the text for execution time humor
          incProgress(
            0.15,
            detail = "Robust regression selected: performing groupwise missForest imputation (this may take a while) ☕"
          )
          imputed_matrix <- groupwise_imputation(log2_matrix, group_labels)

          # Apply downshift for fully missing regions
          global_min <- min(imputed_matrix, na.rm = TRUE)
          set.seed(123)
          imputed_matrix_shifted <- imputed_matrix
          for (i in seq_len(nrow(imputed_matrix_shifted))) {
            na_idx <- which(is.na(imputed_matrix_shifted[i, ]))
            if (length(na_idx) > 0) {
              imputed_matrix_shifted[i, na_idx] <- rnorm(
                length(na_idx),
                mean = global_min - 1,
                sd = 0.3
              )
            }
          }
        } else {
          incProgress(
            0.2,
            detail = "Skipping imputation (running pure LS mode)"
          )
          imputed_matrix_shifted <- log2_matrix
        }

        # 3. Batch Correction
        n_batches <- length(unique(batch_labels))
        if (n_batches > 1) {
          incProgress(0.35, detail = "Calculating ComBat batch correction")
          batch_obj <- tryCatch(
            {
              sva::ComBat(
                dat = imputed_matrix_shifted,
                batch = as.character(batch_labels),
                mod = NULL,
                par.prior = TRUE,
                prior.plots = FALSE
              )
            },
            error = function(e) {
              warning(
                "ComBat failed or invalid NA structure: bypassing ComBat correction: ",
                e$message
              )
              return(imputed_matrix_shifted)
            }
          )
        } else {
          incProgress(0.35, detail = "Single batch detected — skipping ComBat")
          showNotification(
            "Only 1 batch detected: ComBat batch correction was skipped.",
            type = "warning",
            duration = 8
          )
          batch_obj <- imputed_matrix_shifted
        }
        mtx_batch_correct <- as.data.frame(batch_obj)

        # 4. Normalization
        incProgress(0.5, detail = "Normalizing Arrays")
        norm_meth <- input$norm_method
        if (norm_meth == "quantile") {
          limma_mtx <- limma::normalizeBetweenArrays(
            as.matrix(mtx_batch_correct),
            method = "quantile"
          ) %>%
            as.data.frame()
        } else if (norm_meth == "scale") {
          limma_mtx <- limma::normalizeBetweenArrays(
            as.matrix(mtx_batch_correct),
            method = "scale"
          ) %>%
            as.data.frame()
        } else {
          limma_mtx <- limma::normalizeBetweenArrays(
            as.matrix(mtx_batch_correct),
            method = "cyclicloess"
          ) %>%
            as.data.frame()
        }

        # 5. Model Matrix & Contrasts
        incProgress(0.65, detail = "Building Linear Models")
        design_groups <- factor(group_labels)
        design <- model.matrix(~ 0 + design_groups)
        colnames(design) <- str_remove(colnames(design), "design_groups")

        valid_levels <- colnames(design)
        clean_contrasts <- input$contrast_pairs
        contrast_matrix <- limma::makeContrasts(
          contrasts = clean_contrasts,
          levels = design
        )

        # 6. limma Model & eBayes
        fit_model <- limma::lmFit(limma_mtx, design, method = reg_method)
        contrast_fit <- limma::contrasts.fit(fit_model, contrast_matrix)

        use_trend <- as.logical(input$trend_param)

        if (use_trend) {
          # Exclude unestimatable proteins (NA sigma/Amean) to prevent eBayes trend calculation crashes
          valid_idx <- which(
            !is.na(contrast_fit$sigma) & !is.na(contrast_fit$Amean)
          )
          contrast_fit <- contrast_fit[valid_idx, ]
        }

        contrast_fit <- limma::eBayes(contrast_fit, trend = use_trend)

        # 7. Power Analysis
        incProgress(0.8, detail = "Calculating power metrics")
        n_repl <- median(table(group_labels))
        pwr_calc <- pwr::pwr.t.test(
          n = n_repl,
          sig.level = 0.05,
          power = 0.80,
          type = "two.sample"
        )
        power_stats <- extract_power_stats(contrast_fit, pwr_calc)

        # 8. Results Mapping
        incProgress(0.9, detail = "Isolating significance mappings")
        contrast_groups <- tibble(
          comparison = clean_contrasts,
          group_a = str_extract(clean_contrasts, "^.*(?=-)"),
          group_b = str_extract(clean_contrasts, "(?<=-).*$")
        )

        limma_results <- power_stats %>%
          left_join(contrast_groups, by = "comparison") %>%
          rowwise() %>%
          dplyr::mutate(
            imputation_driven = group_na_mask[[group_a]][match(
              Protein,
              group_na_mask$Protein
            )] |
              group_na_mask[[group_b]][match(Protein, group_na_mask$Protein)]
          ) %>%
          ungroup() %>%
          dplyr::mutate(
            status = case_when(
              imputation_driven == TRUE ~ "Not significant",
              logFC > 0 & adj.P.Val <= 0.05 & Is_reliable == TRUE ~ "Increased",
              logFC < 0 & adj.P.Val <= 0.05 & Is_reliable == TRUE ~ "Decreased",
              TRUE ~ "Not significant"
            ),
            status = factor(
              status,
              levels = c("Decreased", "Not significant", "Increased")
            )
          )

        incProgress(1.0, detail = "Complete.")

        list(
          log2_mat = log2_matrix,
          batch_mat = mtx_batch_correct,
          norm_mat = limma_mtx,
          limma_results = limma_results
        )
      })
    })

    ora_results_ev <- eventReactive(input$run_ora, {
      req(limma_results_ev())

      withProgress(
        message = "Executing ORA functional enrichment...",
        value = 0.1,
        {
          lr <- limma_results_ev()$limma_results

          incProgress(0.3, detail = "Mapping Gene Symbols")
          # Attempt multiple strategies to extract gene names from Protein IDs
          id_mapping <- lr %>%
            dplyr::mutate(
              # Strategy 1: UniProt FASTA header with GN= tag
              gene_names = str_extract(Protein, "(?<=GN=)[A-Za-z0-9_]+"),
              # Strategy 2: pipe-separated (sp|ACC|GENE_SPECIES or tr|ACC|GENE_SPECIES)
              gene_names = ifelse(
                is.na(gene_names),
                str_extract(Protein, "(?<=\\|)[A-Za-z0-9_]+(?=_)"),
                gene_names
              ),
              # Strategy 3: use the Protein ID as-is (already a gene symbol)
              gene_names = ifelse(is.na(gene_names), Protein, gene_names)
            ) %>%
            dplyr::filter(!is.na(gene_names) & gene_names != "")

          # Validate organ db
          db_str <- input$org_db
          eval(parse(text = paste0("library(", db_str, ")")))
          current_org_db <- get(db_str)

          ora_results <- id_mapping %>%
            dplyr::filter(status != "Not significant" & Is_reliable == TRUE) %>%
            group_by(status, comparison) %>%
            summarise(proteins_list = list(gene_names), .groups = "drop")

          if (nrow(ora_results) > 0) {
            incProgress(
              0.6,
              detail = "Running clusterProfiler hypergeometric tests... (this may take a while) ☕"
            )
            ora_results <- ora_results %>%
              rowwise() %>%
              dplyr::mutate(
                enrich_result = list(
                  tryCatch(
                    clusterProfiler::enrichGO(
                      gene = proteins_list,
                      OrgDb = current_org_db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2
                    ),
                    error = function(e) {
                      warning("enrichGO failed: ", e$message)
                      NULL
                    }
                  )
                )
              ) %>%
              ungroup()
          } else {
            ora_results <- NULL
          }

          incProgress(1.0, detail = "Rendering Outputs")
          ora_results
        }
      )
    })

    # ── Render Visualizations ────────────────────────────────────────────────
    output$cv_plot <- renderPlot({
      rh(
        function() {
          cv <- cv_results_live()
          ggplot(cv, aes(y = condition, x = CV, fill = condition)) +
            geom_violin(linewidth = 0.2) +
            geom_boxplot(
              width = 0.2,
              fill = "white",
              linewidth = 0.2,
              outlier.alpha = 0.5
            ) +
            geom_vline(xintercept = 20, linetype = "dashed", color = "black") +
            scale_fill_brewer(palette = "Dark2") +
            labs(y = NULL, x = "Coefficient of variation (%)", fill = NULL) +
            theme_linedraw() +
            theme(
              legend.position = "none",
              axis.title = element_text(
                size = 18,
                face = "bold",
                color = "black"
              ),
              axis.text = element_text(
                color = "black",
                size = 14,
                face = "bold"
              )
            )
        },
        "sp_cv"
      )
    })

    output$mv_plot <- renderPlot({
      rh(
        function() {
          mv <- mv_results_live()
          ggplot(mv, aes(x = mean_abundance, y = variance)) +
            geom_point(alpha = 0.3, color = "#2c3e50") +
            geom_smooth(
              method = "loess",
              se = TRUE,
              color = "firebrick",
              linewidth = 1
            ) +
            facet_wrap(~condition) +
            labs(
              x = "Mean log₂(abundance)",
              y = "Protein variance"
            ) +
            theme_bw() +
            theme(
              text = element_text(size = 16),
              strip.text = element_text(face = "bold"),
              strip.background = element_blank(),
              axis.title = element_markdown(face = "bold"),
              axis.text = element_text(color = "black", face = "bold")
            )
        },
        "sp_mv"
      )
    })

    # ── Raw Abundance Boxplots (pre-limma, visible on upload) ───────────────
    output$raw_boxplots <- renderPlot({
      rh(
        function() {
          req(raw_matrix(), meta_edit_df())
          disp_names <- get_display_names()
          raw_log2 <- log2(raw_matrix() + 1)
          raw_df <- raw_log2 |>
            as.data.frame() |>
            tibble::rownames_to_column("protein") |>
            pivot_longer(
              -protein,
              names_to = "sample",
              values_to = "abundance"
            ) |>
            dplyr::mutate(
              condition = meta_edit_df()$Condition[match(
                sample,
                meta_edit_df()$Sample
              )],
              display_name = disp_names[sample]
            )

          ggplot(
            raw_df,
            aes(
              x = display_name,
              y = abundance,
              fill = condition
            )
          ) +
            geom_boxplot(outlier.alpha = 0.3) +
            scale_fill_brewer(palette = "Dark2") +
            labs(
              x = NULL,
              y = "log₂(abundance)"
            ) +
            theme_bw() +
            theme(
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                size = 14,
                face = "bold",
                color = "black"
              ),
              axis.text.y = element_text(
                size = 14,
                face = "bold",
                color = "black"
              ),
              axis.title = element_markdown(size = 16, face = "bold"),
              legend.position = "none",
              text = element_text(size = 16)
            )
        },
        "sp_raw"
      )
    })

    # ── Normalized Abundance Boxplots (post-normalization) ─────────────────────
    output$norm_boxplots <- renderPlot({
      rh(
        function() {
          res <- limma_results_ev()
          disp_names <- get_display_names()
          norm_mtx <- res$norm_mat |>
            as.data.frame() |>
            tibble::rownames_to_column("protein") |>
            pivot_longer(
              -protein,
              names_to = "sample",
              values_to = "abundance"
            ) |>
            dplyr::mutate(
              condition = meta_edit_df()$Condition[match(
                sample,
                meta_edit_df()$Sample
              )],
              display_name = disp_names[sample]
            )

          ggplot(
            norm_mtx,
            aes(
              x = display_name,
              y = abundance,
              fill = condition
            )
          ) +
            geom_boxplot(outlier.alpha = 0.3) +
            scale_fill_brewer(palette = "Dark2") +
            labs(
              x = NULL,
              y = "log₂(abundance)"
            ) +
            theme_bw() +
            theme(
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                size = 14,
                face = "bold",
                color = "black"
              ),
              axis.text.y = element_text(
                size = 14,
                face = "bold",
                color = "black"
              ),
              axis.title = element_markdown(size = 16, face = "bold"),
              legend.position = "none",
              text = element_text(size = 16)
            )
        },
        "sp_norm"
      )
    })

    output$dap_plot <- renderPlot({
      rh(
        function() {
          lr <- limma_results_ev()$limma_results
          label_ids <- input$label_proteins
          lr <- lr %>%
            dplyr::mutate(
              label = ifelse(Protein %in% label_ids, Protein, NA_character_)
            )

          # Count annotations per facet
          count_ann <- lr |>
            dplyr::filter(status != "Not significant") |>
            dplyr::count(comparison, status)

          ggplot(lr, aes(x = AveExpr, y = logFC, color = status)) +
            geom_point(alpha = 0.3) +
            geom_smooth(
              method = "gam",
              se = FALSE,
              color = "darkblue",
              linewidth = 1,
              linetype = "dashed"
            ) +
            ggrepel::geom_text_repel(
              aes(label = label),
              size = 3,
              max.overlaps = 20,
              show.legend = FALSE,
              color = "black"
            ) +
            geom_text(
              data = count_ann,
              aes(label = paste0(status, ": ", n), color = status),
              x = Inf,
              y = Inf,
              hjust = 1.1,
              vjust = ifelse(count_ann$status == "Increased", 1.5, 3),
              size = 4,
              fontface = "bold",
              inherit.aes = FALSE,
              show.legend = FALSE
            ) +
            scale_color_manual(
              values = c(
                "Decreased" = "steelblue",
                "Not significant" = "grey60",
                "Increased" = "firebrick"
              )
            ) +
            guides(
              color = guide_legend(override.aes = list(alpha = 1, size = 3))
            ) +
            facet_wrap(~comparison, ncol = 2) +
            theme_bw() +
            labs(
              x = "log₂ average abundance",
              y = "log₂FC"
            ) +
            theme(
              text = element_text(size = 16),
              strip.text = element_text(face = "bold"),
              strip.background = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              axis.title = element_markdown(face = "bold"),
              axis.text = element_text(color = "black", face = "bold")
            )
        },
        "sp_dap"
      )
    })

    output$volcano_plot <- renderPlot({
      rh(
        function() {
          lr <- limma_results_ev()$limma_results
          label_ids <- input$label_proteins
          lr <- lr %>%
            dplyr::mutate(
              label = ifelse(Protein %in% label_ids, Protein, NA_character_)
            )

          count_ann <- lr |>
            dplyr::filter(status != "Not significant") |>
            dplyr::count(comparison, status)

          ggplot(lr, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
            geom_point(alpha = 0.3) +
            geom_hline(
              yintercept = -log10(0.05),
              linetype = "dashed",
              color = "grey40"
            ) +
            geom_vline(
              xintercept = c(-1, 1),
              linetype = "dotted",
              color = "grey40"
            ) +
            ggrepel::geom_text_repel(
              aes(label = label),
              size = 3,
              max.overlaps = 20,
              min.segment.length = 0,
              show.legend = FALSE,
              color = "black"
            ) +
            geom_text(
              data = count_ann,
              aes(label = paste0(status, ": ", n), color = status),
              x = Inf,
              y = Inf,
              hjust = 1.1,
              vjust = ifelse(count_ann$status == "Increased", 1.5, 3),
              size = 4,
              fontface = "bold",
              inherit.aes = FALSE,
              show.legend = FALSE
            ) +
            scale_color_manual(
              values = c(
                "Decreased" = "steelblue",
                "Not significant" = "grey60",
                "Increased" = "firebrick"
              )
            ) +
            guides(
              color = guide_legend(override.aes = list(alpha = 1, size = 3))
            ) +
            facet_wrap(~comparison, ncol = 2) +
            theme_bw() +
            labs(
              x = "log₂FC",
              y = "-log<sub>10</sub>(adj. p-value)"
            ) +
            theme(
              text = element_text(size = 16),
              strip.text = element_text(face = "bold"),
              strip.background = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              axis.title = element_markdown(face = "bold"),
              axis.text = element_text(color = "black", face = "bold")
            )
        },
        "sp_volcano"
      )
    })

    # Populate the protein label dropdown when limma results are available
    observe({
      lr <- tryCatch(limma_results_ev(), error = function(e) NULL)
      if (!is.null(lr)) {
        proteins <- unique(lr$limma_results$Protein)
        updateSelectizeInput(
          session,
          "label_proteins",
          choices = proteins,
          server = TRUE
        )
      }
    })

    # ── Top 20 Bar Mirror Plot ─────────────────────────────────────────────
    output$barmir_plot <- renderPlot({
      rh(
        function() {
          lr <- limma_results_ev()$limma_results
          sig <- lr %>%
            dplyr::filter(status != "Not significant" & Is_reliable == TRUE)

          if (nrow(sig) == 0) {
            return(NULL)
          }

          top_up <- sig %>%
            dplyr::filter(status == "Increased") %>%
            group_by(comparison) %>%
            dplyr::slice_max(order_by = logFC, n = 20) %>%
            ungroup()

          top_down <- sig %>%
            dplyr::filter(status == "Decreased") %>%
            group_by(comparison) %>%
            dplyr::slice_min(order_by = logFC, n = 20) %>%
            ungroup()

          top_combined <- dplyr::bind_rows(top_up, top_down) %>%
            dplyr::mutate(
              short_name = ifelse(
                nchar(Protein) > 30,
                paste0(substr(Protein, 1, 27), "..."),
                Protein
              )
            )

          ggplot(
            top_combined,
            aes(
              x = reorder(short_name, logFC),
              y = logFC,
              fill = status
            )
          ) +
            geom_col() +
            coord_flip() +
            geom_hline(yintercept = 0, linewidth = 0.5) +
            scale_fill_manual(
              values = c("Decreased" = "steelblue", "Increased" = "firebrick")
            ) +
            facet_wrap(~comparison, scales = "free_y") +
            labs(
              x = NULL,
              y = "log₂FC",
              fill = NULL
            ) +
            theme_bw() +
            theme(
              text = element_text(size = 14),
              strip.text = element_text(face = "bold"),
              strip.background = element_blank(),
              legend.position = "bottom",
              axis.title = element_markdown(face = "bold"),
              axis.text = element_text(color = "black"),
              axis.text.y = element_text(size = 9)
            )
        },
        "sp_barmir"
      )
    })

    # ── P-value distribution histogram ──────────────────────────────────────
    output$pval_histogram <- renderPlot({
      rh(
        function() {
          lr <- limma_results_ev()$limma_results
          ggplot(lr, aes(x = P.Value)) +
            geom_histogram(
              bins = 50,
              fill = "steelblue",
              color = "black",
              alpha = 0.8
            ) +
            geom_vline(
              xintercept = 0.05,
              linetype = "dashed",
              color = "red",
              linewidth = 0.7
            ) +
            annotate(
              "text",
              x = 0.05,
              y = Inf,
              label = "p = 0.05",
              vjust = 2,
              hjust = -0.1,
              color = "red",
              fontface = "bold",
              size = 3.5
            ) +
            facet_wrap(~comparison, ncol = 2) +
            labs(
              title = "Raw p-value distribution",
              x = "p-value",
              y = "Count"
            ) +
            theme_bw() +
            theme(
              text = element_text(size = 14),
              strip.text = element_text(face = "bold"),
              strip.background = element_blank(),
              axis.text = element_text(color = "black", face = "bold"),
              axis.title = element_text(face = "bold"),
              plot.title = element_text(face = "bold", hjust = 0.5)
            )
        },
        "sp_pval_hist"
      )
    })

    # ── Advanced Contrast Correlation ───────────────────────────────────────
    observe({
      lr <- tryCatch(limma_results_ev(), error = function(e) NULL)
      if (!is.null(lr)) {
        contrasts <- unique(lr$limma_results$comparison)
        updateSelectInput(
          session,
          "adv_corr_x",
          choices = contrasts,
          selected = contrasts[1]
        )
        updateSelectInput(
          session,
          "adv_corr_y",
          choices = contrasts,
          selected = if (length(contrasts) >= 2) contrasts[2] else contrasts[1]
        )
      }
    })

    output$adv_corr_plot <- renderPlot({
      rh(
        function() {
          lr <- limma_results_ev()$limma_results
          req(
            input$adv_corr_x,
            input$adv_corr_y,
            input$adv_corr_x != input$adv_corr_y
          )

          dx <- lr |>
            dplyr::filter(comparison == input$adv_corr_x) |>
            dplyr::select(Protein, logFC, status) |>
            dplyr::rename(logFC_x = logFC, status_x = status)
          dy <- lr |>
            dplyr::filter(comparison == input$adv_corr_y) |>
            dplyr::select(Protein, logFC, status) |>
            dplyr::rename(logFC_y = logFC, status_y = status)

          merged <- dplyr::inner_join(dx, dy, by = "Protein") |>
            dplyr::mutate(
              class = dplyr::case_when(
                (status_x == status_y) &
                  (status_x != "Not significant") ~ "Concordant",
                (status_x != "Not significant") &
                  (status_y != "Not significant") &
                  (status_x != status_y) ~ "Inverse",
                xor(
                  status_x == "Not significant",
                  status_y == "Not significant"
                ) ~ "Mismatch",
                TRUE ~ "Not significant"
              )
            )

          rho <- cor(
            merged$logFC_x,
            merged$logFC_y,
            use = "complete.obs",
            method = "spearman"
          )

          ggplot(merged, aes(x = logFC_x, y = logFC_y, color = class)) +
            geom_point(alpha = 0.5) +
            geom_smooth(
              data = merged |> dplyr::filter(class == "Concordant"),
              method = "lm",
              se = FALSE,
              linewidth = 1,
              linetype = "solid"
            ) +
            geom_smooth(
              data = merged |> dplyr::filter(class == "Inverse"),
              method = "lm",
              se = FALSE,
              linewidth = 1,
              linetype = "solid"
            ) +
            geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
            geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
            scale_color_manual(
              values = c(
                "Concordant" = "#1ca957",
                "Inverse" = "#00B3FF",
                "Mismatch" = "#E02121",
                "Not significant" = "#CCCCCCFE"
              )
            ) +
            annotate(
              "text",
              x = Inf,
              y = Inf,
              label = paste0("\u03c1 = ", round(rho, 3)),
              hjust = 1.1,
              vjust = 1.5,
              size = 6,
              fontface = "bold",
              color = "black"
            ) +
            labs(
              title = paste(input$adv_corr_x, "vs", input$adv_corr_y),
              x = paste0("log₂FC: ", input$adv_corr_x),
              y = paste0("log₂FC: ", input$adv_corr_y),
              color = "Classification"
            ) +
            theme_bw() +
            theme(
              text = element_text(size = 14),
              plot.title = element_text(face = "bold", hjust = 0.5),
              axis.text = element_text(color = "black", face = "bold"),
              axis.title = element_markdown(face = "bold"),
              legend.position = "bottom",
              legend.title = element_text(face = "bold", hjust = 0.5),
              legend.title.position = "top"
            )
        },
        "sp_adv_corr"
      )
    })

    output$limma_power_plot <- renderPlot({
      rh(
        function() {
          pw <- limma_results_ev()$limma_results
          ggplot(pw, aes(x = Sigma, y = abs(logFC))) +
            geom_point(aes(color = Is_reliable), alpha = 0.3) +
            geom_line(
              aes(y = Min_Detectable_Log2FC),
              color = "red",
              linetype = "dashed",
              linewidth = 1
            ) +
            labs(
              x = "Residual standard deviation (Sigma)",
              y = "Observed |log₂FC|",
              color = "Significant difference with 80% power"
            ) +
            scale_color_brewer(palette = "Dark2") +
            guides(
              color = guide_legend(override.aes = list(alpha = 1, size = 3))
            ) +
            facet_wrap(~comparison, scales = "free", ncol = 2) +
            theme_bw() +
            theme(
              text = element_text(size = 16),
              legend.position = "bottom",
              axis.title = element_markdown(face = "bold"),
              strip.text = element_text(face = "bold"),
              strip.background = element_blank()
            )
        },
        "sp_pwr"
      )
    })

    output$ora_plot_ui <- renderUI({
      ora <- tryCatch(ora_results_ev(), error = function(e) NULL)
      shinyjs::hide("sp_ora")
      if (is.null(ora) || nrow(ora) == 0) {
        return(tags$p(
          style = "color:#adb5bd;text-align:center;padding:20px;",
          "No significant & reliable proteins found for ORA. Run limma first."
        ))
      }
      n_panels <- nrow(ora)
      dynamic_h <- max(600, n_panels * 350)
      plotOutput(ns("ora_plot"), height = paste0(dynamic_h, "px"))
    })

    output$ora_plot <- renderPlot({
      rh(
        function() {
          ora <- ora_results_ev()
          if (is.null(ora) || nrow(ora) == 0) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "No significant & reliable proteins found for ORA.\nRun limma first and ensure significant proteins exist.",
                  size = 6,
                  color = "grey40",
                  hjust = 0.5
                ) +
                theme_void()
            )
          }
          plots <- lapply(seq_len(nrow(ora)), function(i) {
            row_data <- ora[i, ]
            er <- row_data$enrich_result[[1]]
            # Check that the enrichResult is valid and has actual enriched terms
            if (!is.null(er) && is(er, "enrichResult")) {
              sig_terms <- tryCatch(
                sum(er@result$p.adjust < 0.05, na.rm = TRUE),
                error = function(e) 0
              )
              if (sig_terms > 0) {
                tryCatch(
                  {
                    enrichplot::dotplot(er, font.size = 14) +
                      labs(
                        title = row_data$comparison,
                        subtitle = paste("Status:", row_data$status)
                      ) +
                      theme(
                        plot.title = element_text(face = "bold", size = 16),
                        plot.subtitle = element_text(size = 13),
                        axis.text = element_text(size = 12, color = "black"),
                        axis.title = element_text(size = 14, face = "bold"),
                        legend.position = "right",
                        legend.text = element_text(size = 12)
                      )
                  },
                  error = function(e) {
                    warning(
                      "dotplot rendering failed for ",
                      row_data$comparison,
                      ": ",
                      e$message
                    )
                    NULL
                  }
                )
              } else {
                NULL
              }
            } else {
              NULL
            }
          })
          plots <- Filter(Negate(is.null), plots)
          if (length(plots) == 0) {
            return(
              ggplot() +
                annotate(
                  "text",
                  x = 0.5,
                  y = 0.5,
                  label = "ORA completed but no enriched GO terms were found.\nEnsure your protein identifiers contain valid gene symbols.",
                  size = 6,
                  color = "grey40",
                  hjust = 0.5
                ) +
                theme_void()
            )
          }
          wrap_plots(plots, ncol = 2)
        },
        "sp_ora"
      )
    })

    # ── Downloads ────────────────────────────────────────────────────────────
    output$download_limma <- downloadHandler(
      filename = function() {
        paste0("limma_results_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        readr::write_tsv(limma_results_ev()$limma_results, file)
      }
    )

    output$download_ora <- downloadHandler(
      filename = function() {
        paste0("ORA_results_", Sys.Date(), ".zip")
      },
      content = function(file) {
        ora <- ora_results_ev()
        req(!is.null(ora))
        temp_dir <- tempdir()
        files_to_zip <- c()
        for (i in seq_len(nrow(ora))) {
          row_data <- ora[i, ]
          res_tbl <- row_data$enrich_result[[1]]@result
          fname <- paste0(
            temp_dir,
            "/ORA_",
            row_data$comparison,
            "_",
            row_data$status,
            ".tsv"
          )
          readr::write_tsv(res_tbl, fname)
          files_to_zip <- c(files_to_zip, fname)
        }
        zip::zipr(file, files_to_zip)
      }
    )
    # ── Download All Plots ──────────────────────────────────────────────────
    output$download_all_plots <- downloadHandler(
      filename = function() {
        paste0("PwrQuant_plots_", Sys.Date(), ".zip")
      },
      content = function(file) {
        temp_dir <- tempdir()
        plot_dir <- file.path(temp_dir, "PwrQuant_plots")
        dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

        disp_names <- get_display_names()
        meta <- meta_edit_df()

        # Helper to safely save a ggplot
        safe_save <- function(filename, plot_expr, w = 12, h = 8) {
          tryCatch(
            {
              p <- plot_expr()
              if (!is.null(p)) {
                ggplot2::ggsave(
                  file.path(plot_dir, filename),
                  plot = p,
                  width = w,
                  height = h,
                  dpi = 300,
                  bg = "white"
                )
              }
            },
            error = function(e) {
              warning("Could not save ", filename, ": ", e$message)
            }
          )
        }

        # Helper to safely save a base R plot
        safe_save_base <- function(filename, plot_expr, w = 12, h = 8) {
          tryCatch(
            {
              grDevices::png(
                file.path(plot_dir, filename),
                width = w,
                height = h,
                units = "in",
                res = 300
              )
              plot_expr()
              grDevices::dev.off()
            },
            error = function(e) {
              tryCatch(grDevices::dev.off(), error = function(e2) NULL)
              warning("Could not save ", filename, ": ", e$message)
            }
          )
        }

        # 1. Sparsity plot
        safe_save_base("sparsity_plot.png", function() {
          pm <- as.data.frame(raw_matrix())
          colnames(pm) <- disp_names[colnames(pm)]
          p <- naniar::vis_miss(pm) +
            labs(x = NULL, y = "Number of proteins") +
            theme(
              text = element_text(size = 20),
              axis.text.y = element_text(
                color = "black",
                vjust = 1,
                face = "bold"
              ),
              axis.text.x = element_text(
                angle = 90,
                color = "black",
                face = "bold"
              ),
              axis.title = element_text(color = "black", face = "bold"),
              line = element_blank(),
              panel.background = element_blank()
            )
          print(p)
        })

        # 2. CV plot
        safe_save("cv_plot.png", function() {
          cv <- cv_results_live()
          ggplot(cv, aes(y = condition, x = CV, fill = condition)) +
            geom_violin(linewidth = 0.2) +
            geom_boxplot(
              width = 0.2,
              fill = "white",
              linewidth = 0.2,
              outlier.alpha = 0.5
            ) +
            geom_vline(xintercept = 20, linetype = "dashed", color = "black") +
            scale_fill_brewer(palette = "Dark2") +
            labs(y = NULL, x = "Coefficient of variation (%)", fill = NULL) +
            theme_linedraw() +
            theme(
              legend.position = "none",
              axis.title = element_text(
                size = 18,
                face = "bold",
                color = "black"
              ),
              axis.text = element_text(
                color = "black",
                size = 14,
                face = "bold"
              )
            )
        })

        # 3. Mean-Variance plot
        safe_save("mv_plot.png", function() {
          mv <- mv_results_live()
          ggplot(mv, aes(x = mean_abundance, y = variance)) +
            geom_point(alpha = 0.3, color = "#2c3e50") +
            geom_smooth(
              method = "loess",
              se = TRUE,
              color = "firebrick",
              linewidth = 1
            ) +
            facet_wrap(~condition) +
            labs(
              x = "Mean log₂(abundance)",
              y = "Protein variance"
            ) +
            theme_bw() +
            theme(
              text = element_text(size = 16),
              strip.text = element_text(face = "bold", size = 14),
              strip.background = element_blank(),
              axis.title = element_markdown(face = "bold", size = 16),
              axis.text = element_text(
                color = "black",
                face = "bold",
                size = 14
              )
            )
        })

        # 4. Raw boxplots
        safe_save("raw_boxplots.png", function() {
          raw_log2 <- log2(raw_matrix() + 1)
          raw_df <- raw_log2 |>
            as.data.frame() |>
            tibble::rownames_to_column("protein") |>
            tidyr::pivot_longer(
              -protein,
              names_to = "sample",
              values_to = "abundance"
            ) |>
            dplyr::mutate(
              condition = meta$Condition[match(sample, meta$Sample)],
              display_name = disp_names[sample]
            )
          ggplot(
            raw_df,
            aes(x = display_name, y = abundance, fill = condition)
          ) +
            geom_boxplot(outlier.alpha = 0.3) +
            scale_fill_brewer(palette = "Dark2") +
            labs(x = NULL, y = "log₂(abundance)") +
            theme_bw() +
            theme(
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                size = 14,
                face = "bold",
                color = "black"
              ),
              axis.text.y = element_text(
                size = 14,
                face = "bold",
                color = "black"
              ),
              axis.title = element_markdown(size = 16, face = "bold"),
              legend.position = "none",
              text = element_text(size = 16)
            )
        })

        # 5-9. limma-dependent plots (may not exist yet)
        tryCatch(
          {
            res <- limma_results_ev()
            lr <- res$limma_results

            # 5. Normalized boxplots
            safe_save("norm_boxplots.png", function() {
              batch_mtx <- res$norm_mat |>
                as.data.frame() |>
                tibble::rownames_to_column("protein") |>
                tidyr::pivot_longer(
                  -protein,
                  names_to = "sample",
                  values_to = "abundance"
                ) |>
                dplyr::mutate(
                  condition = meta$Condition[match(sample, meta$Sample)],
                  display_name = disp_names[sample]
                )
              ggplot(
                batch_mtx,
                aes(x = display_name, y = abundance, fill = condition)
              ) +
                geom_boxplot(outlier.alpha = 0.3) +
                scale_fill_brewer(palette = "Dark2") +
                labs(x = NULL, y = "log₂(abundance)") +
                theme_bw() +
                theme(
                  axis.text.x = element_text(
                    angle = 45,
                    hjust = 1,
                    size = 14,
                    face = "bold",
                    color = "black"
                  ),
                  axis.text.y = element_text(
                    size = 14,
                    face = "bold",
                    color = "black"
                  ),
                  axis.title = element_markdown(size = 16, face = "bold"),
                  legend.position = "none",
                  text = element_text(size = 16)
                )
            })

            # 6. MA plot
            safe_save(
              "ma_plot.png",
              function() {
                label_ids <- input$label_proteins
                lr <- lr %>%
                  dplyr::mutate(
                    label = ifelse(
                      Protein %in% label_ids,
                      Protein,
                      NA_character_
                    )
                  )
                ggplot(lr, aes(x = AveExpr, y = logFC, color = status)) +
                  geom_point(alpha = 0.3) +
                  ggrepel::geom_text_repel(
                    aes(label = label),
                    size = 3,
                    max.overlaps = 20,
                    show.legend = FALSE,
                    color = "black"
                  ) +
                  guides(
                    color = guide_legend(
                      override.aes = list(alpha = 1, size = 3)
                    )
                  ) +
                  geom_smooth(
                    method = "gam",
                    se = TRUE,
                    color = "black",
                    linewidth = 1
                  ) +
                  scale_color_manual(
                    values = c(
                      "Decreased" = "steelblue",
                      "Not significant" = "grey60",
                      "Increased" = "firebrick"
                    )
                  ) +
                  facet_wrap(~comparison, ncol = 2) +
                  theme_bw() +
                  labs(x = "log₂ average abundance", y = "log₂FC") +
                  theme(
                    text = element_text(size = 16),
                    strip.text = element_text(face = "bold", size = 14),
                    strip.background = element_blank(),
                    legend.position = "bottom",
                    axis.title = element_markdown(face = "bold", size = 16),
                    axis.text = element_text(
                      color = "black",
                      face = "bold",
                      size = 14
                    )
                  )
              },
              w = 14,
              h = 10
            )

            # 7. Volcano plot
            safe_save(
              "volcano_plot.png",
              function() {
                label_ids <- input$label_proteins
                lr <- lr %>%
                  dplyr::mutate(
                    label = ifelse(
                      Protein %in% label_ids,
                      Protein,
                      NA_character_
                    )
                  )
                ggplot(
                  lr,
                  aes(x = logFC, y = -log10(adj.P.Val), color = status)
                ) +
                  geom_point(alpha = 0.3) +
                  guides(
                    color = guide_legend(
                      override.aes = list(alpha = 1, size = 3)
                    )
                  ) +
                  geom_hline(
                    yintercept = -log10(0.05),
                    linetype = "dashed",
                    color = "grey40"
                  ) +
                  geom_vline(
                    xintercept = c(-1, 1),
                    linetype = "dotted",
                    color = "grey40"
                  ) +
                  ggrepel::geom_text_repel(
                    aes(label = label),
                    size = 3,
                    max.overlaps = 20,
                    min.segment.length = 0,
                    show.legend = FALSE,
                    color = "black"
                  ) +
                  scale_color_manual(
                    values = c(
                      "Decreased" = "steelblue",
                      "Not significant" = "grey60",
                      "Increased" = "firebrick"
                    )
                  ) +
                  facet_wrap(~comparison, ncol = 2) +
                  theme_bw() +
                  labs(
                    x = "log₂FC",
                    y = "-log<sub>10</sub>(adj. p-value)"
                  ) +
                  theme(
                    text = element_text(size = 16),
                    strip.text = element_text(face = "bold", size = 14),
                    strip.background = element_blank(),
                    legend.position = "bottom",
                    axis.title = element_markdown(face = "bold", size = 16),
                    axis.text = element_text(
                      color = "black",
                      face = "bold",
                      size = 14
                    )
                  )
              },
              w = 14,
              h = 10
            )

            # 8. Bar mirror plot
            safe_save(
              "barmir_plot.png",
              function() {
                sig <- lr %>%
                  dplyr::filter(
                    status != "Not significant" & Is_reliable == TRUE
                  )
                if (nrow(sig) == 0) {
                  return(NULL)
                }
                top_up <- sig %>%
                  dplyr::filter(status == "Increased") %>%
                  dplyr::group_by(comparison) %>%
                  dplyr::slice_max(order_by = logFC, n = 20) %>%
                  dplyr::ungroup()
                top_down <- sig %>%
                  dplyr::filter(status == "Decreased") %>%
                  dplyr::group_by(comparison) %>%
                  dplyr::slice_min(order_by = logFC, n = 20) %>%
                  dplyr::ungroup()
                top_combined <- dplyr::bind_rows(top_up, top_down) %>%
                  dplyr::mutate(
                    short_name = ifelse(
                      nchar(Protein) > 30,
                      paste0(substr(Protein, 1, 27), "..."),
                      Protein
                    )
                  )
                ggplot(
                  top_combined,
                  aes(x = reorder(short_name, logFC), y = logFC, fill = status)
                ) +
                  geom_col() +
                  coord_flip() +
                  geom_hline(yintercept = 0, linewidth = 0.5) +
                  scale_fill_manual(
                    values = c(
                      "Decreased" = "steelblue",
                      "Increased" = "firebrick"
                    )
                  ) +
                  facet_wrap(~comparison, scales = "free_y") +
                  labs(x = NULL, y = "log₂FC", fill = NULL) +
                  theme_bw() +
                  theme(
                    text = element_text(size = 14),
                    strip.text = element_text(face = "bold", size = 14),
                    strip.background = element_blank(),
                    legend.position = "bottom",
                    axis.title = element_markdown(face = "bold", size = 14),
                    axis.text = element_text(color = "black", size = 12),
                    axis.text.y = element_text(size = 9)
                  )
              },
              w = 14,
              h = 10
            )

            # 9. Power plot
            safe_save("power_plot.png", function() {
              ggplot(lr, aes(x = Sigma, y = abs(logFC))) +
                geom_point(aes(color = Is_reliable), alpha = 0.3) +
                guides(
                  color = guide_legend(override.aes = list(alpha = 1, size = 3))
                ) +
                geom_line(
                  aes(y = Min_Detectable_Log2FC),
                  color = "red",
                  linetype = "dashed",
                  linewidth = 1
                ) +
                labs(
                  x = "Residual standard deviation (Sigma)",
                  y = "Observed |log₂FC|",
                  color = "80% power"
                ) +
                scale_color_brewer(palette = "Dark2") +
                facet_wrap(~comparison, scales = "free", ncol = 2) +
                theme_bw() +
                theme(
                  text = element_text(size = 16),
                  legend.position = "bottom",
                  axis.title = element_markdown(face = "bold", size = 16),
                  axis.text = element_text(
                    color = "black",
                    face = "bold",
                    size = 14
                  ),
                  strip.text = element_text(face = "bold", size = 14),
                  strip.background = element_blank()
                )
            })
          },
          error = function(e) {
            warning("limma results not available for plot export: ", e$message)
          }
        )

        # 10. ORA plot
        tryCatch(
          {
            ora <- ora_results_ev()
            if (!is.null(ora) && nrow(ora) > 0) {
              plots <- lapply(seq_len(nrow(ora)), function(i) {
                row_data <- ora[i, ]
                er <- row_data$enrich_result[[1]]
                if (!is.null(er) && nrow(er@result) > 0) {
                  enrichplot::dotplot(er, font.size = 14) +
                    labs(
                      title = row_data$comparison,
                      subtitle = paste("Status:", row_data$status)
                    ) +
                    theme(
                      plot.title = element_text(face = "bold", size = 16),
                      legend.position = "right"
                    )
                } else {
                  NULL
                }
              })
              plots <- Filter(Negate(is.null), plots)
              if (length(plots) > 0) {
                p <- patchwork::wrap_plots(plots, ncol = 2)
                ggplot2::ggsave(
                  file.path(plot_dir, "ora_plot.png"),
                  plot = p,
                  width = 16,
                  height = 10,
                  dpi = 300,
                  bg = "white"
                )
              }
            }
          },
          error = function(e) {
            warning("ORA results not available for plot export: ", e$message)
          }
        )

        # Zip all generated PNGs
        png_files <- list.files(
          plot_dir,
          pattern = "\\.png$",
          full.names = TRUE
        )
        if (length(png_files) > 0) {
          zip::zipr(file, png_files)
        } else {
          showNotification("No plots available to download.", type = "warning")
        }
      }
    )

    # ── UpSet — Proteins by Condition ──────────────────────────────────────────────
    upset_sets <- reactive({
      req(raw_matrix(), meta_edit_df())
      mat <- raw_matrix()
      meta <- meta_edit_df()

      # For each condition: union of proteins with at least one non-NA value
      conditions <- unique(meta$Condition)
      sets <- lapply(conditions, function(cond) {
        samples <- meta$Sample[meta$Condition == cond]
        samples <- intersect(samples, colnames(mat))
        if (length(samples) == 0) {
          return(character(0))
        }
        sub <- mat[, samples, drop = FALSE]
        rownames(sub)[apply(sub, 1, function(x) any(!is.na(x)))]
      })
      names(sets) <- conditions
      sets[vapply(sets, length, integer(1)) > 0]
    })

    output$upset_plot_ui <- renderUI({
      req(upset_sets())
      n_cond <- length(upset_sets())
      h <- max(500, 300 + n_cond * 40)
      plotOutput(ns("upset_plot"), height = paste0(h, "px"))
    })

    output$upset_plot <- renderPlot({
      shinyjs::show(id = "sp_upset")
      on.exit(shinyjs::hide(id = "sp_upset"), add = TRUE)
      sets <- upset_sets()
      req(length(sets) >= 2)

      m <- ComplexHeatmap::make_comb_mat(sets)
      ComplexHeatmap::UpSet(
        m,
        comb_order = order(ComplexHeatmap::comb_size(m), decreasing = TRUE),
        top_annotation = ComplexHeatmap::upset_top_annotation(
          m,
          add_numbers = TRUE,
          annotation_name_rot = 0,
          annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
          gp = grid::gpar(fill = "#1B4965")
        ),
        right_annotation = ComplexHeatmap::upset_right_annotation(
          m,
          add_numbers = TRUE,
          annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
          gp = grid::gpar(fill = "#5FA8D3")
        ),
        row_names_gp = grid::gpar(fontsize = 11, fontface = "bold"),
        column_title = "Proteins detected per condition",
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
      )
    })
  })
}
