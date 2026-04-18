## QC4DIANN Shiny Module — QC4DIANN_v1.1.0 module for proteOmni

# ── Sidebar UI (moved to dashboard sidebar via switch in proteOmni.R) ──
QC4DIANN_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "sidebar-section-label", "Data Input"),
    tags$div(
      style = "padding:0 16px;",
      fileInput(
        ns("report"),
        "Parquet Report",
        accept = ".parquet",
        placeholder = "report.parquet"
      ),
      fileInput(
        ns("fasta_file"),
        "Protein FASTA (optional)",
        accept = c(".fasta", ".fa", ".faa"),
        placeholder = "proteome.fasta"
      ),
      numericInput(
        ns("missed_cleavages"),
        "Max Missed Cleavages (FASTA digest)",
        value = 0,
        min = 0,
        max = 3,
        step = 1
      )
    ),
    tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

    tags$div(class = "sidebar-section-label", "QuantUMS Filters"),
    tags$div(
      style = "padding:0 16px;",
      sliderInput(
        ns("PG.MaxLFQ.Quality"),
        "PG MaxLFQ Quality",
        min = 0,
        max = 1,
        value = 0.75,
        step = 0.05
      ),
      sliderInput(
        ns("Empirical.Quality"),
        "Empirical Quality",
        min = 0,
        max = 1,
        value = 0,
        step = 0.05
      )
    ),
    tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

    tags$div(class = "sidebar-section-label", "Interactive Viewer"),
    tags$div(
      style = "padding:0 16px;",
      selectInput(ns("xcol"), "X Sample", choices = NULL),
      selectInput(ns("ycol"), "Y Sample", choices = NULL),
      numericInput(
        ns("efa_factors"),
        "EFA Latent Factors",
        value = 6,
        min = 1,
        max = 20,
        step = 1
      )
    ),
    tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

    tags$div(class = "sidebar-section-label", "Downloads"),
    tags$div(
      style = "padding:0 16px;",
      downloadButton(
        ns("download_matrix"),
        tagList(icon("download"), " Filtered Protein Matrix"),
        class = "dl-btn",
        style = "width:100%;margin-bottom:8px;"
      )
    )
  )
}

# ── Body UI (2-column layout) ───────────────────────────────────────────────
QC4DIANN_ui <- function(id) {
  ns <- NS(id)
  WIDE_HEIGHT <- "700px"

  tagList(
    fluidRow(infoBoxOutput(ns("info_box1"), width = 12)),

    tabsetPanel(
      id = ns("tabs"),
      type = "tabs",

      # ── TAB 1: QC Filters & Distributions (2-column) ──
      tabPanel(
        "QC Filters & Distributions",
        fluidRow(
          column(
            4,
            box(
              title = "QC Plot Controls",
              status = "primary",
              solidHeader = TRUE,
              width = NULL,
              selectInput(
                ns("qc_plot_select"),
                "Select Graphic",
                choices = c(
                  "Reconstruction of XIC" = "plot1",
                  "Density of Ions (m/z vs RT)" = "plot2",
                  "Retention Time Prediction Error" = "plot3",
                  "Charge State Distribution" = "plot4",
                  "Peptide Length Distribution" = "plot5",
                  "Peptides per Sample" = "plot6",
                  "Proteins per Sample" = "plot7",
                  "Cysteine Counts in Peptides" = "plot15",
                  "Sparsity Profile (Missing Values %)" = "plot8",
                  "Missing Values vs Median Abundance" = "plot9",
                  "Abundance Before and After MAD Normalization" = "plot10",
                  "Missed Cleavage Sites" = "plot11",
                  "MS1 Profile Correlation" = "plot12",
                  "QuantUMS Score Distributions" = "plot13",
                  "Gene Quantity Distribution" = "plot14"
                )
              ),
              actionButton(
                ns("run_qc_plot"),
                "Plot Selected Graphic",
                icon = icon("chart-bar"),
                class = "btn-primary",
                style = "width:100%;margin-bottom:8px;"
              ),
              downloadButton(
                ns("download_qc_plot"),
                tagList(icon("download"), " Download (.png)"),
                class = "dl-btn",
                style = "width:100%;"
              )
            )
          ),
          column(
            8,
            box(
              title = "Dynamic QC Plot View",
              status = "primary",
              solidHeader = TRUE,
              width = NULL,
              div(
                class = "plot-wrap",
                tags$div(
                  class = "spinner-overlay",
                  id = ns("sp_qc_main"),
                  icon("spinner", class = "fa-spin")
                ),
                uiOutput(ns("dynamic_qc_plot_ui"))
              )
            )
          )
        )
      ),

      # ── TAB 2: Interactive Viewer ──
      tabPanel(
        "Interactive Viewer",
        fluidRow(
          box(
            title = "Sample Correlation — Non-normalized log₂(Intensity)",
            status = "primary",
            solidHeader = TRUE,
            width = 7,
            plotlyOutput(ns("Corr"), height = WIDE_HEIGHT)
          ),
          tabBox(
            title = "Similarity Metrics",
            side = "right",
            width = 5,
            tabPanel(
              "Cosine Similarity",
              plotOutput(ns("cosine_similarity"), height = "630px")
            ),
            tabPanel(
              "Euclidean Distance",
              plotOutput(ns("euclidean_distance"), height = "630px")
            ),
            tabPanel(
              "Jaccard Similarity",
              plotOutput(ns("jaccard_similarity"), height = "630px")
            )
          )
        ),
        fluidRow(box(
          title = "QuantUMS Score Distribution (3D)",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          plotlyOutput(ns("QuantUMS_dist"), height = WIDE_HEIGHT)
        )),
        fluidRow(
          box(
            title = "Principal Component Analysis (PCA)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp_pca"),
                icon("spinner", class = "fa-spin")
              ),
              plotlyOutput(ns("PCA"), height = WIDE_HEIGHT)
            )
          ),
          box(
            title = "Sample Correlation Matrix",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp_hm"),
                icon("spinner", class = "fa-spin")
              ),
              plotOutput(ns("plot_ggpairs"), height = WIDE_HEIGHT)
            )
          ),
          box(
            title = "Exploratory Factor Analysis (EFA)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("spefa"),
                icon("spinner", class = "fa-spin")
              ),
              plotOutput(ns("plot_efa"), height = WIDE_HEIGHT)
            )
          )
        )
      ),

      # ── TAB 3: Peptide Mapping ──
      tabPanel(
        "Peptide Mapping",
        fluidRow(uiOutput(ns("pep_map_banner"))),
        fluidRow(box(
          title = "Protein Sequence View",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          fluidRow(
            column(
              3,
              selectizeInput(
                ns("selected_protein"),
                "Select Protein:",
                choices = NULL
              )
            ),
            column(
              3,
              selectInput(
                ns("selected_sample_psv"),
                "Select Sample:",
                choices = NULL
              )
            ),
            column(
              3,
              numericInput(
                ns("aa_per_line"),
                "AAs per line:",
                value = 50,
                min = 10,
                max = 200
              )
            ),
            column(
              3,
              checkboxInput(
                ns("show_sequence"),
                "Show Sequence Letters",
                value = TRUE
              )
            )
          ),
          uiOutput(ns("protein_coverage_plot_ui"))
        )),
        fluidRow(box(
          title = "Amino Acid Frequencies (Expected vs Identified)",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_aa"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("plot_aa_freq_ui"))
          )
        )),
        fluidRow(box(
          title = "Proteotypic vs Shared Peptides — Counts",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_pc"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("plot_pep_type_count_ui"))
          )
        )),
        fluidRow(box(
          title = "Proteotypic vs Shared Peptides — Proportions",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_pp"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("plot_pep_type_prop_ui"))
          )
        )),
        fluidRow(box(
          title = "Peptide Mapping Summary Table",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          DTOutput(ns("peptide_table"))
        ))
      ),

      # ── TAB 4: Protease Specificity ──
      tabPanel(
        "Protease Specificity",
        fluidRow(uiOutput(ns("prot_spec_banner"))),
        fluidRow(box(
          title = "Schechter-Berger Sequence Logo (P4–P4') — All Runs Combined",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          plotOutput(ns("seqlogo_all"), height = "420px")
        )),
        fluidRow(box(
          title = "Schechter-Berger Sequence Logo — Per Run",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          uiOutput(ns("seqlogo_perrun_ui"))
        ))
      )
    )
  )
}


# ── Server ──────────────────────────────────────────────────────────────────
QC4DIANN_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # ── Namespace + helpers declared FIRST ──
    ns <- session$ns

    show_spinner <- function(sid) shinyjs::show(id = sid)
    hide_spinner <- function(sid) shinyjs::hide(id = sid)

    pal <- c(
      proteotypic = "#2D6A4F",
      shared = "#E76F51",
      not_in = "#ADB5BD",
      blue1 = "#1B4965",
      blue2 = "#5FA8D3",
      red1 = "#E76F51",
      green1 = "#2D6A4F"
    )

    heatmap_theme <- theme(
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 1),
      legend.position = "bottom",
      legend.key.width = unit(2.5, "cm")
    )

    # ── Dynamic-height helper for Run-faceted plots ─────────────────────────
    qc_plot_h <- reactive({
      d <- data()
      req(d, "Run" %in% names(d))
      n <- length(unique(d$Run))

      grid_rows <- ceiling(sqrt(n))
      max(500L, grid_rows * 180L)
    })

    make_qc_plot_ui <- function(plot_id) {
      renderUI({
        req(data())
        plotOutput(ns(plot_id), height = paste0(qc_plot_h(), "px"))
      })
    }

    output$dynamic_qc_plot_ui <- renderUI({
      req(data())
      plotOutput(ns("dynamic_qc_plot_out"), height = paste0(qc_plot_h(), "px"))
    })

    output$plot_aa_freq_ui <- make_qc_plot_ui("plot_aa_freq")
    output$plot_pep_type_count_ui <- make_qc_plot_ui("plot_pep_type_count")
    output$plot_pep_type_prop_ui <- make_qc_plot_ui("plot_pep_type_prop")

    # ════════════════════════════════════════════════════════════════════════
    # A. DATA REACTIVES
    # ════════════════════════════════════════════════════════════════════════
    data <- reactive({
      req(input$report)
      data_parquet <- arrow::read_parquet(input$report$datapath) %>%
      data.table::as.data.table() # Dataset must be a data.table from data.table package
      
      # Filters
      data_parquet <- data_parquet[
        Lib.PG.Q.Value <= 0.01 &
        Lib.Q.Value <= 0.01 &
        PG.Q.Value <= 0.01 &
        # If input is null, use zero else use the user input
        PG.MaxLFQ.Quality >= ifelse(is.null(input$PG.MaxLFQ.Quality), 0, input$PG.MaxLFQ.Quality) &
        # If input is null, use zero else use the user input
        Empirical.Quality >= ifelse(is.null(input$Empirical.Quality), 0, input$Empirical.Quality)
      ]
      
      # TODO: input$cRAP was not implemented yet #######
      # # If the user select to remove the cRAP, remove it
      # if(input$cRAP) {
      #   data_parquet <- data_parquet[!grepl("cRAP", Protein.Ids)]
      # }
      
      # Mutate
      data_parquet <- data_parquet[, `:=` (
        File.Name = Run,
        peptide_length = nchar(Stripped.Sequence)
      )]
      
      # return(data_parquet) # is implicit, no need to declare it
    })

    proteins <- reactive({
      req(data())
      data() %>%
        dplyr::group_by(Run) %>%
        dplyr::summarise(n_proteins = n_distinct(Protein.Ids))
    })

    peptides_per_run <- reactive({
      req(data())
      data() %>%
        dplyr::group_by(Run) %>%
        dplyr::summarise(n_peptides = n_distinct(Stripped.Sequence))
    })

    unique_genes <- reactive({
      req(data())
      diann::diann_matrix(
        data(),
        id.header = "Protein.Ids",
        quantity.header = "Genes.MaxLFQ.Unique",
        proteotypic.only = FALSE,
        pg.q = 0.01
      )
    })

    observe({
      req(unique_genes())
      cn <- colnames(unique_genes())
      updateSelectInput(session, "xcol", choices = cn)
      updateSelectInput(session, "ycol", choices = cn)
    })

    combined_data <- reactive({
      req(unique_genes())
      raw <- unique_genes() %>%
        log2() %>%
        as.data.frame() %>%
        tidyr::gather("Sample", "Intensity") %>%
        dplyr::mutate(norm = "Raw matrix")
      mad <- unique_genes() %>%
        log2() %>%
        limma::normalizeBetweenArrays(method = "scale") %>%
        as.data.frame() %>%
        tidyr::gather("Sample", "Intensity") %>%
        dplyr::mutate(norm = "MAD normalised")
      dplyr::bind_rows(raw, mad) %>%
        dplyr::mutate(
          norm = factor(norm, levels = c("Raw matrix", "MAD normalised"))
        )
    })

    QuantUMS_scores <- reactive({
      req(input$report)
      data_parquet <- arrow::read_parquet(input$report$datapath) %>%
        data.table::as.data.table() # Dataset must be a data.table from data.table package
      
      # Filters
      data_parquet <- data_parquet[
        Lib.PG.Q.Value <= 0.01 &
        Lib.Q.Value <= 0.01 &
        PG.Q.Value <= 0.01 
      ]
      
      # TODO: input$cRAP was not implemented yet #######
      # # If the user select to remove the cRAP, remove it
      # if(input$cRAP) {
      #   data_parquet <- data_parquet[!grepl("cRAP", Protein.Ids)]
      # }
      
      # Select columns
      data_parquet[, .(
          Run,
          Precursor.Id,
          PG.MaxLFQ.Quality,
          Empirical.Quality,
          Quantity.Quality
      )] %>%
        tidyr::pivot_longer(
          -c(Run, Precursor.Id),
          names_to = "Filter",
          values_to = "Score"
        )
    })

    gene_quantity <- reactive({
      req(data())
      data() %>%
        dplyr::select(Run, PG.MaxLFQ, Genes.MaxLFQ, Genes.MaxLFQ.Unique) %>%
        tidyr::pivot_longer(
          -Run,
          names_to = "protein_metrics",
          values_to = "quantity"
        )
    })

    missing_vs_mean <- reactive({
      req(unique_genes())
      base_df <- function(mat, label) {
        mat %>%
          log2() %>%
          as.data.frame() %>%
          tidyr::gather("Sample", "Intensity") %>%
          dplyr::mutate(missing = is.na(Intensity)) %>%
          dplyr::group_by(Sample) %>%
          dplyr::summarise(
            missing = mean(missing) * 100,
            median_intensity = median(Intensity, na.rm = TRUE)
          ) %>%
          dplyr::mutate(norm = label)
      }
      dplyr::bind_rows(
        base_df(unique_genes(), "Non normalised"),
        base_df(
          limma::normalizeBetweenArrays(log2(unique_genes()), method = "scale"),
          "MAD normalised"
        )
      ) %>%
        dplyr::mutate(
          norm = factor(norm, levels = c("Non normalised", "MAD normalised"))
        )
    })

    MS_corr <- reactive({
      req(input$report)
      data_parquet <- arrow::read_parquet(input$report$datapath) %>%
        data.table::as.data.table() # Dataset must be a data.table from data.table package
      
      # Filters
      data_parquet <- data_parquet[
        Lib.PG.Q.Value <= 0.01 &
        Lib.Q.Value <= 0.01 &
        PG.Q.Value <= 0.01 
      ]    
    
      # Mutate
      data_parquet <- data_parquet[, File.Name := Run]
    })

    PCA_label <- reactive({
      unique_genes() %>%
        log2() %>%
        na.omit() %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Sample")
    })

    pca_data <- reactive({
      unique_genes() %>%
        log2() %>%
        na.omit() %>%
        t() %>%
        prcomp(scale. = TRUE) %>%
        ggplot2::autoplot(data = PCA_label(), colour = "Sample", label = TRUE)
    })

    fasta_data <- reactive({
      req(input$fasta_file)
      seqinr::read.fasta(
        input$fasta_file$datapath,
        seqtype = "AA",
        as.string = TRUE
      )
    })

    digest_fasta <- reactive({
      req(fasta_data())
      mc <- max(0, min(3, as.integer(input$missed_cleavages)))
      seqs <- fasta_data()
      purrr::imap_dfr(seqs, function(seq_obj, prot_id) {
        first_word <- strsplit(prot_id, "\\s+")[[1]][1]
        parts <- strsplit(first_word, "\\|")[[1]]
        clean_id <- if (length(parts) >= 2) parts[2] else parts[1]
        aa <- toupper(as.character(seq_obj))
        pieces <- unlist(strsplit(aa, "(?<=[KR])(?!P)", perl = TRUE))
        pieces <- pieces[nchar(pieces) > 0]
        if (mc == 0) {
          peps <- pieces
        } else {
          all_peps <- character()
          for (i in seq_along(pieces)) {
            for (j in i:min(i + mc, length(pieces))) {
              all_peps <- c(all_peps, paste(pieces[i:j], collapse = ""))
            }
          }
          peps <- all_peps
        }
        peps <- peps[nchar(peps) >= 6 & nchar(peps) <= 52]
        if (length(peps) == 0) {
          return(NULL)
        }
        tibble::tibble(protein_id = clean_id, peptide = peps)
      })
    })

    peptide_type <- reactive({
      req(digest_fasta())
      digest_fasta() %>%
        dplyr::group_by(peptide) %>%
        dplyr::summarise(
          n_proteins = n_distinct(protein_id),
          proteotypic = n_proteins == 1,
          .groups = "drop"
        )
    })

    peptide_mapping <- reactive({
      req(data(), peptide_type())
      data() %>%
        dplyr::select(Run, Stripped.Sequence) %>%
        dplyr::distinct() %>%
        dplyr::left_join(
          peptide_type(),
          by = c("Stripped.Sequence" = "peptide")
        ) %>%
        dplyr::mutate(
          type = dplyr::case_when(
            is.na(proteotypic) ~ "Not in FASTA",
            proteotypic ~ "Proteotypic",
            TRUE ~ "Shared"
          ),
          type = factor(
            type,
            levels = c("Proteotypic", "Shared", "Not in FASTA")
          )
        )
    })

    peptide_summary <- reactive({
      req(peptide_mapping())
      peptide_mapping() %>%
        dplyr::group_by(Run, type) %>%
        dplyr::summarise(n = n(), .groups = "drop")
    })

    # ── Cleavage windows helpers ───────────────────────────────
    extract_one_window <- function(pep, prot_id, full_seqs) {
      if (
        is.null(prot_id) ||
          length(prot_id) == 0 ||
          is.na(prot_id) ||
          nchar(prot_id) == 0
      ) {
        return(rep(NA_character_, 8))
      }
      seq <- full_seqs[[prot_id]]
      if (is.null(seq) || length(seq) == 0 || is.na(seq[1])) {
        return(rep(NA_character_, 8))
      }
      seq <- seq[1]
      pos <- regexpr(pep, seq, fixed = TRUE)
      if (pos == -1L) {
        return(rep(NA_character_, 8))
      }
      start <- as.integer(pos)
      end <- start + nchar(pep) - 1L
      p_side <- substr(seq, max(1L, start - 4L), start - 1L)
      pp_side <- substr(seq, end + 1L, min(nchar(seq), end + 4L))
      p_side <- paste0(strrep("-", 4L - nchar(p_side)), p_side)
      pp_side <- paste0(pp_side, strrep("-", 4L - nchar(pp_side)))
      c(
        substr(p_side, 1, 1),
        substr(p_side, 2, 2),
        substr(p_side, 3, 3),
        substr(p_side, 4, 4),
        substr(pp_side, 1, 1),
        substr(pp_side, 2, 2),
        substr(pp_side, 3, 3),
        substr(pp_side, 4, 4)
      )
    }

    cleavage_windows <- reactive({
      req(data(), fasta_data())
      full_seqs <- sapply(fasta_data(), function(x) toupper(as.character(x)[1]))
      names(full_seqs) <- sapply(names(full_seqs), function(n) {
        fw <- trimws(strsplit(n, "\\s+")[[1]][1])
        p <- strsplit(fw, "\\|")[[1]]
        if (length(p) >= 2) p[2] else p[1]
      })
      identified <- data() %>%
        dplyr::select(Run, Stripped.Sequence, Protein.Ids) %>%
        dplyr::distinct() %>%
        dplyr::mutate(
          prot_id = sapply(Protein.Ids, function(ids) {
            if (is.na(ids) || nchar(ids) == 0) {
              return(NA_character_)
            }
            first <- trimws(strsplit(ids, ";")[[1]][1])
            parts <- strsplit(first, "\\|")[[1]]
            if (length(parts) >= 2) parts[2] else parts[1]
          })
        )
      windows <- mapply(
        extract_one_window,
        pep = identified$Stripped.Sequence,
        prot_id = identified$prot_id,
        MoreArgs = list(full_seqs = full_seqs),
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
      )
      if (!is.matrix(windows)) {
        windows <- matrix(windows, nrow = 8)
      }
      identified %>%
        dplyr::mutate(
          P4 = windows[1, ],
          P3 = windows[2, ],
          P2 = windows[3, ],
          P1 = windows[4, ],
          P1p = windows[5, ],
          P2p = windows[6, ],
          P3p = windows[7, ],
          P4p = windows[8, ]
        ) %>%
        dplyr::filter(!is.na(P1) & P1 != "-" & P1 != "") %>%
        dplyr::select(
          Run,
          peptide = Stripped.Sequence,
          P4,
          P3,
          P2,
          P1,
          P1p,
          P2p,
          P3p,
          P4p
        )
    })

    build_logo_matrix <- function(cw_subset) {
      pos_order <- c("P4", "P3", "P2", "P1", "P1p", "P2p", "P3p", "P4p")
      aa_order <- c(
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
      long <- cw_subset %>%
        tidyr::pivot_longer(
          cols = tidyr::all_of(pos_order),
          names_to = "position",
          values_to = "aa"
        ) %>%
        dplyr::filter(aa %in% aa_order) %>%
        dplyr::mutate(position = factor(position, levels = pos_order))
      mat <- matrix(
        0,
        nrow = length(aa_order),
        ncol = length(pos_order),
        dimnames = list(aa_order, pos_order)
      )
      for (pos in pos_order) {
        counts <- table(factor(
          long$aa[long$position == pos],
          levels = aa_order
        ))
        mat[, pos] <- as.numeric(counts)
      }
      mat
    }

    # ════════════════════════════════════════════════════════════════════════
    # B. INFO + BANNERS
    # ════════════════════════════════════════════════════════════════════════
    output$info_box1 <- renderInfoBox({
      infoBox(
        title = "Active Filters",
        value = HTML(paste0(
          "Lib.PG.Q.Value ≤ 0.01 &bull; Lib.Q.Value ≤ 0.01 &bull; PG.Q.Value ≤ 0.01<br>",
          "PG MaxLFQ Quality ≥ ",
          input$PG.MaxLFQ.Quality,
          " &bull; Empirical Quality ≥ ",
          input$Empirical.Quality
        )),
        icon = icon("filter"),
        color = "blue",
        fill = TRUE
      )
    })

    output$pep_map_banner <- renderUI({
      if (is.null(input$fasta_file)) {
        div(
          class = "plot-placeholder",
          icon("circle-info", style = "margin-right:10px;"),
          "Upload a Protein FASTA file in the sidebar to enable Peptide Mapping analysis."
        )
      }
    })

    output$prot_spec_banner <- renderUI({
      if (is.null(input$fasta_file)) {
        div(
          class = "plot-placeholder",
          icon("circle-info", style = "margin-right:10px;"),
          "Upload a Protein FASTA file in the sidebar to enable Protease Specificity analysis."
        )
      }
    })

    # ════════════════════════════════════════════════════════════════════════
    # C. QC PLOT REACTIVES (plot1–15)  — helpers on top
    # ════════════════════════════════════════════════════════════════════════
    plot1_obj <- reactive({
      req(data())
      data() %>%
        as.data.frame() %>%
        ggplot(aes(x = RT, y = Precursor.Quantity)) +
        geom_line(color = pal["blue1"], alpha = 0.7) +
        labs(x = "Retention time (min)", y = "Precursor quantity") +
        facet_wrap(~Run, scales = "free_y")
    })

    plot2_obj <- reactive({
      req(data())
      data() %>%
        as.data.frame() %>%
        ggplot(aes(x = RT, y = Precursor.Mz)) +
        ggpointdensity::geom_pointdensity(size = 0.25) +
        viridis::scale_color_viridis(option = "plasma") +
        labs(x = "Retention time (min)", y = "Scan range (m/z)", color = NULL) +
        theme(
          legend.position = "bottom",
          legend.key.width = unit(2, "cm"),
          legend.key.height = unit(0.3, "cm")
        ) +
        facet_wrap(~Run)
    })

    plot3_obj <- reactive({
      req(data())
      data() %>%
        as.data.frame() %>%
        ggplot(aes(x = Precursor.Mz, y = RT - Predicted.RT)) +
        ggpointdensity::geom_pointdensity(size = 0.25) +
        viridis::scale_color_viridis(option = "plasma") +
        geom_hline(
          yintercept = c(-1, 0, 1),
          linetype = "dashed",
          color = "black"
        ) +
        labs(x = "Precursor m/z", y = "RT − Predicted RT (min)", color = NULL) +
        theme(
          legend.position = "bottom",
          legend.key.width = unit(2, "cm"),
          legend.key.height = unit(0.3, "cm")
        ) +
        facet_wrap(~Run)
    })

    plot4_obj <- reactive({
      req(data())
      data() %>%
        as.data.frame() %>%
        ggplot(aes(x = Precursor.Charge)) +
        geom_density(fill = pal["blue1"], alpha = 0.7) +
        labs(x = "Precursor charge", y = "Density") +
        facet_wrap(~Run)
    })

    plot5_obj <- reactive({
      req(data())
      data() %>%
        ggplot(aes(x = peptide_length)) +
        geom_histogram(fill = pal["red1"], alpha = 0.8, binwidth = 1) +
        labs(x = "Peptide length (a.a.)", y = "Count") +
        facet_wrap(~Run, scales = "free_y")
    })

    plot6_obj <- reactive({
      req(peptides_per_run())
      peptides_per_run() %>%
        as.data.frame() %>%
        ggplot(aes(y = reorder(Run, n_peptides), x = n_peptides)) +
        geom_bar(stat = "identity", fill = pal["blue1"], alpha = 0.85) +
        geom_text(
          aes(label = scales::comma(n_peptides)),
          color = "white",
          size = 5,
          hjust = 1.08
        ) +
        labs(y = NULL, x = "Number of peptides")
    })

    plot7_obj <- reactive({
      req(proteins())
      proteins() %>%
        as.data.frame() %>%
        ggplot(aes(y = reorder(Run, n_proteins), x = n_proteins)) +
        geom_bar(stat = "identity", fill = pal["green1"], alpha = 0.85) +
        geom_text(
          aes(label = scales::comma(n_proteins)),
          color = "white",
          size = 5,
          hjust = 1.08
        ) +
        labs(y = NULL, x = "Number of proteins")
    })

    plot8_obj <- reactive({
      req(unique_genes())
      unique_genes() %>%
        as.data.frame() %>%
        tidyr::gather("Sample", "Intensity") %>%
        dplyr::mutate(missing = is.na(Intensity)) %>%
        dplyr::group_by(Sample) %>%
        dplyr::summarise(missing = mean(missing) * 100) %>%
        ggplot(aes(x = reorder(Sample, -missing), y = missing)) +
        geom_col(fill = pal["blue1"], alpha = 0.85) +
        geom_text(aes(label = round(missing, 1)), vjust = -0.4, size = 4.5) +
        labs(x = NULL, y = "Missing values (%)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })

    plot9_obj <- reactive({
      req(missing_vs_mean())
      missing_vs_mean() %>%
        ggplot(aes(x = missing, y = median_intensity)) +
        geom_point(size = 3, alpha = 0.8, color = pal["blue1"]) +
        geom_smooth(method = "lm", se = FALSE, color = pal["red1"]) +
        labs(x = "Missing values (%)", y = "Median log₂(intensity)") +
        facet_wrap(~norm)
    })

    plot10_obj <- reactive({
      req(combined_data())
      combined_data() %>%
        as.data.frame() %>%
        ggplot(aes(x = Sample, y = Intensity, fill = norm)) +
        scale_fill_manual(
          values = c(
            "Raw matrix" = pal[["red1"]],
            "MAD normalised" = pal[["blue1"]]
          )
        ) +
        geom_boxplot(alpha = 0.75) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        ) +
        labs(x = NULL, y = "log₂(Intensity)") +
        facet_wrap(~norm)
    })

    plot11_obj <- reactive({
      req(data())
      data() %>%
        as.data.frame() %>%
        dplyr::mutate(
          specificity = dplyr::case_when(
            stringr::str_detect(
              Stripped.Sequence,
              "K$|R$"
            ) ~ "Trypsin C-termini",
            stringr::str_detect(Stripped.Sequence, "E$|D$") ~ "GluC C-termini",
            TRUE ~ "Other C-termini"
          )
        ) %>%
        dplyr::group_by(Run, specificity) %>%
        dplyr::summarise(peptides = n(), .groups = "drop") %>%
        ggplot(aes(y = Run, x = peptides, fill = specificity)) +
        geom_col(alpha = 0.85, position = "stack") +
        scale_fill_manual(
          values = c(
            "Trypsin C-termini" = pal[["red1"]],
            "GluC C-termini" = pal[["blue1"]],
            "Other C-termini" = pal[["green1"]]
          )
        ) +
        labs(y = NULL, x = "Count", fill = "Specificity") +
        theme(legend.position = "top")
    })

    plot12_obj <- reactive({
      req(MS_corr(), input$Empirical.Quality)
      MS_corr() %>%
        as.data.frame() %>%
        dplyr::mutate(
          EQScore_cutoff = ifelse(
            Empirical.Quality >= input$Empirical.Quality,
            "Above threshold",
            "Below threshold"
          )
        ) %>%
        ggplot(aes(x = Ms1.Profile.Corr, fill = EQScore_cutoff)) +
        geom_density(alpha = 0.7) +
        scale_fill_manual(
          values = c(
            "Above threshold" = pal[["blue1"]],
            "Below threshold" = pal[["red1"]]
          )
        ) +
        labs(
          x = "MS1 Profile Correlation",
          y = "Density",
          fill = "Empirical Quality"
        ) +
        theme(legend.position = "top") +
        facet_wrap(~Run)
    })

    plot13_obj <- reactive({
      req(QuantUMS_scores())
      QuantUMS_scores() %>%
        as.data.frame() %>%
        ggplot(aes(x = Score, fill = Filter)) +
        geom_density(alpha = 0.7) +
        scale_fill_manual(
          values = c(
            "PG.MaxLFQ.Quality" = pal[["red1"]],
            "Empirical.Quality" = pal[["blue1"]],
            "Quantity.Quality" = pal[["green1"]]
          )
        ) +
        labs(x = "Score", y = "Density", fill = NULL) +
        theme(legend.position = "bottom") +
        facet_wrap(~Run)
    })

    plot14_obj <- reactive({
      req(gene_quantity())
      gene_quantity() %>%
        as.data.frame() %>%
        ggplot(aes(x = log2(quantity), fill = protein_metrics)) +
        geom_density(alpha = 0.7) +
        scale_fill_manual(
          values = c(
            "Genes.MaxLFQ" = pal[["red1"]],
            "Genes.MaxLFQ.Unique" = pal[["blue1"]],
            "PG.MaxLFQ" = pal[["green1"]]
          )
        ) +
        labs(x = "log₂(Protein Quantity)", y = "Density", fill = NULL) +
        theme(legend.position = "bottom") +
        facet_wrap(~Run)
    })

    plot15_obj <- reactive({
      req(data())
      data() %>%
        as.data.frame() %>%
        dplyr::select(Run, Stripped.Sequence) %>%
        dplyr::distinct() %>%
        dplyr::mutate(
          cysteine_count = stringr::str_count(Stripped.Sequence, "C"),
          cysteine_count_group = ifelse(
            cysteine_count >= 3,
            "3+",
            as.character(cysteine_count)
          )
        ) %>%
        dplyr::count(Run, cysteine_count_group) %>%
        ggplot(aes(x = cysteine_count_group, y = n)) +
        geom_bar(
          stat = "identity",
          position = "dodge",
          fill = pal["blue1"],
          color = "black"
        ) +
        geom_text(aes(label = n), vjust = -0.5, size = 4) +
        labs(x = "Cysteine counts in peptides", y = "Count") +
        facet_wrap(~Run, scales = "free_y")
    })

    # ── Dynamic QC plot dispatcher ──
    current_qc_plot_obj <- eventReactive(input$run_qc_plot, {
      req(input$qc_plot_select)
      show_spinner("sp_qc_main")
      switch(
        input$qc_plot_select,
        "plot1" = plot1_obj(),
        "plot2" = plot2_obj(),
        "plot3" = plot3_obj(),
        "plot4" = plot4_obj(),
        "plot5" = plot5_obj(),
        "plot6" = plot6_obj(),
        "plot7" = plot7_obj(),
        "plot8" = plot8_obj(),
        "plot9" = plot9_obj(),
        "plot10" = plot10_obj(),
        "plot11" = plot11_obj(),
        "plot12" = plot12_obj(),
        "plot13" = plot13_obj(),
        "plot14" = plot14_obj(),
        "plot15" = plot15_obj()
      )
    })

    output$dynamic_qc_plot_out <- renderPlot({
      on.exit(hide_spinner("sp_qc_main"), add = TRUE)
      req(current_qc_plot_obj())
      current_qc_plot_obj()
    })

    # ════════════════════════════════════════════════════════════════════════
    # D. AA FREQUENCY PLOT
    # ════════════════════════════════════════════════════════════════════════
    plot_aa_freq_obj <- reactive({
      req(input$fasta_file, data())
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
      seqs <- fasta_data()
      all_fasta_aa <- paste0(sapply(seqs, as.character), collapse = "")
      all_fasta_aa <- toupper(gsub("[^A-Za-z]", "", all_fasta_aa))
      aa_counts <- table(strsplit(all_fasta_aa, "")[[1]])
      aa_counts <- aa_counts[names(aa_counts) %in% base_aas]
      expected_freq <- as.data.frame(
        aa_counts / sum(aa_counts, na.rm = TRUE) * 100
      )
      colnames(expected_freq) <- c("AA", "Frequency")
      expected_freq$Type <- "Expected (FASTA)"

      obs_list <- lapply(unique(data()$Run), function(r) {
        p_data <- data()[data()$Run == r, ] %>%
          dplyr::select(Stripped.Sequence) %>%
          dplyr::distinct()
        pep_aa <- paste0(p_data$Stripped.Sequence, collapse = "")
        tc <- table(strsplit(pep_aa, "")[[1]])
        tc <- tc[names(tc) %in% base_aas]
        if (sum(tc, na.rm = TRUE) > 0) {
          f <- as.data.frame(tc / sum(tc, na.rm = TRUE) * 100)
          colnames(f) <- c("AA", "Frequency")
          f$Type <- "Identified (Peptides)"
          f$Run <- r
          f
        } else {
          NULL
        }
      })
      obs_df <- do.call(rbind, obs_list)

      exp_repl <- do.call(
        rbind,
        lapply(unique(obs_df$Run), function(r) {
          df <- expected_freq
          df$Run <- r
          df
        })
      )
      comb_df <- rbind(exp_repl, obs_df)
      comb_df$Type <- factor(
        comb_df$Type,
        levels = c("Expected (FASTA)", "Identified (Peptides)")
      )

      ggplot(comb_df, aes(x = AA, y = Frequency, fill = Type)) +
        geom_bar(
          stat = "identity",
          position = "dodge",
          color = "black",
          alpha = 0.85
        ) +
        scale_fill_manual(
          values = c(
            "Expected (FASTA)" = "#ADB5BD",
            "Identified (Peptides)" = pal[["blue1"]]
          )
        ) +
        labs(x = "Amino Acid", y = "Frequency (%)", fill = "") +
        theme(legend.position = "top") +
        facet_wrap(~Run)
    })

    output$plot_aa_freq <- renderPlot({
      on.exit(hide_spinner("sp_aa"), add = TRUE)
      plot_aa_freq_obj()
    })

    # ════════════════════════════════════════════════════════════════════════
    # E. PROTEIN COVERAGE
    # ════════════════════════════════════════════════════════════════════════
    observe({
      req(data(), input$fasta_file)
      identified_proteins <- unique(unlist(strsplit(data()$Protein.Ids, ";")))
      first_parts <- sapply(identified_proteins, function(x) {
        p <- strsplit(x, "\\|")[[1]]
        if (length(p) >= 2) p[2] else p[1]
      })
      valid_proteins <- first_parts[!is.na(first_parts)]
      if (length(valid_proteins) > 0) {
        updateSelectizeInput(
          session,
          "selected_protein",
          choices = c("", sort(valid_proteins)),
          server = TRUE
        )
      }
      runs <- sort(unique(data()$Run))
      updateSelectInput(
        session,
        "selected_sample_psv",
        choices = c("All Samples", runs)
      )
    })

    coverage_plot_data_qc4 <- reactive({
      req(input$selected_protein, input$fasta_file, data())
      if (input$selected_protein == "") {
        return(list(protein_sequence = NULL))
      }
      tp <- input$selected_protein
      tix <- grep(tp, names(fasta_data()), fixed = TRUE)[1]
      if (is.na(tix)) {
        return(list(error = paste("Protein", tp, "not found in FASTA.")))
      }
      pr_seq <- as.character(fasta_data()[[tix]])[1]

      pep_df <- data()[
        stringr::str_detect(data()$Protein.Ids, stringr::fixed(tp)),
      ]
      if (
        !is.null(input$selected_sample_psv) &&
          input$selected_sample_psv != "All Samples"
      ) {
        pep_df <- pep_df[pep_df$Run == input$selected_sample_psv, ]
      }
      if (nrow(pep_df) == 0) {
        return(list(
          error = paste("No peptides for", tp, "in selected sample(s).")
        ))
      }
      agg <- pep_df %>%
        dplyr::group_by(Run, Stripped.Sequence) %>%
        dplyr::summarise(n = n(), .groups = "drop")
      list(protein_sequence = pr_seq, peptides_data = agg)
    })

    output$protein_coverage_plot_ui <- renderUI({
      pd <- coverage_plot_data_qc4()
      if (!is.null(pd$error)) {
        return(div(
          style = "color:red;padding:15px;font-weight:bold;",
          pd$error
        ))
      }
      if (is.null(pd$protein_sequence)) {
        return(NULL)
      }
      apl <- max(
        10,
        min(200, ifelse(is.na(input$aa_per_line), 50, input$aa_per_line))
      )
      runs <- unique(pd$peptides_data$Run)
      lines <- max(1, ceiling(nchar(pd$protein_sequence) / apl))
      ph <- max(400, lines * (40 + 20 * length(runs)) + 100)
      plotOutput(ns("protein_coverage_plot"), height = paste0(ph, "px"))
    })

    output$protein_coverage_plot <- renderPlot({
      pd <- coverage_plot_data_qc4()
      req(is.null(pd$error), pd$protein_sequence)
      apl <- max(
        10,
        min(200, ifelse(is.na(input$aa_per_line), 50, input$aa_per_line))
      )
      pl <- nchar(pd$protein_sequence)
      runs <- sort(unique(pd$peptides_data$Run))
      nr_runs <- length(runs)

      run_masks <- lapply(runs, function(r) {
        mask <- rep(FALSE, pl)
        peps <- pd$peptides_data$Stripped.Sequence[pd$peptides_data$Run == r]
        for (pep in peps) {
          pos <- as.integer(regexpr(pep, pd$protein_sequence, fixed = TRUE))
          if (pos > 0) mask[pos:(pos + nchar(pep) - 1)] <- TRUE
        }
        mask
      })
      names(run_masks) <- runs
      nl <- max(1, ceiling(pl / apl))
      plots <- list()
      run_colors <- scales::hue_pal()(max(1, nr_runs))
      names(run_colors) <- runs

      for (li in seq_len(nl)) {
        sp <- (li - 1) * apl + 1
        ep <- min(sp + apl - 1, pl)
        aas <- strsplit(substr(pd$protein_sequence, sp, ep), "")[[1]]
        ad <- data.frame(
          position = sp:ep,
          aa = aas,
          x_pos = seq_along(aas),
          stringsAsFactors = FALSE
        )
        cov_df <- do.call(
          rbind,
          lapply(seq_along(runs), function(i) {
            r <- runs[i]
            mask <- run_masks[[r]][sp:ep]
            d <- data.frame(
              x_pos = seq_along(aas),
              has_peptide = mask,
              Run = r,
              y_pos = -0.5 - ((i - 1) * 0.9)
            )
            d[d$has_peptide, ]
          })
        )

        p <- ggplot() +
          xlim(0.5, apl + 0.5) +
          ylim(-0.5 - (nr_runs * 0.9), 1.5) +
          theme_void() +
          theme(
            plot.margin = margin(5, 5, 5, 5),
            legend.position = if (li == 1) "right" else "none"
          )
        p <- p +
          geom_tile(
            data = ad,
            aes(x = x_pos, y = 0.5),
            width = 0.9,
            height = 0.9,
            fill = "lightgray",
            color = "black",
            linewidth = 0.2,
            alpha = 0.3
          )
        if (isTRUE(input$show_sequence)) {
          p <- p +
            geom_text(
              data = ad,
              aes(x = x_pos, y = 0.5, label = aa),
              size = 4,
              fontface = "bold"
            )
        }
        pl2 <- ad[ad$position %% 10 == 0, ]
        if (nrow(pl2) > 0) {
          p <- p +
            geom_text(
              data = pl2,
              aes(x = x_pos, y = 1.2, label = position),
              size = 3,
              color = "black",
              fontface = "bold"
            )
        }
        if (!is.null(cov_df) && nrow(cov_df) > 0) {
          p <- p +
            geom_tile(
              data = cov_df,
              aes(x = x_pos, y = y_pos, fill = Run),
              width = 0.9,
              height = 0.8,
              color = "black",
              linewidth = 0.1
            ) +
            scale_fill_manual(values = run_colors, drop = FALSE)
        }
        plots[[li]] <- p +
          labs(fill = "Sample", title = paste0("AA ", sp, "-", ep)) +
          theme(plot.title = element_text(size = 10, face = "bold", hjust = 0))
      }
      tg <- grid::textGrob(
        paste("Protein sequence view (length:", pl, "aa)"),
        gp = grid::gpar(fontsize = 16, fontface = "bold")
      )
      combined <- do.call(gridExtra::arrangeGrob, c(plots, ncol = 1))
      gridExtra::grid.arrange(
        tg,
        combined,
        heights = c(0.5, max(4, nr_runs) * 1.5)
      )
    })

    # ════════════════════════════════════════════════════════════════════════
    # F. INTERACTIVE VIEWER
    # ════════════════════════════════════════════════════════════════════════
    cosine_obj <- reactive({
      req(unique_genes())
      unique_genes() %>%
        log2() %>%
        na.omit() %>%
        lsa::cosine() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Sample") %>%
        tidyr::pivot_longer(
          -Sample,
          names_to = "Match",
          values_to = "value"
        ) %>%
        ggplot(aes(x = Sample, y = Match, fill = value)) +
        geom_tile() +
        viridis::scale_fill_viridis(option = "E") +
        heatmap_theme +
        labs(x = NULL, y = NULL, fill = "Cosine similarity")
    })

    euclidean_obj <- reactive({
      req(unique_genes())
      unique_genes() %>%
        log2() %>%
        t() %>%
        dist(method = "euclidean") %>%
        as.matrix() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Sample") %>%
        tidyr::pivot_longer(
          -Sample,
          names_to = "Match",
          values_to = "value"
        ) %>%
        ggplot(aes(x = Sample, y = Match, fill = value)) +
        geom_tile() +
        viridis::scale_fill_viridis(option = "E") +
        heatmap_theme +
        labs(x = NULL, y = NULL, fill = "Euclidean distance")
    })

    jaccard_obj <- reactive({
      req(unique_genes())
      mat <- unique_genes() %>%
        log2() %>%
        t() %>%
        vegan::vegdist(method = "jaccard", na.rm = TRUE) %>%
        as.matrix()
      as.data.frame(mat) %>%
        tibble::rownames_to_column("Sample") %>%
        tidyr::pivot_longer(
          -Sample,
          names_to = "Match",
          values_to = "value"
        ) %>%
        ggplot(aes(x = Sample, y = Match, fill = value)) +
        geom_tile() +
        viridis::scale_fill_viridis(option = "E") +
        heatmap_theme +
        labs(x = NULL, y = NULL, fill = "Jaccard similarity")
    })

    output$cosine_similarity <- renderPlot({
      cosine_obj()
    })
    output$euclidean_distance <- renderPlot({
      euclidean_obj()
    })
    output$jaccard_similarity <- renderPlot({
      jaccard_obj()
    })

    # ── Sample-correlation scatterplot matrix (ggpairs) ──
    plot_ggpairs_obj <- reactive({
      req(unique_genes())
      d_mat <- unique_genes() %>% log2()
      d_mat[is.na(d_mat)] <- 0
      if (ncol(d_mat) < 2) {
        return(NULL)
      }
      GGally::ggpairs(
        d_mat,
        lower = list(
          continuous = GGally::wrap(
            "points",
            alpha = 0.3,
            size = 0.8,
            color = "#1B4965"
          )
        ),
        diag = list(
          continuous = GGally::wrap(
            "densityDiag",
            fill = "#5FA8D3",
            alpha = 0.8
          )
        ),
        upper = list(
          continuous = GGally::wrap(
            "cor",
            size = 4,
            color = "black",
            stars = FALSE
          )
        )
      ) +
        theme_bw() +
        theme(strip.text = element_text(face = "bold", size = 10))
    })

    output$plot_ggpairs <- renderPlot({
      on.exit(hide_spinner("sp_hm"), add = TRUE)
      plot_ggpairs_obj()
    })

    # ── Exploratory Factor Analysis ──
    plot_efa_obj <- reactive({
      req(unique_genes(), input$efa_factors)
      df <- unique_genes() %>% log2() %>% na.omit() %>% as.data.frame()
      tryCatch(
        {
          fit <- lavaan::efa(data = df, nfactors = input$efa_factors)
          if ("efaList" %in% class(fit)) {
            est <- lavaan::standardizedSolution(fit[[1]])
            loadings <- est[est$op == "=~", ]
            ggplot(loadings, aes(x = lhs, y = est.std)) +
              geom_boxplot(outliers = FALSE, fill = "#5FA8D3", alpha = 0.3) +
              geom_jitter(width = 0.2, color = "#1B4965", alpha = 0.7) +
              geom_hline(
                yintercept = 0.9,
                linetype = "dashed",
                color = "#E76F51"
              ) +
              ggrepel::geom_text_repel(
                data = subset(loadings, abs(est.std) > 0.9),
                aes(label = rhs),
                size = 3,
                color = "red",
                max.overlaps = Inf,
                min.segment.length = 0
              ) +
              theme_bw() +
              labs(
                x = "Latent Factor",
                y = "Standardized Loading",
                title = paste0(
                  "Exploratory Factor Analysis (",
                  input$efa_factors,
                  " Factors)"
                )
              ) +
              theme(text = element_text(size = 14))
          } else {
            ggplot() +
              theme_void() +
              ggtitle("EFA did not return expected model")
          }
        },
        error = function(e) {
          ggplot() +
            theme_void() +
            annotate(
              "text",
              x = 0.5,
              y = 0.5,
              label = paste0("EFA failed:\n", e$message)
            )
        }
      )
    })

    output$plot_efa <- renderPlot({
      on.exit(hide_spinner("spefa"), add = TRUE)
      plot_efa_obj()
    })

    # ── Interactive scatter (Corr) ──
    output$Corr <- renderPlotly({
      req(unique_genes(), input$xcol, input$ycol)
      p <- unique_genes() %>%
        log2() %>%
        as.data.frame() %>%
        ggplot(aes(x = !!sym(input$xcol), y = !!sym(input$ycol))) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE, color = pal["red1"]) +
        labs(
          x = paste0("log₂(", input$xcol, ")"),
          y = paste0("log₂(", input$ycol, ")")
        )
      ggplotly(p)
    })

    # ── 3D QuantUMS score distribution ──
    output$QuantUMS_dist <- renderPlotly({
      req(data())
      data() %>%
        plot_ly(
          x = ~PG.MaxLFQ.Quality,
          y = ~Quantity.Quality,
          z = ~Empirical.Quality,
          color = ~Run,
          alpha = 0.6,
          colors = viridis::viridis(256),
          type = "scatter3d",
          mode = "markers"
        ) %>%
        layout(
          scene = list(
            xaxis = list(title = "PG MaxLFQ Quality"),
            yaxis = list(title = "Quantity Quality"),
            zaxis = list(title = "Empirical Quality")
          )
        )
    })

    # ── PCA ──
    output$PCA <- renderPlotly({
      on.exit(hide_spinner("sp_pca"), add = TRUE)
      req(pca_data())
      ggplotly(pca_data())
    })

    # ════════════════════════════════════════════════════════════════════════
    # G. PEPTIDE MAPPING PLOTS & TABLE
    # ════════════════════════════════════════════════════════════════════════
    plot_pep_count_obj <- reactive({
      req(peptide_summary())
      peptide_summary() %>%
        ggplot(aes(y = Run, x = n, fill = type)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.85) +
        geom_text(
          aes(label = n),
          position = position_stack(vjust = 0.5),
          size = 3,
          color = "white",
          fontface = "bold"
        ) +
        scale_fill_manual(
          values = c(
            "Proteotypic" = pal[["green1"]],
            "Shared" = pal[["red1"]],
            "Not in FASTA" = pal[["not_in"]]
          )
        ) +
        labs(y = NULL, x = "Number of peptides", fill = "Peptide type") +
        theme(legend.position = "top")
    })

    plot_pep_prop_obj <- reactive({
      req(peptide_summary())
      peptide_summary() %>%
        dplyr::group_by(Run) %>%
        dplyr::mutate(prop = n / sum(n) * 100) %>%
        ggplot(aes(y = Run, x = prop, fill = type)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.85) +
        geom_text(
          aes(label = sprintf("%.1f%%", prop)),
          position = position_stack(vjust = 0.5),
          size = 3,
          color = "white",
          fontface = "bold"
        ) +
        scale_fill_manual(
          values = c(
            "Proteotypic" = pal[["green1"]],
            "Shared" = pal[["red1"]],
            "Not in FASTA" = pal[["not_in"]]
          )
        ) +
        labs(y = NULL, x = "Proportion (%)", fill = "Peptide type") +
        theme(legend.position = "top")
    })

    output$plot_pep_type_count <- renderPlot({
      req(input$fasta_file)
      on.exit(hide_spinner("sp_pc"), add = TRUE)
      plot_pep_count_obj()
    })

    output$plot_pep_type_prop <- renderPlot({
      req(input$fasta_file)
      on.exit(hide_spinner("sp_pp"), add = TRUE)
      plot_pep_prop_obj()
    })

    output$peptide_table <- DT::renderDT(
      {
        req(input$fasta_file, peptide_mapping())
        peptide_mapping() %>%
          dplyr::group_by(Run, type) %>%
          dplyr::summarise(n_peptides = n(), .groups = "drop") %>%
          tidyr::pivot_wider(
            names_from = type,
            values_from = n_peptides,
            values_fill = 0
          ) %>%
          dplyr::rename(Sample = Run) %>%
          DT::datatable(
            rownames = FALSE,
            options = list(pageLength = 20, dom = "tp"),
            class = "cell-border stripe hover"
          )
      },
      server = FALSE
    )

    # ════════════════════════════════════════════════════════════════════════
    # H. PROTEASE SPECIFICITY (SeqLogos)
    # ════════════════════════════════════════════════════════════════════════
    seqlogo_x_scale <- scale_x_continuous(
      breaks = 1:8,
      labels = c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
    )

    output$seqlogo_all <- renderPlot({
      req(input$fasta_file)
      cw <- cleavage_windows()
      req(nrow(cw) > 0)
      mat <- build_logo_matrix(cw)
      ggseqlogo::ggseqlogo(mat, method = "prob", seq_type = "aa") +
        seqlogo_x_scale +
        labs(
          title = "All Runs — Schechter-Berger P4–P4' Cleavage Specificity",
          x = "Position relative to cleavage site",
          y = "Probability"
        ) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    })

    output$seqlogo_perrun_ui <- renderUI({
      req(input$fasta_file)
      cw <- cleavage_windows()
      runs <- unique(cw$Run)
      tagList(lapply(seq_along(runs), function(i) {
        fluidRow(box(
          title = paste("Run:", runs[i]),
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          plotOutput(ns(paste0("seqlogo_run_", i)), height = "350px")
        ))
      }))
    })

    observe({
      req(input$fasta_file)
      cw <- cleavage_windows()
      runs <- unique(cw$Run)
      for (i in seq_along(runs)) {
        local({
          rn <- runs[i]
          oid <- paste0("seqlogo_run_", i)
          output[[oid]] <- renderPlot({
            sub_cw <- dplyr::filter(cw, Run == rn)
            if (nrow(sub_cw) < 10) {
              return(NULL)
            }
            ggseqlogo::ggseqlogo(
              build_logo_matrix(sub_cw),
              method = "prob",
              seq_type = "aa"
            ) +
              seqlogo_x_scale +
              labs(x = "Position (P4–P4')", y = "Probability") +
              theme(plot.title = element_text(hjust = 0.5, face = "bold"))
          })
        })
      }
    })

    # ════════════════════════════════════════════════════════════════════════
    # I. DOWNLOADS
    # ════════════════════════════════════════════════════════════════════════
    output$download_matrix <- downloadHandler(
      filename = function() {
        paste0("filtered_protein_matrix_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        readr::write_tsv(
          as.data.frame(unique_genes()) %>%
            tibble::rownames_to_column("protein_id"),
          file
        )
      }
    )

    output$download_qc_plot <- downloadHandler(
      filename = function() {
        paste0("QC4DIANN_", input$qc_plot_select, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        obj <- current_qc_plot_obj()
        req(obj)
        if (inherits(obj, "Heatmap") || inherits(obj, "HeatmapList")) {
          png(file, width = 11, height = 8, units = "in", res = 300)
          ComplexHeatmap::draw(obj)
          dev.off()
        } else {
          ggsave(
            file,
            obj,
            width = 11,
            height = 8,
            bg = "white",
            device = "png",
            dpi = 300
          )
        }
      }
    )
  })
}
