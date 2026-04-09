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
    labs(x = "Retention time (min)", y = "*m/z*", color = "Density") +
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
    labs(x = "*m/z*", y = "Number of data points", color = "Density") +
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

      # ── Data Upload ───────────────────────────────────────────────────────
      tags$div(class = "sidebar-section-label", "Data Upload"),

      fileInput(
        ns("msms_file"),
        "Upload msms.txt",
        accept = c(".txt", ".tsv")
      ),

      fileInput(
        ns("evidence_file"),
        "Upload evidence.txt",
        accept = c(".txt", ".tsv")
      ),

      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── MS/MS Parameters ──────────────────────────────────────────────────
      tags$div(class = "sidebar-section-label", "MS/MS Parameters"),

      selectizeInput(
        ns("peptide_seq"),
        "Select Peptide Sequence",
        choices = NULL,
        options = list(
          placeholder = "Upload msms.txt first…",
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
          style = "width:100%;text-align:left;"
        )
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

    # ── Tab 1: MS/MS Spectrum ─────────────────────────────────────────────
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
    # A. msms.txt REACTIVES
    # ════════════════════════════════════════════════════════════════════

    raw_msms <- reactive({
      req(input$msms_file)
      withProgress(message = "Reading msms.txt…", value = 0.4, {
        df <- read_msms_file(input$msms_file$datapath)
        incProgress(0.6, detail = "Done.")
        df
      })
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
      req(raw_msms(), nchar(input$peptide_seq) > 0)
      show_spinner("sp_spectrum")
      withProgress(message = "Tidying MS/MS data…", value = 0.5, {
        result <- tidy_msms(raw_msms(), input$peptide_seq)
        incProgress(0.5, detail = "Done.")
        result
      })
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
      req(input$evidence_file)
      withProgress(message = "Reading evidence.txt…", value = 0.4, {
        df <- read_evidence_file(input$evidence_file$datapath)
        incProgress(0.6, detail = "Done.")
        df
      })
    })

    # The currently selected evidence plot — recalculate on button click
    current_ev_plot <- eventReactive(input$run_evidence, {
      req(raw_evidence())
      show_spinner("sp_evidence")

      df <- raw_evidence()
      sel <- input$ev_plot_select

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
    })

    # Dynamic height for evidence plots
    ev_plot_height_px <- reactive({
      req(raw_evidence())
      facet_height_px(raw_evidence())
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
      content = function(file) readr::write_tsv(tidy_msms_data(), file)
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
      content = function(file) readr::write_tsv(raw_evidence(), file)
    )
  })
}
