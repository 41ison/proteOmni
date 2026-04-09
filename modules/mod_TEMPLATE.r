# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  proteOmni — Custom Search Engine Module Template                          ║
# ║                                                                            ║
# ║  PURPOSE: A fully documented, copy-and-modify template for adding new      ║
# ║  search engine analysis modules to the proteOmni platform.                 ║
# ║                                                                            ║
# ║  HOW TO USE:                                                               ║
# ║    1. Copy this file and rename it (e.g., mod_MyEngine.r)                  ║
# ║    2. Replace all occurrences of "Template" with your module name           ║
# ║    3. Implement your data processing functions                              ║
# ║    4. Build your UI (sidebar + body) and server logic                       ║
# ║    5. Register the module in proteOmni.r (see Integration section below)   ║
# ║                                                                            ║
# ║  CONVENTIONS:                                                              ║
# ║    • UI functions:  <ModuleName>_sidebar_ui(id)                            ║
# ║                     <ModuleName>_body_ui(id)                               ║
# ║    • Server:        <ModuleName>_server(id)                                ║
# ║    • All inputs/outputs must use ns() for proper namespacing               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


# ── 1. HELPER FUNCTIONS ─────────────────────────────────────────────────────
#
# Define pure, testable functions outside the module. This keeps your server
# logic clean and makes unit testing possible. These functions should NOT
# depend on any Shiny reactive context.
# ─────────────────────────────────────────────────────────────────────────────

#' Filter a proteomics data frame by a minimum score threshold.
#'
#' @param df    A data.frame with at least a "Score" column.
#' @param score Numeric threshold — rows with Score >= this value are kept.
#' @return Filtered data.frame.
filter_by_score <- function(df, score = 0.01) {
  stopifnot("Score" %in% colnames(df))
  dplyr::filter(df, Score >= score)
}


#' Count the number of unique peptides per protein.
#'
#' @param df A data.frame with at least "Protein" and "Peptide" columns.
#' @return A tibble with columns: Protein, n_peptides.
count_peptides_per_protein <- function(df) {
  stopifnot(all(c("Protein", "Peptide") %in% colnames(df)))
  df %>%
    dplyr::distinct(Protein, Peptide) %>%
    dplyr::count(Protein, name = "n_peptides") %>%
    dplyr::arrange(dplyr::desc(n_peptides))
}


#' Compute basic summary statistics per protein.
#'
#' @param df A data.frame with "Protein" and "Intensity" columns.
#' @return A tibble with columns: Protein, mean_intensity, sd_intensity, n.
summarise_protein_intensity <- function(df) {
  stopifnot(all(c("Protein", "Intensity") %in% colnames(df)))
  df %>%
    dplyr::group_by(Protein) %>%
    dplyr::summarise(
      mean_intensity = mean(Intensity, na.rm = TRUE),
      sd_intensity   = sd(Intensity, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
}


# ── 2. SIDEBAR UI ───────────────────────────────────────────────────────────
#
# The sidebar holds all user controls: file upload, analysis parameters,
# action buttons, and download handlers.
#
# IMPORTANT: Every input/output ID MUST be wrapped in ns() to ensure
# proper namespacing when multiple modules coexist in proteOmni.
# ─────────────────────────────────────────────────────────────────────────────

Template_sidebar_ui <- function(id) {
  # Create the module's namespace function

  ns <- NS(id)

  tagList(
    tags$div(
      id = ns("sidebar_content"),

      # ── Section: Data Upload ──────────────────────────────────────────
      tags$div(
        style = "padding:12px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "Data Upload"
      ),
      fileInput(
        ns("data_file"),
        "Upload Search Results (.tsv/.csv)",
        accept = c(".tsv", ".csv", ".txt")
      ),
      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── Section: Analysis Parameters ──────────────────────────────────
      tags$div(
        style = "padding:12px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "Analysis Parameters"
      ),
      numericInput(
        ns("score_threshold"),
        "Minimum Score Threshold",
        value = 0.01,
        min = 0, max = 1, step = 0.01
      ),
      selectInput(
        ns("color_palette"),
        "Plot Color Palette",
        choices = c("Dark2", "Set1", "Set2", "Paired", "Spectral"),
        selected = "Dark2"
      ),
      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── Section: Run & Download ───────────────────────────────────────
      tags$div(
        style = "padding:0 8px;text-align:center;",
        actionButton(
          ns("run_analysis"),
          "Run Analysis",
          class = "btn-primary",
          style = "width:80%;font-weight:bold;margin-top:10px;margin-bottom:10px;"
        )
      ),
      tags$div(
        style = "padding:0 8px;",
        div(
          style = "margin-bottom:8px;",
          downloadButton(
            ns("download_results"),
            "\u2B07 Download Results",
            class = "dl-btn",
            style = "width:100%;text-align:left;"
          )
        )
      )
    )
  )
}


# ── 3. BODY UI ──────────────────────────────────────────────────────────────
#
# The body holds the main content: tabs with data tables and plots.
# Use shinydashboard::box() for consistent card styling with the rest
# of proteOmni. Wrap every plot in a spinner-overlay div for loading
# feedback.
# ─────────────────────────────────────────────────────────────────────────────

Template_body_ui <- function(id) {
  ns <- NS(id)

  tabsetPanel(
    id = ns("tabs"),
    type = "tabs",

    # ── Tab 1: Data Preview ─────────────────────────────────────────────
    tabPanel(
      "Data Preview",
      fluidRow(
        box(
          title = "Uploaded Data",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          DT::dataTableOutput(ns("data_preview"))
        )
      )
    ),

    # ── Tab 2: Intensity Distribution ───────────────────────────────────
    tabPanel(
      "Intensity QC",
      fluidRow(
        box(
          title = "Intensity distributions per sample",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_intensity"),
              icon("spinner", class = "fa-spin")
            ),
            plotOutput(ns("intensity_boxplot"), height = 600)
          )
        )
      )
    ),

    # ── Tab 3: Identification Counts ────────────────────────────────────
    tabPanel(
      "Identifications",
      fluidRow(
        box(
          title = "Top proteins by unique peptide count",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_counts"),
              icon("spinner", class = "fa-spin")
            ),
            plotOutput(ns("protein_bar_chart"), height = 600)
          )
        )
      )
    )
  )
}


# ── 4. SERVER FUNCTION ──────────────────────────────────────────────────────
#
# The moduleServer wraps all reactive logic. Key patterns:
#   • reactive()      — for lazy, cached computations
#   • eventReactive() — for expensive computations triggered by a button
#   • observeEvent()  — for side effects (e.g., updating UI elements)
#   • renderPlot()    — for plot outputs
#
# NAMESPACING: Inside moduleServer, input$X and output$Y are automatically
# namespaced — you do NOT need ns() here. Only use ns() in renderUI()
# calls that generate new dynamic UI elements.
# ─────────────────────────────────────────────────────────────────────────────

Template_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Spinner Helpers (consistent with proteOmni convention) ────────────
    spin_ids <- c("sp_intensity", "sp_counts")
    show_spinners <- function() {
      lapply(spin_ids, function(s) shinyjs::show(id = s))
    }
    hide_spinner <- function(sid) shinyjs::hide(id = sid)

    # rh(): Render-with-Hide. Wraps a render function so the spinner
    # is automatically hidden when the plot finishes rendering.
    rh <- function(expr_fn, sid) {
      on.exit(hide_spinner(sid), add = TRUE)
      expr_fn()
    }

    # ── Data Ingestion ───────────────────────────────────────────────────
    # A reactive that reads the uploaded file and returns a data.frame.
    # This fires automatically whenever the user uploads a new file.
    raw_data <- reactive({
      req(input$data_file)
      show_spinners()

      ext <- tools::file_ext(input$data_file$name)
      if (ext %in% c("tsv", "txt")) {
        readr::read_tsv(input$data_file$datapath, show_col_types = FALSE)
      } else {
        readr::read_csv(input$data_file$datapath, show_col_types = FALSE)
      }
    })

    # ── Filtered Data ────────────────────────────────────────────────────
    # An eventReactive that only recalculates when the user clicks
    # "Run Analysis". This prevents expensive re-computation on every
    # parameter change.
    filtered_data <- eventReactive(input$run_analysis, {
      req(raw_data())

      withProgress(message = "Processing data...", value = 0, {
        incProgress(0.3, detail = "Applying score filter")
        result <- filter_by_score(raw_data(), input$score_threshold)

        incProgress(0.7, detail = "Calculating protein statistics")
        # You can add more processing steps here

        incProgress(1.0, detail = "Done!")
        result
      })
    })

    # ── Data Preview Table ───────────────────────────────────────────────
    output$data_preview <- DT::renderDataTable({
      req(raw_data())
      DT::datatable(
        head(raw_data(), 200),
        rownames = FALSE,
        options = list(
          dom = "frtip",
          pageLength = 15,
          scrollX = TRUE
        ),
        class = "display compact"
      )
    })

    # ── Plot 1: Intensity Boxplot ────────────────────────────────────────
    # Visualises the distribution of log2(Intensity) across samples.
    output$intensity_boxplot <- renderPlot({
      rh(
        function() {
          req(filtered_data())
          df <- filtered_data()

          # Guard: ensure required columns exist
          req("Intensity" %in% colnames(df))

          # If there is a "Sample" column, use it for faceting
          if ("Sample" %in% colnames(df)) {
            plot_df <- df %>%
              dplyr::mutate(log2_intensity = log2(Intensity + 1))

            ggplot(plot_df, aes(x = Sample, y = log2_intensity, fill = Sample)) +
              geom_boxplot(outlier.alpha = 0.3) +
              scale_fill_brewer(palette = input$color_palette) +
              labs(
                x = NULL,
                y = "log<sub>2</sub>(Intensity)"
              ) +
              theme_bw() +
              theme(
                text = element_text(size = 16),
                axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
                axis.title = ggtext::element_markdown(face = "bold"),
                legend.position = "none"
              )
          } else {
            # Fallback: single histogram if no Sample column
            plot_df <- df %>%
              dplyr::mutate(log2_intensity = log2(Intensity + 1))

            ggplot(plot_df, aes(x = log2_intensity)) +
              geom_histogram(bins = 50, fill = "#2c3e50", alpha = 0.7) +
              labs(
                x = "log<sub>2</sub>(Intensity)",
                y = "Count"
              ) +
              theme_bw() +
              theme(
                text = element_text(size = 16),
                axis.title = ggtext::element_markdown(face = "bold")
              )
          }
        },
        "sp_intensity"
      )
    })

    # ── Plot 2: Protein Bar Chart ────────────────────────────────────────
    # Shows the top 20 proteins ranked by unique peptide count.
    output$protein_bar_chart <- renderPlot({
      rh(
        function() {
          req(filtered_data())
          df <- filtered_data()

          # Guard: ensure required columns exist
          req(all(c("Protein", "Peptide") %in% colnames(df)))

          prot_counts <- count_peptides_per_protein(df) %>%
            head(20)

          ggplot(
            prot_counts,
            aes(
              x = reorder(Protein, n_peptides),
              y = n_peptides,
              fill = n_peptides
            )
          ) +
            geom_col() +
            coord_flip() +
            scale_fill_viridis_c(option = "D") +
            labs(
              x = NULL,
              y = "Unique peptides",
              fill = NULL
            ) +
            theme_bw() +
            theme(
              text = element_text(size = 14),
              axis.title = element_text(face = "bold"),
              axis.text = element_text(color = "black"),
              legend.position = "none"
            )
        },
        "sp_counts"
      )
    })

    # ── Download Handler ─────────────────────────────────────────────────
    output$download_results <- downloadHandler(
      filename = function() {
        paste0("template_results_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        readr::write_tsv(filtered_data(), file)
      }
    )
  })
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  5. INTEGRATION INSTRUCTIONS                                               ║
# ║                                                                            ║
# ║  Follow these steps to register your module in the main proteOmni app.     ║
# ║  All changes go into: app/proteOmni.r                                      ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
#
# ── Step 1: Source this module ───────────────────────────────────────────────
#
#   In the "Source modules" section of proteOmni.r, add:
#
#     source("modules/mod_TEMPLATE.r") you can give it any name you want.
#
#
# ── Step 2: Add a landing-page card ──────────────────────────────────────────
#
#   In the body → "Home" tabItem, add a card alongside the existing ones:
#
#     div(
#       class = "tool-card",
#       onclick = "Shiny.setInputValue('tool_select', 'Template', {priority: 'event'})",
#       icon("flask", class = "tool-icon"),      # Choose an appropriate icon
#       h4("Template Engine"),
#       p("Description of what this module does.")
#     )
#
#
# ── Step 3: Add sidebar UI ──────────────────────────────────────────────────
#
#   In the server function, inside the observeEvent(input$tool_select, ...):
#
#     "Template" = {
#       output$dynamic_sidebar <- renderUI({
#         Template_sidebar_ui("template_mod")
#       })
#       updateTabItems(session, "main_tabs", selected = "template_tab")
#     }
#
#
# ── Step 4: Add body tabItem ─────────────────────────────────────────────────
#
#   In the dashboardBody, add a new tabItem:
#
#     tabItem(
#       tabName = "template_tab",
#       Template_body_ui("template_mod")
#     )
#
#
# ── Step 5: Call the server ──────────────────────────────────────────────────
#
#   At the end of the server function:
#
#     Template_server("template_mod")
#
#
# ── Step 6 (Optional): Pass reactive data from another module ────────────────
#
#   If your module needs data from elsewhere (e.g., the PSManalyst matrix),
#   modify the server function signature to accept it:
#
#     Template_server <- function(id, shared_data) {
#       moduleServer(id, function(input, output, session) {
#         # Access the shared reactive:
#         observe({
#           df <- shared_data()
#           # ... use df
#         })
#       })
#     }
#
#   Then in proteOmni.r:
#
#     # Create a shared reactive (e.g., from PSManalyst)
#     psm_matrix <- PSManalyst_server("psm_mod")
#
#     # Pass it to your module
#     Template_server("template_mod", shared_data = psm_matrix)
#
#
# ══════════════════════════════════════════════════════════════════════════════
