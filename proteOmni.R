## proteOmni — Unified Proteomics QC Dashboard
## Modules: PSManalyst (FragPipe - DDA QC), QC4DIANN (DIA-NN - DIA QC),
##          deNovo (Casanovo - de novo MS/MS), PwrQuant (limma + ORA),
##          MaxQuant (MaxQuant - DDA QC), InstaNovo (InstaNovo - de novo MS/MS),
##          EncyclopeDIA (EncyclopeDIA - DIA), Sage (Sage - DDA/DIA QC)

# ── Package installation ─────────────────────────────────────────────────────
options(repos = c(CRAN = "https://cran.rstudio.com/"))
CRAN_packages <- c(
  "shiny",
  "shinydashboard",
  "shinydashboardPlus",
  "shinyjs",
  "fresh",
  "devtools",
  "BiocManager",
  "tidyverse",
  "tidytext",
  "janitor",
  "ggpointdensity",
  "ggtext",
  "ggrepel",
  "ggseqlogo",
  "lsa",
  "vegan",
  "plotly",
  "viridis",
  "ggfortify",
  "seqinr",
  "zip",
  "DT",
  "colourpicker",
  "R6",
  "gridExtra",
  "scales",
  "lavaan",
  "rhandsontable",
  "naniar",
  "patchwork",
  "pwr",
  "missForest",
  "data.table",
  "GGally",
  "ape",
  "ggiraph",
  "ggforce",
  "ggridges"
)

not_inst <- CRAN_packages[
  !(CRAN_packages %in% installed.packages()[, "Package"])
]
if (length(not_inst)) {
  install.packages(not_inst)
}
if (!requireNamespace("diann", quietly = TRUE)) {
  devtools::install_github("https://github.com/vdemichev/diann-rpackage")
}
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", update = FALSE, ask = FALSE)
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings", update = FALSE, ask = FALSE)
}
if (!requireNamespace("sva", quietly = TRUE)) {
  BiocManager::install("sva", update = FALSE, ask = FALSE)
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
}
if (!requireNamespace("GO.db", quietly = TRUE)) {
  BiocManager::install("GO.db", update = FALSE, ask = FALSE)
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot", update = FALSE, ask = FALSE)
}
Sys.setenv(LIBARROW_MINIMAL = "false", ARROW_WITH_ZSTD = "ON")
if (!requireNamespace("arrow", quietly = TRUE)) {
  install.packages("arrow")
}

# ── Libraries ─────────────────────────────────────────────────────────────────
library(shiny)
library(shinydashboard)
library(shinyjs)
library(fresh)
library(diann)
library(arrow)
library(tidyverse)
library(tidytext)
library(janitor)
library(ggpointdensity)
library(ggtext)
library(ggseqlogo)
library(lsa)
library(vegan)
library(plotly)
library(viridis)
library(ggfortify)
library(seqinr)
library(zip)
library(DT)
library(data.table)
library(colourpicker)
library(R6)
library(gridExtra)
library(grid)
library(GGally)
library(Biostrings)
library(limma)
library(lavaan)
library(patchwork)
library(enrichplot)
library(ComplexHeatmap)

# ── Session info & Logging ─────────────────────────────────────────────────
cat("\n========== proteOmni Session Info ==========", "\n")
print(sessionInfo())
cat("============================================\n\n")

.omni_log_path <- file.path(getwd(), "Session_Info_log.txt")
.omni_log_con <- file(.omni_log_path, open = "wt")

writeLines(
  c(
    paste0("proteOmni — Session Log  [", Sys.time(), "]"),
    paste0("R version: ", R.version.string),
    paste0("Platform:  ", R.version$platform),
    paste0("Working directory: ", getwd()),
    "",
    "=== Loaded Package Versions ==="
  ),
  con = .omni_log_con
)

.omni_pkgs <- c(
  "shiny",
  "shinydashboard",
  "tidyverse",
  "ggplot2",
  "limma",
  "sva",
  "clusterProfiler",
  "enrichplot",
  "patchwork",
  "pwr",
  "missForest",
  "naniar",
  "DT",
  "ggrepel",
  "ggtext",
  "GGally",
  "Biostrings",
  "lavaan",
  "arrow",
  "diann"
)
for (.pkg in .omni_pkgs) {
  .ver <- tryCatch(as.character(packageVersion(.pkg)), error = function(e) {
    "not installed"
  })
  writeLines(paste0("  ", .pkg, ": ", .ver), con = .omni_log_con)
}
writeLines(c("", "=== Runtime Log ===", ""), con = .omni_log_con)
flush(.omni_log_con)

.omni_log_env <- new.env(parent = emptyenv())
.omni_log_env$buffer <- character(0)

.omni_buf_append <- function(line) {
  .omni_log_env$buffer <- c(.omni_log_env$buffer, line)
}

.omni_buf_append(paste0("proteOmni \u2014 Console History  [", Sys.time(), "]"))
.omni_buf_append(paste0("R version: ", R.version.string))
.omni_buf_append(paste0("Platform:  ", R.version$platform))
.omni_buf_append(paste0("Working directory: ", getwd()))
.omni_buf_append("")
.omni_buf_append("=== Loaded Package Versions ===")
for (.pkg in .omni_pkgs) {
  .ver <- tryCatch(as.character(packageVersion(.pkg)), error = function(e) {
    "not installed"
  })
  .omni_buf_append(paste0("  ", .pkg, ": ", .ver))
}
.omni_buf_append("")
.omni_buf_append("=== Console Log ===")
.omni_buf_append("")

.omni_log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  line <- paste0(timestamp, " ", conditionMessage(msg))
  tryCatch(
    {
      writeLines(line, con = .omni_log_con)
      flush(.omni_log_con)
    },
    error = function(e) NULL
  )
  .omni_buf_append(line)
  invokeRestart("muffleMessage")
}

.omni_log_warn <- function(w) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  line <- paste0(timestamp, " WARNING: ", conditionMessage(w))
  tryCatch(
    {
      writeLines(line, con = .omni_log_con)
      flush(.omni_log_con)
    },
    error = function(e) NULL
  )
  .omni_buf_append(line)
  invokeRestart("muffleWarning")
}

globalCallingHandlers(message = .omni_log_msg, warning = .omni_log_warn)
message("proteOmni started at ", Sys.time())


# ── Source modules ────────────────────────────────────────────────────────────
source("modules/utils_fasta.r")
source("modules/mod_PSManalyst.r")
source("modules/mod_QC4DIANN.r")
source("modules/mod_PwrQuant.r")
OMNIVIEW_LOADED <- TRUE
source("modules/dash_deNovo.r")
source("modules/mod_InstaNovo.r")
source("modules/mod_EncyclopeDIA.r")
source("modules/mod_Sage.r")
source("modules/mod_MaxQuantMSMS.r")


# ── Global options ────────────────────────────────────────────────────────────
options(shiny.maxRequestSize = 2000 * 1024^2)

theme_set(theme_bw(base_size = 13))
ggplot2::theme_update(
  text = element_text(color = "black", family = "sans"),
  axis.text = element_text(color = "black", face = "bold"),
  axis.title = element_text(color = "black", face = "bold"),
  plot.title = element_text(face = "bold", hjust = 0.5),
  strip.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold", hjust = 0.5),
  legend.title.position = "top",
  panel.grid = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
)


# ── Theme ─────────────────────────────────────────────────────────────────────
omni_theme <- create_theme(
  adminlte_color(light_blue = "#1B4965"),
  adminlte_sidebar(
    width = "290px",
    dark_bg = "#17202a",
    dark_color = "#ecf0f1",
    dark_hover_bg = "#1B4965",
    dark_hover_color = "#ffffff"
  ),
  adminlte_global(content_bg = "#f4f6f9")
)

omni_css <- "
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap');
  body, .content-wrapper, .main-footer, .main-header { font-family: 'Inter', sans-serif; }

  /* ── Button contrast fixes for inactive / default state ─────────────────── */
  .btn-default, .btn-default:link, .btn-default:visited {
    background-color:#e8edf2 !important;
    color:#1a2332 !important;
    border:1px solid #b0bec5 !important;
    font-weight:600 !important;
  }
  .btn-default:hover, .btn-default:focus {
    background-color:#cdd6de !important;
    color:#0d1a26 !important;
    border-color:#90a4ae !important;
  }
  /* fileInput Browse button */
  .btn-file { background-color:#e8edf2 !important; color:#1a2332 !important;
              border:1px solid #b0bec5 !important; font-weight:600 !important; }
  .btn-file:hover { background-color:#cdd6de !important; color:#0d1a26 !important; }
  /* Download links that don't use .dl-btn class */
  .shiny-download-link:not(.dl-btn) {
    color:#1a2332 !important; font-weight:600 !important;
  }

  .skin-blue .main-header .logo { background-color:#1B4965; font-weight:700; font-size:15px; letter-spacing:0.5px; }
  .skin-blue .main-header .navbar { background-color:#1B4965; }
  .skin-blue .main-sidebar { background-color:#17202a; }
  .skin-blue .sidebar-menu > li.active > a,
  .skin-blue .sidebar-menu > li > a:hover { border-left-color:#5FA8D3; background:#1B4965; }

  .box { border-radius:8px; box-shadow:0 2px 12px rgba(0,0,0,0.08); border-top:3px solid #1B4965 !important; }
  .box-header { border-radius:8px 8px 0 0;
                background:linear-gradient(135deg,#1B4965 0%,#2D6A8F 100%) !important; color:#fff !important; }
  .box-title { color:#fff !important; font-weight:600; }

  .dl-btn { display:block; width:90%; margin:6px auto;
             background:linear-gradient(135deg,#1a5632,#2D6A4F);
             color:#fff !important; border:none; border-radius:20px;
             padding:10px 16px; font-size:13px; font-weight:600;
             cursor:pointer; transition:all 0.25s ease; text-align:center;
             box-shadow:0 2px 8px rgba(26,86,50,0.3); text-shadow:0 1px 2px rgba(0,0,0,0.2); }
  .dl-btn:hover { background:linear-gradient(135deg,#145228,#256b45); transform:translateY(-1px);
                   box-shadow:0 4px 12px rgba(26,86,50,0.4); }

  .back-btn .btn { background:linear-gradient(135deg,#2c3e50,#34495e);
                   color:#fff; border:none; border-radius:20px;
                   padding:7px 18px; font-weight:600; font-size:13px;
                   transition:all 0.2s ease; box-shadow:0 2px 6px rgba(0,0,0,0.15); }
  .back-btn .btn:hover { background:linear-gradient(135deg,#1B4965,#5FA8D3); transform:translateY(-1px); }

  .info-box { border-radius:8px; }
  .nav-tabs-custom > .nav-tabs > li.active > a { border-top-color:#1B4965; }
  .plot-placeholder { display:flex; align-items:center; justify-content:center;
                      height:300px; color:#6c757d; font-size:16px;
                      border:2px dashed #dee2e6; border-radius:8px; margin:20px; }
  .plot-wrap { position:relative; min-height:200px; }
  .spinner-overlay { position:absolute; top:0; left:0; width:100%; height:100%;
                     display:flex; align-items:center; justify-content:center;
                     background:rgba(244,246,249,0.85); z-index:10;
                     font-size:2rem; color:#1B4965; }

  .tool-card { cursor:pointer; transition:transform 0.18s,box-shadow 0.18s; border-radius:14px;
               border:2px solid #dee2e6; padding:40px 30px; text-align:center; background:#fff;
               margin:12px; }
  .tool-card:hover { transform:translateY(-5px); box-shadow:0 10px 30px rgba(27,73,101,0.2);
                     border-color:#1B4965; }
  .tool-card .tool-icon { font-size:4rem; margin-bottom:20px; color:#1B4965; }
  .tool-card h3 { font-size:1.8rem; font-weight:700; color:#17202a; margin-bottom:12px; }
  .tool-card p  { font-size:1.15rem; color:#6c757d; line-height:1.5; }
  .tool-badge { display:inline-block; background:#1B4965; color:#fff;
                border-radius:14px; padding:5px 16px; font-size:13px; font-weight:700;
                margin-top:12px; letter-spacing:0.3px; }
  .omni-title    { font-size:3.2rem; font-weight:700; color:#17202a; margin-bottom:10px; }
  .omni-subtitle { font-size:1.35rem; color:#6c757d; margin-bottom:48px; }
  .back-btn { margin-bottom:16px; }

  /* Sidebar section label — used by module sidebars */
  .sidebar-section-label { padding:12px 16px 4px; color:#adb5bd; font-size:11px;
                            font-weight:700; text-transform:uppercase; letter-spacing:1px; }

  /* Inactive button contrast fixes */
  .btn-default { color: #1a1a1a !important; background-color: #f8f9fa !important; border-color: #ced4da !important; }
  .btn-default:hover, .btn-default:focus { color: #1a1a1a !important; background-color: #e2e6ea !important; border-color: #bdc6d0 !important; }
  input[type=file] label { color: #1a1a1a !important; }
  .shiny-download-link:not(.dl-btn) { color: #1a1a1a !important; font-weight: 500; }
"


# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = tags$span(
      icon("knight", style = "margin-right:6px;", lib = "glyphicon"),
      "proteOmni"
    ),
    titleWidth = 220
  ),
  dashboardSidebar(
    width = 290,
    useShinyjs(),
    tags$head(tags$style(HTML(omni_css))),
    use_theme(omni_theme),
    sidebarMenu(
      id = "main_menu",
      tags$div(
        style = "padding:10px 16px 10px 16px; margin-bottom:6px;",
        downloadButton(
          "download_log_history",
          tags$span(
            icon("file-lines", style = "margin-right:4px;"),
            "Download Log History"
          ),
          class = "dl-btn",
          style = "width:100%;text-align:center;font-size:12px;padding:8px 10px;"
        )
      ),
      menuItem("Home", tabName = "home", icon = icon("house")),
      menuItem(
        "PSManalyst",
        tabName = "psm_tab",
        icon = icon("barcode", lib = "glyphicon")
      ),
      menuItem("QC4DIANN", tabName = "qc_tab", icon = icon("circle-nodes")),
      menuItem("PwrQuant", tabName = "pwr_tab", icon = icon("calculator")),
      menuItem("Casanovo", tabName = "dnv_tab", icon = icon("fingerprint")),
      menuItem("InstaNovo", tabName = "ins_tab", icon = icon("brain")),
      menuItem("EncyclopeDIA", tabName = "enc_tab", icon = icon("book-open")),
      menuItem("Sage", tabName = "sag_tab", icon = icon("leaf")),
      menuItem("MaxQuant", tabName = "mqt_tab", icon = icon("table"))
    ),
    tags$hr(style = "border-color:#2d3741;margin:4px 0;"),
    uiOutput("active_sidebar_controls"),
    Fasta_sidebar_ui("global_fasta")
  ),
  dashboardBody(
    tabItems(
      # ── Home ──────────────────────────────────────────────────────────────
      tabItem(
        tabName = "home",
        fluidRow(
          column(
            12,
            align = "center",
            tags$br(),
            tags$h1("proteOmni", class = "omni-title"),
            tags$p(
              "Unified proteomics QC dashboard — select a tool to begin",
              class = "omni-subtitle"
            ),
            fluidRow(
              column(
                4,
                div(
                  class = "tool-card",
                  id = "card_pwr",
                  style = "background-color:#fcf3cf;border:2px solid #f1c40f;",
                  onclick = "Shiny.setInputValue('select_tool','pwr',{priority:'event'})",
                  div(class = "tool-icon", icon("calculator")),
                  tags$h3("PwrQuant"),
                  tags$p("Differential abundance mapping and Power analysis."),
                  tags$span("limma / stats", class = "tool-badge")
                )
              ),
              column(
                4,
                div(
                  class = "tool-card",
                  id = "card_psm",
                  onclick = "Shiny.setInputValue('select_tool','psm',{priority:'event'})",
                  div(class = "tool-icon", icon("barcode", lib = "glyphicon")),
                  tags$h3("PSManalyst"),
                  tags$p(
                    "Visual quality control for FragPipe DDA proteomics results."
                  ),
                  tags$span("FragPipe / DDA", class = "tool-badge")
                )
              ),
              column(
                4,
                div(
                  class = "tool-card",
                  id = "card_qc",
                  onclick = "Shiny.setInputValue('select_tool','qc4',{priority:'event'})",
                  div(class = "tool-icon", icon("circle-nodes")),
                  tags$h3("QC4DIANN"),
                  tags$p("QC report and diagnostics for DIA-NN results."),
                  tags$span("DIA-NN / DIA", class = "tool-badge")
                )
              )
            ),
            fluidRow(
              column(
                4,
                div(
                  class = "tool-card",
                  id = "card_mqt",
                  onclick = "Shiny.setInputValue('select_tool','mqt',{priority:'event'})",
                  div(class = "tool-icon", icon("table")),
                  tags$h3("MaxQuant"),
                  tags$p("QC report and diagnostics for MaxQuant results."),
                  tags$span("MaxQuant / DDA", class = "tool-badge")
                )
              ),
              column(
                4,
                div(
                  class = "tool-card",
                  id = "card_dnv",
                  onclick = "Shiny.setInputValue('select_tool','dnv',{priority:'event'})",
                  div(class = "tool-icon", icon("fingerprint")),
                  tags$h3("Casanovo"),
                  tags$p("Visualise Casanovo de novo sequencing results."),
                  tags$span("de novo", class = "tool-badge")
                )
              ),
              column(
                4,
                div(
                  class = "tool-card",
                  id = "card_ins",
                  onclick = "Shiny.setInputValue('select_tool','ins',{priority:'event'})",
                  div(class = "tool-icon", icon("brain")),
                  tags$h3("InstaNovo"),
                  tags$p("Visualise InstaNovo de novo sequencing results."),
                  tags$span("de novo", class = "tool-badge")
                )
              )
            ),
            fluidRow(
              column(
                4,
                div(
                  class = "tool-card",
                  id = "card_enc",
                  onclick = "Shiny.setInputValue('select_tool','enc',{priority:'event'})",
                  div(class = "tool-icon", icon("book-open")),
                  tags$h3("EncyclopeDIA"),
                  tags$p("Aggregate and explore EncyclopeDIA DIA results."),
                  tags$span("EncyclopeDIA / DIA", class = "tool-badge")
                )
              ),
              column(
                4,
                div(
                  class = "tool-card",
                  id = "card_sag",
                  onclick = "Shiny.setInputValue('select_tool','sag',{priority:'event'})",
                  div(class = "tool-icon", icon("leaf")),
                  tags$h3("Sage"),
                  tags$p(
                    "Visual QC and validation for Sage search engine results."
                  ),
                  tags$span("Sage / DDA / DIA", class = "tool-badge")
                )
              )
            ),
            tags$br(),
            tags$br(),
            tags$hr(
              style = "border-color:#dee2e6;max-width:800px;margin:24px auto;"
            ),
            div(
              style = "max-width:800px;margin:0 auto;text-align:left;color:#6c757d;font-size:13px;",
              tags$h4("Please cite:", style = "font-weight:600;color:#17202a;"),
              tags$ul(
                style = "padding-left:20px;line-height:1.6;",
                tags$li(
                  "Chaves AFA. PSManalyst: A Dashboard for Visual Quality Control of FragPipe Results. J Proteome Res. 2025 Sep 5;24(9):4344-4346. doi: 10.1021/acs.jproteome.5c00557. Epub 2025 Aug 15. PMID: 40815682."
                ),
                tags$li(
                  "Moschem JDC, de Barros BCSC, Serrano SMT, Chaves AFA. Decoding the Impact of Isolation Window Selection and QuantUMS Filtering in DIA-NN for DIA Quantification of Peptides and Proteins. J Proteome Res. 2025 Aug 1;24(8):3860-3873. doi: 10.1021/acs.jproteome.5c00009. Epub 2025 Jul 8. PMID: 40629671."
                )
              )
            )
          )
        )
      ),

      # ── PSManalyst ────────────────────────────────────────────────────────
      tabItem(
        tabName = "psm_tab",
        fluidRow(column(
          12,
          div(
            class = "back-btn",
            actionButton(
              "back_to_home_psm",
              "← Home",
              class = "btn btn-default btn-sm"
            )
          )
        )),
        PSManalyst_ui("psm")
      ),

      # ── QC4DIANN ──────────────────────────────────────────────────────────
      tabItem(
        tabName = "qc_tab",
        fluidRow(column(
          12,
          div(
            class = "back-btn",
            actionButton(
              "back_to_home_qc",
              "← Home",
              class = "btn btn-default btn-sm"
            )
          )
        )),
        QC4DIANN_ui("qc4")
      ),

      # ── Casanovo ─────────────────────────────────────────────────────────
      tabItem(
        tabName = "dnv_tab",
        fluidRow(column(
          12,
          div(
            class = "back-btn",
            actionButton(
              "back_to_home_dnv",
              "← Home",
              class = "btn btn-default btn-sm"
            )
          )
        )),
        deNovo_body_ui("dnv")
      ),

      # ── PwrQuant ──────────────────────────────────────────────────────────
      tabItem(
        tabName = "pwr_tab",
        fluidRow(column(
          12,
          div(
            class = "back-btn",
            actionButton(
              "back_to_home_pwr",
              "← Home",
              class = "btn btn-default btn-sm"
            )
          )
        )),
        PwrQuant_body_ui("pwr")
      ),

      # ── InstaNovo ─────────────────────────────────────────────────────────
      tabItem(
        tabName = "ins_tab",
        fluidRow(column(
          12,
          div(
            class = "back-btn",
            actionButton(
              "back_to_home_ins",
              "← Home",
              class = "btn btn-default btn-sm"
            )
          )
        )),
        InstaNovo_body_ui("ins")
      ),

      # ── EncyclopeDIA ──────────────────────────────────────────────────────
      tabItem(
        tabName = "enc_tab",
        fluidRow(column(
          12,
          div(
            class = "back-btn",
            actionButton(
              "back_to_home_enc",
              "← Home",
              class = "btn btn-default btn-sm"
            )
          )
        )),
        EncyclopeDIA_body_ui("enc")
      ),

      # ── Sage ──────────────────────────────────────────────────────────────
      tabItem(
        tabName = "sag_tab",
        fluidRow(column(
          12,
          div(
            class = "back-btn",
            actionButton(
              "back_to_home_sag",
              "← Home",
              class = "btn btn-default btn-sm"
            )
          )
        )),
        Sage_body_ui("sag")
      ),

      # ── MaxQuant ──────────────────────────────────────────────────────────
      tabItem(
        tabName = "mqt_tab",
        fluidRow(column(
          12,
          div(
            class = "back-btn",
            actionButton(
              "back_to_home_mqt",
              "← Home",
              class = "btn btn-default btn-sm"
            )
          )
        )),
        MaxQuantMSMS_body_ui("mqt")
      )
    )
  )
)


# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  onStop(function() {
    tryCatch(
      {
        writeLines(
          paste0("\n[proteOmni session ended at ", Sys.time(), "]"),
          con = .omni_log_con
        )
        close(.omni_log_con)
      },
      error = function(e) NULL
    )
  })

  # ── Tool card → tab routing ────────────────────────────────────────────────
  observeEvent(input$select_tool, {
    tab <- switch(
      input$select_tool,
      psm = "psm_tab",
      qc4 = "qc_tab",
      dnv = "dnv_tab",
      pwr = "pwr_tab",
      ins = "ins_tab",
      enc = "enc_tab",
      sag = "sag_tab",
      mqt = "mqt_tab"
    )
    if (!is.null(tab)) updateTabItems(session, "main_menu", selected = tab)
  })

  # ── Log download ───────────────────────────────────────────────────────────
  output$download_log_history <- downloadHandler(
    filename = function() {
      paste0(
        "proteOmni_log_history_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".txt"
      )
    },
    content = function(file) {
      file.copy(.omni_log_path, file)
      message("Log history downloaded by user at ", Sys.time())
    },
    contentType = "text/plain"
  )

  # ── Back buttons ──────────────────────────────────────────────────────────
  back_buttons <- list(
    back_to_home_psm = "home",
    back_to_home_qc = "home",
    back_to_home_dnv = "home",
    back_to_home_pwr = "home",
    back_to_home_ins = "home",
    back_to_home_enc = "home",
    back_to_home_sag = "home",
    back_to_home_mqt = "home"
  )
  for (btn in names(back_buttons)) {
    local({
      b <- btn
      observeEvent(input[[b]], {
        updateTabItems(session, "main_menu", selected = "home")
      })
    })
  }

  # ── Sidebar controls ───────────────────────────────────────────────────────
  output$active_sidebar_controls <- renderUI({
    tab <- input$main_menu
    if (is.null(tab)) {
      return(NULL)
    }

    switch(
      tab,
      psm_tab = PSManalyst_sidebar_ui("psm"),
      qc_tab = QC4DIANN_sidebar_ui("qc4"),
      pwr_tab = PwrQuant_sidebar_ui("pwr"),
      dnv_tab = deNovo_sidebar_ui("dnv"),
      ins_tab = InstaNovo_sidebar_ui("ins"),
      enc_tab = EncyclopeDIA_sidebar_ui("enc"),
      sag_tab = Sage_sidebar_ui("sag"),
      mqt_tab = MaxQuantMSMS_sidebar_ui("mqt"),
      tags$div(
        style = "padding:16px;color:#adb5bd;font-size:13px;",
        icon("arrow-left", style = "margin-right:6px;"),
        "Select a tool from the menu above to load its controls."
      )
    )
  })

  # ── Module servers ─────────────────────────────────────────────────────────
  PSManalyst_server("psm")
  QC4DIANN_server("qc4")
  deNovo_server("dnv")
  PwrQuant_server("pwr")

  global_fasta_digest <- Fasta_server("global_fasta")

  InstaNovo_server("ins", global_fasta_digest)
  EncyclopeDIA_server("enc", global_fasta_digest)
  Sage_server("sag", global_fasta_digest)
  MaxQuantMSMS_server("mqt")
}

shinyApp(ui = ui, server = server)
