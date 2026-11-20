## ============================================================
## mod_SeqProfile.r  —  Sequence profiling and in silico digestion
## ============================================================
## Paste a protein / peptide sequence, then use one of two tabs:
##   Prediction Algorithms  pick an algorithm, click Predict
##   In Silico Digestion    pick an enzyme (or a custom rule), click Digest
## All computation lives in utils_fasta.r (section "SEQUENCE PROFILING");
## this file is UI and reactive wiring only.
##
##   SeqProfile_body_ui(id)
##   SeqProfile_server(id)

# Candida albicans candidalysin (UniProt Q07730), used by "Load example sequence".
SEQPROFILE_EXAMPLE <- paste(
  ">sp|Q07730|CECE1_CANAL Extent of cell elongation protein 1 OS=Candida albicans (strain SC5314 / ATCC MYA-2876) OX=237561 GN=ECE1 PE=1 SV=2",
  "MKFSKIACATVFALSSQAAIIHHAPEFNMKRDVAPAAPAAPADQAPTVPAPQEFNTAITK",
  "RSIIGIIMGILGNIPQVIQIIMSIVKAFKGNKREDIDSVVAGIIADMPFVVRAVDTAMTS",
  "VASTKRDGANDDVANAVVRLPEIVARVATGVQQSIENAKRDGVPDVGLNLVANAPRLISN",
  "VFDGVSETVQQAKRDGLEDFLDELLQRLPQLITRSAESALKDSQPVKRDAGSVALSNLIK",
  "KSIETVGIENAAQIVSERDISSLIEEYFGKA",
  sep = "\n"
)

seqprofile_fmt <- function(x, digits = 3) {
  if (length(x) != 1 || is.na(x)) {
    "n/a"
  } else {
    format(round(x, digits), nsmall = 0)
  }
}

# One labelled figure of the summary strip
seqprofile_stat <- function(label, value) {
  tags$div(
    style = "flex:1 1 130px;background:#f4f6f8;border-radius:8px;padding:8px 12px;",
    tags$div(
      style = "font-size:11px;color:#6c757d;text-transform:uppercase;letter-spacing:.5px;",
      label
    ),
    tags$div(style = "font-size:16px;font-weight:700;color:#17202a;", value)
  )
}


# Sequence with the residues covered by in-range peptides highlighted,
# `per_line` residues per row, each row prefixed by its first position.
seqprofile_coverage_map <- function(sequence, covered, per_line = 60L) {
  aa <- strsplit(sequence, "")[[1]]
  n <- length(aa)
  rows <- lapply(seq(1L, n, by = per_line), function(s) {
    idx <- s:min(s + per_line - 1L, n)
    runs <- rle(covered[idx])
    run_end <- cumsum(runs$lengths)
    run_start <- run_end - runs$lengths + 1L
    spans <- Map(
      function(a, b, is_covered) {
        tags$span(
          style = if (is_covered) {
            "background:#bfe3c6;color:#17202a;"
          } else {
            "color:#c0392b;"
          },
          paste(aa[idx[a:b]], collapse = "")
        )
      },
      run_start,
      run_end,
      runs$values
    )
    tags$div(
      tags$span(
        style = "display:inline-block;width:5ch;margin-right:1.5ch;text-align:right;color:#6c757d;user-select:none;",
        s
      ),
      spans
    )
  })
  tags$div(
    style = "font-family:monospace;font-size:13px;line-height:1.7;white-space:nowrap;overflow-x:auto;",
    rows
  )
}


# ── BODY UI ──────────────────────────────────────────────────────────────────
SeqProfile_body_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # ── Sequence input, shared by both tabs ──────────────────────────────────
    fluidRow(
      box(
        title = "Sequence",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        textAreaInput(
          ns("sequence"),
          "Protein or peptide sequence (raw or FASTA)",
          rows = 6,
          width = "100%",
          placeholder = "Paste an amino acid sequence, with or without a >header line"
        ),
        tags$div(
          style = "display:flex;gap:6px;margin-bottom:8px;",
          actionButton(
            ns("load_example"),
            "Load example sequence",
            class = "btn btn-default btn-sm"
          ),
          actionButton(
            ns("clear"),
            "Clear",
            class = "btn btn-default btn-sm"
          )
        ),
        uiOutput(ns("seq_info"))
      )
    ),

    tabsetPanel(
      id = ns("tabs"),
      type = "tabs",

      # ── Tab 1: Prediction algorithms ───────────────────────────────────────
      tabPanel(
        "Prediction Algorithms",
        value = "predict",
        fluidRow(
          box(
            title = "Algorithm",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            # Controls run left to right and wrap; flex-end keeps the inputs
            # on one baseline when a label wraps onto two lines.
            tags$div(
              style = "display:flex;flex-wrap:wrap;align-items:flex-end;gap:0 16px;",
              selectInput(
                ns("algorithm"),
                "Algorithm",
                choices = seqprofile_choices(),
                width = "280px"
              ),
              uiOutput(
                ns("param_ui"),
                style = "display:flex;flex-wrap:wrap;align-items:flex-end;gap:0 16px;"
              ),
              actionButton(
                ns("predict"),
                "Predict",
                icon = icon("play"),
                class = "btn-primary",
                # Matches the bottom margin of the inputs' form-group
                style = "font-weight:bold;min-width:140px;margin-bottom:15px;"
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Profile Plot & Interpretation",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            tags$div(
              style = "text-align:right;margin-bottom:6px;",
              downloadButton(ns("download_plot"), "Download Plot (.png)")
            ),
            plotly::plotlyOutput(ns("profile_plot"), height = 420),
            uiOutput(ns("interpretation"))
          )
        ),
        fluidRow(
          box(
            title = "Summary",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            uiOutput(ns("summary_strip")),
            radioButtons(
              ns("table_view"),
              NULL,
              choices = c(
                "Predicted regions" = "regions",
                "Full per-position profile" = "profile"
              ),
              inline = TRUE
            ),
            DT::dataTableOutput(ns("summary_table"))
          )
        )
      ),

      # ── Tab 2: In silico digestion ─────────────────────────────────────────
      tabPanel(
        "In Silico Digestion",
        value = "digest",
        fluidRow(
          box(
            title = "Enzyme",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            fluidRow(
              # Enzyme and cleavage rule
              column(
                3,
                selectInput(
                  ns("enzyme"),
                  "Enzyme",
                  choices = seqprofile_enzyme_choices(),
                  width = "100%"
                ),
                uiOutput(ns("enzyme_rule")),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'custom'", ns("enzyme")),
                  textInput(
                    ns("custom_residues"),
                    "Cleavage residues (one-letter codes)",
                    placeholder = "e.g. KR",
                    width = "100%"
                  ),
                  radioButtons(
                    ns("custom_terminus"),
                    "Cleave at the",
                    choices = c(
                      "C-terminal side of the residues" = "C",
                      "N-terminal side of the residues" = "N"
                    )
                  ),
                  textInput(
                    ns("custom_blocked"),
                    "No cleavage next to (optional)",
                    placeholder = "e.g. P",
                    width = "100%"
                  )
                )
              ),
              # Missed cleavages and peptide length limits
              column(
                3,
                numericInput(
                  ns("max_missed"),
                  "Missed cleavages",
                  value = 0,
                  min = 0,
                  max = 5,
                  step = 1,
                  width = "100%"
                ),
                fluidRow(
                  column(
                    6,
                    numericInput(
                      ns("min_len"),
                      "Min length (coverage)",
                      value = 6,
                      min = 1,
                      step = 1,
                      width = "100%"
                    )
                  ),
                  column(
                    6,
                    numericInput(
                      ns("max_len"),
                      "Max length (coverage)",
                      value = 30,
                      min = 1,
                      step = 1,
                      width = "100%"
                    )
                  )
                )
              ),
              # m/z window
              column(
                3,
                checkboxInput(
                  ns("mz_filter"),
                  "Count only peptides with a 2+ or 3+ ion in the m/z range",
                  value = TRUE,
                  width = "100%"
                ),
                conditionalPanel(
                  condition = sprintf("input['%s']", ns("mz_filter")),
                  fluidRow(
                    column(
                      6,
                      numericInput(
                        ns("mz_min"),
                        "Min m/z",
                        value = 350,
                        min = 0,
                        step = 50,
                        width = "100%"
                      )
                    ),
                    column(
                      6,
                      numericInput(
                        ns("mz_max"),
                        "Max m/z",
                        value = 1500,
                        min = 0,
                        step = 50,
                        width = "100%"
                      )
                    )
                  )
                )
              ),
              # Modifications and run button
              column(
                3,
                checkboxInput(
                  ns("carbamidomethyl"),
                  "Carbamidomethyl Cys (+57.0215 Da)",
                  value = FALSE,
                  width = "100%"
                ),
                checkboxInput(
                  ns("oxidation"),
                  "Oxidation of Met (+15.9949 Da)",
                  value = FALSE,
                  width = "100%"
                ),
                actionButton(
                  ns("digest"),
                  "Digest",
                  icon = icon("scissors"),
                  class = "btn-primary",
                  style = "width:100%;font-weight:bold;margin-top:6px;"
                )
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Sequence Coverage",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            uiOutput(ns("digest_coverage"))
          )
        ),
        fluidRow(
          box(
            title = "Theoretical Peptides",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            tags$p(
              style = "font-size:12.5px;color:#6c757d;",
              "Click a peptide to view its theoretical b/y spectrum below."
            ),
            DT::dataTableOutput(ns("digest_table"))
          )
        ),
        fluidRow(
          box(
            title = "Theoretical b/y Spectrum",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            radioButtons(
              ns("fragment_charge"),
              "Fragment charge",
              choices = c("1+" = "1", "2+" = "2"),
              inline = TRUE
            ),
            plotly::plotlyOutput(ns("spectrum_plot"), height = 380)
          )
        )
      )
    )
  )
}


# ── SERVER ───────────────────────────────────────────────────────────────────
SeqProfile_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    registry <- seqprofile_registry()

    observeEvent(input$load_example, {
      updateTextAreaInput(session, "sequence", value = SEQPROFILE_EXAMPLE)
    })
    observeEvent(input$clear, {
      updateTextAreaInput(session, "sequence", value = "")
    })

    # ── Live sequence feedback (the only output that follows keystrokes) ────
    output$seq_info <- renderUI({
      cs <- clean_protein_sequence(input$sequence)
      parts <- list(tags$strong(sprintf("Length: %d aa", cs$length)))
      if (cs$n_entries > 1) {
        parts <- c(
          parts,
          list(tags$div(
            style = "color:#b9770e;",
            sprintf(
              "%d FASTA entries pasted; only the first is used.",
              cs$n_entries
            )
          ))
        )
      }
      if (length(cs$illegal_chars) > 0) {
        parts <- c(
          parts,
          list(tags$div(
            style = "color:#c0392b;",
            paste0(
              "Unsupported characters: ",
              paste(cs$illegal_chars, collapse = " "),
              ". Remove them before predicting."
            )
          ))
        )
      } else if (length(cs$invalid_chars) > 0) {
        parts <- c(
          parts,
          list(tags$div(
            style = "color:#b9770e;",
            paste0(
              "Non-canonical residues detected (",
              paste(cs$invalid_chars, collapse = ", "),
              "): windows containing them are not scored."
            )
          ))
        )
      }
      tags$div(style = "font-size:13px;", parts)
    })

    # ── Parameters of the selected algorithm ────────────────────────────────
    output$param_ui <- renderUI({
      spec <- registry[[req(input$algorithm)]]
      d <- spec$defaults
      min_w <- if (!is.null(spec$min_window)) spec$min_window else 2
      tagList(
        numericInput(
          ns("window"),
          "Window size (residues)",
          value = d$window,
          min = min_w,
          step = 1,
          width = "180px"
        ),
        if (!is.null(d$pH)) {
          numericInput(
            ns("ph"),
            "pH",
            value = d$pH,
            min = 0,
            max = 14,
            step = 0.1,
            width = "110px"
          )
        },
        if (!is.null(d$angle)) {
          selectInput(
            ns("angle"),
            "Rotation angle",
            choices = c(
              "100\u00B0 (alpha-helix)" = "100",
              "160\u00B0 (beta-strand)" = "160"
            ),
            selected = as.character(d$angle),
            width = "200px"
          )
        },
        if (!is.null(d$mask)) {
          textInput(
            ns("mask"),
            "Masked anchor positions (1-based, 'C' = C-term)",
            value = d$mask,
            placeholder = "e.g. 1,2,C",
            width = "300px"
          )
        },
        numericInput(
          ns("threshold"),
          "Threshold (leave empty for none)",
          value = d$threshold,
          step = 0.1,
          width = "220px"
        )
      )
    })

    # ── Prediction: runs only when Predict is clicked ───────────────────────
    result <- eventReactive(input$predict, {
      cs <- clean_protein_sequence(input$sequence)
      # jsonlite::validate() masks shiny's in proteOmni.
      shiny::validate(
        shiny::need(
          nzchar(cs$sequence),
          "Paste a sequence (or load the example) and click Predict."
        ),
        shiny::need(
          length(cs$illegal_chars) == 0,
          paste0(
            "The sequence contains unsupported characters: ",
            paste(cs$illegal_chars, collapse = " ")
          )
        )
      )
      if (cs$n_entries > 1) {
        showNotification(
          sprintf(
            "%d FASTA entries pasted; only the first is analysed.",
            cs$n_entries
          ),
          type = "warning"
        )
      }

      res <- tryCatch(
        run_seq_profile(
          cs$sequence,
          input$algorithm,
          list(
            window = input$window,
            pH = input$ph,
            angle = suppressWarnings(as.numeric(input$angle)),
            mask = input$mask,
            threshold = input$threshold
          )
        ),
        seqprofile_error = function(e) e
      )
      if (inherits(res, "condition")) {
        shiny::validate(shiny::need(FALSE, conditionMessage(res)))
      }
      if (length(res$masked_residues) > 0) {
        showNotification(
          sprintf(
            "%d of %d windows contain non-canonical residues (%s) and were not scored.",
            res$global_stats$n_masked,
            res$global_stats$n_windows,
            paste(res$masked_residues, collapse = ", ")
          ),
          type = "warning",
          duration = 8
        )
      }
      res$header <- cs$header
      res
    })

    # ── Plot ────────────────────────────────────────────────────────────────
    output$profile_plot <- plotly::renderPlotly({
      res <- result()
      plotly::ggplotly(plot_seq_profile(res), tooltip = "text") |>
        plotly::config(displaylogo = FALSE)
    })

    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("SeqProfile_", result()$algorithm_id, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot_seq_profile(result()),
          width = 11,
          height = 5,
          bg = "white",
          device = "png"
        )
      }
    )

    # ── Interpretation ──────────────────────────────────────────────────────
    output$interpretation <- renderUI({
      res <- result()
      it <- res$interpretation
      item <- function(label, text) {
        tags$p(style = "margin:0 0 6px;", tags$strong(label), " ", text)
      }
      tags$div(
        style = "margin-top:12px;padding:12px 16px;background:#eaf4fb;border-left:4px solid #5FA8D3;border-radius:6px;font-size:13.5px;",
        tags$h4(
          style = "margin-top:0;font-weight:700;",
          paste("How to read:", res$algorithm_name)
        ),
        item("What it measures.", it$summary),
        item("Reading the plot.", it$reading),
        item("Threshold and window.", it$thresholds),
        item("Caveats.", it$caveats),
        item("Reference.", it$reference)
      )
    })

    # ── Summary strip ───────────────────────────────────────────────────────
    output$summary_strip <- renderUI({
      res <- result()
      gs <- res$global_stats
      gd <- seqprofile_global_descriptors(res$sequence)
      with_pos <- function(score, pos) {
        if (is.na(score)) {
          "n/a"
        } else {
          sprintf("%s (pos %d)", seqprofile_fmt(score), pos)
        }
      }
      row <- function(...) {
        tags$div(
          style = "display:flex;flex-wrap:wrap;gap:8px;margin-bottom:8px;",
          ...
        )
      }
      tagList(
        row(
          seqprofile_stat("Sequence length", paste(gs$seq_length, "aa")),
          seqprofile_stat("Mean score", seqprofile_fmt(gs$mean_score)),
          seqprofile_stat("Min score", with_pos(gs$min_score, gs$min_pos)),
          seqprofile_stat("Max score", with_pos(gs$max_score, gs$max_pos)),
          seqprofile_stat(
            "Windows above threshold",
            if (is.na(gs$pct_above)) {
              "n/a"
            } else {
              paste0(seqprofile_fmt(gs$pct_above, 1), " %")
            }
          ),
          seqprofile_stat("Predicted regions", gs$n_regions)
        ),
        tags$div(
          style = "font-size:11px;color:#6c757d;margin:2px 0 4px;",
          "Whole-sequence descriptors (Peptides package; n/a when the sequence has non-canonical residues)"
        ),
        row(
          seqprofile_stat("MW (Da)", seqprofile_fmt(gd$mw, 2)),
          seqprofile_stat("pI (Lehninger)", seqprofile_fmt(gd$pI, 2)),
          seqprofile_stat("GRAVY", seqprofile_fmt(gd$gravy)),
          seqprofile_stat("Boman index", seqprofile_fmt(gd$boman)),
          seqprofile_stat(
            "Aliphatic index",
            seqprofile_fmt(gd$aliphatic_index, 2)
          ),
          seqprofile_stat(
            "Instability index",
            seqprofile_fmt(gd$instability_index, 2)
          )
        )
      )
    })

    # ── Summary table ───────────────────────────────────────────────────────
    table_data <- reactive({
      res <- result()
      if (identical(input$table_view, "profile")) {
        d <- res$profile
        names(d) <- c(
          "Position",
          "Residue",
          "Window_Sequence",
          "Score",
          "Above_Threshold"
        )
      } else {
        d <- res$regions
        names(d) <- c(
          "Region_ID",
          "Start",
          "End",
          "Length",
          "Sequence",
          "Peak_Position",
          "Peak_Score",
          "Mean_Score",
          "Note"
        )
      }
      num <- vapply(d, is.double, logical(1))
      d[num] <- lapply(d[num], round, 4)
      d
    })

    # server = FALSE so the export buttons write every row, not just the page
    output$summary_table <- DT::renderDataTable(
      {
        DT::datatable(
          table_data(),
          options = list(
            scrollX = TRUE,
            pageLength = 15,
            dom = "Bfrtip",
            buttons = c("copy", "csv", "excel"),
            language = list(
              emptyTable = "No region crosses the threshold (or no threshold is set)."
            )
          ),
          rownames = FALSE,
          extensions = "Buttons",
          class = "cell-border stripe hover"
        )
      },
      server = FALSE
    )

    # ── Digestion: cleavage rule of the selected enzyme ─────────────────────
    output$enzyme_rule <- renderUI({
      spec <- SEQPROFILE_ENZYMES[[req(input$enzyme)]]
      req(spec)
      tags$p(
        style = "font-size:12.5px;color:#6c757d;margin-top:-6px;",
        paste0(
          seqprofile_enzyme_rule(spec$residues, spec$terminus, spec$blocked_by),
          "."
        )
      )
    })

    # ── Digestion: runs only when Digest is clicked ─────────────────────────
    digest_result <- eventReactive(input$digest, {
      cs <- clean_protein_sequence(input$sequence)
      shiny::validate(
        shiny::need(
          nzchar(cs$sequence),
          "Paste a sequence (or load the example) and click Digest."
        ),
        shiny::need(
          length(cs$illegal_chars) == 0,
          paste0(
            "The sequence contains unsupported characters: ",
            paste(cs$illegal_chars, collapse = " ")
          )
        )
      )
      if (cs$n_entries > 1) {
        showNotification(
          sprintf(
            "%d FASTA entries pasted; only the first is digested.",
            cs$n_entries
          ),
          type = "warning"
        )
      }
      res <- tryCatch(
        run_seq_digest(
          cs$sequence,
          input$enzyme,
          list(
            residues = input$custom_residues,
            terminus = input$custom_terminus,
            blocked_by = input$custom_blocked,
            max_missed = input$max_missed,
            min_len = input$min_len,
            max_len = input$max_len,
            mz_min = if (isTRUE(input$mz_filter)) input$mz_min else NA,
            mz_max = if (isTRUE(input$mz_filter)) input$mz_max else NA,
            carbamidomethyl = input$carbamidomethyl,
            oxidation = input$oxidation
          )
        ),
        seqprofile_error = function(e) e
      )
      if (inherits(res, "condition")) {
        shiny::validate(shiny::need(FALSE, conditionMessage(res)))
      }
      res
    })

    output$digest_coverage <- renderUI({
      res <- digest_result()
      st <- res$stats
      p <- res$params
      counted <- sprintf("%d-%d aa", p$min_len, p$max_len)
      if (!is.na(p$mz_min)) {
        counted <- sprintf(
          "%s with a 2+ or 3+ ion at m/z %s-%s",
          counted,
          format(p$mz_min),
          format(p$mz_max)
        )
      }
      mods <- c(
        if (p$carbamidomethyl) "carbamidomethylated Cys" else "unmodified Cys",
        if (p$oxidation) "oxidised Met" else "unmodified Met"
      )
      tagList(
        tags$div(
          style = "display:flex;flex-wrap:wrap;gap:8px;margin-bottom:8px;",
          seqprofile_stat("Enzyme", res$enzyme_name),
          seqprofile_stat("Cleavage sites", st$n_sites),
          seqprofile_stat("Peptides", st$n_peptides),
          seqprofile_stat("Peptides counted", st$n_in_range),
          seqprofile_stat(
            "Sequence coverage",
            sprintf(
              "%s %% (%d / %d aa)",
              seqprofile_fmt(st$coverage_pct, 1),
              st$covered_residues,
              st$seq_length
            )
          )
        ),
        tags$div(
          style = "font-size:12.5px;color:#6c757d;margin-bottom:8px;",
          sprintf(
            "%s. Masses are monoisotopic with %s. Coverage counts residues in at least one peptide of %s (highlighted); residues in red are only in peptides outside these limits.",
            res$rule,
            paste(mods, collapse = " and "),
            counted
          )
        ),
        seqprofile_coverage_map(res$sequence, res$covered)
      )
    })

    digest_table_data <- reactive({
      pep <- digest_result()$peptides
      data.frame(
        Peptide = seq_len(nrow(pep)),
        Start = pep$start,
        End = pep$end,
        Sequence = pep$sequence,
        Residues = pep$length,
        Monoisotopic_Mass = round(pep$mono_mass, 4),
        `m/z (1+)` = round(pep$mz_1, 4),
        `m/z (2+)` = round(pep$mz_2, 4),
        `m/z (3+)` = round(pep$mz_3, 4),
        Missed_Cleavages = pep$missed_cleavages,
        Counted_For_Coverage = ifelse(pep$in_range, "yes", "no"),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })

    output$digest_table <- DT::renderDataTable(
      {
        DT::datatable(
          digest_table_data(),
          options = list(
            scrollX = TRUE,
            pageLength = 15,
            dom = "Bfrtip",
            buttons = c("copy", "csv", "excel")
          ),
          rownames = FALSE,
          selection = "single",
          extensions = "Buttons",
          class = "cell-border stripe hover"
        )
      },
      server = FALSE
    )

    # ── Digestion: b/y spectrum of the clicked peptide ──────────────────────
    spectrum_ions <- reactive({
      res <- digest_result()
      i <- input$digest_table_rows_selected
      shiny::validate(shiny::need(
        length(i) == 1 && i <= nrow(res$peptides),
        "Click a peptide in the table above to view its theoretical b/y spectrum."
      ))
      charge <- suppressWarnings(as.integer(input$fragment_charge))
      if (length(charge) != 1 || is.na(charge)) {
        charge <- 1L
      }
      ions <- tryCatch(
        seqprofile_fragment_ions(
          res$peptides$sequence[i],
          carbamidomethyl = res$params$carbamidomethyl,
          oxidation = res$params$oxidation,
          charge = charge
        ),
        seqprofile_error = function(e) e
      )
      if (inherits(ions, "condition")) {
        shiny::validate(shiny::need(FALSE, conditionMessage(ions)))
      }
      attr(ions, "peptide") <- res$peptides$sequence[i]
      ions
    })

    output$spectrum_plot <- plotly::renderPlotly({
      ions <- spectrum_ions()
      plotly::ggplotly(
        plot_fragment_spectrum(ions, title = attr(ions, "peptide")),
        tooltip = "text"
      ) |>
        plotly::config(displaylogo = FALSE)
    })
  })
}
