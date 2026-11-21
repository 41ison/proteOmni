## ============================================================
## mod_SeqProfile.r  —  Sequence profiling and in silico digestion
## ============================================================
## Paste a protein / peptide sequence (or fetch one from UniProt by accession),
## then use one of two tabs:
##   Prediction Algorithms  pick one to three algorithms, click Predict.
##                          Several algorithms give stacked tracks on a shared
##                          x axis and a consensus-region table; regions of the
##                          Eisenberg moment at 100 degrees open a helical wheel
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
          style = "display:flex;flex-wrap:wrap;align-items:center;gap:6px;margin-bottom:8px;",
          actionButton(
            ns("load_example"),
            "Load example sequence",
            class = "btn btn-default btn-sm"
          ),
          actionButton(
            ns("clear"),
            "Clear",
            class = "btn btn-default btn-sm"
          ),
          tags$span(
            style = "margin-left:14px;font-size:13px;color:#6c757d;",
            "or UniProt accession:"
          ),
          tags$input(
            id = ns("accession"),
            type = "text",
            class = "shiny-input-text form-control input-sm",
            style = "width:130px;",
            placeholder = "e.g. Q07730"
          ),
          actionButton(
            ns("fetch_uniprot"),
            "Fetch from UniProt",
            icon = icon("cloud-arrow-down"),
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
            tags$div(
              style = "display:flex;flex-wrap:wrap;align-items:flex-end;gap:0 16px;",
              selectizeInput(
                ns("algorithm"),
                "Algorithms (up to 3; several are stacked on a shared x axis)",
                choices = seqprofile_choices(),
                selected = names(seqprofile_registry())[1],
                multiple = TRUE,
                options = list(maxItems = 3, plugins = list("remove_button")),
                width = "560px"
              ),
              actionButton(
                ns("predict"),
                "Predict",
                icon = icon("play"),
                class = "btn-primary",
                # Matches the bottom margin of the inputs' form-group
                style = "font-weight:bold;min-width:140px;margin-bottom:15px;"
              )
            ),
            # One row of parameters per selected algorithm
            uiOutput(ns("param_ui"))
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
            # Rendered server-side: the height grows with the number of tracks
            uiOutput(ns("profile_plot_ui")),
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
            uiOutput(ns("summary_algorithm_ui")),
            uiOutput(ns("summary_strip")),
            radioButtons(
              ns("table_view"),
              NULL,
              choices = c(
                "Predicted regions" = "regions",
                "Full per-position profile" = "profile",
                "Consensus regions (all selected algorithms)" = "consensus"
              ),
              inline = TRUE
            ),
            conditionalPanel(
              condition = sprintf(
                "input['%s'] == 'consensus'",
                ns("table_view")
              ),
              numericInput(
                ns("consensus_min_len"),
                "Minimum consensus length (positions)",
                value = 1,
                min = 1,
                step = 1,
                width = "260px"
              )
            ),
            uiOutput(ns("table_note")),
            DT::dataTableOutput(ns("summary_table")),
            uiOutput(ns("wheel_ui"))
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

    # ── UniProt: canonical FASTA of an accession into the text area ─────────
    observeEvent(input$fetch_uniprot, {
      entry <- tryCatch(
        fetch_uniprot_fasta(input$accession),
        seqprofile_error = function(e) e
      )
      if (inherits(entry, "condition")) {
        showNotification(conditionMessage(entry), type = "error", duration = 8)
        return()
      }
      updateTextAreaInput(session, "sequence", value = entry$fasta)
      showNotification(
        sprintf(
          "%s: %s (%d aa) loaded from UniProt.",
          entry$accession,
          entry$protein_name,
          nchar(entry$sequence)
        ),
        type = "message"
      )
    })

    # ── Live sequence feedback (the only output that follows keystrokes) ────
    output$seq_info <- renderUI({
      cs <- clean_protein_sequence(input$sequence)
      parts <- list(tags$strong(sprintf("Length: %d aa", cs$length)))
      if (!is.na(cs$header)) {
        parts <- c(
          parts,
          list(tags$div(style = "color:#6c757d;", cs$header))
        )
      }
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
    # Input ids are "<parameter>_<algorithm id>", e.g. window_hydropathy_kd.
    selected_algorithms <- reactive({
      intersect(input$algorithm, names(registry))
    })

    output$param_ui <- renderUI({
      ids <- selected_algorithms()
      if (length(ids) == 0) {
        return(tags$p(
          style = "color:#6c757d;",
          "Select one to three algorithms."
        ))
      }
      # Adding or removing an algorithm redraws every row, so values already
      # entered are read back instead of being reset to the defaults.
      current <- function(param, id, default) {
        value <- isolate(input[[paste0(param, "_", id)]])
        if (is.null(value)) default else value
      }
      lapply(ids, function(id) {
        spec <- registry[[id]]
        d <- spec$defaults
        min_w <- if (!is.null(spec$min_window)) spec$min_window else 2
        tags$div(
          tags$div(style = "font-weight:700;margin-bottom:2px;", spec$name),
          # Controls run left to right and wrap; flex-end keeps the inputs
          # on one baseline when a label wraps onto two lines.
          tags$div(
            style = "display:flex;flex-wrap:wrap;align-items:flex-end;gap:0 16px;",
            numericInput(
              ns(paste0("window_", id)),
              "Window size (residues)",
              value = current("window", id, d$window),
              min = min_w,
              step = 1,
              width = "180px"
            ),
            if (!is.null(d$pH)) {
              numericInput(
                ns(paste0("ph_", id)),
                "pH",
                value = current("ph", id, d$pH),
                min = 0,
                max = 14,
                step = 0.1,
                width = "110px"
              )
            },
            if (!is.null(d$angle)) {
              selectInput(
                ns(paste0("angle_", id)),
                "Rotation angle",
                choices = c(
                  "100\u00B0 (alpha-helix)" = "100",
                  "160\u00B0 (beta-strand)" = "160"
                ),
                selected = as.character(current("angle", id, d$angle)),
                width = "200px"
              )
            },
            if (!is.null(d$mask)) {
              textInput(
                ns(paste0("mask_", id)),
                "Masked anchor positions (1-based, 'C' = C-term)",
                value = current("mask", id, d$mask),
                placeholder = "e.g. 1,2,C",
                width = "300px"
              )
            },
            numericInput(
              ns(paste0("threshold_", id)),
              "Threshold (leave empty for none)",
              value = current("threshold", id, d$threshold),
              step = 0.1,
              width = "220px"
            )
          )
        )
      })
    })

    # ── Prediction: runs only when Predict is clicked ───────────────────────
    # Named list with one run_seq_profile() result per selected algorithm.
    results <- eventReactive(input$predict, {
      cs <- clean_protein_sequence(input$sequence)
      ids <- selected_algorithms()
      # jsonlite::validate() masks shiny's in proteOmni.
      shiny::validate(
        shiny::need(
          length(ids) >= 1 && length(ids) <= 3,
          "Select one to three algorithms."
        ),
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

      out <- lapply(ids, function(id) {
        param <- function(name) input[[paste0(name, "_", id)]]
        res <- tryCatch(
          run_seq_profile(
            cs$sequence,
            id,
            list(
              window = param("window"),
              pH = param("ph"),
              angle = suppressWarnings(as.numeric(param("angle"))),
              mask = param("mask"),
              threshold = param("threshold")
            )
          ),
          seqprofile_error = function(e) e
        )
        if (inherits(res, "condition")) {
          msg <- conditionMessage(res)
          if (length(ids) > 1) {
            msg <- paste0(registry[[id]]$name, ": ", msg)
          }
          shiny::validate(shiny::need(FALSE, msg))
        }
        res$header <- cs$header
        res
      })
      names(out) <- ids

      # Window sizes differ between algorithms, so report the worst case
      masked <- vapply(out, function(r) r$global_stats$n_masked, numeric(1))
      worst <- out[[which.max(masked)]]
      if (length(worst$masked_residues) > 0) {
        showNotification(
          sprintf(
            "%s%d of %d windows contain non-canonical residues (%s) and were not scored.",
            if (length(out) > 1) "Up to " else "",
            worst$global_stats$n_masked,
            worst$global_stats$n_windows,
            paste(worst$masked_residues, collapse = ", ")
          ),
          type = "warning",
          duration = 8
        )
      }
      out
    })

    # The algorithm shown in the Summary box (the first one unless changed)
    result <- reactive({
      res <- results()
      id <- input$summary_algorithm
      if (is.null(id) || !id %in% names(res)) {
        id <- names(res)[1]
      }
      res[[id]]
    })

    # ── Consensus regions of the last prediction ────────────────────────────
    # Follows the minimum-length input without a new Predict. Returns the
    # condition instead of failing, so the table can show its message while
    # the plot simply draws no consensus bands.
    consensus <- reactive({
      all_res <- results()
      min_len <- input$consensus_min_len
      tryCatch(
        seqprofile_consensus_regions(
          all_res,
          min_length = if (is.null(min_len)) 1L else min_len
        ),
        seqprofile_error = function(e) e
      )
    })
    consensus_bands <- reactive({
      d <- consensus()
      if (is.data.frame(d)) d else NULL
    })

    # ── Plot: one track per algorithm on a shared x axis ────────────────────
    output$profile_plot_ui <- renderUI({
      # Before the first Predict (or on a validation error) keep one track,
      # so the message of profile_plot is still shown.
      n <- tryCatch(length(results()), error = function(e) 1L)
      bands <- tryCatch(consensus_bands(), error = function(e) NULL)
      tagList(
        plotly::plotlyOutput(
          ns("profile_plot"),
          height = if (n > 1) 260 * n else 420
        ),
        if (n > 1) {
          tags$p(
            style = "font-size:12.5px;color:#6c757d;margin-top:6px;",
            if (is.null(bands)) {
              "No consensus bands are drawn: every algorithm needs a threshold and a valid minimum consensus length."
            } else {
              sprintf(
                "Blue bands: regions of each algorithm. Orange bands: the %d consensus region(s), where all selected algorithms cross their thresholds (minimum length set under Summary > Consensus regions).",
                nrow(bands)
              )
            }
          )
        }
      )
    })

    output$profile_plot <- plotly::renderPlotly({
      res <- results()
      bands <- consensus_bands()
      tracks <- lapply(res, function(r) {
        plotly::ggplotly(
          plot_seq_profile(r, consensus = bands),
          tooltip = "text"
        )
      })
      if (length(tracks) == 1) {
        p <- tracks[[1]]
      } else {
        p <- plotly::subplot(
          tracks,
          nrows = length(tracks),
          shareX = TRUE,
          titleY = TRUE,
          margin = 0.035
        ) |>
          plotly::layout(title = list(text = ""))
        # subplot() keeps one title only, so each track is labelled inside
        # its own panel ("y domain", "y2 domain", ... refer to the panels).
        for (i in seq_along(res)) {
          p <- plotly::add_annotations(
            p,
            text = paste0("<b>", res[[i]]$algorithm_name, "</b>"),
            x = 0,
            y = 1,
            xref = "paper",
            yref = paste0("y", if (i > 1) i else "", " domain"),
            xanchor = "left",
            yanchor = "bottom",
            showarrow = FALSE,
            font = list(size = 12)
          )
        }
      }
      plotly::config(p, displaylogo = FALSE)
    })

    output$download_plot <- downloadHandler(
      filename = function() {
        paste0(
          "SeqProfile_",
          paste(names(results()), collapse = "-"),
          "_",
          Sys.Date(),
          ".png"
        )
      },
      content = function(file) {
        plots <- lapply(
          results(),
          plot_seq_profile,
          consensus = consensus_bands()
        )
        ggplot2::ggsave(
          file,
          if (length(plots) == 1) {
            plots[[1]]
          } else {
            patchwork::wrap_plots(plots, ncol = 1)
          },
          width = 11,
          height = if (length(plots) == 1) 5 else 3.5 * length(plots),
          bg = "white",
          device = "png"
        )
      }
    )

    # ── Interpretation ──────────────────────────────────────────────────────
    output$interpretation <- renderUI({
      all_res <- results()
      item <- function(label, text) {
        tags$p(style = "margin:0 0 6px;", tags$strong(label), " ", text)
      }
      lapply(all_res, function(res) {
        it <- res$interpretation
        body <- tagList(
          item("What it measures.", it$summary),
          item("Reading the plot.", it$reading),
          item("Threshold and window.", it$thresholds),
          item("Caveats.", it$caveats),
          item("Reference.", it$reference)
        )
        heading <- paste("How to read:", res$algorithm_name)
        tags$div(
          style = "margin-top:12px;padding:12px 16px;background:#eaf4fb;border-left:4px solid #5FA8D3;border-radius:6px;font-size:13.5px;",
          if (length(all_res) == 1) {
            tagList(
              tags$h4(style = "margin-top:0;font-weight:700;", heading),
              body
            )
          } else {
            # Collapsed, so three guides do not push the summary off screen
            tags$details(
              tags$summary(
                style = "font-size:15px;font-weight:700;cursor:pointer;",
                heading
              ),
              tags$div(style = "margin-top:8px;", body)
            )
          }
        )
      })
    })

    # ── Summary: which algorithm the strip and tables describe ──────────────
    output$summary_algorithm_ui <- renderUI({
      res <- results()
      if (length(res) < 2) {
        return(NULL)
      }
      choices <- stats::setNames(
        names(res),
        vapply(res, function(r) r$algorithm_name, character(1))
      )
      previous <- isolate(input$summary_algorithm)
      selectInput(
        ns("summary_algorithm"),
        "Statistics, regions and profile of",
        choices = choices,
        selected = if (isTRUE(previous %in% choices)) previous else choices[1],
        width = "420px"
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
            if (identical(res$direction, "below")) {
              "Windows below threshold"
            } else {
              "Windows above threshold"
            },
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
      if (identical(input$table_view, "consensus")) {
        all_res <- results()
        shiny::validate(shiny::need(
          length(all_res) >= 2,
          "Select two or three algorithms and click Predict to call consensus regions."
        ))
        d <- consensus()
        if (inherits(d, "condition")) {
          shiny::validate(shiny::need(FALSE, conditionMessage(d)))
        }
        names(d) <- c(
          "Region_ID",
          "Start",
          "End",
          "Length",
          "Sequence",
          paste0(
            "Mean: ",
            vapply(all_res, function(r) r$algorithm_name, character(1))
          )
        )
      } else if (identical(input$table_view, "profile")) {
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

    output$table_note <- renderUI({
      if (!identical(input$table_view, "consensus")) {
        return(NULL)
      }
      tags$p(
        style = "font-size:12.5px;color:#6c757d;",
        "Stretches of window-centre positions where every selected algorithm crosses its own threshold at the same time. Positions an algorithm does not score (sequence ends, masked windows) never count. The minimum length also applies to the orange consensus bands in the plot."
      )
    })

    # The helical wheel applies to regions of the Eisenberg moment at 100 degrees
    wheel_available <- reactive({
      res <- result()
      identical(res$algorithm_id, "hmoment_eisenberg") &&
        isTRUE(res$params$angle == 100) &&
        !identical(input$table_view, "profile") &&
        !identical(input$table_view, "consensus")
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
              emptyTable = if (identical(input$table_view, "consensus")) {
                "No stretch of the minimum length crosses the thresholds of all selected algorithms at once."
              } else {
                "No region crosses the threshold (or no threshold is set)."
              }
            )
          ),
          rownames = FALSE,
          selection = if (wheel_available()) "single" else "none",
          extensions = "Buttons",
          class = "cell-border stripe hover"
        )
      },
      server = FALSE
    )

    # ── Helical wheel of the clicked amphiphilic region ─────────────────────
    output$wheel_ui <- renderUI({
      if (!wheel_available()) {
        return(NULL)
      }
      tags$div(
        style = "margin-top:14px;",
        tags$h4(
          style = "font-weight:700;",
          "Helical wheel (Schiffer-Edmundson projection)"
        ),
        tags$p(
          style = "font-size:12.5px;color:#6c757d;",
          "Click a region in the table to project it down the helix axis (100\u00B0 per residue). Residues are coloured by Eisenberg hydrophobicity and joined in sequence order; the red arrow is the hydrophobic moment and points at the hydrophobic face. Residues 19 onwards sit on outer rings."
        ),
        plotOutput(ns("wheel_plot"), height = 480)
      )
    })

    wheel_data <- reactive({
      req(wheel_available())
      res <- result()
      shiny::validate(shiny::need(
        nrow(res$regions) > 0,
        "No region to project: enter a threshold for the hydrophobic moment and click Predict."
      ))
      i <- input$summary_table_rows_selected
      shiny::validate(shiny::need(
        length(i) == 1 && i <= nrow(res$regions),
        "Click a region in the table above to draw its helical wheel."
      ))
      segment <- seqprofile_wheel_segment(res, i)
      wheel <- tryCatch(
        seqprofile_helical_wheel(segment$sequence, segment$start, angle = 100),
        seqprofile_error = function(e) e
      )
      if (inherits(wheel, "condition")) {
        shiny::validate(shiny::need(FALSE, conditionMessage(wheel)))
      }
      attr(wheel, "segment") <- segment
      attr(wheel, "region_id") <- res$regions$region_id[i]
      wheel
    })

    output$wheel_plot <- renderPlot({
      wheel <- wheel_data()
      segment <- attr(wheel, "segment")
      plot_helical_wheel(
        wheel,
        title = sprintf(
          "Region %d: residues %d-%d",
          attr(wheel, "region_id"),
          segment$start,
          segment$end
        ),
        subtitle = sprintf(
          "%s\nuH of the drawn residues = %.3f%s",
          segment$sequence,
          attr(wheel, "moment")$uH,
          if (segment$truncated) {
            " (long region: 18 residues around its peak)"
          } else {
            ""
          }
        )
      )
    })

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
