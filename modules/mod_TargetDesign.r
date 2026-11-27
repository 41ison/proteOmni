## ============================================================
## mod_TargetDesign.r  —  Targeted Proteomics Assay Design (SRM / MRM / PRM)
## ============================================================
## Surrogate peptide candidate selection, FASTA proteotypicity verification,
## and product-ion transition ranking following CPTAC / Hoofnagle et al.,
## Clin. Chem. 2016;62(1):48-69 (PMC4830481).
##
## Exported functions:
##   TargetDesign_sidebar_ui(id)
##   TargetDesign_body_ui(id)
##   TargetDesign_server(id)

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

TARGETDESIGN_EXAMPLES <- list(
  apoa1 = list(
    label = "Human Apolipoprotein A-I (P02647)",
    fasta = paste(
      ">sp|P02647|APOA1_HUMAN Apolipoprotein A-I OS=Homo sapiens OX=9606 GN=APOA1 PE=1 SV=1",
      "MKAAVLTLAVLFLTGSQARHFWQQDEPPQSPWDRVKDLATVYVDVLKDSGRDYVSQFEGS",
      "ALGKQLNLKLLDNWDSVTSTFSKLREQLGPVTQEFWDNLEKETEGLRQEMSKDLEEVKAK",
      "VQPYLDDFQKKWQEEMELYRQKVEPLRAELQEGARQKLHELQEKLSPLGEEMRDRARAHV",
      "DALRTHLAPYSDELRQRLAARLEALKENGGARLAEYHAKATEHLSTLSEKAKPALEDLRQ",
      "GLLPVLESFKVSFLSALEEYTKKLNTQ",
      sep = "\n"
    )
  ),
  cst3 = list(
    label = "Human Cystatin-C (P01034)",
    fasta = paste(
      ">sp|P01034|CYTC_HUMAN Cystatin-C OS=Homo sapiens OX=9606 GN=CST3 PE=1 SV=1",
      "MAGPLRAPLLLLAILAVALAVSPAAGSSPGKPPRLVGGPMDASVEEEGVRRALDFAVGEY",
      "NKASNDMYHSRALQVVRARKQIVAGVNYFLDVELGRTTCTKTQPNLDNCPFHDQPHLKRK",
      "AFCSFQIYAVPWQGTMTLSKSTCQDA",
      sep = "\n"
    )
  ),
  ece1 = list(
    label = "Candida albicans ECE1 (Q07730)",
    fasta = paste(
      ">sp|Q07730|CECE1_CANAL Extent of cell elongation protein 1 OS=Candida albicans OX=237561 GN=ECE1 PE=1 SV=2",
      "MKFSKIACATVFALSSQAAIIHHAPEFNMKRDVAPAAPAAPADQAPTVPAPQEFNTAITK",
      "RSIIGIIMGILGNIPQVIQIIMSIVKAFKGNKREDIDSVVAGIIADMPFVVRAVDTAMTS",
      "VASTKRDGANDDVANAVVRLPEIVARVATGVQQSIENAKRDGVPDVGLNLVANAPRLISN",
      "VFDGVSETVQQAKRDGLEDFLDELLQRLPQLITRSAESALKDSQPVKRDAGSVALSNLIK",
      "KSIETVGIENAAQIVSERDISSLIEEYFGKA",
      sep = "\n"
    )
  )
)

targetdesign_stat <- function(label, value, accent = "#17202a") {
  tags$div(
    style = "flex:1 1 135px;background:#f4f6f8;border-radius:8px;padding:8px 12px;",
    tags$div(
      style = "font-size:11px;color:#6c757d;text-transform:uppercase;letter-spacing:.5px;",
      label
    ),
    tags$div(
      style = sprintf("font-size:16px;font-weight:700;color:%s;", accent),
      value
    )
  )
}

# Protein sequence map coloured by best candidate tier covering each residue:
# Recommended (green), Conditional (amber), or None / Avoid (muted red).
targetdesign_sequence_map <- function(sequence, tier_state, per_line = 60L) {
  aa <- strsplit(sequence, "")[[1]]
  n <- length(aa)
  if (length(tier_state) != n) {
    tier_state <- rep("None", n)
  }
  style_for <- function(state) {
    switch(
      state,
      Recommended = "background:#bfe3c6;color:#17202a;font-weight:600;",
      Conditional = "background:#fdebd0;color:#7d5a00;",
      "color:#95a5a6;"
    )
  }
  rows <- lapply(seq(1L, n, by = per_line), function(s) {
    idx <- s:min(s + per_line - 1L, n)
    runs <- rle(tier_state[idx])
    run_end <- cumsum(runs$lengths)
    run_start <- run_end - runs$lengths + 1L
    spans <- Map(
      function(a, b, st) {
        tags$span(
          style = style_for(st),
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


# ═══════════════════════════════════════════════════════════════════════════════
# SIDEBAR UI
# ═══════════════════════════════════════════════════════════════════════════════
TargetDesign_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      style = "padding:12px 16px 4px;color:#ffffff;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      icon("bullseye"),
      " TargetDesign (SRM / MRM / PRM)"
    ),

    # ── Proteotypicity FASTA ─────────────────────────────────────────────────
    tags$div(
      style = "padding:4px 16px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      icon("dna"),
      " Proteotypicity Check"
    ),
    tags$div(
      style = "padding:0 16px;",
      fileInput(
        ns("proteome_fasta"),
        "Background Proteome FASTA",
        accept = c(".fasta", ".fa", ".faa", ".fas"),
        placeholder = "proteome.fasta"
      ),
      checkboxInput(
        ns("equate_il"),
        "Equate isobaric Ile / Leu (I = L)",
        value = TRUE
      ),
      uiOutput(ns("fasta_status"))
    ),

    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),

    # ── In Silico Digestion ──────────────────────────────────────────────────
    tags$div(
      style = "padding:4px 16px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      icon("scissors"),
      " In Silico Digestion"
    ),
    tags$div(
      style = "padding:0 16px;",
      selectInput(
        ns("enzyme"),
        "Enzyme",
        choices = seqprofile_enzyme_choices(),
        selected = "trypsin",
        width = "100%"
      ),
      uiOutput(ns("enzyme_rule")),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'custom'", ns("enzyme")),
        textInput(
          ns("custom_residues"),
          "Cleavage residues (one-letter)",
          placeholder = "e.g. KR",
          width = "100%"
        ),
        radioButtons(
          ns("custom_terminus"),
          "Cleave at the",
          choices = c(
            "C-terminal side" = "C",
            "N-terminal side" = "N"
          )
        ),
        textInput(
          ns("custom_blocked"),
          "No cleavage next to (optional)",
          placeholder = "e.g. P",
          width = "100%"
        )
      ),
      numericInput(
        ns("max_missed"),
        "Max missed cleavages",
        value = 0,
        min = 0,
        max = 3,
        step = 1,
        width = "100%"
      )
    ),

    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),

    # ── Candidate Selection Rules ────────────────────────────────────────────
    tags$div(
      style = "padding:4px 16px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      icon("sliders"),
      " Selection Criteria"
    ),
    tags$div(
      style = "padding:0 16px;",
      fluidRow(
        column(
          6,
          numericInput(
            ns("opt_min_len"),
            "Optimal min aa",
            value = 8,
            min = 5,
            max = 30,
            step = 1
          )
        ),
        column(
          6,
          numericInput(
            ns("opt_max_len"),
            "Optimal max aa",
            value = 22,
            min = 6,
            max = 40,
            step = 1
          )
        )
      ),
      fluidRow(
        column(
          6,
          numericInput(
            ns("mz_min"),
            "Min precursor m/z",
            value = 350,
            min = 100,
            step = 50
          )
        ),
        column(
          6,
          numericInput(
            ns("mz_max"),
            "Max precursor m/z",
            value = 1250,
            min = 200,
            step = 50
          )
        )
      ),
      numericInput(
        ns("min_transitions"),
        "Min y-ions > precursor m/z(2+)",
        value = 3,
        min = 1,
        max = 8,
        step = 1
      ),
      checkboxInput(
        ns("carbamidomethyl"),
        "Carbamidomethyl Cys (+57.0215 Da)",
        value = TRUE
      ),
      checkboxInput(
        ns("allow_alkylated_cys"),
        "Allow internal alkylated Cys (Tier 2)",
        value = TRUE
      ),
      checkboxInput(
        ns("penalize_met"),
        "Penalize Met (oxidation risk)",
        value = TRUE
      ),
      checkboxInput(
        ns("penalize_trp"),
        "Penalize Trp (oxidation risk)",
        value = TRUE
      ),
      checkboxInput(
        ns("oxidation"),
        "Static Met oxidation (+15.9949 Da)",
        value = FALSE
      )
    ),

    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),

    # ── Actions & Downloads ──────────────────────────────────────────────────
    tags$div(
      style = "padding:4px 16px 12px;",
      actionButton(
        ns("evaluate"),
        "Evaluate Candidates",
        icon = icon("bullseye"),
        class = "btn-primary",
        style = "width:100%;font-weight:bold;margin-bottom:8px;"
      ),
      downloadButton(
        ns("download_candidates"),
        tagList(icon("download"), " Candidate Table (.tsv)"),
        class = "dl-btn",
        style = "width:100%;text-align:left;margin-bottom:6px;"
      ),
      downloadButton(
        ns("download_transitions"),
        tagList(icon("download"), " Skyline Transitions (.tsv)"),
        class = "dl-btn",
        style = "width:100%;text-align:left;"
      )
    )
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# BODY UI
# ═══════════════════════════════════════════════════════════════════════════════
TargetDesign_body_ui <- function(id) {
  ns <- NS(id)
  example_choices <- stats::setNames(
    names(TARGETDESIGN_EXAMPLES),
    vapply(TARGETDESIGN_EXAMPLES, function(x) x$label, character(1))
  )

  tagList(
    # ── 1. Sequence Input ────────────────────────────────────────────────────
    fluidRow(
      box(
        title = "Target Protein Sequence",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        textAreaInput(
          ns("sequence"),
          "Protein sequence (raw amino acids or FASTA)",
          rows = 6,
          width = "100%",
          placeholder = "Paste an amino acid sequence, with or without a >header line (whitespace, digits, and line breaks are stripped automatically)"
        ),
        tags$div(
          style = "display:flex;flex-wrap:wrap;align-items:center;gap:6px;margin-bottom:8px;",
          tags$div(
            style = "width:255px;margin-bottom:-15px;",
            selectInput(
              ns("example_select"),
              label = NULL,
              choices = example_choices,
              selected = "apoa1",
              width = "100%"
            )
          ),
          actionButton(
            ns("load_example"),
            "Load Example Sequence",
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
            placeholder = "e.g. P02647"
          ),
          actionButton(
            ns("fetch_uniprot"),
            "Fetch from UniProt",
            icon = icon("cloud-arrow-down"),
            class = "btn btn-default btn-sm"
          ),
          actionButton(
            ns("evaluate_inline"),
            "Evaluate Candidates",
            icon = icon("bullseye"),
            class = "btn-primary btn-sm",
            style = "margin-left:auto;font-weight:bold;"
          )
        ),
        uiOutput(ns("seq_info"))
      )
    ),

    # ── 2. Assay Design Summary & Sequence Map ───────────────────────────────
    fluidRow(
      box(
        title = "Assay Design Summary & Candidate Coverage Map",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        uiOutput(ns("design_summary"))
      )
    ),

    # ── 3. Ranked Surrogate Peptide Candidates ───────────────────────────────
    fluidRow(
      box(
        title = "Surrogate Peptide Candidates (SRM / MRM / PRM)",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        div(
          style = "display:flex;flex-wrap:wrap;align-items:center;justify-content:space-between;gap:8px;margin-bottom:6px;",
          radioButtons(
            ns("tier_filter"),
            label = NULL,
            choices = c(
              "Recommended & Conditional (Tier 1–2)" = "viable",
              "Recommended only (Tier 1)" = "recommended",
              "All digested peptides" = "all"
            ),
            selected = "all",
            inline = TRUE
          ),
          tags$span(
            style = "font-size:12.5px;color:#6c757d;",
            "Click any peptide row to inspect its precursor & product transitions below."
          )
        ),
        DT::dataTableOutput(ns("candidates_table"))
      )
    ),

    # ── 4. Transition Explorer & Guidelines ──────────────────────────────────
    fluidRow(
      box(
        title = "Transition Explorer (Precursor → Product Ions)",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        fluidRow(
          column(
            4,
            radioButtons(
              ns("fragment_charge"),
              "Fragment ion charge",
              choices = c("1+ (standard SRM/PRM)" = "1", "2+" = "2"),
              selected = "1",
              inline = TRUE
            )
          ),
          column(
            4,
            numericInput(
              ns("max_transitions"),
              "Transitions per peptide",
              value = 5,
              min = 3,
              max = 10,
              step = 1,
              width = "200px"
            )
          )
        ),
        uiOutput(ns("selected_peptide_banner")),
        plotly::plotlyOutput(ns("transition_spectrum"), height = 390),
        tags$div(
          style = "margin-top:12px;",
          DT::dataTableOutput(ns("transition_table"))
        ),
        tags$details(
          style = "margin-top:14px;padding:10px 14px;background:#eaf4fb;border-left:4px solid #5FA8D3;border-radius:6px;font-size:13px;",
          tags$summary(
            style = "font-weight:700;cursor:pointer;",
            "Targeted Proteomics Selection Rules & Reference (CPTAC / PMC4830481)"
          ),
          tags$div(
            style = "margin-top:8px;line-height:1.55;",
            tags$ul(
              style = "padding-left:18px;margin-bottom:6px;",
              tags$li(
                tags$strong("Proteotypicity:"),
                " Surrogate peptides must map uniquely to the target protein in the organism's digested proteome (equating isobaric Ile/Leu by default)."
              ),
              tags$li(
                tags$strong("Length & Precursor m/z:"),
                " 8–22 aa optimal (7–25 aa acceptable) with a 2+ or 3+ precursor inside the mass analyzer window (default 350–1250 m/z)."
              ),
              tags$li(
                tags$strong("Complete Proteolysis:"),
                " Zero missed cleavages; avoid dibasic ragged ends (KK, KR, RK, RR), internal KP/RP sites, and protein N-terminal peptides (pos 1–2) subject to initiator-Met clipping or N-acetylation."
              ),
              tags$li(
                tags$strong("Pre-analytical & Chemical Stability:"),
                " Avoid N-terminal Q/E or Cys (spontaneous cyclization), DP bonds (acidic hydrolysis), NG/QG motifs (rapid deamidation), N-X-S/T glycosylation sequons, and oxidation-prone Met/Trp when alternatives exist; internal Cys requires quantitative alkylation."
              ),
              tags$li(
                tags$strong("Chromatography & Transitions:"),
                " Moderate hydrophobicity (GRAVY between -2.0 and +1.0) to prevent void-volume elution or vial/column adsorption, C-terminal K/R for heavy SIS standards, and at least 3 singly charged y-ions (y3+) above the 2+ precursor m/z."
              )
            ),
            tags$p(
              style = "margin:0;color:#495057;",
              tags$strong("Reference: "),
              "Hoofnagle AN, Whiteaker JR, Carr SA, Kuhn E, Liu T, Massoni SA, et al. Recommendations for the Generation, Quantification, Storage, and Handling of Peptides Used for Mass Spectrometry-Based Assays. ",
              tags$em("Clin Chem."),
              " 2016;62(1):48-69. PMCID: PMC4830481."
            )
          )
        )
      )
    )
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# SERVER
# ═══════════════════════════════════════════════════════════════════════════════
TargetDesign_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Example loader, Clear, and UniProt REST fetch ────────────────────────
    observeEvent(input$load_example, {
      key <- input$example_select
      if (is.null(key) || !key %in% names(TARGETDESIGN_EXAMPLES)) {
        key <- "apoa1"
      }
      updateTextAreaInput(
        session,
        "sequence",
        value = TARGETDESIGN_EXAMPLES[[key]]$fasta
      )
    })

    observeEvent(input$clear, {
      updateTextAreaInput(session, "sequence", value = "")
    })

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

    # ── Live sanitized sequence feedback ─────────────────────────────────────
    output$seq_info <- renderUI({
      cs <- clean_protein_sequence(input$sequence)
      parts <- list(tags$strong(sprintf("Length: %d aa", cs$length)))
      if (!is.na(cs$header)) {
        parts <- c(
          parts,
          list(tags$div(style = "color:#6c757d;", cs$header))
        )
      }
      if (cs$n_entries > 1L) {
        parts <- c(
          parts,
          list(tags$div(
            style = "color:#b9770e;",
            sprintf(
              "%d FASTA entries pasted; only the first is evaluated.",
              cs$n_entries
            )
          ))
        )
      }
      if (length(cs$illegal_chars) > 0L) {
        parts <- c(
          parts,
          list(tags$div(
            style = "color:#c0392b;font-weight:600;",
            paste0(
              "Unsupported characters detected: ",
              paste(cs$illegal_chars, collapse = " "),
              ". Remove them before evaluating candidates."
            )
          ))
        )
      } else if (length(cs$invalid_chars) > 0L) {
        parts <- c(
          parts,
          list(tags$div(
            style = "color:#b9770e;font-weight:600;",
            icon("triangle-exclamation", style = "margin-right:4px;"),
            paste0(
              "Non-canonical residues detected (",
              paste(cs$invalid_chars, collapse = ", "),
              "): peptides containing B, Z, X, U, O, * or other non-standard codes are flagged Avoid."
            )
          ))
        )
      }
      tags$div(style = "font-size:13px;", parts)
    })

    # ── Enzyme rule caption ──────────────────────────────────────────────────
    output$enzyme_rule <- renderUI({
      enz <- req(input$enzyme)
      if (identical(enz, "custom")) {
        return(NULL)
      }
      spec <- SEQPROFILE_ENZYMES[[enz]]
      req(spec)
      tags$p(
        style = "font-size:11.5px;color:#adb5bd;margin-top:-6px;",
        paste0(
          seqprofile_enzyme_rule(spec$residues, spec$terminus, spec$blocked_by),
          "."
        )
      )
    })

    # ── Background Proteome FASTA parsing & in silico digestion ─────────────
    proteome_seqs <- reactive({
      if (!isTruthy(input$proteome_fasta)) {
        return(NULL)
      }
      notif_id <- ns("proteome_read_notif")
      showNotification(
        "Reading background proteome FASTA...",
        id = notif_id,
        duration = NULL
      )
      on.exit(removeNotification(notif_id), add = TRUE)
      tryCatch(
        read_fasta_custom(input$proteome_fasta$datapath),
        error = function(e) {
          showNotification(
            paste("Failed to read FASTA:", conditionMessage(e)),
            type = "error"
          )
          NULL
        }
      )
    })

    output$fasta_status <- renderUI({
      seqs <- proteome_seqs()
      if (is.null(seqs)) {
        tags$p(
          style = "font-size:11.5px;color:#adb5bd;margin-top:-4px;",
          "No FASTA loaded — checking intra-sequence uniqueness only."
        )
      } else {
        tags$p(
          style = "font-size:11.5px;color:#82e0aa;margin-top:-4px;",
          icon("circle-check"),
          sprintf(
            " %s proteins loaded for proteotypicity check.",
            format(length(seqs), big.mark = ",")
          )
        )
      }
    })

    # Triggered by either the sidebar button or the inline sequence-box button
    eval_trigger <- reactive({
      n1 <- input$evaluate %||% 0L
      n2 <- input$evaluate_inline %||% 0L
      n1 + n2
    })

    proteome_digest <- reactive({
      seqs <- proteome_seqs()
      if (is.null(seqs) || length(seqs) == 0L) {
        return(NULL)
      }
      enz <- input$enzyme %||% "trypsin"
      spec <- if (identical(enz, "custom")) {
        list(
          residues = input$custom_residues,
          terminus = toupper(as.character(input$custom_terminus %||% "C")[1]),
          blocked_by = input$custom_blocked
        )
      } else {
        SEQPROFILE_ENZYMES[[enz]]
      }
      if (is.null(spec)) {
        return(NULL)
      }
      res_cut <- tryCatch(
        seqprofile_parse_residues(spec$residues),
        error = function(e) character(0)
      )
      res_blk <- tryCatch(
        seqprofile_parse_residues(spec$blocked_by),
        error = function(e) character(0)
      )
      if (length(res_cut) == 0L) {
        return(NULL)
      }
      mc <- suppressWarnings(as.integer(input$max_missed %||% 0L))
      if (is.na(mc) || mc < 0L) {
        mc <- 0L
      }

      notif_id <- ns("proteome_digest_notif")
      showNotification(
        sprintf(
          "Digesting %s FASTA proteins in silico...",
          format(length(seqs), big.mark = ",")
        ),
        id = notif_id,
        duration = NULL
      )
      on.exit(removeNotification(notif_id), add = TRUE)

      digest_proteome_fasta(
        seqs,
        residues = res_cut,
        terminus = spec$terminus,
        blocked_by = res_blk,
        max_missed = mc,
        min_len = 5L,
        max_len = 50L
      )
    })

    # ── Candidate evaluation (runs on Evaluate Candidates) ───────────────────
    design_result <- eventReactive(eval_trigger(), {
      req(eval_trigger() > 0L)
      cs <- clean_protein_sequence(input$sequence)
      shiny::validate(
        shiny::need(
          nzchar(cs$sequence),
          "Paste a sequence (or load an example / fetch from UniProt) and click Evaluate Candidates."
        ),
        shiny::need(
          length(cs$illegal_chars) == 0L,
          paste0(
            "The sequence contains unsupported characters: ",
            paste(cs$illegal_chars, collapse = " ")
          )
        )
      )
      if (cs$n_entries > 1L) {
        showNotification(
          sprintf(
            "%d FASTA entries pasted; only the first is evaluated.",
            cs$n_entries
          ),
          type = "warning"
        )
      }

      opt_min <- max(5L, as.integer(input$opt_min_len %||% 8L))
      opt_max <- max(opt_min, as.integer(input$opt_max_len %||% 22L))
      mz_lo <- as.numeric(input$mz_min %||% 350)
      mz_hi <- as.numeric(input$mz_max %||% 1250)
      min_tr <- max(1L, as.integer(input$min_transitions %||% 3L))
      max_tr <- max(min_tr, as.integer(input$max_transitions %||% 5L))

      dig <- tryCatch(
        run_seq_digest(
          cs$sequence,
          input$enzyme %||% "trypsin",
          list(
            residues = input$custom_residues,
            terminus = input$custom_terminus,
            blocked_by = input$custom_blocked,
            max_missed = input$max_missed %||% 0L,
            min_len = opt_min,
            max_len = opt_max,
            mz_min = mz_lo,
            mz_max = mz_hi,
            carbamidomethyl = isTRUE(input$carbamidomethyl),
            oxidation = isTRUE(input$oxidation)
          )
        ),
        seqprofile_error = function(e) e
      )
      if (inherits(dig, "condition")) {
        shiny::validate(shiny::need(FALSE, conditionMessage(dig)))
      }

      prot_db <- proteome_digest()
      cands <- evaluate_targeted_candidates(
        dig,
        digest_df = prot_db,
        opts = list(
          opt_min_len = opt_min,
          opt_max_len = opt_max,
          hard_min_len = min(7L, opt_min),
          hard_max_len = max(25L, opt_max),
          mz_min = mz_lo,
          mz_max = mz_hi,
          min_transitions = min_tr,
          max_transitions = max_tr,
          allow_alkylated_cys = isTRUE(input$allow_alkylated_cys),
          penalize_met = isTRUE(input$penalize_met),
          penalize_trp = isTRUE(input$penalize_trp),
          equate_il = isTRUE(input$equate_il)
        )
      )

      prot_name <- if (!is.na(cs$header) && nzchar(cs$header)) {
        fasta_header_ids(cs$header)
      } else {
        "TargetProtein"
      }

      list(
        header = cs$header,
        protein_name = prot_name,
        digest = dig,
        candidates = cands,
        tier_state = attr(cands, "tier_residue_state"),
        has_fasta = !is.null(prot_db),
        n_fasta_proteins = if (is.null(proteome_seqs())) 0L else length(proteome_seqs())
      )
    })

    # ── Summary Strip & Sequence Candidate Map ───────────────────────────────
    output$design_summary <- renderUI({
      res <- design_result()
      cands <- res$candidates
      n_rec <- sum(cands$tier == "Recommended")
      n_cond <- sum(cands$tier == "Conditional")
      n_avoid <- sum(cands$tier == "Avoid")
      best_score <- if (nrow(cands) > 0L) max(cands$score) else 0L

      tagList(
        tags$div(
          style = "display:flex;flex-wrap:wrap;gap:8px;margin-bottom:10px;",
          targetdesign_stat("Enzyme", res$digest$enzyme_name),
          targetdesign_stat("Digested Peptides", nrow(cands)),
          targetdesign_stat("Tier 1 · Recommended", n_rec, "#1e8449"),
          targetdesign_stat("Tier 2 · Conditional", n_cond, "#b9770e"),
          targetdesign_stat("Tier 3 · Avoid", n_avoid, "#922b21"),
          targetdesign_stat("Best Score", paste0(best_score, " / 100")),
          targetdesign_stat(
            "Proteotypicity",
            if (res$has_fasta) {
              sprintf("FASTA (%s)", format(res$n_fasta_proteins, big.mark = ","))
            } else {
              "Intra-sequence"
            }
          )
        ),
        tags$div(
          style = "font-size:12.5px;color:#6c757d;margin-bottom:8px;",
          tags$span(
            style = "background:#bfe3c6;color:#17202a;padding:1px 6px;border-radius:3px;font-weight:600;",
            "Green"
          ),
          " = covered by Tier 1 (Recommended) peptide; ",
          tags$span(
            style = "background:#fdebd0;color:#7d5a00;padding:1px 6px;border-radius:3px;font-weight:600;",
            "Amber"
          ),
          " = covered by Tier 2 (Conditional) peptide; ",
          tags$span(style = "color:#95a5a6;font-weight:600;", "Grey"),
          " = no recommended/conditional surrogate peptide.",
          if (!res$has_fasta) {
            tags$span(
              style = "color:#b9770e;margin-left:6px;",
              "Upload a background proteome FASTA in the sidebar to verify proteome-wide proteotypicity."
            )
          }
        ),
        targetdesign_sequence_map(res$digest$sequence, res$tier_state)
      )
    })

    # ── Filtered Candidate Table ─────────────────────────────────────────────
    filtered_candidates <- reactive({
      cands <- design_result()$candidates
      mode <- input$tier_filter %||% "all"
      if (identical(mode, "recommended")) {
        cands[cands$tier == "Recommended", , drop = FALSE]
      } else if (identical(mode, "viable")) {
        cands[cands$tier %in% c("Recommended", "Conditional"), , drop = FALSE]
      } else {
        cands
      }
    })

    candidates_display_df <- reactive({
      d <- filtered_candidates()
      data.frame(
        Rank = d$rank,
        Tier = d$tier,
        Score = d$score,
        Flanked_Peptide = d$flanked_sequence,
        Sequence = d$sequence,
        Start = d$start,
        End = d$end,
        Length = d$length,
        `m/z (2+)` = round(d$mz_2, 4),
        `m/z (3+)` = round(d$mz_3, 4),
        GRAVY = round(d$gravy, 3),
        Proteotypicity = d$proteotypicity,
        Mapped_Proteins = ifelse(is.na(d$mapped_proteins), "—", d$mapped_proteins),
        `y > m/z(2+)` = d$n_y_above_mz2,
        Recommended_Transitions = ifelse(
          nzchar(d$top_transitions),
          d$top_transitions,
          "—"
        ),
        Flags_and_Rationale = d$flags,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })

    output$candidates_table <- DT::renderDataTable(
      {
        df <- candidates_display_df()
        dt <- DT::datatable(
          df,
          options = list(
            scrollX = TRUE,
            pageLength = 15,
            dom = "frtip",
            language = list(
              emptyTable = "No peptides match the selected tier filter (switch to 'All digested peptides' to inspect penalized sequences)."
            )
          ),
          rownames = FALSE,
          selection = list(
            mode = "single",
            selected = if (nrow(df) > 0L) 1L else integer(0)
          ),
          class = "cell-border stripe hover"
        )
        if (nrow(df) > 0L) {
          dt <- dt |>
            DT::formatStyle(
              "Tier",
              backgroundColor = DT::styleEqual(
                c("Recommended", "Conditional", "Avoid"),
                c("#d4efdf", "#fcf3cf", "#fadbd8")
              ),
              color = DT::styleEqual(
                c("Recommended", "Conditional", "Avoid"),
                c("#1e8449", "#7d6608", "#922b21")
              ),
              fontWeight = "bold"
            )
        }
        dt
      },
      server = FALSE
    )

    # ── Selected Peptide Transitions & Spectrum ──────────────────────────────
    selected_candidate <- reactive({
      d <- filtered_candidates()
      shiny::validate(shiny::need(
        nrow(d) > 0L,
        "No peptides in the current view. Switch the filter above to inspect candidates."
      ))
      sel <- input$candidates_table_rows_selected
      idx <- if (length(sel) == 1L && sel >= 1L && sel <= nrow(d)) sel else 1L
      d[idx, , drop = FALSE]
    })

    selected_transitions <- reactive({
      cand <- selected_candidate()
      res <- design_result()
      shiny::validate(
        shiny::need(
          cand$length >= 2L && !is.na(cand$mono_mass),
          "Selected peptide is too short or contains non-canonical residues; fragment transitions cannot be computed."
        )
      )
      z <- suppressWarnings(as.integer(input$fragment_charge %||% 1L))
      if (is.na(z) || z < 1L) {
        z <- 1L
      }
      max_tr <- max(1L, as.integer(input$max_transitions %||% 5L))
      ions <- select_targeted_transitions(
        cand$sequence,
        carbamidomethyl = res$digest$params$carbamidomethyl,
        oxidation = res$digest$params$oxidation,
        charge = z,
        max_transitions = max_tr
      )
      attr(ions, "candidate") <- cand
      ions
    })

    output$selected_peptide_banner <- renderUI({
      cand <- selected_candidate()
      badge_bg <- switch(
        cand$tier,
        Recommended = "#d4efdf",
        Conditional = "#fcf3cf",
        "#fadbd8"
      )
      badge_fg <- switch(
        cand$tier,
        Recommended = "#1e8449",
        Conditional = "#7d6608",
        "#922b21"
      )
      tags$div(
        style = "padding:8px 12px;margin-bottom:8px;background:#f8f9fa;border-radius:6px;font-size:13px;",
        tags$span(
          style = sprintf(
            "background:%s;color:%s;padding:2px 8px;border-radius:4px;font-weight:700;margin-right:8px;",
            badge_bg,
            badge_fg
          ),
          sprintf("%s (Score %d)", cand$tier, cand$score)
        ),
        tags$strong(style = "font-family:monospace;font-size:14px;", cand$flanked_sequence),
        sprintf(
          "  |  pos %d–%d (%d aa)  |  [M+2H]2+ = %.4f m/z  |  [M+3H]3+ = %.4f m/z  |  %s",
          cand$start,
          cand$end,
          cand$length,
          cand$mz_2,
          cand$mz_3,
          cand$flags
        )
      )
    })

    output$transition_spectrum <- plotly::renderPlotly({
      ions <- selected_transitions()
      cand <- attr(ions, "candidate")
      p <- plot_targeted_spectrum(
        ions,
        mz_2 = attr(ions, "mz_2"),
        title = sprintf(
          "%s  (dashed line = [M+2H]2+ %.2f m/z)",
          cand$sequence,
          cand$mz_2
        )
      )
      plotly::ggplotly(p, tooltip = "text") |>
        plotly::config(displaylogo = FALSE)
    })

    output$transition_table <- DT::renderDataTable(
      {
        ions <- selected_transitions()
        ord <- order(!ions$recommended, ions$transition_rank, ions$type, -ions$number)
        ions <- ions[ord, , drop = FALSE]
        df <- data.frame(
          Ion = ions$ion,
          Series = ions$type,
          Fragment_Sequence = ions$fragment,
          `Product m/z` = round(ions$mz, 4),
          `Above [M+2H]2+` = ifelse(ions$above_precursor_2, "yes", "no"),
          `Proline Cleavage` = ifelse(ions$proline_enhanced, "yes", "no"),
          Recommended = ifelse(
            ions$recommended,
            paste0("Transition #", ions$transition_rank),
            "—"
          ),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        DT::datatable(
          df,
          options = list(scrollX = TRUE, pageLength = 10, dom = "tip"),
          rownames = FALSE,
          selection = "none",
          class = "cell-border stripe compact"
        )
      },
      server = FALSE
    )

    # ── Downloads ────────────────────────────────────────────────────────────
    output$download_candidates <- downloadHandler(
      filename = function() {
        nm <- tryCatch(design_result()$protein_name, error = function(e) "target")
        nm <- gsub("[^A-Za-z0-9._-]", "_", nm)
        paste0("TargetDesign_candidates_", nm, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        data.table::fwrite(
          candidates_display_df(),
          file = file,
          sep = "\t",
          na = "NA"
        )
      }
    )

    output$download_transitions <- downloadHandler(
      filename = function() {
        nm <- tryCatch(design_result()$protein_name, error = function(e) "target")
        nm <- gsub("[^A-Za-z0-9._-]", "_", nm)
        paste0("TargetDesign_Skyline_transitions_", nm, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        res <- design_result()
        max_tr <- max(1L, as.integer(input$max_transitions %||% 5L))
        tr_df <- build_targeted_transition_list(
          res$candidates,
          protein_name = res$protein_name,
          carbamidomethyl = res$digest$params$carbamidomethyl,
          oxidation = res$digest$params$oxidation,
          tiers_keep = c("Recommended", "Conditional"),
          max_transitions = max_tr,
          mz_min = res$digest$params$mz_min,
          mz_max = res$digest$params$mz_max
        )
        data.table::fwrite(tr_df, file = file, sep = "\t", na = "NA")
      }
    )
  })
}
