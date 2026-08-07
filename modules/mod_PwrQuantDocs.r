# ── PwrQuant Documentation ──────────────────────────────────────────────────
# A static, self-contained user guide for mod_PwrQuant.r: what each control
# does, how to choose it, what the statistics underneath actually claim, and
# how to read the output. UI-only; there is no server component.

pwrquant_docs_css <- "
  .pq-doc { max-width: 1080px; margin: 0 auto; color: #2c3e50; font-size: 14px;
            line-height: 1.65; }
  .pq-doc h3 { font-size: 1.5rem; font-weight: 700; color: #1B4965;
               margin: 26px 0 10px; padding-bottom: 6px;
               border-bottom: 2px solid #e9ecef; }
  .pq-doc h4 { font-size: 1.12rem; font-weight: 700; color: #17202a;
               margin: 20px 0 8px; }
  .pq-doc p  { margin-bottom: 10px; }
  .pq-doc ul, .pq-doc ol { padding-left: 22px; margin-bottom: 12px; }
  .pq-doc li { margin-bottom: 5px; }
  .pq-doc code { background: #f1f3f5; color: #a03030; padding: 1px 5px;
                 border-radius: 4px; font-size: 92%; }
  .pq-doc table { width: 100%; border-collapse: collapse; margin: 12px 0 18px; }
  .pq-doc th { background: #1B4965; color: #fff; text-align: left;
               padding: 8px 10px; font-size: 13px; font-weight: 600; }
  .pq-doc td { padding: 7px 10px; border-bottom: 1px solid #e9ecef;
               vertical-align: top; font-size: 13px; }
  .pq-doc tr:nth-child(even) td { background: #f8f9fa; }
  .pq-callout { border-left: 4px solid #5FA8D3; background: #f0f7fb;
                padding: 10px 14px; margin: 14px 0; border-radius: 0 6px 6px 0; }
  .pq-warn    { border-left: 4px solid #e67e22; background: #fdf3e7;
                padding: 10px 14px; margin: 14px 0; border-radius: 0 6px 6px 0; }
  .pq-formula { background: #f8f9fa; border: 1px solid #e9ecef; border-radius: 6px;
                padding: 10px 14px; margin: 12px 0; font-family: monospace;
                font-size: 13px; color: #17202a; overflow-x: auto; }
  .pq-toc { background: #f8f9fa; border: 1px solid #e9ecef; border-radius: 8px;
            padding: 14px 20px; margin-bottom: 18px; }
  .pq-toc ul { margin-bottom: 0; }
"

pq_note <- function(...) tags$div(class = "pq-callout", ...)
pq_warn <- function(...) tags$div(class = "pq-warn", ...)
pq_eq <- function(...) tags$div(class = "pq-formula", ...)

pq_table <- function(headers, rows) {
  tags$table(
    tags$thead(tags$tr(lapply(headers, tags$th))),
    tags$tbody(lapply(rows, function(r) {
      tags$tr(lapply(r, function(cell) {
        # cells may be raw HTML strings or already-built shiny tags
        tags$td(if (is.character(cell)) HTML(cell) else cell)
      }))
    }))
  )
}

PwrQuantDocs_body_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tags$head(tags$style(HTML(pwrquant_docs_css))),
    fluidRow(box(
      title = "PwrQuant \u2014 user guide and statistical rationale",
      status = "primary",
      solidHeader = TRUE,
      width = 12,
      tabsetPanel(
        id = ns("doc_tabs"),
        type = "pills",

        # ── 1. Overview ─────────────────────────────────────────────────────
        tabPanel(
          "Overview",
          div(
            class = "pq-doc",
            tags$h3("What PwrQuant does"),
            tags$p(
              "PwrQuant takes a protein abundance matrix, fits a linear model ",
              "with ",
              tags$code("limma"),
              ", and reports differentially ",
              "abundant proteins together with an honest assessment of how ",
              "small an effect the experiment could actually have detected. ",
              "It then layers on functional interpretation (GO over-",
              "representation, STRING networks) and a set of diagnostic views."
            ),
            tags$p(
              "The module is deliberately split into a ",
              tags$b("discovery"),
              " path (what is different?) and a ",
              tags$b("sensitivity"),
              " path (what could we have found?). Those two questions are ",
              "answered by different machinery and should not be conflated \u2014 ",
              "see the ",
              tags$b("Sensitivity & power"),
              " tab for why."
            ),

            tags$h3("Input format"),
            tags$p(
              "Upload a ",
              tags$code(".tsv"),
              ", ",
              tags$code(".csv"),
              " or ",
              tags$code(".txt"),
              " file with:"
            ),
            tags$ul(
              tags$li(
                tags$b("First column"),
                " \u2014 protein identifiers. UniProt ",
                "accessions, ",
                tags$code("sp|P12345|NAME_HUMAN"),
                " strings, ",
                "or gene symbols all work; the module extracts a lookup ID for ",
                "enrichment automatically. Semicolon-separated protein groups ",
                "are resolved to their first member."
              ),
              tags$li(
                tags$b("Remaining columns"),
                " \u2014 one per sample, holding ",
                tags$b("raw (non-log) intensities"),
                ". Missing values may be ",
                "blank, ",
                tags$code("NA"),
                ", or ",
                tags$code("0"),
                "."
              )
            ),
            pq_warn(
              tags$b("Do not pre-log your data."),
              " The module applies ",
              tags$code("log2(x + 1)"),
              " itself. Uploading log-space values ",
              "will compress your fold changes into nonsense."
            ),

            tags$h3("Minimum workflow"),
            tags$ol(
              tags$li(
                tags$b("Upload"),
                " the matrix in the sidebar."
              ),
              tags$li(
                tags$b("Map metadata"),
                " on the ",
                tags$i("Metadata Mapping"),
                " tab. Sample names are parsed into conditions automatically, ",
                "but the table is editable \u2014 check it. ",
                tags$code("Condition"),
                " drives the design matrix and ",
                tags$code("Batch"),
                " drives ",
                "ComBat, so errors here propagate everywhere."
              ),
              tags$li(
                tags$b("Inspect sparsity and QC"),
                " before modelling. If the ",
                tags$i("Sparsity"),
                " tab shows one condition is mostly empty, ",
                "no amount of statistics will rescue it."
              ),
              tags$li(
                tags$b("Choose contrasts"),
                " in the sidebar. Contrasts are ",
                "written ",
                tags$code("GroupA-GroupB"),
                ", and the sign of the ",
                "reported log fold change follows that order: positive means ",
                "higher in ",
                tags$code("GroupA"),
                "."
              ),
              tags$li(
                tags$b("Press Start limma."),
                " Everything downstream ",
                "(enrichment, networks, power) depends on this run."
              ),
              tags$li(
                tags$b("Press Run MDD simulation"),
                " to learn what effect size ",
                "your design can resolve."
              )
            ),

            tags$h3("Tab map"),
            pq_table(
              c("Tab", "Purpose", "Needs"),
              list(
                c(
                  "Metadata Mapping",
                  "Assign condition and batch to each sample",
                  "Upload"
                ),
                c(
                  "Sparsity",
                  "Missing-value pattern across the whole matrix",
                  "Upload"
                ),
                c(
                  "Pre-processing QC",
                  "CV, mean&ndash;variance, distributions, PCA, PLS-DA",
                  "Upload / limma"
                ),
                c(
                  "Differential Abundance",
                  "MA, volcano, top hits, p-value histograms",
                  "limma"
                ),
                c("Correlation", "Agreement between two contrasts", "limma"),
                c(
                  "Power Statistics",
                  "Prospective MDD, replicate sweep, conditional sensitivity",
                  "limma (+ MDD run)"
                ),
                c(
                  "Enrichment",
                  "GO over-representation per contrast and direction",
                  "limma + ORA"
                ),
                c(
                  "Interaction Network",
                  "STRING protein&ndash;protein network",
                  "limma + STRING"
                ),
                c(
                  "Selected Proteins",
                  "Per-protein abundance by condition",
                  "Upload"
                ),
                c("UpSet", "Detection overlap between conditions", "Upload"),
                c(
                  "Z-Score Heatmap",
                  "Clustered abundance heatmap",
                  "Upload / limma"
                )
              )
            )
          )
        ),

        # ── 2. Preprocessing ────────────────────────────────────────────────
        tabPanel(
          "Preprocessing",
          div(
            class = "pq-doc",
            tags$h3("Pipeline order"),
            tags$p("Pressing ", tags$b("Start limma"), " runs, in this order:"),
            pq_eq(
              "log2(x + 1)  \u2192  min-valid filter  \u2192  imputation (robust mode only)",
              tags$br(),
              "  \u2192  ComBat  \u2192  normalizeBetweenArrays  \u2192  lmFit  \u2192  contrasts.fit  \u2192  eBayes"
            ),
            pq_warn(
              tags$b("Ordering caveat."),
              " Batch correction runs ",
              tags$i("before"),
              " between-array normalization. That is the ",
              "reverse of the more common convention. If your batches differ ",
              "mainly in overall loading rather than in composition, consider ",
              "setting normalization to ",
              tags$code("quantile"),
              " and ",
              "checking the two PCA panels: batch structure should be visibly ",
              "reduced in the post-processed plot. If it is not, the ",
              "correction did not work and downstream p-values are optimistic."
            ),

            tags$h3("Min. valid values per group (%)"),
            tags$p(
              "Drops any protein that does not have at least this percentage of ",
              "non-missing values in ",
              tags$b("every"),
              " group. This is the ",
              "single most consequential preprocessing choice, because it ",
              "determines how much imputation has to invent."
            ),
            pq_table(
              c("Setting", "Effect", "Use when"),
              list(
                c(
                  "<code>0</code> (default)",
                  "Keeps everything; maximum reliance on imputation",
                  "Exploratory first pass"
                ),
                c(
                  "<code>50</code>",
                  "Half the replicates must be observed per group",
                  "Typical, defensible default for DDA"
                ),
                c(
                  "<code>67&ndash;100</code>",
                  "Near-complete cases only",
                  "When you want conclusions that do not depend on imputation"
                )
              )
            ),
            pq_note(
              "Raising this filter shrinks the number of tests, which ",
              tags$i("raises"),
              " the BH threshold and can increase the number ",
              "of significant hits even though you removed data. That is ",
              "expected behaviour, not a bug."
            ),

            tags$h3("Imputation"),
            pq_warn(
              tags$b(
                "Imputation only runs when Limma Regression Method is set "
              ),
              tags$b(tags$code("robust")),
              ". In ",
              tags$code("ls"),
              " mode the ",
              "matrix keeps its ",
              tags$code("NA"),
              "s and limma fits each ",
              "protein on whatever samples it has. Both are legitimate; they ",
              "are just different analyses, and the choice is currently tied ",
              "to the regression method rather than exposed separately."
            ),
            pq_table(
              c("Method", "Missingness assumption", "Behaviour", "Cost"),
              list(
                c(
                  "<code>knn</code>",
                  "MAR",
                  "Borrows from the 5 most similar proteins within the same group. Shrinks variance, so p-values skew optimistic.",
                  "Fast"
                ),
                c(
                  "<code>minprob</code>",
                  "MNAR",
                  "Draws from the protein's own left tail (mean &minus; 1.8 SD). Appropriate for below-detection dropout, but can manufacture apparent differences.",
                  "Fastest"
                ),
                c(
                  "<code>missforest</code>",
                  "MAR",
                  "Random-forest imputation. Most accurate on structured data.",
                  "Very slow"
                ),
                c(
                  "<code>bpca</code>",
                  "MAR",
                  "Global latent-factor model across all samples \u2014 the only method not run per group.",
                  "Moderate"
                )
              )
            ),
            tags$p(
              "Values still missing after imputation (typically a protein absent ",
              "from an entire group) are filled from ",
              tags$code("N(global min \u2212 1, 0.3)"),
              ". Those proteins are ",
              "flagged and forced to ",
              tags$i("Not significant"),
              " \u2014 see ",
              tags$b("The imputation-driven mask"),
              " on the ",
              tags$i("Linear model"),
              " tab."
            ),
            pq_note(
              "Groupwise imputation (all methods except BPCA) fills each ",
              "condition using only its own samples. This avoids leaking ",
              "signal across the very contrast you are testing, but it also ",
              "means an all-missing group cannot be recovered."
            ),

            tags$h3("Normalization"),
            pq_table(
              c("Method", "What it equalises", "Assumption"),
              list(
                c(
                  "<code>none</code>",
                  "Nothing",
                  "Data are already normalized upstream"
                ),
                c(
                  "<code>scale</code>",
                  "Median intensity per sample",
                  "Differences are pure loading; distribution shapes match"
                ),
                c(
                  "<code>quantile</code>",
                  "The entire distribution",
                  "All samples share one true distribution \u2014 strong, but robust"
                ),
                c(
                  "<code>cyclicloess</code> (default)",
                  "Intensity-dependent bias, pairwise",
                  "Bias varies smoothly with abundance; slowest but gentlest"
                )
              )
            ),

            tags$h3("Batch correction"),
            tags$p(
              "ComBat runs automatically when the metadata table contains more ",
              "than one batch, using ",
              tags$code("par.prior = TRUE"),
              " and no ",
              "covariate model. With a single batch it is skipped and you get a ",
              "notification."
            ),
            pq_warn(
              "ComBat is run with ",
              tags$code("mod = NULL"),
              ", meaning the ",
              "biological design is ",
              tags$i("not"),
              " protected during ",
              "correction. If batch is confounded with condition \u2014 for example ",
              "all treated samples processed on day 2 \u2014 ComBat will remove your ",
              "biological effect along with the batch effect. Check the PCA ",
              "panels and, if confounded, set every sample to one batch and ",
              "report the limitation instead."
            ),

            tags$h3("Reading the QC tab"),
            tags$ul(
              tags$li(
                tags$b("CV distributions"),
                " \u2014 median CV per condition. ",
                "Sharply higher CV in one group usually means a technical ",
                "problem rather than biology."
              ),
              tags$li(
                tags$b("Mean\u2013variance"),
                " \u2014 a visible downward trend is ",
                "normal for proteomics and is exactly what ",
                tags$code("eBayes Trend = TRUE"),
                " exists to model. A flat ",
                "cloud means you can safely turn the trend off."
              ),
              tags$li(
                tags$b("Raw vs normalized boxplots"),
                " \u2014 after ",
                "normalization the medians should line up. If they do not, the ",
                "chosen method is too weak."
              ),
              tags$li(
                tags$b("PCA raw vs post-processed"),
                " \u2014 the honest test of ",
                "whether batch correction worked. Samples should regroup by ",
                "condition, not by batch."
              ),
              tags$li(
                tags$b("PLS-DA"),
                " \u2014 supervised, so it will ",
                tags$i("always"),
                " separate your groups, including on pure ",
                "noise. Use it for visualisation only; never as evidence that ",
                "groups differ."
              )
            )
          )
        ),

        # ── 3. Linear model ─────────────────────────────────────────────────
        tabPanel(
          "Linear model",
          div(
            class = "pq-doc",
            tags$h3("Design and contrasts"),
            tags$p(
              "The design is a means model, ",
              tags$code("~ 0 + condition"),
              ", ",
              "so each coefficient is a group mean and a contrast ",
              tags$code("A-B"),
              " is the difference of two means in log2 space. ",
              "A logFC of 1 is a two-fold increase in A relative to B."
            ),

            tags$h3("Regression method"),
            pq_table(
              c("Setting", "Meaning", "Trade-off"),
              list(
                c(
                  "<code>ls</code> (default)",
                  "Ordinary least squares; NAs kept and handled per protein",
                  "Fastest and unbiased, but a single outlier replicate can drive a hit"
                ),
                c(
                  "<code>robust</code>",
                  "M-estimation, plus imputation of the matrix first",
                  "Downweights outliers, but the imputation step introduces its own assumptions"
                )
              )
            ),

            tags$h3("Empirical Bayes moderation"),
            tags$p(
              "Proteomics experiments have few replicates, so a per-protein ",
              "variance estimate is noisy \u2014 and a protein that happens to get a ",
              "small variance by chance produces a spuriously large t statistic. ",
              tags$code("eBayes()"),
              " fixes this by shrinking each protein's ",
              "variance toward a global prior estimated across all proteins:"
            ),
            pq_eq(
              "s\u00B2_post,g  =  (d\u2080 \u00B7 s\u2080\u00B2  +  d_g \u00B7 s\u00B2_g)  /  (d\u2080 + d_g)"
            ),
            tags$p(
              "where ",
              tags$code("s\u2080\u00B2"),
              " and ",
              tags$code("d\u2080"),
              " are the prior variance and prior degrees of freedom fitted from ",
              "the data, and ",
              tags$code("d_g"),
              " is the protein's residual df. ",
              "The moderated t statistic then uses ",
              tags$code("s_post"),
              " in the denominator and is compared against a t distribution with"
            ),
            pq_eq("df_total  =  d_g  +  d\u2080"),
            pq_note(
              "Those extra ",
              tags$code("d\u2080"),
              " degrees of freedom are real ",
              "and are the reason limma outperforms a per-protein t-test at low ",
              "n. They are also the reason a classical two-sample power ",
              "calculation does not describe this pipeline \u2014 the point taken up ",
              "on the ",
              tags$b("Sensitivity & power"),
              " tab."
            ),

            tags$h4("eBayes Trend"),
            tags$p(
              tags$code("TRUE"),
              " (default) makes the prior variance a smooth ",
              "function of mean abundance instead of a single constant. Use it ",
              "whenever the mean\u2013variance panel shows a trend, which for ",
              "label-free proteomics is almost always. Set ",
              tags$code("FALSE"),
              " only if that panel is flat."
            ),
            pq_note(
              "With the trend enabled, proteins whose sigma or mean abundance ",
              "cannot be estimated are dropped before ",
              tags$code("eBayes()"),
              " to prevent the loess fit from failing. They will be absent from ",
              "the results table."
            ),

            tags$h3("Significance"),
            tags$p("A protein is called Increased or Decreased when all of:"),
            tags$ul(
              tags$li(
                tags$code("adj.P.Val \u2264 Significance FDR"),
                " \u2014 Benjamini\u2013",
                "Hochberg adjusted, computed separately within each contrast"
              ),
              tags$li(
                tags$code("|logFC| \u2265 Minimum |log\u2082FC|"),
                " \u2014 an explicit ",
                "effect-size floor, default 0 (off)"
              ),
              tags$li("the protein is not flagged as imputation-driven")
            ),
            pq_warn(
              tags$b("Two thresholds, both yours."),
              " Earlier versions ",
              "silently ANDed a power-derived cutoff into this rule. That has ",
              "been removed: the only filters are the two you set. If you want ",
              "an effect-size floor, set ",
              tags$code("Minimum |log\u2082FC|"),
              " deliberately \u2014 a common choice is 0.58 (1.5-fold) or 1.0 ",
              "(2-fold) \u2014 and report it."
            ),

            tags$h4("The imputation-driven mask"),
            tags$p(
              "Before modelling, the module records which proteins are entirely ",
              "missing within each condition. If a protein is fully absent from ",
              "either side of a contrast, its apparent fold change is a product ",
              "of the imputation rule rather than of measurement, so it is ",
              "forced to ",
              tags$i("Not significant"),
              " regardless of its ",
              "p-value. Such proteins are often genuinely interesting \u2014 an ",
              "on/off presence pattern \u2014 but they are a qualitative finding, ",
              "not a quantitative one. Inspect them on the ",
              tags$i("UpSet"),
              " tab."
            ),

            tags$h3("Reading the diagnostic plots"),
            tags$h4("P-value histogram"),
            tags$p(
              "This is the fastest way to tell whether the model is sane."
            ),
            pq_table(
              c("Shape", "Interpretation"),
              list(
                c(
                  "Flat, with a spike near 0",
                  "Healthy. The flat part is the null; the spike is signal."
                ),
                c(
                  "Completely flat",
                  "No detectable differential abundance. Not a failure &mdash; a result."
                ),
                c(
                  "Peak near 1 (right-skewed)",
                  "Variance is overestimated. Often over-aggressive imputation or a mis-specified batch."
                ),
                c(
                  "Bimodal or lumpy",
                  "Unmodelled structure &mdash; usually a batch or covariate missing from the design."
                )
              )
            ),
            tags$h4("MA and volcano plots"),
            tags$p(
              "The MA plot shows logFC against mean abundance; a fold change ",
              "that drifts systematically with abundance points at incomplete ",
              "normalization. The volcano plot shows logFC against ",
              tags$code("-log10(adj.P.Val)"),
              ", so the y-axis is already ",
              "multiplicity-corrected."
            ),
            tags$h4("Contrast correlation"),
            tags$p(
              "Compares two contrasts protein by protein and reports Spearman ",
              "\u03C1 within each class: ",
              tags$b("Concordant"),
              " (same ",
              "direction in both), ",
              tags$b("Inverse"),
              " (opposite ",
              "directions), ",
              tags$b("Mismatch"),
              " (significant in only one). ",
              "Spearman is used rather than Pearson because the logFC ",
              "distribution is heavy-tailed."
            )
          )
        ),

        # ── 4. Sensitivity & power ──────────────────────────────────────────
        tabPanel(
          "Sensitivity & power",
          div(
            class = "pq-doc",
            div(
              class = "pq-toc",
              tags$b("On this page"),
              tags$ul(
                tags$li("Why a textbook power calculation does not apply here"),
                tags$li("Conditional MDD \u2014 the per-protein diagnostic"),
                tags$li(
                  "Prospective MDD \u2014 the simulation, and why it is valid"
                ),
                tags$li("The missingness arm"),
                tags$li("The replicate sweep"),
                tags$li("Choosing the parameters"),
                tags$li("Reading the outputs")
              )
            ),

            tags$h3("Why a textbook power calculation does not apply"),
            tags$p(
              "The obvious thing to do is call ",
              tags$code("power.t.test()"),
              " at \u03B1 = 0.05, convert Cohen's d into a log fold change using ",
              "the observed variance, and call the result a minimum detectable ",
              "difference. It is wrong in ",
              "three independent ways."
            ),
            tags$ol(
              tags$li(
                tags$b("The alpha is not the operating threshold."),
                " ",
                "Significance is called on BH-adjusted p-values across ",
                "thousands of proteins. At the BH boundary the largest rejected ",
                "raw p-value is ",
                tags$code("FDR \u00D7 R / m"),
                " \u2014 with 5,000 ",
                "proteins and 1,000 rejections that is 0.01, and with few ",
                "rejections it is far smaller still. Using 0.05 makes the ",
                "calculation anti-conservative."
              ),
              tags$li(
                tags$b("The degrees of freedom are wrong."),
                " A two-sample ",
                "t-test has ",
                tags$code("2n \u2212 2"),
                " df. The moderated ",
                "statistic actually used has ",
                tags$code("2n \u2212 2 + d\u2080"),
                ". ",
                "Pairing a shrunken variance with unshrunken df is internally ",
                "inconsistent, and conservative."
              ),
              tags$li(
                tags$b("It is circular."),
                " Comparing each protein's MDD ",
                "against that same protein's observed fold change is the ",
                "observed-power fallacy. Worse, it is algebraically vacuous: ",
                "since ",
                tags$code("MDD_g = d \u00B7 s_g"),
                " and ",
                tags$code("t_g = logFC_g / (s_g \u00B7 \u221A(2/n))"),
                ", the test ",
                tags$code("|logFC| \u2265 MDD"),
                " reduces to ",
                tags$code("|t| \u2265 d\u00B7\u221A(n/2)"),
                " \u2014 a fixed cutoff on the t ",
                "statistic, with the variance cancelling out entirely. It adds ",
                "no power information; it is just a second, hidden significance ",
                "threshold."
              )
            ),
            pq_note(
              "Note that errors (1) and (2) push in opposite directions, so the ",
              "resulting number can look plausible while being unjustifiable. ",
              "The magnitude isn't the real problem \u2014 the interpretation ",
              "is."
            ),
            tags$p(
              "PwrQuant therefore reports two separate quantities, and labels ",
              "them differently on purpose."
            ),

            tags$h3("Conditional MDD \u2014 the per-protein diagnostic"),
            tags$p(
              "Shown in the ",
              tags$i("Conditional sensitivity"),
              " box, and ",
              "available as ",
              tags$code("Conditional_MDD_Log2FC"),
              " in the ",
              "downloaded results. It fixes the first two errors above but not ",
              "the third, which is why it is labelled a diagnostic."
            ),
            pq_eq(
              "\u03B1_eff   = FDR \u00D7 R / m        (R = rejections, m = tests; falls back to FDR/m when R = 0)",
              tags$br(),
              "df      = d_g + d\u2080",
              tags$br(),
              "ncp     = solve  P(|T(df, ncp)| > t_crit(\u03B1_eff, df)) = target power",
              tags$br(),
              "MDD_g   = ncp \u00D7 s_post,g \u00D7 stdev.unscaled_g"
            ),
            tags$p(
              "The non-centrality parameter is solved numerically from the ",
              "non-central t distribution rather than taken from a normal ",
              "approximation, because the moderated df can still be modest. ",
              "The ",
              tags$code("stdev.unscaled"),
              " factor is limma's own ",
              "design term, so unbalanced groups and non-pairwise contrasts are ",
              "handled without assuming a balanced two-group layout."
            ),
            pq_warn(
              tags$b("What this number is."),
              " \u201CGiven the variance I ",
              "estimated for this protein, and the threshold this run actually ",
              "operated at, an effect would have needed to be at least this ",
              "large.\u201D It is conditional on the observed data, so it describes ",
              "sensitivity ",
              tags$i("after the fact"),
              ". It is ",
              tags$b("not"),
              " used to call significance."
            ),

            tags$h3("Prospective MDD \u2014 the simulation"),
            tags$p(
              "This is the number to quote. It is produced by simulating data ",
              "that never existed, so nothing about it is conditioned on the ",
              "effects you happened to observe."
            ),
            tags$h4("The generative model"),
            tags$p(
              tags$code("eBayes()"),
              " does not merely shrink variances \u2014 it ",
              "fits an explicit hierarchical model in which each protein's true ",
              "variance is a draw from a scaled inverse chi-squared prior:"
            ),
            pq_eq(
              "\u03C3\u00B2_g  ~  d\u2080 \u00B7 s\u2080\u00B2 / \u03C7\u00B2(d\u2080)"
            ),
            tags$p(
              "Sampling from that prior is therefore not an approximation of ",
              "limma's assumptions; it ",
              tags$i("is"),
              " limma's assumptions. ",
              "This is what makes the simulation faithful rather than merely ",
              "plausible."
            ),
            tags$h4("The algorithm"),
            tags$ol(
              tags$li(
                "Draw a variance for every protein from the fitted prior. Under ",
                tags$code("Trend = TRUE"),
                " the prior is a per-protein function ",
                "of abundance, and it is kept as a vector, so the simulated data ",
                "reproduce your mean\u2013variance trend rather than a flat average ",
                "of it. If ",
                tags$code("d\u2080"),
                " is infinite the prior collapses ",
                "to a point mass at ",
                tags$code("s\u2080\u00B2"),
                ", which is the correct ",
                "limit \u2014 no arbitrary finite df is substituted."
              ),
              tags$li(
                "Build a fresh abundance matrix: observed mean abundance plus ",
                "Gaussian noise at those variances."
              ),
              tags$li(
                "Spike a known fold change into a randomly chosen ",
                tags$code("\u03C0\u2081"),
                " fraction of proteins, with random sign. ",
                "The offset is constructed as ",
                tags$code("\u03B2 = lfc \u00B7 w / w'w"),
                " for contrast vector ",
                tags$code("w"),
                ", which satisfies ",
                tags$code("w'\u03B2 = lfc"),
                " for any contrast, not just simple pairwise ones."
              ),
              tags$li(
                "Optionally knock out values and impute \u2014 see below."
              ),
              tags$li(
                "Run the real design and contrast matrix through ",
                tags$code("lmFit \u2192 contrasts.fit \u2192 eBayes \u2192 BH"),
                "."
              ),
              tags$li(
                "Power = the fraction of spiked proteins called at the FDR ",
                "cutoff. Proteins made unestimable by missingness stay in the ",
                "denominator: a protein that cannot be tested was not detected."
              ),
              tags$li(
                "Repeat, average, and interpolate the fold change at which the ",
                "curve crosses the target power."
              )
            ),
            pq_note(
              tags$b("Why this is prospective."),
              " Entirely new random matrices ",
              "are generated and explicit fold changes are spiked in ",
              tags$i("before"),
              " the FDR calculation. The measured quantity is ",
              "the pipeline's ability to recover known truth, which sidesteps ",
              "the data-conditional fallacy completely."
            ),
            pq_warn(
              tags$b("The honest caveat."),
              " The priors ",
              tags$code("s\u2080\u00B2"),
              " and ",
              tags$code("d\u2080"),
              " are still estimated from your data. ",
              "What the simulation removes is the observed-power fallacy, not ",
              "data-dependence as such. Describe the result as a design-level ",
              "sensitivity ",
              tags$i("given the estimated noise structure"),
              " \u2014 not as an independent, pre-registered power analysis."
            ),

            tags$h3("The missingness arm"),
            tags$p(
              "With ",
              tags$b("Include missingness + imputation"),
              " ticked, the ",
              "simulation stops assuming complete data. Dropout is measured from ",
              "your own matrix with a binomial GLM of missing-value counts on ",
              "mean log2 abundance:"
            ),
            pq_eq(
              "logit P(missing)  =  \u03B2\u2080 + \u03B2\u2081 \u00B7 mean_abundance"
            ),
            tags$p(
              "capturing the MNAR structure of proteomics, where low-abundance ",
              "proteins drop out preferentially. That model is replayed onto ",
              "every simulated matrix, the same imputation as the main pipeline ",
              "is applied, and the drop count per protein per group is capped so ",
              "the ",
              tags$code("min-valid"),
              " floor is respected."
            ),
            pq_warn(
              tags$b("Double-counting caveat."),
              " The variance prior was itself ",
              "estimated from already-imputed data. The comparison between the ",
              "missingness and complete-data arms therefore measures the ",
              tags$i("marginal"),
              " cost of a second imputation round, which ",
              "probably understates the total. Treat the gap as a lower bound."
            ),
            pq_note(
              "A useful side output is the realised false discovery proportion, ",
              "reported per contrast. Because the simulation knows ground truth, ",
              "it can check that BH is actually delivering the FDR you asked ",
              "for. A value near or below your cutoff means the machinery is ",
              "calibrated; a value meaningfully above it is evidence that a ",
              "pipeline step \u2014 usually imputation \u2014 is breaking FDR control."
            ),

            tags$h3("The replicate sweep"),
            tags$p(
              tags$b("Run replicate sweep"),
              " rebuilds a balanced design at each ",
              "replicate count you list and reruns the whole simulation, holding ",
              "the variance prior fixed. The result answers the design question ",
              "directly: how many replicates would I need to resolve a given ",
              "fold change?"
            ),
            tags$p(
              "Because only the crossing point is needed, the sweep brackets and ",
              "bisects rather than scanning a grid \u2014 roughly a threefold saving, ",
              "which matters because ",
              tags$code("lmFit"),
              " falls back to a ",
              "per-row fit whenever NAs are present."
            ),
            pq_warn(
              tags$b("Two limits to keep in mind."),
              " First, the prior is held ",
              "fixed across the sweep. That is defensible \u2014 it is a property of ",
              "the assay, not the sample size \u2014 but it was estimated at your ",
              "current residual df, so the small-n end is optimistic about how ",
              "well eBayes would shrink. Second, when the requested replicate ",
              "count equals the min-valid floor, no dropout can be simulated and ",
              "that row will show 0% missingness by construction."
            ),
            tags$p(
              "Expect the curve to follow roughly ",
              tags$code("1/\u221An"),
              ", degrading at the low end where the df penalty bites. Returns ",
              "typically flatten past about eight replicates."
            ),

            tags$h3("Choosing the parameters"),
            pq_table(
              c("Parameter", "Default", "How to choose"),
              list(
                c(
                  "Assumed proportion differential (&pi;<sub>1</sub>)",
                  "0.10",
                  "The fraction of proteins truly changing. Affects the BH threshold and therefore power. Use your observed hit rate as a guide; report the value, and check sensitivity at 0.01 and 0.30 if it is near a decision boundary."
                ),
                c(
                  "Target power",
                  "0.80",
                  "Convention. Raise to 0.90 for confirmatory work; the MDD will rise accordingly."
                ),
                c(
                  "Max log<sub>2</sub>FC tested",
                  "2",
                  "Upper end of the grid. Raise if the summary reports that the target power was never reached."
                ),
                c(
                  "log<sub>2</sub>FC grid step",
                  "0.2",
                  "Resolution of the curve. The crossing region is steep, so 0.1 gives a materially tighter estimate; 0.2 is adequate for a first look."
                ),
                c(
                  "Simulations per grid point",
                  "5",
                  "Monte Carlo replicates. 5 for exploration, 20 for a number you intend to publish. Cost is linear."
                ),
                c(
                  "Random seed",
                  "4817",
                  "Fixed for reproducibility. Change it to confirm your conclusion is not seed-specific."
                ),
                c(
                  "Include missingness + imputation",
                  "on",
                  "Leave on for realism. Turn off to isolate how much sensitivity the missing-data handling costs you."
                ),
                c(
                  "Replicates per group to sweep",
                  "3,4,5,6,8,10,12,15",
                  "Comma-separated, minimum 2. Include your current n so the plot has a reference point."
                )
              )
            ),

            tags$h3("Reading the outputs"),
            tags$h4("Power curve"),
            tags$p(
              "Power against true spiked fold change, with a shaded 95% Monte ",
              "Carlo band. A wide band means you need more simulations. The ",
              "dashed vertical line is the MDD."
            ),
            tags$h4("MDD summary table"),
            pq_table(
              c("Column", "Meaning"),
              list(
                c(
                  "<code>MDD (log2FC)</code>",
                  "Fold change, in log2 units, detectable at the target power"
                ),
                c(
                  "<code>MDD (fold change)</code>",
                  "The same number on a linear scale &mdash; the one to quote in a manuscript"
                ),
                c(
                  "<code>Mean realised FDP</code>",
                  "Ground-truth false discovery proportion; should sit near or below your FDR cutoff"
                ),
                c(
                  "<code>Note</code>",
                  "Flags a target power never reached, or already reached at the smallest fold change tested"
                )
              )
            ),
            tags$h4("Suggestion on how to report the result"),
            pq_note(
              tags$i(
                "\u201CAt 80% power and a 5% BH-FDR, and assuming 10% of proteins ",
                "are differentially abundant, this design resolves a 1.4-fold ",
                "change (0.49 log2). The estimate comes from parametric ",
                "simulation using the empirical Bayes variance prior fitted to ",
                "these data, including the observed intensity-dependent ",
                "missingness.\u201D"
              )
            )
          )
        ),

        # ── 5. Downstream ───────────────────────────────────────────────────
        tabPanel(
          "Enrichment & networks",
          div(
            class = "pq-doc",
            tags$h3("GO over-representation"),
            tags$p(
              "Runs ",
              tags$code("clusterProfiler::enrichGO"),
              " once per ",
              "selected contrast ",
              tags$i("and per direction"),
              ", so ",
              "Increased and Decreased proteins are tested separately rather ",
              "than pooled."
            ),
            pq_table(
              c("Control", "Default", "Guidance"),
              list(
                c(
                  "Organism (OrgDb)",
                  "<i>Homo sapiens</i>",
                  "Selects the Bioconductor annotation package, installed on first use. Also sets the STRING organism."
                ),
                c(
                  "Contrasts to compare",
                  "&mdash;",
                  "Fewer contrasts and a single ontology make the run considerably faster."
                ),
                c(
                  "GO Category",
                  "ALL",
                  "BP for mechanism, CC for localisation, MF for activity. ALL triples the number of tests."
                ),
                c(
                  "GO term FDR cutoff",
                  "0.2",
                  "Deliberately loose. GO terms are heavily correlated, so BH is conservative here; 0.2 is a common exploratory setting. Tighten to 0.05 for claims."
                ),
                c(
                  "Terms shown per contrast",
                  "5",
                  "Display only &mdash; does not affect the statistics or the download."
                )
              )
            ),
            pq_note(
              tags$b("The universe matters."),
              " The background is every protein ",
              "quantified in your experiment, not the whole genome.",
              "A protein that could never have been detected ",
              "cannot be enriched and it means results are not comparable ",
              "across experiments with different depth."
            ),
            tags$p(
              "Identifier type is detected automatically: if at least half the ",
              "IDs look like UniProt accessions the keytype is ",
              tags$code("UNIPROT"),
              ", otherwise ",
              tags$code("SYMBOL"),
              ". ",
              "A coverage diagnostic reports how many of your significant IDs ",
              "were actually annotated \u2014 check it, because low mapping rates ",
              "quietly bias enrichment."
            ),
            tags$h4("Plots"),
            tags$ul(
              tags$li(
                tags$b("Dotplot"),
                " \u2014 point size is GeneRatio (the fraction ",
                "of your list in the term), colour is adjusted p-value."
              ),
              tags$li(
                tags$b("Manhattan plot"),
                " \u2014 all terms laid out with ",
                tags$code("-log10(adj.P)"),
                " on the y-axis, useful for seeing ",
                "whether a handful of terms dominate."
              )
            ),

            tags$h3("STRING networks"),
            pq_table(
              c("Control", "Default", "Guidance"),
              list(
                c("Contrast", "&mdash;", "One contrast at a time."),
                c(
                  "Protein Set",
                  "combined",
                  "<code>combined</code> shows up- and down-regulated proteins together with directional colouring; the other options restrict to one direction."
                ),
                c(
                  "Min. combined score",
                  "400",
                  "STRING's confidence score, 0&ndash;1000. 150 = low, 400 = medium (default), 700 = high, 900 = highest. Raise it if the network is an unreadable hairball."
                )
              )
            ),
            pq_warn(
              "The STRING combined score aggregates evidence channels including ",
              "text mining, so a high score is not proof of physical ",
              "interaction. Well-studied proteins are systematically better ",
              "connected, which makes hub status partly a literature artefact."
            ),
            tags$p(
              "Proteins STRING could not map are listed in a table beneath the ",
              "network. A large unmapped fraction usually means the organism ",
              "selection is wrong or the identifiers are not resolving."
            ),

            tags$h3("Other views"),
            tags$ul(
              tags$li(
                tags$b("UpSet"),
                " \u2014 detection overlap between conditions, ",
                "based on presence/absence rather than on the model. This is ",
                "where on/off proteins excluded by the imputation-driven mask ",
                "become visible. The intersection table is downloadable."
              ),
              tags$li(
                tags$b("Z-Score Heatmap"),
                " \u2014 row-scaled abundances with ",
                "hierarchical clustering. ",
                tags$code("sig_reliable"),
                " restricts to significant proteins; ",
                tags$code("top_n"),
                " takes the most variable N. Note that clustering on ",
                "significant proteins will always look clean, because those ",
                "proteins were selected for separating the groups."
              ),
              tags$li(
                tags$b("Selected Proteins"),
                " \u2014 raw log2 abundance per ",
                "condition for proteins you name. Always worth checking a few ",
                "top hits here: a \u201Csignificant\u201D protein driven by one outlier ",
                "replicate is obvious in this view and invisible in a volcano ",
                "plot."
              )
            ),

            tags$h3("Downloads"),
            pq_table(
              c("Button", "Produces"),
              list(
                c(
                  "Download limma results",
                  "<code>limma_results_DATE.tsv</code> &mdash; full table including logFC, adj.P.Val, status, Sigma, alpha_eff, df_total and Conditional_MDD_Log2FC"
                ),
                c(
                  "Download ORA table",
                  "<code>GO_enrichment_enrichGO_DATE.tsv</code> with contrast and direction columns"
                ),
                c(
                  "Download ORA Plots",
                  "Zip of dotplot and Manhattan plot at 300 dpi"
                ),
                c("Download STRING network", "PNG at 200 dpi"),
                c("Download UpSet table", "Protein-to-intersection mapping"),
                c(
                  "Download all plots",
                  "Zip of every figure in the module at 300 dpi"
                )
              )
            )
          )
        ),

        # ── 6. Troubleshooting ──────────────────────────────────────────────
        tabPanel(
          "Troubleshooting",
          div(
            class = "pq-doc",
            tags$h3("Common problems"),
            pq_table(
              c("Symptom", "Likely cause", "Fix"),
              list(
                c(
                  "No significant proteins anywhere",
                  "Genuinely no effect, or variance inflated by imputation, or batch confounded with condition",
                  "Check the p-value histogram first. Flat means no signal. Right-skewed means inflated variance &mdash; raise the min-valid filter and re-run."
                ),
                c(
                  "Implausibly many significant proteins",
                  "MinProb imputation manufacturing differences, or a normalization failure",
                  "Compare against <code>ls</code> mode with no imputation; inspect the normalized boxplots."
                ),
                c(
                  "MDD simulation reports &ldquo;target power not reached&rdquo;",
                  "The grid stops below the true MDD",
                  "Raise <b>Max log<sub>2</sub>FC tested</b>."
                ),
                c(
                  "MDD is noisy between runs",
                  "Too few simulations per grid point",
                  "Raise <b>Simulations per grid point</b> to 20 and narrow the grid step."
                ),
                c(
                  "Replicate sweep is very slow",
                  "Missingness enabled forces <code>lmFit</code> onto a per-row path",
                  "Expect a few minutes. Reduce the number of replicate counts, or untick missingness for a quick look."
                ),
                c(
                  "Realised FDP well above the FDR cutoff",
                  "A pipeline step is breaking FDR control &mdash; usually imputation",
                  "Re-run the simulation with missingness off to confirm, then reconsider the imputation method."
                ),
                c(
                  "eBayes fails or the results table is short",
                  "Proteins with unestimable variance dropped by the trend filter",
                  "Raise the min-valid filter so every retained protein has enough observations."
                ),
                c(
                  "Most proteins unmapped in STRING or GO",
                  "Wrong organism, or identifiers not resolving",
                  "Check the organism selector and the coverage diagnostic; supply UniProt accessions if possible."
                ),
                c(
                  "PCA shows batches, not conditions, after correction",
                  "ComBat could not separate batch from condition",
                  "Check for confounding. If batch and condition coincide, correction is impossible &mdash; report it as a limitation."
                )
              )
            ),

            tags$h3("A defensible default configuration"),
            pq_table(
              c("Setting", "Value", "Why"),
              list(
                c(
                  "Normalization",
                  "<code>cyclicloess</code>",
                  "Handles intensity-dependent bias"
                ),
                c(
                  "Regression",
                  "<code>ls</code>",
                  "No imputation assumptions; limma handles NAs natively"
                ),
                c(
                  "eBayes Trend",
                  "<code>TRUE</code>",
                  "Proteomics almost always shows a mean&ndash;variance trend"
                ),
                c(
                  "Min. valid per group",
                  "<code>50</code>",
                  "Every protein observed in at least half the replicates of every group"
                ),
                c("Significance FDR", "<code>0.05</code>", "Convention"),
                c(
                  "Minimum |log<sub>2</sub>FC|",
                  "<code>0</code> or <code>0.58</code>",
                  "0 lets FDR decide; 0.58 adds a 1.5-fold floor if you want one"
                ),
                c(
                  "MDD &pi;<sub>1</sub>",
                  "match your observed hit rate",
                  "Affects the BH threshold and hence power"
                ),
                c(
                  "MDD simulations",
                  "<code>20</code>, step <code>0.1</code>",
                  "Tight enough to publish"
                )
              )
            ),

            tags$h3("What to report in a manuscript"),
            tags$ul(
              tags$li(
                "Normalization, imputation (or its absence), and batch correction, in pipeline order"
              ),
              tags$li("The design and the exact contrasts"),
              tags$li(
                "Both significance thresholds: the FDR level and any |log2FC| floor"
              ),
              tags$li("Whether the eBayes trend was used"),
              tags$li(
                HTML(
                  "The prospective MDD with its assumed &pi;<sub>1</sub>, target power and FDR"
                )
              ),
              tags$li(
                "That the MDD came from parametric simulation on the fitted empirical Bayes prior, and that the prior was estimated from the same data"
              ),
              tags$li("The min-valid filter and how many proteins survived it")
            ),

            tags$h3("References"),
            tags$ul(
              tags$li(
                "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. ",
                tags$i(
                  "limma powers differential expression analyses for RNA-sequencing and microarray studies."
                ),
                " Nucleic Acids Research. 2015;43(7):e47."
              ),
              tags$li(
                "Smyth GK. ",
                tags$i(
                  "Linear models and empirical Bayes methods for assessing differential expression in microarray experiments."
                ),
                " Statistical Applications in Genetics and Molecular Biology. 2004;3:Article 3."
              ),
              tags$li(
                "Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK. ",
                tags$i(
                  "Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression."
                ),
                " Annals of Applied Statistics. 2016;10(2):946-963."
              ),
              tags$li(
                "Benjamini Y, Hochberg Y. ",
                tags$i(
                  "Controlling the false discovery rate: a practical and powerful approach to multiple testing."
                ),
                " JRSS B. 1995;57(1):289-300."
              ),
              tags$li(
                "Johnson WE, Li C, Rabinovic A. ",
                tags$i(
                  "Adjusting batch effects in microarray expression data using empirical Bayes methods."
                ),
                " Biostatistics. 2007;8(1):118-127."
              ),
              tags$li(
                "Hoenig JM, Heisey DM. ",
                tags$i(
                  "The abuse of power: the pervasive fallacy of power calculations for data analysis."
                ),
                " The American Statistician. 2001;55(1):19-24."
              ),
              tags$li(
                "Lazar C, Gatto L, Ferro M, Bruley C, Burger T. ",
                tags$i(
                  "Accounting for the multiple natures of missing values in label-free quantitative proteomics data sets."
                ),
                " Journal of Proteome Research. 2016;15(4):1116-1125."
              ),
              tags$li(
                "Wu T, Hu E, Xu S, et al. ",
                tags$i(
                  "clusterProfiler 4.0: A universal enrichment tool for interpreting omics data."
                ),
                " The Innovation. 2021;2(3):100141."
              ),
              tags$li(
                "Szklarczyk D, Kirsch R, Koutrouli M, et al. ",
                tags$i("The STRING database in 2023."),
                " Nucleic Acids Research. 2023;51(D1):D638-D646."
              )
            )
          )
        )
      )
    ))
  )
}
