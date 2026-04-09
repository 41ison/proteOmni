## PSManalyst Shiny Module — wraps PSManalyst logic for proteOmni

# ── Helper functions ─────────────────────────────────────────────────────────

aa_freq <- function(x) table(x) / length(x) * 100

complete_and_reorder_amino_acids <- function(element) {
  aas <- c(
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
  for (a in aas) {
    if (!a %in% names(element)) element[[a]] <- 0
  }
  element[match(aas, names(element))]
}

#' Extract matrix of amino acids at P4-P4' positions from the fingerprint columns.
#'
#' @param data A tibble with columns \code{fingerprint_Nterm}, \code{fingerprint_Cterm}, and \code{sample_name}.
#' @return A named list of 20×8 numeric matrices of amino acid frequencies (%) per position,
#'   one matrix per sample.

extract_matrix_per_sample <- function(data) {
  filtered <- data %>%
    dplyr::filter(
      !is.na(fingerprint_Nterm) &
        !is.na(fingerprint_Cterm) &
        nchar(fingerprint_Nterm) == 8 &
        nchar(fingerprint_Cterm) == 8
    )

  sample_names <- dplyr::group_keys(dplyr::group_by(
    filtered,
    sample_name
  ))$sample_name

  mats <- dplyr::group_by(filtered, sample_name) %>%
    dplyr::group_map(
      ~ {
        fp <- strsplit(
          c(.x$fingerprint_Nterm, .x$fingerprint_Cterm),
          ""
        )

        if (length(fp) == 0) {
          return(NULL)
        }

        mat <- matrix(unlist(fp), ncol = 8, byrow = TRUE)
        mat <- mat[
          !apply(mat, 1, function(x) any(x %in% c("B", "X", "Z", "U"))),
          ,
          drop = FALSE
        ]

        if (nrow(mat) == 0) {
          return(NULL)
        }
        freq_list <- apply(mat, 2, aa_freq)
        new_list <- lapply(freq_list, complete_and_reorder_amino_acids)
        m <- matrix(unlist(new_list), ncol = 8, byrow = FALSE)
        m %>%
          magrittr::set_colnames(c(
            "P4",
            "P3",
            "P2",
            "P1",
            "P1'",
            "P2'",
            "P3'",
            "P4'"
          )) %>%
          magrittr::set_rownames(c(
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
          ))
      }
    )

  result <- stats::setNames(mats, sample_names)
  purrr::compact(result)
}

color_blue_seq <- c(
  "#d4e6f1",
  "#a9cce3",
  "#7fb3d5",
  "#5499c7",
  "#2980b9",
  "#1f618d",
  "#154360"
)

analyze_terminus_cooccurrence <- function(
  df,
  peptide_col,
  show_values = FALSE
) {
  aas <- c(
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

  # ── Per-sample computation ─────────────────────────────────────────────────
  samples <- sort(unique(df$sample_name))

  per_sample <- lapply(stats::setNames(samples, samples), function(s) {
    td <- df %>%
      dplyr::filter(sample_name == s) %>%
      dplyr::filter(
        !is.na(.data[[peptide_col]]) & nchar(.data[[peptide_col]]) > 0
      ) %>%
      dplyr::mutate(
        N_terminus = substr(.data[[peptide_col]], 1, 1),
        C_terminus = substr(
          .data[[peptide_col]],
          nchar(.data[[peptide_col]]),
          nchar(.data[[peptide_col]])
        )
      ) %>%
      dplyr::filter(N_terminus %in% aas & C_terminus %in% aas)

    ct <- table(td$N_terminus, td$C_terminus)
    total <- sum(ct)
    prob <- ct / total

    full <- matrix(
      0,
      nrow = length(aas),
      ncol = length(aas),
      dimnames = list(aas, aas)
    )
    for (i in rownames(ct)) {
      for (j in colnames(ct)) {
        full[i, j] <- prob[i, j]
      }
    }

    hd <- expand.grid(
      N_terminus = rownames(full),
      C_terminus = colnames(full)
    ) %>%
      dplyr::mutate(
        Probability = as.vector(full),
        N_terminus = factor(N_terminus, levels = rownames(full)),
        C_terminus = factor(C_terminus, levels = colnames(full))
      )

    p <- ggplot(hd, aes(x = C_terminus, y = N_terminus, fill = Probability)) +
      geom_tile(color = "white", linewidth = 0.1) +
      scale_fill_gradient2(
        low = "white",
        mid = "#B8E6D9",
        high = "#1BB99A",
        midpoint = max(hd$Probability) / 2,
        name = "Probability"
      ) +
      labs(
        title = s,
        x = "C-terminus AA",
        y = "N-terminus AA"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(
          hjust = 0.5,
          size = 8,
          face = "bold",
          color = "black"
        ),
        axis.text.y = element_text(size = 8, face = "bold", color = "black"),
        axis.title = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 9, face = "bold", hjust = 0.5),
        legend.title.position = "top",
        legend.text = element_text(size = 9, face = "bold"),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA)
      ) +
      coord_fixed()

    if (show_values) {
      p <- p +
        geom_text(
          aes(label = sprintf("%.3f", Probability)),
          size = 1.5,
          color = "black"
        )
    }

    list(
      matrix = full,
      plot = p,
      total_peptides = total,
      summary_stats = list(
        min_prob = min(full[full > 0]),
        max_prob = max(full),
        mean_prob = mean(full[full > 0]),
        n_observed_pairs = sum(full > 0)
      )
    )
  })

  # ── Combined patchwork plot ────────────────────────────────────────────────
  plots <- lapply(per_sample, `[[`, "plot")
  n_samples <- length(plots)
  ncols <- min(3L, n_samples)

  combined_plot <- patchwork::wrap_plots(plots, ncol = ncols) +
    patchwork::plot_annotation(
      title = "N:C-terminus co-occurrence",
      theme = theme(
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
      )
    )

  list(
    per_sample = per_sample,
    plot = combined_plot,
    sample_names = samples
  )
}


# ── Multi-file PSM ingestion helpers ────────────────────────────────────

#' Recursively find all psm.tsv files inside a directory tree.
#'
#' @param root_path  Character. Top-level folder chosen by the user.
#' @return Character vector of full file paths.
find_psm_files <- function(root_path) {
  list.files(
    path = root_path,
    pattern = "^psm\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
}

#' Read one psm.tsv and add a sample_name column derived from the Spectrum column.
#'
#' @param path  Character. Full path to a single psm.tsv file.
#' @return A tibble with an extra \code{sample_name} column.
read_single_psm <- function(path) {
  suppressMessages(readr::read_tsv(path)) %>%
    janitor::clean_names() %>%
    dplyr::mutate(
      sample_name = stringr::str_extract(spectrum, "^[^.]+")
    )
}

#' Read and row-bind all psm.tsv files found under \code{root_path}.
#'
#' @param root_path  Character. Top-level folder.
#' @return A combined tibble with a \code{sample_name} column.
read_all_psm_files <- function(root_path) {
  paths <- find_psm_files(root_path)
  if (length(paths) == 0L) {
    stop("No psm.tsv files found under: ", root_path)
  }
  purrr::map_dfr(paths, read_single_psm)
}


# ── MS/MS spectrum builder (FragPipe PSM columns) ──────────────────────

#' Tidy the combined PSM data for a single peptide sequence.
#'
#' @param data        Combined PSM tibble (output of \code{read_all_psm_files()}).
#' @param peptide_seq Character. Peptide sequence to subset.
#' @return A tidy tibble ready for \code{build_psm_spectrum()}.
tidy_psm_spectrum <- function(data, peptide_seq) {
  data %>%
    dplyr::filter(peptide == peptide_seq) %>%
    dplyr::select(
      sample_name,
      peptide,
      gene,
      charge,
      ions,
      ion_mz,
      ion_int,
      retention
    ) %>%
    dplyr::mutate(
      retention = round(retention / 60, 3),
      ion_list = stringr::str_remove(ions, "\\[\\]"),
      mz_list = stringr::str_remove(ion_mz, "\\[\\]"),
      int_list = stringr::str_remove(ion_int, "\\[\\]"),
      ion_list = stringr::str_replace_all(ion_list, "\\s+", ""),
      mz_list = stringr::str_replace_all(mz_list, "\\s+", ""),
      int_list = stringr::str_replace_all(int_list, "\\s+", ""),
      ion_list = stringr::str_split(ions, ","),
      mz_list = stringr::str_split(ion_mz, ","),
      int_list = stringr::str_split(ion_int, ",")
    ) %>%
    tidyr::unnest(cols = c(ion_list, mz_list, int_list)) %>%
    dplyr::rename(ion = ion_list, mz = mz_list, intensity = int_list) %>%
    dplyr::mutate(
      mz = as.numeric(mz),
      intensity = as.numeric(intensity)
    ) %>%
    dplyr::filter(!is.na(mz), intensity > 0)
}


#' Build an annotated MS/MS fragmentation spectrum ggplot from FragPipe PSM data.
#'
#' @param tidy_data  Tibble from \code{tidy_psm_spectrum()}.
#' @param label_size Numeric. Font size for ion annotations.
#' @return A ggplot object.
build_psm_spectrum <- function(tidy_data, label_size = 3) {
  if (nrow(tidy_data) == 0L) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "Peptide not found in the uploaded data.",
          size = 6,
          colour = "grey50"
        ) +
        theme_void()
    )
  }

  gene_name <- dplyr::first(stats::na.omit(tidy_data$gene))
  peptide <- dplyr::first(tidy_data$peptide)

  ggplot(tidy_data, aes(x = mz, y = intensity)) +
    geom_segment(aes(xend = mz, yend = 0), colour = "grey25", linewidth = 0.4) +
    geom_text(
      aes(
        label = ion,
        color = dplyr::case_when(
          stringr::str_detect(ion, "^y") ~ "y",
          stringr::str_detect(ion, "^b") ~ "b",
          TRUE ~ "other"
        )
      ),
      vjust = -0.8,
      size = label_size,
      fontface = "bold",
      check_overlap = TRUE
    ) +
    scale_color_manual(
      values = c("y" = "#d95f02", "b" = "#1b9e77", "other" = "#7570b3"),
      guide = "none"
    ) +
    facet_wrap(
      ~ sample_name + charge + retention,
      ncol = 2,
      scales = "free",
      labeller = label_both
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15)),
      labels = scales::label_scientific()
    ) +
    scale_x_continuous(
      breaks = tidy_data$mz,
      labels = scales::label_number(accuracy = 0.1)
    ) +
    labs(
      title = paste0(
        "MS/MS Fragmentation: ",
        peptide,
        if (!is.na(gene_name)) paste0(" (", gene_name, ")") else ""
      ),
      x = "*m/z*",
      y = "Intensity"
    ) +
    theme_bw() +
    theme(
      plot.title = ggtext::element_markdown(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold", size = 8),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 8,
        color = "black"
      ),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title = ggtext::element_markdown(size = 12, face = "bold"),
      axis.ticks = element_line(color = "black", linewidth = 0.25)
    )
}


# ── UI ────────────────────────────────────────────────────────────────────────
PSManalyst_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # ── Sidebar content ───────────────────────────────────────────────────────
    tags$div(
      id = ns("sidebar_content"),

      # ── PSM folder upload ─────────────────────────────────────────────────
      tags$div(
        style = "padding:12px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "PSM Viewer"
      ),
      textInput(
        ns("psm_folder"),
        "Path to folder containing psm.tsv files",
        placeholder = "/path/to/fragpipe/output"
      ),
      actionButton(
        ns("load_psm_folder"),
        "Load PSM Files",
        class = "btn-primary btn-sm",
        style = "width:90%;margin:4px auto 8px;display:block;"
      ),
      uiOutput(ns("psm_folder_status")), # feedback: n files found
      sliderInput(
        ns("hyperscore"),
        "PSM hyperscore filter",
        min = 0,
        max = 1000,
        value = 0,
        step = 5
      ),
      sliderInput(
        ns("probability"),
        "PeptideProphet Probability",
        min = 0,
        max = 1,
        value = 0.95,
        step = 0.01
      ),
      selectInput(
        ns("specificity_filter"),
        "Proteolysis fingerprinting specificity",
        choices = c(
          "All" = "all",
          "Fully specific" = "fully_specific",
          "Semi-specific at N-termini" = "semi_n_termini",
          "Semi-specific at C-termini" = "semi_c_termini",
          "Fully semi-specific" = "fully_semi_specific"
        ),
        selected = "all"
      ),
      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── MS/MS Spectrum viewer controls ────────────────────────────────────
      tags$div(
        style = "padding:12px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "MS/MS Spectrum Viewer"
      ),
      selectizeInput(
        ns("spectrum_peptide"),
        "Select Peptide Sequence",
        choices = NULL,
        options = list(
          placeholder = "Load PSM files first…",
          maxOptions = 5000,
          searchField = "value"
        )
      ),
      numericInput(
        ns("spectrum_label_size"),
        "Ion Label Size",
        value = 3,
        min = 1,
        max = 8,
        step = 0.5
      ),
      actionButton(
        ns("run_spectrum"),
        "Plot Spectrum",
        class = "btn-primary btn-sm",
        style = "width:90%;margin:4px auto 8px;display:block;"
      ),
      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── Protein Viewer controls ────────────────────────────────────────────
      tags$div(
        style = "padding:8px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "Protein Viewer"
      ),
      fileInput(
        ns("combined_protein"),
        "Choose the combined_protein.tsv file",
        accept = ".tsv"
      ),
      fileInput(
        ns("fasta_file"),
        "Choose the FASTA file",
        accept = c(".fasta", ".fa", ".fas")
      ),
      selectInput(
        ns("protein_sample_select"),
        "Sample (Protein View)",
        choices = c("All samples" = "all")
      ),
      selectInput(ns("xcol"), "X Sample", choices = NULL),
      selectInput(ns("ycol"), "Y Sample", choices = NULL),
      colourpicker::colourInput(
        ns("plot_color"),
        "Select plot color",
        value = "#5499c7"
      ),
      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── Modification Diagnostic controls ─────────────────────────────────
      tags$div(
        style = "padding:8px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "Mod. Diagnostic"
      ),
      sliderInput(
        ns("rt_tolerance"),
        "\u0394RT Artifact Threshold (min)",
        min = 0.1,
        max = 5,
        value = 0.5,
        step = 0.1
      ),
      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── Downloads ──────────────────────────────────────────────────────────
      tags$div(
        style = "padding:0 8px;",
        div(
          style = "margin-bottom:8px;",
          downloadButton(
            ns("download_all_plots"),
            "⬇ Download all plots",
            class = "dl-btn",
            style = "width:100%;text-align:left;"
          )
        ),
        div(
          downloadButton(
            ns("download_spectrum_plot"),
            "⬇ MS/MS Spectrum (.pdf)",
            class = "dl-btn",
            style = "width:100%;text-align:left;"
          )
        )
      )
    ),

    # ── Body ──────────────────────────────────────────────────────────────────
    tabsetPanel(
      id = ns("tabs"),
      type = "tabs",

      # ── Tab 1: PSM Viewer ─────────────────────────────────────────────────
      tabPanel(
        "PSM Viewer",
        fluidRow(infoBoxOutput(ns("info_box1"), width = 12)),
        fluidRow(
          box(
            title = "Protease fingerprint",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp01"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot01_ui"))
            )
          ),
          box(
            title = "N-termini SeqLogo",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp03"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot03_ui"))
            )
          ),
          box(
            title = "C-termini SeqLogo",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp04"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot04_ui"))
            )
          ),
          box(
            title = "m/z over retention time",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp05"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot05_ui"))
            )
          ),
          box(
            title = "Mass error (ppm)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp06"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot06_ui"))
            )
          ),
          box(
            title = "Peptide length",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp07"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot07_ui"))
            )
          ),
          box(
            title = "GRAVY (Grand Average of Hydropathy)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp08"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot08_ui"))
            )
          ),
          box(
            title = "Isoelectric Point (pI)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp09"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot09_ui"))
            )
          ),
          box(
            title = "Charge state distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp10"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot10_ui"))
            )
          ),
          box(
            title = "Number of missed cleavages",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp11"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot11_ui"))
            )
          ),
          box(
            title = "Uniqueness",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp12"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot12_ui"))
            )
          ),
          box(
            title = "Hyperscore distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp13"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot13_ui"))
            )
          ),
          box(
            title = "Next Score distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp14"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot14_ui"))
            )
          ),
          box(
            title = "PeptideProphet probability",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp15"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot15_ui"))
            )
          ),
          box(
            title = "Expectation (PeptideProphet)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp16"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot16_ui"))
            )
          ),
          box(
            title = "Assigned modifications",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp17"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot17_ui"))
            )
          ),
          box(
            title = "Top 20 proteins with more PSMs",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp18"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot18_ui"))
            )
          ),
          box(
            title = "Co-occurrence probability matrix of N- and C-terminus AA",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp19"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot19_ui"))
            )
          ),
          box(
            title = "Cysteine counts in peptides",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp20"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot20_ui"))
            )
          )
        ),
        box(
          title = "Peptide Yield vs. FDR",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_fdr"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("plot_fdr_curve_ui"))
          )
        )
      ),

      # ── Tab 2: MS/MS Spectrum Viewer ──────────────────────────────────────
      tabPanel(
        title = tagList(icon("chart-bar"), "MS/MS Spectrum Viewer"),
        fluidRow(
          box(
            title = "Annotated MS/MS Fragmentation Spectrum (FragPipe PSMs)",
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
              uiOutput(ns("spectrum_plot_ui"))
            )
          )
        ),
        fluidRow(
          box(
            title = "Tidy fragment ion data for selected peptide",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput(ns("spectrum_tidy_table"))
          )
        )
      ),

      # ── Tab 3: Protein Viewer ─────────────────────────────────────────────
      tabPanel(
        "Protein Viewer",
        fluidRow(infoBoxOutput(ns("info_box2"), width = 12)),
        fluidRow(
          box(
            title = "Protein Sequence View",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            fluidRow(
              column(
                4,
                selectizeInput(
                  ns("selected_protein"),
                  "Select Protein:",
                  choices = NULL
                )
              ),
              column(
                4,
                numericInput(
                  ns("aa_per_line"),
                  "Amino Acids per Line:",
                  value = 50,
                  min = 10,
                  max = 200
                )
              ),
              column(
                4,
                checkboxInput(
                  ns("show_sequence"),
                  "Show Sequence Letters",
                  value = TRUE
                )
              )
            ),
            uiOutput(ns("protein_coverage_plot_ui"))
          ),
          box(
            title = "Amino Acid Frequencies (Expected vs Identified)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            uiOutput(ns("aa_freq_banner")),
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp_aa_psm"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot_aa_freq_ui"))
            )
          ),
          box(
            title = "Protein coverage summary",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp01p"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot01p_ui"))
            )
          ),
          box(
            title = "Number of proteins by organism",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp02p"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot02p_ui"))
            )
          ),
          box(
            title = "Protein existence evidence",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp03p"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot03p_ui"))
            )
          ),
          box(
            title = "Protein probability (ProteinProphet)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp04p"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot04p_ui"))
            )
          ),
          box(
            title = "Top Peptide Probability",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp05p"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot05p_ui"))
            )
          ),
          box(
            title = "Total peptides mapped to the proteins",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp06p"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot06p_ui"))
            )
          ),
          box(
            title = "Razor spectral count",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp07p"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot07p_ui"))
            )
          ),
          box(
            title = "Razor intensity",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp08p"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot08p_ui"))
            )
          ),
          box(
            title = "MaxLFQ intensity distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp10p"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot10p_ui"))
            )
          ),
          box(
            title = "Abundance correlation (MaxLFQ)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("spgg"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("plot_ggpairs_ui"))
            )
          ),
          box(
            title = "Sample correlation — Non-normalized log2(Intensity)",
            status = "primary",
            solidHeader = TRUE,
            height = 600,
            collapsible = FALSE,
            plotlyOutput(ns("plot11p"))
          ),
          tabBox(
            title = "Similarity metrics",
            side = "right",
            height = 600,
            tabPanel("Cosine similarity", uiOutput(ns("cosine_similarity_ui"))),
            tabPanel(
              "Euclidean distance",
              uiOutput(ns("euclidean_distance_ui"))
            ),
            tabPanel("Jaccard similarity", uiOutput(ns("jaccard_similarity_ui")))
          ),
          box(
            title = "Top 20 proteins — MaxLFQ intensity rank (log2)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp_rank_lfq"),
                icon("spinner", class = "fa-spin")
              ),
              plotOutput(ns("plot_rank_lfq"), height = "600px")
            )
          ),
          box(
            title = "Top 20 proteins — Unique spectral count rank",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp_rank_usc"),
                icon("spinner", class = "fa-spin")
              ),
              plotOutput(ns("plot_rank_usc"), height = "600px")
            )
          )
        )
      ),

      # ── Tab 4: Modification Diagnostic ─────────────────────────────────
      tabPanel(
        "Modification Diagnostic",
        fluidRow(
          box(
            title = "RT Shift Profile — Modified vs. Unmodified Peptides",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp_moddiag"),
                icon("spinner", class = "fa-spin")
              ),
              uiOutput(ns("mod_diag_plot_ui"))
            )
          )
        ),
        fluidRow(
          box(
            title = "RT Shift Table",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            DT::dataTableOutput(ns("mod_diag_table"))
          )
        )
      )
    )
  )
}


# ── Server ────────────────────────────────────────────────────────────────────
PSManalyst_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Dynamic height logic
    psm_plot_h <- reactive({
      d <- data()
      req(d, "sample_name" %in% names(d))
      n <- length(unique(d$sample_name))
      max(400L, ceiling(n / 3) * 350L)
    })
    
    make_psm_plot_ui <- function(plot_id) {
      renderUI({ req(data()); plotOutput(ns(plot_id), height = paste0(psm_plot_h(), "px")) })
    }
    
    output$plot01_ui <- make_psm_plot_ui("plot01")
    output$plot03_ui <- make_psm_plot_ui("plot03")
    output$plot04_ui <- make_psm_plot_ui("plot04")
    output$plot05_ui <- make_psm_plot_ui("plot05")
    output$plot06_ui <- make_psm_plot_ui("plot06")
    output$plot07_ui <- make_psm_plot_ui("plot07")
    output$plot08_ui <- make_psm_plot_ui("plot08")
    output$plot09_ui <- make_psm_plot_ui("plot09")
    output$plot10_ui <- make_psm_plot_ui("plot10")
    output$plot11_ui <- make_psm_plot_ui("plot11")
    output$plot12_ui <- make_psm_plot_ui("plot12")
    output$plot13_ui <- make_psm_plot_ui("plot13")
    output$plot14_ui <- make_psm_plot_ui("plot14")
    output$plot15_ui <- make_psm_plot_ui("plot15")
    output$plot16_ui <- make_psm_plot_ui("plot16")
    output$plot17_ui <- make_psm_plot_ui("plot17")
    output$plot18_ui <- make_psm_plot_ui("plot18")
    output$plot19_ui <- make_psm_plot_ui("plot19")
    output$plot20_ui <- make_psm_plot_ui("plot20")
    
    output$plot01p_ui <- make_psm_plot_ui("plot01p")
    output$plot02p_ui <- make_psm_plot_ui("plot02p")
    output$plot03p_ui <- make_psm_plot_ui("plot03p")
    output$plot04p_ui <- make_psm_plot_ui("plot04p")
    output$plot05p_ui <- make_psm_plot_ui("plot05p")
    output$plot06p_ui <- make_psm_plot_ui("plot06p")
    output$plot07p_ui <- make_psm_plot_ui("plot07p")
    output$plot08p_ui <- make_psm_plot_ui("plot08p")
    output$plot10p_ui <- make_psm_plot_ui("plot10p")
    
    output$plot_fdr_curve_ui <- make_psm_plot_ui("plot_fdr_curve")
    output$plot_aa_freq_ui <- make_psm_plot_ui("plot_aa_freq")
    output$plot_ggpairs_ui <- make_psm_plot_ui("plot_ggpairs")
    output$cosine_similarity_ui <- make_psm_plot_ui("cosine_similarity")
    output$euclidean_distance_ui <- make_psm_plot_ui("euclidean_distance")
    output$jaccard_similarity_ui <- make_psm_plot_ui("jaccard_similarity")


    # ── Spinner helpers ───────────────────────────────────────────────────────
    spin_ids_psm <- paste0(
      "sp",
      c(
        "01",
        "03",
        "04",
        "05",
        "06",
        "07",
        "08",
        "09",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "_fdr"
      )
    )
    spin_ids_prot <- c(
      paste0(
        "sp",
        c(
          "01p",
          "02p",
          "03p",
          "04p",
          "05p",
          "06p",
          "07p",
          "08p",
          "09p",
          "10p",
          "gg"
        )
      ),
      "sp_aa_psm"
    )

    show_psm_spinners <- function() {
      lapply(spin_ids_psm, function(s) shinyjs::show(id = s))
    }
    show_prot_spinners <- function() {
      lapply(spin_ids_prot, function(s) shinyjs::show(id = s))
    }
    hide_spinner <- function(sid) shinyjs::hide(id = sid)
    rh <- function(expr_fn, sid) {
      on.exit(hide_spinner(sid), add = TRUE)
      expr_fn()
    }

    # ════════════════════════════════════════════════════════════════════════
    # A. MULTI-FILE PSM INGESTION
    # ════════════════════════════════════════════════════════════════════════

    # Holds the raw combined data (before QC filters)
    raw_psm_data <- reactiveVal(NULL)
    loaded_psm_folder <- reactiveVal(NULL)

    # Load files when the button is clicked
    observeEvent(input$load_psm_folder, {
      folder <- trimws(input$psm_folder)
      req(nchar(folder) > 0)

      if (!dir.exists(folder)) {
        showNotification(
          paste("Folder not found:", folder),
          type = "error",
          duration = 6
        )
        return()
      }

      withProgress(message = "Scanning for psm.tsv files…", value = 0.2, {
        result <- tryCatch(
          {
            df <- read_all_psm_files(folder)
            incProgress(
              0.8,
              detail = paste(
                "Loaded",
                dplyr::n_distinct(df$sample_name),
                "sample(s),",
                nrow(df),
                "PSMs"
              )
            )
            df
          },
          error = function(e) {
            showNotification(conditionMessage(e), type = "error", duration = 8)
            NULL
          }
        )
        raw_psm_data(result)
        loaded_psm_folder(folder)
      })
    })

    # Feedback label below the folder input
    output$psm_folder_status <- renderUI({
      df <- raw_psm_data()
      if (is.null(df)) {
        return(NULL)
      }
      tags$div(
        style = "padding:4px 16px 8px;font-size:12px;color:#1BB99A;font-weight:600;",
        icon("circle-check", style = "margin-right:4px;"),
        sprintf(
          "%d sample(s) loaded · %s PSMs",
          dplyr::n_distinct(df$sample_name),
          format(nrow(df), big.mark = ",")
        )
      )
    })

    # ════════════════════════════════════════════════════════════════════════
    # B. FILTERED PSM DATA
    # ════════════════════════════════════════════════════════════════════════

    data <- reactive({
      req(raw_psm_data())
      show_psm_spinners()

      psm_file <- raw_psm_data() %>%
        dplyr::filter(
          hyperscore >= input$hyperscore,
          probability >= input$probability
        ) %>%
        dplyr::mutate(
          fingerprint_Nterm = dplyr::case_when(
            stringr::str_detect(extended_peptide, "^\\.") ~ NA_character_,
            TRUE ~ substr(extended_peptide, 2, 16)
          ),
          fingerprint_Cterm = dplyr::case_when(
            stringr::str_detect(extended_peptide, "\\.$") ~ NA_character_,
            TRUE ~ stringr::str_sub(extended_peptide, -16, -2)
          ),
          fingerprint_Nterm = stringr::str_extract(
            fingerprint_Nterm,
            ".{4}\\..{4}"
          ),
          fingerprint_Nterm = stringr::str_remove_all(fingerprint_Nterm, "\\."),
          fingerprint_Cterm = stringr::str_extract(
            fingerprint_Cterm,
            ".{4}\\..{4}"
          ),
          fingerprint_Cterm = stringr::str_remove_all(fingerprint_Cterm, "\\."),
          delta_mass_ppm = (observed_m_z - calculated_m_z) /
            calculated_m_z *
            1e6,
          gravy = sapply(peptide, GRAVY),
          isoelectric_point = sapply(peptide, calculate_pI),
          specificity = dplyr::case_when(
            prev_aa %in%
              c("K", "R") &
              stringr::str_sub(peptide, -1) %in% c("K", "R") ~ "Fully specific",
            prev_aa %in%
              c("K", "R") &
              !stringr::str_sub(peptide, -1) %in%
                c("K", "R") ~ "Semi-specific at C-termini",
            !prev_aa %in% c("K", "R") &
              stringr::str_sub(peptide, -1) %in%
                c("K", "R") ~ "Semi-specific at N-termini",
            TRUE ~ "Fully semi-specific"
          )
        ) %>%
        dplyr::relocate(extended_peptide, .before = fingerprint_Nterm) %>%
        dplyr::relocate(specificity, .after = peptide)

      if (input$specificity_filter != "all") {
        spec_map <- c(
          fully_specific = "Fully specific",
          semi_n_termini = "Semi-specific at N-termini",
          semi_c_termini = "Semi-specific at C-termini",
          fully_semi_specific = "Fully semi-specific"
        )
        psm_file <- psm_file %>%
          dplyr::filter(specificity == spec_map[input$specificity_filter])
      }
      psm_file
    })

    frequency_matrix_of_aa <- reactive({
      req(data())
      extract_matrix_per_sample(data())
    })
    cooccurrence_data <- reactive({
      req(data())
      analyze_terminus_cooccurrence(data(), "peptide", show_values = TRUE)
    })

    # ════════════════════════════════════════════════════════════════════════
    # C. MS/MS SPECTRUM VIEWER
    # ════════════════════════════════════════════════════════════════════════

    # Populate peptide selector once PSM data is loaded
    observeEvent(raw_psm_data(), {
      req(raw_psm_data())
      peptides <- sort(unique(raw_psm_data()$peptide))
      updateSelectizeInput(
        session,
        "spectrum_peptide",
        choices = peptides,
        selected = peptides[1],
        server = TRUE
      )
    })

    # Tidy fragment data — recomputed only on button click
    tidy_spectrum_data <- eventReactive(input$run_spectrum, {
      req(raw_psm_data(), nchar(input$spectrum_peptide) > 0)
      shinyjs::show(id = "sp_spectrum")

      withProgress(message = "Building spectrum…", value = 0.5, {
        result <- tidy_psm_spectrum(raw_psm_data(), input$spectrum_peptide)
        incProgress(0.5, detail = "Done.")
        result
      })
    })

    # Number of facets → dynamic height
    n_spectrum_facets <- reactive({
      req(tidy_spectrum_data())
      tidy_spectrum_data() %>%
        dplyr::distinct(sample_name, charge, retention) %>%
        nrow()
    })

    spectrum_height_px <- reactive({
      n_rows <- ceiling(n_spectrum_facets() / 3)
      max(400L, n_rows * 350L)
    })

    # Dynamic-height plot container
    output$spectrum_plot_ui <- renderUI({
      plotOutput(
        ns("spectrum_plot"),
        height = paste0(spectrum_height_px(), "px")
      )
    })

    output$spectrum_plot <- renderPlot({
      on.exit(hide_spinner("sp_spectrum"), add = TRUE)
      req(tidy_spectrum_data())
      build_psm_spectrum(
        tidy_spectrum_data(),
        label_size = input$spectrum_label_size
      )
    })

    # Tidy data table beneath the spectrum
    output$spectrum_tidy_table <- DT::renderDataTable({
      req(tidy_spectrum_data())
      DT::datatable(
        tidy_spectrum_data(),
        rownames = FALSE,
        options = list(dom = "frtip", pageLength = 20, scrollX = TRUE),
        class = "display compact"
      )
    })

    output$download_spectrum_plot <- downloadHandler(
      filename = function() {
        paste0("psm_spectrum_", input$spectrum_peptide, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        p <- build_psm_spectrum(
          tidy_spectrum_data(),
          label_size = input$spectrum_label_size
        )
        n_rows <- ceiling(n_spectrum_facets() / 3)
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

    # ════════════════════════════════════════════════════════════════════════
    # D. INFO BOXES
    # ════════════════════════════════════════════════════════════════════════

    output$info_box1 <- renderInfoBox({
      spec_labels <- c(
        fully_specific = "Fully specific",
        semi_n_termini = "Semi-specific at N-termini",
        semi_c_termini = "Semi-specific at C-termini",
        fully_semi_specific = "Fully semi-specific"
      )
      ft <- paste(
        "Showing PSMs with Hyperscore >=",
        input$hyperscore,
        "and PeptideProphet probability >=",
        input$probability
      )
      if (input$specificity_filter != "all") {
        ft <- paste(
          ft,
          "and",
          spec_labels[input$specificity_filter],
          "peptides"
        )
      }

      # Add sample count when data is available
      df <- raw_psm_data()
      if (!is.null(df)) {
        ft <- paste0(
          ft,
          " | ",
          dplyr::n_distinct(df$sample_name),
          " sample(s), ",
          format(nrow(data()), big.mark = ","),
          " PSMs after filtering"
        )
      }

      infoBox("Filter settings", ft, icon = icon("info"), color = "black")
    })

    output$info_box2 <- renderInfoBox({
      infoBox(
        "protein.tsv files contain FDR-filtered protein results, where each row is an identified protein group",
        icon = icon("info"),
        color = "black"
      )
    })

    # ════════════════════════════════════════════════════════════════════════
    # E. PSM PLOTS
    # ════════════════════════════════════════════════════════════════════════

    seqlogo_scale <- scale_x_continuous(
      breaks = 1:8,
      labels = c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
    )

    output$plot01 <- renderPlot({
      rh(
        function() {
          mat_list <- frequency_matrix_of_aa()
          purrr::imap_dfr(
            mat_list,
            ~ as.data.frame(.x) %>%
              rownames_to_column(var = "residue") %>%
              pivot_longer(
                cols = -residue,
                names_to = "position",
                values_to = "frequency"
              ) %>%
              dplyr::mutate(sample_name = .y)
          ) %>%
            dplyr::mutate(
              position = factor(
                position,
                c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
              ),
              residue = factor(
                residue,
                c(
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
              )
            ) %>%
            ggplot(aes(x = position, y = residue, fill = frequency)) +
            geom_tile(color = "black") +
            scale_fill_gradient(low = "#d4e6f1", high = "#154360") +
            geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
            theme_minimal() +
            facet_wrap(~sample_name, ncol = 3) +
            labs(
              title = "Cleavage Site Specificity",
              x = "Position",
              y = "Amino acid residue",
              fill = "Frequency (%)"
            ) +
            theme(
              plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
              axis.text.x = element_text(
                size = 12,
                face = "bold",
                color = "black"
              ),
              axis.text.y = element_text(
                size = 11,
                face = "bold",
                color = "black"
              ),
              axis.title = element_text(size = 13, face = "bold"),
              strip.background = element_blank(),
              strip.text = element_text(
                size = 14,
                face = "bold",
                color = "black"
              ),
              legend.position = "bottom",
              legend.key.width = unit(2, "cm"),
              legend.key.height = unit(0.25, "cm"),
              legend.title.position = "top",
              panel.grid = element_blank()
            )
        },
        "sp01"
      )
    })

    output$plot03 <- renderPlot({
      rh(
        function() {
          seq_lst <- data() %>%
            as.data.frame() %>%
            dplyr::filter(!is.na(fingerprint_Nterm)) %>%
            dplyr::group_by(sample_name) %>%
            dplyr::summarise(
              seqs = list(fingerprint_Nterm),
              .groups = "drop"
            ) %>%
            tibble::deframe()
          ggseqlogo::ggseqlogo(seq_lst, method = "bits", seq_type = "AA") +
            geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
            geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
            seqlogo_scale +
            labs(
              title = "SeqLogo of the N-termini fingerprint",
              x = "Amino acid position",
              y = "Bits"
            ) +
            theme_bw() +
            theme(
              plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
              axis.text.x = element_text(
                size = 12,
                face = "bold",
                color = "black"
              ),
              axis.text.y = element_text(
                size = 11,
                face = "bold",
                color = "black"
              ),
              axis.title = element_text(size = 13, face = "bold"),
              strip.background = element_blank(),
              strip.text = element_text(
                size = 14,
                face = "bold",
                color = "black"
              )
            )
        },
        "sp03"
      )
    })

    output$plot04 <- renderPlot({
      rh(
        function() {
          seq_lst <- data() %>%
            as.data.frame() %>%
            dplyr::filter(!is.na(fingerprint_Cterm)) %>%
            dplyr::group_by(sample_name) %>%
            dplyr::summarise(
              seqs = list(fingerprint_Cterm),
              .groups = "drop"
            ) %>%
            tibble::deframe()
          ggseqlogo::ggseqlogo(seq_lst, method = "bits", seq_type = "AA") +
            geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
            geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
            seqlogo_scale +
            labs(
              title = "SeqLogo of the C-termini fingerprint",
              x = "Amino acid position",
              y = "Bits"
            ) +
            theme_bw() +
            theme(
              plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
              axis.text.x = element_text(
                size = 12,
                face = "bold",
                color = "black"
              ),
              axis.text.y = element_text(
                size = 11,
                face = "bold",
                color = "black"
              ),
              axis.title = element_text(size = 13, face = "bold"),
              strip.background = element_blank(),
              strip.text = element_text(
                size = 14,
                face = "bold",
                color = "black"
              )
            )
        },
        "sp04"
      )
    })

    output$plot05 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            ggplot(aes(x = retention / 60, y = observed_m_z)) +
            ggpointdensity::geom_pointdensity(size = 0.25) +
            viridis::scale_color_viridis(option = "plasma") +
            labs(
              x = "Retention time (min)",
              y = "Scan range (m/z)",
              color = "Number of Neighborhoods"
            ) +
            facet_wrap(~sample_name, ncol = 3) +
            theme(
              legend.position = "bottom",
              legend.key.width = unit(2, "cm"),
              legend.key.height = unit(0.25, "cm")
            )
        },
        "sp05"
      )
    })

    output$plot06 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            dplyr::filter(abs(delta_mass_ppm) < 100) %>%
            ggplot(aes(x = retention / 60, y = delta_mass_ppm)) +
            geom_point(alpha = 0.1, color = "black", size = 1) +
            geom_hline(
              yintercept = c(10, 0, -10),
              color = "red",
              linetype = "dashed",
              linewidth = 0.2
            ) +
            labs(
              x = "Retention time (min)",
              y = "Mass error (ppm)",
              caption = "ppm = delta m/z / theoretical m/z * 1e6"
            ) +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp06"
      )
    })

    output$plot07 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_density(aes(x = peptide_length), fill = input$plot_color) +
            labs(x = "Peptide Length", y = "Frequency (%)") +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp07"
      )
    })

    output$plot08 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            ggplot(aes(x = gravy, fill = after_stat(x))) +
            geom_histogram(color = "black") +
            scale_fill_viridis_c(name = "GRAVY index", option = "C") +
            labs(
              x = NULL,
              y = "Count",
              caption = "GRAVY: hydropathic character of the sequence"
            ) +
            facet_wrap(~sample_name, ncol = 3) +
            theme(
              legend.position = "bottom",
              legend.key.width = unit(2.5, "cm"),
              legend.key.height = unit(0.25, "cm")
            )
        },
        "sp08"
      )
    })

    output$plot09 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            ggplot(aes(x = isoelectric_point, fill = after_stat(x))) +
            geom_histogram(color = "black") +
            scale_fill_viridis_c(
              name = "Isoelectric Point (pI)",
              option = "C"
            ) +
            labs(x = NULL, y = "Count") +
            facet_wrap(~sample_name, ncol = 3) +
            theme(
              legend.position = "bottom",
              legend.key.width = unit(2.5, "cm"),
              legend.key.height = unit(0.25, "cm")
            )
        },
        "sp09"
      )
    })

    output$plot10 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_bar(
              aes(x = charge),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(x = "Charge state", y = "Count") +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp10"
      )
    })

    output$plot11 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            group_by(sample_name, number_of_missed_cleavages) %>%
            summarise(n = n()) %>%
            ggplot(aes(x = number_of_missed_cleavages, y = n)) +
            geom_bar(
              stat = "identity",
              position = "dodge",
              fill = input$plot_color,
              color = "black"
            ) +
            geom_text(aes(label = n), vjust = -0.5, size = 5) +
            labs(x = "Number of Missed Cleavages", y = "Count") +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp11"
      )
    })

    output$plot12 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            dplyr::mutate(
              uniqueness = ifelse(is_unique == TRUE, "Unique", "Shared")
            ) %>%
            group_by(sample_name, uniqueness) %>%
            summarise(n = n()) %>%
            ggplot(aes(x = uniqueness, y = n)) +
            geom_bar(
              stat = "identity",
              fill = input$plot_color,
              color = "black"
            ) +
            geom_text(aes(label = n), vjust = -0.5, size = 5) +
            labs(x = "Unique peptides", y = "Count") +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp12"
      )
    })

    output$plot13 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = hyperscore),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(
              x = "Hyperscore",
              y = "Count",
              caption = "Higher values indicate greater similarity to theoretical spectra"
            ) +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp13"
      )
    })

    output$plot14 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = nextscore),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(
              x = "Nextscore",
              y = "Count",
              caption = "Second-highest scoring match for the spectrum"
            ) +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp14"
      )
    })

    output$plot15 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = probability),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(x = "PeptideProphet Probability", y = "Count") +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp15"
      )
    })

    output$plot16 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = expectation),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(x = "Expectation value", y = "Count") +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp16"
      )
    })

    output$plot17 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            tidyr::separate_rows(assigned_modifications, sep = ",") %>%
            dplyr::mutate(
              assigned_modifications = stringr::str_remove_all(
                assigned_modifications,
                ".*\\(|\\)"
              ),
              assigned_modifications = ifelse(
                is.na(assigned_modifications),
                "Unassigned",
                assigned_modifications
              )
            ) %>%
            group_by(sample_name, assigned_modifications) %>%
            summarise(n = n()) %>%
            ggplot(aes(y = assigned_modifications, x = n)) +
            geom_col(fill = input$plot_color, color = "black") +
            geom_text(aes(label = n), hjust = -0.1, size = 5) +
            scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
            labs(y = "Assigned Modifications", x = "Count") +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp17"
      )
    })

    output$plot18 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            dplyr::group_by(sample_name, entry_name) %>%
            dplyr::summarise(n_psm = n(), .groups = "drop") %>%
            dplyr::group_by(sample_name) %>%
            dplyr::slice_max(n_psm, n = 20, with_ties = FALSE) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
              entry_name = tidytext::reorder_within(
                entry_name,
                n_psm,
                sample_name
              )
            ) %>%
            ggplot(aes(x = n_psm, y = entry_name)) +
            geom_col(fill = input$plot_color, color = "black") +
            tidytext::scale_y_reordered() +
            labs(x = "Number of PSMs", y = "Protein") +
            facet_wrap(~sample_name, ncol = 3, scales = "free_y")
        },
        "sp18"
      )
    })

    output$plot19 <- renderPlot({
      rh(
        function() {
          cooccurrence_data()[["plot"]]
        },
        "sp19"
      )
    })

    output$plot20 <- renderPlot({
      rh(
        function() {
          data() %>%
            as.data.frame() %>%
            dplyr::mutate(cysteine_count = stringr::str_count(peptide, "C")) %>%
            dplyr::group_by(sample_name, cysteine_count) %>%
            summarise(n = n()) %>%
            ggplot(aes(x = cysteine_count, y = n)) +
            geom_bar(
              stat = "identity",
              fill = input$plot_color,
              color = "black"
            ) +
            geom_text(aes(label = n), vjust = -0.5, size = 5) +
            labs(x = "Cysteine counts in peptides", y = "Count") +
            facet_wrap(~sample_name, ncol = 3)
        },
        "sp20"
      )
    })

    # ════════════════════════════════════════════════════════════════════════
    # E2. SAMPLE SELECTOR FOR PROTEIN VIEW
    # ════════════════════════════════════════════════════════════════════════

    observeEvent(raw_psm_data(), {
      d <- raw_psm_data()
      req(d)
      samples <- unique(d$sample_name)
      updateSelectInput(
        session,
        "protein_sample_select",
        choices = c("All samples" = "all", setNames(samples, samples))
      )
    })

    # ════════════════════════════════════════════════════════════════════════
    # E3. FDR CURVE PLOT
    # ════════════════════════════════════════════════════════════════════════

    fdr_curve_obj <- reactive({
      d <- data()
      req(d, nrow(d) > 0, "probability" %in% names(d))

      q_sorted <- d |>
        dplyr::group_by(sample_name) |>
        dplyr::arrange(probability) |>
        dplyr::mutate(
          qvalue = 1 - probability,
          cumulative_peptides = dplyr::row_number()
        ) |>
        dplyr::ungroup()

      ggplot(q_sorted, aes(x = qvalue, y = cumulative_peptides)) +
        geom_line(color = input$plot_color, linewidth = 1) +
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
          title = "Peptide yield vs. FDR",
          x = "q-value (1 - probability)",
          y = "Cumulative distinct peptides"
        ) +
        facet_wrap(~sample_name, ncol = 3) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text = element_text(face = "bold", color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", face = "bold"),
          panel.border = element_rect(color = "black", fill = NA)
        )
    })

    output$plot_fdr_curve <- renderPlot({
      on.exit(hide_spinner("sp_fdr"), add = TRUE)
      fdr_curve_obj()
    })

    # ════════════════════════════════════════════════════════════════════════
    # E4. MODIFICATION DIAGNOSTIC
    # ════════════════════════════════════════════════════════════════════════

    mod_diag_data <- reactive({
      d <- data()
      req(d, nrow(d) > 0)
      req("assigned_modifications" %in% names(d), "retention" %in% names(d))

      d_mod <- d |>
        dplyr::filter(
          !is.na(assigned_modifications),
          assigned_modifications != ""
        ) |>
        dplyr::select(
          peptide,
          retention,
          assigned_modifications,
          sample_name,
          gene
        ) |>
        dplyr::mutate(retention = retention / 60) |>
        dplyr::rename(rt_mod = retention)

      d_unmod <- d |>
        dplyr::filter(
          is.na(assigned_modifications) | assigned_modifications == ""
        ) |>
        dplyr::group_by(peptide, sample_name) |>
        dplyr::summarise(
          rt_unmod = median(retention / 60, na.rm = TRUE),
          .groups = "drop"
        )

      paired <- dplyr::inner_join(
        d_mod,
        d_unmod,
        by = c("peptide", "sample_name")
      ) |>
        dplyr::mutate(
          delta_rt = abs(rt_mod - rt_unmod),
          classification = dplyr::if_else(
            delta_rt <= input$rt_tolerance,
            "Suspected Artifact",
            "Sample-Derived"
          )
        )
      paired
    })

    mod_diag_plot_obj <- reactive({
      paired <- mod_diag_data()
      req(paired, nrow(paired) > 0)

      top_peptides <- paired |>
        dplyr::group_by(sample_name, peptide) |>
        dplyr::summarise(
          max_delta = max(delta_rt, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::group_by(sample_name) |>
        dplyr::slice_max(order_by = max_delta, n = 30) |>
        dplyr::pull(peptide) |>
        unique()

      plot_data <- paired |> dplyr::filter(peptide %in% top_peptides)
      req(nrow(plot_data) > 0)

      ggplot(plot_data, aes(y = reorder(peptide, delta_rt))) +
        geom_segment(
          aes(
            x = rt_unmod,
            xend = rt_mod,
            yend = reorder(peptide, delta_rt),
            color = classification
          ),
          linewidth = 0.7
        ) +
        geom_point(aes(x = rt_unmod), color = "grey50", size = 2) +
        geom_point(aes(x = rt_mod, color = classification), size = 2.5) +
        scale_color_manual(
          values = c(
            "Suspected Artifact" = "#e74c3c",
            "Sample-Derived" = "#2ecc71"
          ),
          name = "Classification"
        ) +
        labs(
          title = "RT shift profile — modified vs unmodified",
          x = "Retention time (min)",
          y = "Peptide sequence",
          caption = paste0(
            "\u0394RT threshold = ",
            input$rt_tolerance,
            " min | Grey = unmodified, Colored = modified"
          )
        ) +
        facet_wrap(~sample_name, ncol = 2, scales = "free_y") +
        theme_bw() +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 7, face = "bold", color = "black"),
          axis.text.x = element_text(face = "bold", color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", face = "bold"),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "bottom"
        )
    })

    output$mod_diag_plot_ui <- renderUI({
      paired <- tryCatch(mod_diag_data(), error = function(e) NULL)
      hide_spinner("sp_moddiag")
      if (is.null(paired) || nrow(paired) == 0) {
        return(tags$p(
          style = "color:#adb5bd;text-align:center;padding:20px;",
          "No paired modified/unmodified peptides found. Ensure PSM data contains \"assigned_modifications\" and \"retention\" columns."
        ))
      }
      n_samples <- dplyr::n_distinct(paired$sample_name)
      dynamic_h <- max(400, min(n_samples * 500, 2000))
      plotOutput(ns("mod_diag_plot"), height = paste0(dynamic_h, "px"))
    })

    output$mod_diag_plot <- renderPlot({
      mod_diag_plot_obj()
    })

    output$mod_diag_table <- DT::renderDataTable({
      paired <- mod_diag_data()
      req(paired, nrow(paired) > 0)
      tbl <- paired |>
        dplyr::select(
          sample_name,
          peptide,
          assigned_modifications,
          gene,
          rt_unmod,
          rt_mod,
          delta_rt,
          classification
        ) |>
        dplyr::arrange(dplyr::desc(delta_rt))
      DT::datatable(
        tbl,
        rownames = FALSE,
        filter = "top",
        options = list(pageLength = 20, scrollX = TRUE)
      ) |>
        DT::formatRound(
          columns = c("rt_unmod", "rt_mod", "delta_rt"),
          digits = 3
        )
    })

    # ════════════════════════════════════════════════════════════════════════
    # F. PROTEIN DATA & PLOTS
    # ════════════════════════════════════════════════════════════════════════

    protein_data <- reactive({
      folder <- loaded_psm_folder()
      req(folder)
      # Auto-discover protein.tsv from the loaded PSM folder
      prot_file <- list.files(
        folder,
        pattern = "^protein\\.tsv$",
        full.names = TRUE,
        recursive = TRUE
      )
      req(length(prot_file) > 0)
      show_prot_spinners()
      readr::read_tsv(prot_file[1]) %>%
        janitor::clean_names()
    })

    combined_protein_data <- reactive({
      req(input$combined_protein)
      readr::read_tsv(input$combined_protein$datapath) %>%
        janitor::clean_names() %>%
        dplyr::select(protein_id, ends_with("max_lfq_intensity")) %>%
        column_to_rownames("protein_id") %>%
        dplyr::rename_all(~ stringr::str_remove(., "_max_lfq_intensity")) %>%
        log2()
    })

    combined_protein_raw <- reactive({
      req(input$combined_protein)
      raw <- readr::read_tsv(input$combined_protein$datapath) %>%
        janitor::clean_names()
      keep <- c(
        "entry_name",
        grep(
          "max_lfq_intensity$",
          colnames(raw),
          value = TRUE,
          ignore.case = TRUE
        ),
        grep(
          "unique_spectral_count$",
          colnames(raw),
          value = TRUE,
          ignore.case = TRUE
        )
      )
      keep <- keep[!grepl("combined", keep, ignore.case = TRUE)]
      dplyr::select(raw, dplyr::all_of(keep))
    })

    fasta_data <- reactive({
      req(input$fasta_file)
      tryCatch(
        Biostrings::readAAStringSet(input$fasta_file$datapath),
        error = function(e) {
          showNotification("Error reading FASTA file.", type = "error")
          NULL
        }
      )
    })

    observe({
      req(combined_protein_data())
      cn <- colnames(combined_protein_data())
      updateSelectInput(session, "xcol", choices = cn)
      updateSelectInput(session, "ycol", choices = cn)
    })

    observe({
      req(protein_data())
      proteins_in_data <- NULL
      if ("entry_name" %in% colnames(protein_data())) {
        proteins_in_data <- unique(protein_data()$entry_name)
      } else if ("protein_id" %in% colnames(protein_data())) {
        proteins_in_data <- unique(protein_data()$protein_id)
      } else if ("protein" %in% colnames(protein_data())) {
        proteins_in_data <- unique(protein_data()$protein)
      }
      if (!is.null(proteins_in_data) && length(proteins_in_data) > 0) {
        proteins_in_data <- sort(proteins_in_data[!is.na(proteins_in_data)])
        updateSelectizeInput(
          session,
          "selected_protein",
          choices = c("", proteins_in_data),
          server = FALSE,
          options = list(placeholder = "Type to search for a protein")
        )
      }
    })

    coverage_plot_data <- reactive({
      req(input$selected_protein)
      if (input$selected_protein == "") {
        return(list(protein_sequence = NULL))
      }
      if (is.null(input$fasta_file)) {
        return(list(error = "Please upload a FASTA file."))
      }
      if (is.null(raw_psm_data())) {
        return(list(error = "Please load PSM files in the PSM Viewer tab."))
      }
      req(data(), fasta_data())
      tp <- input$selected_protein
      tix <- grep(tp, names(fasta_data()), fixed = TRUE)[1]
      if (is.na(tix)) {
        tix <- grep(tp, names(fasta_data()), ignore.case = TRUE)[1]
      }
      if (is.na(tix)) {
        return(list(error = paste("Protein", tp, "not found in FASTA.")))
      }
      pr_seq <- as.character(fasta_data()[[tix]])
      pd <- data()
      pep_df <- if ("entry_name" %in% colnames(pd)) {
        pd[pd$entry_name == tp, ]
      } else if ("protein_id" %in% colnames(pd)) {
        pd[pd$protein_id == tp, ]
      } else if ("protein" %in% colnames(pd)) {
        pd[pd$protein == tp, ]
      } else {
        return(list(error = "No suitable protein column in PSM data."))
      }
      if (nrow(pep_df) == 0) {
        return(list(error = paste("No peptides for", tp, ".")))
      }
      sc <- if ("peptide" %in% colnames(pep_df)) {
        "peptide"
      } else if ("peptide_sequence" %in% colnames(pep_df)) {
        "peptide_sequence"
      } else {
        "sequence"
      }
      agg <- pep_df %>%
        dplyr::group_by(!!sym(sc)) %>%
        dplyr::summarise(psm_count = n(), .groups = "drop") %>%
        dplyr::rename(sequence = !!sym(sc)) %>%
        dplyr::mutate(
          start = NA_integer_,
          end = NA_integer_,
          sequence = stringr::str_remove_all(
            stringr::str_remove_all(sequence, "\\[.*?\\]"),
            "[^A-Z]"
          )
        )
      list(protein_sequence = pr_seq, peptides_data = agg)
    })

    output$protein_coverage_plot_ui <- renderUI({
      pd <- coverage_plot_data()
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
      ph <- max(400, ceiling(nchar(pd$protein_sequence) / apl) * 80 + 100)
      plotOutput(ns("protein_coverage_plot"), height = paste0(ph, "px"))
    })

    output$protein_coverage_plot <- renderPlot({
      pd <- coverage_plot_data()
      req(is.null(pd$error), pd$protein_sequence)
      apl <- max(
        10,
        min(200, ifelse(is.na(input$aa_per_line), 50, input$aa_per_line))
      )
      pl <- nchar(pd$protein_sequence)
      cov <- rep(0, pl)
      mask <- rep(FALSE, pl)
      for (i in seq_len(nrow(pd$peptides_data))) {
        pep <- pd$peptides_data[i, ]
        pos <- regexpr(pep$sequence, pd$protein_sequence, fixed = TRUE)
        if (pos == -1) {
          next
        }
        s <- as.integer(pos)
        e <- s + nchar(pep$sequence) - 1
        cnt <- ifelse(is.na(pep$psm_count), 1, pep$psm_count)
        cov[s:e] <- cov[s:e] + cnt
        mask[s:e] <- TRUE
      }
      mc <- max(cov[mask], na.rm = TRUE)
      if (!is.finite(mc) || mc == 0) {
        mc <- 1
      }
      nl <- max(1, ceiling(pl / apl))
      plots <- list()
      for (li in seq_len(nl)) {
        sp <- (li - 1) * apl + 1
        ep <- min(sp + apl - 1, pl)
        aas <- strsplit(substr(pd$protein_sequence, sp, ep), "")[[1]]
        lcov <- cov[sp:ep]
        lmask <- mask[sp:ep]
        ad <- data.frame(
          position = sp:ep,
          aa = aas,
          x_pos = seq_along(aas),
          coverage_norm = ifelse(lmask, lcov / mc, NA),
          has_peptide = lmask,
          stringsAsFactors = FALSE
        )
        cvd <- ad[ad$has_peptide, ]
        ucvd <- ad[!ad$has_peptide, ]
        p <- ggplot() +
          xlim(0.5, apl + 0.5) +
          ylim(-0.5, 2.5) +
          theme_void() +
          theme(
            plot.margin = margin(5, 5, 5, 5),
            legend.position = if (li == 1) "right" else "none",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.2, "cm")
          )
        if (nrow(ucvd) > 0) {
          p <- p +
            geom_tile(
              data = ucvd,
              aes(x = x_pos, y = 0.5),
              width = 0.9,
              height = 0.9,
              fill = "lightgray",
              color = "black",
              size = 0.2,
              alpha = 0.7
            )
        }
        if (nrow(cvd) > 0) {
          p <- p +
            geom_tile(
              data = cvd,
              aes(x = x_pos, y = 0.5, fill = coverage_norm),
              width = 0.9,
              height = 0.9,
              color = "black",
              size = 0.2
            ) +
            scale_fill_gradient(
              low = "#e1f5fe",
              high = input$plot_color,
              name = "Coverage (Norm)",
              na.value = "lightgray"
            )
        }
        if (input$show_sequence) {
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
              aes(x = x_pos, y = -0.3, label = position),
              size = 3,
              color = "black",
              fontface = "bold"
            )
        }
        plots[[li]] <- p + labs(y = paste0("AA ", sp, "-", ep))
      }
      tg <- textGrob(
        paste("Protein sequence view (length:", pl, "aa)"),
        gp = gpar(fontsize = 16, fontface = "bold")
      )
      combined <- do.call(arrangeGrob, c(plots, ncol = 1))
      grid.arrange(tg, combined, heights = c(0.5, 8))
    })

    output$aa_freq_banner <- renderUI({
      if (is.null(input$fasta_file)) {
        div(
          style = "display:flex;align-items:center;justify-content:center;height:150px;color:#6c757d;font-size:16px;border:2px dashed #dee2e6;border-radius:8px;margin:20px;",
          icon("circle-info", style = "margin-right:10px;"),
          "Upload a Protein FASTA file in the sidebar to enable amino acid frequencies comparison."
        )
      }
    })

    output$plot_aa_freq <- renderPlot({
      rh(
        function() {
          req(input$fasta_file, fasta_data(), data())
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

          # ── Expected frequencies from FASTA ──────────────────────────────
          seqs <- sapply(fasta_data(), function(x) as.character(x)[1])
          all_fasta_aa <- toupper(gsub(
            "[^A-Za-z]",
            "",
            paste0(seqs, collapse = "")
          ))
          aa_counts_fasta <- table(strsplit(all_fasta_aa, "")[[1]])
          aa_counts_fasta <- aa_counts_fasta[
            names(aa_counts_fasta) %in% base_aas
          ]
          exp_freq <- as.data.frame(
            aa_counts_fasta / sum(aa_counts_fasta, na.rm = TRUE) * 100
          )
          colnames(exp_freq) <- c("AA", "Frequency")
          exp_freq$Type <- "Expected (FASTA)"

          # ── Observed frequencies per sample ──────────────────────────────
          d <- data() %>%
            dplyr::filter(!is.na(peptide), nchar(peptide) > 0) %>%
            dplyr::mutate(
              clean_peptide = gsub("\\[.*?\\]", "", peptide)
            )

          obs_freq <- d %>%
            dplyr::group_by(sample_name) %>%
            dplyr::summarise(
              all_aas = list(
                strsplit(paste0(clean_peptide, collapse = ""), "")[[1]]
              ),
              .groups = "drop"
            ) %>%
            tidyr::unnest(all_aas) %>%
            dplyr::filter(all_aas %in% base_aas) %>%
            dplyr::count(sample_name, AA = all_aas) %>%
            dplyr::group_by(sample_name) %>%
            dplyr::mutate(Frequency = n / sum(n) * 100) %>%
            dplyr::ungroup() %>%
            dplyr::select(sample_name, AA, Frequency) %>%
            dplyr::mutate(Type = "Identified (Peptides)")

          # ── Replicate expected frequencies for every sample ───────────────
          samples <- unique(obs_freq$sample_name)
          exp_all <- purrr::map_dfr(
            samples,
            ~ dplyr::mutate(exp_freq, sample_name = .x)
          )

          comb_df <- dplyr::bind_rows(exp_all, obs_freq) %>%
            dplyr::mutate(
              Type = factor(
                Type,
                levels = c("Expected (FASTA)", "Identified (Peptides)")
              )
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
                "Identified (Peptides)" = input$plot_color
              )
            ) +
            labs(
              title = "Amino acid frequencies (expected vs identified)",
              x = "Amino acid",
              y = "Frequency (%)",
              fill = ""
            ) +
            facet_wrap(~sample_name, ncol = 3) +
            theme_bw() +
            theme(
              legend.position = "top",
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              axis.text = element_text(
                size = 10,
                face = "bold",
                color = "black"
              ),
              axis.title = element_text(size = 11, face = "bold"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold", color = "black")
            )
        },
        "sp_aa_psm"
      )
    })

    output$plot01p <- renderPlot({
      rh(
        function() {
          protein_data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = coverage),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(x = "Protein coverage (%)", y = "Count")
        },
        "sp01p"
      )
    })

    output$plot02p <- renderPlot({
      rh(
        function() {
          protein_data() %>%
            as.data.frame() %>%
            dplyr::count(organism) %>%
            ggplot() +
            geom_bar(
              aes(x = n, y = reorder(organism, n)),
              fill = input$plot_color,
              color = "black",
              stat = "identity"
            ) +
            geom_text(
              aes(x = n, y = reorder(organism, n), label = n),
              hjust = -0.1,
              size = 5
            ) +
            labs(x = "Number of proteins", y = NULL) +
            theme(axis.text.y = element_text(face = "italic"))
        },
        "sp02p"
      )
    })

    output$plot03p <- renderPlot({
      rh(
        function() {
          protein_data() %>%
            as.data.frame() %>%
            dplyr::mutate(
              protein_existence = stringr::str_remove(
                protein_existence,
                ".*\\:"
              ),
              protein_existence = factor(
                protein_existence,
                levels = c(
                  "Experimental evidence at protein level",
                  "Experimental evidence at transcript level",
                  "Protein inferred from homology",
                  "Protein predicted"
                )
              )
            ) %>%
            ggplot() +
            geom_bar(
              aes(y = protein_existence, fill = protein_existence),
              color = "black",
              show.legend = FALSE
            ) +
            scale_fill_manual(
              values = c("#5499c7", "#7fb3d5", "#a9cce3", "#d4e6f1")
            ) +
            geom_text(
              aes(y = protein_existence, label = after_stat(count)),
              stat = "count",
              vjust = -0.5,
              size = 7,
              fontface = "bold"
            ) +
            labs(y = NULL, x = "Count")
        },
        "sp03p"
      )
    })

    output$plot04p <- renderPlot({
      rh(
        function() {
          protein_data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = protein_probability),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(x = "Protein Probability", y = "Count")
        },
        "sp04p"
      )
    })

    output$plot05p <- renderPlot({
      rh(
        function() {
          protein_data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = top_peptide_probability),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(x = "Peptide Probability", y = "Count")
        },
        "sp05p"
      )
    })

    output$plot06p <- renderPlot({
      rh(
        function() {
          protein_data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = total_peptides),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(x = "Total peptides mapped to proteins", y = "Count")
        },
        "sp06p"
      )
    })

    output$plot07p <- renderPlot({
      rh(
        function() {
          protein_data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = razor_spectral_count),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(x = "Razor Spectral Count", y = "Count")
        },
        "sp07p"
      )
    })

    output$plot08p <- renderPlot({
      rh(
        function() {
          protein_data() %>%
            as.data.frame() %>%
            ggplot() +
            geom_histogram(
              aes(x = razor_intensity),
              fill = input$plot_color,
              color = "black"
            ) +
            labs(x = "Razor Intensity", y = "Count")
        },
        "sp08p"
      )
    })

    output$plot10p <- renderPlot({
      rh(
        function() {
          combined_protein_data() %>%
            as.data.frame() %>%
            rownames_to_column(var = "protein_id") %>%
            tidyr::pivot_longer(
              cols = -protein_id,
              names_to = "sample",
              values_to = "maxlfq_intensity"
            ) %>%
            ggplot() +
            geom_violin(
              aes(x = sample, y = maxlfq_intensity),
              fill = input$plot_color,
              alpha = 0.7,
              color = "black"
            ) +
            geom_boxplot(
              aes(x = sample, y = maxlfq_intensity),
              fill = "white",
              outliers = FALSE,
              color = "black",
              width = 0.1,
              show.legend = FALSE
            ) +
            labs(x = NULL, y = "log<sub>2</sub>(MaxLFQ intensity)") +
            theme(
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            )
        },
        "sp10p"
      )
    })

    output$plot_ggpairs <- renderPlot({
      rh(
        function() {
          combined_protein_data() %>%
            as.data.frame() %>%
            GGally::ggpairs(
              lower = list(continuous = wrap("points", alpha = 0.4)),
              diag = list(continuous = "barDiag"),
              upper = list(continuous = "density")
            ) +
            theme_bw() +
            theme(
              strip.background = element_blank(),
              strip.text = element_text(face = "bold")
            )
        },
        "spgg"
      )
    })

    output$plot11p <- renderPlotly({
      p <- combined_protein_data() %>%
        as.data.frame() %>%
        ggplot(aes(x = !!sym(input$xcol), y = !!sym(input$ycol))) +
        geom_point(alpha = 0.7, show.legend = FALSE) +
        geom_smooth(method = "lm", se = FALSE, color = input$plot_color) +
        labs(
          x = paste0("log2(", input$xcol, ")"),
          y = paste0("log2(", input$ycol, ")")
        ) +
        theme_bw()
      ggplotly(p)
    })

    ht <- theme(
      text = element_text(size = 15),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 1),
      legend.position = "bottom",
      legend.key.width = unit(2.5, "cm"),
      legend.key.height = unit(0.25, "cm")
    )

    output$cosine_similarity <- renderPlot({
      combined_protein_data() %>%
        as.matrix() %>%
        na.omit() %>%
        lsa::cosine() %>%
        as.data.frame() %>%
        rownames_to_column(var = "Sample") %>%
        pivot_longer(-Sample, names_to = "Match", values_to = "value") %>%
        ggplot() +
        geom_tile(aes(x = Sample, y = Match, fill = value)) +
        viridis::scale_fill_viridis(option = "E") +
        ht +
        labs(x = NULL, y = NULL, fill = "Cosine similarity")
    })
    output$euclidean_distance <- renderPlot({
      combined_protein_data() %>%
        t() %>%
        dist(method = "euclidean") %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column(var = "Sample") %>%
        pivot_longer(-Sample, names_to = "Match", values_to = "value") %>%
        ggplot() +
        geom_tile(aes(x = Sample, y = Match, fill = value)) +
        viridis::scale_fill_viridis(option = "E") +
        ht +
        labs(x = NULL, y = NULL, fill = "Euclidean distance")
    })
    output$jaccard_similarity <- renderPlot({
      combined_protein_data() %>%
        t() %>%
        vegan::vegdist(method = "jaccard", na.rm = TRUE) %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column(var = "Sample") %>%
        pivot_longer(-Sample, names_to = "Match", values_to = "value") %>%
        ggplot() +
        geom_tile(aes(x = Sample, y = Match, fill = value)) +
        viridis::scale_fill_viridis(option = "E") +
        ht +
        labs(x = NULL, y = NULL, fill = "Jaccard similarity")
    })

    output$plot_rank_lfq <- renderPlot({
      rh(
        function() {
          raw <- combined_protein_raw()
          lfq_cols <- grep(
            "max_lfq_intensity$",
            colnames(raw),
            value = TRUE,
            ignore.case = TRUE
          )
          df <- raw %>%
            dplyr::select(entry_name, dplyr::all_of(lfq_cols)) %>%
            tidyr::pivot_longer(
              cols = -entry_name,
              names_to = "sample_name",
              values_to = "log2_lfq"
            ) %>%
            dplyr::mutate(
              sample_name = stringr::str_remove(
                sample_name,
                "_max_lfq_intensity$"
              ),
              log2_lfq = log2(as.numeric(log2_lfq))
            ) %>%
            dplyr::filter(!is.na(log2_lfq) & is.finite(log2_lfq))

          samples <- sort(unique(df$sample_name))
          plots <- lapply(samples, function(s) {
            d <- df %>%
              dplyr::filter(sample_name == s) %>%
              dplyr::slice_max(log2_lfq, n = 20, with_ties = FALSE) %>%
              dplyr::mutate(
                entry_name = tidytext::reorder_within(
                  entry_name,
                  log2_lfq,
                  sample_name
                )
              )
            ggplot(d, aes(x = log2_lfq, y = entry_name)) +
              geom_col(
                fill = input$plot_color,
                color = "black",
                linewidth = 0.3
              ) +
              tidytext::scale_y_reordered() +
              labs(title = s, x = "log2(MaxLFQ intensity)", y = NULL) +
              theme_bw() +
              theme(
                axis.title = element_text(size = 10, face = "bold"),
                axis.text = element_text(
                  size = 8,
                  face = "bold",
                  color = "black"
                ),
                plot.title = element_text(face = "bold", hjust = 0.5)
              )
          })
          patchwork::wrap_plots(plots, ncol = min(3L, length(plots))) +
            patchwork::plot_annotation(
              title = "Top 20 proteins by MaxLFQ intensity (log2)",
              theme = theme(
                plot.title = element_text(
                  size = 14,
                  face = "bold",
                  hjust = 0.5
                )
              )
            )
        },
        "sp_rank_lfq"
      )
    })

    output$plot_rank_usc <- renderPlot({
      rh(
        function() {
          raw <- combined_protein_raw()
          usc_cols <- grep(
            "unique_spectral_count$",
            colnames(raw),
            value = TRUE,
            ignore.case = TRUE
          )
          df <- raw %>%
            dplyr::select(entry_name, dplyr::all_of(usc_cols)) %>%
            tidyr::pivot_longer(
              cols = -entry_name,
              names_to = "sample_name",
              values_to = "unique_spectral_count"
            ) %>%
            dplyr::mutate(
              sample_name = stringr::str_remove(
                sample_name,
                "_unique_spectral_count$"
              ),
              unique_spectral_count = as.numeric(unique_spectral_count)
            ) %>%
            dplyr::filter(
              !is.na(unique_spectral_count) & unique_spectral_count > 0
            )

          samples <- sort(unique(df$sample_name))
          plots <- lapply(samples, function(s) {
            d <- df %>%
              dplyr::filter(sample_name == s) %>%
              dplyr::slice_max(
                unique_spectral_count,
                n = 20,
                with_ties = FALSE
              ) %>%
              dplyr::mutate(
                entry_name = tidytext::reorder_within(
                  entry_name,
                  unique_spectral_count,
                  sample_name
                )
              )
            ggplot(d, aes(x = unique_spectral_count, y = entry_name)) +
              geom_col(
                fill = input$plot_color,
                color = "black",
                linewidth = 0.3
              ) +
              tidytext::scale_y_reordered() +
              labs(title = s, x = "Unique spectral count", y = NULL) +
              theme_bw() +
              theme(
                axis.text = element_text(
                  size = 8,
                  face = "bold",
                  color = "black"
                ),
                axis.title = element_text(size = 10, face = "bold"),
                plot.title = element_text(face = "bold", hjust = 0.5)
              )
          })
          patchwork::wrap_plots(plots, ncol = min(3L, length(plots))) +
            patchwork::plot_annotation(
              title = "Top 20 proteins by unique spectral count",
              theme = theme(
                plot.title = element_text(
                  size = 14,
                  face = "bold",
                  hjust = 0.5
                )
              )
            )
        },
        "sp_rank_usc"
      )
    })

    # ── Download ──────────────────────────────────────────────────────────────
    output$download_all_plots <- downloadHandler(
      filename = function() {
        paste0("PSM_and_Protein_plots_", Sys.Date(), ".zip")
      },
      content = function(file) {
        td <- tempdir()
        fps <- character()
        save_p <- function(nm, fn, w = 10, h = 6) {
          fp <- file.path(td, nm)
          tryCatch(
            {
              p <- fn()
              if (!is.null(p)) {
                ggsave(fp, p, width = w, height = h, dpi = 300, bg = "white")
                fps <<- c(fps, fp)
              }
            },
            error = function(e) message("Skip ", nm, ": ", conditionMessage(e))
          )
        }
        withProgress(message = "Saving plots...", value = 0, {
          # ── PSM plots ────────────────────────────────────────────────────
          if (!is.null(tryCatch(data(), error = function(e) NULL))) {
            d_psm <- data()
            try({
              # plot01 — Protease fingerprint heatmap
              save_p(
                "PSM_plot01_protease_fingerprint.png",
                function() {
                  purrr::imap_dfr(
                    frequency_matrix_of_aa(),
                    ~ as.data.frame(.x) %>%
                      rownames_to_column(var = "residue") %>%
                      pivot_longer(
                        -residue,
                        names_to = "position",
                        values_to = "frequency"
                      ) %>%
                      dplyr::mutate(sample_name = .y)
                  ) %>%
                    dplyr::mutate(
                      position = factor(
                        position,
                        c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
                      ),
                      residue = factor(
                        residue,
                        c(
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
                      )
                    ) %>%
                    ggplot(aes(x = position, y = residue, fill = frequency)) +
                    geom_tile(color = "black") +
                    scale_fill_gradient(low = "#d4e6f1", high = "#154360") +
                    geom_vline(
                      xintercept = 4.5,
                      color = "black",
                      linetype = "dashed"
                    ) +
                    theme_minimal() +
                    facet_wrap(~sample_name, ncol = 3) +
                    labs(
                      title = "Cleavage Site Specificity",
                      x = "Position",
                      y = "Amino acid residue",
                      fill = "Frequency (%)"
                    )
                },
                w = 14,
                h = 8
              )

              # plot03 — N-term seqLogo
              save_p(
                "PSM_plot03_nterm_seqlogo.png",
                function() {
                  seq_lst <- d_psm %>%
                    dplyr::filter(!is.na(fingerprint_Nterm)) %>%
                    dplyr::group_by(sample_name) %>%
                    dplyr::summarise(
                      seqs = list(fingerprint_Nterm),
                      .groups = "drop"
                    ) %>%
                    tibble::deframe()
                  ggseqlogo::ggseqlogo(
                    seq_lst,
                    method = "bits",
                    seq_type = "AA"
                  ) +
                    seqlogo_scale +
                    theme_bw() +
                    labs(
                      title = "SeqLogo of the N-termini fingerprint",
                      x = "Amino acid position",
                      y = "Bits"
                    )
                },
                w = 12,
                h = 5
              )

              # plot04 — C-term seqLogo
              save_p(
                "PSM_plot04_cterm_seqlogo.png",
                function() {
                  seq_lst <- d_psm %>%
                    dplyr::filter(!is.na(fingerprint_Cterm)) %>%
                    dplyr::group_by(sample_name) %>%
                    dplyr::summarise(
                      seqs = list(fingerprint_Cterm),
                      .groups = "drop"
                    ) %>%
                    tibble::deframe()
                  ggseqlogo::ggseqlogo(
                    seq_lst,
                    method = "bits",
                    seq_type = "AA"
                  ) +
                    seqlogo_scale +
                    theme_bw() +
                    labs(
                      title = "SeqLogo of the C-termini fingerprint",
                      x = "Amino acid position",
                      y = "Bits"
                    )
                },
                w = 12,
                h = 5
              )

              # plot05 — m/z vs retention time
              save_p(
                "PSM_plot05_mz_rt.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    ggplot(aes(x = retention / 60, y = observed_m_z)) +
                    ggpointdensity::geom_pointdensity(size = 0.25) +
                    viridis::scale_color_viridis(option = "plasma") +
                    labs(
                      x = "Retention time (min)",
                      y = "Scan range (m/z)",
                      color = "Neighborhoods"
                    ) +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot06 — mass error
              save_p(
                "PSM_plot06_mass_error.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    dplyr::filter(abs(delta_mass_ppm) < 100) %>%
                    ggplot(aes(x = retention / 60, y = delta_mass_ppm)) +
                    geom_point(alpha = 0.1, color = "black", size = 1) +
                    geom_hline(
                      yintercept = c(10, 0, -10),
                      color = "red",
                      linetype = "dashed",
                      linewidth = 0.2
                    ) +
                    labs(x = "Retention time (min)", y = "Mass error (ppm)") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot07 — peptide length
              save_p(
                "PSM_plot07_peptide_length.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    ggplot() +
                    geom_density(
                      aes(x = peptide_length),
                      fill = input$plot_color,
                      color = "black"
                    ) +
                    labs(x = "Peptide Length", y = "Frequency (%)") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot08 — GRAVY
              save_p(
                "PSM_plot08_gravy.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    ggplot(aes(x = gravy)) +
                    geom_histogram(fill = input$plot_color, color = "black") +
                    labs(x = "GRAVY index", y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot09 — isoelectric point
              save_p(
                "PSM_plot09_isoelectric_point.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    ggplot(aes(x = isoelectric_point)) +
                    geom_histogram(fill = input$plot_color, color = "black") +
                    labs(x = NULL, y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot10 — charge state
              save_p(
                "PSM_plot10_charge_state.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    ggplot() +
                    geom_bar(
                      aes(x = charge),
                      fill = input$plot_color,
                      color = "black"
                    ) +
                    labs(x = "Charge state", y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot11 — missed cleavages
              save_p(
                "PSM_plot11_missed_cleavages.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    dplyr::count(sample_name, number_of_missed_cleavages) %>%
                    dplyr::mutate(
                      number_of_missed_cleavages = factor(
                        number_of_missed_cleavages
                      )
                    ) %>%
                    ggplot(aes(x = number_of_missed_cleavages, y = n)) +
                    geom_col(fill = input$plot_color, color = "black") +
                    labs(x = "Number of Missed Cleavages", y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot12 — unique peptides
              save_p(
                "PSM_plot12_uniqueness.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    dplyr::count(sample_name, is_unique) %>%
                    dplyr::mutate(
                      uniqueness = ifelse(is_unique, "Unique", "Shared")
                    ) %>%
                    ggplot(aes(x = uniqueness, y = n)) +
                    geom_col(fill = input$plot_color, color = "black") +
                    labs(x = "Unique peptides", y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot13 — hyperscore
              save_p(
                "PSM_plot13_hyperscore.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    ggplot() +
                    geom_histogram(
                      aes(x = hyperscore),
                      fill = input$plot_color,
                      color = "black"
                    ) +
                    labs(x = "Hyperscore", y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot14 — nextscore
              save_p(
                "PSM_plot14_nextscore.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    ggplot() +
                    geom_histogram(
                      aes(x = nextscore),
                      fill = input$plot_color,
                      color = "black"
                    ) +
                    labs(x = "Nextscore", y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot15 — PeptideProphet probability
              save_p(
                "PSM_plot15_probability.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    ggplot() +
                    geom_histogram(
                      aes(x = probability),
                      fill = input$plot_color,
                      color = "black"
                    ) +
                    labs(x = "PeptideProphet Probability", y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot16 — expectation value
              save_p(
                "PSM_plot16_expectation.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    ggplot() +
                    geom_histogram(
                      aes(x = expectation),
                      fill = input$plot_color,
                      color = "black"
                    ) +
                    labs(x = "Expectation value", y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot17 — assigned modifications
              save_p(
                "PSM_plot17_assigned_mods.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    tidyr::separate_rows(assigned_modifications, sep = ",") %>%
                    dplyr::mutate(
                      assigned_modifications = stringr::str_remove_all(
                        assigned_modifications,
                        ".*\\(|\\)"
                      ),
                      assigned_modifications = ifelse(
                        is.na(assigned_modifications),
                        "Unassigned",
                        assigned_modifications
                      )
                    ) %>%
                    dplyr::count(sample_name, assigned_modifications) %>%
                    ggplot(aes(y = assigned_modifications, x = n)) +
                    geom_col(fill = input$plot_color, color = "black") +
                    labs(y = "Assigned Modifications", x = "Count") +
                    facet_wrap(~sample_name, ncol = 3, scales = "free_y")
                },
                w = 12,
                h = 8
              )

              # plot18 — top proteins by PSM count
              save_p(
                "PSM_plot18_top_proteins.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    dplyr::group_by(sample_name, entry_name) %>%
                    dplyr::summarise(n_psm = n(), .groups = "drop") %>%
                    dplyr::group_by(sample_name) %>%
                    dplyr::slice_max(n_psm, n = 20, with_ties = FALSE) %>%
                    dplyr::mutate(
                      entry_name = tidytext::reorder_within(
                        entry_name,
                        n_psm,
                        sample_name
                      )
                    ) %>%
                    ggplot(aes(x = n_psm, y = entry_name)) +
                    geom_col(fill = input$plot_color, color = "black") +
                    tidytext::scale_y_reordered() +
                    labs(x = "Number of PSMs", y = "Protein") +
                    facet_wrap(~sample_name, ncol = 3, scales = "free_y")
                },
                w = 14,
                h = 8
              )

              # plot19 — co-occurrence matrix
              save_p(
                "PSM_plot19_cooccurrence.png",
                function() {
                  cooccurrence_data()[["plot"]]
                },
                w = 12,
                h = 10
              )

              # plot20 — cysteine counts
              save_p(
                "PSM_plot20_cysteine.png",
                function() {
                  d_psm %>%
                    as.data.frame() %>%
                    dplyr::mutate(
                      cysteine_count = stringr::str_count(peptide, "C")
                    ) %>%
                    dplyr::group_by(sample_name, cysteine_count) %>%
                    summarise(n = n()) %>%
                    ggplot(aes(x = cysteine_count, y = n)) +
                    geom_bar(
                      stat = "identity",
                      fill = input$plot_color,
                      color = "black"
                    ) +
                    geom_text(aes(label = n), vjust = -0.5, size = 5) +
                    labs(x = "Cysteine counts in peptides", y = "Count") +
                    facet_wrap(~sample_name, ncol = 3)
                },
                w = 12,
                h = 7
              )

              # plot_aa_freq — amino acid frequency vs FASTA background
              if (
                !is.null(input$fasta_file) &&
                  !is.null(tryCatch(fasta_data(), error = function(e) NULL))
              ) {
                save_p(
                  "PSM_plot_aa_frequency.png",
                  function() {
                    fd <- tryCatch(fasta_data(), error = function(e) NULL)
                    if (is.null(fd)) {
                      return(NULL)
                    }
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
                    seqs <- sapply(fd, function(x) as.character(x)[1])
                    all_fasta_aa <- toupper(gsub(
                      "[^A-Za-z]",
                      "",
                      paste0(seqs, collapse = "")
                    ))
                    aa_counts_fasta <- table(strsplit(all_fasta_aa, "")[[1]])
                    aa_counts_fasta <- aa_counts_fasta[
                      names(aa_counts_fasta) %in% base_aas
                    ]
                    exp_freq <- as.data.frame(
                      aa_counts_fasta / sum(aa_counts_fasta, na.rm = TRUE) * 100
                    )
                    colnames(exp_freq) <- c("AA", "Frequency")
                    exp_freq$Type <- "Expected (FASTA)"
                    d_aa <- d_psm %>%
                      dplyr::filter(!is.na(peptide), nchar(peptide) > 0) %>%
                      dplyr::mutate(
                        clean_peptide = gsub("\\[.*?\\]", "", peptide)
                      )
                    obs_freq <- d_aa %>%
                      dplyr::group_by(sample_name) %>%
                      dplyr::summarise(
                        all_aas = list(
                          strsplit(paste0(clean_peptide, collapse = ""), "")[[
                            1
                          ]]
                        ),
                        .groups = "drop"
                      ) %>%
                      tidyr::unnest(all_aas) %>%
                      dplyr::filter(all_aas %in% base_aas) %>%
                      dplyr::count(sample_name, AA = all_aas) %>%
                      dplyr::group_by(sample_name) %>%
                      dplyr::mutate(Frequency = n / sum(n) * 100) %>%
                      dplyr::ungroup() %>%
                      dplyr::select(sample_name, AA, Frequency) %>%
                      dplyr::mutate(Type = "Identified (Peptides)")
                    samples_aa <- unique(obs_freq$sample_name)
                    exp_all <- purrr::map_dfr(
                      samples_aa,
                      ~ dplyr::mutate(exp_freq, sample_name = .x)
                    )
                    comb_df <- dplyr::bind_rows(exp_all, obs_freq) %>%
                      dplyr::mutate(
                        Type = factor(
                          Type,
                          levels = c(
                            "Expected (FASTA)",
                            "Identified (Peptides)"
                          )
                        )
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
                          "Identified (Peptides)" = input$plot_color
                        )
                      ) +
                      labs(
                        title = "Amino acid frequencies (expected vs identified)",
                        x = "Amino acid",
                        y = "Frequency (%)",
                        fill = ""
                      ) +
                      facet_wrap(~sample_name, ncol = 3) +
                      theme_bw() +
                      theme(
                        legend.position = "top",
                        plot.title = element_text(
                          size = 14,
                          face = "bold",
                          hjust = 0.5
                        ),
                        axis.text = element_text(
                          size = 10,
                          face = "bold",
                          color = "black"
                        ),
                        axis.title = element_text(size = 11, face = "bold"),
                        strip.background = element_blank(),
                        strip.text = element_text(
                          face = "bold",
                          color = "black"
                        )
                      )
                  },
                  w = 12,
                  h = 7
                )
              }

              # FDR curve
              tryCatch(
                {
                  save_p(
                    "PSM_FDR_curve.png",
                    function() fdr_curve_obj(),
                    w = 11,
                    h = 8
                  )
                },
                error = function(e) NULL
              )
            })
          }
          incProgress(1 / 3)

          # ── Protein plots ─────────────────────────────────────────────────
          if (!is.null(tryCatch(protein_data(), error = function(e) NULL))) {
            d_prot <- protein_data()
            try({
              save_p("Protein_plot01_coverage.png", function() {
                d_prot %>%
                  as.data.frame() %>%
                  ggplot() +
                  geom_histogram(
                    aes(x = coverage),
                    fill = input$plot_color,
                    color = "black"
                  ) +
                  labs(x = "Protein coverage (%)", y = "Count")
              })
              save_p("Protein_plot02_organisms.png", function() {
                d_prot %>%
                  as.data.frame() %>%
                  dplyr::count(organism) %>%
                  ggplot() +
                  geom_bar(
                    aes(x = n, y = reorder(organism, n)),
                    fill = input$plot_color,
                    color = "black",
                    stat = "identity"
                  ) +
                  geom_text(
                    aes(x = n, y = reorder(organism, n), label = n),
                    hjust = -0.1,
                    size = 5
                  ) +
                  labs(x = "Number of proteins", y = NULL)
              })
              save_p("Protein_plot03_existence.png", function() {
                d_prot %>%
                  as.data.frame() %>%
                  dplyr::mutate(
                    protein_existence = stringr::str_remove(
                      protein_existence,
                      ".*\\:"
                    )
                  ) %>%
                  ggplot() +
                  geom_bar(
                    aes(y = protein_existence, fill = input$plot_color),
                    color = "black",
                    show.legend = FALSE
                  ) +
                  labs(x = "Count", y = "Protein existence")
              })
              save_p("Protein_plot04_probability.png", function() {
                d_prot %>%
                  as.data.frame() %>%
                  ggplot() +
                  geom_histogram(
                    aes(x = protein_probability),
                    fill = input$plot_color,
                    color = "black"
                  ) +
                  labs(x = "Protein probability", y = "Count")
              })
              save_p("Protein_plot05_peptide_prob.png", function() {
                d_prot %>%
                  as.data.frame() %>%
                  ggplot() +
                  geom_histogram(
                    aes(x = top_peptide_probability),
                    fill = input$plot_color,
                    color = "black"
                  ) +
                  labs(x = "Top peptide probability", y = "Count")
              })
              save_p("Protein_plot06_total_peptides.png", function() {
                d_prot %>%
                  as.data.frame() %>%
                  ggplot() +
                  geom_histogram(
                    aes(x = total_peptides),
                    fill = input$plot_color,
                    color = "black"
                  ) +
                  labs(x = "Total peptides", y = "Count")
              })
              save_p("Protein_plot07_spectral_count.png", function() {
                d_prot %>%
                  as.data.frame() %>%
                  ggplot() +
                  geom_histogram(
                    aes(x = razor_spectral_count),
                    fill = input$plot_color,
                    color = "black"
                  ) +
                  labs(x = "Razor Spectral Count", y = "Count")
              })
              save_p("Protein_plot08_razor_intensity.png", function() {
                d_prot %>%
                  as.data.frame() %>%
                  ggplot() +
                  geom_histogram(
                    aes(x = razor_intensity),
                    fill = input$plot_color,
                    color = "black"
                  ) +
                  labs(x = "Razor Intensity", y = "Count")
              })
            })
          }
          incProgress(1 / 3)

          # ── Combined protein plots ────────────────────────────────────────
          if (!is.null(input$combined_protein)) {
            d_comb <- combined_protein_data()
            raw <- combined_protein_raw()
            try({
              save_p(
                "Protein_plot10_maxlfq_distribution.png",
                function() {
                  d_comb %>%
                    as.data.frame() %>%
                    rownames_to_column(var = "protein_id") %>%
                    tidyr::pivot_longer(
                      -protein_id,
                      names_to = "sample",
                      values_to = "log2_lfq"
                    ) %>%
                    ggplot() +
                    geom_violin(
                      aes(x = sample, y = log2_lfq),
                      fill = input$plot_color,
                      color = "black"
                    ) +
                    labs(x = "Sample", y = "log2(MaxLFQ intensity)") +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                },
                w = 10,
                h = 6
              )

              save_p(
                "Protein_plot_ggpairs.png",
                function() {
                  d_comb %>%
                    GGally::ggpairs(
                      lower = list(continuous = wrap("points", alpha = 0.4)),
                      diag = list(continuous = "barDiag"),
                      upper = list(continuous = "density")
                    ) +
                    theme_bw()
                },
                w = 12,
                h = 12
              )

              ht_dl <- theme(
                text = element_text(size = 15),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                legend.position = "bottom",
                legend.key.width = unit(2.5, "cm"),
                legend.key.height = unit(0.25, "cm")
              )
              save_p("Protein_cosine_similarity.png", function() {
                d_comb %>%
                  as.matrix() %>%
                  na.omit() %>%
                  lsa::cosine() %>%
                  as.data.frame() %>%
                  rownames_to_column("Sample") %>%
                  pivot_longer(
                    -Sample,
                    names_to = "Match",
                    values_to = "value"
                  ) %>%
                  ggplot() +
                  geom_tile(aes(x = Sample, y = Match, fill = value)) +
                  viridis::scale_fill_viridis(option = "E") +
                  ht_dl +
                  labs(x = NULL, y = NULL, fill = "Cosine similarity")
              })
              save_p("Protein_euclidean_distance.png", function() {
                d_comb %>%
                  t() %>%
                  dist(method = "euclidean") %>%
                  as.matrix() %>%
                  as.data.frame() %>%
                  rownames_to_column("Sample") %>%
                  pivot_longer(
                    -Sample,
                    names_to = "Match",
                    values_to = "value"
                  ) %>%
                  ggplot() +
                  geom_tile(aes(x = Sample, y = Match, fill = value)) +
                  viridis::scale_fill_viridis(option = "E") +
                  ht_dl +
                  labs(x = NULL, y = NULL, fill = "Euclidean distance")
              })
              save_p("Protein_jaccard_similarity.png", function() {
                d_comb %>%
                  t() %>%
                  vegan::vegdist(method = "jaccard", na.rm = TRUE) %>%
                  as.matrix() %>%
                  as.data.frame() %>%
                  rownames_to_column("Sample") %>%
                  pivot_longer(
                    -Sample,
                    names_to = "Match",
                    values_to = "value"
                  ) %>%
                  ggplot() +
                  geom_tile(aes(x = Sample, y = Match, fill = value)) +
                  viridis::scale_fill_viridis(option = "E") +
                  ht_dl +
                  labs(x = NULL, y = NULL, fill = "Jaccard similarity")
              })

              # Rank plots — reuse same logic as renderPlot
              save_p(
                "Protein_rank_lfq.png",
                function() {
                  lfq_cols <- grep(
                    "max_lfq_intensity$",
                    colnames(raw),
                    value = TRUE,
                    ignore.case = TRUE
                  )
                  df <- raw %>%
                    dplyr::select(entry_name, dplyr::all_of(lfq_cols)) %>%
                    tidyr::pivot_longer(
                      -entry_name,
                      names_to = "sample_name",
                      values_to = "log2_lfq"
                    ) %>%
                    dplyr::mutate(
                      sample_name = stringr::str_remove(
                        sample_name,
                        "_max_lfq_intensity$"
                      ),
                      log2_lfq = log2(as.numeric(log2_lfq))
                    ) %>%
                    dplyr::filter(!is.na(log2_lfq) & is.finite(log2_lfq))
                  plots <- lapply(sort(unique(df$sample_name)), function(s) {
                    d <- df %>%
                      dplyr::filter(sample_name == s) %>%
                      dplyr::slice_max(log2_lfq, n = 20, with_ties = FALSE) %>%
                      dplyr::mutate(
                        entry_name = tidytext::reorder_within(
                          entry_name,
                          log2_lfq,
                          sample_name
                        )
                      )
                    ggplot(d, aes(x = log2_lfq, y = entry_name)) +
                      geom_col(
                        fill = input$plot_color,
                        color = "black",
                        linewidth = 0.3
                      ) +
                      tidytext::scale_y_reordered() +
                      labs(title = s, x = "log2(MaxLFQ intensity)", y = NULL) +
                      theme_bw()
                  })
                  patchwork::wrap_plots(plots, ncol = min(3L, length(plots))) +
                    patchwork::plot_annotation(
                      title = "Top 20 proteins by MaxLFQ intensity (log2)"
                    )
                },
                w = 14,
                h = 8
              )

              save_p(
                "Protein_rank_usc.png",
                function() {
                  usc_cols <- grep(
                    "unique_spectral_count$",
                    colnames(raw),
                    value = TRUE,
                    ignore.case = TRUE
                  )
                  df <- raw %>%
                    dplyr::select(entry_name, dplyr::all_of(usc_cols)) %>%
                    tidyr::pivot_longer(
                      -entry_name,
                      names_to = "sample_name",
                      values_to = "unique_spectral_count"
                    ) %>%
                    dplyr::mutate(
                      sample_name = stringr::str_remove(
                        sample_name,
                        "_unique_spectral_count$"
                      ),
                      unique_spectral_count = as.numeric(unique_spectral_count)
                    ) %>%
                    dplyr::filter(
                      !is.na(unique_spectral_count) & unique_spectral_count > 0
                    )
                  plots <- lapply(sort(unique(df$sample_name)), function(s) {
                    d <- df %>%
                      dplyr::filter(sample_name == s) %>%
                      dplyr::slice_max(
                        unique_spectral_count,
                        n = 20,
                        with_ties = FALSE
                      ) %>%
                      dplyr::mutate(
                        entry_name = tidytext::reorder_within(
                          entry_name,
                          unique_spectral_count,
                          sample_name
                        )
                      )
                    ggplot(d, aes(x = unique_spectral_count, y = entry_name)) +
                      geom_col(
                        fill = input$plot_color,
                        color = "black",
                        linewidth = 0.3
                      ) +
                      tidytext::scale_y_reordered() +
                      labs(title = s, x = "Unique spectral count", y = NULL) +
                      theme_bw()
                  })
                  patchwork::wrap_plots(plots, ncol = min(3L, length(plots))) +
                    patchwork::plot_annotation(
                      title = "Top 20 proteins by unique spectral count"
                    )
                },
                w = 14,
                h = 8
              )
            })
          }

          # ── Modification Diagnostic ───────────────────────────────────────
          try({
            save_p(
              "PSM_Mod_Diagnostic_Dumbbell.png",
              function() mod_diag_plot_obj(),
              w = 14,
              h = 10
            )
          })
          try({
            paired <- mod_diag_data()
            if (!is.null(paired) && nrow(paired) > 0) {
              rt_csv_path <- file.path(td, "PSM_RT_Shift_Table.csv")
              readr::write_csv(paired, rt_csv_path)
              fps <<- c(fps, rt_csv_path)
            }
          })
          incProgress(1 / 3)
        })
        utils::zip(file, files = fps, flags = "-j")
      },
      contentType = "application/zip"
    )
  })
}
