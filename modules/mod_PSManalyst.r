## PSManalyst Shiny Module — PSManalyst module for proteOmni

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
        fp <- strsplit(c(.x$fingerprint_Nterm, .x$fingerprint_Cterm), "")

        if (length(fp) == 0) {
          return(NULL)
        }

        mat <- matrix(unlist(fp), ncol = 8, byrow = TRUE)
        mat <- mat[
          !apply(mat, 1, function(x) any(x %in% c("B", "X", "Z", "U", "O"))),
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

  purrr::compact(stats::setNames(mats, sample_names))
}

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

# Recursively find all psm.tsv files inside a directory tree.
find_psm_files <- function(root_path) {
  list.files(
    path = root_path,
    pattern = "^psm\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
}

# Read one psm.tsv and add a sample_name column derived from the Spectrum column.
read_single_psm <- function(path) {
  df <- suppressMessages(
    readr::read_tsv(path, col_types = readr::cols(.default = "c"))
  ) %>%
    janitor::clean_names()

  if (nrow(df) == 0) {
    empty_line <- as.list(rep(NA_character_, ncol(df)))
    names(empty_line) <- names(df)
    df <- dplyr::bind_rows(df, empty_line) %>%
      dplyr::mutate(sample_name = basename(dirname(path)))
  } else {
    df <- df %>%
      dplyr::mutate(sample_name = stringr::str_extract(spectrum, "^[^.]+"))
  }

  # Experiment folder: links PSM sample names (spectrum-file prefixes) to the
  # per-folder tables (peptide.tsv), whose samples are named after the folder.
  df$experiment <- basename(dirname(path))

  df
}

# Read and row-bind all psm.tsv files found under \code{root_path}.
read_all_psm_files <- function(root_path) {
  paths <- find_psm_files(root_path)
  if (length(paths) == 0L) {
    stop("No psm.tsv files found under: ", root_path)
  }
  purrr::map_dfr(paths, read_single_psm) %>%
    readr::type_convert(na = c("", "NA"))
}


# ── iBAQ from FragPipe peptide.tsv ──────────────────────────────────────
#   iBAQ = summed peptide intensity of a protein in an experiment /
#          number of theoretically observable peptides of that protein
# Same definition as the DIA-NN iBAQ in mod_QC4DIANN; the in-silico digest
# (count_theoretical_peptides, fasta_accessions) lives in modules/utils_fasta.r.

# Read every peptide.tsv under root_path (they sit next to the psm.tsv files).
# One row per peptide and experiment; `experiment` is the parent folder name.
read_all_peptide_files <- function(root_path) {
  paths <- list.files(
    path = root_path,
    pattern = "^peptide\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(paths) == 0L) {
    stop("No peptide.tsv files found under: ", root_path)
  }
  required <- c("peptide", "intensity", "protein", "protein_id")
  optional <- c("entry_name", "gene", "mapped_proteins")

  purrr::map_dfr(paths, function(f) {
    df <- data.table::fread(
      f,
      sep = "\t",
      header = TRUE,
      quote = "",
      colClasses = "character",
      data.table = FALSE
    ) %>%
      janitor::clean_names()
    missing_cols <- setdiff(required, names(df))
    if (length(missing_cols) > 0) {
      stop(
        "Missing required columns in ",
        f,
        ": ",
        paste(missing_cols, collapse = ", ")
      )
    }
    for (col in setdiff(optional, names(df))) {
      df[[col]] <- rep(NA_character_, nrow(df))
    }
    df <- df[, c(required, optional), drop = FALSE]
    df$intensity <- suppressWarnings(as.numeric(df$intensity))
    df$experiment <- rep(basename(dirname(f)), nrow(df))
    df
  })
}

# Named sequences (first word of the FASTA header -> sequence) from an
# AAStringSet. FragPipe's `Protein` column holds that same first word
# (e.g. "sp|P04075|ALDOA_HUMAN", "contam_sp|P02769|ALBU_BOVIN").
fragpipe_fasta_sequences <- function(fasta) {
  seqs <- toupper(gsub("[[:space:]]", "", as.character(fasta)))
  names(seqs) <- sub("\\s.*$", "", trimws(names(fasta)))
  seqs[!duplicated(names(seqs))]
}

# Digestion settings of the search, from the fragger.params FragPipe writes
# into the results folder. NULL when the file is absent; individual fields are
# NA when a key is missing. Only the first enzyme defines the cleavage rule;
# `enzyme2` is its name when a second enzyme was used (otherwise NA).
read_fragger_digest_params <- function(root_path) {
  f <- file.path(root_path, "fragger.params")
  if (!file.exists(f)) {
    found <- list.files(
      root_path,
      pattern = "^fragger\\.params$",
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(found) == 0L) {
      return(NULL)
    }
    f <- found[1]
  }
  lines <- sub("#.*$", "", readLines(f, warn = FALSE))
  kv <- regmatches(
    lines,
    regexec("^\\s*([A-Za-z0-9_]+)\\s*=\\s*(.*?)\\s*$", lines, perl = TRUE)
  )
  kv <- kv[lengths(kv) == 3L]
  vals <- stats::setNames(
    vapply(kv, `[`, character(1), 3L),
    vapply(kv, `[`, character(1), 2L)
  )
  chr <- function(k) if (k %in% names(vals)) vals[[k]] else NA_character_
  int <- function(k) suppressWarnings(as.integer(chr(k)))
  enzyme2 <- chr("search_enzyme_name_2")
  if (is.na(enzyme2) || enzyme2 %in% c("", "null")) {
    enzyme2 <- NA_character_
  }
  list(
    file = f,
    enzyme = chr("search_enzyme_name_1"),
    cut = chr("search_enzyme_cut_1"),
    nocut = chr("search_enzyme_nocut_1"),
    sense = chr("search_enzyme_sense_1"),
    missed_cleavages = int("allowed_missed_cleavage_1"),
    min_len = int("digest_min_length"),
    max_len = int("digest_max_length"),
    enzyme2 = enzyme2
  )
}

# Cleavage-site regex for count_theoretical_peptides() from an MSFragger-style
# enzyme rule: residues that are cut, residues that block the cut, and the side
# of the cut residue that is cleaved ("C" = after it, "N" = before it).
# NULL when no cut residue is given (e.g. a nonspecific search).
protease_regex_from_rule <- function(cut, nocut = "", sense = "C") {
  residues <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) {
      return("")
    }
    gsub("[^A-Z]", "", toupper(x))
  }
  cut <- residues(cut)
  nocut <- residues(nocut)
  if (!nzchar(cut)) {
    return(NULL)
  }
  if (identical(toupper(sense), "N")) {
    # (?<=.) keeps strsplit() from splitting at every character on a pure
    # lookahead, and from emitting an empty first fragment.
    paste0(
      "(?<=.)",
      if (nzchar(nocut)) sprintf("(?<![%s])", nocut),
      sprintf("(?=[%s])", cut)
    )
  } else {
    paste0(
      sprintf("(?<=[%s])", cut),
      if (nzchar(nocut)) sprintf("(?![%s])", nocut)
    )
  }
}

# iBAQ per protein and experiment, long format (untransformed).
#   peptides  : output of read_all_peptide_files().
#   sequences : output of fragpipe_fasta_sequences().
# Each peptide counts towards the protein FragPipe assigned it to (razor
# assignment); unique_only keeps peptides that map to no other protein.
# Proteins with no quantified peptide in an experiment are absent from it.
# Proteins not found in the FASTA are dropped and listed in the
# "unmatched_proteins" attribute. exclude_contaminants drops proteins whose
# FASTA entry name matches contaminant_regex (FragPipe's "contam_" tag); they
# are listed in the "contaminant_proteins" attribute.
compute_ibaq_fragpipe <- function(
  peptides,
  sequences,
  unique_only = FALSE,
  exclude_contaminants = FALSE,
  contaminant_regex = "^contam_",
  protease_regex = "(?<=[KR])(?!P)",
  min_peptide_len = 6,
  max_peptide_len = 30,
  max_missed_cleavages = 0
) {
  contaminants <- character()
  if (isTRUE(exclude_contaminants)) {
    is_contam <- grepl(contaminant_regex, peptides$protein)
    contaminants <- unique(peptides$protein[is_contam])
    peptides <- peptides[!is_contam, , drop = FALSE]
  }
  if (isTRUE(unique_only)) {
    peptides <- peptides[
      is.na(peptides$mapped_proteins) | peptides$mapped_proteins == "",
      ,
      drop = FALSE
    ]
  }

  intensity <- peptides %>%
    dplyr::group_by(sample = experiment, protein) %>%
    dplyr::summarise(
      sum_intensity = sum(intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(sum_intensity > 0)
  if (nrow(intensity) == 0) {
    stop("No peptide with a non-zero intensity; nothing to quantify.")
  }

  quantified <- unique(intensity$protein)
  proteins <- peptides %>%
    dplyr::filter(protein %in% quantified) %>%
    dplyr::distinct(protein, .keep_all = TRUE) %>%
    dplyr::select(protein, protein_id, entry_name, gene)

  # Match on the header's first word, then fall back to the accession for a
  # FASTA whose headers differ from the search database.
  idx <- match(proteins$protein, names(sequences))
  miss <- is.na(idx)
  if (any(miss)) {
    idx[miss] <- match(
      proteins$protein_id[miss],
      fasta_accessions(names(sequences))
    )
  }
  proteins$n_theoretical_peptides <- NA_integer_
  found <- !is.na(idx)
  proteins$n_theoretical_peptides[found] <- unname(count_theoretical_peptides(
    sequences[idx[found]],
    protease_regex = protease_regex,
    min_len = min_peptide_len,
    max_len = max_peptide_len,
    max_missed = max_missed_cleavages
  ))

  out <- intensity %>%
    dplyr::inner_join(proteins[found, , drop = FALSE], by = "protein") %>%
    # Proteins without any observable peptide would divide by zero
    dplyr::filter(
      !is.na(n_theoretical_peptides),
      n_theoretical_peptides > 0
    ) %>%
    dplyr::mutate(iBAQ = sum_intensity / n_theoretical_peptides) %>%
    dplyr::select(
      protein_id,
      entry_name,
      gene,
      n_theoretical_peptides,
      sample,
      iBAQ
    ) %>%
    dplyr::arrange(protein_id, sample)
  attr(out, "unmatched_proteins") <- proteins$protein[!found]
  attr(out, "contaminant_proteins") <- contaminants
  out
}


# ── MS/MS spectrum builder (FragPipe PSM columns) ──────────────────────

# Tidy the combined PSM data for a single peptide sequence.
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


# Build an annotated MS/MS fragmentation spectrum ggplot from FragPipe PSM data.
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


# ── Body UI — organized into 2 columns ─────────────────────────────
PSManalyst_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "sidebar-section-label", "PSManalyst Controls"),
    tags$div(
      style = "padding:0 16px;",
      textInput(
        ns("psm_folder"),
        "Path to folder with psm.tsv files",
        placeholder = "/path/to/fragpipe"
      ),
      actionButton(
        ns("load_psm_folder"),
        "Load PSM Files",
        class = "btn-primary btn-sm",
        style = "width:80%;margin-bottom:6px;"
      ),
      uiOutput(ns("psm_folder_status")),
      selectizeInput(
        ns("sample_select"),
        "Samples (empty = all)",
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "All samples",
          plugins = list("remove_button")
        )
      ),
      tags$hr(style = "border-color:#2d3741;"),
      sliderInput(ns("hyperscore"), "Hyperscore", 0, 1000, 0, 5),
      sliderInput(ns("probability"), "Probability", 0, 1, 0.95, 0.01),
      selectInput(
        ns("specificity_filter"),
        "Proteolysis specificity",
        choices = c(
          "All" = "all",
          "Fully specific" = "fully_specific",
          "Semi N-termini" = "semi_n_termini",
          "Semi C-termini" = "semi_c_termini",
          "Fully semi-specific" = "fully_semi_specific"
        )
      ),
      colourpicker::colourInput(
        ns("plot_color"),
        "Plot color",
        value = "#5499c7"
      ),
      tags$hr(style = "border-color:#2d3741;"),
      fileInput(
        ns("combined_protein"),
        "combined_protein.tsv",
        accept = ".tsv"
      ),
      tags$hr(style = "border-color:#2d3741;"),
      tags$div(
        style = "padding:12px 16px 4px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
        "Abundance Matrix Export"
      ),
      selectInput(
        ns("export_abundance_metric"),
        "Abundance Metric",
        choices = c(
          "Intensity" = "intensity",
          "MaxLFQ Intensity" = "max_lfq_intensity"
        ),
        selected = "max_lfq_intensity"
      ),
      selectInput(
        ns("export_protein_id"),
        "Protein Identifier",
        choices = c(
          "Protein ID" = "protein_id",
          "Entry Name" = "entry_name",
          "Gene" = "gene"
        ),
        selected = "protein_id"
      ),
      downloadButton(
        ns("download_abundance_matrix"),
        tagList(icon("download"), " Export Abundance Matrix"),
        class = "dl-btn",
        style = "width:100%;margin-bottom:8px;"
      ),
      fileInput(
        ns("fasta_file"),
        "FASTA file",
        accept = c(".fasta", ".fa", ".fas")
      ),
      tags$hr(style = "border-color:#2d3741;"),
      downloadButton(
        ns("download_all_plots"),
        tagList(icon("file-zipper"), " Download All Plots (.zip)"),
        class = "dl-btn",
        style = "width:100%;"
      )
    )
  )
}

PSManalyst_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(infoBoxOutput(ns("info_box1"), width = 12)),
    tabsetPanel(
      id = ns("tabs"),
      type = "tabs",

      # ── TAB 1: PSM Viewer (controls on top, full-width plot) ──
      tabPanel(
        "PSM Viewer",
        fluidRow(
          box(
            title = "PSM Plot Controls",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(
                6,
                selectInput(
                  ns("psm_plot_select"),
                  "Select PSM Plot",
                  choices = c(
                    "Protease fingerprint" = "plot01",
                    "N-termini SeqLogo" = "plot03",
                    "C-termini SeqLogo" = "plot04",
                    "m/z over RT" = "plot05",
                    "Mass error (ppm)" = "plot06",
                    "Peptide length" = "plot07",
                    "GRAVY" = "plot08",
                    "Isoelectric Point (pI)" = "plot09",
                    "Charge state" = "plot10",
                    "Missed cleavages" = "plot11",
                    "Uniqueness" = "plot12",
                    "Hyperscore" = "plot13",
                    "Next Score" = "plot14",
                    "PeptideProphet Probability" = "plot15",
                    "Expectation" = "plot16",
                    "Assigned Modifications" = "plot17",
                    "Top 20 proteins" = "plot18",
                    "N:C-terminus matrix" = "plot19",
                    "Cysteine counts" = "plot20",
                    "AA frequency vs FASTA" = "plot_aa_freq",
                    "Peptide Yield vs. FDR" = "plot_fdr_curve"
                  )
                )
              ),
              column(
                3,
                actionButton(
                  ns("run_psm_plot"),
                  "Build Plot",
                  class = "btn-primary",
                  style = "width:100%;margin-top:25px;"
                )
              ),
              column(
                3,
                downloadButton(
                  ns("download_psm_plot"),
                  tagList(icon("download"), " Download (.png)"),
                  class = "dl-btn",
                  style = "width:100%;margin-top:25px;"
                )
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "PSM Plot Output",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp_psm_main"),
                icon("spinner", class = "fa-spin")
              ),
              plotOutput(ns("psm_dynamic_plot_out"), height = "700px")
            )
          )
        )
      ),

      # ── TAB 2: MS/MS Spectrum (controls on top, full-width plot) ──
      tabPanel(
        title = tagList(icon("chart-bar"), "MS/MS Spectrum Viewer"),
        fluidRow(
          box(
            title = "Spectrum Controls",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(
                5,
                selectizeInput(
                  ns("spectrum_peptide"),
                  "Peptide Sequence",
                  choices = NULL,
                  options = list(placeholder = "Load PSM files first…")
                )
              ),
              column(
                3,
                numericInput(
                  ns("spectrum_label_size"),
                  "Ion Label Size",
                  value = 3,
                  min = 1,
                  max = 8,
                  step = 0.5
                )
              ),
              column(
                2,
                actionButton(
                  ns("run_spectrum"),
                  "Plot Spectrum",
                  class = "btn-primary",
                  style = "width:100%;margin-top:25px;"
                )
              ),
              column(
                2,
                downloadButton(
                  ns("download_spectrum_plot"),
                  tagList(icon("download"), " Download (.pdf)"),
                  class = "dl-btn",
                  style = "width:100%;margin-top:25px;"
                )
              )
            )
          )
        ),
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
              uiOutput(ns("spectrum_plot_ui"))
            )
          )
        ),
        fluidRow(
          box(
            title = "Tidy fragment ion data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput(ns("spectrum_tidy_table"))
          )
        )
      ),

      # ── TAB 3: Protein Viewer (controls on top, full-width plot) ──
      tabPanel(
        "Protein Viewer",
        fluidRow(infoBoxOutput(ns("info_box2"), width = 12)),
        fluidRow(
          box(
            title = "Protein Plot Controls",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(
                4,
                selectInput(
                  ns("prot_plot_select"),
                  "Select Protein Plot",
                  choices = c(
                    "Coverage" = "plot01p",
                    "Organisms" = "plot02p",
                    "Protein existence" = "plot03p",
                    "Protein probability" = "plot04p",
                    "Top peptide probability" = "plot05p",
                    "Total peptides" = "plot06p",
                    "Razor spectral count" = "plot07p",
                    "Razor intensity" = "plot08p",
                    "MaxLFQ distribution" = "plot10p",
                    "Top 20 by MaxLFQ (log₂)" = "plot_rank_lfq",
                    "Top 20 by spectral count" = "plot_rank_usc"
                  )
                )
              ),
              column(
                4,
                actionButton(
                  ns("run_prot_plot"),
                  "Build Plot",
                  class = "btn-primary",
                  style = "width:100%;margin-top:25px;"
                )
              ),
              column(
                4,
                downloadButton(
                  ns("download_prot_plot"),
                  tagList(icon("download"), " Download (.png)"),
                  class = "dl-btn",
                  style = "width:100%;margin-top:25px;"
                )
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Protein Plot Output",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            div(
              class = "plot-wrap",
              tags$div(
                class = "spinner-overlay",
                id = ns("sp_prot_main"),
                icon("spinner", class = "fa-spin")
              ),
              plotOutput(ns("prot_dynamic_plot_out"), height = "700px")
            )
          )
        ),
        fluidRow(
          box(
            title = "Sequence Coverage View",
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
                numericInput(
                  ns("aa_per_line"),
                  "AAs per line",
                  value = 50,
                  min = 10,
                  max = 200,
                  step = 10
                )
              ),
              column(
                3,
                checkboxInput(
                  ns("show_sequence"),
                  "Show residues",
                  value = TRUE
                )
              ),
              column(
                3,
                actionButton(
                  ns("plot_seq"),
                  "Update Coverage",
                  style = "margin-top:25px;width:100%;"
                )
              )
            ),
            uiOutput(ns("protein_coverage_plot_ui"))
          )
        )
      ),

      # ── TAB 4: Modification Diagnostic (controls on top, full-width plot) ──
      tabPanel(
        "Modification Diagnostic",
        fluidRow(
          box(
            title = "Diagnostic Controls",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(
                6,
                sliderInput(
                  ns("rt_tolerance"),
                  "ΔRT Artifact Threshold (min)",
                  min = 0.1,
                  max = 10,
                  value = 0.5,
                  step = 0.1
                )
              ),
              column(
                6,
                helpText(
                  style = "margin-top:25px;",
                  "Lower threshold → more strict artifact classification."
                )
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "RT Shift Profile — Modified vs Unmodified",
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
            title = "Paired Peptide Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput(ns("mod_diag_table"))
          )
        )
      ),

      # ── TAB 5: Similarity & Distance ──
      tabPanel(
        "Similarity & Distance",
        fluidRow(
          box(
            title = "Scatterplot matrix",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            uiOutput(ns("plot_ggpairs_ui"))
          ),
          box(
            title = "Cosine Similarity",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            uiOutput(ns("cosine_similarity_ui"))
          )
        ),
        fluidRow(
          box(
            title = "Euclidean Distance",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            uiOutput(ns("euclidean_distance_ui"))
          ),
          box(
            title = "Jaccard Similarity (Protein ID)",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            uiOutput(ns("jaccard_similarity_ui"))
          )
        )
      ),

      # ── TAB 6: iBAQ (from peptide.tsv) ──
      tabPanel(
        "iBAQ",
        fluidRow(uiOutput(ns("ibaq_banner"))),
        fluidRow(box(
          title = "iBAQ Controls",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          fluidRow(
            column(
              2,
              numericInput(
                ns("ibaq_min_len"),
                "Min peptide length",
                value = 6,
                min = 1,
                max = 100,
                step = 1
              )
            ),
            column(
              2,
              numericInput(
                ns("ibaq_max_len"),
                "Max peptide length",
                value = 30,
                min = 1,
                max = 100,
                step = 1
              )
            ),
            column(
              2,
              numericInput(
                ns("ibaq_missed_cleavages"),
                "Missed cleavages",
                value = 0,
                min = 0,
                max = 5,
                step = 1
              )
            ),
            column(
              2,
              textInput(ns("ibaq_cut"), "Cleaved residues", value = "KR")
            ),
            column(
              2,
              textInput(ns("ibaq_nocut"), "Blocked by residue", value = "P")
            ),
            column(
              2,
              selectInput(
                ns("ibaq_sense"),
                "Cleavage side",
                choices = c("C-terminal" = "C", "N-terminal" = "N")
              )
            )
          ),
          uiOutput(ns("ibaq_digest_source")),
          fluidRow(
            column(
              3,
              numericInput(
                ns("ibaq_top_n"),
                "Top N proteins (bar chart)",
                value = 20,
                min = 1,
                max = 100,
                step = 1
              )
            ),
            column(
              3,
              numericInput(
                ns("ibaq_rank_labels"),
                "Labels per facet (ranked abundance)",
                value = 5,
                min = 0,
                max = 50,
                step = 1
              )
            ),
            column(
              2,
              style = "margin-top:25px;",
              checkboxInput(
                ns("ibaq_unique"),
                "Unique peptides only",
                value = FALSE
              )
            ),
            column(
              2,
              style = "margin-top:25px;",
              checkboxInput(
                ns("ibaq_no_contam"),
                "Exclude contaminants",
                value = FALSE
              )
            ),
            column(
              2,
              style = "margin-top:25px;",
              checkboxInput(
                ns("ibaq_log2"),
                "log₂(iBAQ + 1) in table / export",
                value = TRUE
              )
            )
          ),
          actionButton(
            ns("run_ibaq"),
            "Compute iBAQ",
            icon = icon("calculator"),
            class = "btn-primary",
            style = "margin-right:8px;margin-bottom:8px;"
          ),
          downloadButton(
            ns("download_ibaq_table"),
            tagList(icon("download"), " iBAQ Matrix (.tsv)"),
            class = "dl-btn",
            style = "margin-right:8px;margin-bottom:8px;"
          ),
          downloadButton(
            ns("download_all_ibaq_plots"),
            tagList(icon("file-zipper"), " Download All Plots (.zip)"),
            class = "dl-btn",
            style = "margin-bottom:8px;"
          ),
          uiOutput(ns("ibaq_summary"))
        )),
        fluidRow(box(
          title = "iBAQ Distribution per Sample",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          plotOutput(ns("ibaq_dist_plot"), height = "500px")
        )),
        fluidRow(box(
          title = "Dynamic Range — Ranked Protein Abundance",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          uiOutput(ns("ibaq_rank_plot_ui"))
        )),
        fluidRow(box(
          title = "Most Abundant Proteins — Relative iBAQ",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          uiOutput(ns("ibaq_top_plot_ui"))
        )),
        fluidRow(box(
          title = "iBAQ Table",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          DT::dataTableOutput(ns("ibaq_table"))
        ))
      )
    )
  )
}

# ── Server ─────────────────────────────
PSManalyst_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    psm_fns <- new.env()
    prot_fns <- new.env()

    # ── Spinner helpers ──
    hide_spinner <- function(sid) shinyjs::hide(id = sid)
    show_spinner <- function(sid) shinyjs::show(id = sid)

    # ── Dynamic height logic ──
    psm_plot_h <- reactive({
      d <- data()
      req(d, "sample_name" %in% names(d))
      n <- length(unique(d$sample_name))
      max(400L, ceiling(n / 3) * 350L)
    })

    make_psm_plot_ui <- function(plot_id) {
      renderUI({
        req(data())
        plotOutput(ns(plot_id), height = paste0(psm_plot_h(), "px"))
      })
    }

    output$plot_ggpairs_ui <- make_psm_plot_ui("plot_ggpairs")
    output$cosine_similarity_ui <- make_psm_plot_ui("cosine_similarity")
    output$euclidean_distance_ui <- make_psm_plot_ui("euclidean_distance")
    output$jaccard_similarity_ui <- make_psm_plot_ui("jaccard_similarity")

    # ════════════════════════════════════════════════════════════════════════
    # A. MULTI-FILE PSM INGESTION
    # ════════════════════════════════════════════════════════════════════════
    raw_psm_data <- reactiveVal(NULL)
    loaded_psm_folder <- reactiveVal(NULL)

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

    output$psm_folder_status <- renderUI({
      df <- raw_psm_data()
      if (is.null(df)) {
        return(NULL)
      }
      tags$div(
        style = "padding:4px 0 8px;font-size:12px;color:#1BB99A;font-weight:600;",
        icon("circle-check", style = "margin-right:4px;"),
        sprintf(
          "%d sample(s) · %s PSMs",
          dplyr::n_distinct(df$sample_name),
          format(nrow(df), big.mark = ",")
        )
      )
    })

    # ════════════════════════════════════════════════════════════════════════
    # B. FILTERED PSM DATA
    # ════════════════════════════════════════════════════════════════════════
    # ── Global sample selector (sidebar) ──
    # Populated from the raw PSM data; an empty selection means "all samples".
    observeEvent(raw_psm_data(), {
      req(raw_psm_data(), "sample_name" %in% names(raw_psm_data()))
      samples <- sort(unique(raw_psm_data()$sample_name))
      updateSelectizeInput(
        session,
        "sample_select",
        choices = samples,
        selected = character(0)
      )
    })

    selected_samples <- reactive({
      sel <- input$sample_select
      if (is.null(sel) || length(sel) == 0) NULL else sel
    })

    # Restrict any data frame with a sample_name column to the selected samples.
    filter_samples <- function(d) {
      sel <- selected_samples()
      if (is.null(sel) || !"sample_name" %in% names(d)) {
        return(d)
      }
      dplyr::filter(d, sample_name %in% sel)
    }

    # The per-folder tables (protein.tsv, peptide.tsv) name their samples after
    # the experiment folder, while the sidebar selector lists PSM sample names
    # (spectrum-file prefixes): translate the selection. NULL = all.
    selected_experiments <- reactive({
      sel <- selected_samples()
      if (is.null(sel)) {
        return(NULL)
      }
      psm <- raw_psm_data()
      linked <- if ("experiment" %in% names(psm)) {
        psm$experiment[psm$sample_name %in% sel]
      }
      unique(as.character(c(sel, linked)))
    })

    # Same as filter_samples(), for tables whose sample_name is the folder.
    filter_experiments <- function(d) {
      sel <- selected_experiments()
      if (is.null(sel) || !"sample_name" %in% names(d)) {
        return(d)
      }
      dplyr::filter(d, as.character(sample_name) %in% sel)
    }

    data <- reactive({
      req(raw_psm_data())

      psm_file <- raw_psm_data() %>%
        filter_samples() %>%
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
          isoelectric_point = calculate_pI(peptide),
          molecular_weight = calculate_MW(peptide),
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

    tidy_spectrum_data <- eventReactive(input$run_spectrum, {
      req(raw_psm_data(), nchar(input$spectrum_peptide) > 0)
      show_spinner("sp_spectrum")
      withProgress(message = "Building spectrum…", value = 0.5, {
        result <- tidy_psm_spectrum(
          filter_samples(raw_psm_data()),
          input$spectrum_peptide
        )
        incProgress(0.5, detail = "Done.")
        result
      })
    })

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
        "protein.tsv contain FDR-filtered protein results; each row is an identified protein group",
        icon = icon("info"),
        color = "black"
      )
    })

    # ════════════════════════════════════════════════════════════════════════
    # E. PSM PLOT FUNCTIONS
    # ════════════════════════════════════════════════════════════════════════
    seqlogo_scale <- scale_x_continuous(
      breaks = 1:8,
      labels = c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
    )

    psm_fns$plot01 <- function() {
      req(data())
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
          axis.text = element_text(face = "bold", color = "black"),
          axis.title = element_text(size = 13, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold", color = "black"),
          legend.position = "bottom",
          legend.key.width = unit(2, "cm"),
          legend.key.height = unit(0.25, "cm"),
          legend.title.position = "top",
          panel.grid = element_blank()
        )
    }

    psm_fns$plot03 <- function() {
      req(data())
      seq_lst <- data() %>%
        dplyr::filter(!is.na(fingerprint_Nterm)) %>%
        dplyr::group_by(sample_name) %>%
        dplyr::summarise(seqs = list(fingerprint_Nterm), .groups = "drop") %>%
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
          axis.text = element_text(face = "bold", color = "black"),
          axis.title = element_text(size = 13, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold", color = "black")
        )
    }

    psm_fns$plot04 <- function() {
      req(data())
      seq_lst <- data() %>%
        dplyr::filter(!is.na(fingerprint_Cterm)) %>%
        dplyr::group_by(sample_name) %>%
        dplyr::summarise(seqs = list(fingerprint_Cterm), .groups = "drop") %>%
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
          axis.text = element_text(face = "bold", color = "black"),
          axis.title = element_text(size = 13, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold", color = "black")
        )
    }

    psm_fns$plot05 <- function() {
      req(data())
      data() %>%
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
    }

    psm_fns$plot06 <- function() {
      req(data())
      data() %>%
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
    }

    psm_fns$plot07 <- function() {
      req(data())
      data() %>%
        ggplot() +
        geom_density(aes(x = peptide_length), fill = input$plot_color) +
        labs(x = "Peptide Length", y = "Frequency (%)") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot08 <- function() {
      req(data())
      data() %>%
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
    }

    psm_fns$plot09 <- function() {
      req(data())
      plot_2d_gel(
        data(),
        pi_col = "isoelectric_point",
        mw_col = "molecular_weight",
        title = "Virtual 2D Gel",
        facet_col = "sample_name"
      ) +
        theme(
          legend.position = "bottom",
          legend.key.width = unit(2.5, "cm"),
          legend.key.height = unit(0.25, "cm")
        )
    }

    psm_fns$plot10 <- function() {
      req(data())
      data() %>%
        ggplot() +
        geom_bar(aes(x = charge), fill = input$plot_color, color = "black") +
        labs(x = "Charge state", y = "Count") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot11 <- function() {
      req(data())
      data() %>%
        group_by(sample_name, number_of_missed_cleavages) %>%
        summarise(n = n(), .groups = "drop") %>%
        ggplot(aes(x = number_of_missed_cleavages, y = n)) +
        geom_bar(stat = "identity", fill = input$plot_color, color = "black") +
        geom_text(aes(label = n), vjust = -0.5, size = 5) +
        labs(x = "Number of Missed Cleavages", y = "Count") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot12 <- function() {
      req(data())
      data() %>%
        dplyr::mutate(
          uniqueness = ifelse(is_unique == TRUE, "Unique", "Shared")
        ) %>%
        group_by(sample_name, uniqueness) %>%
        summarise(n = n(), .groups = "drop") %>%
        ggplot(aes(x = uniqueness, y = n)) +
        geom_bar(stat = "identity", fill = input$plot_color, color = "black") +
        geom_text(aes(label = n), vjust = -0.5, size = 5) +
        labs(x = "Unique peptides", y = "Count") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot13 <- function() {
      req(data())
      data() %>%
        ggplot() +
        geom_histogram(
          aes(x = hyperscore),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "Hyperscore", y = "Count") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot14 <- function() {
      req(data())
      data() %>%
        ggplot() +
        geom_histogram(
          aes(x = nextscore),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "Nextscore", y = "Count") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot15 <- function() {
      req(data())
      data() %>%
        ggplot() +
        geom_histogram(
          aes(x = probability),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "PeptideProphet Probability", y = "Count") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot16 <- function() {
      req(data())
      data() %>%
        ggplot() +
        geom_histogram(
          aes(x = expectation),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "Expectation value", y = "Count") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot17 <- function() {
      req(data())
      data() %>%
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
        summarise(n = n(), .groups = "drop") %>%
        ggplot(aes(y = assigned_modifications, x = n)) +
        geom_col(fill = input$plot_color, color = "black") +
        geom_text(aes(label = n), hjust = -0.1, size = 5) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(y = "Assigned Modifications", x = "Count") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot18 <- function() {
      req(data())
      data() %>%
        dplyr::group_by(sample_name, entry_name) %>%
        dplyr::summarise(n_psm = n(), .groups = "drop") %>%
        dplyr::group_by(sample_name) %>%
        dplyr::slice_max(n_psm, n = 20, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          entry_name = tidytext::reorder_within(entry_name, n_psm, sample_name)
        ) %>%
        ggplot(aes(x = n_psm, y = entry_name)) +
        geom_col(fill = input$plot_color, color = "black") +
        tidytext::scale_y_reordered() +
        labs(x = "Number of PSMs", y = "Protein") +
        facet_wrap(~sample_name, ncol = 3, scales = "free_y")
    }

    psm_fns$plot19 <- function() {
      req(data())
      cooccurrence_data()[["plot"]]
    }

    psm_fns$plot20 <- function() {
      req(data())
      data() %>%
        dplyr::mutate(cysteine_count = stringr::str_count(peptide, "C")) %>%
        dplyr::group_by(sample_name, cysteine_count) %>%
        summarise(n = n(), .groups = "drop") %>%
        ggplot(aes(x = cysteine_count, y = n)) +
        geom_bar(stat = "identity", fill = input$plot_color, color = "black") +
        geom_text(aes(label = n), vjust = -0.5, size = 5) +
        labs(x = "Cysteine counts in peptides", y = "Count") +
        facet_wrap(~sample_name, ncol = 3)
    }

    psm_fns$plot_aa_freq <- function() {
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
      seqs <- sapply(fasta_data(), function(x) as.character(x)[1])
      all_fasta_aa <- toupper(gsub(
        "[^A-Za-z]",
        "",
        paste0(seqs, collapse = "")
      ))
      aa_counts_fasta <- table(strsplit(all_fasta_aa, "")[[1]])
      aa_counts_fasta <- aa_counts_fasta[names(aa_counts_fasta) %in% base_aas]
      exp_freq <- as.data.frame(
        aa_counts_fasta / sum(aa_counts_fasta, na.rm = TRUE) * 100
      )
      colnames(exp_freq) <- c("AA", "Frequency")
      exp_freq$Type <- "Expected (FASTA)"

      d <- data() %>%
        dplyr::filter(!is.na(peptide), nchar(peptide) > 0) %>%
        dplyr::mutate(clean_peptide = gsub("\\[.*?\\]", "", peptide))

      obs_freq <- d %>%
        dplyr::group_by(sample_name) %>%
        dplyr::summarise(
          all_aas = list(strsplit(paste0(clean_peptide, collapse = ""), "")[[
            1
          ]]),
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
          axis.text = element_text(size = 10, face = "bold", color = "black"),
          axis.title = element_text(size = 11, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", color = "black")
        )
    }

    # ════════════════════════════════════════════════════════════════════════
    # E3. FDR CURVE
    # ════════════════════════════════════════════════════════════════════════
    fdr_curve_obj <- reactive({
      d <- data()
      req(d, nrow(d) > 0, "probability" %in% names(d))

      q_sorted <- d %>%
        dplyr::group_by(sample_name) %>%
        dplyr::arrange(probability) %>%
        dplyr::mutate(
          qvalue = 1 - probability,
          cumulative_peptides = dplyr::row_number()
        ) %>%
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

    # ════════════════════════════════════════════════════════════════════════
    # E4. MODIFICATION DIAGNOSTIC
    # ════════════════════════════════════════════════════════════════════════
    # Sample restriction is inherited from data() via the sidebar selector.
    mod_diag_data <- reactive({
      d <- data()
      req(d, nrow(d) > 0)
      req("assigned_modifications" %in% names(d), "retention" %in% names(d))

      d_mod <- d %>%
        dplyr::filter(
          !is.na(assigned_modifications),
          assigned_modifications != ""
        ) %>%
        dplyr::select(
          peptide,
          retention,
          assigned_modifications,
          sample_name,
          gene
        ) %>%
        dplyr::mutate(retention = retention / 60) %>%
        dplyr::rename(rt_mod = retention)

      d_unmod <- d %>%
        dplyr::filter(
          is.na(assigned_modifications) | assigned_modifications == ""
        ) %>%
        dplyr::group_by(peptide, sample_name) %>%
        dplyr::summarise(
          rt_unmod = median(retention / 60, na.rm = TRUE),
          .groups = "drop"
        )

      dplyr::inner_join(d_mod, d_unmod, by = c("peptide", "sample_name")) %>%
        dplyr::mutate(
          delta_rt = abs(rt_mod - rt_unmod),
          classification = dplyr::if_else(
            delta_rt <= input$rt_tolerance,
            "Suspected Artifact",
            "Sample-Derived"
          )
        )
    })

    mod_diag_plot_obj <- reactive({
      paired <- mod_diag_data()
      req(paired, nrow(paired) > 0)

      top_peptides <- paired %>%
        dplyr::group_by(sample_name, peptide) %>%
        dplyr::summarise(
          max_delta = max(delta_rt, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::group_by(sample_name) %>%
        dplyr::slice_max(order_by = max_delta, n = 30) %>%
        dplyr::pull(peptide) %>%
        unique()

      plot_data <- paired %>% dplyr::filter(peptide %in% top_peptides)
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
            "ΔRT threshold = ",
            input$rt_tolerance,
            " min | Grey = unmodified, Colored = modified"
          )
        ) +
        facet_wrap(~sample_name, ncol = 2, scales = "free") +
        theme_bw() +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 7, face = "bold", color = "black"),
          axis.text.x = element_text(face = "bold", color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", face = "bold"),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "bottom",
          legend.text = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          plot.caption = element_text(size = 12, face = "bold")
        )
    })

    output$mod_diag_plot_ui <- renderUI({
      paired <- tryCatch(mod_diag_data(), error = function(e) NULL)
      hide_spinner("sp_moddiag")
      if (is.null(paired) || nrow(paired) == 0) {
        return(tags$p(
          style = "color:#adb5bd;text-align:center;padding:20px;",
          "No paired modified/unmodified peptides found. Load PSM data first and ensure it contains 'assigned_modifications' and 'retention' columns."
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
      tbl <- paired %>%
        dplyr::select(
          sample_name,
          peptide,
          assigned_modifications,
          gene,
          rt_unmod,
          rt_mod,
          delta_rt,
          classification
        ) %>%
        dplyr::arrange(dplyr::desc(delta_rt))
      DT::datatable(
        tbl,
        rownames = FALSE,
        filter = "top",
        options = list(pageLength = 20, scrollX = TRUE)
      ) %>%
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
      prot_files <- list.files(
        folder,
        pattern = "^protein\\.tsv$",
        full.names = TRUE,
        recursive = TRUE
      )
      req(length(prot_files) > 0)
      # Read every protein.tsv and tag each with a sample_name derived from its
      # parent folder so that plots can be faceted / filtered by sample.
      purrr::map_dfr(prot_files, function(f) {
        data.table::fread(
          f,
          sep = "\t",
          header = TRUE,
          colClasses = "character"
        ) %>%
          janitor::clean_names() %>%
          dplyr::mutate(sample_name = basename(dirname(f)))
      }) %>%
        dplyr::mutate(dplyr::across(
          dplyr::everything(),
          ~ utils::type.convert(.x, as.is = TRUE)
        )) %>%
        filter_experiments()
    })

    combined_protein_data <- reactive({
      req(input$combined_protein)
      data.table::fread(
        input$combined_protein$datapath,
        sep = "\t",
        header = TRUE
      ) %>%
        janitor::clean_names() %>%
        dplyr::select(protein_id, ends_with("max_lfq_intensity")) %>%
        column_to_rownames("protein_id") %>%
        dplyr::rename_all(~ stringr::str_remove(., "_max_lfq_intensity")) %>%
        # MaxLFQ of 0 means "not quantified": treat as NA so log2 does not
        # yield -Inf and downstream metrics only use quantified proteins.
        dplyr::mutate(dplyr::across(everything(), ~ dplyr::na_if(., 0))) %>%
        log2() %>%
        subset_combined_samples()
    })

    # Keep only columns matching the sidebar sample selection. Column names in
    # combined_protein.tsv may not match PSM-derived sample names exactly, so
    # the subset is only applied when at least one selected sample matches.
    subset_combined_samples <- function(m) {
      sel <- selected_samples()
      if (is.null(sel)) {
        return(m)
      }
      keep <- intersect(colnames(m), c(sel, janitor::make_clean_names(sel)))
      if (length(keep) == 0) {
        return(m)
      }
      m[, keep, drop = FALSE]
    }

    combined_protein_raw <- reactive({
      req(input$combined_protein)
      raw <- data.table::fread(
        input$combined_protein$datapath,
        sep = "\t",
        header = TRUE
      ) %>%
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

    # Protein selector from protein_data
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

    # ── Sequence coverage ──
    # Coverage helpers (fasta_index_exact, fasta_seq_normalize,
    # protein_coverage) live in modules/utils_fasta.r.
    coverage_plot_data <- eventReactive(
      input$plot_seq,
      {
        req(input$selected_protein)
        if (input$selected_protein == "") {
          return(list(protein_sequence = NULL))
        }
        if (is.null(input$fasta_file)) {
          return(list(error = "Please upload a FASTA file."))
        }
        if (is.null(raw_psm_data())) {
          return(list(error = "Please load PSM files first."))
        }

        req(data(), fasta_data())
        tp <- input$selected_protein
        tix <- fasta_index_exact(tp, names(fasta_data()))
        if (is.na(tix)) {
          return(list(error = paste("Protein", tp, "not found in FASTA.")))
        }

        pr_seq <- fasta_seq_normalize(fasta_data()[[tix]])
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
            sequence = stringr::str_remove_all(
              stringr::str_remove_all(sequence, "\\[.*?\\]"),
              "[^A-Z]"
            )
          )
        list(protein_sequence = pr_seq, peptides_data = agg)
      },
      ignoreNULL = FALSE
    )

    output$protein_coverage_plot_ui <- renderUI({
      pd <- coverage_plot_data()
      if (is.null(pd) || (is.null(pd$error) && is.null(pd$protein_sequence))) {
        return(div(
          style = "color:#6c757d;padding:15px;text-align:center;",
          "Select a protein and click 'Update Coverage'."
        ))
      }
      if (!is.null(pd$error)) {
        return(div(
          style = "color:red;padding:15px;font-weight:bold;",
          pd$error
        ))
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
      cvg <- protein_coverage(
        pd$protein_sequence,
        pd$peptides_data$sequence,
        pd$peptides_data$psm_count
      )
      cov <- cvg$depth
      mask <- cvg$mask
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
              linewidth = 0.2,
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
              linewidth = 0.2
            ) +
            scale_fill_gradient(
              low = "#e1f5fe",
              high = input$plot_color,
              name = "Coverage (Norm)",
              na.value = "lightgray"
            )
        }
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

    # ════════════════════════════════════════════════════════════════════════
    # F2. PROTEIN PLOTS
    # ════════════════════════════════════════════════════════════════════════

    # Protein data already restricted to the sidebar sample selection.
    prot_plot_data <- function() {
      d <- protein_data()
      req(d, "sample_name" %in% names(d))
      d
    }

    # Facet a protein plot by sample only when more than one sample is shown.
    add_prot_facet <- function(p) {
      if (dplyr::n_distinct(prot_plot_data()$sample_name) > 1) {
        p <- p + facet_wrap(~sample_name, ncol = 3)
      }
      p
    }

    prot_fns$plot01p <- function() {
      p <- prot_plot_data() %>%
        ggplot() +
        geom_histogram(
          aes(x = coverage),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "Protein coverage (%)", y = "Count")
      add_prot_facet(p)
    }

    prot_fns$plot02p <- function() {
      p <- prot_plot_data() %>%
        dplyr::count(sample_name, organism) %>%
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
      add_prot_facet(p)
    }

    prot_fns$plot03p <- function() {
      p <- prot_plot_data() %>%
        dplyr::mutate(
          protein_existence = stringr::str_remove(protein_existence, ".*\\:"),
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
      add_prot_facet(p)
    }

    prot_fns$plot04p <- function() {
      p <- prot_plot_data() %>%
        ggplot() +
        geom_histogram(
          aes(x = protein_probability),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "Protein Probability", y = "Count")
      add_prot_facet(p)
    }

    prot_fns$plot05p <- function() {
      p <- prot_plot_data() %>%
        ggplot() +
        geom_histogram(
          aes(x = top_peptide_probability),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "Peptide Probability", y = "Count")
      add_prot_facet(p)
    }

    prot_fns$plot06p <- function() {
      p <- prot_plot_data() %>%
        ggplot() +
        geom_histogram(
          aes(x = total_peptides),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "Total peptides mapped to proteins", y = "Count")
      add_prot_facet(p)
    }

    prot_fns$plot07p <- function() {
      p <- prot_plot_data() %>%
        ggplot() +
        geom_histogram(
          aes(x = razor_spectral_count),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "Razor Spectral Count", y = "Count")
      add_prot_facet(p)
    }

    prot_fns$plot08p <- function() {
      p <- prot_plot_data() %>%
        ggplot() +
        geom_histogram(
          aes(x = razor_intensity),
          fill = input$plot_color,
          color = "black"
        ) +
        labs(x = "Razor Intensity", y = "Count")
      add_prot_facet(p)
    }

    prot_fns$plot10p <- function() {
      req(combined_protein_data())
      combined_protein_data() %>%
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
        labs(x = NULL, y = "log₂(MaxLFQ intensity)") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }

    prot_fns$plot_rank_lfq <- function() {
      req(combined_protein_raw())
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
          sample_name = stringr::str_remove(sample_name, "_max_lfq_intensity$"),
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
          geom_col(fill = input$plot_color, color = "black", linewidth = 0.3) +
          tidytext::scale_y_reordered() +
          labs(title = s, x = "log₂(MaxLFQ intensity)", y = NULL) +
          theme_bw() +
          theme(
            axis.title = element_text(size = 10, face = "bold"),
            axis.text = element_text(size = 8, face = "bold", color = "black"),
            plot.title = element_text(face = "bold", hjust = 0.5)
          )
      })
      patchwork::wrap_plots(plots, ncol = min(3L, length(plots))) +
        patchwork::plot_annotation(
          title = "Top 20 proteins by MaxLFQ intensity (log₂)",
          theme = theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
        )
    }

    prot_fns$plot_rank_usc <- function() {
      req(combined_protein_raw())
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
        dplyr::filter(!is.na(unique_spectral_count) & unique_spectral_count > 0)

      plots <- lapply(sort(unique(df$sample_name)), function(s) {
        d <- df %>%
          dplyr::filter(sample_name == s) %>%
          dplyr::slice_max(unique_spectral_count, n = 20, with_ties = FALSE) %>%
          dplyr::mutate(
            entry_name = tidytext::reorder_within(
              entry_name,
              unique_spectral_count,
              sample_name
            )
          )
        ggplot(d, aes(x = unique_spectral_count, y = entry_name)) +
          geom_col(fill = input$plot_color, color = "black", linewidth = 0.3) +
          tidytext::scale_y_reordered() +
          labs(title = s, x = "Unique spectral count", y = NULL) +
          theme_bw() +
          theme(
            axis.text = element_text(size = 8, face = "bold", color = "black"),
            axis.title = element_text(size = 10, face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0.5)
          )
      })
      patchwork::wrap_plots(plots, ncol = min(3L, length(plots))) +
        patchwork::plot_annotation(
          title = "Top 20 proteins by unique spectral count",
          theme = theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
        )
    }

    # ── Similarity / distance plots (require combined_protein upload) ──
    # Upper-triangle panel: Pearson r computed only on proteins quantified
    # (non-NA) in BOTH samples of the pair.
    ggpairs_cor_panel <- function(data, mapping, ...) {
      x <- GGally::eval_data_col(data, mapping$x)
      y <- GGally::eval_data_col(data, mapping$y)
      ok <- is.finite(x) & is.finite(y)
      n <- sum(ok)
      r <- if (n > 2) stats::cor(x[ok], y[ok], method = "pearson") else NA_real_
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = sprintf("r = %.3f\nn = %d", r, n),
          size = 5,
          fontface = "bold"
        ) +
        theme_void()
    }

    prot_fns$plot_ggpairs <- function() {
      req(combined_protein_data())
      combined_protein_data() %>%
        GGally::ggpairs(
          lower = list(continuous = wrap("points", alpha = 0.4, na.rm = TRUE)),
          diag = list(continuous = wrap("barDiag", na.rm = TRUE)),
          upper = list(continuous = ggpairs_cor_panel)
        ) +
        theme_bw() +
        theme(
          strip.background = element_blank(),
          strip.text = element_text(face = "bold")
        )
    }

    ht <- theme(
      text = element_text(size = 15),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 1),
      legend.position = "bottom",
      legend.key.width = unit(2.5, "cm"),
      legend.key.height = unit(0.25, "cm")
    )

    # Apply a pairwise metric f(x, y) to every pair of sample columns, using
    # only proteins that are quantified (non-NA) in both samples of the pair.
    pairwise_sample_matrix <- function(m, f) {
      m <- as.matrix(m)
      n <- ncol(m)
      out <- matrix(
        NA_real_,
        n,
        n,
        dimnames = list(colnames(m), colnames(m))
      )
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          ok <- !is.na(m[, i]) & !is.na(m[, j])
          out[i, j] <- if (any(ok)) f(m[ok, i], m[ok, j]) else NA_real_
        }
      }
      out
    }

    matrix_to_long <- function(mat) {
      as.data.frame(mat) %>%
        rownames_to_column(var = "Sample") %>%
        pivot_longer(-Sample, names_to = "Match", values_to = "value")
    }

    similarity_heatmap <- function(long_df, fill_lab) {
      ggplot(long_df, aes(x = Sample, y = Match, fill = value)) +
        geom_tile() +
        viridis::scale_fill_viridis(option = "E") +
        ht +
        labs(x = NULL, y = NULL, fill = fill_lab)
    }

    prot_fns$cosine_similarity <- function() {
      req(combined_protein_data())
      combined_protein_data() %>%
        pairwise_sample_matrix(function(x, y) {
          sum(x * y) / sqrt(sum(x^2) * sum(y^2))
        }) %>%
        matrix_to_long() %>%
        similarity_heatmap("Cosine similarity")
    }

    prot_fns$euclidean_distance <- function() {
      req(combined_protein_data())
      combined_protein_data() %>%
        pairwise_sample_matrix(function(x, y) sqrt(sum((x - y)^2))) %>%
        matrix_to_long() %>%
        similarity_heatmap("Euclidean distance")
    }

    # Jaccard on protein identification (presence/absence):
    # |A ∩ B| / |A ∪ B| where A, B are the sets of quantified proteins.
    prot_fns$jaccard_similarity <- function() {
      req(combined_protein_data())
      pres <- !is.na(as.matrix(combined_protein_data()))
      n <- ncol(pres)
      jac <- matrix(
        NA_real_,
        n,
        n,
        dimnames = list(colnames(pres), colnames(pres))
      )
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          inter <- sum(pres[, i] & pres[, j])
          uni <- sum(pres[, i] | pres[, j])
          jac[i, j] <- if (uni > 0) inter / uni else NA_real_
        }
      }
      jac %>%
        matrix_to_long() %>%
        similarity_heatmap("Jaccard similarity")
    }

    # Render similarity plots (always visible in their tab)
    output$plot_ggpairs <- renderPlot({
      prot_fns$plot_ggpairs()
    })
    output$cosine_similarity <- renderPlot({
      prot_fns$cosine_similarity()
    })
    output$euclidean_distance <- renderPlot({
      prot_fns$euclidean_distance()
    })
    output$jaccard_similarity <- renderPlot({
      prot_fns$jaccard_similarity()
    })

    # ════════════════════════════════════════════════════════════════════════
    # F3. iBAQ  (from the peptide.tsv next to each psm.tsv)
    # ════════════════════════════════════════════════════════════════════════
    peptide_data <- reactive({
      folder <- loaded_psm_folder()
      req(folder)
      read_all_peptide_files(folder)
    })

    # Digest settings of the search become the defaults of the iBAQ controls
    # whenever a folder is loaded; the user can still change them.
    fragger_digest <- reactive({
      folder <- loaded_psm_folder()
      req(folder)
      tryCatch(read_fragger_digest_params(folder), error = function(e) NULL)
    })

    observeEvent(fragger_digest(), {
      p <- fragger_digest()
      if (!is.na(p$min_len)) {
        updateNumericInput(session, "ibaq_min_len", value = p$min_len)
      }
      if (!is.na(p$max_len)) {
        updateNumericInput(session, "ibaq_max_len", value = p$max_len)
      }
      if (!is.na(p$missed_cleavages)) {
        updateNumericInput(
          session,
          "ibaq_missed_cleavages",
          value = p$missed_cleavages
        )
      }
      # Leave the enzyme controls alone when the rule has no cut residue
      # (nonspecific search): there is no digest to mirror.
      if (!is.null(protease_regex_from_rule(p$cut))) {
        updateTextInput(session, "ibaq_cut", value = p$cut)
        updateTextInput(
          session,
          "ibaq_nocut",
          value = if (is.na(p$nocut)) "" else p$nocut
        )
        if (p$sense %in% c("C", "N")) {
          updateSelectInput(session, "ibaq_sense", selected = p$sense)
        }
      }
    })

    output$ibaq_digest_source <- renderUI({
      req(loaded_psm_folder())
      p <- fragger_digest()
      msg <- if (is.null(p)) {
        "No fragger.params found in the loaded folder: digest settings are the app defaults."
      } else {
        paste0(
          "Digest defaults read from fragger.params",
          if (!is.na(p$enzyme)) paste0(" (enzyme: ", p$enzyme, ")"),
          ". Edit the fields above to override them.",
          if (is.null(protease_regex_from_rule(p$cut))) {
            " The search enzyme has no cleavage rule (nonspecific), so the enzyme fields were left unchanged."
          },
          if (!is.na(p$enzyme2)) {
            paste0(
              " A second enzyme (",
              p$enzyme2,
              ") was used in the search and is not included in the digest."
            )
          }
        )
      }
      helpText(msg)
    })

    # Sequences for the in-silico digest: the sidebar FASTA when uploaded,
    # otherwise the protein.fas files FragPipe writes into each experiment
    # folder (sequences of the identified proteins).
    ibaq_sequences <- reactive({
      if (!is.null(input$fasta_file)) {
        fa <- fasta_data()
        req(fa)
        return(list(seqs = fragpipe_fasta_sequences(fa), source = "upload"))
      }
      folder <- loaded_psm_folder()
      req(folder)
      fas <- list.files(
        folder,
        pattern = "^protein\\.fas$",
        full.names = TRUE,
        recursive = TRUE
      )
      if (length(fas) == 0) {
        return(NULL)
      }
      seqs <- unlist(lapply(fas, function(f) {
        fragpipe_fasta_sequences(Biostrings::readAAStringSet(f))
      }))
      list(seqs = seqs[!duplicated(names(seqs))], source = "protein.fas")
    })

    output$ibaq_banner <- renderUI({
      msg <- if (is.null(loaded_psm_folder())) {
        "Load a FragPipe folder in the sidebar to enable iBAQ calculation (peptide.tsv is read from the same subfolders as psm.tsv)."
      } else if (is.null(input$fasta_file)) {
        "No FASTA uploaded: theoretical peptides are counted from the protein.fas files in the FragPipe folder. Upload a FASTA in the sidebar to use it instead."
      }
      if (!is.null(msg)) {
        div(
          class = "plot-placeholder",
          icon("circle-info", style = "margin-right:10px;"),
          msg
        )
      }
    })

    # Untransformed iBAQ in long format. peptide.tsv is already FDR-filtered by
    # FragPipe, so the sidebar Hyperscore / probability / specificity filters
    # (PSM-level) do not apply here; the sample selector does.
    ibaq_long <- eventReactive(input$run_ibaq, {
      if (is.null(loaded_psm_folder())) {
        showNotification(
          "Load a FragPipe folder to compute iBAQ.",
          type = "warning"
        )
        req(FALSE)
      }
      num <- function(x, default) if (isTruthy(x)) as.integer(x) else default
      min_len <- max(1L, num(input$ibaq_min_len, 6L))
      max_len <- max(min_len, num(input$ibaq_max_len, 30L))
      tryCatch(
        withProgress(message = "Computing iBAQ…", value = 0.3, {
          pep <- peptide_data()
          sel <- selected_experiments()
          if (!is.null(sel)) {
            pep <- pep[pep$experiment %in% sel, , drop = FALSE]
          }
          fa <- ibaq_sequences()
          if (is.null(fa)) {
            stop(
              "no protein.fas found in the FragPipe folder; ",
              "upload a FASTA file in the sidebar."
            )
          }
          protease_regex <- protease_regex_from_rule(
            input$ibaq_cut,
            input$ibaq_nocut,
            input$ibaq_sense
          )
          if (is.null(protease_regex)) {
            stop("enter at least one cleaved residue (e.g. KR for trypsin).")
          }
          compute_ibaq_fragpipe(
            pep,
            fa$seqs,
            unique_only = isTRUE(input$ibaq_unique),
            exclude_contaminants = isTRUE(input$ibaq_no_contam),
            protease_regex = protease_regex,
            min_peptide_len = min_len,
            max_peptide_len = max_len,
            max_missed_cleavages = min(
              5L,
              max(0L, num(input$ibaq_missed_cleavages, 0L))
            )
          )
        }),
        error = function(e) {
          showNotification(
            paste("iBAQ error:", conditionMessage(e)),
            type = "error",
            duration = 10
          )
          NULL
        }
      )
    })

    ibaq_top_n <- reactive({
      n <- suppressWarnings(as.integer(input$ibaq_top_n))
      if (length(n) == 0 || is.na(n)) 20L else max(1L, min(100L, n))
    })

    # Labels per facet in the ranked-abundance plot; 0 turns labels off
    ibaq_rank_labels <- reactive({
      n <- suppressWarnings(as.integer(input$ibaq_rank_labels))
      if (length(n) == 0 || is.na(n)) 5L else max(0L, min(50L, n))
    })

    # Per-sample rank and relative iBAQ (riBAQ = share of the sample's total
    # iBAQ), with a gene label that falls back to the protein ID.
    ibaq_per_sample <- reactive({
      d <- ibaq_long()
      req(d, nrow(d) > 0)
      labels <- d[!duplicated(d$protein_id), ]
      lab <- labels$protein_id
      has_gene <- !is.na(labels$gene) & labels$gene != ""
      lab[has_gene] <- labels$gene[has_gene]
      d$label <- make.unique(lab)[match(d$protein_id, labels$protein_id)]
      d %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(
          riBAQ = iBAQ / sum(iBAQ) * 100,
          rank = dplyr::row_number(dplyr::desc(iBAQ))
        ) %>%
        dplyr::ungroup()
    })

    # Facet layout and pixel heights for the per-sample plots
    IBAQ_NCOL <- 3L
    ibaq_rank_h <- reactive({
      n <- dplyr::n_distinct(ibaq_per_sample()$sample)
      max(500L, ceiling(n / IBAQ_NCOL) * 300L)
    })
    ibaq_top_h <- reactive({
      n <- dplyr::n_distinct(ibaq_per_sample()$sample)
      max(500L, ceiling(n / IBAQ_NCOL) * (ibaq_top_n() * 14L + 70L))
    })

    output$ibaq_rank_plot_ui <- renderUI({
      plotOutput(ns("ibaq_rank_plot"), height = paste0(ibaq_rank_h(), "px"))
    })
    output$ibaq_top_plot_ui <- renderUI({
      plotOutput(ns("ibaq_top_plot"), height = paste0(ibaq_top_h(), "px"))
    })

    ibaq_wide <- reactive({
      d <- ibaq_long()
      req(d, nrow(d) > 0)
      if (isTRUE(input$ibaq_log2)) {
        d$iBAQ <- log2(d$iBAQ + 1)
      }
      d %>%
        dplyr::arrange(sample) %>%
        tidyr::pivot_wider(names_from = sample, values_from = iBAQ) %>%
        dplyr::arrange(protein_id)
    })

    output$ibaq_summary <- renderUI({
      d <- ibaq_long()
      req(d)
      n_unmatched <- length(attr(d, "unmatched_proteins"))
      n_contam <- length(attr(d, "contaminant_proteins"))
      tags$p(
        style = "margin:6px 0 0 0;",
        icon("circle-info"),
        sprintf(
          " iBAQ computed for %s proteins across %s samples.",
          scales::comma(dplyr::n_distinct(d$protein_id)),
          dplyr::n_distinct(d$sample)
        ),
        if (n_contam > 0) {
          sprintf(
            " %s contaminant protein(s) were excluded.",
            scales::comma(n_contam)
          )
        },
        if (n_unmatched > 0) {
          sprintf(
            " %s protein(s) were excluded because they are not in the FASTA (e.g. contaminants, or a FASTA that differs from the search database).",
            scales::comma(n_unmatched)
          )
        }
      )
    })

    ibaq_dist_obj <- reactive({
      d <- ibaq_long()
      req(d, nrow(d) > 0)
      ggplot(d, aes(x = sample, y = log2(iBAQ + 1))) +
        geom_boxplot(fill = input$plot_color, alpha = 0.75) +
        labs(x = NULL, y = "log₂(iBAQ + 1)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })

    ibaq_rank_obj <- reactive({
      s <- ibaq_per_sample()
      ggplot(s, aes(x = rank, y = log2(iBAQ + 1))) +
        geom_point(color = input$plot_color, alpha = 0.6, size = 0.8) +
        ggrepel::geom_text_repel(
          data = dplyr::filter(s, rank <= ibaq_rank_labels()),
          aes(label = label),
          size = 3,
          max.overlaps = Inf,
          min.segment.length = 0
        ) +
        labs(x = "Protein rank", y = "log₂(iBAQ + 1)") +
        facet_wrap(~sample, ncol = IBAQ_NCOL)
    })

    ibaq_top_obj <- reactive({
      top <- ibaq_per_sample() %>%
        dplyr::filter(rank <= ibaq_top_n())
      ggplot(
        top,
        aes(y = tidytext::reorder_within(label, riBAQ, sample), x = riBAQ)
      ) +
        geom_col(fill = input$plot_color, color = "black", alpha = 0.85) +
        tidytext::scale_y_reordered() +
        labs(y = NULL, x = "Relative iBAQ (% of sample total)") +
        facet_wrap(~sample, ncol = IBAQ_NCOL, scales = "free_y")
    })

    output$ibaq_dist_plot <- renderPlot({
      ibaq_dist_obj()
    })
    output$ibaq_rank_plot <- renderPlot({
      ibaq_rank_obj()
    })
    output$ibaq_top_plot <- renderPlot({
      ibaq_top_obj()
    })

    output$ibaq_table <- DT::renderDataTable(
      {
        tbl <- ibaq_wide()
        DT::datatable(
          tbl,
          rownames = FALSE,
          filter = "top",
          options = list(pageLength = 20, scrollX = TRUE)
        ) %>%
          DT::formatRound(
            columns = unique(ibaq_long()$sample),
            digits = 2
          )
      },
      server = TRUE
    )

    output$download_ibaq_table <- downloadHandler(
      filename = function() {
        paste0(
          if (isTRUE(input$ibaq_log2)) "log2_" else "",
          "iBAQ_matrix_",
          Sys.Date(),
          ".tsv"
        )
      },
      content = function(file) {
        data.table::fwrite(ibaq_wide(), file = file, sep = "\t", na = "NA")
      }
    )

    output$download_all_ibaq_plots <- downloadHandler(
      filename = function() {
        paste0("PSManalyst_iBAQ_plots_", Sys.Date(), ".zip")
      },
      content = function(file) {
        td <- tempdir()
        fps <- character()
        save_p <- function(nm, fn, w = 10, h = 8) {
          fp <- file.path(td, nm)
          tryCatch(
            {
              ggsave(
                fp,
                fn(),
                width = w,
                height = h,
                dpi = 300,
                bg = "white",
                limitsize = FALSE
              )
              fps <<- c(fps, fp)
            },
            error = function(e) message("Skip ", nm, ": ", conditionMessage(e))
          )
        }

        withProgress(message = "Saving plots...", value = 0, {
          save_p("iBAQ_distribution.png", function() ibaq_dist_obj(), w = 12)
          incProgress(1 / 3)
          # Same proportions as on screen (px / 80), capped for ggsave
          save_p(
            "iBAQ_ranked_abundance.png",
            function() ibaq_rank_obj(),
            w = 16,
            h = min(45, ibaq_rank_h() / 80)
          )
          incProgress(1 / 3)
          save_p(
            "iBAQ_top_proteins.png",
            function() ibaq_top_obj(),
            w = 16,
            h = min(45, ibaq_top_h() / 80)
          )
          incProgress(1 / 3)
        })

        req(length(fps) > 0)
        utils::zip(file, files = fps, flags = "-j")
      },
      contentType = "application/zip"
    )

    # ════════════════════════════════════════════════════════════════════════
    # G. DISPATCHERS
    # ════════════════════════════════════════════════════════════════════════
    current_psm_plot <- eventReactive(
      input$run_psm_plot,
      {
        req(input$psm_plot_select)
        show_spinner("sp_psm_main")
        p <- NULL
        if (input$psm_plot_select == "plot_fdr_curve") {
          p <- fdr_curve_obj()
        } else if (exists(input$psm_plot_select, envir = psm_fns)) {
          p <- get(input$psm_plot_select, envir = psm_fns)()
        }
        p
      },
      ignoreNULL = FALSE
    )

    output$psm_dynamic_plot_out <- renderPlot({
      on.exit(hide_spinner("sp_psm_main"), add = TRUE)
      req(current_psm_plot())
      current_psm_plot()
    })

    current_prot_plot <- eventReactive(
      input$run_prot_plot,
      {
        req(input$prot_plot_select)
        show_spinner("sp_prot_main")
        p <- NULL
        if (exists(input$prot_plot_select, envir = prot_fns)) {
          p <- get(input$prot_plot_select, envir = prot_fns)()
        }
        p
      },
      ignoreNULL = FALSE
    )

    output$prot_dynamic_plot_out <- renderPlot({
      on.exit(hide_spinner("sp_prot_main"), add = TRUE)
      req(current_prot_plot())
      current_prot_plot()
    })

    # ════════════════════════════════════════════════════════════════════════
    # H. DOWNLOADS
    # ════════════════════════════════════════════════════════════════════════
    output$download_psm_plot <- downloadHandler(
      filename = function() {
        paste0("PSM_", input$psm_plot_select, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        p <- current_psm_plot()
        req(p)
        ggsave(
          file,
          p,
          width = 12,
          height = 8,
          bg = "white",
          device = "png",
          dpi = 300
        )
      }
    )

    output$download_prot_plot <- downloadHandler(
      filename = function() {
        paste0("Protein_", input$prot_plot_select, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        p <- current_prot_plot()
        req(p)
        ggsave(
          file,
          p,
          width = 12,
          height = 8,
          bg = "white",
          device = "png",
          dpi = 300
        )
      }
    )

    # ── Download all plots as ZIP ──
    output$download_all_plots <- downloadHandler(
      filename = function() paste0("PSManalyst_all_plots_", Sys.Date(), ".zip"),
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
          # PSM plots
          if (!is.null(tryCatch(data(), error = function(e) NULL))) {
            psm_plot_ids <- c(
              "plot01",
              "plot03",
              "plot04",
              "plot05",
              "plot06",
              "plot07",
              "plot08",
              "plot09",
              "plot10",
              "plot11",
              "plot12",
              "plot13",
              "plot14",
              "plot15",
              "plot16",
              "plot17",
              "plot18",
              "plot19",
              "plot20"
            )
            for (pid in psm_plot_ids) {
              if (exists(pid, envir = psm_fns)) {
                fn <- get(pid, envir = psm_fns)
                save_p(paste0("PSM_", pid, ".png"), fn, w = 12, h = 8)
              }
            }
            # FDR
            try(save_p(
              "PSM_FDR_curve.png",
              function() fdr_curve_obj(),
              w = 12,
              h = 8
            ))
            # AA freq (requires fasta)
            if (!is.null(input$fasta_file)) {
              try(save_p(
                "PSM_plot_aa_frequency.png",
                function() psm_fns$plot_aa_freq(),
                w = 12,
                h = 8
              ))
            }
            # Mod diagnostic
            try(save_p(
              "PSM_Mod_Diagnostic.png",
              function() mod_diag_plot_obj(),
              w = 14,
              h = 10
            ))
            try({
              paired <- mod_diag_data()
              if (!is.null(paired) && nrow(paired) > 0) {
                rt_csv <- file.path(td, "PSM_RT_Shift_Table.csv")
                data.table::fwrite(paired, file = rt_csv, sep = ",", na = "NA")
                fps <<- c(fps, rt_csv)
              }
            })
          }
          incProgress(1 / 3)

          # Protein plots (from protein.tsv)
          if (!is.null(tryCatch(protein_data(), error = function(e) NULL))) {
            prot_ids <- c(
              "plot01p",
              "plot02p",
              "plot03p",
              "plot04p",
              "plot05p",
              "plot06p",
              "plot07p",
              "plot08p"
            )
            for (pid in prot_ids) {
              if (exists(pid, envir = prot_fns)) {
                fn <- get(pid, envir = prot_fns)
                save_p(paste0("Protein_", pid, ".png"), fn, w = 10, h = 8)
              }
            }
          }
          incProgress(1 / 3)

          # Combined protein plots (from combined_protein.tsv)
          if (!is.null(input$combined_protein)) {
            try(save_p(
              "Protein_plot10p_maxlfq.png",
              function() prot_fns$plot10p(),
              w = 10,
              h = 8
            ))
            try(save_p(
              "Protein_ggpairs.png",
              function() prot_fns$plot_ggpairs(),
              w = 14,
              h = 14
            ))
            try(save_p(
              "Protein_cosine_similarity.png",
              function() prot_fns$cosine_similarity(),
              w = 10,
              h = 10
            ))
            try(save_p(
              "Protein_euclidean_distance.png",
              function() prot_fns$euclidean_distance(),
              w = 10,
              h = 10
            ))
            try(save_p(
              "Protein_jaccard_similarity.png",
              function() prot_fns$jaccard_similarity(),
              w = 10,
              h = 10
            ))
            try(save_p(
              "Protein_rank_lfq.png",
              function() prot_fns$plot_rank_lfq(),
              w = 14,
              h = 10
            ))
            try(save_p(
              "Protein_rank_usc.png",
              function() prot_fns$plot_rank_usc(),
              w = 14,
              h = 10
            ))
          }
          incProgress(1 / 3)
        })

        utils::zip(file, files = fps, flags = "-j")
      },
      contentType = "application/zip"
    )

    # ── Abundance Matrix Export ──────────────────────────────────────────────
    output$download_abundance_matrix <- downloadHandler(
      filename = function() {
        metric <- input$export_abundance_metric %||% "max_lfq_intensity"
        paste0("abundance_matrix_", metric, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        req(input$combined_protein)
        metric <- input$export_abundance_metric %||% "max_lfq_intensity"
        id_col <- input$export_protein_id %||% "protein_id"

        raw <- data.table::fread(
          input$combined_protein$datapath,
          sep = "\t",
          header = TRUE
        ) %>%
          janitor::clean_names()

        # Validate that the chosen identifier column exists
        if (!id_col %in% colnames(raw)) {
          showNotification(
            paste0("Column '", id_col, "' not found in combined_protein.tsv."),
            type = "error",
            duration = 8
          )
          return()
        }

        # Select abundance columns based on the user chosen metric
        if (metric == "max_lfq_intensity") {
          abundance_cols <- grep(
            "max_lfq_intensity$",
            colnames(raw),
            value = TRUE,
            ignore.case = TRUE
          )
          abundance_cols <- abundance_cols[
            !grepl("combined", abundance_cols, ignore.case = TRUE)
          ]
          suffix_pattern <- "_max_lfq_intensity$"
        } else {
          all_intensity <- grep(
            "_intensity$",
            colnames(raw),
            value = TRUE,
            ignore.case = TRUE
          )
          max_lfq <- grep(
            "max_lfq_intensity$",
            colnames(raw),
            value = TRUE,
            ignore.case = TRUE
          )
          abundance_cols <- setdiff(all_intensity, max_lfq)
          abundance_cols <- abundance_cols[
            !grepl("combined", abundance_cols, ignore.case = TRUE)
          ]
          suffix_pattern <- "_intensity$"
        }

        if (length(abundance_cols) == 0) {
          showNotification(
            paste0("No '", metric, "' columns found in the data."),
            type = "error",
            duration = 8
          )
          return()
        }

        # Build matrix
        export_df <- raw %>%
          dplyr::select(dplyr::all_of(c(id_col, abundance_cols)))

        export_df <- export_df %>%
          dplyr::mutate(
            dplyr::across(
              dplyr::all_of(abundance_cols),
              ~ dplyr::if_else(.x == 0, NA_real_, as.numeric(.x))
            )
          )

        clean_names <- stringr::str_remove(
          abundance_cols,
          stringr::regex(suffix_pattern, ignore_case = TRUE)
        )
        colnames(export_df) <- c(id_col, clean_names)

        data.table::fwrite(export_df, file = file, sep = "\t", na = "NA")
      }
    )
  })
}
