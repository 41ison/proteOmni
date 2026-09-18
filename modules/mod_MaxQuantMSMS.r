# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  proteOmni — MaxQuant Evidence / Peptides / MS/MS Scans QC Module         ║
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


# ── MaxQuant folder discovery ─────────────────────────────────────────────────

#' Locate the MaxQuant output files under a user-supplied folder.
#'
#' Accepts either the \code{combined} folder or the \code{combined/txt}
#' subfolder; searches recursively and prefers matches inside a \code{txt}
#' directory when more than one exists.
#'
#' For the heavy tables (\code{msmsScans}, \code{evidence}) the \code{.parquet}
#' file produced by \code{ensure_mq_parquet()} takes precedence; the
#' \code{.txt} is returned only when no parquet exists yet, and is then used
#' solely to build the parquet.
#'
#' @param path Character. Folder to search.
#' @return Named list with elements \code{msmsScans}, \code{evidence},
#'   \code{summary}, \code{proteinGroups}, \code{peptides} (full path or
#'   \code{NA}). msms.txt is deliberately not used (too large).
MQ_FILE_TARGETS <- c(
  msmsScans = "msmsScans",
  evidence = "evidence",
  summary = "summary",
  proteinGroups = "proteinGroups",
  peptides = "peptides"
)

find_mq_files <- function(path) {
  targets <- MQ_FILE_TARGETS
  path <- path.expand(trimws(path))
  if (!nzchar(path) || !dir.exists(path)) {
    return(setNames(
      as.list(rep(NA_character_, length(targets))),
      names(targets)
    ))
  }
  lapply(targets, function(f) {
    exts <- if (f %in% MQ_PARQUET_TABLES) "(txt|parquet)" else "txt"
    hits <- list.files(
      path,
      pattern = paste0("^", f, "\\.", exts, "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(hits) == 0) {
      return(NA_character_)
    }
    in_txt <- grepl("/txt/", hits, fixed = TRUE)
    if (any(in_txt)) {
      hits <- hits[in_txt]
    }
    # Keep a single directory, then prefer the .parquet within it
    hits <- hits[dirname(hits) == dirname(hits[1])]
    is_pq <- grepl("\\.parquet$", hits)
    if (any(is_pq)) hits[is_pq][1] else hits[1]
  })
}


# ── Parquet cache for heavy tables ────────────────────────────────────────────

#' Tables that are cached as parquet next to the original .txt.
MQ_PARQUET_TABLES <- c("msmsScans", "evidence")

#' Path of the parquet sibling of a MaxQuant .txt file.
mq_parquet_path <- function(txt_path) {
  sub("\\.txt$", ".parquet", txt_path)
}

#' Make sure a parquet copy of a MaxQuant .txt exists.
#'
#' If the parquet already exists it is used as-is. Otherwise the full .txt is
#' converted once (all columns), the in-memory txt table is released
#' immediately, and the parquet is written next to the original. The file is
#' written to a temp name and renamed so a failed conversion never leaves a
#' truncated parquet behind.
#'
#' The parquet is the only format the app works with for these tables, so a
#' conversion failure (e.g. read-only folder) is an error rather than a
#' fallback to the .txt.
#'
#' @param txt_path Character. Full path to the .txt file.
#' @param log Function \code{(msg, level)} used to report progress.
#' @return Parquet path.
ensure_mq_parquet <- function(txt_path, log = function(msg, level) NULL) {
  pq_path <- mq_parquet_path(txt_path)
  if (file.exists(pq_path)) {
    return(pq_path)
  }
  if (!file.exists(txt_path)) {
    stop("Neither ", basename(pq_path), " nor ", basename(txt_path), " found.")
  }

  log(
    sprintf(
      "Converting %s to parquet (one-time; %s MB)…",
      basename(txt_path),
      round(file.info(txt_path)$size / 1024^2, 1)
    ),
    "info"
  )
  tmp <- tempfile(
    pattern = paste0(".", basename(pq_path), "-"),
    tmpdir = dirname(pq_path),
    fileext = ".tmp"
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

  t0 <- Sys.time()
  tryCatch(
    {
      dt <- data.table::fread(txt_path)
      arrow::write_parquet(dt, tmp)
      # Release the txt table before anything else is loaded
      rm(dt)
      invisible(gc())
      if (!file.rename(tmp, pq_path)) {
        stop("could not move the temporary file into place")
      }
    },
    error = function(e) {
      stop(
        "Could not create ",
        basename(pq_path),
        " in ",
        dirname(pq_path),
        ": ",
        conditionMessage(e),
        ". The folder must be writable so proteOmni can store the parquet ",
        "next to the original .txt.",
        call. = FALSE
      )
    }
  )
  log(
    sprintf(
      "%s written (%s MB) in %s s",
      basename(pq_path),
      round(file.info(pq_path)$size / 1024^2, 1),
      round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    ),
    "ok"
  )
  pq_path
}

#' Read selected columns of a heavy MaxQuant table from its parquet file.
#'
#' The .txt is never read for analysis: when a .txt path is given it is first
#' converted with \code{ensure_mq_parquet()} and the parquet is read instead.
#'
#' @param path   Character. Path to a .parquet file, or to the .txt when no
#'   parquet exists yet.
#' @param select Character. Columns to read; missing ones are skipped.
#' @param log    Function \code{(msg, level)} for progress messages.
#' @return A data.table.
read_mq_table <- function(path, select, log = function(msg, level) NULL) {
  pq_path <- if (grepl("\\.parquet$", path)) {
    path
  } else {
    ensure_mq_parquet(path, log)
  }

  available <- arrow::ParquetFileReader$create(pq_path)$GetSchema()$names
  cols <- intersect(select, available)
  missing <- setdiff(select, available)
  if (length(missing)) {
    log(
      sprintf(
        "%s: column(s) not found, skipping: %s",
        basename(pq_path),
        paste(missing, collapse = ", ")
      ),
      "warn"
    )
  }
  arrow::read_parquet(pq_path, col_select = tidyselect::all_of(cols)) |>
    data.table::as.data.table()
}


# ── summary.txt helpers ───────────────────────────────────────────────────────

#' Read the MaxQuant summary.txt file, dropping columns that are entirely empty.
#' @param path Character. Full path to summary.txt.
#' @return A data.table with one row per raw file (plus the MaxQuant Total row).
read_summary_file <- function(path) {
  dt <- data.table::fread(path)
  empty <- vapply(
    dt,
    function(x) all(is.na(x) | (is.character(x) & !nzchar(x))),
    logical(1)
  )
  dt[, .SD, .SDcols = !empty]
}


# ── proteinGroups.txt helpers ─────────────────────────────────────────────────

#' Quantity types that can be exported from proteinGroups.txt.
#' Names are display labels, values are the column prefixes.
MQ_QUANT_TYPES <- c(
  "Intensity" = "Intensity",
  "LFQ intensity" = "LFQ intensity",
  "iBAQ" = "iBAQ"
)

#' Read and filter proteinGroups.txt.
#'
#' Removes reverse hits, potential contaminants and groups only identified by
#' site, and adds a \code{Protein ID} column holding the first accession of
#' each protein group.
#'
#' @param path Character. Full path to proteinGroups.txt.
#' @return Filtered data.table with all original columns plus \code{Protein ID}.
read_protein_groups <- function(path) {
  pg <- data.table::fread(path)

  flagged <- function(col) {
    if (!col %in% names(pg)) {
      return(rep(FALSE, nrow(pg)))
    }
    x <- pg[[col]]
    !is.na(x) & as.character(x) == "+"
  }
  keep <- !(flagged("Reverse") |
    flagged("Potential contaminant") |
    flagged("Only identified by site"))
  pg <- pg[keep]
  pg[, `Protein ID` := sub(";.*$", "", `Protein IDs`)]
  pg
}


#' Per-sample column names for a quantity type.
#'
#' Sample-level columns are "<quant> <sample>"; the bare "<quant>" column is
#' the summed total and "iBAQ peptides" is a count, so both are excluded.
quant_sample_cols <- function(pg, quant) {
  cols <- grep(paste0("^", quant, " .+"), names(pg), value = TRUE)
  setdiff(cols, paste0(quant, " peptides"))
}


#' Quantity types present (with per-sample columns) in a proteinGroups table.
#' @return Named character vector, subset of \code{MQ_QUANT_TYPES}.
available_quant_types <- function(pg) {
  present <- vapply(
    MQ_QUANT_TYPES,
    function(q) length(quant_sample_cols(pg, q)) > 0,
    logical(1)
  )
  MQ_QUANT_TYPES[present]
}


#' Build a protein abundance matrix from a filtered proteinGroups table.
#'
#' @param pg         data.table from \code{read_protein_groups()}.
#' @param quant      Character. One of \code{MQ_QUANT_TYPES}.
#' @param zero_to_na Logical. Replace 0 (MaxQuant's "not quantified") with NA.
#' @param log2_transform Logical. Apply log2 to the abundance columns. Zeros
#'   are always converted to NA first, since log2(0) is -Inf.
#' @return A data.table with a \code{Protein ID} column followed by one
#'   abundance column per sample.
build_protein_matrix <- function(
  pg,
  quant = "Intensity",
  zero_to_na = FALSE,
  log2_transform = FALSE
) {
  cols <- quant_sample_cols(pg, quant)
  if (length(cols) == 0) {
    stop("No per-sample '", quant, "' columns found in proteinGroups.txt.")
  }
  out <- cbind(
    data.table::data.table(`Protein ID` = pg[["Protein ID"]]),
    pg[, .SD, .SDcols = cols]
  )
  data.table::setnames(out, cols, sub(paste0("^", quant, " "), "", cols))
  val_cols <- setdiff(names(out), "Protein ID")

  if (zero_to_na || log2_transform) {
    out[,
      (val_cols) := lapply(.SD, function(x) {
        x <- as.numeric(x)
        x[!is.na(x) & x == 0] <- NA_real_
        x
      }),
      .SDcols = val_cols
    ]
  }
  if (log2_transform) {
    out[, (val_cols) := lapply(.SD, log2), .SDcols = val_cols]
  }
  out
}


# ── peptides.txt helpers ──────────────────────────────────────────────────────

#' Read and filter the MaxQuant peptides.txt file.
#'
#' One row per non-redundant peptide sequence. Reverse hits and potential
#' contaminants are removed. Small enough to be read directly with fread.
#'
#' @param path Character. Full path to peptides.txt.
#' @return A data.table.
read_peptides_file <- function(path) {
  pep <- data.table::fread(path)
  flagged <- function(col) {
    if (!col %in% names(pep)) {
      return(rep(FALSE, nrow(pep)))
    }
    x <- pep[[col]]
    !is.na(x) & as.character(x) == "+"
  }
  pep <- pep[!(flagged("Reverse") | flagged("Potential contaminant"))]
  if ("Missed cleavages" %in% names(pep)) {
    pep[, `Missed cleavages` := as.character(`Missed cleavages`)]
  }
  pep
}


#' Experiment names present in peptides.txt ("Experiment <name>" columns).
#' @return Character vector (possibly empty).
pep_experiments <- function(pep) {
  cols <- grep("^Experiment .+", names(pep), value = TRUE)
  sub("^Experiment ", "", cols)
}


#' Peptide summary table shown in the "Peptides" tab.
pep_summary_table <- function(pep) {
  wanted <- c(
    "Sequence",
    "Gene names",
    "Leading razor protein",
    "Length",
    "Missed cleavages",
    "Charges",
    "MS/MS Count",
    "Score",
    "PEP",
    "Unique (Groups)",
    "Unique (Proteins)",
    "Intensity"
  )
  pep[, .SD, .SDcols = intersect(wanted, names(pep))] |>
    dplyr::arrange(dplyr::desc(`MS/MS Count`))
}


#' Long table of per-experiment peptide counts (rows where the peptide was
#' identified in that experiment).
pep_experiment_long <- function(pep) {
  exps <- pep_experiments(pep)
  if (length(exps) == 0) {
    return(NULL)
  }
  pep |>
    dplyr::select(Sequence, dplyr::all_of(paste0("Experiment ", exps))) |>
    tidyr::pivot_longer(
      -Sequence,
      names_to = "Experiment",
      names_prefix = "Experiment ",
      values_to = "n_msms"
    ) |>
    dplyr::filter(!is.na(n_msms), n_msms > 0)
}


# ── Individual peptides.txt plot builders ─────────────────────────────────────

plot_pep_msms_count <- function(pep) {
  cap <- 20L
  pep |>
    dplyr::mutate(
      n = pmin(`MS/MS Count`, cap),
      n = factor(
        ifelse(n >= cap, paste0(cap, "+"), as.character(n)),
        levels = c(as.character(0:(cap - 1)), paste0(cap, "+"))
      )
    ) |>
    ggplot(aes(x = n)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    labs(x = "MS/MS spectra per peptide", y = "Number of peptides") +
    theme_mq()
}

plot_pep_score <- function(pep) {
  ggplot(pep, aes(x = Score)) +
    geom_density(
      fill = "#1b9e77",
      alpha = 0.8,
      color = "white",
      linewidth = 0.25
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(x = "Andromeda score", y = "Density") +
    theme_mq()
}

plot_pep_pep <- function(pep) {
  pep |>
    dplyr::filter(!is.na(PEP), PEP > 0) |>
    ggplot(aes(x = PEP)) +
    geom_density(
      fill = "#1b9e77",
      alpha = 0.8,
      color = "white",
      linewidth = 0.25
    ) +
    scale_x_log10() +
    labs(x = "Posterior error probability (PEP, log scale)", y = "Density") +
    theme_mq()
}

plot_pep_charges <- function(pep) {
  pep |>
    dplyr::select(Sequence, Charges) |>
    tidyr::separate_rows(Charges, sep = ";") |>
    dplyr::filter(nzchar(Charges)) |>
    ggplot(aes(x = Charges)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      vjust = -0.3,
      size = 5,
      fontface = "bold"
    ) +
    scale_x_discrete(limits = as.character(1:6)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x = "Charge state observed for the peptide",
      y = "Number of peptides"
    ) +
    theme_mq()
}

plot_pep_unique <- function(pep) {
  pep |>
    dplyr::select(
      Sequence,
      dplyr::any_of(c("Unique (Groups)", "Unique (Proteins)"))
    ) |>
    tidyr::pivot_longer(-Sequence, names_to = "Level", values_to = "Unique") |>
    ggplot(aes(x = Unique)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      vjust = -0.3,
      size = 5,
      fontface = "bold"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    facet_wrap(~Level) +
    labs(x = "Unique peptide", y = "Number of peptides") +
    theme_mq()
}

plot_pep_per_experiment <- function(pep) {
  long <- pep_experiment_long(pep)
  if (is.null(long)) {
    stop("peptides.txt has no 'Experiment <name>' columns.")
  }
  long |>
    dplyr::count(Experiment, name = "n_peptides") |>
    ggplot(aes(x = Experiment, y = n_peptides)) +
    geom_col(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = n_peptides),
      vjust = -0.3,
      size = 5,
      fontface = "bold"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Experiment", y = "Peptides identified") +
    theme_mq() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
}

plot_pep_experiment_overlap <- function(pep) {
  long <- pep_experiment_long(pep)
  if (is.null(long)) {
    stop("peptides.txt has no 'Experiment <name>' columns.")
  }
  n_exp <- length(pep_experiments(pep))
  long |>
    dplyr::count(Sequence, name = "n_experiments") |>
    dplyr::mutate(n_experiments = factor(n_experiments, levels = 1:n_exp)) |>
    ggplot(aes(x = n_experiments)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      vjust = -0.3,
      size = 5,
      fontface = "bold"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x = "Number of experiments in which the peptide was identified",
      y = "Number of peptides"
    ) +
    theme_mq()
}

plot_pep_intensity <- function(pep) {
  exps <- pep_experiments(pep)
  cols <- intersect(paste0("Intensity ", exps), names(pep))
  if (length(cols) == 0) {
    stop("peptides.txt has no per-experiment 'Intensity <name>' columns.")
  }
  pep |>
    dplyr::select(Sequence, dplyr::all_of(cols)) |>
    tidyr::pivot_longer(
      -Sequence,
      names_to = "Experiment",
      names_prefix = "Intensity ",
      values_to = "Intensity"
    ) |>
    # 0 means "not quantified" in MaxQuant
    dplyr::filter(!is.na(Intensity), Intensity > 0) |>
    ggplot(aes(x = log10(Intensity), fill = Experiment)) +
    geom_density(alpha = 0.5, color = "white", linewidth = 0.25) +
    scale_fill_viridis_d(option = "D") +
    labs(x = "log10 peptide intensity", y = "Density") +
    theme_mq()
}

plot_pep_length <- function(pep) {
  ggplot(pep, aes(x = Length)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(
      x = "Peptide length (number of amino acids)",
      y = "Number of peptides"
    ) +
    theme_mq()
}

plot_pep_missed_cleavages <- function(pep) {
  ggplot(pep, aes(x = `Missed cleavages`)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      vjust = -0.3,
      size = 5,
      fontface = "bold"
    ) +
    scale_x_discrete(limits = as.character(0:5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Number of missed cleavages", y = "Number of peptides") +
    theme_mq()
}

plot_pep_terminal_aa <- function(pep) {
  pep |>
    dplyr::select(
      Sequence,
      dplyr::any_of(c(
        "First amino acid",
        "Last amino acid",
        "Amino acid after"
      ))
    ) |>
    tidyr::pivot_longer(-Sequence, names_to = "Position", values_to = "aa") |>
    dplyr::filter(!is.na(aa), nzchar(aa)) |>
    ggplot(aes(x = aa)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    facet_wrap(~Position, ncol = 1, scales = "free_y") +
    labs(x = "Amino acid", y = "Number of peptides") +
    theme_mq()
}

plot_pep_per_protein <- function(pep) {
  cap <- 30L
  pep |>
    dplyr::filter(
      !is.na(`Leading razor protein`),
      nzchar(`Leading razor protein`)
    ) |>
    dplyr::count(`Leading razor protein`, name = "n_peptides") |>
    dplyr::mutate(
      n = pmin(n_peptides, cap),
      n = factor(
        ifelse(n >= cap, paste0(cap, "+"), as.character(n)),
        levels = c(as.character(1:(cap - 1)), paste0(cap, "+"))
      )
    ) |>
    ggplot(aes(x = n)) +
    geom_bar(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    labs(x = "Peptides per leading razor protein", y = "Number of proteins") +
    theme_mq()
}

#' UpSet-style plot of peptide sharing between experiments.
#'
#' Top panel: number of peptides per exact experiment combination
#' (intersection). Bottom panel: dot matrix marking which experiments make
#' up each combination. Built with ggplot2 + patchwork only.
plot_pep_upset <- function(pep, max_sets = 40L) {
  long <- pep_experiment_long(pep)
  if (is.null(long)) {
    stop("peptides.txt has no 'Experiment <name>' columns.")
  }
  exps <- pep_experiments(pep)

  combos <- long |>
    dplyr::group_by(Sequence) |>
    dplyr::summarise(
      members = list(sort(unique(Experiment))),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      combo = vapply(members, paste, character(1), collapse = " & ")
    ) |>
    dplyr::count(combo, members, name = "n_peptides") |>
    dplyr::arrange(dplyr::desc(n_peptides)) |>
    dplyr::slice_head(n = max_sets) |>
    dplyr::mutate(combo = factor(combo, levels = combo))

  matrix_df <- combos |>
    dplyr::select(combo, members) |>
    tidyr::unnest(members) |>
    dplyr::rename(Experiment = members) |>
    dplyr::mutate(present = TRUE) |>
    tidyr::complete(
      combo,
      Experiment = exps,
      fill = list(present = FALSE)
    ) |>
    dplyr::mutate(Experiment = factor(Experiment, levels = rev(exps)))

  set_sizes <- long |>
    dplyr::count(Experiment, name = "n") |>
    dplyr::mutate(Experiment = factor(Experiment, levels = rev(exps)))

  p_bars <- ggplot(combos, aes(x = combo, y = n_peptides)) +
    geom_col(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = n_peptides),
      vjust = -0.3,
      size = 4,
      fontface = "bold"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Peptides in intersection") +
    theme_mq() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank()
    )

  p_matrix <- ggplot(matrix_df, aes(x = combo, y = Experiment)) +
    geom_point(aes(color = present), size = 4) +
    geom_line(
      data = dplyr::filter(matrix_df, present),
      aes(group = combo),
      color = "black",
      linewidth = 0.8
    ) +
    scale_color_manual(
      values = c(`TRUE` = "black", `FALSE` = "grey85"),
      guide = "none"
    ) +
    labs(x = "Experiment combination", y = NULL) +
    theme_mq() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )

  p_sets <- ggplot(set_sizes, aes(x = n, y = Experiment)) +
    geom_col(fill = "grey50", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = n),
      hjust = 1.1,
      size = 3.5,
      color = "white",
      fontface = "bold"
    ) +
    scale_x_reverse(
      expand = expansion(mult = c(0.05, 0)),
      breaks = scales::pretty_breaks(n = 3)
    ) +
    labs(x = "Peptides per experiment", y = NULL) +
    theme_mq() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  patchwork::wrap_plots(
    patchwork::plot_spacer(),
    p_bars,
    p_sets,
    p_matrix,
    ncol = 2,
    widths = c(1, 4),
    heights = c(3, 1)
  )
}


#' Peptides QC plot registry: key -> (label, builder, height in px).
PEP_PLOTS <- list(
  msms_count = list("MS/MS Spectra per Peptide", plot_pep_msms_count, 550L),
  per_experiment = list(
    "Peptides per Experiment",
    plot_pep_per_experiment,
    550L
  ),
  experiment_overlap = list(
    "Experiment Overlap (count)",
    plot_pep_experiment_overlap,
    550L
  ),
  upset = list("Experiment Overlap (UpSet)", plot_pep_upset, 750L),
  intensity = list(
    "Peptide Intensity per Experiment",
    plot_pep_intensity,
    550L
  ),
  charges = list("Charge States", plot_pep_charges, 550L),
  unique = list("Unique Peptides", plot_pep_unique, 550L),
  score = list("Andromeda Score", plot_pep_score, 550L),
  pep = list("PEP Distribution", plot_pep_pep, 550L),
  length = list("Peptide Length", plot_pep_length, 550L),
  missed_cleavages = list("Missed Cleavages", plot_pep_missed_cleavages, 550L),
  terminal_aa = list("Terminal Amino Acids", plot_pep_terminal_aa, 900L),
  per_protein = list("Peptides per Protein", plot_pep_per_protein, 550L)
)

plot_registry_choices <- function(registry) {
  setNames(names(registry), vapply(registry, `[[`, character(1), 1))
}


# ── msmsScans.txt helpers ─────────────────────────────────────────────────────
#
# msmsScans.txt has one row per MS/MS scan acquired (identified or not), so it
# is the right table for instrument-level QC: identification rate, TIC, ion
# injection time, precursor sampling. It replaces msms.txt, which was too
# large to be practical. Reads through the parquet cache.

#' Read the MaxQuant msmsScans.txt file (selected columns).
#' @param path Character. Full path to msmsScans.txt (or .parquet).
#' @param log  Function \code{(msg, level)} for progress messages.
#' @return A data.table with an \code{Identified} factor column.
read_msmsscans_file <- function(path, log = function(msg, level) NULL) {
  dt <- read_mq_table(
    path,
    log = log,
    select = c(
      "Raw file",
      "Scan number",
      "Retention time",
      "Ion injection time",
      "Total ion current",
      "Base peak intensity",
      "Identified",
      "Sequence",
      "Length",
      "Filtered peaks",
      "m/z",
      "Charge",
      "Scan event number",
      "Precursor intensity",
      "Precursor apex fraction",
      "Score",
      "PEP",
      "Modifications"
    )
  )
  if ("Identified" %in% names(dt)) {
    dt[,
      Identified := factor(
        ifelse(
          !is.na(Identified) & Identified == "+",
          "Identified",
          "Not identified"
        ),
        levels = c("Identified", "Not identified")
      )
    ]
  }
  if ("Charge" %in% names(dt)) {
    dt[, Charge := as.character(Charge)]
  }
  dt
}


# ── Individual msmsScans plot builders ────────────────────────────────────────

SCAN_ID_COLORS <- c("Identified" = "#1b9e77", "Not identified" = "#d95f02")

plot_scan_id_rate <- function(sc) {
  sc |>
    dplyr::count(`Raw file`, Identified) |>
    dplyr::group_by(`Raw file`) |>
    dplyr::mutate(frac = n / sum(n)) |>
    dplyr::ungroup() |>
    ggplot(aes(y = `Raw file`, x = n, fill = Identified)) +
    geom_col(color = "white", linewidth = 0.25, alpha = 0.9) +
    geom_text(
      aes(
        label = sprintf("%s (%.1f%%)", format(n, big.mark = ","), 100 * frac)
      ),
      position = position_stack(vjust = 0.5),
      size = 4,
      fontface = "bold",
      color = "white"
    ) +
    scale_fill_manual(values = SCAN_ID_COLORS) +
    labs(x = "MS/MS scans", y = NULL, fill = NULL) +
    theme_mq()
}

plot_scan_rt_hist <- function(sc) {
  ggplot(sc, aes(x = `Retention time`, fill = Identified)) +
    geom_histogram(
      binwidth = 1,
      color = "white",
      linewidth = 0.1,
      alpha = 0.9
    ) +
    scale_fill_manual(values = SCAN_ID_COLORS) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(
      x = "Retention time (min)",
      y = "MS/MS scans per minute",
      fill = NULL
    ) +
    theme_mq()
}

plot_scan_id_rate_rt <- function(sc) {
  sc |>
    dplyr::mutate(rt_bin = floor(`Retention time`)) |>
    dplyr::group_by(`Raw file`, rt_bin) |>
    dplyr::summarise(
      id_rate = mean(Identified == "Identified"),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    ggplot(aes(x = rt_bin, y = id_rate)) +
    geom_line(color = "#1b9e77", linewidth = 0.7) +
    geom_point(aes(size = n), color = "#1b9e77", alpha = 0.6) +
    scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
    scale_size_area(max_size = 3) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(
      x = "Retention time (min)",
      y = "Identification rate",
      size = "Scans"
    ) +
    theme_mq()
}

plot_scan_tic_rt <- function(sc) {
  ggplot(sc, aes(x = `Retention time`, y = log10(`Total ion current`))) +
    ggpointdensity::geom_pointdensity(
      method = "kde2d",
      adjust = 3,
      size = 0.6
    ) +
    scale_color_viridis_c(option = "D", direction = -1) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(
      x = "Retention time (min)",
      y = "log10 total ion current",
      color = "Density"
    ) +
    theme_mq()
}

plot_scan_injection_time <- function(sc) {
  ggplot(sc, aes(x = `Ion injection time`, fill = Identified)) +
    geom_histogram(bins = 50, color = "white", linewidth = 0.1, alpha = 0.9) +
    scale_fill_manual(values = SCAN_ID_COLORS) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Ion injection time (ms)", y = "MS/MS scans", fill = NULL) +
    theme_mq()
}

plot_scan_injection_time_rt <- function(sc) {
  ggplot(sc, aes(x = `Retention time`, y = `Ion injection time`)) +
    ggpointdensity::geom_pointdensity(
      method = "kde2d",
      adjust = 3,
      size = 0.6
    ) +
    scale_color_viridis_c(option = "D", direction = -1) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(
      x = "Retention time (min)",
      y = "Ion injection time (ms)",
      color = "Density"
    ) +
    theme_mq()
}

plot_scan_base_peak <- function(sc) {
  ggplot(sc, aes(x = log10(`Base peak intensity`), fill = Identified)) +
    geom_density(alpha = 0.6, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = SCAN_ID_COLORS) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "log10 base peak intensity", y = "Density", fill = NULL) +
    theme_mq()
}

plot_scan_precursor_intensity <- function(sc) {
  sc |>
    dplyr::filter(!is.na(`Precursor intensity`), `Precursor intensity` > 0) |>
    ggplot(aes(x = log10(`Precursor intensity`), fill = Identified)) +
    geom_density(alpha = 0.6, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = SCAN_ID_COLORS) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "log10 precursor intensity", y = "Density", fill = NULL) +
    theme_mq()
}

plot_scan_apex_fraction <- function(sc) {
  sc |>
    dplyr::filter(!is.na(`Precursor apex fraction`)) |>
    ggplot(aes(x = `Precursor apex fraction`, fill = Identified)) +
    geom_density(alpha = 0.6, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = SCAN_ID_COLORS) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(
      x = "Precursor apex fraction (1 = sampled at elution apex)",
      y = "Density",
      fill = NULL
    ) +
    theme_mq()
}

plot_scan_charge <- function(sc) {
  ggplot(sc, aes(x = Charge, fill = Identified)) +
    geom_bar(color = "white", linewidth = 0.25, alpha = 0.9) +
    scale_fill_manual(values = SCAN_ID_COLORS) +
    scale_x_discrete(limits = as.character(0:6)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(
      x = "Precursor charge (0 = undetermined)",
      y = "MS/MS scans",
      fill = NULL
    ) +
    theme_mq()
}

plot_scan_event <- function(sc) {
  sc |>
    dplyr::group_by(`Raw file`, `Scan event number`) |>
    dplyr::summarise(
      id_rate = mean(Identified == "Identified"),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    ggplot(aes(x = `Scan event number`, y = id_rate)) +
    geom_col(fill = "#1b9e77", alpha = 0.8, color = "white", linewidth = 0.25) +
    geom_text(aes(label = n), vjust = -0.3, size = 3.5, fontface = "bold") +
    scale_y_continuous(
      labels = scales::label_percent(),
      expand = expansion(mult = c(0, 0.15))
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(
      x = "Scan event number (position in TopN cycle)",
      y = "Identification rate (label: number of scans)"
    ) +
    theme_mq()
}

plot_scan_filtered_peaks <- function(sc) {
  ggplot(sc, aes(x = `Filtered peaks`, fill = Identified)) +
    geom_density(alpha = 0.6, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = SCAN_ID_COLORS) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Filtered peaks per MS/MS spectrum", y = "Density", fill = NULL) +
    theme_mq()
}

plot_scan_score <- function(sc) {
  sc |>
    dplyr::filter(!is.na(Score)) |>
    ggplot(aes(x = Score, fill = Identified)) +
    geom_density(alpha = 0.6, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = SCAN_ID_COLORS) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Andromeda score", y = "Density", fill = NULL) +
    theme_mq()
}

plot_scan_mz_rt <- function(sc) {
  ggplot(sc, aes(x = `Retention time`, y = `m/z`)) +
    geom_point(aes(color = Identified), size = 0.4, alpha = 0.4) +
    scale_color_manual(values = SCAN_ID_COLORS) +
    facet_wrap(~`Raw file`, ncol = 3) +
    labs(x = "Retention time (min)", y = "Precursor m/z", color = NULL) +
    theme_mq() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

#' MS/MS scans QC plot registry: key -> (label, builder). All plots are
#' faceted by raw file, so the height comes from \code{facet_height_px()}.
SCAN_PLOTS <- list(
  id_rate = list("Identification Rate per Raw File", plot_scan_id_rate),
  rt_hist = list("MS/MS Scans over Retention Time", plot_scan_rt_hist),
  id_rate_rt = list(
    "Identification Rate over Retention Time",
    plot_scan_id_rate_rt
  ),
  mz_rt = list("Precursor m/z vs Retention Time", plot_scan_mz_rt),
  tic_rt = list("Total Ion Current vs Retention Time", plot_scan_tic_rt),
  injection_time = list("Ion Injection Time", plot_scan_injection_time),
  injection_time_rt = list(
    "Ion Injection Time vs Retention Time",
    plot_scan_injection_time_rt
  ),
  base_peak = list("Base Peak Intensity", plot_scan_base_peak),
  precursor_intensity = list(
    "Precursor Intensity",
    plot_scan_precursor_intensity
  ),
  apex_fraction = list("Precursor Apex Fraction", plot_scan_apex_fraction),
  charge = list("Precursor Charge", plot_scan_charge),
  scan_event = list(
    "Identification Rate by Scan Event (TopN)",
    plot_scan_event
  ),
  filtered_peaks = list(
    "Filtered Peaks per Spectrum",
    plot_scan_filtered_peaks
  ),
  score = list("Andromeda Score", plot_scan_score)
)


# ── evidence.txt helpers ──────────────────────────────────────────────────────

#' Read the MaxQuant evidence.txt file.
#'
#' Reads through the parquet cache (see \code{read_mq_table()}).
#'
#' @param path Character. Full path to evidence.txt (or evidence.parquet).
#' @param log  Function \code{(msg, level)} for progress messages.
#' @return A data.table.
read_evidence_file <- function(path, log = function(msg, level) NULL) {
  dt <- read_mq_table(
    path,
    log = log,
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
    labs(x = "Retention time (min)", y = "m/z", color = "Density") +
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
    labs(x = "m/z", y = "Number of data points", color = "Density") +
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

      # ── Data Input ────────────────────────────────────────────────────────
      tags$div(class = "sidebar-section-label", "MaxQuant Output Folder"),

      tags$div(
        style = "padding:0 8px;",
        textInput(
          ns("mq_folder"),
          "Path to MaxQuant 'combined' folder",
          value = "",
          placeholder = "/path/to/combined"
        ),
        tags$p(
          style = "color:#adb5bd;font-size:11px;margin-top:-6px;",
          "msmsScans.txt, evidence.txt, peptides.txt, summary.txt and ",
          "proteinGroups.txt are located automatically in the txt/ ",
          "subfolder. msmsScans and evidence are read from .parquet; on ",
          "first load the .txt files are converted once and the .parquet ",
          "saved next to them. msms.txt is not used."
        ),
        tags$div(
          style = "text-align:center;",
          actionButton(
            ns("load_files"),
            "Load MaxQuant Files",
            class = "btn-primary",
            style = "width:80%;font-weight:bold;margin-bottom:6px;"
          )
        ),
        uiOutput(ns("file_status"))
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
        selected = "mz_rt"
      ),

      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── Peptides QC Parameters ────────────────────────────────────────────
      tags$div(class = "sidebar-section-label", "Peptides QC Parameters"),

      selectInput(
        ns("pep_plot_select"),
        "Select Peptides Plot",
        choices = plot_registry_choices(PEP_PLOTS),
        selected = "msms_count"
      ),

      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── MS/MS Scans QC Parameters ─────────────────────────────────────────
      tags$div(class = "sidebar-section-label", "MS/MS Scans QC Parameters"),

      selectInput(
        ns("scan_plot_select"),
        "Select Scans Plot",
        choices = plot_registry_choices(SCAN_PLOTS),
        selected = "id_rate"
      ),

      tags$hr(style = "border-color:#2d3741;margin:4px 0;"),

      # ── Actions ───────────────────────────────────────────────────────────
      tags$div(
        style = "padding:0 8px;text-align:center;",
        actionButton(
          ns("run_evidence"),
          "Plot Evidence QC",
          class = "btn-primary",
          style = "width:80%;font-weight:bold;margin-top:2px;margin-bottom:6px;"
        ),
        actionButton(
          ns("run_peptides"),
          "Plot Peptides QC",
          class = "btn-primary",
          style = "width:80%;font-weight:bold;margin-top:2px;margin-bottom:6px;"
        ),
        actionButton(
          ns("run_scans"),
          "Plot MS/MS Scans QC",
          class = "btn-primary",
          style = "width:80%;font-weight:bold;margin-top:2px;margin-bottom:10px;"
        )
      ),

      tags$div(
        style = "padding:0 8px;",
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
          style = "width:100%;text-align:left;margin-bottom:6px;"
        ),
        downloadButton(
          ns("download_pep_plot"),
          "\u2B07 Peptides Plot (.pdf)",
          class = "dl-btn",
          style = "width:100%;text-align:left;margin-bottom:6px;"
        ),
        downloadButton(
          ns("download_pep_data"),
          "\u2B07 Peptides Data (.tsv)",
          class = "dl-btn",
          style = "width:100%;text-align:left;margin-bottom:6px;"
        ),
        downloadButton(
          ns("download_scan_plot"),
          "\u2B07 MS/MS Scans Plot (.pdf)",
          class = "dl-btn",
          style = "width:100%;text-align:left;margin-bottom:6px;"
        ),
        downloadButton(
          ns("download_scan_data"),
          "\u2B07 MS/MS Scans Data (.tsv)",
          class = "dl-btn",
          style = "width:100%;text-align:left;margin-bottom:6px;"
        ),
        tags$hr(style = "border-color:#2d3741;margin:8px 0 4px 0;"),
        selectInput(
          ns("pg_quant"),
          "Protein abundance values",
          choices = MQ_QUANT_TYPES,
          selected = "Intensity"
        ),
        checkboxInput(
          ns("pg_zero_na"),
          "Replace zeros with NA",
          value = FALSE
        ),
        checkboxInput(
          ns("pg_log2"),
          "log2-transform values (zeros become NA)",
          value = FALSE
        ),
        downloadButton(
          ns("download_protein_matrix"),
          "\u2B07 Protein Abundance Matrix (.tsv)",
          class = "dl-btn",
          style = "width:100%;text-align:left;"
        )
      )
    ),
    tags$hr(style = "border-color:#2d3741;margin:4px 0;"),
    tags$div(class = "sidebar-section-label", "Modification Diagnostic"),
    tags$div(
      style = "padding:0 16px;",
      sliderInput(
        ns("rt_tolerance"),
        "ΔRT artifact threshold (min)",
        min = 0,
        max = 5,
        value = 0.5,
        step = 0.1
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

    # ── Tab 1: Run Summary (summary.txt) ─────────────────────────────────
    tabPanel(
      title = tagList(icon("clipboard-list"), "Run Summary"),
      fluidRow(
        box(
          title = "Processing Status",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          uiOutput(ns("status_log_ui"))
        )
      ),
      fluidRow(
        box(
          title = "MaxQuant run summary (summary.txt)",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          DT::dataTableOutput(ns("run_summary_table"))
        )
      ),
      fluidRow(
        box(
          title = "Protein abundance matrix preview (proteinGroups.txt)",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          DT::dataTableOutput(ns("protein_matrix_table"))
        )
      )
    ),

    # ── Tab 2: Evidence QC ────────────────────────────────────────────────
    tabPanel(
      title = tagList(icon("bar"), "Evidence QC"),
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

    # ── Tab: Peptides QC (peptides.txt) ──────────────────────────────────
    tabPanel(
      title = tagList(icon("chart-column"), "Peptides QC"),
      fluidRow(
        box(
          title = uiOutput(ns("pep_plot_title")),
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_peptides"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("peptides_plot_ui"))
          )
        )
      )
    ),

    # ── Tab: MS/MS Scans QC (msmsScans.txt) ──────────────────────────────
    tabPanel(
      title = tagList(icon("wave-square"), "MS/MS Scans QC"),
      fluidRow(
        box(
          title = uiOutput(ns("scan_plot_title")),
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(
            class = "plot-wrap",
            tags$div(
              class = "spinner-overlay",
              id = ns("sp_scans"),
              icon("spinner", class = "fa-spin")
            ),
            uiOutput(ns("scans_plot_ui"))
          )
        )
      ),
      fluidRow(
        box(
          title = "MS/MS scans table preview (first 500 rows)",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          DT::dataTableOutput(ns("scans_table"))
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

    # ── Tab 5: Peptides Summary (peptides.txt) ────────────────────────────
    tabPanel(
      title = tagList(icon("list"), "Peptides Summary"),
      fluidRow(
        box(
          title = "Identified peptides (peptides.txt; reverse hits and contaminants removed)",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          DT::dataTableOutput(ns("summary_table"))
        )
      )
    ),

    # ── Tab 6: Modification Diagnostic ───────────────────────────────────────
    tabPanel(
      title = tagList(icon("flask"), "Modification Diagnostic"),
      fluidRow(
        box(
          title = "RT Shift Profile — Modified vs Unmodified Peptides",
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
          title = "Modification Diagnostic Table",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          DT::dataTableOutput(ns("mod_diag_table"))
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
    # 0. PROCESSING LOG + FILE LOADING
    # ════════════════════════════════════════════════════════════════════

    # Timestamped step log shown in the sidebar and the Run Summary tab so
    # the user can always tell whether the app is working or idle.
    status_log <- reactiveVal(
      data.frame(
        time = character(0),
        level = character(0),
        msg = character(0),
        stringsAsFactors = FALSE
      )
    )

    log_step <- function(
      msg,
      level = c("info", "ok", "warn", "error"),
      notify = TRUE
    ) {
      level <- match.arg(level)
      status_log(rbind(
        status_log(),
        data.frame(
          time = format(Sys.time(), "%H:%M:%S"),
          level = level,
          msg = msg,
          stringsAsFactors = FALSE
        )
      ))
      if (notify) {
        showNotification(
          msg,
          type = switch(
            level,
            info = "message",
            ok = "message",
            warn = "warning",
            error = "error"
          ),
          duration = if (level == "error") NULL else 4
        )
      }
    }

    mq_files <- reactiveVal(NULL)
    raw_scans_rv <- reactiveVal(NULL)
    raw_evidence_rv <- reactiveVal(NULL)
    raw_summary_rv <- reactiveVal(NULL)
    protein_groups_rv <- reactiveVal(NULL)
    peptides_rv <- reactiveVal(NULL)

    # Read one file with progress + log; returns NULL on failure
    load_one <- function(label, path, reader, step, n_steps) {
      if (is.na(path)) {
        log_step(
          paste0(label, " not found — related tabs will stay empty."),
          "warn"
        )
        return(NULL)
      }
      size_mb <- round(file.info(path)$size / 1024^2, 1)
      log_step(
        sprintf("Reading %s (%s MB)…", label, size_mb),
        "info",
        notify = FALSE
      )
      setProgress(
        value = (step - 1) / n_steps,
        message = sprintf("Step %d/%d: reading %s", step, n_steps, label),
        detail = sprintf("%s MB — large files may take a minute", size_mb)
      )
      t0 <- Sys.time()
      res <- tryCatch(reader(path), error = function(e) {
        log_step(
          sprintf("Failed to read %s: %s", label, conditionMessage(e)),
          "error"
        )
        NULL
      })
      if (!is.null(res)) {
        secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
        msg <- if (is.character(res)) {
          # Lazy table: only the parquet path is kept
          sprintf(
            "%s ready for on-demand queries (%s) in %s s",
            label,
            basename(res),
            secs
          )
        } else {
          sprintf(
            "%s loaded: %s rows in %s s",
            label,
            format(nrow(res), big.mark = ","),
            secs
          )
        }
        log_step(msg, "ok")
      }
      res
    }

    observeEvent(input$load_files, {
      shinyjs::disable("load_files")
      on.exit(shinyjs::enable("load_files"), add = TRUE)

      # Reset previous state
      raw_scans_rv(NULL)
      raw_evidence_rv(NULL)
      raw_summary_rv(NULL)
      protein_groups_rv(NULL)
      peptides_rv(NULL)
      status_log(status_log()[0, ])

      folder <- trimws(input$mq_folder)
      if (!nzchar(folder)) {
        log_step(
          "Please enter the path to the MaxQuant 'combined' folder.",
          "error"
        )
        return()
      }
      if (!dir.exists(path.expand(folder))) {
        log_step(paste0("Folder not found: ", folder), "error")
        return()
      }

      log_step(
        paste0("Scanning ", folder, " for MaxQuant files…"),
        "info",
        notify = FALSE
      )
      files <- find_mq_files(folder)
      mq_files(files)

      found <- names(files)[!is.na(unlist(files))]
      missing <- setdiff(names(files), found)
      log_step(
        sprintf(
          "Found %d/%d files%s",
          length(found),
          length(files),
          if (length(missing)) {
            paste0(" (missing: ", paste(missing, collapse = ", "), ")")
          } else {
            ""
          }
        ),
        if (length(missing)) "warn" else "ok"
      )
      if (length(found) == 0) {
        return()
      }

      withProgress(message = "Loading MaxQuant files", value = 0, {
        n_steps <- 5
        # Small files first so the Run Summary tab populates quickly
        raw_summary_rv(load_one(
          "summary.txt",
          files$summary,
          read_summary_file,
          1,
          n_steps
        ))
        protein_groups_rv(load_one(
          "proteinGroups.txt",
          files$proteinGroups,
          read_protein_groups,
          2,
          n_steps
        ))
        peptides_rv(load_one(
          "peptides.txt",
          files$peptides,
          read_peptides_file,
          3,
          n_steps
        ))
        # Offer only the quantity types actually present in this run
        if (!is.null(protein_groups_rv())) {
          avail <- available_quant_types(protein_groups_rv())
          if (length(avail) == 0) {
            log_step(
              "proteinGroups.txt has no per-sample Intensity / LFQ / iBAQ columns.",
              "warn"
            )
          } else {
            log_step(
              paste0(
                "Abundance values available: ",
                paste(names(avail), collapse = ", ")
              ),
              "info",
              notify = FALSE
            )
          }
          updateSelectInput(
            session,
            "pg_quant",
            choices = avail,
            selected = if ("Intensity" %in% avail) "Intensity" else avail[1]
          )
        }
        # Heavy tables go through the parquet cache; conversion is logged
        cache_log <- function(msg, level) log_step(msg, level, notify = FALSE)
        raw_evidence_rv(load_one(
          "evidence.txt",
          files$evidence,
          function(p) read_evidence_file(p, log = cache_log),
          4,
          n_steps
        ))
        raw_scans_rv(load_one(
          "msmsScans.txt",
          files$msmsScans,
          function(p) read_msmsscans_file(p, log = cache_log),
          5,
          n_steps
        ))
        setProgress(1, message = "Loading complete", detail = "")
      })

      log_step("All available files processed. Ready for analysis.", "ok")
    })

    # Sidebar status block
    output$file_status <- renderUI({
      files <- mq_files()
      if (is.null(files)) {
        return(tags$p(
          style = "color:#adb5bd;font-size:11px;",
          icon("circle-info"),
          " No folder loaded yet."
        ))
      }
      loaded <- list(
        summary = raw_summary_rv(),
        proteinGroups = protein_groups_rv(),
        peptides = peptides_rv(),
        evidence = raw_evidence_rv(),
        msmsScans = raw_scans_rv()
      )
      row <- function(nm) {
        if (is.na(files[[nm]])) {
          tags$li(
            style = "color:#e74c3c;",
            icon("xmark"),
            " ",
            nm,
            if (nm %in% MQ_PARQUET_TABLES) {
              ".parquet / .txt missing"
            } else {
              ".txt missing"
            }
          )
        } else if (is.null(loaded[[nm]])) {
          tags$li(
            style = "color:#f39c12;",
            icon("spinner", class = "fa-spin"),
            " ",
            basename(files[[nm]]),
            " found — loading…"
          )
        } else {
          # Show the file actually used (parquet cache when available)
          shown <- files[[nm]]
          if (
            nm %in% MQ_PARQUET_TABLES && file.exists(mq_parquet_path(shown))
          ) {
            shown <- mq_parquet_path(shown)
          }
          detail <- if (is.character(loaded[[nm]])) {
            " (queried on demand)"
          } else {
            paste0(" (", format(nrow(loaded[[nm]]), big.mark = ","), " rows)")
          }
          tags$li(
            style = "color:#2ecc71;",
            icon("check"),
            " ",
            basename(shown),
            detail
          )
        }
      }
      tags$ul(
        style = "list-style:none;padding-left:4px;font-size:11px;margin-top:4px;",
        lapply(
          c("summary", "proteinGroups", "peptides", "evidence", "msmsScans"),
          row
        )
      )
    })

    # Full step log in the Run Summary tab
    output$status_log_ui <- renderUI({
      lg <- status_log()
      if (nrow(lg) == 0) {
        return(tags$p(
          style = "color:#adb5bd;",
          icon("circle-info"),
          " Enter the MaxQuant 'combined' folder path in the sidebar and click ",
          tags$b("Load MaxQuant Files"),
          ". Progress will be reported here."
        ))
      }
      colors <- c(
        info = "#3498db",
        ok = "#2ecc71",
        warn = "#f39c12",
        error = "#e74c3c"
      )
      icons <- c(
        info = "circle-info",
        ok = "check",
        warn = "triangle-exclamation",
        error = "xmark"
      )
      tags$ul(
        style = "list-style:none;padding-left:0;font-family:monospace;font-size:12px;margin:0;",
        lapply(rev(seq_len(nrow(lg))), function(i) {
          tags$li(
            style = paste0("color:", colors[lg$level[i]], ";"),
            tags$span(style = "color:#7f8c8d;", lg$time[i]),
            "  ",
            icon(icons[lg$level[i]]),
            " ",
            lg$msg[i]
          )
        })
      )
    })

    # ════════════════════════════════════════════════════════════════════
    # B. evidence.txt REACTIVES
    # ════════════════════════════════════════════════════════════════════

    raw_evidence <- reactive({
      req(raw_evidence_rv())
    })

    # The currently selected evidence plot — recalculate on button click
    current_ev_plot <- eventReactive(input$run_evidence, {
      if (is.null(raw_evidence_rv())) {
        log_step(
          "evidence.txt is not loaded — load the MaxQuant folder first.",
          "warn"
        )
      }
      req(raw_evidence())
      show_spinner("sp_evidence")

      df <- raw_evidence()
      sel <- input$ev_plot_select
      log_step(
        paste0("Building evidence QC plot: ", sel, "…"),
        "info",
        notify = FALSE
      )

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
      log_step("Evidence QC plot built — rendering…", "ok", notify = FALSE)
      p
    })

    # Dynamic height for evidence plots
    ev_plot_height_px <- reactive({
      req(raw_evidence())
      facet_height_px(raw_evidence())
    })

    # ════════════════════════════════════════════════════════════════════
    # B2. peptides.txt REACTIVES
    # ════════════════════════════════════════════════════════════════════

    peptides <- reactive({
      req(peptides_rv())
    })

    current_pep_plot <- eventReactive(input$run_peptides, {
      if (is.null(peptides_rv())) {
        log_step(
          "peptides.txt is not loaded — load the MaxQuant folder first.",
          "warn"
        )
      }
      req(peptides())
      show_spinner("sp_peptides")

      sel <- input$pep_plot_select
      log_step(
        paste0("Building peptides QC plot: ", sel, "…"),
        "info",
        notify = FALSE
      )
      p <- tryCatch(
        PEP_PLOTS[[sel]][[2]](peptides()),
        error = function(e) {
          log_step(
            sprintf("Peptides plot '%s' failed: %s", sel, conditionMessage(e)),
            "error"
          )
          hide_spinner("sp_peptides")
          NULL
        }
      )
      req(p)
      log_step("Peptides QC plot built — rendering…", "ok", notify = FALSE)
      p
    })

    pep_plot_height_px <- reactive({
      PEP_PLOTS[[input$pep_plot_select]][[3]]
    })

    # ════════════════════════════════════════════════════════════════════
    # B3. msmsScans.txt REACTIVES
    # ════════════════════════════════════════════════════════════════════

    raw_scans <- reactive({
      req(raw_scans_rv())
    })

    current_scan_plot <- eventReactive(input$run_scans, {
      if (is.null(raw_scans_rv())) {
        log_step(
          "msmsScans.txt is not loaded — load the MaxQuant folder first.",
          "warn"
        )
      }
      req(raw_scans())
      show_spinner("sp_scans")

      sel <- input$scan_plot_select
      log_step(
        paste0("Building MS/MS scans QC plot: ", sel, "…"),
        "info",
        notify = FALSE
      )
      p <- tryCatch(
        SCAN_PLOTS[[sel]][[2]](raw_scans()),
        error = function(e) {
          log_step(
            sprintf("Scans plot '%s' failed: %s", sel, conditionMessage(e)),
            "error"
          )
          hide_spinner("sp_scans")
          NULL
        }
      )
      req(p)
      log_step("MS/MS scans QC plot built — rendering…", "ok", notify = FALSE)
      p
    })

    scan_plot_height_px <- reactive({
      req(raw_scans())
      if (identical(input$scan_plot_select, "id_rate")) {
        max(300L, 60L * dplyr::n_distinct(raw_scans()$`Raw file`) + 150L)
      } else {
        facet_height_px(raw_scans())
      }
    })

    # ════════════════════════════════════════════════════════════════════
    # C0. OUTPUTS — Run Summary / Protein matrix
    # ════════════════════════════════════════════════════════════════════

    output$run_summary_table <- DT::renderDataTable({
      req(raw_summary_rv())
      DT::datatable(
        raw_summary_rv(),
        rownames = FALSE,
        extensions = "FixedColumns",
        options = list(
          dom = "frtip",
          pageLength = 25,
          scrollX = TRUE,
          fixedColumns = list(leftColumns = 1)
        ),
        class = "display compact"
      )
    })

    # Abundance matrix for the currently selected quantity type
    protein_matrix <- reactive({
      req(protein_groups_rv(), input$pg_quant)
      req(input$pg_quant %in% available_quant_types(protein_groups_rv()))
      build_protein_matrix(
        protein_groups_rv(),
        input$pg_quant,
        zero_to_na = isTRUE(input$pg_zero_na),
        log2_transform = isTRUE(input$pg_log2)
      )
    })

    # Human-readable description of the current matrix settings
    protein_matrix_label <- reactive({
      paste0(
        input$pg_quant,
        if (isTRUE(input$pg_log2)) {
          " (log2, zeros -> NA)"
        } else if (isTRUE(input$pg_zero_na)) {
          " (zeros -> NA)"
        } else {
          ""
        }
      )
    })

    # Keep the two checkboxes consistent: log2 implies zeros -> NA
    observeEvent(
      input$pg_log2,
      {
        if (isTRUE(input$pg_log2)) {
          updateCheckboxInput(session, "pg_zero_na", value = TRUE)
          shinyjs::disable("pg_zero_na")
        } else {
          shinyjs::enable("pg_zero_na")
        }
      },
      ignoreInit = TRUE
    )

    output$protein_matrix_table <- DT::renderDataTable({
      pm <- protein_matrix()
      DT::datatable(
        head(pm, 500),
        rownames = FALSE,
        caption = paste0(
          "Values: ",
          protein_matrix_label(),
          " — ",
          format(nrow(pm), big.mark = ","),
          " protein groups × ",
          ncol(pm) - 1,
          " samples (first 500 rows shown)"
        ),
        options = list(dom = "frtip", pageLength = 20, scrollX = TRUE),
        class = "display compact"
      ) |>
        DT::formatSignif(
          columns = setdiff(names(pm), "Protein ID"),
          digits = 4
        )
    })

    # Peptide-level summary from peptides.txt (replaces the former msms.txt
    # sequence count, which required the whole msms table in memory)
    output$summary_table <- DT::renderDataTable({
      req(peptides())
      tbl <- pep_summary_table(peptides())
      DT::datatable(
        tbl,
        rownames = FALSE,
        filter = "top",
        options = list(dom = "frtip", pageLength = 25, scrollX = TRUE),
        class = "display compact"
      ) |>
        DT::formatSignif(
          columns = intersect(c("Score", "PEP", "Intensity"), names(tbl)),
          digits = 4
        )
    })

    # ════════════════════════════════════════════════════════════════════
    # D2. OUTPUTS — Peptides QC
    # ════════════════════════════════════════════════════════════════════

    output$pep_plot_title <- renderUI({
      PEP_PLOTS[[input$pep_plot_select]][[1]]
    })

    output$peptides_plot_ui <- renderUI({
      plotOutput(
        ns("peptides_plot"),
        height = paste0(pep_plot_height_px(), "px")
      )
    })

    output$peptides_plot <- renderPlot({
      on.exit(hide_spinner("sp_peptides"), add = TRUE)
      req(current_pep_plot())
      current_pep_plot()
    })

    # ════════════════════════════════════════════════════════════════════
    # D3. OUTPUTS — MS/MS Scans QC
    # ════════════════════════════════════════════════════════════════════

    output$scan_plot_title <- renderUI({
      SCAN_PLOTS[[input$scan_plot_select]][[1]]
    })

    output$scans_plot_ui <- renderUI({
      plotOutput(
        ns("scans_plot"),
        height = paste0(scan_plot_height_px(), "px")
      )
    })

    output$scans_plot <- renderPlot({
      on.exit(hide_spinner("sp_scans"), add = TRUE)
      req(current_scan_plot())
      current_scan_plot()
    })

    output$scans_table <- DT::renderDataTable({
      req(raw_scans())
      DT::datatable(
        head(raw_scans(), 500),
        rownames = FALSE,
        options = list(dom = "frtip", pageLength = 20, scrollX = TRUE),
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
      content = function(file) {
        data.table::fwrite(raw_evidence(), file = file, sep = ",", na = "NA")
      }
    )

    # Peptides plot PDF
    output$download_pep_plot <- downloadHandler(
      filename = function() {
        paste0("peptides_", input$pep_plot_select, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggplot2::ggsave(
          file,
          plot = current_pep_plot(),
          device = "pdf",
          width = 10,
          # ~90 px per inch keeps the PDF proportions close to the screen
          height = max(6, pep_plot_height_px() / 90),
          units = "in"
        )
      }
    )

    # MS/MS scans plot PDF
    output$download_scan_plot <- downloadHandler(
      filename = function() {
        paste0("msmsScans_", input$scan_plot_select, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        n_rows <- ceiling(dplyr::n_distinct(raw_scans()$`Raw file`) / 3)
        ggplot2::ggsave(
          file,
          plot = current_scan_plot(),
          device = "pdf",
          width = 16,
          height = max(5, n_rows * 4),
          units = "in"
        )
      }
    )

    # MS/MS scans TSV
    output$download_scan_data <- downloadHandler(
      filename = function() {
        paste0("msmsScans_data_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        data.table::fwrite(raw_scans(), file = file, sep = "\t", na = "NA")
      }
    )

    # Peptides table TSV (filtered peptides.txt)
    output$download_pep_data <- downloadHandler(
      filename = function() {
        paste0("peptides_data_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        data.table::fwrite(peptides(), file = file, sep = "\t", na = "NA")
      }
    )

    # Protein abundance matrix TSV
    output$download_protein_matrix <- downloadHandler(
      filename = function() {
        tag <- gsub("[^A-Za-z0-9]+", "_", input$pg_quant)
        if (isTRUE(input$pg_log2)) {
          tag <- paste0(tag, "_log2")
        } else if (isTRUE(input$pg_zero_na)) {
          tag <- paste0(tag, "_zeroNA")
        }
        paste0("protein_abundance_matrix_", tag, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        pm <- protein_matrix()
        data.table::fwrite(pm, file = file, sep = "\t", na = "NA")
        log_step(
          sprintf(
            "Protein abundance matrix (%s) downloaded: %s proteins x %d samples.",
            protein_matrix_label(),
            format(nrow(pm), big.mark = ","),
            ncol(pm) - 1
          ),
          "ok"
        )
      }
    )

    # ════════════════════════════════════════════════════════════════════════
    # MODIFICATION DIAGNOSTIC
    # ════════════════════════════════════════════════════════════════════════
    mod_diag_data <- reactive({
      d <- raw_evidence()
      req(d, nrow(d) > 0)
      req(all(
        c("Sequence", "Modifications", "Retention time", "Raw file") %in%
          names(d)
      ))

      d_mod <- d |>
        dplyr::filter(
          !is.na(Modifications),
          !Modifications %in% c("", "Unmodified")
        ) |>
        dplyr::select(
          `Sequence`,
          `Modifications`,
          `Retention time`,
          `Raw file`
        ) |>
        dplyr::rename(
          peptide = `Sequence`,
          mod_label = `Modifications`,
          rt_mod = `Retention time`,
          sample_name = `Raw file`
        )

      d_unmod <- d |>
        dplyr::filter(
          is.na(Modifications) | Modifications %in% c("", "Unmodified")
        ) |>
        dplyr::group_by(`Sequence`, `Raw file`) |>
        dplyr::summarise(
          rt_unmod = median(`Retention time`, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(peptide = `Sequence`, sample_name = `Raw file`)

      dplyr::inner_join(d_mod, d_unmod, by = c("peptide", "sample_name")) |>
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
          title = "RT shift profile \u2014 modified vs unmodified",
          x = "Retention time (min)",
          y = "Peptide sequence",
          caption = paste0(
            "\u0394RT threshold: ",
            input$rt_tolerance,
            " min | Artifact if |\u0394RT| \u2264 threshold"
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
      shinyjs::hide("sp_moddiag")
      if (is.null(paired) || nrow(paired) == 0) {
        return(tags$p(
          style = "color:#adb5bd;text-align:center;padding:20px;",
          "No paired modified/unmodified peptides found. Load the MaxQuant folder (evidence.txt) and ensure it contains 'Sequence', 'Modifications', and 'Retention time' columns."
        ))
      }
      n_samples <- dplyr::n_distinct(paired$sample_name)
      dynamic_h <- max(400L, min(n_samples * 500L, 2000L))
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
          mod_label,
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
  })
}
