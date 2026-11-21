## ============================================================
## utils_fasta.r  —  FASTA parsing and proteotypic peptides
## ============================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Reading FASTA purely in R to avoid heavy Bioconductor dependencies
read_fasta_custom <- function(file) {
  lines <- readLines(file)
  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0) {
    return(character(0))
  }

  headers <- sub("^>", "", lines[header_idx])
  # Entry ID: the accession of UniProt-style headers ("sp|P02769|ALBU_BOVIN",
  # also with a prefix such as "contam_sp|..."), otherwise the first word up
  # to any '|'. The database tag has to be the whole first field: testing for
  # "sp"/"tr" anywhere in it reads ">trypsin_pig Porcine trypsin" as
  # UniProt-style and names the entry "Porcine".
  first_word <- sub("\\s.*$", "", trimws(headers))
  prot_ids <- vapply(
    strsplit(first_word, "|", fixed = TRUE),
    function(x) {
      if (length(x) >= 2 && grepl("(^|_)(sp|tr)$", x[1])) x[2] else x[1]
    },
    character(1)
  )

  seqs <- character(length(header_idx))
  starts <- header_idx + 1
  ends <- c(header_idx[-1] - 1, length(lines))

  for (i in seq_along(header_idx)) {
    if (starts[i] <= ends[i]) {
      seqs[i] <- paste(lines[starts[i]:ends[i]], collapse = "")
    }
  }
  names(seqs) <- prot_ids
  return(seqs)
}

# In silico digestion
in_silico_digest <- function(fasta_seqs, max_missed = 2) {
  pep_list <- lapply(fasta_seqs, function(seq) {
    if (is.na(seq) || nchar(seq) == 0) {
      return(character(0))
    }
    cleaved <- strsplit(seq, "(?<=[KR])(?!P)", perl = TRUE)[[1]]

    peps <- cleaved
    if (max_missed > 0 && length(cleaved) > 1) {
      for (m in 1:max_missed) {
        if (length(cleaved) > m) {
          for (i in 1:(length(cleaved) - m)) {
            peps <- c(peps, paste(cleaved[i:(i + m)], collapse = ""))
          }
        }
      }
    }
    # Standard peptide lengths for typical proteomics
    peps <- peps[nchar(peps) >= 6 & nchar(peps) <= 50]
    return(peps)
  })

  len <- sapply(pep_list, length)
  prot_ids <- names(fasta_seqs)

  df <- data.frame(
    protein = rep(prot_ids, len),
    peptide = unlist(pep_list),
    stringsAsFactors = FALSE
  )

  df |> distinct(peptide, protein)
}

# Classifying observed peptides against the defined FASTA search space
classify_peptides <- function(detected_peptides, digest_df) {
  mapping <- digest_df |>
    filter(peptide %in% detected_peptides) |>
    group_by(peptide) |>
    summarise(
      mapped_proteins = paste(unique(protein), collapse = ";"),
      n_proteins = n_distinct(protein),
      .groups = "drop"
    )

  out <- data.frame(
    peptide = unique(detected_peptides),
    stringsAsFactors = FALSE
  )
  out <- left_join(out, mapping, by = "peptide")
  out$classification <- ifelse(
    is.na(out$n_proteins) | out$n_proteins == 0,
    "Unmapped",
    ifelse(out$n_proteins == 1, "Proteotypic", "Shared")
  )
  out
}

# Constructing the PwrQuant-compatible matrix based strictly on proteotypic sum/top3
compute_protein_abundance <- function(
  df,
  pep_col,
  abundance_col,
  class_col,
  mapped_prot_col,
  method = "sum"
) {
  req_cols <- c("filename", pep_col, abundance_col, class_col, mapped_prot_col)

  # Ensure all necessary columns exist (filename defaults to unknown if missing)
  if (!"filename" %in% names(df)) {
    df$filename <- "Unknown"
  }
  missing_cols <- setdiff(req_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }

  df_sub <- df |>
    filter(!!sym(class_col) == "Proteotypic") |>
    filter(!is.na(!!sym(abundance_col)))

  if (nrow(df_sub) == 0) {
    return(NULL)
  }

  # Top3 or Sum calculation
  if (method == "sum") {
    mat <- df_sub |>
      group_by(!!sym(mapped_prot_col), filename) |>
      summarise(
        Abundance = sum(!!sym(abundance_col), na.rm = TRUE),
        .groups = "drop"
      ) |>
      pivot_wider(
        names_from = filename,
        values_from = Abundance,
        values_fill = NA
      )
  } else if (method == "top3") {
    mat <- df_sub |>
      group_by(!!sym(mapped_prot_col), filename) |>
      slice_max(order_by = !!sym(abundance_col), n = 3, with_ties = FALSE) |>
      summarise(
        Abundance = sum(!!sym(abundance_col), na.rm = TRUE),
        .groups = "drop"
      ) |>
      pivot_wider(
        names_from = filename,
        values_from = Abundance,
        values_fill = NA
      )
  }

  return(mat)
}

# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL FASTA MODULE
# ═══════════════════════════════════════════════════════════════════════════════

Fasta_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$hr(style = "border-color:#2d3741;margin:6px 0;"),
    tags$div(
      style = "padding:4px 16px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
      icon("dna"),
      " Global FASTA Mapping"
    ),
    fileInput(
      ns("fasta_file"),
      "Protein FASTA File",
      accept = c(".fasta", ".fa")
    ),
    numericInput(
      ns("missed_cleavages"),
      "Max missed cleavages",
      value = 2,
      min = 0,
      max = 5
    )
  )
}

Fasta_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactive({
      if (!isTruthy(input$fasta_file)) {
        return(NULL)
      }

      showNotification(
        "Reading and digesting FASTA globally...",
        id = "fasta_global_notif",
        duration = NULL
      )
      seqs <- tryCatch(
        read_fasta_custom(input$fasta_file$datapath),
        error = function(e) {
          showNotification("Failed to read FASTA", type = "error")
          NULL
        }
      )
      if (is.null(seqs)) {
        removeNotification("fasta_global_notif")
        return(NULL)
      }

      df <- in_silico_digest(seqs, max_missed = input$missed_cleavages)
      removeNotification("fasta_global_notif")
      df
    })
  })
}

# ── Global Peptide Property Utilities ──────────────────────────────────────────
GRAVY <- function(sequence) {
  hi <- c(
    A = 1.8,
    R = -4.5,
    N = -3.5,
    D = -3.5,
    C = 2.5,
    Q = -3.5,
    E = -3.5,
    G = -0.4,
    H = -3.2,
    I = 4.5,
    L = 3.8,
    K = -3.9,
    M = 1.9,
    F = 2.8,
    P = -1.6,
    S = -0.8,
    T = -0.7,
    W = -0.9,
    Y = -1.3,
    V = 4.2
  )
  mean(sapply(strsplit(sequence, NULL)[[1]], function(aa) hi[aa]), na.rm = TRUE)
}

# Vectorized isoelectric point
calculate_pI <- function(sequence) {
  seqs <- toupper(as.character(sequence))
  n <- length(seqs)
  if (n == 0) {
    return(numeric(0))
  }

  pk_nterm <- 9.69
  pk_cterm <- 2.34
  pk_acid <- c(D = 3.86, E = 4.25, C = 8.33, Y = 10.07) # HA -> A- + H+
  pk_base <- c(H = 6.00, K = 10.53, R = 12.48) # BH+ -> B + H+

  count_res <- function(letter) {
    nchar(gsub(paste0("[^", letter, "]"), "", seqs))
  }

  acid_counts <- vapply(names(pk_acid), count_res, numeric(n))
  base_counts <- vapply(names(pk_base), count_res, numeric(n))
  if (n == 1) {
    dim(acid_counts) <- c(1, length(pk_acid))
    dim(base_counts) <- c(1, length(pk_base))
  }

  n_standard <- nchar(gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", seqs))
  valid <- n_standard > 0 & !is.na(seqs)
  net_charge <- function(pH) {
    pos <- 1 /
      (1 + 10^(pH - pk_nterm)) +
      as.vector(base_counts %*% (1 / (1 + 10^(pH - pk_base))))
    neg <- 1 /
      (1 + 10^(pk_cterm - pH)) +
      as.vector(acid_counts %*% (1 / (1 + 10^(pk_acid - pH))))
    pos - neg
  }

  grid <- seq(0, 14, by = 0.1)
  charges <- vapply(grid, net_charge, numeric(n))
  if (n == 1) {
    dim(charges) <- c(1, length(grid))
  }

  sign_change <- charges[, -1, drop = FALSE] *
    charges[, -ncol(charges), drop = FALSE] <=
    0
  first_cross <- apply(sign_change, 1, function(x) which(x)[1])
  pI <- numeric(n)
  no_cross <- is.na(first_cross)
  if (any(no_cross)) {
    pI[no_cross] <- grid[apply(
      abs(charges[no_cross, , drop = FALSE]),
      1,
      which.min
    )]
  }

  has_cross <- !no_cross
  if (any(has_cross)) {
    lo <- grid[first_cross[has_cross]]
    hi <- grid[first_cross[has_cross] + 1]
    ac <- acid_counts[has_cross, , drop = FALSE]
    bc <- base_counts[has_cross, , drop = FALSE]
    nc_sub <- function(pH) {
      pos <- 1 /
        (1 + 10^(pH - pk_nterm)) +
        rowSums(bc * (1 / (1 + 10^(outer(pH, pk_base, "-")))))
      neg <- 1 /
        (1 + 10^(pk_cterm - pH)) +
        rowSums(ac * (1 / (1 + 10^(-outer(pH, pk_acid, "-")))))
      pos - neg
    }
    for (i in seq_len(60)) {
      mid <- (lo + hi) / 2
      c_mid <- nc_sub(mid)
      lo <- ifelse(c_mid > 0, mid, lo)
      hi <- ifelse(c_mid > 0, hi, mid)
    }
    pI[has_cross] <- (lo + hi) / 2
  }

  pI[!valid] <- NA_real_
  pI
}

# Vectorized average molecular weight (Da) from amino acid composition.
calculate_MW <- function(sequence) {
  seqs <- toupper(as.character(sequence))
  n <- length(seqs)
  if (n == 0) {
    return(numeric(0))
  }
  res_mass <- c(
    A = 71.0788,
    R = 156.1875,
    N = 114.1038,
    D = 115.0886,
    C = 103.1388,
    E = 129.1155,
    Q = 128.1307,
    G = 57.0519,
    H = 137.1411,
    I = 113.1594,
    L = 113.1594,
    K = 128.1741,
    M = 131.1926,
    F = 147.1766,
    P = 97.1167,
    S = 87.0782,
    T = 101.1051,
    W = 186.2132,
    Y = 163.1760,
    V = 99.1326
  )
  water <- 18.01528
  total <- rep(0, n)
  has_res <- rep(FALSE, n)
  for (aa in names(res_mass)) {
    k <- nchar(gsub(paste0("[^", aa, "]"), "", seqs))
    total <- total + k * res_mass[[aa]]
    has_res <- has_res | k > 0
  }
  mw <- total + water
  mw[!has_res | is.na(seqs)] <- NA_real_
  mw
}

point_density_2d <- function(x, y, n = 100) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & is.finite(y)
  if (
    sum(ok) < 3 ||
      length(unique(x[ok])) < 2 ||
      length(unique(y[ok])) < 2 ||
      !requireNamespace("MASS", quietly = TRUE)
  ) {
    out[ok] <- 1
    return(out)
  }
  dens <- MASS::kde2d(x[ok], y[ok], n = n)
  ix <- findInterval(x[ok], dens$x, all.inside = TRUE)
  iy <- findInterval(y[ok], dens$y, all.inside = TRUE)
  out[ok] <- dens$z[cbind(ix, iy)]
  out
}

plot_2d_gel <- function(
  d,
  pi_col = "pI",
  mw_col = "MW",
  title = "Virtual 2D gel (pI vs MW)",
  facet_col = NULL,
  ncol = 3
) {
  d <- d[is.finite(d[[pi_col]]) & is.finite(d[[mw_col]]), , drop = FALSE]
  if (nrow(d) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No pI / MW data available"
        ) +
        ggplot2::theme_void()
    )
  }
  faceted <- !is.null(facet_col) && facet_col %in% names(d)
  if (faceted) {
    d$.density <- stats::ave(
      seq_len(nrow(d)),
      d[[facet_col]],
      FUN = function(i) point_density_2d(d[[pi_col]][i], d[[mw_col]][i])
    )
  } else {
    d$.density <- point_density_2d(d[[pi_col]], d[[mw_col]])
  }
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = .data[[pi_col]],
      y = .data[[mw_col]],
      color = .data[[".density"]]
    )
  ) +
    ggplot2::geom_point(size = 1.5, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(name = "Spot density", option = "D") +
    ggplot2::labs(
      title = title,
      x = "Isoelectric point (pI)",
      y = "Molecular weight (Da)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.title.position = "top",
      legend.text = element_text(size = 9, face = "bold"),
      legend.key.height = unit(0.25, "cm"),
      legend.key.width = unit(1, "cm"),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA)
    )
  if (faceted) {
    p <- p + ggplot2::facet_wrap(stats::reformulate(facet_col), ncol = ncol)
  }
  p
}

# ── Peptide-to-protein coverage ───────────────────────────────────────────────
# Shared by mod_QC4DIANN (section E) and mod_PSManalyst (sequence coverage).

# Normalize a FASTA sequence for matching: one string, no internal whitespace,
# upper case. Some tools emit lower-case FASTA files, which would otherwise
# match zero peptides and report 0% coverage with no error.
fasta_seq_normalize <- function(seq) {
  toupper(gsub("[[:space:]]", "", as.character(seq)[1]))
}

# The '|'-delimited tokens of a FASTA header's identifier field, e.g.
# "sp|P02769|ALBU_BOVIN Albumin OS=..." -> c("sp", "P02769", "ALBU_BOVIN").
fasta_header_tokens <- function(headers) {
  lapply(headers, function(h) {
    strsplit(trimws(strsplit(h, "\\s+")[[1]][1]), "|", fixed = TRUE)[[1]]
  })
}

# Index of the FASTA entry whose accession or entry name equals `target`.
# A substring grep selects the wrong entry (searching "P0231" matches the
# header of "P02310"), so require a whole-token match.
fasta_index_exact <- function(target, headers) {
  toks <- fasta_header_tokens(headers)
  hit <- which(vapply(toks, function(t) target %in% t, logical(1)))
  if (length(hit) == 0) {
    hit <- which(vapply(
      toks,
      function(t) any(tolower(t) == tolower(target)),
      logical(1)
    ))
  }
  if (length(hit) == 0) NA_integer_ else hit[1]
}

# Whether a search-engine protein field ("sp|P1|A_HUMAN;sp|P2|B_HUMAN")
# contains `target` as a complete accession or entry name rather than as a
# substring of some other identifier.
protein_ids_contain <- function(ids, target) {
  vapply(
    ids,
    function(x) {
      if (is.na(x)) {
        return(FALSE)
      }
      toks <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
      toks <- unlist(strsplit(trimws(toks), "|", fixed = TRUE))
      target %in% toks
    },
    logical(1),
    USE.NAMES = FALSE
  )
}

# 1-based start positions of every non-overlapping occurrence of `pep`, so
# peptides in repeat regions are covered wherever they occur, not only once.
peptide_starts <- function(pep, prot_seq) {
  if (length(pep) != 1 || is.na(pep) || nchar(pep) == 0) {
    return(integer(0))
  }
  hits <- as.integer(gregexpr(
    toupper(pep),
    toupper(prot_seq),
    fixed = TRUE
  )[[1]])
  hits[hits > 0L]
}

# Per-residue coverage depth and mask for a set of peptides against one
# protein sequence. `counts` weights each peptide (e.g. PSM counts); NULL
# weights every peptide equally. Returns the number of peptides that did not
# map anywhere in the sequence.
protein_coverage <- function(prot_seq, peptides, counts = NULL) {
  pl <- nchar(prot_seq)
  depth <- numeric(pl)
  mask <- rep(FALSE, pl)
  if (is.null(counts)) {
    counts <- rep(1, length(peptides))
  }
  counts <- ifelse(is.na(counts), 1, counts)
  n_unmapped <- 0L
  for (i in seq_along(peptides)) {
    starts <- peptide_starts(peptides[i], prot_seq)
    if (length(starts) == 0) {
      n_unmapped <- n_unmapped + 1L
      next
    }
    for (s in starts) {
      e <- s + nchar(peptides[i]) - 1L
      if (e > pl) {
        next
      }
      depth[s:e] <- depth[s:e] + counts[i]
      mask[s:e] <- TRUE
    }
  }
  list(depth = depth, mask = mask, n_unmapped = n_unmapped)
}

# Function to extract the protein abundace matrix from the output of the DIA-NN pipeline

extract_diann_matrix <- function(
  x,
  id.header = "Genes",
  quantity.header = "Genes.MaxLFQ.Unique",
  proteotypic.only = FALSE,
  q = 0.01,
  protein.q = 0.01,
  pg.q = 0.01,
  gg.q = 0.01
) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)

  # 1. Filter by proteotypic status if requested
  if (proteotypic.only) {
    df <- df[which(df[["Proteotypic"]] != 0), , drop = FALSE]
  }

  # 2. Filter by validity, positive quantity, and Q-value thresholds
  keep <- which(
    !is.na(df[[id.header]]) &
      df[[id.header]] != "" &
      !is.na(df[[quantity.header]]) &
      df[[quantity.header]] > 0 &
      df[["Q.Value"]] <= q &
      df[["Protein.Q.Value"]] <= protein.q &
      df[["PG.Q.Value"]] <= pg.q &
      df[["GG.Q.Value"]] <= gg.q
  )

  df <- unique(df[keep, c("Run", id.header, quantity.header), drop = FALSE])

  if (nrow(df) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  # 3. Check for duplicate (Run + ID) combinations
  is_duplicated <- any(duplicated(paste0(df[["Run"]], ":", df[[id.header]])))
  ids <- sort(unique(df[[id.header]]))
  files <- sort(unique(df[["Run"]]))

  # 4. Pivot to wide matrix (Rows = IDs, Columns = Run)
  if (is_duplicated) {
    warning(
      "Multiple quantities per id: the maximum of these will be calculated"
    )

    row_fac <- factor(df[[id.header]], levels = ids)
    col_fac <- factor(df[["Run"]], levels = files)
    mat <- tapply(df[[quantity.header]], list(row_fac, col_fac), FUN = max)
    attr(mat, "dimnames") <- list(ids, files)
  } else {
    mat <- matrix(
      NA_real_,
      nrow = length(ids),
      ncol = length(files),
      dimnames = list(ids, files)
    )
    row_idx <- match(df[[id.header]], ids)
    col_idx <- match(df[["Run"]], files)

    mat[cbind(row_idx, col_idx)] <- df[[quantity.header]]
  }

  return(mat)
}

# ── iBAQ from a DIA-NN report ─────────────────────────────────────────────────
# Ported from the iBAQ_for_diann reference script (iBAQ_from_DIANN.R) so the
# app does not need it as a package dependency.
#   iBAQ = summed precursor intensity of a protein in a run /
#          number of theoretically observable peptides of that protein
# Schwanhäusser et al. (2011) Nature 473:337-342.

# Accession for each FASTA header: the second '|' field of UniProt-style
# headers ("sp|P02769|ALBU_BOVIN ..."), otherwise the first word.
fasta_accessions <- function(headers) {
  first <- sub("\\s.*$", "", trimws(sub("^>", "", headers)))
  uniprot <- grepl("^[a-z]{2}\\|[^|]+\\|", first)
  first[uniprot] <- sub("^[a-z]{2}\\|([^|]+)\\|.*$", "\\1", first[uniprot])
  first
}

# Named character vector (accession -> sequence) from a FASTA path or from an
# already-parsed FASTA (seqinr::read.fasta() list, read_fasta_custom() vector).
ibaq_fasta_sequences <- function(fasta) {
  if (is.character(fasta) && length(fasta) == 1 && is.null(names(fasta))) {
    fasta <- read_fasta_custom(fasta)
  }
  seqs <- vapply(fasta, fasta_seq_normalize, character(1), USE.NAMES = FALSE)
  names(seqs) <- fasta_accessions(names(fasta))
  seqs[!duplicated(names(seqs))]
}

# Number of in-silico peptides within [min_len, max_len] for each sequence.
# Peptide lengths with missed cleavages are sums of consecutive fragment
# lengths, so nothing has to be pasted together. NA for empty sequences.
count_theoretical_peptides <- function(
  sequences,
  protease_regex = "(?<=[KR])(?!P)",
  min_len = 6,
  max_len = 30,
  max_missed = 0
) {
  fragments <- strsplit(as.character(sequences), protease_regex, perl = TRUE)
  counts <- vapply(
    fragments,
    function(fr) {
      len <- nchar(fr)
      len <- len[!is.na(len) & len > 0]
      n <- length(len)
      if (n == 0) {
        return(NA_integer_)
      }
      ends <- c(0L, cumsum(len))
      total <- 0L
      for (mc in 0:max_missed) {
        if (n - mc < 1) {
          break
        }
        pep_len <- ends[(mc + 2):(n + 1)] - ends[seq_len(n - mc)]
        total <- total + sum(pep_len >= min_len & pep_len <= max_len)
      }
      as.integer(total)
    },
    integer(1)
  )
  names(counts) <- names(sequences)
  counts
}

build_theoretical_peptide_counts <- function(sequences, ...) {
  counts <- count_theoretical_peptides(sequences, ...)
  tibble::tibble(
    protein_id = names(counts),
    n_theoretical_peptides = unname(counts)
  )
}

# ── iBAQ helpers shared by the search-engine modules (FragPipe, Sage) ─────────

# First word of each FASTA header: the protein identifier that FragPipe and
# Sage report (e.g. "sp|P04075|ALDOA_HUMAN", "contam_sp|P02769|ALBU_BOVIN").
fasta_header_ids <- function(headers) {
  sub("\\s.*$", "", trimws(sub("^>", "", headers)))
}

# Gene name from the "GN=" field of UniProt-style headers; NA when absent.
fasta_header_genes <- function(headers) {
  gene <- rep(NA_character_, length(headers))
  has <- grepl("\\bGN=\\S+", headers)
  gene[has] <- sub("^.*\\bGN=(\\S+).*$", "\\1", headers[has])
  gene
}

# Named sequences (header's first word -> sequence) from an AAStringSet.
fasta_sequences_by_header_id <- function(fasta) {
  seqs <- toupper(gsub("[[:space:]]", "", as.character(fasta)))
  names(seqs) <- fasta_header_ids(names(fasta))
  seqs[!duplicated(names(seqs))]
}

# Index into `sequences` for each protein: match on the header's first word,
# then fall back to the accession for a FASTA whose headers differ from the
# search database. NA when the protein is not in the FASTA.
match_fasta_sequences <- function(protein, accession, sequences) {
  idx <- match(protein, names(sequences))
  miss <- is.na(idx)
  if (any(miss)) {
    idx[miss] <- match(accession[miss], fasta_accessions(names(sequences)))
  }
  idx
}

# Cleavage-site regex for count_theoretical_peptides() from an enzyme rule:
# residues that are cut, residues that block the cut, and the side of the cut
# residue that is cleaved ("C" = after it, "N" = before it).
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

# Quality filters for a DIA-NN report. A NULL cutoff (or a column absent from
# the report) skips that filter.
filter_diann_report <- function(
  report,
  q_value_cutoff = 0.01,
  pg_q_value_cutoff = 0.01,
  lib_pg_q_value_cutoff = 0.01,
  proteotypic_only = TRUE,
  lfq_quality_cutoff = 0.50
) {
  keep <- rep(TRUE, nrow(report))
  at_most <- function(col, cutoff) {
    if (!is.null(cutoff) && col %in% names(report)) {
      keep <<- keep & !is.na(report[[col]]) & report[[col]] <= cutoff
    }
  }
  at_most("Q.Value", q_value_cutoff)
  at_most("PG.Q.Value", pg_q_value_cutoff)
  at_most("Lib.PG.Q.Value", lib_pg_q_value_cutoff)
  if (isTRUE(proteotypic_only) && "Proteotypic" %in% names(report)) {
    keep <- keep & !is.na(report$Proteotypic) & report$Proteotypic == 1
  }
  if (
    !is.null(lfq_quality_cutoff) &&
      lfq_quality_cutoff > 0 &&
      "PG.MaxLFQ.Quality" %in% names(report)
  ) {
    keep <- keep &
      !is.na(report$PG.MaxLFQ.Quality) &
      report$PG.MaxLFQ.Quality >= lfq_quality_cutoff
  }
  report[keep, , drop = FALSE]
}

# iBAQ per protein and run.
#   report    : DIA-NN report data frame (e.g. from arrow::read_parquet()).
#   sequences : named sequences from ibaq_fasta_sequences().
# Returns a tibble, wide (one column per run) or long (sample, iBAQ). Protein
# groups whose ID is not a FASTA accession (e.g. "P1;P2" groups, or a FASTA
# that differs from the search database) are dropped; they are listed in the
# "unmatched_proteins" attribute of the result.
compute_ibaq <- function(
  report,
  sequences,
  intensity_col = "Precursor.Normalised",
  q_value_cutoff = 0.01,
  pg_q_value_cutoff = 0.01,
  lib_pg_q_value_cutoff = 0.01,
  proteotypic_only = TRUE,
  lfq_quality_cutoff = 0.50,
  protease_regex = "(?<=[KR])(?!P)",
  min_peptide_len = 6,
  max_peptide_len = 30,
  max_missed_cleavages = 0,
  log2_transform = TRUE,
  output_long = FALSE,
  id_col = "Protein.Group"
) {
  report <- as.data.frame(report, stringsAsFactors = FALSE)
  missing_cols <- setdiff(c(id_col, "Run", intensity_col), names(report))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in the DIA-NN report: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  report <- filter_diann_report(
    report,
    q_value_cutoff = q_value_cutoff,
    pg_q_value_cutoff = pg_q_value_cutoff,
    lib_pg_q_value_cutoff = lib_pg_q_value_cutoff,
    proteotypic_only = proteotypic_only,
    lfq_quality_cutoff = lfq_quality_cutoff
  )
  if (nrow(report) == 0) {
    stop("No rows passed the quality filters; relax the filter thresholds.")
  }

  # id_col may itself be an annotation column (e.g. id_col = "Genes"); it is
  # already reported as protein_id, so it must not be selected twice.
  annot_cols <- setdiff(
    intersect(c("Protein.Names", "Genes"), names(report)),
    id_col
  )
  annot <- report[!duplicated(report[[id_col]]), c(id_col, annot_cols)]
  names(annot)[1] <- "protein_id"

  intensity <- data.frame(
    Run = report$Run,
    protein_id = report[[id_col]],
    intensity = report[[intensity_col]]
  ) |>
    group_by(Run, protein_id) |>
    summarise(sum_intensity = sum(intensity, na.rm = TRUE), .groups = "drop")

  theo <- build_theoretical_peptide_counts(
    sequences,
    protease_regex = protease_regex,
    min_len = min_peptide_len,
    max_len = max_peptide_len,
    max_missed = max_missed_cleavages
  )
  unmatched <- setdiff(unique(intensity$protein_id), theo$protein_id)

  ibaq <- intensity |>
    inner_join(theo, by = "protein_id") |>
    left_join(annot, by = "protein_id") |>
    # Proteins without any observable peptide would divide by zero
    filter(!is.na(n_theoretical_peptides), n_theoretical_peptides > 0) |>
    mutate(iBAQ = sum_intensity / n_theoretical_peptides)
  if (isTRUE(log2_transform)) {
    ibaq$iBAQ <- log2(ibaq$iBAQ + 1)
  }

  id_cols <- c("protein_id", annot_cols, "n_theoretical_peptides")
  out <- ibaq |>
    select(all_of(c(id_cols, "Run", "iBAQ")))
  out <- if (isTRUE(output_long)) {
    out |> rename(sample = Run) |> arrange(protein_id, sample)
  } else {
    out |>
      arrange(Run) |>
      pivot_wider(
        id_cols = all_of(id_cols),
        names_from = Run,
        values_from = iBAQ
      ) |>
      arrange(protein_id)
  }
  attr(out, "unmatched_proteins") <- unmatched
  out
}

# ═══════════════════════════════════════════════════════════════════════════════
# SEQUENCE PROFILING  (backend of modules/mod_SeqProfile.r)
# ═══════════════════════════════════════════════════════════════════════════════
# Sliding-window profiles of a single protein / peptide sequence.
#
# SCALE PROVENANCE: every score in this section is computed by a function of
# the CRAN package Peptides (called as Peptides::fn(), never attached, because
# it exports names such as charge(), pI() and mw()). No amino-acid scale value,
# pKa constant or weighting factor is defined here. The only numbers below are
# user-editable parameter defaults (window size, pH, rotation angle, threshold).
#
# Not implemented, because Peptides does not provide them and their tables
# would have to be entered by hand: Karplus-Schulz flexibility, von Heijne
# transmembrane, Emini surface probability, Chou-Fasman, Jameson-Wolf
# antigenicity, Sette MHC-II. Add one entry to seqprofile_registry() per
# algorithm when a verified source is available.

SEQPROFILE_AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

# Signals an error the Shiny module can show to the user as a validation
# message instead of a stack trace.
seqprofile_abort <- function(message, class = NULL) {
  stop(structure(
    class = c(class, "seqprofile_error", "error", "condition"),
    list(message = message, call = NULL)
  ))
}

#' Clean pasted sequence text (raw residues or FASTA).
#'
#' Header lines (">...") are dropped, whitespace and digits removed, and the
#' residues upper-cased. When several FASTA entries are pasted only the first
#' is kept; `n_entries` reports how many were found.
#'
#' @param raw_text Character vector, as typed or pasted by the user.
#' @return list(header, sequence, length, invalid_chars, is_valid, n_entries,
#'   illegal_chars). `invalid_chars` are all characters outside the 20
#'   canonical residues (B, Z, X, U, O, * ...); `illegal_chars` is the subset
#'   that is not even a letter or "*", which makes `is_valid` FALSE.
clean_protein_sequence <- function(raw_text) {
  if (is.null(raw_text) || length(raw_text) == 0) {
    raw_text <- ""
  }
  raw_text <- paste(as.character(raw_text), collapse = "\n")
  lines <- trimws(strsplit(raw_text, "\r\n|\r|\n")[[1]])
  lines <- lines[nzchar(lines)]
  is_header <- startsWith(lines, ">")
  # Entry 0 is any text placed before the first header
  entry <- cumsum(is_header)
  n_entries <- length(unique(entry))

  header <- NA_character_
  body <- character(0)
  if (length(lines) > 0) {
    first <- entry == entry[1]
    if (any(is_header & first)) {
      header <- sub("^>\\s*", "", lines[is_header & first][1])
    }
    body <- lines[!is_header & first]
  }

  sequence <- toupper(gsub(
    "[[:space:][:digit:]]",
    "",
    paste(body, collapse = "")
  ))
  chars <- unique(strsplit(sequence, "")[[1]])
  invalid_chars <- sort(setdiff(chars, SEQPROFILE_AA))
  illegal_chars <- invalid_chars[!grepl("^[A-Z*]$", invalid_chars)]

  list(
    header = header,
    sequence = sequence,
    length = nchar(sequence),
    invalid_chars = invalid_chars,
    is_valid = nchar(sequence) > 0 && length(illegal_chars) == 0,
    n_entries = n_entries,
    illegal_chars = illegal_chars
  )
}

#' Sliding windows along a sequence.
#'
#' The reported position is the window centre, `start + (window - 1) %/% 2`.
#' For an even window this is the residue just left of the midpoint (window 6
#' starting at residue 1 is reported at position 3).
#'
#' @return data.frame(position, residue, window_seq), one row per window.
seqprofile_windows <- function(sequence, window) {
  n <- nchar(sequence)
  if (n < window) {
    seqprofile_abort(
      sprintf(
        "The sequence has %d residue(s), shorter than the window of %d. Use a longer sequence or a smaller window.",
        n,
        window
      ),
      "seqprofile_short_sequence"
    )
  }
  starts <- seq_len(n - window + 1L)
  centre <- starts + (window - 1L) %/% 2L
  data.frame(
    position = as.integer(centre),
    residue = substring(sequence, centre, centre),
    window_seq = substring(sequence, starts, starts + window - 1L),
    stringsAsFactors = FALSE
  )
}

#' Group contiguous positions crossing a threshold into regions.
#'
#' @param profile_df data.frame with `position`, `residue` and `score`, one row
#'   per consecutive position.
#' @param threshold Numeric cutoff; `NA` returns no regions.
#' @param direction "above" keeps `score >= threshold`, "below" keeps
#'   `score <= threshold`. Positions with an `NA` score end a region.
#' @return data.frame(region_id, start, end, length, sequence, peak_pos,
#'   peak_score, mean_score, note). `sequence` is the run of window-centre
#'   residues from `start` to `end`.
extract_profile_regions <- function(
  profile_df,
  threshold,
  direction = "above"
) {
  direction <- match.arg(direction, c("above", "below"))
  empty <- data.frame(
    region_id = integer(0),
    start = integer(0),
    end = integer(0),
    length = integer(0),
    sequence = character(0),
    peak_pos = integer(0),
    peak_score = numeric(0),
    mean_score = numeric(0),
    note = character(0),
    stringsAsFactors = FALSE
  )
  if (nrow(profile_df) == 0 || length(threshold) != 1 || is.na(threshold)) {
    return(empty)
  }
  score <- profile_df$score
  hit <- !is.na(score) &
    (if (direction == "above") score >= threshold else score <= threshold)
  if (!any(hit)) {
    return(empty)
  }
  runs <- rle(hit)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  keep <- which(runs$values)
  pick <- if (direction == "above") which.max else which.min

  rows <- lapply(seq_along(keep), function(k) {
    idx <- starts[keep[k]]:ends[keep[k]]
    peak <- idx[pick(score[idx])]
    data.frame(
      region_id = k,
      start = profile_df$position[idx[1]],
      end = profile_df$position[idx[length(idx)]],
      length = length(idx),
      sequence = paste(profile_df$residue[idx], collapse = ""),
      peak_pos = profile_df$position[peak],
      peak_score = score[peak],
      mean_score = mean(score[idx]),
      note = "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ── Immunological Propensity Scales ───────────────

# Kolaskar & Tongaonkar (1990) FEBS Lett 276(1-2):172-174.
# Source: IEDB B-Cell Analysis Resource (tools.iedb.org/bcell/help/).
SEQPROFILE_SCALE_KOLASKAR <- c(
  A = 1.064,
  C = 1.412,
  D = 0.866,
  E = 0.851,
  F = 1.091,
  G = 0.874,
  H = 1.105,
  I = 1.152,
  K = 0.930,
  L = 1.250,
  M = 0.826,
  N = 0.776,
  P = 1.064,
  Q = 1.015,
  R = 0.873,
  S = 1.012,
  T = 0.909,
  V = 1.383,
  W = 0.893,
  Y = 1.161
)

# Emini et al. (1985) J Virol 55(3):836-839 (fractional surface probability d_n).
# Source: IEDB B-Cell Analysis Resource (tools.iedb.org/bcell/help/).
SEQPROFILE_SCALE_EMINI <- c(
  A = 0.49,
  C = 0.26,
  D = 0.81,
  E = 0.84,
  F = 0.42,
  G = 0.48,
  H = 0.66,
  I = 0.34,
  K = 0.97,
  L = 0.40,
  M = 0.48,
  N = 0.78,
  P = 0.75,
  Q = 0.84,
  R = 0.95,
  S = 0.65,
  T = 0.70,
  V = 0.36,
  W = 0.51,
  Y = 0.76
)

# Parker, Guo & Hodges (1986) Biochemistry 25(19):5425-5432 (HPLC hydrophilicity).
# Source: IEDB B-Cell Analysis Resource (tools.iedb.org/bcell/help/).
SEQPROFILE_SCALE_PARKER <- c(
  A = 2.1,
  C = 1.4,
  D = 10.0,
  E = 7.8,
  F = -9.2,
  G = 5.7,
  H = 2.1,
  I = -8.0,
  K = 5.7,
  L = -9.2,
  M = -4.2,
  N = 7.0,
  P = 2.1,
  Q = 6.0,
  R = 4.2,
  S = 6.5,
  T = 5.2,
  V = -3.7,
  W = -10.0,
  Y = -1.9
)

# Chou & Fasman (1978) Adv Enzymol Relat Areas Mol Biol 47:45-148 (beta-turn P_t).
# Source: IEDB B-Cell Analysis Resource (tools.iedb.org/bcell/help/).
SEQPROFILE_SCALE_CHOU_FASMAN_TURN <- c(
  A = 0.66,
  C = 1.19,
  D = 1.46,
  E = 0.74,
  F = 0.60,
  G = 1.56,
  H = 0.95,
  I = 0.47,
  K = 1.01,
  L = 0.59,
  M = 0.60,
  N = 1.56,
  P = 1.52,
  Q = 0.98,
  R = 0.95,
  S = 1.43,
  T = 0.96,
  V = 0.50,
  W = 0.96,
  Y = 1.14
)

# Calis et al. (2013) PLoS Comput Biol 9(10):e1003266 (Table 2 & IEDB Class I
# Immunogenicity). Logarithmic amino-acid enrichment scores and nonamer
# position weights (Kullback-Leibler divergence; P1, P2, P9 = 0).
SEQPROFILE_SCALE_CALIS <- c(
  A = 0.127,
  C = -0.175,
  D = 0.072,
  E = 0.325,
  F = 0.380,
  G = 0.110,
  H = 0.105,
  I = 0.432,
  K = -0.700,
  L = -0.036,
  M = -0.570,
  N = -0.021,
  P = -0.036,
  Q = -0.376,
  R = 0.168,
  S = -0.537,
  T = 0.126,
  V = 0.134,
  W = 0.719,
  Y = -0.012
)

SEQPROFILE_WEIGHTS_CALIS_9MER <- c(
  0.00,
  0.00,
  0.10,
  0.31,
  0.30,
  0.29,
  0.26,
  0.18,
  0.00
)

# ── Immunological Scoring Helpers ─────────────────────────────────────

#' Mean amino acid scale value per window (Kolaskar, Parker, Chou-Fasman).
seqprofile_scale_mean <- function(window_seq, scale) {
  if (length(window_seq) == 0L) {
    return(numeric(0))
  }
  vapply(
    strsplit(toupper(window_seq), "", fixed = TRUE),
    function(aa) {
      if (length(aa) == 0L) {
        return(NA_real_)
      }
      vals <- unname(scale[aa])
      if (anyNA(vals)) NA_real_ else mean(vals)
    },
    numeric(1),
    USE.NAMES = FALSE
  )
}

#' Kolaskar & Tongaonkar (1990) antigenicity score per window.
score_kolaskar_tongaonkar <- function(window_seq) {
  seqprofile_scale_mean(window_seq, SEQPROFILE_SCALE_KOLASKAR)
}

#' Emini et al. (1985) surface accessibility score per window:
#' S_n = prod(d_i) * (0.37)^(-w)
score_emini_surface <- function(window_seq) {
  if (length(window_seq) == 0L) {
    return(numeric(0))
  }
  vapply(
    strsplit(toupper(window_seq), "", fixed = TRUE),
    function(aa) {
      w <- length(aa)
      if (w == 0L) {
        return(NA_real_)
      }
      vals <- unname(SEQPROFILE_SCALE_EMINI[aa])
      if (anyNA(vals)) NA_real_ else prod(vals) * (0.37)^(-w)
    },
    numeric(1),
    USE.NAMES = FALSE
  )
}

#' Parker et al. (1986) HPLC hydrophilicity score per window.
score_parker_hydrophilicity <- function(window_seq) {
  seqprofile_scale_mean(window_seq, SEQPROFILE_SCALE_PARKER)
}

#' Chou & Fasman (1978) beta-turn propensity score per window.
score_chou_fasman_turn <- function(window_seq) {
  seqprofile_scale_mean(window_seq, SEQPROFILE_SCALE_CHOU_FASMAN_TURN)
}

#' Parse a comma/space-separated mask specification (1-based indices or "C" for
#' C-terminus) into integer positions for a peptide of length `pep_len`.
parse_calis_mask <- function(mask, pep_len) {
  if (is.null(mask) || length(mask) == 0L) {
    return(integer(0))
  }
  if (is.numeric(mask)) {
    idx <- as.integer(mask[!is.na(mask)])
  } else {
    raw <- trimws(paste(as.character(mask), collapse = ","))
    if (!nzchar(raw)) {
      return(integer(0))
    }
    tokens <- toupper(strsplit(raw, "[,;[:space:]]+")[[1]])
    tokens <- tokens[nzchar(tokens)]
    idx <- vapply(
      tokens,
      function(tok) {
        if (tok %in% c("C", "CTERM", "C-TERM")) {
          as.integer(pep_len)
        } else if (grepl("^[0-9]+$", tok)) {
          as.integer(tok)
        } else {
          seqprofile_abort(
            sprintf(
              "Invalid mask token '%s'. Use 1-based positions and/or 'C' (e.g. '1,2,C').",
              tok
            ),
            "seqprofile_bad_param"
          )
        }
      },
      integer(1),
      USE.NAMES = FALSE
    )
  }
  unique(idx[idx >= 1L & idx <= pep_len])
}

#' Calis et al. (2013) MHC-I T-cell immunogenicity score for peptides of length >= 9.
#'
#' Replicates IEDB Class I Immunogenicity (`predict_immunogenicity.py`):
#' nonamers use the published 9-position weights; longer peptides insert the
#' central P5 weight (0.30) for each extra internal residue beyond 9.
score_calis_immunogenicity <- function(window_seq, mask = "1,2,C") {
  if (length(window_seq) == 0L) {
    return(numeric(0))
  }
  vapply(
    strsplit(toupper(window_seq), "", fixed = TRUE),
    function(aa) {
      L <- length(aa)
      if (L < 9L) {
        seqprofile_abort(
          "Calis MHC-I immunogenicity requires a window/peptide length of at least 9 residues.",
          "seqprofile_bad_param"
        )
      }
      vals <- unname(SEQPROFILE_SCALE_CALIS[aa])
      if (anyNA(vals)) {
        return(NA_real_)
      }
      w9 <- SEQPROFILE_WEIGHTS_CALIS_9MER
      weights <- if (L == 9L) {
        w9
      } else {
        c(w9[1:5], rep(w9[5], L - 9L), w9[6:9])
      }
      mask_idx <- parse_calis_mask(mask, L)
      if (length(mask_idx) > 0L) {
        weights[mask_idx] <- 0
      }
      sum(vals * weights)
    },
    numeric(1),
    USE.NAMES = FALSE
  )
}

# ── Classic Biophysical & IDR Scales / Helpers ───────────────────────────────

# Karplus & Schulz (1985) Naturwissenschaften 72:212-213.
# Source: AAindex KARP850102 (Flexibility parameter for one rigid neighbor;
# GenomeNet https://www.genome.jp/dbget-bin/www_bget?aaindex:KARP850102).
SEQPROFILE_SCALE_KS_FLEX <- c(
  A = 0.946, R = 1.028, N = 1.006, D = 1.089, C = 0.878,
  Q = 1.025, E = 1.036, G = 1.042, H = 0.952, I = 0.892,
  L = 0.961, K = 1.082, M = 0.862, F = 0.912, P = 1.085,
  S = 1.048, T = 1.051, W = 0.917, Y = 0.930, V = 0.927
)

#' Triangular smoothing weights for Karplus-Schulz flexibility.
#' For the default 7-residue window this gives c(0.25, 0.5, 0.75, 1, 0.75, 0.5, 0.25) / 4.
ks_triangular_weights <- function(w) {
  mid <- (w + 1) / 2
  dist <- abs(seq_len(w) - mid)
  raw <- 1 - dist / (mid)
  raw / sum(raw)
}

#' Karplus & Schulz (1985) backbone flexibility score per window.
score_karplus_schulz <- function(window_seq) {
  if (length(window_seq) == 0L) {
    return(numeric(0))
  }
  vapply(
    strsplit(toupper(window_seq), "", fixed = TRUE),
    function(aa) {
      w <- length(aa)
      if (w == 0L) {
        return(NA_real_)
      }
      vals <- unname(SEQPROFILE_SCALE_KS_FLEX[aa])
      if (anyNA(vals)) {
        return(NA_real_)
      }
      sum(vals * ks_triangular_weights(w))
    },
    numeric(1),
    USE.NAMES = FALSE
  )
}

#' Trapezoidal weights for von Heijne (1992) TopPred.
#' For window = 21: 5-residue flanks ramping linearly + 11-residue flat core.
toppred_trapezoid_weights <- function(w) {
  if (w == 21L) {
    flank <- 5L
    core <- 11L
  } else {
    flank <- max(1L, as.integer(floor(w / 4)))
    core <- max(1L, w - 2L * flank)
  }
  ramp <- seq_len(flank) / (flank + 1)
  raw <- c(ramp, rep(1, core), rev(ramp))
  raw / sum(raw)
}

#' von Heijne (1992) TopPred transmembrane hydrophobicity profile using the
#' Engelman-Steitz-Goldman (GES) scale from Peptides and a trapezoidal window.
score_von_heijne_tm <- function(window_seq) {
  if (length(window_seq) == 0L) {
    return(numeric(0))
  }
  # Ensure positive = hydrophobic (GES convention in TopPred)
  ges_vals <- vapply(
    SEQPROFILE_AA,
    function(a) as.numeric(Peptides::hydrophobicity(a, scale = "Engelman")),
    numeric(1)
  )
  if (ges_vals[["I"]] < 0) {
    ges_vals <- -ges_vals
  }

  vapply(
    strsplit(toupper(window_seq), "", fixed = TRUE),
    function(aa) {
      w <- length(aa)
      if (w == 0L) {
        return(NA_real_)
      }
      vals <- unname(ges_vals[aa])
      if (anyNA(vals)) {
        return(NA_real_)
      }
      sum(vals * toppred_trapezoid_weights(w))
    },
    numeric(1),
    USE.NAMES = FALSE
  )
}

#' Annotate von Heijne TopPred TM regions with certainty (putative >= 0.6,
#' certain >= 1.0) and positive-inside topology bias (K+R in 15-aa flanks).
annotate_von_heijne_regions <- function(regions, sequence, params) {
  if (nrow(regions) == 0L) {
    return(character(0))
  }
  n <- nchar(sequence)
  half_w <- (params$window - 1L) %/% 2L
  flank_len <- 15L

  count_kr <- function(s, e) {
    if (s > e || e < 1L || s > n) {
      return(0L)
    }
    subseq <- substring(sequence, max(1L, s), min(n, e))
    nchar(gsub("[^KR]", "", subseq))
  }

  vapply(
    seq_len(nrow(regions)),
    function(i) {
      cert <- if (regions$peak_score[i] >= 1.0) "Certain TM" else "Putative TM"
      tm_start <- max(1L, regions$start[i] - half_w)
      tm_end <- min(n, regions$end[i] + half_w)
      kr_n <- count_kr(tm_start - flank_len, tm_start - 1L)
      kr_c <- count_kr(tm_end + 1L, tm_end + flank_len)
      orient <- if (kr_n > kr_c) {
        "N-in/C-out"
      } else if (kr_c > kr_n) {
        "N-out/C-in"
      } else {
        "Balanced flanks"
      }
      sprintf("%s | %s (N-flank K+R=%d, C-flank K+R=%d)", cert, orient, kr_n, kr_c)
    },
    character(1)
  )
}

#' Jameson & Wolf (1988) composite Antigenic Index per window.
#'
#' Combines discretized Hopp-Woods hydrophilicity (weight 0.30), Emini surface
#' accessibility (0.15), Karplus-Schulz flexibility (0.15), and Chou-Fasman
#' beta-turn propensity (0.40, carrying the combined secondary-structure turn
#' weight because Garnier-Robson's 1,360-parameter matrix is not in Peptides).
score_jameson_wolf <- function(window_seq) {
  if (length(window_seq) == 0L) {
    return(numeric(0))
  }
  h_raw <- as.numeric(Peptides::hydrophobicity(window_seq, scale = "HoppWoods"))
  # Emini uses a 6-residue hexapeptide; evaluate inner 6 aa when window >= 6
  w_len <- nchar(window_seq[1])
  emini_sub <- if (w_len > 6L) {
    offset <- (w_len - 6L) %/% 2L
    substring(window_seq, 1L + offset, 6L + offset)
  } else {
    window_seq
  }
  s_raw <- score_emini_surface(emini_sub)
  f_raw <- score_karplus_schulz(window_seq)
  t_raw <- score_chou_fasman_turn(window_seq)

  # Table I discretization (Jameson & Wolf 1988, CABIOS 4:181-186)
  h_disc <- ifelse(
    h_raw > 0.5, 2,
    ifelse(h_raw > 0.0, 1, ifelse(h_raw > -0.4, -1, -2))
  )
  s_disc <- ifelse(s_raw >= 1.0, 1, 0)
  f_disc <- ifelse(f_raw >= 1.0, 1, 0)
  t_disc <- ifelse(t_raw >= 1.2, 2, ifelse(t_raw >= 1.0, 1, 0))

  0.30 * h_disc + 0.15 * s_disc + 0.15 * f_disc + 0.40 * t_disc
}

#' FoldIndex (Prilusky et al., 2005; Uversky et al., 2000) sliding-window IDR proxy.
#'
#' I_F = 2.785 * <H_scaled> - |<R>| - 1.151
#' where <H_scaled> is mean Kyte-Doolittle hydropathy rescaled to [0, 1]
#' via (KD + 4.5) / 9, and |<R>| is absolute mean net charge per residue.
#' Negative scores predict intrinsically disordered (unfolded) regions.
score_foldindex <- function(window_seq, pH = 7.4) {
  if (length(window_seq) == 0L) {
    return(numeric(0))
  }
  kd <- as.numeric(Peptides::hydrophobicity(window_seq, scale = "KyteDoolittle"))
  h_scaled <- (kd + 4.5) / 9.0
  abs_charge <- abs(
    as.numeric(Peptides::charge(window_seq, pH = pH, pKscale = "Lehninger")) /
      nchar(window_seq)
  )
  2.785 * h_scaled - abs_charge - 1.151
}

#' Registry of sequence-profiling algorithms.
#'
#' One entry per algorithm. `defaults` lists exactly the parameters the
#' algorithm takes (the module renders a control for each); `score_fn` receives
#' the window subsequences (canonical residues only) and the resolved
#' parameters and returns one score per window.
seqprofile_registry <- function() {
  list(
    flexibility_ks = list(
      id = "flexibility_ks",
      name = "Backbone Flexibility (Karplus-Schulz)",
      category = "Structural / Membrane",
      y_label = "Flexibility index (Bnorm)",
      defaults = list(window = 7L, threshold = 1.0),
      min_window = 3L,
      direction = "above",
      baseline = 1.0,
      needs_peptides = FALSE,
      region_note = "Flexible loop / segment",
      score_fn = function(window_seq, p) {
        score_karplus_schulz(window_seq)
      },
      interpretation = list(
        summary = "Triangularly weighted sliding-window backbone flexibility derived from normalized crystallographic alpha-carbon B-factors (AAindex KARP850102).",
        reading = "Values above 1.0 indicate mobile, solvent-exposed loops or linkers; values below 1.0 mark rigid secondary-structure cores.",
        thresholds = "The default 7-residue window uses Karplus-Schulz triangular weights (0.25, 0.50, 0.75, 1.00, 0.75, 0.50, 0.25) with a 1.0 baseline.",
        caveats = "Uses the single-table KARP850102 (one-rigid-neighbor) scale rather than full 3-table neighbor switching. High flexibility often overlaps with surface loops and linear B-cell epitopes.",
        reference = "Karplus PA, Schulz GE (1985) Naturwissenschaften 72:212-213; AAindex KARP850102."
      )
    ),
    tm_von_heijne = list(
      id = "tm_von_heijne",
      name = "Transmembrane Helices (von Heijne TopPred)",
      category = "Structural / Membrane",
      y_label = "Trapezoidal GES hydropathy (kcal/mol)",
      defaults = list(window = 21L, threshold = 0.6),
      min_window = 11L,
      direction = "above",
      baseline = 1.0,
      region_note = "Predicted TM helix",
      region_note_fn = annotate_von_heijne_regions,
      score_fn = function(window_seq, p) {
        score_von_heijne_tm(window_seq)
      },
      interpretation = list(
        summary = "TopPred sliding trapezoidal window (11-residue core + 5-residue flanking wedges = 21 aa) over the Engelman-Steitz-Goldman (GES) hydrophobicity scale, paired with positive-inside flanking K+R charge bias.",
        reading = "High positive peaks indicate nonpolar transbilayer alpha-helices. The regions table classifies each hit as Certain (>= 1.0) or Putative (>= 0.6) and reports N- vs C-flank K+R counts (positive-inside rule: the more positively charged loop tends to face the cytoplasm).",
        thresholds = "TopPred defines 0.6 as the lower (putative TM) cutoff (dashed red line) and 1.0 as the upper (certain TM) cutoff (solid grey line) over a 21-residue trapezoidal window.",
        caveats = "N-terminal hydrophobic signal peptides frequently exceed the 1.0 cutoff and can mimic a first transmembrane helix. Orientation calls are heuristic for single helices.",
        reference = "von Heijne G (1992) J Mol Biol 225:487-494; Engelman DM, Steitz TA, Goldman A (1986) Annu Rev Biophys Biophys Chem 15:321-353. GES scale via Peptides::hydrophobicity(scale = \"Engelman\")."
      )
    ),
    antigenicity_jw = list(
      id = "antigenicity_jw",
      name = "Antigenic Index (Jameson-Wolf)",
      category = "Immunogenicity / Epitopes",
      y_label = "Antigenic Index (AI)",
      defaults = list(window = 7L, threshold = 0.0),
      min_window = 6L,
      direction = "above",
      baseline = 0.0,
      region_note = "Predicted antigenic peak",
      score_fn = function(window_seq, p) {
        score_jameson_wolf(window_seq)
      },
      interpretation = list(
        summary = "Composite linear epitope index combining discretized Hopp-Woods hydrophilicity (0.30), Emini surface probability (0.15), Karplus-Schulz flexibility (0.15), and Chou-Fasman beta-turn propensity (0.40).",
        reading = "Positive peaks highlight hydrophilic, flexible, surface-exposed turns likely to act as linear B-cell antigenic determinants; negative regions correspond to buried or hydrophobic cores.",
        thresholds = "Default window is 7 residues (evaluating Emini over the inner 6 residues) with a 0.0 baseline. Raising the threshold to 0.5-1.0 isolates only the strongest composite peaks.",
        caveats = "The original 1988 algorithm splits the 0.40 secondary-structure turn weight evenly between Chou-Fasman (0.20) and Garnier-Robson (0.20); here all 0.40 is assigned to Chou-Fasman beta-turn propensity, and peak-broadening smoothing is omitted so raw scores stay transparent.",
        reference = "Jameson BA, Wolf H (1988) Comput Appl Biosci 4(1):181-186; Wolf H et al. (1988) Comput Appl Biosci 4(1):187-191."
      )
    ),
    idr_foldindex = list(
      id = "idr_foldindex",
      name = "Intrinsically Disordered Regions (FoldIndex / Uversky CH)",
      category = "Structural / Membrane",
      y_label = "FoldIndex (IF; < 0 = disordered)",
      defaults = list(window = 51L, pH = 7.4, threshold = 0.0),
      min_window = 11L,
      direction = "below",
      baseline = 0.0,
      region_note = "Predicted intrinsically disordered region (IDR)",
      score_fn = function(window_seq, p) {
        score_foldindex(window_seq, pH = p$pH)
      },
      interpretation = list(
        summary = "Sliding-window FoldIndex (IF = 2.785 * <H_scaled> - |<R>| - 1.151) derived from the Uversky charge-hydropathy boundary, combining rescaled Kyte-Doolittle hydropathy and absolute Lehninger net charge per residue.",
        reading = "Negative values (< 0, shaded) predict intrinsically disordered / unfolded regions characterized by low mean hydrophobicity and high net charge; positive values (> 0) predict ordered, folded domains.",
        thresholds = "Prilusky et al. (2005) recommend a 51-residue window and a 0.0 fold-boundary cutoff. For short peptides (< 51 aa), reduce the window size (e.g. 15-21 aa).",
        caveats = "FoldIndex is a fast biophysical proxy based on net charge and hydropathy; neutral polar IDRs (e.g. Q/N-rich prion-like domains) with low net charge can score near or above 0.",
        reference = "Prilusky J et al. (2005) Bioinformatics 21(16):3435-3438; Uversky VN, Gillespie JR, Fink AL (2000) Proteins 41(3):415-427. Computed with Peptides::hydrophobicity() and Peptides::charge()."
      )
    ),
    charge_lehninger = list(
      id = "charge_lehninger",
      name = "Charge Density (Lehninger)",
      category = "Physicochemical",
      y_label = "Net charge per residue",
      defaults = list(window = 9L, pH = 7.4, threshold = 0),
      direction = "above",
      baseline = 0,
      region_note = "Net positive segment",
      score_fn = function(window_seq, p) {
        Peptides::charge(window_seq, pH = p$pH, pKscale = "Lehninger") /
          nchar(window_seq)
      },
      interpretation = list(
        summary = "Net charge of each window at the chosen pH, from the Henderson-Hasselbalch equation with the Lehninger pKa set, divided by the window length.",
        reading = "Positive values mark basic (Lys/Arg/His-rich) stretches, negative values acidic (Asp/Glu-rich) stretches; values near 0 are neutral or charge-balanced.",
        thresholds = "The default window of 9 residues gives a local view; the 0 line separates net positive from net negative segments. Shaded regions are net positive at the chosen pH.",
        caveats = "Each window is scored as a free peptide, so one N-terminal amine and one C-terminal carboxyl are included in every window. pKa values are for free amino acids and ignore the local environment in a folded protein.",
        reference = "Moore DS (1985) Biochem Educ 13:10-11; pKa set from Nelson DL, Cox MM, Lehninger Principles of Biochemistry, 4th ed. (2004). Computed with Peptides::charge()."
      )
    ),
    hydropathy_kd = list(
      id = "hydropathy_kd",
      name = "Hydropathy (Kyte-Doolittle)",
      category = "Physicochemical",
      y_label = "Mean hydropathy (Kyte-Doolittle)",
      defaults = list(window = 9L, threshold = 0),
      direction = "above",
      baseline = 0,
      region_note = "Hydrophobic segment",
      score_fn = function(window_seq, p) {
        Peptides::hydrophobicity(window_seq, scale = "KyteDoolittle")
      },
      interpretation = list(
        summary = "Mean Kyte-Doolittle hydropathy of the residues in each window.",
        reading = "Positive = hydrophobic (likely buried or membrane-embedded), negative = hydrophilic (likely solvent-exposed). Note the sign is opposite to Hopp-Woods.",
        thresholds = "A window of 9 suits surface / interior patterns of globular proteins (threshold 0). To look for membrane-spanning segments use a window of 19 and a threshold of 1.6, as proposed by Kyte and Doolittle.",
        caveats = "A hydropathy peak is suggestive, not proof, of a transmembrane helix; signal peptides also score high.",
        reference = "Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132. Computed with Peptides::hydrophobicity(scale = \"KyteDoolittle\")."
      )
    ),
    hydrophilicity_hw = list(
      id = "hydrophilicity_hw",
      name = "Hydrophilicity (Hopp-Woods)",
      category = "Physicochemical",
      y_label = "Mean hydrophilicity (Hopp-Woods)",
      defaults = list(window = 6L, threshold = 0),
      direction = "above",
      baseline = 0,
      region_note = "Hydrophilic segment",
      score_fn = function(window_seq, p) {
        Peptides::hydrophobicity(window_seq, scale = "HoppWoods")
      },
      interpretation = list(
        summary = "Mean Hopp-Woods hydrophilicity of the residues in each window.",
        reading = "Positive = hydrophilic, likely surface-exposed; negative = hydrophobic. Note the sign is opposite to Kyte-Doolittle.",
        thresholds = "Hopp and Woods used a 6-residue window and looked for the highest peaks as candidate antigenic determinants; the 0 line separates net hydrophilic from net hydrophobic windows.",
        caveats = "With an even window the score is reported at the residue just left of the window midpoint. Hydrophilicity alone is a weak predictor of epitopes.",
        reference = "Hopp TP, Woods KR (1981) Proc Natl Acad Sci USA 78:3824-3828. Computed with Peptides::hydrophobicity(scale = \"HoppWoods\")."
      )
    ),
    hmoment_eisenberg = list(
      id = "hmoment_eisenberg",
      name = "Amphiphilicity (Eisenberg hydrophobic moment)",
      category = "Structural / Membrane",
      y_label = "Hydrophobic moment (uH)",
      defaults = list(window = 11L, angle = 100, threshold = NA_real_),
      direction = "above",
      baseline = NA_real_,
      region_note = "Amphiphilic segment",
      score_fn = function(window_seq, p) {
        Peptides::hmoment(window_seq, angle = p$angle, window = p$window)
      },
      interpretation = list(
        summary = "Hydrophobic moment of each window: how strongly hydrophobic residues cluster on one face of a periodic structure with the chosen rotation angle per residue.",
        reading = "Higher values = more amphiphilic (one hydrophobic and one polar face). Values are never negative. Low values mean hydrophobicity is spread evenly around the axis.",
        thresholds = "The default window is 11 residues. An angle of 100 degrees assumes an alpha-helix, 160 degrees a beta-strand. No cutoff is preset; enter a threshold to shade regions above it.",
        caveats = "The result is only meaningful if the segment really adopts the assumed periodicity; a high moment at 100 degrees does not show that the segment is helical.",
        reference = "Eisenberg D, Weiss RM, Terwilliger TC (1984) Proc Natl Acad Sci USA 81:140-144. Computed with Peptides::hmoment()."
      )
    ),
    antigenicity_kolaskar = list(
      id = "antigenicity_kolaskar",
      name = "Antigenicity (Kolaskar-Tongaonkar)",
      category = "Immunogenicity / Epitopes",
      y_label = "Antigenic propensity",
      defaults = list(window = 7L, threshold = 1.0),
      direction = "above",
      baseline = 1.0,
      needs_peptides = FALSE,
      region_note = "Predicted antigenic determinant",
      score_fn = function(window_seq, p) {
        score_kolaskar_tongaonkar(window_seq)
      },
      interpretation = list(
        summary = "Semi-empirical antigenic propensity averaged over each window, derived from amino acid frequencies in experimentally known linear B-cell epitopes versus protein surfaces.",
        reading = "Values above 1.0 indicate residues enriched in antigenic determinants (notably Cys, Val, Leu, Tyr, Ile). Higher peaks mark stronger linear epitope candidates.",
        thresholds = "The default 7-residue window and 1.0 parity line match the Kolaskar-Tongaonkar convention; IEDB also commonly uses the protein's own mean score (shown in the summary strip) as a sequence-specific cutoff.",
        caveats = "Hydrophobic residues such as Val, Leu, and Cys carry high propensities, so buried hydrophobic cores can score high. Cross-check peaks with Emini surface accessibility or Parker hydrophilicity.",
        reference = "Kolaskar AS, Tongaonkar PC (1990) FEBS Lett 276(1-2):172-174; scale values from IEDB B-Cell Analysis Resource."
      )
    ),
    surface_emini = list(
      id = "surface_emini",
      name = "Surface Accessibility (Emini)",
      category = "Structural / Membrane",
      y_label = "Surface probability (Sn)",
      defaults = list(window = 6L, threshold = 1.0),
      direction = "above",
      baseline = 1.0,
      needs_peptides = FALSE,
      region_note = "Surface-exposed segment",
      score_fn = function(window_seq, p) {
        score_emini_surface(window_seq)
      },
      interpretation = list(
        summary = "Multiplicative surface probability Sn = prod(d_i) * (0.37)^(-w) of each window, estimating the likelihood that a stretch is exposed on the protein surface.",
        reading = "Scores above 1.0 indicate an increased probability of surface exposure relative to an average hexapeptide. Because the formula is multiplicative, a single strongly buried residue sharply suppresses the window score.",
        thresholds = "Emini et al. defined a 6-residue window with a threshold of 1.0. With an even window, the score is assigned to the residue just left of the midpoint (position 3 of residues 1-6).",
        caveats = "Scores can span several orders of magnitude in highly charged loops. Surface exposure is necessary for antibody binding but does not guarantee immunogenicity on its own.",
        reference = "Emini EA, Hughes JV, Perlow DS, Boger J (1985) J Virol 55(3):836-839; scale values from IEDB B-Cell Analysis Resource."
      )
    ),
    hydrophilicity_parker = list(
      id = "hydrophilicity_parker",
      name = "Hydrophilicity (Parker HPLC)",
      category = "Physicochemical",
      y_label = "Mean hydrophilicity (Parker)",
      defaults = list(window = 7L, threshold = 0),
      direction = "above",
      baseline = 0,
      needs_peptides = FALSE,
      region_note = "Hydrophilic / surface loop",
      score_fn = function(window_seq, p) {
        score_parker_hydrophilicity(window_seq)
      },
      interpretation = list(
        summary = "Mean hydrophilicity per window based on reversed-phase HPLC peptide retention times (scale ranges from -10.0 for Trp to +10.0 for Asp).",
        reading = "Positive peaks mark hydrophilic, solvent-exposed stretches often associated with linear B-cell epitopes; negative valleys correspond to hydrophobic segments.",
        thresholds = "Parker et al. used a 7-residue window with the score assigned to the central (4th) residue. The 0 line separates net hydrophilic from net hydrophobic windows; you can raise the threshold to the sequence mean to isolate top peaks.",
        caveats = "Hydrophilicity correlates with solvent accessibility rather than antibody binding directly; combine with beta-turn and surface accessibility profiles.",
        reference = "Parker JMR, Guo D, Hodges RS (1986) Biochemistry 25(19):5425-5432; scale values from IEDB B-Cell Analysis Resource."
      )
    ),
    betaturn_chou_fasman = list(
      id = "betaturn_chou_fasman",
      name = "Beta-Turn Propensity (Chou-Fasman)",
      category = "Structural / Membrane",
      y_label = "Mean beta-turn propensity (Pt)",
      defaults = list(window = 7L, threshold = 1.0),
      direction = "above",
      baseline = 1.0,
      needs_peptides = FALSE,
      region_note = "Predicted beta-turn region",
      score_fn = function(window_seq, p) {
        score_chou_fasman_turn(window_seq)
      },
      interpretation = list(
        summary = "Mean Chou-Fasman beta-turn conformational parameter (Pt) over each window, highlighting regions enriched in turn-forming residues (Gly, Asn, Pro, Asp, Ser).",
        reading = "Values above 1.0 indicate above-average beta-turn propensity. Exposed beta-turns frequently coincide with linear B-cell epitopes.",
        thresholds = "The IEDB sliding-window convention uses a 7-residue window. A threshold of 1.0 marks windows above random-coil expectation; the sequence mean can be used for stricter peak picking.",
        caveats = "This profile evaluates single-residue turn propensity only, not full 4-residue Chou-Fasman positional bend rules or competing alpha-helix/beta-sheet propensities.",
        reference = "Chou PY, Fasman GD (1978) Adv Enzymol Relat Areas Mol Biol 47:45-148; scale values from IEDB B-Cell Analysis Resource."
      )
    ),
    immunogenicity_calis = list(
      id = "immunogenicity_calis",
      name = "MHC-I T-Cell Immunogenicity (Calis)",
      category = "Immunogenicity / Epitopes",
      y_label = "Immunogenicity score",
      defaults = list(window = 9L, mask = "1,2,C", threshold = 0),
      min_window = 9L,
      direction = "above",
      baseline = 0,
      needs_peptides = FALSE,
      region_note = "Predicted immunogenic MHC-I peptide",
      score_fn = function(window_seq, p) {
        score_calis_immunogenicity(window_seq, mask = p$mask)
      },
      interpretation = list(
        summary = "Position-weighted sum of amino acid log-enrichment scores across each candidate MHC class I peptide (9-mer or longer), emphasizing TCR-facing residues (positions 4-6) enriched in large/aromatic amino acids (Trp, Ile, Phe, Glu).",
        reading = "Higher positive scores indicate a greater likelihood of T-cell recognition once presented on MHC-I; negative scores indicate depletion of immunogenic features.",
        thresholds = "Default window is 9 residues (standard MHC-I peptide length) with anchor positions 1, 2, and C-terminus masked ('1,2,C') and threshold 0. Longer windows (>=10 aa) replicate IEDB's central-weight insertion.",
        caveats = "Calis et al. models TCR recognition of peptides already bound to MHC-I; it does not predict HLA binding affinity or proteasomal cleavage. Combine top windows with an MHC-I binding predictor.",
        reference = "Calis JJA et al. (2013) PLoS Comput Biol 9(10):e1003266 (IEDB Class I Immunogenicity)."
      )
    )
  )
}

# Display order of the algorithm groups in the selectInput().
SEQPROFILE_CATEGORIES <- c(
  "Physicochemical",
  "Structural / Membrane",
  "Immunogenicity / Epitopes"
)

# Choices for a grouped selectInput(): list(category = list(name = id, ...)).
# Groups follow SEQPROFILE_CATEGORIES (any other category is appended);
# algorithms keep their registry order within a group.
seqprofile_choices <- function() {
  reg <- seqprofile_registry()
  cats <- vapply(reg, function(x) x$category, character(1))
  lvls <- union(SEQPROFILE_CATEGORIES, cats)
  lapply(split(reg, factor(cats, levels = lvls), drop = TRUE), function(group) {
    as.list(stats::setNames(
      vapply(group, function(x) x$id, character(1)),
      vapply(group, function(x) x$name, character(1))
    ))
  })
}

#' Run one profiling algorithm on a cleaned sequence.
#'
#' Windows containing a residue outside the 20 canonical amino acids get an
#' `NA` score (Peptides only warns and still returns a number for them).
#'
#' @param sequence Cleaned sequence, e.g. `clean_protein_sequence()$sequence`.
#' @param algorithm_id A name of `seqprofile_registry()`.
#' @param params Named list overriding the algorithm's defaults; unknown names
#'   and NULL entries are ignored.
#' @return Named list: algorithm_id, algorithm_name, sequence, params, profile,
#'   regions, global_stats, interpretation, plus y_label, direction, baseline
#'   and masked_residues.
run_seq_profile <- function(sequence, algorithm_id, params = list()) {
  spec <- seqprofile_registry()[[algorithm_id]]
  if (is.null(spec)) {
    seqprofile_abort(
      paste0("Unknown algorithm: ", algorithm_id),
      "seqprofile_bad_param"
    )
  }
  if (
    !identical(spec$needs_peptides, FALSE) &&
      !requireNamespace("Peptides", quietly = TRUE)
  ) {
    seqprofile_abort(
      "The 'Peptides' package is not installed. Start proteOmni through run.R so it is installed from proteOmni.lock.",
      "seqprofile_missing_package"
    )
  }
  sequence <- toupper(as.character(sequence)[1])
  if (is.na(sequence) || !nzchar(sequence)) {
    seqprofile_abort("The sequence is empty.", "seqprofile_empty")
  }

  p <- spec$defaults
  for (nm in intersect(names(params), names(p))) {
    value <- params[[nm]]
    if (length(value) == 1) p[[nm]] <- value
  }
  min_w <- if (!is.null(spec$min_window)) as.integer(spec$min_window) else 2L
  if (is.na(p$window) || p$window < min_w || p$window != round(p$window)) {
    seqprofile_abort(
      sprintf("Window size must be a whole number of at least %d.", min_w),
      "seqprofile_bad_param"
    )
  }
  p$window <- as.integer(p$window)
  if (!is.null(p$pH) && (is.na(p$pH) || p$pH < 0 || p$pH > 14)) {
    seqprofile_abort("pH must be between 0 and 14.", "seqprofile_bad_param")
  }
  if (!is.null(p$angle) && is.na(p$angle)) {
    seqprofile_abort("Rotation angle is missing.", "seqprofile_bad_param")
  }
  p$threshold <- as.numeric(p$threshold)

  profile <- seqprofile_windows(sequence, p$window)
  scorable <- !grepl("[^ACDEFGHIKLMNPQRSTVWY]", profile$window_seq)
  profile$score <- NA_real_
  if (any(scorable)) {
    profile$score[scorable] <- as.numeric(
      spec$score_fn(profile$window_seq[scorable], p)
    )
  }
  profile$above_thr <- if (is.na(p$threshold)) {
    NA
  } else if (spec$direction == "above") {
    profile$score >= p$threshold
  } else {
    profile$score <= p$threshold
  }

  regions <- extract_profile_regions(profile, p$threshold, spec$direction)
  regions$note <- if (!is.null(spec$region_note_fn) && nrow(regions) > 0L) {
    spec$region_note_fn(regions, sequence, p)
  } else {
    rep(spec$region_note, nrow(regions))
  }

  has_score <- any(scorable)
  i_min <- if (has_score) which.min(profile$score) else NA_integer_
  i_max <- if (has_score) which.max(profile$score) else NA_integer_
  global_stats <- list(
    seq_length = nchar(sequence),
    mean_score = if (has_score) mean(profile$score, na.rm = TRUE) else NA_real_,
    sd_score = if (has_score) {
      stats::sd(profile$score, na.rm = TRUE)
    } else {
      NA_real_
    },
    min_score = if (has_score) profile$score[i_min] else NA_real_,
    max_score = if (has_score) profile$score[i_max] else NA_real_,
    n_regions = nrow(regions),
    # Share of scored windows, not of all residues
    pct_above = if (has_score && !is.na(p$threshold)) {
      100 * mean(profile$above_thr, na.rm = TRUE)
    } else {
      NA_real_
    },
    min_pos = if (has_score) profile$position[i_min] else NA_integer_,
    max_pos = if (has_score) profile$position[i_max] else NA_integer_,
    n_windows = nrow(profile),
    n_masked = sum(!scorable)
  )

  list(
    algorithm_id = spec$id,
    algorithm_name = spec$name,
    sequence = sequence,
    params = p,
    profile = profile,
    regions = regions,
    global_stats = global_stats,
    interpretation = spec$interpretation,
    y_label = spec$y_label,
    direction = spec$direction,
    baseline = spec$baseline,
    masked_residues = sort(setdiff(
      unique(strsplit(sequence, "")[[1]]),
      SEQPROFILE_AA
    ))
  )
}

#' Whole-sequence descriptors from Peptides.
#'
#' @return Named list (mw, pI, gravy, boman, aliphatic_index,
#'   instability_index); all `NA` when the sequence holds a non-canonical
#'   residue, because Peptides does not score those reliably.
seqprofile_global_descriptors <- function(sequence) {
  nms <- c(
    "mw",
    "pI",
    "gravy",
    "boman",
    "aliphatic_index",
    "instability_index"
  )
  sequence <- toupper(as.character(sequence)[1])
  if (
    is.na(sequence) ||
      !nzchar(sequence) ||
      grepl("[^ACDEFGHIKLMNPQRSTVWY]", sequence) ||
      !requireNamespace("Peptides", quietly = TRUE)
  ) {
    return(stats::setNames(as.list(rep(NA_real_, length(nms))), nms))
  }
  list(
    mw = Peptides::mw(sequence),
    pI = Peptides::pI(sequence, pKscale = "Lehninger"),
    gravy = Peptides::hydrophobicity(sequence, scale = "KyteDoolittle"),
    boman = Peptides::boman(sequence),
    aliphatic_index = Peptides::aIndex(sequence),
    instability_index = Peptides::instaIndex(sequence)
  )
}

#' Profile plot for a `run_seq_profile()` result.
#'
#' Line of score vs. window-centre position, a dashed threshold line, a solid
#' baseline where the algorithm has one, and shaded bands over the predicted
#' regions. The points carry a `text` aesthetic for plotly hover.
#'
#' @param consensus Optional `seqprofile_consensus_regions()` output; its
#'   intervals are shaded in orange on top of the algorithm's own regions, so
#'   the same bands line up across stacked tracks.
plot_seq_profile <- function(result, consensus = NULL) {
  prof <- result$profile
  if (all(is.na(prof$score))) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No window could be scored (non-canonical residues in every window)"
        ) +
        ggplot2::theme_void()
    )
  }
  thr <- result$params$threshold
  base <- result$baseline
  prof$hover <- sprintf(
    "Position: %d<br>Residue: %s<br>Window: %s<br>Score: %.3f",
    prof$position,
    prof$residue,
    prof$window_seq,
    prof$score
  )
  # Finite band limits: plotly does not draw rectangles with infinite bounds
  y_rng <- range(c(prof$score, thr, base), na.rm = TRUE)
  pad <- diff(y_rng) * 0.05
  if (pad == 0) {
    pad <- 0.5
  }
  y_rng <- y_rng + c(-pad, pad)

  p <- ggplot2::ggplot(prof, ggplot2::aes(x = position, y = score))
  if (nrow(result$regions) > 0) {
    bands <- data.frame(
      xmin = result$regions$start - 0.5,
      xmax = result$regions$end + 0.5,
      ymin = y_rng[1],
      ymax = y_rng[2]
    )
    p <- p +
      ggplot2::geom_rect(
        data = bands,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "#5FA8D3",
        alpha = 0.25
      )
  }
  if (!is.null(consensus) && nrow(consensus) > 0) {
    consensus_bands <- data.frame(
      xmin = consensus$start - 0.5,
      xmax = consensus$end + 0.5,
      ymin = y_rng[1],
      ymax = y_rng[2]
    )
    p <- p +
      ggplot2::geom_rect(
        data = consensus_bands,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "#E8A33D",
        alpha = 0.45
      )
  }
  if (!is.na(base)) {
    p <- p + ggplot2::geom_hline(yintercept = base, colour = "grey60")
  }
  if (!is.na(thr)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = thr,
        linetype = "dashed",
        colour = "#c0392b"
      )
  }
  p +
    ggplot2::geom_line(colour = "#1B4965", na.rm = TRUE) +
    # `text` is not a ggplot2 aesthetic; it is read by plotly::ggplotly()
    suppressWarnings(ggplot2::geom_point(
      ggplot2::aes(text = hover),
      colour = "#1B4965",
      size = 0.7,
      na.rm = TRUE
    )) +
    ggplot2::labs(
      title = result$algorithm_name,
      x = "Residue position (window centre)",
      y = result$y_label
    )
}

# ── UniProt retrieval ────────────────────────────────────────────────────────

# UniProtKB accession format (uniprot.org/help/accession_numbers), with an
# optional isoform suffix such as P05067-2.
SEQPROFILE_UNIPROT_ACCESSION <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$"

# Trims and upper-cases an accession; aborts when it is not a UniProtKB one, so
# nothing typed by the user reaches the request URL unchecked.
seqprofile_clean_accession <- function(accession) {
  acc <- toupper(trimws(paste(as.character(accession), collapse = "")))
  if (!grepl(SEQPROFILE_UNIPROT_ACCESSION, acc)) {
    seqprofile_abort(
      "Enter a UniProtKB accession such as Q07730 or P05067-2.",
      "seqprofile_bad_accession"
    )
  }
  acc
}

#' Split the text of a single-entry UniProt FASTA file.
#'
#' @return list(fasta, header, protein_name, sequence). `protein_name` is the
#'   header text between the identifier and the first OS= / OX= / GN= / PE= /
#'   SV= field.
parse_uniprot_fasta <- function(text) {
  text <- trimws(paste(as.character(text), collapse = "\n"))
  cs <- clean_protein_sequence(text)
  if (!startsWith(text, ">") || !nzchar(cs$sequence)) {
    seqprofile_abort(
      "UniProt returned no FASTA sequence for this accession (the entry may be obsolete).",
      "seqprofile_fetch_failed"
    )
  }
  description <- sub("^\\S+\\s*", "", cs$header)
  list(
    fasta = gsub("\r\n|\r", "\n", text),
    header = cs$header,
    protein_name = trimws(sub("\\s+(OS|OX|GN|PE|SV)=.*$", "", description)),
    sequence = cs$sequence
  )
}

#' Download one entry from the UniProt REST API as FASTA.
#'
#' Requests rest.uniprot.org/uniprotkb/{accession}.fasta, which returns the
#' canonical sequence (or the named isoform when the accession has a suffix).
#'
#' @return `parse_uniprot_fasta()` output plus `accession`.
fetch_uniprot_fasta <- function(accession, timeout = 15) {
  acc <- seqprofile_clean_accession(accession)
  if (!requireNamespace("httr", quietly = TRUE)) {
    seqprofile_abort(
      "The 'httr' package is not installed. Start proteOmni through run.R so it is installed from proteOmni.lock.",
      "seqprofile_missing_package"
    )
  }
  url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", acc)
  resp <- tryCatch(
    httr::GET(url, httr::timeout(timeout)),
    error = function(e) e
  )
  if (inherits(resp, "error")) {
    seqprofile_abort(
      paste0("Could not reach UniProt: ", conditionMessage(resp)),
      "seqprofile_fetch_failed"
    )
  }
  status <- httr::status_code(resp)
  if (status %in% c(400L, 404L)) {
    seqprofile_abort(
      sprintf("UniProt has no entry for accession %s.", acc),
      "seqprofile_fetch_failed"
    )
  }
  if (status != 200L) {
    seqprofile_abort(
      sprintf("The UniProt request failed (HTTP %d).", status),
      "seqprofile_fetch_failed"
    )
  }
  out <- parse_uniprot_fasta(
    httr::content(resp, as = "text", encoding = "UTF-8")
  )
  out$accession <- acc
  out
}

# ── Consensus regions across several profiles ────────────────────────────────

#' Intervals where every profile crosses its own threshold.
#'
#' A position counts when it is the window centre of a window that passes the
#' threshold in every result. Positions that one algorithm does not score (the
#' sequence ends, where its window does not fit, or masked windows) never count.
#'
#' @param results List of `run_seq_profile()` results for the same sequence.
#' @param min_length Shortest stretch (in positions) reported as a region;
#'   regions are numbered after this filter.
#' @return data.frame(region_id, start, end, length, sequence, mean_<id>...),
#'   with one mean-score column per algorithm.
seqprofile_consensus_regions <- function(results, min_length = 1L) {
  if (length(results) < 2) {
    seqprofile_abort(
      "Consensus regions need at least two algorithms.",
      "seqprofile_bad_param"
    )
  }
  if (
    length(min_length) != 1 ||
      is.na(min_length) ||
      min_length < 1 ||
      min_length != round(min_length)
  ) {
    seqprofile_abort(
      "The minimum consensus length must be a whole number of at least 1.",
      "seqprofile_bad_param"
    )
  }
  no_thr <- vapply(results, function(r) is.na(r$params$threshold), logical(1))
  if (any(no_thr)) {
    seqprofile_abort(
      sprintf(
        "Consensus regions need a threshold for every algorithm. None is set for: %s.",
        paste(
          vapply(results[no_thr], function(r) r$algorithm_name, character(1)),
          collapse = ", "
        )
      ),
      "seqprofile_bad_param"
    )
  }
  sequence <- results[[1]]$sequence
  same_seq <- vapply(
    results,
    function(r) identical(r$sequence, sequence),
    logical(1)
  )
  if (!all(same_seq)) {
    seqprofile_abort(
      "Consensus regions need profiles of the same sequence.",
      "seqprofile_bad_param"
    )
  }

  n <- nchar(sequence)
  k <- length(results)
  hit <- matrix(FALSE, nrow = n, ncol = k)
  score <- matrix(NA_real_, nrow = n, ncol = k)
  for (j in seq_len(k)) {
    prof <- results[[j]]$profile
    hit[prof$position, j] <- prof$above_thr %in% TRUE
    score[prof$position, j] <- prof$score
  }
  runs <- rle(rowSums(hit) == k)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  keep <- which(runs$values & runs$lengths >= min_length)

  out <- data.frame(
    region_id = seq_along(keep),
    start = starts[keep],
    end = ends[keep],
    length = runs$lengths[keep],
    # substring() rejects zero-length positions when no region is found
    sequence = if (length(keep) > 0) {
      substring(sequence, starts[keep], ends[keep])
    } else {
      character(0)
    },
    stringsAsFactors = FALSE
  )
  for (j in seq_len(k)) {
    out[[paste0("mean_", results[[j]]$algorithm_id)]] <- vapply(
      keep,
      function(r) mean(score[starts[r]:ends[r], j]),
      numeric(1)
    )
  }
  out
}

# ── Helical wheel (Schiffer-Edmundson projection) ────────────────────────────

#' Residues to draw for one predicted region of a profile.
#'
#' A region lists window centres, so the stretch that produced it also includes
#' the flanks of its first and last window. Stretches longer than `max_len` are
#' cut to `wheel_len` residues around the region's peak position.
#'
#' @param region Row number in `result$regions`.
#' @return list(start, end, sequence, truncated).
seqprofile_wheel_segment <- function(
  result,
  region,
  max_len = 36L,
  wheel_len = 18L
) {
  regions <- result$regions
  if (
    length(region) != 1 || is.na(region) || region < 1 || region > nrow(regions)
  ) {
    seqprofile_abort("Unknown region.", "seqprofile_bad_param")
  }
  window <- result$params$window
  left <- (window - 1L) %/% 2L
  lo <- regions$start[region] - left
  hi <- regions$end[region] + (window - 1L - left)
  truncated <- (hi - lo + 1L) > max_len
  if (truncated) {
    centred <- regions$peak_pos[region] - (wheel_len - 1L) %/% 2L
    lo <- max(lo, min(centred, hi - wheel_len + 1L))
    hi <- lo + wheel_len - 1L
  }
  list(
    start = as.integer(lo),
    end = as.integer(hi),
    sequence = substring(result$sequence, lo, hi),
    truncated = truncated
  )
}

#' Helical wheel coordinates of a segment.
#'
#' Residue i is placed `angle` degrees clockwise from residue i - 1, starting at
#' the top. At 100 degrees residue 19 falls on residue 1, so every `per_ring`
#' residues move one ring outwards. Hydrophobicity is the Eisenberg consensus
#' scale read from Peptides, the scale Peptides::hmoment() uses.
#'
#' @param start Sequence position of the first residue of `segment`.
#' @return data.frame(position, residue, angle_deg, x, y, hydrophobicity) with
#'   attribute `moment` = list(x, y, uH): the mean hydrophobicity vector of the
#'   whole segment and its length (the hydrophobic moment).
seqprofile_helical_wheel <- function(
  segment,
  start = 1L,
  angle = 100,
  per_ring = 18L
) {
  segment <- toupper(as.character(segment)[1])
  if (
    is.na(segment) ||
      nchar(segment) < 3 ||
      grepl("[^ACDEFGHIKLMNPQRSTVWY]", segment)
  ) {
    seqprofile_abort(
      "A helical wheel needs at least 3 canonical residues.",
      "seqprofile_bad_param"
    )
  }
  if (!requireNamespace("Peptides", quietly = TRUE)) {
    seqprofile_abort(
      "The 'Peptides' package is not installed. Start proteOmni through run.R so it is installed from proteOmni.lock.",
      "seqprofile_missing_package"
    )
  }
  aa <- strsplit(segment, "")[[1]]
  i <- seq_along(aa) - 1L
  theta <- i * angle * pi / 180
  h <- as.numeric(Peptides::hydrophobicity(aa, scale = "Eisenberg"))
  radius <- 1 + 0.3 * (i %/% per_ring)
  wheel <- data.frame(
    position = as.integer(start) + i,
    residue = aa,
    angle_deg = (i * angle) %% 360,
    x = radius * sin(theta),
    y = radius * cos(theta),
    hydrophobicity = h,
    stringsAsFactors = FALSE
  )
  mx <- mean(h * sin(theta))
  my <- mean(h * cos(theta))
  attr(wheel, "moment") <- list(x = mx, y = my, uH = sqrt(mx^2 + my^2))
  wheel
}

#' Helical wheel plot for a `seqprofile_helical_wheel()` result.
#'
#' Residues are filled by Eisenberg hydrophobicity and joined in sequence
#' order; the arrow is the hydrophobic moment vector (length = uH), pointing at
#' the hydrophobic face.
plot_helical_wheel <- function(wheel, title = NULL, subtitle = NULL) {
  moment <- attr(wheel, "moment")
  lim <- max(abs(c(wheel$x, wheel$y))) + 0.25
  ggplot2::ggplot(wheel, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_path(colour = "grey75", linewidth = 0.4) +
    ggplot2::geom_point(
      ggplot2::aes(fill = hydrophobicity),
      shape = 21,
      size = 10,
      colour = "grey30"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(residue, "\n", position)),
      size = 3,
      lineheight = 0.85
    ) +
    ggplot2::annotate(
      "segment",
      x = 0,
      y = 0,
      xend = moment$x,
      yend = moment$y,
      colour = "#c0392b",
      linewidth = 1,
      arrow = grid::arrow(length = grid::unit(0.12, "inches"))
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#5FA8D3",
      mid = "white",
      high = "#E8A33D",
      midpoint = 0,
      name = "Eisenberg\nhydrophobicity"
    ) +
    ggplot2::coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

# ── In silico digestion of one sequence (Digestion tab of mod_SeqProfile.r) ──

# Monoisotopic residue masses in Da (amino acid minus water), unmodified.
SEQPROFILE_MONO_MASS <- c(
  G = 57.021464,
  A = 71.037114,
  S = 87.032028,
  P = 97.052764,
  V = 99.068414,
  T = 101.047679,
  C = 103.009185,
  L = 113.084064,
  I = 113.084064,
  N = 114.042927,
  D = 115.026943,
  Q = 128.058578,
  K = 128.094963,
  E = 129.042593,
  M = 131.040485,
  H = 137.058912,
  F = 147.068414,
  R = 156.101111,
  Y = 163.063329,
  W = 186.079313
)
SEQPROFILE_WATER_MONO <- 18.010565
SEQPROFILE_PROTON_MONO <- 1.007276467
# Carbamidomethylation of Cys (iodoacetamide alkylation), C2H3NO
SEQPROFILE_CARBAMIDOMETHYL_MONO <- 57.021464
# Oxidation of Met, +O
SEQPROFILE_OXIDATION_MONO <- 15.994915

# Preset cleavage rules. `terminus` is the side of `residues` that is cut;
# `blocked_by` holds residues that prevent cleavage when they sit on the other
# side of the bond (e.g. trypsin does not cut K|P).
SEQPROFILE_ENZYMES <- list(
  trypsin = list(
    name = "Trypsin",
    residues = "KR",
    terminus = "C",
    blocked_by = "P"
  ),
  gluc = list(name = "GluC", residues = "E", terminus = "C", blocked_by = "P"),
  lysn = list(name = "LysN", residues = "K", terminus = "N", blocked_by = ""),
  aspn = list(name = "AspN", residues = "D", terminus = "N", blocked_by = ""),
  chymotrypsin = list(
    name = "Chymotrypsin",
    residues = "FWYL",
    terminus = "C",
    blocked_by = "P"
  ),
  elastase = list(
    name = "Elastase",
    residues = "AGILV",
    terminus = "C",
    blocked_by = "P"
  )
)

# Choices for the enzyme selectInput(); "custom" is handled by the module.
seqprofile_enzyme_choices <- function() {
  c(
    stats::setNames(
      names(SEQPROFILE_ENZYMES),
      vapply(SEQPROFILE_ENZYMES, function(x) x$name, character(1))
    ),
    "Custom rule" = "custom"
  )
}

# "KR", "k, r" or c("K", "R") -> c("K", "R"); aborts on anything but letters.
seqprofile_parse_residues <- function(x, what = "Cleavage residues") {
  if (is.null(x) || all(is.na(x))) {
    x <- ""
  }
  x <- toupper(gsub("[[:space:],;/|]", "", paste(x, collapse = "")))
  if (!grepl("^[A-Z]*$", x)) {
    seqprofile_abort(
      paste0(what, " must be one-letter amino acid codes, e.g. KR."),
      "seqprofile_bad_param"
    )
  }
  unique(strsplit(x, "")[[1]])
}

# Plain-language description of a cleavage rule.
seqprofile_enzyme_rule <- function(residues, terminus, blocked_by = "") {
  residues <- seqprofile_parse_residues(residues)
  blocked_by <- seqprofile_parse_residues(blocked_by)
  rule <- sprintf(
    "Cleaves %s-terminal to %s",
    terminus,
    paste(residues, collapse = "/")
  )
  if (length(blocked_by) > 0) {
    rule <- sprintf(
      "%s, except when %s by %s",
      rule,
      if (terminus == "C") "followed" else "preceded",
      paste(blocked_by, collapse = "/")
    )
  }
  rule
}

# Residue masses with the chosen modifications applied to every Cys / Met.
seqprofile_residue_masses <- function(
  carbamidomethyl = FALSE,
  oxidation = FALSE
) {
  masses <- SEQPROFILE_MONO_MASS
  if (isTRUE(carbamidomethyl)) {
    masses[["C"]] <- masses[["C"]] + SEQPROFILE_CARBAMIDOMETHYL_MONO
  }
  if (isTRUE(oxidation)) {
    masses[["M"]] <- masses[["M"]] + SEQPROFILE_OXIDATION_MONO
  }
  masses
}

#' Monoisotopic mass of peptides (neutral, [M]).
#'
#' @param masses Residue masses; the default is unmodified, see
#'   `seqprofile_residue_masses()` for modified Cys / Met.
#' @return Numeric vector; `NA` for a peptide with a non-canonical residue.
seqprofile_mono_mass <- function(peptides, masses = SEQPROFILE_MONO_MASS) {
  vapply(
    strsplit(as.character(peptides), ""),
    function(aa) {
      if (length(aa) == 0) {
        return(NA_real_)
      }
      sum(masses[aa]) + SEQPROFILE_WATER_MONO
    },
    numeric(1)
  )
}

#' Theoretical b and y fragment ions of one peptide.
#'
#' b(i) holds the first i residues, y(i) the last i residues plus water; both
#' carry `charge` protons. Ions run from 1 to length - 1.
#'
#' @return data.frame(ion, type, number, fragment, mz, intensity), with
#'   `intensity` fixed at 1.
seqprofile_fragment_ions <- function(
  peptide,
  carbamidomethyl = FALSE,
  oxidation = FALSE,
  charge = 1L
) {
  peptide <- toupper(as.character(peptide)[1])
  aa <- strsplit(peptide, "")[[1]]
  n <- length(aa)
  if (n < 2) {
    seqprofile_abort(
      "A peptide needs at least 2 residues to have b/y fragment ions.",
      "seqprofile_short_sequence"
    )
  }
  res_mass <- seqprofile_residue_masses(carbamidomethyl, oxidation)[aa]
  if (anyNA(res_mass)) {
    seqprofile_abort(
      "Fragment ions cannot be computed for a peptide with non-canonical residues.",
      "seqprofile_bad_param"
    )
  }
  if (length(charge) != 1 || is.na(charge) || charge < 1) {
    seqprofile_abort(
      "Fragment charge must be 1 or more.",
      "seqprofile_bad_param"
    )
  }
  i <- seq_len(n - 1L)
  b <- cumsum(res_mass)[i]
  y <- cumsum(rev(res_mass))[i] + SEQPROFILE_WATER_MONO
  suffix <- if (charge > 1) paste0(" ", charge, "+") else ""
  data.frame(
    ion = paste0(rep(c("b", "y"), each = n - 1L), c(i, i), suffix),
    type = rep(c("b", "y"), each = n - 1L),
    number = c(i, i),
    fragment = c(substring(peptide, 1L, i), substring(peptide, n - i + 1L, n)),
    mz = (unname(c(b, y)) + charge * SEQPROFILE_PROTON_MONO) / charge,
    intensity = 1,
    stringsAsFactors = FALSE
  )
}

#' Stick spectrum for a `seqprofile_fragment_ions()` result.
#'
#' The points carry a `text` aesthetic for plotly hover.
plot_fragment_spectrum <- function(ions, title = NULL) {
  ions$hover <- sprintf(
    "%s (%s)<br>m/z: %.4f",
    ions$ion,
    ions$fragment,
    ions$mz
  )
  # b and y labels at different heights so near-isobaric pairs stay readable
  ions$label_y <- ions$intensity + ifelse(ions$type == "b", 0.04, 0.09)
  ggplot2::ggplot(ions, ggplot2::aes(x = mz, colour = type)) +
    ggplot2::geom_segment(ggplot2::aes(xend = mz, y = 0, yend = intensity)) +
    # `text` is not a ggplot2 aesthetic; it is read by plotly::ggplotly()
    suppressWarnings(ggplot2::geom_point(
      ggplot2::aes(y = intensity, text = hover),
      size = 0.8
    )) +
    ggplot2::geom_text(
      ggplot2::aes(y = label_y, label = ion),
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(b = "#1B4965", y = "#c0392b"),
      name = "Ion series"
    ) +
    ggplot2::labs(title = title, x = "m/z", y = "Relative intensity")
}

#' Cleave one sequence with a residue-specific rule.
#'
#' @param residues,blocked_by Character vectors of one-letter codes.
#' @param terminus "C" cuts after `residues`, "N" cuts before them.
#' @param max_missed Peptides spanning up to this many uncut sites are added.
#' @return data.frame(start, end, length, sequence, missed_cleavages,
#'   mono_mass), ordered by start then end.
digest_single_sequence <- function(
  sequence,
  residues,
  terminus = "C",
  blocked_by = character(0),
  max_missed = 0L
) {
  aa <- strsplit(sequence, "")[[1]]
  n <- length(aa)
  # Bond i joins residues i and i + 1
  i <- seq_len(max(n - 1L, 0L))
  cut <- if (terminus == "C") {
    aa[i] %in% residues & !(aa[i + 1L] %in% blocked_by)
  } else {
    aa[i + 1L] %in% residues & !(aa[i] %in% blocked_by)
  }
  bonds <- i[cut]
  starts <- c(1L, bonds + 1L)
  ends <- c(bonds, n)
  k <- length(starts)

  pep <- do.call(
    rbind,
    lapply(0:min(max_missed, k - 1L), function(m) {
      idx <- seq_len(k - m)
      data.frame(
        start = starts[idx],
        end = ends[idx + m],
        missed_cleavages = rep(as.integer(m), length(idx))
      )
    })
  )
  pep <- pep[order(pep$start, pep$end), ]
  rownames(pep) <- NULL
  pep$length <- pep$end - pep$start + 1L
  pep$sequence <- substring(sequence, pep$start, pep$end)
  pep$mono_mass <- seqprofile_mono_mass(pep$sequence)
  pep[c("start", "end", "length", "sequence", "missed_cleavages", "mono_mass")]
}

#' Theoretical digestion of a cleaned sequence and its sequence coverage.
#'
#' Coverage counts the residues that lie in at least one peptide whose length
#' is within `min_len`..`max_len` (a complete digest always covers 100 % of
#' the sequence, so only the peptides in a detectable size range are counted).
#' When `mz_min` and `mz_max` are given, a peptide must also have its 2+ or
#' 3+ ion inside that m/z range.
#'
#' @param sequence Cleaned sequence, e.g. `clean_protein_sequence()$sequence`.
#' @param enzyme_id A name of `SEQPROFILE_ENZYMES`, or "custom".
#' @param params Named list: `max_missed`, `min_len`, `max_len`, `mz_min`,
#'   `mz_max` (both `NA` or absent for no m/z filter), `carbamidomethyl` (TRUE
#'   adds 57.021464 Da per Cys), `oxidation` (TRUE adds 15.994915 Da per Met),
#'   and for the custom rule `residues`, `terminus` ("C" / "N") and
#'   `blocked_by`.
#' @return Named list: enzyme_id, enzyme_name, rule, sequence, params,
#'   peptides (with `mz_1`, `mz_2`, `mz_3` for [M+zH]z+ and an `in_range`
#'   column), covered (one logical per residue) and stats.
run_seq_digest <- function(sequence, enzyme_id = "trypsin", params = list()) {
  sequence <- toupper(as.character(sequence)[1])
  if (is.na(sequence) || !nzchar(sequence)) {
    seqprofile_abort("The sequence is empty.", "seqprofile_empty")
  }

  if (identical(enzyme_id, "custom")) {
    spec <- list(
      name = "Custom rule",
      residues = params$residues,
      terminus = toupper(as.character(params$terminus)[1]),
      blocked_by = params$blocked_by
    )
  } else {
    spec <- SEQPROFILE_ENZYMES[[enzyme_id]]
    if (is.null(spec)) {
      seqprofile_abort(
        paste0("Unknown enzyme: ", enzyme_id),
        "seqprofile_bad_param"
      )
    }
  }
  residues <- seqprofile_parse_residues(spec$residues)
  blocked_by <- seqprofile_parse_residues(spec$blocked_by, "Blocking residues")
  if (length(residues) == 0) {
    seqprofile_abort(
      "Enter at least one cleavage residue for the custom rule.",
      "seqprofile_bad_param"
    )
  }
  if (is.na(spec$terminus) || !spec$terminus %in% c("C", "N")) {
    seqprofile_abort(
      "The cleavage side must be 'C' or 'N'.",
      "seqprofile_bad_param"
    )
  }

  p <- list(max_missed = 0L, min_len = 6L, max_len = 30L)
  for (nm in intersect(names(params), names(p))) {
    value <- params[[nm]]
    if (length(value) == 1) p[[nm]] <- value
  }
  whole <- function(x, min) !is.na(x) && x >= min && x == round(x)
  if (!whole(p$max_missed, 0)) {
    seqprofile_abort(
      "Missed cleavages must be a whole number of 0 or more.",
      "seqprofile_bad_param"
    )
  }
  if (!whole(p$min_len, 1) || !whole(p$max_len, 1) || p$min_len > p$max_len) {
    seqprofile_abort(
      "Peptide length limits must be whole numbers with minimum <= maximum.",
      "seqprofile_bad_param"
    )
  }
  p <- lapply(p, as.integer)
  p$carbamidomethyl <- isTRUE(params$carbamidomethyl)
  p$oxidation <- isTRUE(params$oxidation)
  one_number <- function(x) {
    if (length(x) == 1) suppressWarnings(as.numeric(x)) else NA_real_
  }
  p$mz_min <- one_number(params$mz_min)
  p$mz_max <- one_number(params$mz_max)
  use_mz <- !is.na(p$mz_min) || !is.na(p$mz_max)
  if (
    use_mz &&
      (is.na(p$mz_min) ||
        is.na(p$mz_max) ||
        p$mz_min < 0 ||
        p$mz_min >= p$mz_max)
  ) {
    seqprofile_abort(
      "Enter both m/z limits, with minimum < maximum, or switch the m/z filter off.",
      "seqprofile_bad_param"
    )
  }

  peptides <- digest_single_sequence(
    sequence,
    residues,
    spec$terminus,
    blocked_by,
    p$max_missed
  )
  peptides$mono_mass <- seqprofile_mono_mass(
    peptides$sequence,
    seqprofile_residue_masses(p$carbamidomethyl, p$oxidation)
  )
  for (z in 1:3) {
    peptides[[paste0("mz_", z)]] <- (peptides$mono_mass +
      z * SEQPROFILE_PROTON_MONO) /
      z
  }
  in_mz <- function(mz) !is.na(mz) & mz >= p$mz_min & mz <= p$mz_max
  peptides$in_range <- peptides$length >= p$min_len &
    peptides$length <= p$max_len &
    (!use_mz | in_mz(peptides$mz_2) | in_mz(peptides$mz_3))

  n <- nchar(sequence)
  kept <- peptides[peptides$in_range, ]
  depth <- cumsum(
    tabulate(kept$start, n + 1L) - tabulate(kept$end + 1L, n + 1L)
  )
  covered <- depth[seq_len(n)] > 0

  list(
    enzyme_id = enzyme_id,
    enzyme_name = spec$name,
    rule = seqprofile_enzyme_rule(residues, spec$terminus, blocked_by),
    sequence = sequence,
    params = p,
    peptides = peptides,
    covered = covered,
    stats = list(
      seq_length = n,
      n_sites = sum(peptides$missed_cleavages == 0) - 1L,
      n_peptides = nrow(peptides),
      n_in_range = sum(peptides$in_range),
      covered_residues = sum(covered),
      coverage_pct = 100 * mean(covered)
    )
  )
}
