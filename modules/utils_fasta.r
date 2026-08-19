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
  prot_ids <- sapply(strsplit(headers, "\\s+|\\|"), function(x) {
    if (length(x) >= 2 && grepl("sp|tr", x[1])) {
      return(x[2])
    }
    return(x[1])
  })

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
  pk_base <- c(H = 6.00, K = 9.74, R = 12.48) # BH+ -> B + H+

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
