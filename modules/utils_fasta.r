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
  if(length(header_idx) == 0) return(character(0))
  
  headers <- sub("^>", "", lines[header_idx])
  prot_ids <- sapply(strsplit(headers, "\\s+|\\|"), function(x) {
     if (length(x) >= 2 && grepl("sp|tr", x[1])) return(x[2])
     return(x[1])
  })
  
  seqs <- character(length(header_idx))
  starts <- header_idx + 1
  ends <- c(header_idx[-1] - 1, length(lines))
  
  for(i in seq_along(header_idx)) {
    if(starts[i] <= ends[i]) {
      seqs[i] <- paste(lines[starts[i]:ends[i]], collapse="")
    }
  }
  names(seqs) <- prot_ids
  return(seqs)
}

# In silico digestion
in_silico_digest <- function(fasta_seqs, max_missed = 2) {
  pep_list <- lapply(fasta_seqs, function(seq) {
    if(is.na(seq) || nchar(seq) == 0) return(character(0))
    cleaved <- strsplit(seq, "(?<=[KR])(?!P)", perl = TRUE)[[1]]
    
    peps <- cleaved
    if(max_missed > 0 && length(cleaved) > 1) {
      for(m in 1:max_missed) {
        if(length(cleaved) > m) {
          for(i in 1:(length(cleaved)-m)) {
            peps <- c(peps, paste(cleaved[i:(i+m)], collapse=""))
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
      mapped_proteins = paste(unique(protein), collapse=";"),
      n_proteins = n_distinct(protein),
      .groups = "drop"
    )
    
  out <- data.frame(peptide = unique(detected_peptides), stringsAsFactors=FALSE)
  out <- left_join(out, mapping, by="peptide")
  out$classification <- ifelse(is.na(out$n_proteins) | out$n_proteins == 0, "Unmapped",
                        ifelse(out$n_proteins == 1, "Proteotypic", "Shared"))
  out
}

# Constructing the PwrQuant-compatible matrix based strictly on proteotypic sum/top3
compute_protein_abundance <- function(df, pep_col, abundance_col, class_col, mapped_prot_col, method="sum") {
  req_cols <- c("filename", pep_col, abundance_col, class_col, mapped_prot_col)
  
  # Ensure all necessary columns exist (filename defaults to unknown if missing)
  if (!"filename" %in% names(df)) df$filename <- "Unknown"
  missing_cols <- setdiff(req_cols, names(df))
  if(length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse=", ")))
  
  df_sub <- df |> 
    filter(!!sym(class_col) == "Proteotypic") |>
    filter(!is.na(!!sym(abundance_col)))
    
  if(nrow(df_sub) == 0) return(NULL)
  
  # Top3 or Sum calculation
  if(method == "sum") {
    mat <- df_sub |>
      group_by(!!sym(mapped_prot_col), filename) |>
      summarise(Abundance = sum(!!sym(abundance_col), na.rm=TRUE), .groups="drop") |>
      pivot_wider(names_from = filename, values_from = Abundance, values_fill = NA)
  } else if(method == "top3") {
    mat <- df_sub |>
      group_by(!!sym(mapped_prot_col), filename) |>
      slice_max(order_by = !!sym(abundance_col), n = 3, with_ties = FALSE) |>
      summarise(Abundance = sum(!!sym(abundance_col), na.rm=TRUE), .groups="drop") |>
      pivot_wider(names_from = filename, values_from = Abundance, values_fill = NA)
  }
  
  return(mat)
}

# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL FASTA MODULE
# ═══════════════════════════════════════════════════════════════════════════════

Fasta_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$hr(style="border-color:#2d3741;margin:6px 0;"),
    tags$div(style="padding:4px 16px;color:#adb5bd;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:1px;",
             icon("dna"), " Global FASTA Mapping"),
    fileInput(ns("fasta_file"), "Protein FASTA File", accept = c(".fasta", ".fa")),
    numericInput(ns("missed_cleavages"), "Max missed cleavages", value = 2, min = 0, max = 5)
  )
}

Fasta_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactive({
      if (!isTruthy(input$fasta_file)) return(NULL)
      
      showNotification("Reading and digesting FASTA globally...", id="fasta_global_notif", duration=NULL)
      seqs <- tryCatch(read_fasta_custom(input$fasta_file$datapath), error = function(e){
        showNotification("Failed to read FASTA", type="error")
        NULL
      })
      if(is.null(seqs)) {
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
  hi <- c(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,
          I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
  mean(sapply(strsplit(sequence,NULL)[[1]], function(aa) hi[aa]), na.rm=TRUE)
}

calculate_pI <- function(sequence) {
  pk <- c(A=2.34,R=12.48,N=10.76,D=3.86,C=8.33,Q=10.76,E=4.25,G=2.34,H=6.00,
          I=6.04,L=6.04,K=9.74,M=5.74,F=5.48,P=1.99,S=2.21,T=2.15,W=9.39,Y=10.07,V=6.02)
  mean(sapply(strsplit(sequence,NULL)[[1]], function(aa) pk[aa]), na.rm=TRUE)
}
