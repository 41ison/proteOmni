## proteOmni — regenerate proteOmni.lock (MAINTAINERS ONLY)
##
## End users never run this; bootstrap.R reads the lockfile it produces.
##
## Usage, from the app/ directory:
##
##   Rscript tools/update_lockfile.R <snapshot-date> [bioconductor-version]
##   Rscript tools/update_lockfile.R 2026-08-11 3.23
##
## Versions are resolved from a dated Posit Public Package Manager (P3M)
## snapshot, so every locked version is guaranteed to be downloadable — on
## Windows and macOS as pre-built binaries. GitHub-only packages (none at
## present, see `github_direct`) are pinned to the commit currently installed
## in your library.
##
## After writing the lock, the script compares it with your own library. Before
## committing a new lock, run `source("bootstrap.R")` in a fresh R session (it
## will bring your library in line with the lock) and test the app.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript tools/update_lockfile.R <YYYY-MM-DD> [bioc-version]")
}
snapshot <- as.character(as.Date(args[[1]]))
bioc_version <- if (length(args) >= 2) {
  args[[2]]
} else {
  as.character(BiocManager::version())
}

# ── Direct dependencies ─────────────────────────────────────────────────────
# Edit these lists when proteOmni starts or stops using a package. Indirect
# dependencies are resolved automatically.
cran_bioc_direct <- c(
  # CRAN
  "shiny",
  "shinydashboard",
  "shinyjs",
  "fresh",
  "BiocManager",
  "tidyverse",
  "tidytext",
  "janitor",
  "data.table",
  "ggpointdensity",
  "ggtext",
  "ggrepel",
  "ggseqlogo",
  "GGally",
  "ggsci",
  "ggfortify",
  "plotly",
  "viridis",
  "RColorBrewer",
  "scales",
  "patchwork",
  "seqinr",
  "zip",
  "DT",
  "colourpicker",
  "R6",
  "gridExtra",
  "lavaan",
  "naniar",
  "missForest",
  "httr",
  "jsonlite",
  "arrow",
  # Bioconductor
  "limma",
  "edgeR",
  "Biostrings",
  "sva",
  "impute",
  "pcaMethods",
  "ComplexHeatmap",
  "STRINGdb",
  "clusterProfiler",
  "enrichplot",
  "AnnotationDbi",
  "mixOmics",
  "org.Hs.eg.db"
)
# GitHub-only packages, as `pkg = c(user = "owner", repo = "repository")`.
# None at present.
github_direct <- list()

# OrgDb packages installed on demand by the PwrQuant enrichment module. They
# are locked (so on-demand installs are pinned too) but marked Optional, so
# bootstrap.R does not install them up front.
pwrquant_src <- readLines("modules/mod_PwrQuant.r", warn = FALSE)
optional <- unique(gsub(
  '"',
  "",
  unlist(regmatches(
    pwrquant_src,
    gregexpr('"org\\.[A-Za-z0-9]+(\\.[A-Za-z0-9]+)*\\.db"', pwrquant_src)
  ))
))
optional <- setdiff(optional, cran_bioc_direct)

# ── Snapshot repositories ───────────────────────────────────────────────────
p3m <- "https://packagemanager.posit.co"
repos <- c(
  CRAN = sprintf("%s/cran/%s", p3m, snapshot),
  BioCsoft = sprintf(
    "%s/bioconductor/%s/packages/%s/bioc",
    p3m,
    snapshot,
    bioc_version
  ),
  BioCann = sprintf(
    "%s/bioconductor/%s/packages/%s/data/annotation",
    p3m,
    snapshot,
    bioc_version
  )
)
db <- available.packages(repos = repos, type = "source")
db_repo <- vapply(
  db[, "Repository"],
  function(u) names(repos)[startsWith(u, repos)][1],
  character(1)
)

# Base and recommended packages ship with R itself and are not pinned.
r_pkgs <- c(
  rownames(installed.packages(.Library, priority = "base")),
  rownames(installed.packages(.Library, priority = "recommended"))
)

# ── GitHub packages: pin the installed commit ───────────────────────────────
dep_names <- function(x) {
  x <- unlist(strsplit(x[!is.na(x)], ","))
  x <- trimws(sub("\\(.*$", "", x))
  x[nzchar(x) & x != "R"]
}
gh_records <- list()
gh_deps <- character(0)
for (pkg in names(github_direct)) {
  gh <- github_direct[[pkg]]
  sha <- tryCatch(
    utils::packageDescription(pkg, fields = "RemoteSha"),
    warning = function(w) NA
  )
  if (is.na(sha)) {
    stop(
      pkg,
      " is not installed from GitHub; install it first from ",
      "https://github.com/",
      gh[["user"]],
      "/",
      gh[["repo"]]
    )
  }
  desc <- read.dcf(url(sprintf(
    "https://raw.githubusercontent.com/%s/%s/%s/DESCRIPTION",
    gh[["user"]],
    gh[["repo"]],
    sha
  )))
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(desc))
  gh_deps <- c(gh_deps, dep_names(desc[1, fields]))
  gh_records[[pkg]] <- data.frame(
    Package = pkg,
    Version = desc[1, "Version"],
    Source = "GitHub",
    Repository = NA,
    RemoteUsername = gh[["user"]],
    RemoteRepo = gh[["repo"]],
    RemoteSha = sha,
    Direct = "yes",
    Optional = "no"
  )
}

# ── Resolve the dependency closure ──────────────────────────────────────────
closure <- function(roots) {
  deps <- tools::package_dependencies(
    roots,
    db = db,
    which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
  )
  setdiff(unique(c(roots, unlist(deps))), r_pkgs)
}
required_roots <- setdiff(c(cran_bioc_direct, gh_deps), r_pkgs)
absent <- setdiff(required_roots, rownames(db))
if (length(absent)) {
  stop(
    "Not available in the ",
    snapshot,
    " snapshot: ",
    paste(absent, collapse = ", ")
  )
}
absent_opt <- setdiff(optional, rownames(db))
if (length(absent_opt)) {
  warning(
    "Optional packages not in snapshot (left unpinned): ",
    paste(absent_opt, collapse = ", ")
  )
  optional <- setdiff(optional, absent_opt)
}
required <- closure(required_roots)
optional_only <- setdiff(closure(optional), required)

all_pkgs <- c(required, optional_only)
repo_records <- data.frame(
  Package = all_pkgs,
  Version = db[all_pkgs, "Version"],
  Source = ifelse(db_repo[all_pkgs] == "CRAN", "CRAN", "Bioconductor"),
  Repository = db_repo[all_pkgs],
  RemoteUsername = NA,
  RemoteRepo = NA,
  RemoteSha = NA,
  Direct = ifelse(all_pkgs %in% cran_bioc_direct, "yes", "no"),
  Optional = ifelse(all_pkgs %in% optional_only, "yes", "no")
)
records <- rbind(repo_records, do.call(rbind, gh_records))
records <- records[order(tolower(records$Package)), ]
rownames(records) <- NULL

# ── Write the lockfile ──────────────────────────────────────────────────────
header <- data.frame(
  Lockfile = "proteOmni",
  LockfileVersion = "1",
  Generated = format(Sys.Date()),
  Note = "Generated by tools/update_lockfile.R. Do not edit by hand.",
  R = paste(R.version$major, R.version$minor, sep = "."),
  Bioconductor = bioc_version,
  Snapshot = snapshot,
  CRAN = repos[["CRAN"]],
  BioCsoft = repos[["BioCsoft"]],
  BioCann = repos[["BioCann"]]
)
lock_path <- "proteOmni.lock"
write.dcf(header, lock_path)
cat("\n", file = lock_path, append = TRUE)
write.dcf(records, lock_path, append = TRUE)

message(sprintf(
  "Wrote %s: %d required + %d optional packages (snapshot %s, Bioconductor %s).",
  lock_path,
  sum(records$Optional == "no"),
  sum(records$Optional == "yes"),
  snapshot,
  bioc_version
))

# ── Compare with the local library ──────────────────────────────────────────
ip <- installed.packages(fields = "RemoteSha")
ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]
req <- records[records$Optional == "no", ]
local_v <- unname(ip[match(req$Package, ip[, "Package"]), "Version"])
local_sha <- unname(ip[match(req$Package, ip[, "Package"]), "RemoteSha"])
diff <- is.na(local_v) |
  local_v != req$Version |
  (req$Source == "GitHub" & (is.na(local_sha) | local_sha != req$RemoteSha))
diff[is.na(diff)] <- TRUE
if (any(diff)) {
  message(
    "\nYour library differs from the new lock for ",
    sum(diff),
    " package(s). Run source('bootstrap.R') in a fresh R session and re-test ",
    "the app before committing:"
  )
  print(data.frame(
    Package = req$Package[diff],
    local = local_v[diff],
    locked = req$Version[diff]
  ))
} else {
  message("Local library matches the lockfile.")
}
