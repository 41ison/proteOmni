## proteOmni — dependency bootstrap
##
## Runs with BASE R ONLY. Nothing here may reference an add-on package before
## it has been installed, so this file can execute on a completely fresh R
## installation (no shiny, no BiocManager, no pak).
##
## Package versions are pinned in proteOmni.lock (next to this file). Every
## package whose installed version differs from the lock — missing, older or
## newer — is (re)installed into the user's R library at the locked version,
## from a dated Posit Package Manager snapshot. Regenerate the lock with
## tools/update_lockfile.R; never edit it by hand.
##
## Sourced by run.R before Shiny is loaded, and defensively by proteOmni.r so
## that launching the app directly from an IDE still works.

if (!isTRUE(getOption("proteOmni.bootstrapped"))) {
  # ── Locate the lockfile ──────────────────────────────────────────────────
  .omni_dir <- getwd()
  for (.i in rev(seq_len(sys.nframe()))) {
    .f <- sys.frame(.i)$ofile
    if (!is.null(.f) && basename(.f) == "bootstrap.R") {
      .omni_dir <- dirname(normalizePath(.f))
      break
    }
  }
  .omni_lock_path <- file.path(.omni_dir, "proteOmni.lock")
  if (!file.exists(.omni_lock_path)) {
    stop(
      "proteOmni.lock not found in ",
      .omni_dir,
      ". Start proteOmni from its app folder.",
      call. = FALSE
    )
  }
  .omni_lock <- read.dcf(.omni_lock_path, all = TRUE)
  .omni_hdr <- .omni_lock[!is.na(.omni_lock$Lockfile), , drop = FALSE][1, ]
  .omni_lock <- .omni_lock[!is.na(.omni_lock$Package), , drop = FALSE]
  .omni_req <- .omni_lock[.omni_lock$Optional != "yes", , drop = FALSE]
  .omni_repos <- c(
    CRAN = .omni_hdr$CRAN,
    BioCsoft = .omni_hdr$BioCsoft,
    BioCann = .omni_hdr$BioCann
  )
  # Used by on-demand installs (e.g. OrgDb packages in PwrQuant) so they come
  # from the same snapshot.
  options(
    proteOmni.repos = .omni_repos,
    proteOmni.snapshot = .omni_hdr$Snapshot
  )

  .omni_r <- paste(
    R.version$major,
    sub("\\..*$", "", R.version$minor),
    sep = "."
  )
  .omni_lock_r <- sub("^([0-9]+\\.[0-9]+).*$", "\\1", .omni_hdr$R)
  if (!identical(.omni_r, .omni_lock_r)) {
    message(
      "[proteOmni] WARNING: proteOmni.lock was built for R ",
      .omni_lock_r,
      ".x (Bioconductor ",
      .omni_hdr$Bioconductor,
      ") but you are running R ",
      .omni_r,
      ". Installation may need to compile from source or fail. ",
      "Installing R ",
      .omni_lock_r,
      ".x is recommended."
    )
  }

  # ── Writable library ─────────────────────────────────────────────────────
  # Rscript does not offer to create the personal library the way interactive R
  # does, so a fresh install has nowhere writable to put packages.
  .omni_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(.omni_lib)) {
    .omni_lib <- file.path(
      path.expand("~"),
      "R",
      paste0(R.version$platform, "-library"),
      .omni_r
    )
  }
  .omni_lib <- strsplit(.omni_lib, .Platform$path.sep)[[1]][1]
  if (!dir.exists(.omni_lib)) {
    dir.create(.omni_lib, recursive = TRUE, showWarnings = FALSE)
  }
  if (dir.exists(.omni_lib) && file.access(.omni_lib, 2) == 0) {
    .libPaths(c(.omni_lib, .libPaths()))
  }
  if (file.access(.libPaths()[1], 2) != 0) {
    stop(
      "No writable R library found. Set R_LIBS_USER to a writable directory ",
      "and restart proteOmni.",
      call. = FALSE
    )
  }

  # ── Compare installed versions with the lock ─────────────────────────────
  # The version that counts is the one R would load: the first match along
  # .libPaths().
  .omni_outdated <- function() {
    ip <- utils::installed.packages(fields = "RemoteSha", noCache = TRUE)
    ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]
    idx <- match(.omni_req$Package, ip[, "Package"])
    v <- unname(ip[idx, "Version"])
    sha <- unname(ip[idx, "RemoteSha"])
    bad <- is.na(v) |
      v != .omni_req$Version |
      (.omni_req$Source == "GitHub" &
        (is.na(sha) | sha != .omni_req$RemoteSha))
    bad[is.na(bad)] <- TRUE
    out <- .omni_req[bad, , drop = FALSE]
    out$Installed <- v[bad]
    out
  }

  .omni_todo <- .omni_outdated()

  # ── CRAN + Bioconductor (from the dated snapshot) ────────────────────────
  .omni_repo_pkgs <- .omni_todo$Package[.omni_todo$Source != "GitHub"]
  if (length(.omni_repo_pkgs)) {
    message(
      "[proteOmni] Installing ",
      length(.omni_repo_pkgs),
      " package(s) at the versions pinned in proteOmni.lock (snapshot ",
      .omni_hdr$Snapshot,
      "). This can take a while on the first launch ..."
    )
    # arrow needs these set at install time to get a full-featured libarrow
    # when it is built from source (e.g. on Linux).
    Sys.setenv(LIBARROW_MINIMAL = "false", ARROW_WITH_ZSTD = "ON")
    tryCatch(
      utils::install.packages(
        .omni_repo_pkgs,
        lib = .libPaths()[1],
        repos = .omni_repos,
        # Every dependency is itself in the lock, so it is either already at
        # the locked version or part of this same install.
        dependencies = FALSE,
        Ncpus = max(1L, min(4L, parallel::detectCores(), na.rm = TRUE))
      ),
      error = function(e) {
        message(
          "[proteOmni] Package installation failed: ",
          conditionMessage(e)
        )
      }
    )
  }

  # ── GitHub (pinned commit, installed with base R) ────────────────────────
  for (.i in which(.omni_todo$Source == "GitHub")) {
    .r <- .omni_todo[.i, ]
    message(
      "[proteOmni] Installing ",
      .r$Package,
      " from GitHub (",
      .r$RemoteUsername,
      "/",
      .r$RemoteRepo,
      "@",
      substr(.r$RemoteSha, 1, 7),
      ") ..."
    )
    tryCatch(
      {
        .tmp <- tempfile("proteOmni-gh-")
        dir.create(.tmp)
        .tar <- file.path(.tmp, "src.tar.gz")
        utils::download.file(
          sprintf(
            "https://github.com/%s/%s/archive/%s.tar.gz",
            .r$RemoteUsername,
            .r$RemoteRepo,
            .r$RemoteSha
          ),
          .tar,
          mode = "wb",
          quiet = TRUE
        )
        utils::untar(.tar, exdir = .tmp)
        .src <- setdiff(list.dirs(.tmp, recursive = FALSE), .tmp)[1]
        # Some repos commit build artefacts (diann ships src/*.o built on
        # another platform); R CMD build would strip them, so do the same.
        unlink(list.files(
          file.path(.src, "src"),
          pattern = "\\.(o|so|dll|dylib)$",
          full.names = TRUE
        ))
        # Record the commit in DESCRIPTION so the next launch can verify it.
        .desc_path <- file.path(.src, "DESCRIPTION")
        .desc <- read.dcf(.desc_path)
        .remote <- c(
          RemoteType = "github",
          RemoteHost = "api.github.com",
          RemoteUsername = .r$RemoteUsername,
          RemoteRepo = .r$RemoteRepo,
          RemoteSha = .r$RemoteSha
        )
        .desc <- .desc[, setdiff(colnames(.desc), names(.remote)), drop = FALSE]
        .desc <- cbind(.desc, t(.remote))
        write.dcf(.desc, .desc_path)
        utils::install.packages(
          .src,
          lib = .libPaths()[1],
          repos = NULL,
          type = "source"
        )
      },
      error = function(e) {
        message(
          "[proteOmni] Failed to install ",
          .r$Package,
          ": ",
          conditionMessage(e)
        )
      },
      finally = unlink(.tmp, recursive = TRUE)
    )
  }

  # ── Verify and report ────────────────────────────────────────────────────
  if (nrow(.omni_todo)) {
    .omni_todo <- .omni_outdated()
  }
  if (nrow(.omni_todo)) {
    message(
      "\n[proteOmni] These packages do not match proteOmni.lock:\n",
      paste0(
        "  - ",
        .omni_todo$Package,
        ": installed ",
        ifelse(is.na(.omni_todo$Installed), "none", .omni_todo$Installed),
        ", locked ",
        .omni_todo$Version,
        collapse = "\n"
      ),
      "\nproteOmni may fail to start or behave differently. See the messages ",
      "above for the cause (commonly no internet connection, a missing ",
      "compiler/system library, or - on Windows - the package being in use by ",
      "another open R session).\n"
    )
  } else {
    message(
      "[proteOmni] All dependencies match proteOmni.lock (snapshot ",
      .omni_hdr$Snapshot,
      ")."
    )
  }

  rm(
    list = intersect(
      ls(all.names = TRUE),
      c(
        ".omni_dir",
        ".omni_lock_path",
        ".omni_lock",
        ".omni_hdr",
        ".omni_req",
        ".omni_repos",
        ".omni_r",
        ".omni_lock_r",
        ".omni_lib",
        ".omni_outdated",
        ".omni_todo",
        ".omni_repo_pkgs",
        ".i",
        ".f",
        ".r",
        ".tmp",
        ".tar",
        ".src",
        ".desc_path",
        ".desc",
        ".remote"
      )
    )
  )
  options(proteOmni.bootstrapped = TRUE)
}
