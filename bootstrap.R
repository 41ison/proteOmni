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
## This file stops with an error, before the app is loaded, when
##   * the running R minor version is not the one the lock was built for
##     (override: environment variable PROTEOMNI_ALLOW_R_MISMATCH=1), or
##   * a required package is still not installed after the retries.
## A package that is installed at a version other than the locked one only
## produces a warning message.
##
## Sourced by run.R before Shiny is loaded, and defensively by proteOmni.r so
## that launching the app directly from an IDE still works.

if (!isTRUE(getOption("proteOmni.bootstrapped"))) {
  local({
    # ── Locate the lockfile ────────────────────────────────────────────────
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
    # write.dcf() drops fields that are NA in every record, so a lock without
    # GitHub packages has no Remote* columns. They must exist: comparing
    # against a NULL column yields a zero-length result, which made every
    # package look up to date.
    for (.f in c("RemoteUsername", "RemoteRepo", "RemoteSha")) {
      if (is.null(.omni_lock[[.f]])) {
        .omni_lock[[.f]] <- NA_character_
      }
    }
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

    # ── R version ──────────────────────────────────────────────────────────
    # Bioconductor releases are tied to an R minor version: the snapshot has
    # no Bioconductor binaries for any other R, so a mismatch means compiling
    # dozens of packages from source against an R they were not built for.
    .omni_r <- paste(
      R.version$major,
      sub("\\..*$", "", R.version$minor),
      sep = "."
    )
    .omni_lock_r <- sub("^([0-9]+\\.[0-9]+).*$", "\\1", .omni_hdr$R)
    if (!identical(.omni_r, .omni_lock_r)) {
      .omni_r_msg <- paste0(
        "proteOmni.lock was built for R ",
        .omni_lock_r,
        ".x (Bioconductor ",
        .omni_hdr$Bioconductor,
        ") but you are running R ",
        .omni_r,
        "."
      )
      .omni_allow <- tolower(Sys.getenv("PROTEOMNI_ALLOW_R_MISMATCH"))
      if (.omni_allow %in% c("1", "true", "yes")) {
        message(
          "[proteOmni] WARNING: ",
          .omni_r_msg,
          " Continuing because PROTEOMNI_ALLOW_R_MISMATCH is set; packages ",
          "may have to compile from source or fail to install."
        )
      } else {
        stop(
          .omni_r_msg,
          "\nInstall R ",
          .omni_lock_r,
          ".x from https://cran.r-project.org and start proteOmni again.",
          "\n(Unsupported: set the environment variable ",
          "PROTEOMNI_ALLOW_R_MISMATCH=1 to try anyway.)",
          call. = FALSE
        )
      }
    }

    # ── Writable library ───────────────────────────────────────────────────
    # Rscript does not offer to create the personal library the way
    # interactive R does, so a fresh install has nowhere writable to put
    # packages.
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

    # ── Compare installed versions with the lock ───────────────────────────
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
      stopifnot(length(bad) == nrow(.omni_req))
      out <- .omni_req[bad, , drop = FALSE]
      out$Installed <- v[bad]
      out
    }

    # install.packages() reports most failures (download errors, compilation
    # failures, "not available") as warnings, which Rscript only prints when
    # the script ends. Show them as they happen and return them so the final
    # report can say why a package is missing.
    .omni_install <- function(pkgs) {
      why <- character(0)
      note <- function(cond) {
        # First line only: the rest is R's generic advice about CRAN.
        msg <- sub("\n.*$", "", trimws(conditionMessage(cond)))
        why <<- c(why, msg)
        message("[proteOmni] ", msg)
      }
      withCallingHandlers(
        tryCatch(
          utils::install.packages(
            pkgs,
            lib = .libPaths()[1],
            repos = .omni_repos,
            # Every dependency is itself in the lock, so it is either already
            # at the locked version or part of this same bootstrap.
            dependencies = FALSE,
            Ncpus = max(1L, min(4L, parallel::detectCores(), na.rm = TRUE))
          ),
          error = note
        ),
        warning = function(w) {
          note(w)
          invokeRestart("muffleWarning")
        }
      )
      unique(why)
    }

    .omni_todo <- .omni_outdated()
    .omni_why <- list()

    # ── CRAN + Bioconductor (from the dated snapshot) ──────────────────────
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
      # R's default of 60 s per file is too short for the larger downloads
      # (arrow, BH, GO.db, org.Hs.eg.db) on a slow connection. Left raised so
      # the on-demand OrgDb installs in PwrQuant benefit too.
      options(timeout = max(600, getOption("timeout"), na.rm = TRUE))

      .omni_install(.omni_repo_pkgs)

      # Retry whatever is still wrong one package at a time, so a single
      # failure cannot take others down with it. Repeat while a pass makes
      # progress: a source package can only build once its dependencies have.
      .omni_n_left <- length(.omni_repo_pkgs)
      for (.pass in 1:5) {
        .omni_left <- .omni_outdated()
        .omni_left <- .omni_left$Package[.omni_left$Source != "GitHub"]
        if (!length(.omni_left)) {
          break
        }
        # A large batch in which nothing at all succeeded is a systemic
        # failure (no network, blocked host); retrying each package would only
        # repeat it hundreds of times.
        .omni_stuck <- if (.pass == 1) {
          length(.omni_left) == .omni_n_left && .omni_n_left > 10
        } else {
          length(.omni_left) >= .omni_n_left
        }
        if (.omni_stuck) {
          break
        }
        .omni_n_left <- length(.omni_left)
        message(
          "[proteOmni] Retrying ",
          length(.omni_left),
          " package(s) one at a time: ",
          paste(.omni_left, collapse = ", ")
        )
        for (.p in .omni_left) {
          .omni_why[[.p]] <- .omni_install(.p)
        }
      }
    }

    # ── GitHub (pinned commit, installed with base R) ──────────────────────
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
          # Some repos commit build artefacts built on another platform;
          # R CMD build would strip them, so do the same.
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
          .desc <- .desc[,
            setdiff(colnames(.desc), names(.remote)),
            drop = FALSE
          ]
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
          .omni_why[[.r$Package]] <<- conditionMessage(e)
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

    # ── Verify and report ──────────────────────────────────────────────────
    if (nrow(.omni_todo)) {
      .omni_todo <- .omni_outdated()
    }
    .omni_absent <- .omni_todo[is.na(.omni_todo$Installed), , drop = FALSE]
    .omni_differ <- .omni_todo[!is.na(.omni_todo$Installed), , drop = FALSE]

    if (nrow(.omni_differ)) {
      message(
        "\n[proteOmni] WARNING: these packages are installed at a version ",
        "other than the one in proteOmni.lock:\n",
        paste0(
          "  - ",
          .omni_differ$Package,
          ": installed ",
          .omni_differ$Installed,
          ", locked ",
          .omni_differ$Version,
          collapse = "\n"
        ),
        "\nproteOmni will start but may behave differently. On Windows this ",
        "usually means the package is in use by another open R session: ",
        "close it and start proteOmni again.\n"
      )
    }

    if (nrow(.omni_absent)) {
      .omni_reason <- vapply(
        .omni_absent$Package,
        function(p) {
          w <- .omni_why[[p]]
          if (length(w)) paste(w, collapse = "; ") else "no error reported"
        },
        character(1)
      )
      stop(
        "proteOmni cannot start: ",
        nrow(.omni_absent),
        " required package(s) could not be installed.\n",
        paste0(
          "  - ",
          .omni_absent$Package,
          " ",
          .omni_absent$Version,
          ": ",
          .omni_reason,
          collapse = "\n"
        ),
        "\nCommon causes: no or unstable internet connection, a proxy or ",
        "firewall blocking packagemanager.posit.co, or a missing ",
        "compiler/system library (Linux; Rtools on Windows).",
        "\nFix the cause and start proteOmni again; only the packages still ",
        "missing are installed.",
        call. = FALSE
      )
    }

    if (!nrow(.omni_differ)) {
      message(
        "[proteOmni] All dependencies match proteOmni.lock (snapshot ",
        .omni_hdr$Snapshot,
        ")."
      )
    }

    options(proteOmni.bootstrapped = TRUE)
  })
}
