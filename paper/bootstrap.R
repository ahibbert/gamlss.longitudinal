jss_bootstrap <- function(root, restore = TRUE, install_source = TRUE) {
  lock <- file.path(root, "paper", "renv.lock")
  version <- paste(R.version$major, strsplit(R.version$minor, "[.]", fixed = FALSE)[[1L]][1L], sep = ".")
  library_override <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_LIBRARY", unset = "")
  lock_key <- substr(unname(tools::md5sum(lock)), 1L, 12L)
  lib <- if (nzchar(library_override)) {
    library_override
  } else {
    file.path(
      tools::R_user_dir("gamlss.longitudinal", which = "cache"),
      "paper-replication", paste0("lock-", lock_key), paste0("R-", version)
    )
  }
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  inherited_libs <- .libPaths()
  .libPaths(unique(c(lib, inherited_libs)))
  # A restored reviewer run sees only the pinned paper library. The explicit
  # --no-restore developer shortcut may reuse already installed dependencies,
  # including inside the child R process used by R CMD INSTALL.
  user_libs <- if (restore) lib else paste(unique(c(lib, inherited_libs)), collapse = .Platform$path.sep)
  Sys.setenv(R_LIBS_USER = user_libs)

  if (restore) {
    if (!requireNamespace("renv", quietly = TRUE)) {
      install.packages("renv", lib = lib, repos = "https://cloud.r-project.org")
    }
    renv::restore(project = root, lockfile = lock, library = lib, prompt = FALSE)
  }
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Dependency bootstrap did not provide 'targets'. Rerun without --no-restore.", call. = FALSE)
  }

  installed <- FALSE
  if (install_source) {
    cmd <- file.path(R.home("bin"), "R")
    status <- system2(cmd, c("CMD", "INSTALL", "--no-multiarch", "--with-keep.source", "-l", shQuote(lib), shQuote(root)))
    if (!identical(status, 0L)) stop("Failed to install the checked-out package source.", call. = FALSE)
    installed <- TRUE
  } else if (!requireNamespace("gamlss.longitudinal", quietly = TRUE)) {
    stop("--no-install requires an available gamlss.longitudinal package.", call. = FALSE)
  }
  list(library = normalizePath(lib, winslash = "/", mustWork = TRUE), restored = restore, source_installed = installed)
}
