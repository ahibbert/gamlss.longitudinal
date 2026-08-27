#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  at <- match(flag, args)
  if (is.na(at)) return(default)
  if (at == length(args) || startsWith(args[[at + 1L]], "--")) stop("Missing value after ", flag, call. = FALSE)
  args[[at + 1L]]
}

profile <- value_after("--profile", Sys.getenv("GAMLSS_LONGITUDINAL_JSS_PROFILE", "smoke"))
if (identical(profile, "expanded")) {
  warning("Profile 'expanded' is deprecated; using 'paper'.", call. = FALSE)
  profile <- "paper"
}
profile <- match.arg(profile, c("smoke", "paper", "full"))
workers <- as.integer(value_after("--workers", Sys.getenv("GAMLSS_LONGITUDINAL_JSS_WORKERS", "1")))
if (is.na(workers) || workers < 1L) stop("--workers must be a positive integer.", call. = FALSE)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else file.path("paper", "replicate.R")
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
setwd(root)
source(file.path("paper", "bootstrap.R"))
bootstrap <- jss_bootstrap(root, restore = !"--no-restore" %in% args, install_source = !"--no-install" %in% args)

Sys.setenv(
  GAMLSS_LONGITUDINAL_JSS_PROFILE = profile,
  GAMLSS_LONGITUDINAL_JSS_ROOT = root,
  GAMLSS_LONGITUDINAL_JSS_WORKERS = workers
)
store <- file.path("paper", "_targets", profile)
message("Running public JSS replication profile '", profile, "' with ", workers, " worker(s).")
started <- Sys.time()
if (workers > 1L && "tar_make_future" %in% getNamespaceExports("targets")) {
  targets::tar_make_future(script = file.path("paper", "_targets.R"), store = store, workers = workers)
} else {
  targets::tar_make(script = file.path("paper", "_targets.R"), store = store)
}

source(file.path("paper", "R", "replication-helpers.R"))
settings <- jss_settings(create = TRUE)
jss_write_run_metadata(settings, started, bootstrap, store)
jss_check_target_warnings(store, strict = identical(tolower(Sys.getenv("GAMLSS_LONGITUDINAL_JSS_STRICT_WARNINGS")), "true"))
invisible(jss_write_provenance_hashes(settings))
message("Replication complete. Outputs: ", settings$out_dir)
