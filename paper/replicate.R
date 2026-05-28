#!/usr/bin/env Rscript

profile <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_PROFILE", unset = "smoke")
profile <- match.arg(profile, c("smoke", "expanded"))

args_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
script_dir <- if (!is.na(args_file)) dirname(args_file) else file.path(getwd(), "paper")
root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(root)

if (!requireNamespace("targets", quietly = TRUE)) {
  stop(
    "The full JSS replication workflow requires the optional package 'targets'. ",
    "Install it with install.packages('targets') and rerun paper/replicate.R.",
    call. = FALSE
  )
}

Sys.setenv(GAMLSS_LONGITUDINAL_JSS_PROFILE = profile)
Sys.setenv(GAMLSS_LONGITUDINAL_JSS_ROOT = root)

message("Running JSS replication profile: ", profile)
targets::tar_make(script = file.path("paper", "_targets.R"))
