#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(Sys.getenv("JSS_MAIN_RECOVERY_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
input_root <- if (length(args) >= 1L) args[[1L]] else file.path(root, "paper", "data", "public-derived")
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(input_root, "main-recovery")
resume <- tolower(Sys.getenv("JSS_MAIN_RECOVERY_RESUME", unset = "true")) %in% c("true", "t", "1", "yes")
conf.level <- as.numeric(Sys.getenv("JSS_MAIN_RECOVERY_CONF_LEVEL", unset = "0.95"))
if (!is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
  stop("JSS_MAIN_RECOVERY_CONF_LEVEL must be strictly between zero and one.", call. = FALSE)
}

source(file.path(root, "paper", "R", "main-recovery-evidence.R"), local = FALSE)
paths <- jss_build_main_recovery_evidence(
  input_root = input_root,
  output_dir = output_dir,
  resume = resume,
  conf.level = conf.level,
  repo_root = root
)

cat(
  if (isTRUE(attr(paths, "resumed"))) "Main-recovery evidence is current; reused checkpoint.\n" else "Main-recovery evidence regenerated.\n",
  "Output directory: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n",
  "Artifacts: ", length(paths), "\n",
  sep = ""
)
