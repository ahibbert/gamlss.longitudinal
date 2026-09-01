#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else if (i == length(args)) stop("Missing value after ", flag, call. = FALSE) else args[[i + 1L]] }
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
repo <- value_after("--paper-repo", "")
if (!nzchar(repo)) stop("Supply --paper-repo <clean clone>.", call. = FALSE)
repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
profile <- value_after("--profile", "paper")
if (identical(profile, "expanded")) profile <- "paper"
profile <- match.arg(profile, c("paper", "full"))
apply <- "--apply" %in% args
if (apply && "--dry-run" %in% args) stop("Choose either --dry-run or --apply.", call. = FALSE)

status <- system2("git", c("-C", shQuote(repo), "status", "--porcelain"), stdout = TRUE, stderr = TRUE)
if (length(status)) stop("Paper repository must be clean before publishing assets.", call. = FALSE)
manifest_path <- file.path(root, "results", "jss-replication", profile, "manifest.csv")
if (!file.exists(manifest_path)) stop("Run the requested replication profile first: ", manifest_path, call. = FALSE)
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
rows <- manifest$access == "public" & manifest$publication_status == "active" & nzchar(manifest$output_path)
publish <- manifest[rows, , drop = FALSE]

actions <- lapply(seq_len(nrow(publish)), function(i) {
  src <- publish$output_path[[i]]
  if (!grepl("^([A-Za-z]:[/\\\\]|/)", src)) src <- file.path(dirname(manifest_path), src)
  dest_rel <- publish$manuscript_path[[i]]
  if (identical(publish$artifact_type[[i]], "table") && !grepl("[.]tex$", dest_rel, ignore.case = TRUE)) dest_rel <- paste0(dest_rel, ".tex")
  dest_rel <- gsub("\\\\", "/", dest_rel)
  if (!grepl("^(charts|tables)/[^/].+", dest_rel) || grepl("(^|/)\\.\\.(/|$)", dest_rel)) {
    stop("Publisher destination is outside the charts/tables allowlist: ", dest_rel, call. = FALSE)
  }
  dest <- file.path(repo, dest_rel)
  if (!file.exists(src)) stop("Missing allowlisted generated asset: ", src, call. = FALSE)
  state <- if (!file.exists(dest)) "ADD" else if (identical(readBin(src, "raw", file.info(src)$size), readBin(dest, "raw", file.info(dest)$size))) "UNCHANGED" else "UPDATE"
  cat(sprintf("%-9s %s <- %s\n", state, dest_rel, src))
  if (apply && state != "UNCHANGED") {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src, dest, overwrite = TRUE)) stop("Failed to publish ", dest_rel, call. = FALSE)
  }
  data.frame(artifact_id = publish$artifact_id[[i]], action = state, destination = dest_rel, stringsAsFactors = FALSE)
})
report <- do.call(rbind, actions)
cat(if (apply) "Assets staged in the paper clone; no commit or push was performed.\n" else "Dry run only; no files were changed.\n")
invisible(report)
