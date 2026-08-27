#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[[i + 1L]] }
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
repo <- value_after("--paper-repo", "")
temporary <- !nzchar(repo)
if (temporary) {
  repo <- tempfile("working-paper-gamlss-long-")
  status <- system2("git", c("clone", "--depth", "1", "https://github.com/ahibbert/working-paper-gamlss-long.git", shQuote(repo)))
  if (!identical(status, 0L)) stop("Could not clone the canonical paper repository. Supply --paper-repo <clone>.", call. = FALSE)
}
repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
tex <- readLines(file.path(repo, "main.tex"), warn = FALSE)
joined <- paste(tex, collapse = "\n")
matches <- regmatches(joined, gregexpr("(?:includegraphics(?:\\[[^]]*\\])?\\{|input\\{)([^}]+)", joined, perl = TRUE))[[1L]]
refs <- sub("^.*\\{", "", matches)
refs <- refs[grepl("^(charts/|tables/|orcid[.]pdf|lipid - term plot[.]png)", refs)]
inline <- unique(sub("^.*\\\\label\\{", "", regmatches(joined, gregexpr("\\\\label\\{tab:[^}]+", joined, perl = TRUE))[[1L]]))
refs <- c(refs, paste0("inline:", inline))
manifest <- utils::read.csv(file.path(root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
key <- function(x) sub("[.](png|pdf|tex)$", "", gsub("\\\\", "/", x), ignore.case = TRUE)
idx <- match(key(refs), key(manifest$manuscript_path))
audit <- data.frame(manuscript_path = refs, classified = !is.na(idx), artifact_id = manifest$artifact_id[idx],
  access = manifest$access[idx], publication_status = manifest$publication_status[idx], stringsAsFactors = FALSE)
out <- file.path(root, "paper", "manuscript-asset-audit.csv")
utils::write.csv(audit, out, row.names = FALSE)
cat(sprintf("Audited %d manuscript assets: %d classified, %d unclassified.\n", nrow(audit), sum(audit$classified), sum(!audit$classified)))
if (any(!audit$classified)) stop("Unclassified manuscript assets: ", paste(audit$manuscript_path[!audit$classified], collapse = ", "), call. = FALSE)
