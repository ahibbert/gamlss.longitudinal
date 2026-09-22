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
asset_locations <- gregexpr("(?:includegraphics(?:\\[[^]]*\\])?\\{|input\\{)([^}]+)", joined, perl = TRUE)[[1L]]
matches <- regmatches(joined, list(asset_locations))[[1L]]
refs <- sub("^.*\\{", "", matches)
keep <- grepl("^(charts/|tables/|orcid[.]pdf)", refs)
refs <- refs[keep]
asset_locations <- asset_locations[keep]
asset_labels <- vapply(asset_locations, function(position) {
  prefix <- substr(joined, 1L, position)
  begin <- gregexpr("\\\\begin\\{figure\\*?\\}", prefix, perl = TRUE)[[1L]]
  end <- gregexpr("\\\\end\\{figure\\*?\\}", prefix, perl = TRUE)[[1L]]
  begin <- begin[begin > 0L]; end <- end[end > 0L]
  if (!length(begin) || (length(end) && max(end) > max(begin))) return("")
  suffix <- substr(joined, position, nchar(joined))
  close <- regexpr("\\\\end\\{figure\\*?\\}", suffix, perl = TRUE)[[1L]]
  if (close < 0L) return("")
  block <- substr(joined, max(begin), position + close - 1L)
  label <- regmatches(block, regexpr("\\\\label\\{[^}]+\\}", block, perl = TRUE))
  if (!length(label)) "" else sub("\\}$", "", sub("^\\\\label\\{", "", label))
}, character(1L))
inline <- unique(sub("^.*\\\\label\\{", "", regmatches(joined, gregexpr("\\\\label\\{tab:[^}]+", joined, perl = TRUE))[[1L]]))
refs <- c(refs, paste0("inline:", inline))
asset_labels <- c(asset_labels, inline)
manifest <- utils::read.csv(file.path(root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
key <- function(x) sub("[.](png|pdf|tex)$", "", gsub("\\\\", "/", x), ignore.case = TRUE)
idx <- match(key(refs), key(manifest$manuscript_path))
label_idx <- match(asset_labels, manifest$manuscript_label)
private_label_keys <- ifelse(nzchar(asset_labels), paste0(
  "private-label-sha256:",
  vapply(asset_labels, function(label) digest::digest(
    charToRaw(label), algo = "sha256", serialize = FALSE
  ), character(1L))
), "")
label_idx[is.na(label_idx)] <- match(
  private_label_keys[is.na(label_idx)], manifest$manuscript_label
)
idx[is.na(idx)] <- label_idx[is.na(idx)]
audit <- data.frame(manuscript_path = refs, classified = !is.na(idx), artifact_id = manifest$artifact_id[idx],
  manuscript_label = asset_labels, access = manifest$access[idx],
  publication_status = manifest$publication_status[idx], actual_sha256 = NA_character_,
  reference_sha256 = manifest$approved_sha256[idx], hash_verified = NA,
  stringsAsFactors = FALSE)

inline_hash <- function(label) {
  label_line <- which(grepl(paste0("\\label{", label, "}"), tex, fixed = TRUE))
  if (length(label_line) != 1L) return(NA_character_)
  before <- seq_len(label_line)
  starts <- before[
    grepl("\\begin{table", tex[before], fixed = TRUE) |
      grepl("\\begin{longtable", tex[before], fixed = TRUE)
  ]
  after <- seq.int(label_line, length(tex))
  ends <- after[
    grepl("\\end{table", tex[after], fixed = TRUE) |
      grepl("\\end{longtable", tex[after], fixed = TRUE)
  ]
  if (!length(starts) || !length(ends)) return(NA_character_)
  digest::digest(charToRaw(paste(tex[max(starts):min(ends)], collapse = "\n")),
    algo = "sha256", serialize = FALSE)
}
private <- which(audit$classified & audit$access == "private")
if (length(private)) {
  actual <- vapply(private, function(i) {
    if (startsWith(refs[[i]], "inline:")) {
      inline_hash(sub("^inline:", "", refs[[i]]))
    } else {
      path <- file.path(repo, refs[[i]])
      if (!file.exists(path)) NA_character_ else digest::digest(file = path, algo = "sha256", serialize = FALSE)
    }
  }, character(1L))
  expected <- manifest$approved_sha256[idx[private]]
  audit$actual_sha256[private] <- actual
  audit$hash_verified[private] <- !is.na(actual) & nzchar(expected) & actual == expected
}
cat(sprintf("Audited %d manuscript assets: %d classified, %d unclassified.\n", nrow(audit), sum(audit$classified), sum(!audit$classified)))
# Do not copy withdrawn clinical labels or asset names from the private paper
# repository into the public checkout. Neutral manifest identifiers and hashes
# preserve the audit result without disclosing those strings.
if (length(private)) {
  audit$manuscript_path[private] <- manifest$manuscript_path[idx[private]]
  audit$manuscript_label[private] <- manifest$manuscript_label[idx[private]]
}
out <- file.path(root, "paper", "manuscript-asset-audit.csv")
utils::write.csv(audit, out, row.names = FALSE)
if (any(!audit$classified)) {
  stop("Unclassified manuscript assets: ", paste(audit$manuscript_path[!audit$classified], collapse = ", "), call. = FALSE)
}
retired <- audit$classified & audit$publication_status == "retired"
if (any(retired)) {
  stop(
    "Retired artifacts are still referenced by the manuscript: ",
    paste(audit$artifact_id[retired], collapse = ", "),
    call. = FALSE
  )
}
if (length(private) && any(!audit$hash_verified[private])) {
  stop("Private static manuscript artifact hash mismatch: ",
    paste(audit$artifact_id[private][!audit$hash_verified[private]], collapse = ", "), call. = FALSE)
}
