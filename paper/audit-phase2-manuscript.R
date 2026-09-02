#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing required argument: ", flag, call. = FALSE)
  args[[i + 1L]]
}

paper_repo <- value_after("--paper-repo")
main_tex <- file.path(paper_repo, "main.tex")
if (!file.exists(main_tex)) stop("Canonical paper repository has no main.tex: ", paper_repo, call. = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else file.path("paper", "audit-phase2-manuscript.R")
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "paper", "R", "replication-helpers.R"))
source(file.path(root, "paper", "R", "phase2-paper-evidence.R"))

jss_phase2_validate_canonical_source(root, paper_repo)
claims <- jss_phase2_claim_register(root)
marker_audit <- jss_phase2_validate_claim_markers(main_tex, claims, allow_pending = TRUE)
scope_audit <- jss_phase2_validate_manuscript_scope(root, main_tex, allow_pending = TRUE)
asset_audit <- jss_phase2_validate_manuscript_assets(root, main_tex, allow_pending = TRUE)
known <- jss_phase2_known_canonical_issues(root, main_tex)
prose_error <- tryCatch({ jss_phase2_lint_empirical_prose(root, main_tex, allow_pending = TRUE); "" },
  error = function(e) conditionMessage(e))

cat("Canonical Phase 2 manuscript audit\n")
cat("  main.tex: ", normalizePath(main_tex, winslash = "/"), "\n", sep = "")
cat("  pending claim markers: ", sum(marker_audit$manuscript_state == "pending_overleaf"), "\n", sep = "")
cat("  pending/missing scope anchors: ", sum(scope_audit$publication_state == "pending_overleaf" & !scope_audit$anchor_present), "\n", sep = "")
cat("  pending artifacts: ", sum(asset_audit$canonical_state == "pending"), "\n", sep = "")
active_failures <- known$status == "fail" & known$excerpt_present
cat("  known stale empirical blocks present: ", sum(active_failures), "\n", sep = "")
cat("  unregistered empirical prose: ", if (nzchar(prose_error)) "present" else "none", "\n", sep = "")
if (any(active_failures)) {
  for (i in which(active_failures)) cat("    - ", known$issue_id[[i]], ": ", known$required_action[[i]], "\n", sep = "")
}
if (any(active_failures) || nzchar(prose_error) || any(marker_audit$manuscript_state %in% c("pending_overleaf", "pending_evidence")) ||
    any(scope_audit$publication_state == "pending_overleaf") || any(asset_audit$canonical_state == "pending")) {
  quit(save = "no", status = 1L)
}
