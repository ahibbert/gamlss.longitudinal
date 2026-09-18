#!/usr/bin/env Rscript

root <- file.path("paper", "reviews", "actg-shortlist", "time1-polish-screen")
checkpoint_files <- list.files(file.path(root, "checkpoints"), full.names = TRUE)

diagnostics <- do.call(rbind, lapply(checkpoint_files, function(path) {
  object <- readRDS(path)
  if (is.null(object$diagnostics)) return(NULL)
  transform(object$diagnostics$by_time, label = object$metadata$label)
}))

diagnostics <- diagnostics[
  order(diagnostics$label, diagnostics$time),
  c("label", "time", "n", "pit_mean", "pit_sd", "pit_cvm", "worm_rms",
    "worm_intercept", "worm_linear", "worm_quadratic", "worm_cubic")
]
write.csv(
  diagnostics,
  file.path(root, "all-candidate-diagnostics-by-time.csv"),
  row.names = FALSE
)

print(diagnostics[, c("label", "time", "pit_sd", "pit_cvm", "worm_rms")],
      row.names = FALSE)
