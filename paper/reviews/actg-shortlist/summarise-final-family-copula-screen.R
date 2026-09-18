#!/usr/bin/env Rscript

root <- file.path(
  "paper", "reviews", "actg-shortlist", "final-family-copula-screen"
)
checkpoint_files <- list.files(
  file.path(root, "checkpoints"), pattern = "[.]rds$", full.names = TRUE
)

objects <- lapply(checkpoint_files, readRDS)
objects <- objects[vapply(objects, function(x) {
  is.list(x) && !is.null(x$metadata$label)
}, logical(1L))]

result_rows <- do.call(rbind, lapply(objects, function(object) {
  data.frame(
    label = object$metadata$label,
    margin = object$metadata$margin,
    copula = object$metadata$copula,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    EDF = object$EDF,
    converged = object$converged,
    elapsed_sec = object$elapsed_sec,
    stop_reason = object$stop_reason,
    stringsAsFactors = FALSE
  )
}))
reference_AIC <- result_rows$AIC[result_rows$label == "zisichel-g"]
reference_BIC <- result_rows$BIC[result_rows$label == "zisichel-g"]
result_rows$delta_AIC_from_zisichel_g <- result_rows$AIC - reference_AIC
result_rows$delta_BIC_from_zisichel_g <- result_rows$BIC - reference_BIC
result_rows <- result_rows[order(result_rows$AIC), ]
write.csv(
  result_rows,
  file.path(root, "completed-family-copula-results.csv"),
  row.names = FALSE
)

diagnostics <- do.call(rbind, lapply(objects, function(object) {
  if (is.null(object$diagnostics)) return(NULL)
  transform(object$diagnostics$by_time, label = object$metadata$label)
}))
diagnostics <- diagnostics[
  order(diagnostics$label, diagnostics$time),
  c("label", "time", "n", "pit_mean", "pit_sd", "pit_cvm", "lower_10",
    "upper_10", "worm_rms", "worm_intercept", "worm_linear",
    "worm_quadratic", "worm_cubic")
]
write.csv(
  diagnostics,
  file.path(root, "completed-diagnostics-by-time.csv"),
  row.names = FALSE
)

print(result_rows, row.names = FALSE)
