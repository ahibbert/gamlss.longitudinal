#!/usr/bin/env Rscript

script_path <- normalizePath(
  file.path("paper", "reviews", "actg-shortlist", "build-final-family-screen-report.R"),
  mustWork = TRUE
)
review_root <- dirname(script_path)
package_root <- normalizePath(file.path(review_root, "..", "..", ".."))
results_root <- file.path(review_root, "final-family-copula-screen")
figure_dir <- file.path(results_root, "figures")
object_dir <- file.path(results_root, "objects")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

devtools::load_all(package_root, quiet = TRUE)

checkpoint_paths <- c(
  `ZISICHEL--Gumbel` = file.path(results_root, "checkpoints", "reference-zisichel-g.rds"),
  `ZIBNB--Gumbel` = file.path(results_root, "checkpoints", "zibnb-g-74b87e02c9e4.rds"),
  `ZIPIG--Gumbel` = file.path(results_root, "checkpoints", "zipig-g-f27eacee6cb8.rds")
)
if (!all(file.exists(checkpoint_paths))) stop("Required Gumbel checkpoints are missing.")
objects <- lapply(checkpoint_paths, readRDS)

slug <- function(x) tolower(gsub("[^A-Za-z0-9]+", "-", x))
for (label in names(objects)) {
  object <- objects[[label]]
  fit <- object$fit
  name <- slug(label)
  saveRDS(fit, file.path(object_dir, paste0(name, ".rds")), compress = "gzip")

  plots <- list(
    pit = pithist(fit, randomize = TRUE, seed = 20260916L, by_time = TRUE),
    worm = wormplot(fit, randomize = TRUE, seed = 20260916L, by_time = TRUE),
    rootogram = rootogram(fit, by_time = TRUE)
  )
  for (plot_name in names(plots)) {
    ggplot2::ggsave(
      file.path(figure_dir, paste0(name, "-", plot_name, ".png")),
      plots[[plot_name]], width = 14, height = 10, dpi = 160
    )
  }
}

results <- read.csv(file.path(results_root, "completed-family-copula-results.csv"))
diagnostics <- read.csv(file.path(results_root, "completed-diagnostics-by-time.csv"))
time1 <- diagnostics[diagnostics$time == 1, c(
  "label", "pit_sd", "pit_cvm", "worm_rms"
)]
comparison <- merge(results, time1, by = "label", all.x = TRUE)
comparison <- comparison[order(comparison$AIC), ]

display_label <- c(
  `zisichel-g` = "ZISICHEL--Gumbel",
  `zibnb-g` = "ZIBNB--Gumbel",
  `zipig-g` = "ZIPIG--Gumbel"
)
comparison$model <- unname(display_label[comparison$label])

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
rows <- paste(vapply(seq_len(nrow(comparison)), function(index) {
  row <- comparison[index, ]
  paste0(
    "<tr><td>", html_escape(row$model), "</td>",
    "<td>", sprintf("%.2f", row$AIC), "</td>",
    "<td>", sprintf("%.2f", row$BIC), "</td>",
    "<td>", sprintf("%+.2f", row$delta_AIC_from_zisichel_g), "</td>",
    "<td>", sprintf("%+.2f", row$delta_BIC_from_zisichel_g), "</td>",
    "<td>", sprintf("%.4f", row$pit_sd), "</td>",
    "<td>", sprintf("%.4f", row$pit_cvm), "</td>",
    "<td>", sprintf("%.4f", row$worm_rms), "</td></tr>"
  )
}, character(1L)), collapse = "\n")

plot_sections <- paste(vapply(names(objects), function(label) {
  name <- slug(label)
  paste0(
    "<h3>", html_escape(label), "</h3><div class='plots'>",
    "<figure><img src='figures/", name, "-pit.png'><figcaption>PIT by visit</figcaption></figure>",
    "<figure><img src='figures/", name, "-worm.png'><figcaption>Worm plot by visit</figcaption></figure>",
    "<figure><img src='figures/", name, "-rootogram.png'><figcaption>Rootogram by visit</figcaption></figure>",
    "</div>"
  )
}, character(1L)), collapse = "\n")

html <- paste0(
  "<!doctype html><html><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>ACTG final margin sensitivity screen</title><style>",
  "body{font-family:Arial,sans-serif;max-width:1500px;margin:2rem auto;padding:0 1rem;color:#20262e;line-height:1.45}",
  "table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccd2d9;padding:.5rem;text-align:left}th{background:#eef2f6}",
  ".note{background:#fff7df;border-left:4px solid #d99b00;padding:1rem;margin:1rem 0}.good{background:#edf8f1;border-left:4px solid #2d8a50;padding:1rem;margin:1rem 0}",
  ".plots{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:1rem}figure{margin:0}img{width:100%;border:1px solid #d7dce2}figcaption{text-align:center;color:#53606d}",
  "@media(max-width:950px){.plots{grid-template-columns:1fr}}</style></head><body>",
  "<h1>ACTG final margin sensitivity screen</h1>",
  "<p>Full ACTG data; selected covariate structure; Gumbel dependence; 20-seed randomized-PIT summaries.</p>",
  "<div class='good'><strong>Recommendation:</strong> retain ZISICHEL--Gumbel. ZIBNB improves AIC/BIC modestly but has substantially poorer PIT, worm, and rootogram calibration. ZIPIG is inferior by information criteria and calibration.</div>",
  "<div class='note'><strong>Student-t sensitivity:</strong> ZIPIG--Student-t terminated its R worker. ZIBNB--Student-t remained active beyond the two-hour cap without returning a checkpoint and was stopped. These are not valid fitted alternatives.</div>",
  "<table><thead><tr><th>Model</th><th>AIC</th><th>BIC</th><th>&Delta;AIC vs ZISICHEL</th><th>&Delta;BIC</th><th>Visit-1 PIT SD</th><th>Visit-1 CvM</th><th>Visit-1 worm RMS</th></tr></thead><tbody>",
  rows, "</tbody></table><p><small>Uniform PIT target SD = 0.2887; lower CvM and worm RMS are better.</small></p>",
  plot_sections,
  "<h2>Files</h2><ul><li><a href='completed-family-copula-results.csv'>Model comparison CSV</a></li>",
  "<li><a href='completed-diagnostics-by-time.csv'>All visit-specific diagnostics</a></li>",
  "<li><a href='../run-final-family-copula-screen.R'>Checkpointed fitting script</a></li></ul>",
  "</body></html>"
)
writeLines(html, file.path(results_root, "index.html"), useBytes = TRUE)
