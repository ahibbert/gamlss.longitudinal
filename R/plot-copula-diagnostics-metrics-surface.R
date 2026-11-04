.copula_v2_surface_metrics <- function(emp_density, fit_density, overlap_probs = c(0.7, 0.85, 0.95)) {
  emp <- as.numeric(emp_density)

  fit <- as.numeric(fit_density)

  ok <- is.finite(emp) & is.finite(fit)

  emp <- emp[ok]

  fit <- fit[ok]

  if (length(emp) == 0) {
    return(list(summary = data.frame(), overlap = data.frame()))
  }

  # Scale both surfaces to unit mass before computing distance metrics.

  emp <- pmax(emp, 0)

  fit <- pmax(fit, 0)

  emp <- emp / max(sum(emp), .Machine$double.eps)

  fit <- fit / max(sum(fit), .Machine$double.eps)

  summary_df <- data.frame(
    rmse = sqrt(mean((fit - emp)^2)),
    mae = mean(abs(fit - emp)),
    surface_cor = suppressWarnings(stats::cor(emp, fit, use = "complete.obs")),
    stringsAsFactors = FALSE
  )

  overlap_df <- do.call(rbind, lapply(overlap_probs, function(p) {
    thr_emp <- stats::quantile(emp, probs = p, na.rm = TRUE, type = 7)

    thr_fit <- stats::quantile(fit, probs = p, na.rm = TRUE, type = 7)

    mask_emp <- emp >= thr_emp

    mask_fit <- fit >= thr_fit

    union_n <- sum(mask_emp | mask_fit)

    iou <- if (union_n == 0) NA_real_ else sum(mask_emp & mask_fit) / union_n

    data.frame(level_prob = p, contour_iou = iou, stringsAsFactors = FALSE)
  }))

  list(summary = summary_df, overlap = overlap_df)
}
