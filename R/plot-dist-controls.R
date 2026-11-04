#' Normalize plot_dist() argument controls
#'
#' @noRd
.plot_dist_resolve_controls <- function(fit,
                                        offdiag_scale,
                                        transform,
                                        overlay,
                                        copula_dist,
                                        margin_dist) {
  if (!is.null(fit) && !inherits(fit, "gamlss.longitudinal")) {
    stop("'fit' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  offdiag_scale <- match.arg(offdiag_scale, c("pseudo", "response"))
  transform <- match.arg(transform, c("normal", "uniform"))
  overlay_choices <- c("none", "margin", "copula", "model")
  overlay <- if (is.null(overlay)) {
    if (!is.null(fit)) "model" else "none"
  } else {
    match.arg(overlay, overlay_choices)
  }
  if (overlay == "model" && is.null(fit)) {
    stop("'fit' is required when overlay = 'model'.", call. = FALSE)
  }
  if (overlay == "copula" && is.null(copula_dist)) {
    stop("'copula_dist' is required when overlay = 'copula'.", call. = FALSE)
  }
  if (overlay %in% c("copula", "model") && offdiag_scale != "pseudo") {
    warning("Copula overlays are only drawn when offdiag_scale = 'pseudo'.", call. = FALSE)
  }
  if (overlay == "margin" && is.null(margin_dist) && is.null(fit)) {
    stop("'margin_dist' is required when overlay = 'margin'.", call. = FALSE)
  }

  list(
    offdiag_scale = offdiag_scale,
    transform = transform,
    overlay = overlay
  )
}

#' Resolve optional fitted-model and selected-copula plot overlays
#'
#' @noRd
.plot_dist_overlay_state <- function(fit, dataset, num_margins, copula_dist) {
  fit_diag_data <- NULL
  fit_pair_data <- NULL
  fit_copula_spec <- NULL
  if (!is.null(fit)) {
    fit_overlay_data <- .plot_dist_fit_overlay_data(fit, dataset, num_margins)
    fit_diag_data <- fit_overlay_data$diag_data
    fit_pair_data <- fit_overlay_data$pair_data
    fit_copula_spec <- fit_overlay_data$copula_spec
  }
  copula_spec <- if (!is.null(copula_dist)) .plot_copula_selection_spec(copula_dist) else NULL

  list(
    fit_diag_data = fit_diag_data,
    fit_pair_data = fit_pair_data,
    fit_copula_spec = fit_copula_spec,
    copula_spec = copula_spec
  )
}
