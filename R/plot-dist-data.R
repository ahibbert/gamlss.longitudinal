#' @keywords internal
#' @noRd
.plot_dist_fit_var <- function(fit, data_names, role, explicit = NULL) {
  if (!is.null(explicit)) {
    return(explicit)
  }
  candidates <- character(0)
  if (!is.null(fit)) {
    stored <- fit[[paste0(role, "_var")]]
    candidates <- c(candidates, stored)
    if (!is.null(fit$var_map)) {
      candidates <- c(candidates, names(fit$var_map)[fit$var_map == role])
    }
  }
  candidates <- c(candidates, role)
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit <- candidates[candidates %in% data_names]
  if (length(hit) > 0L) {
    return(hit[1L])
  }
  NULL
}

#' @keywords internal
#' @noRd
.plot_dist_time_values <- function(time) {
  if (is.factor(time)) {
    lev <- levels(time)
    lev[lev %in% as.character(unique(time))]
  } else {
    values <- unique(time)
    if (is.numeric(values) || is.integer(values)) sort(values) else sort(as.character(values))
  }
}

#' @keywords internal
#' @noRd
.plot_dist_normalise_data <- function(dataset, fit, subject_var, time_var, response_var) {
  if (is.null(dataset)) {
    if (!is.null(fit) && !is.null(fit$dataset)) {
      dataset <- fit$dataset
    } else {
      stop("Supply 'dataset', or supply a fitted 'fit' with stored data.", call. = FALSE)
    }
  }
  dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)

  if (!is.null(fit)) {
    subject_var <- .plot_dist_fit_var(fit, names(dataset), "subject", subject_var)
    time_var <- .plot_dist_fit_var(fit, names(dataset), "time", time_var)
    response_var <- .plot_dist_fit_var(fit, names(dataset), "response", response_var)
  }

  missing_inputs <- c(
    subject_var = is.null(subject_var),
    time_var = is.null(time_var),
    response_var = is.null(response_var)
  )
  if (any(missing_inputs)) {
    stop(
      "Raw-data plotting requires 'subject_var', 'time_var', and 'response_var'.",
      call. = FALSE
    )
  }

  missing_cols <- setdiff(c(subject_var, time_var, response_var), names(dataset))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in 'dataset': ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  data.frame(
    subject = dataset[[subject_var]],
    time = dataset[[time_var]],
    response = dataset[[response_var]],
    stringsAsFactors = FALSE
  )
}

#' Build per-time response and pseudo-observation slices for distribution plots
#'
#' @noRd
.plot_dist_margin_slices <- function(plot_data, time_values) {
  margin_data <- list()
  margin_pseudo <- list()

  for (i in seq_along(time_values)) {
    margin_data[[i]] <- plot_data[
      as.character(plot_data$time) == as.character(time_values[i]),
      c("subject", "response")
    ]

    r <- rank(margin_data[[i]]$response, ties.method = "average", na.last = "keep")
    n_obs <- sum(!is.na(margin_data[[i]]$response))
    u <- r / (n_obs + 1)
    margin_pseudo[[i]] <- data.frame(subject = margin_data[[i]]$subject, u = u)
  }

  list(margin_data = margin_data, margin_pseudo = margin_pseudo)
}

#' Build the off-diagonal panel data used by plot_dist()
#'
#' @noRd
.plot_dist_offdiag_data <- function(
    i,
    j,
    margin_data,
    margin_pseudo,
    offdiag_scale,
    transform) {
  if (offdiag_scale == "pseudo") {
    input_data <- merge(
      margin_pseudo[[i]],
      margin_pseudo[[j]],
      by = "subject",
      suffixes = c(".i", ".j"),
      all = FALSE
    )
    input_data <- input_data[complete.cases(input_data$u.i, input_data$u.j), c("u.i", "u.j")]
    names(input_data) <- c("X1", "X2")
    x_lab <- latex2exp::TeX(paste("$U_", i, "$"))
    y_lab <- latex2exp::TeX(paste("$U_", j, "$"))
    if (identical(transform, "normal")) {
      input_data$X1 <- stats::qnorm(.copula_v2_clamp01(input_data$X1))
      input_data$X2 <- stats::qnorm(.copula_v2_clamp01(input_data$X2))
      x_lab <- latex2exp::TeX(paste0("$\\Phi^{-1}(U_{", i, "})$"))
      y_lab <- latex2exp::TeX(paste0("$\\Phi^{-1}(U_{", j, "})$"))
    }
  } else {
    input_data <- merge(
      margin_data[[i]],
      margin_data[[j]],
      by = "subject",
      suffixes = c(".i", ".j"),
      all = FALSE
    )
    input_data <- input_data[complete.cases(input_data$response.i, input_data$response.j), c("response.i", "response.j")]
    names(input_data) <- c("X1", "X2")
    x_lab <- latex2exp::TeX(paste("$Y_", i, "$"))
    y_lab <- latex2exp::TeX(paste("$Y_", j, "$"))
  }

  list(input_data = input_data, x_lab = x_lab, y_lab = y_lab)
}
