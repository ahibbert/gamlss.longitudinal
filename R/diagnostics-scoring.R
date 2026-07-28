.gl_crps_sample <- function(y, draws) {
  mean(abs(draws - y)) - 0.5 * mean(abs(outer(draws, draws, "-")))
}

#' @export

proscore.gamlss.longitudinal <- function(object, type = c("logs", "crps", "mae", "mse", "dss"), crps_grid = 25, ...) {
  type <- match.arg(type, several.ok = TRUE)

  diag_data <- .gl_diag_data(object)

  y <- diag_data$response

  params <- diag_data$params

  mu_hat <- diag_data$mu_hat

  sigma_hat <- diag_data$sigma_hat

  out <- setNames(numeric(length(type)), type)

  density_hat <- .gl_call_family_fun("d", diag_data$family, y, params)

  if ("logs" %in% type) {
    out["logs"] <- mean(-log(pmax(density_hat, .Machine$double.eps)), na.rm = TRUE)
  }

  if ("mae" %in% type) {
    out["mae"] <- mean(abs(y - mu_hat), na.rm = TRUE)
  }

  if ("mse" %in% type) {
    out["mse"] <- mean((y - mu_hat)^2, na.rm = TRUE)
  }

  if ("dss" %in% type) {
    out["dss"] <- mean(log(sigma_hat^2) + ((y - mu_hat)^2 / sigma_hat^2), na.rm = TRUE)
  }

  if ("crps" %in% type) {
    p_grid <- seq_len(crps_grid) / (crps_grid + 1)

    sample_mat <- vapply(p_grid, function(prob) {
      .gl_call_family_fun("q", diag_data$family, prob, params)
    }, numeric(length(y)))

    if (is.null(dim(sample_mat))) {
      sample_mat <- matrix(sample_mat, ncol = 1)
    }

    out["crps"] <- mean(vapply(seq_len(nrow(sample_mat)), function(i) {
      .gl_crps_sample(y[i], sample_mat[i, ])
    }, numeric(1)), na.rm = TRUE)
  }

  out
}

#' @export

procast.gamlss.longitudinal <- function(object, type = c("quantile", "cdf", "density"), at = c(0.025, 0.5, 0.975), newdata = NULL, ...) {
  type <- match.arg(type)

  require_response <- type %in% c("cdf", "density")

  diag_data <- .gl_fitted_distribution(object, newdata = newdata, require_response = require_response)

  y <- diag_data$response

  params <- diag_data$params

  if (type == "quantile") {
    quantile_df <- data.frame(response = y)

    for (prob in at) {
      quantile_df[[paste0("q", gsub("^0\\.", "", format(prob, trim = TRUE)))]] <- .gl_call_family_fun("q", diag_data$family, prob, params)
    }

    return(quantile_df)
  }

  if (type == "cdf") {
    return(data.frame(response = y, cdf = .gl_call_family_fun("p", diag_data$family, y, params)))
  }

  data.frame(response = y, density = .gl_call_family_fun("d", diag_data$family, y, params))
}
