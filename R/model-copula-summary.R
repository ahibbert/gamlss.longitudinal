#' Summarise fitted copula parameters by time
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param lags Integer lags to assess, measured in ordered time steps.
#' @param stat Character summary statistic for fitted values, one of "mean" or "median".
#'
#' @return A data frame with fitted theta and tau summaries by time.
#' @export
copula_time_summary <- function(object, lags = 1, stat = c("mean", "median")) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.")
  }

  stat <- match.arg(stat)
  copula_info <- get_copula_dist(object$copula_dist)
  has_zeta <- "zeta" %in% copula_info$parameters

  fit_data <- .copula_v2_fit_data(object)
  pair_data <- .copula_v2_pair_data(fit_data, lags = lags)

  agg_fun <- .gl_copula_summary_agg_fun(stat)
  time_summary <- .gl_copula_time_table(fit_data, agg_fun, has_zeta)
  pair_summary <- .gl_copula_pair_table(pair_data, agg_fun, has_zeta)
  tidy_data <- .gl_copula_drop_zeta_if_absent(fit_data, pair_data, has_zeta)
  fit_data <- tidy_data$fit_data
  pair_data <- tidy_data$pair_data

  out <- list(
    time_summary = time_summary,
    pair_summary = pair_summary,
    fit_data = fit_data,
    pair_data = pair_data
  )
  class(out) <- "copula_time_summary"
  out
}
