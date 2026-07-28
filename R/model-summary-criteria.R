#' Build summary log-likelihood, information criteria, and EDF tables
#'
#' @noRd
.gl_summary_fit_criteria <- function(object, n_obs) {
  loglik_vec <- c(marginal = NA_real_, copula = NA_real_, joint = NA_real_)
  if (!is.null(object$calc_lik_out_end) && !is.null(object$calc_lik_out_end$log_lik)) {
    ll_in <- object$calc_lik_out_end$log_lik
    for (nm in c("marginal", "copula", "joint")) {
      if (nm %in% names(ll_in)) {
        loglik_vec[nm] <- as.numeric(ll_in[[nm]])
      }
    }
  }
  loglik_joint <- as.numeric(loglik_vec["joint"])

  p_cop <- object$par[grepl("theta", names(object$par)) | grepl("zeta", names(object$par))]
  p_mar <- object$par[!(grepl("theta", names(object$par)) | grepl("zeta", names(object$par)))]

  df_s_total <- 0
  df_s_cop_total <- 0
  df_s_margin_total <- 0
  if (!is.null(object$df_s) && length(object$df_s) > 0) {
    for (par_name in names(object$df_s)) {
      df_val <- suppressWarnings(sum(as.numeric(unlist(object$df_s[[par_name]])), na.rm = TRUE))
      if (!is.finite(df_val)) df_val <- 0

      if (par_name %in% c("theta", "zeta")) {
        df_s_cop_total <- df_s_cop_total + df_val
      } else {
        df_s_margin_total <- df_s_margin_total + df_val
      }
      df_s_total <- df_s_total + df_val
    }
  }

  edf_vec <- c(
    marginal = length(p_mar) + df_s_margin_total,
    copula = length(p_cop) + df_s_cop_total,
    joint = length(object$par) + df_s_total
  )

  aic_vec <- -2 * loglik_vec + 2 * edf_vec
  bic_vec <- -2 * loglik_vec + log(max(1, n_obs)) * edf_vec
  model_selection <- rbind(LogLik = loglik_vec, AIC = aic_vec, BIC = bic_vec, EDF = edf_vec)

  list(
    loglik_vec = loglik_vec,
    loglik_joint = loglik_joint,
    edf_vec = edf_vec,
    aic_vec = aic_vec,
    bic_vec = bic_vec,
    model_selection = model_selection
  )
}
