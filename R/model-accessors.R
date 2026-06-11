#' Log-likelihood for a fitted longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return Named log-likelihood components.
#' @export
logLik.gamlss.longitudinal=function(object, ...) {
  return(object$calc_lik_out$log_lik)
}

#' Coefficients for a fitted longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return Named coefficient vector.
#' @export
coef.gamlss.longitudinal=function(object, ...) {
  return(object$par)
}

#' Access components of a fitted longitudinal GAMLSS-copula model
#'
#' These S3 methods expose the standard regression components expected by
#' model-auditing workflows: formulas, terms, observation counts, model frames,
#' fitted values, and residuals. They return the expanded internal data by
#' default because `gamlss_longitudinal()` represents structurally missing
#' subject-time combinations explicitly.
#'
#' @param x,object,formula A fitted `gamlss.longitudinal` object.
#' @param parameter Distributional parameter to extract. Defaults to `"mu"`.
#' @param internal Logical; return internally translated formulas when `TRUE`.
#' @param type Observation-count, model-frame, or residual type.
#' @param finite Logical; restrict fitted values or residuals to finite observed
#'   responses and finite fitted parameters.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return A formula, terms object, integer count, data frame, or numeric vector.
#' @name gamlss_longitudinal_accessors
NULL

#' @rdname gamlss_longitudinal_accessors
#' @export
formula.gamlss.longitudinal <- function(x, parameter = c("mu", "sigma", "nu", "tau", "theta", "zeta"), internal = FALSE, ...) {
  parameter <- match.arg(parameter)
  formulas <- if (isTRUE(internal)) x$formulas_int else x$formulas
  fml <- formulas[[parameter]]
  if (inherits(fml, "formula")) {
    return(fml)
  }
  stats::as.formula(fml)
}

#' @rdname gamlss_longitudinal_accessors
#' @export
terms.gamlss.longitudinal <- function(x, parameter = c("mu", "sigma", "nu", "tau", "theta", "zeta"), internal = FALSE, ...) {
  stats::terms(formula.gamlss.longitudinal(x, parameter = parameter, internal = internal, ...))
}

#' @rdname gamlss_longitudinal_accessors
#' @export
nobs.gamlss.longitudinal <- function(object, type = c("observed", "expanded", "submitted"), ...) {
  type <- match.arg(type)
  if (identical(type, "submitted") && !is.null(object$dataset_original)) {
    return(nrow(object$dataset_original))
  }
  if (identical(type, "expanded") && !is.null(object$dataset)) {
    return(nrow(object$dataset))
  }
  sum(is.finite(object$response))
}

#' @rdname gamlss_longitudinal_accessors
#' @export
model.frame.gamlss.longitudinal <- function(formula, type = c("expanded", "observed", "submitted"), ...) {
  object <- formula
  type <- match.arg(type)
  out <- if (identical(type, "submitted") && !is.null(object$dataset_original)) {
    as.data.frame(object$dataset_original, stringsAsFactors = FALSE)
  } else {
    as.data.frame(object$dataset, stringsAsFactors = FALSE)
  }
  if (identical(type, "observed") && "response" %in% names(out)) {
    out <- out[is.finite(out$response), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

.gl_fitted_parameter_values <- function(object, parameter = "mu") {
  copula_link <- get_copula_dist(object$copula_dist)$copula_link
  eta_out <- calc_eta(
    par_cov = object$par,
    mm = object$model_matrix,
    margin_dist = object$margin_dist,
    copula_link = copula_link,
    par_s = object$par_s
  )
  params <- eta_out$eta_inv[names(object$margin_dist$parameters)]
  if (!parameter %in% names(params)) {
    stop("Parameter '", parameter, "' is not available in the fitted margin.", call. = FALSE)
  }
  as.numeric(params[[parameter]])
}

#' @rdname gamlss_longitudinal_accessors
#' @export
fitted.gamlss.longitudinal <- function(object, parameter = "mu", finite = FALSE, ...) {
  fit <- .gl_fitted_parameter_values(object, parameter = parameter)
  if (isTRUE(finite)) {
    keep <- is.finite(object$response) & is.finite(fit)
    fit <- fit[keep]
  }
  fit
}

#' @rdname gamlss_longitudinal_accessors
#' @export
residuals.gamlss.longitudinal <- function(object, type = c("response", "pearson", "quantile"), finite = TRUE, ...) {
  type <- match.arg(type)
  mu <- .gl_fitted_parameter_values(object, parameter = "mu")
  y <- as.numeric(object$response)
  keep <- is.finite(y) & is.finite(mu)

  if (identical(type, "response")) {
    out <- y - mu
  } else if (identical(type, "pearson")) {
    sigma <- if ("sigma" %in% names(object$margin_dist$parameters)) {
      .gl_fitted_parameter_values(object, parameter = "sigma")
    } else {
      rep(stats::sd(y, na.rm = TRUE), length(y))
    }
    sigma <- pmax(as.numeric(sigma), .Machine$double.eps)
    out <- (y - mu) / sigma
    keep <- keep & is.finite(sigma)
  } else {
    copula_link <- get_copula_dist(object$copula_dist)$copula_link
    eta_out <- calc_eta(
      par_cov = object$par,
      mm = object$model_matrix,
      margin_dist = object$margin_dist,
      copula_link = copula_link,
      par_s = object$par_s
    )
    params <- eta_out$eta_inv[names(object$margin_dist$parameters)]
    pit <- .gl_call_family_fun("p", object$margin_dist$family[1], y, params)
    pit <- pmin(pmax(pit, .Machine$double.eps), 1 - .Machine$double.eps)
    out <- stats::qnorm(pit)
    keep <- keep & is.finite(out)
  }

  if (isTRUE(finite)) {
    out <- out[keep]
  }
  out
}

#' Summarise a fitted longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param include_vcov Logical; include coefficient uncertainty summaries when
#'   possible.
#' @param numderiv Logical; use the numerical Hessian path when computing
#'   variance-covariance output.
#' @param ci_level Confidence level for coefficient intervals.
#' @param ... Additional arguments passed to `vcov.gamlss.longitudinal()`.
#'
#' @return A `summary.gamlss.longitudinal` object with model metadata,
#'   coefficient summaries, smooth-term summaries, and optional vcov output.
#' @export

summary.gamlss.longitudinal = function(
  object,
  include_vcov = TRUE,
  numderiv = FALSE,
  ci_level = 0.95,
  ...
) {
  if(!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  n_obs = length(object$response)
  n_subjects = length(unique(object$response_subject))
  n_timepoints = length(unique(object$response_margin))

  n_fixed = length(object$par)
  n_smooth_terms = 0
  if(!is.null(object$par_s)) {
    n_smooth_terms = sum(vapply(object$par_s, length, integer(1)))
  }

  edf_smooth = NA_real_
  if(!is.null(object$df_s) && length(object$df_s) > 0) {
    df_vals = suppressWarnings(as.numeric(unlist(object$df_s, use.names = FALSE)))
    df_vals = df_vals[is.finite(df_vals)]
    if(length(df_vals) > 0) {
      edf_smooth = sum(df_vals)
    }
  }

  loglik_vec = c(marginal = NA_real_, copula = NA_real_, joint = NA_real_)
  if(!is.null(object$calc_lik_out_end) && !is.null(object$calc_lik_out_end$log_lik)) {
    ll_in = object$calc_lik_out_end$log_lik
    for(nm in c("marginal", "copula", "joint")) {
      if(nm %in% names(ll_in)) {
        loglik_vec[nm] = as.numeric(ll_in[[nm]])
      }
    }
  }
  loglik_joint = as.numeric(loglik_vec["joint"])

  p_cop = object$par[grepl("theta", names(object$par)) | grepl("zeta", names(object$par))]
  p_mar = object$par[!(grepl("theta", names(object$par)) | grepl("zeta", names(object$par)))]

  df_s_total = 0
  df_s_cop_total = 0
  df_s_margin_total = 0
  if(!is.null(object$df_s) && length(object$df_s) > 0) {
    for(par_name in names(object$df_s)) {
      df_val = suppressWarnings(sum(as.numeric(unlist(object$df_s[[par_name]])), na.rm = TRUE))
      if(!is.finite(df_val)) df_val = 0

      if(par_name %in% c("theta", "zeta")) {
        df_s_cop_total = df_s_cop_total + df_val
      } else {
        df_s_margin_total = df_s_margin_total + df_val
      }
      df_s_total = df_s_total + df_val
    }
  }

  edf_vec = c(
    marginal = length(p_mar) + df_s_margin_total,
    copula = length(p_cop) + df_s_cop_total,
    joint = length(object$par) + df_s_total
  )

  aic_vec = -2 * loglik_vec + 2 * edf_vec
  bic_vec = -2 * loglik_vec + log(max(1, n_obs)) * edf_vec
  model_selection = rbind(LogLik = loglik_vec, AIC = aic_vec, BIC = bic_vec, EDF = edf_vec)

  coef_tbl = data.frame(
    term = names(object$par),
    estimate = as.numeric(object$par),
    std_error = NA_real_,
    p_value = NA_real_,
    signif = NA_character_,
    stringsAsFactors = FALSE
  )

  coef_tbl$.original_order = seq_len(nrow(coef_tbl))
  coef_tbl$parameter = sub("\\..*$", "", coef_tbl$term)

  param_order = c("mu", "sigma", "nu", "tau", "theta", "zeta")
  coef_tbl$.param_rank = match(coef_tbl$parameter, param_order)
  coef_tbl$.param_rank[is.na(coef_tbl$.param_rank)] = length(param_order) + 1L

  vcov_out = NULL
  if(isTRUE(include_vcov)) {
    vcov_out = .resolve_vcov(
      object = object,
      numderiv = numderiv,
      extra_args = list(...)
    )

    if(!is.null(vcov_out$vcov) && !is.null(vcov_out$vcov$overall)) {
      V = vcov_out$vcov$overall
      se = NULL
      if(!is.null(vcov_out$se) && !is.null(vcov_out$se$overall)) {
        se = as.numeric(vcov_out$se$overall)
        se_names = names(vcov_out$se$overall)
      } else {
        se = sqrt(pmax(0, diag(V)))
        se_names = names(diag(V))
      }

      if(is.null(se_names) && !is.null(rownames(V)) && length(rownames(V)) == length(se)) {
        se_names = rownames(V)
      }

      if(!is.null(se_names)) {
        names(se) = se_names
      }

      if(!is.null(names(se))) {
        idx = match(coef_tbl$term, names(se))
        coef_tbl$std_error = se[idx]
      } else if(length(se) == nrow(coef_tbl)) {
        coef_tbl$std_error = se
      }

      z_abs = abs(coef_tbl$estimate / coef_tbl$std_error)
      coef_tbl$p_value = 2 * stats::pnorm(z_abs, lower.tail = FALSE)
      coef_tbl$signif = ifelse(
        is.na(coef_tbl$p_value),
        NA_character_,
        ifelse(coef_tbl$p_value < 0.001, "***",
               ifelse(coef_tbl$p_value < 0.01, "**",
                      ifelse(coef_tbl$p_value < 0.05, "*",
                             ifelse(coef_tbl$p_value < 0.1, ".", " "))))
      )
    }
  }

  coef_tbl = coef_tbl[order(coef_tbl$.param_rank, coef_tbl$.original_order), , drop = FALSE]
  rownames(coef_tbl) = NULL

  out = list(
    model = list(
      margin_dist = if(!is.null(object$margin_dist$family[1])) as.character(object$margin_dist$family[1]) else NA_character_,
      copula_dist = object$copula_dist,
      n_obs = n_obs,
      n_subjects = n_subjects,
      n_timepoints = n_timepoints,
      n_fixed = n_fixed,
      n_smooth_terms = n_smooth_terms,
      edf_smooth = edf_smooth
    ),
    fit = list(
      logLik = loglik_joint,
      AIC = as.numeric(aic_vec["joint"]),
      BIC = as.numeric(bic_vec["joint"]),
      ci_level = ci_level,
      vcov_included = isTRUE(include_vcov),
      vcov_numderiv = isTRUE(numderiv),
      vcov_method = vcov_out$method %||% object$vcov_meta$method_used %||% object$vcov_meta$method %||% NA_character_,
      vcov_method_requested = vcov_out$method_requested %||% object$vcov_meta$method %||% NA_character_,
      hessian_diagnostics = vcov_out$hessian_diagnostics %||% NULL,
      model_selection = model_selection
    ),
    smooth_terms = {
      st = list()
      if(!is.null(object$par_s) && length(object$par_s) > 0) {
        for(par_name in names(object$par_s)) {
          if(length(object$par_s[[par_name]]) == 0) next
          for(s_name in names(object$par_s[[par_name]])) {
            smooth_edf = NA_real_
            if(!is.null(object$df_s) && !is.null(object$df_s[[par_name]]) && s_name %in% names(object$df_s[[par_name]])) {
              smooth_edf = suppressWarnings(as.numeric(object$df_s[[par_name]][[s_name]])[1])
            }
            st[[length(st) + 1]] = data.frame(
              parameter = par_name,
              smooth_term = s_name,
              edf = smooth_edf,
              stringsAsFactors = FALSE
            )
          }
        }
      }
      if(length(st) == 0) {
        data.frame(parameter = character(0), smooth_term = character(0), edf = numeric(0), stringsAsFactors = FALSE)
      } else {
        do.call(rbind, st)
      }
    },
    coefficients = within(coef_tbl, {
      .original_order = NULL
      .param_rank = NULL
    }),
    vcov = vcov_out
  )
  class(out) = "summary.gamlss.longitudinal"
  out
}

.copula_v2_gaussian_limit_warning <- function(object, threshold = 7.5) {
  if (is.null(object) || is.null(object$model$copula_dist)) {
    return(NULL)
  }

  copula_dist <- object$model$copula_dist
  if (!identical(copula_dist, "t")) {
    return(NULL)
  }

  coef_tbl <- object$coefficients
  if (is.null(coef_tbl) || nrow(coef_tbl) == 0 || !("parameter" %in% names(coef_tbl))) {
    return(NULL)
  }

  zeta_rows <- coef_tbl[coef_tbl$parameter == "zeta" & is.finite(coef_tbl$estimate), , drop = FALSE]
  if (nrow(zeta_rows) == 0) {
    return(NULL)
  }

  zeta_link <- zeta_rows$estimate
  near_limit <- is.finite(zeta_link) & zeta_link >= threshold
  if (!any(near_limit)) {
    return(NULL)
  }

  zeta_nat <- exp(zeta_link[near_limit]) + 2
  zeta_label <- paste0(formatC(zeta_link[near_limit], format = "f", digits = 2), collapse = ", ")
  df_label <- paste0(formatC(zeta_nat, format = "f", digits = 1), collapse = ", ")

  paste0(
    "WARNING: t-copula zeta is near the Gaussian-limit regime (link-scale zeta = ",
    zeta_label,
    "; implied degrees of freedom ",
    df_label,
    "). The t-copula may be collapsing toward Gaussian dependence."
  )
}

#' @export
print.summary.gamlss.longitudinal = function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nGAMLSS Longitudinal Model Summary\n")
  cat("--------------------------------\n")
  cat("Margin distribution:", x$model$margin_dist, "\n")
  cat("Copula distribution:", x$model$copula_dist, "\n")

  copula_warning <- .copula_v2_gaussian_limit_warning(x)
  if (!is.null(copula_warning)) {
    cat(copula_warning, "\n")
  }

  cat("Observations:", x$model$n_obs,
      " | Subjects:", x$model$n_subjects,
      " | Time points:", x$model$n_timepoints, "\n")
  cat("Fixed coefficients:", x$model$n_fixed,
      " | Smooth terms:", x$model$n_smooth_terms,
      " | Smooth EDF:", format(round(x$model$edf_smooth, digits), nsmall = 2), "\n")
  if (isTRUE(x$fit$vcov_included)) {
    cat(
      "VCOV:", x$fit$vcov_method %||% "unknown",
      " | Requested:", x$fit$vcov_method_requested %||% "unknown",
      "\n"
    )
    hd <- x$fit$hessian_diagnostics
    if (!is.null(hd) && is.finite(hd$condition_number %||% NA_real_)) {
      cat("Hessian condition number:", formatC(hd$condition_number, digits = digits, format = "fg"), "\n")
    }
  }

  cat("\nFixed coefficients:\n")
  cat("--------------------\n")
  coef_tbl = x$coefficients
  p_value_raw = coef_tbl$p_value
  coef_tbl$estimate = round(coef_tbl$estimate, digits)
  coef_tbl$std_error = round(coef_tbl$std_error, digits)
  coef_tbl$p_value = round(coef_tbl$p_value, digits + 1)

  fmt_num = function(v, d) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d))
  fmt_p_value = function(v, v_raw, d) {
    ifelse(
      is.na(v_raw),
      "NA",
      ifelse(v_raw > 0 & v_raw < 10^(-d), paste0("<", formatC(10^(-d), format = "f", digits = d)), fmt_num(v, d))
    )
  }
  coef_disp = data.frame(
    term = as.character(coef_tbl$term),
    estimate = fmt_num(coef_tbl$estimate, digits),
    std_error = fmt_num(coef_tbl$std_error, digits),
    p_value = fmt_p_value(coef_tbl$p_value, p_value_raw, digits + 1),
    signif = ifelse(is.na(coef_tbl$signif), "", as.character(coef_tbl$signif)),
    parameter = as.character(coef_tbl$parameter),
    stringsAsFactors = FALSE
  )

  w_term = max(nchar("term"), nchar(coef_disp$term, type = "width"), na.rm = TRUE)
  w_est = max(nchar("estimate"), nchar(coef_disp$estimate, type = "width"), na.rm = TRUE)
  w_se = max(nchar("std_error"), nchar(coef_disp$std_error, type = "width"), na.rm = TRUE)
  w_p = max(nchar("p_value"), nchar(coef_disp$p_value, type = "width"), na.rm = TRUE)
  w_sig = max(nchar("signif"), nchar(coef_disp$signif, type = "width"), na.rm = TRUE)

  format_row = function(term, estimate, std_error, p_value, signif) {
    sprintf(
      "%-*s  %*s  %*s  %*s  %-*s",
      w_term, term,
      w_est, estimate,
      w_se, std_error,
      w_p, p_value,
      w_sig, signif
    )
  }

  print_coef_block = function(block, prefix = "    ") {
    hdr = format_row("term", "estimate", "std_error", "p_value", "signif")
    cat(prefix, hdr, "\n", sep = "")
    for(ii in seq_len(nrow(block))) {
      row_txt = format_row(
        block$term[ii],
        block$estimate[ii],
        block$std_error[ii],
        block$p_value[ii],
        block$signif[ii]
      )
      cat(prefix, row_txt, "\n", sep = "")
    }
  }

  param_order = c("mu", "sigma", "nu", "tau", "theta", "zeta")
  params_present = unique(coef_tbl$parameter)
  params_print = c(param_order[param_order %in% params_present], setdiff(params_present, param_order))

  for(k in seq_along(params_print)) {
    p = params_print[k]
    block = coef_disp[coef_disp$parameter == p, c("term", "estimate", "std_error", "p_value", "signif"), drop = FALSE]
    cat(sprintf("  [%s]\n", p))
    print_coef_block(block, prefix = "    ")
    if(k < length(params_print)) cat("  --------------------\n")
  }

  cat("  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")


  cat("\nSmooth terms:\n")
  cat("--------------------\n")
  if(!is.null(x$smooth_terms) && nrow(x$smooth_terms) > 0) {
    smooth_disp = x$smooth_terms
    if("edf" %in% names(smooth_disp)) {
      smooth_disp$edf = round(smooth_disp$edf, digits)
    }
    print(smooth_disp, row.names = FALSE)
    cat("Use plot(object) to visualize smooth and fixed terms with confidence bands.\n")
  } else {
    cat("None\n")
  }

  cat("\nModel Selection Criteria:\n")
  cat("--------------------\n")
  if(!is.null(x$fit$model_selection)) {
    print(round(x$fit$model_selection, digits))
  } else {
    fit_tbl = data.frame(
      metric = c("logLik", "AIC", "BIC"),
      value = c(x$fit$logLik, x$fit$AIC, x$fit$BIC),
      stringsAsFactors = FALSE
    )
    print(fit_tbl, row.names = FALSE, digits = digits)
  }
  cat("--------------------------------\n")

  invisible(x)
}

#' Plot all smooth terms with confidence bands
#'
#' This utility plots every smooth term in a fitted `gamlss.longitudinal` object
#' and computes pointwise confidence bands using the smooth coefficient
#' covariance matrices returned by `vcov.gamlss.longitudinal()`.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param vcov_obj Optional output from `vcov(object, ...)`. If `NULL`, this is
#' computed internally with the analytical vcov path.
#' @param data Optional data frame containing original covariates used for the
#' x-axis variable of each smooth.
#' @param ci_level Confidence level for pointwise intervals.
#' @param ncol Number of plot columns (defaults to 2 or fewer if needed).
#' @param ci_col Color for confidence bands.
#' @param fit_col Color for fitted smooth line.
#' @param ci_lty Line type for confidence bands.
#' @param fit_lwd Line width for fitted smooth.
#' @param sort_x Logical; sort points by x before plotting lines.
#' @param even_grid Logical; if TRUE, plot smooths on an evenly spaced x-grid
#' built over observed x-range.
#' @param grid_n Number of grid points when `even_grid = TRUE`.
#' @param fallback_to_index Logical; if x variable cannot be inferred, plot
#' against row index.
#' @param setup_mfrow Logical; if TRUE (default), configure par(mfrow) inside
#' this function. Set FALSE when caller configures layout.
#' @param show_legend Logical; if TRUE, draw a small legend in each panel.
#'
#' @return Invisibly returns a nested list with x, fitted values, standard
#' errors, and confidence limits for each smooth term.
#' @export
