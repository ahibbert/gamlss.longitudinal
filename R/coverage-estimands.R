.coverage_capture_conditions <- function(expr) {

  warnings <- character(0)

  value <- withCallingHandlers(

    tryCatch(expr, error = function(e) e),

    warning = function(w) {

      warnings <<- c(warnings, conditionMessage(w))

      invokeRestart("muffleWarning")

    }

  )

  list(value = value, warnings = warnings)

}


#' @keywords internal

#' @noRd

.coverage_taxonomy <- function(success, error, warnings, elapsed, max_elapsed_sec) {

  if (isTRUE(success)) return("ok")

  text <- paste(c(error, warnings), collapse = " ")

  if (grepl("max_elapsed_sec|elapsed", text, ignore.case = TRUE) || elapsed >= max_elapsed_sec) {

    return("convergence timeout")

  }

  if (grepl("starting|start", text, ignore.case = TRUE)) return("start failure")

  if (grepl("non-finite|NA/NaN|Inf|infinite", text, ignore.case = TRUE)) return("non-finite likelihood")

  if (grepl("boundary|outside|domain", text, ignore.case = TRUE)) return("boundary estimate")

  if (grepl("derivative|hessian|gradient", text, ignore.case = TRUE)) return("numerical derivative issue")

  "fit error"

}


#' @keywords internal

#' @noRd

.coverage_natural_estimates <- function(fit) {

  if (!inherits(fit, "gamlss.longitudinal")) return(stats::setNames(numeric(0), character(0)))

  margin_dist <- fit$margin_dist

  copula_link <- get_copula_dist(fit$copula_dist)$copula_link

  par <- fit$par

  out <- numeric(0)

  for (par_name in c(names(margin_dist$parameters), get_copula_dist(fit$copula_dist)$parameters)) {

    coef_name <- paste0(par_name, ".intercept")

    if (!coef_name %in% names(par)) next

    linkinv <- if (par_name %in% names(margin_dist$parameters)) {

      margin_dist[[paste0(par_name, ".linkinv")]]

    } else {

      copula_link[[paste0(par_name, ".linkinv")]]

    }

    if (is.null(linkinv)) next

    out[par_name] <- as.numeric(linkinv(unname(par[coef_name]))[1])

  }

  out

}


#' @keywords internal

#' @noRd

.coverage_truth_summary <- function(dat, copula) {

  truth_cols <- grep("^true_", names(dat), value = TRUE)

  out <- vapply(truth_cols, function(nm) mean(dat[[nm]], na.rm = TRUE), numeric(1))

  names(out) <- sub("^true_", "", names(out))

  if ("theta" %in% names(out)) {

    par2 <- if ("true_zeta" %in% names(dat)) dat$true_zeta else 0

    valid <- is.finite(dat$true_theta) & is.finite(par2)

    out["copula_tau"] <- if (any(valid)) {

      mean(.copula_par_to_tau(copula, dat$true_theta[valid], par2[valid]), na.rm = TRUE)

    } else {

      NA_real_

    }

  }

  out

}


#' @keywords internal

#' @noRd

.coverage_truth_adjacent_start <- function(dat, family, copula, design) {

  if (!identical(design, "intercept")) {

    return(NA)

  }


  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  copula_link <- get_copula_dist(copula)$copula_link

  truth <- .coverage_truth_summary(dat, copula)

  out <- numeric(0)


  for (par_name in names(margin_dist$parameters)) {

    if (!par_name %in% names(truth)) next

    linkfun <- margin_dist[[paste0(par_name, ".linkfun")]]

    if (is.null(linkfun) || !is.finite(truth[[par_name]])) next

    out[paste0(par_name, ".intercept")] <- as.numeric(linkfun(truth[[par_name]]))[1]

  }


  if ("theta" %in% names(truth) && is.finite(truth[["theta"]])) {

    out["theta.intercept"] <- as.numeric(copula_link$theta.linkfun(truth[["theta"]]))[1]

  }

  if (identical(copula, "t") && "zeta" %in% names(truth) && is.finite(truth[["zeta"]])) {

    out["zeta.intercept"] <- as.numeric(copula_link$zeta.linkfun(truth[["zeta"]]))[1]

  }


  if (length(out) == 0L) NA else out

}


#' @keywords internal

#' @noRd

.coverage_true_eta_coefficients <- function(dat, family, copula, design) {

  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  rows <- list()


  add_margin_rows <- function(par_name) {

    truth_col <- paste0("true_", par_name)

    if (!truth_col %in% names(dat)) return(NULL)

    linkfun <- margin_dist[[paste0(par_name, ".linkfun")]]

    if (is.null(linkfun)) return(NULL)

    ok <- is.finite(dat[[truth_col]])

    if (!any(ok)) return(NULL)

    eta <- as.numeric(linkfun(dat[[truth_col]][ok]))

    if (!all(is.finite(eta))) return(NULL)

    varies_with_x <- identical(design, "covariate") ||

      (identical(design, "scale") && identical(par_name, "sigma"))

    if (isTRUE(varies_with_x) && "x" %in% names(dat)) {

      fit <- stats::lm(eta ~ dat$x[ok])

      cf <- stats::coef(fit)

      data.frame(

        parameter = par_name,

        term = c("intercept", "x"),

        true_eta = as.numeric(c(cf[[1]], cf[[2]])),

        stringsAsFactors = FALSE

      )

    } else {

      data.frame(

        parameter = par_name,

        term = "intercept",

        true_eta = mean(eta),

        stringsAsFactors = FALSE

      )

    }

  }


  for (par_name in names(margin_dist$parameters)) {

    rows[[length(rows) + 1L]] <- add_margin_rows(par_name)

  }


  if ("true_theta" %in% names(dat)) {

    theta_ok <- is.finite(dat$true_theta)

    if (any(theta_ok)) {

      theta_link <- get_copula_dist(copula)$copula_link$theta.linkfun

      theta_eta <- as.numeric(theta_link(dat$true_theta[theta_ok]))

      if (identical(design, "time_dependence") && "time" %in% names(dat)) {

        fit <- stats::lm(theta_eta ~ dat$time[theta_ok])

        cf <- stats::coef(fit)

        rows[[length(rows) + 1L]] <- data.frame(

          parameter = "theta",

          term = c("intercept", "time_covariate"),

          true_eta = as.numeric(c(cf[[1]], cf[[2]])),

          stringsAsFactors = FALSE

        )

      } else {

        rows[[length(rows) + 1L]] <- data.frame(

          parameter = "theta",

          term = "intercept",

          true_eta = mean(theta_eta),

          stringsAsFactors = FALSE

        )

      }

    }

  }


  if (identical(copula, "t") && "true_zeta" %in% names(dat)) {

    zeta_ok <- is.finite(dat$true_zeta)

    if (any(zeta_ok)) {

      zeta_link <- get_copula_dist(copula)$copula_link$zeta.linkfun

      rows[[length(rows) + 1L]] <- data.frame(

        parameter = "zeta",

        term = "intercept",

        true_eta = mean(as.numeric(zeta_link(dat$true_zeta[zeta_ok]))),

        stringsAsFactors = FALSE

      )

    }

  }


  rows <- rows[!vapply(rows, is.null, logical(1))]

  if (length(rows) == 0L) {

    return(data.frame(parameter = character(0), term = character(0), true_eta = numeric(0)))

  }

  do.call(rbind, rows)

}


#' @keywords internal

#' @noRd

.coverage_longitudinal_eta_estimates <- function(fit) {

  if (!inherits(fit, "gamlss.longitudinal")) {

    return(data.frame(parameter = character(0), term = character(0), estimate_eta = numeric(0)))

  }

  par <- fit$par

  pieces <- strsplit(names(par), ".", fixed = TRUE)

  data.frame(

    parameter = vapply(pieces, `[`, character(1), 1L),

    term = vapply(pieces, function(x) paste(x[-1L], collapse = "."), character(1)),

    estimate_eta = as.numeric(par),

    stringsAsFactors = FALSE

  )

}


#' @keywords internal

#' @noRd

.coverage_gamlss2_eta_estimates <- function(fit) {

  if (inherits(fit, "error") || is.null(fit)) {

    return(data.frame(parameter = character(0), term = character(0), estimate_eta = numeric(0)))

  }

  cf <- tryCatch(stats::coef(fit), error = function(e) numeric(0))

  if (length(cf) == 0L) {

    return(data.frame(parameter = character(0), term = character(0), estimate_eta = numeric(0)))

  }

  nms <- names(cf)

  parameter <- sub("^([^.]+)\\.p\\..*$", "\\1", nms)

  term <- sub("^[^.]+\\.p\\.", "", nms)

  term <- ifelse(term %in% c("(Intercept)", "intercept"), "intercept", term)

  data.frame(

    parameter = parameter,

    term = term,

    estimate_eta = as.numeric(cf),

    stringsAsFactors = FALSE

  )

}


#' @keywords internal

#' @noRd

.coverage_eta_error_class <- function(term, abs_error) {

  if (!is.finite(abs_error)) return("missing")

  if (identical(term, "intercept")) {

    if (abs_error <= 0.25) return("acceptable")

    if (abs_error <= 0.50) return("review")

    return("concern")

  }

  if (abs_error <= 0.15) return("acceptable")

  if (abs_error <= 0.35) return("review")

  "concern"

}


#' @keywords internal

#' @noRd

.coverage_parameter_results <- function(estimates, truth) {

  out <- merge(truth, estimates, by = c("parameter", "term"), all = TRUE)

  out$eta_error <- out$estimate_eta - out$true_eta

  out$abs_eta_error <- abs(out$eta_error)

  out$eta_error_class <- vapply(seq_len(nrow(out)), function(i) {

    .coverage_eta_error_class(out$term[[i]], out$abs_eta_error[[i]])

  }, character(1))

  out[order(out$parameter, out$term), , drop = FALSE]

}


#' @keywords internal

#' @noRd

.coverage_smooth_eta_recovery <- function(fit, dat, copula) {

  empty <- c(

    smooth_eta_rmse = NA_real_,

    smooth_eta_max_abs_error = NA_real_,

    smooth_eta_n = 0

  )

  if (!inherits(fit, "gamlss.longitudinal")) {

    return(empty)

  }

  eta_out <- tryCatch(

    calc_eta(

      par_cov = fit$par,

      mm = fit$model_matrix,

      margin_dist = fit$margin_dist,

      copula_link = get_copula_dist(fit$copula_dist)$copula_link,

      par_s = fit$par_s

    ),

    error = function(e) NULL

  )

  if (is.null(eta_out) || is.null(eta_out$eta)) {

    return(empty)

  }


  errors <- numeric(0)

  add_errors <- function(parameter, truth) {

    if (!parameter %in% names(eta_out$eta)) return()

    estimate <- as.numeric(eta_out$eta[[parameter]])

    truth <- as.numeric(truth)

    ok <- is.finite(truth)

    if (!any(ok)) return()

    truth <- truth[ok]

    if (length(estimate) == length(ok)) {

      estimate <- estimate[ok]

    } else if (length(estimate) != length(truth)) {

      return()

    }

    valid <- is.finite(estimate) & is.finite(truth)

    if (any(valid)) {

      errors <<- c(errors, estimate[valid] - truth[valid])

    }

  }


  for (parameter in names(fit$margin_dist$parameters)) {

    truth_col <- paste0("true_", parameter)

    if (!truth_col %in% names(dat)) next

    linkfun <- fit$margin_dist[[paste0(parameter, ".linkfun")]]

    if (!is.function(linkfun)) next

    truth_eta <- tryCatch(as.numeric(linkfun(dat[[truth_col]])), error = function(e) NULL)

    if (!is.null(truth_eta)) add_errors(parameter, truth_eta)

  }


  if ("true_theta" %in% names(dat)) {

    theta_link <- get_copula_dist(copula)$copula_link$theta.linkfun

    truth_eta <- tryCatch(as.numeric(theta_link(dat$true_theta)), error = function(e) NULL)

    if (!is.null(truth_eta)) add_errors("theta", truth_eta)

  }

  if (identical(copula, "t") && "true_zeta" %in% names(dat)) {

    zeta_link <- get_copula_dist(copula)$copula_link$zeta.linkfun

    truth_eta <- tryCatch(as.numeric(zeta_link(dat$true_zeta)), error = function(e) NULL)

    if (!is.null(truth_eta)) add_errors("zeta", truth_eta)

  }


  if (length(errors) == 0L) {

    return(empty)

  }

  c(

    smooth_eta_rmse = sqrt(mean(errors^2, na.rm = TRUE)),

    smooth_eta_max_abs_error = max(abs(errors), na.rm = TRUE),

    smooth_eta_n = length(errors)

  )

}


#' @keywords internal

#' @noRd

.coverage_gamlss_eta_estimates <- function(fit, margin_dist) {

  if (inherits(fit, "error") || is.null(fit)) {

    return(data.frame(parameter = character(0), term = character(0), estimate_eta = numeric(0)))

  }

  rows <- list()

  for (parameter in names(margin_dist$parameters)) {

    cf <- tryCatch(stats::coef(fit, what = parameter), error = function(e) numeric(0))

    if (length(cf) == 0L) next

    term <- names(cf)

    term <- ifelse(term %in% c("(Intercept)", "intercept"), "intercept", term)

    rows[[length(rows) + 1L]] <- data.frame(

      parameter = parameter,

      term = term,

      estimate_eta = as.numeric(cf),

      stringsAsFactors = FALSE

    )

  }

  if (length(rows) == 0L) {

    return(data.frame(parameter = character(0), term = character(0), estimate_eta = numeric(0)))

  }

  do.call(rbind, rows)

}


#' @keywords internal

#' @noRd
