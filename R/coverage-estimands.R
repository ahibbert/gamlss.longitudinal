# Consolidated support module for reviewer navigation.
# Source blocks are grouped mechanically from the files named below.

# ---- coverage-case-grid.R ----

#' @keywords internal

#' @noRd

.coverage_make_case_grid <- function(
    families = NULL,
    copulas = .coverage_supported_copulas(),
    methods = .coverage_default_methods(),
    designs = "intercept",
    include_mixed = FALSE) {
  catalog <- .coverage_family_catalog(include_mixed = include_mixed)

  if (is.null(families)) {
    families <- catalog$family[catalog$supported]
  }

  unsupported <- catalog[match(families, catalog$family), , drop = FALSE]

  unsupported <- unsupported[is.na(unsupported$family) | !unsupported$supported, , drop = FALSE]

  if (nrow(unsupported) > 0L) {
    bad <- ifelse(is.na(unsupported$family), families[is.na(match(families, catalog$family))], unsupported$family)

    stop("Unsupported coverage family/families: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  bad_methods <- setdiff(methods, .coverage_supported_methods())

  if (length(bad_methods) > 0L) {
    stop("Unsupported coverage method(s): ", paste(bad_methods, collapse = ", "), call. = FALSE)
  }

  bad_designs <- setdiff(designs, .coverage_supported_designs())

  if (length(bad_designs) > 0L) {
    stop("Unsupported coverage design(s): ", paste(bad_designs, collapse = ", "), call. = FALSE)
  }

  grid <- expand.grid(
    family = families,
    copula = copulas,
    method = methods,
    design = designs,
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )

  grid$case_id <- seq_len(nrow(grid))

  grid
}

# ---- coverage-catalog.R ----

#' @keywords internal

#' @noRd

.coverage_supported_copulas <- function() {
  c("N", "C", "F", "G", "J", "t")
}


#' @keywords internal

#' @noRd

.coverage_supported_methods <- function() {
  c("gamlss", "gamlss2", "rs_separate", "rs_joint", "cg", "gee", "glmm", "gam", "glmmTMB")
}


#' @keywords internal

#' @noRd

.coverage_supported_designs <- function() {
  c("intercept", "covariate", "scale", "time_dependence", "smooth")
}


#' @keywords internal

#' @noRd

.coverage_default_methods <- function() {
  c("gamlss", "rs_separate", "rs_joint", "cg")
}


#' @keywords internal

#' @noRd

.coverage_attach_namespace <- function(package) {
  if (!paste0("package:", package) %in% search()) {
    suppressPackageStartupMessages(attachNamespace(package))
  }

  invisible(TRUE)
}


#' @keywords internal

#' @noRd

.coverage_family_overrides <- function() {
  data.frame(
    family = c(
      "BI", "BB", "DBI", "ZABB", "ZABI", "ZIBB", "ZIBI",
      "LG",
      "MN3", "MN4", "MN5"
    ),
    supported = FALSE,
    unsupported_reason = c(
      rep("requires denominator/bounded-binomial response support not yet represented in the coverage likelihood calls", 7L),
      "logarithmic-series family needs family-specific starting/support handling before it can be all-method comparable",
      rep("ordinal/multinomial response support needs a family-specific simulation and likelihood path", 3L)
    ),
    stringsAsFactors = FALSE
  )
}


#' @keywords internal

#' @noRd

.coverage_family_catalog <- function(include_mixed = FALSE) {
  requireNamespace("gamlss.dist", quietly = TRUE)


  objects <- getNamespaceExports("gamlss.dist")

  candidates <- objects[grepl("^[A-Z][A-Z0-9]+$", objects)]

  rows <- lapply(candidates, function(family) {
    family_obj <- tryCatch(do.call(get(family, envir = asNamespace("gamlss.dist")), list()), error = function(e) NULL)

    if (is.null(family_obj) || is.null(family_obj$family) || is.null(family_obj$parameters)) {
      return(NULL)
    }

    family_name <- as.character(family_obj$family[1])

    q_name <- paste0("q", family_name)

    p_name <- paste0("p", family_name)

    d_name <- paste0("d", family_name)

    has_qpd <- all(vapply(c(q_name, p_name, d_name), exists, logical(1),
      envir = asNamespace("gamlss.dist"), inherits = FALSE
    ))

    type <- paste(as.character(family_obj$type), collapse = " ")

    mixed <- grepl("mixed", type, ignore.case = TRUE)

    supported <- has_qpd && (!mixed || isTRUE(include_mixed))

    unsupported_reason <- if (!has_qpd) {
      "missing q/p/d function"
    } else if (mixed && !isTRUE(include_mixed)) {
      "mixed support families are excluded from the default coverage grid"
    } else {
      NA_character_
    }

    data.frame(
      family = family_name,
      type = type,
      parameters = paste(names(family_obj$parameters), collapse = ","),
      supported = supported,
      unsupported_reason = unsupported_reason,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

  out <- out[!duplicated(out$family), , drop = FALSE]

  overrides <- .coverage_family_overrides()

  override_idx <- match(overrides$family, out$family)

  override_idx <- override_idx[!is.na(override_idx)]

  if (length(override_idx) > 0L) {
    matched <- match(out$family[override_idx], overrides$family)

    out$supported[override_idx] <- overrides$supported[matched]

    out$unsupported_reason[override_idx] <- overrides$unsupported_reason[matched]
  }

  out[order(out$family), , drop = FALSE]
}

# ---- coverage-conditions.R ----

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
  if (isTRUE(success)) {
    return("ok")
  }

  text <- paste(c(error, warnings), collapse = " ")

  if (grepl("max_elapsed_sec|elapsed", text, ignore.case = TRUE) || elapsed >= max_elapsed_sec) {
    return("convergence timeout")
  }

  if (grepl("starting|start", text, ignore.case = TRUE)) {
    return("start failure")
  }

  if (grepl("non-finite|NA/NaN|Inf|infinite", text, ignore.case = TRUE)) {
    return("non-finite likelihood")
  }

  if (grepl("boundary|outside|domain", text, ignore.case = TRUE)) {
    return("boundary estimate")
  }

  if (grepl("derivative|hessian|gradient", text, ignore.case = TRUE)) {
    return("numerical derivative issue")
  }

  "fit error"
}

# ---- coverage-parameter-defaults.R ----

#' @keywords internal

#' @noRd

.coverage_safe_margin_param_overrides <- function(family_name) {
  overrides <- list(
    BI = list(mu = 0.5),
    LG = list(mu = 0.5),
    NO = list(mu = 0, sigma = 1),
    NO2 = list(mu = 0, sigma = 1),
    EXP = list(mu = 1),
    GA = list(mu = 2, sigma = 0.5),
    IG = list(mu = 2, sigma = 0.4),
    IGAMMA = list(mu = 2, sigma = 0.5),
    WEI = list(mu = 2, sigma = 0.8),
    WEI2 = list(mu = 2, sigma = 0.8),
    WEI3 = list(mu = 2, sigma = 0.8),
    PO = list(mu = 3),
    GEOM = list(mu = 2),
    YULE = list(mu = 3, sigma = 2),
    ZIPF = list(mu = 2, sigma = 2),
    NBI = list(mu = 4, sigma = 0.5),
    NBII = list(mu = 4, sigma = 0.5),
    DEL = list(mu = 4, sigma = 0.6, nu = 0.5),
    PIG = list(mu = 4, sigma = 0.5),
    ZIP = list(mu = 4, sigma = 0.2),
    ZIP2 = list(mu = 4, sigma = 0.2),
    ZAP = list(mu = 4, sigma = 0.2, nu = 0.5),
    ZINBI = list(mu = 4, sigma = 0.5, nu = 0.2),
    ZANBI = list(mu = 4, sigma = 0.5, nu = 0.2),
    ZAPIG = list(mu = 4, sigma = 0.5, nu = 0.2),
    ZASICHEL = list(mu = 4, sigma = 0.5, nu = 0.2),
    ZAZIPF = list(mu = 0.5, sigma = 0.2),
    BNB = list(mu = 4, sigma = 0.5, nu = 2),
    RGE = list(mu = 2, sigma = 0.3, nu = 1),
    SI = list(mu = 0.5, sigma = 0.05, nu = -0.5),
    ZABNB = list(mu = 4, sigma = 0.5, nu = 2, tau = 0.2),
    ZIBNB = list(mu = 4, sigma = 0.5, nu = 2, tau = 0.2),
    ZINBF = list(mu = 4, sigma = 0.5, nu = 2, tau = 0.2)
  )

  overrides[[family_name]] %||% list()
}

#' @keywords internal

#' @noRd

.coverage_default_margin_params <- function(margin_dist) {
  family_name <- as.character(margin_dist$family[1])

  qfun <- get(paste0("q", family_name), envir = asNamespace("gamlss.dist"), inherits = FALSE)

  q_formals <- formals(qfun)

  par_names <- names(margin_dist$parameters)

  out <- list()

  for (par_name in par_names) {
    value <- q_formals[[par_name]]

    value <- tryCatch(eval(value, envir = baseenv()), error = function(e) NA_real_)

    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      value <- switch(par_name,
        mu = if (.is_discrete_margin(margin_dist)) 3 else 1,
        sigma = 0.7,
        nu = 0.3,
        tau = 2,
        0
      )
    }

    out[[par_name]] <- as.numeric(value)
  }


  safe_overrides <- .coverage_safe_margin_param_overrides(family_name)

  for (par_name in intersect(names(safe_overrides), names(out))) {
    out[[par_name]] <- safe_overrides[[par_name]]
  }

  if ("sigma" %in% names(out) && .is_discrete_margin(margin_dist)) out$sigma <- min(max(out$sigma, 0.3), 1.2)

  if ("nu" %in% names(out) && family_name %in% c("DEL", "ZIP", "ZAP", "ZINBI")) out$nu <- min(max(out$nu, 0.25), 0.75)

  out
}


#' @keywords internal

#' @noRd

.coverage_copula_params <- function(copula, dependence = c("moderate", "near_independent", "strong")) {
  dependence <- match.arg(dependence)

  tau <- switch(dependence,
    near_independent = 0.05,
    moderate = 0.25,
    strong = 0.55
  )

  rho <- switch(dependence,
    near_independent = 0.05,
    moderate = 0.3,
    strong = 0.65
  )


  if (copula %in% c("N")) {
    list(theta = rho)
  } else if (copula %in% c("t")) {
    list(theta = rho, zeta = 5)
  } else {
    list(tau = tau)
  }
}


#' @keywords internal

#' @noRd

.coverage_simulation_u_bounds <- function(family) {
  fragile <- c(
    "GT", "SEP", "SEP1", "SEP2", "SEP3", "SEP4",
    "SHASH", "SST", "ST1", "ST2", "ST3", "ST4", "ST5",
    "TF", "TF2"
  )

  if (family %in% fragile) c(1e-4, 1 - 1e-4) else NULL
}

# ---- coverage-smooth-designs.R ----

#' @keywords internal

#' @noRd

.coverage_smooth_eta_component <- function(data, amplitude = 0.18) {
  x <- if ("x" %in% names(data)) {
    as.numeric(data$x)
  } else {
    sim_rescale01(as.numeric(data$.sim_subject_index)) - 0.5
  }

  x_scaled <- sim_rescale01(x)

  wave <- sin(2 * pi * x_scaled) + 0.5 * cos(4 * pi * x_scaled)

  wave <- wave - mean(wave, na.rm = TRUE)

  amplitude * wave
}


#' @keywords internal

#' @noRd

.coverage_make_smooth_param <- function(linkfun, linkinv, base_value, amplitude = 0.18) {
  force(linkfun)

  force(linkinv)

  force(base_value)

  force(amplitude)

  function(data) {
    base_eta <- as.numeric(linkfun(base_value))[1L]

    out <- linkinv(base_eta + .coverage_smooth_eta_component(data, amplitude = amplitude))

    as.numeric(out)
  }
}


#' @keywords internal

#' @noRd

.coverage_time_varying_copula_params <- function(copula) {
  if (copula %in% c("N")) {
    return(list(theta = function(edge_data) {
      time_scaled <- sim_rescale01(edge_data$time_left)

      0.1 + 0.45 * time_scaled
    }))
  }

  if (copula %in% c("t")) {
    return(list(
      theta = function(edge_data) {
        time_scaled <- sim_rescale01(edge_data$time_left)

        0.1 + 0.45 * time_scaled
      },
      zeta = 5
    ))
  }

  list(tau = function(edge_data) {
    time_scaled <- sim_rescale01(edge_data$time_left)

    0.1 + 0.3 * time_scaled
  })
}


#' @keywords internal

#' @noRd

# ---- coverage-truth-estimands.R ----

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

    if (!truth_col %in% names(dat)) {
      return(NULL)
    }

    linkfun <- margin_dist[[paste0(par_name, ".linkfun")]]

    if (is.null(linkfun)) {
      return(NULL)
    }

    ok <- is.finite(dat[[truth_col]])

    if (!any(ok)) {
      return(NULL)
    }

    eta <- as.numeric(linkfun(dat[[truth_col]][ok]))

    if (!all(is.finite(eta))) {
      return(NULL)
    }

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

# ---- coverage-eta-estimates.R ----

#' @keywords internal

#' @noRd

.coverage_natural_estimates <- function(fit) {
  if (!inherits(fit, "gamlss.longitudinal")) {
    return(stats::setNames(numeric(0), character(0)))
  }

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

# ---- coverage-eta-recovery.R ----

#' @keywords internal

#' @noRd

.coverage_eta_error_class <- function(term, abs_error) {
  if (!is.finite(abs_error)) {
    return("missing")
  }

  if (identical(term, "intercept")) {
    if (abs_error <= 0.25) {
      return("acceptable")
    }

    if (abs_error <= 0.50) {
      return("review")
    }

    return("concern")
  }

  if (abs_error <= 0.15) {
    return("acceptable")
  }

  if (abs_error <= 0.35) {
    return("review")
  }

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
    if (!parameter %in% names(eta_out$eta)) {
      return()
    }

    estimate <- as.numeric(eta_out$eta[[parameter]])

    truth <- as.numeric(truth)

    ok <- is.finite(truth)

    if (!any(ok)) {
      return()
    }

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

# ---- coverage-review-metrics.R ----

#' @keywords internal
#' @noRd
.coverage_margin_class <- function(method, gap_pct) {
  if (!is.finite(gap_pct)) {
    return("missing")
  }

  if (method %in% c("gamlss", "gamlss2")) {
    return("reference")
  }

  if (identical(method, "rs_separate")) {
    if (gap_pct <= 1) {
      return("excellent")
    }

    if (gap_pct <= 3) {
      return("acceptable")
    }

    return("review")
  }

  if (gap_pct <= 3) {
    return("excellent")
  }

  if (gap_pct <= 7.5) {
    return("acceptable")
  }

  "review"
}


#' @keywords internal
#' @noRd
.coverage_joint_class <- function(method, delta_pct) {
  if (identical(method, "gamlss2")) {
    return("not_applicable")
  }

  if (!is.finite(delta_pct)) {
    return("missing")
  }

  if (delta_pct >= -1) {
    return("acceptable")
  }

  "review"
}


#' @keywords internal
#' @noRd
.coverage_add_review_metrics <- function(results) {
  results$gamlss2_marginal_loglik <- NA_real_

  results$margin_gap_pct_vs_gamlss2 <- NA_real_

  results$reference_marginal_method <- NA_character_

  results$reference_marginal_loglik <- NA_real_

  results$margin_gap_pct_vs_reference <- NA_real_

  results$margin_review_class <- "missing"

  results$rs_separate_joint_loglik <- NA_real_

  results$joint_delta_vs_rs_separate <- NA_real_

  results$joint_delta_pct_vs_rs_separate <- NA_real_

  results$joint_review_class <- "missing"


  group_cols <- c("family", "copula", "design", "n_subject", "n_time", "dependence", "missingness", "start_mode")

  group_key <- interaction(results[group_cols], drop = TRUE, lex.order = TRUE)


  for (key in unique(group_key)) {
    idx <- which(group_key == key)

    group <- results[idx, , drop = FALSE]

    gamlss2_ll <- group$marginal_loglik[group$method == "gamlss2" & group$success][1]

    gamlss_ll <- group$marginal_loglik[group$method == "gamlss" & group$success][1]

    gamlss2_fit_method <- group$marginal_fit_method[group$method == "gamlss2" & group$success][1]

    gamlss_fit_method <- group$marginal_fit_method[group$method == "gamlss" & group$success][1]

    rs_joint_ll <- group$joint_loglik[group$method == "rs_separate" & group$success][1]


    if (length(gamlss2_ll) == 1L && is.finite(gamlss2_ll)) {
      results$gamlss2_marginal_loglik[idx] <- gamlss2_ll

      results$margin_gap_pct_vs_gamlss2[idx] <- 100 * (gamlss2_ll - results$marginal_loglik[idx]) / pmax(abs(gamlss2_ll), 1e-8)
    }

    reference_ll <- NA_real_

    reference_method <- NA_character_

    if (length(gamlss2_ll) == 1L && is.finite(gamlss2_ll)) {
      reference_ll <- gamlss2_ll

      reference_method <- if (length(gamlss2_fit_method) == 1L && !is.na(gamlss2_fit_method)) gamlss2_fit_method else "gamlss2"
    } else if (length(gamlss_ll) == 1L && is.finite(gamlss_ll)) {
      reference_ll <- gamlss_ll

      reference_method <- if (length(gamlss_fit_method) == 1L && !is.na(gamlss_fit_method)) gamlss_fit_method else "gamlss"
    }

    if (is.finite(reference_ll)) {
      results$reference_marginal_method[idx] <- reference_method

      results$reference_marginal_loglik[idx] <- reference_ll

      results$margin_gap_pct_vs_reference[idx] <- 100 * (reference_ll - results$marginal_loglik[idx]) / pmax(abs(reference_ll), 1e-8)
    }

    if (length(rs_joint_ll) == 1L && is.finite(rs_joint_ll)) {
      results$rs_separate_joint_loglik[idx] <- rs_joint_ll

      results$joint_delta_vs_rs_separate[idx] <- results$joint_loglik[idx] - rs_joint_ll

      results$joint_delta_pct_vs_rs_separate[idx] <- 100 * results$joint_delta_vs_rs_separate[idx] / pmax(abs(rs_joint_ll), 1e-8)
    }
  }


  results$margin_review_class <- vapply(seq_len(nrow(results)), function(i) {
    .coverage_margin_class(results$method[[i]], results$margin_gap_pct_vs_reference[[i]])
  }, character(1))

  results$joint_review_class <- vapply(seq_len(nrow(results)), function(i) {
    .coverage_joint_class(results$method[[i]], results$joint_delta_pct_vs_rs_separate[[i]])
  }, character(1))

  results
}


#' @keywords internal
#' @noRd
.coverage_runtime_summary <- function(results) {
  stats::aggregate(

    elapsed_sec ~ family + copula + method + design + success + failure_type,
    data = results,
    FUN = function(x) c(n = length(x), median = stats::median(x), max = max(x))
  )
}
