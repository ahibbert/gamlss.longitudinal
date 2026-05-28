#' @keywords internal
#' @noRd
.coverage_supported_copulas <- function() {
  c("N", "C", "F", "G", "J", "t")
}

#' @keywords internal
#' @noRd
.coverage_supported_methods <- function() {
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
    ZAZIPF = list(mu = 0.5, sigma = 0.2)
  )
  overrides[[family_name]] %||% list()
}

#' @keywords internal
#' @noRd
.coverage_make_case_grid <- function(
  families = NULL,
  copulas = .coverage_supported_copulas(),
  methods = .coverage_supported_methods(),
  designs = "intercept",
  include_mixed = FALSE
) {
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
.coverage_simulate_case <- function(
  family,
  copula,
  design = c("intercept", "covariate", "smooth"),
  n = 80,
  times = 1:3,
  seed = 1,
  dependence = "moderate"
) {
  design <- match.arg(design)
  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
  margin_params <- .coverage_default_margin_params(margin_dist)
  covariates <- NULL

  if (design %in% c("covariate", "smooth")) {
    covariates <- function(base) {
      x <- sim_rescale01(as.numeric(base$.sim_subject_index))
      data.frame(x = x - mean(x), stringsAsFactors = FALSE)
    }
    if ("mu" %in% names(margin_params)) {
      base_mu <- margin_params$mu
      margin_params$mu <- function(data) {
        x <- sim_rescale01(as.numeric(data$.sim_subject_index))
        x <- x - mean(x)
        if (base_mu > 0) {
          pmax(base_mu * exp(0.15 * x), .Machine$double.eps)
        } else {
          base_mu + 0.15 * x
        }
      }
    }
  }

  simulate_longitudinal_dataset(
    n = n,
    times = times,
    margin_dist = margin_dist,
    copula_dist = copula,
    margin_params = margin_params,
    copula_params = .coverage_copula_params(copula, dependence = dependence),
    covariates = covariates,
    seed = seed,
    include_truth = TRUE
  )
}

#' @keywords internal
#' @noRd
.coverage_apply_missingness <- function(dat, missingness = c("none", "mcar", "drop_rows"), prop = 0.05) {
  missingness <- match.arg(missingness)
  if (identical(missingness, "none")) {
    return(dat)
  }

  n_subject <- length(unique(dat$subject))
  n_time <- length(unique(dat$time))
  n_target <- max(1L, floor(nrow(dat) * prop))

  if (identical(missingness, "mcar")) {
    eligible <- which(dat$time != min(dat$time))
    idx <- eligible[seq_len(min(length(eligible), n_target))]
    dat$response[idx] <- NA_real_
    return(dat)
  }

  subject_index <- as.integer(dat$subject)
  time_values <- sort(unique(dat$time))
  drop_idx <- subject_index <= max(1L, floor(n_subject * prop)) & dat$time == time_values[min(2L, n_time)]
  dat[!drop_idx, , drop = FALSE]
}

#' @keywords internal
#' @noRd
.coverage_fit_formulas <- function(design) {
  if (identical(design, "covariate")) {
    list(mu = response ~ x, sigma = ~1, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)
  } else if (identical(design, "smooth")) {
    list(mu = response ~ s(x, bs = "ps"), sigma = ~1, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)
  } else {
    list(mu = response ~ 1, sigma = ~1, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)
  }
}

#' @keywords internal
#' @noRd
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
    if (identical(design, "covariate") && "x" %in% names(dat)) {
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
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = "theta",
        term = "intercept",
        true_eta = mean(as.numeric(theta_link(dat$true_theta[theta_ok]))),
        stringsAsFactors = FALSE
      )
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
.coverage_fit_gamlss <- function(dat, family, copula, design, max_inner_iter = 8, max_elapsed_sec = Inf) {
  .coverage_attach_namespace("gamlss")
  .coverage_attach_namespace("gamlss.dist")
  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
  formulas <- .coverage_fit_formulas(design)
  fml <- formulas$mu
  fit_vars <- unique(unlist(lapply(formulas[names(formulas) %in% names(margin_dist$parameters)], all.vars)))
  fit_vars <- fit_vars[fit_vars %in% names(dat)]
  fit_dat <- dat[, fit_vars, drop = FALSE]
  start <- Sys.time()
  captured <- .coverage_capture_conditions({
    fit <- NULL
    invisible(utils::capture.output({
      fit <- gamlss::gamlss(
        formula = fml,
        sigma.formula = formulas$sigma,
        nu.formula = formulas$nu,
        tau.formula = formulas$tau,
        data = fit_dat,
        family = margin_dist,
        control = gamlss::gamlss.control(n.cyc = max_inner_iter, trace = FALSE)
      )
    }))
    fit
  })
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  success <- !inherits(captured$value, "error") && is.finite(as.numeric(stats::logLik(captured$value)))
  converged <- success
  if (success && !is.null(captured$value$converged)) {
    converged <- isTRUE(captured$value$converged)
  }
  out <- data.frame(
    method = "gamlss",
    success = success,
    converged = converged,
    failure_type = .coverage_taxonomy(success, if (inherits(captured$value, "error")) conditionMessage(captured$value) else NA_character_, captured$warnings, elapsed, max_elapsed_sec),
    error = if (inherits(captured$value, "error")) conditionMessage(captured$value) else NA_character_,
    warnings = paste(unique(captured$warnings), collapse = " | "),
    marginal_loglik = if (success) as.numeric(stats::logLik(captured$value)) else NA_real_,
    copula_loglik = NA_real_,
    joint_loglik = NA_real_,
    elapsed_sec = elapsed,
    max_abs_param_error = NA_real_,
    max_rel_param_error = NA_real_,
    fitted_copula_tau = NA_real_,
    true_copula_tau = NA_real_,
    stringsAsFactors = FALSE
  )
  truth_eta <- .coverage_true_eta_coefficients(dat, family, copula, design)
  truth_eta <- truth_eta[truth_eta$parameter %in% names(margin_dist$parameters), , drop = FALSE]
  attr(out, "parameter_results") <- .coverage_parameter_results(
    .coverage_gamlss_eta_estimates(captured$value, margin_dist),
    truth_eta
  )
  out
}

#' @keywords internal
#' @noRd
.coverage_fit_gamlss2 <- function(dat, family, copula, design, max_elapsed_sec = Inf) {
  if (!requireNamespace("gamlss2", quietly = TRUE)) {
    stop("Package 'gamlss2' is required for method = 'gamlss2'.", call. = FALSE)
  }
  gamlss2_fit <- getExportedValue("gamlss2", "gamlss2")
  gamlss2_control <- getExportedValue("gamlss2", "gamlss2_control")
  .coverage_attach_namespace("gamlss.dist")
  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
  fml <- if (identical(design, "intercept")) response ~ 1 else response ~ x
  start <- Sys.time()
  captured <- .coverage_capture_conditions({
    fit <- gamlss2_fit(
      fml,
      data = dat,
      family = margin_dist,
      control = gamlss2_control(trace = FALSE)
    )
    fit
  })
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  success <- !inherits(captured$value, "error") && is.finite(as.numeric(stats::logLik(captured$value)))
  out <- data.frame(
    method = "gamlss2",
    success = success,
    converged = success,
    failure_type = .coverage_taxonomy(success, if (inherits(captured$value, "error")) conditionMessage(captured$value) else NA_character_, captured$warnings, elapsed, max_elapsed_sec),
    error = if (inherits(captured$value, "error")) conditionMessage(captured$value) else NA_character_,
    warnings = paste(unique(captured$warnings), collapse = " | "),
    marginal_loglik = if (success) as.numeric(stats::logLik(captured$value)) else NA_real_,
    copula_loglik = NA_real_,
    joint_loglik = NA_real_,
    elapsed_sec = elapsed,
    max_abs_param_error = NA_real_,
    max_rel_param_error = NA_real_,
    fitted_copula_tau = NA_real_,
    true_copula_tau = NA_real_,
    stringsAsFactors = FALSE
  )
  truth_eta <- .coverage_true_eta_coefficients(dat, family, copula, design)
  truth_eta <- truth_eta[truth_eta$parameter %in% names(margin_dist$parameters), , drop = FALSE]
  attr(out, "parameter_results") <- .coverage_parameter_results(
    .coverage_gamlss2_eta_estimates(captured$value),
    truth_eta
  )
  out
}

#' @keywords internal
#' @noRd
.coverage_fit_longitudinal <- function(
  dat,
  family,
  copula,
  method,
  design,
  max_outer_iter = 8,
  max_inner_iter = 8,
  max_elapsed_sec = 20,
  start_from = NA
) {
  .coverage_attach_namespace("gamlss")
  .coverage_attach_namespace("gamlss.dist")
  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
  formulas <- .coverage_fit_formulas(design)
  fit_method <- if (identical(method, "cg")) "CG" else "RS"
  include_dlcopdpar <- identical(method, "rs_joint") || identical(method, "cg")
  start <- Sys.time()
  captured <- .coverage_capture_conditions({
    fit <- NULL
    invisible(utils::capture.output({
      fit <- gamlss.longitudinal(
        dataset = dat,
        margin_dist = margin_dist,
        copula_dist = copula,
        time_var = "time",
        subject_var = "subject",
        mu.formula = formulas$mu,
        sigma.formula = formulas$sigma,
        nu.formula = formulas$nu,
        tau.formula = formulas$tau,
        theta.formula = formulas$theta,
        zeta.formula = formulas$zeta,
        include_dlcopdpar = include_dlcopdpar,
        method = fit_method,
        start_from = start_from,
        warm_start_joint = isTRUE(include_dlcopdpar) && all(is.na(start_from)),
        warm_start_joint_iter = min(5L, max_outer_iter),
        max_outer_iter = max_outer_iter,
        max_inner_iter = max_inner_iter,
        outer_stop_crit = 1e-4,
        inner_stop_crit = 1e-4,
        max_elapsed_sec = max_elapsed_sec,
        compute_vcov = FALSE,
        verbose = 0
      )
    }))
    fit
  })
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  fit <- captured$value
  is_fit <- inherits(fit, "gamlss.longitudinal")
  loglik <- if (is_fit) fit$calc_lik_out_end$log_lik else c(marginal = NA_real_, copula = NA_real_, joint = NA_real_)
  estimates <- if (is_fit) .coverage_natural_estimates(fit) else stats::setNames(numeric(0), character(0))
  truth <- .coverage_truth_summary(dat, copula)
  common <- intersect(names(estimates), names(truth))
  abs_error <- if (length(common)) max(abs(estimates[common] - truth[common]), na.rm = TRUE) else NA_real_
  rel_error <- if (length(common)) {
    denom <- pmax(abs(truth[common]), 1e-8)
    max(abs(estimates[common] - truth[common]) / denom, na.rm = TRUE)
  } else {
    NA_real_
  }
  theta_tau <- if (is_fit && "theta" %in% names(estimates)) {
    zeta <- if ("zeta" %in% names(estimates)) estimates[["zeta"]] else 0
    .copula_par_to_tau(copula, estimates[["theta"]], zeta)
  } else {
    NA_real_
  }
  success <- is_fit && all(is.finite(as.numeric(loglik))) && !inherits(fit, "error")
  out <- data.frame(
    method = method,
    success = success,
    converged = if (is_fit && !is.null(fit$convergence$converged)) isTRUE(fit$convergence$converged) else FALSE,
    failure_type = .coverage_taxonomy(success, if (inherits(fit, "error")) conditionMessage(fit) else NA_character_, captured$warnings, elapsed, max_elapsed_sec),
    error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
    warnings = paste(unique(captured$warnings), collapse = " | "),
    marginal_loglik = as.numeric(loglik["marginal"]),
    copula_loglik = as.numeric(loglik["copula"]),
    joint_loglik = as.numeric(loglik["joint"]),
    elapsed_sec = elapsed,
    max_abs_param_error = abs_error,
    max_rel_param_error = rel_error,
    fitted_copula_tau = theta_tau,
    true_copula_tau = if ("copula_tau" %in% names(truth)) truth[["copula_tau"]] else NA_real_,
    stringsAsFactors = FALSE
  )
  attr(out, "parameter_results") <- .coverage_parameter_results(
    .coverage_longitudinal_eta_estimates(fit),
    .coverage_true_eta_coefficients(dat, family, copula, design)
  )
  out
}

#' @keywords internal
#' @noRd
.coverage_run_case <- function(
  case,
  n = 80,
  times = 1:3,
  seed = 1,
  max_outer_iter = 8,
  max_inner_iter = 8,
  max_elapsed_sec = 20,
  dependence = "moderate",
  missingness = "none",
  start_mode = c("default", "truth_adjacent")
) {
  start_mode <- match.arg(start_mode)
  dat <- .coverage_simulate_case(
    family = case$family,
    copula = case$copula,
    design = case$design,
    n = n,
    times = times,
    seed = seed,
    dependence = dependence
  )
  dat <- .coverage_apply_missingness(dat, missingness = missingness)
  start_from <- if (identical(start_mode, "truth_adjacent") && !case$method %in% c("gamlss", "gamlss2")) {
    .coverage_truth_adjacent_start(dat, case$family, case$copula, case$design)
  } else {
    NA
  }
  fit_row <- if (identical(case$method, "gamlss")) {
    .coverage_fit_gamlss(
      dat,
      case$family,
      case$copula,
      case$design,
      max_inner_iter = max_inner_iter,
      max_elapsed_sec = max_elapsed_sec
    )
  } else if (identical(case$method, "gamlss2")) {
    .coverage_fit_gamlss2(dat, case$family, case$copula, case$design, max_elapsed_sec = max_elapsed_sec)
  } else {
    .coverage_fit_longitudinal(
      dat,
      family = case$family,
      copula = case$copula,
      method = case$method,
      design = case$design,
      max_outer_iter = max_outer_iter,
      max_inner_iter = max_inner_iter,
      max_elapsed_sec = max_elapsed_sec,
      start_from = start_from
    )
  }
  context <- data.frame(
    case_id = case$case_id %||% NA_integer_,
    family = case$family,
    copula = case$copula,
    design = case$design,
    n_subject = n,
    n_time = length(times),
    dependence = dependence,
    missingness = missingness,
    start_mode = start_mode,
    stringsAsFactors = FALSE
  )
  out <- cbind(
    context,
    fit_row,
    row.names = NULL
  )
  parameter_results <- attr(fit_row, "parameter_results")
  if (!is.null(parameter_results) && nrow(parameter_results) > 0L) {
    attr(out, "parameter_results") <- cbind(
      context[rep(1L, nrow(parameter_results)), , drop = FALSE],
      data.frame(
        method = fit_row$method,
        success = fit_row$success,
        elapsed_sec = fit_row$elapsed_sec,
        stringsAsFactors = FALSE
      ),
      parameter_results,
      row.names = NULL
    )
  }
  out
}

#' @keywords internal
#' @noRd
.coverage_run_grid <- function(
  grid,
  n = 80,
  times = 1:3,
  seed = 1,
  max_outer_iter = 8,
  max_inner_iter = 8,
  max_elapsed_sec = 20,
  dependence = "moderate",
  missingness = "none",
  start_mode = "default"
) {
  rows <- vector("list", nrow(grid))
  parameter_rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    case_result <- .coverage_run_case(
      grid[i, , drop = FALSE],
      n = n,
      times = times,
      seed = seed + i,
      max_outer_iter = max_outer_iter,
      max_inner_iter = max_inner_iter,
      max_elapsed_sec = max_elapsed_sec,
      dependence = dependence,
      missingness = missingness,
      start_mode = start_mode
    )
    parameter_rows[[i]] <- attr(case_result, "parameter_results")
    rows[[i]] <- case_result
  }
  results <- .coverage_add_review_metrics(do.call(rbind, rows))
  parameter_rows <- parameter_rows[!vapply(parameter_rows, is.null, logical(1))]
  attr(results, "parameter_results") <- if (length(parameter_rows) == 0L) {
    data.frame()
  } else {
    do.call(rbind, parameter_rows)
  }
  attr(results, "runtime_summary") <- .coverage_runtime_summary(results)
  results
}

#' @keywords internal
#' @noRd
.coverage_margin_class <- function(method, gap_pct) {
  if (!is.finite(gap_pct)) return("missing")
  if (method %in% c("gamlss", "gamlss2")) return("reference")
  if (identical(method, "rs_separate")) {
    if (gap_pct <= 1) return("excellent")
    if (gap_pct <= 3) return("acceptable")
    return("review")
  }
  if (gap_pct <= 3) return("excellent")
  if (gap_pct <= 7.5) return("acceptable")
  "review"
}

#' @keywords internal
#' @noRd
.coverage_joint_class <- function(method, delta_pct) {
  if (identical(method, "gamlss2")) return("not_applicable")
  if (!is.finite(delta_pct)) return("missing")
  if (delta_pct >= -1) return("acceptable")
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
    rs_joint_ll <- group$joint_loglik[group$method == "rs_separate" & group$success][1]

    if (length(gamlss2_ll) == 1L && is.finite(gamlss2_ll)) {
      results$gamlss2_marginal_loglik[idx] <- gamlss2_ll
      results$margin_gap_pct_vs_gamlss2[idx] <- 100 * (gamlss2_ll - results$marginal_loglik[idx]) / pmax(abs(gamlss2_ll), 1e-8)
    }
    reference_ll <- NA_real_
    reference_method <- NA_character_
    if (length(gamlss2_ll) == 1L && is.finite(gamlss2_ll)) {
      reference_ll <- gamlss2_ll
      reference_method <- "gamlss2"
    } else if (length(gamlss_ll) == 1L && is.finite(gamlss_ll)) {
      reference_ll <- gamlss_ll
      reference_method <- "gamlss"
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

#' Run opt-in distribution/copula/method coverage simulations
#'
#' @param families Optional family names. Defaults to all supported non-mixed
#'   `gamlss.dist` families with `q`, `p`, and `d` functions.
#' @param copulas Copula codes.
#' @param methods Fit methods. Defaults to `"gamlss"`, `"rs_separate"`,
#'   `"rs_joint"`, and `"cg"`. Method `"gamlss2"` is also supported when the
#'   optional non-CRAN package is installed and explicitly requested.
#' @param designs Simulation designs: `"intercept"`, `"covariate"`, or
#'   `"smooth"`.
#' @param include_mixed Logical; include mixed-support `gamlss.dist` families in
#'   the candidate family grid.
#' @param output_dir Directory for CSV/RDS outputs.
#' @param write_results Write result files when `TRUE`.
#' @param ... Passed to the grid runner, e.g. `n`, `times`, `max_outer_iter`,
#'   `dependence`, `missingness`, and `start_mode`.
#'
#' @return A data frame of per-fit results. The return value also carries
#'   `"parameter_results"` and `"runtime_summary"` attributes.
#' @export
run_coverage_simulations <- function(
  families = NULL,
  copulas = c("N", "C", "F", "G", "J", "t"),
  methods = .coverage_supported_methods(),
  designs = "intercept",
  include_mixed = FALSE,
  output_dir = file.path("results", "coverage_simulations"),
  write_results = TRUE,
  ...
) {
  grid <- .coverage_make_case_grid(
    families = families,
    copulas = copulas,
    methods = methods,
    designs = designs,
    include_mixed = include_mixed
  )
  results <- .coverage_run_grid(grid, ...)
  if (isTRUE(write_results)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    parameter_results <- attr(results, "parameter_results")
    runtime_summary <- attr(results, "runtime_summary")
    utils::write.csv(results, file.path(output_dir, paste0("coverage_results_", stamp, ".csv")), row.names = FALSE)
    saveRDS(results, file.path(output_dir, paste0("coverage_results_", stamp, ".rds")))
    if (!is.null(parameter_results) && nrow(parameter_results) > 0L) {
      utils::write.csv(
        parameter_results,
        file.path(output_dir, paste0("coverage_parameter_results_", stamp, ".csv")),
        row.names = FALSE
      )
    }
    if (!is.null(runtime_summary) && nrow(runtime_summary) > 0L) {
      utils::write.csv(
        runtime_summary,
        file.path(output_dir, paste0("coverage_runtime_summary_", stamp, ".csv")),
        row.names = FALSE
      )
    }
  }
  results
}
