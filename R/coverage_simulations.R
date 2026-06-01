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
.coverage_make_case_grid <- function(
  families = NULL,
  copulas = .coverage_supported_copulas(),
  methods = .coverage_default_methods(),
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
.coverage_simulate_case <- function(
  family,
  copula,
  design = .coverage_supported_designs(),
  n = 80,
  times = 1:3,
  seed = 1,
  dependence = "moderate"
) {
  design <- match.arg(design)
  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
  margin_params <- .coverage_default_margin_params(margin_dist)
  covariates <- NULL

  if (design %in% c("covariate", "scale", "smooth")) {
    covariates <- function(base) {
      x <- sim_rescale01(as.numeric(base$.sim_subject_index))
      data.frame(x = x - mean(x), stringsAsFactors = FALSE)
    }
    if (identical(design, "covariate") && "mu" %in% names(margin_params)) {
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
    if (identical(design, "smooth")) {
      for (par_name in names(margin_params)) {
        linkfun <- margin_dist[[paste0(par_name, ".linkfun")]]
        linkinv <- margin_dist[[paste0(par_name, ".linkinv")]]
        if (!is.function(linkfun) || !is.function(linkinv)) next
        margin_params[[par_name]] <- .coverage_make_smooth_param(
          linkfun,
          linkinv,
          margin_params[[par_name]],
          amplitude = if (identical(par_name, "mu")) 0.22 else 0.14
        )
      }
    }
    if (identical(design, "scale") && "sigma" %in% names(margin_params)) {
      base_sigma <- margin_params$sigma
      margin_params$sigma <- function(data) {
        x <- sim_rescale01(as.numeric(data$.sim_subject_index))
        x <- x - mean(x)
        pmax(base_sigma * exp(0.45 * x), .Machine$double.eps)
      }
    }
  }

  copula_params <- if (identical(design, "time_dependence")) {
    .coverage_time_varying_copula_params(copula)
  } else if (identical(design, "smooth")) {
    copula_link <- get_copula_dist(copula)$copula_link
    base <- .coverage_copula_params(copula, dependence = dependence)
    if ("tau" %in% names(base)) {
      list(tau = function(edge_data) {
        eta <- stats::qlogis(base$tau) + .coverage_smooth_eta_component(edge_data, amplitude = 0.12)
        stats::plogis(eta)
      })
    } else {
      base$theta <- .coverage_make_smooth_param(copula_link$theta.linkfun, copula_link$theta.linkinv, base$theta, amplitude = 0.12)
      base
    }
  } else {
    .coverage_copula_params(copula, dependence = dependence)
  }

  simulate_longitudinal_dataset(
    n = n,
    times = times,
    margin_dist = margin_dist,
    copula_dist = copula,
    margin_params = margin_params,
    copula_params = copula_params,
    covariates = covariates,
    seed = seed,
    include_truth = TRUE,
    u_bounds = .coverage_simulation_u_bounds(family)
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
  } else if (identical(design, "scale")) {
    list(mu = response ~ 1, sigma = ~x, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)
  } else if (identical(design, "time_dependence")) {
    list(mu = response ~ 1, sigma = ~1, nu = ~1, tau = ~1, theta = ~time, zeta = ~1)
  } else if (identical(design, "smooth")) {
    list(
      mu = response ~ s(x, bs = "ps", k = 6),
      sigma = ~s(x, bs = "ps", k = 6),
      nu = ~s(x, bs = "ps", k = 6),
      tau = ~s(x, bs = "ps", k = 6),
      theta = ~s(x, bs = "ps", k = 6),
      zeta = ~s(x, bs = "ps", k = 6)
    )
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
  success <- !inherits(captured$value, "error") &&
    is.finite(as.numeric(stats::logLik(captured$value)))
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
    marginal_fit_method = "gamlss",
    marginal_fallback_used = FALSE,
    marginal_fallback_error = NA_character_,
    fit_attempt = "gamlss",
    fit_attempt_count = 1L,
    fit_attempt_trace = "gamlss",
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
  .coverage_attach_namespace("gamlss.dist")
  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
  fml <- if (design %in% c("intercept", "scale", "time_dependence")) response ~ 1 else response ~ x
  fit_vars <- unique(all.vars(fml))
  fit_vars <- fit_vars[fit_vars %in% names(dat)]
  fit_dat <- dat[, fit_vars, drop = FALSE]
  start <- Sys.time()
  fallback_used <- FALSE
  fallback_error <- NA_character_
  gamlss2_available <- requireNamespace("gamlss2", quietly = TRUE)
  captured <- if (gamlss2_available) {
    gamlss2_fit <- getExportedValue("gamlss2", "gamlss2")
    gamlss2_control <- getExportedValue("gamlss2", "gamlss2_control")
    .coverage_capture_conditions({
      fit <- gamlss2_fit(
        fml,
        data = fit_dat,
        family = margin_dist,
        control = gamlss2_control(trace = FALSE)
      )
      fit
    })
  } else {
    list(value = simpleError("Package 'gamlss2' is required for method = 'gamlss2'."), warnings = character(0))
  }

  success <- !inherits(captured$value, "error") &&
    is.finite(as.numeric(stats::logLik(captured$value)))
  if (!success) {
    fallback_error <- if (inherits(captured$value, "error")) conditionMessage(captured$value) else "gamlss2 returned non-finite logLik"
    fallback_used <- TRUE
    gamlss_captured <- .coverage_capture_conditions({
      fit <- NULL
      invisible(utils::capture.output({
        fit <- gamlss::gamlss(
          formula = fml,
          data = fit_dat,
          family = margin_dist,
          control = gamlss::gamlss.control(n.cyc = 8, trace = FALSE)
        )
      }))
      fit
    })
    if (!inherits(gamlss_captured$value, "error") && is.finite(as.numeric(stats::logLik(gamlss_captured$value)))) {
      captured <- gamlss_captured
      success <- TRUE
    } else {
      captured$warnings <- c(captured$warnings, gamlss_captured$warnings)
      captured$value <- if (inherits(gamlss_captured$value, "error")) gamlss_captured$value else captured$value
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  success <- success && (!is.finite(max_elapsed_sec) || elapsed <= max_elapsed_sec)
  actual_method <- if (success && isTRUE(fallback_used)) "gamlss" else "gamlss2"
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
    marginal_fit_method = actual_method,
    marginal_fallback_used = fallback_used && success,
    marginal_fallback_error = fallback_error,
    fit_attempt = actual_method,
    fit_attempt_count = if (isTRUE(fallback_used)) 2L else 1L,
    fit_attempt_trace = if (isTRUE(fallback_used)) "gamlss2 > gamlss" else "gamlss2",
    stringsAsFactors = FALSE
  )
  truth_eta <- .coverage_true_eta_coefficients(dat, family, copula, design)
  truth_eta <- truth_eta[truth_eta$parameter %in% names(margin_dist$parameters), , drop = FALSE]
  estimates <- if (identical(actual_method, "gamlss")) {
    .coverage_gamlss_eta_estimates(captured$value, margin_dist)
  } else {
    .coverage_gamlss2_eta_estimates(captured$value)
  }
  attr(out, "parameter_results") <- .coverage_parameter_results(estimates, truth_eta)
  out
}

#' @keywords internal
#' @noRd
.coverage_standard_family <- function(family) {
  if (family %in% c("NO", "NO2")) {
    return(stats::gaussian())
  }
  if (family %in% c("PO", "ZIP", "ZIP2")) {
    return(stats::poisson())
  }
  if (family %in% c("GA", "EXP")) {
    return(stats::Gamma(link = "log"))
  }
  NULL
}

#' @keywords internal
#' @noRd
.coverage_standard_formula <- function(design, comparator) {
  if (design %in% c("intercept", "scale", "time_dependence")) {
    return(response ~ 1)
  }
  if (identical(design, "smooth") && identical(comparator, "gam")) {
    return(response ~ s(x, bs = "ps"))
  }
  response ~ x
}

#' @keywords internal
#' @noRd
.coverage_true_margin_distribution <- function(dat, family, p = 0.9) {
  empty <- list(q = rep(NA_real_, nrow(dat)), cdf = rep(NA_real_, nrow(dat)))
  if (is.null(family)) {
    return(empty)
  }
  margin_dist <- tryCatch(do.call(get(family, envir = asNamespace("gamlss.dist")), list()), error = function(e) NULL)
  if (is.null(margin_dist) || is.null(margin_dist$family)) {
    return(empty)
  }
  family_name <- as.character(margin_dist$family[1])
  qfun <- tryCatch(get(paste0("q", family_name), envir = asNamespace("gamlss.dist"), inherits = FALSE), error = function(e) NULL)
  pfun <- tryCatch(get(paste0("p", family_name), envir = asNamespace("gamlss.dist"), inherits = FALSE), error = function(e) NULL)
  if (is.null(qfun) || is.null(pfun)) {
    return(empty)
  }
  par_names <- names(margin_dist$parameters)
  true_cols <- paste0("true_", par_names)
  if (!all(true_cols %in% names(dat))) {
    return(empty)
  }
  par_args <- stats::setNames(lapply(true_cols, function(nm) dat[[nm]]), par_names)
  q <- tryCatch(do.call(qfun, c(list(p = p), par_args)), error = function(e) rep(NA_real_, nrow(dat)))
  cdf <- tryCatch(do.call(pfun, c(list(q = q), par_args)), error = function(e) rep(NA_real_, nrow(dat)))
  list(q = as.numeric(q), cdf = as.numeric(cdf))
}

#' @keywords internal
#' @noRd
.coverage_comparator_dispersion <- function(y, fitted, family) {
  ok <- is.finite(y) & is.finite(fitted)
  if (sum(ok) < 3L) {
    return(NA_real_)
  }
  if (identical(family$family, "gaussian")) {
    return(stats::sd(y[ok] - fitted[ok]))
  }
  if (identical(family$family, "Gamma")) {
    mu <- pmax(fitted[ok], .Machine$double.eps)
    pearson <- (y[ok] - mu) / mu
    return(sum(pearson^2, na.rm = TRUE) / max(1L, sum(ok) - 1L))
  }
  NA_real_
}

#' @keywords internal
#' @noRd
.coverage_comparator_distribution <- function(y, fitted, family, q, p = 0.9) {
  n <- length(fitted)
  empty <- list(q = rep(NA_real_, n), cdf_at_q = rep(NA_real_, n))
  ok <- is.finite(y) & is.finite(fitted) & is.finite(q)
  if (!any(ok)) {
    return(empty)
  }

  if (identical(family$family, "gaussian")) {
    sigma_hat <- .coverage_comparator_dispersion(y, fitted, family)
    if (!is.finite(sigma_hat) || sigma_hat <= 0) {
      return(empty)
    }
    out_q <- fitted + stats::qnorm(p) * sigma_hat
    cdf <- stats::pnorm(q, mean = fitted, sd = sigma_hat)
  } else if (identical(family$family, "poisson")) {
    lambda <- pmax(fitted, .Machine$double.eps)
    out_q <- stats::qpois(p, lambda = lambda)
    cdf <- stats::ppois(q, lambda = lambda)
  } else if (identical(family$family, "Gamma")) {
    dispersion <- .coverage_comparator_dispersion(y, fitted, family)
    if (!is.finite(dispersion) || dispersion <= 0) {
      return(empty)
    }
    mu <- pmax(fitted, .Machine$double.eps)
    shape <- 1 / dispersion
    scale <- mu * dispersion
    out_q <- stats::qgamma(p, shape = shape, scale = scale)
    cdf <- stats::pgamma(q, shape = shape, scale = scale)
  } else {
    return(empty)
  }

  list(
    q = as.numeric(out_q),
    cdf_at_q = pmin(pmax(as.numeric(cdf), 0), 1)
  )
}

#' @keywords internal
#' @noRd
.coverage_benchmark_truth_metrics <- function(dat, fitted, family, gamlss_family = NULL) {
  n <- nrow(dat)
  empty <- c(
    benchmark_mean_bias = NA_real_,
    benchmark_mean_mae = NA_real_,
    benchmark_mean_rmse = NA_real_,
    benchmark_q90_mae = NA_real_,
    benchmark_upper_tail_error_90 = NA_real_,
    benchmark_interval_coverage_95 = NA_real_,
    benchmark_pit_ks_p_value = NA_real_,
    benchmark_pit_mean_abs_error = NA_real_,
    benchmark_tail_error_lower_05 = NA_real_,
    benchmark_tail_error_upper_05 = NA_real_
  )
  if (length(fitted) != n) {
    return(empty)
  }

  out <- empty
  if ("true_mu" %in% names(dat)) {
    ok_truth <- is.finite(dat$true_mu) & is.finite(fitted)
    if (any(ok_truth)) {
      err_mu <- fitted[ok_truth] - dat$true_mu[ok_truth]
      out["benchmark_mean_bias"] <- mean(err_mu)
      out["benchmark_mean_mae"] <- mean(abs(err_mu))
      out["benchmark_mean_rmse"] <- sqrt(mean(err_mu^2))
    }
  }

  truth_dist <- .coverage_true_margin_distribution(dat, gamlss_family, p = 0.9)
  comparator_dist <- .coverage_comparator_distribution(dat$response, fitted, family, truth_dist$q, p = 0.9)
  ok_q90 <- is.finite(truth_dist$q) & is.finite(comparator_dist$q)
  if (any(ok_q90)) {
    out["benchmark_q90_mae"] <- mean(abs(comparator_dist$q[ok_q90] - truth_dist$q[ok_q90]), na.rm = TRUE)
  }
  ok_tail <- is.finite(truth_dist$cdf) & is.finite(comparator_dist$cdf_at_q)
  if (any(ok_tail)) {
    true_upper_tail <- 1 - truth_dist$cdf[ok_tail]
    comparator_upper_tail <- 1 - comparator_dist$cdf_at_q[ok_tail]
    out["benchmark_upper_tail_error_90"] <- mean(comparator_upper_tail - true_upper_tail, na.rm = TRUE)
  }

  ok_obs <- is.finite(dat$response) & is.finite(fitted)
  if (!identical(family$family, "gaussian") || sum(ok_obs) < 3L) {
    return(out)
  }
  sigma_hat <- stats::sd(dat$response[ok_obs] - fitted[ok_obs])
  if (!is.finite(sigma_hat) || sigma_hat <= 0) {
    return(out)
  }
  lower <- fitted[ok_obs] + stats::qnorm(0.025) * sigma_hat
  upper <- fitted[ok_obs] + stats::qnorm(0.975) * sigma_hat
  out["benchmark_interval_coverage_95"] <- mean(dat$response[ok_obs] >= lower & dat$response[ok_obs] <= upper)
  pit <- stats::pnorm(dat$response[ok_obs], mean = fitted[ok_obs], sd = sigma_hat)
  pit <- pmin(pmax(pit, 0), 1)
  out["benchmark_pit_ks_p_value"] <- tryCatch(
    suppressWarnings(stats::ks.test(pit, "punif")$p.value),
    error = function(e) NA_real_
  )
  out["benchmark_pit_mean_abs_error"] <- abs(mean(pit, na.rm = TRUE) - 0.5)
  out["benchmark_tail_error_lower_05"] <- mean(pit <= 0.05, na.rm = TRUE) - 0.05
  out["benchmark_tail_error_upper_05"] <- mean(pit >= 0.95, na.rm = TRUE) - 0.05
  out
}

#' @keywords internal
#' @noRd
.coverage_benchmark_gamlss_metrics <- function(dat, fit, family) {
  n <- nrow(dat)
  empty <- c(
    benchmark_mae = NA_real_,
    benchmark_rmse = NA_real_,
    benchmark_mean_bias = NA_real_,
    benchmark_mean_mae = NA_real_,
    benchmark_mean_rmse = NA_real_,
    benchmark_q90_mae = NA_real_,
    benchmark_upper_tail_error_90 = NA_real_,
    benchmark_interval_coverage_95 = NA_real_,
    benchmark_pit_ks_p_value = NA_real_,
    benchmark_pit_mean_abs_error = NA_real_,
    benchmark_tail_error_lower_05 = NA_real_,
    benchmark_tail_error_upper_05 = NA_real_
  )
  if (!inherits(fit, "gamlss.longitudinal")) {
    return(empty)
  }

  tryCatch({
  diag_data <- tryCatch(.gl_fitted_distribution(fit, newdata = NULL, require_response = FALSE), error = function(e) NULL)
  if (is.null(diag_data) || length(diag_data$keep_index) == 0L) {
    return(empty)
  }
  idx <- diag_data$keep_index
  idx <- idx[idx >= 1L & idx <= n]
  if (length(idx) == 0L) {
    return(empty)
  }

  params <- diag_data$params
  mu_hat <- if ("mu" %in% names(params)) as.numeric(params$mu) else as.numeric(params[[1L]])
  fitted <- rep(NA_real_, n)
  fitted[idx] <- mu_hat[seq_along(idx)]

  out <- empty
  ok_obs <- is.finite(dat$response) & is.finite(fitted)
  if (any(ok_obs)) {
    err_obs <- fitted[ok_obs] - dat$response[ok_obs]
    out["benchmark_mae"] <- mean(abs(err_obs))
    out["benchmark_rmse"] <- sqrt(mean(err_obs^2))
  }
  if ("true_mu" %in% names(dat)) {
    ok_truth <- is.finite(dat$true_mu) & is.finite(fitted)
    if (any(ok_truth)) {
      err_mu <- fitted[ok_truth] - dat$true_mu[ok_truth]
      out["benchmark_mean_bias"] <- mean(err_mu)
      out["benchmark_mean_mae"] <- mean(abs(err_mu))
      out["benchmark_mean_rmse"] <- sqrt(mean(err_mu^2))
    }
  }

  truth_dist <- .coverage_true_margin_distribution(dat, family, p = 0.9)
  fitted_q90 <- rep(NA_real_, n)
  fitted_cdf_at_true_q90 <- rep(NA_real_, n)
  fitted_q025 <- rep(NA_real_, n)
  fitted_q975 <- rep(NA_real_, n)
  fitted_pit <- rep(NA_real_, n)

  fitted_q90[idx] <- tryCatch(
    as.numeric(.gl_call_family_fun("q", diag_data$family, 0.9, params))[seq_along(idx)],
    error = function(e) rep(NA_real_, length(idx))
  )
  fitted_cdf_at_true_q90[idx] <- tryCatch(
    as.numeric(.gl_call_family_fun("p", diag_data$family, truth_dist$q[idx], params))[seq_along(idx)],
    error = function(e) rep(NA_real_, length(idx))
  )
  fitted_q025[idx] <- tryCatch(
    as.numeric(.gl_call_family_fun("q", diag_data$family, 0.025, params))[seq_along(idx)],
    error = function(e) rep(NA_real_, length(idx))
  )
  fitted_q975[idx] <- tryCatch(
    as.numeric(.gl_call_family_fun("q", diag_data$family, 0.975, params))[seq_along(idx)],
    error = function(e) rep(NA_real_, length(idx))
  )
  fitted_pit[idx] <- tryCatch(
    as.numeric(.gl_call_family_fun("p", diag_data$family, dat$response[idx], params))[seq_along(idx)],
    error = function(e) rep(NA_real_, length(idx))
  )

  ok_q90 <- is.finite(truth_dist$q) & is.finite(fitted_q90)
  if (any(ok_q90)) {
    out["benchmark_q90_mae"] <- mean(abs(fitted_q90[ok_q90] - truth_dist$q[ok_q90]), na.rm = TRUE)
  }
  ok_tail <- is.finite(truth_dist$cdf) & is.finite(fitted_cdf_at_true_q90)
  if (any(ok_tail)) {
    true_upper_tail <- 1 - truth_dist$cdf[ok_tail]
    fitted_upper_tail <- 1 - fitted_cdf_at_true_q90[ok_tail]
    out["benchmark_upper_tail_error_90"] <- mean(fitted_upper_tail - true_upper_tail, na.rm = TRUE)
  }
  ok_interval <- is.finite(dat$response) & is.finite(fitted_q025) & is.finite(fitted_q975)
  if (any(ok_interval)) {
    out["benchmark_interval_coverage_95"] <- mean(
      dat$response[ok_interval] >= fitted_q025[ok_interval] &
        dat$response[ok_interval] <= fitted_q975[ok_interval]
    )
  }
  ok_pit <- is.finite(fitted_pit)
  if (sum(ok_pit) >= 3L) {
    pit <- pmin(pmax(fitted_pit[ok_pit], 0), 1)
    out["benchmark_pit_ks_p_value"] <- tryCatch(
      suppressWarnings(stats::ks.test(pit, "punif")$p.value),
      error = function(e) NA_real_
    )
    out["benchmark_pit_mean_abs_error"] <- abs(mean(pit, na.rm = TRUE) - 0.5)
    out["benchmark_tail_error_lower_05"] <- mean(pit <= 0.05, na.rm = TRUE) - 0.05
    out["benchmark_tail_error_upper_05"] <- mean(pit >= 0.95, na.rm = TRUE) - 0.05
  }

  out
  }, error = function(e) empty)
}

#' @keywords internal
#' @noRd
.coverage_fit_standard_comparator <- function(dat, family, copula, design, method, max_elapsed_sec = Inf) {
  comparator_family <- .coverage_standard_family(family)
  start <- Sys.time()
  truth <- .coverage_truth_summary(dat, copula)
  true_tau <- truth["copula_tau"] %||% NA_real_
  empty_row <- function(success, failure_type, elapsed, error = NA_character_, warning = NA_character_) {
    data.frame(
      method = method,
      success = success,
      failure_type = failure_type,
      elapsed_sec = elapsed,
      max_abs_error = NA_real_,
      max_rel_error = NA_real_,
      fitted_copula_tau = NA_real_,
      true_copula_tau = unname(true_tau),
      marginal_loglik = NA_real_,
      copula_loglik = NA_real_,
      joint_loglik = NA_real_,
      marginal_fit_method = method,
      fit_attempt = method,
      fit_attempt_trace = method,
      benchmark_comparator = method,
      benchmark_class = NA_character_,
      benchmark_estimator = NA_character_,
      benchmark_mae = NA_real_,
      benchmark_rmse = NA_real_,
      benchmark_mean_bias = NA_real_,
      benchmark_mean_mae = NA_real_,
      benchmark_mean_rmse = NA_real_,
      benchmark_q90_mae = NA_real_,
      benchmark_upper_tail_error_90 = NA_real_,
      benchmark_interval_coverage_95 = NA_real_,
      benchmark_pit_ks_p_value = NA_real_,
      benchmark_pit_mean_abs_error = NA_real_,
      benchmark_tail_error_lower_05 = NA_real_,
      benchmark_tail_error_upper_05 = NA_real_,
      benchmark_error = error,
      benchmark_warning = warning,
      stringsAsFactors = FALSE
    )
  }

  if (is.null(comparator_family)) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    return(empty_row(
      success = FALSE,
      failure_type = "unsupported comparator family",
      elapsed = elapsed,
      error = paste0("No standard comparator family mapping for GAMLSS family '", family, "'.")
    ))
  }

  dat_fit <- dat
  dat_fit$subject <- factor(dat_fit$subject)
  formula <- .coverage_standard_formula(design, method)
  captured <- .coverage_capture_conditions({
    benchmark_standard_models(
      data = dat_fit,
      formula = formula,
      subject_var = "subject",
      family = comparator_family,
      comparators = method,
      correlation = "exchangeable"
    )
  })
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  if (inherits(captured$value, "error")) {
    return(empty_row(
      success = FALSE,
      failure_type = .coverage_taxonomy(FALSE, conditionMessage(captured$value), captured$warnings, elapsed, max_elapsed_sec),
      elapsed = elapsed,
      error = conditionMessage(captured$value),
      warning = if (length(captured$warnings)) paste(unique(captured$warnings), collapse = " | ") else NA_character_
    ))
  }

  row <- captured$value$results[1L, , drop = FALSE]
  success <- isTRUE(row$success) && (!is.finite(max_elapsed_sec) || elapsed <= max_elapsed_sec)
  failure_type <- if (success) {
    "ok"
  } else if (!isTRUE(row$available)) {
    "comparator unavailable"
  } else {
    .coverage_taxonomy(FALSE, row$error, row$warning, elapsed, max_elapsed_sec)
  }
  fit <- captured$value$fits[[method]]
  fitted <- if (success && !is.null(fit)) {
    .benchmark_predict_response(fit, dat_fit)
  } else {
    rep(NA_real_, nrow(dat_fit))
  }
  truth_metrics <- .coverage_benchmark_truth_metrics(dat_fit, fitted, comparator_family, gamlss_family = family)
  data.frame(
    method = method,
    success = success,
    failure_type = failure_type,
    elapsed_sec = elapsed,
    max_abs_error = row$mae,
    max_rel_error = NA_real_,
    fitted_copula_tau = NA_real_,
    true_copula_tau = unname(true_tau),
    marginal_loglik = row$logLik,
    copula_loglik = NA_real_,
    joint_loglik = row$logLik,
    marginal_fit_method = row$estimator,
    fit_attempt = method,
    fit_attempt_trace = row$estimator,
    benchmark_comparator = row$comparator,
    benchmark_class = row$comparator_class,
    benchmark_estimator = row$estimator,
    benchmark_mae = row$mae,
    benchmark_rmse = row$rmse,
    benchmark_mean_bias = unname(truth_metrics[["benchmark_mean_bias"]]),
    benchmark_mean_mae = unname(truth_metrics[["benchmark_mean_mae"]]),
    benchmark_mean_rmse = unname(truth_metrics[["benchmark_mean_rmse"]]),
    benchmark_q90_mae = unname(truth_metrics[["benchmark_q90_mae"]]),
    benchmark_upper_tail_error_90 = unname(truth_metrics[["benchmark_upper_tail_error_90"]]),
    benchmark_interval_coverage_95 = unname(truth_metrics[["benchmark_interval_coverage_95"]]),
    benchmark_pit_ks_p_value = unname(truth_metrics[["benchmark_pit_ks_p_value"]]),
    benchmark_pit_mean_abs_error = unname(truth_metrics[["benchmark_pit_mean_abs_error"]]),
    benchmark_tail_error_lower_05 = unname(truth_metrics[["benchmark_tail_error_lower_05"]]),
    benchmark_tail_error_upper_05 = unname(truth_metrics[["benchmark_tail_error_upper_05"]]),
    benchmark_error = row$error,
    benchmark_warning = row$warning,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.coverage_longitudinal_start_fit <- function(
  dat,
  family,
  copula,
  design,
  max_outer_iter = 5,
  max_inner_iter = 8,
  max_elapsed_sec = 20
) {
  .coverage_attach_namespace("gamlss")
  .coverage_attach_namespace("gamlss.dist")
  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
  formulas <- .coverage_fit_formulas(design)
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
      include_dlcopdpar = FALSE,
      method = "RS",
      start_from = NA,
      warm_start_joint = FALSE,
      warm_start_joint_iter = 0L,
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
}

#' @keywords internal
#' @noRd
.coverage_longitudinal_attempt <- function(
  dat,
  family,
  copula,
  method,
  design,
  max_outer_iter = 8,
  max_inner_iter = 8,
  max_elapsed_sec = 20,
  start_from = NA,
  warm_start_joint = TRUE,
  start_step_size = 0.5,
  cg_max_delta = 0.5,
  attempt_label = "default"
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
        warm_start_joint = isTRUE(warm_start_joint) && isTRUE(include_dlcopdpar) && all(is.na(start_from)),
        warm_start_joint_iter = min(5L, max_outer_iter),
        start_step_size = start_step_size,
        cg_max_delta = cg_max_delta,
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
  list(
    value = captured$value,
    warnings = captured$warnings,
    elapsed = elapsed,
    attempt_label = attempt_label
  )
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
  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
  attempts <- list()
  start_from_supplied <- !all(is.na(start_from))

  if (identical(method, "cg") && !start_from_supplied) {
    rs_start <- tryCatch(
      .coverage_longitudinal_start_fit(
        dat,
        family = family,
        copula = copula,
        design = design,
        max_outer_iter = min(5L, max_outer_iter),
        max_inner_iter = max_inner_iter,
        max_elapsed_sec = max_elapsed_sec
      )$par,
      error = function(e) NA
    )
    attempts[[length(attempts) + 1L]] <- .coverage_longitudinal_attempt(
      dat, family, copula, method, design,
      max_outer_iter = max_outer_iter,
      max_inner_iter = max_inner_iter,
      max_elapsed_sec = max_elapsed_sec,
      start_from = rs_start,
      warm_start_joint = FALSE,
      attempt_label = "rs_separate_start"
    )
  }

  attempts[[length(attempts) + 1L]] <- .coverage_longitudinal_attempt(
    dat, family, copula, method, design,
    max_outer_iter = max_outer_iter,
    max_inner_iter = max_inner_iter,
    max_elapsed_sec = max_elapsed_sec,
    start_from = start_from,
    warm_start_joint = TRUE,
    attempt_label = if (start_from_supplied) "explicit_start" else "default"
  )

  if (identical(method, "cg")) {
    attempts[[length(attempts) + 1L]] <- .coverage_longitudinal_attempt(
      dat, family, copula, method, design,
      max_outer_iter = max_outer_iter,
      max_inner_iter = max_inner_iter,
      max_elapsed_sec = max_elapsed_sec,
      start_from = NA,
      warm_start_joint = FALSE,
      attempt_label = "cold_start"
    )
    interior_start <- .coverage_truth_adjacent_start(dat, family, copula, design)
    attempts[[length(attempts) + 1L]] <- .coverage_longitudinal_attempt(
      dat, family, copula, method, design,
      max_outer_iter = max_outer_iter,
      max_inner_iter = max_inner_iter,
      max_elapsed_sec = max_elapsed_sec,
      start_from = interior_start,
      warm_start_joint = FALSE,
      start_step_size = 0.2,
      cg_max_delta = 0.2,
      attempt_label = "damped_interior_start"
    )
  }

  successful <- vapply(attempts, function(x) {
    fit <- x$value
    is_fit <- inherits(fit, "gamlss.longitudinal")
    if (!is_fit) return(FALSE)
    loglik <- fit$calc_lik_out_end$log_lik
    all(is.finite(as.numeric(loglik)))
  }, logical(1))
  chosen_idx <- if (any(successful)) which(successful)[1L] else length(attempts)
  chosen <- attempts[[chosen_idx]]
  elapsed <- sum(vapply(attempts[seq_len(chosen_idx)], `[[`, numeric(1), "elapsed"))
  fit <- chosen$value
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
  success <- is_fit &&
    all(is.finite(as.numeric(loglik))) &&
    !inherits(fit, "error") &&
    (!is.finite(max_elapsed_sec) || elapsed <= max_elapsed_sec)
  truth_metrics <- if (success) {
    .coverage_benchmark_gamlss_metrics(dat, fit, family)
  } else {
    .coverage_benchmark_gamlss_metrics(dat, NULL, family)
  }
  smooth_metrics <- if (success && identical(design, "smooth")) {
    .coverage_smooth_eta_recovery(fit, dat, copula)
  } else {
    c(smooth_eta_rmse = NA_real_, smooth_eta_max_abs_error = NA_real_, smooth_eta_n = 0)
  }
  out <- data.frame(
    method = method,
    success = success,
    converged = if (is_fit && !is.null(fit$convergence$converged)) isTRUE(fit$convergence$converged) else FALSE,
    failure_type = .coverage_taxonomy(success, if (inherits(fit, "error")) conditionMessage(fit) else NA_character_, chosen$warnings, elapsed, max_elapsed_sec),
    error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
    warnings = paste(unique(unlist(lapply(attempts[seq_len(chosen_idx)], `[[`, "warnings"))), collapse = " | "),
    marginal_loglik = as.numeric(loglik["marginal"]),
    copula_loglik = as.numeric(loglik["copula"]),
    joint_loglik = as.numeric(loglik["joint"]),
    elapsed_sec = elapsed,
    max_abs_param_error = abs_error,
    max_rel_param_error = rel_error,
    fitted_copula_tau = theta_tau,
    true_copula_tau = if ("copula_tau" %in% names(truth)) truth[["copula_tau"]] else NA_real_,
    marginal_fit_method = method,
    marginal_fallback_used = FALSE,
    marginal_fallback_error = NA_character_,
    fit_attempt = chosen$attempt_label,
    fit_attempt_count = chosen_idx,
    fit_attempt_trace = paste(vapply(attempts[seq_len(chosen_idx)], `[[`, character(1), "attempt_label"), collapse = " > "),
    smooth_eta_rmse = unname(smooth_metrics[["smooth_eta_rmse"]]),
    smooth_eta_max_abs_error = unname(smooth_metrics[["smooth_eta_max_abs_error"]]),
    smooth_eta_n = unname(smooth_metrics[["smooth_eta_n"]]),
    benchmark_comparator = "gamlss.longitudinal",
    benchmark_class = "gamlss_longitudinal",
    benchmark_estimator = paste0("gamlss.longitudinal::", method),
    benchmark_mae = unname(truth_metrics[["benchmark_mae"]]),
    benchmark_rmse = unname(truth_metrics[["benchmark_rmse"]]),
    benchmark_mean_bias = unname(truth_metrics[["benchmark_mean_bias"]]),
    benchmark_mean_mae = unname(truth_metrics[["benchmark_mean_mae"]]),
    benchmark_mean_rmse = unname(truth_metrics[["benchmark_mean_rmse"]]),
    benchmark_q90_mae = unname(truth_metrics[["benchmark_q90_mae"]]),
    benchmark_upper_tail_error_90 = unname(truth_metrics[["benchmark_upper_tail_error_90"]]),
    benchmark_interval_coverage_95 = unname(truth_metrics[["benchmark_interval_coverage_95"]]),
    benchmark_pit_ks_p_value = unname(truth_metrics[["benchmark_pit_ks_p_value"]]),
    benchmark_pit_mean_abs_error = unname(truth_metrics[["benchmark_pit_mean_abs_error"]]),
    benchmark_tail_error_lower_05 = unname(truth_metrics[["benchmark_tail_error_lower_05"]]),
    benchmark_tail_error_upper_05 = unname(truth_metrics[["benchmark_tail_error_upper_05"]]),
    benchmark_error = NA_character_,
    benchmark_warning = NA_character_,
    stringsAsFactors = FALSE
  )
  parameter_results <- .coverage_parameter_results(
    .coverage_longitudinal_eta_estimates(fit),
    .coverage_true_eta_coefficients(dat, family, copula, design)
  )
  theta_time <- parameter_results[
    parameter_results$parameter == "theta" &
      parameter_results$term == "time_covariate",
    ,
    drop = FALSE
  ]
  out$benchmark_theta_time_abs_error <- if (nrow(theta_time) == 1L) {
    theta_time$abs_eta_error[[1L]]
  } else {
    NA_real_
  }
  attr(out, "parameter_results") <- parameter_results
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
  method_max_outer_iter = NULL,
  method_max_inner_iter = NULL,
  method_max_elapsed_sec = NULL,
  dependence = "moderate",
  missingness = "none",
  start_mode = c("default", "truth_adjacent")
) {
  start_mode <- match.arg(start_mode)
  method_value <- as.character(case$method)
  pick_method_control <- function(default, override) {
    if (is.null(override)) return(default)
    if (is.list(override)) {
      value <- override[[method_value]]
    } else {
      value <- override[[method_value]]
    }
    if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) default else value[[1L]]
  }
  max_outer_iter <- pick_method_control(max_outer_iter, method_max_outer_iter)
  max_inner_iter <- pick_method_control(max_inner_iter, method_max_inner_iter)
  max_elapsed_sec <- pick_method_control(max_elapsed_sec, method_max_elapsed_sec)
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
  } else if (case$method %in% c("gee", "glmm", "gam", "glmmTMB")) {
    .coverage_fit_standard_comparator(
      dat,
      family = case$family,
      copula = case$copula,
      design = case$design,
      method = case$method,
      max_elapsed_sec = max_elapsed_sec
    )
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
  method_max_outer_iter = NULL,
  method_max_inner_iter = NULL,
  method_max_elapsed_sec = NULL,
  dependence = "moderate",
  missingness = "none",
  start_mode = "default"
) {
  row_bind_fill <- function(rows) {
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0L) {
      return(data.frame())
    }
    all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
    rows <- lapply(rows, function(x) {
      missing <- setdiff(all_names, names(x))
      for (nm in missing) {
        x[[nm]] <- NA
      }
      x[all_names]
    })
    do.call(rbind, rows)
  }
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
      method_max_outer_iter = method_max_outer_iter,
      method_max_inner_iter = method_max_inner_iter,
      method_max_elapsed_sec = method_max_elapsed_sec,
      dependence = dependence,
      missingness = missingness,
      start_mode = start_mode
    )
    parameter_rows[[i]] <- attr(case_result, "parameter_results")
    rows[[i]] <- case_result
  }
  results <- .coverage_add_review_metrics(row_bind_fill(rows))
  parameter_rows <- parameter_rows[!vapply(parameter_rows, is.null, logical(1))]
  attr(results, "parameter_results") <- if (length(parameter_rows) == 0L) {
    data.frame()
  } else {
    row_bind_fill(parameter_rows)
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

#' Run opt-in distribution/copula/method coverage simulations
#'
#' @param families Optional family names. Defaults to all supported non-mixed
#'   `gamlss.dist` families with `q`, `p`, and `d` functions.
#' @param copulas Copula codes.
#' @param methods Fit methods. Defaults to `"gamlss"`, `"rs_separate"`,
#'   `"rs_joint"`, and `"cg"`. Method `"gamlss2"` is also supported when the
#'   optional non-CRAN package is installed and explicitly requested. Standard
#'   comparator methods `"gee"`, `"glmm"`, `"gam"`, and `"glmmTMB"` are
#'   available for families that map to common mean-model families.
#' @param designs Simulation designs: `"intercept"`, `"covariate"`,
#'   `"scale"`, `"time_dependence"`, or `"smooth"`.
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
  methods = .coverage_default_methods(),
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
