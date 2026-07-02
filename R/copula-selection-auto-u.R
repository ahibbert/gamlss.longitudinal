.select_copula_auto_u <- function(
    data,
    response_var,
    margin_dist,
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    time_var,
    time_intercepts) {
  if (is.null(response_var) || !is.character(response_var) || length(response_var) != 1L) {
    stop(
      "Provide u_var, u, u1/u2, object, or response_var so select_copula() can create pseudo-observations.",
      call. = FALSE
    )
  }
  if (!response_var %in% names(data)) {
    stop("response_var='", response_var, "' not found in 'data'.", call. = FALSE)
  }

  margin_selection <- NULL
  if (inherits(margin_dist, "margin_selection")) {
    margin_selection <- margin_dist
    if (isTRUE(attr(margin_selection, "time_intercepts"))) {
      time_intercepts <- TRUE
      selected_time_var <- attr(margin_selection, "time_var")
      if (!is.null(selected_time_var)) {
        time_var <- selected_time_var
      }
    }
    margin_dist <- best_fit_family(margin_selection)
    if (is.null(margin_dist)) {
      stop("The supplied select_margin() result did not contain a usable marginal family.", call. = FALSE)
    }
  } else if (is.null(margin_dist)) {
    warning(
      "'margin_dist' was not supplied; select_copula() is selecting a temporary marginal distribution with select_margin().",
      call. = FALSE
    )
    margin_selection <- suppressWarnings(suppressMessages(
      select_margin(
        data,
        response_var = response_var,
        time_var = if (isTRUE(time_intercepts)) time_var else NULL,
        time_intercepts = time_intercepts,
        trace = FALSE
      )
    ))
    margin_dist <- best_fit_family(margin_selection)
    if (is.null(margin_dist)) {
      stop("select_margin() did not return a usable marginal family.", call. = FALSE)
    }
  }

  if (isTRUE(time_intercepts)) {
    if (is.null(time_var) || !is.character(time_var) || length(time_var) != 1L) {
      stop("'time_var' must be a single column name when 'time_intercepts = TRUE'.", call. = FALSE)
    }
    if (!time_var %in% names(data)) {
      stop("time_var='", time_var, "' not found in 'data'.", call. = FALSE)
    }
    mu.formula <- .select_copula_time_response_formula(mu.formula, response_var, time_var)
    sigma.formula <- .select_copula_time_rhs_formula(sigma.formula, time_var)
    nu.formula <- .select_copula_time_rhs_formula(nu.formula, time_var)
    tau.formula <- .select_copula_time_rhs_formula(tau.formula, time_var)
  }

  pfun <- tryCatch(
    get(paste0("p", margin_dist$family[1]), envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE),
    error = function(e) NULL
  )
  if (is.null(pfun)) {
    stop("'margin_dist' must provide a gamlss.dist CDF function.", call. = FALSE)
  }

  mu.formula <- .select_copula_response_formula(mu.formula, response_var)
  args <- list(
    formula = mu.formula,
    family = margin_dist,
    data = data,
    trace = FALSE,
    control = gamlss::gamlss.control(n.cyc = 100, trace = FALSE)
  )
  parameter_names <- names(margin_dist$parameters)
  if ("sigma" %in% parameter_names) args$sigma.fo <- .select_copula_rhs_formula(sigma.formula)
  if ("nu" %in% parameter_names) args$nu.fo <- .select_copula_rhs_formula(nu.formula)
  if ("tau" %in% parameter_names) args$tau.fo <- .select_copula_rhs_formula(tau.formula)

  formula_vars <- unique(unlist(
    lapply(
      c(list(args$formula), args[c("sigma.fo", "nu.fo", "tau.fo")]),
      function(formula) if (is.null(formula)) character(0) else all.vars(stats::as.formula(formula))
    ),
    use.names = FALSE
  ))
  missing_vars <- setdiff(formula_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      "Temporary marginal formula variable(s) not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  keep <- is.finite(data[[response_var]]) & stats::complete.cases(data[, formula_vars, drop = FALSE])
  if (sum(keep) < 3L) {
    stop("Need at least three complete response rows to create pseudo-observations.", call. = FALSE)
  }
  data_fit <- data[keep, , drop = FALSE]
  args$data <- data_fit[, formula_vars, drop = FALSE]

  margin_fit <- tryCatch(
    do.call(gamlss::gamlss, args),
    error = function(e) {
      stop("Temporary marginal model for copula selection failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (identical(margin_fit$converged, FALSE)) {
    warning(
      "Temporary marginal model for copula selection did not report convergence; ",
      "pseudo-observations will still be used for screening. Review the selected fit before final modelling.",
      call. = FALSE
    )
  }

  cdf_args <- list(q = data_fit[[response_var]], y = data_fit[[response_var]], x = data_fit[[response_var]])
  for (parameter in parameter_names) {
    cdf_args[[parameter]] <- stats::fitted(margin_fit, what = parameter)
  }
  fixed_unlinked_values <- attr(margin_dist, "fixed_unlinked_values")
  if (length(fixed_unlinked_values) > 0L) {
    for (parameter in names(fixed_unlinked_values)) {
      value <- fixed_unlinked_values[[parameter]]
      if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
        cdf_args[[parameter]] <- rep(value, nrow(data_fit))
      }
    }
  }
  cdf_args <- cdf_args[names(cdf_args) %in% formalArgs(pfun)]
  u <- tryCatch(
    do.call(pfun, cdf_args),
    error = function(e) {
      stop("Could not compute pseudo-observations from the temporary marginal model: ", conditionMessage(e), call. = FALSE)
    }
  )
  u_full <- rep(NA_real_, nrow(data))
  u_full[keep] <- .copula_clamp01(u)
  data[[".u_auto"]] <- u_full
  list(
    data = data,
    u_var = ".u_auto",
    margin_selection = margin_selection,
    source = if (is.null(margin_selection)) "margin_dist" else "select_margin"
  )
}
