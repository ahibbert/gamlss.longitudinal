# Public application candidates for the worked application. This module is
# deliberately outside the default targets graph
# until one candidate has been selected for the paper and vignette.

jss_public_application_names <- function() {
  c("patents", "pbc", "steps")
}

jss_load_package_dataset <- function(object) {
  environment <- new.env(parent = emptyenv())
  loaded <- try(
    utils::data(list = object, package = "gamlss.longitudinal", envir = environment),
    silent = TRUE
  )
  if (!inherits(loaded, "try-error") && exists(object, envir = environment, inherits = FALSE)) {
    return(get(object, envir = environment, inherits = FALSE))
  }

  source_file <- file.path("data", paste0(object, ".rda"))
  if (!file.exists(source_file)) {
    stop("Cannot locate package dataset '", object, "'.", call. = FALSE)
  }
  load(source_file, envir = environment)
  get(object, envir = environment, inherits = FALSE)
}

jss_public_application_spec <- function(name = jss_public_application_names()) {
  name <- match.arg(name, jss_public_application_names())

  switch(name,
    patents = list(
      name = "patents",
      title = "US patents and R&D",
      data_object = "patents_panel",
      response = "patents",
      subject = "firm",
      time = "year",
      margin = "NBI",
      copula = "C",
      missingness = "error",
      margin_candidates = c("PO", "NBI", "DEL"),
      copula_candidates = c("N", "C", "F", "G", "J", "t"),
      selection_cycles = 200L,
      margin_selection_note = paste(
        "DEL ranks first in the lightweight margin screen, but its exact",
        "joint rectangle fit does not complete within the routine application",
        "budget. NBI is retained from the selected shortlist for full fitting."
      ),
      copula_selection_note = paste(
        "Gumbel ranks first in the NBI copula screen, but the full exact",
        "rectangle fit exceeds the routine application budget. Clayton is",
        "retained from the selected shortlist for full fitting."
      ),
      fit_control = list(
        outer_tol = 0.01, inner_tol = 1e-5, max_outer_iter = 60L,
        max_inner_iter = 150L, max_elapsed_sec = 1200
      ),
      formulas = list(
        mu = patents ~ year_factor + log_rd + scientific + capital_z,
        sigma = ~ scientific + capital_z,
        nu = ~ 1,
        tau = ~ 1,
        theta = ~ capital_z,
        zeta = ~ 1
      ),
      comparisons = list(theta_null = ~ 1, sigma_null = ~ 1)
    ),
    pbc = list(
      name = "pbc",
      title = "Mayo PBC prothrombin",
      data_object = "pbc_prothrombin",
      response = "prothrombin",
      subject = "subject",
      time = "visit",
      margin = "GG",
      copula = "F",
      missingness = "error",
      margin_candidates = c("GG", "LOGNO", "GA", "NO"),
      copula_candidates = c("N", "C", "F", "G", "J", "t"),
      selection_cycles = 500L,
      fit_control = list(
        outer_tol = 1e-4, inner_tol = 1e-6, max_outer_iter = 80L,
        max_inner_iter = 200L, max_elapsed_sec = 900
      ),
      formulas = list(
        mu = prothrombin ~ splines::ns(years, df = 4) + baseline_stage +
          baseline_age_z + sex + treatment + baseline_log_bilirubin_z,
        sigma = ~ baseline_log_bilirubin_z,
        nu = ~ 1,
        tau = ~ 1,
        theta = ~ baseline_stage_z,
        zeta = ~ 1
      ),
      comparisons = list(theta_null = ~ 1, sigma_null = ~ 1)
    ),
    steps = list(
      name = "steps",
      title = "Vietnamese adolescent daily steps",
      data_object = "vietnam_steps",
      response = "steps",
      subject = "subject",
      time = "day",
      margin = "NBI",
      copula = "G",
      missingness = "segment",
      margin_candidates = c("PO", "NBI", "ZINBI"),
      copula_candidates = c("N", "C", "F", "G", "J", "t"),
      selection_cycles = 200L,
      fit_control = list(
        outer_tol = 1e-3, inner_tol = 1e-4, max_outer_iter = 50L,
        max_inner_iter = 100L, max_elapsed_sec = 1200
      ),
      formulas = list(
        mu = steps ~ day_name + bmi_z + sex + paqc_z,
        sigma = ~ 1,
        nu = ~ 1,
        tau = ~ 1,
        theta = ~ sex,
        zeta = ~ 1
      ),
      comparisons = list(theta_null = ~ 1, sigma_candidate = ~ sex)
    )
  )
}

jss_prepare_public_application <- function(name = jss_public_application_names(), data = NULL) {
  spec <- jss_public_application_spec(name)
  if (is.null(data)) data <- jss_load_package_dataset(spec$data_object)
  data <- as.data.frame(data)

  if (identical(spec$name, "patents")) {
    data$year_factor <- factor(data$year)
    data$log_rd <- log(data$rd)
    data$capital_z <- as.numeric(scale(log(data$capital_1972)))
    if (any(!is.finite(data$log_rd)) || any(!is.finite(data$capital_z))) {
      stop("Patents predictors must be finite after log transformation.", call. = FALSE)
    }
  } else if (identical(spec$name, "pbc")) {
    visits <- table(data$subject)
    retained <- as.integer(names(visits[visits >= 4L]))
    data <- data[data$subject %in% retained, ]
    data <- data[order(data$subject, data$visit), ]
    data$baseline_age_z <- as.numeric(scale(data$baseline_age))
    data$baseline_log_bilirubin_z <- as.numeric(scale(log(data$baseline_bilirubin)))
    stage_number <- as.numeric(data$baseline_stage)
    data$baseline_stage_z <- as.numeric(scale(stage_number))
  } else {
    data <- data[order(data$subject, data$day), ]
    data$bmi_z <- as.numeric(scale(data$bmi))
    data$paqc_z <- as.numeric(scale(data$paqc))
  }

  if (anyDuplicated(data[c(spec$subject, spec$time)])) {
    stop("Candidate data must contain one row per subject/time combination.", call. = FALSE)
  }
  data
}

jss_public_application_contract <- function(name = jss_public_application_names()) {
  spec <- jss_public_application_spec(name)
  data <- jss_prepare_public_application(name)
  response <- data[[spec$response]]
  observed <- is.finite(response)
  visits <- table(data[[spec$subject]])

  list(
    specification = spec,
    subjects = length(visits),
    observations = nrow(data),
    observed_responses = sum(observed),
    missing_responses = sum(!observed),
    visits = stats::quantile(as.numeric(visits), c(0, 0.25, 0.5, 0.75, 1)),
    response_range = range(response[observed]),
    integer_response = all(response[observed] == floor(response[observed]))
  )
}

jss_margin_family <- function(code) {
  getExportedValue("gamlss.dist", code)()
}

jss_select_public_application_distributions <- function(
    name = jss_public_application_names(),
    data = NULL,
    margin_families = NULL,
    copula_families = NULL,
    criterion = "AIC") {
  spec <- jss_public_application_spec(name)
  if (is.null(data)) data <- jss_prepare_public_application(name)
  if (is.null(margin_families)) margin_families <- spec$margin_candidates
  if (is.null(copula_families)) copula_families <- spec$copula_candidates

  margin_selection <- gamlss.longitudinal::select_margin(
    data,
    response_var = spec$response,
    time_var = spec$time,
    time_intercepts = TRUE,
    time_n_cycles = spec$selection_cycles,
    families = margin_families,
    trace = FALSE
  )
  if (!nrow(margin_selection) || !isTRUE(margin_selection$converged[[1L]])) {
    stop(
      "The selected marginal screen did not produce a converged best fit for ",
      spec$name,
      ".",
      call. = FALSE
    )
  }

  screen_selected_margin <- gamlss.longitudinal::best_fit(margin_selection)$family_name
  selected_margin <- spec$margin
  if (!selected_margin %in% margin_selection$family) {
    stop(
      "The analysis margin '", selected_margin,
      "' was not retained by select_margin() for ", spec$name, ".",
      call. = FALSE
    )
  }

  copula_selection <- gamlss.longitudinal::select_copula(
    data = data,
    response_var = spec$response,
    margin_dist = jss_margin_family(selected_margin),
    subject_var = spec$subject,
    time_var = spec$time,
    families = copula_families,
    criterion = criterion,
    time_intercepts = TRUE,
    margin_n_cycles = spec$selection_cycles,
    copula_time_intercepts = TRUE
  )
  screen_selected_copula <- gamlss.longitudinal::best_fit_family(copula_selection)
  selected_copula <- spec$copula
  if (!selected_copula %in% copula_selection$family) {
    stop(
      "The analysis copula '", selected_copula,
      "' was not retained by select_copula() for ", spec$name, ".",
      call. = FALSE
    )
  }

  list(
    margin = margin_selection,
    copula = copula_selection,
    screen_selected_margin = screen_selected_margin,
    selected_margin = selected_margin,
    screen_selected_copula = screen_selected_copula,
    selected_copula = selected_copula,
    margin_selection_note = if (is.null(spec$margin_selection_note)) {
      "Top-ranked margin retained."
    } else {
      spec$margin_selection_note
    },
    copula_selection_note = if (is.null(spec$copula_selection_note)) {
      "Top-ranked registered copula retained."
    } else {
      spec$copula_selection_note
    },
    criterion = criterion,
    contract = list(
      registry_version = unique(
        gamlss.longitudinal::longitudinal_capabilities("routes")$registry_version
      ),
      margin_time_intercepts = TRUE,
      copula_time_intercepts = TRUE,
      margin_n_cycles = spec$selection_cycles,
      copula_screen_version = attr(copula_selection, "screen_version"),
      copula_screen_likelihood = attr(copula_selection, "screen_likelihood")
    )
  )
}

jss_screen_public_application_margins <- function(
    name = jss_public_application_names(),
    data = NULL,
    n_cycles = 100L) {
  spec <- jss_public_application_spec(name)
  if (is.null(data)) data <- jss_prepare_public_application(name)
  gamlss.longitudinal::select_margin(
    data,
    response_var = spec$response,
    time_var = spec$time,
    time_intercepts = TRUE,
    time_n_cycles = max(as.integer(n_cycles), spec$selection_cycles),
    families = spec$margin_candidates,
    trace = FALSE
  )
}

jss_patents_refined_margin_screen <- function(data = NULL, n_cycles = 200L) {
  if (is.null(data)) data <- jss_prepare_public_application("patents")
  specifications <- list(
    NBI_base = list(family = "NBI", sigma = ~ scientific + capital_z),
    NBI_rd_dispersion = list(
      family = "NBI", sigma = ~ scientific + capital_z + log_rd
    ),
    DEL_base = list(
      family = "DEL", sigma = ~ scientific + capital_z, nu = ~ 1
    ),
    DEL_rd_dispersion = list(
      family = "DEL", sigma = ~ scientific + capital_z + log_rd, nu = ~ 1
    ),
    DEL_shape_covariates = list(
      family = "DEL", sigma = ~ scientific + capital_z,
      nu = ~ scientific + capital_z
    )
  )
  mu_formula <- jss_public_application_spec("patents")$formulas$mu
  needed <- unique(c(
    all.vars(mu_formula),
    unlist(lapply(specifications, function(x) {
      unlist(lapply(x[intersect(names(x), c("sigma", "nu"))], all.vars))
    }))
  ))
  data <- data[stats::complete.cases(data[needed]), needed, drop = FALSE]

  rows <- lapply(names(specifications), function(model) {
    specification <- specifications[[model]]
    warnings <- character()
    started <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      withCallingHandlers(
        gamlss::gamlss(
          formula = mu_formula,
          sigma.formula = specification$sigma,
          nu.formula = if (is.null(specification$nu)) ~ 1 else specification$nu,
          family = jss_margin_family(specification$family),
          data = data,
          control = gamlss::gamlss.control(
            n.cyc = as.integer(n_cycles), trace = FALSE
          )
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = identity
    )
    elapsed <- proc.time()[["elapsed"]] - started
    if (inherits(fit, "error")) {
      return(data.frame(
        model = model, family = specification$family, converged = FALSE,
        df = NA_real_, logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
        elapsed_sec = elapsed, warning = "", error = conditionMessage(fit)
      ))
    }
    likelihood <- stats::logLik(fit)
    data.frame(
      model = model,
      family = specification$family,
      converged = isTRUE(fit$converged),
      df = attr(likelihood, "df"),
      logLik = as.numeric(likelihood),
      AIC = stats::AIC(fit),
      BIC = gamlss::GAIC(fit, k = log(nrow(data))),
      elapsed_sec = elapsed,
      warning = paste(unique(warnings), collapse = " | "),
      error = ""
    )
  })
  result <- do.call(rbind, rows)
  result[order(result$AIC), , drop = FALSE]
}

jss_patents_joint_refinement_specs <- function() {
  base <- jss_public_application_spec("patents")$formulas
  list(
    NBI_base = list(margin = "NBI", formulas = base),
    NBI_rd_dispersion = list(
      margin = "NBI",
      formulas = utils::modifyList(base, list(
        sigma = ~ scientific + capital_z + log_rd
      ))
    ),
    DEL_base = list(margin = "DEL", formulas = base),
    DEL_rd_dispersion = list(
      margin = "DEL",
      formulas = utils::modifyList(base, list(
        sigma = ~ scientific + capital_z + log_rd
      ))
    ),
    DEL_shape_covariates = list(
      margin = "DEL",
      formulas = utils::modifyList(base, list(
        nu = ~ scientific + capital_z
      ))
    )
  )
}

jss_patents_rd_dispersion_start <- function(fit) {
  coefficients <- stats::coef(fit)
  insertion <- match("sigma.capital_z", names(coefficients))
  if (is.na(insertion)) {
    stop("The patents base fit does not contain sigma.capital_z.", call. = FALSE)
  }
  start <- append(coefficients, values = 0, after = insertion)
  names(start)[insertion + 1L] <- "sigma.log_rd"
  start
}

jss_patents_nbi_rd_start <- jss_patents_rd_dispersion_start

jss_fit_public_application <- function(
    name = jss_public_application_names(),
    data = NULL,
    formulas = NULL,
    margin = NULL,
    copula = NULL,
    compute_vcov = FALSE,
    outer_tol = NULL,
    inner_tol = NULL,
    max_outer_iter = NULL,
    max_inner_iter = NULL,
    max_elapsed_sec = NULL,
    start_from = NA,
    verbose = 0) {
  spec <- jss_public_application_spec(name)
  if (is.null(data)) data <- jss_prepare_public_application(name)
  if (is.null(formulas)) formulas <- spec$formulas
  if (is.null(margin)) margin <- spec$margin
  if (is.null(copula)) copula <- spec$copula
  if (is.null(outer_tol)) outer_tol <- spec$fit_control$outer_tol
  if (is.null(max_outer_iter)) max_outer_iter <- spec$fit_control$max_outer_iter
  if (is.null(max_inner_iter)) max_inner_iter <- spec$fit_control$max_inner_iter
  if (is.null(max_elapsed_sec)) max_elapsed_sec <- spec$fit_control$max_elapsed_sec
  user_start <- !all(is.na(start_from))
  if (is.null(inner_tol)) inner_tol <- spec$fit_control$inner_tol
  use_warm_start <- !user_start && !identical(copula, "t")

  control <- gamlss.longitudinal::gamlss_longitudinal_control(
    outer_tol = outer_tol,
    max_outer_iter = as.integer(max_outer_iter),
    max_elapsed_sec = max_elapsed_sec,
    rs = list(
      inner_tol = inner_tol,
      max_inner_iter = as.integer(max_inner_iter),
      warm_start_joint = use_warm_start,
      warm_start_joint_iter = if (use_warm_start) 5L else 0L
    )
  )
  segmented <- identical(spec$missingness, "segment")

  fit <- gamlss.longitudinal::gamlss_longitudinal(
    dataset = data,
    margin_dist = jss_margin_family(margin),
    copula_dist = copula,
    subject_var = spec$subject,
    time_var = spec$time,
    missingness = spec$missingness,
    mu.formula = formulas$mu,
    sigma.formula = formulas$sigma,
    nu.formula = formulas$nu,
    tau.formula = formulas$tau,
    theta.formula = formulas$theta,
    zeta.formula = formulas$zeta,
    method = "RS",
    start_from = start_from,
    optimizer_control = control,
    compute_vcov = compute_vcov && !segmented,
    vcov_method = "analytical",
    verbose = verbose
  )
  if (compute_vcov && segmented && isTRUE(fit$convergence$converged)) {
    fit$vcov <- stats::vcov(
      fit,
      method = "sandwich",
      details = TRUE,
      sandwich_bread_method = "numderiv",
      progress = FALSE
    )
    fit$vcov_meta$precomputed <- TRUE
    fit$vcov_meta$numderiv <- FALSE
    fit$vcov_meta$method <- "sandwich"
    fit$vcov_meta$method_used <- fit$vcov$method
  }
  fit$application_candidate <- spec$name
  fit$application_selection <- list(
    margin = margin,
    copula = copula
  )
  fit
}

jss_fit_public_application_comparator <- function(
    fit,
    component = c("theta_null", "sigma_null", "sigma_candidate"),
    ...) {
  component <- match.arg(component)
  name <- fit$application_candidate
  if (is.null(name) || !name %in% jss_public_application_names()) {
    stop("Supply comparators through a named candidate fit result.", call. = FALSE)
  }
  spec <- jss_public_application_spec(name)
  formula <- spec$comparisons[[component]]
  if (is.null(formula)) stop("Comparator is not defined for this candidate.", call. = FALSE)
  formulas <- spec$formulas
  if (grepl("^theta", component)) formulas$theta <- formula else formulas$sigma <- formula
  jss_fit_public_application(
    name,
    formulas = formulas,
    margin = fit$margin_dist$family[[1L]],
    copula = fit$copula_dist,
    ...
  )
}

jss_public_application_outputs <- function(fit) {
  capture_output <- function(code) {
    tryCatch(
      code(),
      error = function(e) structure(
        list(error = conditionMessage(e)),
        class = "candidate_output_error"
      )
    )
  }
  vcov_method <- if (identical(fit$missingness$requested, "segment")) {
    "sandwich"
  } else {
    "analytical"
  }
  diagnostic_data <- fit$dataset
  time_var <- fit$var_map[[fit$time_var]]
  response_var <- fit$var_map[[fit$response_var]]
  if (is.null(time_var)) time_var <- fit$time_var
  if (is.null(response_var)) response_var <- fit$response_var
  if (!is.null(diagnostic_data) &&
      time_var %in% names(diagnostic_data) &&
      response_var %in% names(diagnostic_data)) {
    observed_by_time <- tapply(
      is.finite(diagnostic_data[[response_var]]),
      diagnostic_data[[time_var]],
      sum
    )
    diagnostic_times <- names(sort(observed_by_time, decreasing = TRUE))[seq_len(
      min(4L, length(observed_by_time))
    )]
    diagnostic_data <- diagnostic_data[
      as.character(diagnostic_data[[time_var]]) %in% diagnostic_times,
      ,
      drop = FALSE
    ]
  }
  vcov_output <- if (
      isTRUE(fit$vcov_meta$precomputed) &&
      identical(fit$vcov_meta$method, vcov_method) &&
      !is.null(fit$vcov)) {
    fit$vcov
  } else {
    capture_output(function() {
      vcov_args <- list(
        object = fit,
        method = vcov_method,
        details = TRUE,
        progress = FALSE
      )
      if (identical(vcov_method, "sandwich")) {
        vcov_args$sandwich_bread_method <- "numderiv"
      }
      do.call(stats::vcov, vcov_args)
    })
  }
  fit_for_outputs <- fit
  if (!inherits(vcov_output, "candidate_output_error")) {
    fit_for_outputs$vcov <- vcov_output
    fit_for_outputs$vcov_meta$precomputed <- TRUE
    fit_for_outputs$vcov_meta$numderiv <- FALSE
    fit_for_outputs$vcov_meta$method <- vcov_method
  }

  list(
    variance_covariance = vcov_output,
    model_specification = capture_output(function() {
      gamlss.longitudinal::model_spec(fit_for_outputs)
    }),
    model_check = capture_output(function() {
      gamlss.longitudinal::check_model(fit_for_outputs)
    }),
    coefficient_table = capture_output(function() {
      gamlss.longitudinal::publication_table(
        fit_for_outputs,
        vcov_method = vcov_method
      )
    }),
    margin_diagnostics = capture_output(function() {
      gamlss.longitudinal::plot_dist(
        dataset = diagnostic_data,
        fit = fit_for_outputs
      )
    }),
    copula_diagnostics = capture_output(function() {
      gamlss.longitudinal::plot_copula_diagnostics(fit_for_outputs)
    }),
    term_plot = capture_output(function() {
      if (identical(vcov_method, "sandwich") &&
          !inherits(vcov_output, "candidate_output_error")) {
        gamlss.longitudinal::plot_fixed_terms(
          fit_for_outputs,
          vcov_obj = vcov_output
        )
      } else {
        gamlss.longitudinal::plot_terms(fit_for_outputs)
      }
    })
  )
}
