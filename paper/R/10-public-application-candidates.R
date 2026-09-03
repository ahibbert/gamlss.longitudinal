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
      copula_candidates = "C",
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
      copula = "N",
      missingness = "error",
      margin_candidates = c("GG", "LOGNO", "GA", "NO"),
      copula_candidates = "N",
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
      copula = "C",
      missingness = "segment",
      margin_candidates = c("PO", "NBI", "ZINBI"),
      copula_candidates = "C",
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

jss_screen_public_application_margins <- function(
    name = jss_public_application_names(),
    data = NULL,
    n_cycles = 100L) {
  spec <- jss_public_application_spec(name)
  if (is.null(data)) data <- jss_prepare_public_application(name)
  needed <- unique(all.vars(spec$formulas$mu))
  data <- data[is.finite(data[[spec$response]]), needed, drop = FALSE]
  data <- data[stats::complete.cases(data), , drop = FALSE]

  rows <- lapply(spec$margin_candidates, function(code) {
    family <- jss_margin_family(code)
    arguments <- list(
      formula = spec$formulas$mu,
      family = family,
      data = data,
      control = gamlss::gamlss.control(n.cyc = as.integer(n_cycles), trace = FALSE)
    )
    if ("sigma" %in% family$parameters) arguments$sigma.formula <- ~ 1
    if ("nu" %in% family$parameters) arguments$nu.formula <- ~ 1
    if ("tau" %in% family$parameters) arguments$tau.formula <- ~ 1

    warnings <- character()
    started <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      withCallingHandlers(
        do.call(gamlss::gamlss, arguments),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - started

    if (inherits(fit, "error")) {
      return(data.frame(
        family = code, converged = FALSE, logLik = NA_real_, df = NA_real_,
        AIC = NA_real_, BIC = NA_real_, elapsed_sec = elapsed,
        warning = "", error = conditionMessage(fit)
      ))
    }
    likelihood <- logLik(fit)
    data.frame(
      family = code,
      converged = isTRUE(fit$converged),
      logLik = as.numeric(likelihood),
      df = attr(likelihood, "df"),
      AIC = AIC(fit),
      BIC = gamlss::GAIC(fit, k = log(nrow(data))),
      elapsed_sec = elapsed,
      warning = paste(unique(warnings), collapse = " | "),
      error = ""
    )
  })
  result <- do.call(rbind, rows)
  result[order(result$AIC), ]
}

jss_fit_public_application <- function(
    name = jss_public_application_names(),
    data = NULL,
    formulas = NULL,
    compute_vcov = FALSE,
    max_outer_iter = 40L,
    max_inner_iter = 50L,
    max_elapsed_sec = 600,
    start_from = NA,
    verbose = 0) {
  spec <- jss_public_application_spec(name)
  if (is.null(data)) data <- jss_prepare_public_application(name)
  if (is.null(formulas)) formulas <- spec$formulas
  user_start <- !all(is.na(start_from))

  control <- gamlss.longitudinal::gamlss_longitudinal_control(
    outer_tol = 1e-3,
    max_outer_iter = as.integer(max_outer_iter),
    max_elapsed_sec = max_elapsed_sec,
    rs = list(
      inner_tol = 1e-3,
      max_inner_iter = as.integer(max_inner_iter),
      warm_start_joint = !user_start,
      warm_start_joint_iter = if (user_start) 0L else 5L
    )
  )

  fit <- gamlss.longitudinal::gamlss_longitudinal(
    dataset = data,
    margin_dist = jss_margin_family(spec$margin),
    copula_dist = spec$copula,
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
    compute_vcov = compute_vcov,
    verbose = verbose
  )
  fit$application_candidate <- spec$name
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
  jss_fit_public_application(name, formulas = formulas, ...)
}

jss_public_application_outputs <- function(fit) {
  coefficient_output <- if (isTRUE(fit$vcov$precomputed)) {
    gamlss.longitudinal::publication_table(fit)
  } else {
    data.frame(
      term = names(stats::coef(fit)),
      estimate = unname(stats::coef(fit)),
      inference = "Not computed in candidate fit",
      stringsAsFactors = FALSE
    )
  }
  list(
    model_specification = gamlss.longitudinal::model_spec(fit),
    model_check = gamlss.longitudinal::check_model(fit),
    coefficient_table = coefficient_output,
    margin_diagnostics = gamlss.longitudinal::plot_dist(fit),
    copula_diagnostics = gamlss.longitudinal::plot_copula_diagnostics(fit),
    term_plot = gamlss.longitudinal::plot_terms(fit)
  )
}
