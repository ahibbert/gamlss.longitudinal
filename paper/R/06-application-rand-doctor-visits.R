jss_run_06_rand_doctor_visits <- function(settings) {
  paths <- jss_rand_paths(settings)
  data <- jss_rand_load_data(settings)
  prepared <- jss_rand_prepare_data(data)
  sample_data <- jss_rand_sample_data(prepared$data, settings)
  fit_budget <- jss_rand_fit_budget(settings)

  audit <- jss_rand_data_audit(prepared$data, sample_data, settings)
  utils::write.csv(audit, paths$data, row.names = FALSE)

  family_screen <- jss_rand_family_screen(sample_data)
  utils::write.csv(family_screen$table, paths$family_screen, row.names = FALSE)
  selected_family <- jss_rand_select_margin_family(family_screen$table)

  copula_screen <- jss_rand_copula_screen(sample_data, fit_budget, settings$profile, selected_family)
  utils::write.csv(copula_screen, paths$copula_screen, row.names = FALSE)
  selected_copula <- jss_rand_select_copula(copula_screen)

  theta_screen <- jss_rand_theta_screen(sample_data, fit_budget, settings$profile, selected_family, selected_copula)
  theta_screen$table <- jss_rand_add_theta_evidence(theta_screen$table)
  utils::write.csv(theta_screen$table, paths$theta_screen, row.names = FALSE)

  primary <- jss_rand_select_primary_fit(theta_screen$fits, theta_screen$table)
  contrast <- jss_rand_select_covariate_contrast(theta_screen$fits, theta_screen$table)
  final_fit <- primary$fit
  contrast_fit <- contrast$fit
  intercept_fit <- theta_screen$fits[["theta_intercept"]]

  coefficient_table <- jss_rand_coefficient_table(final_fit, primary$model)
  utils::write.csv(coefficient_table, paths$coefficients, row.names = FALSE)

  se_comparison <- jss_rand_standard_error_comparison(
    sample_data = sample_data,
    marginal_se_table = family_screen$se_tables[[selected_family]],
    intercept_fit = intercept_fit,
    final_fit = contrast_fit
  )
  utils::write.csv(se_comparison, paths$se_comparison, row.names = FALSE)

  convergence <- jss_rand_convergence_table(
    theta_fits = theta_screen$fits,
    selected_model = primary$model
  )
  utils::write.csv(convergence, paths$convergence, row.names = FALSE)

  workflow <- jss_rand_write_workflow_artifacts(
    full_data = prepared$data,
    sample_data = sample_data,
    family_screen = family_screen,
    copula_screen = copula_screen,
    theta_screen = theta_screen,
    selected_family = selected_family,
    selected_copula = selected_copula,
    final_fit = final_fit,
    contrast_fit = contrast_fit,
    intercept_fit = intercept_fit,
    paths = paths,
    settings = settings
  )

  status <- jss_rand_status_table(
    settings = settings,
    paths = paths,
    audit = audit,
    family_screen = family_screen$table,
    copula_screen = copula_screen,
    theta_screen = theta_screen$table,
    selected_model = primary$model,
    contrast_model = contrast$model,
    selected_family = selected_family,
    selected_copula = selected_copula,
    final_fit = final_fit
  )
  utils::write.csv(status, paths$table, row.names = FALSE)

  jss_rand_write_distribution_figure(sample_data, final_fit, paths$distribution_figure)
  jss_rand_write_dependence_figure(sample_data, intercept_fit, contrast_fit, paths$figure)
  jss_rand_write_theta_effect_figure(sample_data, contrast_fit, contrast$model, paths$theta_figure)
  jss_rand_write_se_figure(se_comparison, paths$se_figure)
  jss_rand_write_latex_note(
    paths = paths,
    audit = audit,
    theta_screen = theta_screen$table,
    primary_model = primary$model,
    contrast_model = contrast$model,
    final_fit = final_fit,
    se_comparison = se_comparison
  )

  list(
    module_id = jss_rand_module_id(),
    title = "RAND Health and aging doctor-visits application",
    status = "current",
    data = paths$data,
    tables = c(
      paths$table,
      paths$family_screen,
      paths$copula_screen,
      paths$theta_screen,
      paths$coefficients,
      paths$se_comparison,
      paths$convergence,
      workflow$tables
    ),
    figures = c(
      workflow$figures,
      paths$distribution_figure,
      paths$figure,
      paths$theta_figure,
      paths$se_figure
    ),
    notes = paths$latex_note
  )
}

jss_rand_module_id <- function() {
  "06-application-rand-doctor-visits"
}

jss_rand_paths <- function(settings) {
  module_id <- jss_rand_module_id()
  list(
    data = file.path(settings$data_dir, paste0(module_id, "-data.csv")),
    table = file.path(settings$tables_dir, paste0(module_id, "-status.csv")),
    family_screen = file.path(settings$tables_dir, paste0(module_id, "-family-screen.csv")),
    copula_screen = file.path(settings$tables_dir, paste0(module_id, "-copula-screen.csv")),
    theta_screen = file.path(settings$tables_dir, paste0(module_id, "-theta-screen.csv")),
    coefficients = file.path(settings$tables_dir, paste0(module_id, "-coefficients.csv")),
    se_comparison = file.path(settings$tables_dir, paste0(module_id, "-standard-error-comparison.csv")),
    convergence = file.path(settings$tables_dir, paste0(module_id, "-convergence.csv")),
    data_shape = file.path(settings$tables_dir, paste0(module_id, "-data-shape.csv")),
    missingness_response = file.path(settings$tables_dir, paste0(module_id, "-missingness-response.csv")),
    missingness_predictors = file.path(settings$tables_dir, paste0(module_id, "-missingness-predictors.csv")),
    missingness_terms = file.path(settings$tables_dir, paste0(module_id, "-missingness-terms.csv")),
    time_moments = file.path(settings$tables_dir, paste0(module_id, "-time-moments.csv")),
    pairwise_correlation = file.path(settings$tables_dir, paste0(module_id, "-pairwise-correlation.csv")),
    check_basic = file.path(settings$tables_dir, paste0(module_id, "-check-model-basic.csv")),
    check_all = file.path(settings$tables_dir, paste0(module_id, "-check-model-all.csv")),
    check_scores = file.path(settings$tables_dir, paste0(module_id, "-check-model-scores.csv")),
    check_tail = file.path(settings$tables_dir, paste0(module_id, "-check-model-tail.csv")),
    check_residual_dependence = file.path(settings$tables_dir, paste0(module_id, "-check-model-residual-dependence.csv")),
    copula_diagnostic_quartiles = file.path(settings$tables_dir, paste0(module_id, "-copula-diagnostic-quartiles.csv")),
    copula_diagnostic_kendall = file.path(settings$tables_dir, paste0(module_id, "-copula-diagnostic-kendall.csv")),
    copula_diagnostic_tail = file.path(settings$tables_dir, paste0(module_id, "-copula-diagnostic-tail.csv")),
    copula_diagnostic_conditional_tail = file.path(settings$tables_dir, paste0(module_id, "-copula-diagnostic-conditional-tail.csv")),
    copula_diagnostic_residual_lag = file.path(settings$tables_dir, paste0(module_id, "-copula-diagnostic-residual-lag.csv")),
    benchmark_status = file.path(settings$tables_dir, paste0(module_id, "-benchmark-status.csv")),
    benchmark_results = file.path(settings$tables_dir, paste0(module_id, "-benchmark-results.csv")),
    benchmark_coefficients = file.path(settings$tables_dir, paste0(module_id, "-benchmark-coefficients.csv")),
    prediction = file.path(settings$tables_dir, paste0(module_id, "-prediction-example.csv")),
    marginal_effects = file.path(settings$tables_dir, paste0(module_id, "-marginal-effects.csv")),
    simulation_summary = file.path(settings$tables_dir, paste0(module_id, "-simulation-summary.csv")),
    workflow_status = file.path(settings$tables_dir, paste0(module_id, "-workflow-status.csv")),
    distribution_matrix_figure = file.path(settings$figures_dir, paste0(module_id, "-distribution-matrix.png")),
    margin_overlay_figure = file.path(settings$figures_dir, paste0(module_id, "-margin-overlay.png")),
    margin_by_time_figure = file.path(settings$figures_dir, paste0(module_id, "-margin-by-time.png")),
    screened_copula_figure = file.path(settings$figures_dir, paste0(module_id, "-screened-copula-fit.png")),
    term_plot_figure = file.path(settings$figures_dir, paste0(module_id, "-term-plots.png")),
    marginal_diagnostics_figure = file.path(settings$figures_dir, paste0(module_id, "-marginal-diagnostics.png")),
    copula_diagnostics_figure = file.path(settings$figures_dir, paste0(module_id, "-copula-diagnostics.png")),
    simulation_figure = file.path(settings$figures_dir, paste0(module_id, "-simulation-check.png")),
    distribution_figure = file.path(settings$figures_dir, paste0(module_id, "-distribution-fit.png")),
    figure = file.path(settings$figures_dir, paste0(module_id, "-dependence-fit.png")),
    theta_figure = file.path(settings$figures_dir, paste0(module_id, "-theta-effects.png")),
    se_figure = file.path(settings$figures_dir, paste0(module_id, "-standard-errors.png")),
    latex_note = file.path(settings$root, "paper", "notes", paste0(module_id, "-results.tex"))
  )
}

jss_rand_load_data <- function(settings) {
  local_path <- file.path(settings$root, "results", "rand_hrs_doctor_visits", "rand_hrs_doctor_visits_long.rds")
  env_path <- Sys.getenv("GAMLSS_LONGITUDINAL_RAND_DATA", unset = "")
  path <- if (file.exists(local_path)) local_path else env_path
  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "RAND doctor-visits data were not found. Expected local RDS at ",
      local_path,
      " or a valid GAMLSS_LONGITUDINAL_RAND_DATA path.",
      call. = FALSE
    )
  }
  if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    readRDS(path)
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE)
  }
}

jss_rand_prepare_data <- function(data) {
  required <- c(
    "hhidpn", "wave", "survey_year", "doctor_visits", "hospital_stay",
    "self_rated_health", "hypertension", "diabetes", "cancer",
    "heart_problem", "stroke", "arthritis", "sex", "education", "age",
    "wave_factor"
  )
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("RAND data are missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  data$subject <- as.character(data$hhidpn)
  data$time <- data$wave
  wave_range <- range(data$wave, na.rm = TRUE)
  data$wave_scaled <- if (diff(wave_range) > 0) {
    (data$wave - wave_range[1]) / diff(wave_range)
  } else {
    0
  }
  data$age_scaled <- as.numeric(scale(data$age))

  data$wave_factor <- factor(data$wave_factor, levels = sort(unique(data$wave_factor)))
  data$self_rated_health <- factor(
    data$self_rated_health,
    levels = c("Excellent", "Very good", "Good", "Fair", "Poor"),
    ordered = FALSE
  )
  for (nm in c("hospital_stay", "hypertension", "diabetes", "cancer", "heart_problem", "stroke", "arthritis")) {
    data[[nm]] <- factor(data[[nm]], levels = c("No", "Yes"), ordered = FALSE)
  }
  data$sex <- factor(data$sex, levels = c("Male", "Female"), ordered = FALSE)
  data$education <- factor(data$education, levels = c(
    "Less than high school",
    "GED",
    "High school graduate",
    "Some college",
    "College and above"
  ), ordered = FALSE)
  data$doctor_visits <- as.integer(data$doctor_visits)

  list(data = data)
}

jss_rand_model_covariates <- function(include_broad = TRUE) {
  base <- c(
    "wave_factor", "self_rated_health", "hospital_stay", "hypertension",
    "cancer", "education", "sex", "arthritis"
  )
  if (!isTRUE(include_broad)) {
    return(base)
  }
  unique(c(base, "diabetes", "heart_problem", "stroke", "age_scaled", "wave_scaled"))
}

jss_rand_sample_data <- function(data, settings) {
  sample_n <- if (identical(settings$profile, "expanded")) 200L else 40L
  sample_n <- as.integer(Sys.getenv("GAMLSS_LONGITUDINAL_RAND_SAMPLE_N", unset = sample_n))
  sample_n <- max(1L, sample_n)

  covariates <- jss_rand_model_covariates(include_broad = identical(settings$profile, "expanded"))
  row_complete <- stats::complete.cases(data[, covariates, drop = FALSE])
  subject_complete <- tapply(row_complete, data$subject, all)
  observed_counts <- tapply(!is.na(data$doctor_visits) & row_complete, data$subject, sum)
  eligible <- sort(names(subject_complete)[subject_complete & observed_counts >= 8L])
  if (length(eligible) == 0L) {
    stop("No RAND subjects met the complete-covariate sampling rule.", call. = FALSE)
  }

  selected <- jss_rand_deterministic_subject_sample(eligible, sample_n)
  out <- data[data$subject %in% selected, , drop = FALSE]
  out <- out[order(out$subject, out$time), , drop = FALSE]
  out <- droplevels(out)
  attr(out, "eligible_subjects") <- length(eligible)
  attr(out, "requested_subjects") <- sample_n
  attr(out, "sampling_rule") <- "deterministic equally spaced eligible subject IDs; complete model covariates; at least 8 observed doctor-visit responses"
  rownames(out) <- NULL
  out
}

jss_rand_deterministic_subject_sample <- function(eligible, n) {
  if (length(eligible) <= n) {
    return(eligible)
  }
  idx <- unique(pmax(1L, pmin(length(eligible), round(seq(1, length(eligible), length.out = n)))))
  if (length(idx) < n) {
    idx <- sort(unique(c(idx, seq_len(length(eligible)))))
    idx <- idx[seq_len(n)]
  }
  eligible[idx]
}

jss_rand_fit_budget <- function(settings) {
  if (identical(settings$profile, "expanded")) {
    return(list(max_outer_iter = 60L, max_inner_iter = 60L))
  }
  list(max_outer_iter = 12L, max_inner_iter = 12L)
}

jss_rand_mu_formula <- function() {
  doctor_visits ~ wave_factor + self_rated_health + hospital_stay + hypertension + cancer + education
}

jss_rand_sigma_formula <- function() {
  ~hypertension + self_rated_health + sex
}

jss_rand_fit_longitudinal <- function(data, copula_dist = "C", theta_formula = ~1,
                                      margin_family = "NBI",
                                      fit_budget = jss_rand_fit_budget(list(profile = "smoke")),
                                      compute_vcov = FALSE) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(gamlss.longitudinal::gamlss_longitudinal(
    dataset = data,
    margin_dist = jss_rand_margin_dist(margin_family),
    copula_dist = copula_dist,
    time_var = "time",
    subject_var = "subject",
    mu.formula = jss_rand_mu_formula(),
    sigma.formula = jss_rand_sigma_formula(),
    theta.formula = theta_formula,
    include_dlcopdpar = FALSE,
    max_outer_iter = fit_budget$max_outer_iter,
    max_inner_iter = fit_budget$max_inner_iter,
    verbose = 0,
    compute_vcov = compute_vcov
  )))
  out
}

jss_rand_margin_dist <- function(margin_family) {
  switch(margin_family,
    NBI = gamlss.dist::NBI(),
    NBII = gamlss.dist::NBII(),
    ZIP = gamlss.dist::ZIP(),
    PO = gamlss.dist::PO(),
    stop("Unsupported RAND longitudinal margin family: ", margin_family, call. = FALSE)
  )
}

jss_rand_data_audit <- function(data, sample_data, settings) {
  wave_rows <- do.call(rbind, lapply(split(data, data$wave), function(x) {
    y <- x$doctor_visits
    y_obs <- y[is.finite(y)]
    data.frame(
      scope = "full_data",
      profile = settings$profile,
      wave = x$wave[[1L]],
      survey_year = x$survey_year[[1L]],
      n_rows = nrow(x),
      n_observed = length(y_obs),
      missing_response_pct = 100 * mean(!is.finite(y)),
      mean = mean(y_obs, na.rm = TRUE),
      variance = stats::var(y_obs, na.rm = TRUE),
      zero_pct = 100 * mean(y_obs == 0),
      skewness = jss_rand_skewness(y_obs),
      n_subjects = length(unique(data$subject)),
      eligible_subjects = attr(sample_data, "eligible_subjects") %||% NA_integer_,
      sampling_rule = attr(sample_data, "sampling_rule") %||% NA_character_,
      sample_subjects = length(unique(sample_data$subject)),
      sample_rows = nrow(sample_data),
      sample_observed = sum(is.finite(sample_data$doctor_visits)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(wave_rows) <- NULL
  wave_rows
}

jss_rand_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || stats::sd(x) == 0) {
    return(NA_real_)
  }
  mean((x - mean(x))^3) / stats::sd(x)^3
}

jss_rand_family_screen <- function(sample_data) {
  model_vars <- unique(c(
    all.vars(jss_rand_mu_formula()),
    all.vars(jss_rand_sigma_formula())
  ))
  observed <- sample_data[is.finite(sample_data$doctor_visits), model_vars, drop = FALSE]
  observed <- stats::na.omit(observed)
  families <- c("NBI", "NBII", "ZIP", "PO")
  rows <- lapply(families, function(fam) {
    fit <- tryCatch(
      {
        out <- NULL
        utils::capture.output(out <- jss_rand_fit_gamlss_family(fam, observed))
        out
      },
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(list(
        row = data.frame(
          family = fam,
          aic = NA_real_,
          bic = NA_real_,
          converged = FALSE,
          error = conditionMessage(fit),
          stringsAsFactors = FALSE
        ),
        fit = NULL,
        se_table = jss_rand_empty_se_table("marginal_gamlss", conditionMessage(fit))
      ))
    }
    list(
      row = data.frame(
        family = fam,
        aic = stats::AIC(fit),
        bic = stats::BIC(fit),
        converged = isTRUE(fit$converged),
        error = NA_character_,
        stringsAsFactors = FALSE
      ),
      fit = fit,
      se_table = jss_rand_gamlss_se_table(fit, model = "marginal_gamlss")
    )
  })
  out <- do.call(rbind, lapply(rows, `[[`, "row"))
  out <- out[order(out$aic), , drop = FALSE]
  rownames(out) <- NULL
  fits <- stats::setNames(lapply(families, function(fam) rows[[match(fam, families)]]$fit), families)
  se_tables <- stats::setNames(lapply(families, function(fam) rows[[match(fam, families)]]$se_table), families)
  list(table = out, fits = fits, se_tables = se_tables)
}

jss_rand_fit_gamlss_family <- function(fam, observed) {
  control <- gamlss::gamlss.control(n.cyc = 80, trace = FALSE)
  switch(fam,
    NBI = gamlss::gamlss(
      formula = jss_rand_mu_formula(),
      sigma.formula = jss_rand_sigma_formula(),
      family = gamlss.dist::NBI(),
      data = observed,
      trace = FALSE,
      control = control
    ),
    NBII = gamlss::gamlss(
      formula = jss_rand_mu_formula(),
      sigma.formula = jss_rand_sigma_formula(),
      family = gamlss.dist::NBII(),
      data = observed,
      trace = FALSE,
      control = control
    ),
    ZIP = gamlss::gamlss(
      formula = jss_rand_mu_formula(),
      sigma.formula = jss_rand_sigma_formula(),
      family = gamlss.dist::ZIP(),
      data = observed,
      trace = FALSE,
      control = control
    ),
    PO = gamlss::gamlss(
      formula = jss_rand_mu_formula(),
      sigma.formula = jss_rand_sigma_formula(),
      family = gamlss.dist::PO(),
      data = observed,
      trace = FALSE,
      control = control
    ),
    stop("Unsupported RAND family: ", fam, call. = FALSE)
  )
}

jss_rand_select_copula <- function(copula_screen) {
  usable <- copula_screen[copula_screen$converged & is.finite(copula_screen$aic_joint), , drop = FALSE]
  if (nrow(usable) == 0L) {
    return("C")
  }
  usable$copula_code[[which.min(usable$aic_joint)]]
}

jss_rand_copula_screen <- function(sample_data, fit_budget, profile, selected_family) {
  copulas <- if (identical(profile, "expanded")) {
    c(C = "Clayton", N = "Gaussian", F = "Frank")
  } else {
    c(C = "Clayton")
  }
  rows <- lapply(names(copulas), function(code) {
    fit <- tryCatch(
      jss_rand_fit_longitudinal(sample_data, copula_dist = code, theta_formula = ~1, margin_family = selected_family, fit_budget = fit_budget),
      error = function(e) e
    )
    jss_rand_fit_screen_row(fit, model = paste0("copula_", code), family = selected_family, copula_code = code, copula_name = copulas[[code]], theta_formula = "~1")
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$aic_joint), , drop = FALSE]
  rownames(out) <- NULL
  out
}

jss_rand_theta_candidates <- function(profile) {
  candidates <- list(
    theta_intercept = ~1,
    theta_wave = ~wave_scaled,
    theta_arthritis = ~arthritis,
    theta_wave_arthritis = ~wave_scaled + arthritis
  )
  if (identical(profile, "expanded")) {
    candidates <- c(candidates, list(
      theta_diabetes = ~diabetes,
      theta_hospital_stay = ~hospital_stay,
      theta_stroke = ~stroke,
      theta_age = ~age_scaled,
      theta_self_rated_health = ~self_rated_health,
      theta_heart_problem = ~heart_problem
    ))
  }
  candidates
}

jss_rand_theta_screen <- function(sample_data, fit_budget, profile, selected_family, selected_copula) {
  candidates <- jss_rand_theta_candidates(profile)
  fits <- list()
  rows <- lapply(names(candidates), function(model) {
    fit <- tryCatch(
      jss_rand_fit_longitudinal(sample_data, copula_dist = selected_copula, theta_formula = candidates[[model]], margin_family = selected_family, fit_budget = fit_budget),
      error = function(e) e
    )
    if (!inherits(fit, "error")) {
      fits[[model]] <<- fit
    }
    jss_rand_fit_screen_row(
      fit,
      model = model,
      family = selected_family,
      copula_code = selected_copula,
      copula_name = jss_rand_copula_name(selected_copula),
      theta_formula = jss_rand_formula_label(candidates[[model]])
    )
  })
  table <- do.call(rbind, rows)
  table <- table[order(table$aic_joint), , drop = FALSE]
  rownames(table) <- NULL
  list(table = table, fits = fits)
}

jss_rand_fit_screen_row <- function(fit, model, family, copula_code, copula_name, theta_formula) {
  if (inherits(fit, "error")) {
    return(data.frame(
      model = model,
      family = family,
      copula_code = copula_code,
      copula_name = copula_name,
      theta_formula = theta_formula,
      loglik_joint = NA_real_,
      aic_joint = NA_real_,
      bic_joint = NA_real_,
      edf_joint = NA_real_,
      converged = FALSE,
      hit_outer_limit = NA,
      outer_iterations = NA_real_,
      max_outer_iter = NA_real_,
      error = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }
  crit <- jss_rand_fit_criteria(fit)
  data.frame(
    model = model,
    family = family,
    copula_code = copula_code,
    copula_name = copula_name,
    theta_formula = theta_formula,
    loglik_joint = crit["LogLik", "joint"],
    aic_joint = crit["AIC", "joint"],
    bic_joint = crit["BIC", "joint"],
    edf_joint = crit["EDF", "joint"],
    converged = isTRUE(fit$convergence$converged),
    hit_outer_limit = isTRUE(fit$convergence$hit_outer_limit),
    outer_iterations = fit$convergence$outer_iterations %||% NA_real_,
    max_outer_iter = fit$convergence$max_outer_iter %||% NA_real_,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

jss_rand_fit_criteria <- function(fit) {
  summary(fit, include_vcov = FALSE)$fit$model_selection
}

jss_rand_formula_label <- function(formula) {
  paste(deparse(formula), collapse = " ")
}

jss_rand_copula_name <- function(code) {
  names <- c(C = "Clayton", N = "Gaussian", F = "Frank", G = "Gumbel", J = "Joe", t = "t")
  names[[code]] %||% code
}

jss_rand_add_theta_evidence <- function(theta_table) {
  intercept <- theta_table[theta_table$model == "theta_intercept", , drop = FALSE]
  theta_table$delta_loglik_vs_intercept <- NA_real_
  theta_table$delta_aic_vs_intercept <- NA_real_
  theta_table$lr_stat_vs_intercept <- NA_real_
  theta_table$lr_df_vs_intercept <- NA_real_
  theta_table$lr_p_value_vs_intercept <- NA_real_
  theta_table$evidence_label <- NA_character_
  if (nrow(intercept) != 1L || !is.finite(intercept$loglik_joint) || !is.finite(intercept$edf_joint)) {
    return(theta_table)
  }
  theta_table$delta_loglik_vs_intercept <- theta_table$loglik_joint - intercept$loglik_joint[[1L]]
  theta_table$delta_aic_vs_intercept <- theta_table$aic_joint - intercept$aic_joint[[1L]]
  theta_table$lr_stat_vs_intercept <- 2 * theta_table$delta_loglik_vs_intercept
  theta_table$lr_df_vs_intercept <- pmax(0, theta_table$edf_joint - intercept$edf_joint[[1L]])
  has_lr <- is.finite(theta_table$lr_stat_vs_intercept) & theta_table$lr_df_vs_intercept > 0
  theta_table$lr_p_value_vs_intercept[has_lr] <- stats::pchisq(
    theta_table$lr_stat_vs_intercept[has_lr],
    df = theta_table$lr_df_vs_intercept[has_lr],
    lower.tail = FALSE
  )
  theta_table$evidence_label <- ifelse(
    theta_table$model == "theta_intercept",
    "reference",
    ifelse(theta_table$delta_aic_vs_intercept < -2, "aic_support",
      ifelse(theta_table$delta_loglik_vs_intercept > 0, "likelihood_gain_only", "no_improvement")
    )
  )
  theta_table
}

jss_rand_select_primary_fit <- function(fits, theta_table) {
  usable <- theta_table[theta_table$converged & is.finite(theta_table$aic_joint), , drop = FALSE]
  if (nrow(usable) > 0L) {
    row <- usable[which.min(usable$aic_joint), , drop = FALSE]
  } else {
    row <- theta_table[which.min(theta_table$aic_joint), , drop = FALSE]
  }
  model <- row$model[[1L]]
  list(model = model, fit = fits[[model]])
}

jss_rand_select_covariate_contrast <- function(fits, theta_table) {
  usable <- theta_table[
    theta_table$converged &
      is.finite(theta_table$loglik_joint) &
      theta_table$model != "theta_intercept",
    ,
    drop = FALSE
  ]
  if (nrow(usable) == 0L) {
    return(jss_rand_select_primary_fit(fits, theta_table))
  }
  preferred <- usable[usable$model == "theta_wave_arthritis", , drop = FALSE]
  if (nrow(preferred) == 1L && is.finite(preferred$delta_loglik_vs_intercept) && preferred$delta_loglik_vs_intercept > 0) {
    row <- preferred
  } else {
    row <- usable[which.max(usable$delta_loglik_vs_intercept), , drop = FALSE]
  }
  model <- row$model[[1L]]
  list(model = model, fit = fits[[model]])
}

jss_rand_coefficient_table <- function(fit, selected_model) {
  coef <- summary(fit, include_vcov = FALSE)$coefficients
  coef$selected_theta_model <- selected_model
  coef
}

jss_rand_standard_error_comparison <- function(sample_data, marginal_se_table, intercept_fit, final_fit) {
  marginal <- marginal_se_table
  intercept <- jss_rand_longitudinal_se_table(intercept_fit, model = "intercept_copula")
  final <- jss_rand_longitudinal_se_table(final_fit, model = "covariate_copula")
  out <- rbind(marginal, intercept, final)
  out <- out[out$parameter %in% c("mu", "sigma"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

jss_rand_gamlss_se_table <- function(fit, model) {
  if (is.null(fit)) {
    return(jss_rand_empty_se_table(model, "marginal fit unavailable"))
  }
  sm <- tryCatch({
    out <- NULL
    utils::capture.output(out <- suppressWarnings(summary(fit)))
    out
  }, error = function(e) NULL)
  mu <- tryCatch(stats::coef(fit, what = "mu"), error = function(e) NULL)
  sigma <- tryCatch(stats::coef(fit, what = "sigma"), error = function(e) NULL)
  if (is.null(sm) || is.null(mu) || is.null(sigma)) {
    return(jss_rand_empty_se_table(model, "summary coefficients unavailable"))
  }
  mu_idx <- seq_along(mu)
  sigma_idx <- length(mu) + seq_along(sigma)
  rbind(
    jss_rand_gamlss_se_rows(model, "mu", rownames(sm)[mu_idx], sm[mu_idx, , drop = FALSE]),
    jss_rand_gamlss_se_rows(model, "sigma", rownames(sm)[sigma_idx], sm[sigma_idx, , drop = FALSE])
  )
}

jss_rand_gamlss_se_rows <- function(model, parameter, terms, table) {
  se <- as.numeric(table[, "Std. Error"])
  data.frame(
    model = model,
    parameter = parameter,
    term = terms,
    estimate = as.numeric(table[, "Estimate"]),
    std_error = se,
    std_error_status = ifelse(is.finite(se), "available", "unavailable"),
    stringsAsFactors = FALSE
  )
}

jss_rand_longitudinal_se_table <- function(fit, model) {
  coef <- summary(fit, include_vcov = FALSE)$coefficients
  out <- data.frame(
    model = model,
    parameter = coef$parameter,
    term = sub("^[^.]+[.]", "", coef$term),
    estimate = coef$estimate,
    std_error = NA_real_,
    std_error_status = "not_computed_for_exact_discrete_smoke_fit",
    stringsAsFactors = FALSE
  )
  out
}

jss_rand_empty_se_table <- function(model, status) {
  data.frame(
    model = character(0),
    parameter = character(0),
    term = character(0),
    estimate = numeric(0),
    std_error = numeric(0),
    std_error_status = character(0),
    stringsAsFactors = FALSE
  )
}

jss_rand_convergence_table <- function(theta_fits, selected_model) {
  rows <- lapply(names(theta_fits), function(model) {
    fit <- theta_fits[[model]]
    data.frame(
      model = model,
      selected = identical(model, selected_model),
      converged = isTRUE(fit$convergence$converged),
      hit_outer_limit = isTRUE(fit$convergence$hit_outer_limit),
      hit_max_stall = isTRUE(fit$convergence$hit_max_stall),
      stop_reason = fit$convergence$stop_reason %||% NA_character_,
      outer_iterations = fit$convergence$outer_iterations %||% NA_real_,
      max_outer_iter = fit$convergence$max_outer_iter %||% NA_real_,
      outer_log_lik_change = fit$convergence$outer_log_lik_change %||% NA_real_,
      outer_stop_crit = fit$convergence$outer_stop_crit %||% NA_real_,
      method = fit$convergence$method %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_rand_write_workflow_artifacts <- function(full_data, sample_data, family_screen,
                                              copula_screen, theta_screen,
                                              selected_family, selected_copula,
                                              final_fit, contrast_fit, intercept_fit,
                                              paths, settings) {
  status <- list()

  jss_rand_write_csv(jss_rand_data_shape_table(full_data, sample_data, settings), paths$data_shape)
  status <- c(status, list(jss_rand_artifact_status("data_shape", "table", paths$data_shape, "ok")))

  missingness <- jss_rand_safe_eval("missingness_check", function() {
    gamlss.longitudinal::check_missingness(
      full_data,
      response_var = "doctor_visits",
      time_var = "time",
      subject_var = "subject",
      predictors = c(
        "time", "hospital_stay", "self_rated_health", "hypertension", "diabetes",
        "cancer", "heart_problem", "stroke", "arthritis", "sex", "education", "age_scaled"
      )
    )
  })
  if (missingness$ok) {
    jss_rand_write_csv(missingness$value$response, paths$missingness_response)
    jss_rand_write_csv(missingness$value$predictors, paths$missingness_predictors)
    jss_rand_write_csv(missingness$value$terms, paths$missingness_terms)
  } else {
    jss_rand_write_csv(jss_rand_error_table(missingness$error), paths$missingness_response)
    jss_rand_write_csv(jss_rand_error_table(missingness$error), paths$missingness_predictors)
    jss_rand_write_csv(jss_rand_error_table(missingness$error), paths$missingness_terms)
  }
  status <- c(status, list(jss_rand_artifact_status("missingness_check", "table", paths$missingness_response, missingness$status, missingness$error)))

  jss_rand_write_csv(jss_rand_time_moments(sample_data), paths$time_moments)
  jss_rand_write_csv(jss_rand_pairwise_correlation(sample_data), paths$pairwise_correlation)
  status <- c(status, list(
    jss_rand_artifact_status("time_moments", "table", paths$time_moments, "ok"),
    jss_rand_artifact_status("pairwise_correlation", "table", paths$pairwise_correlation, "ok")
  ))

  check <- jss_rand_safe_eval("check_model", function() {
    suppressWarnings(gamlss.longitudinal::check_model(final_fit, include_vcov = FALSE, include_plots = FALSE))
  })
  if (check$ok) {
    jss_rand_write_csv(check$value$basic_checks, paths$check_basic)
    jss_rand_write_csv(check$value$checks, paths$check_all)
    jss_rand_write_csv(check$value$scores, paths$check_scores)
    jss_rand_write_csv(check$value$tail, paths$check_tail)
    jss_rand_write_csv(check$value$residual_dependence, paths$check_residual_dependence)
  } else {
    for (path in c(paths$check_basic, paths$check_all, paths$check_scores, paths$check_tail, paths$check_residual_dependence)) {
      jss_rand_write_csv(jss_rand_error_table(check$error), path)
    }
  }
  status <- c(status, list(jss_rand_artifact_status("check_model", "table", paths$check_basic, check$status, check$error)))

  plot_data <- jss_rand_representative_wave_data(sample_data, max_waves = 6L)
  fig_status <- list(
    jss_rand_save_plot(paths$distribution_matrix_figure, 9, 9, function() {
      gamlss.longitudinal::plot_dist(
        plot_data,
        margin_dist = jss_rand_margin_dist(selected_family),
        subject_var = "subject",
        time_var = "time",
        response_var = "doctor_visits",
        offdiag_scale = "pseudo",
        transform = "normal"
      )
    }, title = "Response and Pairwise Dependence Matrix"),
    jss_rand_save_plot(paths$margin_overlay_figure, 8, 4.8, function() {
      gamlss.longitudinal::plot_margin_fit(final_fit, bins = 30, plot = FALSE)$plot
    }, title = "Final Margin Overlay"),
    jss_rand_save_plot(paths$margin_by_time_figure, 11, 8, function() {
      gamlss.longitudinal::plot_margin_fit(final_fit, bins = 18, by_time = TRUE, plot = FALSE)$plot
    }, title = "Final Margin Overlay by Wave"),
    jss_rand_save_plot(paths$screened_copula_figure, 11, 8, function() {
      gamlss.longitudinal::plot_copula_fit(fit = final_fit, by_time = TRUE, bins = 14, plot = FALSE)$plot
    }, title = "Screened Copula Fit by Adjacent Wave"),
    jss_rand_save_device_plot(paths$term_plot_figure, 12, 8, function() {
      gamlss.longitudinal::plot_terms(final_fit, data = sample_data, ncol = 3)
    }, title = "Final Model Term Plots"),
    jss_rand_save_plot(paths$marginal_diagnostics_figure, 10, 8, function() {
      out <- graphics::plot(final_fit, data = sample_data, randomize = TRUE)
      out$dashboard
    }, title = "Marginal Diagnostics"),
    jss_rand_save_plot(paths$copula_diagnostics_figure, 13, 13, function() {
      diag <- gamlss.longitudinal::plot_copula_diagnostics(final_fit, data = sample_data, plot = FALSE, dashboard_ncol = 3)
      jss_rand_write_copula_diagnostic_tables(diag, paths)
      diag$dashboard
    }, title = "Copula Diagnostics")
  )
  status <- c(status, fig_status)

  if (!file.exists(paths$copula_diagnostic_quartiles)) {
    for (path in c(
      paths$copula_diagnostic_quartiles,
      paths$copula_diagnostic_kendall,
      paths$copula_diagnostic_tail,
      paths$copula_diagnostic_conditional_tail,
      paths$copula_diagnostic_residual_lag
    )) {
      jss_rand_write_csv(jss_rand_error_table("copula diagnostics unavailable"), path)
    }
  }

  benchmark <- jss_rand_write_benchmark_tables(sample_data, final_fit, paths, settings)
  prediction <- jss_rand_write_prediction_tables(sample_data, final_fit, paths)
  simulation <- jss_rand_write_simulation_tables_and_figure(sample_data, final_fit, paths)
  status <- c(status, benchmark$status, prediction$status, simulation$status)

  workflow_status <- do.call(rbind, status)
  workflow_status$profile <- settings$profile
  workflow_status$selected_family <- selected_family
  workflow_status$selected_copula <- jss_rand_copula_name(selected_copula)
  workflow_status$primary_theta_model <- jss_rand_select_primary_fit(theta_screen$fits, theta_screen$table)$model
  jss_rand_write_csv(workflow_status, paths$workflow_status)

  list(
    tables = c(
      paths$data_shape,
      paths$missingness_response,
      paths$missingness_predictors,
      paths$missingness_terms,
      paths$time_moments,
      paths$pairwise_correlation,
      paths$check_basic,
      paths$check_all,
      paths$check_scores,
      paths$check_tail,
      paths$check_residual_dependence,
      paths$copula_diagnostic_quartiles,
      paths$copula_diagnostic_kendall,
      paths$copula_diagnostic_tail,
      paths$copula_diagnostic_conditional_tail,
      paths$copula_diagnostic_residual_lag,
      paths$benchmark_status,
      paths$benchmark_results,
      paths$benchmark_coefficients,
      paths$prediction,
      paths$marginal_effects,
      paths$simulation_summary,
      paths$workflow_status
    ),
    figures = c(
      paths$distribution_matrix_figure,
      paths$margin_overlay_figure,
      paths$margin_by_time_figure,
      paths$screened_copula_figure,
      paths$term_plot_figure,
      paths$marginal_diagnostics_figure,
      paths$copula_diagnostics_figure,
      paths$simulation_figure
    )
  )
}

jss_rand_data_shape_table <- function(full_data, sample_data, settings) {
  subject_time <- paste(full_data$subject, full_data$time, sep = "::")
  sample_subject_time <- paste(sample_data$subject, sample_data$time, sep = "::")
  data.frame(
    profile = settings$profile,
    full_rows = nrow(full_data),
    full_subjects = length(unique(full_data$subject)),
    full_waves = length(unique(full_data$wave)),
    full_duplicate_subject_time_rows = sum(duplicated(subject_time)),
    sample_rows = nrow(sample_data),
    sample_subjects = length(unique(sample_data$subject)),
    sample_waves = length(unique(sample_data$wave)),
    sample_duplicate_subject_time_rows = sum(duplicated(sample_subject_time)),
    response_column = "doctor_visits",
    subject_column = "subject",
    time_column = "time",
    time_ordered_by = "numeric RAND wave",
    stringsAsFactors = FALSE
  )
}

jss_rand_time_moments <- function(data) {
  do.call(rbind, lapply(split(data, data$wave), function(x) {
    y <- x$doctor_visits
    y <- y[is.finite(y)]
    data.frame(
      wave = x$wave[[1L]],
      survey_year = x$survey_year[[1L]],
      n = length(y),
      mean = mean(y),
      variance = stats::var(y),
      skewness = jss_rand_skewness(y),
      kurtosis = jss_rand_kurtosis(y),
      zero_pct = 100 * mean(y == 0),
      stringsAsFactors = FALSE
    )
  }))
}

jss_rand_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L || stats::sd(x) == 0) {
    return(NA_real_)
  }
  mean((x - mean(x))^4) / stats::sd(x)^4
}

jss_rand_pairwise_correlation <- function(data) {
  wide <- stats::reshape(
    data[, c("subject", "wave", "doctor_visits")],
    idvar = "subject",
    timevar = "wave",
    direction = "wide"
  )
  response_cols <- grep("^doctor_visits[.]", names(wide), value = TRUE)
  waves <- sub("^doctor_visits[.]", "", response_cols)
  rows <- list()
  for (i in seq_along(response_cols)) {
    for (j in seq_along(response_cols)) {
      if (j <= i) next
      x <- wide[[response_cols[[i]]]]
      y <- wide[[response_cols[[j]]]]
      keep <- is.finite(x) & is.finite(y)
      rows[[length(rows) + 1L]] <- data.frame(
        wave_left = as.numeric(waves[[i]]),
        wave_right = as.numeric(waves[[j]]),
        lag = as.numeric(waves[[j]]) - as.numeric(waves[[i]]),
        n_pairs = sum(keep),
        pearson = if (sum(keep) >= 3L) suppressWarnings(stats::cor(x[keep], y[keep], method = "pearson")) else NA_real_,
        kendall = if (sum(keep) >= 3L) suppressWarnings(stats::cor(x[keep], y[keep], method = "kendall")) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

jss_rand_representative_wave_data <- function(data, max_waves = 6L) {
  waves <- sort(unique(data$wave))
  if (length(waves) <= max_waves) {
    return(data)
  }
  keep <- unique(round(seq(1, length(waves), length.out = max_waves)))
  droplevels(data[data$wave %in% waves[keep], , drop = FALSE])
}

jss_rand_write_copula_diagnostic_tables <- function(diag, paths) {
  jss_rand_write_csv(diag$quartile_summary %||% data.frame(), paths$copula_diagnostic_quartiles)
  jss_rand_write_csv(diag$kendall_summary %||% data.frame(), paths$copula_diagnostic_kendall)
  jss_rand_write_csv(diag$tail_summary %||% data.frame(), paths$copula_diagnostic_tail)
  jss_rand_write_csv(diag$conditional_tail_summary %||% data.frame(), paths$copula_diagnostic_conditional_tail)
  jss_rand_write_csv(diag$residual_lag_summary %||% data.frame(), paths$copula_diagnostic_residual_lag)
}

jss_rand_write_benchmark_tables <- function(sample_data, final_fit, paths, settings) {
  observed <- sample_data[is.finite(sample_data$doctor_visits), , drop = FALSE]
  observed <- observed[stats::complete.cases(observed[, c("doctor_visits", "wave_factor", "self_rated_health", "hospital_stay", "hypertension", "cancer", "education"), drop = FALSE]), , drop = FALSE]
  comparators <- if (identical(settings$profile, "expanded")) c("glm", "gee", "glmm", "gam", "gamm") else c("glm", "gam")
  status <- gamlss.longitudinal::benchmark_comparator_status()
  status$requested <- status$comparator %in% comparators
  jss_rand_write_csv(status, paths$benchmark_status)
  bench <- jss_rand_safe_eval("benchmark_standard_models", function() {
    suppressWarnings(gamlss.longitudinal::benchmark_standard_models(
      data = observed,
      formula = doctor_visits ~ wave_factor + self_rated_health + hospital_stay + hypertension + cancer + education,
      subject_var = "subject",
      family = stats::poisson(),
      comparators = comparators,
      fit = final_fit,
      fit_name = "gamlss_longitudinal",
      distributional_metrics = TRUE
    ))
  })
  if (bench$ok) {
    jss_rand_write_csv(bench$value$results, paths$benchmark_results)
    jss_rand_write_csv(bench$value$coefficients, paths$benchmark_coefficients)
  } else {
    jss_rand_write_csv(jss_rand_error_table(bench$error), paths$benchmark_results)
    jss_rand_write_csv(jss_rand_error_table(bench$error), paths$benchmark_coefficients)
  }
  list(status = list(jss_rand_artifact_status("benchmark_standard_models", "table", paths$benchmark_results, bench$status, bench$error)))
}

jss_rand_write_prediction_tables <- function(sample_data, final_fit, paths) {
  newdata <- jss_rand_reference_prediction_grid(sample_data)
  pred <- jss_rand_safe_eval("predict_mean_confidence", function() {
    suppressWarnings(stats::predict(
      final_fit,
      newdata = newdata,
      type = "mean",
      se.fit = TRUE,
      interval = "confidence"
    ))
  })
  if (!pred$ok) {
    pred <- jss_rand_safe_eval("predict_mean", function() {
      stats::predict(final_fit, newdata = newdata, type = "mean")
    })
  }
  if (pred$ok) {
    out <- if (is.data.frame(pred$value)) pred$value else data.frame(fit = as.numeric(pred$value))
    out <- cbind(newdata[, c("wave", "wave_scaled", "survey_year"), drop = FALSE], out)
    jss_rand_write_csv(out, paths$prediction)
  } else {
    jss_rand_write_csv(jss_rand_error_table(pred$error), paths$prediction)
  }

  effects <- jss_rand_safe_eval("marginal_effects", function() {
    suppressWarnings(gamlss.longitudinal::marginal_effects(
      final_fit,
      newdata = sample_data[is.finite(sample_data$doctor_visits), , drop = FALSE],
      variable = "hospital_stay",
      values = levels(sample_data$hospital_stay),
      parameter = "mu",
      se.fit = FALSE
    ))
  })
  if (effects$ok) {
    jss_rand_write_csv(effects$value, paths$marginal_effects)
  } else {
    jss_rand_write_csv(jss_rand_error_table(effects$error), paths$marginal_effects)
  }
  list(status = list(
    jss_rand_artifact_status("prediction", "table", paths$prediction, pred$status, pred$error),
    jss_rand_artifact_status("marginal_effects", "table", paths$marginal_effects, effects$status, effects$error)
  ))
}

jss_rand_reference_prediction_grid <- function(sample_data) {
  waves <- sort(unique(sample_data$wave))
  modal <- function(x) {
    ux <- unique(x[!is.na(x)])
    ux[which.max(tabulate(match(x, ux)))]
  }
  grid <- data.frame(
    wave = waves,
    time = waves,
    survey_year = as.numeric(vapply(waves, function(w) modal(sample_data$survey_year[sample_data$wave == w]), numeric(1))),
    stringsAsFactors = FALSE
  )
  grid$wave_factor <- factor(as.character(grid$wave), levels = levels(sample_data$wave_factor))
  grid$wave_scaled <- (grid$wave - min(sample_data$wave)) / (max(sample_data$wave) - min(sample_data$wave))
  grid$age_scaled <- 0
  for (nm in c("hospital_stay", "hypertension", "diabetes", "cancer", "heart_problem", "stroke", "arthritis", "sex", "education", "self_rated_health")) {
    val <- modal(sample_data[[nm]])
    grid[[nm]] <- factor(val, levels = levels(sample_data[[nm]]))
  }
  grid$subject <- "reference"
  grid$doctor_visits <- NA_integer_
  grid
}

jss_rand_write_simulation_tables_and_figure <- function(sample_data, final_fit, paths) {
  sim <- jss_rand_safe_eval("simulate", function() {
    stats::simulate(final_fit, nsim = 25, seed = 20260528)
  })
  if (!sim$ok) {
    jss_rand_write_csv(jss_rand_error_table(sim$error), paths$simulation_summary)
    jss_rand_write_placeholder_plot(paths$simulation_figure, "Simulation Check", sim$error)
    return(list(status = list(jss_rand_artifact_status("simulate", "table/figure", paths$simulation_summary, sim$status, sim$error))))
  }
  sim_data <- cbind(sample_data[, c("wave", "doctor_visits"), drop = FALSE], sim$value)
  rows <- lapply(sort(unique(sim_data$wave)), function(w) {
    x <- sim_data[sim_data$wave == w, , drop = FALSE]
    obs <- x$doctor_visits[is.finite(x$doctor_visits)]
    sims <- as.matrix(x[grep("^sim_", names(x)), drop = FALSE])
    data.frame(
      wave = w,
      observed_mean = mean(obs),
      observed_variance = stats::var(obs),
      simulated_mean = mean(sims, na.rm = TRUE),
      simulated_variance = stats::var(as.numeric(sims), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  jss_rand_write_csv(summary, paths$simulation_summary)
  p <- ggplot2::ggplot(summary, ggplot2::aes(x = wave)) +
    ggplot2::geom_line(ggplot2::aes(y = observed_mean, color = "observed"), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = observed_mean, color = "observed"), size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = simulated_mean, color = "simulated"), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = simulated_mean, color = "simulated"), size = 2) +
    ggplot2::labs(x = "RAND wave", y = "Doctor visits", title = "Observed and simulated mean doctor visits") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
  ggplot2::ggsave(paths$simulation_figure, p, width = 7.5, height = 4.8, dpi = 320, bg = "white")
  list(status = list(jss_rand_artifact_status("simulate", "table/figure", paths$simulation_summary, "ok")))
}

jss_rand_safe_eval <- function(name, expr) {
  tryCatch(
    list(name = name, ok = TRUE, status = "ok", value = expr(), error = NA_character_),
    error = function(e) list(name = name, ok = FALSE, status = "failed", value = NULL, error = conditionMessage(e))
  )
}

jss_rand_save_plot <- function(path, width, height, expr, title) {
  res <- jss_rand_safe_eval(title, function() {
    plot <- expr()
    ggplot2::ggsave(path, plot, width = width, height = height, dpi = 320, bg = "white")
    plot
  })
  if (!res$ok) {
    jss_rand_write_placeholder_plot(path, title, res$error)
  }
  jss_rand_artifact_status(title, "figure", path, res$status, res$error)
}

jss_rand_save_device_plot <- function(path, width, height, expr, title) {
  res <- jss_rand_safe_eval(title, function() {
    grDevices::png(path, width = width, height = height, units = "in", res = 320, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    expr()
    invisible(TRUE)
  })
  if (!res$ok) {
    jss_rand_write_placeholder_plot(path, title, res$error)
  }
  jss_rand_artifact_status(title, "figure", path, res$status, res$error)
}

jss_rand_write_placeholder_plot <- function(path, title, message) {
  p <- ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) +
    ggplot2::geom_text(ggplot2::aes(label = message), size = 3.2) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void()
  ggplot2::ggsave(path, p, width = 7, height = 4.5, dpi = 220, bg = "white")
}

jss_rand_write_csv <- function(x, path) {
  if (is.null(x) || !is.data.frame(x) || nrow(x) == 0L) {
    x <- data.frame(note = "no rows", stringsAsFactors = FALSE)
  }
  utils::write.csv(x, path, row.names = FALSE)
}

jss_rand_error_table <- function(message) {
  data.frame(status = "failed", error = message, stringsAsFactors = FALSE)
}

jss_rand_artifact_status <- function(step, type, path, status, error = NA_character_) {
  data.frame(
    workflow_step = step,
    artifact_type = type,
    path = path,
    status = status,
    error = error,
    exists = file.exists(path),
    bytes = if (file.exists(path)) file.info(path)$size else NA_real_,
    stringsAsFactors = FALSE
  )
}

jss_rand_status_table <- function(settings, paths, audit, family_screen, copula_screen,
                                  theta_screen, selected_model, contrast_model, selected_family, selected_copula, final_fit) {
  selected_row <- theta_screen[theta_screen$model == selected_model, , drop = FALSE]
  contrast_row <- theta_screen[theta_screen$model == contrast_model, , drop = FALSE]
  best_row <- theta_screen[which.min(theta_screen$aic_joint), , drop = FALSE]
  data.frame(
    module_id = jss_rand_module_id(),
    profile = settings$profile,
    analysis_state = "current",
    sampling_rule = unique(audit$sampling_rule)[1],
    eligible_subjects = unique(audit$eligible_subjects)[1],
    n_sample_subjects = unique(audit$sample_subjects)[1],
    n_sample_rows = unique(audit$sample_rows)[1],
    n_sample_observed = unique(audit$sample_observed)[1],
    selected_family = selected_family,
    selected_copula = jss_rand_copula_name(selected_copula),
    selected_copula_code = selected_copula,
    primary_theta_model = selected_model,
    primary_theta_formula = selected_row$theta_formula[[1L]],
    primary_aic_joint = selected_row$aic_joint[[1L]],
    covariate_contrast_model = contrast_model,
    covariate_contrast_formula = contrast_row$theta_formula[[1L]],
    covariate_delta_loglik_vs_intercept = contrast_row$delta_loglik_vs_intercept[[1L]],
    covariate_delta_aic_vs_intercept = contrast_row$delta_aic_vs_intercept[[1L]],
    covariate_lr_p_value_vs_intercept = contrast_row$lr_p_value_vs_intercept[[1L]],
    best_aic_theta_model = best_row$model[[1L]],
    best_aic_theta_formula = best_row$theta_formula[[1L]],
    best_aic_joint = best_row$aic_joint[[1L]],
    primary_converged = isTRUE(final_fit$convergence$converged),
    primary_hit_outer_limit = isTRUE(final_fit$convergence$hit_outer_limit),
    primary_outer_iterations = final_fit$convergence$outer_iterations %||% NA_real_,
    canonical_data = paths$data,
    canonical_table = paths$table,
    canonical_figure = paths$figure,
    stringsAsFactors = FALSE
  )
}

jss_rand_empirical_adjacent_tau <- function(data) {
  waves <- sort(unique(data$wave))
  rows <- lapply(seq_len(length(waves) - 1L), function(i) {
    w1 <- waves[[i]]
    w2 <- waves[[i + 1L]]
    left <- data[data$wave == w1, c("subject", "doctor_visits"), drop = FALSE]
    right <- data[data$wave == w2, c("subject", "doctor_visits"), drop = FALSE]
    names(left)[2] <- "y1"
    names(right)[2] <- "y2"
    pair <- merge(left, right, by = "subject", all = FALSE)
    pair <- pair[is.finite(pair$y1) & is.finite(pair$y2), , drop = FALSE]
    data.frame(
      wave = w1,
      next_wave = w2,
      time_mid = (w1 + w2) / 2,
      n_pairs = nrow(pair),
      empirical_tau = if (nrow(pair) >= 4L) suppressWarnings(stats::cor(pair$y1, pair$y2, method = "kendall")) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

jss_rand_theta_grid <- function(data, fit, model_name) {
  waves <- sort(unique(data$wave))
  grid <- expand.grid(
    wave = waves,
    arthritis = levels(data$arthritis),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$wave_scaled <- (grid$wave - min(data$wave)) / (max(data$wave) - min(data$wave))
  grid$arthritis <- factor(grid$arthritis, levels = levels(data$arthritis))
  grid$theta <- jss_rand_predict_theta(fit, grid)
  grid$tau <- jss_rand_theta_to_tau(fit$copula_dist, grid$theta)
  grid$model <- model_name
  grid
}

jss_rand_predict_theta <- function(fit, newdata) {
  theta_par <- fit$par[grepl("^theta[.]", names(fit$par))]
  if (length(theta_par) == 0L) {
    return(rep(NA_real_, nrow(newdata)))
  }
  formula <- fit$formulas_int$theta %||% fit$formulas$theta
  X <- stats::model.matrix(formula, newdata)
  col_key <- paste0("theta.", colnames(X))
  col_key <- sub("theta.[(]Intercept[)]", "theta.intercept", col_key, fixed = FALSE)
  beta <- theta_par[col_key]
  beta[is.na(beta)] <- 0
  eta <- as.numeric(X %*% beta)
  copula_info <- getFromNamespace("get_copula_dist", "gamlss.longitudinal")(fit$copula_dist)
  copula_info$copula_link$theta.linkinv(eta)
}

jss_rand_theta_to_tau <- function(copula_dist, theta) {
  tau_fun <- getFromNamespace(".copula_par_to_tau", "gamlss.longitudinal")
  suppressWarnings(tau_fun(copula_dist, theta))
}

jss_rand_write_distribution_figure <- function(sample_data, fit, path) {
  observed <- sample_data[is.finite(sample_data$doctor_visits), , drop = FALSE]
  fitted_mean <- tryCatch(stats::predict(fit, type = "mean"), error = function(e) rep(NA_real_, nrow(sample_data)))
  sample_data$fitted_mean <- fitted_mean
  fitted_by_wave <- stats::aggregate(fitted_mean ~ wave, data = sample_data, FUN = mean, na.rm = TRUE)
  p <- ggplot2::ggplot(observed, ggplot2::aes(x = factor(wave), y = doctor_visits)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.18, width = 0.65, fill = "#d9e8ef", color = "#355c68") +
    ggplot2::geom_point(
      data = fitted_by_wave,
      ggplot2::aes(x = factor(wave), y = fitted_mean),
      inherit.aes = FALSE,
      color = "#b23a48",
      size = 2.2
    ) +
    ggplot2::coord_cartesian(ylim = stats::quantile(observed$doctor_visits, c(0, 0.98), na.rm = TRUE)) +
    ggplot2::labs(x = "RAND wave", y = "Doctor visits", title = "Observed doctor visits and fitted mean by wave") +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(path, p, width = 8, height = 4.8, dpi = 320, bg = "white")
}

jss_rand_write_dependence_figure <- function(sample_data, intercept_fit, final_fit, path) {
  empirical <- jss_rand_empirical_adjacent_tau(sample_data)
  intercept_tau <- jss_rand_theta_grid(sample_data, intercept_fit, "intercept copula")
  final_tau <- jss_rand_theta_grid(sample_data, final_fit, "covariate copula")
  fitted <- rbind(intercept_tau, final_tau)
  fitted <- stats::aggregate(tau ~ wave + model, data = fitted, FUN = mean, na.rm = TRUE)
  p <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = empirical,
      ggplot2::aes(x = factor(wave), y = empirical_tau),
      width = 0.55,
      fill = "#c8c6a7",
      alpha = 0.85
    ) +
    ggplot2::geom_line(
      data = fitted,
      ggplot2::aes(x = factor(wave), y = tau, color = model, group = model),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = fitted,
      ggplot2::aes(x = factor(wave), y = tau, color = model),
      size = 2
    ) +
    ggplot2::labs(x = "Left wave in adjacent pair", y = "Kendall's tau", title = "Adjacent-wave dependence: empirical and fitted") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
  ggplot2::ggsave(path, p, width = 8, height = 4.8, dpi = 320, bg = "white")
}

jss_rand_write_theta_effect_figure <- function(sample_data, fit, selected_model, path) {
  grid <- jss_rand_theta_grid(sample_data, fit, selected_model)
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = wave, y = tau, color = arthritis, group = arthritis)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "RAND wave", y = "Fitted Kendall's tau", title = "Covariate effect on fitted dependence") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
  ggplot2::ggsave(path, p, width = 7.5, height = 4.8, dpi = 320, bg = "white")
}

jss_rand_write_se_figure <- function(se_comparison, path) {
  plot_data <- se_comparison[is.finite(se_comparison$std_error), , drop = FALSE]
  if (nrow(plot_data) == 0L) {
    plot_data <- data.frame(
      model = "not available",
      term = "standard errors",
      std_error = 0,
      stringsAsFactors = FALSE
    )
  }
  plot_data <- plot_data[plot_data$parameter == "mu", , drop = FALSE]
  plot_data <- plot_data[seq_len(min(nrow(plot_data), 18L)), , drop = FALSE]
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stats::reorder(term, std_error), y = std_error, fill = model)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Standard error", title = "Available marginal GAMLSS standard errors") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
  ggplot2::ggsave(path, p, width = 8, height = 5.5, dpi = 320, bg = "white")
}

jss_rand_write_latex_note <- function(paths, audit, theta_screen, primary_model, contrast_model, final_fit, se_comparison) {
  primary <- theta_screen[theta_screen$model == primary_model, , drop = FALSE]
  contrast <- theta_screen[theta_screen$model == contrast_model, , drop = FALSE]
  intercept <- theta_screen[theta_screen$model == "theta_intercept", , drop = FALSE]
  check_basic <- tryCatch(utils::read.csv(paths$check_basic, stringsAsFactors = FALSE), error = function(e) data.frame())
  check_result <- if (nrow(check_basic) > 0L && "status" %in% names(check_basic)) {
    paste(check_basic$area, check_basic$status, sep = "=", collapse = "; ")
  } else {
    "check_model output unavailable"
  }
  delta_aic <- if (nrow(intercept) == 1L && nrow(contrast) == 1L) contrast$aic_joint - intercept$aic_joint else NA_real_
  delta_loglik <- if (nrow(intercept) == 1L && nrow(contrast) == 1L) contrast$loglik_joint - intercept$loglik_joint else NA_real_
  lr_p <- if (nrow(contrast) == 1L) contrast$lr_p_value_vs_intercept[[1L]] else NA_real_
  evidence <- if (is.finite(delta_aic) && delta_aic < -2) {
    "supported by AIC"
  } else if (is.finite(delta_loglik) && delta_loglik > 0) {
    "a likelihood-improving exploratory contrast rather than the AIC-preferred model"
  } else {
    "not supported in the smoke sample"
  }
  lines <- c(
    "% Auto-generated by paper/R/06-application-rand-doctor-visits.R",
    "\\paragraph{RAND doctor-visits application.}",
    sprintf(
      paste0(
        "The RAND Health and Retirement Study doctor-visits analysis used a deterministic ",
        "%s-profile sample of %d subjects from %d eligible complete-covariate subjects ",
        "(%d observed responses)."
      ),
      unique(audit$profile)[1],
      unique(audit$sample_subjects)[1],
      unique(audit$eligible_subjects)[1],
      unique(audit$sample_observed)[1]
    ),
    sprintf(
      paste0(
        "The primary model selected by joint AIC used a %s marginal distribution with a %s copula ",
        "and theta formula \\texttt{%s} (joint AIC %.2f)."
      ),
      primary$family[[1L]],
      primary$copula_name[[1L]],
      primary$theta_formula[[1L]],
      primary$aic_joint[[1L]]
    ),
    sprintf(
      paste0(
        "The pre-specified covariate-dependence contrast \\texttt{%s} was %s: ",
        "joint log likelihood changed by %.2f, AIC changed by %.2f, and the nested LR-screen p-value was %.3f."
      ),
      contrast$theta_formula[[1L]],
      evidence,
      delta_loglik,
      delta_aic,
      lr_p
    ),
    sprintf(
      paste0(
        "The primary refit reported convergence as \\texttt{%s} after %s outer iterations, ",
        "using convergence fields read directly from the fitted object."
      ),
      as.character(isTRUE(final_fit$convergence$converged)),
      as.character(final_fit$convergence$outer_iterations %||% NA)
    ),
    sprintf(
      paste0(
        "Following the native workflow, the analysis also records data shape, missingness, ",
        "time-stratified moments, pairwise dependence, screening overlays, \\texttt{check\\_model()}, ",
        "term plots, marginal diagnostics, copula diagnostics, predictions, marginal effects, ",
        "model-based simulations, and benchmark-standard-model comparisons. The compact ",
        "\\texttt{check\\_model()} statuses were: %s."
      ),
      check_result
    ),
    paste0(
      "The associated figures report the observed doctor-visit distribution, adjacent-wave ",
      "Kendall dependence, fitted covariate effects on the selected copula dependence parameter, ",
      "and available marginal-parameter standard errors. Exact-discrete longitudinal ",
      "standard errors use analytical rectangle-likelihood Hessians when available ",
      "and are marked unavailable only when Hessian estimation is singular or intentionally skipped."
    )
  )
  writeLines(lines, paths$latex_note, useBytes = TRUE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}
jss_rand_select_margin_family <- function(family_screen) {
  usable <- family_screen[family_screen$converged & is.finite(family_screen$aic), , drop = FALSE]
  if (nrow(usable) == 0L) {
    return("NBI")
  }
  usable$family[[which.min(usable$aic)]]
}
