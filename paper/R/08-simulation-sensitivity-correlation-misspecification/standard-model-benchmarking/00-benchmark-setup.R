# Shared helpers for the JSS correlation-misspecification standard-model benchmark.
#
# These scripts are intentionally scoped to the module 08 paper workflow and
# results/jss-exploratory/08-simulation-sensitivity-correlation-misspecification.
# Package source files are not modified by this workflow.

bmk_find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(file.path(current, "R"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find repository root from ", start, call. = FALSE)
    }
    current <- parent
  }
}

bmk_repo_root <- bmk_find_repo_root()
bmk_script_dir <- file.path(
  bmk_repo_root,
  "paper",
  "R",
  "08-simulation-sensitivity-correlation-misspecification",
  "standard-model-benchmarking"
)
bmk_output_root <- file.path(
  bmk_repo_root,
  "results",
  "jss-exploratory",
  "08-simulation-sensitivity-correlation-misspecification",
  "standard-model-benchmarking"
)

bmk_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  if (length(value) == 0L || is.na(value) || !nzchar(value)) default else value
}

bmk_env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(bmk_env(name, as.character(default))))
  if (length(value) == 0L || is.na(value) || !is.finite(value)) default else value
}

bmk_env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(bmk_env(name, as.character(default))))
  if (length(value) == 0L || is.na(value) || !is.finite(value)) default else value
}

bmk_env_flag <- function(name, default = FALSE) {
  value <- tolower(bmk_env(name, if (isTRUE(default)) "true" else "false"))
  value %in% c("1", "true", "t", "yes", "y")
}

bmk_env_vector <- function(name, default = character()) {
  value <- bmk_env(name, "")
  if (!nzchar(value)) return(default)
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

bmk_env_int_vector <- function(name, default = integer()) {
  value <- bmk_env_vector(name, character())
  if (length(value) == 0L) return(default)
  out <- suppressWarnings(as.integer(value))
  out <- out[is.finite(out) & out > 0L]
  if (length(out) == 0L) default else unique(out)
}

bmk_timestamp <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

bmk_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

bmk_stable_library <- function() {
  lib <- file.path(bmk_output_root, "_stable-library")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  lib
}

bmk_default_source <- function() {
  tarballs <- Sys.glob(file.path(bmk_repo_root, "gamlss.longitudinal_*.tar.gz"))
  if (length(tarballs) > 0L) "tarball" else "installed"
}

bmk_load_package <- function() {
  source <- bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE", bmk_default_source())
  source <- match.arg(source, c("tarball", "installed", "github", "local"))

  if (identical(source, "tarball")) {
    tarball <- bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_TARBALL", "")
    if (!nzchar(tarball)) {
      tarballs <- Sys.glob(file.path(bmk_repo_root, "gamlss.longitudinal_*.tar.gz"))
      if (length(tarballs) == 0L) {
        stop("No package tarball found; set GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE=installed or github.", call. = FALSE)
      }
      info <- file.info(tarballs)
      tarball <- tarballs[order(info$mtime, decreasing = TRUE)][1L]
    }
    lib <- bmk_stable_library()
    .libPaths(c(lib, .libPaths()))
    utils::install.packages(tarball, lib = lib, repos = NULL, type = "source", quiet = TRUE)
  } else if (identical(source, "github")) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      stop("Package 'remotes' is required for GitHub benchmark source.", call. = FALSE)
    }
    lib <- bmk_stable_library()
    .libPaths(c(lib, .libPaths()))
    repo <- bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_GITHUB_REPO", "ahibbert/gamlss.longitudinal")
    ref <- bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_GITHUB_REF", "")
    remotes::install_github(repo, ref = if (nzchar(ref)) ref else "HEAD", lib = lib, upgrade = "never", quiet = TRUE)
  } else if (identical(source, "local")) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Package 'pkgload' is required for local source loading.", call. = FALSE)
    }
    pkgload::load_all(bmk_repo_root, quiet = TRUE)
    return(invisible(source))
  }

  suppressPackageStartupMessages(library(gamlss.longitudinal))
  invisible(source)
}

bmk_require_namespaces <- function(packages, strict = TRUE) {
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (isTRUE(strict) && any(!available)) {
    stop(
      "Missing required package(s): ",
      paste(packages[!available], collapse = ", "),
      call. = FALSE
    )
  }
  data.frame(package = packages, available = available, stringsAsFactors = FALSE)
}

bmk_family_specs <- function(include_binary = TRUE) {
  specs <- list(
    gaussian = list(
      family = "gaussian",
      gamlss_family = "NO",
      label = "Gaussian",
      margin_dist = function() gamlss.dist::NO(),
      standard_family = stats::gaussian(),
      linkinv = identity,
      intercept = 0.2,
      slope = 0.6,
      sigma = 1
    ),
    gamma = list(
      family = "gamma",
      gamlss_family = "GA",
      label = "Gamma",
      margin_dist = function() gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
      standard_family = stats::Gamma(link = "log"),
      linkinv = exp,
      intercept = log(2),
      slope = 0.45,
      sigma = 0.45
    ),
    poisson = list(
      family = "poisson",
      gamlss_family = "PO",
      label = "Poisson",
      margin_dist = function() gamlss.dist::PO(mu.link = "log"),
      standard_family = stats::poisson(),
      linkinv = exp,
      intercept = log(3),
      slope = 0.35,
      sigma = NA_real_
    )
  )

  if (isTRUE(include_binary)) {
    specs$binary <- list(
      family = "binary",
      gamlss_family = "BI",
      label = "Binary",
      margin_dist = function() gamlss.dist::BI(mu.link = "logit"),
      standard_family = stats::binomial(),
      linkinv = stats::plogis,
      intercept = stats::qlogis(0.35),
      slope = 0.75,
      sigma = NA_real_
    )
  }

  specs
}

bmk_scenario_specs <- function() {
  list(
    external_exchangeable_moderate = list(
      scenario = "external_exchangeable_moderate",
      generator = "external",
      correlation = "exchangeable",
      correlation_level = "moderate",
      rho = 0.45,
      tau_base = NA_real_,
      n_time = 4L,
      theta_formula = ~1,
      gee_correlations = c("exchangeable", "ar1", "unstructured")
    ),
    external_exchangeable_high = list(
      scenario = "external_exchangeable_high",
      generator = "external",
      correlation = "exchangeable",
      correlation_level = "high",
      rho = 0.75,
      tau_base = NA_real_,
      n_time = 4L,
      theta_formula = ~1,
      gee_correlations = c("exchangeable", "ar1", "unstructured")
    ),
    external_ar1_moderate = list(
      scenario = "external_ar1_moderate",
      generator = "external",
      correlation = "ar1",
      correlation_level = "moderate",
      rho = 0.55,
      tau_base = NA_real_,
      n_time = 4L,
      theta_formula = ~1,
      gee_correlations = c("ar1", "exchangeable", "unstructured")
    ),
    external_ar1_high = list(
      scenario = "external_ar1_high",
      generator = "external",
      correlation = "ar1",
      correlation_level = "high",
      rho = 0.78,
      tau_base = NA_real_,
      n_time = 4L,
      theta_formula = ~1,
      gee_correlations = c("ar1", "exchangeable", "unstructured")
    ),
    internal_time_varying_high = list(
      scenario = "internal_time_varying_high",
      generator = "internal",
      correlation = "time_varying_adjacent",
      correlation_level = "moderate_high",
      rho = NA_real_,
      tau_base = NA_real_,
      tau_edges = c(0.35, 0.55, 0.70),
      n_time = 4L,
      theta_formula = ~time,
      gee_correlations = c("ar1", "exchangeable", "unstructured")
    ),
    internal_covariate_dependent_high = list(
      scenario = "internal_covariate_dependent_high",
      generator = "internal",
      correlation = "covariate_dependent_adjacent",
      correlation_level = "moderate_high",
      rho = NA_real_,
      tau_base = 0.55,
      tau_effect = 1.1,
      n_time = 4L,
      theta_formula = ~x,
      gee_correlations = c("exchangeable", "ar1", "unstructured")
    )
  )
}

bmk_resize_scenario_time <- function(scenario, n_time) {
  scenario$n_time <- as.integer(n_time)
  if (identical(scenario$correlation, "time_varying_adjacent")) {
    scenario$tau_edges <- seq(0.35, 0.70, length.out = max(1L, scenario$n_time - 1L))
  }
  scenario
}

bmk_expand_scenarios_for_timepoints <- function(scenarios, timepoints) {
  timepoints <- as.integer(timepoints)
  timepoints <- timepoints[is.finite(timepoints) & timepoints > 1L]
  if (length(timepoints) == 0L) return(scenarios)
  out <- list()
  for (nm in names(scenarios)) {
    for (tt in timepoints) {
      scenario <- bmk_resize_scenario_time(scenarios[[nm]], tt)
      out_name <- paste0(nm, "_t", tt)
      scenario$scenario_base <- scenario$scenario
      scenario$scenario <- out_name
      out[[out_name]] <- scenario
    }
  }
  out
}

bmk_select_named <- function(x, keep) {
  if (length(keep) == 0L) return(x)
  missing <- setdiff(keep, names(x))
  if (length(missing) > 0L) {
    stop("Unknown name(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x[keep]
}

bmk_design_covariates <- function(base) {
  subject_index <- as.numeric(base$.sim_subject_index)
  x <- gamlss.longitudinal::sim_rescale01(subject_index) - 0.5
  z_null <- as.numeric((subject_index %% 2L) == 0L)
  z_null <- z_null - mean(z_null)
  data.frame(x = x, z_null = z_null, stringsAsFactors = FALSE)
}

bmk_eta_mu <- function(data, spec) {
  spec$intercept + spec$slope * data$x
}

bmk_margin_params <- function(spec) {
  params <- list(mu = function(data) spec$linkinv(bmk_eta_mu(data, spec)))
  if (identical(spec$gamlss_family, "NO")) {
    params$sigma <- spec$sigma
  } else if (identical(spec$gamlss_family, "GA")) {
    params$sigma <- spec$sigma
  }
  params
}

bmk_add_truth_columns <- function(dat, spec) {
  dat$true_eta_mu <- bmk_eta_mu(dat, spec)
  dat$true_beta_intercept <- spec$intercept
  dat$true_beta_x <- spec$slope
  dat$true_beta_z_null <- 0
  dat
}

bmk_internal_copula_params <- function(scenario) {
  if (identical(scenario$correlation, "time_varying_adjacent")) {
    return(list(tau = scenario$tau_edges))
  }
  if (identical(scenario$correlation, "covariate_dependent_adjacent")) {
    return(list(tau = function(edge_data) {
      stats::plogis(stats::qlogis(scenario$tau_base) + scenario$tau_effect * edge_data$x)
    }))
  }
  list(tau = scenario$tau_base)
}

bmk_simulate_internal <- function(spec, scenario, n, seed) {
  dat <- gamlss.longitudinal::simulate_longitudinal_dataset(
    n = n,
    times = seq_len(scenario$n_time),
    margin_dist = spec$margin_dist(),
    copula_dist = "N",
    margin_params = bmk_margin_params(spec),
    copula_params = bmk_internal_copula_params(scenario),
    covariates = bmk_design_covariates,
    seed = seed,
    include_truth = TRUE,
    u_bounds = if (spec$gamlss_family %in% c("PO", "BI")) c(1e-8, 1 - 1e-8) else NULL
  )
  bmk_add_truth_columns(dat, spec)
}

bmk_external_correlation_matrix <- function(n_time, structure, rho) {
  idx <- seq_len(n_time)
  if (identical(structure, "exchangeable")) {
    out <- matrix(rho, n_time, n_time)
    diag(out) <- 1
    return(out)
  }
  if (identical(structure, "ar1")) {
    return(rho ^ abs(outer(idx, idx, "-")))
  }
  stop("Unsupported external correlation structure: ", structure, call. = FALSE)
}

bmk_simulate_external <- function(spec, scenario, n, seed) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for the external simulator.", call. = FALSE)
  }
  set.seed(seed)
  n_time <- scenario$n_time
  subject <- factor(rep(seq_len(n), each = n_time))
  time <- rep(seq_len(n_time), times = n)
  base <- data.frame(
    .sim_subject_index = rep(seq_len(n), each = n_time),
    .sim_time_index = time,
    subject = subject,
    time = time,
    stringsAsFactors = FALSE
  )
  covars <- bmk_design_covariates(base)
  dat <- cbind(base[c("subject", "time")], covars)
  eta <- bmk_eta_mu(dat, spec)
  mu <- spec$linkinv(eta)

  R <- bmk_external_correlation_matrix(n_time, scenario$correlation, scenario$rho)
  z <- mvtnorm::rmvnorm(n, sigma = R)
  u <- as.numeric(t(stats::pnorm(z)))
  qfun <- get(paste0("q", spec$gamlss_family), envir = asNamespace("gamlss.dist"))
  if (identical(spec$gamlss_family, "NO")) {
    response <- qfun(u, mu = mu, sigma = spec$sigma)
    dat$true_sigma <- spec$sigma
  } else if (identical(spec$gamlss_family, "GA")) {
    response <- qfun(u, mu = mu, sigma = spec$sigma)
    dat$true_sigma <- spec$sigma
  } else if (identical(spec$gamlss_family, "PO")) {
    response <- qfun(u, mu = mu)
  } else if (identical(spec$gamlss_family, "BI")) {
    response <- qfun(u, bd = 1, mu = mu)
  } else {
    stop("Unsupported family for external simulator: ", spec$gamlss_family, call. = FALSE)
  }
  dat$response <- as.numeric(response)
  dat$u <- u
  dat$true_mu <- as.numeric(mu)
  dat$true_theta <- NA_real_
  dat$true_zeta <- NA_real_
  right_edge <- dat$time > min(dat$time)
  dat$true_theta[right_edge] <- scenario$rho
  dat$true_zeta[right_edge] <- 0
  bmk_add_truth_columns(dat[c("subject", "time", "response", "x", "z_null", setdiff(names(dat), c("subject", "time", "response", "x", "z_null")))], spec)
}

bmk_simulate_dataset <- function(spec, scenario, n, seed) {
  if (identical(scenario$generator, "internal")) {
    bmk_simulate_internal(spec, scenario, n = n, seed = seed)
  } else {
    bmk_simulate_external(spec, scenario, n = n, seed = seed)
  }
}

bmk_fit_rs_joint <- function(dat, spec, scenario, max_elapsed_sec = Inf, max_outer_iter = 100L, max_inner_iter = 100L) {
  gamlss.longitudinal::gamlss_longitudinal(
    dataset = dat,
    margin_dist = spec$margin_dist(),
    copula_dist = "N",
    time_var = "time",
    subject_var = "subject",
    mu.formula = response ~ x + z_null,
    sigma.formula = ~1,
    nu.formula = ~1,
    tau.formula = ~1,
    theta.formula = scenario$theta_formula,
    zeta.formula = ~1,
    include_dlcopdpar = TRUE,
    method = "RS",
    start_from = NA,
    warm_start_joint = TRUE,
    compute_vcov = TRUE,
    vcov_method = "analytical",
    max_elapsed_sec = max_elapsed_sec,
    max_outer_iter = max_outer_iter,
    max_inner_iter = max_inner_iter,
    verbose = 0
  )
}

bmk_capture_fit <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

bmk_primary_fit_row <- function(fit_result, elapsed, scenario, spec, rep_id, seed) {
  success <- inherits(fit_result$value, "gamlss.longitudinal")
  converged <- if (success && !is.null(fit_result$value$convergence$converged)) {
    isTRUE(fit_result$value$convergence$converged)
  } else {
    NA
  }
  loglik <- if (success && !is.null(fit_result$value$calc_lik_out_end$log_lik)) {
    fit_result$value$calc_lik_out_end$log_lik
  } else {
    c(marginal = NA_real_, copula = NA_real_, joint = NA_real_)
  }
  data.frame(
    generator = scenario$generator,
    scenario = scenario$scenario,
    correlation = scenario$correlation,
    correlation_level = scenario$correlation_level,
    n_time = scenario$n_time,
    family = spec$family,
    gamlss_family = spec$gamlss_family,
    rep = rep_id,
    seed = seed,
    method = "rs_joint",
    success = success,
    converged = converged,
    elapsed_sec = elapsed,
    marginal_loglik = as.numeric(loglik["marginal"]),
    copula_loglik = as.numeric(loglik["copula"]),
    joint_loglik = as.numeric(loglik["joint"]),
    error = if (success) NA_character_ else conditionMessage(fit_result$value),
    warning = if (length(fit_result$warnings)) paste(fit_result$warnings, collapse = " | ") else NA_character_,
    stringsAsFactors = FALSE
  )
}

bmk_run_primary_fit <- function(dat, spec, scenario, rep_id, seed, max_elapsed_sec, max_outer_iter, max_inner_iter) {
  start <- Sys.time()
  primary_timeout <- bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_PRIMARY_TIMEOUT_SEC", Inf)
  if (is.finite(primary_timeout) && primary_timeout > 0) {
    child <- tryCatch(
      bmk_run_in_callr(
        fun = function(repo_root, package_source, dat, spec, scenario, rep_id, seed, max_elapsed_sec, max_outer_iter, max_inner_iter) {
          setwd(repo_root)
          Sys.setenv(GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE = package_source)
          source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking", "00-benchmark-setup.R"))
          bmk_load_package()
          bmk_require_namespaces(c("gamlss.dist", "mvtnorm"), strict = TRUE)
          start <- Sys.time()
          fit_result <- bmk_capture_fit(
            bmk_fit_rs_joint(
              dat = dat,
              spec = spec,
              scenario = scenario,
              max_elapsed_sec = max_elapsed_sec,
              max_outer_iter = max_outer_iter,
              max_inner_iter = max_inner_iter
            )
          )
          elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
          list(
            fit = if (inherits(fit_result$value, "gamlss.longitudinal")) fit_result$value else NULL,
            status = bmk_primary_fit_row(fit_result, elapsed, scenario, spec, rep_id, seed)
          )
        },
        args = list(
          repo_root = bmk_repo_root,
          package_source = bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE", bmk_default_source()),
          dat = dat,
          spec = spec,
          scenario = scenario,
          rep_id = rep_id,
          seed = seed,
          max_elapsed_sec = max_elapsed_sec,
          max_outer_iter = max_outer_iter,
          max_inner_iter = max_inner_iter
        ),
        timeout = primary_timeout
      ),
      error = function(e) e
    )
    if (!inherits(child, "error")) return(child)
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    return(list(
      fit = NULL,
      status = data.frame(
        generator = scenario$generator,
        scenario = scenario$scenario,
        correlation = scenario$correlation,
        correlation_level = scenario$correlation_level,
        n_time = scenario$n_time,
        family = spec$family,
        gamlss_family = spec$gamlss_family,
        rep = rep_id,
        seed = seed,
        method = "rs_joint",
        success = FALSE,
        converged = NA,
        elapsed_sec = elapsed,
        marginal_loglik = NA_real_,
        copula_loglik = NA_real_,
        joint_loglik = NA_real_,
        error = conditionMessage(child),
        warning = NA_character_,
        stringsAsFactors = FALSE
      )
    ))
  }
  fit_result <- bmk_capture_fit(
    bmk_with_elapsed_limit(
      bmk_fit_rs_joint(
        dat = dat,
        spec = spec,
        scenario = scenario,
        max_elapsed_sec = max_elapsed_sec,
        max_outer_iter = max_outer_iter,
        max_inner_iter = max_inner_iter
      ),
      elapsed_sec = primary_timeout
    )
  )
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  list(
    fit = if (inherits(fit_result$value, "gamlss.longitudinal")) fit_result$value else NULL,
    status = bmk_primary_fit_row(fit_result, elapsed, scenario, spec, rep_id, seed)
  )
}

bmk_comparator_methods <- function() {
  methods <- bmk_env_vector("GAMLSS_LONGITUDINAL_BENCHMARK_COMPARATORS", "glm")
  methods <- setdiff(methods, "gee")
  unique(methods)
}

bmk_prefix_results <- function(results, scenario, spec, rep_id, seed, suffix = NULL) {
  if (nrow(results) == 0L) return(results)
  if (!is.null(suffix)) {
    results$method <- paste0(results$method, suffix)
    results$comparator <- paste0(results$comparator, suffix)
  }
  cbind(
    data.frame(
      generator = scenario$generator,
      scenario = scenario$scenario,
      correlation = scenario$correlation,
      correlation_level = scenario$correlation_level,
      n_time = scenario$n_time,
      family = spec$family,
      gamlss_family = spec$gamlss_family,
      rep = rep_id,
      seed = seed,
      stringsAsFactors = FALSE
    ),
    results,
    row.names = NULL
  )
}

bmk_prefix_coefficients <- function(coefs, scenario, spec, rep_id, seed, suffix = NULL) {
  if (is.null(coefs) || nrow(coefs) == 0L) return(data.frame())
  if (!is.null(suffix)) {
    coefs$method <- paste0(coefs$method, suffix)
  }
  cbind(
    data.frame(
      generator = scenario$generator,
      scenario = scenario$scenario,
      correlation = scenario$correlation,
      correlation_level = scenario$correlation_level,
      n_time = scenario$n_time,
      family = spec$family,
      gamlss_family = spec$gamlss_family,
      rep = rep_id,
      seed = seed,
      stringsAsFactors = FALSE
    ),
    coefs,
    row.names = NULL
  )
}

bmk_bind_rows_fill <- function(...) {
  pieces <- list(...)
  pieces <- unlist(pieces, recursive = FALSE)
  pieces <- pieces[vapply(pieces, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (length(pieces) == 0L) return(data.frame())
  cols <- unique(unlist(lapply(pieces, names), use.names = FALSE))
  pieces <- lapply(pieces, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, pieces)
}

bmk_with_elapsed_limit <- function(expr, elapsed_sec = Inf) {
  elapsed_sec <- suppressWarnings(as.numeric(elapsed_sec)[1L])
  if (!is.finite(elapsed_sec) || elapsed_sec <= 0) return(force(expr))
  old <- setTimeLimit(elapsed = elapsed_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  force(expr)
}

bmk_run_in_callr <- function(fun, args, timeout) {
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("Package 'callr' is required for process-level benchmark timeouts.", call. = FALSE)
  }
  callr::r(
    func = fun,
    args = args,
    timeout = timeout,
    show = FALSE,
    spinner = FALSE
  )
}

bmk_model_df <- function(fit, method) {
  if (is.null(fit)) {
    return(c(mean_df = NA_real_, dependence_df = NA_real_, total_df = NA_real_))
  }
  if (inherits(fit, "gamlss.longitudinal")) {
    fixed_df <- length(fit[["par"]])
    smooth_df <- 0
    if (!is.null(fit[["df_s"]])) {
      smooth_df <- sum(unlist(fit[["df_s"]], use.names = FALSE), na.rm = TRUE)
    }
    dependence_df <- sum(grepl("^theta\\.|^zeta\\.", names(fit[["par"]])))
    return(c(
      mean_df = fixed_df - dependence_df + smooth_df,
      dependence_df = dependence_df,
      total_df = fixed_df + smooth_df
    ))
  }
  coef_df <- length(stats::coef(fit))
  alpha <- if (inherits(fit, "geeglm") && !is.null(fit[["geese"]])) fit[["geese"]][["alpha"]] else numeric()
  dependence_df <- length(alpha)
  c(mean_df = coef_df, dependence_df = dependence_df, total_df = coef_df + dependence_df)
}

bmk_fit_complexity_row <- function(fit, method, scenario, spec, rep_id, seed) {
  df <- bmk_model_df(fit, method)
  data.frame(
    generator = scenario$generator,
    scenario = scenario$scenario,
    correlation = scenario$correlation,
    correlation_level = scenario$correlation_level,
    n_time = scenario$n_time,
    family = spec$family,
    gamlss_family = spec$gamlss_family,
    rep = rep_id,
    seed = seed,
    method = method,
    mean_df = unname(df[["mean_df"]]),
    dependence_df = unname(df[["dependence_df"]]),
    total_df = unname(df[["total_df"]]),
    stringsAsFactors = FALSE
  )
}

bmk_true_pair_dependence <- function(dat, scenario) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  time_values <- sort(unique(dat$time))
  time_lookup <- setNames(seq_along(time_values), as.character(time_values))
  dat$time_idx <- unname(time_lookup[as.character(dat$time)])
  pair_list <- list()
  idx <- 1L

  for (subject_id in unique(dat$subject)) {
    subject_rows <- dat[dat$subject == subject_id, , drop = FALSE]
    subject_rows <- subject_rows[order(subject_rows$time_idx), , drop = FALSE]
    if (nrow(subject_rows) < 2L) next

    for (left_pos in seq_len(nrow(subject_rows) - 1L)) {
      for (right_pos in seq.int(left_pos + 1L, nrow(subject_rows))) {
        left_idx <- subject_rows$time_idx[left_pos]
        right_idx <- subject_rows$time_idx[right_pos]
        lag <- right_idx - left_idx

        true_theta <- NA_real_
        if (identical(scenario$generator, "external") && is.finite(scenario$rho)) {
          true_theta <- if (identical(scenario$correlation, "ar1")) {
            scenario$rho ^ lag
          } else {
            scenario$rho
          }
        } else if ("true_theta" %in% names(subject_rows)) {
          edge_rows <- subject_rows$time_idx > left_idx & subject_rows$time_idx <= right_idx
          edge_theta <- as.numeric(subject_rows$true_theta[edge_rows])
          if (length(edge_theta) > 0L && all(is.finite(edge_theta))) {
            true_theta <- prod(edge_theta)
          }
        }

        pair_list[[idx]] <- data.frame(
          subject = subject_id,
          time_left = as.character(subject_rows$time[left_pos]),
          time_right = as.character(subject_rows$time[right_pos]),
          time_left_idx = left_idx,
          time_right_idx = right_idx,
          lag = lag,
          true_theta = true_theta,
          true_zeta = 0,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  if (length(pair_list) == 0L) return(data.frame())
  out <- do.call(rbind, pair_list)
  out$true_tau <- (2 / pi) * asin(pmax(-1, pmin(1, out$true_theta)))
  out
}

bmk_gee_pair_alpha <- function(fit, truth, n_time, corstr) {
  if (is.null(fit) || !inherits(fit, "geeglm") || is.null(fit[["geese"]])) return(rep(NA_real_, nrow(truth)))
  alpha <- suppressWarnings(as.numeric(fit[["geese"]][["alpha"]]))
  if (length(alpha) == 0L) return(rep(NA_real_, nrow(truth)))
  if (identical(corstr, "exchangeable")) {
    return(rep(alpha[[1L]], nrow(truth)))
  }
  if (identical(corstr, "ar1")) {
    return(alpha[[1L]] ^ truth$lag)
  }
  if (identical(corstr, "unstructured")) {
    pair_grid <- utils::combn(seq_len(n_time), 2L)
    pair_key <- paste(pair_grid[1L, ], pair_grid[2L, ], sep = "::")
    truth_key <- paste(truth$time_left_idx, truth$time_right_idx, sep = "::")
    pair_index <- match(truth_key, pair_key)
    out <- rep(NA_real_, nrow(truth))
    valid <- is.finite(pair_index) & pair_index <= length(alpha)
    out[valid] <- alpha[pair_index[valid]]
    return(out)
  }
  rep(NA_real_, nrow(truth))
}

bmk_fitted_pair_theta_gamlss <- function(fit, truth) {
  dep <- tryCatch(gamlss.longitudinal::copula_time_summary(fit), error = function(e) NULL)
  if (is.null(dep) || !is.data.frame(dep$fit_data) || nrow(dep$fit_data) == 0L) {
    return(rep(NA_real_, nrow(truth)))
  }

  fit_data <- dep$fit_data
  time_values <- sort(unique(fit_data$time))
  time_lookup <- setNames(seq_along(time_values), as.character(time_values))
  fit_data$time_idx <- unname(time_lookup[as.character(fit_data$time)])

  out <- rep(NA_real_, nrow(truth))
  for (i in seq_len(nrow(truth))) {
    subject_rows <- fit_data[as.character(fit_data$subject) == as.character(truth$subject[i]), , drop = FALSE]
    if (nrow(subject_rows) == 0L) next
    edge_rows <- subject_rows$time_idx >= truth$time_left_idx[i] & subject_rows$time_idx < truth$time_right_idx[i]
    edge_theta <- as.numeric(subject_rows$theta_fit[edge_rows])
    if (length(edge_theta) > 0L && all(is.finite(edge_theta))) {
      out[i] <- prod(edge_theta)
    }
  }
  out
}

bmk_dependence_recovery_row <- function(dat, fit, method, scenario, spec, rep_id, seed, corstr = NA_character_) {
  truth <- bmk_true_pair_dependence(dat, scenario)
  if (nrow(truth) == 0L) return(data.frame())

  fitted_theta <- rep(NA_real_, nrow(truth))
  fitted_tau <- rep(NA_real_, nrow(truth))
  estimand <- "none"

  if (inherits(fit, "gamlss.longitudinal")) {
    fitted_theta <- bmk_fitted_pair_theta_gamlss(fit, truth)
    if (any(is.finite(fitted_theta))) {
      fitted_tau <- (2 / pi) * asin(pmax(-1, pmin(1, fitted_theta)))
      estimand <- "all_pair_fitted_copula"
    }
  } else if (inherits(fit, "geeglm")) {
    fitted_theta <- bmk_gee_pair_alpha(fit, truth, scenario$n_time, corstr)
    fitted_tau <- (2 / pi) * asin(pmax(-1, pmin(1, fitted_theta)))
    estimand <- "all_pair_working_correlation"
  }

  theta_err <- fitted_theta - truth$true_theta
  tau_err <- fitted_tau - truth$true_tau
  data.frame(
    generator = scenario$generator,
    scenario = scenario$scenario,
    correlation = scenario$correlation,
    correlation_level = scenario$correlation_level,
    n_time = scenario$n_time,
    family = spec$family,
    gamlss_family = spec$gamlss_family,
    rep = rep_id,
    seed = seed,
    method = method,
    dependence_estimand = estimand,
    dependence_scope = "all_pairs",
    dependence_n = sum(is.finite(theta_err) | is.finite(tau_err)),
    theta_mean_error = mean(theta_err, na.rm = TRUE),
    theta_mae = mean(abs(theta_err), na.rm = TRUE),
    theta_rmse = sqrt(mean(theta_err^2, na.rm = TRUE)),
    tau_mean_error = mean(tau_err, na.rm = TRUE),
    tau_mae = mean(abs(tau_err), na.rm = TRUE),
    tau_rmse = sqrt(mean(tau_err^2, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

bmk_failed_gee_result <- function(dat, scenario, spec, rep_id, seed, corstr, elapsed, error, warning = NA_character_) {
  out <- data.frame(
    method = paste0("gee_", corstr),
    comparator = paste0("gee_", corstr),
    comparator_class = "gee",
    estimator = "geepack::geeglm",
    package = "geepack",
    available = requireNamespace("geepack", quietly = TRUE),
    success = FALSE,
    elapsed_sec = elapsed,
    nobs = nrow(dat),
    logLik = NA_real_,
    AIC = NA_real_,
    mae = NA_real_,
    rmse = NA_real_,
    benchmark_mae = NA_real_,
    benchmark_rmse = NA_real_,
    as.list(gamlss.longitudinal:::.benchmark_distribution_metric_empty()),
    warning = warning,
    error = error,
    stringsAsFactors = FALSE
  )
  bmk_prefix_results(out, scenario, spec, rep_id, seed)
}

bmk_run_gee_benchmark <- function(dat, spec, corstr, interval_level, timeout = Inf) {
  if (!is.finite(timeout) || timeout <= 0) {
    return(gamlss.longitudinal::benchmark_standard_models(
      data = dat,
      formula = response ~ x + z_null,
      subject_var = "subject",
      family = spec$standard_family,
      comparators = "gee",
      correlation = corstr,
      truth_family = spec$gamlss_family,
      interval_level = interval_level,
      waves = dat$time
    ))
  }
  bmk_run_in_callr(
    fun = function(repo_root, package_source, dat, spec, corstr, interval_level) {
      setwd(repo_root)
      Sys.setenv(GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE = package_source)
      source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking", "00-benchmark-setup.R"))
      bmk_load_package()
      gamlss.longitudinal::benchmark_standard_models(
        data = dat,
        formula = response ~ x + z_null,
        subject_var = "subject",
        family = spec$standard_family,
        comparators = "gee",
        correlation = corstr,
        truth_family = spec$gamlss_family,
        interval_level = interval_level,
        waves = dat$time
      )
    },
    args = list(
      repo_root = bmk_repo_root,
      package_source = bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE", bmk_default_source()),
      dat = dat,
      spec = spec,
      corstr = corstr,
      interval_level = interval_level
    ),
    timeout = timeout
  )
}

bmk_truth_for_term <- function(term, spec) {
  term <- as.character(term)
  out <- rep(NA_real_, length(term))
  out[term %in% c("intercept", "(Intercept)", "Intercept")] <- spec$intercept
  out[term == "x"] <- spec$slope
  out[term == "z_null"] <- 0
  out
}

bmk_annotate_coefficients <- function(coefs, spec, level = 0.95) {
  if (nrow(coefs) == 0L) return(coefs)
  coefs$truth <- bmk_truth_for_term(coefs$term, spec)
  z <- stats::qnorm((1 + level) / 2)
  if (!"conf.low" %in% names(coefs)) coefs$conf.low <- coefs$estimate - z * coefs$std_error
  if (!"conf.high" %in% names(coefs)) coefs$conf.high <- coefs$estimate + z * coefs$std_error
  coefs$bias <- coefs$estimate - coefs$truth
  coefs$ci_width <- coefs$conf.high - coefs$conf.low
  coefs$ci_covers_truth <- is.finite(coefs$truth) & is.finite(coefs$conf.low) &
    is.finite(coefs$conf.high) & coefs$conf.low <= coefs$truth & coefs$conf.high >= coefs$truth
  coefs$p_value <- 2 * stats::pnorm(abs(coefs$estimate / coefs$std_error), lower.tail = FALSE)
  coefs$false_positive <- is.finite(coefs$truth) & abs(coefs$truth) < 1e-12 &
    is.finite(coefs$p_value) & coefs$p_value < 0.05
  coefs
}

bmk_run_benchmark_standard_models <- function(dat, primary_fit, spec, scenario, rep_id, seed, interval_level = 0.95, primary_elapsed_sec = NA_real_) {
  main <- gamlss.longitudinal::benchmark_standard_models(
    data = dat,
    formula = response ~ x + z_null,
    subject_var = "subject",
    family = spec$standard_family,
    comparators = bmk_comparator_methods(),
    correlation = scenario$gee_correlations[[1L]],
    fit = primary_fit,
    fit_name = "rs_joint",
    truth_family = spec$gamlss_family,
    interval_level = interval_level
  )

  if ("elapsed_sec" %in% names(main$results)) {
    main$results$elapsed_sec[main$results$method == "rs_joint"] <- primary_elapsed_sec
  }

  result_rows <- list(bmk_prefix_results(main$results, scenario, spec, rep_id, seed))
  coef_rows <- list(bmk_prefix_coefficients(main$coefficients$long, scenario, spec, rep_id, seed))
  complexity_rows <- list()
  dependence_rows <- list()

  for (method in names(main$fits)) {
    complexity_rows[[length(complexity_rows) + 1L]] <- bmk_fit_complexity_row(main$fits[[method]], method, scenario, spec, rep_id, seed)
    dependence_rows[[length(dependence_rows) + 1L]] <- bmk_dependence_recovery_row(dat, main$fits[[method]], method, scenario, spec, rep_id, seed)
  }

  for (corstr in scenario$gee_correlations) {
    timeout <- if (identical(corstr, "unstructured")) {
      bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_GEE_UNSTRUCTURED_TIMEOUT_SEC", Inf)
    } else {
      bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_GEE_TIMEOUT_SEC", Inf)
    }
    start_gee <- Sys.time()
    gee_capture <- bmk_capture_fit(
      bmk_run_gee_benchmark(dat, spec, corstr, interval_level, timeout = timeout)
    )
    if (inherits(gee_capture$value, "error")) {
      elapsed <- as.numeric(difftime(Sys.time(), start_gee, units = "secs"))
      result_rows[[length(result_rows) + 1L]] <- bmk_failed_gee_result(
        dat = dat,
        scenario = scenario,
        spec = spec,
        rep_id = rep_id,
        seed = seed,
        corstr = corstr,
        elapsed = elapsed,
        error = conditionMessage(gee_capture$value),
        warning = if (length(gee_capture$warnings)) paste(gee_capture$warnings, collapse = " | ") else NA_character_
      )
      next
    }
    gee <- gee_capture$value
    suffix <- paste0("_", corstr)
    result_rows[[length(result_rows) + 1L]] <- bmk_prefix_results(gee$results, scenario, spec, rep_id, seed, suffix = suffix)
    coef_rows[[length(coef_rows) + 1L]] <- bmk_prefix_coefficients(gee$coefficients$long, scenario, spec, rep_id, seed, suffix = suffix)
    gee_fit <- gee$fits$gee
    gee_method <- paste0("gee_", corstr)
    complexity_rows[[length(complexity_rows) + 1L]] <- bmk_fit_complexity_row(gee_fit, gee_method, scenario, spec, rep_id, seed)
    dependence_rows[[length(dependence_rows) + 1L]] <- bmk_dependence_recovery_row(dat, gee_fit, gee_method, scenario, spec, rep_id, seed, corstr = corstr)
  }

  results <- bmk_bind_rows_fill(result_rows)
  coefficients <- bmk_bind_rows_fill(coef_rows)
  complexity <- bmk_bind_rows_fill(complexity_rows)
  dependence <- bmk_bind_rows_fill(dependence_rows)

  duplicated_primary_gee <- results$method == "gee"
  if (any(duplicated_primary_gee)) {
    results <- results[!duplicated_primary_gee, , drop = FALSE]
  }
  if (nrow(coefficients) > 0L && any(coefficients$method == "gee")) {
    coefficients <- coefficients[coefficients$method != "gee", , drop = FALSE]
  }

  coefficients <- bmk_annotate_coefficients(coefficients, spec, level = interval_level)
  list(results = results, coefficients = coefficients, complexity = complexity, dependence = dependence)
}

bmk_metric_columns <- function(x) {
  grep("^(benchmark_|elapsed_sec$|AIC$|logLik$)", names(x), value = TRUE)
}

bmk_available_methods <- function() {
  status <- gamlss.longitudinal::benchmark_comparator_status()
  status$available[match("glmmTMB", status$comparator)] <-
    requireNamespace("glmmTMB", quietly = TRUE)
  status
}

bmk_session_info <- function(source) {
  data.frame(
    item = c(
      "timestamp",
      "repo_root",
      "package_source",
      "package_version",
      "R_version",
      "output_root"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      bmk_repo_root,
      source,
      as.character(utils::packageVersion("gamlss.longitudinal")),
      as.character(getRversion()),
      bmk_output_root
    ),
    stringsAsFactors = FALSE
  )
}
