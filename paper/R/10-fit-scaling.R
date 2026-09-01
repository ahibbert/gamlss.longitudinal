jss_scaling_module_id <- function() "10-fit-scaling"

jss_scaling_schema_version <- function() 3L

jss_scaling_seed_scheme_version <- function() "scenario-key-v1"

jss_scaling_hash <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for fit-scaling checkpoint fingerprints.", call. = FALSE)
  }
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

jss_scaling_seed <- function(base_seed, scenario_id, attempt_id, stream = "fit") {
  key <- paste(
    jss_scaling_seed_scheme_version(), as.integer(base_seed), scenario_id,
    as.integer(attempt_id), stream, sep = "|"
  )
  bytes <- utf8ToInt(enc2utf8(key))
  # An integer-only rolling hash is stable across R versions and platforms.
  value <- 104729
  for (byte in bytes) value <- (value * 131 + byte) %% 2147483646
  as.integer(value + 1)
}

jss_scaling_precision_rule <- function(profile = "full") {
  profile <- match.arg(profile, c("smoke", "paper", "full"))
  if (identical(profile, "smoke")) {
    return(list(
      min_attempts = 2L, initial_max_attempts = 3L, extension_attempts = 1L,
      max_attempts = 4L, target_relative_mcse = 0.50, bootstrap_draws = 99L,
      warmup_attempts = 1L
    ))
  }
  list(
    min_attempts = 5L, initial_max_attempts = 15L, extension_attempts = 5L,
    max_attempts = 45L, target_relative_mcse = 0.10, bootstrap_draws = 999L,
    warmup_attempts = 1L
  )
}

jss_scaling_design <- function(profile = "full") {
  profile <- match.arg(profile, c("smoke", "paper", "full"))
  base <- if (identical(profile, "smoke")) {
    list(n_subjects = 20L, n_visits = 3L, smooth_k = 4L)
  } else {
    list(n_subjects = 100L, n_visits = 4L, smooth_k = 5L)
  }
  high <- if (identical(profile, "smoke")) {
    list(n_subjects = 30L, n_visits = 4L, smooth_k = 5L)
  } else {
    list(n_subjects = 250L, n_visits = 8L, smooth_k = 10L)
  }
  contrasts <- data.frame(
    scenario_suffix = c("base", "subjects", "visits", "smooth_basis"),
    changed_factor = c("none", "n_subjects", "n_visits", "smooth_k"),
    n_subjects = c(base$n_subjects, high$n_subjects, base$n_subjects, base$n_subjects),
    n_visits = c(base$n_visits, base$n_visits, high$n_visits, base$n_visits),
    smooth_k = c(base$smooth_k, base$smooth_k, base$smooth_k, high$smooth_k),
    stringsAsFactors = FALSE
  )
  families <- data.frame(
    model_class = c("continuous", "discrete"),
    family = c("NO", "NBI"),
    copula = c("N", "C"),
    stringsAsFactors = FALSE
  )
  design <- merge(families, contrasts, by = NULL, sort = FALSE)
  design$scenario_id <- paste(design$model_class, design$scenario_suffix, sep = "-")
  design$profile_design <- if (identical(profile, "smoke")) "smoke" else "authoritative"
  design$primary_comparison <- design$scenario_suffix != "base"
  design[, c(
    "scenario_id", "model_class", "family", "copula", "changed_factor",
    "n_subjects", "n_visits", "smooth_k", "profile_design", "primary_comparison"
  )]
}

jss_scaling_validate_design <- function(design) {
  required <- c(
    "scenario_id", "model_class", "family", "copula", "changed_factor",
    "n_subjects", "n_visits", "smooth_k", "profile_design", "primary_comparison"
  )
  missing <- setdiff(required, names(design))
  if (length(missing)) stop("Scaling design is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(design$scenario_id)) stop("Scaling scenario IDs must be unique.", call. = FALSE)
  if (nrow(design) != 8L) stop("Scaling design must contain exactly the registered 8 cells.", call. = FALSE)
  if (any(!design$changed_factor %in% c("none", "n_subjects", "n_visits", "smooth_k"))) {
    stop("Scaling design has an unregistered changed factor.", call. = FALSE)
  }
  numeric_fields <- c("n_subjects", "n_visits", "smooth_k")
  if (any(vapply(design[numeric_fields], function(x) any(!is.finite(x) | x < 1 | x != as.integer(x)), logical(1)))) {
    stop("Scaling dimensions must be positive integers.", call. = FALSE)
  }
  for (model_class in unique(design$model_class)) {
    block <- design[design$model_class == model_class, , drop = FALSE]
    base <- block[block$changed_factor == "none", , drop = FALSE]
    if (nrow(base) != 1L) stop("Each model class needs exactly one base scenario.", call. = FALSE)
    for (i in which(block$changed_factor != "none")) {
      changed <- numeric_fields[as.numeric(block[i, numeric_fields]) != as.numeric(base[1L, numeric_fields])]
      if (!identical(changed, block$changed_factor[[i]])) {
        stop("Scenario ", block$scenario_id[[i]], " does not change exactly its named factor.", call. = FALSE)
      }
    }
  }
  profile_design <- unique(as.character(design$profile_design))
  if (length(profile_design) != 1L || !profile_design %in% c("smoke", "authoritative")) {
    stop("Scaling design must declare one registered profile_design.", call. = FALSE)
  }
  expected <- jss_scaling_design(if (identical(profile_design, "smoke")) "smoke" else "full")
  compare_fields <- names(expected)
  actual <- design[order(design$scenario_id), compare_fields, drop = FALSE]
  expected <- expected[order(expected$scenario_id), compare_fields, drop = FALSE]
  rownames(actual) <- rownames(expected) <- NULL
  for (field in compare_fields) {
    if (!identical(as.character(actual[[field]]), as.character(expected[[field]]))) {
      stop(
        "Scaling design does not match the registered 8-cell design in field ",
        field, ". Family, copula, dimensions, and completeness are fixed.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

jss_scaling_normalize_design <- function(design) {
  jss_scaling_validate_design(design)
  design <- design[order(design$scenario_id), names(jss_scaling_design("full")), drop = FALSE]
  rownames(design) <- NULL
  design
}

jss_scaling_configuration <- function(profile = "full") {
  list(
    profile = match.arg(profile, c("smoke", "paper", "full")),
    precision_rule = jss_scaling_precision_rule(profile),
    simulation = list(tau = 0.25, nbi_sigma = 0.35, no_sigma = 0.80),
    fit = list(
      method = "RS", outer_tol = 1e-3, max_outer_iter = 50L,
      max_elapsed_sec = 180, inner_tol = 1e-3, max_inner_iter = 50L,
      compute_vcov = FALSE, include_dlcopdpar = TRUE
    )
  )
}

jss_scaling_code_fingerprint <- function() {
  functions <- c(
    "jss_scaling_design", "jss_scaling_precision_rule", "jss_scaling_seed",
    "jss_scaling_configuration", "jss_scaling_simulate", "jss_scaling_fit_attempt",
    "jss_scaling_converged", "jss_scaling_precision_met", "jss_scaling_bootstrap_median",
    "jss_scaling_timing_environment", "jss_scaling_checkpoint_spec",
    "jss_scaling_attempt_order", "jss_scaling_checkpoint_valid", "jss_scaling_run"
  )
  jss_scaling_hash(lapply(functions, function(name) paste(deparse(body(get(name, mode = "function"))), collapse = "\n")))
}

jss_scaling_package_versions <- function() {
  packages <- c("gamlss.longitudinal", "gamlss", "gamlss.dist", "mgcv", "digest")
  versions <- vapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NA_character_
  }, character(1))
  stats::setNames(versions, packages)
}

jss_scaling_git_value <- function(root, field = c("sha", "dirty")) {
  field <- match.arg(field)
  helper <- if (identical(field, "sha")) "jss_git_sha" else "jss_git_dirty"
  if (exists(helper, mode = "function")) return(as.character(get(helper, mode = "function")(root)))
  args <- if (identical(field, "sha")) {
    c("-C", shQuote(root), "rev-parse", "HEAD")
  } else {
    c("-C", shQuote(root), "status", "--porcelain")
  }
  value <- tryCatch(system2("git", args, stdout = TRUE, stderr = FALSE), error = function(e) character())
  if (identical(field, "sha")) {
    if (length(value)) value[[1L]] else NA_character_
  } else {
    if (length(value)) "dirty" else "clean"
  }
}

jss_scaling_timing_environment <- function(
    root = getwd(), execution_context = "standalone_or_unverified") {
  info <- Sys.info()
  versions <- jss_scaling_package_versions()
  list(
    timing_git_sha = jss_scaling_git_value(root, "sha"),
    timing_git_dirty = jss_scaling_git_value(root, "dirty"),
    timing_r_version = R.version.string,
    timing_platform = R.version$platform,
    timing_os = paste(info[["sysname"]], info[["release"]], info[["machine"]]),
    timing_cpu_model = jss_scaling_cpu_model(),
    timing_execution_context = as.character(execution_context),
    timing_dependency_versions = paste(names(versions), versions, sep = "=", collapse = ";")
  )
}

jss_scaling_checkpoint_spec <- function(
    design, profile, base_seed, root = getwd(),
    execution_context = "standalone_or_unverified") {
  normalized <- jss_scaling_normalize_design(design)
  configuration <- jss_scaling_configuration(profile)
  versions <- jss_scaling_package_versions()
  timing_environment <- jss_scaling_timing_environment(root, execution_context)
  spec <- list(
    checkpoint_schema_version = jss_scaling_schema_version(),
    design_fingerprint = jss_scaling_hash(normalized),
    configuration_fingerprint = jss_scaling_hash(configuration),
    code_fingerprint = jss_scaling_code_fingerprint(),
    package_fingerprint = jss_scaling_hash(versions),
    timing_environment_fingerprint = jss_scaling_hash(timing_environment),
    seed_scheme_version = jss_scaling_seed_scheme_version(),
    base_seed = as.integer(base_seed),
    checkpoint_fingerprint = jss_scaling_hash(list(
      schema = jss_scaling_schema_version(), design = normalized,
      configuration = configuration, code = jss_scaling_code_fingerprint(),
      packages = versions, timing_environment = timing_environment,
      seed_scheme = jss_scaling_seed_scheme_version(),
      base_seed = as.integer(base_seed)
    ))
  )
  c(spec, timing_environment)
}

jss_scaling_bootstrap_median <- function(x, draws = 999L, seed = 20260528L) {
  x <- sort(as.numeric(x[is.finite(x)]))
  if (!length(x)) return(c(median = NA_real_, mcse = NA_real_, lower = NA_real_, upper = NA_real_))
  if (length(x) == 1L) return(c(median = x, mcse = NA_real_, lower = x, upper = x))
  set.seed(seed)
  medians <- replicate(as.integer(draws), stats::median(sample(x, length(x), replace = TRUE)))
  c(
    median = stats::median(x),
    mcse = stats::sd(medians),
    lower = unname(stats::quantile(medians, 0.025, names = FALSE)),
    upper = unname(stats::quantile(medians, 0.975, names = FALSE))
  )
}

jss_scaling_precision_met <- function(elapsed, rule, seed = 20260528L) {
  elapsed <- elapsed[is.finite(elapsed)]
  if (length(elapsed) < rule$min_attempts) return(FALSE)
  estimate <- jss_scaling_bootstrap_median(elapsed, rule$bootstrap_draws, seed)
  if (!is.finite(estimate[["mcse"]]) || estimate[["median"]] <= 0) return(FALSE)
  estimate[["mcse"]] / estimate[["median"]] <= rule$target_relative_mcse
}

jss_scaling_simulate <- function(scenario, seed) {
  n_subjects <- as.integer(scenario$n_subjects[[1L]])
  n_visits <- as.integer(scenario$n_visits[[1L]])
  set.seed(seed)
  subject <- data.frame(x = stats::rnorm(n_subjects), stringsAsFactors = FALSE)
  covariates <- subject[rep(seq_len(n_subjects), each = n_visits), , drop = FALSE]
  covariates$time_scaled <- rep(seq(0, 1, length.out = n_visits), times = n_subjects)
  covariates$s1 <- stats::runif(nrow(covariates), -1, 1)
  family <- scenario$family[[1L]]
  margin <- do.call(get(family, envir = asNamespace("gamlss.dist"), mode = "function"), list())
  mu_fun <- if (identical(family, "NBI")) {
    function(data) exp(1 + 0.20 * data$x + 0.15 * data$time_scaled + 0.25 * sin(pi * data$s1))
  } else {
    function(data) 1 + 0.30 * data$x + 0.20 * data$time_scaled + 0.40 * sin(pi * data$s1)
  }
  sigma_value <- if (identical(family, "NBI")) 0.35 else 0.80
  gamlss.longitudinal::simulate_longitudinal_dataset(
    n = n_subjects,
    times = seq_len(n_visits),
    margin_dist = margin,
    copula_dist = scenario$copula[[1L]],
    margin_params = list(mu = mu_fun, sigma = sigma_value),
    copula_params = list(tau = 0.25),
    covariates = covariates,
    seed = seed,
    subject_var = "subject",
    time_var = "time",
    response_var = "response",
    include_truth = FALSE
  )
}

jss_scaling_converged <- function(fit) {
  candidates <- list(fit$converged, fit$optimizer$converged, fit$convergence$converged)
  for (value in candidates) if (length(value) == 1L && !is.na(value)) return(isTRUE(value))
  FALSE
}

jss_scaling_fit_attempt <- function(
    scenario, attempt_id, seed, warmup = FALSE, checkpoint_spec = NULL,
    execution_order = NA_integer_) {
  warnings <- character()
  error <- NA_character_
  fit <- NULL
  dat <- NULL
  formula <- NULL
  margin <- NULL
  stage <- "simulation"
  simulation_elapsed <- system.time({
    dat <- withCallingHandlers(
      tryCatch(
        jss_scaling_simulate(scenario, seed),
        error = function(e) {
          error <<- conditionMessage(e)
          NULL
        }
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })[["elapsed"]]
  setup_elapsed <- NA_real_
  fit_elapsed <- NA_real_
  if (is.null(dat) && is.na(error)) error <- "Simulation did not return a dataset."
  if (!is.null(dat)) {
    stage <- "setup"
    setup_elapsed <- system.time({
      withCallingHandlers(
        tryCatch({
        formula <- stats::as.formula(sprintf(
          "response ~ x + time_scaled + s(s1, k = %d)",
          as.integer(scenario$smooth_k[[1L]])
        ))
        margin <- do.call(
          get(scenario$family[[1L]], envir = asNamespace("gamlss.dist"), mode = "function"),
          list()
        )
        }, error = function(e) {
          error <<- conditionMessage(e)
          NULL
        }),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    })[["elapsed"]]
  }
  if (!is.null(formula) && !is.null(margin) && is.na(error)) {
    stage <- "fit"
    fit_elapsed <- system.time({
      fit <- withCallingHandlers(
        tryCatch(gamlss.longitudinal::gamlss_longitudinal(
          dataset = dat,
          margin_dist = margin,
          copula_dist = scenario$copula[[1L]],
          time_var = "time",
          subject_var = "subject",
          mu.formula = formula,
          sigma.formula = ~1,
          theta.formula = ~1,
          include_dlcopdpar = TRUE,
          method = "RS",
          optimizer_control = gamlss.longitudinal::gamlss_longitudinal_control(
            outer_tol = 1e-3,
            max_outer_iter = 50L,
            max_elapsed_sec = 180,
            rs = list(inner_tol = 1e-3, max_inner_iter = 50L)
          ),
          compute_vcov = FALSE,
          verbose = 0,
          plot_results = FALSE
        ), error = function(e) {
            error <<- conditionMessage(e)
            NULL
          }
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    })[["elapsed"]]
  }
  fitted <- inherits(fit, "gamlss.longitudinal")
  converged <- fitted && jss_scaling_converged(fit)
  failure_reason <- if (converged) {
    "none"
  } else if (!is.na(error) && identical(stage, "simulation")) {
    "simulation_error"
  } else if (!is.na(error) && identical(stage, "setup")) {
    "setup_error"
  } else if (!is.na(error)) {
    "fit_error"
  } else if (fitted) {
    "optimizer_nonconvergence"
  } else {
    "fit_not_returned"
  }
  row <- data.frame(
    scenario_id = scenario$scenario_id[[1L]],
    attempt_id = as.integer(attempt_id),
    execution_order = as.integer(execution_order),
    seed = as.integer(seed),
    model_class = scenario$model_class[[1L]],
    family = scenario$family[[1L]],
    copula = scenario$copula[[1L]],
    changed_factor = scenario$changed_factor[[1L]],
    n_subjects = as.integer(scenario$n_subjects[[1L]]),
    n_visits = as.integer(scenario$n_visits[[1L]]),
    smooth_k = as.integer(scenario$smooth_k[[1L]]),
    n_observations = if (is.null(dat)) NA_integer_ else nrow(dat),
    attempted = TRUE,
    warmup = isTRUE(warmup),
    fit_returned = fitted,
    converged = converged,
    retained = converged && !isTRUE(warmup),
    elapsed_sec = as.numeric(fit_elapsed),
    fit_elapsed_sec = as.numeric(fit_elapsed),
    simulation_elapsed_sec = as.numeric(simulation_elapsed),
    setup_elapsed_sec = as.numeric(setup_elapsed),
    total_elapsed_sec = sum(c(simulation_elapsed, setup_elapsed, fit_elapsed), na.rm = TRUE),
    timing_scope = "fit_only",
    exclusion_reason = if (isTRUE(warmup)) "warmup" else if (!converged) "fit_failure" else "none",
    failure_reason = failure_reason,
    warning_count = length(unique(warnings)),
    warning_messages = paste(unique(warnings), collapse = " | "),
    error_message = error,
    stringsAsFactors = FALSE
  )
  if (!is.null(checkpoint_spec)) {
    for (field in names(checkpoint_spec)) row[[field]] <- checkpoint_spec[[field]]
  }
  row
}

jss_scaling_bind_rows <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_scaling_atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  if (file.exists(path) && !file.remove(path)) stop("Could not replace scaling checkpoint: ", path, call. = FALSE)
  if (!file.rename(temporary, path)) stop("Could not publish scaling checkpoint: ", path, call. = FALSE)
  invisible(path)
}

jss_scaling_checkpoint_valid <- function(existing, scenario, checkpoint_spec) {
  if (!nrow(existing)) return(TRUE)
  required <- c(
    "scenario_id", "attempt_id", "execution_order", "seed", "model_class", "family", "copula",
    "changed_factor", "n_subjects", "n_visits", "smooth_k", "warmup", "retained",
    "elapsed_sec", "fit_elapsed_sec", "simulation_elapsed_sec", "setup_elapsed_sec",
    "total_elapsed_sec", "timing_scope",
    names(checkpoint_spec)
  )
  if (length(setdiff(required, names(existing)))) return(FALSE)
  for (field in names(checkpoint_spec)) {
    if (!all(as.character(existing[[field]]) == as.character(checkpoint_spec[[field]]))) return(FALSE)
  }
  fixed <- c("scenario_id", "model_class", "family", "copula", "changed_factor", "n_subjects", "n_visits", "smooth_k")
  for (field in fixed) {
    if (!all(as.character(existing[[field]]) == as.character(scenario[[field]][[1L]]))) return(FALSE)
  }
  warmup <- existing$warmup %in% TRUE
  if (any(!is.finite(existing$execution_order)) || anyDuplicated(existing$execution_order)) return(FALSE)
  if (!all(existing$timing_scope == "fit_only")) return(FALSE)
  same_fit_time <- (is.na(existing$elapsed_sec) & is.na(existing$fit_elapsed_sec)) |
    (is.finite(existing$elapsed_sec) & is.finite(existing$fit_elapsed_sec) &
      abs(existing$elapsed_sec - existing$fit_elapsed_sec) <= 1e-12)
  if (!all(same_fit_time)) return(FALSE)
  if (anyDuplicated(paste(existing$attempt_id, warmup, sep = "|"))) return(FALSE)
  timed_ids <- sort(as.integer(existing$attempt_id[!warmup]))
  if (length(timed_ids) && !identical(timed_ids, seq_len(length(timed_ids)))) return(FALSE)
  warmup_ids <- sort(as.integer(existing$attempt_id[warmup]), decreasing = TRUE)
  if (length(warmup_ids) && !identical(warmup_ids, 1L - seq_along(warmup_ids))) return(FALSE)
  streams <- ifelse(warmup, "warmup", "fit")
  expected_seeds <- mapply(
    jss_scaling_seed,
    base_seed = checkpoint_spec$base_seed,
    scenario_id = as.character(existing$scenario_id),
    attempt_id = as.integer(existing$attempt_id),
    stream = streams
  )
  identical(as.integer(existing$seed), as.integer(expected_seeds))
}

jss_scaling_checkpoint_contract_valid <- function(attempts, design, profile = "full") {
  if (!nrow(attempts) || !"base_seed" %in% names(attempts)) return(FALSE)
  base_seeds <- unique(suppressWarnings(as.integer(attempts$base_seed)))
  if (length(base_seeds) != 1L || !is.finite(base_seeds)) return(FALSE)
  spec <- tryCatch(
    jss_scaling_checkpoint_spec(design, profile, base_seeds[[1L]]),
    error = function(e) NULL
  )
  if (is.null(spec)) return(FALSE)
  all(vapply(seq_len(nrow(design)), function(i) {
    scenario <- design[i, , drop = FALSE]
    existing <- attempts[attempts$scenario_id == scenario$scenario_id[[1L]], , drop = FALSE]
    nrow(existing) > 0L && jss_scaling_checkpoint_valid(existing, scenario, spec)
  }, logical(1)))
}

jss_scaling_recorded_contract_valid <- function(attempts, design) {
  metadata <- c(
    "checkpoint_schema_version", "design_fingerprint", "configuration_fingerprint",
    "code_fingerprint", "package_fingerprint", "timing_environment_fingerprint",
    "seed_scheme_version", "base_seed", "checkpoint_fingerprint",
    "timing_git_sha", "timing_git_dirty", "timing_r_version", "timing_platform",
    "timing_os", "timing_cpu_model", "timing_execution_context",
    "timing_dependency_versions"
  )
  fixed <- c("scenario_id", "model_class", "family", "copula", "changed_factor", "n_subjects", "n_visits", "smooth_k")
  required <- c(
    metadata, fixed, "attempt_id", "execution_order", "seed", "warmup",
    "elapsed_sec", "fit_elapsed_sec", "simulation_elapsed_sec", "setup_elapsed_sec",
    "total_elapsed_sec", "timing_scope"
  )
  if (!nrow(attempts) || length(setdiff(required, names(attempts)))) return(FALSE)
  values <- lapply(metadata, function(field) unique(as.character(attempts[[field]])))
  if (any(vapply(values, function(value) length(value) != 1L || is.na(value) || !nzchar(value), logical(1)))) {
    return(FALSE)
  }
  if (!identical(as.integer(values[[match("checkpoint_schema_version", metadata)]]), jss_scaling_schema_version())) return(FALSE)
  if (!identical(values[[match("design_fingerprint", metadata)]], jss_scaling_hash(jss_scaling_normalize_design(design)))) return(FALSE)
  if (!identical(values[[match("seed_scheme_version", metadata)]], jss_scaling_seed_scheme_version())) return(FALSE)
  base_seed <- suppressWarnings(as.integer(values[[match("base_seed", metadata)]]))
  if (length(base_seed) != 1L || !is.finite(base_seed)) return(FALSE)
  timing_fields <- c(
    "timing_git_sha", "timing_git_dirty", "timing_r_version", "timing_platform",
    "timing_os", "timing_cpu_model", "timing_execution_context",
    "timing_dependency_versions"
  )
  recorded_environment <- stats::setNames(
    lapply(timing_fields, function(field) values[[match(field, metadata)]]),
    timing_fields
  )
  if (!identical(
    values[[match("timing_environment_fingerprint", metadata)]],
    jss_scaling_hash(recorded_environment)
  )) return(FALSE)
  execution_order <- sort(as.integer(attempts$execution_order))
  if (!identical(execution_order, seq_len(nrow(attempts)))) return(FALSE)
  if (!all(attempts$timing_scope == "fit_only")) return(FALSE)
  same_fit_time <- (is.na(attempts$elapsed_sec) & is.na(attempts$fit_elapsed_sec)) |
    (is.finite(attempts$elapsed_sec) & is.finite(attempts$fit_elapsed_sec) &
      abs(attempts$elapsed_sec - attempts$fit_elapsed_sec) <= 1e-12)
  if (!all(same_fit_time)) return(FALSE)
  all(vapply(seq_len(nrow(design)), function(i) {
    scenario <- design[i, , drop = FALSE]
    existing <- attempts[attempts$scenario_id == scenario$scenario_id[[1L]], , drop = FALSE]
    if (!nrow(existing)) return(FALSE)
    if (any(vapply(fixed, function(field) {
      !all(as.character(existing[[field]]) == as.character(scenario[[field]][[1L]]))
    }, logical(1)))) return(FALSE)
    warmup <- existing$warmup %in% TRUE
    timed_ids <- sort(as.integer(existing$attempt_id[!warmup]))
    warmup_ids <- sort(as.integer(existing$attempt_id[warmup]), decreasing = TRUE)
    if (length(timed_ids) && !identical(timed_ids, seq_len(length(timed_ids)))) return(FALSE)
    if (length(warmup_ids) && !identical(warmup_ids, 1L - seq_along(warmup_ids))) return(FALSE)
    expected_seeds <- mapply(
      jss_scaling_seed, base_seed = base_seed,
      scenario_id = as.character(existing$scenario_id),
      attempt_id = as.integer(existing$attempt_id),
      stream = ifelse(warmup, "warmup", "fit")
    )
    identical(as.integer(existing$seed), as.integer(expected_seeds))
  }, logical(1)))
}

jss_scaling_attempt_order <- function(scenario_ids, round_id, base_seed, stream = "order") {
  keys <- vapply(
    scenario_ids,
    function(id) jss_scaling_seed(base_seed, id, round_id, stream),
    integer(1)
  )
  scenario_ids[order(keys, scenario_ids)]
}

jss_scaling_run <- function(
    design, checkpoint_dir, profile = "full", base_seed = 20260528L,
    root = getwd(), execution_context = "standalone_or_unverified") {
  jss_scaling_validate_design(design)
  rule <- jss_scaling_precision_rule(profile)
  checkpoint_spec <- jss_scaling_checkpoint_spec(
    design, profile, base_seed, root = root, execution_context = execution_context
  )
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  design <- design[order(design$scenario_id), , drop = FALSE]
  scenario_ids <- as.character(design$scenario_id)
  checkpoints <- stats::setNames(vector("list", nrow(design)), scenario_ids)
  for (i in seq_len(nrow(design))) {
    scenario <- design[i, , drop = FALSE]
    path <- file.path(checkpoint_dir, paste0(scenario$scenario_id, ".csv"))
    existing <- if (file.exists(path)) utils::read.csv(path, stringsAsFactors = FALSE) else data.frame()
    if (!jss_scaling_checkpoint_valid(existing, scenario, checkpoint_spec)) {
      stop(
        "Stale or incompatible fit-scaling checkpoint: ", path,
        ". Remove it explicitly before rerunning; checkpoints are bound to schema, design, configuration, code, packages, and seeds.",
        call. = FALSE
      )
    }
    checkpoints[[scenario$scenario_id]] <- existing
  }

  existing_orders <- unlist(lapply(checkpoints, function(x) x$execution_order), use.names = FALSE)
  execution_order <- if (length(existing_orders)) max(as.integer(existing_orders)) else 0L

  warmup_order <- jss_scaling_attempt_order(scenario_ids, 0L, base_seed, "warmup-order")
  for (scenario_id in warmup_order) {
    scenario <- design[design$scenario_id == scenario_id, , drop = FALSE]
    existing <- checkpoints[[scenario_id]]
    completed <- sum(existing$warmup %in% TRUE)
    if (completed < rule$warmup_attempts) {
      for (warmup_index in seq.int(completed + 1L, rule$warmup_attempts)) {
        attempt_id <- 1L - warmup_index
        execution_order <- execution_order + 1L
        row <- jss_scaling_fit_attempt(
          scenario, attempt_id,
          jss_scaling_seed(base_seed, scenario_id, attempt_id, "warmup"),
          warmup = TRUE, checkpoint_spec = checkpoint_spec,
          execution_order = execution_order
        )
        existing <- jss_scaling_bind_rows(list(existing, row))
        checkpoints[[scenario_id]] <- existing
        jss_scaling_atomic_write_csv(existing, file.path(checkpoint_dir, paste0(scenario_id, ".csv")))
      }
    }
  }

  caps <- unique(c(
    seq.int(rule$initial_max_attempts, rule$max_attempts, by = rule$extension_attempts),
    rule$max_attempts
  ))
  for (current_cap in caps) for (round_id in seq_len(current_cap)) {
    active <- scenario_ids[vapply(scenario_ids, function(scenario_id) {
      existing <- checkpoints[[scenario_id]]
      timed <- !existing$warmup %in% TRUE
      elapsed <- existing$elapsed_sec[timed & existing$retained %in% TRUE]
      n_timed <- sum(timed)
      n_timed < min(round_id, current_cap) && !jss_scaling_precision_met(
        elapsed, rule, jss_scaling_seed(base_seed, scenario_id, 0L, "precision")
      )
    }, logical(1))]
    if (!length(active)) next
    active <- jss_scaling_attempt_order(active, round_id, base_seed, "round-order")
    for (scenario_id in active) {
      scenario <- design[design$scenario_id == scenario_id, , drop = FALSE]
      existing <- checkpoints[[scenario_id]]
      next_attempt <- sum(!existing$warmup %in% TRUE) + 1L
      execution_order <- execution_order + 1L
      row <- jss_scaling_fit_attempt(
        scenario, next_attempt,
        jss_scaling_seed(base_seed, scenario_id, next_attempt, "fit"),
        warmup = FALSE, checkpoint_spec = checkpoint_spec,
        execution_order = execution_order
      )
      existing <- jss_scaling_bind_rows(list(existing, row))
      checkpoints[[scenario_id]] <- existing
      jss_scaling_atomic_write_csv(existing, file.path(checkpoint_dir, paste0(scenario_id, ".csv")))
    }
  }
  out <- jss_scaling_bind_rows(checkpoints)
  out <- out[order(out$execution_order), , drop = FALSE]
  rownames(out) <- NULL
  out
}

jss_scaling_summary <- function(attempts, rule = jss_scaling_precision_rule("full"), seed = 20260528L) {
  required <- c("scenario_id", "attempt_id", "attempted", "converged", "retained", "elapsed_sec", "failure_reason")
  missing <- setdiff(required, names(attempts))
  if (length(missing)) stop("Scaling attempts are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  split_rows <- split(attempts, attempts$scenario_id)
  rows <- lapply(names(split_rows), function(scenario_id) {
    x <- split_rows[[scenario_id]]
    warmup <- if ("warmup" %in% names(x)) x$warmup %in% TRUE else rep(FALSE, nrow(x))
    retained <- x$retained %in% TRUE & !warmup
    precision_seed <- jss_scaling_seed(seed, scenario_id, 0L, "precision")
    timing <- jss_scaling_bootstrap_median(x$elapsed_sec[retained], rule$bootstrap_draws, precision_seed)
    q <- if (any(retained)) stats::quantile(x$elapsed_sec[retained], c(0.25, 0.75), names = FALSE) else c(NA_real_, NA_real_)
    first_index <- if (any(!warmup)) which(!warmup)[[1L]] else 1L
    first <- x[first_index, intersect(c("scenario_id", "model_class", "family", "copula", "changed_factor", "n_subjects", "n_visits", "smooth_k", "n_observations"), names(x)), drop = FALSE]
    timed <- !warmup
    precision_met <- jss_scaling_precision_met(x$elapsed_sec[retained], rule, precision_seed)
    attempted <- sum(timed)
    cbind(
      first,
      data.frame(
        warmup_attempted = sum(warmup),
        attempted = attempted,
        converged = sum(x$converged %in% TRUE & timed),
        retained = sum(retained),
        failed = sum(timed & !retained),
        failure_reasons = paste(names(table(x$failure_reason[timed & !retained])), as.integer(table(x$failure_reason[timed & !retained])), sep = ":", collapse = ";"),
        median_elapsed_sec = unname(timing[["median"]]),
        q1_elapsed_sec = q[[1L]],
        q3_elapsed_sec = q[[2L]],
        median_mcse_sec = unname(timing[["mcse"]]),
        median_ci_lower_sec = unname(timing[["lower"]]),
        median_ci_upper_sec = unname(timing[["upper"]]),
        target_relative_mcse = rule$target_relative_mcse,
        initial_attempt_cap = if (!is.null(rule$initial_max_attempts)) rule$initial_max_attempts else rule$max_attempts,
        absolute_attempt_cap = rule$max_attempts,
        extension_used = attempted > if (!is.null(rule$initial_max_attempts)) rule$initial_max_attempts else rule$max_attempts,
        precision_met = precision_met,
        cell_eligible = precision_met && sum(retained) >= rule$min_attempts,
        stringsAsFactors = FALSE
      )
    )
  })
  jss_scaling_bind_rows(rows)
}

jss_scaling_contrasts <- function(attempts, draws = 999L, seed = 20260528L) {
  required <- c("scenario_id", "model_class", "changed_factor", "retained", "elapsed_sec")
  missing <- setdiff(required, names(attempts))
  if (length(missing)) stop("Scaling attempts are missing contrast fields: ", paste(missing, collapse = ", "), call. = FALSE)
  rows <- list()
  out_index <- 0L
  for (model_class in sort(unique(attempts$model_class))) {
    block <- attempts[attempts$model_class == model_class & attempts$retained %in% TRUE, , drop = FALSE]
    base <- sort(block$elapsed_sec[block$changed_factor == "none" & is.finite(block$elapsed_sec)])
    for (factor in c("n_subjects", "n_visits", "smooth_k")) {
      changed <- sort(block$elapsed_sec[block$changed_factor == factor & is.finite(block$elapsed_sec)])
      if (!length(base) || !length(changed)) next
      set.seed(jss_scaling_seed(seed, paste(model_class, factor, sep = "-"), 0L, "contrast"))
      bootstrap_base <- replicate(as.integer(draws), stats::median(sample(base, length(base), replace = TRUE)))
      bootstrap_changed <- replicate(as.integer(draws), stats::median(sample(changed, length(changed), replace = TRUE)))
      differences <- bootstrap_changed - bootstrap_base
      ratios <- bootstrap_changed / bootstrap_base
      out_index <- out_index + 1L
      scenario_id <- unique(block$scenario_id[block$changed_factor == factor])
      rows[[out_index]] <- data.frame(
        model_class = model_class,
        scenario_id = scenario_id[[1L]],
        changed_factor = factor,
        base_attempts = length(base),
        comparison_attempts = length(changed),
        base_median_sec = stats::median(base),
        comparison_median_sec = stats::median(changed),
        median_difference_sec = stats::median(changed) - stats::median(base),
        difference_mcse_sec = stats::sd(differences),
        difference_ci_lower_sec = unname(stats::quantile(differences, 0.025, names = FALSE)),
        difference_ci_upper_sec = unname(stats::quantile(differences, 0.975, names = FALSE)),
        median_runtime_ratio = stats::median(changed) / stats::median(base),
        ratio_mcse = stats::sd(ratios),
        ratio_ci_lower = unname(stats::quantile(ratios, 0.025, names = FALSE)),
        ratio_ci_upper = unname(stats::quantile(ratios, 0.975, names = FALSE)),
        stringsAsFactors = FALSE
      )
    }
  }
  jss_scaling_bind_rows(rows)
}

jss_scaling_reconciliation_tolerance <- function() {
  # CSV round trips preserve these statistics much more closely than this bound.
  list(absolute = 1e-10, relative = 1e-8)
}

jss_scaling_tables_agree <- function(actual, expected, keys) {
  if (!is.data.frame(actual) || !is.data.frame(expected)) return(FALSE)
  if (!setequal(names(actual), names(expected)) || length(names(actual)) != length(names(expected))) return(FALSE)
  if (!all(keys %in% names(actual))) return(FALSE)
  make_key <- function(x) do.call(paste, c(lapply(x[keys], as.character), sep = "\r"))
  actual_key <- make_key(actual)
  expected_key <- make_key(expected)
  if (anyDuplicated(actual_key) || anyDuplicated(expected_key) || !setequal(actual_key, expected_key)) return(FALSE)
  actual <- actual[match(expected_key, actual_key), names(expected), drop = FALSE]
  tolerance <- jss_scaling_reconciliation_tolerance()
  all(vapply(names(expected), function(field) {
    observed <- actual[[field]]
    reference <- expected[[field]]
    if (is.numeric(reference) || is.integer(reference)) {
      observed <- suppressWarnings(as.numeric(observed))
      reference <- as.numeric(reference)
      same_na <- is.na(observed) == is.na(reference)
      finite <- is.finite(observed) & is.finite(reference)
      close <- rep(FALSE, length(reference))
      close[same_na & is.na(reference)] <- TRUE
      close[finite] <- abs(observed[finite] - reference[finite]) <=
        tolerance$absolute + tolerance$relative * abs(reference[finite])
      return(all(same_na & close))
    }
    if (is.logical(reference)) return(identical(as.logical(observed), reference))
    observed <- as.character(observed)
    reference <- as.character(reference)
    observed[is.na(observed)] <- ""
    reference[is.na(reference)] <- ""
    identical(observed, reference)
  }, logical(1)))
}

jss_scaling_attempt_provenance_agrees <- function(attempts, hardware) {
  hardware_fields <- c(
    timing_git_sha = "git_sha", timing_git_dirty = "git_dirty",
    timing_r_version = "r_version", timing_platform = "platform",
    timing_os = "os", timing_cpu_model = "cpu_model",
    timing_execution_context = "execution_context",
    timing_dependency_versions = "dependency_versions"
  )
  required_attempts <- c(names(hardware_fields), "timing_environment_fingerprint")
  if (is.null(hardware) || !is.data.frame(hardware) || nrow(hardware) != 1L ||
      length(setdiff(unname(hardware_fields), names(hardware))) ||
      length(setdiff(required_attempts, names(attempts)))) return(FALSE)
  environment <- stats::setNames(
    lapply(unname(hardware_fields), function(field) as.character(hardware[[field]][[1L]])),
    names(hardware_fields)
  )
  if (!all(vapply(names(hardware_fields), function(field) {
    value <- unique(as.character(attempts[[field]]))
    length(value) == 1L && identical(value, environment[[field]])
  }, logical(1)))) return(FALSE)
  fingerprints <- unique(as.character(attempts$timing_environment_fingerprint))
  length(fingerprints) == 1L && identical(fingerprints, jss_scaling_hash(environment))
}

jss_scaling_publication_status <- function(
    design, attempts, summary, contrasts,
    rule = jss_scaling_precision_rule("full"), hardware = NULL) {
  reasons <- character()
  design_error <- tryCatch({
    jss_scaling_validate_design(design)
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(design_error)) reasons <- c(reasons, paste0("invalid_design: ", design_error))
  expected_ids <- if (is.null(design_error)) sort(as.character(design$scenario_id)) else character()
  actual_attempt_ids <- if ("scenario_id" %in% names(attempts)) sort(unique(as.character(attempts$scenario_id))) else character()
  if (!identical(actual_attempt_ids, expected_ids)) reasons <- c(reasons, "attempt_scenario_incomplete")
  if (is.null(design_error)) {
    if (!jss_scaling_recorded_contract_valid(attempts, design)) {
      reasons <- c(reasons, "checkpoint_contract_invalid")
    }
  }
  if (!jss_scaling_attempt_provenance_agrees(attempts, hardware)) {
    reasons <- c(reasons, "timing_provenance_reconciliation_failed")
  }
  base_seed <- if ("base_seed" %in% names(attempts)) unique(suppressWarnings(as.integer(attempts$base_seed))) else integer()
  recomputed <- if (length(base_seed) == 1L && is.finite(base_seed)) {
    tryCatch(list(
      summary = jss_scaling_summary(attempts, rule, base_seed),
      contrasts = jss_scaling_contrasts(attempts, rule$bootstrap_draws, base_seed)
    ), error = function(e) NULL)
  } else {
    NULL
  }
  if (is.null(recomputed)) {
    reasons <- c(reasons, "attempt_recomputation_failed")
  } else {
    if (!jss_scaling_tables_agree(summary, recomputed$summary, "scenario_id")) {
      reasons <- c(reasons, "summary_reconciliation_failed")
    }
    if (!jss_scaling_tables_agree(
      contrasts, recomputed$contrasts,
      c("model_class", "scenario_id", "changed_factor")
    )) {
      reasons <- c(reasons, "contrast_reconciliation_failed")
    }
  }
  actual_summary_ids <- if ("scenario_id" %in% names(summary)) sort(as.character(summary$scenario_id)) else character()
  if (!identical(actual_summary_ids, expected_ids) || anyDuplicated(actual_summary_ids)) {
    reasons <- c(reasons, "summary_scenario_incomplete")
  }
  reconciled_summary <- if (is.null(recomputed)) data.frame() else recomputed$summary
  eligible <- if ("cell_eligible" %in% names(reconciled_summary)) {
    reconciled_summary$cell_eligible %in% TRUE
  } else if (all(c("precision_met", "retained") %in% names(reconciled_summary))) {
    reconciled_summary$precision_met %in% TRUE & reconciled_summary$retained >= rule$min_attempts
  } else {
    rep(FALSE, nrow(summary))
  }
  if (length(eligible) != 8L || any(!eligible)) reasons <- c(reasons, "cell_precision_or_success_ineligible")
  if (all(c("scenario_id", "warmup") %in% names(attempts))) {
    warmups <- table(factor(attempts$scenario_id[attempts$warmup %in% TRUE], levels = expected_ids))
    if (length(warmups) != 8L || any(warmups != rule$warmup_attempts)) reasons <- c(reasons, "warmup_control_incomplete")
    if (any(attempts$retained[attempts$warmup %in% TRUE] %in% TRUE)) reasons <- c(reasons, "warmup_retained")
  } else {
    reasons <- c(reasons, "warmup_control_missing")
  }
  expected_contrasts <- if (is.null(design_error)) {
    comparison <- design[design$changed_factor != "none", , drop = FALSE]
    paste(comparison$model_class, comparison$scenario_id, comparison$changed_factor, sep = "|")
  } else {
    character()
  }
  actual_contrasts <- if (all(c("model_class", "scenario_id", "changed_factor") %in% names(contrasts))) {
    paste(contrasts$model_class, contrasts$scenario_id, contrasts$changed_factor, sep = "|")
  } else {
    character()
  }
  if (!setequal(actual_contrasts, expected_contrasts) || length(actual_contrasts) != 6L || anyDuplicated(actual_contrasts)) {
    reasons <- c(reasons, "contrast_cells_incomplete")
  }
  contrast_metrics <- c(
    "median_difference_sec", "difference_mcse_sec", "difference_ci_lower_sec",
    "difference_ci_upper_sec", "median_runtime_ratio", "ratio_mcse",
    "ratio_ci_lower", "ratio_ci_upper"
  )
  if (!all(contrast_metrics %in% names(contrasts)) ||
      nrow(contrasts) != 6L ||
      any(!vapply(contrasts[intersect(contrast_metrics, names(contrasts))], function(x) all(is.finite(x)), logical(1)))) {
    reasons <- c(reasons, "contrast_uncertainty_incomplete")
  }
  reasons <- unique(reasons)
  data.frame(
    publication_eligible = !length(reasons),
    reasons = if (length(reasons)) paste(reasons, collapse = ";") else "none",
    stringsAsFactors = FALSE
  )
}

jss_scaling_assert_publication_eligible <- function(...) {
  status <- jss_scaling_publication_status(...)
  if (!isTRUE(status$publication_eligible[[1L]])) {
    stop("Fit-scaling evidence is not publication eligible: ", status$reasons[[1L]], call. = FALSE)
  }
  invisible(status)
}

jss_scaling_cpu_model <- function() {
  candidate <- Sys.getenv("PROCESSOR_IDENTIFIER", unset = "")
  if (nzchar(candidate)) return(candidate)
  if (file.exists("/proc/cpuinfo")) {
    lines <- readLines("/proc/cpuinfo", warn = FALSE)
    model <- sub("^[^:]+:[[:space:]]*", "", lines[grepl("^model name[[:space:]]*:", lines)])
    if (length(model)) return(model[[1L]])
  }
  unname(Sys.info()[["machine"]])
}

jss_scaling_provenance_status <- function(hardware, attempts = NULL) {
  reasons <- character()
  if (!all(c("git_dirty", "execution_context") %in% names(hardware)) || nrow(hardware) != 1L) {
    reasons <- c(reasons, "provenance_fields_missing")
  } else {
    dirty <- tolower(as.character(hardware$git_dirty[[1L]]))
    if (!dirty %in% c("false", "clean", "0")) reasons <- c(reasons, "worktree_not_clean")
    if (!identical(as.character(hardware$execution_context[[1L]]), "targets")) {
      reasons <- c(reasons, "not_target_integrated")
    }
  }
  if (is.null(attempts)) {
    reasons <- c(reasons, "attempt_provenance_missing")
  } else {
    if (!jss_scaling_attempt_provenance_agrees(attempts, hardware)) {
      reasons <- c(reasons, "attempt_hardware_provenance_mismatch")
    }
    if (!all(c("timing_git_dirty", "timing_execution_context") %in% names(attempts))) {
      reasons <- c(reasons, "attempt_provenance_fields_missing")
    } else {
      attempt_dirty <- unique(tolower(as.character(attempts$timing_git_dirty)))
      attempt_context <- unique(as.character(attempts$timing_execution_context))
      if (length(attempt_dirty) != 1L || !attempt_dirty %in% c("false", "clean", "0")) {
        reasons <- c(reasons, "attempt_worktree_not_clean")
      }
      if (length(attempt_context) != 1L || !identical(attempt_context, "targets")) {
        reasons <- c(reasons, "attempt_not_target_integrated")
      }
    }
  }
  data.frame(
    production_eligible = !length(reasons),
    reasons = if (length(reasons)) paste(reasons, collapse = ";") else "none",
    stringsAsFactors = FALSE
  )
}

jss_scaling_assert_production_provenance <- function(hardware, attempts) {
  status <- jss_scaling_provenance_status(hardware, attempts)
  if (!isTRUE(status$production_eligible[[1L]])) {
    stop("Fit-scaling evidence lacks production provenance: ", status$reasons[[1L]], call. = FALSE)
  }
  invisible(status)
}

jss_scaling_hardware <- function(settings, attempts) {
  info <- Sys.info()
  versions <- jss_scaling_package_versions()
  execution_context <- if (!is.null(settings$target_integrated) && isTRUE(settings$target_integrated)) "targets" else "standalone_or_unverified"
  environment <- jss_scaling_timing_environment(settings$root, execution_context)
  data.frame(
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    git_sha = environment$timing_git_sha,
    git_dirty = environment$timing_git_dirty,
    r_version = environment$timing_r_version,
    platform = environment$timing_platform,
    os = environment$timing_os,
    cpu_model = environment$timing_cpu_model,
    logical_cores = parallel::detectCores(logical = TRUE),
    workers = settings$workers,
    package_version = as.character(utils::packageVersion("gamlss.longitudinal")),
    dependency_versions = environment$timing_dependency_versions,
    timing_environment_fingerprint = jss_scaling_hash(environment),
    checkpoint_schema_version = jss_scaling_schema_version(),
    seed_scheme_version = jss_scaling_seed_scheme_version(),
    execution_context = environment$timing_execution_context,
    timing_scope = "fit_only",
    timing_unit = "fit-only elapsed seconds",
    attempted_fits = nrow(attempts),
    stringsAsFactors = FALSE
  )
}

jss_run_10_fit_scaling <- function(settings) {
  module_id <- jss_scaling_module_id()
  tracked <- file.path(settings$public_data_dir, "fit-scaling")
  output <- list(
    scenarios = file.path(settings$data_dir, paste0(module_id, "-scenarios.csv")),
    attempts = file.path(settings$data_dir, paste0(module_id, "-attempts.csv")),
    summary = file.path(settings$tables_dir, paste0(module_id, "-summary.csv")),
    contrasts = file.path(settings$tables_dir, paste0(module_id, "-contrasts.csv")),
    hardware = file.path(settings$tables_dir, paste0(module_id, "-hardware.csv"))
  )
  if (identical(settings$profile, "paper")) {
    source <- file.path(tracked, c("scenarios.csv", "attempts.csv", "summary.csv", "contrasts.csv", "hardware.csv"))
    if (any(!file.exists(source))) stop("Tracked fit-scaling evidence is incomplete. Run the full profile and approve its public-derived bundle.", call. = FALSE)
    tracked_design <- utils::read.csv(source[[1L]], stringsAsFactors = FALSE)
    tracked_attempts <- utils::read.csv(source[[2L]], stringsAsFactors = FALSE)
    tracked_summary <- utils::read.csv(source[[3L]], stringsAsFactors = FALSE)
    tracked_contrasts <- utils::read.csv(source[[4L]], stringsAsFactors = FALSE)
    tracked_hardware <- utils::read.csv(source[[5L]], stringsAsFactors = FALSE)
    jss_scaling_assert_publication_eligible(
      tracked_design, tracked_attempts, tracked_summary, tracked_contrasts,
      jss_scaling_precision_rule("full"), hardware = tracked_hardware
    )
    jss_scaling_assert_production_provenance(tracked_hardware, tracked_attempts)
    invisible(lapply(dirname(unlist(output)), dir.create, recursive = TRUE, showWarnings = FALSE))
    copied <- file.copy(source, unlist(output), overwrite = TRUE)
    if (any(!copied)) stop("Could not stage tracked fit-scaling evidence.", call. = FALSE)
  } else {
    design <- jss_scaling_design(settings$profile)
    attempts <- jss_scaling_run(
      design,
      checkpoint_dir = file.path(settings$out_dir, "checkpoints", module_id),
      profile = settings$profile,
      base_seed = settings$seed,
      root = settings$root,
      execution_context = if (!is.null(settings$target_integrated) && isTRUE(settings$target_integrated)) "targets" else "standalone_or_unverified"
    )
    summary <- jss_scaling_summary(attempts, jss_scaling_precision_rule(settings$profile), settings$seed)
    contrasts <- jss_scaling_contrasts(attempts, jss_scaling_precision_rule(settings$profile)$bootstrap_draws, settings$seed)
    hardware <- jss_scaling_hardware(settings, attempts)
    if (identical(settings$profile, "full")) {
      jss_scaling_assert_publication_eligible(
        design, attempts, summary, contrasts, jss_scaling_precision_rule("full"),
        hardware = hardware
      )
    }
    jss_scaling_atomic_write_csv(design, output$scenarios)
    jss_scaling_atomic_write_csv(attempts, output$attempts)
    jss_scaling_atomic_write_csv(summary, output$summary)
    jss_scaling_atomic_write_csv(contrasts, output$contrasts)
    jss_scaling_atomic_write_csv(hardware, output$hardware)
  }
  output_hardware <- utils::read.csv(output$hardware, stringsAsFactors = FALSE)
  output_attempts <- utils::read.csv(output$attempts, stringsAsFactors = FALSE)
  provenance <- jss_scaling_provenance_status(output_hardware, output_attempts)
  list(
    module_id = module_id,
    status = if (isTRUE(provenance$production_eligible[[1L]])) "current" else "candidate_not_production_approved",
    data = c(output$scenarios, output$attempts),
    tables = c(output$summary, output$contrasts, output$hardware),
    figures = character(),
    notes = "Measured one-factor-at-a-time fit scaling for representative continuous and discrete models."
  )
}
