#!/usr/bin/env Rscript

# Backward formula screening for the full ACTG 193A dataset.
# This deliberately remains outside the manuscript targets graph.
#
# At each step, all valid one-term deletions are fitted. The deletion with the
# lowest information criterion is accepted if it improves on the current
# model. Independent candidates within a step can run in parallel, and every
# candidate is checkpointed so interrupted runs can resume.
#
# Full run (defaults shown):
# Rscript paper/reviews/2026-09-05-medical-candidates/run-screening-v2.R \
#   --criterion=AIC --hierarchy=strong --workers=4 --parameters=mu
# Smooth-aware staged run:
# Rscript paper/reviews/2026-09-05-medical-candidates/run-screening-v2.R \
#   --stage=mu-smooth --smooth-k=8 --workers=4
# Quick configuration check without fitting:
# Rscript paper/reviews/2026-09-05-medical-candidates/run-screening-v2.R \
#   --dry-run=true

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = FALSE)
} else {
  file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "run-screening-v2.R"
  )
}
package_root <- normalizePath(
  file.path(dirname(script_path), "..", "..", ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(package_root, "DESCRIPTION"))) {
  stop("Could not locate the gamlss.longitudinal package root from this script.")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools with install.packages('devtools').")
}
devtools::load_all(package_root, quiet = TRUE)

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

argument_flag <- function(name, default) {
  value <- tolower(argument_value(name, if (default) "true" else "false"))
  if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
    stop("--", name, " must be true or false.")
  }
  value %in% c("true", "1", "yes")
}

positive_integer <- function(value, name, allow_inf = FALSE) {
  if (allow_inf && identical(tolower(value), "inf")) return(Inf)
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L) {
    stop("--", name, " must be a positive integer", if (allow_inf) " or Inf" else "", ".")
  }
  parsed
}

physical_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(physical_cores)) physical_cores <- 2L
default_workers <- max(1L, min(4L, as.integer(physical_cores) - 1L))

criterion <- toupper(argument_value("criterion", "AIC"))
if (!criterion %in% c("AIC", "BIC")) stop("--criterion must be AIC or BIC.")
hierarchy <- tolower(argument_value("hierarchy", "strong"))
if (!hierarchy %in% c("strong", "none")) stop("--hierarchy must be strong or none.")
workers <- positive_integer(
  argument_value("workers", Sys.getenv("N_CORES", unset = as.character(default_workers))),
  "workers"
)
max_steps <- positive_integer(argument_value("max-steps", "Inf"), "max-steps", TRUE)
initial_inner_iter <- positive_integer(
  argument_value("initial-inner-iter", "2"), "initial-inner-iter"
)
max_elapsed_sec <- suppressWarnings(as.numeric(argument_value("max-elapsed-sec", "1190")))
if (length(max_elapsed_sec) != 1L || is.na(max_elapsed_sec) || max_elapsed_sec <= 0) {
  stop("--max-elapsed-sec must be a positive number.")
}
pilot_low_lambda_max <- suppressWarnings(as.numeric(
  argument_value("pilot-low-lambda-max", "1")
))
if (length(pilot_low_lambda_max) != 1L || is.na(pilot_low_lambda_max) ||
    pilot_low_lambda_max < 0) {
  stop("--pilot-low-lambda-max must be a non-negative number.")
}
pilot_smooth_projection_rmse_min <- suppressWarnings(as.numeric(
  argument_value("pilot-smooth-projection-rmse-min", "0.25")
))
if (length(pilot_smooth_projection_rmse_min) != 1L ||
    is.na(pilot_smooth_projection_rmse_min) ||
    pilot_smooth_projection_rmse_min < 0) {
  stop("--pilot-smooth-projection-rmse-min must be a non-negative number.")
}
pilot_delta_threshold <- suppressWarnings(as.numeric(
  argument_value("pilot-delta-threshold", "100")
))
if (length(pilot_delta_threshold) != 1L || is.na(pilot_delta_threshold) ||
    pilot_delta_threshold < 0) {
  stop("--pilot-delta-threshold must be a non-negative number.")
}
pilot_fixed_projection_rmse_min <- suppressWarnings(as.numeric(
  argument_value("pilot-fixed-projection-rmse-min", "5")
))
if (length(pilot_fixed_projection_rmse_min) != 1L ||
    is.na(pilot_fixed_projection_rmse_min) ||
    pilot_fixed_projection_rmse_min < 0) {
  stop("--pilot-fixed-projection-rmse-min must be a non-negative number.")
}
pilot_fixed_parent_terms <- positive_integer(
  argument_value("pilot-fixed-parent-terms", "3"),
  "pilot-fixed-parent-terms"
)
pilot_outer_iter <- positive_integer(
  argument_value("pilot-outer-iter", "2"), "pilot-outer-iter"
)
pilot_bfgs_max_iter <- positive_integer(
  argument_value("pilot-bfgs-max-iter", "20"), "pilot-bfgs-max-iter"
)
low_lambda_initial_inner_iter <- positive_integer(
  argument_value("low-lambda-initial-inner-iter", "5"),
  "low-lambda-initial-inner-iter"
)
resume <- argument_flag("resume", TRUE)
dry_run <- argument_flag("dry-run", FALSE)
subject_argument <- tolower(argument_value("subjects", "all"))
stage <- tolower(argument_value("stage", "mu-linear"))
margin_family <- toupper(argument_value("margin", "ZIPIG"))
if (!margin_family %in% c("ZIPIG", "ZISICHEL", "ZIBNB")) {
  stop("--margin must be ZIPIG, ZISICHEL, or ZIBNB.")
}
copula_argument <- tolower(argument_value("copula", "G"))
copula_code <- switch(copula_argument, g = "G", gumbel = "G", t = "t", student = "t",
  `student-t` = "t", stop("--copula must be G/Gumbel or t/Student-t."))
copula_label <- if (identical(copula_code, "G")) "Gumbel (G)" else "Student-t (t)"
has_tau <- margin_family %in% c("ZISICHEL", "ZIBNB")
has_zeta <- identical(copula_code, "t")
valid_stages <- c("mu-linear", "mu-smooth", "sigma", "nu", "tau", "theta", "zeta", "refine")
if (!stage %in% valid_stages) {
  stop("--stage must be one of: ", paste(valid_stages, collapse = ", "), ".")
}
if (identical(stage, "tau") && !has_tau) {
  stop("--stage=tau is unavailable for margin ", margin_family, ".")
}
if (identical(stage, "zeta") && !has_zeta) {
  stop("--stage=zeta requires --copula=t.")
}
smooth_k <- positive_integer(argument_value("smooth-k", "8"), "smooth-k")
parameters <- trimws(strsplit(
  tolower(argument_value("parameters", "mu")), ",", fixed = TRUE
)[[1L]])
valid_parameters <- c("mu", "sigma", "nu", if (has_tau) "tau", "theta", if (has_zeta) "zeta")
if (!length(parameters) || any(!parameters %in% valid_parameters)) {
  stop("--parameters must contain one or more of: ", paste(valid_parameters, collapse = ","), ".")
}
if (identical(stage, "refine") && length(parameters) != 1L) {
  stop("--stage=refine requires exactly one value in --parameters.")
}
input_formulas_argument <- argument_value("input-formulas", "")
initial_fit_argument <- argument_value("initial-fit", "")
output_stem_argument <- argument_value("output-stem", "")
results_subdir <- argument_value("results-subdir", "")
if (nzchar(output_stem_argument) && !grepl("^[A-Za-z0-9._-]+$", output_stem_argument)) {
  stop("--output-stem may contain only letters, numbers, periods, underscores, and hyphens.")
}
if (nzchar(results_subdir) && !grepl("^[A-Za-z0-9._-]+$", results_subdir)) {
  stop("--results-subdir may contain only letters, numbers, periods, underscores, and hyphens.")
}

output_dir <- file.path(
  package_root, "paper", "reviews", "2026-09-05-medical-candidates", "results",
  results_subdir
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load_package_data <- function(name) {
  data_environment <- new.env(parent = emptyenv())
  utils::data(
    list = name, package = "gamlss.longitudinal", envir = data_environment
  )
  get(name, envir = data_environment, inherits = FALSE)
}

read_actg <- function() {
  data <- load_package_data("actg_cd4")
  data$week_factor <- factor(data$week)
  data$baseline_log_z <- as.numeric(scale(log1p(data$baseline_cd4)))
  data$age_z <- as.numeric(scale(data$age))
  data
}

sample_patients <- function(data, id_var, n, seed = 1L) {
  subject_rows <- data[!duplicated(data[[id_var]]), c(id_var, "treatment")]
  groups <- split(subject_rows[[id_var]], subject_rows$treatment)
  set.seed(seed)
  allocation <- rep(n %/% length(groups), length(groups))
  allocation[seq_len(n %% length(groups))] <-
    allocation[seq_len(n %% length(groups))] + 1L
  ids <- unlist(Map(
    function(group_ids, size) sample(group_ids, min(size, length(group_ids))),
    groups,
    allocation
  ), use.names = FALSE)
  data[data[[id_var]] %in% ids, , drop = FALSE]
}

actg_data <- read_actg()
if (!identical(subject_argument, "all")) {
  n_subjects <- positive_integer(subject_argument, "subjects")
  total_subjects <- length(unique(actg_data$subject))
  if (n_subjects > total_subjects) stop("--subjects cannot exceed ", total_subjects, ".")
  actg_data <- sample_patients(actg_data, "subject", n_subjects, seed = 20260909L)
}
stage_label <- if (identical(stage, "mu-linear")) "" else paste0("-", stage)
family_stem <- if (identical(margin_family, "ZIPIG")) "actg" else paste0("actg-", tolower(margin_family))
if (has_zeta) family_stem <- paste0(family_stem, "-t")
default_run_stem <- if (identical(subject_argument, "all")) {
  paste0(family_stem, "-backward-selection", stage_label)
} else {
  paste0(
    family_stem, "-backward-selection", stage_label, "-n",
    length(unique(actg_data$subject))
  )
}
run_stem <- if (nzchar(output_stem_argument)) output_stem_argument else default_run_stem
checkpoint_dir <- file.path(output_dir, paste0(run_stem, "-checkpoints"))
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

# Do not reuse fits made against a different version of the package source.
package_source_files <- c(
  file.path(package_root, "DESCRIPTION"),
  list.files(file.path(package_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
)
package_source_fingerprint <- tools::md5sum(package_source_files)

linear_mu_maximum <- cd4 ~ treatment + age + sex + baseline_cd4 +
  week_factor + treatment:week_factor
smooth_mu_maximum <- stats::as.formula(sprintf(
  paste0(
    "cd4 ~ treatment + s(age, bs = 'ps', k = %d) + sex + ",
    "s(baseline_cd4, bs = 'ps', k = %d) + week_factor + ",
    "treatment:week_factor"
  ),
  smooth_k,
  smooth_k
))
smooth_parameter_maximum <- stats::as.formula(sprintf(
  paste0(
    "~ treatment + s(age, bs = 'ps', k = %d) + sex + ",
    "s(baseline_cd4, bs = 'ps', k = %d) + week_factor + ",
    "treatment:week_factor"
  ),
  smooth_k,
  smooth_k
))

stage_parameter <- switch(
  stage,
  `mu-linear` = "mu",
  `mu-smooth` = "mu",
  sigma = "sigma",
  nu = "nu",
  tau = "tau",
  theta = "theta",
  zeta = "zeta",
  refine = NULL
)
if (!stage %in% c("mu-linear", "refine")) parameters <- stage_parameter

prior_stage <- switch(
  stage,
  sigma = "mu-smooth",
  nu = "sigma",
  tau = "nu",
  theta = if (has_tau) "tau" else "nu",
  zeta = "theta",
  NULL
)
if (identical(stage, "refine")) {
  if (!nzchar(input_formulas_argument) || !file.exists(input_formulas_argument)) {
    stop("--stage=refine requires an existing --input-formulas RDS file.")
  }
  maximum_formulas <- readRDS(input_formulas_argument)
  required_formulas <- c("mu", "sigma", "nu", if (has_tau) "tau", "theta", if (has_zeta) "zeta")
  if (!is.list(maximum_formulas) || !all(required_formulas %in% names(maximum_formulas))) {
    stop("The refinement input must contain mu, sigma, nu, and theta formulas.")
  }
  active_parameter <- parameters[[1L]]
  maximum_formulas[[active_parameter]] <- if (identical(active_parameter, "mu")) {
    smooth_mu_maximum
  } else {
    smooth_parameter_maximum
  }
} else if (is.null(prior_stage)) {
  maximum_formulas <- list(
    mu = if (identical(stage, "mu-smooth")) smooth_mu_maximum else linear_mu_maximum,
    sigma = ~ 1,
    nu = ~ 1,
    tau = if (has_tau) ~ 1 else NULL,
    theta = ~ 1,
    zeta = if (has_zeta) ~ 1 else NULL
  )
  maximum_formulas <- maximum_formulas[!vapply(maximum_formulas, is.null, logical(1L))]
} else {
  prior_stem <- if (identical(subject_argument, "all")) {
    paste0(family_stem, "-backward-selection-", prior_stage)
  } else {
    paste0(
      family_stem, "-backward-selection-", prior_stage, "-n",
      length(unique(actg_data$subject))
    )
  }
  prior_formula_path <- file.path(
    output_dir, paste0(prior_stem, "-final-formulas.rds")
  )
  if (!file.exists(prior_formula_path)) {
    stop(
      "Stage '", stage, "' requires the completed '", prior_stage,
      "' formulas at: ", prior_formula_path
    )
  }
  maximum_formulas <- readRDS(prior_formula_path)
  maximum_formulas[[stage_parameter]] <- smooth_parameter_maximum
}

optimizer_control <- gamlss_longitudinal_control(
  outer_tol = 1e-4,
  max_outer_iter = 50L,
  max_elapsed_sec = max_elapsed_sec,
  rs = list(
    inner_tol = 1e-5,
    initial_inner_iter = initial_inner_iter
  )
)
optimizer_control <- gamlss.longitudinal:::.gl_apply_route_optimizer_defaults(
  control = optimizer_control,
  method = "RS",
  margin_dist = getExportedValue("gamlss.dist", margin_family)(),
  copula_dist = copula_code,
  optimization = "joint",
  purpose = "screening"
)

formula_text <- function(formula) {
  paste(deparse(formula, width.cutoff = 500L), collapse = " ")
}

formula_signature <- function(formulas) {
  vapply(formulas, formula_text, character(1L))
}

object_md5 <- function(object) {
  path <- tempfile("actg-screen-signature-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 2)
  unname(tools::md5sum(path))
}

checkpoint_metadata_without_elapsed_cap <- function(metadata) {
  value <- metadata
  if (is.list(value$optimizer_control) && is.list(value$optimizer_control$shared)) {
    value$optimizer_control$shared$max_elapsed_sec <- NULL
  }
  value
}

checkpoint_metadata_compatible <- function(stored, requested) {
  if (identical(stored, requested)) return(TRUE)
  if (!is.list(stored) || !is.list(requested)) return(FALSE)
  stored_cap <- stored$optimizer_control$shared$max_elapsed_sec
  requested_cap <- requested$optimizer_control$shared$max_elapsed_sec
  length(stored_cap) == 1L && is.finite(stored_cap) &&
    length(requested_cap) == 1L && is.finite(requested_cap) &&
    as.numeric(stored_cap) <= as.numeric(requested_cap) &&
    identical(
      checkpoint_metadata_without_elapsed_cap(stored),
      checkpoint_metadata_without_elapsed_cap(requested)
    )
}

safe_name <- function(value) {
  value <- gsub("[^A-Za-z0-9]+", "-", value)
  tolower(gsub("(^-+|-+$)", "", value))
}

remove_formula_term <- function(formula, term) {
  formula_terms <- stats::terms(formula)
  labels <- attr(formula_terms, "term.labels")
  response <- if (attr(formula_terms, "response") == 1L) {
    deparse(formula[[2L]], width.cutoff = 500L)
  } else {
    NULL
  }
  stats::reformulate(
    labels[labels != term],
    response = response,
    intercept = attr(formula_terms, "intercept") == 1L,
    env = environment(formula)
  )
}

replace_formula_term <- function(formula, term, replacement) {
  formula_terms <- stats::terms(formula)
  labels <- attr(formula_terms, "term.labels")
  labels[labels == term] <- replacement
  labels <- unique(labels)
  response <- if (attr(formula_terms, "response") == 1L) {
    deparse(formula[[2L]], width.cutoff = 500L)
  } else {
    NULL
  }
  stats::reformulate(
    labels,
    response = response,
    intercept = attr(formula_terms, "intercept") == 1L,
    env = environment(formula)
  )
}

droppable_terms <- function(formula, hierarchy) {
  labels <- attr(stats::terms(formula), "term.labels")
  if (!length(labels) || identical(hierarchy, "none")) return(labels)
  interactions <- labels[grepl(":", labels, fixed = TRUE)]
  protected <- unique(unlist(strsplit(interactions, ":", fixed = TRUE), use.names = FALSE))
  setdiff(labels, protected)
}

smooth_linear_term <- function(term) {
  if (!grepl("^\\s*s\\(", term)) return(NULL)
  call <- tryCatch(parse(text = term)[[1L]], error = function(error) NULL)
  if (is.null(call) || length(call) < 2L) return(NULL)
  paste(deparse(call[[2L]], width.cutoff = 500L), collapse = " ")
}

formula_moves <- function(formula, hierarchy) {
  drops <- droppable_terms(formula, hierarchy)
  moves <- lapply(drops, function(term) {
    list(label = term, formula = remove_formula_term(formula, term), type = "drop")
  })
  labels <- attr(stats::terms(formula), "term.labels")
  smooth_terms <- labels[grepl("^\\s*s\\(", labels)]
  for (term in smooth_terms) {
    replacement <- smooth_linear_term(term)
    if (is.null(replacement) || replacement %in% labels) next
    moves[[length(moves) + 1L]] <- list(
      label = paste0(term, " -> ", replacement),
      formula = replace_formula_term(formula, term, replacement),
      type = "simplify"
    )
  }
  moves
}

start_for_formulas <- function(fit, formulas, data) {
  if (!inherits(fit, "gamlss.longitudinal")) return(NULL)
  # The package maps retained fixed coefficients, compatible smooth basis
  # coefficients, and smoothing penalties by name. Passing the fitted parent
  # also bypasses the generic five-outer separate warm phase.
  fit
}

screening_start_strategy <- function(
    copula,
    step,
    screening_stage,
    parameter,
    dropped_term,
    change_type,
    project_smooth_drop = FALSE) {
  # In the ACTG Student-t routes, inheriting the fitted zeta surface after
  # deleting a zeta interaction reproducibly enters a poorer local optimum.
  # Reset only that screened parameter; every other parameter, plus the
  # fitted smoothing penalties, remains warm-started from the parent.
  reset_zeta_interaction <- identical(copula, "t") &&
    identical(parameter, "zeta") &&
    identical(change_type, "drop") &&
    length(dropped_term) == 1L && !is.na(dropped_term) &&
    grepl(":", dropped_term, fixed = TRUE)
  # Once the ACTG zeta model has been reduced to treatment + sex + baseline,
  # deleting treatment from the fitted-parent state can wander for more than
  # an hour and fail in another parameter block. The same candidate reaches a
  # converged, clearly noncompetitive endpoint in under 20 minutes from a
  # zero-zeta start. Keep this exception terminal and term-specific.
  reset_terminal_zeta_treatment <- identical(copula, "t") &&
    as.integer(step) >= 5L && identical(parameter, "zeta") &&
    identical(change_type, "drop") && identical(dropped_term, "treatment")
  project_refinement_smooth_drop <- isTRUE(project_smooth_drop) &&
    identical(screening_stage, "refine") &&
    identical(change_type, "drop") && length(dropped_term) == 1L &&
    !is.na(dropped_term) && grepl("^\\s*s\\(", dropped_term)
  if (reset_zeta_interaction || reset_terminal_zeta_treatment) {
    "zero_screened_parameter"
  } else if (project_refinement_smooth_drop) {
    "project_parent_predictor"
  } else {
    "fitted_parent"
  }
}

zero_screened_parameter_start <- function(start, parameter) {
  if (!inherits(start, "gamlss.longitudinal")) return(start)
  fixed_names <- startsWith(names(start$par), paste0(parameter, "."))
  start$par[fixed_names] <- 0
  if (parameter %in% names(start$par_s)) {
    start$par_s[[parameter]] <- lapply(
      start$par_s[[parameter]],
      function(value) value * 0
    )
  }
  start
}

project_smooth_drop_start <- function(start, parameter, dropped_term) {
  if (!inherits(start, "gamlss.longitudinal") ||
      is.na(dropped_term) || !nzchar(dropped_term)) {
    return(start)
  }
  term <- trimws(sub("\\s*->.*$", "", dropped_term))
  fixed_design <- start$model_matrix$x[[parameter]]
  fixed_names <- paste(parameter, colnames(fixed_design), sep = ".")
  smooth_design <- start$model_matrix$s[[parameter]]
  smooth_coefficients <- start$par_s[[parameter]]
  if (is.null(fixed_design) || !all(fixed_names %in% names(start$par)) ||
      is.null(smooth_design) || is.null(smooth_coefficients) ||
      !term %in% names(smooth_coefficients)) {
    return(start)
  }
  fixed_design <- as.matrix(fixed_design)
  eta <- as.numeric(fixed_design %*% as.numeric(start$par[fixed_names]))
  smooth_names <- intersect(names(smooth_design), names(smooth_coefficients))
  for (name in smooth_names) {
    basis <- as.matrix(smooth_design[[name]])
    coefficient <- as.numeric(smooth_coefficients[[name]])
    if (nrow(basis) != length(eta) || ncol(basis) != length(coefficient)) {
      return(start)
    }
    eta <- eta + as.numeric(basis %*% coefficient)
  }
  retained_smooths <- setdiff(smooth_names, term)
  retained_designs <- lapply(retained_smooths, function(name) {
    as.matrix(smooth_design[[name]])
  })
  candidate_design <- do.call(cbind, c(list(fixed_design), retained_designs))
  if (!ncol(candidate_design) || any(!is.finite(candidate_design)) ||
      any(!is.finite(eta))) {
    return(start)
  }
  projected <- tryCatch(stats::lm.fit(candidate_design, eta), error = identity)
  if (inherits(projected, "error") ||
      length(projected$coefficients) != ncol(candidate_design) ||
      any(!is.finite(projected$coefficients))) {
    return(start)
  }
  offset <- ncol(fixed_design)
  start$par[fixed_names] <- projected$coefficients[seq_len(offset)]
  for (name in retained_smooths) {
    width <- ncol(as.matrix(smooth_design[[name]]))
    index <- offset + seq_len(width)
    value <- projected$coefficients[index]
    names(value) <- names(start$par_s[[parameter]][[name]])
    start$par_s[[parameter]][[name]] <- value
    offset <- offset + width
  }
  start
}

start_signature <- function(start) {
  if (!inherits(start, "gamlss.longitudinal")) return(object_md5(start))
  object_md5(list(
    par = start$par,
    par_s = start$par_s,
    lambda_s = start$lambda_s,
    margin = start$margin_dist$family,
    copula = start$copula_dist,
    formulas = formula_signature(start$formulas)
  ))
}

empty_metrics <- function() {
  c(logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, EDF = NA_real_)
}

candidate_checkpoint_paths <- function(task) {
  prior <- if (length(task$checkpoint_pattern) == 1L) {
    list.files(
      dirname(task$checkpoint),
      pattern = task$checkpoint_pattern,
      full.names = TRUE
    )
  } else {
    character()
  }
  unique(c(task$checkpoint, prior))
}

fit_screening_candidate <- function(task) {
  checkpoint_paths <- if (isTRUE(resume)) {
    candidate_checkpoint_paths(task)
  } else {
    task$checkpoint
  }
  for (checkpoint_path in checkpoint_paths[file.exists(checkpoint_paths)]) {
    cached <- tryCatch(readRDS(checkpoint_path), error = identity)
    cached_has_metrics <- is.list(cached) && is.numeric(cached$metrics) &&
      all(names(empty_metrics()) %in% names(cached$metrics))
    stored_cap <- if (is.list(cached) && is.list(cached$metadata)) {
      cached$metadata$optimizer_control$shared$max_elapsed_sec
    } else {
      NA_real_
    }
    requested_cap <- task$metadata$optimizer_control$shared$max_elapsed_sec
    looser_elapsed_cap <- is.list(cached) && length(stored_cap) == 1L &&
      is.finite(stored_cap) && length(requested_cap) == 1L &&
      is.finite(requested_cap) && requested_cap > stored_cap
    elapsed_timeout <- is.list(cached) && length(cached$error) == 1L &&
      !is.na(cached$error) && grepl("exceeded max_elapsed_sec", cached$error, fixed = TRUE)
    cache_is_reusable <- is.list(cached) &&
      checkpoint_metadata_compatible(cached$metadata, task$metadata) &&
      inherits(cached$fit, "gamlss.longitudinal") && cached_has_metrics &&
      (isTRUE(cached$converged) || identical(cached$converged, FALSE)) &&
      !(looser_elapsed_cap && elapsed_timeout)
    if (cache_is_reusable) {
      cached$from_checkpoint <- TRUE
      cached$checkpoint <- checkpoint_path
      # The complete fit is already durable in the checkpoint. Returning it
      # through a PSOCK socket makes the coordinator deserialize every
      # candidate fit in a step, even though it only needs scalar metrics to
      # rank candidates. Reload only the accepted checkpoint below.
      cached$fit <- NULL
      return(cached)
    }
  }

  optional_formulas <- c(
    if (has_tau) list(tau.formula = task$formulas$tau) else list(),
    if (has_zeta) list(zeta.formula = task$formulas$zeta) else list()
  )
  run_fit <- function(control) {
    started <- proc.time()[["elapsed"]]
    captured_warnings <- character()
    fit <- tryCatch(
      withCallingHandlers(
        do.call(gamlss_longitudinal, c(list(
          data = actg_data,
          time_var = "week_factor",
          subject_var = "subject",
          margin_dist = getExportedValue("gamlss.dist", margin_family)(),
          copula_dist = copula_code,
          mu.formula = task$formulas$mu,
          sigma.formula = task$formulas$sigma,
          nu.formula = task$formulas$nu,
          theta.formula = task$formulas$theta,
          missingness = "segment",
          method = "RS",
          optimization = "joint",
          optimizer_control = control,
          start = task$start,
          compute_vcov = FALSE,
          verbose = 0
        ), optional_formulas)),
        warning = function(warning) {
          captured_warnings <<- c(captured_warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      ),
      error = identity
    )
    list(
      fit = fit,
      elapsed_sec = proc.time()[["elapsed"]] - started,
      warnings = unique(captured_warnings)
    )
  }

  result_from_attempt <- function(attempt, elapsed_sec = attempt$elapsed_sec) {
    fit <- attempt$fit
    returned_fit <- inherits(fit, "gamlss.longitudinal")
    metrics <- empty_metrics()
    if (returned_fit && !is.null(fit$model_selection)) {
      metrics <- c(
        logLik = as.numeric(fit$model_selection["LogLik", "joint"]),
        AIC = as.numeric(fit$model_selection["AIC", "joint"]),
        BIC = as.numeric(fit$model_selection["BIC", "joint"]),
        EDF = as.numeric(fit$model_selection["EDF", "joint"])
      )
    }
    list(
      metadata = task$metadata,
      fit = if (returned_fit) fit else NULL,
      metrics = metrics,
      converged = returned_fit && isTRUE(fit$convergence$converged),
      elapsed_sec = elapsed_sec,
      outer_iterations = if (returned_fit) fit$convergence$outer_iterations else NA_integer_,
      stop_reason = if (returned_fit) fit$convergence$stop_reason else "error",
      warnings = attempt$warnings,
      error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
      checkpoint = task$checkpoint,
      from_checkpoint = FALSE
    )
  }

  candidate_control <- optimizer_control
  inner_override <- task$metadata$rs_initial_inner_iter_override
  if (length(inner_override) == 1L && is.finite(inner_override)) {
    candidate_control$rs$initial_inner_iter <- as.integer(inner_override)
  }

  pilot <- task$metadata$dominance_pilot
  if (is.list(pilot)) {
    pilot_control <- candidate_control
    pilot_control$shared$max_outer_iter <- pilot$max_outer_iter
    pilot_control$rs$bfgs_max_iter <- pilot$bfgs_max_iter
    pilot_attempt <- run_fit(pilot_control)
    pilot_result <- result_from_attempt(pilot_attempt)
    pilot_delta <- pilot_result$metrics[[pilot$criterion]] - pilot$parent_value
    dominated <- !isTRUE(pilot_result$converged) &&
      identical(pilot_result$stop_reason, "max_iterations") &&
      is.finite(pilot_delta) && pilot_delta > pilot$delta_threshold
    if (dominated) {
      pilot_result$dominance_pilot_rejected <- TRUE
      result <- pilot_result
    } else {
      full_control <- candidate_control
      remaining <- full_control$shared$max_elapsed_sec - pilot_attempt$elapsed_sec
      if (!is.finite(remaining) || remaining <= 0) {
        result <- pilot_result
      } else {
        full_control$shared$max_elapsed_sec <- remaining
        full_attempt <- run_fit(full_control)
        result <- result_from_attempt(
          full_attempt,
          elapsed_sec = pilot_attempt$elapsed_sec + full_attempt$elapsed_sec
        )
        result$dominance_pilot_rejected <- FALSE
      }
    }
  } else {
    result <- result_from_attempt(run_fit(candidate_control))
  }
  temporary <- tempfile("candidate-", tmpdir = dirname(task$checkpoint), fileext = ".rds")
  saveRDS(result, temporary, compress = "gzip")
  saved <- if (file.exists(task$checkpoint)) {
    file.copy(temporary, task$checkpoint, overwrite = TRUE)
  } else {
    file.rename(temporary, task$checkpoint)
  }
  if (file.exists(temporary)) unlink(temporary)
  if (!saved) {
    stop("Could not save checkpoint: ", task$checkpoint)
  }
  result$fit <- NULL
  result
}

load_checkpoint_fit <- function(result, task) {
  checkpoint <- result$checkpoint %||% task$checkpoint
  complete <- tryCatch(readRDS(checkpoint), error = identity)
  valid <- is.list(complete) &&
    checkpoint_metadata_compatible(complete$metadata, task$metadata) &&
    isTRUE(complete$converged) &&
    inherits(complete$fit, "gamlss.longitudinal")
  if (!valid) {
    detail <- if (inherits(complete, "error")) conditionMessage(complete) else {
      "checkpoint did not contain the expected converged fitted object"
    }
    stop("Could not reload accepted checkpoint ", checkpoint, ": ", detail)
  }
  complete$from_checkpoint <- isTRUE(result$from_checkpoint)
  complete
}

smooth_penalty_for_term <- function(start, parameter, dropped_term) {
  if (!inherits(start, "gamlss.longitudinal") || is.null(start$lambda_s) ||
      is.na(dropped_term) || !nzchar(dropped_term)) {
    return(NA_real_)
  }
  term <- trimws(sub("\\s*->.*$", "", dropped_term))
  if (!grepl("^s\\(", term)) return(NA_real_)
  values <- unlist(start$lambda_s, recursive = TRUE, use.names = TRUE)
  if (!length(values) || is.null(names(values))) return(NA_real_)
  prefix <- paste0(parameter, ".", term)
  matched <- startsWith(names(values), prefix)
  values <- suppressWarnings(as.numeric(values[matched]))
  values <- values[is.finite(values)]
  if (length(values) == 1L) values[[1L]] else NA_real_
}

fixed_drop_projection_rmse <- function(start, parameter, dropped_term) {
  if (!inherits(start, "gamlss.longitudinal") ||
      is.na(dropped_term) || !nzchar(dropped_term)) {
    return(NA_real_)
  }
  term <- trimws(sub("\\s*->.*$", "", dropped_term))
  if (grepl("^s\\(", term)) return(NA_real_)
  formula <- start$formulas[[parameter]]
  x <- start$model_matrix$x[[parameter]]
  if (is.null(formula) || is.null(x)) return(NA_real_)
  x <- as.matrix(x)
  assignment <- attr(start$model_matrix$x[[parameter]], "assign")
  formula_terms <- attr(stats::terms(formula), "term.labels")
  term_index <- match(term, formula_terms)
  if (is.na(term_index) || length(assignment) != ncol(x) ||
      !any(assignment == term_index)) {
    return(NA_real_)
  }
  fixed_names <- paste(parameter, colnames(x), sep = ".")
  if (!all(fixed_names %in% names(start$par))) return(NA_real_)
  eta <- as.numeric(x %*% as.numeric(start$par[fixed_names]))
  smooth_design <- start$model_matrix$s[[parameter]] %||% list()
  smooth_coefficients <- start$par_s[[parameter]] %||% list()
  retained_smooths <- intersect(names(smooth_design), names(smooth_coefficients))
  for (name in retained_smooths) {
    basis <- as.matrix(smooth_design[[name]])
    coefficient <- as.numeric(smooth_coefficients[[name]])
    if (nrow(basis) != length(eta) || ncol(basis) != length(coefficient)) {
      return(NA_real_)
    }
    eta <- eta + as.numeric(basis %*% coefficient)
  }
  retained_fixed <- assignment != term_index
  candidate_design <- x[, retained_fixed, drop = FALSE]
  if (length(retained_smooths)) {
    candidate_design <- do.call(cbind, c(
      list(candidate_design),
      lapply(smooth_design[retained_smooths], as.matrix)
    ))
  }
  if (!ncol(candidate_design) || !all(is.finite(candidate_design)) ||
      !all(is.finite(eta))) {
    return(NA_real_)
  }
  projected <- tryCatch(
    stats::lm.fit(candidate_design, eta),
    error = identity
  )
  if (inherits(projected, "error") ||
      length(projected$fitted.values) != length(eta) ||
      any(!is.finite(projected$fitted.values))) {
    return(NA_real_)
  }
  sqrt(mean((eta - projected$fitted.values)^2))
}

smooth_drop_projection_rmse <- function(start, parameter, dropped_term) {
  if (!inherits(start, "gamlss.longitudinal") ||
      is.na(dropped_term) || !nzchar(dropped_term)) {
    return(NA_real_)
  }
  term <- trimws(sub("\\s*->.*$", "", dropped_term))
  if (!grepl("^s\\(", term)) return(NA_real_)
  fixed_design <- start$model_matrix$x[[parameter]]
  fixed_names <- paste(parameter, colnames(fixed_design), sep = ".")
  smooth_design <- start$model_matrix$s[[parameter]]
  smooth_coefficients <- start$par_s[[parameter]]
  if (is.null(fixed_design) || !all(fixed_names %in% names(start$par)) ||
      is.null(smooth_design) || is.null(smooth_coefficients) ||
      !term %in% names(smooth_coefficients)) {
    return(NA_real_)
  }
  fixed_design <- as.matrix(fixed_design)
  eta <- as.numeric(fixed_design %*% as.numeric(start$par[fixed_names]))
  smooth_names <- intersect(names(smooth_design), names(smooth_coefficients))
  for (name in smooth_names) {
    basis <- as.matrix(smooth_design[[name]])
    coefficient <- as.numeric(smooth_coefficients[[name]])
    if (nrow(basis) != length(eta) || ncol(basis) != length(coefficient)) {
      return(NA_real_)
    }
    eta <- eta + as.numeric(basis %*% coefficient)
  }
  retained_smooths <- setdiff(smooth_names, term)
  candidate_design <- do.call(cbind, c(
    list(fixed_design),
    lapply(retained_smooths, function(name) as.matrix(smooth_design[[name]]))
  ))
  if (!ncol(candidate_design) || any(!is.finite(candidate_design)) ||
      any(!is.finite(eta))) {
    return(NA_real_)
  }
  projected <- tryCatch(
    stats::lm.fit(candidate_design, eta),
    error = identity
  )
  if (inherits(projected, "error") ||
      length(projected$fitted.values) != length(eta) ||
      any(!is.finite(projected$fitted.values))) {
    return(NA_real_)
  }
  sqrt(mean((eta - projected$fitted.values)^2))
}

dominance_pilot_for_task <- function(
    step,
    parameter,
    dropped_term,
    change_type,
    start) {
  if (!identical(stage, "refine")) return(NULL)
  if (!change_type %in% c("drop", "simplify")) return(NULL)
  lambda <- smooth_penalty_for_term(start, parameter, dropped_term)
  low_lambda_smooth <- is.finite(lambda) && lambda <= pilot_low_lambda_max
  smooth_projection_rmse <- if (identical(change_type, "drop")) {
    smooth_drop_projection_rmse(start, parameter, dropped_term)
  } else {
    NA_real_
  }
  high_projection_smooth <- is.finite(smooth_projection_rmse) &&
    smooth_projection_rmse >= pilot_smooth_projection_rmse_min
  projection_rmse <- if (identical(change_type, "drop")) {
    fixed_drop_projection_rmse(start, parameter, dropped_term)
  } else {
    NA_real_
  }
  parent_term_count <- tryCatch(
    length(attr(stats::terms(start$formulas[[parameter]]), "term.labels")),
    error = function(error) NA_integer_
  )
  high_projection_fixed <- is.finite(projection_rmse) &&
    projection_rmse >= pilot_fixed_projection_rmse_min &&
    is.finite(parent_term_count) &&
    parent_term_count == pilot_fixed_parent_terms
  if (!low_lambda_smooth && !high_projection_smooth &&
      !high_projection_fixed) return(NULL)
  parent_value <- tryCatch(
    as.numeric(start$model_selection[criterion, "joint"]),
    error = function(error) NA_real_
  )
  if (!is.finite(parent_value)) return(NULL)
  list(
    criterion = criterion,
    parent_value = parent_value,
    delta_threshold = pilot_delta_threshold,
    low_lambda_max = pilot_low_lambda_max,
    smooth_lambda = lambda,
    smooth_projection_rmse_min = pilot_smooth_projection_rmse_min,
    smooth_projection_rmse = smooth_projection_rmse,
    fixed_projection_rmse_min = pilot_fixed_projection_rmse_min,
    fixed_projection_rmse = projection_rmse,
    fixed_parent_terms_required = pilot_fixed_parent_terms,
    fixed_parent_terms_observed = parent_term_count,
    trigger = paste(c(
      if (low_lambda_smooth) "low_lambda_smooth",
      if (high_projection_smooth) "high_projection_smooth",
      if (high_projection_fixed) "high_projection_fixed"
    ), collapse = "+"),
    max_outer_iter = pilot_outer_iter,
    bfgs_max_iter = pilot_bfgs_max_iter
  )
}

low_lambda_inner_override_for_task <- function(
    step,
    parameter,
    dropped_term,
    change_type,
    start) {
  if (!identical(stage, "refine") || as.integer(step) <= 1L) return(NULL)
  if (!identical(change_type, "drop")) return(NULL)
  lambda <- smooth_penalty_for_term(start, parameter, dropped_term)
  if (!is.finite(lambda) || lambda > pilot_low_lambda_max) return(NULL)
  max(as.integer(initial_inner_iter), as.integer(low_lambda_initial_inner_iter))
}

make_task <- function(
    step,
    parameter,
    dropped_term,
    formulas,
    start,
    change_type = "drop") {
  dominance_pilot <- dominance_pilot_for_task(
    step, parameter, dropped_term, change_type, start
  )
  rs_initial_inner_iter_override <- low_lambda_inner_override_for_task(
    step, parameter, dropped_term, change_type, start
  )
  project_piloted_smooth_drop <- identical(stage, "refine") &&
    identical(change_type, "drop") && length(dropped_term) == 1L &&
    !is.na(dropped_term) && grepl("^\\s*s\\(", dropped_term) &&
    (is.list(dominance_pilot) ||
      length(rs_initial_inner_iter_override) == 1L)
  start_strategy <- screening_start_strategy(
    copula_code, step, stage, parameter, dropped_term, change_type,
    project_smooth_drop = project_piloted_smooth_drop
  )
  if (identical(start_strategy, "zero_screened_parameter")) {
    start <- zero_screened_parameter_start(start, parameter)
  } else if (identical(start_strategy, "project_parent_predictor")) {
    start <- project_smooth_drop_start(start, parameter, dropped_term)
  }
  metadata <- list(
    step = as.integer(step), parameter = parameter, dropped_term = dropped_term,
    change_type = change_type,
    formulas = formula_signature(formulas), start_signature = start_signature(start),
    margin = margin_family, copula = copula_code, missingness = "segment",
    optimizer_control = optimizer_control,
    package_source_fingerprint = package_source_fingerprint,
    r_version = R.version.string,
    gamlss_version = as.character(utils::packageVersion("gamlss")),
    gamlss_dist_version = as.character(utils::packageVersion("gamlss.dist")),
    n_rows = nrow(actg_data),
    n_subjects = length(unique(actg_data$subject)),
    n_observed = sum(is.finite(actg_data$cd4))
  )
  if (is.list(dominance_pilot)) {
    metadata$dominance_pilot <- dominance_pilot
  }
  if (length(rs_initial_inner_iter_override) == 1L) {
    metadata$rs_initial_inner_iter_override <- rs_initial_inner_iter_override
  }
  # Omit the ordinary strategy to preserve all compatible fitted-parent
  # checkpoints. The exceptional strategy is part of the key so the known
  # poor-basin checkpoint cannot be reused.
  if (!identical(start_strategy, "fitted_parent")) {
    metadata$start_strategy <- start_strategy
  }
  signature <- substr(object_md5(metadata), 1L, 12L)
  label <- if (is.na(dropped_term)) "full" else paste0("drop-", safe_name(dropped_term))
  checkpoint_stem <- sprintf("step-%02d-%s-%s", step, parameter, label)
  list(
    step = as.integer(step), parameter = parameter, dropped_term = dropped_term,
    change_type = change_type, formulas = formulas, start = start, metadata = metadata,
    checkpoint = file.path(
      checkpoint_dir,
      paste0(checkpoint_stem, "-", signature, ".rds")
    ),
    checkpoint_pattern = paste0("^", checkpoint_stem, "-[[:xdigit:]]{12}\\.rds$")
  )
}

candidate_row <- function(
    result,
    task,
    accepted = FALSE,
    parent_metrics = NULL,
    criterion_rank = NA_integer_) {
  if (is.null(parent_metrics)) parent_metrics <- empty_metrics()
  data.frame(
    model_id = tools::file_path_sans_ext(basename(task$checkpoint)),
    step = task$step,
    parameter = task$parameter,
    dropped_term = task$dropped_term,
    change_type = task$change_type,
    candidate_formula = formula_text(task$formulas[[task$parameter]]),
    mu_formula = formula_text(task$formulas$mu),
    sigma_formula = formula_text(task$formulas$sigma),
    nu_formula = formula_text(task$formulas$nu),
    tau_formula = if (has_tau) formula_text(task$formulas$tau) else NA_character_,
    theta_formula = formula_text(task$formulas$theta),
    zeta_formula = if (has_zeta) formula_text(task$formulas$zeta) else NA_character_,
    logLik = unname(result$metrics[["logLik"]]),
    AIC = unname(result$metrics[["AIC"]]),
    BIC = unname(result$metrics[["BIC"]]),
    EDF = unname(result$metrics[["EDF"]]),
    delta_AIC_from_parent = unname(result$metrics[["AIC"]] - parent_metrics[["AIC"]]),
    delta_BIC_from_parent = unname(result$metrics[["BIC"]] - parent_metrics[["BIC"]]),
    selection_criterion = criterion,
    criterion_rank_within_step = criterion_rank,
    converged = result$converged,
    accepted = accepted,
    elapsed_sec = result$elapsed_sec,
    outer_iterations = result$outer_iterations,
    stop_reason = result$stop_reason,
    dominance_pilot_rejected = isTRUE(result$dominance_pilot_rejected),
    from_checkpoint = result$from_checkpoint,
    warnings = paste(result$warnings, collapse = " | "),
    error = result$error,
    checkpoint = result$checkpoint,
    stringsAsFactors = FALSE
  )
}

write_rows <- function(rows, path) {
  output <- if (length(rows)) do.call(rbind, rows) else data.frame()
  temporary <- tempfile("results-", tmpdir = dirname(path), fileext = ".csv")
  utils::write.csv(output, temporary, row.names = FALSE, na = "")
  copied <- file.copy(temporary, path, overwrite = TRUE)
  unlink(temporary)
  if (!copied) stop("Could not write results: ", path)
  invisible(output)
}

initial_fit <- NULL
if (nzchar(initial_fit_argument)) {
  if (!file.exists(initial_fit_argument)) {
    stop("--initial-fit does not exist: ", initial_fit_argument)
  }
  initial_container <- readRDS(initial_fit_argument)
  initial_fit <- if (inherits(initial_container, "gamlss.longitudinal")) {
    initial_container
  } else if (is.list(initial_container) &&
             inherits(initial_container$fit, "gamlss.longitudinal")) {
    initial_container$fit
  } else if (is.list(initial_container)) {
    initial_container$value
  } else {
    NULL
  }
  if (!inherits(initial_fit, "gamlss.longitudinal")) {
    stop("--initial-fit must contain a gamlss.longitudinal fitted object.")
  }
  if (!identical(initial_fit$margin_dist$family[[1L]], margin_family) ||
      !identical(initial_fit$copula_dist, copula_code)) {
    stop("--initial-fit must use the requested margin and copula route.")
  }
}

n_subjects <- length(unique(actg_data$subject))
n_observed <- sum(is.finite(actg_data$cd4))
cat(
  "ACTG backward screening\n",
  "  subjects: ", n_subjects, " (", nrow(actg_data), " rows; ", n_observed,
  " observed responses)\n",
  "  family/copula: ", margin_family, "/", copula_label, "\n",
  "  criterion/hierarchy: ", criterion, "/", hierarchy, "\n",
  "  stage/smooth k: ", stage, "/", smooth_k, "\n",
  "  parameters: ", paste(parameters, collapse = ", "), "\n",
  "  workers: ", workers, "\n",
  "  dominance pilot: refinement smooth lambda <= ", pilot_low_lambda_max,
  " or smooth-drop projection RMSE >= ", pilot_smooth_projection_rmse_min,
  "; sparse fixed-drop projection RMSE >= ", pilot_fixed_projection_rmse_min,
  " at ", pilot_fixed_parent_terms, " parent terms",
  ", delta ", criterion, " > ", pilot_delta_threshold,
  ", outer/BFGS = ", pilot_outer_iter, "/", pilot_bfgs_max_iter, "\n",
  "  later low-lambda initial inner iterations: ",
  low_lambda_initial_inner_iter, "\n",
  "  checkpoints: ", checkpoint_dir, "\n",
  sep = ""
)

if (dry_run) {
  cat("\nDry run: no models were fitted.\n")
  print(maximum_formulas)
  quit(save = "no", status = 0L)
}

cluster <- NULL
if (workers > 1L) {
  cluster <- parallel::makePSOCKcluster(
    workers,
    outfile = file.path(output_dir, paste0(run_stem, "-workers.log")),
    setup_strategy = "parallel"
  )
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(root) {
    setwd(root)
    Sys.setenv(
      OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
    )
    devtools::load_all(root, quiet = TRUE)
    TRUE
  }, package_root)
  parallel::clusterExport(
    cluster,
    c(
      "actg_data", "optimizer_control", "resume", "empty_metrics",
      "candidate_checkpoint_paths", "checkpoint_metadata_compatible",
      "checkpoint_metadata_without_elapsed_cap", "fit_screening_candidate",
      "margin_family", "copula_code", "has_tau", "has_zeta"
    ),
    envir = environment()
  )
}

dispatch_tasks <- function(tasks) {
  if (!length(tasks)) return(list())
  if (is.null(cluster)) {
    lapply(tasks, fit_screening_candidate)
  } else {
    parallel::parLapplyLB(cluster, tasks, fit_screening_candidate)
  }
}

run_tasks <- function(tasks) {
  isolated <- which(vapply(tasks, function(task) {
    length(task$metadata$rs_initial_inner_iter_override) == 1L
  }, logical(1L)))
  if (!length(isolated)) return(dispatch_tasks(tasks))

  results <- vector("list", length(tasks))
  parallel_indices <- setdiff(seq_along(tasks), isolated)
  if (length(parallel_indices)) {
    results[parallel_indices] <- dispatch_tasks(tasks[parallel_indices])
  }
  for (index in isolated) {
    results[[index]] <- dispatch_tasks(tasks[index])[[1L]]
  }
  results
}

candidate_rows <- list()
path_rows <- list()
current_formulas <- maximum_formulas
full_task <- make_task(
  0L, parameters[[1L]], NA_character_, current_formulas, initial_fit, "maximum"
)
cat("\nStep 0: fitting the maximum model...\n")
current_summary <- run_tasks(list(full_task))[[1L]]
if (!isTRUE(current_summary$converged)) {
  stop(
    "The maximum model did not converge (", current_summary$stop_reason, "). ",
    if (!is.na(current_summary$error)) current_summary$error else "Inspect its checkpoint and warnings."
  )
}
current_result <- load_checkpoint_fit(current_summary, full_task)
candidate_rows[[1L]] <- candidate_row(
  current_result, full_task, TRUE, criterion_rank = 1L
)
path_rows[[1L]] <- candidate_row(
  current_result, full_task, TRUE, criterion_rank = 1L
)
write_rows(
  candidate_rows,
  file.path(output_dir, paste0(run_stem, "-model-selection-table.csv"))
)
write_rows(path_rows, file.path(output_dir, paste0(run_stem, "-path.csv")))
rm(current_summary, full_task, initial_fit)
gc(verbose = FALSE)

global_step <- 0L
for (parameter in parameters) {
  parameter_step <- 0L
  repeat {
    if (parameter_step >= max_steps) break
    moves <- formula_moves(current_formulas[[parameter]], hierarchy)
    if (!length(moves)) break
    global_step <- global_step + 1L
    parameter_step <- parameter_step + 1L
    tasks <- lapply(moves, function(move) {
      formulas <- current_formulas
      formulas[[parameter]] <- move$formula
      make_task(
        global_step, parameter, move$label, formulas,
        start_for_formulas(current_result$fit, formulas, actg_data),
        move$type
      )
    })
    cat(
      "\nStep ", global_step, " (", parameter, "): fitting ", length(tasks),
      " backward move", if (length(tasks) == 1L) "" else "s", "...\n",
      sep = ""
    )
    results <- run_tasks(tasks)
    parent_metrics <- current_result$metrics
    values <- vapply(results, function(result) {
      if (isTRUE(result$converged)) result$metrics[[criterion]] else Inf
    }, numeric(1L))
    within_step_rank <- rep(NA_integer_, length(values))
    finite_values <- is.finite(values)
    within_step_rank[finite_values] <- as.integer(rank(
      values[finite_values], ties.method = "min"
    ))
    best <- which.min(values)
    current_value <- parent_metrics[[criterion]]
    improved <- length(best) && is.finite(values[[best]]) && values[[best]] < current_value
    for (index in seq_along(results)) {
      candidate_rows[[length(candidate_rows) + 1L]] <- candidate_row(
        results[[index]], tasks[[index]], improved && index == best,
        parent_metrics = parent_metrics,
        criterion_rank = within_step_rank[[index]]
      )
    }
    if (improved) {
      current_result <- load_checkpoint_fit(results[[best]], tasks[[best]])
      current_formulas <- tasks[[best]]$formulas
      path_rows[[length(path_rows) + 1L]] <- candidate_row(
        current_result, tasks[[best]], TRUE,
        parent_metrics = parent_metrics,
        criterion_rank = within_step_rank[[best]]
      )
      cat(
        "Accepted ", tasks[[best]]$change_type, " of ", tasks[[best]]$dropped_term,
        ": ", criterion, " ",
        format(current_value, digits = 10), " -> ", format(values[[best]], digits = 10),
        "\n", sep = ""
      )
    } else {
      cat("No deletion improved ", criterion, "; ", parameter, " selection stopped.\n", sep = "")
    }
    write_rows(
      candidate_rows,
      file.path(output_dir, paste0(run_stem, "-model-selection-table.csv"))
    )
    write_rows(path_rows, file.path(output_dir, paste0(run_stem, "-path.csv")))
    rm(results, tasks)
    gc(verbose = FALSE)
    if (!improved) break
  }
}

final_fit_path <- file.path(output_dir, paste0(run_stem, "-final-fit.rds"))
saveRDS(current_result$fit, final_fit_path, compress = "gzip")
saveRDS(
  current_formulas,
  file.path(output_dir, paste0(run_stem, "-final-formulas.rds"))
)
model_table <- write_rows(
  candidate_rows,
  file.path(output_dir, paste0(run_stem, "-model-selection-table.csv"))
)
path_table <- write_rows(
  path_rows,
  file.path(output_dir, paste0(run_stem, "-path.csv"))
)
run_manifest_path <- file.path(output_dir, paste0(run_stem, "-run-manifest.rds"))
saveRDS(list(
  created_at = Sys.time(),
  selection_criterion = criterion,
  stage = stage,
  smooth_k = smooth_k,
  hierarchy = hierarchy,
  parameters = parameters,
  margin_family = margin_family,
  copula_family = copula_label,
  missingness = "segment",
  optimizer_control = optimizer_control,
  data_summary = list(
    rows = nrow(actg_data),
    subjects = length(unique(actg_data$subject)),
    observed_responses = sum(is.finite(actg_data$cd4))
  ),
  maximum_formulas = maximum_formulas,
  input_formulas = if (nzchar(input_formulas_argument)) normalizePath(input_formulas_argument) else NULL,
  initial_fit = if (nzchar(initial_fit_argument)) normalizePath(initial_fit_argument) else NULL,
  selected_formulas = current_formulas,
  model_table = model_table,
  accepted_path = path_table,
  model_checkpoints = stats::setNames(model_table$checkpoint, model_table$model_id),
  final_fit_path = final_fit_path,
  package_source_fingerprint = package_source_fingerprint
), run_manifest_path, compress = "gzip")

cat("\nSelected formulas:\n")
print(current_formulas)
cat(
  "\nFinal ", criterion, ": ", format(current_result$metrics[[criterion]], digits = 10),
  "\nSaved final fit: ", final_fit_path,
  "\nSaved selection path: ",
  file.path(output_dir, paste0(run_stem, "-path.csv")), "\n",
  "Saved model-selection table: ",
  file.path(output_dir, paste0(run_stem, "-model-selection-table.csv")), "\n",
  "Saved run manifest: ", run_manifest_path, "\n",
  sep = ""
)
