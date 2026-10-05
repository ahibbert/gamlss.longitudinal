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
max_elapsed_sec <- suppressWarnings(as.numeric(argument_value("max-elapsed-sec", "1200")))
if (length(max_elapsed_sec) != 1L || is.na(max_elapsed_sec) || max_elapsed_sec <= 0) {
  stop("--max-elapsed-sec must be a positive number.")
}
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

route_refinement <- switch(
  paste(margin_family, copula_code, sep = "/"),
  `ZISICHEL/G` = list(after = 10L, iterations = 100L, refinements = 1L,
                       rescue = FALSE, repolish = 50L),
  `ZIBNB/G` = list(after = 2L, iterations = 170L, refinements = 1L,
                    rescue = FALSE, repolish = 50L),
  `ZISICHEL/t` = list(after = 1L, iterations = 120L, refinements = 1L,
                       rescue = FALSE, repolish = 50L),
  `ZIBNB/t` = list(after = 1L, iterations = 40L, refinements = 2L,
                    rescue = TRUE, repolish = 40L),
  list(after = 10L, iterations = 100L, refinements = 1L,
       rescue = FALSE, repolish = 50L)
)
target_route <- margin_family %in% c("ZISICHEL", "ZIBNB")
optimizer_control <- gamlss_longitudinal_control(
  outer_tol = 1e-4,
  max_outer_iter = 50L,
  max_elapsed_sec = max_elapsed_sec,
  rs = list(
    inner_tol = 1e-5,
    max_inner_iter = 20L,
    initial_inner_iter = 5L,
    start_step_size = 2,
    warm_start_joint = TRUE,
    warm_start_joint_iter = if (target_route) 2L else 5L,
    use_backtracking = TRUE,
    backtracking_max_halves = 2L,
    adaptive_step_size = FALSE,
    update_lambda = TRUE,
    lambda_update_every = 5L,
    bfgs_refine = target_route,
    bfgs_refine_after = route_refinement$after,
    bfgs_max_refinements = route_refinement$refinements,
    bfgs_max_iter = route_refinement$iterations,
    bfgs_reltol = 1e-12,
    bfgs_coordinate_rescue = route_refinement$rescue,
    bfgs_repolish_iter = route_refinement$repolish,
    discrete_score_method = "analytical"
  )
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

fit_screening_candidate <- function(task) {
  if (isTRUE(resume) && file.exists(task$checkpoint)) {
    cached <- tryCatch(readRDS(task$checkpoint), error = identity)
    cache_is_reusable <- is.list(cached) &&
      identical(cached$metadata, task$metadata) &&
      isTRUE(cached$converged) &&
      inherits(cached$fit, "gamlss.longitudinal")
    if (cache_is_reusable) {
      cached$from_checkpoint <- TRUE
      return(cached)
    }
  }

  started <- proc.time()[["elapsed"]]
  warnings <- character()
  optional_formulas <- c(
    if (has_tau) list(tau.formula = task$formulas$tau) else list(),
    if (has_zeta) list(zeta.formula = task$formulas$zeta) else list()
  )
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
        optimizer_control = optimizer_control,
        start = task$start,
        compute_vcov = FALSE,
        verbose = 0
      ), optional_formulas)),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )
  elapsed <- proc.time()[["elapsed"]] - started
  returned_fit <- inherits(fit, "gamlss.longitudinal")
  converged <- returned_fit && isTRUE(fit$convergence$converged)
  metrics <- empty_metrics()
  if (returned_fit && !is.null(fit$model_selection)) {
    metrics <- c(
      logLik = as.numeric(fit$model_selection["LogLik", "joint"]),
      AIC = as.numeric(fit$model_selection["AIC", "joint"]),
      BIC = as.numeric(fit$model_selection["BIC", "joint"]),
      EDF = as.numeric(fit$model_selection["EDF", "joint"])
    )
  }
  result <- list(
    metadata = task$metadata,
    fit = if (returned_fit) fit else NULL,
    metrics = metrics,
    converged = converged,
    elapsed_sec = elapsed,
    outer_iterations = if (returned_fit) fit$convergence$outer_iterations else NA_integer_,
    stop_reason = if (returned_fit) fit$convergence$stop_reason else "error",
    warnings = unique(warnings),
    error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
    checkpoint = task$checkpoint,
    from_checkpoint = FALSE
  )
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
  result
}

make_task <- function(
    step,
    parameter,
    dropped_term,
    formulas,
    start,
    change_type = "drop") {
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
  signature <- substr(object_md5(metadata), 1L, 12L)
  label <- if (is.na(dropped_term)) "full" else paste0("drop-", safe_name(dropped_term))
  list(
    step = as.integer(step), parameter = parameter, dropped_term = dropped_term,
    change_type = change_type, formulas = formulas, start = start, metadata = metadata,
    checkpoint = file.path(
      checkpoint_dir,
      sprintf("step-%02d-%s-%s-%s.rds", step, parameter, label, signature)
    )
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
      "fit_screening_candidate", "margin_family", "copula_code", "has_tau", "has_zeta"
    ),
    envir = environment()
  )
}

run_tasks <- function(tasks) {
  if (is.null(cluster)) {
    lapply(tasks, fit_screening_candidate)
  } else {
    parallel::parLapplyLB(cluster, tasks, fit_screening_candidate)
  }
}

candidate_rows <- list()
path_rows <- list()
current_formulas <- maximum_formulas
full_task <- make_task(
  0L, parameters[[1L]], NA_character_, current_formulas, initial_fit, "maximum"
)
cat("\nStep 0: fitting the maximum model...\n")
current_result <- run_tasks(list(full_task))[[1L]]
if (!isTRUE(current_result$converged)) {
  stop(
    "The maximum model did not converge (", current_result$stop_reason, "). ",
    if (!is.na(current_result$error)) current_result$error else "Inspect its checkpoint and warnings."
  )
}
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
      current_result <- results[[best]]
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
