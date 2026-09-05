#!/usr/bin/env Rscript

# Exploratory comparison of four public medical longitudinal datasets.
# This is deliberately outside the manuscript targets graph while the paper
# application is being selected.

if (file.exists("DESCRIPTION") &&
    any(grepl("^Package: gamlss.longitudinal$", readLines("DESCRIPTION", n = 2L)))) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Install devtools to run this analysis from the package checkout.")
  }
  devtools::load_all(quiet = TRUE)
} else if (!requireNamespace("gamlss.longitudinal", quietly = TRUE)) {
  stop("Install gamlss.longitudinal or run this script from its package checkout.")
}

output_dir <- file.path("paper", "reviews", "2026-09-05-medical-candidates", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

family_object <- function(code) getExportedValue("gamlss.dist", code)()

read_tlc <- function() {
  wide <- utils::read.table(
    "https://content.sph.harvard.edu/fitzmaur/ala/tlc.dat",
    col.names = c("subject", "treatment_code", "week0", "week1", "week4", "week6")
  )
  long <- stats::reshape(
    wide,
    varying = c("week0", "week1", "week4", "week6"),
    v.names = "lead",
    timevar = "week",
    times = c(0, 1, 4, 6),
    idvar = "subject",
    direction = "long"
  )
  long <- long[order(long$subject, long$week), ]
  rownames(long) <- NULL
  long$treatment <- factor(
    long$treatment_code,
    levels = c("A", "P"),
    labels = c("succimer", "placebo")
  )
  long$week_factor <- factor(long$week)
  long$week_placebo <- long$week * as.integer(long$treatment == "placebo")
  long
}

read_actg <- function() {
  if (!requireNamespace("bcmixed", quietly = TRUE)) {
    stop("The exploratory ACTG analysis requires the CRAN package 'bcmixed'.")
  }
  environment <- new.env(parent = emptyenv())
  utils::data("aidscd4", package = "bcmixed", envir = environment)
  data <- get("aidscd4", envir = environment)
  data <- data[order(data$id, data$weekc), ]
  data$treatment <- factor(data$treatment)
  data$week_factor <- factor(data$weekc)
  data$baseline_log_z <- as.numeric(scale(log1p(data$cd4.bl)))
  data$age_z <- as.numeric(scale(data$age))
  data$sex <- factor(data$sex, levels = c(0, 1), labels = c("female", "male"))
  data
}

read_skin <- function() {
  lines <- readLines(
    "https://content.sph.harvard.edu/fitzmaur/ala/skin.txt",
    warn = FALSE
  )
  lines <- lines[grepl("^\\s*[0-9]+\\s+", lines)]
  data <- utils::read.table(
    text = lines,
    col.names = c(
      "subject", "center", "age", "skin_type", "sex", "exposure",
      "cancers", "treatment_code", "year"
    )
  )
  data <- data[order(data$subject, data$year), ]
  data$treatment <- factor(
    data$treatment_code,
    levels = c(0, 1),
    labels = c("placebo", "beta_carotene")
  )
  data$sex <- factor(data$sex, levels = c(0, 1), labels = c("female", "male"))
  data$skin_type <- factor(data$skin_type, levels = c(0, 1), labels = c("other", "burns"))
  data$year_factor <- factor(data$year)
  data$center_factor <- factor(data$center)
  data$exposure_z <- as.numeric(scale(log1p(data$exposure)))
  data$age_z <- as.numeric(scale(data$age))
  data
}

read_epilepsy <- function() {
  environment <- new.env(parent = emptyenv())
  utils::data("epil", package = "MASS", envir = environment)
  data <- get("epil", envir = environment)
  data <- data[order(data$subject, data$period), ]
  data$period_factor <- factor(data$period)
  data$treatment <- data$trt
  data$baseline_log_z <- as.numeric(scale(log1p(data$base)))
  data$age_z <- as.numeric(scale(data$age))
  data
}

candidate_specs <- list(
  tlc = list(
    read = read_tlc,
    response = "lead", subject = "subject", time = "week",
    margins = c("NO", "GA", "LOGNO", "GG"),
    mu = lead ~ week_factor + treatment + week_placebo,
    parameter_candidates = list(
      constant = list(sigma = ~ 1, nu = ~ 1),
      treatment_sigma = list(sigma = ~ treatment, nu = ~ 1),
      treatment_nu = list(sigma = ~ 1, nu = ~ treatment),
      treatment_sigma_nu = list(sigma = ~ treatment, nu = ~ treatment)
    )
  ),
  actg = list(
    read = read_actg,
    response = "cd4", subject = "id", time = "weekc",
    margins = c("PO", "NBI", "DEL", "ZINBI"),
    mu = cd4 ~ week_factor + treatment + baseline_log_z + age_z + sex,
    parameter_candidates = list(
      constant = list(sigma = ~ 1, nu = ~ 1),
      treatment_sigma = list(sigma = ~ treatment, nu = ~ 1),
      baseline_sigma = list(sigma = ~ baseline_log_z, nu = ~ 1),
      treatment_baseline_sigma = list(sigma = ~ treatment + baseline_log_z, nu = ~ 1),
      treatment_nu = list(sigma = ~ 1, nu = ~ treatment),
      baseline_nu = list(sigma = ~ 1, nu = ~ baseline_log_z)
    )
  ),
  skin = list(
    read = read_skin,
    response = "cancers", subject = "subject", time = "year",
    margins = c("PO", "NBI", "DEL", "ZINBI"),
    mu = cancers ~ year_factor + treatment + exposure_z + age_z + sex + skin_type,
    parameter_candidates = list(
      constant = list(sigma = ~ 1, nu = ~ 1),
      treatment_sigma = list(sigma = ~ treatment, nu = ~ 1),
      exposure_sigma = list(sigma = ~ exposure_z, nu = ~ 1),
      treatment_exposure_sigma = list(sigma = ~ treatment + exposure_z, nu = ~ 1),
      age_sigma = list(sigma = ~ age_z, nu = ~ 1),
      sex_sigma = list(sigma = ~ sex, nu = ~ 1),
      skin_type_sigma = list(sigma = ~ skin_type, nu = ~ 1),
      center_sigma = list(sigma = ~ center_factor, nu = ~ 1),
      treatment_nu = list(sigma = ~ 1, nu = ~ treatment),
      exposure_nu = list(sigma = ~ 1, nu = ~ exposure_z),
      age_nu = list(sigma = ~ 1, nu = ~ age_z),
      sex_nu = list(sigma = ~ 1, nu = ~ sex),
      skin_type_nu = list(sigma = ~ 1, nu = ~ skin_type),
      center_nu = list(sigma = ~ 1, nu = ~ center_factor)
    )
  ),
  epilepsy = list(
    read = read_epilepsy,
    response = "y", subject = "subject", time = "period",
    margins = c("PO", "NBI", "DEL", "ZINBI"),
    mu = y ~ period_factor + treatment + baseline_log_z + age_z +
      treatment:baseline_log_z,
    parameter_candidates = list(
      constant = list(sigma = ~ 1, nu = ~ 1),
      treatment_sigma = list(sigma = ~ treatment, nu = ~ 1),
      baseline_sigma = list(sigma = ~ baseline_log_z, nu = ~ 1),
      treatment_baseline_sigma = list(sigma = ~ treatment + baseline_log_z, nu = ~ 1),
      treatment_nu = list(sigma = ~ 1, nu = ~ treatment),
      baseline_nu = list(sigma = ~ 1, nu = ~ baseline_log_z)
    )
  )
)

capture_fit <- function(expression) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(
      expression,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )
  list(
    value = value,
    warnings = unique(warnings),
    elapsed = proc.time()[["elapsed"]] - started
  )
}

screen_margin_parameters <- function(data, spec, family_code) {
  family <- family_object(family_code)
  parameter_names <- names(family$parameters)
  needed <- unique(c(
    all.vars(spec$mu),
    unlist(lapply(spec$parameter_candidates, function(x) {
      unlist(lapply(x, all.vars), use.names = FALSE)
    }), use.names = FALSE)
  ))
  analysis_data <- data[stats::complete.cases(data[needed]), , drop = FALSE]

  rows <- lapply(names(spec$parameter_candidates), function(model) {
    parameter_spec <- spec$parameter_candidates[[model]]
    args <- list(
      formula = spec$mu,
      family = family,
      data = analysis_data,
      trace = FALSE,
      control = gamlss::gamlss.control(n.cyc = 300L, trace = FALSE)
    )
    if ("sigma" %in% parameter_names) args$sigma.formula <- parameter_spec$sigma
    if ("nu" %in% parameter_names) args$nu.formula <- parameter_spec$nu
    attempt <- capture_fit(do.call(gamlss::gamlss, args))
    failed <- inherits(attempt$value, "error")
    data.frame(
      model = model,
      family = family_code,
      n = nrow(analysis_data),
      converged = if (failed) FALSE else isTRUE(attempt$value$converged),
      df = if (failed) NA_real_ else attempt$value$df.fit,
      logLik = if (failed) NA_real_ else as.numeric(stats::logLik(attempt$value)),
      AIC = if (failed) NA_real_ else stats::AIC(attempt$value),
      BIC = if (failed) NA_real_ else stats::BIC(attempt$value),
      elapsed_sec = attempt$elapsed,
      warnings = paste(attempt$warnings, collapse = " | "),
      error = if (failed) conditionMessage(attempt$value) else "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$AIC, na.last = TRUE), ]
}

dataset_summary <- function(data, spec, name) {
  response <- data[[spec$response]]
  observed <- is.finite(response)
  visits <- table(data[[spec$subject]][observed])
  adjacent <- do.call(rbind, lapply(split(data, data[[spec$subject]]), function(subject_data) {
    subject_data <- subject_data[order(subject_data[[spec$time]]), ]
    if (nrow(subject_data) < 2L) return(NULL)
    data.frame(
      left = subject_data[[spec$response]][-nrow(subject_data)],
      right = subject_data[[spec$response]][-1L]
    )
  }))
  data.frame(
    dataset = name,
    subjects = length(unique(data[[spec$subject]])),
    rows = nrow(data),
    observed_responses = sum(observed),
    missing_responses = sum(!observed),
    median_observed_visits = stats::median(as.numeric(visits)),
    max_observed_visits = max(visits),
    adjacent_spearman = stats::cor(
      adjacent$left, adjacent$right, method = "spearman", use = "complete.obs"
    ),
    response_mean = mean(response[observed]),
    response_variance = stats::var(response[observed]),
    zero_fraction = mean(response[observed] == 0),
    stringsAsFactors = FALSE
  )
}

run_one <- function(name) {
  message("Screening ", name)
  spec <- candidate_specs[[name]]
  data <- spec$read()
  summary <- dataset_summary(data, spec, name)

  cache_path <- file.path(output_dir, paste0(name, "-screening.rds"))
  cached <- if (file.exists(cache_path)) readRDS(cache_path) else NULL
  cache_current <- !is.null(cached) &&
    identical(attr(cached$copula, "screen_version"), "2026.1")
  if (cache_current) {
    message("Reusing current margin and copula screens for ", name)
    margin <- cached$margin
    copula <- cached$copula
  } else {
    margin <- suppressMessages(gamlss.longitudinal::select_margin(
      data,
      response_var = spec$response,
      time_var = spec$time,
      time_intercepts = TRUE,
      families = spec$margins,
      time_n_cycles = 300L,
      trace = FALSE
    ))
    copula <- suppressMessages(gamlss.longitudinal::select_copula(
      data = data,
      response_var = spec$response,
      margin_dist = margin,
      subject_var = spec$subject,
      time_var = spec$time,
      families = c("N", "C", "F", "G", "J", "t"),
      time_intercepts = TRUE,
      margin_n_cycles = 300L,
      pit_seed = 20260905L
    ))
  }
  selected_family <- gamlss.longitudinal::best_fit(margin)$family_name

  parameter_screen <- screen_margin_parameters(data, spec, selected_family)
  utils::write.csv(as.data.frame(margin), file.path(output_dir, paste0(name, "-margin.csv")), row.names = FALSE)
  utils::write.csv(as.data.frame(copula), file.path(output_dir, paste0(name, "-copula.csv")), row.names = FALSE)
  utils::write.csv(parameter_screen, file.path(output_dir, paste0(name, "-margin-covariates.csv")), row.names = FALSE)
  saveRDS(
    list(
      summary = summary,
      margin = margin,
      copula = copula,
      margin_covariates = parameter_screen,
      selected_margin = selected_family,
      selected_copula = gamlss.longitudinal::best_fit_family(copula)
    ),
    file.path(output_dir, paste0(name, "-screening.rds"))
  )
  summary$selected_margin <- selected_family
  summary$selected_copula <- gamlss.longitudinal::best_fit_family(copula)
  summary$copula_screen_likelihood <- attr(copula, "screen_likelihood")
  summary
}

if (!identical(Sys.getenv("MEDICAL_CANDIDATE_DEFINE_ONLY"), "true")) {
  args <- commandArgs(trailingOnly = TRUE)
  requested <- if (length(args)) args else names(candidate_specs)
  unknown <- setdiff(requested, names(candidate_specs))
  if (length(unknown)) stop("Unknown dataset(s): ", paste(unknown, collapse = ", "))

  summaries <- lapply(requested, run_one)
  summary <- do.call(rbind, summaries)
  utils::write.csv(summary, file.path(output_dir, "dataset-summary.csv"), row.names = FALSE)
  print(summary)
}
