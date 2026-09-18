#!/usr/bin/env Rscript

# Crash-safe manual reproducer for count-margin/Student-t fits.
# Run this file in a disposable Rscript process. It deliberately does not catch
# native faults; the parent process should record the worker exit status.

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

family <- toupper(argument_value("family", "ZIPIG"))
copula <- argument_value("copula", "t")
n_subjects <- as.integer(argument_value("subjects", "25"))
max_outer <- as.integer(argument_value("outer", "1"))
max_inner <- as.integer(argument_value("inner", "1"))
max_elapsed <- as.numeric(argument_value("timeout", "300"))
score_method <- argument_value("score", "analytical")
backend <- argument_value("backend", "native")
stage <- argument_value("stage", "fit")
output <- argument_value("output", tempfile(fileext = ".rds"))

mark <- function(label, fields = list()) {
  values <- c(
    list(
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"),
      pid = Sys.getpid(),
      stage = label
    ),
    fields
  )
  cat(paste(names(values), unlist(values), sep = "=", collapse = " "), "\n")
  flush.console()
}

mark("startup", list(
  family = family,
  copula = copula,
  subjects = n_subjects,
  outer = max_outer,
  inner = max_inner,
  score = score_method,
  backend = backend,
  requested_stage = stage
))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)
options(gamlss.longitudinal.copula_backend = backend)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath("tests/manual/repro-count-t-stability.R", mustWork = TRUE)
}
package_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload.")
pkgload::load_all(package_root, quiet = TRUE)

data_env <- new.env(parent = emptyenv())
utils::data("actg_cd4", package = "gamlss.longitudinal", envir = data_env)
dat <- get("actg_cd4", envir = data_env, inherits = FALSE)
dat$week_factor <- factor(dat$week)
ids <- sort(unique(dat$subject))
if (is.finite(n_subjects) && n_subjects < length(ids)) {
  ids <- ids[seq_len(n_subjects)]
  dat <- dat[dat$subject %in% ids, , drop = FALSE]
}
mark("data_ready", list(rows = nrow(dat), observed = sum(is.finite(dat$cd4))))

margin_dist <- get(family, envir = asNamespace("gamlss.dist"))()
parameter_names <- names(margin_dist$parameters)[margin_dist$parameters %in% TRUE]
mark("margin_constructed", list(parameters = paste(parameter_names, collapse = ",")))

if (identical(stage, "margin")) {
  margin_data <- dat[is.finite(dat$cd4), , drop = FALSE]
  margin_fit <- gamlss::gamlss(
    cd4 ~ 1,
    family = margin_dist,
    data = margin_data,
    trace = FALSE,
    control = gamlss::gamlss.control(n.cyc = 20)
  )
  params <- lapply(parameter_names, function(parameter) {
    stats::fitted(margin_fit, what = parameter)
  })
  names(params) <- parameter_names
  call_args <- c(list(q = margin_data$cd4), params)
  mark("before_margin_cdf")
  probability <- do.call(
    get(paste0("p", family), envir = asNamespace("gamlss.dist")),
    call_args
  )
  mark("after_margin_cdf", list(finite = sum(is.finite(probability))))
  saveRDS(list(fit = margin_fit, probability = probability), output)
  quit(save = "no", status = 0L)
}

control <- gamlss_longitudinal_control(
  max_outer_iter = max_outer,
  max_elapsed_sec = max_elapsed,
  rs = list(
    max_inner_iter = max_inner,
    warm_start_joint = FALSE,
    discrete_score_method = score_method
  )
)

mark("before_fit")
result <- select_joint_distribution(
  data = dat,
  time_var = "week_factor",
  subject_var = "subject",
  response_var = "cd4",
  margin_families = family,
  copula_families = copula,
  trace = TRUE,
  progress = TRUE,
  keep_fits = TRUE,
  fit_args = list(
    optimizer_control = control,
    verbose = 1
  )
)
mark("after_fit", list(rows = nrow(result)))
saveRDS(result, output)
mark("saved", list(output = normalizePath(output, mustWork = TRUE)))
