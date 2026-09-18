.joint_selection_resolve_isolation <- function(mode) {
  mode <- match.arg(mode, c("auto", "always", "never"))
  available <- requireNamespace("callr", quietly = TRUE)
  if (identical(mode, "always") && !available) {
    stop(
      "Candidate isolation requires package 'callr'.",
      call. = FALSE
    )
  }
  if (identical(mode, "auto") && !available) {
    warning(
      "Package 'callr' is not installed; joint candidates will run in the current R process.",
      call. = FALSE
    )
  }
  identical(mode, "always") || (identical(mode, "auto") && available)
}

.joint_selection_fit_one_isolated <- function(
    data,
    response_var,
    time_var,
    subject_var,
    margin_family,
    copula_family,
    fit_args,
    time_intercepts,
    copula_time_intercepts,
    trace,
    keep_fit,
    timeout) {
  started <- Sys.time()
  package_path <- path.package("gamlss.longitudinal")
  backend <- getOption("gamlss.longitudinal.copula_backend", "native")
  worker_args <- list(
    data = data,
    response_var = response_var,
    time_var = time_var,
    subject_var = subject_var,
    margin_family = margin_family,
    copula_family = copula_family,
    fit_args = fit_args,
    time_intercepts = time_intercepts,
    copula_time_intercepts = copula_time_intercepts,
    trace = trace,
    keep_fit = keep_fit
  )

  result <- tryCatch(
    callr::r(
      function(package_path, backend, worker_args) {
        if (file.exists(file.path(package_path, "R", "joint-distribution-selection-fit.R")) &&
            requireNamespace("pkgload", quietly = TRUE)) {
          suppressMessages(pkgload::load_all(package_path, quiet = TRUE))
        } else {
          loadNamespace("gamlss.longitudinal")
        }
        options(gamlss.longitudinal.copula_backend = backend)
        namespace <- asNamespace("gamlss.longitudinal")
        worker_args$margin_dist <- get(
          ".margin_family_object", envir = namespace, inherits = FALSE
        )(worker_args$margin_family)
        keep_fit <- worker_args$keep_fit
        worker_args$keep_fit <- NULL
        answer <- do.call(
          get(".joint_selection_fit_one", envir = namespace, inherits = FALSE),
          worker_args
        )
        if (!isTRUE(keep_fit)) answer$fit <- NULL
        answer
      },
      args = list(
        package_path = package_path,
        backend = backend,
        worker_args = worker_args
      ),
      libpath = .libPaths(),
      show = isTRUE(trace),
      timeout = timeout,
      user_profile = FALSE,
      system_profile = FALSE,
      env = c(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1"
      )
    ),
    error = identity
  )

  if (!inherits(result, "error")) return(result)

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  condition_text <- conditionMessage(result)
  timed_out <- inherits(result, "system_command_timeout_error") ||
    inherits(result, "callr_timeout_error") ||
    grepl("timed out|timeout", condition_text, ignore.case = TRUE)
  reason <- if (timed_out) "candidate_timeout" else "worker_failure"
  message <- if (timed_out) {
    paste0("Candidate worker exceeded the ", timeout, " second timeout.")
  } else {
    paste0("Candidate worker terminated: ", condition_text)
  }
  list(
    row = .joint_selection_failed_row(
      margin_family = margin_family,
      copula_family = copula_family,
      elapsed = elapsed,
      warnings = character(),
      error = message,
      stop_reason = reason
    ),
    fit = NULL
  )
}

.joint_selection_checkpoint_signature <- function(
    data, response_var, time_var, subject_var, combinations, fit_args,
    time_intercepts, copula_time_intercepts, keep_fits) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop(
      "Checkpointed joint selection requires the suggested package 'digest'.",
      call. = FALSE
    )
  }
  digest::digest(
    list(
      data = data[, c(response_var, time_var, subject_var), drop = FALSE],
      combinations = combinations,
      fit_args = fit_args,
      time_intercepts = time_intercepts,
      copula_time_intercepts = copula_time_intercepts,
      keep_fits = keep_fits
    ),
    algo = "xxhash64",
    serialize = TRUE
  )
}

.joint_selection_checkpoint_path <- function(directory, index, margin, copula) {
  safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "-", x)
  file.path(
    directory,
    sprintf("%03d-%s-%s.rds", index, safe(margin), safe(copula))
  )
}

.joint_selection_read_checkpoint <- function(path, signature) {
  if (!file.exists(path)) return(NULL)
  value <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(value) || !identical(value$signature, signature) ||
      !is.list(value$result) || is.null(value$result$row)) {
    return(NULL)
  }
  if (!is.na(value$result$row$error[[1L]])) return(NULL)
  value$result
}

.joint_selection_write_checkpoint <- function(path, signature, result) {
  temporary <- tempfile(tmpdir = dirname(path), fileext = ".rds")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(list(signature = signature, result = result), temporary)
  if (!file.rename(temporary, path)) {
    if (!file.copy(temporary, path, overwrite = TRUE)) {
      stop("Could not write joint-selection checkpoint: ", path, call. = FALSE)
    }
    unlink(temporary)
  }
  invisible(path)
}
