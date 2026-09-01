.gl_validate_bootstrap_args <- function(object, R, fit_args) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  R <- as.integer(R)
  if (length(R) != 1L || is.na(R) || R < 1L) {
    stop("'R' must be a positive integer.", call. = FALSE)
  }

  if (!is.list(fit_args) || (length(fit_args) > 0L && (is.null(names(fit_args)) || any(names(fit_args) == "")))) {
    stop("'fit_args' must be a named list.", call. = FALSE)
  }

  list(R = R, fit_args = fit_args)
}

.gl_normalize_bootstrap_simulation_args <- function(dots) {
  if ("simulation_type" %in% names(dots)) {
    simulation_type <- dots$simulation_type

    if (!identical(simulation_type, "copula")) {
      stop(
        "'simulation_type' is no longer supported; bootstrap_inference() simulates from the fitted copula model.",
        call. = FALSE
      )
    }

    dots$simulation_type <- NULL
  }

  dots
}

.gl_cluster_bootstrap_dataset <- function(object) {
  if (is.null(object$dataset) || !"subject" %in% names(object$dataset)) {
    stop(
      "Cluster bootstrap requires a fitted object with a stored internal dataset and subject column.",
      call. = FALSE
    )
  }

  dat <- object$dataset
  subjects <- unique(dat$subject)
  draws <- sample(subjects, length(subjects), replace = TRUE)
  pieces <- vector("list", length(draws))

  for (i in seq_along(draws)) {
    piece <- dat[dat$subject == draws[[i]], , drop = FALSE]
    piece$subject <- i
    pieces[[i]] <- piece
  }

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

.gl_bootstrap_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }

  as.numeric(stats::quantile(x, probs = prob, names = FALSE))
}

.gl_bootstrap_summary <- function(estimates, terms_use, boot_coef, level) {
  alpha <- 1 - level

  data.frame(
    term = terms_use,
    estimate = as.numeric(estimates[terms_use]),
    bootstrap_mean = colMeans(boot_coef, na.rm = TRUE),
    bootstrap_se = apply(boot_coef, 2, stats::sd, na.rm = TRUE),
    conf.low = apply(boot_coef, 2, .gl_bootstrap_quantile, prob = alpha / 2),
    conf.high = apply(boot_coef, 2, .gl_bootstrap_quantile, prob = 1 - alpha / 2),
    reps = colSums(is.finite(boot_coef)),
    stringsAsFactors = FALSE
  )
}

.gl_bootstrap_result <- function(summary,
                                 boot_coef,
                                 errors,
                                 R,
                                 level,
                                 fits = NULL,
                                 simulation_type = "copula") {
  successful <- sum(stats::complete.cases(boot_coef))
  contract_id <- if (identical(simulation_type, "cluster")) {
    "bootstrap_cluster_fixed"
  } else {
    "bootstrap_parametric_fixed"
  }
  min_reps <- if ("reps" %in% names(summary) && nrow(summary) > 0L) min(summary$reps) else 0L
  contract <- .gl_inference_contract(
    contract_id,
    coefficient_names = colnames(boot_coef),
    method = simulation_type,
    validity_status = if (min_reps >= 2L) "available_with_reported_failures" else "insufficient_successful_replicates",
    failure_states = unique(stats::na.omit(errors))
  )
  contract$computational_minimum <- 2L
  contract$inferential_adequacy <- "not_assessed_from_replicate_count"
  contract$adequacy_note <- paste0(
    "Two successful replicates are the computational minimum for a finite ",
    "standard deviation, not evidence of adequate bootstrap inference."
  )
  out <- list(
    summary = summary,
    replicates = as.data.frame(boot_coef, stringsAsFactors = FALSE),
    errors = errors,
    R = R,
    successful_replicates = successful,
    failed_replicates = R - successful,
    level = level,
    simulation_type = simulation_type,
    fits = fits,
    inference_contract = contract
  )

  class(out) <- "gamlss_longitudinal_bootstrap"
  attr(out$summary, "inference_contract") <- contract
  attr(out$replicates, "inference_contract") <- contract
  attr(out, "inference_contract") <- contract
  out
}
