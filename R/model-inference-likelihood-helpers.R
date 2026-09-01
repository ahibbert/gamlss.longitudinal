#' Normalize and validate likelihood-comparison model inputs
#'
#' @param models List collected from `...`.
#' @param labels Optional labels captured from the unevaluated call.
#' @return Named list of fitted models.
#' @noRd
.gl_likelihood_compare_models <- function(models, labels = NULL) {
  if (length(models) == 1L && is.list(models[[1L]]) && !inherits(models[[1L]], "gamlss.longitudinal")) {
    models <- models[[1L]]
  }

  if (length(models) < 2L) {
    stop("At least two fitted models are required.", call. = FALSE)
  }

  ok <- vapply(models, inherits, logical(1), what = "gamlss.longitudinal")

  if (!all(ok)) {
    stop("All inputs must be fitted 'gamlss.longitudinal' objects.", call. = FALSE)
  }

  model_names <- names(models)

  if (is.null(model_names)) {
    model_names <- rep("", length(models))
  }

  if (!is.null(labels) && length(labels) == length(models)) {
    missing_names <- model_names == ""

    model_names[missing_names] <- labels[missing_names]
  }

  if (any(model_names == "")) {
    model_names[model_names == ""] <- paste0("model_", which(model_names == ""))
  }

  names(models) <- model_names
  models
}

#' Capture likelihood-comparison labels from the user call
#'
#' @param expr Unevaluated `list(...)` call from `likelihood_compare()`.
#' @return Character vector of labels, with empty strings for non-symbol inputs.
#' @noRd
.gl_likelihood_compare_call_labels <- function(expr) {
  args <- as.list(expr)[-1L]

  if (length(args) == 1L && is.call(args[[1L]]) && identical(args[[1L]][[1L]], as.name("list"))) {
    args <- as.list(args[[1L]])[-1L]
  }

  labels <- rep("", length(args))

  arg_names <- names(args)

  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  for (i in seq_along(args)) {
    if (nzchar(arg_names[[i]])) {
      labels[[i]] <- arg_names[[i]]
    } else if (is.symbol(args[[i]])) {
      labels[[i]] <- as.character(args[[i]])
    }
  }

  labels
}

#' Build the likelihood-comparison result table
#'
#' @param models Named list of fitted `gamlss.longitudinal` objects.
#' @param sort Deprecated logical; AIC sorting is not allowed for a sequential
#'   likelihood-ratio table.
#' @return Classed likelihood-comparison data frame.
#' @noRd
.gl_likelihood_compare_table <- function(models, sort = FALSE) {
  if (isTRUE(sort)) {
    stop(
      "'sort = TRUE' is not supported: supply models in the intended reduced-to-full order and inspect AIC separately.",
      call. = FALSE
    )
  }
  labels <- names(models)

  n_obs <- vapply(models, function(x) sum(is.finite(x$response)), integer(1))

  df <- vapply(models, .gl_model_edf, numeric(1))

  loglik <- vapply(models, .gl_joint_loglik, numeric(1))

  aic <- -2 * loglik + 2 * df

  bic <- -2 * loglik + log(pmax(1, n_obs)) * df

  delta_df <- c(NA_real_, diff(df))

  lr <- c(NA_real_, 2 * diff(loglik))

  negative_lr <- which(is.finite(lr) & lr < -sqrt(.Machine$double.eps))
  if (length(negative_lr)) {
    stop(.gl_inference_unavailable(list(
      status = "unavailable",
      failure_codes = "negative_likelihood_ratio",
      message = paste0(
        "Likelihood-ratio inference is unavailable because a fuller model has ",
        "a lower maximized log likelihood than its preceding model."
      ),
      rows = negative_lr,
      likelihood_ratio = lr[negative_lr]
    )))
  }
  lr[is.finite(lr) & lr < 0] <- 0

  p_value <- rep(NA_real_, length(models))
  reference_status <- c("baseline", rep("unverified", length(models) - 1L))
  reference_failure <- c(NA_character_, rep(NA_character_, length(models) - 1L))

  for (i in seq.int(2L, length(models))) {
    checks <- .gl_likelihood_reference_checks(models[[i - 1L]], models[[i]])
    failures <- checks$failure_codes
    if (!is.finite(loglik[[i]]) || !is.finite(loglik[[i - 1L]])) {
      failures <- c(failures, "nonfinite_loglik")
    }
    if (!is.finite(delta_df[[i]]) || delta_df[[i]] <= 0) {
      failures <- c(failures, "nonpositive_df_increment")
    }
    failures <- unique(failures)
    verified <- length(failures) == 0L && isTRUE(checks$sample_identical) &&
      isTRUE(checks$objective_identical) && isTRUE(checks$nested)
    if (verified) {
      reference_status[[i]] <- "reference_available"
      p_value[[i]] <- stats::pchisq(lr[[i]], df = delta_df[[i]], lower.tail = FALSE)
    } else {
      reference_status[[i]] <- if (any(grepl("mismatch|nonpositive|nonfinite", failures))) {
        "unavailable"
      } else {
        "unverified"
      }
      reference_failure[[i]] <- paste(failures, collapse = ";")
    }
  }

  out <- data.frame(
    model = labels,
    n_obs = as.integer(n_obs),
    df = as.numeric(df),
    logLik = as.numeric(loglik),
    AIC = as.numeric(aic),
    BIC = as.numeric(bic),
    delta_df = as.numeric(delta_df),
    LR_statistic = as.numeric(lr),
    p_value = as.numeric(p_value),
    reference_status = reference_status,
    reference_failure = reference_failure,
    stringsAsFactors = FALSE
  )

  class(out) <- c("gamlss_longitudinal_likelihood_compare", "data.frame")

  failures <- unique(reference_failure[!is.na(reference_failure) & nzchar(reference_failure)])
  contract <- .gl_inference_contract(
    "likelihood_ratio_nested",
    validity_status = if (all(reference_status[-1L] == "reference_available")) {
      "reference_available"
    } else if (any(reference_status[-1L] == "reference_available")) {
      "partially_available"
    } else {
      "unverified_or_unavailable"
    },
    failure_states = failures
  )
  contract$comparison_order <- labels
  contract$reference_status <- reference_status
  attr(out, "inference_contract") <- contract

  out
}

.gl_likelihood_sample_signature <- function(object) {
  response <- object$response
  if (is.null(response)) return(NULL)
  required <- c("response_subject", "response_margin")
  if (any(vapply(required, function(field) {
    is.null(object[[field]]) || length(object[[field]]) != length(response)
  }, logical(1)))) return(NULL)
  observed <- is.finite(response)
  components <- list(observed = observed, response = response[observed])
  for (field in c("response_subject", "response_margin", "response_time")) {
    value <- object[[field]]
    if (!is.null(value) && length(value) == length(response)) {
      components[[field]] <- value[observed]
    }
  }
  components
}

.gl_likelihood_objective_signature <- function(object) {
  margin <- object$margin_dist$family %||% NULL
  copula <- object$copula_dist %||% NULL
  if (is.null(margin) || is.null(copula)) return(NULL)
  list(
    margin = as.character(margin)[1L],
    copula = as.character(copula)[1L],
    include_dlcopdpar = object$include_dlcopdpar %||% NA
  )
}

.gl_likelihood_nested_design <- function(reduced, full) {
  reduced_x <- reduced$model_matrix$x %||% NULL
  full_x <- full$model_matrix$x %||% NULL
  if (is.null(reduced_x) || is.null(full_x)) return(NA)
  parameters <- union(names(reduced_x), names(full_x))
  fixed_nested <- vapply(parameters, function(parameter) {
    reduced_matrix <- reduced_x[[parameter]] %||% NULL
    full_matrix <- full_x[[parameter]] %||% NULL
    if (is.null(reduced_matrix) || ncol(reduced_matrix) == 0L) return(TRUE)
    if (is.null(full_matrix) || nrow(reduced_matrix) != nrow(full_matrix)) return(FALSE)
    reduced_cols <- colnames(reduced_matrix)
    full_cols <- colnames(full_matrix)
    if (is.null(reduced_cols) || is.null(full_cols) || !all(reduced_cols %in% full_cols)) {
      return(FALSE)
    }
    identical(
      unname(as.matrix(reduced_matrix[, reduced_cols, drop = FALSE])),
      unname(as.matrix(full_matrix[, reduced_cols, drop = FALSE]))
    )
  }, logical(1))
  reduced_s <- reduced$model_matrix$s %||% list()
  full_s <- full$model_matrix$s %||% list()
  isTRUE(all(fixed_nested)) && identical(reduced_s, full_s)
}

.gl_likelihood_reference_checks <- function(reduced, full) {
  reduced_sample <- .gl_likelihood_sample_signature(reduced)
  full_sample <- .gl_likelihood_sample_signature(full)
  sample_identical <- if (is.null(reduced_sample) || is.null(full_sample)) NA else {
    identical(reduced_sample, full_sample)
  }
  reduced_objective <- .gl_likelihood_objective_signature(reduced)
  full_objective <- .gl_likelihood_objective_signature(full)
  objective_identical <- if (is.null(reduced_objective) || is.null(full_objective)) NA else {
    identical(reduced_objective, full_objective)
  }
  nested <- .gl_likelihood_nested_design(reduced, full)
  failures <- character()
  if (isFALSE(sample_identical)) failures <- c(failures, "observed_sample_mismatch")
  if (is.na(sample_identical)) failures <- c(failures, "observed_sample_unverified")
  if (isFALSE(objective_identical)) failures <- c(failures, "likelihood_objective_mismatch")
  if (is.na(objective_identical)) failures <- c(failures, "likelihood_objective_unverified")
  if (isFALSE(nested)) failures <- c(failures, "model_nesting_not_verified")
  if (is.na(nested)) failures <- c(failures, "model_nesting_unverified")
  list(
    sample_identical = sample_identical,
    objective_identical = objective_identical,
    nested = nested,
    failure_codes = failures
  )
}
