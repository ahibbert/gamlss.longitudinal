.gl_prediction_frame <- function(object, newdata = NULL, require_response = FALSE) {
  diag_data <- .gl_fitted_distribution(object, newdata = newdata, require_response = require_response)

  out <- data.frame(
    subject = diag_data$subject,
    time = diag_data$time,
    response = diag_data$response,
    stringsAsFactors = FALSE
  )

  for (par_name in names(diag_data$params)) {
    out[[par_name]] <- as.numeric(diag_data$params[[par_name]])
  }

  out
}

.gl_quantile_columns <- function(params, family, probs) {
  out <- list()

  for (prob in probs) {
    label <- paste0("q", gsub("[^0-9]+", "", format(prob, trim = TRUE, scientific = FALSE)))

    if (identical(label, "q")) {
      label <- paste0("q", seq_along(out) + 1L)
    }

    out[[label]] <- .gl_call_family_fun("q", family, prob, params)
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}

.gl_mu_values <- function(params) {
  if ("mu" %in% names(params)) {
    return(as.numeric(params$mu))
  }

  as.numeric(params[[1L]])
}

.gl_distribution_mean <- function(params, family, probs = seq(0.001, 0.999, length.out = 199L)) {
  family <- as.character(family)[1L]

  if (identical(family, "NO") && "mu" %in% names(params)) {
    return(as.numeric(params$mu))
  }

  if (identical(family, "GA") && "mu" %in% names(params)) {
    return(as.numeric(params$mu))
  }

  if (identical(family, "LOGNO") && all(c("mu", "sigma") %in% names(params))) {
    return(exp(as.numeric(params$mu) + 0.5 * as.numeric(params$sigma)^2))
  }

  probs <- as.numeric(probs)

  probs <- probs[is.finite(probs) & probs > 0 & probs < 1]

  if (length(probs) < 3L) {
    stop("'probs' must contain at least three finite probabilities inside (0, 1).", call. = FALSE)
  }

  q_mat <- vapply(probs, function(prob) {
    as.numeric(.gl_call_family_fun("q", family, prob, params))
  }, numeric(length(params[[1L]])))

  as.numeric(rowMeans(q_mat, na.rm = TRUE))
}

.gl_prediction_values <- function(type, params, family) {
  switch(type,
    response = .gl_mu_values(params),
    mu = .gl_mu_values(params),
    mean = .gl_distribution_mean(params, family),
    median = as.numeric(.gl_call_family_fun("q", family, 0.5, params)),
    stop("Unsupported prediction value type.", call. = FALSE)
  )
}

.gl_prediction_eval_values <- function(value, label, type, diag_data, require_response) {
  if (is.null(value)) {
    if (!require_response) {
      stop("'", label, "' is required for type = '", type, "' when no response is available.", call. = FALSE)
    }

    return(as.numeric(diag_data$response))
  }

  value <- as.numeric(value)

  n <- length(diag_data$subject)

  if (length(value) == 1L) {
    return(rep(value, n))
  }

  if (length(value) != n) {
    stop("'", label, "' must be length 1 or match the number of prediction rows.", call. = FALSE)
  }

  value
}
