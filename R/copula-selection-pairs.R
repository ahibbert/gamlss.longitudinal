.select_copula_pairs <- function(
    data,
    object,
    u1,
    u2,
    u,
    u_var,
    response_var,
    margin_dist,
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    subject_var,
    time_var,
    lags,
    time_intercepts = FALSE) {
  if (!is.null(u1) || !is.null(u2)) {
    if (is.null(u1) || is.null(u2)) {
      stop("Both u1 and u2 must be supplied for direct pseudo-observation pairs.", call. = FALSE)
    }
    n <- max(length(u1), length(u2))
    return(data.frame(
      u1 = rep(.copula_clamp01(u1), length.out = n),
      u2 = rep(.copula_clamp01(u2), length.out = n)
    ))
  }

  if (!is.null(object)) {
    u <- .select_copula_u_from_fit(object)
    data <- data.frame(
      .subject = object$response_subject,
      .time = object$response_margin,
      .u = u
    )
    subject_var <- ".subject"
    time_var <- ".time"
    u_var <- ".u"
  } else {
    if (is.null(data)) {
      stop("Provide either object, u1/u2, or data with u/u_var.", call. = FALSE)
    }
    data <- as.data.frame(data)
    margin_selection_time_var <- .select_copula_margin_selection_time_var(margin_dist)
    if (!is.null(margin_selection_time_var) && !time_var %in% names(data) && margin_selection_time_var %in% names(data)) {
      time_var <- margin_selection_time_var
    }
    if (!is.null(u)) {
      if (length(u) != nrow(data)) {
        stop("u must have one value per row of data.", call. = FALSE)
      }
      data[[".u"]] <- u
      u_var <- ".u"
    }
    if (is.null(u_var)) {
      auto <- .select_copula_auto_u(
        data = data,
        response_var = response_var,
        margin_dist = margin_dist,
        mu.formula = mu.formula,
        sigma.formula = sigma.formula,
        nu.formula = nu.formula,
        tau.formula = tau.formula,
        time_var = time_var,
        time_intercepts = time_intercepts
      )
      data <- auto$data
      u_var <- auto$u_var
    }
  }

  if (!all(c(subject_var, time_var, u_var) %in% names(data))) {
    stop("data must contain subject_var, time_var, and u_var columns.", call. = FALSE)
  }
  pairs <- .select_copula_adjacent_pairs(data, subject_var = subject_var, time_var = time_var, u_var = u_var, lags = lags)
  if (exists("auto", inherits = FALSE)) {
    attr(pairs, "margin_selection") <- auto$margin_selection
    attr(pairs, "pseudo_observation_source") <- auto$source
  } else if (!is.null(object)) {
    attr(pairs, "pseudo_observation_source") <- "fitted_object"
  } else if (!is.null(u) || identical(u_var, ".u")) {
    attr(pairs, "pseudo_observation_source") <- "u"
  } else {
    attr(pairs, "pseudo_observation_source") <- u_var
  }
  pairs
}



.select_copula_adjacent_pairs <- function(data, subject_var, time_var, u_var, lags) {
  ord <- order(data[[subject_var]], data[[time_var]])
  data <- data[ord, , drop = FALSE]
  subjects <- unique(data[[subject_var]])
  out <- vector("list", length(subjects) * length(lags))
  k <- 0L

  for (subject in subjects) {
    subject_data <- data[data[[subject_var]] == subject, , drop = FALSE]
    subject_data <- subject_data[order(subject_data[[time_var]]), , drop = FALSE]
    n_time <- nrow(subject_data)
    for (lag in lags) {
      if (n_time <= lag) next
      left <- seq_len(n_time - lag)
      right <- left + lag
      k <- k + 1L
      out[[k]] <- data.frame(
        u1 = .copula_clamp01(subject_data[[u_var]][left]),
        u2 = .copula_clamp01(subject_data[[u_var]][right]),
        copula_time = paste(subject_data[[time_var]][left], subject_data[[time_var]][right], sep = "->")
      )
    }
  }

  if (k == 0L) {
    return(data.frame(u1 = numeric(), u2 = numeric()))
  }
  do.call(rbind, out[seq_len(k)])
}
