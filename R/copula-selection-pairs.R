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
    time_intercepts = FALSE,
    margin_n_cycles = 100L,
    pit_seed = 1L,
    margin_args = list()) {
  upper_var <- NULL
  lower_var <- NULL
  log_mass_var <- NULL
  log_survival_var <- NULL
  screen_likelihood <- "copula_density"
  randomized <- FALSE
  used_pit_seed <- NA_integer_

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
    fitted_pit <- .select_copula_pit_from_fit(object, pit_seed = pit_seed)
    data <- data.frame(
      .subject = fitted_pit$subject,
      .time = fitted_pit$time,
      .u = fitted_pit$pit
    )
    if (!is.null(fitted_pit$upper)) data$.u_upper <- fitted_pit$upper
    if (!is.null(fitted_pit$lower)) data$.u_lower <- fitted_pit$lower
    if (!is.null(fitted_pit$log_mass)) data$.log_mass <- fitted_pit$log_mass
    if (!is.null(fitted_pit$log_survival)) data$.log_survival <- fitted_pit$log_survival
    subject_var <- ".subject"
    time_var <- ".time"
    u_var <- ".u"
    upper_var <- if (!is.null(fitted_pit$upper)) ".u_upper" else NULL
    lower_var <- if (!is.null(fitted_pit$lower)) ".u_lower" else NULL
    log_mass_var <- if (!is.null(fitted_pit$log_mass)) ".log_mass" else NULL
    log_survival_var <- if (!is.null(fitted_pit$log_survival)) ".log_survival" else NULL
    screen_likelihood <- fitted_pit$likelihood
    randomized <- fitted_pit$randomized
    used_pit_seed <- fitted_pit$pit_seed
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
        time_intercepts = time_intercepts,
        margin_n_cycles = margin_n_cycles,
        pit_seed = pit_seed,
        margin_args = margin_args
      )
      data <- auto$data
      u_var <- auto$u_var
      upper_var <- auto$upper_var
      lower_var <- auto$lower_var
      log_mass_var <- auto$log_mass_var
      log_survival_var <- auto$log_survival_var
      screen_likelihood <- auto$likelihood
      randomized <- auto$randomized
      used_pit_seed <- auto$pit_seed
    }
  }

  if (!all(c(subject_var, time_var, u_var) %in% names(data))) {
    stop("data must contain subject_var, time_var, and u_var columns.", call. = FALSE)
  }
  pairs <- .select_copula_adjacent_pairs(
    data,
    subject_var = subject_var,
    time_var = time_var,
    u_var = u_var,
    lags = lags,
    upper_var = upper_var,
    lower_var = lower_var,
    log_mass_var = log_mass_var,
    log_survival_var = log_survival_var
  )
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
  attr(pairs, "screen_likelihood") <- screen_likelihood
  attr(pairs, "randomized_pit") <- isTRUE(randomized)
  attr(pairs, "pit_seed") <- used_pit_seed
  pairs
}

.select_copula_adjacent_pairs <- function(
    data,
    subject_var,
    time_var,
    u_var,
    lags,
    upper_var = NULL,
    lower_var = NULL,
    log_mass_var = NULL,
    log_survival_var = NULL) {
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
      pair <- data.frame(
        u1 = .copula_clamp01(subject_data[[u_var]][left]),
        u2 = .copula_clamp01(subject_data[[u_var]][right]),
        copula_time = paste(subject_data[[time_var]][left], subject_data[[time_var]][right], sep = "->")
      )
      if (!is.null(upper_var) && !is.null(lower_var)) {
        pair$upper1 <- subject_data[[upper_var]][left]
        pair$upper2 <- subject_data[[upper_var]][right]
        pair$lower1 <- subject_data[[lower_var]][left]
        pair$lower2 <- subject_data[[lower_var]][right]
        if (!is.null(log_mass_var)) {
          pair$log_mass1 <- subject_data[[log_mass_var]][left]
          pair$log_mass2 <- subject_data[[log_mass_var]][right]
        }
        if (!is.null(log_survival_var)) {
          pair$log_survival1 <- subject_data[[log_survival_var]][left]
          pair$log_survival2 <- subject_data[[log_survival_var]][right]
        }
      }
      out[[k]] <- pair
    }
  }

  if (k == 0L) {
    empty <- data.frame(u1 = numeric(), u2 = numeric(), copula_time = character())
    if (!is.null(upper_var) && !is.null(lower_var)) {
      for (name in c("upper1", "upper2", "lower1", "lower2")) empty[[name]] <- numeric()
    }
    if (!is.null(log_mass_var)) {
      empty$log_mass1 <- numeric()
      empty$log_mass2 <- numeric()
    }
    if (!is.null(log_survival_var)) {
      empty$log_survival1 <- numeric()
      empty$log_survival2 <- numeric()
    }
    return(empty)
  }
  do.call(rbind, out[seq_len(k)])
}
