jss_actg_load_data <- function() {
  data_environment <- new.env(parent = emptyenv())
  utils::data(
    list = "actg_cd4",
    package = "gamlss.longitudinal",
    envir = data_environment
  )
  data <- get("actg_cd4", envir = data_environment, inherits = FALSE)
  data$week_factor <- factor(data$week, levels = sort(unique(data$week)))
  data$visit1 <- as.integer(data$week == min(data$week))
  data <- data[order(data$subject, data$week), , drop = FALSE]
  rownames(data) <- NULL
  data
}

jss_actg_capture <- function(expr) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- withCallingHandlers(
    expr,
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(
    value = value,
    warnings = unique(warnings),
    elapsed_sec = unname(proc.time()[["elapsed"]] - started)
  )
}

jss_actg_subject_patterns <- function(data) {
  rows <- lapply(split(data, data$subject), function(subject_data) {
    observed <- is.finite(subject_data$cd4)
    observed_index <- which(observed)
    intermittent <- length(observed_index) >= 2L &&
      any(!observed[min(observed_index):max(observed_index)])
    leading <- length(observed_index) && min(observed_index) > 1L
    terminal <- length(observed_index) && max(observed_index) < length(observed)
    segments <- if (any(observed)) sum(rle(observed)$values) else 0L
    data.frame(
      subject = subject_data$subject[[1L]],
      scheduled_rows = nrow(subject_data),
      observed_responses = sum(observed),
      observation_pattern = paste(as.integer(observed), collapse = ""),
      leading_missing = leading,
      intermittent_gap = intermittent,
      terminal_missing = terminal,
      observed_segments = segments,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_actg_sample_skewness <- function(x) {
  x <- x[is.finite(x)]
  centre <- x - mean(x)
  mean(centre^3) / mean(centre^2)^1.5
}

jss_actg_sample_excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  centre <- x - mean(x)
  mean(centre^4) / mean(centre^2)^2 - 3
}

jss_actg_pseudo_uniform <- function(x) {
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  out[keep] <- rank(x[keep], ties.method = "average") / (sum(keep) + 1)
  out
}

jss_actg_empirical_tail <- function(u1, u2, q = 0.05,
    tail = c("upper", "lower")) {
  tail <- match.arg(tail)
  keep <- is.finite(u1) & is.finite(u2)
  u1 <- u1[keep]
  u2 <- u2[keep]
  if (length(u1) < 2L) return(NA_real_)
  if (identical(tail, "upper")) {
    mean(u1 > 1 - q & u2 > 1 - q) / q
  } else {
    mean(u1 <= q & u2 <= q) / q
  }
}

jss_actg_data_summary <- function(data) {
  time_summary <- do.call(rbind, lapply(
    split(data, data$week_factor),
    function(part) {
      response <- part$cd4[is.finite(part$cd4)]
      data.frame(
        week = part$week[[1L]],
        n = length(response),
        n_zero = sum(response == 0),
        zero_prop = mean(response == 0),
        mean = mean(response),
        variance = stats::var(response),
        skewness = jss_actg_sample_skewness(response),
        excess_kurtosis = jss_actg_sample_excess_kurtosis(response),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(time_summary) <- NULL

  patterns <- jss_actg_subject_patterns(data)
  n_subjects <- nrow(patterns)
  n_observed <- sum(is.finite(data$cd4))
  cohort_summary <- data.frame(
    characteristic = c(
      "Participants", "Scheduled post-baseline observations",
      "Observed CD4 counts", "Missing CD4 counts",
      "Participants with all four CD4 counts",
      "Participants with at least one missing CD4 count",
      "Participants with a leading or intermittent gap",
      "Contiguous observed response segments",
      "Participants with leading unobserved visits",
      "Participants with an intermittent gap",
      "Participants with terminal missing observations"
    ),
    value = c(
      n_subjects, nrow(data), n_observed, sum(!is.finite(data$cd4)),
      sum(patterns$observed_responses == patterns$scheduled_rows),
      sum(patterns$observed_responses < patterns$scheduled_rows),
      sum(patterns$leading_missing | patterns$intermittent_gap),
      sum(patterns$observed_segments),
      sum(patterns$leading_missing), sum(patterns$intermittent_gap),
      sum(patterns$terminal_missing)
    ),
    denominator = c(n_subjects, nrow(data), nrow(data), nrow(data), rep(n_subjects, 7L)),
    stringsAsFactors = FALSE
  )

  wide <- reshape(
    data[c("subject", "week", "cd4")],
    idvar = "subject", timevar = "week", direction = "wide"
  )
  wide <- wide[, -1L, drop = FALSE]
  names(wide) <- sub("^cd4[.]", "week_", names(wide))
  pearson <- stats::cor(wide, use = "pairwise.complete.obs", method = "pearson")
  kendall <- stats::cor(wide, use = "pairwise.complete.obs", method = "kendall")
  pseudo <- as.data.frame(lapply(wide, jss_actg_pseudo_uniform))
  weeks <- names(pseudo)
  upper <- lower <- matrix(
    NA_real_, length(weeks), length(weeks), dimnames = list(weeks, weeks)
  )
  for (i in seq_along(weeks)) {
    for (j in seq_along(weeks)) {
      if (i == j) {
        upper[i, j] <- lower[i, j] <- 1
      } else {
        upper[i, j] <- jss_actg_empirical_tail(
          pseudo[[i]], pseudo[[j]], tail = "upper"
        )
        lower[i, j] <- jss_actg_empirical_tail(
          pseudo[[i]], pseudo[[j]], tail = "lower"
        )
      }
    }
  }

  missingness <- jss_actg_capture(check_missingness(
    data = data,
    time_var = "week_factor",
    subject_var = "subject",
    response_var = "cd4"
  ))
  list(
    time_summary = time_summary,
    cohort_summary = cohort_summary,
    subject_patterns = patterns,
    pearson = pearson,
    kendall = kendall,
    upper_tail = upper,
    lower_tail = lower,
    missingness = missingness$value,
    missingness_warnings = missingness$warnings,
    missingness_elapsed_sec = missingness$elapsed_sec
  )
}

jss_actg_optimizer_control <- function(
    max_elapsed_sec = 1190,
    margin_family = "ZIBNB",
    copula_family = "t",
    screening = TRUE) {
  margin_family <- toupper(margin_family)
  copula_family <- if (tolower(copula_family) == "t") "t" else toupper(copula_family)
  route <- paste(margin_family, copula_family, sep = "/")
  refinement <- switch(
    route,
    `ZISICHEL/G` = list(after = if (isTRUE(screening)) 2L else 10L,
                         iterations = 100L, refinements = 1L,
                         scale = 25, rescue = FALSE, repolish = 50L),
    `ZIBNB/G` = list(after = 2L, iterations = 170L, refinements = 1L,
                      scale = 25, rescue = FALSE, repolish = 50L),
    `ZISICHEL/t` = list(after = 1L, iterations = 60L, refinements = 1L,
                         scale = 5, rescue = FALSE, repolish = 50L),
    `ZIBNB/t` = list(after = 1L, iterations = 40L, refinements = 2L,
                      scale = 25, rescue = TRUE, repolish = 40L),
    list(after = 10L, iterations = 100L, refinements = 1L,
         scale = 25, rescue = FALSE, repolish = 50L)
  )
  gamlss_longitudinal_control(
    outer_tol = 1e-4,
    max_outer_iter = 50L,
    max_elapsed_sec = max_elapsed_sec,
    rs = list(
      inner_tol = 1e-5,
      max_inner_iter = 20L,
      initial_inner_iter = 2L,
      start_step_size = 2,
      warm_start_joint = TRUE,
      warm_start_joint_iter = 2L,
      use_backtracking = TRUE,
      backtracking_max_halves = 2L,
      adaptive_step_size = FALSE,
      update_lambda = TRUE,
      lambda_update_every = 5L,
      bfgs_refine = margin_family %in% c("ZISICHEL", "ZIBNB"),
      bfgs_refine_after = refinement$after,
      bfgs_max_refinements = refinement$refinements,
      bfgs_max_iter = refinement$iterations,
      bfgs_reltol = 1e-12,
      bfgs_parscale_max = refinement$scale,
      bfgs_coordinate_rescue = refinement$rescue,
      bfgs_repolish_iter = refinement$repolish,
      discrete_score_method = "analytical"
    )
  )
}

jss_actg_screen <- function(data, settings) {
  screen_margins <- c("ZIBNB", "ZIPIG", "ZISICHEL")
  screen_copulas <- c("t", "G")
  route_controls <- stats::setNames(
    lapply(
      as.vector(outer(screen_margins, screen_copulas, paste, sep = "/")),
      function(route) {
        fields <- strsplit(route, "/", fixed = TRUE)[[1L]]
        jss_actg_optimizer_control(
          margin_family = fields[[1L]], copula_family = fields[[2L]],
          screening = TRUE
        )
      }
    ),
    as.vector(outer(screen_margins, screen_copulas, paste, sep = "/"))
  )
  pooled <- jss_actg_capture(select_margin(
    data = data,
    time_var = "week_factor",
    subject_var = "subject",
    response_var = "cd4",
    verbose = 0
  ))
  time_varying <- jss_actg_capture(select_margin(
    data = data,
    time_var = "week_factor",
    subject_var = "subject",
    response_var = "cd4",
    time_intercepts = TRUE,
    verbose = 0
  ))
  copula <- jss_actg_capture(select_copula(
    data = data,
    time_var = "week_factor",
    subject_var = "subject",
    response_var = "cd4",
    margin_dist = gamlss.dist::ZIBNB(),
    pit_seed = 20260922L,
    verbose = 0
  ))
  checkpoint_dir <- file.path(settings$data_dir, "actg-joint-screen-checkpoints")
  joint <- jss_actg_capture(select_joint_distribution(
    data = data,
    time_var = "week_factor",
    subject_var = "subject",
    response_var = "cd4",
    margin_families = screen_margins,
    copula_families = screen_copulas,
    progress = FALSE,
    checkpoint_dir = checkpoint_dir,
    resume = TRUE,
    fit_args = list(
      missingness = "segment",
      method = "RS",
      optimization = "joint",
      optimizer_controls_by_route = route_controls,
      compute_vcov = FALSE,
      verbose = 0
    )
  ))
  list(
    pooled = pooled$value,
    time_varying = time_varying$value,
    copula = copula$value,
    joint = joint$value,
    warnings = list(
      pooled = pooled$warnings,
      time_varying = time_varying$warnings,
      copula = copula$warnings,
      joint = joint$warnings
    ),
    timings = c(
      pooled_margin = pooled$elapsed_sec,
      time_varying_margin = time_varying$elapsed_sec,
      copula = copula$elapsed_sec,
      joint = joint$elapsed_sec
    )
  )
}

jss_actg_load_screen_ledger <- function(settings) {
  input_dir <- file.path(settings$public_data_dir, "actg")
  read_screen <- function(name) {
    path <- file.path(input_dir, name)
    if (!file.exists(path)) {
      stop("Missing tracked ACTG screening input: ", path, call. = FALSE)
    }
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  list(
    pooled = read_screen("screen-margin-pooled.csv"),
    time_varying = read_screen("screen-margin-time.csv"),
    copula = read_screen("screen-copula.csv"),
    joint = read_screen("screen-joint.csv"),
    warnings = list(
      pooled = character(), time_varying = character(),
      copula = character(), joint = character()
    ),
    timings = c(
      pooled_margin = 0, time_varying_margin = 0,
      copula = 0, joint = 0
    ),
    source = "tracked_public_screen_ledger"
  )
}

jss_actg_fit_primary <- function(data) {
  fit <- jss_actg_capture(gamlss_longitudinal(
    data = data,
    time_var = "week_factor",
    subject_var = "subject",
    margin_dist = gamlss.dist::ZISICHEL(),
    copula_dist = "G",
    mu.formula = cd4 ~ treatment + age +
      s(baseline_cd4, bs = "ps", k = 8) +
      week_factor + baseline_cd4:visit1,
    sigma.formula = ~ treatment + baseline_cd4,
    nu.formula = ~ visit1,
    tau.formula = ~ treatment + sex + baseline_cd4,
    theta.formula = ~ treatment +
      s(baseline_cd4, bs = "ps", k = 8) + week_factor,
    zeta.formula = ~ 1,
    missingness = "segment",
    method = "RS",
    optimization = "joint",
    optimizer_control = jss_actg_optimizer_control(
      margin_family = "ZISICHEL", copula_family = "G", screening = FALSE
    ),
    verbose = 0
  ))
  if (!inherits(fit$value, "gamlss.longitudinal") ||
      !isTRUE(fit$value$convergence$converged)) {
    stop("The prespecified ACTG model did not converge.", call. = FALSE)
  }
  fit
}

jss_actg_cvm_uniform <- function(x) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  1 / (12 * n) + sum((x - (2 * seq_len(n) - 1) / (2 * n))^2)
}

jss_actg_worm_rms <- function(pit) {
  z <- sort(stats::qnorm(pmin(pmax(pit, .Machine$double.eps), 1 - .Machine$double.eps)))
  theoretical <- stats::qnorm((seq_along(z) - 0.5) / length(z))
  keep <- abs(theoretical) <= 2.5
  sqrt(mean((z[keep] - theoretical[keep])^2))
}

jss_actg_fit_calibration <- function(fit, seeds = 20260916L + 0:19) {
  observed <- is.finite(fit$dataset$response)
  times <- fit$dataset$time[observed]
  raw <- do.call(rbind, lapply(seeds, function(seed) {
    pit <- pithist(fit, randomize = TRUE, seed = seed, plot = FALSE)$pit
    if (length(pit) != length(times)) stop("ACTG PIT alignment failure.", call. = FALSE)
    do.call(rbind, lapply(sort(unique(times)), function(time) {
      value <- pit[times == time]
      data.frame(
        seed = seed,
        time = time,
        n = length(value),
        pit_sd = stats::sd(value),
        pit_cvm = jss_actg_cvm_uniform(value),
        worm_rms = jss_actg_worm_rms(value),
        stringsAsFactors = FALSE
      )
    }))
  }))
  aggregate(
    raw[c("pit_sd", "pit_cvm", "worm_rms")],
    by = list(time = raw$time), FUN = mean
  )
}

jss_actg_fit_candidate <- function(data, family, formulas) {
  fit <- jss_actg_capture(gamlss_longitudinal(
    data = data,
    time_var = "week_factor",
    subject_var = "subject",
    margin_dist = get(family, envir = asNamespace("gamlss.dist"))(),
    copula_dist = "G",
    mu.formula = formulas$mu,
    sigma.formula = formulas$sigma,
    nu.formula = formulas$nu,
    tau.formula = formulas$tau,
    theta.formula = formulas$theta,
    zeta.formula = formulas$zeta,
    missingness = "segment",
    method = "RS",
    optimization = "joint",
    optimizer_control = jss_actg_optimizer_control(
      margin_family = family, copula_family = "G", screening = FALSE
    ),
    compute_vcov = FALSE,
    verbose = 0
  ))
  if (!inherits(fit$value, "gamlss.longitudinal") ||
      !isTRUE(fit$value$convergence$converged)) {
    stop("ACTG family-comparison fit did not converge for ", family, ".", call. = FALSE)
  }
  fit
}

jss_actg_family_comparison <- function(data, primary_fit) {
  reference <- primary_fit$value
  formulas <- reference$formulas
  zibnb <- jss_actg_fit_candidate(data, "ZIBNB", formulas)
  zipig_formulas <- formulas
  zipig_formulas$nu <- formulas$tau
  zipig_formulas$tau <- ~ 1
  zipig <- jss_actg_fit_candidate(data, "ZIPIG", zipig_formulas)
  fits <- list(ZISICHEL = reference, ZIBNB = zibnb$value, ZIPIG = zipig$value)
  calibration <- lapply(fits, jss_actg_fit_calibration)
  table <- do.call(rbind, lapply(names(fits), function(family) {
    fit <- fits[[family]]
    diagnostics <- calibration[[family]]
    data.frame(
      margin = family,
      logLik = as.numeric(logLik(fit)),
      AIC = AIC(fit),
      BIC = BIC(fit),
      EDF = as.numeric(fit$model_selection["EDF", "joint"]),
      mean_pit_sd = mean(diagnostics$pit_sd),
      mean_pit_cvm = mean(diagnostics$pit_cvm),
      mean_worm_rms = mean(diagnostics$worm_rms),
      converged = isTRUE(fit$convergence$converged),
      stringsAsFactors = FALSE
    )
  }))
  table$delta_AIC <- table$AIC - min(table$AIC)
  table <- table[order(table$AIC), , drop = FALSE]
  rownames(table) <- NULL
  list(
    table = table,
    calibration = calibration,
    fits = fits,
    warnings = list(
      primary = primary_fit$warnings,
      ZIBNB = zibnb$warnings,
      ZIPIG = zipig$warnings
    ),
    timings = c(
      primary = primary_fit$elapsed_sec,
      ZIBNB = zibnb$elapsed_sec,
      ZIPIG = zipig$elapsed_sec
    )
  )
}

jss_actg_fit_sensitivity <- function(data) {
  patterns <- jss_actg_subject_patterns(data)
  excluded <- patterns$subject[patterns$intermittent_gap]
  sensitivity_data <- data[!data$subject %in% excluded, , drop = FALSE]
  fit <- jss_actg_capture(gamlss_longitudinal(
    data = sensitivity_data,
    time_var = "week_factor",
    subject_var = "subject",
    margin_dist = gamlss.dist::ZISICHEL(),
    copula_dist = "G",
    mu.formula = cd4 ~ treatment + age +
      s(baseline_cd4, bs = "ps", k = 8) +
      week_factor + baseline_cd4:visit1,
    sigma.formula = ~ treatment + baseline_cd4,
    nu.formula = ~ visit1,
    tau.formula = ~ treatment + sex + baseline_cd4,
    theta.formula = ~ treatment +
      s(baseline_cd4, bs = "ps", k = 8) + week_factor,
    zeta.formula = ~ 1,
    missingness = "segment",
    method = "RS",
    optimization = "joint",
    optimizer_control = jss_actg_optimizer_control(
      margin_family = "ZISICHEL", copula_family = "G", screening = FALSE
    ),
    compute_vcov = FALSE,
    verbose = 0
  ))
  if (!inherits(fit$value, "gamlss.longitudinal") ||
      !isTRUE(fit$value$convergence$converged)) {
    stop("The ACTG intermittent-gap sensitivity model did not converge.", call. = FALSE)
  }
  list(
    fit = fit$value,
    warnings = fit$warnings,
    elapsed_sec = fit$elapsed_sec,
    excluded_subjects = excluded,
    retained_subjects = length(unique(sensitivity_data$subject)),
    retained_observations = sum(is.finite(sensitivity_data$cd4))
  )
}

jss_actg_dependence_stats <- function(y, panel, week, tail_q = 0.05) {
  order_index <- order(match(panel$subject, unique(panel$subject)), panel$time)
  id <- panel$subject[order_index]
  time <- as.integer(panel$time[order_index])
  visit <- as.character(week[order_index])
  y <- y[order_index]
  left <- seq_len(length(y) - 1L)
  right <- left + 1L
  keep <- id[left] == id[right] & time[right] == time[left] + 1L &
    is.finite(y[left]) & is.finite(y[right])
  left <- left[keep]
  right <- right[keep]
  pair <- paste0(visit[left], "--", visit[right])
  do.call(rbind, lapply(split(seq_along(left), pair), function(index) {
    x <- y[left[index]]
    z <- y[right[index]]
    n <- length(x)
    u <- rank(x, ties.method = "average") / (n + 1)
    v <- rank(z, ties.method = "average") / (n + 1)
    data.frame(
      visit_pair = pair[index[[1L]]],
      metric = c("kendall_tau", "lower_tail_dependence", "upper_tail_dependence"),
      value = c(
        stats::cor(x, z, method = "kendall"),
        mean(u <= tail_q & v <= tail_q) / tail_q,
        mean(u >= 1 - tail_q & v >= 1 - tail_q) / tail_q
      ),
      stringsAsFactors = FALSE
    )
  }))
}

jss_actg_marginal_stats <- function(y, week) {
  do.call(rbind, lapply(split(seq_along(y), week), function(index) {
    value <- y[index]
    value <- value[is.finite(value)]
    data.frame(
      week = as.character(week[index[[1L]]]),
      metric = c("mean", "q10", "median", "q90", "zero_proportion"),
      value = c(
        mean(value), stats::quantile(value, 0.10, names = FALSE),
        stats::median(value), stats::quantile(value, 0.90, names = FALSE),
        mean(value == 0)
      ),
      stringsAsFactors = FALSE
    )
  }))
}

jss_actg_compare_simulations <- function(fit, simulations) {
  panel <- model.frame(fit, type = "expanded")
  if (nrow(panel) != nrow(simulations)) stop("ACTG simulation alignment failure.", call. = FALSE)
  week <- if ("time_covariate" %in% names(panel)) panel$time_covariate else panel$time
  values <- cbind(observed = panel$response, as.matrix(simulations))
  marginal_raw <- do.call(rbind, lapply(seq_len(ncol(values)), function(index) {
    transform(
      jss_actg_marginal_stats(values[, index], week),
      draw = colnames(values)[index]
    )
  }))
  dependence_raw <- do.call(rbind, lapply(seq_len(ncol(values)), function(index) {
    transform(
      jss_actg_dependence_stats(values[, index], panel, week),
      draw = colnames(values)[index]
    )
  }))
  summarise <- function(raw, keys) {
    groups <- split(
      seq_len(nrow(raw)),
      interaction(raw[keys], drop = TRUE, lex.order = TRUE)
    )
    out <- lapply(groups, function(index) {
      observed <- raw$value[index][raw$draw[index] == "observed"][[1L]]
      simulated <- raw$value[index][raw$draw[index] != "observed"]
      row <- raw[index[[1L]], keys, drop = FALSE]
      row$observed <- observed
      row$simulated_median <- stats::median(simulated, na.rm = TRUE)
      row$simulated_lower <- stats::quantile(
        simulated, 0.025, na.rm = TRUE, names = FALSE
      )
      row$simulated_upper <- stats::quantile(
        simulated, 0.975, na.rm = TRUE, names = FALSE
      )
      row
    })
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    out
  }
  list(
    marginal = summarise(marginal_raw, c("week", "metric")),
    dependence = summarise(dependence_raw, c("visit_pair", "metric"))
  )
}

jss_actg_postfit <- function(data, primary_fit) {
  fit <- primary_fit$value
  inference <- jss_actg_capture(summary(
    fit, include_vcov = TRUE, vcov_method = "sandwich"
  ))
  diagnostics <- jss_actg_capture(check_model(fit))
  predictions <- jss_actg_capture(publication_table(
    fit,
    table = "predictions",
    newdata = data,
    by = c("treatment", "week_factor"),
    probs = c(0.1, 0.5, 0.9),
    output = "data.frame"
  ))
  trajectories <- jss_actg_capture(simulate(
    fit, nsim = 500, seed = 20260922L
  ))
  comparison <- jss_actg_compare_simulations(fit, trajectories$value)

  observed_at_24 <- data$subject[data$week == 24 & is.finite(data$cd4)]
  future_week32 <- data[
    data$week == 32 & !is.finite(data$cd4) & data$subject %in% observed_at_24,
    , drop = FALSE
  ]
  continuation <- jss_actg_capture(simulate(
    fit,
    nsim = 500,
    seed = 20260922L,
    newdata = future_week32,
    condition_on = "fitted"
  ))
  continuation_matrix <- as.matrix(continuation$value)
  continuation_summary <- data.frame(
    subject = future_week32$subject,
    week = future_week32$week,
    simulated_mean = rowMeans(continuation_matrix, na.rm = TRUE),
    simulated_median = apply(
      continuation_matrix, 1L, stats::median, na.rm = TRUE
    ),
    simulated_q10 = apply(
      continuation_matrix, 1L, stats::quantile, probs = 0.1, na.rm = TRUE
    ),
    simulated_q90 = apply(
      continuation_matrix, 1L, stats::quantile, probs = 0.9, na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  list(
    summary = inference$value,
    check_model = diagnostics$value,
    predictions = as.data.frame(predictions$value, stringsAsFactors = FALSE),
    simulation_comparison = comparison,
    continuation_summary = continuation_summary,
    warnings = list(
      inference = inference$warnings,
      diagnostics = diagnostics$warnings,
      predictions = predictions$warnings,
      trajectories = trajectories$warnings,
      continuation = continuation$warnings
    ),
    timings = c(
      inference = inference$elapsed_sec,
      diagnostics = diagnostics$elapsed_sec,
      predictions = predictions$elapsed_sec,
      trajectories = trajectories$elapsed_sec,
      continuation = continuation$elapsed_sec
    )
  )
}

jss_actg_matrix_frame <- function(x, name) {
  out <- as.data.frame(as.table(x), stringsAsFactors = FALSE)
  names(out) <- c("week_1", "week_2", name)
  out
}

jss_actg_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

jss_actg_latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
  gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
}

jss_actg_fmt <- function(x, digits = 3L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "--")
}

jss_actg_write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  path
}

jss_actg_write_cohort_tex <- function(summary, path) {
  value <- summary$value
  denominator <- summary$denominator
  display <- ifelse(
    summary$characteristic %in% c("Observed CD4 counts", "Missing CD4 counts") |
      grepl("^Participants with", summary$characteristic),
    sprintf("%s (%.1f\\%%)", format(value, big.mark = ","), 100 * value / denominator),
    format(value, big.mark = ",")
  )
  main <- seq_len(8L)
  lines <- c(
    "% Generated by the public ACTG replication workflow.",
    "\\begin{table}[!h]", "\\centering", "\\small",
    "\\caption{ACTG analysis cohort and response availability.}",
    "\\label{tab:actg-cohort}",
    "\\begin{tabular}{@{}lr@{}}", "\\toprule",
    "Characteristic & Summary \\\\ ", "\\midrule",
    paste0(
      jss_actg_latex_escape(summary$characteristic[main]),
      " & ", display[main], " \\\\"
    ),
    "\\bottomrule", "\\end{tabular}", "\\begin{minipage}{0.92\\linewidth}",
    "\\footnotesize",
    sprintf(
      paste0(
        "A total of %s participants had one or more leading unobserved visits, ",
        "%s had an intermittent gap, and %s had terminal missing observations. ",
        "These categories can overlap."
      ),
      format(value[[9L]], big.mark = ","),
      format(value[[10L]], big.mark = ","),
      format(value[[11L]], big.mark = ",")
    ),
    "\\end{minipage}", "\\end{table}"
  )
  jss_actg_write_lines(lines, path)
}

jss_actg_treatment_labels <- function(x) {
  labels <- c(
    alternating_zidovudine_didanosine = "Alternating ZDV/ddI",
    zidovudine_zalcitabine = "ZDV + ddC",
    zidovudine_didanosine = "ZDV + ddI",
    zidovudine_didanosine_nevirapine = "ZDV + ddI + NVP"
  )
  out <- unname(labels[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

jss_actg_write_family_tex <- function(table, path) {
  rows <- apply(table, 1L, function(row) paste0(
    row[["margin"]], " & ", jss_actg_fmt(as.numeric(row[["AIC"]]), 2L), " & ",
    jss_actg_fmt(as.numeric(row[["delta_AIC"]]), 2L), " & ",
    jss_actg_fmt(as.numeric(row[["mean_pit_sd"]]), 3L), " & ",
    jss_actg_fmt(as.numeric(row[["mean_pit_cvm"]]), 3L), " & ",
    jss_actg_fmt(as.numeric(row[["mean_worm_rms"]]), 3L), " \\\\"
  ))
  jss_actg_write_lines(c(
    "% Generated by the public ACTG replication workflow.",
    "\\begin{table}[!h]", "\\centering", "\\small",
    paste0(
      "\\caption{Joint family comparison under the candidate Gumbel model ",
      "and common predictor structure.}"
    ),
    "\\label{tab:joint-family-screen}",
    "\\begin{tabular}{lrrrrr}", "\\toprule",
    "Margin & AIC & $\\Delta$AIC & Mean PIT SD & Mean PIT CvM & Mean worm RMS \\\\ ",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
    paste0(
      "\\begin{minipage}{\\linewidth}\\footnotesize Lower AIC, PIT ",
      "Cram{\\'e}r--von Mises statistic, and worm-plot RMS indicate better fit. ",
      "The target PIT standard deviation is $1/\\sqrt{12}=0.2887$.",
      "\\end{minipage}"
    ),
    "\\end{table}"
  ), path)
}

jss_actg_write_model_specification_tex <- function(root, path) {
  formulas <- utils::read.csv(
    file.path(
      root, "paper", "data", "public-derived", "actg",
      "final-formulas.csv"
    ),
    stringsAsFactors = FALSE
  )
  expected <- c("mu", "sigma", "nu", "tau", "theta")
  formulas <- formulas[match(expected, formulas$parameter), , drop = FALSE]
  if (anyNA(formulas$parameter)) {
    stop("The compact ACTG final-formula ledger is incomplete.", call. = FALSE)
  }
  feature <- c(
    "Mean (non-zero)", "Dispersion", "Shape",
    "Zero-inflation", "Adjacent-visit dependence"
  )
  link <- c(
    "$\\eta_{\\mu}=\\log(\\mu)$",
    "$\\eta_{\\sigma}=\\log(\\sigma)$",
    "$\\eta_{\\nu}=\\nu$",
    "$\\eta_{\\tau}=\\log\\!\\{\\tau/(1-\\tau)\\}$",
    "$\\eta_{\\theta}=\\log(\\theta-1)$"
  )
  symbol <- c("$\\mu$", "$\\sigma$", "$\\nu$", "$\\tau$", "$\\theta$")
  rows <- vapply(seq_along(expected), function(i) paste0(
    symbol[[i]], " & ", feature[[i]], " & ", link[[i]], " & ",
    jss_actg_abbreviate_formula(formulas$formula[[i]]), " \\\\"
  ), character(1L))
  jss_actg_write_lines(c(
    "% Generated from the compact public ACTG final-formula ledger.",
    "\\begin{table}[h]", "\\centering", "\\small",
    "\\caption{Final ZISICHEL--Gumbel model specification for the ACTG analysis.}",
    "\\label{tab:actg-final-model}",
    "\\begin{tabularx}{\\linewidth}{@{}c L{2.4cm} L{2.7cm} Y@{}}",
    "\\toprule",
    "Parameter & Modelled feature & Link & Linear predictor \\\\ ",
    "\\midrule", rows, "\\bottomrule", "\\end{tabularx}",
    paste0(
      "\\begin{minipage}{\\linewidth}\\footnotesize All predictors include an ",
      "intercept. The notation $s(\\cdot)$ denotes a penalized spline with ",
      "maximum basis dimension $k=8$. For the Gumbel copula, ",
      "$\\theta=1+\\exp(\\eta_{\\theta})$.\\end{minipage}"
    ),
    "\\end{table}"
  ), path)
}

jss_actg_write_treatment_tex <- function(coefficients, path) {
  keep <- coefficients$parameter %in% c("mu", "sigma", "tau", "theta") &
    grepl("treatment", coefficients$term, fixed = TRUE)
  x <- coefficients[keep, , drop = FALSE]
  treatment <- sub("^[^.]+[.]treatment", "", x$term)
  treatment <- jss_actg_treatment_labels(treatment)
  symbols <- c(
    mu = "$\\mu$", sigma = "$\\sigma$",
    tau = "$\\tau$", theta = "$\\theta$"
  )
  rows <- character(nrow(x))
  previous <- ""
  for (i in seq_len(nrow(x))) {
    parameter <- x$parameter[[i]]
    parameter_cell <- if (!identical(parameter, previous)) {
      symbols[[parameter]]
    } else {
      ""
    }
    ci <- x$estimate[[i]] + c(-1, 1) * stats::qnorm(0.975) *
      x$robust_std_error[[i]]
    p <- if (x$p_value[[i]] < 0.001) {
      "$<0.001$"
    } else {
      jss_actg_fmt(x$p_value[[i]], 3L)
    }
    rows[[i]] <- paste0(
      parameter_cell, " & ", treatment[[i]], " & ",
      jss_actg_fmt(x$estimate[[i]], 3L), " & ",
      jss_actg_fmt(x$model_std_error[[i]], 3L), " & ",
      jss_actg_fmt(x$robust_std_error[[i]], 3L),
      " & [", jss_actg_fmt(ci[[1L]], 3L), ", ",
      jss_actg_fmt(ci[[2L]], 3L), "] & ", p, " \\\\"
    )
    previous <- parameter
  }
  jss_actg_write_lines(c(
    "% Generated by the public ACTG replication workflow.",
    "\\begin{table}", "\\centering", "\\scriptsize",
    paste0(
      "\\caption{Treatment contrasts and standard-error sensitivity in the ",
      "final ZISICHEL--Gumbel model.}"
    ),
    "\\label{tab:actg-treatment-effects}",
    "\\begin{tabular}{clrrrrr}", "\\toprule",
    paste0(
      "Parameter & Treatment vs alternating ZDV/ddI & Estimate & Model SE & ",
      "Cluster SE & Cluster 95\\% CI & $p$-value \\\\ "
    ),
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
    paste0(
      "\\begin{minipage}{\\linewidth}\\footnotesize Estimates are on their ",
      "respective link scales. Cluster SEs use a participant-clustered sandwich ",
      "estimator; confidence intervals and $p$-values use the cluster SEs.",
      "\\end{minipage}"
    ),
    "\\end{table}"
  ), path)
}

jss_actg_write_prediction_tex <- function(predictions, path) {
  names_lower <- tolower(gsub("[^A-Za-z0-9]+", "_", names(predictions)))
  names(predictions) <- gsub("(^_|_$)", "", names_lower)
  treatment_column <- intersect("treatment", names(predictions))
  week_column <- intersect(c("week_factor", "week"), names(predictions))
  if (length(treatment_column) != 1L || length(week_column) != 1L) {
    stop("ACTG prediction grouping columns changed.", call. = FALSE)
  }
  treatment_column <- treatment_column[[1L]]
  week_column <- week_column[[1L]]
  predictions[[treatment_column]] <- jss_actg_treatment_labels(
    predictions[[treatment_column]]
  )
  required <- c(
    treatment_column, week_column, "n", "mean", "mu",
    "median", "q01", "q09"
  )
  if (!all(required %in% names(predictions))) {
    stop(
      "ACTG prediction table schema changed: ",
      paste(names(predictions), collapse = ", "),
      call. = FALSE
    )
  }
  rows <- vapply(seq_len(nrow(predictions)), function(i) paste0(
    jss_actg_latex_escape(predictions[[treatment_column]][[i]]), " & ",
    predictions[[week_column]][[i]], " & ", predictions$n[[i]], " & ",
    paste(
      jss_actg_fmt(
        as.numeric(predictions[i, c("mean", "mu", "median", "q01", "q09")]),
        3L
      ),
      collapse = " & "
    ),
    " \\\\"
  ), character(1L))
  jss_actg_write_lines(c(
    "% Generated by the public ACTG replication workflow.",
    "\\begin{table}[!ht]", "\\centering", "\\scriptsize",
    "\\caption{Fitted marginal CD4 summaries by treatment and follow-up week.}",
    "\\label{tab:actg-fitted-predictions}",
    "\\begin{tabular}{lrrrrrrr}", "\\toprule",
    paste0(
      "Treatment & Week & $n$ & Mean & $\\mu$ & Median & ",
      "$q_{0.10}$ & $q_{0.90}$ \\\\ "
    ),
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
    paste0(
      "\\begin{minipage}{\\linewidth}\\footnotesize Mean is the fitted marginal ",
      "response mean, while $\\mu$ is the fitted mean parameter of the Sichel ",
      "component before zero inflation.\\end{minipage}"
    ),
    "\\end{table}"
  ), path)
}

jss_actg_write_simulation_tex <- function(comparison, path) {
  marginal <- comparison$marginal[
    comparison$marginal$metric %in% c("mean", "zero_proportion"), ,
    drop = FALSE
  ]
  marginal <- marginal[
    order(as.numeric(marginal$week), match(
      marginal$metric, c("mean", "zero_proportion")
    )), ,
    drop = FALSE
  ]
  wide <- reshape(
    marginal, idvar = "week", timevar = "metric", direction = "wide"
  )
  marginal_rows <- vapply(seq_len(nrow(wide)), function(i) paste0(
    wide$week[[i]], " & ", jss_actg_fmt(wide$observed.mean[[i]], 1L), " & ",
    "\\makecell{", jss_actg_fmt(wide$simulated_median.mean[[i]], 1L),
    "\\\\{[", jss_actg_fmt(wide$simulated_lower.mean[[i]], 1L), ", ",
    jss_actg_fmt(wide$simulated_upper.mean[[i]], 1L), "]}} & ",
    jss_actg_fmt(100 * wide$observed.zero_proportion[[i]], 2L), "\\% & ",
    "\\makecell{",
    jss_actg_fmt(100 * wide$simulated_median.zero_proportion[[i]], 2L),
    "\\%\\\\{[",
    jss_actg_fmt(100 * wide$simulated_lower.zero_proportion[[i]], 2L), ", ",
    jss_actg_fmt(100 * wide$simulated_upper.zero_proportion[[i]], 2L),
    "]\\%}} \\\\"
  ), character(1L))
  dependence <- comparison$dependence
  dep_wide <- reshape(
    dependence, idvar = "visit_pair", timevar = "metric", direction = "wide"
  )
  dep_rows <- vapply(seq_len(nrow(dep_wide)), function(i) {
    cell <- function(metric) paste0(
      jss_actg_fmt(dep_wide[[paste0("observed.", metric)]][[i]], 3L), " & ",
      "\\makecell{",
      jss_actg_fmt(dep_wide[[paste0("simulated_median.", metric)]][[i]], 3L),
      "\\\\{[",
      jss_actg_fmt(dep_wide[[paste0("simulated_lower.", metric)]][[i]], 3L),
      ", ",
      jss_actg_fmt(dep_wide[[paste0("simulated_upper.", metric)]][[i]], 3L),
      "]}}"
    )
    paste0(
      dep_wide$visit_pair[[i]], " & ", cell("kendall_tau"), " & ",
      cell("lower_tail_dependence"), " & ",
      cell("upper_tail_dependence"), " \\\\"
    )
  }, character(1L))
  jss_actg_write_lines(c(
    "% Generated by the public ACTG replication workflow.",
    "\\begin{table}[!ht]", "\\centering", "\\scriptsize",
    "\\caption{Observed and model-simulated marginal and dependence summaries.}",
    "\\label{tab:actg-simulation-comparison}",
    "\\textit{Panel A: Selected marginal summaries}",
    "\\begin{tabular}{@{}lcccc@{}}", "\\toprule",
    paste0(
      "& \\multicolumn{2}{c}{Mean CD4} & ",
      "\\multicolumn{2}{c}{Zero proportion} \\\\ "
    ),
    "Week & Observed & Simulated & Observed & Simulated \\\\ ",
    "\\midrule", marginal_rows, "\\bottomrule", "\\end{tabular}",
    "\\vspace{0.8em}", "\\textit{Panel B: Adjacent-visit dependence}",
    "\\begin{tabular}{@{}lcccccc@{}}", "\\toprule",
    paste0(
      "& \\multicolumn{2}{c}{Kendall's $\\tau$} & ",
      "\\multicolumn{2}{c}{$\\widehat{\\lambda}_{L}(0.05)$} & ",
      "\\multicolumn{2}{c}{$\\widehat{\\lambda}_{U}(0.05)$} \\\\ "
    ),
    paste0(
      "Visit pair & Observed & Simulated & Observed & Simulated & ",
      "Observed & Simulated \\\\ "
    ),
    "\\midrule", dep_rows, "\\bottomrule", "\\end{tabular}",
    paste0(
      "\\begin{minipage}{\\linewidth}\\footnotesize Simulated entries show the ",
      "median and 95\\% simulation interval across 500 fitted joint trajectories. ",
      "The tail quantities are finite-threshold co-occurrence ratios.",
      "\\end{minipage}"
    ),
    "\\end{table}"
  ), path)
}

jss_actg_abbreviate_formula <- function(x) {
  x <- trimws(sub("^.*~", "", x))
  x <- gsub(
    's\\(baseline_cd4, bs = "ps", k = 8\\)', "s(base)", x
  )
  x <- gsub('s\\(age, bs = "ps", k = 8\\)', "s(age)", x)
  x <- gsub("baseline_cd4", "base", x, fixed = TRUE)
  x <- gsub("week_factor", "week", x, fixed = TRUE)
  x <- gsub("treatment", "trt", x, fixed = TRUE)
  x <- gsub("visit1", "v1", x, fixed = TRUE)
  paste0("\\texttt{", jss_actg_latex_escape(x), "}")
}

jss_actg_write_selection_tex <- function(root, path) {
  input_dir <- file.path(root, "paper", "data", "public-derived", "actg")
  accepted <- utils::read.csv(
    file.path(input_dir, "covariate-selection-accepted.csv"),
    stringsAsFactors = FALSE
  )
  targeted <- utils::read.csv(
    file.path(input_dir, "covariate-selection-targeted.csv"),
    stringsAsFactors = FALSE
  )
  final <- utils::read.csv(
    file.path(input_dir, "final-formulas.csv"),
    stringsAsFactors = FALSE
  )
  stage_order <- unique(accepted$stage)
  stage_labels <- c(
    "Initial $\\mu$", "Initial $\\sigma$", "Initial $\\nu$", "Initial $\\theta$",
    "Cycle 1 $\\mu$", "Cycle 1 $\\sigma$", "Cycle 1 $\\nu$", "Cycle 1 $\\theta$",
    "Cycle 2 $\\mu$", "Cycle 2 $\\sigma$", "Cycle 2 $\\nu$", "Cycle 2 $\\theta$"
  )
  names(stage_labels) <- stage_order
  accepted_rows <- character()
  previous <- ""
  for (i in seq_len(nrow(accepted))) {
    stage <- accepted$stage[[i]]
    stage_cell <- if (!identical(stage, previous)) stage_labels[[stage]] else ""
    change <- accepted$dropped_term[[i]]
    change <- if (!nzchar(change)) {
      "Starting model"
    } else {
      paste0(
        if (identical(accepted$change_type[[i]], "simplify")) {
          "Simplify "
        } else {
          "Drop "
        },
        "\\texttt{", jss_actg_latex_escape(change), "}"
      )
    }
    accepted_rows <- c(accepted_rows, paste0(
      stage_cell, " & ", accepted$step[[i]], " & ", change, " & ",
      jss_actg_abbreviate_formula(accepted$active_formula[[i]]), " & ",
      formatC(accepted$AIC[[i]], format = "f", digits = 2L, big.mark = ","),
      " & ",
      formatC(accepted$BIC[[i]], format = "f", digits = 2L, big.mark = ","),
      " & -- \\\\"
    ))
    previous <- stage
  }
  candidate_labels <- c(
    "preserved-model" = "No visit-1 addition",
    "mu_baseline_visit1" = "$\\mu$: base $\\times$ v1",
    "sigma_visit1" = "$\\sigma$: v1",
    "sigma_baseline_visit1" = "$\\sigma$: v1 + base $\\times$ v1",
    "nu_visit1" = "$\\nu$: v1",
    "tau_visit1" = "$\\tau$: v1",
    "tau_baseline_visit1" = "$\\tau$: v1 + base $\\times$ v1",
    "mu_baseline_visit1+nu_visit1" =
      "$\\mu$: base $\\times$ v1; $\\nu$: v1",
    "mu_baseline_visit1+sigma_visit1" =
      "$\\mu$: base $\\times$ v1; $\\sigma$: v1",
    "mu_baseline_visit1+sigma_baseline_visit1" =
      "$\\mu$: base $\\times$ v1; $\\sigma$: v1 + base $\\times$ v1",
    "mu_baseline_visit1+tau_visit1" =
      "$\\mu$: base $\\times$ v1; $\\tau$: v1",
    "mu_baseline_visit1+tau_baseline_visit1" =
      "$\\mu$: base $\\times$ v1; $\\tau$: v1 + base $\\times$ v1",
    "nu_week_factor_baseline_visit1" =
      "$\\nu$: week + base + base $\\times$ v1",
    "nu_visit1_baseline" = "$\\nu$: base + v1 + base $\\times$ v1",
    "nu_week_factor" = "$\\nu$: week",
    "nu_week_factor_baseline" = "$\\nu$: week + base"
  )
  target_labels <- c(
    "single-addition" = "T1: one addition",
    "second-addition" = "T2: second addition",
    "nu-refinement" = "T3: $\\nu$ refinement"
  )
  targeted <- targeted[
    order(match(targeted$stage, names(target_labels)), targeted$AIC), ,
    drop = FALSE
  ]
  targeted$rank <- ave(
    targeted$AIC, targeted$stage,
    FUN = function(x) rank(x, ties.method = "first")
  )
  target_rows <- character()
  previous <- ""
  for (i in seq_len(nrow(targeted))) {
    stage <- targeted$stage[[i]]
    label <- candidate_labels[[targeted$label[[i]]]]
    if (isTRUE(targeted$selected[[i]])) {
      label <- paste0(label, "$^{\\dagger}$")
    }
    if (identical(
      targeted$label[[i]], "mu_baseline_visit1+nu_visit1"
    )) {
      label <- paste0("\\textbf{", label, "$^{*}$}")
    }
    target_rows <- c(target_rows, paste0(
      if (!identical(stage, previous)) target_labels[[stage]] else "",
      " & ", targeted$rank[[i]], " & ", label, " & -- & ",
      formatC(targeted$AIC[[i]], format = "f", digits = 2L, big.mark = ","),
      " & ",
      formatC(targeted$BIC[[i]], format = "f", digits = 2L, big.mark = ","),
      " & ",
      formatC(targeted$time1_pit_cvm[[i]], format = "f", digits = 3L),
      " \\\\"
    ))
    previous <- stage
  }
  final_symbols <- c(
    mu = "$\\mu$", sigma = "$\\sigma$", nu = "$\\nu$",
    tau = "$\\tau$", theta = "$\\theta$"
  )
  selected <- targeted[
    targeted$label == "mu_baseline_visit1+nu_visit1", ,
    drop = FALSE
  ][1L, , drop = FALSE]
  final_rows <- vapply(seq_len(nrow(final)), function(i) paste0(
    if (i == 1L) "Final retained" else "",
    " & -- & ", final_symbols[[final$parameter[[i]]]], " & ",
    jss_actg_abbreviate_formula(final$formula[[i]]), " & ",
    if (i == 1L) {
      formatC(selected$AIC, format = "f", digits = 2L, big.mark = ",")
    } else {
      "--"
    },
    " & ",
    if (i == 1L) {
      formatC(selected$BIC, format = "f", digits = 2L, big.mark = ",")
    } else {
      "--"
    },
    " & ",
    if (i == 1L) {
      formatC(selected$time1_pit_cvm, format = "f", digits = 3L)
    } else {
      "--"
    },
    " \\\\"
  ), character(1L))
  jss_actg_write_lines(c(
    "% Generated from the compact public ACTG selection ledger.",
    "\\begingroup", "\\scriptsize", "\\setlength{\\tabcolsep}{2.6pt}",
    paste0(
      "\\begin{longtable}{@{}p{1.75cm}p{0.90cm}p{3.15cm}",
      "p{4.70cm}rrr@{}}"
    ),
    paste0(
      "\\caption{ACTG covariate-selection path for the retained ",
      "ZISICHEL--Gumbel model.}",
      "\\label{tab:actg-covariate-selection}\\\\"
    ),
    "\\toprule",
    paste0(
      "Stage & No. & Accepted change or candidate & Active predictor & ",
      "AIC & BIC & $W_1$ \\\\ "
    ),
    "\\midrule", "\\endfirsthead", "\\toprule",
    paste0(
      "Stage & No. & Accepted change or candidate & Active predictor & ",
      "AIC & BIC & $W_1$ \\\\ "
    ),
    "\\midrule", "\\endhead",
    paste0(
      "\\multicolumn{7}{@{}l}{\\textit{Panel A: accepted greedy backward ",
      "and cyclic-refinement steps}} \\\\ "
    ),
    accepted_rows, "\\midrule",
    paste0(
      "\\multicolumn{7}{@{}l}{\\textit{Panel B: targeted ZISICHEL ",
      "visit-1 diagnostic refinement}} \\\\ "
    ),
    target_rows, "\\midrule",
    paste0(
      "\\multicolumn{7}{@{}l}{\\textit{Panel C: final retained predictors}}",
      " \\\\ "
    ),
    final_rows, "\\bottomrule", "\\end{longtable}",
    paste0(
      "\\footnotesize The ledger records the accepted AIC-guided search and ",
      "targeted visit-1 refinement. $W_1$ is the visit-1 randomized-PIT ",
      "Cram{\\'e}r--von Mises statistic. $^{\\dagger}$ marks the minimum-AIC ",
      "candidate within a targeted stage and $^{*}$ the retained specification."
    ),
    "\\endgroup"
  ), path)
}

jss_actg_coef_vector <- function(fit) {
  summary(fit, include_vcov = FALSE)$coefficients[
    c("term", "estimate", "parameter")
  ]
}

jss_actg_sensitivity_comparison <- function(primary_fit, sensitivity) {
  primary <- jss_actg_coef_vector(primary_fit$value)
  secondary <- jss_actg_coef_vector(sensitivity$fit)
  names(primary)[names(primary) == "estimate"] <- "primary"
  names(secondary)[names(secondary) == "estimate"] <- "sensitivity"
  comparison <- merge(
    primary, secondary,
    by = c("term", "parameter"), all = TRUE, sort = FALSE
  )
  comparison$same_direction <- sign(comparison$primary) ==
    sign(comparison$sensitivity)
  comparison$relative_change_percent <- 100 *
    (comparison$sensitivity - comparison$primary) / abs(comparison$primary)
  comparison
}

jss_actg_save_plots <- function(data, primary_fit, sensitivity, settings) {
  fit <- primary_fit$value
  plot_warnings <- character()
  capture_plot <- function(expr) {
    result <- jss_actg_capture(expr)
    plot_warnings <<- c(plot_warnings, result$warnings)
    result$value
  }
  paths <- file.path(settings$figures_dir, c(
    "final-model-plot_dist.png",
    "nbi_plot_adjusted.png",
    "zisichel_plot_adjusted.png",
    "gumbel_copula_fit.png",
    "t_copula_fit.png",
    "plot_margin_trt.png",
    "plot_copula_trt.png",
    "final_model_copula_diagnostics.png",
    "final_model_plot.png",
    "plot_terms.png",
    "ACTG no intermittent missingness model.png"
  ))
  names(paths) <- c(
    "dist", "nbi", "zisichel", "gumbel", "t", "margin_trt",
    "copula_trt", "copula_diag", "model_diag", "terms", "sensitivity"
  )

  dist <- capture_plot(plot_dist(
    data,
    margin_dist = gamlss.dist::ZISICHEL(),
    time_var = "week_factor",
    subject_var = "subject",
    response_var = "cd4"
  ))
  ggplot2::ggsave(
    paths[["dist"]], dist,
    width = 12, height = 11, dpi = 240, bg = "white"
  )

  margin_plot <- function(family) {
    capture_plot(plot_margin_fit(
      data = data,
      time_var = "week_factor",
      response_var = "cd4",
      margin_dist = family,
      time_intercepts = TRUE,
      by_time = TRUE,
      bins = 100,
      separate_zero = TRUE
    ))$plot +
      ggplot2::coord_cartesian(xlim = c(0, 150)) +
      ggplot2::facet_wrap(~split_group, nrow = 1)
  }
  ggplot2::ggsave(
    paths[["nbi"]], margin_plot(gamlss.dist::NBI()),
    width = 12, height = 4, dpi = 240, bg = "white"
  )
  ggplot2::ggsave(
    paths[["zisichel"]], margin_plot(gamlss.dist::ZISICHEL()),
    width = 12, height = 4, dpi = 240, bg = "white"
  )

  copula_plot <- function(family) {
    capture_plot(plot_copula_fit(
      data = data,
      time_var = "week_factor",
      response_var = "cd4",
      margin_dist = gamlss.dist::ZISICHEL(),
      copula_dist = family,
      by_time = TRUE
    ))$plot
  }
  ggplot2::ggsave(
    paths[["gumbel"]], copula_plot("G"),
    width = 12, height = 4, dpi = 240, bg = "white"
  )
  ggplot2::ggsave(
    paths[["t"]], copula_plot("t"),
    width = 12, height = 4, dpi = 240, bg = "white"
  )

  margin_trt <- capture_plot(plot_margin_fit(
    fit = fit, data = data, by = "treatment", bins = 100
  ))$plot +
    ggplot2::coord_cartesian(xlim = c(0, 150), ylim = c(0, 0.3)) +
    ggplot2::facet_wrap(~split_group, nrow = 1)
  ggplot2::ggsave(
    paths[["margin_trt"]], margin_trt,
    width = 12, height = 4, dpi = 240, bg = "white"
  )
  copula_trt <- capture_plot(plot_copula_fit(
    fit = fit, data = data, by = "treatment"
  ))$plot +
    ggplot2::facet_wrap(~split_group, nrow = 1)
  ggplot2::ggsave(
    paths[["copula_trt"]], copula_trt,
    width = 12, height = 4, dpi = 240, bg = "white"
  )

  diagnostics <- capture_plot(plot(fit))
  copula_diagnostics <- capture_plot(
    plot_copula_diagnostics(fit, type = "all", ncol = 3)
  )
  terms <- capture_plot(plot_terms(fit))
  ggplot2::ggsave(
    paths[["model_diag"]], diagnostics$dashboard,
    width = 12, height = 9, dpi = 240, bg = "white"
  )
  ggplot2::ggsave(
    paths[["copula_diag"]], copula_diagnostics$dashboard,
    width = 12, height = 12, dpi = 240, bg = "white"
  )
  term_dashboard <- ggpubr::ggarrange(
    plotlist = terms$dashboard$plotlist,
    ncol = terms$dashboard$ncol,
    nrow = terms$dashboard$nrow
  )
  ggplot2::ggsave(
    paths[["terms"]], term_dashboard,
    width = 12, height = 10, dpi = 240, bg = "white"
  )

  comparison <- jss_actg_sensitivity_comparison(primary_fit, sensitivity)
  plot_data <- rbind(
    data.frame(
      parameter = comparison$parameter,
      term = comparison$term,
      analysis = "Primary",
      estimate = comparison$primary
    ),
    data.frame(
      parameter = comparison$parameter,
      term = comparison$term,
      analysis = "No intermittent gaps",
      estimate = comparison$sensitivity
    )
  )
  plot_data <- plot_data[
    !grepl("intercept$", plot_data$term, ignore.case = TRUE), ,
    drop = FALSE
  ]
  plot_data$label <- sub("^[^.]+[.]", "", plot_data$term)
  plot_data$label <- factor(
    plot_data$label, levels = rev(unique(plot_data$label))
  )
  sensitivity_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = estimate, y = label, colour = analysis)
  ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", colour = "grey60"
    ) +
    ggplot2::geom_point(
      size = 2.2,
      position = ggplot2::position_dodge(width = 0.45)
    ) +
    ggplot2::facet_grid(
      parameter ~ ., scales = "free_y", space = "free_y"
    ) +
    ggplot2::labs(
      x = "Coefficient estimate (link scale)", y = NULL, colour = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank()
    )
  ggplot2::ggsave(
    paths[["sensitivity"]], sensitivity_plot,
    width = 10, height = 11, dpi = 240, bg = "white"
  )
  output_paths <- unname(paths)
  attr(output_paths, "warnings") <- unique(plot_warnings)
  output_paths
}

jss_run_05_actg <- function(settings, data, summaries, screening, primary_fit,
    family_comparison, sensitivity, postfit) {
  final_model_path <- file.path(settings$data_dir, "actg-final-model.rds")
  saveRDS(primary_fit$value, final_model_path, compress = "gzip")
  data_paths <- c(
    jss_actg_write_csv(
      summaries$time_summary,
      file.path(settings$data_dir, "actg-time-summary.csv")
    ),
    jss_actg_write_csv(
      summaries$cohort_summary,
      file.path(settings$data_dir, "actg-cohort-summary.csv")
    ),
    jss_actg_write_csv(
      summaries$subject_patterns,
      file.path(settings$data_dir, "actg-subject-patterns.csv")
    ),
    jss_actg_write_csv(
      jss_actg_matrix_frame(summaries$pearson, "pearson"),
      file.path(settings$data_dir, "actg-pearson-correlation.csv")
    ),
    jss_actg_write_csv(
      jss_actg_matrix_frame(summaries$kendall, "kendall_tau"),
      file.path(settings$data_dir, "actg-kendall-correlation.csv")
    ),
    jss_actg_write_csv(
      jss_actg_matrix_frame(summaries$upper_tail, "upper_tail"),
      file.path(settings$data_dir, "actg-upper-tail.csv")
    ),
    jss_actg_write_csv(
      jss_actg_matrix_frame(summaries$lower_tail, "lower_tail"),
      file.path(settings$data_dir, "actg-lower-tail.csv")
    ),
    jss_actg_write_csv(
      as.data.frame(screening$pooled),
      file.path(settings$data_dir, "actg-margin-screen-pooled.csv")
    ),
    jss_actg_write_csv(
      as.data.frame(screening$time_varying),
      file.path(settings$data_dir, "actg-margin-screen-time.csv")
    ),
    jss_actg_write_csv(
      as.data.frame(screening$copula),
      file.path(settings$data_dir, "actg-copula-screen.csv")
    ),
    jss_actg_write_csv(
      as.data.frame(screening$joint),
      file.path(settings$data_dir, "actg-joint-screen.csv")
    ),
    jss_actg_write_csv(
      family_comparison$table,
      file.path(settings$data_dir, "actg-final-family-comparison.csv")
    ),
    jss_actg_write_csv(
      postfit$summary$coefficients,
      file.path(settings$data_dir, "actg-coefficients.csv")
    ),
    jss_actg_write_csv(
      postfit$predictions,
      file.path(settings$data_dir, "actg-predictions.csv")
    ),
    jss_actg_write_csv(
      postfit$simulation_comparison$marginal,
      file.path(settings$data_dir, "actg-simulation-marginal.csv")
    ),
    jss_actg_write_csv(
      postfit$simulation_comparison$dependence,
      file.path(settings$data_dir, "actg-simulation-dependence.csv")
    ),
    jss_actg_write_csv(
      postfit$continuation_summary,
      file.path(settings$data_dir, "actg-continuation-summary.csv")
    ),
    jss_actg_write_csv(
      jss_actg_sensitivity_comparison(primary_fit, sensitivity),
      file.path(settings$data_dir, "actg-intermittent-gap-sensitivity.csv")
    )
  )
  data_paths <- c(
    data_paths,
    jss_actg_write_lines(
      capture.output(print(summaries$missingness)),
      file.path(settings$data_dir, "actg-missingness-output.txt")
    ),
    jss_actg_write_lines(
      capture.output(print(postfit$summary)),
      file.path(settings$data_dir, "actg-model-summary.txt")
    ),
    jss_actg_write_lines(
      capture.output(print(postfit$check_model)),
      file.path(settings$data_dir, "actg-check-model.txt")
    )
  )

  figure_paths <- jss_actg_save_plots(
    data, primary_fit, sensitivity, settings
  )
  family_warnings <- family_comparison$warnings
  family_warnings$primary <- NULL
  warning_groups <- c(
    list(
      primary_fit = primary_fit$warnings,
      sensitivity = sensitivity$warnings,
      plots = if (is.null(attr(figure_paths, "warnings"))) {
        character()
      } else {
        attr(figure_paths, "warnings")
      }
    ),
    screening$warnings,
    family_warnings,
    postfit$warnings
  )
  warning_table <- do.call(rbind, lapply(names(warning_groups), function(stage) {
    value <- warning_groups[[stage]]
    if (!length(value)) return(NULL)
    data.frame(stage = stage, warning = value, stringsAsFactors = FALSE)
  }))
  if (is.null(warning_table)) {
    warning_table <- data.frame(stage = character(), warning = character())
  }
  warning_code <- function(value) {
    if (grepl("segmented likelihood", value, fixed = TRUE)) {
      "segmented_likelihood_notice"
    } else if (grepl("Analytical Hessian vcov failed", value, fixed = TRUE)) {
      "analytical_hessian_fallback"
    } else if (grepl("fitted_gradient_too_large", value, fixed = TRUE)) {
      "fitted_gradient_review"
    } else {
      "unexpected"
    }
  }
  warning_table$warning_code <- vapply(
    warning_table$warning, warning_code, character(1L)
  )
  warning_table$expected <- warning_table$warning_code != "unexpected"
  if (any(!warning_table$expected)) {
    stop(
      "Unexpected warning(s) occurred in the ACTG workflow: ",
      paste(unique(warning_table$warning[!warning_table$expected]), collapse = " | "),
      call. = FALSE
    )
  }
  timings <- c(
    screening$timings,
    family_comparison$timings,
    sensitivity = sensitivity$elapsed_sec,
    postfit$timings
  )
  audit_paths <- c(
    jss_actg_write_csv(
      warning_table,
      file.path(settings$logs_dir, "actg-warnings.csv")
    ),
    jss_actg_write_csv(
      data.frame(stage = names(timings), elapsed_sec = unname(timings)),
      file.path(settings$logs_dir, "actg-timings.csv")
    )
  )

  table_paths <- c(
    jss_actg_write_cohort_tex(
      summaries$cohort_summary,
      file.path(settings$tables_dir, "actg-cohort.tex")
    ),
    jss_actg_write_family_tex(
      family_comparison$table,
      file.path(settings$tables_dir, "actg-joint-family-screen.tex")
    ),
    jss_actg_write_model_specification_tex(
      settings$root,
      file.path(settings$tables_dir, "actg-final-model.tex")
    ),
    jss_actg_write_treatment_tex(
      postfit$summary$coefficients,
      file.path(settings$tables_dir, "actg-treatment-effects.tex")
    ),
    jss_actg_write_prediction_tex(
      postfit$predictions,
      file.path(settings$tables_dir, "actg-fitted-predictions.tex")
    ),
    jss_actg_write_simulation_tex(
      postfit$simulation_comparison,
      file.path(settings$tables_dir, "actg-simulation-comparison.tex")
    ),
    jss_actg_write_selection_tex(
      settings$root,
      file.path(settings$tables_dir, "actg-covariate-selection.tex")
    )
  )
  list(
    module_id = "05-actg-workflow",
    title = "Reproducible ACTG workflow",
    status = "regenerated",
    data = c(final_model_path, data_paths, audit_paths),
    tables = table_paths,
    figures = figure_paths,
    notes = paste0(
      "Public package-shipped ACTG data; final model and intermittent-gap ",
      "sensitivity refitted."
    )
  )
}
