.build_margin_eval_cache <- function(margin_dist, calc_d2 = FALSE) {
  if (calc_d2 == TRUE) {
    to_include <- grepl("dld", names(margin_dist)) | grepl("d2ld", names(margin_dist))
  } else {
    to_include <- grepl("dld", names(margin_dist))
  }

  margin_deriv_names <- names(margin_dist)[to_include]
  margin_deriv_cache <- lapply(margin_deriv_names, function(deriv_name) {
    FUN <- margin_dist[[deriv_name]]
    list(name = deriv_name, FUN = FUN, args = formalArgs(FUN))
  })

  margin_pFUN <- get(
    paste("p", margin_dist$family[1], sep = ""),
    envir = asNamespace("gamlss.dist"),
    mode = "function",
    inherits = FALSE
  )
  margin_dFUN <- get(
    paste("d", margin_dist$family[1], sep = ""),
    envir = asNamespace("gamlss.dist"),
    mode = "function",
    inherits = FALSE
  )

  list(
    calc_d2 = calc_d2,
    margin_deriv_cache = margin_deriv_cache,
    margin_pFUN = margin_pFUN,
    margin_p_args = formalArgs(margin_pFUN),
    margin_dFUN = margin_dFUN,
    margin_d_args = formalArgs(margin_dFUN),
    force_scalar_family = .gl_margin_requires_scalar_calls(margin_dist$family[1]),
    family_call_cache = new.env(parent = emptyenv())
  )
}

.gl_margin_requires_scalar_calls <- function(family) {
  as.character(family)[1L] %in% c("ZANBI", "ZAP", "ZAPIG", "ZABI")
}

.call_margin_family_cached <- function(
    FUN,
    args,
    FUN_args,
    cacheable = FALSE,
    cache_env = NULL,
    cache_prefix = "",
    force_scalar = FALSE) {
  call_args <- args[FUN_args]
  safe_call <- function(call_args, fallback_n = NULL) {
    if (is.null(fallback_n)) {
      fallback_n <- max(c(1L, vapply(call_args, length, integer(1))), na.rm = TRUE)
    }
    if (isTRUE(force_scalar) && fallback_n > 1L) {
      lengths <- vapply(call_args, length, integer(1))
      if (!all(lengths %in% c(1L, fallback_n))) return(rep(NA_real_, fallback_n))
      expanded <- lapply(call_args, rep, length.out = fallback_n)
      return(vapply(seq_len(fallback_n), function(ii) {
        scalar_args <- lapply(expanded, `[`, ii)
        value <- tryCatch(do.call(FUN, args = scalar_args), error = function(e) NA_real_)
        as.numeric(value)[1L]
      }, numeric(1)))
    }
    tryCatch(
      do.call(FUN, args = call_args),
      error = function(e) rep(NA_real_, fallback_n)
    )
  }
  if (!isTRUE(cacheable) || length(call_args) == 0L) {
    return(safe_call(call_args))
  }

  arg_lengths <- vapply(call_args, length, integer(1))
  n <- max(arg_lengths)
  if (n <= 1L) {
    return(safe_call(call_args, fallback_n = n))
  }

  if (any(!arg_lengths %in% c(1L, n))) {
    return(safe_call(call_args, fallback_n = n))
  }
  if (!all(vapply(call_args, function(x) is.atomic(x) && !is.character(x), logical(1)))) {
    return(safe_call(call_args, fallback_n = n))
  }

  expanded <- lapply(call_args, rep, length.out = n)
  key_parts <- lapply(expanded, function(x) {
    if (is.numeric(x) || is.integer(x)) {
      format(signif(as.numeric(x), 15), scientific = FALSE, trim = TRUE)
    } else {
      as.character(x)
    }
  })
  key <- do.call(paste, c(key_parts, sep = "\r"))
  unique_idx <- !duplicated(key)
  if (sum(unique_idx) > 0.8 * n) {
    return(safe_call(call_args, fallback_n = n))
  }

  if (!is.null(cache_env)) {
    cached <- rep(NA_real_, n)
    hit <- logical(n)
    lookup_key <- paste0(cache_prefix, "\r", key)
    for (ii in seq_len(n)) {
      if (exists(lookup_key[[ii]], envir = cache_env, inherits = FALSE)) {
        cached[[ii]] <- get(lookup_key[[ii]], envir = cache_env, inherits = FALSE)
        hit[[ii]] <- TRUE
      }
    }
    if (all(hit)) {
      return(cached)
    }
    need <- !hit
    need_unique <- !duplicated(key[need])
    need_args <- lapply(expanded, function(x) x[need][need_unique])
    names(need_args) <- names(call_args)
    need_val <- safe_call(need_args, fallback_n = sum(need_unique))
    if (length(need_val) != sum(need_unique)) {
      return(safe_call(call_args, fallback_n = n))
    }
    need_lookup <- paste0(cache_prefix, "\r", key[need][need_unique])
    for (jj in seq_along(need_lookup)) {
      assign(need_lookup[[jj]], need_val[[jj]], envir = cache_env)
    }
    cached[need] <- need_val[match(key[need], key[need][need_unique])]
    return(cached)
  }

  unique_args <- lapply(expanded, `[`, unique_idx)
  names(unique_args) <- names(call_args)
  unique_val <- safe_call(unique_args, fallback_n = sum(unique_idx))
  if (length(unique_val) != sum(unique_idx)) {
    return(safe_call(call_args, fallback_n = n))
  }
  unique_val[match(key, key[unique_idx])]
}

.fast_delaporte_cdf <- function(q, mu, sigma, nu, max_terms = 5e5) {
  n <- max(length(q), length(mu), length(sigma), length(nu))
  q <- rep(q, length.out = n)
  mu <- rep(mu, length.out = n)
  sigma <- rep(sigma, length.out = n)
  nu <- rep(nu, length.out = n)
  if (any(!is.finite(q)) || any(q != floor(q)) ||
      any(!is.finite(mu) | mu <= 0) ||
      any(!is.finite(sigma) | sigma <= 0) ||
      any(!is.finite(nu) | nu <= 0 | nu >= 1)) {
    return(NULL)
  }

  out <- numeric(n)
  nonnegative <- q >= 0
  poisson_limit <- nonnegative & sigma < 1e-4
  out[poisson_limit] <- stats::ppois(q[poisson_limit], mu[poisson_limit])
  regular <- which(nonnegative & !poisson_limit)
  if (!length(regular)) return(out)

  terms <- q[regular] + 1L
  chunk <- floor((cumsum(terms) - 1) / max_terms)
  for (rows in split(regular, chunk)) {
    counts <- q[rows] + 1L
    row_id <- rep.int(seq_along(rows), counts)
    k <- sequence(counts) - 1L
    lambda <- mu[rows] * nu[rows]
    nb_mean <- mu[rows] * (1 - nu[rows])
    nb_size <- 1 / sigma[rows]
    probability <- stats::dpois(k, lambda[row_id]) * stats::pnbinom(
      q[rows][row_id] - k,
      size = nb_size[row_id],
      mu = nb_mean[row_id]
    )
    out[rows] <- as.numeric(rowsum(probability, row_id, reorder = FALSE))
  }
  pmin(pmax(out, 0), 1)
}

.fast_delaporte_log_survival <- function(q, mu, sigma, nu, max_terms = 5e5) {
  n <- max(length(q), length(mu), length(sigma), length(nu))
  q <- rep(q, length.out = n)
  mu <- rep(mu, length.out = n)
  sigma <- rep(sigma, length.out = n)
  nu <- rep(nu, length.out = n)
  if (any(!is.finite(q)) || any(q < 0) || any(q != floor(q)) ||
      any(!is.finite(mu) | mu <= 0) ||
      any(!is.finite(sigma) | sigma <= 0) ||
      any(!is.finite(nu) | nu <= 0 | nu >= 1)) {
    return(NULL)
  }

  out <- rep(NA_real_, n)
  poisson_limit <- sigma < 1e-4
  out[poisson_limit] <- stats::ppois(
    q[poisson_limit], mu[poisson_limit], lower.tail = FALSE, log.p = TRUE
  )
  regular <- which(!poisson_limit)
  if (!length(regular)) return(out)

  terms <- q[regular] + 1L
  chunk <- floor((cumsum(terms) - 1) / max_terms)
  for (rows in split(regular, chunk)) {
    counts <- q[rows] + 1L
    row_id <- rep.int(seq_along(rows), counts)
    k <- sequence(counts) - 1L
    lambda <- mu[rows] * nu[rows]
    nb_mean <- mu[rows] * (1 - nu[rows])
    nb_size <- 1 / sigma[rows]
    log_terms <- stats::dpois(k, lambda[row_id], log = TRUE) +
      stats::pnbinom(
        q[rows][row_id] - k,
        size = nb_size[row_id],
        mu = nb_mean[row_id],
        lower.tail = FALSE,
        log.p = TRUE
      )
    maximum <- as.numeric(tapply(log_terms, row_id, max))
    log_sum <- rep(-Inf, length(rows))
    finite_maximum <- is.finite(maximum)
    if (any(finite_maximum)) {
      scaled <- exp(log_terms - maximum[row_id])
      scaled[!is.finite(scaled)] <- 0
      log_sum[finite_maximum] <- maximum[finite_maximum] + log(
        as.numeric(rowsum(scaled, row_id, reorder = FALSE))[finite_maximum]
      )
    }
    poisson_tail <- stats::ppois(
      q[rows], lambda, lower.tail = FALSE, log.p = TRUE
    )
    combined_maximum <- pmax(log_sum, poisson_tail)
    both_zero <- is.infinite(combined_maximum) & combined_maximum < 0
    combined <- rep(-Inf, length(rows))
    combined[!both_zero] <- combined_maximum[!both_zero] + log(
      exp(log_sum[!both_zero] - combined_maximum[!both_zero]) +
        exp(poisson_tail[!both_zero] - combined_maximum[!both_zero])
    )
    out[rows] <- combined
  }
  out
}

.fast_nbi_probability <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  n <- max(length(q), length(mu), length(sigma))
  q <- rep(q, length.out = n)
  mu <- rep(mu, length.out = n)
  sigma <- rep(sigma, length.out = n)
  if (any(!is.finite(mu) | mu <= 0) || any(!is.finite(sigma) | sigma <= 0)) {
    return(NULL)
  }
  out <- numeric(n)
  poisson_limit <- sigma <= 1e-4
  out[poisson_limit] <- stats::ppois(
    q[poisson_limit], mu[poisson_limit], lower.tail = lower.tail, log.p = log.p
  )
  regular <- !poisson_limit
  out[regular] <- stats::pnbinom(
    q[regular], size = 1 / sigma[regular], mu = mu[regular],
    lower.tail = lower.tail, log.p = log.p
  )
  out[q < 0] <- 0
  out
}

.fast_nbi_density <- function(x, mu, sigma, log = FALSE) {
  n <- max(length(x), length(mu), length(sigma))
  x <- rep(x, length.out = n)
  mu <- rep(mu, length.out = n)
  sigma <- rep(sigma, length.out = n)
  if (any(!is.finite(mu) | mu <= 0) || any(!is.finite(sigma) | sigma <= 0)) {
    return(NULL)
  }
  out <- numeric(n)
  poisson_limit <- sigma <= 1e-4
  out[poisson_limit] <- stats::dpois(
    x[poisson_limit], mu[poisson_limit], log = log
  )
  regular <- !poisson_limit
  out[regular] <- stats::dnbinom(
    x[regular], size = 1 / sigma[regular], mu = mu[regular], log = log
  )
  out[x < 0] <- 0
  out
}

.fast_bnb_probability <- function(
    q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  n <- max(length(q), length(mu), length(sigma), length(nu))
  q <- rep(q, length.out = n)
  mu <- rep(mu, length.out = n)
  sigma <- rep(sigma, length.out = n)
  nu <- rep(nu, length.out = n)
  if (any(!is.finite(q) | q != floor(q)) ||
      any(!is.finite(mu) | mu <= 0) ||
      any(!is.finite(sigma) | sigma <= 0) ||
      any(!is.finite(nu) | nu <= 0)) {
    return(NULL)
  }

  q_integer <- as.integer(q)
  probability <- numeric(n)
  active <- which(q_integer >= 0L)
  if (length(active)) {
    shape_m <- 1 / sigma + 1
    shape_n <- mu * nu / sigma
    shape_k <- 1 / nu
    mass <- numeric(n)
    mass[active] <- exp(
      lbeta(shape_n[active], shape_m[active] + shape_k[active]) -
        lbeta(shape_n[active], shape_m[active])
    )
    probability[active] <- mass[active]
    max_q <- max(q_integer[active])
    if (max_q >= 1L) {
      for (x in 0:(max_q - 1L)) {
        next_active <- active[q_integer[active] >= x + 1L]
        if (!length(next_active)) break
        mass[next_active] <- mass[next_active] *
          ((x + shape_n[next_active]) /
             (x + shape_n[next_active] + shape_m[next_active] + shape_k[next_active])) *
          ((x + shape_k[next_active]) / (x + 1))
        probability[next_active] <- probability[next_active] + mass[next_active]
      }
    }
  }
  probability <- pmin(pmax(probability, 0), 1)
  result <- if (isTRUE(lower.tail)) probability else 1 - probability
  if (isTRUE(log.p)) result <- log(result)
  result[q_integer < 0L] <- 0
  result
}

.fast_bnb_density <- function(x, mu, sigma, nu, log = FALSE) {
  n <- max(length(x), length(mu), length(sigma), length(nu))
  x <- rep(x, length.out = n)
  mu <- rep(mu, length.out = n)
  sigma <- rep(sigma, length.out = n)
  nu <- rep(nu, length.out = n)
  if (any(!is.finite(mu) | mu <= 0) ||
      any(!is.finite(sigma) | sigma <= 0) ||
      any(!is.finite(nu) | nu <= 0)) {
    return(NULL)
  }
  shape_m <- 1 / sigma + 1
  shape_n <- mu * nu / sigma
  shape_k <- 1 / nu
  log_mass <- lbeta(x + shape_n, shape_m + shape_k) -
    lbeta(shape_n, shape_m) - lgamma(x + 1) - lgamma(shape_k) +
    lgamma(x + shape_k)
  result <- if (isTRUE(log)) log_mass else exp(log_mass)
  result[x < 0] <- 0
  result
}

.fast_zibnb_probability <- function(
    q, mu, sigma, nu, tau, lower.tail = TRUE, log.p = FALSE) {
  n <- max(length(q), length(mu), length(sigma), length(nu), length(tau))
  tau <- rep(tau, length.out = n)
  if (any(!is.finite(tau) | tau <= 0 | tau >= 1)) return(NULL)
  base <- .fast_bnb_probability(
    q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE
  )
  if (is.null(base)) return(NULL)
  probability <- tau + (1 - tau) * base
  result <- if (isTRUE(lower.tail)) probability else 1 - probability
  if (isTRUE(log.p)) result <- log(result)
  result[rep(q, length.out = n) < 0] <- 0
  result
}

.fast_zibnb_density <- function(x, mu, sigma, nu, tau, log = FALSE) {
  n <- max(length(x), length(mu), length(sigma), length(nu), length(tau))
  x <- rep(x, length.out = n)
  tau <- rep(tau, length.out = n)
  if (any(!is.finite(tau) | tau <= 0 | tau >= 1)) return(NULL)
  base_log_mass <- .fast_bnb_density(x, mu, sigma, nu, log = TRUE)
  if (is.null(base_log_mass)) return(NULL)
  log_mass <- log1p(-tau) + base_log_mass
  at_zero <- x == 0
  log_mass[at_zero] <- log(
    tau[at_zero] + (1 - tau[at_zero]) * exp(base_log_mass[at_zero])
  )
  result <- if (isTRUE(log)) log_mass else exp(log_mass)
  result[x < 0] <- 0
  result
}

.fast_zero_modified_probability <- function(
    q, mu, zero_probability, base = c("poisson", "nbi"), sigma = NULL,
    hurdle = FALSE, lower.tail = TRUE, log.p = FALSE) {
  base <- match.arg(base)
  n <- max(length(q), length(mu), length(zero_probability), length(sigma))
  q <- rep(q, length.out = n)
  mu <- rep(mu, length.out = n)
  zero_probability <- rep(zero_probability, length.out = n)
  if (any(!is.finite(mu) | mu <= 0) ||
      any(!is.finite(zero_probability) |
          zero_probability <= 0 | zero_probability >= 1)) {
    return(NULL)
  }
  if (identical(base, "nbi")) {
    sigma <- rep(sigma, length.out = n)
    if (any(!is.finite(sigma) | sigma <= 0)) return(NULL)
    base_log_survival <- .fast_nbi_probability(
      q, mu, sigma, lower.tail = FALSE, log.p = TRUE
    )
  } else {
    base_log_survival <- stats::ppois(
      q, mu, lower.tail = FALSE, log.p = TRUE
    )
  }
  if (is.null(base_log_survival)) return(NULL)

  log_survival <- log1p(-zero_probability) + base_log_survival
  if (isTRUE(hurdle)) {
    base_log_positive <- if (identical(base, "poisson")) {
      stats::ppois(0, mu, lower.tail = FALSE, log.p = TRUE)
    } else {
      .fast_nbi_probability(
        rep(0, n), mu, sigma, lower.tail = FALSE, log.p = TRUE
      )
    }
    log_survival <- log_survival - base_log_positive
  }
  out <- if (isTRUE(lower.tail)) {
    lower <- -expm1(log_survival)
    if (isTRUE(log.p)) log(lower) else lower
  } else if (isTRUE(log.p)) {
    log_survival
  } else {
    exp(log_survival)
  }
  out[q < 0] <- 0
  out
}

.fast_zero_modified_density <- function(
    x, mu, zero_probability, base = c("poisson", "nbi"), sigma = NULL,
    hurdle = FALSE, log = FALSE) {
  base <- match.arg(base)
  n <- max(length(x), length(mu), length(zero_probability), length(sigma))
  x <- rep(x, length.out = n)
  mu <- rep(mu, length.out = n)
  zero_probability <- rep(zero_probability, length.out = n)
  if (any(!is.finite(mu) | mu <= 0) ||
      any(!is.finite(zero_probability) |
          zero_probability <= 0 | zero_probability >= 1)) {
    return(NULL)
  }
  if (identical(base, "nbi")) {
    sigma <- rep(sigma, length.out = n)
    base_log_mass <- .fast_nbi_density(x, mu, sigma, log = TRUE)
    base_log_zero <- .fast_nbi_density(rep(0, n), mu, sigma, log = TRUE)
  } else {
    base_log_mass <- stats::dpois(x, mu, log = TRUE)
    base_log_zero <- stats::dpois(0, mu, log = TRUE)
  }
  if (is.null(base_log_mass) || is.null(base_log_zero)) return(NULL)

  log_mass <- log1p(-zero_probability) + base_log_mass
  if (isTRUE(hurdle)) {
    log_mass <- log_mass - log1p(-exp(base_log_zero))
    log_mass[x == 0] <- log(zero_probability[x == 0])
  } else {
    at_zero <- x == 0
    log_mass[at_zero] <- log(
      zero_probability[at_zero] +
        (1 - zero_probability[at_zero]) * exp(base_log_zero[at_zero])
    )
  }
  out <- if (isTRUE(log)) log_mass else exp(log_mass)
  out[x < 0] <- 0
  out
}

.call_fast_count_family <- function(prefix, family_name, args) {
  if (!"mu" %in% names(args)) return(NULL)

  x <- args$x %||% args$q %||% args$y
  mu <- args$mu
  if (is.null(x) || is.null(mu)) {
    return(NULL)
  }

  if (identical(family_name, "PO")) {
    if (prefix %in% c("d", "logd")) {
      return(stats::dpois(x, lambda = mu, log = identical(prefix, "logd")))
    }
    if (identical(prefix, "p")) {
      return(stats::ppois(
        x,
        lambda = mu,
        lower.tail = args$lower.tail %||% TRUE,
        log.p = args$log.p %||% FALSE
      ))
    }
  }

  if (identical(family_name, "NBI") && "sigma" %in% names(args)) {
    if (prefix %in% c("d", "logd")) {
      return(.fast_nbi_density(
        x, mu, args$sigma, log = identical(prefix, "logd")
      ))
    }
    if (identical(prefix, "p")) {
      return(.fast_nbi_probability(
        x, mu, args$sigma,
        lower.tail = args$lower.tail %||% TRUE,
        log.p = args$log.p %||% FALSE
      ))
    }
  }

  if (identical(family_name, "BNB") &&
      all(c("sigma", "nu") %in% names(args))) {
    if (prefix %in% c("d", "logd")) {
      return(.fast_bnb_density(
        x, mu, args$sigma, args$nu, log = identical(prefix, "logd")
      ))
    }
    if (identical(prefix, "p")) {
      return(.fast_bnb_probability(
        x, mu, args$sigma, args$nu,
        lower.tail = args$lower.tail %||% TRUE,
        log.p = args$log.p %||% FALSE
      ))
    }
  }

  if (identical(family_name, "ZIBNB") &&
      isTRUE(getOption("gamlss.longitudinal.fast_zibnb_family", TRUE)) &&
      all(c("sigma", "nu", "tau") %in% names(args))) {
    if (prefix %in% c("d", "logd")) {
      return(.fast_zibnb_density(
        x, mu, args$sigma, args$nu, args$tau,
        log = identical(prefix, "logd")
      ))
    }
    if (identical(prefix, "p")) {
      return(.fast_zibnb_probability(
        x, mu, args$sigma, args$nu, args$tau,
        lower.tail = args$lower.tail %||% TRUE,
        log.p = args$log.p %||% FALSE
      ))
    }
  }

  if (family_name %in% c("ZIP", "ZAP") && "sigma" %in% names(args)) {
    if (prefix %in% c("d", "logd")) {
      return(.fast_zero_modified_density(
        x, mu, args$sigma,
        base = "poisson",
        hurdle = identical(family_name, "ZAP"),
        log = identical(prefix, "logd")
      ))
    }
    if (identical(prefix, "p")) {
      return(.fast_zero_modified_probability(
        x, mu, args$sigma,
        base = "poisson",
        hurdle = identical(family_name, "ZAP"),
        lower.tail = args$lower.tail %||% TRUE,
        log.p = args$log.p %||% FALSE
      ))
    }
  }

  if (identical(family_name, "ZINBI") &&
      all(c("sigma", "nu") %in% names(args))) {
    if (prefix %in% c("d", "logd")) {
      return(.fast_zero_modified_density(
        x, mu, args$nu,
        base = "nbi", sigma = args$sigma,
        log = identical(prefix, "logd")
      ))
    }
    if (identical(prefix, "p")) {
      return(.fast_zero_modified_probability(
        x, mu, args$nu,
        base = "nbi", sigma = args$sigma,
        lower.tail = args$lower.tail %||% TRUE,
        log.p = args$log.p %||% FALSE
      ))
    }
  }

  if (identical(family_name, "DEL") && identical(prefix, "p") &&
      isTRUE(args$lower.tail %||% TRUE) && !isTRUE(args$log.p %||% FALSE) &&
      all(c("sigma", "nu") %in% names(args))) {
    probability <- .fast_delaporte_cdf(x, mu, args$sigma, args$nu)
    if (is.null(probability)) return(NULL)
    return(probability)
  }
  if (identical(family_name, "DEL") && identical(prefix, "p") &&
      !isTRUE(args$lower.tail %||% TRUE) && isTRUE(args$log.p %||% FALSE) &&
      all(c("sigma", "nu") %in% names(args))) {
    return(.fast_delaporte_log_survival(x, mu, args$sigma, args$nu))
  }
  NULL
}
