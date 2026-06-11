#' Rescale a numeric vector to the unit interval
#'
#' @param x Numeric vector.
#' @return Numeric vector scaled to `[0, 1]`. Constant finite inputs return 0.
#' @export
sim_rescale01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng))) {
    return(rep(NA_real_, length(x)))
  }
  width <- diff(rng)
  if (width <= 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / width
}

#' Simple deterministic smooth shapes for simulation truth
#'
#' These helpers are intended for readable simulation specifications. They work
#' best with inputs already scaled to `[0, 1]`, for example by
#' [sim_rescale01()].
#'
#' @param x Numeric input, usually scaled to `[0, 1]`.
#' @param slope Linear slope.
#' @param center Logical; if `TRUE`, center `x` around 0.5 for the linear shape.
#' @param amplitude Effect amplitude.
#' @param period Period of the sinusoid in units of `x`.
#' @param phase Phase shift for the sinusoid.
#' @param location Center/location of the bump, sigmoid, or U-shape.
#' @param width Width of the Gaussian bump.
#' @param steepness Steepness of the sigmoid transition.
#'
#' @return Numeric vector of the same length as `x`.
#' @name sim_smooth_shapes
NULL

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_linear <- function(x, slope = 1, center = TRUE) {
  x <- as.numeric(x)
  if (isTRUE(center)) x <- x - 0.5
  slope * x
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_sin <- function(x, amplitude = 1, period = 1, phase = 0) {
  x <- as.numeric(x)
  amplitude * sin(2 * pi * (x / period + phase))
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_bump <- function(x, amplitude = 1, location = 0.5, width = 0.15) {
  x <- as.numeric(x)
  amplitude * exp(-0.5 * ((x - location) / width)^2)
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_sigmoid <- function(x, amplitude = 1, location = 0.5, steepness = 10) {
  x <- as.numeric(x)
  amplitude * (stats::plogis(steepness * (x - location)) - 0.5)
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_u <- function(x, amplitude = 1, location = 0.5) {
  x <- as.numeric(x)
  amplitude * ((x - location)^2 - mean((x - location)^2, na.rm = TRUE))
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_wiggle <- function(x, amplitude = 1) {
  x <- as.numeric(x)
  amplitude * (0.65 * sin(2 * pi * x) + 0.35 * sin(6 * pi * x + 0.4))
}

#' Look up named effects for a factor-like variable
#'
#' @param x Factor, character, or numeric vector.
#' @param effects Named numeric vector of effects.
#' @param reference Optional reference level to prepend with effect 0 when it is
#'   not already present in `effects`.
#'
#' @return Numeric vector of effects aligned with `x`.
#' @export
sim_factor_effect <- function(x, effects, reference = NULL) {
  if (is.null(names(effects)) || any(!nzchar(names(effects)))) {
    stop("effects must be a named numeric vector.", call. = FALSE)
  }
  effects <- as.numeric(effects) |>
    stats::setNames(names(effects))
  if (!is.null(reference) && !reference %in% names(effects)) {
    effects <- c(stats::setNames(0, reference), effects)
  }
  x_chr <- as.character(x)
  out <- effects[x_chr]
  if (anyNA(out)) {
    missing_levels <- unique(x_chr[is.na(out)])
    stop(
      "Missing factor effect for level(s): ",
      paste(missing_levels, collapse = ", "),
      ". Available effect names: ",
      paste(names(effects), collapse = ", "),
      ". Supply one named effect per observed level, or use reference to set ",
      "a zero effect for a missing reference level.",
      call. = FALSE
    )
  }
  as.numeric(out)
}

#' Generate covariates for longitudinal simulations
#'
#' This helper is designed to be used inside the `covariates` argument of
#' [simulate_longitudinal_dataset()]. Subject-level variables are evaluated once
#' per subject and repeated across time; observation-level variables are
#' evaluated once per row of the long-format design.
#'
#' @param data Base long-format design data supplied by
#'   [simulate_longitudinal_dataset()] to a covariate callback.
#' @param subject Named list of subject-level covariate specifications.
#' @param observation Named list of observation-level covariate specifications.
#'
#' @return A data frame of covariates aligned row-for-row with `data`.
#' @export
simulate_longitudinal_covariates <- function(data, subject = list(), observation = list()) {
  if (!is.data.frame(data) || !".sim_subject_index" %in% names(data)) {
    stop(
      "data must be the base design supplied by simulate_longitudinal_dataset().",
      call. = FALSE
    )
  }
  out <- data.frame(row.names = seq_len(nrow(data)))
  subject_rows <- !duplicated(data$.sim_subject_index)
  subject_data <- data[subject_rows, , drop = FALSE]

  for (nm in names(subject)) {
    values <- .sim_eval_covariate_spec(subject[[nm]], subject_data, nm)
    if (length(values) == 1L) {
      values <- rep(values, nrow(subject_data))
    }
    if (length(values) != nrow(subject_data)) {
      stop("Subject covariate '", nm, "' must return one value per subject.", call. = FALSE)
    }
    out[[nm]] <- values[match(data$.sim_subject_index, subject_data$.sim_subject_index)]
  }

  for (nm in names(observation)) {
    values <- .sim_eval_covariate_spec(observation[[nm]], data, nm)
    if (length(values) == 1L) {
      values <- rep(values, nrow(data))
    }
    if (length(values) != nrow(data)) {
      stop("Observation covariate '", nm, "' must return one value per row.", call. = FALSE)
    }
    out[[nm]] <- values
  }

  out
}

.sim_eval_covariate_spec <- function(spec, data, label) {
  if (is.function(spec)) {
    return(spec(data))
  }
  spec
}

#' Simulate longitudinal GAMLSS-copula data
#'
#' Simulates a balanced longitudinal dataset from a GAMLSS marginal distribution
#' and a first-order copula dependence model between consecutive time points.
#' Parameter specifications may be constants, time-indexed vectors,
#' subject-by-time matrices, or functions of the generated design data.
#'
#' @param n Number of subjects.
#' @param times Vector of observed time values.
#' @param margin_dist A `gamlss.dist` family object, for example `NO()` or
#'   `GA(mu.link = "log", sigma.link = "log")`.
#' @param copula_dist Copula family code. Supported codes are `"N"`, `"C"`,
#'   `"F"`, `"G"`, `"J"`, and `"t"`.
#' @param margin_params Named list of marginal distribution parameter
#'   specifications. Names should match the quantile function arguments, such
#'   as `mu`, `sigma`, `nu`, and `tau`.
#' @param copula_params Named list of copula parameter specifications. Use
#'   `theta` or `par` for the primary copula parameter, `tau` to specify
#'   Kendall's tau instead, and `zeta` or `par2` for the t-copula degrees of
#'   freedom.
#' @param covariates Optional data frame or function returning covariates. A
#'   data frame may have either `n` rows for subject-level covariates or
#'   `n * length(times)` rows for long-format covariates. A function is called
#'   with the base long-format data. It may return only new covariate columns
#'   or the input data with new columns added; columns already present in the
#'   base design are ignored.
#' @param seed Optional random seed.
#' @param subject_var,time_var,response_var Column names for the subject, time,
#'   and response variables.
#' @param include_truth If `TRUE`, include simulated uniforms and true
#'   parameter columns.
#' @param u_bounds Optional length-two numeric vector giving lower and upper
#'   bounds used to clamp simulated uniforms before applying the marginal
#'   quantile function. The default `NULL` leaves uniforms unchanged.
#'
#' @return A long-format data frame with one row per subject-time observation.
#' @export
simulate_longitudinal_dataset <- function(
  n = 100,
  times = seq_len(3),
  margin_dist = gamlss.dist::NO(),
  copula_dist = "N",
  margin_params = list(mu = 0, sigma = 1),
  copula_params = list(theta = 0),
  covariates = NULL,
  seed = NULL,
  subject_var = "subject",
  time_var = "time",
  response_var = "response",
  include_truth = TRUE,
  u_bounds = NULL
) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  n <- .sim_check_count(n, "n")
  if (length(times) < 2L) {
    stop("times must contain at least two values.", call. = FALSE)
  }
  if (anyDuplicated(times)) {
    stop("times must not contain duplicate values.", call. = FALSE)
  }
  copula_dist <- .copula_family_code(copula_dist)

  long <- .sim_base_long_data(n, times, subject_var, time_var)
  long <- .sim_add_covariates(long, covariates, n, length(times), subject_var)

  u_mat <- .sim_copula_uniform_matrix(
    n = n,
    n_time = length(times),
    family = copula_dist,
    copula_params = copula_params,
    long_data = long,
    subject_var = subject_var,
    time_var = time_var
  )

  long$.sim_u <- as.vector(t(u_mat))
  if (!is.null(u_bounds)) {
    if (!is.numeric(u_bounds) || length(u_bounds) != 2L ||
        any(!is.finite(u_bounds)) || u_bounds[1L] < 0 ||
        u_bounds[2L] > 1 || u_bounds[1L] >= u_bounds[2L]) {
      stop("u_bounds must be NULL or a finite increasing length-two vector inside [0, 1].", call. = FALSE)
    }
    long$.sim_u <- pmin(pmax(long$.sim_u, u_bounds[1L]), u_bounds[2L])
  }
  qfun <- .sim_margin_quantile_function(margin_dist)
  qargs <- list(p = long$.sim_u)
  for (param_name in names(margin_params)) {
    qargs[[param_name]] <- .sim_eval_long_param(
      margin_params[[param_name]],
      long,
      n = n,
      n_time = length(times),
      label = paste0("margin_params$", param_name)
    )
  }
  qargs <- qargs[names(qargs) %in% formalArgs(qfun)]
  long[[response_var]] <- do.call(qfun, qargs)

  if (include_truth) {
    long$u <- long$.sim_u
    for (param_name in names(margin_params)) {
      long[[paste0("true_", param_name)]] <- .sim_eval_long_param(
        margin_params[[param_name]],
        long,
        n = n,
        n_time = length(times),
        label = paste0("margin_params$", param_name)
      )
    }
    edge_truth <- attr(u_mat, "edge_truth")
    if (!is.null(edge_truth)) {
      long$true_theta <- NA_real_
      long$true_zeta <- NA_real_
      edge_rows <- long$.sim_time_index > 1L
      long$true_theta[edge_rows] <- as.vector(t(edge_truth$theta))
      long$true_zeta[edge_rows] <- as.vector(t(edge_truth$zeta))
    }
  }

  long$.sim_u <- NULL
  long$.sim_time_index <- NULL
  long$.sim_subject_index <- NULL
  front_cols <- c(subject_var, time_var, response_var)
  long <- long[c(front_cols, setdiff(names(long), front_cols))]
  rownames(long) <- NULL
  long
}

.sim_check_count <- function(x, label) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(label, " must be a single positive integer.", call. = FALSE)
  }
  x <- as.integer(x)
  if (!is.finite(x) || x < 1L) {
    stop(label, " must be a single positive integer.", call. = FALSE)
  }
  x
}

.sim_base_long_data <- function(n, times, subject_var, time_var) {
  data.frame(
    .sim_subject_index = rep(seq_len(n), each = length(times)),
    .sim_time_index = rep(seq_along(times), times = n),
    stringsAsFactors = FALSE
  ) |>
    within({
      subject <- factor(.sim_subject_index)
      time <- times[.sim_time_index]
    }) |>
    .sim_rename_base_columns(subject_var, time_var)
}

.sim_rename_base_columns <- function(x, subject_var, time_var) {
  names(x)[names(x) == "subject"] <- subject_var
  names(x)[names(x) == "time"] <- time_var
  x
}

.sim_add_covariates <- function(long, covariates, n, n_time, subject_var) {
  if (is.null(covariates)) {
    return(long)
  }
  if (is.function(covariates)) {
    covariates <- covariates(long)
  }
  covariates <- as.data.frame(covariates, stringsAsFactors = FALSE)
  if (nrow(covariates) == n) {
    covariates <- covariates[long$.sim_subject_index, , drop = FALSE]
  } else if (nrow(covariates) != n * n_time) {
    stop(
      "covariates must have either n rows or n * length(times) rows.",
      call. = FALSE
    )
  }
  overlap <- intersect(names(long), names(covariates))
  if (length(overlap) > 0L) {
    covariates <- covariates[, setdiff(names(covariates), overlap), drop = FALSE]
  }
  if (ncol(covariates) == 0L) return(long)
  cbind(long, covariates)
}

.sim_margin_quantile_function <- function(margin_dist) {
  family_name <- margin_dist$family[1]
  if (!is.character(family_name) || length(family_name) != 1L || is.na(family_name)) {
    stop("margin_dist must be a gamlss.dist family object.", call. = FALSE)
  }
  qfun_name <- paste0("q", family_name)
  if (!exists(qfun_name, envir = asNamespace("gamlss.dist"), inherits = FALSE)) {
    stop("Could not find gamlss.dist quantile function ", qfun_name, "().", call. = FALSE)
  }
  get(qfun_name, envir = asNamespace("gamlss.dist"), inherits = FALSE)
}

.sim_eval_long_param <- function(spec, data, n, n_time, label) {
  if (is.function(spec)) {
    spec <- spec(data)
  }
  .sim_expand_param(spec, n = n, n_time = n_time, index = data$.sim_time_index, label = label)
}

.sim_eval_edge_param <- function(spec, edge_data, n, n_edge, label) {
  if (is.function(spec)) {
    spec <- spec(edge_data)
  }
  .sim_expand_param(spec, n = n, n_time = n_edge, index = edge_data$.sim_edge_index, label = label)
}

.sim_expand_param <- function(spec, n, n_time, index, label) {
  if (is.matrix(spec)) {
    if (!identical(dim(spec), c(n, n_time))) {
      stop(label, " matrix must have dimensions n by number of time positions.", call. = FALSE)
    }
    return(as.vector(t(spec)))
  }
  spec <- as.numeric(spec)
  if (length(spec) == 1L) {
    return(rep(spec, length(index)))
  }
  if (length(spec) == n_time) {
    return(spec[index])
  }
  if (length(spec) == length(index)) {
    return(spec)
  }
  stop(
    label,
    " must be a scalar, a time-position vector, a full-length vector, ",
    "a subject-by-time matrix, or a function returning one of those.",
    call. = FALSE
  )
}

.sim_copula_uniform_matrix <- function(n, n_time, family, copula_params, long_data, subject_var, time_var) {
  u <- matrix(stats::runif(n * n_time), nrow = n, ncol = n_time)
  if (n_time == 1L) {
    return(u)
  }

  edge_data <- .sim_edge_data(long_data, n_time, subject_var, time_var)
  theta <- .sim_resolve_copula_theta(copula_params, edge_data, n, n_time - 1L, family)
  zeta <- .sim_resolve_copula_zeta(copula_params, edge_data, n, n_time - 1L, family)

  theta_mat <- matrix(theta, nrow = n, ncol = n_time - 1L, byrow = TRUE)
  zeta_mat <- matrix(zeta, nrow = n, ncol = n_time - 1L, byrow = TRUE)

  for (time_index in 2:n_time) {
    target <- stats::runif(n)
    edge_index <- time_index - 1L
    for (subject_index in seq_len(n)) {
      u[subject_index, time_index] <- .sim_invert_hfunc1(
        u1 = u[subject_index, time_index - 1L],
        target = target[subject_index],
        family = family,
        par = theta_mat[subject_index, edge_index],
        par2 = zeta_mat[subject_index, edge_index]
      )
    }
  }

  attr(u, "edge_truth") <- list(theta = theta_mat, zeta = zeta_mat)
  u
}

.sim_edge_data <- function(long_data, n_time, subject_var, time_var) {
  edge_data <- long_data[long_data$.sim_time_index < n_time, , drop = FALSE]
  edge_data$.sim_edge_index <- edge_data$.sim_time_index
  edge_data$time_left <- edge_data[[time_var]]
  edge_data$time_right <- long_data[[time_var]][long_data$.sim_time_index > 1L]
  edge_data
}

.sim_resolve_copula_theta <- function(copula_params, edge_data, n, n_edge, family) {
  if (!is.null(copula_params$tau)) {
    tau <- .sim_eval_edge_param(copula_params$tau, edge_data, n, n_edge, "copula_params$tau")
    return(.copula_tau_to_par(family, tau))
  }
  theta_spec <- copula_params$theta
  if (is.null(theta_spec)) {
    theta_spec <- copula_params$par
  }
  if (is.null(theta_spec)) {
    stop("copula_params must include theta, par, or tau.", call. = FALSE)
  }
  .sim_eval_edge_param(theta_spec, edge_data, n, n_edge, "copula_params$theta")
}

.sim_resolve_copula_zeta <- function(copula_params, edge_data, n, n_edge, family) {
  zeta_spec <- copula_params$zeta
  if (is.null(zeta_spec)) {
    zeta_spec <- copula_params$par2
  }
  if (is.null(zeta_spec)) {
    zeta_spec <- if (identical(family, "t")) 4 else 0
  }
  .sim_eval_edge_param(zeta_spec, edge_data, n, n_edge, "copula_params$zeta")
}

.sim_invert_hfunc1 <- function(u1, target, family, par, par2) {
  eps <- sqrt(.Machine$double.eps)
  target <- min(max(target, eps), 1 - eps)
  objective <- function(u2) {
    .copula_hfunc1(u1, u2, family = family, par = par, par2 = par2) - target
  }
  lower_value <- objective(eps)
  upper_value <- objective(1 - eps)
  if (!is.finite(lower_value) || !is.finite(upper_value)) {
    stop("Non-finite copula h-function encountered during simulation.", call. = FALSE)
  }
  if (lower_value >= 0) {
    return(eps)
  }
  if (upper_value <= 0) {
    return(1 - eps)
  }
  stats::uniroot(objective, interval = c(eps, 1 - eps), tol = 1e-10)$root
}
