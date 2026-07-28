#' Simulate responses from a longitudinal GAMLSS-copula model

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param nsim Number of simulated response columns.

#' @param seed Optional random seed.

#' @param newdata Optional new data. When supplied, simulation is unconditional:

#'   subject identifiers define independent trajectory groups, and the fitted

#'   marginal and copula formulas are evaluated on the supplied rows. Responses

#'   may be missing.

#' @param ... Additional arguments reserved for future methods.

#'

#' @return A data frame with one column per simulation.

#' @importFrom stats simulate

#' @export

simulate.gamlss.longitudinal <- function(
    object,
    nsim = 1,
    seed = NULL,
    newdata = NULL,
    ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  nsim <- as.integer(nsim)

  if (!is.finite(nsim) || nsim < 1L) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
  }

  dots <- list(...)

  if ("type" %in% names(dots)) {
    type <- dots$type

    if (!identical(type, "copula")) {
      stop("'type' is no longer supported; simulate() always uses the fitted copula model.", call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }

    on.exit(
      {
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      },
      add = TRUE
    )

    set.seed(seed)
  }

  if (is.null(newdata)) {
    diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)

    sim_mat <- .gl_simulate_copula_matrix(object, diag_data, nsim)
  } else {
    sim_data <- .gl_simulation_newdata(object, newdata)

    diag_data <- sim_data$diag_data

    sim_mat <- .gl_simulate_copula_matrix(

      object,
      diag_data,
      nsim,
      fit_data = sim_data$fit_data,
      time_levels = sim_data$time_levels
    )
  }

  if (is.null(newdata) && nrow(sim_mat) != length(object$response)) {
    full_sim_mat <- matrix(NA_real_, nrow = length(object$response), ncol = nsim)

    full_sim_mat[diag_data$keep_index, ] <- sim_mat

    sim_mat <- full_sim_mat
  }

  out <- as.data.frame(sim_mat, stringsAsFactors = FALSE)

  names(out) <- paste0("sim_", seq_len(nsim))

  out
}
