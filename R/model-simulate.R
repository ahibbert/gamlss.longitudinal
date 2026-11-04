#' Simulate responses from a longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param nsim Number of simulated response columns.
#' @param seed Optional random seed.
#' @param newdata Optional new data. With the default `condition_on = NULL`,
#'   simulation is unconditional: subject identifiers define independent
#'   trajectory groups, and the fitted marginal and copula formulas are
#'   evaluated on the supplied rows. Responses may be missing.
#' @param condition_on Conditioning history for forward simulation. The default,
#'   `NULL`, retains unconditional simulation. Use `"fitted"` to condition each
#'   subject's first immediately adjacent future visit on the most recent usable
#'   response in the fitted model, or supply a data frame containing custom
#'   response histories. Conditional simulation requires `newdata`, which must
#'   contain future rows only. A missing scheduled visit starts a new independent
#'   segment; subsequent contiguous future visits remain dependent recursively.
#' @param ... Additional arguments reserved for future methods.
#'
#' @details Conditional simulation uses the package's global scheduled-visit
#'   ordering rather than numeric time distance. It conditions only on the most
#'   recent usable response under the fitted first-order adjacent-pair model; it
#'   is not a latent-state filter over an incomplete history. Continuous history
#'   responses are transformed by their fitted marginal CDF. For a discrete
#'   response, a latent uniform is randomized over its fitted CDF interval and
#'   retained while future values are generated. Marginal predictions remain the
#'   responsibility of [predict.gamlss.longitudinal()]; analytic conditional
#'   prediction is not performed here.
#'
#' @return A data frame with one column per simulation. Conditional simulations
#'   contain future rows only and have a `conditioning` attribute recording the
#'   history source and each conditioned or unconditional segment start.
#' @importFrom stats simulate
#' @export
simulate.gamlss.longitudinal <- function(
    object,
    nsim = 1,
    seed = NULL,
    newdata = NULL,
    condition_on = NULL,
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

  conditional <- !is.null(condition_on)
  if (conditional && is.null(newdata)) {
    stop("'newdata' is required when 'condition_on' is used.", call. = FALSE)
  }

  if (!conditional && is.null(newdata)) {
    diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)
    sim_mat <- .gl_simulate_copula_matrix(object, diag_data, nsim)
  } else if (!conditional) {
    sim_data <- .gl_simulation_newdata(object, newdata)
    diag_data <- sim_data$diag_data
    sim_mat <- .gl_simulate_copula_matrix(
      object,
      diag_data,
      nsim,
      fit_data = sim_data$fit_data,
      time_levels = sim_data$time_levels
    )
  } else {
    setup <- .gl_conditional_simulation_setup(
      object = object,
      newdata = newdata,
      condition_on = condition_on
    )
    simulated <- .gl_simulate_conditional_copula_matrix(
      object = object,
      setup = setup,
      nsim = nsim
    )
    sim_mat <- simulated$values
  }

  if (!conditional && is.null(newdata) && nrow(sim_mat) != length(object$response)) {
    full_sim_mat <- matrix(NA_real_, nrow = length(object$response), ncol = nsim)
    full_sim_mat[diag_data$keep_index, ] <- sim_mat
    sim_mat <- full_sim_mat
  }

  out <- as.data.frame(sim_mat, stringsAsFactors = FALSE)
  names(out) <- paste0("sim_", seq_len(nsim))

  if (conditional) {
    attr(out, "conditioning") <- simulated$metadata
    .gl_warn_conditional_simulation_gaps(simulated$metadata)
  }

  out
}
