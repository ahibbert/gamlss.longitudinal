.gl_simulation_time_levels <- function(time, reference_time = NULL) {

  if (is.factor(time) || is.factor(reference_time)) {

    lev <- unique(c(levels(reference_time), levels(time)))

    return(lev[!is.na(lev)])

  }

  u <- unique(c(reference_time, time))

  if (is.numeric(u) || is.integer(u)) sort(u) else sort(as.character(u))

}


.gl_simulation_newdata <- function(object, newdata) {

  copula_spec <- get_copula_dist(object$copula_dist)

  nd <- .gl_prepare_newdata_internal(object, newdata, require_response = FALSE)

  nd$.gl_sim_row_id <- seq_len(nrow(nd))

  time_for_grid <- if ("time_covariate" %in% names(nd) && is.factor(nd$time_covariate)) {

    nd$time_covariate

  } else {

    nd$time

  }


  time_levels <- .gl_simulation_time_levels(time_for_grid, reference_time = object$response_margin)

  time_lookup <- stats::setNames(seq_along(time_levels), as.character(time_levels))

  nd$.gl_sim_time_idx <- unname(time_lookup[as.character(time_for_grid)])

  if (any(!is.finite(nd$.gl_sim_time_idx))) {

    stop("Could not map newdata time values to an ordered simulation grid.", call. = FALSE)

  }

  if (anyDuplicated(paste(nd$subject, as.character(time_for_grid), sep = "\r"))) {

    stop("newdata must contain at most one row per subject/time combination.", call. = FALSE)

  }


  nd_eval <- nd[order(nd$subject, nd$.gl_sim_time_idx, nd$.gl_sim_row_id), , drop = FALSE]

  mm_use <- do.call(

    create_model_matrices,

    list(

      mu.formula = object$formulas_int$mu,

      sigma.formula = object$formulas_int$sigma,

      nu.formula = object$formulas_int$nu,

      tau.formula = object$formulas_int$tau,

      theta.formula = object$formulas_int$theta,

      zeta.formula = object$formulas_int$zeta,

      margin.family = object$margin_dist,

      copula.family = object$copula_dist,

      copula.link = copula_spec$copula_link,

      dataset = nd_eval,

      quiet_gamlss2 = TRUE,

      preserve_factor_levels = TRUE

    )

  )

  mm_use <- .gl_align_model_matrix_columns(mm_use, object$model_matrix)

  eta_out <- calc_eta(

    par_cov = object$par,

    mm = mm_use,

    margin_dist = object$margin_dist,

    copula_link = copula_spec$copula_link,

    par_s = object$par_s

  )


  margin_params <- eta_out$eta_inv[names(object$margin_dist$parameters)]

  common_n <- min(

    nrow(nd_eval),

    if (length(margin_params) > 0) min(vapply(margin_params, length, integer(1))) else nrow(nd_eval)

  )

  if (!is.finite(common_n) || common_n < 1L) {

    stop("No newdata rows are available for simulation.", call. = FALSE)

  }


  nd_eval <- nd_eval[seq_len(common_n), , drop = FALSE]

  margin_params <- lapply(margin_params, function(x) x[seq_len(common_n)])

  keep <- rep(TRUE, common_n)

  for (par_name in names(margin_params)) {

    keep <- keep & is.finite(margin_params[[par_name]])

  }

  if (!all(keep)) {

    stop("newdata produced non-finite fitted marginal parameters.", call. = FALSE)

  }


  align_to_rows <- function(param_vec) {

    if (length(param_vec) == 0L) {

      return(rep(NA_real_, common_n))

    }

    if (length(param_vec) == common_n) {

      return(as.numeric(param_vec))

    }

    left_rows <- which(nd_eval$.gl_sim_time_idx < length(time_levels))

    if (length(param_vec) == length(left_rows)) {

      out <- rep(NA_real_, common_n)

      out[left_rows] <- as.numeric(param_vec)

      return(out)

    }

    rep(as.numeric(param_vec), length.out = common_n)

  }


  theta_fit <- align_to_rows(eta_out$eta_inv$theta %||% numeric(0))

  zeta_fit <- align_to_rows(eta_out$eta_inv$zeta %||% numeric(0))

  ord <- order(nd_eval$.gl_sim_row_id)

  margin_params <- lapply(margin_params, function(x) as.numeric(x[ord]))


  list(

    diag_data = list(

      response = nd_eval$response[ord],

      params = margin_params,

      family = object$margin_dist$family[1],

      subject = nd_eval$subject[ord],

      time = nd_eval$time[ord],

      keep_index = nd_eval$.gl_sim_row_id[ord]

    ),

    time_levels = time_levels,

    fit_data = data.frame(

      subject = nd_eval$subject[ord],

      time = time_for_grid[nd_eval$.gl_sim_row_id[ord]],

      theta_fit = theta_fit[ord],

      zeta_fit = zeta_fit[ord],

      stringsAsFactors = FALSE

    )

  )

}


.gl_simulate_copula_matrix <- function(object, diag_data, nsim, fit_data = NULL, time_levels = NULL) {

  if (is.null(fit_data)) {

    fit_data <- .copula_v2_fit_data(object)

  }

  n <- nrow(fit_data)

  qfun <- get(paste0("q", diag_data$family), envir = asNamespace("gamlss.dist"), inherits = FALSE)

  out <- matrix(NA_real_, nrow = n, ncol = nsim)


  if (is.null(time_levels)) {

    time_levels <- .gl_simulation_time_levels(fit_data$time)

  }

  time_lookup <- stats::setNames(seq_along(time_levels), as.character(time_levels))

  fit_data$.time_idx <- unname(time_lookup[as.character(fit_data$time)])

  fit_data$.row_id <- seq_len(n)


  fit_data_ordered <- fit_data[order(fit_data$subject, fit_data$.time_idx), , drop = FALSE]

  split_rows <- split(fit_data_ordered, fit_data_ordered$subject)


  for (j in seq_len(nsim)) {

    u <- rep(NA_real_, n)

    for (subject_rows in split_rows) {

      row_ids <- subject_rows$.row_id

      if (length(row_ids) == 0L) next

      u[row_ids[[1L]]] <- stats::runif(1L)

      if (length(row_ids) > 1L) {

        for (k in 2:length(row_ids)) {

          left_row <- row_ids[[k - 1L]]

          current_row <- row_ids[[k]]

          if (fit_data$.time_idx[[current_row]] != fit_data$.time_idx[[left_row]] + 1L) {

            u[current_row] <- stats::runif(1L)

            next

          }

          theta <- fit_data$theta_fit[[left_row]]

          zeta <- fit_data$zeta_fit[[left_row]]

          target <- stats::runif(1L)

          if (!is.finite(theta)) {

            u[current_row] <- target

          } else {

            u[current_row] <- .sim_invert_hfunc1(

              u1 = u[left_row],

              target = target,

              family = object$copula_dist,

              par = theta,

              par2 = if (is.finite(zeta)) zeta else 0

            )

          }

        }

      }

    }

    args <- c(list(p = u), diag_data$params)

    args <- args[names(args) %in% formalArgs(qfun)]

    out[, j] <- do.call(qfun, args)

  }

  out

}


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

  ...

) {

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


