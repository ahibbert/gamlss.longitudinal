#' Check whether a stopping criterion is set to automatic (NULL or NA)
#'
#' @noRd
.gl_is_auto_stop_crit <- function(x) {
  is.null(x) || (length(x) == 1 && is.na(x))
}


#' Validate a positive and finite stopping criterion
#'
#' @noRd
.gl_validate_stop_crit <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    stop(name, " must be a single positive finite number, or NA/NULL for automatic selection.")
  }
  as.numeric(x)
}

#' Set stopping criteria for inner and outer loops to user specified values OR
#' automatically determine them based on the initial likelihood value and dataset size.
#' 
#' The calculation is based on the initial joint log-likelihood, with the idea that 
#' the outer stopping criterion should be a small fraction of the initial log-likelihood (or dataset size), 
#' and the inner stopping criterion should be a smaller fraction of that. 
#' The function also checks that the stopping criteria are not too small or too large by setting reasonable bounds.
#'
#' @noRd
.gl_resolve_stop_criteria <- function(
    inner_stop_crit,
    outer_stop_crit,
    cg_grad_tol,
    cg_step_tol,
    method,
    par_cov,
    par_s,
    mm,
    margin_dist,
    copula_dist,
    copula_link,
    dataset,
    pair_cache,
    margin_eval_cache,
    verbose = 1) {
  if (.gl_is_auto_stop_crit(inner_stop_crit) || .gl_is_auto_stop_crit(outer_stop_crit)) {
    eta_init_out <- calc_eta(par_cov, mm, margin_dist, copula_link, par_s = par_s)
    eta_init <- eta_init_out$eta_inv
    calc_lik_init <- calc_likelihood_minimal(
      eta_init,
      mm = mm$x,
      margin_dist,
      copula_dist,
      calc_d2 = FALSE,
      response = dataset$response,
      response_margin = dataset$time,
      response_subject = dataset$subject,
      pair_cache = pair_cache,
      margin_eval_cache = margin_eval_cache
    )

    init_joint_ll <- as.numeric(calc_lik_init$log_lik["joint"])
    if (!is.finite(init_joint_ll)) {
      init_joint_ll <- 0
    }

    scale_base <- max(1, abs(init_joint_ll), nrow(dataset))
    auto_outer_stop_crit <- min(0.05, max(1e-4, 1e-6 * scale_base))
    if (identical(method, "CG")) {
      auto_outer_stop_crit <- auto_outer_stop_crit / 10
    }
    auto_inner_stop_crit <- min(0.01, max(1e-5, auto_outer_stop_crit / 5))

    if (.gl_is_auto_stop_crit(outer_stop_crit)) {
      outer_stop_crit <- auto_outer_stop_crit
    } else {
      outer_stop_crit <- .gl_validate_stop_crit(outer_stop_crit, "outer_stop_crit")
    }

    if (.gl_is_auto_stop_crit(inner_stop_crit)) {
      inner_stop_crit <- auto_inner_stop_crit
    } else {
      inner_stop_crit <- .gl_validate_stop_crit(inner_stop_crit, "inner_stop_crit")
    }

    if (verbose > 0) {
      cat(
        "\nUsing stop criteria:",
        "inner_stop_crit=", format(inner_stop_crit, digits = 6),
        "| outer_stop_crit=", format(outer_stop_crit, digits = 6), "\n"
      )
    }
  } else {
    inner_stop_crit <- .gl_validate_stop_crit(inner_stop_crit, "inner_stop_crit")
    outer_stop_crit <- .gl_validate_stop_crit(outer_stop_crit, "outer_stop_crit")
  }

  cg_grad_tol_eff <- if (.gl_is_auto_stop_crit(cg_grad_tol)) {
    max(1e-3, 10 * outer_stop_crit)
  } else {
    .gl_validate_stop_crit(cg_grad_tol, "cg_grad_tol")
  }

  cg_step_tol_eff <- if (.gl_is_auto_stop_crit(cg_step_tol)) {
    max(1e-5, 0.1 * outer_stop_crit)
  } else {
    .gl_validate_stop_crit(cg_step_tol, "cg_step_tol")
  }

  list(
    inner_stop_crit = inner_stop_crit,
    outer_stop_crit = outer_stop_crit,
    cg_grad_tol_eff = cg_grad_tol_eff,
    cg_step_tol_eff = cg_step_tol_eff
  )
}
