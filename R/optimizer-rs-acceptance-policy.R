#' Accept, backtrack, or reject one RS backfitting proposal
#'
#' Requires a valid finite likelihood and a nondecreasing supplied objective.
#' The registered RS caller supplies the conditional penalized objective at a
#' frozen lambda. Backtracking halves the step; disabling it permits one try,
#' not a downhill acceptance. Raw likelihood values remain separate diagnostics.
#' 
#' If no proposal is accepted, the current parameter values are retained and the step is rejected.
#'
#' @noRd
.gl_rs_accept_backfitting_step <- function(
    proposed_results,
    current_results,
    nominal_step_size,
    use_backtracking,
    backtracking_max_halves,
    proposal_fn,
    outer_iteration,
    inner_iteration,
    global_inner_iteration,
    parameter,
    elapsed_sec = NA_real_,
    objective_fn = function(result) as.numeric(result$calc_lik_out_end$log_lik["joint"])) {
  start_joint_loglik <- as.numeric(current_results$calc_lik_out_end$log_lik["joint"])
  proposed_joint_loglik <- as.numeric(proposed_results$calc_lik_out_end$log_lik["joint"])
  finite_scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  start_objective <- objective_fn(current_results)
  proposed_objective <- objective_fn(proposed_results)
  if (!finite_scalar(start_joint_loglik) || !finite_scalar(start_objective) ||
      identical(current_results$calc_lik_out_end$valid, FALSE)) {
    stop("Invalid RS baseline likelihood or conditional objective.", call. = FALSE)
  }

  accepted_results <- proposed_results
  accepted_step_size <- nominal_step_size
  step_rejected <- FALSE
  backtracking_attempts_used <- 0L
  max_backtrack <- 0L
  backtracking_applied <- FALSE
  proposed_invalid <- !finite_scalar(proposed_joint_loglik) || !finite_scalar(proposed_objective) ||
    identical(proposed_results$calc_lik_out_end$valid, FALSE)
  proposed_decrease <- !proposed_invalid && proposed_objective < start_objective

  if (isTRUE(use_backtracking) &&
    is.finite(start_joint_loglik) &&
    (isTRUE(proposed_invalid) || proposed_decrease)) {
    backtracking_applied <- TRUE
    max_backtrack <- as.integer(backtracking_max_halves)
    trial_step <- nominal_step_size
    accepted <- FALSE

    for (bt in seq_len(max_backtrack)) {
      backtracking_attempts_used <- bt
      trial_step <- trial_step / 2
      trial_results <- proposal_fn(trial_step)
      trial_joint_loglik <- as.numeric(trial_results$calc_lik_out_end$log_lik["joint"])
      trial_objective <- objective_fn(trial_results)
      trial_valid <- !identical(trial_results$calc_lik_out_end$valid, FALSE)

      if (isTRUE(trial_valid) && finite_scalar(trial_joint_loglik) &&
          finite_scalar(trial_objective) && trial_objective >= start_objective) {
        accepted_results <- trial_results
        accepted_step_size <- trial_step
        accepted <- TRUE
        break
      }
    }

    if (!accepted) {
      accepted_results <- current_results
      accepted_step_size <- 0
      step_rejected <- TRUE
    }
  } else if (isTRUE(proposed_invalid) || proposed_decrease) {
    accepted_results <- current_results
    accepted_step_size <- 0
    step_rejected <- TRUE
  }

  accepted_joint_loglik <- as.numeric(accepted_results$calc_lik_out_end$log_lik["joint"])
  accepted_objective <- objective_fn(accepted_results)
  trace_row <- data.frame(
    outer_iteration = as.integer(outer_iteration),
    inner_iteration = as.integer(inner_iteration),
    global_inner_iteration = as.integer(global_inner_iteration),
    parameter = parameter,
    start_logLik = as.numeric(start_joint_loglik),
    proposed_logLik = as.numeric(proposed_joint_loglik),
    accepted_logLik = as.numeric(accepted_joint_loglik),
    proposed_change = as.numeric(proposed_joint_loglik - start_joint_loglik),
    accepted_change = as.numeric(accepted_joint_loglik - start_joint_loglik),
    nominal_step_size = as.numeric(nominal_step_size),
    accepted_step_size = as.numeric(accepted_step_size),
    backtracking_attempts = as.integer(backtracking_attempts_used),
    max_backtracking_attempts = as.integer(max_backtrack),
    rejected = isTRUE(step_rejected),
    proposed_likelihood_valid = !isTRUE(proposed_invalid),
    start_objective = start_objective,
    proposed_objective = proposed_objective,
    accepted_objective = accepted_objective,
    proposed_objective_change = proposed_objective - start_objective,
    accepted_objective_change = accepted_objective - start_objective,
    elapsed_sec = as.numeric(elapsed_sec),
    stringsAsFactors = FALSE
  )

  list(
    accepted_results = accepted_results,
    accepted_step_size = accepted_step_size,
    rejected = step_rejected,
    backtracking_applied = backtracking_applied,
    backtracking_attempts = backtracking_attempts_used,
    max_backtracking_attempts = max_backtrack,
    start_joint_loglik = start_joint_loglik,
    proposed_joint_loglik = proposed_joint_loglik,
    accepted_joint_loglik = accepted_joint_loglik,
    start_objective = start_objective,
    proposed_objective = proposed_objective,
    accepted_objective = accepted_objective,
    change_objective = accepted_objective - start_objective,
    trace_row = trace_row
  )
}

#' Report RS acceptance and backtracking diagnostics
#'
#' @noRd
.gl_report_rs_acceptance <- function(
    par_name,
    step_size,
    use_backtracking,
    rs_acceptance,
    calc_lik_out_end,
    verbose,
    cat_fn = cat,
    print_fn = print) {
  if (isTRUE(rs_acceptance$backtracking_applied) && verbose > 1) {
    cat_fn(paste0(
      "\nBacktracking applied for ", par_name,
      ": step_size ", signif(step_size, 4),
      " -> ", signif(rs_acceptance$accepted_step_size, 4),
      " (halves tried=", rs_acceptance$backtracking_attempts,
      "/", rs_acceptance$max_backtracking_attempts, ")",
      "\n"
    ))
  }

  if (par_name == "theta" && verbose > 2) {
    cat_fn(paste0(
      "\nTheta step diagnostics: start=", signif(rs_acceptance$start_joint_loglik, 8),
      ", proposed=", signif(rs_acceptance$proposed_joint_loglik, 8),
      ", accepted=", signif(rs_acceptance$accepted_joint_loglik, 8),
      ", backtracking=", if (isTRUE(use_backtracking)) "on" else "off",
      ", step=", signif(step_size, 4),
      ", accepted_step=", signif(rs_acceptance$accepted_step_size, 4),
      ", halves_tried=", rs_acceptance$backtracking_attempts,
      "/", rs_acceptance$max_backtracking_attempts,
      ", rejected=", if (isTRUE(rs_acceptance$rejected)) "yes" else "no",
      "\n"
    ))
  }

  if (verbose > 2) {
    cat_fn("\nLogLik:\n")
    print_fn(calc_lik_out_end$log_lik)
  }

  invisible(TRUE)
}

#' Apply one accepted RS step to loop bookkeeping state
#'
#' @noRd
.gl_apply_rs_acceptance_state <- function(
    rs_acceptance,
    rs_block_trace,
    calc_lik_out,
    run_counter,
    outer_run_counter,
    inner_run_counter) {
  accepted_results <- rs_acceptance$accepted_results
  rs_block_trace[[length(rs_block_trace) + 1L]] <- rs_acceptance$trace_row
  calc_lik_out_end <- accepted_results$calc_lik_out_end
  change_log_lik <- calc_lik_out_end$log_lik["joint"] - calc_lik_out$log_lik["joint"]

  list(
    par_cov = accepted_results$par_cov,
    par_s = accepted_results$par_s,
    calc_lik_out_end = calc_lik_out_end,
    df_s = accepted_results$df_s,
    rs_block_trace = rs_block_trace,
    change_log_lik = change_log_lik,
    change_objective = rs_acceptance$change_objective %||% change_log_lik,
    coefficient_step_status = if (isTRUE(rs_acceptance$rejected)) "rejected" else "accepted",
    run_counter = run_counter + 1,
    outer_run_counter = outer_run_counter + 1,
    inner_run_counter = inner_run_counter + 1
  )
}
