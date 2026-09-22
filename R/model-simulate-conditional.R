#' Prepare conditional forward simulation inputs
#'
#' @noRd
.gl_conditional_simulation_setup <- function(object, newdata, condition_on) {
  if (!(identical(condition_on, "fitted") || is.data.frame(condition_on))) {
    stop(
      "'condition_on' must be NULL, \"fitted\", or a data frame containing response history.",
      call. = FALSE
    )
  }

  future <- .gl_prepare_conditional_panel(
    object,
    newdata,
    context = "'newdata'",
    require_response = FALSE
  )
  future$response <- NA_real_
  future$.gl_conditional_future_row <- seq_len(nrow(future))

  # Validate family-specific inputs before history and future are combined, so
  # missing future predictors cannot be hidden by row binding.
  future_only <- .gl_simulation_newdata(object, future)
  schedule <- future_only$time_levels

  if (identical(condition_on, "fitted")) {
    history_source <- "fitted"
    if (is.null(object$dataset)) {
      stop(
        "Fitted response history is unavailable in this model object; supply a history data frame instead.",
        call. = FALSE
      )
    }
    history_all <- .gl_prepare_newdata_internal(object, object$dataset, require_response = FALSE)
  } else {
    history_source <- "data.frame"
    history_all <- .gl_prepare_conditional_panel(
      object,
      condition_on,
      context = "'condition_on'",
      require_response = FALSE
    )
  }

  history_subjects <- unique(as.character(history_all$subject))
  history <- history_all[is.finite(history_all$response), , drop = FALSE]

  if (nrow(history)) {
    history_time <- .gl_conditional_time_values(history)
    history_idx <- match(as.character(history_time), as.character(schedule))
    if (anyNA(history_idx)) {
      stop(
        "Usable history contains time values outside the fitted and supplied future schedule.",
        call. = FALSE
      )
    }
    history$.gl_conditional_time_idx <- history_idx
    history$.gl_conditional_subject_key <- as.character(history$subject)
    history <- .gl_latest_conditional_history(history)
  } else {
    history$.gl_conditional_time_idx <- integer(0)
    history$.gl_conditional_subject_key <- character(0)
  }

  future_time <- .gl_conditional_time_values(future)
  future_idx <- match(as.character(future_time), as.character(schedule))
  if (anyNA(future_idx)) {
    stop("Could not map future time values to the simulation schedule.", call. = FALSE)
  }
  future$.gl_conditional_time_idx <- future_idx
  future$.gl_conditional_subject_key <- as.character(future$subject)

  .gl_validate_future_after_history(future, history)

  history$.gl_conditional_is_history <- TRUE
  future$.gl_conditional_is_history <- FALSE
  history$.gl_conditional_future_row <- NA_integer_

  combined <- .gl_bind_conditional_rows(history, future)
  evaluated <- .gl_simulation_newdata(object, combined)

  n_history <- nrow(history)
  history_rows <- if (n_history) seq_len(n_history) else integer(0)
  future_rows <- n_history + seq_len(nrow(future))

  history_eval <- .gl_conditional_history_evaluation(
    object,
    evaluated,
    history_rows,
    history
  )

  future_diag <- evaluated$diag_data
  future_diag$response <- future_diag$response[future_rows]
  future_diag$subject <- future_diag$subject[future_rows]
  future_diag$time <- future_diag$time[future_rows]
  future_diag$params <- .gl_subset_family_params(future_diag$params, future_rows)
  future_diag$keep_index <- seq_len(nrow(future))

  future_fit <- evaluated$fit_data[future_rows, , drop = FALSE]
  future_fit$.time_idx <- future$.gl_conditional_time_idx
  future_fit$.future_row <- future$.gl_conditional_future_row
  future_fit$.subject_key <- future$.gl_conditional_subject_key

  plan <- .gl_conditional_segment_plan(
    future_fit = future_fit,
    history = history_eval,
    history_subjects = history_subjects,
    schedule = schedule
  )

  list(
    future_diag = future_diag,
    future_fit = future_fit,
    history = history_eval,
    plan = plan,
    time_levels = schedule,
    history_source = history_source
  )
}

#' Prepare and validate one conditional-simulation panel
#'
#' @noRd
.gl_prepare_conditional_panel <- function(object, data, context, require_response) {
  if (!is.data.frame(data)) {
    stop(context, " must be a data frame.", call. = FALSE)
  }
  if (nrow(data) < 1L) {
    stop(context, " must contain at least one row.", call. = FALSE)
  }

  translated_names <- .gl_translate_newdata_names(
    as.data.frame(data, stringsAsFactors = FALSE),
    object$var_map
  )
  missing_keys <- setdiff(c("subject", "time"), names(translated_names))
  if (length(missing_keys)) {
    stop(context, " must contain the fitted subject and time columns.", call. = FALSE)
  }

  nd <- .gl_prepare_newdata_internal(object, data, require_response = require_response)
  required <- unique(unlist(lapply(object$formulas_int, function(fml) {
    if (is.null(fml)) character(0) else all.vars(stats::as.formula(fml))
  }), use.names = FALSE))
  required_predictors <- setdiff(required, "response")
  missing_predictors <- setdiff(required_predictors, names(nd))
  if (length(missing_predictors)) {
    stop(
      context, " is missing model predictor column(s): ",
      paste(missing_predictors, collapse = ", "), ".",
      call. = FALSE
    )
  }
  bad_predictors <- required_predictors[vapply(
    required_predictors,
    function(nm) any(is.na(nd[[nm]])),
    logical(1)
  )]
  if (length(bad_predictors)) {
    stop(
      context, " contains missing values in required predictor column(s): ",
      paste(bad_predictors, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (anyNA(nd$subject) || anyNA(nd$time)) {
    stop(context, " cannot contain missing subject or time values.", call. = FALSE)
  }
  if (anyDuplicated(paste(
    as.character(nd$subject),
    as.character(.gl_conditional_time_values(nd)),
    sep = "\r"
  ))) {
    stop(context, " must contain at most one row per subject/time combination.", call. = FALSE)
  }
  nd
}

#' Return the user-scale scheduled time values from prepared data
#'
#' @noRd
.gl_conditional_time_values <- function(data) {
  if ("time_covariate" %in% names(data) && is.factor(data$time_covariate)) {
    return(data$time_covariate)
  }
  data$time
}

#' Select the most recent usable history response for each subject
#'
#' @noRd
.gl_latest_conditional_history <- function(history) {
  ord <- order(
    history$.gl_conditional_subject_key,
    history$.gl_conditional_time_idx,
    seq_len(nrow(history))
  )
  history <- history[ord, , drop = FALSE]
  keep <- !duplicated(history$.gl_conditional_subject_key, fromLast = TRUE)
  history[keep, , drop = FALSE]
}

#' Require future rows to follow the most recent usable response
#'
#' @noRd
.gl_validate_future_after_history <- function(future, history) {
  if (!nrow(history)) return(invisible(TRUE))
  first_future <- tapply(
    future$.gl_conditional_time_idx,
    future$.gl_conditional_subject_key,
    min
  )
  history_idx <- stats::setNames(
    history$.gl_conditional_time_idx,
    history$.gl_conditional_subject_key
  )
  common <- intersect(names(first_future), names(history_idx))
  invalid <- common[first_future[common] <= history_idx[common]]
  if (length(invalid)) {
    stop(
      "Conditional 'newdata' must contain only visits after each subject's most recent usable history response. Problem subject(s): ",
      paste(utils::head(invalid, 10L), collapse = ", "),
      if (length(invalid) > 10L) ", ..." else "",
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Bind prepared history and future rows while preserving model columns
#'
#' @noRd
.gl_bind_conditional_rows <- function(history, future) {
  all_names <- union(names(history), names(future))
  add_missing <- function(data) {
    for (nm in setdiff(all_names, names(data))) data[[nm]] <- NA
    data[, all_names, drop = FALSE]
  }
  rbind(add_missing(history), add_missing(future))
}

#' Evaluate fitted history CDF intervals and outgoing copula parameters
#'
#' @noRd
.gl_conditional_history_evaluation <- function(object, evaluated, rows, history) {
  if (!length(rows)) {
    return(data.frame(
      subject = character(), time = character(), time_idx = integer(),
      response = numeric(), cdf_lower = numeric(), cdf_upper = numeric(),
      theta_fit = numeric(), zeta_fit = numeric(), stringsAsFactors = FALSE
    ))
  }
  params <- .gl_subset_family_params(evaluated$diag_data$params, rows)
  response <- history$response
  upper <- .gl_call_family_fun("p", evaluated$diag_data$family, response, params)
  if (.is_discrete_margin(object$margin_dist)) {
    lower <- .gl_call_family_fun("p", evaluated$diag_data$family, response - 1, params)
  } else {
    lower <- upper
  }
  eps <- sqrt(.Machine$double.eps)
  lower <- pmin(pmax(as.numeric(lower), eps), 1 - eps)
  upper <- pmin(pmax(as.numeric(upper), eps), 1 - eps)
  if (any(!is.finite(lower) | !is.finite(upper) | upper < lower)) {
    stop("Could not obtain valid fitted CDF values for the conditioning history.", call. = FALSE)
  }

  data.frame(
    subject = history$.gl_conditional_subject_key,
    time = as.character(.gl_conditional_time_values(history)),
    time_idx = history$.gl_conditional_time_idx,
    response = response,
    cdf_lower = lower,
    cdf_upper = upper,
    theta_fit = evaluated$fit_data$theta_fit[rows],
    zeta_fit = evaluated$fit_data$zeta_fit[rows],
    stringsAsFactors = FALSE
  )
}

#' Classify conditional and independent future segment starts
#'
#' @noRd
.gl_conditional_segment_plan <- function(future_fit, history, history_subjects, schedule) {
  starts <- vector("list", 0L)
  groups <- split(seq_len(nrow(future_fit)), future_fit$.subject_key, drop = TRUE)
  history_row <- stats::setNames(seq_len(nrow(history)), history$subject)
  start_number <- 0L

  for (subject_key in names(groups)) {
    rows <- groups[[subject_key]]
    rows <- rows[order(future_fit$.time_idx[rows], future_fit$.future_row[rows])]
    hrow <- unname(history_row[subject_key])

    for (k in seq_along(rows)) {
      row <- rows[[k]]
      if (k == 1L) {
        if (length(hrow) && !is.na(hrow)) {
          adjacent <- future_fit$.time_idx[row] == history$time_idx[hrow] + 1L
          start_type <- if (adjacent) "conditioned_history" else "unconditional_gap"
          reason <- if (adjacent) "adjacent_to_history" else "nonadjacent_to_history"
          previous_time <- history$time[hrow]
        } else {
          known_subject <- subject_key %in% history_subjects
          start_type <- if (known_subject) "unconditional_no_usable_history" else "unconditional_new_subject"
          reason <- if (known_subject) "no_usable_history" else "new_subject"
          previous_time <- NA_character_
        }
      } else {
        previous <- rows[[k - 1L]]
        adjacent <- future_fit$.time_idx[row] == future_fit$.time_idx[previous] + 1L
        if (adjacent) next
        start_type <- "unconditional_gap"
        reason <- "omitted_scheduled_visit"
        previous_time <- as.character(future_fit$time[previous])
      }

      start_number <- start_number + 1L
      starts[[start_number]] <- data.frame(
        subject = subject_key,
        future_row = future_fit$.future_row[row],
        time = as.character(future_fit$time[row]),
        previous_time = previous_time,
        start_type = start_type,
        reason = reason,
        stringsAsFactors = FALSE
      )
    }
  }

  starts <- if (length(starts)) do.call(rbind, starts) else data.frame(
    subject = character(), future_row = integer(), time = character(),
    previous_time = character(), start_type = character(), reason = character()
  )
  rownames(starts) <- NULL

  list(groups = groups, history_row = history_row, starts = starts, schedule = schedule)
}

#' Simulate future trajectories conditionally under the adjacent-pair model
#'
#' @noRd
.gl_simulate_conditional_copula_matrix <- function(object, setup, nsim) {
  future_fit <- setup$future_fit
  history <- setup$history
  n <- nrow(future_fit)
  qfun <- get(
    paste0("q", setup$future_diag$family),
    envir = asNamespace("gamlss.dist"),
    inherits = FALSE
  )
  out <- matrix(NA_real_, nrow = n, ncol = nsim)
  eps <- sqrt(.Machine$double.eps)

  for (j in seq_len(nsim)) {
    u <- rep(NA_real_, n)

    for (subject_key in names(setup$plan$groups)) {
      rows <- setup$plan$groups[[subject_key]]
      rows <- rows[order(future_fit$.time_idx[rows], future_fit$.future_row[rows])]
      hrow <- unname(setup$plan$history_row[subject_key])

      for (k in seq_along(rows)) {
        row <- rows[[k]]
        previous_u <- NA_real_
        theta <- NA_real_
        zeta <- NA_real_
        adjacent <- FALSE

        if (k == 1L && length(hrow) && !is.na(hrow)) {
          adjacent <- future_fit$.time_idx[row] == history$time_idx[hrow] + 1L
          if (adjacent) {
            previous_u <- if (history$cdf_upper[hrow] > history$cdf_lower[hrow]) {
              stats::runif(1L, history$cdf_lower[hrow], history$cdf_upper[hrow])
            } else {
              history$cdf_upper[hrow]
            }
            theta <- history$theta_fit[hrow]
            zeta <- history$zeta_fit[hrow]
          }
        } else if (k > 1L) {
          previous <- rows[[k - 1L]]
          adjacent <- future_fit$.time_idx[row] == future_fit$.time_idx[previous] + 1L
          if (adjacent) {
            previous_u <- u[previous]
            theta <- future_fit$theta_fit[previous]
            zeta <- future_fit$zeta_fit[previous]
          }
        }

        target <- stats::runif(1L)
        if (!adjacent || !is.finite(theta)) {
          u[row] <- target
        } else {
          u[row] <- .sim_invert_hfunc1(
            u1 = min(max(previous_u, eps), 1 - eps),
            target = target,
            family = object$copula_dist,
            par = theta,
            par2 = if (is.finite(zeta)) zeta else 0
          )
        }
      }
    }

    family_args <- attr(setup$future_diag$params, "family_args", exact = TRUE) %||% list()
    args <- c(list(p = u), family_args, setup$future_diag$params)
    args <- args[names(args) %in% formalArgs(qfun)]
    out[, j] <- do.call(qfun, args)
  }

  metadata <- .gl_conditional_simulation_metadata(setup, nsim)
  list(values = out, metadata = metadata)
}

#' Assemble conditioning provenance for simulated trajectories
#'
#' @noRd
.gl_conditional_simulation_metadata <- function(setup, nsim) {
  starts <- setup$plan$starts
  count <- function(type) as.integer(sum(starts$start_type == type))
  out <- list(
    mode = "conditional_forward",
    history_source = setup$history_source,
    conditioning_scope = "most_recent_usable_response_first_order",
    gap_assumption = "independent_segments",
    n_simulations = as.integer(nsim),
    n_future_rows = nrow(setup$future_fit),
    n_subjects = length(unique(setup$future_fit$.subject_key)),
    counts = c(
      conditioned_starts = count("conditioned_history"),
      unconditional_new_subject_starts = count("unconditional_new_subject"),
      unconditional_no_usable_history_starts = count("unconditional_no_usable_history"),
      segmented_gap_starts = count("unconditional_gap")
    ),
    starts = starts
  )
  class(out) <- "gamlss_longitudinal_conditioning_metadata"
  out
}

#' Warn once when conditional simulation encounters scheduled-visit gaps
#'
#' @noRd
.gl_warn_conditional_simulation_gaps <- function(metadata) {
  n_gap <- unname(metadata$counts[["segmented_gap_starts"]])
  if (!is.finite(n_gap) || n_gap < 1L) return(invisible(FALSE))

  warning(structure(
    list(
      message = paste0(
        "Conditional simulation started ", n_gap,
        " future segment", if (n_gap == 1L) "" else "s",
        " unconditionally because scheduled visits were omitted between the usable history or preceding simulated visit and the next future row. ",
        "Under the segmented-independence assumption, dependence is not carried across these gaps."
      ),
      call = NULL,
      n_gap_starts = n_gap,
      starts = metadata$starts[metadata$starts$start_type == "unconditional_gap", , drop = FALSE]
    ),
    class = c(
      "gamlss.longitudinal_conditional_gap_warning",
      "warning",
      "condition"
    )
  ))
  invisible(TRUE)
}
