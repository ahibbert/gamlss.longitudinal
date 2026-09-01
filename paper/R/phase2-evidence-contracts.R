jss_phase2_attempt_key <- function(x) {
  key <- intersect(c("study_id", "scenario_id", "method", "attempt_id"), names(x))
  if (!all(c("scenario_id", "attempt_id") %in% key)) {
    stop("Attempt data require scenario_id and attempt_id.", call. = FALSE)
  }
  key
}

jss_phase2_validate_attempts <- function(attempts) {
  required <- c("scenario_id", "attempt_id", "seed", "attempted", "converged", "retained", "failure_reason")
  missing <- setdiff(required, names(attempts))
  if (length(missing)) stop("Attempt data are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(attempts)) stop("Attempt data must contain at least one row.", call. = FALSE)
  if (anyNA(attempts$scenario_id) || any(!nzchar(as.character(attempts$scenario_id)))) {
    stop("Every attempt requires a scenario_id.", call. = FALSE)
  }
  if (anyNA(attempts$attempt_id) || anyNA(attempts$seed)) stop("Every attempt requires an attempt ID and seed.", call. = FALSE)
  if (any(attempts$attempted %in% FALSE) || anyNA(attempts$attempted)) stop("An attempt row must have attempted = TRUE.", call. = FALSE)
  if (any(attempts$retained %in% TRUE & !(attempts$converged %in% TRUE))) {
    stop("A retained attempt must be converged.", call. = FALSE)
  }
  failed <- !(attempts$retained %in% TRUE)
  if (any(is.na(attempts$failure_reason[failed]) | !nzchar(as.character(attempts$failure_reason[failed])) | attempts$failure_reason[failed] == "none")) {
    stop("Every non-retained attempt requires a named failure reason.", call. = FALSE)
  }
  key <- jss_phase2_attempt_key(attempts)
  signature <- do.call(paste, c(attempts[key], sep = "\r"))
  if (anyDuplicated(signature)) stop("Attempt keys are not unique.", call. = FALSE)
  invisible(TRUE)
}

jss_phase2_denominators <- function(attempts) {
  jss_phase2_validate_attempts(attempts)
  groups <- intersect(c("study_id", "scenario_id", "method"), names(attempts))
  split_key <- interaction(attempts[groups], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(attempts, split_key), function(x) {
    failed <- !(x$retained %in% TRUE)
    first <- x[1L, groups, drop = FALSE]
    reasons <- sort(table(as.character(x$failure_reason[failed])), decreasing = TRUE)
    cbind(
      first,
      data.frame(
        attempted = nrow(x),
        converged = sum(x$converged %in% TRUE),
        retained = sum(x$retained %in% TRUE),
        failed = sum(failed),
        failure_reasons = paste(names(reasons), as.integer(reasons), sep = ":", collapse = ";"),
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_phase2_wilson_interval <- function(successes, attempts, level = 0.95) {
  if (length(successes) != length(attempts) || any(attempts < 1) || any(successes < 0 | successes > attempts)) {
    stop("Invalid binomial counts.", call. = FALSE)
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- successes / attempts
  denom <- 1 + z^2 / attempts
  centre <- (p + z^2 / (2 * attempts)) / denom
  half <- z * sqrt(p * (1 - p) / attempts + z^2 / (4 * attempts^2)) / denom
  data.frame(estimate = p, mcse = sqrt(p * (1 - p) / attempts), lower = pmax(0, centre - half), upper = pmin(1, centre + half))
}

jss_phase2_mean_interval <- function(x, level = 0.95) {
  x <- as.numeric(x[is.finite(x)])
  n <- length(x)
  if (!n) return(data.frame(n = 0L, estimate = NA_real_, mcse = NA_real_, lower = NA_real_, upper = NA_real_))
  estimate <- mean(x)
  mcse <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
  half <- if (n > 1L) stats::qt(1 - (1 - level) / 2, df = n - 1L) * mcse else NA_real_
  data.frame(n = n, estimate = estimate, mcse = mcse, lower = estimate - half, upper = estimate + half)
}

jss_phase2_paired_difference <- function(metrics, lhs_method, rhs_method, metric_col = "value", level = 0.95) {
  required <- c("scenario_id", "attempt_id", "method", metric_col)
  missing <- setdiff(required, names(metrics))
  if (length(missing)) stop("Paired metrics are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  key <- intersect(c("study_id", "scenario_id", "attempt_id", "seed", "metric"), names(metrics))
  lhs <- metrics[metrics$method == lhs_method, c(key, metric_col), drop = FALSE]
  rhs <- metrics[metrics$method == rhs_method, c(key, metric_col), drop = FALSE]
  names(lhs)[ncol(lhs)] <- "lhs"
  names(rhs)[ncol(rhs)] <- "rhs"
  paired <- merge(lhs, rhs, by = key, all = FALSE, sort = TRUE)
  paired$difference <- as.numeric(paired$lhs) - as.numeric(paired$rhs)
  groups <- intersect(c("study_id", "scenario_id", "metric"), names(paired))
  split_key <- interaction(paired[groups], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(paired, split_key), function(x) {
    interval <- jss_phase2_mean_interval(x$difference, level)
    cbind(
      x[1L, groups, drop = FALSE],
      data.frame(
        lhs_method = lhs_method,
        rhs_method = rhs_method,
        paired_attempts = interval$n,
        estimate = interval$estimate,
        mcse = interval$mcse,
        lower = interval$lower,
        upper = interval$upper,
        sign_probability = mean(x$difference > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  })
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_phase2_reconcile_scenarios <- function(attempts, scenarios, design_fields) {
  jss_phase2_validate_attempts(attempts)
  required <- c("scenario_id", design_fields)
  missing <- setdiff(required, names(scenarios))
  if (length(missing)) stop("Scenario table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(scenarios$scenario_id)) stop("Scenario table has duplicate scenario IDs.", call. = FALSE)
  if (!setequal(unique(as.character(attempts$scenario_id)), as.character(scenarios$scenario_id))) {
    stop("Attempt and scenario IDs do not agree.", call. = FALSE)
  }
  for (scenario_id in scenarios$scenario_id) {
    a <- attempts[attempts$scenario_id == scenario_id, , drop = FALSE]
    s <- scenarios[scenarios$scenario_id == scenario_id, , drop = FALSE]
    for (field in design_fields) {
      if (!field %in% names(a)) stop("Attempt data do not contain design field: ", field, call. = FALSE)
      values <- unique(as.character(a[[field]]))
      if (length(values) != 1L || !identical(values, as.character(s[[field]][[1L]]))) {
        stop("Scenario metadata mismatch for ", scenario_id, " / ", field, ".", call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

jss_phase2_reconcile_denominators <- function(attempts, summary) {
  observed <- jss_phase2_denominators(attempts)
  key <- intersect(c("study_id", "scenario_id", "method"), names(observed))
  required <- c(key, "attempted", "converged", "retained", "failed")
  missing <- setdiff(required, names(summary))
  if (length(missing)) stop("Summary is missing denominator fields: ", paste(missing, collapse = ", "), call. = FALSE)
  comparison <- merge(observed[required], summary[required], by = key, all = TRUE, suffixes = c("_attempt", "_summary"))
  for (field in c("attempted", "converged", "retained", "failed")) {
    if (anyNA(comparison[[paste0(field, "_attempt")]]) || anyNA(comparison[[paste0(field, "_summary")]]) ||
        any(comparison[[paste0(field, "_attempt")]] != comparison[[paste0(field, "_summary")]])) {
      stop("Summary denominator mismatch for ", field, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

jss_phase2_validate_claims <- function(claims, effects) {
  claim_fields <- c("claim_id", "study_id", "scenario_id", "metric", "lhs_method", "rhs_method", "expected_direction", "attempt_source")
  missing <- setdiff(claim_fields, names(claims))
  if (length(missing)) stop("Claim register is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(claims$claim_id)) stop("Claim IDs must be unique.", call. = FALSE)
  if (any(!claims$expected_direction %in% c("positive", "negative", "no_direction"))) {
    stop("Claims use an unregistered direction.", call. = FALSE)
  }
  if (any(!nzchar(as.character(claims$attempt_source)))) stop("Every claim must name its attempt-level source.", call. = FALSE)
  effect_key <- c("study_id", "scenario_id", "metric", "lhs_method", "rhs_method")
  missing_effect <- setdiff(c(effect_key, "estimate", "mcse", "lower", "upper", "paired_attempts"), names(effects))
  if (length(missing_effect)) stop("Effect table is missing: ", paste(missing_effect, collapse = ", "), call. = FALSE)
  matched <- merge(claims, effects, by = effect_key, all.x = TRUE, sort = FALSE)
  if (nrow(matched) != nrow(claims) || anyNA(matched$estimate) || any(matched$paired_attempts < 1L)) {
    stop("Every empirical claim must map to exactly one estimated effect.", call. = FALSE)
  }
  invalid_interval <- !is.finite(matched$estimate) | !is.finite(matched$mcse) | matched$mcse < 0 |
    !is.finite(matched$lower) | !is.finite(matched$upper) |
    matched$lower > matched$estimate | matched$estimate > matched$upper
  if (any(invalid_interval)) stop("One or more claims has invalid Monte Carlo uncertainty.", call. = FALSE)
  unsupported <- (matched$expected_direction == "positive" & matched$lower <= 0) |
    (matched$expected_direction == "negative" & matched$upper >= 0)
  if (any(unsupported)) stop("One or more narrative directions lacks confidence-interval support from attempt-level evidence.", call. = FALSE)
  matched[, c("claim_id", effect_key, "expected_direction", "attempt_source", "paired_attempts", "estimate", "mcse", "lower", "upper"), drop = FALSE]
}
