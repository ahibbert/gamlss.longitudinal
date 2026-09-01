#' Wald tests for fixed coefficients

#'

#' `wald_test()` provides a small reporting-friendly hypothesis-test surface for

#' fitted `gamlss.longitudinal` models. It uses the same variance-covariance

#' route as [summary.gamlss.longitudinal()] and [confint.gamlss.longitudinal()],

#' so numerical-Hessian tests should be reported as approximate.

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param terms Optional coefficient names, formula-term names such as

#'   `"mu.treatment"`, coefficient-name prefixes, or numeric indices. When `L`

#'   is `NULL`, these select coefficients for individual tests or a joint test.

#' @param L Optional contrast matrix. Columns must either be named with

#'   coefficient names or have one column per fixed coefficient in model order.

#' @param rhs Null-hypothesis value. Either a scalar or one value per tested row.

#' @param joint Logical; when `TRUE`, test selected `terms` jointly. Contrast

#'   matrices supplied through `L` are always tested jointly.

#' @param method Variance-covariance method passed to [vcov.gamlss.longitudinal()].

#' @param ... Additional arguments passed to [vcov.gamlss.longitudinal()].

#'

#' @return An object of class `gamlss_longitudinal_wald_test`.

#' @export

wald_test <- function(
    object,
    terms = NULL,
    L = NULL,
    rhs = 0,
    joint = FALSE,
    method = "analytical",
    ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  estimates <- stats::coef(object)

  if (length(estimates) == 0L) {
    stop("No fixed coefficients are available to test.", call. = FALSE)
  }

  vc <- .resolve_vcov(object, extra_args = list(method = method, ...))
  .gl_require_available_inference(vc)

  V <- vc$vcov$overall

  if (is.null(V) || !is.matrix(V)) {
    stop("A fixed-effect variance-covariance matrix is required for Wald tests.", call. = FALSE)
  }

  V_names <- colnames(V) %||% rownames(V)

  if (is.null(V_names)) {
    V_names <- names(estimates)
  }

  idx_v <- match(names(estimates), V_names)

  if (any(is.na(idx_v))) {
    stop("Variance-covariance matrix names do not match fitted coefficients.", call. = FALSE)
  }

  V <- V[idx_v, idx_v, drop = FALSE]

  rownames(V) <- colnames(V) <- names(estimates)

  contrast <- .gl_wald_contrast_matrix(
    L = L,
    terms = terms,
    estimates = estimates,
    object = object,
    joint = joint
  )
  L <- contrast$L
  joint <- contrast$joint

  rhs <- rep(rhs, length.out = nrow(L))

  estimate <- as.numeric(L %*% estimates)

  diff <- estimate - rhs

  LVL <- L %*% V %*% t(L)

  if (isTRUE(joint)) {
    stat <- tryCatch(

      as.numeric(t(diff) %*% solve(LVL, diff)),
      error = function(e) NA_real_
    )

    out <- data.frame(
      hypothesis = paste(rownames(L), collapse = ", "),
      df = nrow(L),
      statistic = stat,
      p_value = stats::pchisq(stat, df = nrow(L), lower.tail = FALSE),
      method = vc$method %||% method,
      stringsAsFactors = FALSE
    )
  } else {
    se <- .gl_sqrt_derived_variance(
      diag(LVL), "Wald contrast covariance", allow_zero = FALSE
    )

    z <- diff / se

    out <- data.frame(
      term = rownames(L),
      estimate = estimate,
      rhs = rhs,
      std_error = se,
      statistic = z,
      p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      method = vc$method %||% method,
      stringsAsFactors = FALSE
    )
  }

  attr(out, "joint") <- isTRUE(joint)

  attr(out, "method_requested") <- method

  attr(out, "vcov_method") <- vc$method %||% method

  class(out) <- c("gamlss_longitudinal_wald_test", "data.frame")

  out
}

#' @export

print.gamlss_longitudinal_wald_test <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nWald Test for gamlss.longitudinal\n")

  cat("---------------------------------\n")

  cat("Test type:", if (isTRUE(attr(x, "joint"))) "joint" else "individual", "\n")

  cat("VCOV method:", attr(x, "vcov_method") %||% "unknown", "\n\n")

  print.data.frame(x, digits = digits, row.names = FALSE, ...)

  invisible(x)
}
