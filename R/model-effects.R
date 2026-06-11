#' Counterfactual marginal effects from fitted distributional parameters

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param newdata Data used as the counterfactual baseline.

#' @param variable Single variable to vary.

#' @param values Values to assign to `variable`. Defaults to observed factor

#'   levels for factors/characters or the 25th, 50th, and 75th percentiles for

#'   numeric variables.

#' @param parameter Distributional parameter to summarize, usually `"mu"`.

#' @param reference Optional reference value. Defaults to the first value.

#' @param se.fit Logical; when `TRUE` and `parameter = "mu"`, attach

#'   approximate delta-method standard errors for response-scale averages.

#' @param level Confidence level used when `se.fit = TRUE`.

#' @param vcov_method Variance-covariance method passed to [vcov.gamlss.longitudinal()].

#' @param ... Additional arguments passed to [predict.gamlss.longitudinal()].

#'

#' @return A data frame with average fitted parameter values and contrasts.

#' @importFrom stats predict

#' @export

marginal_effects <- function(

  object,

  newdata,

  variable,

  values = NULL,

  parameter = "mu",

  reference = NULL,

  se.fit = FALSE,

  level = 0.95,

  vcov_method = "analytical",

  ...

) {

  if (!inherits(object, "gamlss.longitudinal")) {

    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)

  }

  if (missing(newdata) || is.null(newdata)) {

    stop("'newdata' is required for counterfactual marginal effects.", call. = FALSE)

  }

  newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)

  if (!is.character(variable) || length(variable) != 1L || !variable %in% names(newdata)) {

    stop("'variable' must be a single column name in 'newdata'.", call. = FALSE)

  }


  x <- newdata[[variable]]

  if (is.null(values)) {

    values <- if (is.factor(x)) {

      levels(x)

    } else if (is.character(x)) {

      sort(unique(x))

    } else {

      as.numeric(stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE))

    }

  }

  if (length(values) == 0L) {

    stop("'values' must contain at least one counterfactual value.", call. = FALSE)

  }

  if (is.null(reference)) {

    reference <- values[[1L]]

  }


  baseline_n <- nrow(newdata)

  add_factor_calibration_rows <- function(nd, template) {

    factor_cols <- names(nd)[vapply(nd, is.factor, logical(1))]

    extras <- list()

    for (fc in factor_cols) {

      levs <- levels(nd[[fc]])

      if (length(levs) < 2L) next

      present <- unique(as.character(nd[[fc]]))

      missing <- setdiff(levs, present)

      if (length(missing) == 0L) next

      for (lev in missing) {

        row <- template[1L, , drop = FALSE]

        row[[fc]] <- factor(lev, levels = levs, ordered = is.ordered(nd[[fc]]))

        for (other_fc in factor_cols) {

          if (other_fc == fc) next

          row[[other_fc]] <- factor(

            as.character(row[[other_fc]]),

            levels = levels(nd[[other_fc]]),

            ordered = is.ordered(nd[[other_fc]])

          )

        }

        extras[[length(extras) + 1L]] <- row

      }

    }

    if (length(extras) == 0L) {

      return(nd)

    }

    rbind(nd, do.call(rbind, extras))

  }


  rows <- lapply(values, function(value) {

    nd <- newdata

    if (is.factor(nd[[variable]])) {

      nd[[variable]] <- factor(

        rep(value, nrow(nd)),

        levels = levels(nd[[variable]]),

        ordered = is.ordered(nd[[variable]])

      )

    } else {

      nd[[variable]] <- value

    }

    nd_pred <- add_factor_calibration_rows(nd, newdata)

    if (isTRUE(se.fit) && identical(parameter, "mu")) {

      pred <- predict(

        object,

        newdata = nd_pred,

        type = "response",

        se.fit = TRUE,

        interval = "none",

        vcov_method = vcov_method,

        ...

      )

      pred <- pred[seq_len(baseline_n), , drop = FALSE]

      estimate <- mean(pred$fit, na.rm = TRUE)

      n_finite <- sum(is.finite(pred$se.fit))

      std_error <- if (n_finite > 0L) {

        sqrt(sum(pred$se.fit^2, na.rm = TRUE)) / n_finite

      } else {

        NA_real_

      }

    } else {

      pred <- predict(object, newdata = nd_pred, type = "parameters", ...)

      pred <- pred[seq_len(baseline_n), , drop = FALSE]

      if (!parameter %in% names(pred)) {

        stop("Parameter '", parameter, "' is not available in model predictions.", call. = FALSE)

      }

      estimate <- mean(pred[[parameter]], na.rm = TRUE)

      std_error <- NA_real_

    }

    data.frame(

      variable = variable,

      value = as.character(value),

      parameter = parameter,

      estimate = estimate,

      std_error = std_error,

      stringsAsFactors = FALSE

    )

  })

  out <- do.call(rbind, rows)

  ref_idx <- match(as.character(reference), out$value)

  if (is.na(ref_idx)) {

    stop("'reference' must be one of the counterfactual values.", call. = FALSE)

  }

  out$reference <- out$value[[ref_idx]]

  out$contrast <- out$estimate - out$estimate[[ref_idx]]

  if (isTRUE(se.fit) && any(is.finite(out$std_error))) {

    alpha <- 1 - level

    z <- stats::qnorm(1 - alpha / 2)

    out$conf.low <- out$estimate - z * out$std_error

    out$conf.high <- out$estimate + z * out$std_error

  }

  rownames(out) <- NULL

  out

}
