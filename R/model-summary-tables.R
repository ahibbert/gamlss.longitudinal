#' Build the summary coefficient table
#'
#' @noRd
.gl_summary_coefficient_table <- function(object, vcov_out = NULL) {
  coef_tbl <- data.frame(
    term = names(object$par),
    estimate = as.numeric(object$par),
    std_error = NA_real_,
    p_value = NA_real_,
    signif = NA_character_,
    stringsAsFactors = FALSE
  )

  coef_tbl$.original_order <- seq_len(nrow(coef_tbl))
  coef_tbl$parameter <- sub("\\..*$", "", coef_tbl$term)

  param_order <- c("mu", "sigma", "nu", "tau", "theta", "zeta")
  coef_tbl$.param_rank <- match(coef_tbl$parameter, param_order)
  coef_tbl$.param_rank[is.na(coef_tbl$.param_rank)] <- length(param_order) + 1L

  if (!is.null(vcov_out) && !is.null(vcov_out$vcov) && !is.null(vcov_out$vcov$overall)) {
    V <- vcov_out$vcov$overall
    se <- NULL
    if (!is.null(vcov_out$se) && !is.null(vcov_out$se$overall)) {
      se <- as.numeric(vcov_out$se$overall)
      se_names <- names(vcov_out$se$overall)
    } else {
      se <- sqrt(pmax(0, diag(V)))
      se_names <- names(diag(V))
    }

    if (is.null(se_names) && !is.null(rownames(V)) && length(rownames(V)) == length(se)) {
      se_names <- rownames(V)
    }

    if (!is.null(se_names)) {
      names(se) <- se_names
    }

    if (!is.null(names(se))) {
      idx <- match(coef_tbl$term, names(se))
      coef_tbl$std_error <- se[idx]
    } else if (length(se) == nrow(coef_tbl)) {
      coef_tbl$std_error <- se
    }

    z_abs <- abs(coef_tbl$estimate / coef_tbl$std_error)
    coef_tbl$p_value <- 2 * stats::pnorm(z_abs, lower.tail = FALSE)
    coef_tbl$signif <- ifelse(
      is.na(coef_tbl$p_value),
      NA_character_,
      ifelse(coef_tbl$p_value < 0.001, "***",
        ifelse(coef_tbl$p_value < 0.01, "**",
          ifelse(coef_tbl$p_value < 0.05, "*",
            ifelse(coef_tbl$p_value < 0.1, ".", " ")
          )
        )
      )
    )
  }

  coef_tbl <- coef_tbl[order(coef_tbl$.param_rank, coef_tbl$.original_order), , drop = FALSE]
  rownames(coef_tbl) <- NULL

  within(coef_tbl, {
    .original_order <- NULL
    .param_rank <- NULL
  })
}

#' Build the smooth-term table for summary output
#'
#' @noRd
.gl_summary_smooth_terms <- function(object) {
  st <- list()
  if (!is.null(object$par_s) && length(object$par_s) > 0) {
    for (par_name in names(object$par_s)) {
      if (length(object$par_s[[par_name]]) == 0) next
      for (s_name in names(object$par_s[[par_name]])) {
        smooth_edf <- NA_real_
        if (!is.null(object$df_s) && !is.null(object$df_s[[par_name]]) && s_name %in% names(object$df_s[[par_name]])) {
          smooth_edf <- suppressWarnings(as.numeric(object$df_s[[par_name]][[s_name]])[1])
        }
        st[[length(st) + 1]] <- data.frame(
          parameter = par_name,
          smooth_term = s_name,
          edf = smooth_edf,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(st) == 0) {
    data.frame(parameter = character(0), smooth_term = character(0), edf = numeric(0), stringsAsFactors = FALSE)
  } else {
    do.call(rbind, st)
  }
}
