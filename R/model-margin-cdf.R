#' @keywords internal

#' @noRd

calc_F_x <- function(eta_inv, mm, margin_dist, response) {
  # Setup input matrix of response and parameters

  # margin_names=unique(response_margin)

  # num_margins=length(margin_names)

  margin_deriv_input <- list()

  margin_deriv_input[["y"]] <- response

  margin_deriv_input[["q"]] <- response

  margin_deriv_input[["x"]] <- response

  for (par_name in names(mm)) {
    if (par_name %in% c("mu", "sigma", "nu", "tau")) {
      margin_deriv_input[[par_name]] <- eta_inv[[par_name]]
    }
  }

  fixed_unlinked_values <- attr(margin_dist, "fixed_unlinked_values")

  if (length(fixed_unlinked_values) > 0L) {
    for (par_name in names(fixed_unlinked_values)) {
      value <- fixed_unlinked_values[[par_name]]

      if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
        margin_deriv_input[[par_name]] <- rep(value, length(response))
      }
    }
  }

  margin_deriv_input <- c(
    margin_deriv_input,
    .gl_margin_fixed_family_args(margin_dist, length(response))
  )

  negative_response <- is.finite(response) & response < 0

  if (.is_discrete_margin(margin_dist) && any(negative_response)) {
    margin_p <- rep(NA_real_, length(response))

    margin_p[negative_response] <- 0

    valid_response <- !negative_response

    if (any(valid_response)) {
      call_input <- lapply(margin_deriv_input, function(value) {
        if (length(value) == length(response)) value[valid_response] else value
      })

      margin_pFUN <- get(

        paste("p", margin_dist$family[1], sep = ""),
        envir = asNamespace("gamlss.dist"),
        mode = "function",
        inherits = FALSE
      )

      FUN_args <- names(call_input)[names(call_input) %in% formalArgs(margin_pFUN)]

      margin_p[valid_response] <- tryCatch(

        do.call(margin_pFUN, args = call_input[FUN_args]),
        error = function(e) rep(NA_real_, sum(valid_response))
      )
    }

    return(margin_p)
  }

  margin_pFUN <- get(

    paste("p", margin_dist$family[1], sep = ""),
    envir = asNamespace("gamlss.dist"),
    mode = "function",
    inherits = FALSE
  )

  FUN <- margin_pFUN

  FUN_args <- names(margin_deriv_input)[names(margin_deriv_input) %in% formalArgs(FUN)]

  margin_p <- tryCatch(

    do.call(FUN, args = margin_deriv_input[FUN_args]),
    error = function(e) rep(NA_real_, length(response))
  )

  return(margin_p)
}
