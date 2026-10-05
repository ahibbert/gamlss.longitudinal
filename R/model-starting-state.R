#' Generate starting values for model covariates
#'
#' Either uses start_from or generates starting values from the data and family and expands this 
#' into `par_cov`, `par_s`, `df_s`, and `lambda_s` which are used by the optimiser loop.
#'
#' @noRd
.gl_build_initial_parameter_state <- function(
    start_from,
    warm_start_par_s = NULL,
    warm_start_lambda_s = NULL,
    mm,
    margin_dist,
    copula_dist,
    dataset,
    lambda_start = NA) {
  if (is.list(start_from) && !is.null(start_from$par)) {
    start_from <- start_from$par
  }

  if (all(is.na(start_from))) {
    par_eta <- get_starting_values(copula_dist, margin_dist, dataset = dataset, eta_transform = TRUE)
    par_cov <- as.numeric(vector())

    for (par_name in names(mm$x)) {
      par_cov_single <- as.numeric(vector(length = length(colnames(mm$x[[par_name]]))))
      names(par_cov_single) <- paste(par_name, colnames(mm$x[[par_name]]), sep = ".")
      par_cov_single[1] <- par_eta[par_name]
      if (length(par_cov_single) > 1) {
        par_cov_single[2:length(par_cov_single)] <- 0
      }
      par_cov <- c(par_cov, par_cov_single)
    }
  } else {
    expected <- unlist(Map(function(par_name, design) {
      paste(par_name, colnames(design), sep = ".")
    }, names(mm$x), mm$x), use.names = FALSE)
    mapped_fitted_state <- !is.null(warm_start_par_s) || !is.null(warm_start_lambda_s)
    if (!is.null(names(start_from)) && all(expected %in% names(start_from))) {
      par_cov <- stats::setNames(as.numeric(start_from[expected]), expected)
    } else if (!is.null(names(start_from)) && mapped_fitted_state &&
               any(expected %in% names(start_from))) {
      # Formula screening can introduce a fixed coefficient that was not in
      # the fitted parent (for example, replacing s(age) by age or expanding a
      # previously intercept-only parameter). Retain every compatible fitted
      # coefficient and start genuinely new effects at the neutral value zero.
      # Raw user vectors remain strict so misspelled coefficient names are not
      # silently accepted.
      par_cov <- stats::setNames(numeric(length(expected)), expected)
      retained <- intersect(expected, names(start_from))
      par_cov[retained] <- as.numeric(start_from[retained])
    } else if (length(start_from) == length(expected)) {
      par_cov <- stats::setNames(as.numeric(start_from), expected)
    } else {
      stop(
        "The supplied start does not contain every fixed coefficient required by the model.",
        call. = FALSE
      )
    }
  }

  par_s <- list()
  df_s <- list()
  lambda_s <- list()

  for (par_name in names(mm$x)) {
    par_s[[par_name]] <- list()
    df_s[[par_name]] <- list()
    lambda_s[[par_name]] <- list()

    for (s_name in names(mm$s[[par_name]])) {
      B <- mm$s[[par_name]][[s_name]]
      par_s[[par_name]][[s_name]] <- c(par_s[[par_name]][[s_name]], rep(0, ncol(B)))
      names(par_s[[par_name]][[s_name]]) <- paste(par_name, s_name, 1:ncol(B), sep = ".")
      df_s[[par_name]][[s_name]] <- 0

      S_init <- attr(B, "penalty")
      if (is.na(lambda_start)) {
        if (!is.null(S_init) && is.matrix(S_init) && sum(diag(S_init)) > 0) {
          lambda_s[[par_name]][[s_name]] <- sum(diag(t(B) %*% B)) / sum(diag(S_init))
        } else {
          lambda_s[[par_name]][[s_name]] <- 10
        }
      } else {
        lambda_s[[par_name]][[s_name]] <- lambda_start
      }

      names(df_s[[par_name]][[s_name]]) <- names(lambda_s[[par_name]][[s_name]]) <- s_name
    }
  }

  if (!is.null(warm_start_par_s)) {
    for (par_name in intersect(names(par_s), names(warm_start_par_s))) {
      if (length(par_s[[par_name]]) == 0 || length(warm_start_par_s[[par_name]]) == 0) next
      for (s_name in intersect(names(par_s[[par_name]]), names(warm_start_par_s[[par_name]]))) {
        warm_beta <- warm_start_par_s[[par_name]][[s_name]]
        if (length(warm_beta) == length(par_s[[par_name]][[s_name]])) {
          par_s[[par_name]][[s_name]] <- warm_beta
        }
      }
    }
  }

  if (!is.null(warm_start_lambda_s)) {
    for (par_name in intersect(names(lambda_s), names(warm_start_lambda_s))) {
      if (!length(lambda_s[[par_name]]) || !length(warm_start_lambda_s[[par_name]])) next
      for (s_name in intersect(names(lambda_s[[par_name]]), names(warm_start_lambda_s[[par_name]]))) {
        warm_lambda <- warm_start_lambda_s[[par_name]][[s_name]]
        if (length(warm_lambda) == 1L && is.finite(warm_lambda) && warm_lambda > 0) {
          lambda_s[[par_name]][[s_name]] <- warm_lambda
        }
      }
    }
  }

  list(
    par_cov = par_cov,
    par_s = par_s,
    df_s = df_s,
    lambda_s = lambda_s
  )
}
