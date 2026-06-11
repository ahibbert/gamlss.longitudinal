.coverage_simulate_case <- function(

  family,

  copula,

  design = .coverage_supported_designs(),

  n = 80,

  times = 1:3,

  seed = 1,

  dependence = "moderate"

) {

  design <- match.arg(design)

  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  margin_params <- .coverage_default_margin_params(margin_dist)

  covariates <- NULL


  if (design %in% c("covariate", "scale", "smooth")) {

    covariates <- function(base) {

      x <- sim_rescale01(as.numeric(base$.sim_subject_index))

      data.frame(x = x - mean(x), stringsAsFactors = FALSE)

    }

    if (identical(design, "covariate") && "mu" %in% names(margin_params)) {

      base_mu <- margin_params$mu

      margin_params$mu <- function(data) {

        x <- sim_rescale01(as.numeric(data$.sim_subject_index))

        x <- x - mean(x)

        if (base_mu > 0) {

          pmax(base_mu * exp(0.15 * x), .Machine$double.eps)

        } else {

          base_mu + 0.15 * x

        }

      }

    }

    if (identical(design, "smooth")) {

      for (par_name in names(margin_params)) {

        linkfun <- margin_dist[[paste0(par_name, ".linkfun")]]

        linkinv <- margin_dist[[paste0(par_name, ".linkinv")]]

        if (!is.function(linkfun) || !is.function(linkinv)) next

        margin_params[[par_name]] <- .coverage_make_smooth_param(

          linkfun,

          linkinv,

          margin_params[[par_name]],

          amplitude = if (identical(par_name, "mu")) 0.22 else 0.14

        )

      }

    }

    if (identical(design, "scale") && "sigma" %in% names(margin_params)) {

      base_sigma <- margin_params$sigma

      margin_params$sigma <- function(data) {

        x <- sim_rescale01(as.numeric(data$.sim_subject_index))

        x <- x - mean(x)

        pmax(base_sigma * exp(0.45 * x), .Machine$double.eps)

      }

    }

  }


  copula_params <- if (identical(design, "time_dependence")) {

    .coverage_time_varying_copula_params(copula)

  } else if (identical(design, "smooth")) {

    copula_link <- get_copula_dist(copula)$copula_link

    base <- .coverage_copula_params(copula, dependence = dependence)

    if ("tau" %in% names(base)) {

      list(tau = function(edge_data) {

        eta <- stats::qlogis(base$tau) + .coverage_smooth_eta_component(edge_data, amplitude = 0.12)

        stats::plogis(eta)

      })

    } else {

      base$theta <- .coverage_make_smooth_param(copula_link$theta.linkfun, copula_link$theta.linkinv, base$theta, amplitude = 0.12)

      base

    }

  } else {

    .coverage_copula_params(copula, dependence = dependence)

  }


  simulate_longitudinal_dataset(

    n = n,

    times = times,

    margin_dist = margin_dist,

    copula_dist = copula,

    margin_params = margin_params,

    copula_params = copula_params,

    covariates = covariates,

    seed = seed,

    include_truth = TRUE,

    u_bounds = .coverage_simulation_u_bounds(family)

  )

}


#' @keywords internal

#' @noRd

.coverage_apply_missingness <- function(dat, missingness = c("none", "mcar", "drop_rows"), prop = 0.05) {

  missingness <- match.arg(missingness)

  if (identical(missingness, "none")) {

    return(dat)

  }


  n_subject <- length(unique(dat$subject))

  n_time <- length(unique(dat$time))

  n_target <- max(1L, floor(nrow(dat) * prop))


  if (identical(missingness, "mcar")) {

    eligible <- which(dat$time != min(dat$time))

    idx <- eligible[seq_len(min(length(eligible), n_target))]

    dat$response[idx] <- NA_real_

    return(dat)

  }


  subject_index <- as.integer(dat$subject)

  time_values <- sort(unique(dat$time))

  drop_idx <- subject_index <= max(1L, floor(n_subject * prop)) & dat$time == time_values[min(2L, n_time)]

  dat[!drop_idx, , drop = FALSE]

}


#' @keywords internal

#' @noRd

.coverage_fit_formulas <- function(design) {

  if (identical(design, "covariate")) {

    list(mu = response ~ x, sigma = ~1, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)

  } else if (identical(design, "scale")) {

    list(mu = response ~ 1, sigma = ~x, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)

  } else if (identical(design, "time_dependence")) {

    list(mu = response ~ 1, sigma = ~1, nu = ~1, tau = ~1, theta = ~time, zeta = ~1)

  } else if (identical(design, "smooth")) {

    list(

      mu = response ~ s(x, bs = "ps", k = 6),

      sigma = ~s(x, bs = "ps", k = 6),

      nu = ~s(x, bs = "ps", k = 6),

      tau = ~s(x, bs = "ps", k = 6),

      theta = ~s(x, bs = "ps", k = 6),

      zeta = ~s(x, bs = "ps", k = 6)

    )

  } else {

    list(mu = response ~ 1, sigma = ~1, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)

  }

}


#' @keywords internal

#' @noRd
