calc_true_SE_numderiv_only_covariates <- function(object, par, mm, margin_dist, response, testing = FALSE, response_margin = NA, response_subject = NA, h = .0001, progress = interactive()) {
  # object=fit; par=NA; numderiv=TRUE; sep_d2=TRUE

  response <- object$response
  response_margin <- object$response_margin
  response_subject <- object$response_subject

  # margin_names=unique(object$response_margin)
  # num_margins=length(margin_names)

  # se_out=object$par*0;
  margin_dist <- object$margin_dist
  copula_dist <- object$copula_dist
  copula_link <- get_copula_dist(copula_dist)$copula_link
  mm <- object$model_matrix

  par_cov <- object$par
  par_s <- object$par_s

  input_par <- par_cov
  progress <- isTRUE(progress)

  adj_fac <- h
  par_names <- names(input_par)
  nd_impact <- rep(0, length(par_names))
  names(nd_impact) <- par_names # [names(eta_inv) %in% c("mu","sigma","nu","tau")]

  # Reuse fixed copula pairing metadata across all numerical derivative evaluations.
  pair_cache <- build_copula_pair_cache(response, response_margin, response_subject)
  eval_joint_loglik <- function(par_vec) {
    eta_out <- calc_eta(par_cov = par_vec, mm = mm, margin_dist, copula_link, par_s)
    eta_inv <- eta_out$eta_inv
    calc_likelihood_minimal(
      eta_inv,
      mm = mm$x,
      margin_dist,
      copula_dist,
      calc_d2 = FALSE,
      response = response,
      response_margin = response_margin,
      response_subject = response_subject,
      pair_cache = pair_cache,
      calc_margin_deriv = FALSE
    )$log_lik["joint"]
  }

  base_loglik <- eval_joint_loglik(input_par)

  if (progress) cat("Calculating numerical first derivates for Hessian matrix...\n")
  pb_first <- NULL
  if (progress) {
    pb_first <- utils::txtProgressBar(min = 0, max = length(par_names), style = 3)
    on.exit(close(pb_first), add = TRUE)
  }
  first_counter <- 0L
  for (par_names_nd in par_names) {
    # print(par_names_nd)
    change <- rep(0, 3)
    i <- 1
    for (adj in c(-1 * adj_fac, adj_fac)) {
      par_cov <- input_par
      par_cov[[par_names_nd]] <- par_cov[[par_names_nd]] + adj

      change[i] <- eval_joint_loglik(par_cov)
      i <- i + 1
    }
    change[3] <- base_loglik
    nd_impact[par_names_nd] <- (change[2] + change[1] - 2 * change[3]) / (adj_fac^2)

    first_counter <- first_counter + 1L
    if (progress) {
      utils::setTxtProgressBar(pb_first, first_counter)
    }

    # print(c(change,nd_impact[eta_par_names_nd]))
  }
  if (progress) cat("\n")

  if (progress) cat("Calculating numerical second derivates for Hessian matrix... this may take a while\n")
  p <- length(par_names)
  nd_cross <- matrix(0, nrow = p, ncol = p)
  colnames(nd_cross) <- rownames(nd_cross) <- par_names
  second_total <- if (p > 1) choose(p, 2) else 0
  pb_second <- NULL
  if (progress && second_total > 0) {
    pb_second <- utils::txtProgressBar(min = 0, max = second_total, style = 3)
    on.exit(close(pb_second), add = TRUE)
  }
  second_counter <- 0L
  if (p > 1) {
    for (i in 1:(p - 1)) {
      name1 <- par_names[i]
      for (j in (i + 1):p) {
        name2 <- par_names[j]
        cross_sum <- 0
        for (adj1 in c(-1 * adj_fac, adj_fac)) {
          for (adj2 in c(-1 * adj_fac, adj_fac)) {
            par <- input_par
            par[[name1]] <- par[[name1]] + adj1
            par[[name2]] <- par[[name2]] + adj2

            change <- eval_joint_loglik(par)
            cross_sum <- cross_sum + change * if (adj1 == adj2) {
              1
            } else {
              -1
            }
          }
        }
        nd_cross[name1, name2] <- cross_sum
        nd_cross[name2, name1] <- cross_sum

        second_counter <- second_counter + 1L
        if (progress && !is.null(pb_second)) {
          utils::setTxtProgressBar(pb_second, second_counter)
        }
      }
    }
  }
  if (progress && !is.null(pb_second)) cat("\n")
  nd_cross <- nd_cross / (4 * (adj_fac^2))

  nd2 <- (diag(nd_impact) + nd_cross)

  return(nd2)
}
