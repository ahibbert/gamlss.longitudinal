#' @keywords internal
#' @noRd
score_function_v2 <- function(eta, dldpar, d2ldpar, dpardeta, response = NA, phi = 1, step_size = 1, verbose = FALSE, crit_wk = 0.0000001) {
  u_k <- dldeta <- dldpar * dpardeta
  f_k <- d2ldpar
  w_k <- -f_k * (dpardeta * dpardeta)

  # Stop if weights are too small
  w_k[abs(w_k) < crit_wk] <- 1
  u_k[abs(w_k) < crit_wk] <- 0

  w_k[abs(u_k) < crit_wk] <- 1
  u_k[abs(u_k) < crit_wk] <- 0

  z_k <- (1 - phi) * eta + phi * (eta + step_size * (u_k / w_k))

  if (verbose == TRUE) {
    steps_mean <- round(rbind(
      colMeans(as.matrix(eta)),
      colMeans(as.matrix(dldpar - dlcopdpar)),
      colMeans(as.matrix(dlcopdpar)),
      colMeans(as.matrix(dpardeta)),
      colMeans(as.matrix(dpardeta * dpardeta)),
      colMeans(as.matrix(f_k)),
      colMeans(as.matrix(w_k)),
      colMeans(as.matrix(u_k)),
      colMeans(as.matrix(u_k / w_k)),
      colMeans(as.matrix(z_k))
    ), 8)
    rownames(steps_mean) <- c("eta", "dldpar", "dlcopdpar", "dpardeta", "dpardeta2", "f_k", "w_k", "u_k", "(1/w_k)*u_k", "z_k")
    print(steps_mean)
  }
  return_list <- list(colMeans(as.matrix(z_k)), as.matrix(u_k), as.matrix(f_k), as.matrix(w_k), as.matrix(z_k))
  names(return_list) <- c("par", "u_k", "f_k", "w_k", "z_k")
  return(return_list)
}
