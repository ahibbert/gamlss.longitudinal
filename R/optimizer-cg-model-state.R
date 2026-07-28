#' Build augmented CG model matrices and coefficient vector
#'
#' Moves smooth basis columns into the fixed matrix representation used by the
#' CG optimizer and appends smooth coefficients to the optimization vector.
#'
#' @noRd
.gl_build_cg_model <- function(mm, par_cov, par_s) {
  mm_cg <- mm
  beta <- par_cov
  for (pn in names(mm$x)) {
    if (length(mm$s[[pn]]) > 0) {
      for (sn in names(mm$s[[pn]])) {
        B <- mm$s[[pn]][[sn]]
        b <- par_s[[pn]][[sn]]
        colnames(B) <- sub(paste0("^", pn, "\\."), "", names(b))
        mm_cg$x[[pn]] <- cbind(mm_cg$x[[pn]], B)
        beta <- c(beta, b)
      }
    }
    mm_cg$s[[pn]] <- list()
  }
  list(mm = mm_cg, beta = beta)
}

#' Unpack CG coefficient vector into fixed and smooth parameter state
#'
#' @noRd
.gl_unpack_cg_beta <- function(beta_vec, par_cov_template, par_s_template) {
  par_cov_new <- beta_vec[names(par_cov_template)]
  par_s_new <- par_s_template
  for (pn in names(par_s_new)) {
    if (length(par_s_new[[pn]]) == 0) next
    for (sn in names(par_s_new[[pn]])) {
      b_names <- names(par_s_new[[pn]][[sn]])
      par_s_new[[pn]][[sn]] <- beta_vec[b_names]
    }
  }
  list(par_cov = par_cov_new, par_s = par_s_new)
}
