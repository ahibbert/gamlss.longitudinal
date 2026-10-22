#' Choose a finite-difference step on the natural parameter scale
#'
#' Identity-linked parameters use the base step; other links scale by parameter
#' magnitude to reduce avoidable cancellation.
#'
#' @noRd
.natural_fd_step <- function(par, par_name, margin_dist, h) {
  link_name <- margin_dist[[paste(par_name, "link", sep = ".")]]

  if (is.character(link_name) && identical(link_name[1], "identity")) {
    return(rep(h, length(par)))
  }

  step <- h * pmax(abs(par), 1)

  # Keep symmetric natural-scale perturbations strictly inside the parameter
  # domain. This matters for fitted zero-inflation probabilities that can be
  # much closer to zero or one than the nominal Hessian step. An invalid
  # observed-row evaluation after this contraction remains a computational
  # error; it is not replaced by a zero derivative.
  if (is.character(link_name) && identical(link_name[1], "log")) {
    step <- pmin(step, 0.25 * par)
  } else if (is.character(link_name) && identical(link_name[1], "logit")) {
    step <- pmin(step, 0.25 * pmin(par, 1 - par))
  }

  if (any(!is.finite(step) | step <= 0)) {
    stop(structure(
      list(
        message = paste0(
          "Could not construct an interior finite-difference step for margin parameter '",
          par_name, "'."
        ),
        call = NULL,
        rows = which(!is.finite(step) | step <= 0),
        context = "Natural-scale Hessian finite-difference step"
      ),
      class = c("gamlss_longitudinal_score_evaluation_error", "error", "condition")
    ))
  }

  step
}
