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

  h * pmax(abs(par), 1)
}
