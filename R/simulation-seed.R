#' Temporarily set and restore simulation random seed
#'
#' @noRd
.sim_set_seed_restore <- function(seed) {
  if (is.null(seed)) {
    return(function() invisible(NULL))
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  set.seed(seed)

  function() {
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}
