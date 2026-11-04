#' Prepare bootstrap RNG seed restoration
#'
#' @param seed Optional seed supplied to `bootstrap_inference()`.
#' @return A zero-argument restore function, or `NULL` when `seed` is `NULL`.
#' @noRd
.gl_prepare_bootstrap_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
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
  }
}
