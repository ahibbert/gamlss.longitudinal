.sim_margin_quantile_function <- function(margin_dist) {
  family_name <- margin_dist$family[1]
  if (!is.character(family_name) || length(family_name) != 1L || is.na(family_name)) {
    stop("margin_dist must be a gamlss.dist family object.", call. = FALSE)
  }
  qfun_name <- paste0("q", family_name)
  if (!exists(qfun_name, envir = asNamespace("gamlss.dist"), inherits = FALSE)) {
    stop("Could not find gamlss.dist quantile function ", qfun_name, "().", call. = FALSE)
  }
  get(qfun_name, envir = asNamespace("gamlss.dist"), inherits = FALSE)
}

.sim_eval_long_param <- function(spec, data, n, n_time, label) {
  if (is.function(spec)) {
    spec <- spec(data)
  }
  .sim_expand_param(spec, n = n, n_time = n_time, index = data$.sim_time_index, label = label)
}

.sim_eval_edge_param <- function(spec, edge_data, n, n_edge, label) {
  if (is.function(spec)) {
    spec <- spec(edge_data)
  }
  .sim_expand_param(spec, n = n, n_time = n_edge, index = edge_data$.sim_edge_index, label = label)
}

.sim_expand_param <- function(spec, n, n_time, index, label) {
  if (is.matrix(spec)) {
    if (!identical(dim(spec), c(n, n_time))) {
      stop(label, " matrix must have dimensions n by number of time positions.", call. = FALSE)
    }
    return(as.vector(t(spec)))
  }
  spec <- as.numeric(spec)
  if (length(spec) == 1L) {
    return(rep(spec, length(index)))
  }
  if (length(spec) == n_time) {
    return(spec[index])
  }
  if (length(spec) == length(index)) {
    return(spec)
  }
  stop(
    label,
    " must be a scalar, a time-position vector, a full-length vector, ",
    "a subject-by-time matrix, or a function returning one of those.",
    call. = FALSE
  )
}
