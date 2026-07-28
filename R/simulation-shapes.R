#' Rescale a numeric vector to the unit interval
#'
#' @param x Numeric vector.
#' @return Numeric vector scaled to `[0, 1]`. Constant finite inputs return 0.
#' @export
sim_rescale01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng))) {
    return(rep(NA_real_, length(x)))
  }
  width <- diff(rng)
  if (width <= 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / width
}

#' Simple deterministic smooth shapes for simulation truth
#'
#' These helpers are intended for readable simulation specifications. They work
#' best with inputs already scaled to `[0, 1]`, for example by
#' [sim_rescale01()].
#'
#' @param x Numeric input, usually scaled to `[0, 1]`.
#' @param slope Linear slope.
#' @param center Logical; if `TRUE`, center `x` around 0.5 for the linear shape.
#' @param amplitude Effect amplitude.
#' @param period Period of the sinusoid in units of `x`.
#' @param phase Phase shift for the sinusoid.
#' @param location Center/location of the bump, sigmoid, or U-shape.
#' @param width Width of the Gaussian bump.
#' @param steepness Steepness of the sigmoid transition.
#'
#' @return Numeric vector of the same length as `x`.
#' @name sim_smooth_shapes
NULL

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_linear <- function(x, slope = 1, center = TRUE) {
  x <- as.numeric(x)
  if (isTRUE(center)) x <- x - 0.5
  slope * x
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_sin <- function(x, amplitude = 1, period = 1, phase = 0) {
  x <- as.numeric(x)
  amplitude * sin(2 * pi * (x / period + phase))
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_bump <- function(x, amplitude = 1, location = 0.5, width = 0.15) {
  x <- as.numeric(x)
  amplitude * exp(-0.5 * ((x - location) / width)^2)
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_sigmoid <- function(x, amplitude = 1, location = 0.5, steepness = 10) {
  x <- as.numeric(x)
  amplitude * (stats::plogis(steepness * (x - location)) - 0.5)
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_u <- function(x, amplitude = 1, location = 0.5) {
  x <- as.numeric(x)
  amplitude * ((x - location)^2 - mean((x - location)^2, na.rm = TRUE))
}

#' @rdname sim_smooth_shapes
#' @export
sim_smooth_wiggle <- function(x, amplitude = 1) {
  x <- as.numeric(x)
  amplitude * (0.65 * sin(2 * pi * x) + 0.35 * sin(6 * pi * x + 0.4))
}

#' Look up named effects for a factor-like variable
#'
#' @param x Factor, character, or numeric vector.
#' @param effects Named numeric vector of effects.
#' @param reference Optional reference level to prepend with effect 0 when it is
#'   not already present in `effects`.
#'
#' @return Numeric vector of effects aligned with `x`.
#' @export
sim_factor_effect <- function(x, effects, reference = NULL) {
  if (is.null(names(effects)) || any(!nzchar(names(effects)))) {
    stop("effects must be a named numeric vector.", call. = FALSE)
  }
  effects <- as.numeric(effects) |>
    stats::setNames(names(effects))
  if (!is.null(reference) && !reference %in% names(effects)) {
    effects <- c(stats::setNames(0, reference), effects)
  }
  x_chr <- as.character(x)
  out <- effects[x_chr]
  if (anyNA(out)) {
    missing_levels <- unique(x_chr[is.na(out)])
    stop(
      "Missing factor effect for level(s): ",
      paste(missing_levels, collapse = ", "),
      ". Available effect names: ",
      paste(names(effects), collapse = ", "),
      ". Supply one named effect per observed level, or use reference to set ",
      "a zero effect for a missing reference level.",
      call. = FALSE
    )
  }
  as.numeric(out)
}
