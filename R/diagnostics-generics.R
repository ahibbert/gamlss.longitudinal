# Silence NSE checks for ggplot mappings.

utils::globalVariables(c("theoretical", "observed", "detrended", "midpoint", "root_diff", "band_lower", "band_upper", "time", "split_group", "expected"))

#' Diagnostic generics for fitted longitudinal GAMLSS-copula models

#'

#' These generics dispatch to methods for `gamlss.longitudinal` objects and

#' return either diagnostic plot objects, scoring summaries, or forecast data.

#'

#' @name diagnostics

#' @aliases pithist qqrplot wormplot rootogram proscore procast

NULL

#' @export

pithist <- function(object, ...) {
  UseMethod("pithist")
}

#' @export

qqrplot <- function(object, ...) {
  UseMethod("qqrplot")
}

#' @export

wormplot <- function(object, ...) {
  UseMethod("wormplot")
}

#' @export

rootogram <- function(object, ...) {
  UseMethod("rootogram")
}

#' @export

proscore <- function(object, ...) {
  UseMethod("proscore")
}

#' @export

procast <- function(object, ...) {
  UseMethod("procast")
}
