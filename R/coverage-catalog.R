#' @keywords internal

#' @noRd

.coverage_supported_copulas <- function() {

  c("N", "C", "F", "G", "J", "t")

}


#' @keywords internal

#' @noRd

.coverage_supported_methods <- function() {

  c("gamlss", "gamlss2", "rs_separate", "rs_joint", "cg", "gee", "glmm", "gam", "glmmTMB")

}


#' @keywords internal

#' @noRd

.coverage_supported_designs <- function() {

  c("intercept", "covariate", "scale", "time_dependence", "smooth")

}


#' @keywords internal

#' @noRd

.coverage_default_methods <- function() {

  c("gamlss", "rs_separate", "rs_joint", "cg")

}


#' @keywords internal

#' @noRd

.coverage_attach_namespace <- function(package) {

  if (!paste0("package:", package) %in% search()) {

    suppressPackageStartupMessages(attachNamespace(package))

  }

  invisible(TRUE)

}


#' @keywords internal

#' @noRd

.coverage_family_overrides <- function() {

  data.frame(

    family = c(

      "BI", "BB", "DBI", "ZABB", "ZABI", "ZIBB", "ZIBI",

      "LG",

      "MN3", "MN4", "MN5"

    ),

    supported = FALSE,

    unsupported_reason = c(

      rep("requires denominator/bounded-binomial response support not yet represented in the coverage likelihood calls", 7L),

      "logarithmic-series family needs family-specific starting/support handling before it can be all-method comparable",

      rep("ordinal/multinomial response support needs a family-specific simulation and likelihood path", 3L)

    ),

    stringsAsFactors = FALSE

  )

}


#' @keywords internal

#' @noRd

.coverage_family_catalog <- function(include_mixed = FALSE) {

  requireNamespace("gamlss.dist", quietly = TRUE)


  objects <- getNamespaceExports("gamlss.dist")

  candidates <- objects[grepl("^[A-Z][A-Z0-9]+$", objects)]

  rows <- lapply(candidates, function(family) {

    family_obj <- tryCatch(do.call(get(family, envir = asNamespace("gamlss.dist")), list()), error = function(e) NULL)

    if (is.null(family_obj) || is.null(family_obj$family) || is.null(family_obj$parameters)) {

      return(NULL)

    }

    family_name <- as.character(family_obj$family[1])

    q_name <- paste0("q", family_name)

    p_name <- paste0("p", family_name)

    d_name <- paste0("d", family_name)

    has_qpd <- all(vapply(c(q_name, p_name, d_name), exists, logical(1),

      envir = asNamespace("gamlss.dist"), inherits = FALSE

    ))

    type <- paste(as.character(family_obj$type), collapse = " ")

    mixed <- grepl("mixed", type, ignore.case = TRUE)

    supported <- has_qpd && (!mixed || isTRUE(include_mixed))

    unsupported_reason <- if (!has_qpd) {

      "missing q/p/d function"

    } else if (mixed && !isTRUE(include_mixed)) {

      "mixed support families are excluded from the default coverage grid"

    } else {

      NA_character_

    }

    data.frame(

      family = family_name,

      type = type,

      parameters = paste(names(family_obj$parameters), collapse = ","),

      supported = supported,

      unsupported_reason = unsupported_reason,

      stringsAsFactors = FALSE

    )

  })

  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

  out <- out[!duplicated(out$family), , drop = FALSE]

  overrides <- .coverage_family_overrides()

  override_idx <- match(overrides$family, out$family)

  override_idx <- override_idx[!is.na(override_idx)]

  if (length(override_idx) > 0L) {

    matched <- match(out$family[override_idx], overrides$family)

    out$supported[override_idx] <- overrides$supported[matched]

    out$unsupported_reason[override_idx] <- overrides$unsupported_reason[matched]

  }

  out[order(out$family), , drop = FALSE]

}


#' @keywords internal

#' @noRd

.coverage_safe_margin_param_overrides <- function(family_name) {

  overrides <- list(

    BI = list(mu = 0.5),

    LG = list(mu = 0.5),

    NO = list(mu = 0, sigma = 1),

    NO2 = list(mu = 0, sigma = 1),

    EXP = list(mu = 1),

    GA = list(mu = 2, sigma = 0.5),

    IG = list(mu = 2, sigma = 0.4),

    IGAMMA = list(mu = 2, sigma = 0.5),

    WEI = list(mu = 2, sigma = 0.8),

    WEI2 = list(mu = 2, sigma = 0.8),

    WEI3 = list(mu = 2, sigma = 0.8),

    PO = list(mu = 3),

    GEOM = list(mu = 2),

    YULE = list(mu = 3, sigma = 2),

    ZIPF = list(mu = 2, sigma = 2),

    NBI = list(mu = 4, sigma = 0.5),

    NBII = list(mu = 4, sigma = 0.5),

    DEL = list(mu = 4, sigma = 0.6, nu = 0.5),

    PIG = list(mu = 4, sigma = 0.5),

    ZIP = list(mu = 4, sigma = 0.2),

    ZIP2 = list(mu = 4, sigma = 0.2),

    ZAP = list(mu = 4, sigma = 0.2, nu = 0.5),

    ZINBI = list(mu = 4, sigma = 0.5, nu = 0.2),

    ZANBI = list(mu = 4, sigma = 0.5, nu = 0.2),

    ZAPIG = list(mu = 4, sigma = 0.5, nu = 0.2),

    ZASICHEL = list(mu = 4, sigma = 0.5, nu = 0.2),

    ZAZIPF = list(mu = 0.5, sigma = 0.2),

    BNB = list(mu = 4, sigma = 0.5, nu = 2),

    RGE = list(mu = 2, sigma = 0.3, nu = 1),

    SI = list(mu = 0.5, sigma = 0.05, nu = -0.5),

    ZABNB = list(mu = 4, sigma = 0.5, nu = 2, tau = 0.2),

    ZIBNB = list(mu = 4, sigma = 0.5, nu = 2, tau = 0.2),

    ZINBF = list(mu = 4, sigma = 0.5, nu = 2, tau = 0.2)

  )

  overrides[[family_name]] %||% list()

}


#' @keywords internal

#' @noRd

.coverage_make_case_grid <- function(

  families = NULL,

  copulas = .coverage_supported_copulas(),

  methods = .coverage_default_methods(),

  designs = "intercept",

  include_mixed = FALSE

) {

  catalog <- .coverage_family_catalog(include_mixed = include_mixed)

  if (is.null(families)) {

    families <- catalog$family[catalog$supported]

  }

  unsupported <- catalog[match(families, catalog$family), , drop = FALSE]

  unsupported <- unsupported[is.na(unsupported$family) | !unsupported$supported, , drop = FALSE]

  if (nrow(unsupported) > 0L) {

    bad <- ifelse(is.na(unsupported$family), families[is.na(match(families, catalog$family))], unsupported$family)

    stop("Unsupported coverage family/families: ", paste(bad, collapse = ", "), call. = FALSE)

  }

  bad_methods <- setdiff(methods, .coverage_supported_methods())

  if (length(bad_methods) > 0L) {

    stop("Unsupported coverage method(s): ", paste(bad_methods, collapse = ", "), call. = FALSE)

  }

  bad_designs <- setdiff(designs, .coverage_supported_designs())

  if (length(bad_designs) > 0L) {

    stop("Unsupported coverage design(s): ", paste(bad_designs, collapse = ", "), call. = FALSE)

  }

  grid <- expand.grid(

    family = families,

    copula = copulas,

    method = methods,

    design = designs,

    stringsAsFactors = FALSE,

    KEEP.OUT.ATTRS = FALSE

  )

  grid$case_id <- seq_len(nrow(grid))

  grid

}


#' @keywords internal

#' @noRd

.coverage_default_margin_params <- function(margin_dist) {

  family_name <- as.character(margin_dist$family[1])

  qfun <- get(paste0("q", family_name), envir = asNamespace("gamlss.dist"), inherits = FALSE)

  q_formals <- formals(qfun)

  par_names <- names(margin_dist$parameters)

  out <- list()

  for (par_name in par_names) {

    value <- q_formals[[par_name]]

    value <- tryCatch(eval(value, envir = baseenv()), error = function(e) NA_real_)

    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {

      value <- switch(par_name,

        mu = if (.is_discrete_margin(margin_dist)) 3 else 1,

        sigma = 0.7,

        nu = 0.3,

        tau = 2,

        0

      )

    }

    out[[par_name]] <- as.numeric(value)

  }


  safe_overrides <- .coverage_safe_margin_param_overrides(family_name)

  for (par_name in intersect(names(safe_overrides), names(out))) {

    out[[par_name]] <- safe_overrides[[par_name]]

  }

  if ("sigma" %in% names(out) && .is_discrete_margin(margin_dist)) out$sigma <- min(max(out$sigma, 0.3), 1.2)

  if ("nu" %in% names(out) && family_name %in% c("DEL", "ZIP", "ZAP", "ZINBI")) out$nu <- min(max(out$nu, 0.25), 0.75)

  out

}


#' @keywords internal

#' @noRd

.coverage_copula_params <- function(copula, dependence = c("moderate", "near_independent", "strong")) {

  dependence <- match.arg(dependence)

  tau <- switch(dependence,

    near_independent = 0.05,

    moderate = 0.25,

    strong = 0.55

  )

  rho <- switch(dependence,

    near_independent = 0.05,

    moderate = 0.3,

    strong = 0.65

  )


  if (copula %in% c("N")) {

    list(theta = rho)

  } else if (copula %in% c("t")) {

    list(theta = rho, zeta = 5)

  } else {

    list(tau = tau)

  }

}


#' @keywords internal

#' @noRd

.coverage_simulation_u_bounds <- function(family) {

  fragile <- c(

    "GT", "SEP", "SEP1", "SEP2", "SEP3", "SEP4",

    "SHASH", "SST", "ST1", "ST2", "ST3", "ST4", "ST5",

    "TF", "TF2"

  )

  if (family %in% fragile) c(1e-4, 1 - 1e-4) else NULL

}


#' @keywords internal

#' @noRd

.coverage_smooth_eta_component <- function(data, amplitude = 0.18) {

  x <- if ("x" %in% names(data)) {

    as.numeric(data$x)

  } else {

    sim_rescale01(as.numeric(data$.sim_subject_index)) - 0.5

  }

  x_scaled <- sim_rescale01(x)

  wave <- sin(2 * pi * x_scaled) + 0.5 * cos(4 * pi * x_scaled)

  wave <- wave - mean(wave, na.rm = TRUE)

  amplitude * wave

}


#' @keywords internal

#' @noRd

.coverage_make_smooth_param <- function(linkfun, linkinv, base_value, amplitude = 0.18) {

  force(linkfun)

  force(linkinv)

  force(base_value)

  force(amplitude)

  function(data) {

    base_eta <- as.numeric(linkfun(base_value))[1L]

    out <- linkinv(base_eta + .coverage_smooth_eta_component(data, amplitude = amplitude))

    as.numeric(out)

  }

}


#' @keywords internal

#' @noRd

.coverage_time_varying_copula_params <- function(copula) {

  if (copula %in% c("N")) {

    return(list(theta = function(edge_data) {

      time_scaled <- sim_rescale01(edge_data$time_left)

      0.1 + 0.45 * time_scaled

    }))

  }

  if (copula %in% c("t")) {

    return(list(

      theta = function(edge_data) {

        time_scaled <- sim_rescale01(edge_data$time_left)

        0.1 + 0.45 * time_scaled

      },

      zeta = 5

    ))

  }

  list(tau = function(edge_data) {

    time_scaled <- sim_rescale01(edge_data$time_left)

    0.1 + 0.3 * time_scaled

  })

}


#' @keywords internal

#' @noRd
