.build_margin_eval_cache <- function(margin_dist, calc_d2 = FALSE) {
  if (calc_d2 == TRUE) {
    to_include <- grepl("dld", names(margin_dist)) | grepl("d2ld", names(margin_dist))
  } else {
    to_include <- grepl("dld", names(margin_dist))
  }

  margin_deriv_names <- names(margin_dist)[to_include]
  margin_deriv_cache <- lapply(margin_deriv_names, function(deriv_name) {
    FUN <- margin_dist[[deriv_name]]
    list(name = deriv_name, FUN = FUN, args = formalArgs(FUN))
  })

  margin_pFUN <- get(
    paste("p", margin_dist$family[1], sep = ""),
    envir = asNamespace("gamlss.dist"),
    mode = "function",
    inherits = FALSE
  )
  margin_dFUN <- get(
    paste("d", margin_dist$family[1], sep = ""),
    envir = asNamespace("gamlss.dist"),
    mode = "function",
    inherits = FALSE
  )

  list(
    calc_d2 = calc_d2,
    margin_deriv_cache = margin_deriv_cache,
    margin_pFUN = margin_pFUN,
    margin_p_args = formalArgs(margin_pFUN),
    margin_dFUN = margin_dFUN,
    margin_d_args = formalArgs(margin_dFUN),
    family_call_cache = new.env(parent = emptyenv())
  )
}

.call_margin_family_cached <- function(FUN, args, FUN_args, cacheable = FALSE, cache_env = NULL, cache_prefix = "") {
  call_args <- args[FUN_args]
  safe_call <- function(call_args, fallback_n = NULL) {
    if (is.null(fallback_n)) {
      fallback_n <- max(c(1L, vapply(call_args, length, integer(1))), na.rm = TRUE)
    }
    tryCatch(
      do.call(FUN, args = call_args),
      error = function(e) rep(NA_real_, fallback_n)
    )
  }
  if (!isTRUE(cacheable) || length(call_args) == 0L) {
    return(safe_call(call_args))
  }

  arg_lengths <- vapply(call_args, length, integer(1))
  n <- max(arg_lengths)
  if (n <= 1L) {
    return(safe_call(call_args, fallback_n = n))
  }

  if (any(!arg_lengths %in% c(1L, n))) {
    return(safe_call(call_args, fallback_n = n))
  }
  if (!all(vapply(call_args, function(x) is.atomic(x) && !is.character(x), logical(1)))) {
    return(safe_call(call_args, fallback_n = n))
  }

  expanded <- lapply(call_args, rep, length.out = n)
  key_parts <- lapply(expanded, function(x) {
    if (is.numeric(x) || is.integer(x)) {
      format(signif(as.numeric(x), 15), scientific = FALSE, trim = TRUE)
    } else {
      as.character(x)
    }
  })
  key <- do.call(paste, c(key_parts, sep = "\r"))
  unique_idx <- !duplicated(key)
  if (sum(unique_idx) > 0.8 * n) {
    return(safe_call(call_args, fallback_n = n))
  }

  if (!is.null(cache_env)) {
    cached <- rep(NA_real_, n)
    hit <- logical(n)
    lookup_key <- paste0(cache_prefix, "\r", key)
    for (ii in seq_len(n)) {
      if (exists(lookup_key[[ii]], envir = cache_env, inherits = FALSE)) {
        cached[[ii]] <- get(lookup_key[[ii]], envir = cache_env, inherits = FALSE)
        hit[[ii]] <- TRUE
      }
    }
    if (all(hit)) {
      return(cached)
    }
    need <- !hit
    need_unique <- !duplicated(key[need])
    need_args <- lapply(expanded, function(x) x[need][need_unique])
    names(need_args) <- names(call_args)
    need_val <- safe_call(need_args, fallback_n = sum(need_unique))
    if (length(need_val) != sum(need_unique)) {
      return(safe_call(call_args, fallback_n = n))
    }
    need_lookup <- paste0(cache_prefix, "\r", key[need][need_unique])
    for (jj in seq_along(need_lookup)) {
      assign(need_lookup[[jj]], need_val[[jj]], envir = cache_env)
    }
    cached[need] <- need_val[match(key[need], key[need][need_unique])]
    return(cached)
  }

  unique_args <- lapply(expanded, `[`, unique_idx)
  names(unique_args) <- names(call_args)
  unique_val <- safe_call(unique_args, fallback_n = sum(unique_idx))
  if (length(unique_val) != sum(unique_idx)) {
    return(safe_call(call_args, fallback_n = n))
  }
  unique_val[match(key, key[unique_idx])]
}

.call_fast_count_family <- function(prefix, family_name, args) {
  if (!identical(family_name, "PO") || !"mu" %in% names(args)) {
    return(NULL)
  }
  x <- args$x %||% args$q %||% args$y
  mu <- args$mu
  if (is.null(x) || is.null(mu)) {
    return(NULL)
  }
  if (identical(prefix, "d")) {
    return(stats::dpois(x, lambda = mu))
  }
  if (identical(prefix, "p")) {
    return(stats::ppois(x, lambda = mu))
  }
  NULL
}
