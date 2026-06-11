#' Create model matrices for model fitting

#'

#' This function takes the forumlas for each parameter mu,sigma,nu,tau,theta,zeta

#' and creates a list of model matrices mm with items mm$x and mm$s for

#' fixed and smooth terms respectively, with each of those lists being lists of each parameter

#' and their respective model matrices

#' @param mu.formula Formula for the mean parameter of the marginal distribution

#' @param sigma.formula Formula for the sigma parameter of the marginal distribution

#' @param nu.formula Formula for the nu parameter of the marginal distribution

#' @param tau.formula Formula for the tau parameter of the marginal distribution

#' @param theta.formula Formula for the theta parameter of the copula distribution

#' @param zeta.formula Formula for the zeta parameter of the copula distribution

#' @param margin.family Marginal distribution specified as a gamlss family object,

#' e.g. GA(), NO(), PO(), NBI(), etc.

#' @param copula.family Copula distribution code, one of "N", "C", "F", "G", "J", or "t".

#' @param copula.link List of link functions for the copula parameters

#' @return Returns a list mm with items mm$x and mm$s for fixed and smooth terms respectively,

#' with each of those lists being lists of each parameter and their respective model matrices

#'

#' @keywords internal

#' @noRd

create_model_matrices<-function(

    mu.formula = ("response ~ 1"),

    sigma.formula = ("1"),

    nu.formula = ("1"),

    tau.formula = ("1"),

    theta.formula=("1"),

    zeta.formula=("1"),

    margin.family=NO(),

    copula.family="N",

    copula.link=NA,

    dataset=NA,

    quiet_gamlss2=TRUE,

    preserve_factor_levels=FALSE

) {


  dataset_mm <- dataset

  normalize_ordered_factor <- function(col) {

    if (!is.factor(col)) return(col)

    if (!is.ordered(col)) return(col)


    levs <- levels(col)

    col_nom <- factor(as.character(col), levels = levs, ordered = FALSE)

    if (length(levs) > 1) {

      contr <- contr.treatment(length(levs))

      colnames(contr) <- levs[-1]

      contrasts(col_nom) <- contr

    }

    col_nom

  }


  mode_value <- function(x) {

    x_non_na <- x[!is.na(x)]

    if (length(x_non_na) == 0) return(NA)

    tab <- table(x_non_na)

    names(tab)[which.max(tab)]

  }


  # Build model matrices from an NA-free proxy dataset.

  # This does not alter likelihood calculations (which still use original dataset).

  for (nm in names(dataset_mm)) {

    if (!any(is.na(dataset_mm[[nm]]))) next

    if (nm %in% c("time", "subject")) next


    col <- dataset_mm[[nm]]

    if (is.numeric(col) || is.integer(col)) {

      obs <- col[!is.na(col)]

      fill_val <- if (length(obs) > 0) mean(obs) else 0

      col[is.na(col)] <- fill_val

      dataset_mm[[nm]] <- col

    } else if (is.factor(col)) {

      col <- normalize_ordered_factor(col)

      fill_val <- mode_value(col)

      if (is.na(fill_val)) {

        fill_val <- if (length(levels(col)) > 0) levels(col)[1] else "missing"

      }

      col_chr <- as.character(col)

      col_chr[is.na(col_chr)] <- fill_val

      dataset_mm[[nm]] <- factor(col_chr, levels = levels(col), ordered = FALSE)

    } else {

      fill_val <- mode_value(col)

      if (is.na(fill_val)) fill_val <- "missing"

      col[is.na(col)] <- fill_val

      dataset_mm[[nm]] <- col

    }

  }


  if ("response" %in% names(dataset_mm) && any(is.na(dataset_mm$response))) {

    obs_resp <- dataset_mm$response[!is.na(dataset_mm$response)]

    fill_val <- if (length(obs_resp) > 0) mean(obs_resp) else 0

    dataset_mm$response[is.na(dataset_mm$response)] <- fill_val

  }


  normalize_time_covariate_colnames <- function(nms) {

    if (length(nms) == 0) return(nms)


    out <- nms

    suffix_map <- c(L = "1", Q = "2", C = "3")


    for (sx in names(suffix_map)) {

      out <- gsub(

        paste0("time_covariate\\.", sx, "\\b"),

        paste0("time_covariate.", suffix_map[[sx]]),

        out,

        perl = TRUE

      )

    }


    # contr.poly names can appear as ^4, ^5, ...; normalize to .4, .5, ...

    out <- gsub("time_covariate\\^([0-9]+)", "time_covariate.\\1", out, perl = TRUE)


    out

  }


  sanitize_for_gamlss2 <- function(data_in, fml) {

    vars_needed <- unique(all.vars(stats::as.formula(fml)))

    vars_needed <- vars_needed[vars_needed %in% names(data_in)]

    data_out <- data_in[, vars_needed, drop = FALSE]


    for (nm in names(data_out)) {

      col <- data_out[[nm]]

      if (is.factor(col)) {

        col <- normalize_ordered_factor(col)

        if (!isTRUE(preserve_factor_levels)) {

          col <- droplevels(col)

        }

        data_out[[nm]] <- col

      } else if (is.numeric(col) || is.integer(col)) {

        col[!is.finite(col)] <- NA

        if (any(is.na(col))) {

          obs <- col[!is.na(col)]

          fill_val <- if (length(obs) > 0) mean(obs) else 0

          col[is.na(col)] <- fill_val

        }

        data_out[[nm]] <- col

      } else {

        if (any(is.na(col))) {

          x_non_na <- col[!is.na(col)]

          fill_val <- if (length(x_non_na) > 0) {

            tab <- table(x_non_na)

            names(tab)[which.max(tab)]

          } else {

            "missing"

          }

          col[is.na(col)] <- fill_val

        }

        data_out[[nm]] <- col

      }

    }


    data_out

  }


  to_response_formula <- function(fml, response_name = "response") {

    if (inherits(fml, "formula")) {

      rhs_txt <- if (length(fml) == 3L) {

        paste(deparse(fml[[3]]), collapse = " ")

      } else {

        paste(deparse(fml[[2]]), collapse = " ")

      }

    } else if (is.character(fml) && length(fml) == 1L) {

      txt <- trimws(fml)

      if (grepl("~", txt, fixed = TRUE)) {

        parts <- strsplit(txt, "~", fixed = TRUE)[[1]]

        rhs_txt <- trimws(parts[length(parts)])

      } else {

        rhs_txt <- txt

      }

    } else {

      stop("Invalid formula input: ", deparse(fml))

    }


    as.formula(paste(response_name, "~", rhs_txt), env = parent.frame())

  }


  if(copula.family %in% c("t", "T", "Student")){two_par_cop=TRUE} else {two_par_cop=FALSE}

  included_parameters <- c(names(margin.family$parameters), if(two_par_cop) c("theta","zeta") else c("theta"))


  formulas=list()

  for (parameter in included_parameters) {

    formulas[[parameter]]=get(paste(parameter,"formula",sep="."))

  }


  if (!requireNamespace("mgcv", quietly = TRUE)) {

    stop("Package 'mgcv' is required to construct smooth-term model matrices.")

  }


  formulas[["mu"]] <- as.formula(mu.formula)

  for (parameter in included_parameters[2:length(included_parameters)]) {

    formulas[[parameter]] <- to_response_formula(formulas[[parameter]], response_name = "response")

  }


  mm_x=list()

  mm_s=list()


  smooth_eval_env <- new.env(parent = baseenv())

  smooth_eval_env$s <- mgcv::s


  for(parameter in included_parameters) {

    data_for_par <- if(parameter %in% c("theta","zeta")) {

      dataset_mm[dataset_mm$time %in% unique(dataset_mm$time)[1:(length(unique(dataset_mm$time))-1)], , drop = FALSE]

    } else {

      dataset_mm

    }


    data_for_par <- sanitize_for_gamlss2(data_for_par, formulas[[parameter]])

    formula_terms <- stats::terms(formulas[[parameter]])

    has_intercept <- as.integer(attr(formula_terms, "intercept")) == 1L

    term_labels <- attr(formula_terms, "term.labels")


    # Keep all non-smooth RHS terms so model.matrix can expand interactions

    # like `time*gender` into main effects + interaction columns.

    fixed_terms <- term_labels[!grepl("^\\s*s\\(", term_labels)]


    if(length(fixed_terms) > 0 || has_intercept) {

      fixed_formula <- if (length(fixed_terms) == 0L && has_intercept) {

        stats::as.formula("~ 1")

      } else {

        stats::reformulate(termlabels = fixed_terms, intercept = has_intercept)

      }

      X_fixed <- stats::model.matrix(fixed_formula, data = data_for_par)

      fixed_assign <- attr(X_fixed, "assign")

      fixed_term_labels <- attr(stats::terms(fixed_formula), "term.labels")

      colnames(X_fixed) <- sub("^\\(Intercept\\)$", "intercept", colnames(X_fixed))

      colnames(X_fixed) <- normalize_time_covariate_colnames(colnames(X_fixed))

      mm_x[[parameter]] <- as.data.frame(X_fixed, check.names = FALSE)

      attr(mm_x[[parameter]], "assign") <- fixed_assign

      attr(mm_x[[parameter]], "term.labels") <- fixed_term_labels

    } else {

      mm_x[[parameter]] <- data.frame(row.names = seq_len(nrow(data_for_par)))

    }


    smooth_terms <- term_labels[grepl("^\\s*s\\(", term_labels)]

    if(length(smooth_terms) == 0) {

      mm_s[[parameter]] <- NULL

    } else {

      mm_s[[parameter]] <- list()

      for (s_label in smooth_terms) {

        s_txt <- trimws(s_label)

        s_call <- tryCatch(parse(text = s_txt)[[1]], error = function(e) NULL)

        s_obj <- eval(parse(text = s_txt), envir = smooth_eval_env)

        s_con <- mgcv::smoothCon(s_obj, data = data_for_par, knots = NULL, absorb.cons = TRUE)

        if (length(s_con) > 0 && !is.null(s_con[[1]]$X)) {

          B_s <- s_con[[1]]$X

          # Store the basis-specific penalty matrix returned by smoothCon so the

          # optimizer can use it instead of a generic second-difference fallback.

          if (!is.null(s_con[[1]]$S) && length(s_con[[1]]$S) > 0) {

            attr(B_s, "penalty") <- s_con[[1]]$S[[1]]

          }

          if (!is.null(s_call) && length(s_call) >= 2) {

            x_expr <- s_call[[2]]

            x_var <- trimws(gsub("`", "", paste(deparse(x_expr), collapse = " "), fixed = TRUE))

            x_value <- tryCatch(eval(x_expr, envir = data_for_par, enclos = parent.frame()), error = function(e) NULL)

            if (!is.null(x_value) && length(x_value) == nrow(B_s)) {

              attr(B_s, "smooth_x") <- as.numeric(x_value)

              attr(B_s, "smooth_var") <- x_var

            } else if (nzchar(x_var)) {

              attr(B_s, "smooth_var") <- x_var

            }

          }

          mm_s[[parameter]][[s_txt]] <- B_s

        }

      }

      if(length(mm_s[[parameter]]) == 0) {

        mm_s[[parameter]] <- NULL

      }

    }

  }



  mm=list(mm_x,mm_s)

  names(mm)=c("x","s")

  return(mm)

}


#' Calculate eta, eta inverse and eta derivative based on the given parameters and model matrices

#'

#' This function calculates the linear predictors (eta) for each parameter

#' based on the given covariate parameters (par_cov) and model matrices (mm).

#' It also computes the inverse link function (eta_inv) and the derivative of the link function (eta_dr)

#' for each parameter using the specified marginal distribution and copula link functions.

#'

#' @param par_cov A named vector of covariate parameters for each model term.

#' @param mm A list containing model matrices for fixed effects (mm$x) and smooth terms (mm$s).

#' @param margin_dist A list of functions for the marginal distribution, including link inverse and derivative functions.

#' @param copula_link A list of functions for the copula link, including link inverse and derivative functions.

#' @param par_s A list of smooth term parameters for each model parameter (optional).

#'

#' @return A list containing:

#' \item{eta}{A list of linear predictors for each parameter.}

#' \item{eta_inv}{A list of inverse link function values for each parameter.}

#' \item{eta_dr}{A list of derivatives of the link function for each parameter.}

#'

#' @keywords internal

#' @noRd

calc_eta=function(par_cov,mm,margin_dist,copula_link,par_s=NA) {

  eta=list()

  #par_s=list(); par_s[["mu"]]=list(); par_s[["sigma"]]=list(); par_s[["nu"]]=list(); par_s[["tau"]]=list()

  #par_s[["mu"]][["s(age)"]]=matrix(rep(1,ncol(mm$s[["mu"]][["s(age)"]])),ncol=1)

  for (par_name in names(mm$x)) {

    par_cov_single=par_cov[grepl(par_name,names(par_cov))]

    mm_temp=mm$x[[par_name]]

    #If there are no smooth terms for the parameter then just do standard calculation

    if(all(is.na(par_s[[par_name]]))) {

      eta[[par_name]]=rowSums(mm_temp * matrix(rep(par_cov_single,each=nrow(mm_temp)),ncol=length(par_cov_single),dimnames=list(NULL,c(names(par_cov_single)))))

    } else {

      eta[[par_name]]=

        rowSums(mm_temp * matrix(rep(par_cov_single,each=nrow(mm_temp)),ncol=length(par_cov_single),dimnames=list(NULL,c(names(par_cov_single)))))

      for (s_name in names(mm$s[[par_name]])) {

        eta[[par_name]]=eta[[par_name]] + mm$s[[par_name]][[s_name]] %*% par_s[[par_name]][[s_name]]

      }

    }

  }

  #Get link transforms (eta) and derivatives w.r.t to link for parameters

  eta_dr=eta_inv=list()

  for (par_name in names(mm$x)) {

    if(par_name %in% c("mu","sigma","nu","tau")) {

      eta_inv[[par_name]]=margin_dist[[paste(par_name,".linkinv",sep="")]](eta[[par_name]])

      eta_dr[[par_name]]=margin_dist[[paste(par_name,".dr",sep="")]](eta[[par_name]])

    }

    if(par_name %in% c("theta","zeta")) {

      eta_inv[[par_name]]=copula_link[[paste(par_name,".linkinv",sep="")]](eta[[par_name]])

      eta_dr[[par_name]]=copula_link[[paste(par_name,".dr",sep="")]](eta[[par_name]])

    }

  }

  return(list(eta=eta,eta_inv=eta_inv,eta_dr=eta_dr))

}


.calc_eta_rs_cached <- function(

  rs_design_cache,

  par_cov,

  par_s,

  margin_dist,

  copula_link,

  update_only = NULL,

  eta_out = NULL

) {

  if (is.null(eta_out)) {

    eta_out <- list(eta = list(), eta_inv = list(), eta_dr = list())

  }


  par_names <- names(rs_design_cache)

  if (!is.null(update_only)) {

    par_names <- intersect(update_only, par_names)

  }


  for (par_name in par_names) {

    design_info <- rs_design_cache[[par_name]]

    X <- design_info$X

    beta <- numeric(ncol(X))

    names(beta) <- colnames(X)


    fixed_names <- design_info$fixed_names

    beta[fixed_names] <- par_cov[fixed_names]


    if (length(par_s[[par_name]]) > 0) {

      smooth_beta <- unlist(par_s[[par_name]], use.names = FALSE)

      smooth_names <- setdiff(colnames(X), fixed_names)

      if (length(smooth_beta) != length(smooth_names)) {

        stop("Smooth coefficient length does not match cached design columns for ", par_name, ".", call. = FALSE)

      }

      beta[smooth_names] <- smooth_beta

    }


    eta_vec <- as.numeric(X %*% beta)

    eta_out$eta[[par_name]] <- eta_vec


    if (par_name %in% c("mu", "sigma", "nu", "tau")) {

      eta_out$eta_inv[[par_name]] <- margin_dist[[paste(par_name, ".linkinv", sep = "")]](eta_vec)

      eta_out$eta_dr[[par_name]] <- margin_dist[[paste(par_name, ".dr", sep = "")]](eta_vec)

    } else if (par_name %in% c("theta", "zeta")) {

      eta_out$eta_inv[[par_name]] <- copula_link[[paste(par_name, ".linkinv", sep = "")]](eta_vec)

      eta_out$eta_dr[[par_name]] <- copula_link[[paste(par_name, ".dr", sep = "")]](eta_vec)

    }

  }


  eta_out

}


