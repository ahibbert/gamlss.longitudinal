#' Variance-covariance matrix for a fitted longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param par Optional parameter list for evaluating uncertainty away from the
#'   fitted coefficients.
#' @param sep_d2 Logical legacy argument retained for compatibility.
#' @param numderiv Logical; use the numerical Hessian path.
#' @param method Character; Hessian method to use. `"analytical"` (default)
#'   uses the semi-analytical Hessian from `R/hessian-analytical.R`.
#'   `"numderiv"` uses full finite-difference numerical second derivatives as a
#'   slower reference path. The legacy `numderiv` logical argument is still
#'   accepted and maps to `method = "numderiv"` when `TRUE`.
#' @param progress Logical; show progress bars for slow Hessian calculations.
#' @param h Numeric finite-difference step used by the analytical Hessian
#'   helper.
#' @param ... Additional arguments, currently unused.
#'
#' @return A list containing variance-covariance matrices and standard errors.
#' @export
vcov.gamlss.longitudinal=function(object,par=NA,sep_d2=TRUE,numderiv=FALSE,

                                   method=c("analytical","numderiv","analytical_only"),

                                   progress=interactive(), h=1e-4, ...) {


  #object=fit; par=NA; numderiv=TRUE; sep_d2=TRUE


  method <- match.arg(method)

  method_requested <- method

  method_used <- method

  # Legacy: numderiv=TRUE overrides method selection

  if (isTRUE(numderiv)) method <- "numderiv"

  method_requested <- method

  method_used <- method


  progress = isTRUE(progress)


  include_dlcopdpar=TRUE

  response=object$response

  response_margin=object$response_margin

  response_subject=object$response_subject


  margin_names=unique(object$response_margin)

  num_margins=length(margin_names)


  #se_out=object$par*0;

  margin_dist=object$margin_dist; copula_dist=object$copula_dist; copula_link=get_copula_dist(copula_dist)$copula_link

  mm=object$model_matrix


  if(all(is.na(par))) {

    par_cov=object$par

    par_s=object$par_s

  } else {

    par_cov=par$par

    par_s=par$par_s

  }


  eta_out=calc_eta(par_cov,mm,margin_dist,copula_link,par_s=par_s)

  eta_inv=eta_out$eta_inv; eta_dr=eta_out$eta_dr; eta=eta_out$eta; eta_dr=eta_out$eta_dr


  if (identical(method, "analytical") &&

      identical(as.character(margin_dist$family[1]), "GG") &&

      "nu" %in% names(eta_inv)) {

    nu_abs_min <- suppressWarnings(min(abs(as.numeric(eta_inv$nu[is.finite(eta_inv$nu)])), na.rm = TRUE))

    if (is.finite(nu_abs_min) && nu_abs_min < 0.06) {

      warning(

        sprintf(

          paste(

            "Analytical Hessian for GG may be numerically unstable because fitted",

            "nu is close to 0 (min |nu| = %.4g); falling back to numerical Hessian."

          ),

          nu_abs_min

        ),

        call. = FALSE

      )

      method <- "numderiv"

      method_used <- "numderiv"

    }

  }


  #if(!all(is.na(par))) {response=eta_inv[["mu"]]}

  calc_lik_out=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=TRUE

            ,response=response,response_margin=response_margin,response_subject = response_subject)


  Fx_1_2=calc_lik_out$Fx_1_2; margin_p=calc_lik_out$margin_p; margin_d=calc_lik_out$margin_d; copula_d=calc_lik_out$copula_d


  if (method %in% c("analytical", "analytical_only") &&

      identical(calc_lik_out$likelihood_type, "discrete_rectangle")) {

    zero_fraction <- mean(response == 0, na.rm = TRUE)

    msg <- paste0(

      "Analytical Hessian is not yet implemented for exact discrete rectangle likelihoods; ",

      "falling back to numerical Hessian."

    )

    if (is.finite(zero_fraction) && zero_fraction >= 0.35) {

      msg <- paste0(

        "Analytical Hessian for zero-heavy discrete margins may be numerically delicate; ",

        msg

      )

    }

    if (identical(method, "analytical_only")) {

      stop(msg, call. = FALSE)

    }

    warning(msg, call. = FALSE)

    method <- "numderiv"

    method_used <- "numderiv"

  }


  ###Calculate derivaties: margin and copula d1 and d2

  margin_derivatives=calc_lik_out$margin_deriv

  copula_derivatives=calc_copula_derivatives(eta_inv, Fx_1_2, copula_dist,calc_d2 = TRUE)


  #Calculate numerical derivatives for F(x) - also used as a reference for overall likelihood derivative

  nd_impact_F=calc_Fx_derivatives(eta_inv,mm$x,margin_dist,response)

  nd_impact_F2=calc_Fx2_derivatives(eta_inv,mm$x,margin_dist,response)


  solve_hessian_vcov <- function(H) {

    vc <- -solve(H)

    se <- sqrt(abs(diag(solve(H))))

    if(!is.matrix(vc) || any(!is.finite(vc)) || any(!is.finite(se))) {

      stop("Hessian inversion produced non-finite variance-covariance values.", call. = FALSE)

    }

    eig <- tryCatch(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)

    list(

      vcov = vc,

      se = se,

      hessian_diagnostics = list(

        condition_number = tryCatch(kappa(H), error = function(e) NA_real_),

        min_abs_eigen = if (any(is.finite(eig))) min(abs(eig[is.finite(eig)])) else NA_real_,

        max_abs_eigen = if (any(is.finite(eig))) max(abs(eig[is.finite(eig)])) else NA_real_

      )

    )

  }


  if (method == "numderiv") {

    #nd2_joint_lik=calc_true_SE_numderiv_only(eta_inv,mm,margin_dist,response,testing=TRUE,response_margin,response_subject)

    hessian_nd=calc_true_SE_numderiv_only_covariates(object=object,par=par_cov,mm=mm$x,margin_dist=margin_dist,response=response,testing=FALSE,response_margin=response_margin,response_subject=response_subject,progress=progress)

  } else if (method %in% c("analytical", "analytical_only")) {

    # Source the analytical Hessian helpers if not already loaded.

    if (!exists("calc_analytical_hessian", mode = "function")) {

      # Try relative to this file, then working directory.

      hessian_files <- c(
        "hessian-setup.R",
        "hessian-margin-cdf.R",
        "hessian-copula.R",
        "hessian-assembly.R",
        "hessian-analytical.R"
      )

      candidate_dirs <- c(
        dirname(attr(body(vcov.gamlss.longitudinal), "srcfile")$filename %||% ""),
        "R",
        file.path(getwd(), "R")
      )

      loaded <- FALSE

      for (dir in candidate_dirs) {

        candidate_paths <- file.path(dir, hessian_files)

        if (all(file.exists(candidate_paths))) {
          for (cp in candidate_paths) source(cp, local = FALSE)
          loaded <- TRUE
          break
        }

      }

      if (!loaded) stop("Cannot locate hessian helper files. Source them manually or ensure the working directory is the package root.")

    }

    analytical_hessian <- tryCatch(

      calc_analytical_hessian(object, progress = progress, h = h),

      error = function(e) structure(list(error = conditionMessage(e)), class = "vcov_hessian_error")

    )

    analytical_vcov <- if(inherits(analytical_hessian, "vcov_hessian_error")) {

      analytical_hessian

    } else {

      tryCatch(

        solve_hessian_vcov(analytical_hessian),

        error = function(e) structure(list(error = conditionMessage(e)), class = "vcov_hessian_error")

      )

    }

    if(inherits(analytical_vcov, "vcov_hessian_error")) {

      if(identical(method, "analytical_only")) {

        stop("Analytical Hessian vcov failed: ", analytical_vcov$error, call. = FALSE)

      }

      warning(

        "Analytical Hessian vcov failed; falling back to numerical Hessian. Reason: ",

        analytical_vcov$error,

        call. = FALSE

      )

      method <- "numderiv"

      method_used <- "numderiv"

      hessian_nd=calc_true_SE_numderiv_only_covariates(object=object,par=par_cov,mm=mm$x,margin_dist=margin_dist,response=response,testing=FALSE,response_margin=response_margin,response_subject=response_subject,progress=progress)

    } else {

      method_used <- "analytical"

      hessian_nd <- analytical_hessian

      vcov_final <- analytical_vcov$vcov

      se_final <- analytical_vcov$se

      hessian_diagnostics <- analytical_vcov$hessian_diagnostics

    }

  } else {


    #######to delete############

    nd_impact_C2=calc_Fx2_derivatives(eta_inv,mm,margin_dist,response,testing=TRUE,response_margin,response_subject)[[2]]


    ### MARGIN LIKELIHOOD DERIVATIVES

    margin_deriv_subnames=c("m","d","v","t")

    names(margin_deriv_subnames)=c("mu","sigma","nu","tau")

    margin_par=names(mm)[names(mm) %in% c("mu","sigma","nu","tau")]

    order_margin=cbind(object$response_margin,object$response_subject)

    colnames(order_margin)=c("time","subject")


    dcdu1=copula_derivatives$dcdu1; dcdu2=copula_derivatives$dcdu2;d2cdu12=copula_derivatives$d2cdu12;d2cdu22=copula_derivatives$d2cdu22


    #####################For each parameter...

    d1_cop=d2_cop=matrix(0,nrow=length(response),ncol=length(margin_par))

    colnames(d1_cop)=colnames(d2_cop)=margin_par

    pair_cache_diag <- build_copula_pair_cache(response, response_margin, response_subject)

    for (par_name in margin_par) {

      if(object$include_dlcopdpar==TRUE | include_dlcopdpar==TRUE) {


        order_copula=data.frame()

        for (i in 1:(num_margins-1)) {

          order_copula=rbind(order_copula,cbind(order_margin[response_margin == margin_names[i],c("time","subject")],order_margin[response_margin == margin_names[i+1],c("time","subject")]))

        }

        colnames(order_copula)=c("time1","subject1","time2","subject2")


        margin_deriv_1=matrix(0,ncol=length(margin_par),nrow=length(response))

        colnames(margin_deriv_1)=paste("dld",margin_par,sep="")

        margin_deriv_1[,paste("dld",par_name,sep="")]=margin_derivatives[grepl("dld",names(margin_derivatives))][[which(margin_par==par_name)]]


        #COPULA DERIVS WITH RESPECT TO


        mu=eta_inv[["mu"]]

        F_nd=nd_impact_F[[par_name]]

        F_nd2=nd_impact_F2[[par_name]]

        c_nd2=nd_impact_C2[[par_name]]


        margin_components=cbind(order_margin,response,margin_p,margin_d,margin_deriv_1,mu,F_nd,F_nd2)

        margin_components_Ft_plus=margin_components

        margin_components_Ft_plus[,"time"]=normalize_lag_time(margin_components_Ft_plus[,"time"])

        margin_plus=merge(margin_components,margin_components_Ft_plus,by=c("time","subject"),all.x=TRUE)


        copula_components=cbind(

          order_copula,

          row_id1=pair_cache_diag$row_id1,

          row_id2=pair_cache_diag$row_id2,

          dcdu1,

          dcdu2,

          copula_d,

          d2cdu12,

          d2cdu22,

          c_nd2

        )

        copula_merged=merge(copula_components,margin_plus,by.x=c("time1","subject1"),by.y=c("time","subject"),all.x=TRUE)


        #Calculate copula derivative with respect to marginal parameters

        input=copula_merged

        d1_cop[,par_name]=calc_deriv_copula_wrt_margin(input,margin_par,par_name,calc_d2=FALSE)[,which(margin_par==par_name)]


        #OK so let's calcute the numerical d2lcopdpar and pass it through input

        d2_cop[,par_name]=calc_deriv_copula_wrt_margin(input,margin_par,par_name,calc_d2=TRUE)[,which(margin_par==par_name)]

      }

    }


    ###########Need d1 and d2 for score function


    m_d1_names=names(margin_derivatives)[grepl("dld",names(margin_derivatives))]

    c_d1_names=names(copula_derivatives)[grepl("dld",names(copula_derivatives))]


    m_d2_names=names(margin_derivatives)[grepl("d2ld",names(margin_derivatives))]

    c_d2_names=names(copula_derivatives)[grepl("d2ld",names(copula_derivatives))]


    d1_all=list(); d2_all=list()


    i=1

    for(par_name in c(m_d1_names)) {

      d1_all[[par_name]]=c(margin_derivatives[[par_name]])

      if((object$include_dlcopdpar==TRUE | include_dlcopdpar==TRUE)) {

        d1_all[[par_name]]=c(margin_derivatives[[par_name]]+d1_cop[,i])

        i=i+1

      }

    }


    for(par_name in c(c_d1_names)) {

      d1_all[[par_name]]=c(copula_derivatives[[par_name]])

    }


    names(d1_all)=names(mm)


    i=1

    for(par_name in c(m_d2_names)) {

      d2_all[[par_name]]=c(margin_derivatives[[par_name]])

      if((object$include_dlcopdpar==TRUE | include_dlcopdpar==TRUE) & endsWith(par_name,"2")) {

        d2_all[[par_name]]=c(margin_derivatives[[par_name]]+d2_cop[,i]*(if(sep_d2==TRUE) {0} else {1}))

        i=i+1

      }

    }

    for(par_name in c(c_d2_names)) {

      d2_all[[par_name]]=c(copula_derivatives[[par_name]])

    }


    d2_all_mean=rep(0,length=length(d2_all))

    names(d2_all_mean)=names(d2_all)

    for (deriv_name in names(d2_all)) {

      d2_all_mean[deriv_name]=mean(d2_all[[deriv_name]])

    }


    d2_mat_diag=d2_all_mean[endsWith(names(d2_all_mean),"2")]

    d2_mat_cross=d2_all_mean[!endsWith(names(d2_all_mean),"2")]

    d2_mat=matrix(nrow=length(eta),ncol=length(eta))

    #print(d2_mat);print(names(eta))

    colnames(d2_mat)=rownames(d2_mat)=names(eta)


    copula_deriv_subnames=c("th","z")

    names(copula_deriv_subnames)=c("theta","zeta")

    all_names=c(margin_deriv_subnames,copula_deriv_subnames)

    sub_names_in=all_names[names(eta)]

    print(sub_names_in)


    for (row_name in rownames(d2_mat)) {

      for (col_name in colnames(d2_mat)) {

        if(!row_name==col_name) {

          deriv_name_temp=paste("d2ld",sub_names_in[row_name],"d",sub_names_in[col_name],sep="")


          if(is.na(d2_all_mean[deriv_name_temp])) {

            deriv_name_temp=paste("d2ld",sub_names_in[col_name],"d",sub_names_in[row_name],sep="")

          }


          deriv_val_temp=d2_all_mean[deriv_name_temp]

          d2_mat[row_name,col_name]=deriv_val_temp

        }

      }

    }


    #cop_row=(grepl("theta",rownames(d2_mat))|grepl("zeta",rownames(d2_mat)))

    #d2_mat[!cop_row,!cop_row][upper.tri(d2_mat[!cop_row,!cop_row])]=d2_mat_cross

    diag(d2_mat)=d2_mat_diag


    d2_mat[is.na(d2_mat)]=0


  }


  if (method %in% c("numderiv", "analytical", "analytical_only")) {

    if(!exists("vcov_final", inherits = FALSE) || !exists("se_final", inherits = FALSE)) {

      vcov_solved <- solve_hessian_vcov(hessian_nd)

      vcov_final <- vcov_solved$vcov

      se_final <- vcov_solved$se

      hessian_diagnostics <- vcov_solved$hessian_diagnostics

    }

  } else {

    vcov_final = -(solve((d2_mat)))/(length(response))

    se_final=sqrt(abs(diag(vcov_final)))

  }


  ###########TESTING APPROACH FOR ESTIMATING SMOOTHER VARIANCE


  # Method 1: Bayesian/Mixed Model Approach for Smoother Variance

  # Calculate variance-covariance matrix for smooth terms using: Var(beta) = (X'WX + lambda*P)^(-1) * sigma^2


  smooth_vcov_list = list()

  smooth_se_list = list()


  # Extract residual variance estimate (using reciprocal of mean weights as proxy for sigma^2)

  if(!is.null(object$weights) && length(object$weights) > 0 && is.numeric(object$weights)) {

    sigma2_est = 1 / mean(object$weights, na.rm = TRUE)

  } else {

    # Fallback: estimate from residuals if weights not available

    fitted_response = eta_inv[["mu"]]

    residuals = response - fitted_response

    sigma2_est = var(residuals, na.rm = TRUE)

  }


  # Ensure scalar numeric to avoid deprecated array recycling warnings.

  sigma2_est = as.numeric(sigma2_est)[1]

  if(!is.finite(sigma2_est)) {

    sigma2_est = 1

  }


  # Process each parameter that has smooth terms

  for(par_name in names(object$par_s)) {

    if(length(object$par_s[[par_name]]) > 0) {


      smooth_vcov_list[[par_name]] = list()

      smooth_se_list[[par_name]] = list()


      # Process each smooth term for this parameter

      for(s_name in names(object$par_s[[par_name]])) {


        # Get the B-spline basis matrix

        B = object$model_matrix$s[[par_name]][[s_name]]


        # Get the smoothing parameter

        lambda = object$lambda_s[[par_name]][[s_name]]


        # Use the mgcv-generated penalty stored on the basis matrix; fall back to

        # a generic second-difference penalty only when unavailable.

        k = ncol(B)

        pen_attr = attr(B, "penalty")

        if (!is.null(pen_attr) && is.matrix(pen_attr) &&

            nrow(pen_attr) == k && ncol(pen_attr) == k) {

          P = pen_attr

        } else if (k > 2) {

          D2 = diff(diag(k), differences = 2)

          P = t(D2) %*% D2

        } else {

          P = diag(k)

        }


        # Get per-parameter IRLS working weights. object$weights is a named list

        # keyed by parameter name; fall back to unit weights if not available.

        w_par = object$weights[[par_name]]

        if (!is.null(w_par) && is.numeric(w_par) && length(w_par) == nrow(B)) {

          w_diag = as.vector(w_par)

        } else {

          w_diag = rep(1, nrow(B))

        }

        W = diag(w_diag)


        # Per-parameter sigma2: scale consistent with IRLS, 1/mean(w)

        sigma2_par = if (all(w_diag > 0)) 1 / mean(w_diag) else sigma2_est


        # Calculate the penalized precision matrix: X'WX + lambda*P

        XWX = t(B) %*% W %*% B

        penalized_precision = XWX + lambda * P


        # Variance-covariance matrix for this smooth: (X'WX + lambda*P)^(-1) * sigma^2

        tryCatch({

          smooth_vcov = solve(penalized_precision) * sigma2_par

          smooth_se = sqrt((diag(smooth_vcov)))


          # Store results

          smooth_vcov_list[[par_name]][[s_name]] = smooth_vcov

          smooth_se_list[[par_name]][[s_name]] = smooth_se


          # Also calculate the smoother matrix for fitted values variance

          # A = X(X'WX + lambda*P)^(-1)X'W

          smoother_matrix = B %*% solve(penalized_precision) %*% t(B) %*% W

          fitted_se = sqrt(abs(as.vector(diag(smoother_matrix))) * sigma2_par)


          cat(sprintf("\nSmooth term variance estimates for %s:%s\n", par_name, s_name))

          cat(sprintf("  Basis coefficients SE: min=%.4f, max=%.4f, mean=%.4f\n",

                     min(smooth_se), max(smooth_se), mean(smooth_se)))

          cat(sprintf("  Fitted values SE: min=%.4f, max=%.4f, mean=%.4f\n",

                     min(fitted_se), max(fitted_se), mean(fitted_se)))

          cat(sprintf("  Effective DF: %.2f (trace of smoother matrix)\n",

                     sum(diag(smoother_matrix))))

          cat(sprintf("  Smoothing parameter lambda: %.4f\n", lambda))


        }, error = function(e) {

          warning(sprintf("Could not calculate variance for smooth %s:%s - %s",

                         par_name, s_name, e$message))

          smooth_vcov_list[[par_name]][[s_name]] = NULL

          smooth_se_list[[par_name]][[s_name]] = NULL

        })

      }

    }

  }


  # Add smooth variance results to return list

  vcov_final_with_smooth = list(

    overall = vcov_final,

    smooth_vcov = smooth_vcov_list,

    smooth_se = smooth_se_list

  )


  se_final_with_smooth = list(

    overall = se_final,

    smooth_se = smooth_se_list

  )


  return(list(

    vcov=vcov_final_with_smooth,

    se=se_final_with_smooth,

    method=method_used,

    method_requested=method_requested,

    hessian_diagnostics=hessian_diagnostics %||% NULL

  ))


}


.can_use_cached_vcov <- function(object, numderiv = FALSE, method = NULL, extra_args = list()) {

  if(!inherits(object, "gamlss.longitudinal")) return(FALSE)

  extra_args_cache <- extra_args

  extra_args_cache$method <- NULL

  if(!is.null(extra_args_cache) && length(extra_args_cache) > 0) return(FALSE)

  if(is.null(object$vcov) || !is.list(object$vcov)) return(FALSE)

  if(is.null(object$vcov$vcov) || is.null(object$vcov$vcov$overall)) return(FALSE)


  if(!is.null(object$vcov_meta) && !is.null(object$vcov_meta$numderiv)) {

    numderiv_ok <- identical(isTRUE(object$vcov_meta$numderiv), isTRUE(numderiv))

    method_ok <- is.null(method) ||

      is.null(object$vcov_meta$method) ||

      identical(as.character(object$vcov_meta$method)[1], as.character(method)[1])

    return(numderiv_ok && method_ok)

  }


  TRUE

}


.resolve_vcov <- function(object, numderiv = FALSE, extra_args = list()) {

  vcov_method <- extra_args$method %||% if (isTRUE(numderiv)) "numderiv" else "analytical"

  if(.can_use_cached_vcov(object, numderiv = numderiv, method = vcov_method, extra_args = extra_args)) {

    return(object$vcov)

  }


  if(is.null(extra_args$method)) {

    extra_args$method <- vcov_method

  }


  do.call(

    vcov.gamlss.longitudinal,

    c(list(object = object, numderiv = numderiv), extra_args)

  )

}


#' Summarize a fitted gamlss.longitudinal model

#'

#' Creates a compact summary of key model diagnostics and coefficient estimates

#' for a fitted `gamlss.longitudinal` object. By default it wraps

#' `vcov.gamlss.longitudinal()` to provide standard errors and confidence

#' intervals for fixed effects.

#'

#' @param object A fitted object of class `gamlss.longitudinal`.

#' @param include_vcov Logical; if `TRUE`, compute and include variance-covariance

#'   output via `vcov.gamlss.longitudinal()`.

#' @param numderiv Logical passed to `vcov.gamlss.longitudinal()`.

#' @param ci_level Confidence level for coefficient intervals.

#' @param ... Additional arguments passed to `vcov.gamlss.longitudinal()`.

#'

#' @return An object of class `summary.gamlss.longitudinal` containing:

#' - model dimensions and parameter counts,

#' - likelihood and information criteria,

#' - fixed-effect coefficient table,

#' - smooth-term table with effective degrees of freedom,

#' - optional `vcov` output.

#' @export
