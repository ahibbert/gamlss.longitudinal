###########NEW SIMPLIFIED FUNCTIONS

#' Fit a longitudinal joint regression model
#'
#' This function fits a longitudinal model to a dataset with gamlss margins
#' and copula fit to dependence. Any linear or factor covariates can be fit
#' to any parameters of the copula or margin distributions. The model is fit
#' using RS() optimisation and the joint likelihood by default. Select
#' use dlcopdpar=FALSE to fit separately optimised models for the margin and
#' copula likelihoods which can be quicker with a slight loss to overall fit.
#'
#' @param margin_dist Marginal distribution specified as a gamlss family object,
#' e.g. GA(), NO(), PO(), NBI(), etc.
#' @param copula_dist Copula distribution specified as a string in VineCopula style
#'  e.g. "t", "C" etc.
#' @param mu.formula Formula for the mean parameter of the marginal distribution
#' @param sigma.formula Formula for the sigma parameter of the marginal distribution
#' @param nu.formula Formula for the nu parameter of the marginal distribution
#' @param tau.formula Formula for the tau parameter of the marginal distribution
#' @param theta.formula Formula for the theta parameter of the copula distribution
#' @param zeta.formula Formula for the zeta parameter of the copula distribution
#' @param include_dlcopdpar Include the derivative of the copula likelihood with respect
#' to the margin parameters in the joint likelihood.
#' @param inner_stop_crit Stopping criterion for the inner loop
#' @param outer_stop_crit Stopping criterion for the outer loop
#' @param start_step_size Initial step size for the backfitting algorithm
#' @param step_adjustment Step size adjustment factor
#' @param max_steps Maximum number of times for reducing the step size
#' @param start_from Starting values for the parameters if needed
#' @param verbose Level of output to the console 3 = ALL, 0 = Minimal
#' @param plot_results Plot the results of the optimisation
#' @param true_val True values for the parameters if known for plotting
#' @param method Optimisation method to use, RS() is the default
#' @param max_outer_iter Maximum number of outer iterations
#' @param max_inner_iter Maximum number of inner iterations
#' @param use_Rcpp Use Rcpp for matrix operations
#'
#' @return Returns matrix of log likleihoods, margin parameters and copula parameters
#' by iteration of the optimisation.
#' @export
gamlss.longitudinal=function(dataset,
                        margin_dist,
                        copula_dist,
                        time_var=NA,
                        subject_var=NA,
                        mu.formula = ("response ~ 1"),
                        sigma.formula = ("1"),
                        nu.formula = ("1"),
                        tau.formula = ("1"),
                        theta.formula=("1"),
                        zeta.formula=("1"),
                        include_dlcopdpar=FALSE,
                        inner_stop_crit=.1,
                        outer_stop_crit=.1,
                        start_step_size=.5,
                        step_adjustment=.5,
                        max_steps=5,
                        start_from=NA,
                        verbose=3,
                        plot_results=FALSE,
                        true_val=NA,
                        method="RS",
                        max_outer_iter=20,
                        max_inner_iter=20,
                        use_Rcpp=FALSE,
                        lambda_start=5,
                        lambda_penalty_K=2
                      )
{
  fit_start_time <- Sys.time()

  ##################### DATA CHECKS AND VALIDATION #####################

  # Save original dataset
  dataset_original <- dataset

  # Force plain data.frame (safe for tibble/data.table too)
  if (!is.data.frame(dataset)) {
    dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)
  } else {
    dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)
  }

  # Validate and prepare input data
  if(all(is.na(time_var)) || all(is.na(subject_var))) {
    stop("ERROR: Required input variables not specified.\n",
         "Please specify:\n",
         "  - time_var: column name for time/margin variable (e.g., 'time')\n",
         "  - subject_var: column name for subject ID variable (e.g., 'subject')\n",
         "Example: gamlss.longitudinal(..., time_var='time', subject_var='subject')")
  }

  # Validate dataset contains required columns
  if(!time_var %in% colnames(dataset)) {
    stop("ERROR: time_var='", time_var, "' not found in dataset.\n",
         "Available columns: ", paste(colnames(dataset), collapse=", "))
  }
  if(!subject_var %in% colnames(dataset)) {
    stop("ERROR: subject_var='", subject_var, "' not found in dataset.\n",
         "Available columns: ", paste(colnames(dataset), collapse=", "))
  }

  # Extract response variable name from mu formula
  mu_formula_obj <- as.formula(mu.formula)
  response_var <- all.vars(mu_formula_obj)[1]
  if (verbose > 1) print(paste("Identified response variable:", response_var))

  # Validate response variable exists
  if(!response_var %in% names(dataset)) {
    stop("ERROR: response variable '", response_var, "' not found in dataset.\n",
         "Available columns: ", paste(names(dataset), collapse=", "))
  }

  # Rename columns for internal use
  names(dataset)[names(dataset) == time_var] <- "time"
  names(dataset)[names(dataset) == subject_var] <- "subject"
  names(dataset)[names(dataset) == response_var] <- "response"

  if (any(is.na(dataset$time)) || any(is.na(dataset$subject))) {
    stop("ERROR: time and subject variables cannot contain NA values.")
  }

  # Robust formula normalizer:
  # - accepts formula objects
  # - accepts strings without "~" (e.g. "time + s(age)")
  # - for mu, adds response on LHS if missing
  normalize_formula <- function(fml, response_name = "response", require_lhs = FALSE) {
    if (inherits(fml, "formula")) {
      return(fml)
    }
    if (!is.character(fml) || length(fml) != 1 || is.na(fml) || nchar(trimws(fml)) == 0) {
      stop("ERROR: Invalid formula input: ", deparse(fml))
    }

    txt <- trimws(fml)

    if (!grepl("~", txt, fixed = TRUE)) {
      txt <- if (require_lhs) paste0(response_name, " ~ ", txt) else paste0("~ ", txt)
    } else if (require_lhs) {
      parts <- strsplit(txt, "~", fixed = TRUE)[[1]]
      lhs <- trimws(parts[1])
      rhs <- trimws(parts[2])
      if (nchar(lhs) == 0) txt <- paste0(response_name, " ~ ", rhs)
    }

    as.formula(txt, env = parent.frame())
  }

  # Variable-name translation from user names -> internal names
  translate_formula_vars <- function(fml, var_map, response_name = "response", require_lhs = FALSE) {
    f_obj <- normalize_formula(fml, response_name = response_name, require_lhs = require_lhs)
    f_txt <- paste(deparse(f_obj), collapse = " ")

    for (old_name in names(var_map)) {
      new_name <- var_map[[old_name]]
      if (!identical(old_name, new_name)) {
        old_esc <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", old_name)
        f_txt <- gsub(paste0("`", old_esc, "`"), new_name, f_txt, perl = TRUE)
        f_txt <- gsub(paste0("\\b", old_esc, "\\b"), new_name, f_txt, perl = TRUE)
      }
    }

    as.formula(f_txt, env = environment(f_obj))
  }

  var_map <- c()
  var_map[[time_var]] <- "time"
  var_map[[subject_var]] <- "subject"
  var_map[[response_var]] <- "response"

  mu.formula.int    <- translate_formula_vars(mu.formula,    var_map, response_name = "response", require_lhs = TRUE)
  sigma.formula.int <- translate_formula_vars(sigma.formula, var_map, response_name = "response", require_lhs = FALSE)
  nu.formula.int    <- translate_formula_vars(nu.formula,    var_map, response_name = "response", require_lhs = FALSE)
  tau.formula.int   <- translate_formula_vars(tau.formula,   var_map, response_name = "response", require_lhs = FALSE)
  theta.formula.int <- translate_formula_vars(theta.formula, var_map, response_name = "response", require_lhs = FALSE)
  zeta.formula.int  <- translate_formula_vars(zeta.formula,  var_map, response_name = "response", require_lhs = FALSE)

  if(verbose > 1) {
    cat("Input validation successful.\n")
    cat("Data dimensions:", nrow(dataset), "x", ncol(dataset), "\n")
    cat("Response variable:", response_var, "-> renamed to 'response'\n")
    cat("Time variable:", time_var, "-> renamed to 'time'\n")
    cat("Subject variable:", subject_var, "-> renamed to 'subject'\n")
    cat("Time points:", length(unique(dataset$time)), "\n")
    cat("Subjects:", length(unique(dataset$subject)), "\n")
  }

  # Validate that all subject/time combinations are unique
  subject_time_combo <- paste(dataset$subject, dataset$time, sep="_")
  if(length(subject_time_combo) != length(unique(subject_time_combo))) {
    duplicate_combos <- subject_time_combo[duplicated(subject_time_combo)]
    stop("ERROR: Duplicate subject/time combinations found.\n",
         "Each subject must have exactly one observation per time point.\n",
         "Duplicate combinations (first 10): ", 
         paste(unique(duplicate_combos)[1:min(10, length(unique(duplicate_combos)))], collapse=", "))
  }

  if(verbose > 1) {
    cat("Subject/time uniqueness check passed.\n")
    cat("Unique subject/time combinations:", length(unique(subject_time_combo)), "\n\n")
  }

  # Missingness summary by time and consecutive time pairs.
  time_levels <- sort(unique(dataset$time))
  n_time_levels <- length(time_levels)

  miss_by_time <- do.call(rbind, lapply(time_levels, function(ti) {
    idx <- dataset$time == ti
    n_total <- sum(idx)
    n_na <- sum(is.na(dataset$response[idx]))
    c(time = ti, n_total = n_total, n_na_response = n_na, n_observed_response = n_total - n_na)
  }))
  miss_by_time <- as.data.frame(miss_by_time)
  rownames(miss_by_time) <- NULL

  pair_summary <- data.frame(
    time1 = numeric(0),
    time2 = numeric(0),
    complete_pairs = integer(0),
    total_pairs = integer(0),
    total_observations = integer(0)
  )

  if (n_time_levels > 1) {
    for (i in seq_len(n_time_levels - 1)) {
      t1 <- time_levels[i]
      t2 <- time_levels[i + 1]
      d1 <- dataset[dataset$time == t1, c("subject", "response")]
      d2 <- dataset[dataset$time == t2, c("subject", "response")]
      names(d1) <- c("subject", "response_t1")
      names(d2) <- c("subject", "response_t2")
      merged_pair <- merge(d1, d2, by = "subject", all = FALSE)

      total_pairs <- nrow(merged_pair)
      complete_pairs <- sum(!is.na(merged_pair$response_t1) & !is.na(merged_pair$response_t2))

      pair_summary <- rbind(
        pair_summary,
        data.frame(
          time1 = t1,
          time2 = t2,
          complete_pairs = complete_pairs,
          total_pairs = total_pairs,
          total_observations = nrow(dataset)
        )
      )
    }
  }

  if (verbose > 0) {
    cat("Missingness Summary (by time):\n")
    print(miss_by_time)
    cat("\nConsecutive Pair Completeness:\n")
    if (nrow(pair_summary) > 0) {
      print(pair_summary)
    } else {
      cat("No consecutive time pairs available.\n")
    }
    cat("\n")
  }

  # Hard stop if any margin is 100% missing.
  margin_all_missing <- miss_by_time$n_observed_response == 0
  if (any(margin_all_missing)) {
    bad_times <- miss_by_time$time[margin_all_missing]
    stop(
      "ERROR: 100% missing response values detected for margin time point(s): ",
      paste(bad_times, collapse = ", "),
      "\nModel fitting stopped because at least one margin has no observed outcomes."
    )
  }

  # Hard stop if any consecutive copula pair has 0 complete pairs while pairs exist.
  if (nrow(pair_summary) > 0) {
    pair_all_missing <- pair_summary$total_pairs > 0 & pair_summary$complete_pairs == 0
    if (any(pair_all_missing)) {
      bad_pairs <- apply(pair_summary[pair_all_missing, c("time1", "time2"), drop = FALSE], 1, function(x) {
        paste0("(", x[1], ",", x[2], ")")
      })
      stop(
        "ERROR: 100% missing complete copula pairs detected for consecutive time pair(s): ",
        paste(bad_pairs, collapse = ", "),
        "\nModel fitting stopped because at least one copula pair contributes no complete observations."
      )
    }
  }

  #Setup model matrix from given formulas
  copula_link=get_copula_dist(copula_dist)$copula_link
  mm=suppressWarnings(create_model_matrices(
    mu.formula.int,
    sigma.formula.int,
    nu.formula.int,
    tau.formula.int,
    theta.formula.int,
    zeta.formula.int,
    margin.family = margin_dist,
    copula.family = copula_dist,
    copula.link = copula_link,
    dataset = dataset
  ))

  #Create vector of starting covariate values, currently starting at zero before first fit with the intercept as the mean
  if(all(is.na(start_from))) {
    par_eta=get_starting_values(copula_dist,margin_dist,dataset=dataset,eta_transform=TRUE)
    par_cov=as.numeric(vector())
    for (par_name in names(mm$x)) {
      par_cov_single=as.numeric(vector(length=length(colnames(mm$x[[par_name]]))))
      names(par_cov_single)=paste(par_name,colnames(mm$x[[par_name]]),sep=".")
      par_cov_single[1]=par_eta[par_name]
      if(length(par_cov_single)>1) {
        par_cov_single[2:length(par_cov_single)]=0
      }
      par_cov=c(par_cov,par_cov_single)
    }
  } else {
    par_cov=start_from
  }
  par_s=list()
  df_s=list()
  lambda_s=list()
  #names(par_s)=names(df_s)=names(lambda_s)=names(mm$x)
  for (par_name in names(mm$x)) {
    par_s[[par_name]]=list()
    df_s[[par_name]]=list()
    lambda_s[[par_name]]=list()
    for (s_name in names(mm$s[[par_name]])) {
      B=mm$s[[par_name]][[s_name]]
      par_s[[par_name]][[s_name]]=c(par_s[[par_name]][[s_name]],rep(0,ncol(B)))
      names(par_s[[par_name]][[s_name]])=paste(par_name,s_name,1:ncol(B),sep=".")
      df_s[[par_name]][[s_name]]=0; lambda_s[[par_name]][[s_name]]=lambda_start
      names(df_s[[par_name]][[s_name]])=names(lambda_s[[par_name]][[s_name]])=s_name
   }
  }
  #Starting parameters for fixed parameters: par_cov

  #Parameters used in optimisation loops
  first_outer_run=TRUE
  outer_log_lik_change=outer_start_log_lik=outer_end_log_lik=0
  log_lik_history=matrix(ncol=3,nrow=0)
  par_history=matrix(ncol=length(par_cov),nrow=0); colnames(par_history)=names(par_cov)
  outer_run_counter=1; outer_only_run_counter=1
  step_size=start_step_size
  weights_final=list()

  #OUTER ITERATION (MAIN LOOP)
  while ((first_outer_run==TRUE | (abs(outer_log_lik_change)>outer_stop_crit)) & outer_only_run_counter < max_outer_iter) {

    cat(paste("\nOUTER ITERATION:",outer_only_run_counter))
    first_outer_run=TRUE

    # RUN INNER ITERATION FOR EACH PARAMETER
    for (par_name in names(mm$x)) {

      if(verbose > 2) {
        cat(paste("\nINNER ITERATION: Parameter:",par_name))
      }

      first_inner_run=TRUE; change_log_lik=0; beta_change_inner=99
      run_counter=1
      inner_run_counter=1

      # INNER ITERATION (GLIM)
      while ( (first_inner_run==TRUE | abs(change_log_lik)>inner_stop_crit) & inner_run_counter<max_inner_iter) { #

        timer=c()
        timer_start=Sys.time()

        first_inner_run=FALSE

        eta_out=calc_eta(par_cov,mm,margin_dist,copula_link,par_s=if(first_inner_run) {NA} else {par_s})
        eta=eta_out$eta; eta_dr=eta_out$eta_dr; eta_inv=eta_out$eta_inv

        # Guard against silent row dropping in matrix construction when response has NAs.
        n_resp <- length(dataset$response)
        margin_params <- intersect(names(mm$x), c("mu", "sigma", "nu", "tau"))
        bad_lengths <- margin_params[sapply(margin_params, function(pn) length(eta_inv[[pn]]) != n_resp)]
        if (length(bad_lengths) > 0) {
          detail <- paste(sapply(bad_lengths, function(pn) {
            paste0(pn, "=", length(eta_inv[[pn]]), " vs response=", n_resp)
          }), collapse = ", ")
          stop(
            "ERROR: Parameter vector lengths do not match response length. ",
            detail,
            ".\nThis usually indicates model-matrix rows were dropped (often due to NA handling)."
          )
        }

        calc_lik_out=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=FALSE
          ,response=dataset$response,response_margin=(dataset$time),response_subject = dataset$subject)
        log_lik=calc_lik_out$log_lik; margin_d=calc_lik_out$margin_d; margin_p=calc_lik_out$margin_p; 
        margin_deriv=calc_lik_out$margin_deriv; copula_d=calc_lik_out$copula_d; copula_p=calc_lik_out$copula_p; 
        Fx_1_2=calc_lik_out$Fx_1_2;order_copula=calc_lik_out$order_copula

        if(first_outer_run==TRUE) {
          outer_start_log_lik=log_lik["joint"]; first_outer_run=FALSE
        }

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Calc Lik")

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Numerical Derivatives")

        #Capturing log lik and parameter estimates for each iteration
        log_lik_history=rbind(log_lik_history,calc_lik_out$log_lik)
        par_history=rbind(par_history,par_cov)

        #Fixing extreme values if they exist, though they shouldn't
        Fx_1_2[Fx_1_2>1]=1;Fx_1_2[Fx_1_2<0]=0

        ########CALCULATE COPULA DERIVATIVES
        copula_derivatives=calc_copula_derivatives(
          eta_inv,
          Fx_1_2,
          copula_dist,
          par1 = calc_lik_out$copula_par1,
          par2 = calc_lik_out$copula_par2,
          pair_complete = calc_lik_out$pair_complete
        )
        dldth=copula_derivatives$dldth; dcdth=copula_derivatives$dcdth; dcdu1=copula_derivatives$dcdu1; dcdu2=copula_derivatives$dcdu2
        if("zeta" %in% names(eta_inv)) {dldz=copula_derivatives$dldz; dcdz=copula_derivatives$dcdz}

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Copula Derivatives")

        ### Calculate copula derivatives w.r.t margin parameters
        if(!par_name %in% c("mu","sigma","nu","tau")) {
          if(par_name == "theta") {
            d1=as.matrix(dldth)
            colnames(d1)="dldtheta"
            #d2=d2ldth2
          } else if(par_name == "zeta") {
            d1=as.matrix(dldz)
            colnames(d1)="dldzeta"
            #d2=d2ldz2
          } else {
            stop("Unexpected copula parameter in optimisation: ", par_name)
          }
        } else {

          ### MARGIN LIKELIHOOD DERIVATIVES
          margin_deriv_subnames=c("m","d","v","t")
          names(margin_deriv_subnames)=c("mu","sigma","nu","tau")
          margin_par=names(mm$x)[names(mm$x) %in% c("mu","sigma","nu","tau")]
          response=dataset$response

          #Extract margin calculations for F(x), f(x), response and derivatives at time 1 and time 2, join to copula values for time 1 and time 2
          margin_deriv_1=matrix(0,ncol=length(margin_par),nrow=length(response))
          colnames(margin_deriv_1)=paste("dld",margin_par,sep="")
          #margin_deriv_2=margin_deriv_2cross
          margin_deriv_1[,paste("dld",par_name,sep="")]=margin_deriv[grepl("dld",names(margin_deriv))][[which(margin_par==par_name)]]
          #margin_deriv_2[,i]=margin_deriv[grepl("d2ld",names(margin_deriv))&endsWith(names(margin_deriv),"2")][[i]]

          d1=as.matrix(margin_deriv[grepl(paste("dld",margin_deriv_subnames[par_name],sep=""),names(margin_deriv))][[1]])
          colnames(d1)=paste("dld",par_name,sep="")

          if(include_dlcopdpar==TRUE) {

            #Calculate numerical derivatives for F(x) - also used as a reference for overall likelihood derivative
            nd_impact_F=calc_Fx_derivatives(eta_inv,mm$x,margin_dist,response=dataset$response)

            ### COPULA LIKELIHOOD DERIVATIVES
            order_margin=dataset[,c("time","subject")]
            mu=eta_inv[["mu"]]
            F_nd=nd_impact_F[[par_name]]

            margin_components=cbind(order_margin,response,margin_p,margin_d,margin_deriv_1,mu,F_nd)
            margin_components_Ft_plus=margin_components
            margin_components_Ft_plus$time=margin_components_Ft_plus$time-1
            margin_plus=merge(margin_components,margin_components_Ft_plus,by=c("time","subject"),all.x=TRUE)

            copula_components=cbind(order_copula,dcdu1,dcdu2,copula_d)

            copula_merged=merge(copula_components,margin_plus,by.x=c("time1","subject1"),by.y=c("time","subject"),all.x=TRUE)

            #Calculate copula derivative with respect to marginal parameters
            input=copula_merged
            d1_cop=calc_deriv_copula_wrt_margin(input,margin_par,par_name)
            d1_m=d1
            d1=d1_m+d1_cop
            #d1=d1*0+(nd_impact[par_name]/nrow(d1))
          }

          timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
          names(timer)[length(timer)]=paste("Margin Derivatives")

          if(verbose>=4) {print(timer)}

          d1=d1[,grepl(par_name,colnames(d1))]

          if(include_dlcopdpar==FALSE) {d1_cop=d1*0; d1_m=d1}

          #nd=round(c(nd_impact[par_name],nd_impact_m[par_name],nd_impact_c[par_name],sum(d1),sum(d1_m),sum(d1_cop*0.5)),2)
          #names(nd)=c("joint_nd","marginal_nd","copula_nd","joint_calc","margin_calc","copula_calc")
          #print(nd)
        }

        ### INNER ITERATION / BACKFITTING STEP

        # 1. Calculate y_k, w_k

        ########### FIRST ITERATION CALCULATES B_k without smooths
        score=score_function_v2(eta=eta[[par_name]],dldpar=d1,d2ldpar=-(d1*d1),dpardeta=eta_dr[[par_name]])

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Backfitting")

        # Setup model matrices
        X=as.matrix(mm$x[[par_name]])
        W=diag(as.vector(score$w_k))

        z_k=score$z_k
        for (s_name in names(mm$s[[par_name]])) {
          B=mm$s[[par_name]][[s_name]]
          #par_s[[par_name]][[s_name]]=c(par_s[[par_name]][[s_name]],rep(10,ncol(B)))
          colnames(B)=paste(par_name,s_name,1:ncol(B),sep=".")
          X=cbind(X,B)
        }
        if(length(par_s[[par_name]])==0) {
          paste("No smooths found for parameter; running basic IRLS",par_name)
          beta_start=c(par_cov[paste(paste(par_name,sep=" "),colnames(mm$x[[par_name]]),sep=".")])
        } else {
            ############# UNPENALISED VERSION
          temp_par_s_unlisted=unlist(par_s[[par_name]],use.names=TRUE)
          names(temp_par_s_unlisted)=colnames(X)[(ncol(mm$x[[par_name]])+1):ncol(X)]
          beta_start=c(par_cov[paste(paste(par_name,sep=" "),colnames(mm$x[[par_name]]),sep=".")],temp_par_s_unlisted)
        }

        backfitting_iteration <- function(
          par_s,
          par_cov,
          beta_start,
          lambda_s,
          first_inner_run,
          K,
          margin_dist,
          copula_dist,
          dataset,
          mm,
          copula_link,
          df_s,
          step_size,
          par_name
        ) {
          
          ############# Backfitting with penalisation
          pen_mat=matrix(0,nrow=ncol(X),ncol=ncol(X))
          if(length(par_s[[par_name]])>0) {
            start_idx=ncol(mm$x[[par_name]])+1
            for (s_name in names(mm$s[[par_name]])) {
              #print(lambda_s[[par_name]])
              B=mm$s[[par_name]][[s_name]]
              n_B=ncol(B)
              idx=start_idx:(start_idx+n_B-1)
              D=diff(diag(n_B),differences=2)
              S=t(D)%*%D
              G=par_s[[par_name]][[s_name]]
              #pen_val=lambda%*% t(G)%*%S%*%G
              if(first_inner_run==TRUE) {
                pen_mat[idx,idx]=S*0
              } else {
                pen_mat[idx,idx]=S*lambda_s[[par_name]][[s_name]]
              }
              df_s[[par_name]][[s_name]]=sum(diag(B%*%ginv(t(B)%*%B+pen_mat[idx,idx])%*%t(B)))
              start_idx=start_idx+n_B
            }
          }

          beta_update=as.vector(ginv(t(X)%*%W%*%X + pen_mat)%*%t(X)%*%W%*%z_k)
          beta_change_inner=beta_update-beta_start
          beta_new=beta_start*(1-step_size) + (step_size)*(beta_update)
          beta_new

          temp_par_cov_new=beta_new[grepl(par_name,names(beta_new))]
          par_cov_new=c(beta_new[names(par_cov)[grepl(par_name,names(par_cov))]],par_cov[!names(par_cov) %in% names(temp_par_cov_new)])
          #par_cov_new[names(beta)]=beta
          temp_par_s_new=beta_new[!names(beta_new) %in% names(par_cov_new)]
          #Select all beta_new which have names corresponding to s_name
          par_s_new=par_s
          for(s_name in names(par_s[[par_name]])) {
            par_s_new[[par_name]][[s_name]]=temp_par_s_new[grepl(s_name,names(temp_par_s_new),fixed=TRUE)]
          }

          eta_out=calc_eta(par_cov_new,mm,margin_dist,copula_link,par_s=par_s_new)

          eta=eta_out$eta; eta_dr=eta_out$eta_dr; eta_inv=eta_out$eta_inv
          par_cov=par_cov_new
          par_s=par_s_new

          calc_lik_out_end=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=FALSE
            ,response=dataset$response,response_margin=(dataset$time),response_subject = dataset$subject)      
          

          #print(sum(unlist(df_s[[par_name]])))
          GAIC_lambda_k=-2*calc_lik_out_end$log_lik["joint"]+ K*sum(unlist(df_s[[par_name]]))

          #print(paste("K*DF_S",K*sum(unlist(df_s[[par_name]]))))
          
          return_list=list(par_cov,par_s,calc_lik_out_end,GAIC_lambda_k,df_s)
          names(return_list)=c("par_cov","par_s","calc_lik_out_end","GAIC_lambda_k","df_s")
          return(return_list)
        }
        
        optim_lambda <- function(lambda_val,smooth_name,
        par_s,par_cov, beta_start, lambda_s, first_inner_run=FALSE,K=K,
                margin_dist, copula_dist, dataset, mm, copula_link,df_s,step_size,par_name) {
          lambda_s_temp=lambda_s
          lambda_s_temp[[par_name]][[smooth_name]]=lambda_val
          backfitting_iteration_results=backfitting_iteration(par_s=par_s,par_cov=par_cov, beta_start=beta_start, lambda_s=lambda_s_temp, first_inner_run=FALSE,K=K,
                margin_dist=margin_dist, copula_dist=copula_dist, dataset=dataset, mm=mm, copula_link=copula_link
                ,df_s=df_s,step_size=step_size,par_name=par_name)
          
          # Debug output
          loglik <- backfitting_iteration_results$calc_lik_out_end$log_lik["joint"]
          df_total <- sum(unlist(backfitting_iteration_results$df_s[[par_name]]))
          gaic_val <- backfitting_iteration_results$GAIC_lambda_k
          #print(sprintf("λ=%.3f | LogLik=%.2f | DF=%.2f | GAIC=%.2f\n", 
          #           lambda_val, loglik, df_total, gaic_val))

          return(backfitting_iteration_results$GAIC_lambda_k)
        }

        K=lambda_penalty_K
        num_smooths=length(lambda_s[[par_name]])
        if(num_smooths==0|outer_only_run_counter==1) {
          backfitting_iteration_results=backfitting_iteration(par_s=par_s,par_cov=par_cov, beta_start=beta_start, lambda_s=lambda_s, first_inner_run=TRUE,K=K,
                margin_dist=margin_dist, copula_dist=copula_dist, dataset=dataset, mm=mm, copula_link=copula_link
                ,df_s=df_s,step_size=step_size,par_name=par_name)
        } else {
          for (smooth_name in names(lambda_s[[par_name]])) {
           #Optimize lambda for each smooth
           if(inner_run_counter==1) {
             cat(paste("\nOptimising smoothing parameter for",par_name,"-",smooth_name))
              optim_lambda_out=optim(par=lambda_s[[par_name]][[smooth_name]],fn=optim_lambda,
                smooth_name=smooth_name,
                par_s=par_s,par_cov=par_cov, beta_start=beta_start, lambda_s=lambda_s, first_inner_run=FALSE,K=K,
                  margin_dist=margin_dist, copula_dist=copula_dist, dataset=dataset, mm=mm, copula_link=copula_link
                  ,df_s=df_s,step_size=step_size,par_name=par_name,
                  method="L-BFGS-B",lower=1,upper=1e3,control = list(factr=1,pgtol=.1)
              )
              lambda_s[[par_name]][[smooth_name]]=optim_lambda_out$par
              if(verbose>2) {
                print(paste("Chosen lambda:" ,round(lambda_s[[par_name]][[smooth_name]],2), "| Penalty K =", K))
              }
            } #end if inner_run_counter
          } #end for smooth_name
        } #end if num_smooths
        
        backfitting_iteration_results=backfitting_iteration(par_s=par_s,par_cov=par_cov, beta_start=beta_start, lambda_s=lambda_s, first_inner_run=FALSE,K=K,
          margin_dist=margin_dist, copula_dist=copula_dist, dataset=dataset, mm=mm, copula_link=copula_link,df_s=df_s,step_size=step_size,par_name=par_name)
        par_cov=backfitting_iteration_results$par_cov
        par_s=backfitting_iteration_results$par_s
        calc_lik_out_end=backfitting_iteration_results$calc_lik_out_end
        df_s=backfitting_iteration_results$df_s

        if (verbose>2) {
          cat("\nLogLik:\n")
          print(calc_lik_out_end$log_lik)
        }

        if(plot_results==TRUE) {
          plot_count=3+length(par_cov)
          sides=round(sqrt(plot_count))

          par(mfrow=c(sides+1,sides))
          plot(log_lik_history[,3],type="l",main="LogLik - Overall")
          plot(log_lik_history[,1],type="l",main="LogLik - Margin")
          plot(log_lik_history[,2],type="l",main="LogLik - Copula")

          for(i in 1:length(colnames(par_history))) {


            if(!all(is.na(true_val))) {
              plot(par_history[,i],type="l",main=colnames(par_history)[i],xlab="Iteration",ylab="Parameter estimate",ylim=range(c(par_history[,i],true_val[i])))
              abline(h=true_val[i],col="red")
            } else {
              plot(par_history[,i],type="l",main=colnames(par_history)[i],xlab="Iteration",ylab="Parameter estimate")
            }
          }
        }

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Plotting")
        #print(timer)

        change_log_lik=calc_lik_out_end$log_lik["joint"]-calc_lik_out$log_lik["joint"]

        run_counter=run_counter+1
        outer_run_counter=outer_run_counter+1
        inner_run_counter=inner_run_counter+1

      }
      weights_final[[par_name]]=score$w_k
    }

    step_size = (step_adjustment^min(outer_only_run_counter,max_steps))*start_step_size
    outer_only_run_counter=outer_only_run_counter+1
    outer_end_log_lik=calc_lik_out_end$log_lik["joint"]
    outer_log_lik_change=outer_end_log_lik-outer_start_log_lik

    out_temp=c(outer_start_log_lik,outer_end_log_lik,outer_log_lik_change)
    names(out_temp) = c("Start LogLik","End LogLik","Change")
    cat("\n")
    print(out_temp)


    if(abs(outer_log_lik_change)<=outer_stop_crit) {
      print(c(outer_end_log_lik-outer_start_log_lik))
      cat("\nOUTER CONVERGED")
    }
    
  }

  cat("\n\n############ MODEL FIT ############\n")
  cat(paste("\nMargin distribution:",margin_dist$family[2]))
  cat(paste("\nCopula distribution:",copula_dist))
  cat("\n")
  cat(paste("\nParameter count:",length(par_cov)))
  cat(paste("\nObservations:",nrow(dataset)))
  cat(paste("\nMargins:",length(unique(dataset$time))))
  cat("\n")
  total_fit_time <- as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
  cat(paste("\nTotal time (seconds):",round(total_fit_time,2)))
  cat("\n\n")
  par_mat_out_temp=t(t((par_cov)))
  colnames(par_mat_out_temp) = c("estimate")
  print(par_mat_out_temp)
  cat("\n")
  cat("Model Selection Criteria:")
  cat("\n")

  p_cop=par_cov[grepl("theta",names(par_cov))|grepl("zeta",names(par_cov))]
  p_mar=par_cov[!(grepl("theta",names(par_cov))|grepl("zeta",names(par_cov)))]

  df_s_total=df_s_cop_total=df_s_margin_total=0
  for(par_name in names(par_s)) {
    if(par_name %in% c("theta","zeta")) 
      df_s_cop_total=df_s_cop_total+sum(unlist(df_s[[par_name]]))
    else {
      df_s_margin_total=df_s_margin_total+sum(unlist(df_s[[par_name]]))
    }
    df_s_total=df_s_total+sum(unlist(df_s[[par_name]]))
  }

  aics=rbind(t(calc_lik_out_end$log_lik),
             t(-calc_lik_out_end$log_lik*2)+2*c(length(p_mar)+df_s_margin_total,length(p_cop)+df_s_cop_total,length(par_cov)+df_s_total),
             t(-calc_lik_out_end$log_lik*2)+c(length(p_mar)+df_s_margin_total,length(p_cop)+df_s_cop_total,length(par_cov)+df_s_total)*log(nrow(dataset)),
             t(c(length(p_mar)+df_s_margin_total,length(p_cop)+df_s_cop_total,length(par_cov)+df_s_total))
             )

  rownames(aics)=c("LogLik","AIC","BIC","EDF")
  print(aics)

  cat("\n####################################\n")

  return_list=list(par_cov,log_lik_history,par_history,calc_lik_out_end,mm,margin_dist,copula_dist,include_dlcopdpar,dataset$response,dataset$time,dataset$subject,par_s,lambda_s,df_s,weights_final)
  names(return_list)=c("par","log_lik_history","par_history","calc_lik_out_end","model_matrix","margin_dist","copula_dist","include_dlcopdpar","response","response_margin","response_subject","par_s","lambda_s","df_s","weights")
  class(return_list)="gamlss.longitudinal"
  return(return_list)
}

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
#' @param copula.family Copula distribution specified as a string in VineCopula style
#' e.g. "t", "C" etc.
#' @param copula.link List of link functions for the copula parameters
#' @return Returns a list mm with items mm$x and mm$s for fixed and smooth terms respectively,
#' with each of those lists being lists of each parameter and their respective model matrices
#' 
#' @export
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
    quiet_gamlss2=TRUE
) {

  dataset_mm <- dataset
  if ("response" %in% names(dataset_mm) && any(is.na(dataset_mm$response))) {
    obs_resp <- dataset_mm$response[!is.na(dataset_mm$response)]
    fill_val <- if (length(obs_resp) > 0) mean(obs_resp) else 0
    dataset_mm$response[is.na(dataset_mm$response)] <- fill_val
  }

  run_gamlss2 <- function(...) {
    if (isTRUE(quiet_gamlss2)) {
      fit <- NULL
      invisible(utils::capture.output({
        fit <- suppressMessages(suppressWarnings(gamlss2(...)))
      }, type = "output"))
      return(fit)
    }
    gamlss2(...)
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

  if(copula.family %in% c("t")){two_par_cop=TRUE} else {two_par_cop=FALSE}
  included_parameters <- c(names(margin.family$parameters), if(two_par_cop) c("theta","zeta") else c("theta"))

  formulas=list()
  for (parameter in included_parameters) {
    formulas[[parameter]]=get(paste(parameter,"formula",sep="."))
  }

  m_temp=list()
  m_temp[["mu"]] <- run_gamlss2(
    formula = as.formula(mu.formula),
    family = NO(),
    data = dataset_mm,
    control = gamlss2_control(maxit = 1)
  )

  for (parameter in included_parameters[2:length(included_parameters)]) {
    formulas[[parameter]] <- to_response_formula(formulas[[parameter]], response_name = "response")
    m_temp[[parameter]] <- run_gamlss2(
      formula = formulas[[parameter]],
      family = NO(), # margin.family
      data = if(parameter %in% c("theta","zeta")) {
        dataset_mm[dataset_mm$time %in% unique(dataset_mm$time)[1:(length(unique(dataset_mm$time))-1)], ]
      } else dataset_mm,
      control = gamlss2_control(maxit = 1)
    )
  }

  # print(m_temp)

  mm_x=list()
  mm_s=list()
  for(parameter in included_parameters) {
    m=m_temp[[parameter]]
    formula_terms = stats::terms(formulas[[parameter]])
    has_intercept = as.integer(attr(formula_terms, "intercept")) == 1L

    # Use formula term labels (not xterms) so factor bases like "gender"
    # are retained even when xterms contains expanded names like "gender1".
    fixed_terms = attr(formula_terms, "term.labels")
    fixed_terms = fixed_terms[!grepl("^\\s*s\\(", fixed_terms)]
    fixed_terms = intersect(fixed_terms, names(m$model))

    if(length(fixed_terms) > 0 || has_intercept) {
      fixed_formula = stats::reformulate(termlabels = fixed_terms, intercept = has_intercept)
      X_fixed = stats::model.matrix(fixed_formula, data = m$model)
      colnames(X_fixed) = sub("^\\(Intercept\\)$", "intercept", colnames(X_fixed))
      mm_x[[parameter]] = as.data.frame(X_fixed, check.names = FALSE)
    } else {
      mm_x[[parameter]] = data.frame(row.names = seq_len(nrow(m$model)))
    }

    if(length(m$sterms$mu)==0) {
      mm_s[[parameter]]=NULL
    } else {
      mm_s[[parameter]]=list()
      for (s in m$sterms$mu) {
        mm_s[[parameter]][[s]]=m$specials[[s]]$X
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
#' @export
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

#' Calculate the likelihood components for the joint model
#' 
#' This function calculates the marginal and copula log likelihoods and components
#' for the joint model based on the inverse link function values (eta_inv),
#' model matrices (mm), marginal distribution functions (margin_dist),
#' and copula distribution (copula_dist). It computes the marginal density,
#' copula density, and joint log-likelihood.
#' 
#' @param eta_inv A list of inverse link function values for each parameter.
#' @param mm A list containing model matrices for fixed effects (mm$x).
#' @param margin_dist A list of functions for the marginal distribution, including density and distribution functions.
#' @param copula_dist A string specifying the copula distribution (e.g., "C", "t").
#' @param calc_d2 A logical indicating whether to calculate second derivatives (default is FALSE).
#' @param response A numeric vector of response values.
#' @param response_margin A numeric vector indicating the margin (time) for each response.
#' @param response_subject A numeric vector indicating the subject for each response.
#' @return A list containing:
#' \item{log_lik}{A named vector with marginal, copula, and joint log-likelihoods.}
#' \item{margin_d}{A numeric vector of marginal densities.}
#' \item{copula_d}{A numeric vector of copula densities.}
#' \item{margin_p}{A numeric vector of marginal distribution function values.}
#' \item{Fx_1_2}{A matrix of marginal distribution function values for pairs of margins.}
#' \item{order_copula}{A matrix indicating the order of margins and subjects for copula calculations.}
#' \item{margin_deriv}{A list of marginal derivatives.}
#' \item{order_copula}{A matrix indicating the order of margins and subjects for copula calculations.}
#' 
#' @export
calc_likelihood_minimal <- function(eta_inv,mm,margin_dist,copula_dist,calc_d2=FALSE,response,response_margin,response_subject,penalize_smooth=FALSE,par_s=NA) {
  #Setup input matrix of response and parameters
  #response=dataset$response; response_subject=dataset$subject; response_margin=dataset$time; dataset=NA
  margin_names=sort(unique(response_margin))
  num_margins=length(margin_names)
  n_obs=length(response)

  obs_response=!is.na(response)

  order_margin=cbind(response_margin,response_subject)
  colnames(order_margin)=c("time","subject")

  margin_deriv_input=list()
  margin_deriv_input[["y"]]=response
  margin_deriv_input[["q"]]=response
  margin_deriv_input[["x"]]=response
  for (par_name in names(mm)) {
    if (par_name %in% c("mu","sigma","nu","tau")) {
      margin_deriv_input[[par_name]]=eta_inv[[par_name]]
    }
  }

  #Calculate all derivatives

  ################## MARGIN DERIVATIVES
  if(calc_d2==TRUE) {
    to_include=grepl("dld",names(margin_dist))|grepl("d2ld",names(margin_dist))
  } else {
    to_include=grepl("dld",names(margin_dist))
  }

  margin_deriv_names=names(margin_dist)[to_include]
  margin_deriv=list()
  for (deriv_name in margin_deriv_names) {
    FUN=margin_dist[[deriv_name]]
    FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
    deriv_val=do.call(FUN,args=margin_deriv_input[FUN_args])
    if(length(deriv_val)==n_obs) {
      deriv_val[!obs_response]=0
      deriv_val[!is.finite(deriv_val)]=0
    }
    margin_deriv[[deriv_name]]=deriv_val
  }

  margin_pFUN=eval(parse( text=paste("p",margin_dist$family[1],sep="") ))
  FUN=margin_pFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_p=do.call(FUN,args=margin_deriv_input[FUN_args])
  margin_p[!obs_response]=NA
  margin_p[!is.finite(margin_p)]=NA

  margin_dFUN=eval(parse( text=paste("d",margin_dist$family[1],sep="") ))
  FUN=margin_dFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_d=do.call(FUN,args=margin_deriv_input[FUN_args])
  margin_d[!obs_response]=NA
  margin_d[!is.finite(margin_d) | margin_d<=0]=NA

  ################COPULA DERIVATIVES
  #First calculate margin F(x1), F(x2) as inputs to copula

  pair_df_all=list()
  base_df=data.frame(
    row_id=seq_len(n_obs),
    time=response_margin,
    subject=response_subject,
    u=margin_p,
    observed=obs_response,
    stringsAsFactors = FALSE
  )

  if(num_margins>1) {
    for (i in seq_len(num_margins-1)) {
      t1=margin_names[i]
      t2=margin_names[i+1]

      left=base_df[base_df$time==t1,c("row_id","subject","time","u","observed")]
      right=base_df[base_df$time==t2,c("row_id","subject","time","u","observed")]
      names(left)=c("row_id1","subject","time1","u1","observed1")
      names(right)=c("row_id2","subject","time2","u2","observed2")

      pair_i=merge(left,right,by="subject",all=FALSE)
      if(nrow(pair_i)>0) {
        pair_df_all[[length(pair_df_all)+1]]=pair_i
      }
    }
  }

  if(length(pair_df_all)==0) {
    pair_df=data.frame(
      subject=response_subject[0],
      row_id1=integer(0),
      time1=response_margin[0],
      u1=numeric(0),
      observed1=logical(0),
      row_id2=integer(0),
      time2=response_margin[0],
      u2=numeric(0),
      observed2=logical(0)
    )
  } else {
    pair_df=do.call(rbind,pair_df_all)
  }

  order_copula=as.matrix(pair_df[,c("time1","subject","time2","subject")])
  colnames(order_copula)=c("time1","subject1","time2","subject2")

  Fx_1_2=as.matrix(pair_df[,c("u1","u2")])
  if(nrow(Fx_1_2)==0) {
    Fx_1_2=matrix(numeric(0),ncol=2)
  }
  colnames(Fx_1_2)=c("u1","u2")

  pair_complete=rep(FALSE,nrow(pair_df))
  if(nrow(pair_df)>0) {
    pair_complete=pair_df$observed1 & pair_df$observed2 & is.finite(pair_df$u1) & is.finite(pair_df$u2)
  }

  par1=rep(NA_real_,nrow(pair_df))
  par2=rep(NA_real_,nrow(pair_df))
  if(nrow(pair_df)>0) {
    theta_rows = which(response_margin %in% margin_names[seq_len(max(1, num_margins-1))])
    theta_index_map=rep(NA_integer_,n_obs)
    theta_index_map[theta_rows]=seq_along(theta_rows)
    theta_idx=theta_index_map[pair_df$row_id1]

    par1=eta_inv[["theta"]][theta_idx]
    if("zeta" %in% names(eta_inv)) {
      par2=eta_inv[["zeta"]][theta_idx]
    } else {
      par2=rep(0,length(par1))
    }
  }

  pair_complete=pair_complete & is.finite(par1) & is.finite(par2)

  Fx_eval=Fx_1_2
  if(nrow(Fx_eval)>0) {
    Fx_eval[!is.finite(Fx_eval)]=0.5
    Fx_eval[Fx_eval>1]=1
    Fx_eval[Fx_eval<0]=0
  }

  par1_eval=par1
  par2_eval=par2
  par1_eval[!is.finite(par1_eval)]=0
  par2_eval[!is.finite(par2_eval)]=0

  if(copula_dist=="C") {
    par1_eval[par1_eval>=28]=27.9
  }

  if(length(par1_eval)==0) {
    copula_d=numeric(0)
  } else {
    copula_d=BiCopPDF(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval)
  }
  if(length(copula_d)>0) {
    copula_d[!is.finite(copula_d) | copula_d<=0]=1
    copula_d[!pair_complete]=1
  }

  ########COMBINE MARGINS AND COPULA DERVIATIVES

  margin_loglik_terms=log(margin_d[!is.na(margin_d)])
  margin_loglik_terms=margin_loglik_terms[is.finite(margin_loglik_terms)]
  copula_loglik_terms=log(copula_d[pair_complete])
  copula_loglik_terms=copula_loglik_terms[is.finite(copula_loglik_terms)]

  log_lik=c(sum(margin_loglik_terms),sum(copula_loglik_terms),sum(margin_loglik_terms)+sum(copula_loglik_terms))
  names(log_lik)=c("marginal","copula","joint")

  copula_p=rep(NA_real_,length(copula_d))

  return_list=list(log_lik,margin_d,copula_d,margin_p,copula_p,Fx_1_2,order_copula,margin_deriv,pair_complete,par1,par2)
  names(return_list)=c("log_lik","margin_d","copula_d","margin_p","copula_p","Fx_1_2","order_copula","margin_deriv","pair_complete","copula_par1","copula_par2")
  return(return_list)
}

#' 
#' 
#' @export
score_function_v2 <- function(eta,dldpar,d2ldpar,dpardeta,response=NA,phi=1,step_size=1,verbose=FALSE,crit_wk=0.0000001) {

  u_k=dldeta = dldpar * dpardeta
  f_k=d2ldpar
  w_k=-f_k*(dpardeta*dpardeta)

  #Stop if weights are too small
  w_k[abs(w_k)<crit_wk]=1
  u_k[abs(w_k)<crit_wk]=0

  w_k[abs(u_k)<crit_wk]=1
  u_k[abs(u_k)<crit_wk]=0

  z_k=(1-phi)*eta+phi*(eta+step_size*(u_k/w_k))

  if(verbose==TRUE) {
    steps_mean=round(rbind(colMeans(as.matrix(eta))
                           ,colMeans(as.matrix(dldpar-dlcopdpar))
                           ,colMeans(as.matrix(dlcopdpar))
                           ,colMeans(as.matrix(dpardeta))
                           ,colMeans(as.matrix(dpardeta*dpardeta))
                           ,colMeans(as.matrix(f_k))
                           ,colMeans(as.matrix(w_k))
                           ,colMeans(as.matrix(u_k))
                           ,colMeans(as.matrix(u_k/w_k))
                           ,colMeans(as.matrix(z_k))
    ),8)
    rownames(steps_mean)=c("eta","dldpar","dlcopdpar","dpardeta","dpardeta2","f_k","w_k","u_k","(1/w_k)*u_k","z_k")
    print(steps_mean)
  }
  return_list=list(colMeans(as.matrix(z_k)),as.matrix(u_k),as.matrix(f_k),as.matrix(w_k),as.matrix(z_k))
  names(return_list)=c("par","u_k","f_k","w_k","z_k")
  return(return_list)
}

calc_deriv_copula_wrt_margin = function(input,margin_par,par_name,calc_d2=FALSE) {

    #Calculate copula derivative with respect to marginal parameters
    #input=copula_merged

    if(calc_d2==FALSE) {

      dlcopdpar=matrix(0,nrow=nrow(input),ncol=length(margin_par))
      i=1
      for (inner_par_name in margin_par) {

        if(inner_par_name==par_name) {
          #Take parameters from input for clarity
          dc_tplus_du_t=input[,"dcdu1"]
          dc_tplus_du_tplus=input[,"dcdu2"]
          l_t=input[,paste(paste("dld",inner_par_name,sep=""),".x",sep="")]
          l_t_plus=input[,paste(paste("dld",inner_par_name,sep=""),".y",sep="")]
          x_t=input[,"response.x"]
          x_t_plus=input[,"response.y"]
          f_t=input[,"margin_d.x"]
          f_t_plus=input[,"margin_d.y"]
          c_tplus=input[,"copula_d"]
          mu_t=input[,"mu.x"]
          mu_t_plus=input[,"mu.y"]

          F_nd_t=input[,"F_nd.x"]
          F_nd_t_plus=input[,"F_nd.y"]

          du_t_dmu=F_nd_t
          du_t_plus_dmu=F_nd_t_plus

          dc_plus_dt_dmu=dc_tplus_du_t * du_t_dmu
          dc_plus_dt_plus_dmu=dc_tplus_du_tplus * du_t_plus_dmu
          dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)]=0
          dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)]=0
          dcdmu_tplus=((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
          dcdmu_tplus[is.nan(dcdmu_tplus)|is.na(dcdmu_tplus)]=0

          dlcopdpar[,i]=dcdmu_tplus

        }
        i=i+1
      }
      colnames(dlcopdpar)=paste("dlcopd",margin_par,sep="")

      par_dlcopdpar=dlcopdpar[,paste("dlcopd",margin_par,sep="")]
      merged_dlcopdpar=merge(cbind(input[,c("time1","time2","subject1","subject2")],par_dlcopdpar),cbind(input[,c("time1","time2","subject1","subject2")],par_dlcopdpar),by.x=c("time2","subject2"),by.y=c("time1","subject1"),all=TRUE)
      merged_dlcopdpar[is.na(merged_dlcopdpar)]=0

      x_comp=grepl("dlcopd",colnames(merged_dlcopdpar))&grepl(".x",colnames(merged_dlcopdpar))
      y_comp=grepl("dlcopd",colnames(merged_dlcopdpar))&grepl(".y",colnames(merged_dlcopdpar))

      d1_cop=0.5*(merged_dlcopdpar[,x_comp]+merged_dlcopdpar[,y_comp])

      return(d1_cop)

    } else {

      d2lcopdpar2=matrix(0,nrow=nrow(input),ncol=length(margin_par))
      i=1
      for (inner_par_name in margin_par) {

        if(inner_par_name==par_name) {
          #Take parameters from input for clarity
          dc_tplus_du_t=input[,"dcdu1"]
          dc_tplus_du_tplus=input[,"dcdu2"]
          #l_t=input[,paste(paste("dld",inner_par_name,sep=""),".x",sep="")]
          #l_t_plus=input[,paste(paste("dld",inner_par_name,sep=""),".y",sep="")]
          #x_t=input[,"response.x"]
          #x_t_plus=input[,"response.y"]
          #f_t=input[,"margin_d.x"]
          #f_t_plus=input[,"margin_d.y"]
          c_tplus=input[,"copula_d"]
          mu_t=input[,"mu.x"]
          mu_t_plus=input[,"mu.y"]

          F_nd_t=input[,"F_nd.x"]
          F_nd_t_plus=input[,"F_nd.y"]

          du_t_dmu=F_nd_t
          du_t_plus_dmu=F_nd_t_plus

          dc_plus_dt_dmu=dc_tplus_du_t * du_t_dmu
          dc_plus_dt_plus_dmu=dc_tplus_du_tplus * du_t_plus_dmu
          dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)]=0
          dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)]=0
          dcdmu_tplus=((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
          dcdmu_tplus[is.nan(dcdmu_tplus)|is.na(dcdmu_tplus)]=0

          #dlcopdpar[,i]=dcdmu_tplus

          #######NOW FOR SECOND DERIVATIVE OF COPULA TERM

          F_nd2=input[,"F_nd2.x"]
          F_nd2_plus=input[,"F_nd2.y"]

          d2u_t_dmu2=F_nd2
          d2u_t_plus_dmu2=F_nd2_plus

          d2cdu_t2=input[,"d2cdu12"]
          d2cdu_t_plus2=input[,"d2cdu22"]
          d2cdu_t2[is.nan(d2cdu_t2)]=0
          d2cdu_t_plus2[is.nan(d2cdu_t_plus2)]=0

          d2cdmu2=  d2cdu_t2*du_t_dmu^2 +
                    dc_tplus_du_t * d2u_t_dmu2 +
                    d2cdu_t_plus2*du_t_plus_dmu^2 +
                    dc_tplus_du_tplus * d2u_t_plus_dmu2

          d2lcdmu2=as.matrix((d2cdmu2*c_tplus-(dcdmu_tplus^2))/(c_tplus^2))
          d2lcdmu2=input[,"c_nd2"]

          d2lcopdpar2[,i]=d2lcdmu2

        }
        i=i+1
      }
      colnames(d2lcopdpar2)=paste("d2lcopd",margin_par,sep="")

      par_d2lcopdpar=d2lcopdpar2[,paste("d2lcopd",margin_par,sep="")]
      merged_d2lcopdpar=merge(cbind(input[,c("time1","time2","subject1","subject2")],par_d2lcopdpar)
                              ,cbind(input[,c("time1","time2","subject1","subject2")],par_d2lcopdpar)
                              ,by.x=c("time2","subject2"),by.y=c("time1","subject1"),all=TRUE)
      merged_d2lcopdpar[is.na(merged_d2lcopdpar)]=0

      x_comp=grepl("d2lcopd",colnames(merged_d2lcopdpar))&grepl(".x",colnames(merged_d2lcopdpar))
      y_comp=grepl("d2lcopd",colnames(merged_d2lcopdpar))&grepl(".y",colnames(merged_d2lcopdpar))

      d2_cop=0.5*(merged_d2lcopdpar[,x_comp]+merged_d2lcopdpar[,y_comp])

      #plot(d2lcopdpar2[,paste("d2lcopd",par_name,sep="")],input[,"c_nd2"],main="d2",ylab="numerical")

      return(d2_cop)
    }


}

calc_copula_derivatives = function(eta_inv, Fx_1_2, copula_dist, calc_d2=FALSE, calc_d2_marginal=FALSE, par1=NULL, par2=NULL, pair_complete=NULL) {

  if(is.null(par1)) {
    par1=eta_inv[["theta"]]
  }

  if(is.null(par2)) {
    if("zeta" %in% names(eta_inv)) {
      par2=eta_inv[["zeta"]]
    } else {
      par2=eta_inv[["theta"]]*0
    }
  }

  if(is.null(pair_complete)) {
    pair_complete=rep(TRUE,length(par1))
  }

  if(length(par1)==0) {
    if("zeta" %in% names(eta_inv)) {
      if(calc_d2==TRUE) {
        return(list(dldth=numeric(0),dcdth=numeric(0),dldz=numeric(0),dcdz=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0),d2ldth2=numeric(0),d2ldz2=numeric(0),d2ldthdz=numeric(0),d2cdu12=numeric(0),d2cdu22=numeric(0)))
      }
      return(list(dldth=numeric(0),dcdth=numeric(0),dldz=numeric(0),dcdz=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0)))
    }
    if(calc_d2==TRUE) {
      return(list(dldth=numeric(0),dcdth=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0),d2ldth2=numeric(0),d2cdu12=numeric(0),d2cdu22=numeric(0)))
    }
    return(list(dldth=numeric(0),dcdth=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0)))
  }

  Fx_eval=as.matrix(Fx_1_2)
  Fx_eval[!is.finite(Fx_eval)]=0.5
  Fx_eval[Fx_eval>1]=1
  Fx_eval[Fx_eval<0]=0

  par1_eval=par1
  par2_eval=par2
  par1_eval[!is.finite(par1_eval)]=0
  par2_eval[!is.finite(par2_eval)]=0

  if(copula_dist=="C") {
    par1_eval[par1_eval>=28]=27.9
  }

  copula_d=BiCopPDF(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval)
  copula_d[!is.finite(copula_d) | copula_d<=0]=1
  copula_d[!pair_complete]=1

  dldth=BiCopDeriv(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="par",log=TRUE)
  dcdth=BiCopDeriv(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="par",log=FALSE)

  if(calc_d2==TRUE) {
    d2cdth=BiCopDeriv2(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="par")
    d2ldth2=(1/(copula_d^2))*(copula_d*d2cdth-dcdth^2)
  }

  if("zeta" %in% names(eta_inv)) {
    dldz=BiCopDeriv(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="par2",log=TRUE)
    dcdz=BiCopDeriv(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="par2",log=FALSE)

    if(calc_d2==TRUE) {
      d2cdz=BiCopDeriv2(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="par2")
      d2ldz2=(1/(copula_d^2))*(copula_d*d2cdz-dcdz^2)

      d2cdthdz=BiCopDeriv2(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="par1par2")
      d2ldthdz=(d2cdthdz*copula_d-dcdth*dcdz)/(copula_d^2)
    }

  }
  dcdu1=BiCopDeriv(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="u1",log=FALSE)
  dcdu2=BiCopDeriv(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="u2",log=FALSE)

  dldth[!is.finite(dldth)]=0; if(calc_d2==TRUE) {d2ldth2[!is.finite(d2ldth2)]=0  }
  dldth[!pair_complete]=0
  dcdth[!pair_complete]=0
  dcdu1[!pair_complete]=0
  dcdu2[!pair_complete]=0

  if(calc_d2==TRUE) {
    d2cdu12=BiCopDeriv2(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="u1")
    d2cdu22=BiCopDeriv2(Fx_eval[,1],Fx_eval[,2],family = as.numeric(BiCopName(copula_dist)),par=par1_eval,par2=par2_eval,deriv="u2")
    d2cdu12[!is.finite(d2cdu12)]=0
    d2cdu22[!is.finite(d2cdu22)]=0
    d2cdu12[!pair_complete]=0
    d2cdu22[!pair_complete]=0
  }

  if("zeta" %in% names(eta_inv)) {
    dldz[!is.finite(dldz)]=0
    dcdz[!is.finite(dcdz)]=0
    dldz[!pair_complete]=0
    dcdz[!pair_complete]=0
    if(calc_d2==TRUE) {
      d2ldz2[!is.finite(d2ldz2)]=0
      d2ldthdz[!is.finite(d2ldthdz)]=0
      d2ldz2[!pair_complete]=0
      d2ldthdz[!pair_complete]=0
    }
  }

  if(calc_d2==TRUE) {
    d2ldth2[!pair_complete]=0
  }

  ############# RETURN LIST

  if("zeta" %in% names(eta_inv)) {
    if(calc_d2==TRUE) {
      return_list=list(dldth,dcdth,dldz,dcdz,dcdu1,dcdu2,d2ldth2,d2ldz2,d2ldthdz, d2cdu12, d2cdu22)
      names(return_list)=c("dldth","dcdth","dldz","dcdz","dcdu1","dcdu2","d2ldth2","d2ldz2","d2ldthdz","d2cdu12","d2cdu22")
    } else {
      return_list=list(dldth,dcdth,dldz,dcdz,dcdu1,dcdu2)
      names(return_list)=c("dldth","dcdth","dldz","dcdz","dcdu1","dcdu2")
    }
  } else {
    if(calc_d2==TRUE) {
      return_list=list(dldth,dcdth,dcdu1,dcdu2,d2ldth2, d2cdu12, d2cdu22)
      names(return_list)=c("dldth","dcdth","dcdu1","dcdu2","d2ldth2","d2cdu12","d2cdu22")
    } else {
      return_list=list(dldth,dcdth,dcdu1,dcdu2)
      names(return_list)=c("dldth","dcdth","dcdu1","dcdu2")
    }
  }

  return(return_list)
}

calc_Fx_derivatives = function(eta_inv, mm, margin_dist,response) {
  # Allow callers to pass full model matrix object; we only need fixed-effect blocks.
  if (is.list(mm) && all(c("x", "s") %in% names(mm))) {
    mm = mm$x
  }

  nd_impact=nd_impact_m=nd_impact_c=rep(0,length(names(eta_inv)))
  nd_impact_F=list() ###### WE DON"T NEED TO BE CALCULATING THIS FOR ALL PARAMETERS EVERY TIME

  margin_par_names=names(eta_inv)[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (eta_par_names_nd in margin_par_names) {
    adj_fac=.0001
    change=change_m=change_c=c(0,0)
    change_F=matrix(0,nrow=length(eta_inv[[1]]),ncol=2)
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      #change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["joint"]
      #change_m[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["marginal"]
      #change_c[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["copula"]
      change_F[,i]=calc_F_x(eta_inv_adj,mm,margin_dist,response)#calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
      #print(change_F[,i]==calc_F_x(eta_inv_adj,mm,margin_dist))
      i=i+1
    }
    #nd_impact[eta_par_names_nd]=(change[2]-change[1])/(2*adj_fac)
    #nd_impact_m[eta_par_names_nd]=(change_m[2]-change_m[1])/(2*adj_fac)
    #nd_impact_c[eta_par_names_nd]=(change_c[2]-change_c[1])/(2*adj_fac)
    nd_impact_F[[eta_par_names_nd]]=(change_F[,2]-change_F[,1])/(2*adj_fac)
  }
  return(nd_impact_F)
}

calc_Fx2_derivatives = function(eta_inv, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA) {
  nd_impact=nd_impact_m=nd_impact_c=rep(0,length(names(eta_inv)))
  nd_impact_F=nd_impact_c=list() ###### WE DON"T NEED TO BE CALCULATING THIS FOR ALL PARAMETERS EVERY TIME

  margin_par_names=names(eta_inv)[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (eta_par_names_nd in margin_par_names) {
    adj_fac=.0001
    change_F=matrix(0,nrow=length(eta_inv[[1]]),ncol=3)
    change_c=matrix(0,nrow=length(eta_inv[["theta"]]),ncol=3)
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      #change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["joint"]
      #change_m[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["marginal"]
      if(testing==TRUE) {
        change_c[,i]=log(calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2=FALSE,response,response_margin,response_subject)$copula_d)
      }
      change_F[,i]=calc_F_x(eta_inv_adj,mm,margin_dist,response)#calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
      #print(change_F[,i]==calc_F_x(eta_inv_adj,mm,margin_dist))
      i=i+1
    }
    change_F[,3]=calc_F_x(eta_inv,mm,margin_dist,response)
    if(testing==TRUE) {
      change_c[,3]=log(calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2=FALSE,response,response_margin,response_subject)$copula_d)
    }
    #nd_impact[eta_par_names_nd]=(change[2]-change[1])/(2*adj_fac)
    #nd_impact_m[eta_par_names_nd]=(change_m[2]-change_m[1])/(2*adj_fac)
    nd_impact_c[[eta_par_names_nd]]=(change_c[,2]+change_c[,1]-2*change_c[,3])/(adj_fac^2)
    nd_impact_F[[eta_par_names_nd]]=(change_F[,2]+change_F[,1]-2*change_F[,3])/(adj_fac^2)
  }
  if(testing==FALSE) {
    return(nd_impact_F)
  } else {
    return(list(nd_impact_F,nd_impact_c))
  }
}

calc_true_SE_numderiv_only = function(eta_inv, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA) {

  adj_fac=.001
  nd_impact=rep(0,length(names(eta_inv)))
  names(nd_impact)=margin_par_names=names(eta_inv)#[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (eta_par_names_nd in margin_par_names) {

    change=rep(0,length(names(eta_inv)))
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
      i=i+1
    }
    change[3]=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
    nd_impact[eta_par_names_nd]=(change[2]+change[1]-2*change[3])/(adj_fac^2)

    #print(c(change,nd_impact[eta_par_names_nd]))
  }

  nd_cross=matrix(0,nrow=length(names(eta_inv)),ncol=length(names(eta_inv)))
  colnames(nd_cross)=rownames(nd_cross)=names(eta_inv)
  for (name1 in margin_par_names) {
    for (name2 in margin_par_names) {
      for(adj1 in c(-1*adj_fac,adj_fac)) {
        for(adj2 in c(-1*adj_fac,adj_fac)) {
          if(name1!=name2) {
            eta_inv_adj=eta_inv
            eta_inv_adj[[name1]]=eta_inv_adj[[name1]]+adj1
            eta_inv_adj[[name2]]=eta_inv_adj[[name2]]+adj2
            change=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
            nd_cross[name1,name2]=nd_cross[name1,name2]+change*if(adj1==adj2){1} else {-1}
          }
        }
      }
    }
  }
  nd_cross=nd_cross/(4*(adj_fac^2))

  nd2=(diag(nd_impact)+nd_cross)

  return(nd2)
}

calc_true_SE_numderiv_only_rowwise = function(eta_inv, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA) {

  adj_fac=.00001
  nd_impact_m=nd_impact_c=list()

  for (eta_par_names_nd in margin_par_names) {

    change_m=change_c=list()
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj

      calc_lik_temp=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)

      change_m[[i]]=calc_lik_temp$margin_d
      change_c[[i]]=calc_lik_temp$copula_d
      i=i+1
    }
    calc_lik_temp=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)
    change_m[[3]]=calc_lik_temp$margin_d
    change_c[[3]]=calc_lik_temp$copula_d
    nd_impact_m[[eta_par_names_nd]]=(change_m[[3]]+change_m[[1]]-2*change_m[[3]])/(adj_fac^2)
    nd_impact_c[[eta_par_names_nd]]=(change_c[[3]]+change_c[[1]]-2*change_c[[3]])/(adj_fac^2)
  }
  names(nd_impact_m)=names(nd_impact_c)=margin_par_names=names(eta_inv)

  nd_cross_m
  colnames(nd_cross)=rownames(nd_cross)=names(eta_inv)
  for (name1 in margin_par_names) {
    for (name2 in margin_par_names) {
      for(adj1 in c(-1*adj_fac,adj_fac)) {
        for(adj2 in c(-1*adj_fac,adj_fac)) {
          if(name1!=name2) {
            eta_inv_adj=eta_inv
            eta_inv_adj[[name1]]=eta_inv_adj[[name1]]+adj1
            eta_inv_adj[[name2]]=eta_inv_adj[[name2]]+adj2
            change=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
            nd_cross[name1,name2]=nd_cross[name1,name2]+change*if(adj1==adj2){1} else {-1}
          }
        }
      }
    }
  }
  nd_cross=nd_cross/(4*(adj_fac^2))

  nd2=(diag(nd_impact)+nd_cross)

  return(nd2)
}

calc_true_SE_numderiv_only_covariates = function(object, par, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA,h=.0001) {

  #object=fit; par=NA; numderiv=TRUE; sep_d2=TRUE

  response=object$response
  response_margin=object$response_margin
  response_subject=object$response_subject

  #margin_names=unique(object$response_margin)
  #num_margins=length(margin_names)

  #se_out=object$par*0;
  margin_dist=object$margin_dist; copula_dist=object$copula_dist; copula_link=get_copula_dist(copula_dist)$copula_link
  mm=object$model_matrix

  par_cov=object$par
  par_s=object$par_s

  input_par=par_cov

  adj_fac=h
  nd_impact=rep(0,length(names(input_par)))
  names(nd_impact)=names(input_par)#[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  cat("Calculating numerical first derivates for Hessian matrix...\n")
  for (par_names_nd in names(input_par)) {

    #print(par_names_nd)
    change=rep(0,length(names(input_par)))
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {

      par_cov=input_par
      par_cov[[par_names_nd]]=par_cov[[par_names_nd]]+adj

      eta_out=calc_eta(par_cov=par_cov,mm=mm,margin_dist,copula_link,par_s)
      eta_inv=eta_out$eta_inv#; eta_dr=eta_out$eta_dr; eta=eta_out$eta; eta_dr=eta_out$eta_dr
      change[i]=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=TRUE
        ,response=response,response_margin=response_margin,response_subject = response_subject)$log_lik["joint"]
      i=i+1
    }
    par_cov=input_par
    eta_out=calc_eta(par_cov=par_cov,mm=mm,margin_dist,copula_link,par_s)
    eta_inv=eta_out$eta_inv#; eta_dr=eta_out$eta_dr; eta=eta_out$eta; eta_dr=eta_out$eta_dr

        change[3]=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=TRUE
          ,response=response,response_margin=response_margin,response_subject = response_subject)$log_lik["joint"]
    nd_impact[par_names_nd]=(change[2]+change[1]-2*change[3])/(adj_fac^2)

    #print(c(change,nd_impact[eta_par_names_nd]))
  }

  cat("Calculating numerical second derivates for Hessian matrix... this may take a while\n")
  nd_cross=matrix(0,nrow=length(names(input_par)),ncol=length(names(input_par)))
  colnames(nd_cross)=rownames(nd_cross)=names(input_par)
  for (name1 in names(input_par)) {
    for (name2 in names(input_par)) {
      for(adj1 in c(-1*adj_fac,adj_fac)) {
        for(adj2 in c(-1*adj_fac,adj_fac)) {
          if(name1!=name2) {

            par=input_par
            par[[name1]]=par[[name1]]+adj1
            par[[name2]]=par[[name2]]+adj2

            eta_out=calc_eta(par_cov=par,mm=mm,margin_dist,copula_link,par_s)
            eta_inv=eta_out$eta_inv#; eta_dr=eta_out$eta_dr; eta=eta_out$eta; eta_dr=eta_out$eta_dr

            change=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=TRUE
            ,response=response,response_margin=response_margin,response_subject = response_subject)$log_lik["joint"]
            nd_cross[name1,name2]=nd_cross[name1,name2]+change*if(adj1==adj2){1} else {-1}
          }
        }
      }
    }
  }
  nd_cross=nd_cross/(4*(adj_fac^2))

  nd2=(diag(nd_impact)+nd_cross)

  return(nd2)
}

#' This function returns the log likelihood for a fitted gamlss.longitudinal object
#' @export
logLik.gamlss.longitudinal=function(object) {
  return(object$calc_lik_out$log_lik)
}

# This function returns the coefficients for a fitted gamlss.longitudinal object
#' @export
coef.gamlss.longitudinal=function(object) {
  return(object$par)
}

# This function returns the variance-covariance matrix for a given gamlss longitudinal object
#' @export
vcov.gamlss.longitudinal=function(object,par=NA,sep_d2=TRUE,numderiv=FALSE) {

  #object=fit; par=NA; numderiv=TRUE; sep_d2=TRUE
  
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

  #if(!all(is.na(par))) {response=eta_inv[["mu"]]}
  calc_lik_out=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=TRUE
            ,response=response,response_margin=response_margin,response_subject = response_subject)

  Fx_1_2=calc_lik_out$Fx_1_2; margin_p=calc_lik_out$margin_p; margin_d=calc_lik_out$margin_d; copula_d=calc_lik_out$copula_d

  ###Calculate derivaties: margin and copula d1 and d2
  margin_derivatives=calc_lik_out$margin_deriv
  copula_derivatives=calc_copula_derivatives(eta_inv, Fx_1_2, copula_dist,calc_d2 = TRUE)

  #Calculate numerical derivatives for F(x) - also used as a reference for overall likelihood derivative
  nd_impact_F=calc_Fx_derivatives(eta_inv,mm$x,margin_dist,response)
  nd_impact_F2=calc_Fx2_derivatives(eta_inv,mm$x,margin_dist,response)

  if(numderiv==TRUE) {
    #nd2_joint_lik=calc_true_SE_numderiv_only(eta_inv,mm,margin_dist,response,testing=TRUE,response_margin,response_subject)
    hessian_nd=calc_true_SE_numderiv_only_covariates(object=object,par=par_cov,mm=mm$x,margin_dist=margin_dist,response=response,testing=FALSE,response_margin=response_margin,response_subject=response_subject)
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
    for (par_name in margin_par) {
      if(object$include_dlcopdpar==TRUE | include_dlcopdpar==TRUE) {

        order_copula=matrix(ncol=4,nrow=0)
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
        margin_components_Ft_plus[,"time"]=margin_components_Ft_plus[,"time"]-1
        margin_plus=merge(margin_components,margin_components_Ft_plus,by=c("time","subject"),all.x=TRUE)

        copula_components=cbind(order_copula,dcdu1,dcdu2,copula_d,d2cdu12,d2cdu22,c_nd2)
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

  if(numderiv==TRUE) {
    vcov_final=-solve(hessian_nd)
    se_final=sqrt(abs(diag(solve(hessian_nd))))
  } else {
    vcov_final = -(solve((d2_mat)))/(length(response))
    se_final=sqrt(abs(diag(vcov_final)))
  }

  ###########TESTING APPROACH FOR ESTIMATING SMOOTHER VARIANCE
  
  # Method 1: Bayesian/Mixed Model Approach for Smoother Variance
  # Calculate variance-covariance matrix for smooth terms using: Var(β) = (X'WX + λP)^(-1) * σ²
  
  smooth_vcov_list = list()
  smooth_se_list = list()
  
  # Extract residual variance estimate (using reciprocal of mean weights as proxy for σ²)
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
        
        # Create difference penalty matrix (2nd order differences)
        # This assumes B-splines with difference penalty
        k = ncol(B)  # number of basis functions
        if(k > 2) {
          # Create second difference penalty matrix
          D2 = diff(diag(k), differences = 2)
          P = t(D2) %*% D2
        } else {
          # For very small basis, use identity penalty
          P = diag(k)
        }
        
        # Get working weights (diagonal of W matrix)
        # Use weights from the score function calculation
        if(!is.null(object$weights) && length(object$weights) == nrow(B)) {
          w_diag = object$weights
        } else {
          # Fallback: use unit weights
          w_diag = rep(1, nrow(B))
        }
        W = diag(w_diag)
        
        # Calculate the penalized precision matrix: X'WX + λP
        XWX = t(B) %*% W %*% B
        penalized_precision = XWX + lambda * P
        
        # Variance-covariance matrix for this smooth: (X'WX + λP)^(-1) * σ²
        tryCatch({
          smooth_vcov = solve(penalized_precision) * sigma2_est
          smooth_se = sqrt((diag(smooth_vcov)))
          
          # Store results
          smooth_vcov_list[[par_name]][[s_name]] = smooth_vcov
          smooth_se_list[[par_name]][[s_name]] = smooth_se
          
          # Also calculate the smoother matrix for fitted values variance
          # A = X(X'WX + λP)^(-1)X'W
          smoother_matrix = B %*% solve(penalized_precision) %*% t(B) %*% W
          fitted_se = sqrt(abs(as.vector(diag(smoother_matrix))) * sigma2_est)
          
          cat(sprintf("\nSmooth term variance estimates for %s:%s\n", par_name, s_name))
          cat(sprintf("  Basis coefficients SE: min=%.4f, max=%.4f, mean=%.4f\n", 
                     min(smooth_se), max(smooth_se), mean(smooth_se)))
          cat(sprintf("  Fitted values SE: min=%.4f, max=%.4f, mean=%.4f\n", 
                     min(fitted_se), max(fitted_se), mean(fitted_se)))
          cat(sprintf("  Effective DF: %.2f (trace of smoother matrix)\n", 
                     sum(diag(smoother_matrix))))
          cat(sprintf("  Smoothing parameter λ: %.4f\n", lambda))
          
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

  return(list(vcov=vcov_final_with_smooth, se=se_final_with_smooth))

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
#' - optional `vcov` output.
#' @export
summary.gamlss.longitudinal = function(
  object,
  include_vcov = TRUE,
  numderiv = TRUE,
  ci_level = 0.95,
  ...
) {
  if(!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  n_obs = length(object$response)
  n_subjects = length(unique(object$response_subject))
  n_timepoints = length(unique(object$response_margin))

  n_fixed = length(object$par)
  n_smooth_terms = 0
  if(!is.null(object$par_s)) {
    n_smooth_terms = sum(vapply(object$par_s, length, integer(1)))
  }

  edf_smooth = NA_real_
  if(!is.null(object$df_s) && length(object$df_s) > 0) {
    df_vals = suppressWarnings(as.numeric(unlist(object$df_s, use.names = FALSE)))
    df_vals = df_vals[is.finite(df_vals)]
    if(length(df_vals) > 0) {
      edf_smooth = sum(df_vals)
    }
  }

  loglik_joint = NA_real_
  if(!is.null(object$calc_lik_out_end) && !is.null(object$calc_lik_out_end$log_lik)) {
    if("joint" %in% names(object$calc_lik_out_end$log_lik)) {
      loglik_joint = as.numeric(object$calc_lik_out_end$log_lik["joint"])
    }
  }

  k_total = n_fixed + ifelse(is.finite(edf_smooth), edf_smooth, 0)
  aic = if(is.finite(loglik_joint)) -2 * loglik_joint + 2 * k_total else NA_real_
  bic = if(is.finite(loglik_joint)) -2 * loglik_joint + log(max(1, n_obs)) * k_total else NA_real_

  coef_tbl = data.frame(
    term = names(object$par),
    estimate = as.numeric(object$par),
    std_error = NA_real_,
    p_value = NA_real_,
    signif = NA_character_,
    stringsAsFactors = FALSE
  )

  coef_tbl$.original_order = seq_len(nrow(coef_tbl))
  coef_tbl$parameter = sub("\\..*$", "", coef_tbl$term)

  param_order = c("mu", "sigma", "nu", "tau", "theta", "zeta")
  coef_tbl$.param_rank = match(coef_tbl$parameter, param_order)
  coef_tbl$.param_rank[is.na(coef_tbl$.param_rank)] = length(param_order) + 1L

  vcov_out = NULL
  if(isTRUE(include_vcov)) {
    cat("Calculating variance-covariance matrix...\n")
    vcov_out = NULL
    invisible(utils::capture.output({
      vcov_out = do.call(
        vcov.gamlss.longitudinal,
        c(list(object = object, numderiv = numderiv), list(...))
      )
    }, type = "output"))

    if(!is.null(vcov_out$vcov) && !is.null(vcov_out$vcov$overall)) {
      V = vcov_out$vcov$overall
      se = NULL
      if(!is.null(vcov_out$se) && !is.null(vcov_out$se$overall)) {
        se = as.numeric(vcov_out$se$overall)
        se_names = names(vcov_out$se$overall)
      } else {
        se = sqrt(pmax(0, diag(V)))
        se_names = names(diag(V))
      }

      if(is.null(se_names) && !is.null(rownames(V)) && length(rownames(V)) == length(se)) {
        se_names = rownames(V)
      }

      if(!is.null(se_names)) {
        names(se) = se_names
      }

      if(!is.null(names(se))) {
        idx = match(coef_tbl$term, names(se))
        coef_tbl$std_error = se[idx]
      } else if(length(se) == nrow(coef_tbl)) {
        coef_tbl$std_error = se
      }

      z_abs = abs(coef_tbl$estimate / coef_tbl$std_error)
      coef_tbl$p_value = 2 * stats::pnorm(z_abs, lower.tail = FALSE)
      coef_tbl$signif = ifelse(
        is.na(coef_tbl$p_value),
        NA_character_,
        ifelse(coef_tbl$p_value < 0.001, "***",
               ifelse(coef_tbl$p_value < 0.01, "**",
                      ifelse(coef_tbl$p_value < 0.05, "*",
                             ifelse(coef_tbl$p_value < 0.1, ".", " "))))
      )
    }
  }

  coef_tbl = coef_tbl[order(coef_tbl$.param_rank, coef_tbl$.original_order), , drop = FALSE]
  rownames(coef_tbl) = NULL

  out = list(
    model = list(
      margin_dist = if(!is.null(object$margin_dist$family[1])) as.character(object$margin_dist$family[1]) else NA_character_,
      copula_dist = object$copula_dist,
      n_obs = n_obs,
      n_subjects = n_subjects,
      n_timepoints = n_timepoints,
      n_fixed = n_fixed,
      n_smooth_terms = n_smooth_terms,
      edf_smooth = edf_smooth
    ),
    fit = list(
      logLik = loglik_joint,
      AIC = aic,
      BIC = bic,
      ci_level = ci_level,
      vcov_included = isTRUE(include_vcov),
      vcov_numderiv = isTRUE(numderiv)
    ),
    smooth_terms = {
      st = list()
      if(!is.null(object$par_s) && length(object$par_s) > 0) {
        for(par_name in names(object$par_s)) {
          if(length(object$par_s[[par_name]]) == 0) next
          for(s_name in names(object$par_s[[par_name]])) {
            st[[length(st) + 1]] = data.frame(
              parameter = par_name,
              smooth_term = s_name,
              stringsAsFactors = FALSE
            )
          }
        }
      }
      if(length(st) == 0) {
        data.frame(parameter = character(0), smooth_term = character(0), stringsAsFactors = FALSE)
      } else {
        do.call(rbind, st)
      }
    },
    coefficients = within(coef_tbl, {
      .original_order = NULL
      .param_rank = NULL
    }),
    vcov = vcov_out
  )
  class(out) = "summary.gamlss.longitudinal"
  out
}

#' @export
print.summary.gamlss.longitudinal = function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nGAMLSS Longitudinal Model Summary\n")
  cat("--------------------------------\n")
  cat("Margin distribution:", x$model$margin_dist, "\n")
  cat("Copula distribution:", x$model$copula_dist, "\n")
  cat("Observations:", x$model$n_obs,
      " | Subjects:", x$model$n_subjects,
      " | Time points:", x$model$n_timepoints, "\n")
  cat("Fixed coefficients:", x$model$n_fixed,
      " | Smooth terms:", x$model$n_smooth_terms,
      " | Smooth EDF:", format(round(x$model$edf_smooth, digits), nsmall = 2), "\n")

  cat("\nFixed coefficients:\n")
  cat("--------------------\n")
  coef_tbl = x$coefficients
  coef_tbl$estimate = round(coef_tbl$estimate, digits)
  coef_tbl$std_error = round(coef_tbl$std_error, digits)
  coef_tbl$p_value = round(coef_tbl$p_value, digits + 1)

  fmt_num = function(v, d) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d))
  coef_disp = data.frame(
    term = as.character(coef_tbl$term),
    estimate = fmt_num(coef_tbl$estimate, digits),
    std_error = fmt_num(coef_tbl$std_error, digits),
    p_value = fmt_num(coef_tbl$p_value, digits + 1),
    signif = ifelse(is.na(coef_tbl$signif), "", as.character(coef_tbl$signif)),
    parameter = as.character(coef_tbl$parameter),
    stringsAsFactors = FALSE
  )

  w_term = max(nchar("term"), nchar(coef_disp$term, type = "width"), na.rm = TRUE)
  w_est = max(nchar("estimate"), nchar(coef_disp$estimate, type = "width"), na.rm = TRUE)
  w_se = max(nchar("std_error"), nchar(coef_disp$std_error, type = "width"), na.rm = TRUE)
  w_p = max(nchar("p_value"), nchar(coef_disp$p_value, type = "width"), na.rm = TRUE)
  w_sig = max(nchar("signif"), nchar(coef_disp$signif, type = "width"), na.rm = TRUE)

  format_row = function(term, estimate, std_error, p_value, signif) {
    sprintf(
      "%-*s  %*s  %*s  %*s  %-*s",
      w_term, term,
      w_est, estimate,
      w_se, std_error,
      w_p, p_value,
      w_sig, signif
    )
  }

  print_coef_block = function(block, prefix = "    ") {
    hdr = format_row("term", "estimate", "std_error", "p_value", "signif")
    cat(prefix, hdr, "\n", sep = "")
    for(ii in seq_len(nrow(block))) {
      row_txt = format_row(
        block$term[ii],
        block$estimate[ii],
        block$std_error[ii],
        block$p_value[ii],
        block$signif[ii]
      )
      cat(prefix, row_txt, "\n", sep = "")
    }
  }

  param_order = c("mu", "sigma", "nu", "tau", "theta", "zeta")
  params_present = unique(coef_tbl$parameter)
  params_print = c(param_order[param_order %in% params_present], setdiff(params_present, param_order))

  for(k in seq_along(params_print)) {
    p = params_print[k]
    block = coef_disp[coef_disp$parameter == p, c("term", "estimate", "std_error", "p_value", "signif"), drop = FALSE]
    cat(sprintf("  [%s]\n", p))
    print_coef_block(block, prefix = "    ")
    if(k < length(params_print)) cat("  --------------------\n")
  }

  cat("  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  
  cat("\nSmooth terms:\n")
  cat("--------------------\n")
  if(!is.null(x$smooth_terms) && nrow(x$smooth_terms) > 0) {
    print(x$smooth_terms, row.names = FALSE)
    cat("Use plot(object) to visualize smooth and fixed terms with confidence bands.\n")
  } else {
    cat("None\n")
  }

  cat("\nFit statistics:\n")
  cat("--------------------\n")
  fit_tbl = data.frame(
    metric = c("logLik", "AIC", "BIC"),
    value = c(x$fit$logLik, x$fit$AIC, x$fit$BIC),
    stringsAsFactors = FALSE
  )
  print(fit_tbl, row.names = FALSE, digits = digits)
  cat("--------------------------------\n")

  invisible(x)
}

#' Plot all smooth terms with confidence bands
#'
#' This utility plots every smooth term in a fitted `gamlss.longitudinal` object
#' and computes pointwise confidence bands using the smooth coefficient
#' covariance matrices returned by `vcov.gamlss.longitudinal()`.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param vcov_obj Optional output from `vcov(object, ...)`. If `NULL`, this is
#' computed internally with `vcov(object, numderiv = TRUE)`.
#' @param data Optional data frame containing original covariates used for the
#' x-axis variable of each smooth.
#' @param ci_level Confidence level for pointwise intervals.
#' @param ncol Number of plot columns (defaults to 2 or fewer if needed).
#' @param ci_col Color for confidence bands.
#' @param fit_col Color for fitted smooth line.
#' @param ci_lty Line type for confidence bands.
#' @param fit_lwd Line width for fitted smooth.
#' @param sort_x Logical; sort points by x before plotting lines.
#' @param fallback_to_index Logical; if x variable cannot be inferred, plot
#' against row index.
#' @param setup_mfrow Logical; if TRUE (default), configure par(mfrow) inside
#' this function. Set FALSE when caller configures layout.
#' @param show_legend Logical; if TRUE, draw a small legend in each panel.
#'
#' @return Invisibly returns a nested list with x, fitted values, standard
#' errors, and confidence limits for each smooth term.
#' @export
plot_smooth_terms = function(
  object,
  vcov_obj = NULL,
  data = NULL,
  ci_level = 0.95,
  ncol = NULL,
  ci_col = "red",
  fit_col = "black",
  ci_lty = 2,
  fit_lwd = 2,
  sort_x = TRUE,
  fallback_to_index = TRUE,
  setup_mfrow = TRUE,
  show_legend = TRUE
) {
  if(!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  if(is.null(vcov_obj)) {
    vcov_obj = vcov(object, numderiv = TRUE)
  }

  if(!is.list(vcov_obj) || is.null(vcov_obj$vcov)) {
    stop("'vcov_obj' must be the output of vcov.gamlss.longitudinal().")
  }

  smooth_vcov_list = vcov_obj$vcov$smooth_vcov
  smooth_se_list = vcov_obj$vcov$smooth_se

  extract_smooth_var = function(s_name) {
    out = sub("^s\\(([^,\\)]+).*$", "\\1", s_name)
    out = trimws(out)
    out = gsub("`", "", out, fixed = TRUE)
    out
  }

  get_x_for_smooth = function(par_name, s_name, B) {
    x_var = extract_smooth_var(s_name)
    x = NULL

    if(!is.null(data) && is.data.frame(data)) {
      data_names = names(data)
      # Prefer exact match, then case-insensitive match, then make.names match.
      idx_exact = which(data_names == x_var)
      idx_ci = which(tolower(data_names) == tolower(x_var))
      idx_mn = which(make.names(data_names) == make.names(x_var))
      idx = c(idx_exact, idx_ci, idx_mn)
      idx = idx[!duplicated(idx)]
      if(length(idx) > 0) {
        matched_name = data_names[idx[1]]
        x = data[[matched_name]]
        x_var = matched_name
      }
    }

    if(is.null(x) && !is.null(object$model_matrix$x[[par_name]]) && x_var %in% colnames(object$model_matrix$x[[par_name]])) {
      x = object$model_matrix$x[[par_name]][, x_var]
    }

    if(is.null(x) && fallback_to_index) {
      x = seq_len(nrow(B))
      x_var = "index"
      warning("Falling back to index for smooth term '", s_name,
              "' because covariate was not found in supplied data/model matrix.")
    }

    if(is.null(x)) {
      stop("Could not infer x-axis for smooth term '", s_name, "'. Provide 'data' with the smooth covariate columns.")
    }

    if(length(x) != nrow(B)) {
      if(length(x) > nrow(B)) {
        x = x[seq_len(nrow(B))]
      } else {
        stop("Length mismatch for smooth term '", s_name, "': length(x)=", length(x), " but nrow(B)=", nrow(B), ".")
      }
    }

    list(x = as.numeric(x), x_var = x_var)
  }

  z = qnorm((1 + ci_level) / 2)
  smooth_index = list()

  for(par_name in names(object$par_s)) {
    if(length(object$par_s[[par_name]]) == 0) next
    for(s_name in names(object$par_s[[par_name]])) {
      B = object$model_matrix$s[[par_name]][[s_name]]
      beta_s = object$par_s[[par_name]][[s_name]]
      if(is.null(B) || is.null(beta_s)) next

      smooth_index[[length(smooth_index) + 1]] = list(par_name = par_name, s_name = s_name)
    }
  }

  n_plots = length(smooth_index)
  if(n_plots == 0) {
    warning("No smooth terms found to plot.")
    return(invisible(list()))
  }

  if(is.null(ncol)) {
    ncol = min(2, n_plots)
  }
  nrow = ceiling(n_plots / ncol)

  if(setup_mfrow) {
    old_par = par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(nrow, ncol))
  }

  out = list()
  for(i in seq_len(n_plots)) {
    par_name = smooth_index[[i]]$par_name
    s_name = smooth_index[[i]]$s_name

    B = object$model_matrix$s[[par_name]][[s_name]]
    beta_s = object$par_s[[par_name]][[s_name]]
    x_info = get_x_for_smooth(par_name, s_name, B)
    x = x_info$x

    fitted_smooth = as.numeric(B %*% beta_s)

    smooth_vcov = NULL
    smooth_se = NULL
    if(!is.null(smooth_vcov_list) && !is.null(smooth_vcov_list[[par_name]])) {
      smooth_vcov = smooth_vcov_list[[par_name]][[s_name]]
    }
    if(!is.null(smooth_se_list) && !is.null(smooth_se_list[[par_name]])) {
      smooth_se = smooth_se_list[[par_name]][[s_name]]
    }

    if(!is.null(smooth_vcov) && all(dim(smooth_vcov) == c(ncol(B), ncol(B)))) {
      smooth_fit_se = sqrt(pmax(0, diag(B %*% smooth_vcov %*% t(B))))
    } else if(!is.null(smooth_se) && length(smooth_se) == ncol(B)) {
      beta_var_diag = as.numeric(smooth_se)^2
      smooth_fit_se = sqrt(pmax(0, rowSums((B^2) * rep(beta_var_diag, each = nrow(B)))))
    } else {
      smooth_fit_se = rep(NA_real_, nrow(B))
    }

    ci_lower = fitted_smooth - z * smooth_fit_se
    ci_upper = fitted_smooth + z * smooth_fit_se
    ord = if(sort_x) order(x) else seq_along(x)

    main_title = paste(par_name, s_name, sep = ": ")
    ylab_text = paste("smooth(", x_info$x_var, ")", sep = "")

    y_vals = c(fitted_smooth, ci_lower, ci_upper)
    y_vals = y_vals[is.finite(y_vals)]
    if(length(y_vals) > 0) {
      y_rng = range(y_vals)
      y_pad = 0.05 * max(1e-8, diff(y_rng))
      y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
    } else {
      y_lim = NULL
    }

    plot(x[ord], fitted_smooth[ord], type = "l", lwd = fit_lwd, col = fit_col,
         main = main_title, xlab = x_info$x_var, ylab = ylab_text, ylim = y_lim)

    if(any(is.finite(smooth_fit_se))) {
      lines(x[ord], ci_lower[ord], col = ci_col, lty = ci_lty)
      lines(x[ord], ci_upper[ord], col = ci_col, lty = ci_lty)
    }

    if(show_legend) {
      legend("topright",
             legend = c("fit", paste0(round(ci_level * 100), "% CI")),
             col = c(fit_col, ci_col),
             lty = c(1, ci_lty),
             lwd = c(fit_lwd, 1),
             bty = "n",
             cex = 0.8)
    }

    if(is.null(out[[par_name]])) out[[par_name]] = list()
    out[[par_name]][[s_name]] = list(
      x = x,
      fitted = fitted_smooth,
      se = smooth_fit_se,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )
  }

  invisible(out)
}

#' Plot all fixed terms with confidence bands
#'
#' This utility plots fixed-effect term contributions for a fitted
#' `gamlss.longitudinal` object using coefficient uncertainty from
#' `vcov.gamlss.longitudinal()`.
#'
#' For each fixed-effect design-matrix column $x_j$, it plots
#' $x_j \hat\beta_j$ with pointwise confidence bands
#' $x_j \hat\beta_j \pm z_{\alpha/2}\sqrt{x_j^2 \mathrm{Var}(\hat\beta_j)}$.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param vcov_obj Optional output from `vcov(object, ...)`. If `NULL`, this is
#' computed internally with `vcov(object, numderiv = TRUE)`.
#' @param ci_level Confidence level for pointwise intervals.
#' @param ncol Number of plot columns (defaults to 2 or fewer if needed).
#' @param include_intercept Logical; include intercept columns in plots.
#' @param ci_col Color for confidence bands.
#' @param fit_col Color for fitted fixed-term line.
#' @param ci_lty Line type for confidence bands.
#' @param fit_lwd Line width for fitted fixed-term line.
#' @param sort_x Logical; sort x-values before drawing lines.
#' @param fallback_to_index Logical; if x has one unique value, use index on x-axis.
#' @param setup_mfrow Logical; if TRUE (default), configure par(mfrow) inside
#' this function. Set FALSE when caller configures layout.
#' @param data Optional data frame used to detect factor columns and show
#' factor levels on x-axis for categorical fixed terms.
#' @param factor_pch Point symbol for factor-level estimates.
#' @param factor_cex Point size for factor-level estimates.
#' @param show_legend Logical; if TRUE, draw a small legend in each panel.
#'
#' @return Invisibly returns a nested list with x, fitted values, standard
#' errors, and confidence limits for each fixed term.
#' @export
plot_fixed_terms = function(
  object,
  vcov_obj = NULL,
  ci_level = 0.95,
  ncol = NULL,
  include_intercept = FALSE,
  ci_col = "red",
  fit_col = "black",
  ci_lty = 2,
  fit_lwd = 2,
  sort_x = TRUE,
  fallback_to_index = TRUE,
  setup_mfrow = TRUE,
  data = NULL,
  factor_pch = 16,
  factor_cex = 1.2,
  show_legend = TRUE
) {
  if(!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  if(is.null(vcov_obj)) {
    vcov_obj = vcov(object, numderiv = TRUE)
  }

  if(!is.list(vcov_obj) || is.null(vcov_obj$vcov) || is.null(vcov_obj$vcov$overall)) {
    stop("'vcov_obj' must be the output of vcov.gamlss.longitudinal() with vcov$overall present.")
  }

  V = vcov_obj$vcov$overall
  if(is.null(rownames(V)) || is.null(colnames(V))) {
    stop("vcov$overall must have row and column names matching fixed coefficients.")
  }

  z = qnorm((1 + ci_level) / 2)
  build_factor_groups = function(X, data) {
    groups = list()
    if(is.null(data) || !is.data.frame(data) || is.null(X) || ncol(X) == 0) return(groups)

    x_cols = colnames(X)
    for(var_name in names(data)) {
      v = data[[var_name]]
      if(!is.factor(v)) next
      levs = levels(v)
      if(length(levs) < 2) next

      level_col_map = list()
      matched_cols = character(0)
      for(lev in levs[-1]) {
        candidates = unique(c(
          paste0(var_name, lev),
          paste0(var_name, make.names(lev)),
          paste0(var_name, "_", lev),
          paste0(var_name, "_", make.names(lev))
        ))
        hit = candidates[candidates %in% x_cols]
        if(length(hit) > 0) {
          level_col_map[[lev]] = hit[1]
          matched_cols = c(matched_cols, hit[1])
        }
      }

      if(length(matched_cols) == 0 && length(levs) == 2 && var_name %in% x_cols) {
        level_col_map[[levs[2]]] = var_name
        matched_cols = var_name
      }

      if(length(level_col_map) > 0) {
        groups[[var_name]] = list(
          var_name = var_name,
          levels = levs,
          ref_level = levs[1],
          level_col_map = level_col_map,
          matched_cols = unique(matched_cols)
        )
      }
    }
    groups
  }

  plot_specs = list()
  for(par_name in names(object$model_matrix$x)) {
    X = object$model_matrix$x[[par_name]]
    if(is.null(X) || ncol(X) == 0) next

    factor_groups = build_factor_groups(X, data)
    grouped_cols = unique(unlist(lapply(factor_groups, function(g) g$matched_cols), use.names = FALSE))
    if(length(grouped_cols) == 0) grouped_cols = character(0)

    for(var_name in names(factor_groups)) {
      fg = factor_groups[[var_name]]
      has_valid_coef = FALSE
      for(lev in names(fg$level_col_map)) {
        coef_name = paste(par_name, fg$level_col_map[[lev]], sep = ".")
        if(coef_name %in% names(object$par) && coef_name %in% rownames(V) && coef_name %in% colnames(V)) {
          has_valid_coef = TRUE
          break
        }
      }
      if(has_valid_coef) {
        plot_specs[[length(plot_specs) + 1]] = list(
          type = "factor",
          par_name = par_name,
          var_name = var_name,
          group = fg
        )
      }
    }

    for(col_name in colnames(X)) {
      if(col_name %in% grouped_cols) next
      if(!include_intercept && col_name == "intercept") next
      coef_name = paste(par_name, col_name, sep = ".")
      if(!coef_name %in% names(object$par)) next
      if(!coef_name %in% rownames(V) || !coef_name %in% colnames(V)) next

      plot_specs[[length(plot_specs) + 1]] = list(
        type = "continuous",
        par_name = par_name,
        col_name = col_name,
        coef_name = coef_name
      )
    }
  }

  n_plots = length(plot_specs)
  if(n_plots == 0) {
    warning("No fixed terms found to plot with matching vcov entries.")
    return(invisible(list()))
  }

  if(is.null(ncol)) {
    ncol = min(2, n_plots)
  }
  nrow = ceiling(n_plots / ncol)

  if(setup_mfrow) {
    old_par = par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(nrow, ncol))
  }

  out = list()
  for(i in seq_len(n_plots)) {
    spec = plot_specs[[i]]
    par_name = spec$par_name

    if(identical(spec$type, "factor")) {
      fg = spec$group
      levs = fg$levels
      x_plot = seq_along(levs)
      fitted_term = rep(NA_real_, length(levs))
      term_se = rep(NA_real_, length(levs))

      for(j in seq_along(levs)) {
        lev = levs[j]
        if(identical(lev, fg$ref_level)) {
          fitted_term[j] = 0
          term_se[j] = 0
        } else if(lev %in% names(fg$level_col_map)) {
          col_name_lev = fg$level_col_map[[lev]]
          coef_name_lev = paste(par_name, col_name_lev, sep = ".")
          if(coef_name_lev %in% names(object$par) && coef_name_lev %in% rownames(V) && coef_name_lev %in% colnames(V)) {
            fitted_term[j] = as.numeric(object$par[coef_name_lev])
            term_se[j] = sqrt(pmax(0, as.numeric(V[coef_name_lev, coef_name_lev])))
          }
        }
      }

      keep = is.finite(fitted_term) & is.finite(term_se)
      ci_lower = fitted_term - z * term_se
      ci_upper = fitted_term + z * term_se

      y_vals = c(fitted_term[keep], ci_lower[keep], ci_upper[keep])
      if(length(y_vals) > 0) {
        y_rng = range(y_vals)
        y_pad = 0.05 * max(1e-8, diff(y_rng))
        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
      } else {
        y_lim = NULL
      }

      plot(x_plot[keep], fitted_term[keep],
           type = "p",
           pch = factor_pch,
           cex = factor_cex,
           col = fit_col,
           xaxt = "n",
           main = paste(par_name, fg$var_name, sep = ": "),
           xlab = fg$var_name,
           ylab = paste("fixed contribution:", paste(par_name, fg$var_name, sep = ".")),
           ylim = y_lim)
      axis(1, at = x_plot, labels = as.character(levs))
      arrows(x_plot[keep], ci_lower[keep], x_plot[keep], ci_upper[keep],
             angle = 90, code = 3, length = 0.05, col = ci_col, lty = ci_lty)
      abline(h = 0, col = "grey70", lty = 3)

      if(show_legend) {
        legend("topright",
               legend = c("estimate", paste0(round(ci_level * 100), "% CI")),
               col = c(fit_col, ci_col),
               pch = c(factor_pch, NA),
               lty = c(NA, ci_lty),
               bty = "n",
               cex = 0.8)
      }

      if(is.null(out[[par_name]])) out[[par_name]] = list()
      out[[par_name]][[fg$var_name]] = list(
        coefficient = paste(par_name, fg$var_name, sep = "."),
        x = x_plot,
        levels = levs,
        fitted = fitted_term,
        se = term_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper
      )
    } else {
      col_name = spec$col_name
      coef_name = spec$coef_name
      X = object$model_matrix$x[[par_name]]
      x_raw = as.numeric(X[, col_name])
      beta_hat = as.numeric(object$par[coef_name])
      var_beta = as.numeric(V[coef_name, coef_name])

      fitted_term = x_raw * beta_hat
      term_se = sqrt(pmax(0, (x_raw^2) * var_beta))
      ci_lower = fitted_term - z * term_se
      ci_upper = fitted_term + z * term_se

      if(length(unique(x_raw)) <= 1 && fallback_to_index) {
        x_plot = seq_along(x_raw)
        xlab_text = paste(col_name, "(index)")
        ord = seq_along(x_plot)
      } else {
        x_plot = x_raw
        xlab_text = col_name
        ord = if(sort_x) order(x_plot) else seq_along(x_plot)
      }

      main_title = paste(par_name, col_name, sep = ": ")
      ylab_text = paste("fixed contribution:", coef_name)

      y_vals = c(fitted_term, ci_lower, ci_upper)
      y_vals = y_vals[is.finite(y_vals)]
      if(length(y_vals) > 0) {
        y_rng = range(y_vals)
        y_pad = 0.05 * max(1e-8, diff(y_rng))
        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
      } else {
        y_lim = NULL
      }

      plot(x_plot[ord], fitted_term[ord], type = "l", lwd = fit_lwd, col = fit_col,
           main = main_title, xlab = xlab_text, ylab = ylab_text, ylim = y_lim)
      lines(x_plot[ord], ci_lower[ord], col = ci_col, lty = ci_lty)
      lines(x_plot[ord], ci_upper[ord], col = ci_col, lty = ci_lty)

      if(show_legend) {
        legend("topright",
               legend = c("fit", paste0(round(ci_level * 100), "% CI")),
               col = c(fit_col, ci_col),
               lty = c(1, ci_lty),
               lwd = c(fit_lwd, 1),
               bty = "n",
               cex = 0.8)
      }
      if(is.null(out[[par_name]])) out[[par_name]] = list()
      out[[par_name]][[col_name]] = list(
        coefficient = coef_name,
        x = x_plot,
        fitted = fitted_term,
        se = term_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper
      )
    }
  }

  invisible(out)
}

#' Plot method for fitted gamlss.longitudinal objects
#'
#' This S3 plot method generates diagnostic plots for a fitted
#' `gamlss.longitudinal` object, including both fixed-effect terms and
#' smooth terms with confidence bands. If the total number of plots
#' exceeds `max_plots_per_page`, the user is prompted to advance pages.
#'
#' @param x A fitted `gamlss.longitudinal` object.
#' @param y Unused; included for S3 generic compatibility.
#' @param ci_level Confidence level for pointwise intervals (default: 0.95).
#' @param max_plots_per_page Maximum number of plots before pagination prompt
#' (default: 6). Set to 0 or negative to disable pagination.
#' @param ncol Number of columns in each plot frame (default: 2).
#' @param include_intercept Logical; include intercept terms in fixed plots
#' (default: FALSE).
#' @param data Optional data frame containing original covariates used to infer
#' smooth-term x-axis variables. If omitted, smooth terms may fall back to index.
#' @param ci_col Color for confidence bands (default: "red").
#' @param fit_col Color for fitted lines (default: "black").
#' @param show_legend Logical; if TRUE, draw a legend in each panel.
#' @param ... Additional arguments passed to plotting functions (unused).
#'
#' @return Invisibly returns a list with elements:
#' - `smooth_terms`: output from `plot_smooth_terms()`
#' - `fixed_terms`: output from `plot_fixed_terms()`
#'
#' @export
plot.gamlss.longitudinal = function(
  x,
  y,
  data = NULL,
  ci_level = 0.95,
  max_plots_per_page = 6,
  ncol = 2,
  include_intercept = FALSE,
  ci_col = "red",
  fit_col = "black",
  show_legend = TRUE,
  ...
) {
  cat("\n=== Plotting gamlss.longitudinal object ===\n")
  cat("Computing variance-covariance matrix...\n\n")

  vcov_obj = vcov(x, numderiv = TRUE)

  count_plot_terms = function(obj) {
    n_smooth = 0
    n_fixed = 0

    for(par_name in names(obj$par_s)) {
      if(length(obj$par_s[[par_name]]) > 0) {
        n_smooth = n_smooth + length(obj$par_s[[par_name]])
      }
    }

    for(par_name in names(obj$model_matrix$x)) {
      X = obj$model_matrix$x[[par_name]]
      if(!is.null(X) && ncol(X) > 0) {
        n_cols = ncol(X)
        if(!include_intercept && "intercept" %in% colnames(X)) {
          n_cols = n_cols - 1
        }
        coef_names = paste(par_name, colnames(X), sep = ".")
        coef_names = coef_names[!(colnames(X) == "intercept" & !include_intercept)]
        n_valid = sum(coef_names %in% names(obj$par))
        n_fixed = n_fixed + n_valid
      }
    }

    list(smooth = n_smooth, fixed = n_fixed, total = n_smooth + n_fixed)
  }

  counts = count_plot_terms(x)
  cat(sprintf("Found %d smooth terms and %d fixed terms (total: %d plots).\n\n",
              counts$smooth, counts$fixed, counts$total))

  smooth_results = fixed_results = list()

  if(counts$total == 0) {
    warning("No plots to display.")
    return(invisible(list(smooth_terms = list(), fixed_terms = list())))
  }

  nrow_total = ceiling(counts$total / ncol)
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nrow_total, ncol))

  cat("Plotting all terms...\n")
  if(counts$smooth > 0) {
    smooth_results = plot_smooth_terms(
      object = x,
      vcov_obj = vcov_obj,
      data = data,
      ci_level = ci_level,
      ncol = ncol,
      ci_col = ci_col,
      fit_col = fit_col,
      setup_mfrow = FALSE,
      show_legend = show_legend
    )
  }

  if(counts$fixed > 0) {
    fixed_results = plot_fixed_terms(
      object = x,
      vcov_obj = vcov_obj,
      ci_level = ci_level,
      ncol = ncol,
      include_intercept = include_intercept,
      ci_col = ci_col,
      fit_col = fit_col,
      setup_mfrow = FALSE,
      data = data,
      show_legend = show_legend
    )
  }

  cat("\n")
  invisible(list(
    smooth_terms = smooth_results,
    fixed_terms = fixed_results
  ))
}

#' @export
eta_to_par=function(eta,margin_dist,copula_dist) {
  par=eta*0
  for (par_name in names(eta)) {
    if(par_name %in% names(margin_dist$parameters)) {
      FUN = eval(parse(text=paste(paste(paste("margin_dist$",par_name,sep=""),"linkinv",sep="."))))
      par[par_name]=FUN(eta[par_name])
    }
    if(par_name %in% names(copula_dist$parameters)) {
      FUN = eval(parse(text=paste(paste(paste("copula_dist$",par_name,sep=""),"linkinv",sep="."))))
      par[par_name]=FUN(eta[par_name])
    }
  }
  return(par)
}
#' @export
get_starting_values = function(copula_dist,margin_dist,dataset,eta_transform=FALSE) {

  margin_names=unique(dataset$time)
  num_margins=length(margin_names)

  tau_start=cor(dataset[dataset$time%in%(margin_names[1:(num_margins-1)]),"response"]
                ,dataset[dataset$time%in%(margin_names[2:(num_margins)]),"response"],method="kendall",use="complete.obs")
  if(!is.finite(tau_start)) {
    warning("Non-finite Kendall tau in get_starting_values(); using tau = 0 for copula initialisation.")
    tau_start=0
  }
  tau_start=max(min(tau_start,0.9999),-0.9999)

  cop_par=BiCopTau2Par(family=as.numeric(BiCopName(copula_dist)),tau=tau_start)
  names(cop_par)=get_copula_dist(copula_dist)$parameters

  if(margin_dist$family[1]=="GA" | margin_dist$family[1]=="EXP") {
    margin_par=c(
      mean(dataset$response)
      , sd(dataset$response)/mean(dataset$response)
      , skewness(dataset$response)
      , kurtosis(dataset$response)
    )
  } else if (margin_dist$family[1]=="NO") {
    margin_par=c(
      mean(dataset$response)
      , sd(dataset$response)
    )
  } else if (margin_dist$family[1]=="PO") {
    margin_par=c(
      mean(dataset$response)
    )
  } else if (margin_dist$family[1]=="NBI") {
    margin_par=c(
      mean(dataset$response),
      sd(dataset$response)/mean(dataset$response)
    )
  } else {
    cat("Fitting initial GAMLSS model for margin to obtain starting values...\n")
    # Deliberately low-iteration startup fit; silence expected convergence warnings.
    start_fit=suppressWarnings(suppressMessages(
      gamlss(dataset$response~1, family=margin_dist,method=RS(1))
    ))
    margin_par=unlist(coefAll(start_fit))
    names(margin_par)=names(margin_dist$parameters)
    #margin_par=eta_to_par(margin_par_temp,margin_dist,get_copula_dist(copula_dist))
    eta_transform=FALSE
  }

  names(margin_par)=names(margin_dist$parameters)
  margin_par=margin_par[!is.na(names(margin_par))]

  if(eta_transform==TRUE) {
    margin_par_eta=margin_par
    cop_par_eta=cop_par

    for (par_name in names(margin_par)) {
      FUN = eval(parse(text=paste(paste(paste("margin_dist$",par_name,sep=""),"linkfun",sep="."))))
      margin_par_eta[par_name]=FUN(margin_par[par_name])
    }

    for (par_name in names(cop_par)) {
      cop_par_eta[par_name]=get_copula_dist(copula_dist)$copula_link[[paste(par_name,".linkfun",sep="")]](cop_par[par_name])
    }

    return_list=c(margin_par_eta,cop_par_eta)
  } else {
    return_list=c(margin_par,cop_par)
  }

  return(return_list)
}
#' @export
par_to_eta = function(par,copula_dist,margin_dist) {

  margin_par=par[names(margin_dist$parameters)]
  names(margin_par)=names(margin_dist$parameters)

  cop_par=par[get_copula_dist(copula_dist)$parameters]
  names(cop_par)=get_copula_dist(copula_dist)$parameters

    margin_par_eta=margin_par
    cop_par_eta=cop_par

    for (par_name in names(margin_par)) {
      FUN = eval(parse(text=paste(paste(paste("margin_dist$",par_name,sep=""),"linkfun",sep="."))))
      margin_par_eta[par_name]=FUN(margin_par[par_name])
    }

    for (par_name in names(cop_par)) {
      cop_par_eta[par_name]=get_copula_dist(copula_dist)$copula_link[[paste(par_name,".linkfun",sep="")]](cop_par[par_name])
      names(cop_par_eta)=names(cop_par)
    }

    return_list=c(margin_par_eta,cop_par_eta)

  return(return_list)
}

#' @export
fit_jointreg_nocov <- function(input_par,margin_dist,copula_dist,data
                               , use_dlcopdpar=TRUE, verbose=TRUE, plot_results=TRUE
                               , crit_lik_change=0.05, start_step_size=.5, step_adjustment=.9, max_steps=5
                               , true_val = NA) {

  log_lik_history=matrix(ncol=3+2,nrow=0)
  par_history=matrix(ncol=length(input_par)+2,nrow=0)

  ### Run fit for separate and joint optimisation
    copula_deriv=if(use_dlcopdpar==TRUE){1}else{0}
    ### CORE ITERATION
    change=1;log_lik_start=0;log_lik_change=1000;run_counter=1;step_size=start_step_size;
    while (abs(log_lik_change)>crit_lik_change) {
      step_size=step_size*(step_adjustment^min(max_steps,run_counter))
      par_history=rbind(par_history,c(copula_deriv,run_counter,input_par))

      #Run optimisation
      outer_optim_output=optim_outer(par=input_par,dataset,margin_dist,copula_dist,use_dlcopdpar=use_dlcopdpar,verbose=FALSE,step_size=step_size)

      #Capture outputs
      input_par=outer_optim_output$par_end
      change=sum(outer_optim_output$par_change)
      #print(outer_optim_output$log_lik)

      log_lik=outer_optim_output$log_lik["joint"]
      log_lik_change=log_lik-log_lik_start
      log_lik_start=log_lik

      #Capture changes in parameters
      log_lik_history=rbind(log_lik_history,c(copula_deriv,run_counter,outer_optim_output$log_lik))
      run_counter=run_counter+1

    }
    par_history=rbind(par_history,c(copula_deriv,run_counter,input_par))

    outer_optim_output=optim_outer(par=input_par,dataset,margin_dist,copula_dist,use_dlcopdpar=use_dlcopdpar,verbose=FALSE,step_size=step_size)
    log_lik=outer_optim_output$log_lik["joint"]
    log_lik_change=log_lik-log_lik_start
    log_lik_start=log_lik
    log_lik_history=rbind(log_lik_history,c(copula_deriv,run_counter,outer_optim_output$log_lik))

    colnames(log_lik_history)[1:2]=colnames(par_history)[1:2]=c("use_dlcopdpar","run_counter")

  #Plot likelihood and parameters
  if(plot_results==TRUE) {

    plot.new()
    par_count=round(sqrt((ncol(par_history)+1)),0)+1
    par(mfrow=c(par_count,par_count))

    for (i in colnames(log_lik_history)[3:5]) {
      plot( log_lik_history[,i],xlab="Iteration",ylab="LogLik",main=i,type = "l",col="blue",xlim=c(1,max(log_lik_history[,"run_counter"])))
      #lines(log_lik_history[log_lik_history[,"use_dlcopdpar"]==0,i],xlab="LogLik",ylab="Iteration",main=i,type = "l",col="red",xlim=c(1,max(log_lik_history[,"run_counter"])),ylim=range(log_lik_history[,i]))
      #legend("bottomright",c("Joint","Separate"), lwd=c(5,2), col=c("blue","red"))

    }

    for (i in 1:(ncol(par_history)-2)) {
      #lines(par_nodlcop[,i+2],col="red",type="l")
      if (!all(is.na(true_val))) {
        plot(par_history[,i+2],col="blue",type="l",main=colnames(par_history)[i+2],ylab="Parameter estimate",ylim=range(c(par_history[,i+2],true_val[i])))
        abline(h=true_val[i])
      } else {
        plot(par_history[,i+2],col="blue",type="l",main=colnames(par_history)[i+2],ylab="Parameter estimate")
      }
      #legend("bottomright",c("Joint","Separate"), lwd=c(5,2), col=c("blue","red"))
    }
  }

  return_list=list(par_history,log_lik_history)
  names(return_list)=c("par_history","log_lik_history")
  return(return_list)
}




#' @export
calc_F_x <- function(eta_inv,mm,margin_dist,response) {
  #Setup input matrix of response and parameters
  #margin_names=unique(response_margin)
  #num_margins=length(margin_names)

  margin_deriv_input=list()
  margin_deriv_input[["y"]]=response
  margin_deriv_input[["q"]]=response
  margin_deriv_input[["x"]]=response
  for (par_name in names(mm)) {
    if (par_name %in% c("mu","sigma","nu","tau")) {
      margin_deriv_input[[par_name]]=eta_inv[[par_name]]
    }
  }

  margin_pFUN=eval(parse( text=paste("p",margin_dist$family[1],sep="") ))
  FUN=margin_pFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_p=do.call(FUN,args=margin_deriv_input[FUN_args])

  return(margin_p)
}

#' @export
get_copula_dist=function(copula_dist) {

  if(copula_dist=="C" | copula_dist=="Clayton") {
    copula_link=list(log,exp,dloginv=exp); two_par_cop=FALSE
    copula_dist=BiCopName(copula_dist)
    parameters=c("theta")
  }
  else if(copula_dist=="N" | copula_dist=="Normal") {
    copula_link=list(logit,logit_inv,dlogit_inv); two_par_cop=FALSE
    copula_dist=BiCopName(copula_dist)
    parameters=c("theta")
  } else if(copula_dist=="t" | copula_dist=="Normal") {
    copula_link=list(logit,logit_inv,dlogit_inv,log_2plus,log_2plus_inv,dlog_2plus_inv); two_par_cop=TRUE
    copula_dist=BiCopName(copula_dist)
    parameters=c("theta","zeta")
  } else {
    stop("ERROR: COPULA DIST LINK FUNCTIONS NOT YET IMPLEMENTED: DEFINE MANUALLY WITH copula_link ARGUMENT.")
  }

  if(two_par_cop) {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")} else {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr")}

  return_list=list()
  return_list[["copula_link"]]=copula_link
  return_list[["copula_dist"]]=copula_dist
  return_list[["parameters"]]=parameters

  return(return_list)
}
#' @export
plotDist <- function (dataset,dist) {

  num_margins=length(unique(dataset[,"time"]))

  margin_data=list()
  margin_unif=list()
  margin_fit=list()

  for (i in 1:num_margins) {
    margin_data[[i]]<-dataset[dataset[,"time"]==i,"response"]
    margin_fit[[i]]<-gamlss(margin_data[[i]]~1,family=dist)
    margin_unif[[i]]<-(margin_fit[[i]]$residuals)
  }

  ##plot.new()
  #par(mfrow=c(1,num_margins))

  #for (i in 1:num_margins) {histDist(margin_data[[i]],family=dist,xlab=TeX(paste("$Y_",i,"$")),main=paste("Histogram of margin",i,"and fitted",dist))}
  #invisible(readline(prompt="Press [enter] to continue"))

  plots=list()

  z=1
  for (i in 1:(num_margins)) {
    for (j in 1:(num_margins)) {
      if(i==j) {
        input_data=data.frame(margin_data[[i]])
        colnames(input_data)<-"X1"
        p <- ggplot(input_data, aes(x=X1)) +
          geom_histogram() +
          labs(x = TeX(paste("$Y_",i,"$")))
      }
      if(i!=j) {
        input_data=data.frame(cbind(margin_unif[[i]],margin_unif[[j]]))
        p=ggplot(data=input_data,aes(x=X1,y=X2)) +
          #geom_point(size=0.25,color="black") +
          geom_density_2d(contour_var="density",bins=10,color="black") +
          scale_fill_brewer() +
          labs(x = TeX(paste("$Y_",i,"$")), y=TeX(paste("$Y_",j,"$")),fill="density")
      }

      plots[[z]]=p
      z=z+1
    }
  }
  ggarrange(plotlist=plots,ncol=num_margins,nrow=num_margins)

}
#' @export
create_longitudinal_dataset <- function(response,covariates,labels=NA) {
  num_time_points=ncol(response)
  if(num_time_points <=1) {print('Not enough time points')}

  dataset<-matrix(data=NA,ncol=2+length(covariates),nrow=0)
  subject<-as.factor(seq(1:nrow(response)))

  for (t in 1:ncol(response)) {

    dataset_temp<-cbind(subject,t,response[,t])

    for (i in 1:length(covariates)) {
      if (ncol(covariates[[i]]) == 1 ) {
        covariate_for_time=covariates[[i]]
      } else {
        covariate_for_time=covariates[[i]][,t]
      }
      dataset_temp<-cbind(dataset_temp,covariate_for_time)
    }

    ###Add dataset temp to full table
    dataset <- rbind(dataset,dataset_temp)
  }

  if(!all(is.na(labels))) {
    colnames(dataset) <- labels
  }

  dataset=dataset[order(dataset$time,dataset$subject),] ###NOTE THIS WILL BREAK GLMM
  rownames(dataset)=1:nrow(dataset)

  return(dataset)
}

#' @export
loadDataset <- function(simOption=5,plot_dist=FALSE,n=100,d=3,copula_dist=NA, margin_dist,copula.link=NA,par.copula,par.margin,covariates_input=NA) {

  if (simOption==1) {
    load("Data/rand_mvt.rds")
    head(rand_mvt)

    # Basic data setup
    response = rand_mvt[,4:18]#[,4:18](4+2) ####Currently limiting to just 5 margins for simplicity
    covariates=list()
    covariates[[1]] = as.data.frame(rand_mvt[,19]) #Age 19:33 - changed to age at start to avoid correlation with time
    covariates[[2]] = as.data.frame(rand_mvt[,34:48]) #Time 34:48
    covariates[[3]] = as.data.frame(rand_mvt[,3]) #Gender

    # Setup data as longitudinal file
    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
  }
  else if (simOption==2) {

    # set up D-vine copula model with mixed pair-copulas
    d <- 3
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(2, 2, 0)
    par <- c(logit_inv(.8), logit_inv(.8), logit_inv(.8))
    par2 <- c(log_2plus_inv(2.1),log_2plus_inv(2.1),log_2plus_inv(2.1))

    # transform to R-vine matrix notation
    RVM <- D2RVine(order, family, par, par2)
    contour(RVM)

    t=d
    copsim=RVineSim(n*t,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {
      margin[covariates[[3]]==0,i]=qZISICHEL(copsim[,i],mu=exp(0.3+0.2*i),sigma=exp(0.3+0.2*i),nu=-0.8,tau=0.05)[covariates[[3]]==0]#Update to i*mu/sigma as needed
      margin[covariates[[3]]==1,i]=qZISICHEL(copsim[,i],mu=exp(0.3+0.2*i+0.1),sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.05)[covariates[[3]]==1]#Update to i*mu/sigma as needed
    }

    response = as.data.frame(margin)

    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

  }
  else if (simOption==3) {

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==d-1){
      par=c(copula.link$theta.linkinv(par.copula),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(copula.link$theta.linkinv(par.copula[1:(length(par.copula)/2)]), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(copula.link$zeta.linkinv(par.copula[(length(par.copula)/2+1):(length(par.copula))]),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- D2RVine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=RVineSim(n,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }
  else if (simOption==4) {

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==d-1){
      par=c(rep(par.copula,d-1),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(rep(par.copula["theta"],d-1), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"],d-1),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- D2RVine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=RVineSim(n,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }
  else if (simOption==5) {

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==d-1){
      par=c(rep(par.copula,d-1),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(rep(par.copula["theta"],d-1), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"],d-1),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- D2RVine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=RVineSim(n,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }
  else if (simOption==6) {

    t=d
    margin_sim=matrix(0,ncol=d,nrow=n)

    for (i in 1:d) {
      margin_sim[,i]=rnorm(n,1,3)
    }

    W=rnorm(n,0,3)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin_sim_out=matrix(0,ncol=d,nrow=n)
    for (i in 1:d) {
      margin_sim_out[,i]=margin_sim[,i]+W
    }

    response=as.data.frame(margin_sim_out)

  }
  else if (simOption==7) {

    t=d
    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))

    #par.copula=c(.3); names(par.copula)=c("theta")
    theta_intercept=unlist(par.copula["theta"])
    theta_out=theta_intercept+matrix(rep(covariates_input$theta.time*1:(d-1),n),ncol=d-1,byrow=TRUE) + 
    matrix(rep(as.matrix(covariates_input$theta.age*((covariates[[1]]-50)/100)^2),d-1),ncol=d-1)

    theta_inv=copula_input$copula_link$theta.linkinv(theta_out)

    if(length(par.copula)==d-1){
      par=c((par.copula),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(par.copula[1:(length(par.copula)/2)], rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(par.copula[(length(par.copula)/2+1):(length(par.copula))],rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notatio
    
    RVM=list()
    
    for (i in 1:n) {
      RVM[[i]] = D2RVine(order, c(rep(copula.family,length(theta_inv[i,])),rep(0,dd-(length(theta_inv[i,])))), par=c(theta_inv[i,],rep(0,dd-(length(theta_inv[i,])))), par2=c(theta_inv[i,],rep(0,dd-(length(theta_inv[i,])))))
    }
    #RVM <- D2RVine(order, rep(family[1],nrow(theta_inv)), theta_inv, theta_inv*0)
    #contour(RVM)

    copsim=RVineSim(n,RVM)


    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }
  else if (simOption==8) {
    #Multivariate Gamma

    U=margin_sim=matrix(0,ncol=d,nrow=n)

    a=.25;b=1.75;mu=rep(1,d)
    W=rbeta(n,shape1=a,shape2=b)

    for (i in 1:d) {
      U[,i]=rgamma(n,shape=a+b,rate=1/mu[i])
      margin_sim[,i]=U[,i]*W
    }

    #Fake covariates
    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=d,nrow=n))*(1:d)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    response=margin_sim
  } else if (simOption==9) { ########TIME VARIANT MU

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==1){
      par=c(rep(par.copula,d-1),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(rep(par.copula["theta"],d-1), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"],d-1),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- D2RVine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=RVineSim(n,RVM)
    

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  } else if (simOption==9) { ########TIME VARIANT MU

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==1){
      par=c(rep(par.copula,d-1),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(rep(par.copula["theta"],d-1), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"],d-1),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- D2RVine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=RVineSim(n,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }  else if (simOption==10) { ########TIME VARIANT SIGMA AND MU

    print("WARNING: SIMULATION IS IMPLEMENTED WITHOUT LINK FUNCTIONS FOR COVARIATES SO MAY NOT BE APPROPRIATE FOR ALL PARAMETER RANGES, ERRORS LIKELY")

    t=d

    # Setup covariates from covariates_input

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,2),0)) #Gender

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    mu_out=par.margin[1]+matrix(rep(covariates_input$mu.time*1:(d),n),ncol=d,byrow=TRUE) + 
      matrix(rep(as.matrix(covariates_input$mu.age*((covariates[[1]]-50)/100)^2),d),ncol=d) + 
      matrix(rep(as.matrix(covariates_input$mu.gender*(covariates[[3]])),d),ncol=d)
    sigma_out=exp(par.margin[2]+matrix(rep(covariates_input$sigma.time*1:(d),n),ncol=d,byrow=TRUE) + 
      matrix(rep(as.matrix(covariates_input$sigma.age*((covariates[[1]]-50)/100)^2),d),ncol=d)) +
      matrix(rep(as.matrix(covariates_input$sigma.gender*(covariates[[3]])),d),ncol=d)
    nu_out=par.margin[3]+matrix(rep(covariates_input$nu.time*1:(d),n),ncol=d,byrow=TRUE) + 
      matrix(rep(as.matrix(covariates_input$nu.age*((covariates[[1]]-50)/100)^2),d),ncol=d) + 
      matrix(rep(as.matrix(covariates_input$nu.gender*(covariates[[3]])),d),ncol=d)
    tau_out=(par.margin[4]+matrix(rep(covariates_input$tau.time*1:(d),n),ncol=d,byrow=TRUE) +
      matrix(rep(as.matrix(covariates_input$tau.age*((covariates[[1]]-50)/100)^2),d),ncol=d)) +
      matrix(rep(as.matrix(covariates_input$tau.gender*(covariates[[3]])),d),ncol=d)
    theta_out=par.copula[1]+matrix(rep(covariates_input$theta.time*1:(d-1),n),ncol=d-1,byrow=TRUE) + 
      matrix(rep(as.matrix(covariates_input$theta.age*((covariates[[1]]-50)/100)^2),d-1),ncol=d-1) +
      matrix(rep(as.matrix(covariates_input$theta.gender*(covariates[[3]])),d-1),ncol=d-1)
    zeta_out=par.copula[2]+matrix(rep(covariates_input$zeta.time*1:(d-1),n),ncol=d-1,byrow=TRUE) + 
      matrix(rep(as.matrix(covariates_input$zeta.age*((covariates[[1]]-50)/100)^2),d-1),ncol=d-1) + 
      matrix(rep(as.matrix(covariates_input$zeta.gender*(covariates[[3]])),d-1),ncol=d-1)
    
    #Print the ranges for all the variables in one simple output
    print(paste("MU RANGE: ", round(range(mu_out),2), " | SIGMA RANGE: ", round(range(sigma_out),2), " | NU RANGE: ", round(range(nu_out),2), " | TAU RANGE: ", round(range(tau_out),2), " | THETA RANGE: ", round(range(theta_out),2), " | ZETA RANGE: ", round(range(zeta_out),2)))

    #Define margin distribution
    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))

    # OK now for each row in theta_out and zeta_out, we need to create a new RVM and simulate from it, then apply the qFUN with the appropriate parameters to get the margin values for that row. This is going to be computationally intensive but should work.

    # row-specific copula simulation from theta_out / zeta_out
    copsim <- matrix(NA_real_, nrow = n, ncol = d)
    for (r in 1:n) {
      par_r  <- c(as.numeric(theta_out[r, ]), rep(0, dd - (d - 1)))
      par2_r <- c(if (length(par.copula) > 1 && all(is.finite(zeta_out[r, ]))) as.numeric(zeta_out[r, ]) else rep(0, d - 1),
                  rep(0, dd - (d - 1)))
      copsim[r, ] <- as.numeric(RVineSim(1, D2RVine(order, family, par_r, par2_r)))
    }

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=mu_out[,i],sigma=sigma_out[,i],nu=nu_out[,i],tau=tau_out[,i])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }

  dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

  if(plot_dist==TRUE) {plotDist(dataset,margin_dist)}

  return(dataset)
}
#' @export
bvt_norm_true_SE_B0_Bt <- function(sigma_x,sigma_y,rho,n,d){

  #sigma_x=2
  #sigma_y=2
  #rho=.75

  #((1)/(1-rho^2))
  #(1/(sigma_x^2))
  #(1/(sigma_y^2))
  #(1/(sigma_x))
  #(1/(sigma_y))
  #(rho/(sigma_x*sigma_y))

  hessian=((-1)/(1-(rho^2)))*matrix(c( (1/(sigma_x^2)) +(1/(sigma_y^2)) - 2*(rho/(sigma_x*sigma_y)),
                                       ((rho/(sigma_x*sigma_y))-(1/(sigma_y^2))),
                                       ((rho/(sigma_x*sigma_y))-(1/(sigma_y^2))),
                                       (1/(sigma_y^2))),nrow=2)

  vcov_matrix=-solve(hessian)

  true_SE=(diag(vcov_matrix))
  names(true_SE)=c("B0","Bt")

  return(true_SE)
}

########## ARCHIVE ###########

#Given a parameter vector starting values par = (mu,sigma,nu,tau,theta,zeta), return best fit parameters
#' @export
optim_outer <- function(par,dataset,margin_dist,copula_dist,
                        step_size=0.1,verbose=TRUE,use_dlcopdpar=TRUE) {

  #print("THIS FUNCTION ASSUMES RESPONSE IS ORDERED AS TIME, SUBJECT | PAR INPUT MUST BE NAMED")

  copula_input=get_copula_dist(copula_dist)
  copula_number=copula_input$copula_dist
  copula_link=copula_input$copula_link

  num_margins=length(unique(dataset$time))
  margin_names=unique(dataset$time)
  response=dataset$response

  #Set up parameter vector so names are consistent with the distributions

  if(all(is.null(names(par))|is.na(names(par)))) {print("ERROR: par vector must be named"); break}
  margin_par=par[names(par)%in%c("mu","sigma","nu","tau")]
  copula_par=par[!names(par)%in%c("mu","sigma","nu","tau")]

  ##### Calculate all relevant derivatives / CG method with first and second derivatives

  ### Calculate margin derivatives w.r.t. margin parameters

  #Get names for margin derivatives from margin_dist
  margin_deriv_names=names(margin_dist)[grepl("dld",names(margin_dist))|grepl("d2ld",names(margin_dist))]

  #Get link transforms (eta) and derivatives w.r.t to link for parameters
  par_eta_dr=par_eta=par*0
  for (par_name in names(par)) {
    if(par_name %in% names(margin_par)) {
      par_eta[par_name]=margin_dist[[paste(par_name,".linkfun",sep="")]](par[par_name])
      par_eta_dr[par_name]=margin_dist[[paste(par_name,".dr",sep="")]](par_eta[par_name])
    }
    if(par_name %in% names(copula_par)) {
      par_eta[par_name]=copula_link[[paste(par_name,".linkfun",sep="")]](par[par_name])
      par_eta_dr[par_name]=copula_link[[paste(par_name,".dr",sep="")]](par_eta[par_name])
    }
  }

  #Setup input matrix of response and parameters
  margin_deriv_input=list()
  margin_deriv_input[["y"]]=response
  margin_deriv_input[["q"]]=response
  margin_deriv_input[["x"]]=response
  for (par_name in names(margin_par)) {
    margin_deriv_input[[par_name]]=rep(margin_par[par_name],length(response))
  }

  #Calculate all derivatives
  margin_deriv=list()
  for (deriv_name in margin_deriv_names) {
    FUN=margin_dist[[deriv_name]]
    FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
    margin_deriv[[deriv_name]]=do.call(FUN,args=margin_deriv_input[FUN_args])
    margin_deriv[[deriv_name]][!is.finite(margin_deriv[[deriv_name]])]=0
  }

  margin_pFUN=eval(parse( text=paste("p",margin_dist$family[1],sep="") ))
  FUN=margin_pFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_p=do.call(FUN,args=margin_deriv_input[FUN_args])

  margin_dFUN=eval(parse( text=paste("d",margin_dist$family[1],sep="") ))
  FUN=margin_dFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_d=do.call(FUN,args=margin_deriv_input[FUN_args])

  ### Calculate copula derivatives w.r.t. copula parameters

  #First calculate margin F(x1), F(x2) as inputs to copula

  Fx_1_2=matrix(ncol=2,nrow=0)
  order_copula=matrix(ncol=4,nrow=0)
  for (i in 1:(num_margins-1)) {
    Fx_1_2=rbind(Fx_1_2,cbind(margin_p[dataset$time == margin_names[i]],margin_p[dataset$time == margin_names[i+1]]))
    order_copula=rbind(order_copula,cbind(dataset[dataset$time == margin_names[i],c("time","subject")],dataset[dataset$time == margin_names[i+1],c("time","subject")]))
  }
  names(order_copula)=c("time1","subject1","time2","subject2")

  par1=copula_par["theta"]
  if(is.na(copula_par["zeta"])) {par2=0} else {par2=copula_par["zeta"]}

  #Handling extreme values
  Fx_1_2[Fx_1_2>1]=1;Fx_1_2[Fx_1_2<0]=0

  if(copula_number==3) {
    if(par1>28){par1=28}
  }

  copula_d=BiCopPDF(  Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2)
  dldth=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par",log=TRUE)
  dcdth=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par",log=FALSE)
  d2cdth=BiCopDeriv2( Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par")
  d2ldth2=(1/(copula_d^2))*(copula_d*d2cdth-dcdth^2)
  if(!is.na(copula_par["zeta"])) {
    dldz=BiCopDeriv(    Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par2",log=TRUE)
    dcdz=BiCopDeriv(    Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par2",log=FALSE)
    d2cdz=BiCopDeriv2(  Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par2")
    d2ldz2=(1/(copula_d^2))*(copula_d*d2cdz-dcdz^2)

    d2cdthdz=BiCopDeriv2(  Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par1par2")
    d2ldthdz=(d2cdthdz*copula_d-dcdth*dcdz)/(copula_d^2)
  }
  dcdu1=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u1",log=FALSE)
  dcdu2=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u2",log=FALSE)

  d2cdu12=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u1",log=FALSE)
  d2cdu22=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u2",log=FALSE)

  d2ldth2[!is.finite(d2ldth2)]=0

  ### Calculate copula derivatives w.r.t margin parameters

  #Extract margin calculations for F(x), f(x), response and derivatives at time 1 and time 2, join to copula values for time 1 and time 2
  margin_deriv_1=margin_deriv_2=margin_deriv_2cross=matrix(ncol=length(margin_par),nrow=length(response))
  for (i in 1:length(margin_par)) {
    margin_deriv_1[,i]=margin_deriv[grepl("dld",names(margin_deriv))][[i]]
    margin_deriv_2[,i]=margin_deriv[grepl("d2ld",names(margin_deriv))&endsWith(names(margin_deriv),"2")][[i]]
  }
  colnames(margin_deriv_1)=paste("dld",names(margin_par),sep="")
  colnames(margin_deriv_2)=paste(paste("d2ld",names(margin_par),sep=""),"2",sep="")

  #colnames(margin_deriv_2)=paste("d2ld",names(margin_par),sep="")

  order_margin=dataset[,c("time","subject")]
  margin_components=cbind(order_margin,response,margin_p,margin_d,margin_deriv_1,margin_deriv_2)
  margin_components_Ft_plus=margin_components
  margin_components_Ft_plus$time=margin_components_Ft_plus$time-1
  margin_plus=merge(margin_components,margin_components_Ft_plus,by=c("time","subject"),all.x=TRUE)

  copula_components=cbind(order_copula,dcdu1,dcdu2,copula_d,d2cdu12,d2cdu22)
  copula_merged=merge(copula_components,margin_plus,by.x=c("time1","subject1"),by.y=c("time","subject"),all.x=TRUE)

  #Calculate copula derivative with respect to marginal parameters
  input=copula_merged
  dlcopdpar=matrix(0,nrow=nrow(input),ncol=length(margin_par))
  d2lcopdpar2=matrix(0,nrow=nrow(input),ncol=length(margin_par))

  i=1
  for (par_name in names(margin_par)) {

    #Take parameters from input for clarity
    dc_tplus_du_t=input[,"dcdu1"]
    dc_tplus_du_tplus=input[,"dcdu2"]
    l_t=input[,paste(paste("dld",par_name,sep=""),".x",sep="")]
    l_t_plus=input[,paste(paste("dld",par_name,sep=""),".y",sep="")]
    x_t=input[,"response.x"]
    x_t_plus=input[,"response.y"]
    f_t=input[,"margin_d.x"]
    f_t_plus=input[,"margin_d.y"]
    du_t_dmu=x_t*f_t*l_t
    du_t_plus_dmu=x_t_plus*f_t_plus*l_t_plus
    c_tplus=input[,"copula_d"]

    du_t_dmu=x_t*f_t*l_t
    du_t_plus_dmu=x_t_plus*f_t_plus*l_t_plus

    dc_plus_dt_dmu=dc_tplus_du_t * du_t_dmu
    dc_plus_dt_plus_dmu=dc_tplus_du_tplus * du_t_plus_dmu
    dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)]=0
    dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)]=0
    dcdmu_tplus=((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
    dcdmu_tplus[is.nan(dcdmu_tplus)|is.na(dcdmu_tplus)]=0

    dlcopdpar[,i]=dcdmu_tplus


    #######NOW FOR SECOND DERIVATIVE OF COPULA TERM

    l2_t=input[,paste(paste(paste("d2ld",par_name,sep=""),"2",sep=""),".x",sep="")]
    l2_tplus=input[,paste(paste(paste("d2ld",par_name,sep=""),"2",sep=""),".y",sep="")]

    df_t_dmu=f_t*l_t
    df_t_plus_dmu=f_t_plus*l_t_plus

    d2f_t_dmu=df_t_dmu*l_t + f_t*l2_t
    d2f_t_plus_dmu=df_t_plus_dmu*l_t_plus + f_t_plus*l2_tplus

    d2u_t_dmu2=x_t*d2f_t_dmu
    d2u_t_plus_dmu2=x_t_plus*d2f_t_plus_dmu

    d2cdu_t2=input[,"d2cdu12"]
    d2cdu_t_plus2=input[,"d2cdu22"]
    d2cdu_t2[is.nan(d2cdu_t2)]=0
    d2cdu_t_plus2[is.nan(d2cdu_t_plus2)]=0

    d2cdmu2=d2cdu_t2*du_t_dmu^2 + dc_tplus_du_t * d2u_t_dmu2 + d2cdu_t_plus2*du_t_plus_dmu^2 + dc_tplus_du_tplus * d2u_t_plus_dmu2

    d2lcdmu2=as.matrix((d2cdmu2*c_tplus-(dcdmu_tplus^2))/(c_tplus^2))

    d2lcopdpar2[,i]=d2lcdmu2
    #num_deriv=margin_copula_merged_2[,"num_dlcopdpar_ordered.Ft"]
    #num_deriv_nolog=margin_copula_merged_2[,"num_dlcopdpar_nolog_ordered.Ft"]

    i=i+1
  }
  colnames(dlcopdpar)=paste("dlcopd",names(margin_par),sep="")
  colnames(d2lcopdpar2)=paste(paste("d2lcd",names(margin_par),sep=""),"2",sep="")

  dlcopdpar[!is.finite(dlcopdpar)]=0
  d2lcopdpar2[!is.finite(d2lcopdpar2)]=0

  #### Define score and hessian

  score=par*0
  hessian=matrix(0,nrow=length(par),ncol=length(par))
  colnames(hessian)=names(par);rownames(hessian)=names(par)
  names(score)=names(par)

  margin_deriv_sum=vector()
  for (i in 1:length(margin_deriv)) {
    margin_deriv[[i]][!is.finite(margin_deriv[[i]])]=0
    margin_deriv_sum[i]=sum(margin_deriv[[i]])
  }
  names(margin_deriv_sum)=names(margin_deriv)

  margin_d1=margin_deriv_sum[grepl("dld",names(margin_deriv))]
  margin_d2=margin_deriv_sum[grepl("d2ld",names(margin_deriv))&endsWith(names(margin_deriv),"2")]
  margin_d2d=margin_deriv_sum[grepl("d2ld",names(margin_deriv))&!endsWith(names(margin_deriv),"2")]

  if(is.na(copula_par["zeta"])) {
    copula_d1=sum(dldth)
    copula_d2=sum(d2ldth2)
  } else {
    copula_d1=colSums(cbind(dldth,dldz))
    copula_d2=colSums(cbind(d2ldth2,d2ldz2))
  }
  margin_d1_dlcopdpar=margin_d1+if(use_dlcopdpar==TRUE){ colSums(dlcopdpar)} else {colSums(dlcopdpar)*0}
  margin_d2_dlcopdpar=margin_d2+if(use_dlcopdpar==TRUE){ colSums(d2lcopdpar2)*0} else {colSums(d2lcopdpar2)*0}
  score=c(margin_d1_dlcopdpar,copula_d1)

  ###CALCULATING HESSIAN USING D2
  diag(hessian)=c(margin_d2_dlcopdpar,copula_d2)
  hessian[1:length(margin_par),1:length(margin_par)][upper.tri(hessian[1:length(margin_par),1:length(margin_par)])]=margin_d2d
  hessian[1:length(margin_par),1:length(margin_par)][lower.tri(hessian[1:length(margin_par),1:length(margin_par)])]=margin_d2d

  #Why isn't d2 for copula negative?
  copula_hess=hessian[(length(margin_par)+1):(length(margin_par)+length(copula_par)),(length(margin_par)+1):(length(margin_par)+length(copula_par))]
  if(!is.na(copula_par["zeta"])) {
    copula_hess[upper.tri(copula_hess)]=sum(d2ldthdz)
    copula_hess[lower.tri(copula_hess)]=sum(d2ldthdz)
  }
  hessian[(length(margin_par)+1):(length(margin_par)+length(copula_par)),(length(margin_par)+1):(length(margin_par)+length(copula_par))]=copula_hess


  ###STILL NEED TO CALCULATE d2 for marginal parameters with respect to copula likelihood and add to hessian values

  #par_end=par-(solve(-hessian)%*%(score))

  #weights_eta=diag((1/(score_eta^2)))
  #weights=-diag(score*score)

  #score=score
  #weights=-solve(hessian)

  weights_eta=-solve(hessian*par_eta_dr*par_eta_dr)

  #weights_eta=diag(1/(score*score*par_eta_dr*par_eta_dr))
  score_eta=score*par_eta_dr
  #par_end=par*(1-step_size) + step_size*(par+par_change)

  par_end=par*0
  names(par_end)=names(par_eta_end)=names(par)
  #Get end paraemters re-transformed
  #for (par_name in names(par)) {
  #  if(par_name %in% names(margin_par)) {
  #    par_end[par_name]=margin_dist[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
  #  }
  #  if(par_name %in% names(copula_par)) {
  #    par_end[par_name]=copula_link[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
  #  }
  #}
  ###If calculating for eta
  for (par_name in names(par)) {
    if(par_name %in% names(margin_par)) {
      par_end[par_name]=margin_dist[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
    }
    if(par_name %in% names(copula_par)) {
      par_end[par_name]=copula_link[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
    }
  }


  sum_log_margin_p=sum(log(margin_d)[is.finite(log(margin_d))])
  sum_log_copula_d=sum(log(copula_d)[is.finite(log(copula_d))])

  log_lik=c(sum_log_copula_d,sum_log_margin_p,sum_log_copula_d+sum_log_margin_p)
  names(log_lik)=c("copula","margin","joint")

  if(verbose==TRUE) {
    print("Start Parameters")
    print(par)
    print("End Parameters:")
    print(par_end)
    print("Score:")
    print(score)
    print("Hessian:")
    print(hessian)
    print("Weights:")
    print(weights_eta)

    print(log_lik)
  }

  return(list(score=score,hessian=hessian,par_end=par_end,par_eta_end=par_eta_end,par_start=par,log_lik=log_lik))
}
