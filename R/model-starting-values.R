#' @keywords internal

#' @noRd

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


.select_t_copula_zeta_start <- function(dataset, margin_dist, copula_dist, margin_par, theta_start) {

  # Grid is ordered low-to-high: return the first candidate that yields a finite

  # joint log-likelihood. Starting as low as possible avoids the optimizer being

  # trapped at high df (Gaussian limit) where the link-scale step sizes collapse.

  zeta_grid <- c(2.05, 2.2, 2.5, 3, 4, 5, 8, 12, 20, 35)

  fallback_zeta <- 3


  param_names <- c(names(margin_dist$parameters), get_copula_dist(copula_dist)$parameters)

  mm_stub <- as.list(setNames(rep(1, length(param_names)), param_names))

  pair_cache <- build_copula_pair_cache(dataset$response, dataset$time, dataset$subject)


  base_eta_inv <- c(as.list(margin_par), list(theta = as.numeric(theta_start)[1]))


  for (candidate_zeta in zeta_grid) {

    eta_inv <- base_eta_inv

    eta_inv$zeta <- candidate_zeta


    candidate_fit <- tryCatch(

      calc_likelihood_minimal(

        eta_inv = eta_inv,

        mm = mm_stub,

        margin_dist = margin_dist,

        copula_dist = copula_dist,

        calc_d2 = FALSE,

        response = dataset$response,

        response_margin = dataset$time,

        response_subject = dataset$subject,

        pair_cache = pair_cache

      ),

      error = function(e) NULL

    )


    if (!is.null(candidate_fit) && is.finite(candidate_fit$log_lik["joint"])) {

      return(candidate_zeta)

    }

  }


  fallback_zeta

}


#' @keywords internal

#' @noRd

get_starting_values = function(copula_dist,margin_dist,dataset,eta_transform=FALSE) {


  margin_dist <- .normalise_margin_dist_links(margin_dist)


  margin_names=unique(dataset$time)

  num_margins=length(margin_names)

  finite_response <- dataset$response[is.finite(dataset$response)]

  margin_par_already_eta <- FALSE

  moment_skewness <- function(x) {

    x <- x[is.finite(x)]

    if (length(x) < 3) return(0)

    s <- stats::sd(x)

    if (!is.finite(s) || s <= 0) return(0)

    mean(((x - mean(x)) / s)^3)

  }

  moment_kurtosis <- function(x) {

    x <- x[is.finite(x)]

    if (length(x) < 4) return(3)

    s <- stats::sd(x)

    if (!is.finite(s) || s <= 0) return(3)

    mean(((x - mean(x)) / s)^4)

  }


  tau_start=cor(dataset[dataset$time%in%(margin_names[1:(num_margins-1)]),"response"]

                ,dataset[dataset$time%in%(margin_names[2:(num_margins)]),"response"],method="kendall",use="complete.obs")

  if(!is.finite(tau_start)) {

    warning("Non-finite Kendall tau in get_starting_values(); using tau = 0 for copula initialisation.")

    tau_start=0

  }

  tau_start=max(min(tau_start,0.9999),-0.9999)


  copula_spec=get_copula_dist(copula_dist)

  theta_start=.copula_tau_to_par(

    family=copula_dist,

    tau=tau_start

  )


  if(margin_dist$family[1]=="GA" | margin_dist$family[1]=="EXP") {

    margin_par=c(

      mean(finite_response)

      , stats::sd(finite_response)/mean(finite_response)

      , moment_skewness(finite_response)

      , moment_kurtosis(finite_response)

    )

  } else if (margin_dist$family[1]=="NO") {

    margin_par=c(

      mean(finite_response)

      , stats::sd(finite_response)

    )

  } else if (margin_dist$family[1]=="PO") {

    margin_par=c(

      mean(finite_response)

    )

  } else if (margin_dist$family[1]=="NBI") {

    margin_par=c(

      mean(finite_response),

      stats::sd(finite_response)/mean(finite_response)

    )

  } else {

    cat("Fitting initial GAMLSS model for margin to obtain starting values...\n")

    # Deliberately low-iteration startup fit; silence expected convergence warnings.

    start_fit <- tryCatch({

      fit <- NULL

      invisible(utils::capture.output({

        fit <- suppressWarnings(suppressMessages(

          gamlss(

            dataset$response ~ 1,

            family = margin_dist,

            control = gamlss::gamlss.control(n.cyc = 5, trace = FALSE)

          )

        ))

      }))

      fit

    }, error = function(e) NULL)

    margin_par <- vapply(names(margin_dist$parameters), function(parameter) {

      cf <- tryCatch(stats::coef(start_fit, what = parameter), error = function(e) numeric(0))

      if (length(cf) > 0L && is.finite(cf[[1L]])) {

        return(as.numeric(cf[[1L]]))

      }

      qfun <- get(paste0("q", margin_dist$family[1]), envir = asNamespace("gamlss.dist"), inherits = FALSE)

      value <- tryCatch(eval(formals(qfun)[[parameter]], envir = baseenv()), error = function(e) NA_real_)

      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {

        value <- switch(parameter,

          mu = if (.is_discrete_margin(margin_dist)) 3 else mean(finite_response),

          sigma = stats::sd(finite_response) / max(abs(mean(finite_response)), 1e-8),

          nu = 0.5,

          tau = 2,

          0

        )

      }

      linkfun <- margin_dist[[paste0(parameter, ".linkfun")]]

      if (is.null(linkfun)) {

        return(as.numeric(value))

      }

      as.numeric(linkfun(value))[1L]

    }, numeric(1))

    names(margin_par)=names(margin_dist$parameters)

    margin_par_already_eta <- TRUE

  }


  names(margin_par)=names(margin_dist$parameters)

  margin_par=margin_par[!is.na(names(margin_par))]


  if("zeta" %in% copula_spec$parameters) {

    # .copula_tau_to_par() returns only theta for t-copula; select zeta by a small grid search.

    zeta_start <- .select_t_copula_zeta_start(

      dataset = dataset,

      margin_dist = margin_dist,

      copula_dist = copula_dist,

      margin_par = margin_par,

      theta_start = theta_start

    )

    cop_par=c(theta=as.numeric(theta_start)[1], zeta=as.numeric(zeta_start))

  } else {

    cop_par=c(theta=as.numeric(theta_start)[1])

  }


  if(eta_transform==TRUE) {

    margin_par_eta=margin_par

    cop_par_eta=cop_par


    if(!isTRUE(margin_par_already_eta)) {

      for (par_name in names(margin_par)) {

        FUN = eval(parse(text=paste(paste(paste("margin_dist$",par_name,sep=""),"linkfun",sep="."))))

        margin_par_eta[par_name]=FUN(margin_par[par_name])

      }

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

#' @keywords internal

#' @noRd

par_to_eta = function(par,copula_dist,margin_dist) {


  margin_dist <- .normalise_margin_dist_links(margin_dist)


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


#' @keywords internal

#' @noRd

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





#' @keywords internal

#' @noRd

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

  fixed_unlinked_values <- attr(margin_dist, "fixed_unlinked_values")

  if (length(fixed_unlinked_values) > 0L) {

    for (par_name in names(fixed_unlinked_values)) {

      value <- fixed_unlinked_values[[par_name]]

      if (is.numeric(value) && length(value) == 1L && is.finite(value)) {

        margin_deriv_input[[par_name]] <- rep(value, length(response))

      }

    }

  }


  negative_response <- is.finite(response) & response < 0

  if (.is_discrete_margin(margin_dist) && any(negative_response)) {

    margin_p <- rep(NA_real_, length(response))

    margin_p[negative_response] <- 0

    valid_response <- !negative_response

    if (any(valid_response)) {

      call_input <- lapply(margin_deriv_input, function(value) {

        if (length(value) == length(response)) value[valid_response] else value

      })

      margin_pFUN <- get(

        paste("p", margin_dist$family[1], sep = ""),

        envir = asNamespace("gamlss.dist"),

        mode = "function",

        inherits = FALSE

      )

      FUN_args <- names(call_input)[names(call_input) %in% formalArgs(margin_pFUN)]

      margin_p[valid_response] <- tryCatch(

        do.call(margin_pFUN, args = call_input[FUN_args]),

        error = function(e) rep(NA_real_, sum(valid_response))

      )

    }

    return(margin_p)

  }


  margin_pFUN <- get(

    paste("p", margin_dist$family[1], sep = ""),

    envir = asNamespace("gamlss.dist"),

    mode = "function",

    inherits = FALSE

  )

  FUN=margin_pFUN

  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]

  margin_p <- tryCatch(

    do.call(FUN, args = margin_deriv_input[FUN_args]),

    error = function(e) rep(NA_real_, length(response))

  )


  return(margin_p)

}


#' @keywords internal

#' @noRd

get_copula_dist=function(copula_dist) {


  copula_dist <- .copula_family_code(copula_dist)


  if(copula_dist=="C") {

    copula_link=list(log,exp,dloginv=exp); two_par_cop=FALSE

    parameters=c("theta")

  }

  else if(copula_dist=="F") {

    copula_link=list(identity,identity,function(x) rep(1, length(x))); two_par_cop=FALSE

    parameters=c("theta")

  }

  else if(copula_dist=="J") {

    copula_link=list(log_1plus,log_1plus_inv,dlog_1plus_inv); two_par_cop=FALSE

    parameters=c("theta")

  }

  else if(copula_dist=="G") {

    copula_link=list(gumbel_linkfun,gumbel_linkinv,dgumbel_linkinv); two_par_cop=FALSE

    parameters=c("theta")

  }

  else if(copula_dist=="N") {

    copula_link=list(fisher_z,fisher_z_inv,dfisher_z_inv); two_par_cop=FALSE

    parameters=c("theta")

  } else if(copula_dist=="t") {

    copula_link=list(fisher_z,fisher_z_inv,dfisher_z_inv,log_2plus,log_2plus_inv,dlog_2plus_inv); two_par_cop=TRUE

    parameters=c("theta","zeta")

  } else {

    stop("ERROR: COPULA DIST LINK FUNCTIONS NOT YET IMPLEMENTED.")

  }


  if(two_par_cop) {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")} else {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr")}


  return_list=list()

  return_list[["copula_link"]]=copula_link

  return_list[["copula_dist"]]=copula_dist

  return_list[["parameters"]]=parameters


  return(return_list)

}


