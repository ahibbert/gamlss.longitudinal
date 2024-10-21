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
fit_jointreg=function(  dataset,
                        margin_dist,
                        copula_dist,
                        mu.formula = ("response ~ 1"),
                        sigma.formula = ("~ 1"),
                        nu.formula = ("~ 1"),
                        tau.formula = ("~ 1"),
                        theta.formula=("~1"),
                        zeta.formula=("~1"),
                        include_dlcopdpar=TRUE,
                        inner_stop_crit=.01,
                        outer_stop_crit=.0001,
                        start_step_size=.5,
                        step_adjustment=.5,
                        max_steps=5,
                        start_from=NA,
                        verbose=TRUE,
                        plot_results=TRUE,
                        true_val=NA,
                        method="RS",
                        max_outer_iter=100,
                        max_inner_iter=20,
                        use_Rcpp=FALSE
                      )
{

  #Setup model matrix from given formulas
  copula_link=get_copula_dist(copula_dist)$copula_link
  mm=suppressWarnings(create_model_matrices(
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    theta.formula,
    zeta.formula,
    margin.family=margin_dist,copula.family=copula_dist,copula.link=copula_link
  ))

  #Create vector of starting covariate values, currently starting at zero before first fit with the intercept as the mean

  if(all(is.na(start_from))) {
    par_eta=get_starting_values(copula_dist,margin_dist,dataset,eta_transform=TRUE)
  } else {
    par_eta=start_from
  }

  par_cov=as.numeric(vector())
  for (par_name in names(mm)) {
    par_cov_single=as.numeric(vector(length=length(colnames(mm[[par_name]]))))
    names(par_cov_single)=paste(par_name,colnames(mm[[par_name]]),sep=".")
    par_cov_single[1]=par_eta[par_name]
    if(length(par_cov_single)>1) {
      par_cov_single[2:length(par_cov_single)]=0
    }
    par_cov=c(par_cov,par_cov_single)
  }

  #Parameters used in optimisation loops
  first_outer_run=TRUE
  outer_log_lik_change=outer_start_log_lik=outer_end_log_lik=0
  log_lik_history=matrix(ncol=3,nrow=0)
  par_history=matrix(ncol=length(par_cov),nrow=0); colnames(par_history)=names(par_cov)
  outer_run_counter=1; outer_only_run_counter=1
  step_size=start_step_size

  #OUTER ITERATION (MAIN LOOP)
  while ((first_outer_run==TRUE | (abs(outer_log_lik_change)>outer_stop_crit)) & outer_only_run_counter < max_outer_iter) {

    cat(paste("\nOUTER ITERATION:",outer_only_run_counter))
    first_outer_run=TRUE

    # RUN INNER ITERATION FOR EACH PARAMETER
    for (par_name in names(mm)) {

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
        eta_out=calc_eta(par_cov,mm,margin_dist,copula_link)
        eta=eta_out$eta; eta_dr=eta_out$eta_dr; eta_inv=eta_out$eta_inv

        calc_lik_out=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2=FALSE,response=dataset$response,margin_names=unique(dataset$time))
        log_lik=calc_lik_out$log_lik; margin_d=calc_lik_out$margin_d; margin_p=calc_lik_out$margin_p; margin_deriv=calc_lik_out$margin_deriv; copula_d=calc_lik_out$copula_d; copula_p=calc_lik_out$copula_p; Fx_1_2=calc_lik_out$Fx_1_2;order_copula=calc_lik_out$order_copula

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
        copula_derivatives=calc_copula_derivatives(eta_inv, Fx_1_2, copula_dist)
        dldth=copula_derivatives$dldth; dcdth=copula_derivatives$dcdth; dcdu1=copula_derivatives$dcdu1; dcdu2=copula_derivatives$dcdu2
        if(!"zeta" %in% names(eta_inv)) {dldz=copula_derivatives$dldz; dcdz=copula_derivatives$dcdz}

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Copula Derivatives")

        ### Calculate copula derivatives w.r.t margin parameters
        if(!par_name %in% c("mu","sigma","nu","tau")) {
          if(!"zeta" %in% names(eta_inv)) {
            d1=dldth
            #d2=d2ldth2
          } else {
            d1=cbind(dldth,dldz)
            #d2=cbind(d2ldth2,d2ldz2)
          }
        } else {

          ### MARGIN LIKELIHOOD DERIVATIVES
          margin_deriv_subnames=c("m","d","v","t")
          names(margin_deriv_subnames)=c("mu","sigma","nu","tau")
          margin_par=names(mm)[names(mm) %in% c("mu","sigma","nu","tau")]
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
            nd_impact_F=calc_Fx_derivatives(eta_inv,mm,margin_dist)

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


        ### BACKFITTING STEP
        score=score_function_v2(eta=eta[[par_name]],dldpar=d1,d2ldpar=-(d1*d1),dpardeta=eta_dr[[par_name]])

        X=as.matrix(mm[[par_name]])
        W=diag(as.vector(score$w_k))
        z_k=score$z_k
        beta_start=par_cov[paste(paste(par_name,sep=" "),colnames(mm[[par_name]]),sep=".")]

        if(use_Rcpp==TRUE) {
          sourceCpp("test.cpp")
          t_x_w=eigenMapMatMult(t(X),W)
          t_x_w_x=eigenMapMatMult(t_x_w,X)

          inv_t_x_w_x=tryCatch( chol2inv(chol(t_x_w_x)),
                                error=function(e) {
                                                    return(solve(t_x_w_x))
                                                  }
          )


          update_no_z_k=eigenMapMatMult(inv_t_x_w_x,t_x_w)
          beta_update=eigenMapMatMult(update_no_z_k,z_k)
        } else {
          beta_update=as.vector(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z_k)
        }

        beta_change_inner=beta_update-beta_start
        beta=beta_start*(1-step_size) + (step_size)*(beta_update)
        names(beta)=paste(paste(par_name,sep=" "),colnames(mm[[par_name]]),sep=".")

        par_cov_new=par_cov
        par_cov_new[names(beta)]=beta

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Backfitting")

        eta_out=calc_eta(par_cov_new,mm,margin_dist,copula_link)

        eta=eta_out$eta; eta_dr=eta_out$eta_dr; eta_inv=eta_out$eta_inv
        par_cov=par_cov_new

        calc_lik_out_end=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2=FALSE,response=dataset$response,margin_names=unique(dataset$time))

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
  cat(paste("\nTotal time:",round(max(timer),2)))
  cat("\n\n")
  par_mat_out_temp=t(t((par_cov)))
  colnames(par_mat_out_temp) = c("estimate")
  print(par_mat_out_temp)
  cat("\n")
  cat("Model Selection Criteria:")
  cat("\n")

  p_cop=par_cov[grepl("theta",names(par_cov))|grepl("zeta",names(par_cov))]
  p_mar=par_cov[!(grepl("theta",names(par_cov))|grepl("zeta",names(par_cov)))]

  aics=rbind(t(calc_lik_out_end$log_lik),
             t(calc_lik_out_end$log_lik)-2*c(length(p_mar),length(p_cop),length(par_cov)),
             t(calc_lik_out_end$log_lik)-c(length(p_mar),length(p_cop),length(par_cov))*log(nrow(dataset)))

  rownames(aics)=c("LogLik","AIC","BIC")
  print(aics)

  cat("\n####################################\n")

  return_list=list(par_cov,log_lik_history,par_history,calc_lik_out_end,mm,margin_dist,copula_dist,include_dlcopdpar,dataset$response)
  names(return_list)=c("par","log_lik_history","par_history","calc_lik_out_end","model_matrix","margin_dist","copula_dist","include_dlcopdpar","response")
  class(return_list)="gamlss.longitudinal"
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

          #if(margin_dist$family[1]=="EXP") {
          #  du_t_dmu_approx=du_t_dmu
          #  du_t_plus_dmu_approx=du_t_plus_dmu
          #   du_t_dmu=-x_t*exp(-x_t*mu_t)
          #   du_t_plus_dmu=-x_t_plus*exp(-x_t_plus*mu_t_plus)
          #
          #   par(mfrow=c(2,2))
          #   plot(du_t_dmu_approx,du_t_dmu)
          #   plot(du_t_plus_dmu_approx,du_t_plus_dmu)
          #   plot(F_nd_t,du_t_dmu)
          #   plot(F_nd_t_plus,du_t_plus_dmu)
          # }

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

          d2cdmu2=d2cdu_t2*du_t_dmu^2 + dc_tplus_du_t * d2u_t_dmu2 + d2cdu_t_plus2*du_t_plus_dmu^2 + dc_tplus_du_tplus * d2u_t_plus_dmu2

          d2lcdmu2=as.matrix((d2cdmu2*c_tplus-(dcdmu_tplus^2))/(c_tplus^2))

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

      return(d2_cop)
    }


}

calc_copula_derivatives = function(eta_inv, Fx_1_2, copula_dist, calc_d2=FALSE, calc_d2_marginal=FALSE) {

  par1=eta_inv[["theta"]]

  if("zeta" %in% names(eta_inv)) {
    par2=eta_inv[["zeta"]]
  } else {
    par2=eta_inv[["theta"]]*0
  }

  if(copula_dist=="C") {
    par1[par1>=28]=27.9
  }

  copula_d=BiCopPDF(   Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2)

  if(!"zeta" %in% names(eta_inv)) {par2=eta_inv[["theta"]]*0} else {par2=eta_inv[["zeta"]]}
  dldth=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par",log=TRUE)
  dcdth=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par",log=FALSE)

  if(calc_d2==TRUE) {
    d2cdth=BiCopDeriv2( Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par")
    d2ldth2=(1/(copula_d^2))*(copula_d*d2cdth-dcdth^2)
  }

  if("zeta" %in% names(eta_inv)) {
    dldz=BiCopDeriv(    Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par2",log=TRUE)
    dcdz=BiCopDeriv(    Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par2",log=FALSE)

    if(calc_d2==TRUE) {
      d2cdz=BiCopDeriv2(  Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par2")
      d2ldz2=(1/(copula_d^2))*(copula_d*d2cdz-dcdz^2)

      d2cdthdz=BiCopDeriv2(  Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par1par2")
      d2ldthdz=(d2cdthdz*copula_d-dcdth*dcdz)/(copula_d^2)
    }

  }
  dcdu1=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="u1",log=FALSE)
  dcdu2=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="u2",log=FALSE)

  dldth[!is.finite(dldth)]=0; if(calc_d2==TRUE) {d2ldth2[!is.finite(d2ldth2)]=0  }

  if(calc_d2==TRUE) {
    d2cdu12=BiCopDeriv2( Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="u1")
    d2cdu22=BiCopDeriv2( Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="u2")
  }

  if("zeta" %in% names(eta_inv)) {
    dldz[!is.finite(dldz)]=0
    if(calc_d2==TRUE) {
      d2ldz2[!is.finite(d2ldz2)]=0
      d2ldthdz[!is.finite(d2ldthdz)]=0
    }
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

calc_Fx_derivatives = function(eta_inv, mm, margin_dist) {
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
      change_F[,i]=calc_F_x(eta_inv_adj,mm,margin_dist)#calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
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

calc_Fx2_derivatives = function(eta_inv, mm, margin_dist) {
  nd_impact=nd_impact_m=nd_impact_c=rep(0,length(names(eta_inv)))
  nd_impact_F=list() ###### WE DON"T NEED TO BE CALCULATING THIS FOR ALL PARAMETERS EVERY TIME

  margin_par_names=names(eta_inv)[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (eta_par_names_nd in margin_par_names) {
    adj_fac=.0001
    change=change_m=change_c=c(0,0)
    change_F=matrix(0,nrow=length(eta_inv[[1]]),ncol=3)
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      #change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["joint"]
      #change_m[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["marginal"]
      #change_c[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["copula"]
      change_F[,i]=calc_F_x(eta_inv_adj,mm,margin_dist)#calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
      #print(change_F[,i]==calc_F_x(eta_inv_adj,mm,margin_dist))
      i=i+1
    }
    change_F[,3]=calc_F_x(eta_inv,mm,margin_dist)
    #nd_impact[eta_par_names_nd]=(change[2]-change[1])/(2*adj_fac)
    #nd_impact_m[eta_par_names_nd]=(change_m[2]-change_m[1])/(2*adj_fac)
    #nd_impact_c[eta_par_names_nd]=(change_c[2]-change_c[1])/(2*adj_fac)
    nd_impact_F[[eta_par_names_nd]]=(change_F[,2]+change_F[,1]-2*change_F[,3])/(adj_fac^2)
  }
  return(nd_impact_F)
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
vcov.gamlss.longitudinal=function(object) {

  #object=w_dl
  #object=out[[i]][[j]]

  include_dlcopdpar=TRUE

  margin_names=unique(dataset$time)
  num_margins=length(margin_names)

  se_out=object$par*0;
  margin_dist=no_dl$margin_dist; copula_dist=object$copula_dist; copula_link=get_copula_dist(copula_dist)$copula_link
  mm=object$model_matrix

  eta_out=calc_eta(par_cov=object$par,mm=mm,margin_dist,copula_link)
  eta_inv=eta_out$eta_inv; eta_dr=eta_out$eta_dr; eta=eta_out$eta
  calc_lik_out=calc_likelihood_minimal(eta_inv=eta_inv,mm=object$model_matrix,margin_dist,copula_dist,calc_d2=TRUE,response=object$response,margin_names)
  Fx_1_2=calc_lik_out$Fx_1_2; margin_p=calc_lik_out$margin_p; margin_d=calc_lik_out$margin_d; copula_d=calc_lik_out$copula_d

  ###Calculate derivaties: margin and copula d1 and d2
  margin_derivatives=calc_lik_out$margin_deriv
  copula_derivatives=calc_copula_derivatives(eta_inv, Fx_1_2, copula_dist,calc_d2 = TRUE)
  nd_impact_F=calc_Fx_derivatives(eta_inv,mm,margin_dist)

  #Calculate numerical derivatives for F(x) - also used as a reference for overall likelihood derivative
  nd_impact_F=calc_Fx_derivatives(eta_inv,mm,margin_dist)
  nd_impact_F2=calc_Fx2_derivatives(eta_inv,mm,margin_dist)

  ### MARGIN LIKELIHOOD DERIVATIVES
  margin_deriv_subnames=c("m","d","v","t")
  names(margin_deriv_subnames)=c("mu","sigma","nu","tau")
  margin_par=names(mm)[names(mm) %in% c("mu","sigma","nu","tau")]
  response=dataset$response
  order_margin=dataset[,c("time","subject")]

  dcdu1=copula_derivatives$dcdu1; dcdu2=copula_derivatives$dcdu2;d2cdu12=copula_derivatives$d2cdu12;d2cdu22=copula_derivatives$d2cdu22

  ######################For each parameter...

  d1_cop=d2_cop=matrix(0,nrow=nrow(dataset),ncol=length(margin_par))
  colnames(d1_cop)=colnames(d2_cop)=margin_par
  for (par_name in margin_par) {
    if(object$include_dlcopdpar==TRUE | include_dlcopdpar==TRUE) {

      order_copula=matrix(ncol=4,nrow=0)
      for (i in 1:(num_margins-1)) {
        order_copula=rbind(order_copula,cbind(dataset[dataset$time == margin_names[i],c("time","subject")],dataset[dataset$time == margin_names[i+1],c("time","subject")]))
      }
      names(order_copula)=c("time1","subject1","time2","subject2")

      margin_deriv_1=matrix(0,ncol=length(margin_par),nrow=length(response))
      colnames(margin_deriv_1)=paste("dld",margin_par,sep="")
      margin_deriv_1[,paste("dld",par_name,sep="")]=margin_derivatives[grepl("dld",names(margin_derivatives))][[which(margin_par==par_name)]]

      #COPULA DERIVS WITH RESPECT TO

      mu=eta_inv[["mu"]]
      F_nd=nd_impact_F[[par_name]]
      F_nd2=nd_impact_F2[[par_name]]

      margin_components=cbind(order_margin,response,margin_p,margin_d,margin_deriv_1,mu,F_nd,F_nd2)
      margin_components_Ft_plus=margin_components
      margin_components_Ft_plus$time=margin_components_Ft_plus$time-1
      margin_plus=merge(margin_components,margin_components_Ft_plus,by=c("time","subject"),all.x=TRUE)

      copula_components=cbind(order_copula,dcdu1,dcdu2,copula_d,d2cdu12,d2cdu22)
      copula_merged=merge(copula_components,margin_plus,by.x=c("time1","subject1"),by.y=c("time","subject"),all.x=TRUE)

      #Calculate copula derivative with respect to marginal parameters
      input=copula_merged
      d1_cop[,par_name]=calc_deriv_copula_wrt_margin(input,margin_par,par_name,calc_d2=FALSE)[,which(margin_par==par_name)]
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
      d2_all[[par_name]]=c(margin_derivatives[[par_name]]+d2_cop[,i])
      i=i+1
    }
  }
  for(par_name in c(c_d2_names)) {
    d2_all[[par_name]]=c(copula_derivatives[[par_name]])
  }

  d2_all_diag=d2_all[endsWith(names(d2_all),"2")]

  names(d2_all_diag)=names(mm)

  #Score function does 1 parameter at a time
  vcov_mat=list()
  for (par_name in names(mm)) {

    d1=d1_all[[par_name]]
    d2=d2_all_diag[[par_name]]
    score=score_function_v2(eta=eta[[par_name]],dldpar=d1,d2ldpar=d2,dpardeta=eta_dr[[par_name]])
    X=as.matrix(mm[[par_name]])
    W=diag(as.vector(score$w_k))

    vcov_mat[[par_name]]=diag(solve(t(X)%*%W%*%X))
    #print(vcov_mat[[par_name]])
    #se_out[par_name]=diag(vcov_mat[[par_name]])
  }

  ########CALCULATE COPULA DERIVATIVES

  #dldth=copula_derivatives$dldth; dcdth=copula_derivatives$dcdth; dcdu1=copula_derivatives$dcdu1; dcdu2=copula_derivatives$dcdu2
  #if(!"zeta" %in% names(eta_inv)) {dldz=copula_derivatives$dldz; dcdz=copula_derivatives$dcdz}

###NEXT STEPS:
  ### USE REAL SECOND DERIVATIVES
  ### USE INCORPORATE CROSS-PARAMETER DERIVATIVES LIKE dldmudtheta

  return(vcov_mat)
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

  cop_par=BiCopTau2Par(family=as.numeric(BiCopName(copula_dist))
                       ,tau=cor(dataset[dataset$time%in%(margin_names[1:(num_margins-1)]),"response"]
                                ,dataset[dataset$time%in%(margin_names[2:(num_margins)]),"response"],method="kendall"))
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
    print("ERROR: MARGIN DISTRIBUTION STARTING VALUES NOT DEFINED, STARTING FROM GAMLSS START")
    start_fit=gamlss(dataset$response~1, family=margin_dist,method=RS(1))
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

      input_list=list(p=copsim[,i],mu=exp(par.margin[1]+par.margin[2]*i),sigma=exp(1),nu=1,tau=0.1)
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      input_list=list(p=copsim[,i],mu=exp(par.margin[1]+par.margin[2]*i+par.margin[3]),sigma=exp(1),nu=1,tau=0.1)
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed
    }

    response = as.data.frame(margin)

    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

  }
  else if (simOption==4) {

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

      input_list=list(p=copsim[,i],mu=exp(par.margin[1]),sigma=exp(par.margin[2]),nu=par.margin[3],tau=logit_inv(par.margin[4]))
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      input_list=list(p=copsim[,i],mu=exp(par.margin[1]),sigma=exp(par.margin[2]),nu=par.margin[3],tau=logit_inv(par.margin[4]))
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed
    }

    response = as.data.frame(margin)

    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

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

      input_list=list(p=copsim[,i],mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
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

    par.copula=par.copula["theta"]+covariates_input$theta.time*1:(d-1)#+covariates_input$theta.age*covariates[[1]]+covariates_input$theta.gender*covariates[[3]]
    par.copula=copula_input$copula_link$theta.linkinv(par.copula)

    if(length(par.copula)==d-1){
      par=c((par.copula),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(par.copula[1:(length(par.copula)/2)], rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(par.copula[(length(par.copula)/2+1):(length(par.copula))],rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- D2RVine(order, family, par, par2)
    #contour(RVM)

    copsim=RVineSim(n,RVM)

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=   if("mu.linkfun" %in% names(margin_dist))    {margin_dist$mu.linkinv(   rep(par.margin[1],n) +covariates_input$mu.time*covariates[[2]][,i]    + covariates_input$mu.age*covariates[[1]]    + covariates_input$mu.gender*covariates[[3]])} else {rep(0,n)}
                      ,sigma=if("sigma.linkfun" %in% names(margin_dist)) {margin_dist$sigma.linkinv(rep(par.margin[2],n) +covariates_input$sigma.time*covariates[[2]][,i] + covariates_input$sigma.age*covariates[[1]] + covariates_input$sigma.gender*covariates[[3]])} else {rep(0,n)}
                      ,nu=   if("nu.linkfun" %in% names(margin_dist))    {margin_dist$nu.linkinv(   rep(par.margin[3],n) +covariates_input$nu.time*covariates[[2]][,i]    + covariates_input$nu.age*covariates[[1]]    + covariates_input$nu.gender*covariates[[3]])} else {rep(0,n)}
                      ,tau=  if("tau.linkfun" %in% names(margin_dist))   {margin_dist$tau.linkinv(  rep(par.margin[4],n) +covariates_input$tau.time*covariates[[2]][,i]   + covariates_input$tau.age*covariates[[1]]   + covariates_input$tau.gender*covariates[[3]])} else {rep(0,n)}
                      )
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]

      qFunOutput_1=matrix(0,ncol=1,nrow=n)
      for ( j in 1:n) {
        input_list_one_item=list()
        for (item in names(input_list)) {
          input_list_one_item=append(input_list_one_item,as.matrix(input_list[[item]])[j,])
        }
        names(input_list_one_item)=names(input_list)
        qFunOutput_1[j,]=do.call(qFUN,args=(input_list_one_item[args]))
      }

      #qFunOutput_1=do.call(qFUN,args=(input_list[args]))
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

    a=.25;b=2;mu=1:d
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
  }

  dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

  if(plot_dist==TRUE) {plotDist(dataset,margin_dist)}

  return(dataset)
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

  par_eta_end=(1-step_size)*par_eta+step_size*(par_eta+weights_eta%*%score_eta)

  par_end=par*0
  names(par_end)=names(par_eta_end)=names(par)
  #Get end paraemters re-transformed
  #for (par_name in names(par)) {
  #  if(par_name %in% names(margin_par)) {
  #    par_eta_end[par_name]=margin_dist[[paste(par_name,".linkfun",sep="")]](par_end[par_name])
  #  }
  #  if(par_name %in% names(copula_par)) {
  #    par_eta_end[par_name]=copula_link[[paste(par_name,".linkfun",sep="")]](par_end[par_name])
  #  }
  #}
  ###If calculating for eta
  for (par_name in names(par)) {
   if(par_name %in% names(margin_par)) {
     par_end[par_name]=margin_dist[[paste(par_name,".linkinv",sep="")]](par_eta_end[par_name])
   }
   if(par_name %in% names(copula_par)) {
     par_end[par_name]=copula_link[[paste(par_name,".linkinv",sep="")]](par_eta_end[par_name])
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

create_model_matrices<-function(
    mu.formula = ("response ~ 1"),#mu.formula = formula("response ~ as.factor(time)+as.factor(gender)+age")
    sigma.formula = ("~ 1"),#sigma.formula = formula("~ as.factor(time)+age")
    nu.formula = ("~ 1"),#nu.formula = formula("~ as.factor(time)+as.factor(gender)")
    tau.formula = ("~ 1"),#tau.formula = formula("~ age")
    theta.formula=("~1"),#theta.formula=formula("response~as.factor(gender)")
    zeta.formula=("~1"),#zeta.formula=formula("response~1")
    margin.family=NO(),
    copula.family="N",
    copula.link=NA,
    start="nofit"
) {

  if(copula.family %in% c("t")){two_par_cop=TRUE} else {two_par_cop=FALSE}

  #Turn text formula inputs into formulas
  mu.formula = formula(mu.formula)#mu.formula = formula("response ~ as.factor(time)+as.factor(gender)+age")
  sigma.formula = formula(sigma.formula)#sigma.formula = formula("~ as.factor(time)+age")
  nu.formula = formula(nu.formula)#nu.formula = formula("~ as.factor(time)+as.factor(gender)")
  tau.formula = formula(tau.formula)#tau.formula = formula("~ age")
  theta.formula=formula(paste("response",theta.formula,sep=""))#theta.formula=formula("response~as.factor(gender)")
  zeta.formula=formula(paste("response",zeta.formula,sep=""))#zeta.formula=formula("response~1")

  #Fit gamlss marginal models
  invisible(capture.output(gamlss_model <- gamlss(formula = mu.formula
                         , sigma.formula = sigma.formula
                         , nu.formula = nu.formula
                         , tau.formula = tau.formula
                         , family=margin.family,data=dataset,method = RS(1))))

  #fitted_copulas=fitted_margins=list()
  #fitted_margins[[1]]=gamlss_model

  #Extract margin model matrix
  mm=list()
  for (parameter in gamlss_model$parameters) {
    mm[[parameter]]= model.matrix(gamlss_model,what=parameter)
  }

  #Fit starting copula with just an intercept for each parameter

  #if(start=="fit") {
  #  gamlss_unif_resid=pnorm(gamlss_model$residuals)
  #  timepoints=unique(dataset$time)
  #  copula_response=matrix(ncol=2,nrow=0)
  #  for (i in timepoints[1:(length(timepoints)-1)]) {
  #    copula_response=rbind(copula_response,cbind(gamlss_unif_resid[dataset$time==i],gamlss_unif_resid[dataset$time==i+1]))
  #  }
  #  fitted_copula=BiCopEst(copula_response[,1],copula_response[,2],family=BiCopName(copula.family))
  #  start_par<-c(extract_parameters_from_fitted(fitted_margins,fitted_copulas=list(fitted_copula),copula_link=copula.link))
  #} else {
  #  cop_par=log(c(BiCopTau2Par(family=BiCopName(copula.family),tau=
  #                               cor( dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],"response"],
  #                                    dataset[dataset$time %in% unique(dataset$time)[2:(length(unique(dataset$time)))],"response"],
  #                                    method="kendall"))))
  #  if(two_par_cop) {cop_par=c(cop_par,2);names(cop_par)=c("theta","zeta")}else{names(cop_par)="theta"}
  #  start_par<-c(extract_parameters_from_fitted(fitted_margins,fitted_copulas=NA,copula_link=copula.link),cop_par)
  #  fitted_copula=NA
  #}

  #Get copula model matrix
  #mm_cop=generate_cop_model_matrix(dataset=dataset,formula=theta.formula,zeta.formula=zeta.formula,time="time")
  invisible(capture.output(mm[["theta"]]<-model.matrix(gamlss(formula=theta.formula
                                                              ,data=(dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],]),method=RS(1)))))
  if(two_par_cop) {
    invisible(capture.output(mm[["zeta"]]<-model.matrix(gamlss(formula=zeta.formula
                                                              ,data=(dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],]),method=RS(1)))))
  }

  #print("Extracting start parameters and model matrices")
  #Create parameter vector from model matrix for copulas
  #temp_cop_parameter_names=c()
  #for (parameter in names(mm_cop)) {
  #  temp_cop_parameter_names=c(temp_cop_parameter_names,paste(parameter,colnames(mm_cop[[parameter]]),sep="."))
  #}

  #Setting up starting parameter vector with correct factors from mm_mar and mm_cop
  #theta_par_loc=grepl("theta",names(start_par))
  #if(two_par_cop) {
  #  zeta_par_loc=grepl("zeta",names(start_par))
  #}

  #temp_cop_start_par=vector(length=length(temp_cop_parameter_names))
  #names(temp_cop_start_par)=temp_cop_parameter_names

  #temp_cop_start_par[grepl("Intercept",names(temp_cop_start_par))&grepl("theta",names(temp_cop_start_par))]=start_par[theta_par_loc][1] #Theta starting value
  #if(two_par_cop) {
  #  temp_cop_start_par[grepl("Intercept",names(temp_cop_start_par))&grepl("zeta",names(temp_cop_start_par))]=start_par[zeta_par_loc][1] #Theta starting value
  #  start_par=c(start_par[!(theta_par_loc | zeta_par_loc)],temp_cop_start_par)
  #} else {
  #  start_par=c(start_par[!(theta_par_loc)],temp_cop_start_par)
  #}

  return(mm)
}
#' @export
calc_eta=function(par_cov,mm,margin_dist,copula_link) {
  eta=list()
  for (par_name in names(mm)) {
    par_cov_single=par_cov[grepl(par_name,names(par_cov))]
    mm_temp=mm[[par_name]]
    eta[[par_name]]=rowSums(mm_temp * matrix(rep(par_cov_single,each=nrow(mm_temp)),ncol=length(par_cov_single),dimnames=list(NULL,c(names(par_cov_single)))))
  }
  #Get link transforms (eta) and derivatives w.r.t to link for parameters
  eta_dr=eta_inv=list()
  for (par_name in names(mm)) {
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
#' @export
calc_F_x <- function(eta_inv,mm,margin_dist) {
  #Setup input matrix of response and parameters
  response=dataset$response
  num_margins=length(unique(dataset$time))
  margin_names=unique(dataset$time)

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
calc_likelihood_minimal <- function(eta_inv,mm,margin_dist,copula_dist,calc_d2=FALSE,response,margin_names) {
  #Setup input matrix of response and parameters
  #response=dataset$response

  num_margins=length(unique(dataset$time))

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
    margin_deriv[[deriv_name]]=do.call(FUN,args=margin_deriv_input[FUN_args])
  }

  margin_pFUN=eval(parse( text=paste("p",margin_dist$family[1],sep="") ))
  FUN=margin_pFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_p=do.call(FUN,args=margin_deriv_input[FUN_args])

  margin_dFUN=eval(parse( text=paste("d",margin_dist$family[1],sep="") ))
  FUN=margin_dFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_d=do.call(FUN,args=margin_deriv_input[FUN_args])

  #First calculate margin F(x1), F(x2) as inputs to copula

  Fx_1_2=matrix(ncol=2,nrow=0)
  order_copula=matrix(ncol=4,nrow=0)
  for (i in 1:(num_margins-1)) {
    Fx_1_2=rbind(Fx_1_2,cbind(margin_p[dataset$time == margin_names[i]],margin_p[dataset$time == margin_names[i+1]]))
    order_copula=rbind(order_copula,cbind(dataset[dataset$time == margin_names[i],c("time","subject")],dataset[dataset$time == margin_names[i+1],c("time","subject")]))
  }
  names(order_copula)=c("time1","subject1","time2","subject2")

  Fx_1_2[Fx_1_2>1]=1;Fx_1_2[Fx_1_2<0]=0

  par1=eta_inv[["theta"]]
  if(!"zeta" %in% names(eta_inv)) {par2=eta_inv[["theta"]]*0} else {par2=eta_inv[["zeta"]]}

  if(copula_dist=="C") {
    par1[par1>=28]=27.9
  }

  copula_d=BiCopPDF(  Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=as.vector(par1),par2=as.vector(par2))

  log_lik=c(sum(log(margin_d)),sum(log(copula_d)),sum(log(margin_d))+sum(log(copula_d)))
  names(log_lik)=c("marginal","copula","joint")
  return_list=list(log_lik,margin_d,copula_d,margin_p,Fx_1_2,order_copula,margin_deriv,order_copula)
  names(return_list)=c("log_lik","margin_d","copula_d","margin_p","Fx_1_2","order_copula","margin_deriv","order_copula")
  return(return_list)
}
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

      input_list=list(p=copsim[,i],mu=exp(par.margin[1]+par.margin[2]*i),sigma=exp(1),nu=1,tau=0.1)
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      input_list=list(p=copsim[,i],mu=exp(par.margin[1]+par.margin[2]*i+par.margin[3]),sigma=exp(1),nu=1,tau=0.1)
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed
    }

    response = as.data.frame(margin)

    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

  }
  else if (simOption==4) {

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

      input_list=list(p=copsim[,i],mu=exp(par.margin[1]),sigma=exp(par.margin[2]),nu=par.margin[3],tau=logit_inv(par.margin[4]))
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      input_list=list(p=copsim[,i],mu=exp(par.margin[1]),sigma=exp(par.margin[2]),nu=par.margin[3],tau=logit_inv(par.margin[4]))
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed
    }

    response = as.data.frame(margin)

    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

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

      input_list=list(p=copsim[,i],mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
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

    par.copula=par.copula["theta"]+covariates_input$theta.time*1:(d-1)#+covariates_input$theta.age*covariates[[1]]+covariates_input$theta.gender*covariates[[3]]
    par.copula=copula_input$copula_link$theta.linkinv(par.copula)

    if(length(par.copula)==d-1){
      par=c((par.copula),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(par.copula[1:(length(par.copula)/2)], rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(par.copula[(length(par.copula)/2+1):(length(par.copula))],rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- D2RVine(order, family, par, par2)
    #contour(RVM)

    copsim=RVineSim(n,RVM)

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=   if("mu.linkfun" %in% names(margin_dist))    {margin_dist$mu.linkinv(   rep(par.margin[1],n) +covariates_input$mu.time*covariates[[2]][,i]    + covariates_input$mu.age*covariates[[1]]    + covariates_input$mu.gender*covariates[[3]])} else {rep(0,n)}
                      ,sigma=if("sigma.linkfun" %in% names(margin_dist)) {margin_dist$sigma.linkinv(rep(par.margin[2],n) +covariates_input$sigma.time*covariates[[2]][,i] + covariates_input$sigma.age*covariates[[1]] + covariates_input$sigma.gender*covariates[[3]])} else {rep(0,n)}
                      ,nu=   if("nu.linkfun" %in% names(margin_dist))    {margin_dist$nu.linkinv(   rep(par.margin[3],n) +covariates_input$nu.time*covariates[[2]][,i]    + covariates_input$nu.age*covariates[[1]]    + covariates_input$nu.gender*covariates[[3]])} else {rep(0,n)}
                      ,tau=  if("tau.linkfun" %in% names(margin_dist))   {margin_dist$tau.linkinv(  rep(par.margin[4],n) +covariates_input$tau.time*covariates[[2]][,i]   + covariates_input$tau.age*covariates[[1]]   + covariates_input$tau.gender*covariates[[3]])} else {rep(0,n)}
      )
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]

      qFunOutput_1=matrix(0,ncol=1,nrow=n)
      for ( j in 1:n) {
        input_list_one_item=list()
        for (item in names(input_list)) {
          input_list_one_item=append(input_list_one_item,as.matrix(input_list[[item]])[j,])
        }
        names(input_list_one_item)=names(input_list)
        qFunOutput_1[j,]=do.call(qFUN,args=(input_list_one_item[args]))
      }

      #qFunOutput_1=do.call(qFUN,args=(input_list[args]))
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
  }

  dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

  if(plot_dist==TRUE) {plotDist(dataset,margin_dist)}

  return(dataset)
}
