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
gamlss.longitudinal=function(  dataset,
                        margin_dist,
                        copula_dist,
                        mu.formula = ("response ~ 1"),
                        sigma.formula = ("1"),
                        nu.formula = ("1"),
                        tau.formula = ("1"),
                        theta.formula=("1"),
                        zeta.formula=("1"),
                        include_dlcopdpar=FALSE,
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
                        max_outer_iter=20,
                        max_inner_iter=20,
                        use_Rcpp=FALSE,
                        lambda_start=5,
                        lambda_penalty_K=2
                      )
{
  #include_dlcopdpar=FALSE;      inner_stop_crit=.01;                        outer_stop_crit=.0001;                        start_step_size=.5;                        step_adjustment=.5;
  #                      max_steps=5;                        start_from=NA;                        verbose=TRUE;                        plot_results=TRUE;
  #                      true_val=NA;                        method="RS";
  #                      max_outer_iter=20;                        max_inner_iter=20;                        use_Rcpp=FALSE;
  #                      lambda_start=1;                        lambda_penalty_K=4
  #mu.formula=mu_formula;  sigma.formula=sigma_formula;  nu.formula=nu_formula;  tau.formula=tau_formula;  theta.formula=theta_formula;  zeta.formula=zeta_formula;

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
            nd_impact_F=calc_Fx_derivatives(eta_inv,mm,margin_dist,response=dataset$response)

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
          print(sprintf("λ=%.3f | LogLik=%.2f | DF=%.2f | GAIC=%.2f\n", 
                     lambda_val, loglik, df_total, gaic_val))

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
              print(paste("Chosen lambda:" ,round(lambda_s[[par_name]][[smooth_name]],2), "| Penalty K =", K))
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

  return_list=list(par_cov,log_lik_history,par_history,calc_lik_out_end,mm,margin_dist,copula_dist,include_dlcopdpar,dataset$response,dataset$time,dataset$subject,par_s)
  names(return_list)=c("par","log_lik_history","par_history","calc_lik_out_end","model_matrix","margin_dist","copula_dist","include_dlcopdpar","response","response_margin","response_subject","par_s")
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
    mu.formula = ("response ~ 1"),#mu.formula = formula("response ~ as.factor(time)+as.factor(gender)+age")
    sigma.formula = ("1"),#sigma.formula = formula("~ as.factor(time)+age")
    nu.formula = ("1"),#nu.formula = formula("~ as.factor(time)+as.factor(gender)")
    tau.formula = ("1"),#tau.formula = formula("~ age")
    theta.formula=("1"),#theta.formula=formula("response~as.factor(gender)")
    zeta.formula=("1"),#zeta.formula=formula("response~1")
    margin.family=NO(),
    copula.family="N",
    copula.link=NA
) {

  if(copula.family %in% c("t")){two_par_cop=TRUE} else {two_par_cop=FALSE}
  included_parameters <- c(names(margin.family$parameters), if(two_par_cop) c("theta","zeta") else c("theta"))

  formulas=list()
  for (parameter in included_parameters) {
    formulas[[parameter]]=get(paste(parameter,"formula",sep="."))
  }

  m_temp=list()
  m_temp[["mu"]]=gamlss2(formula=as.formula(mu.formula), family=NO(), data=dataset,  control=gamlss2_control(maxit=1))
  for (parameter in included_parameters[2:length(included_parameters)]) {
    formulas[[parameter]]<-as.formula(paste(as.formula(mu.formula)[[2]],formulas[[parameter]],sep="~"))
    print(formulas[[parameter]])
    m_temp[[parameter]]=gamlss2(formula=formulas[[parameter]], family=NO() #margin.family
      , data=if(parameter %in% c("theta","zeta")) (dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],]) else dataset,  control=gamlss2_control(maxit=1))
  }

  mm_x=list()
  mm_s=list()
  for(parameter in included_parameters) {
    m=m_temp[[parameter]]
    mm_x[[parameter]]=m$model[m$xterms$mu[m$xterms$mu!="(Intercept)"]]
    if(m$model[m$xterms$mu=="(Intercept)"] %>% length()>0) {
      mm_x[[parameter]]=cbind(intercept=1,mm_x[[parameter]])
    }
    if(length(m$sterms$mu)==0) { mm_s[[parameter]]=NULL} else {
      for (s in m$sterms$mu) {
        print(s)
        mm_s[[parameter]]=list()
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
  margin_names=unique(response_margin)
  num_margins=length(margin_names)

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

  ################COPULA DERIVATIVES
  #First calculate margin F(x1), F(x2) as inputs to copula

  Fx_1_2=matrix(ncol=2,nrow=0)
  order_copula=matrix(ncol=4,nrow=0)
  for (i in 1:(num_margins-1)) {
    Fx_1_2=rbind(Fx_1_2,cbind(margin_p[response_margin == margin_names[i]],margin_p[response_margin == margin_names[i+1]]))
    order_copula=rbind(order_copula,cbind(order_margin[response_margin == margin_names[i],c("time","subject")],order_margin[response_margin == margin_names[i+1],c("time","subject")]))
  }
  colnames(order_copula)=c("time1","subject1","time2","subject2")

  Fx_1_2[Fx_1_2>1]=1;Fx_1_2[Fx_1_2<0]=0

  par1=eta_inv[["theta"]]
  if(!"zeta" %in% names(eta_inv)) {par2=eta_inv[["theta"]]*0} else {par2=eta_inv[["zeta"]]}

  if(copula_dist=="C") {
    par1[par1>=28]=27.9
  }

  copula_d=BiCopPDF(  Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=as.vector(par1),par2=as.vector(par2))

  ########COMBINE MARGINS AND COPULA DERVIATIVES

  log_lik=c(sum(log(margin_d)),sum(log(copula_d)),sum(log(margin_d))+sum(log(copula_d)))
  names(log_lik)=c("marginal","copula","joint")

  return_list=list(log_lik,margin_d,copula_d,margin_p,Fx_1_2,order_copula,margin_deriv,order_copula)
  names(return_list)=c("log_lik","margin_d","copula_d","margin_p","Fx_1_2","order_copula","margin_deriv","order_copula")
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

calc_Fx_derivatives = function(eta_inv, mm, margin_dist,response) {
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

  #object=no_dl; dataset=NA; par=NA; numderiv=TRUE; sep_d2=TRUE

  #include_dlcopdpar=TRUE
  response=object$response
  response_margin=object$response_margin
  response_subject=object$response_subject

  #margin_names=unique(object$response_margin)
  #num_margins=length(margin_names)

  #se_out=object$par*0;
  margin_dist=object$margin_dist; copula_dist=object$copula_dist; copula_link=get_copula_dist(copula_dist)$copula_link
  mm=object$model_matrix

  input_par=par

  adj_fac=h
  nd_impact=rep(0,length(names(input_par)))
  names(nd_impact)=names(input_par)#[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (par_names_nd in names(input_par)) {

    change=rep(0,length(names(input_par)))
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {

      par=input_par
      par[[par_names_nd]]=par[[par_names_nd]]+adj

      eta_out=calc_eta(par_cov=par,mm=mm,margin_dist,copula_link)
      eta_inv=eta_out$eta_inv#; eta_dr=eta_out$eta_dr; eta=eta_out$eta; eta_dr=eta_out$eta_dr
      change[i]=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
      i=i+1
    }
    par=input_par
    eta_out=calc_eta(par_cov=par,mm=mm,margin_dist,copula_link)
    eta_inv=eta_out$eta_inv#; eta_dr=eta_out$eta_dr; eta=eta_out$eta; eta_dr=eta_out$eta_dr

    change[3]=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
    nd_impact[par_names_nd]=(change[2]+change[1]-2*change[3])/(adj_fac^2)

    #print(c(change,nd_impact[eta_par_names_nd]))
  }

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

            eta_out=calc_eta(par_cov=par,mm=mm,margin_dist,copula_link)
            eta_inv=eta_out$eta_inv#; eta_dr=eta_out$eta_dr; eta=eta_out$eta; eta_dr=eta_out$eta_dr

            change=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
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

  #object=no_dl; dataset=NA; par=NA; numderiv=TRUE; sep_d2=TRUE
  #object=w_dl; dataset=NA; par=NA; numderiv=TRUE; sep_d2=TRUE
  #object=out[[i]][[j]]

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
    input_par=object$par
  } else {
    input_par=par
  }

  eta_out=calc_eta(par_cov=input_par,mm=mm,margin_dist,copula_link)
  eta_inv=eta_out$eta_inv; eta_dr=eta_out$eta_dr; eta=eta_out$eta; eta_dr=eta_out$eta_dr

  #if(!all(is.na(par))) {response=eta_inv[["mu"]]}
  calc_lik_out=calc_likelihood_minimal(eta_inv=eta_inv,mm=mm,margin_dist,copula_dist,calc_d2=TRUE,response=response,response_margin=object$response_margin,response_subject=object$response_subject)
  Fx_1_2=calc_lik_out$Fx_1_2; margin_p=calc_lik_out$margin_p; margin_d=calc_lik_out$margin_d; copula_d=calc_lik_out$copula_d

  ###Calculate derivaties: margin and copula d1 and d2
  margin_derivatives=calc_lik_out$margin_deriv

  copula_derivatives=calc_copula_derivatives(eta_inv, Fx_1_2, copula_dist,calc_d2 = TRUE)

  nd_impact_F=calc_Fx_derivatives(eta_inv,mm,margin_dist,response)

  #Calculate numerical derivatives for F(x) - also used as a reference for overall likelihood derivative
  nd_impact_F=calc_Fx_derivatives(eta_inv,mm,margin_dist,response)
  nd_impact_F2=calc_Fx2_derivatives(eta_inv,mm,margin_dist,response)

  if(numderiv==TRUE) {
    #nd2_joint_lik=calc_true_SE_numderiv_only(eta_inv,mm,margin_dist,response,testing=TRUE,response_margin,response_subject)
    hessian_nd=calc_true_SE_numderiv_only_covariates(object,par=input_par, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA)
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

  # rownames(vcov_final)=colnames(vcov_final)=names(eta)
  #
  # vcov_covariates=list()
  # vcov_diag=diag(vcov_final)
  # names(vcov_diag)=rownames(vcov_final)
  # for (par_name in names(eta)) {
  #
  #   vcov_diag[par_name]
  #   eta_dr[[par_name]]
  #
  #   #var_yi=rep(vcov_diag[par_name],length(eta_dr[[par_name]]))
  #
  #   deriv_name_temp=paste("d2ld",sub_names_in[par_name],"2",sep="")
  #
  #   #print(d2_all[[deriv_name_temp]])
  #
  #   #W=diag((eta_dr[[par_name]]^2)/var_yi)
  #   W=diag((eta_dr[[par_name]]^2)*-d2_all[[deriv_name_temp]])
  #
  #   X=as.matrix(mm[[par_name]])
  #
  #   #if(par_name=="sigma") {print(eta_dr[[par_name]]);print(d2_all[[deriv_name_temp]])}
  #
  #   print(t(X)%*%W%*%X)
  #   print(solve(t(X)%*%W%*%X))
  #
  #   vcov_covariates[[par_name]]=diag((solve(t(X)%*%W%*%X)))
  #   #print(vcov_covariates)
  #
  # }

  return(list(vcov_final,se_final))

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
                      ,mu=   if("mu.linkfun" %in% names(margin_dist))    {margin_dist$mu.linkinv(   rep(par.margin[1],n) +covariates_input$mu.time*covariates[[2]][,i]    + ((covariates_input$mu.age    ))*((covariates[[1]]-50)/100)^2    + covariates_input$mu.gender*covariates[[3]])} else {rep(0,n)}
                      ,sigma=if("sigma.linkfun" %in% names(margin_dist)) {margin_dist$sigma.linkinv(rep(par.margin[2],n) +covariates_input$sigma.time*covariates[[2]][,i] + ((covariates_input$sigma.age ))*((covariates[[1]]-50)/100)^3 + covariates_input$sigma.gender*covariates[[3]])} else {rep(0,n)}
                      ,nu=   if("nu.linkfun" %in% names(margin_dist))    {margin_dist$nu.linkinv(   rep(par.margin[3],n) +covariates_input$nu.time*covariates[[2]][,i]    + ((covariates_input$nu.age    ))*(((covariates[[1]]-50)/100)^2)    + covariates_input$nu.gender*covariates[[3]])} else {rep(0,n)}
                      ,tau=  if("tau.linkfun" %in% names(margin_dist))   {margin_dist$tau.linkinv(  rep(par.margin[4],n) +covariates_input$tau.time*covariates[[2]][,i]   + ((covariates_input$tau.age   ))*(((covariates[[1]]-50)/100)^2)   + covariates_input$tau.gender*covariates[[3]])} else {rep(0,n)}
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
                      ,mu=margin_dist$mu.linkinv(margin_dist$mu.linkfun(par.margin[1])+(i-1))
                      ,sigma=par.margin[2]#margin_dist$sigma.linkinv(margin_dist$sigma.linkfun(par.margin[2])+(i-1))
                      ,nu=par.margin[3]
                      ,tau=par.margin[4])
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
                      ,mu=margin_dist$mu.linkinv(margin_dist$mu.linkfun(par.margin[1])+(i-1))
                      ,sigma=par.margin[2]#margin_dist$sigma.linkinv(margin_dist$sigma.linkfun(par.margin[2])+(i-1))
                      ,nu=par.margin[3]
                      ,tau=par.margin[4])
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
                      ,mu=margin_dist$mu.linkinv(margin_dist$mu.linkfun(par.margin[1])+(i-1))
                      ,sigma=margin_dist$sigma.linkinv(margin_dist$sigma.linkfun(par.margin[2])+(((i-1)/10)))
                      ,nu=par.margin[3]
                      ,tau=par.margin[4])
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
