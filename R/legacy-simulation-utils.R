#' Legacy simulation and research-era helpers
#'
#' These helpers are retained for compatibility and internal historical support.
#' They are not part of the core reviewer path for fitting, prediction, or diagnostics.
#'
#' @noRd
NULL

#' @keywords internal
#' @noRd
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

#' @keywords internal
#' @noRd
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
    RVM <- .copula_dvine(order, family, par, par2)
    contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n*t,RVM)

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

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)

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

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)

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

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)

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
      RVM[[i]] = .copula_dvine(order, c(rep(copula.family,length(theta_inv[i,])),rep(0,dd-(length(theta_inv[i,])))), par=c(theta_inv[i,],rep(0,dd-(length(theta_inv[i,])))), par2=c(theta_inv[i,],rep(0,dd-(length(theta_inv[i,])))))
    }
    #RVM <- .copula_dvine(order, rep(family[1],nrow(theta_inv)), theta_inv, theta_inv*0)
    #contour(RVM)

    copsim=.copula_rvine_sim(n,RVM)


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

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)


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

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)

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

    print("WARNING: SIMULATION MAPS MARGIN AND COPULA PARAMETERS THROUGH LINK-INVERSE FUNCTIONS AFTER ADDING COVARIATE EFFECTS.")

    t=d

    # Setup covariates from covariates_input

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,2),0)) #Gender

    # Build treatment-coded factor effects for simulation inputs.
    # Supported coefficient formats for each *.time / *.gender entry:
    # - scalar: same effect for all non-reference levels
    # - vector length (L-1): explicit non-reference effects
    # - vector length L: full per-level effects
    resolve_factor_effect <- function(levels, coef_input, label) {
      lvl <- as.character(levels)
      n_lvl <- length(lvl)

      if (is.null(coef_input) || any(is.na(coef_input))) {
        return(stats::setNames(rep(0, n_lvl), lvl))
      }

      coef_vec <- as.numeric(coef_input)
      if (length(coef_vec) == 1) {
        out <- c(0, rep(coef_vec, max(0, n_lvl - 1)))
      } else if (length(coef_vec) == (n_lvl - 1)) {
        out <- c(0, coef_vec)
      } else if (length(coef_vec) == n_lvl) {
        out <- coef_vec
      } else {
        stop(
          "simOption 10 factor effect '", label, "' has invalid length ", length(coef_vec),
          ". Expected 1, ", n_lvl - 1, ", or ", n_lvl, " for levels: ",
          paste(lvl, collapse = ", "),
          "."
        )
      }

      stats::setNames(out, lvl)
    }

    make_time_factor_component <- function(coef_input, n_cols, label) {
      levels <- as.character(seq_len(n_cols))
      level_effects <- resolve_factor_effect(levels, coef_input, label)
      matrix(rep(level_effects[levels], n), ncol = n_cols, byrow = TRUE)
    }

    make_gender_factor_component <- function(coef_input, n_cols, label) {
      gender_vals <- as.character(as.vector(covariates[[3]][, 1]))
      levels <- sort(unique(gender_vals))
      level_effects <- resolve_factor_effect(levels, coef_input, label)
      subj_effect <- as.numeric(level_effects[gender_vals])
      matrix(rep(subj_effect, n_cols), ncol = n_cols)
    }

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    apply_margin_link <- function(par_name, par_value, eta_component) {
      if (par_name %in% names(margin_dist$parameters)) {
        linkfun_name <- paste0(par_name, ".linkfun")
        linkinv_name <- paste0(par_name, ".linkinv")
        par_eta_base <- eval(parse(text=paste0("margin_dist$", linkfun_name)))(par_value)
        par_eta <- par_eta_base + eta_component
        return(eval(parse(text=paste0("margin_dist$", linkinv_name)))(par_eta))
      }

      return(NULL)
    }

    mu_out = NULL
    sigma_out = NULL
    nu_out = NULL
    tau_out = NULL

    if ("mu" %in% names(margin_dist$parameters)) {
      mu_eta = make_time_factor_component(covariates_input$mu.time, d, "mu.time") +
        matrix(rep(as.matrix(covariates_input$mu.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$mu.gender, d, "mu.gender")
      mu_out = apply_margin_link("mu", par.margin[1], mu_eta)
    }

    if ("sigma" %in% names(margin_dist$parameters)) {
      sigma_eta = make_time_factor_component(covariates_input$sigma.time, d, "sigma.time") +
        matrix(rep(as.matrix(covariates_input$sigma.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$sigma.gender, d, "sigma.gender")
      sigma_out = apply_margin_link("sigma", par.margin[2], sigma_eta)
    }

    if ("nu" %in% names(margin_dist$parameters)) {
      nu_eta = make_time_factor_component(covariates_input$nu.time, d, "nu.time") +
        matrix(rep(as.matrix(covariates_input$nu.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$nu.gender, d, "nu.gender")
      nu_out = apply_margin_link("nu", par.margin[3], nu_eta)
    }

    if ("tau" %in% names(margin_dist$parameters)) {
      tau_eta = make_time_factor_component(covariates_input$tau.time, d, "tau.time") +
        matrix(rep(as.matrix(covariates_input$tau.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$tau.gender, d, "tau.gender")
      tau_out = apply_margin_link("tau", par.margin[4], tau_eta)
    }
    theta_eta_out=par.copula[1]+make_time_factor_component(covariates_input$theta.time, d - 1, "theta.time") +
      matrix(rep(as.matrix(covariates_input$theta.age*((covariates[[1]]-50)/100)^2),d-1),ncol=d-1) +
      make_gender_factor_component(covariates_input$theta.gender, d - 1, "theta.gender")
    theta_out = copula_input$copula_link$theta.linkinv(theta_eta_out)

    if ("zeta" %in% copula_input$parameters) {
      zeta_eta_out=par.copula[2]+make_time_factor_component(covariates_input$zeta.time, d - 1, "zeta.time") +
        matrix(rep(as.matrix(covariates_input$zeta.age*((covariates[[1]]-50)/100)^2),d-1),ncol=d-1) +
        make_gender_factor_component(covariates_input$zeta.gender, d - 1, "zeta.gender")
      zeta_out = copula_input$copula_link$zeta.linkinv(zeta_eta_out)
    } else {
      zeta_out = matrix(0, nrow = n, ncol = d - 1)
    }

    if (!is.null(mu_out) && any(!is.finite(mu_out))) {
      stop("simOption 10 generated non-finite mu values after link inverse transformation.")
    }
    if (!is.null(sigma_out) && any(!is.finite(sigma_out))) {
      stop("simOption 10 generated non-finite sigma values after link inverse transformation.")
    }
    if (!is.null(nu_out) && any(!is.finite(nu_out))) {
      stop("simOption 10 generated non-finite nu values after link inverse transformation.")
    }
    if (!is.null(tau_out) && any(!is.finite(tau_out))) {
      stop("simOption 10 generated non-finite tau values after link inverse transformation.")
    }

    if (any(!is.finite(theta_out))) {
      stop("simOption 10 generated non-finite theta values after link inverse transformation.")
    }
    if ("zeta" %in% copula_input$parameters && any(!is.finite(zeta_out))) {
      stop("simOption 10 generated non-finite zeta values after link inverse transformation.")
    }

    # Print parameter ranges for quick simulation diagnostics.
    range_str <- function(label, x) {
      if (is.null(x)) return(paste0(label, ": [NA, NA]"))
      sprintf("%s: [%.2f, %.2f]", label, min(x), max(x))
    }
    margin_range_msg <- paste(
      c(
        range_str("MU", mu_out),
        range_str("SIGMA", sigma_out),
        range_str("NU", nu_out),
        range_str("TAU", tau_out)
      ),
      collapse = " | "
    )
    copula_range_msg <- paste0(
      sprintf("THETA: [%.2f, %.2f]", min(theta_out), max(theta_out)),
      if ("zeta" %in% copula_input$parameters) sprintf(" | ZETA: [%.2f, %.2f]", min(zeta_out), max(zeta_out)) else ""
    )
    print(paste("MARGIN RANGES ->", margin_range_msg, "| COPULA RANGES ->", copula_range_msg))

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
      par2_r <- c(
        if ("zeta" %in% copula_input$parameters) as.numeric(zeta_out[r, ]) else rep(0, d - 1),
        rep(0, dd - (d - 1))
      )
      copsim[r, ] <- as.numeric(.copula_rvine_sim(1, .copula_dvine(order, family, par_r, par2_r)))
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

        if(plot_dist==TRUE) {
          plot_dist(
            dataset,
            margin_dist = margin_dist,
            subject_var = "subject",
            time_var = "time",
            response_var = "response"
          )
        }

  return(dataset)
}
#' @keywords internal
#' @noRd
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
#' @keywords internal
#' @noRd
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

  if(all(is.null(names(par))|is.na(names(par)))) {stop("ERROR: par vector must be named")}
  margin_par=par[names(par)%in%c("mu","sigma","nu","tau")]
  copula_par=par[!names(par)%in%c("mu","sigma","nu","tau")]

  ##### Calculate all relevant derivatives / CG method with first and second derivatives

  ### Calculate margin derivatives w.r.t. margin parameters

  #Get names for margin derivatives from margin_dist
            n_par <- length(eta[[par_name]])
            d1_full=matrix(0,nrow=n_par,ncol=1)

  #Get link transforms (eta) and derivatives w.r.t to link for parameters
              if(n_par == length(dataset$response)) {
                par_idx <- row_id1
              } else {
                margin_names = sort(unique(dataset$time))
                theta_rows = which(dataset$time %in% margin_names[seq_len(max(1, length(margin_names)-1))])
                theta_index_map=rep(NA_integer_,length(dataset$response))
                theta_index_map[theta_rows]=seq_along(theta_rows)
                par_idx <- theta_index_map[row_id1]
              }

              valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= nrow(d1_full)
              d1_full[par_idx[valid_idx],1] <- dldth[valid_idx]
    if(par_name %in% names(margin_par)) {
      par_eta[par_name]=margin_dist[[paste(par_name,".linkfun",sep="")]](par[par_name])
      par_eta_dr[par_name]=margin_dist[[paste(par_name,".dr",sep="")]](par_eta[par_name])
    }
    if(par_name %in% names(copula_par)) {
            n_par <- length(eta[[par_name]])
            d1_full=matrix(0,nrow=n_par,ncol=1)
      par_eta_dr[par_name]=copula_link[[paste(par_name,".dr",sep="")]](par_eta[par_name])
    }
              if(n_par == length(dataset$response)) {
                par_idx <- row_id1
              } else {
                margin_names = sort(unique(dataset$time))
                theta_rows = which(dataset$time %in% margin_names[seq_len(max(1, length(margin_names)-1))])
                theta_index_map=rep(NA_integer_,length(dataset$response))
                theta_index_map[theta_rows]=seq_along(theta_rows)
                par_idx <- theta_index_map[row_id1]
              }

              valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= nrow(d1_full)
              d1_full[par_idx[valid_idx],1] <- dldz[valid_idx]
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
  order_copula=data.frame()
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

  copula_d=.copula_pdf(  Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2)
  dldth=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par",log=TRUE)
  dcdth=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par",log=FALSE)
  d2cdth=.copula_deriv2( Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par")
  d2ldth2=(1/(copula_d^2))*(copula_d*d2cdth-dcdth^2)
  if(!is.na(copula_par["zeta"])) {
    dldz=.copula_deriv(    Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par2",log=TRUE)
    dcdz=.copula_deriv(    Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par2",log=FALSE)
    d2cdz=.copula_deriv2(  Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par2")
    d2ldz2=(1/(copula_d^2))*(copula_d*d2cdz-dcdz^2)

    d2cdthdz=.copula_deriv2(  Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par1par2")
    d2ldthdz=(d2cdthdz*copula_d-dcdth*dcdz)/(copula_d^2)
  }
  dcdu1=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u1",log=FALSE)
  dcdu2=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u2",log=FALSE)

  d2cdu12=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u1",log=FALSE)
  d2cdu22=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u2",log=FALSE)

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
  margin_components_Ft_plus$time=normalize_lag_time(margin_components_Ft_plus$time)
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

# Load analytical Hessian helpers when sourcing this file directly (development workflow).
local({
  hessian_files <- c(
    "hessian-setup.R",
    "hessian-margin-cdf.R",
    "hessian-copula.R",
    "hessian-assembly.R",
    "hessian-analytical.R"
  )
  candidate_dirs <- c("R", file.path(getwd(), "R"))
  for (dir in candidate_dirs) {
    candidates <- file.path(dir, hessian_files)
    if (all(file.exists(candidates))) {
      for (p in candidates) source(p, local = FALSE)
      break
    }
  }
})
