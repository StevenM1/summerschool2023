### PMWG ---

jointLL <- function(pars, dadms){
  parPreFixs <- unique(gsub("[-].*", "", names(pars)))
  i <- 0
  total_ll <- 0
  for(dadm in dadms){
    if(is.data.frame(dadm)){
      i <- i + 1
      parPrefix <- parPreFixs[i]
      currentPars <- pars[grep(paste0(parPrefix, "-"), names(pars))]
      names(currentPars) <- gsub(".*[-]", "", names(currentPars))
      total_ll <- total_ll +  attr(dadm, "model")$log_likelihood(currentPars, dadm)
    }
  }
  return(total_ll)
}


dm_list <- function(dadm)
  # Makes data model into subjects list for use by likelihood  
  # Assumes each subject has the same design.
{
  
  sub_design <- function(designs,isin)
    lapply(designs,function(x) {
      attr(x,"expand") <- attr(x,"expand")[isin]
      x
    })
  
  
  model <- attr(dadm,"model")
  p_names <- attr(dadm,"p_names")
  sampled_p_names <- attr(dadm,"sampled_p_names")
  designs <- attr(dadm,"designs")
  expand <- attr(dadm,"expand")
  s_expand <- attr(dadm,"s_expand")
  unique_nort <- attr(dadm,"unique_nort") 
  expand_nort <- attr(dadm,"expand_nort") 
  unique_nortR <- attr(dadm,"unique_nortR") 
  expand_nortR <- attr(dadm,"expand_nortR") 
  ok_trials <- attr(dadm,"ok_trials")
  ok_dadm_winner <- attr(dadm,"ok_dadm_winner")
  ok_dadm_looser <- attr(dadm,"ok_dadm_looser")
  ok_da_winner <- attr(dadm,"ok_da_winner")
  ok_da_looser <- attr(dadm,"ok_da_looser")
  expand_uc <- attr(dadm,"expand_uc") 
  expand_lc <- attr(dadm,"expand_lc") 
  
  # winner on expanded dadm  
  expand_winner <- attr(dadm,"expand_winner")
  # subjects for first level of lR in expanded dadm
  slR1=dadm$subjects[expand][dadm$lR[expand]==levels(dadm$lR)[[1]]]
  
  dl <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                 levels(dadm$subjects))
  for (i in levels(dadm$subjects)) {
    isin <- dadm$subjects==i         # dadm
    isin1 <- s_expand==i             # da
    isin2 <- attr(dadm,"s_data")==i  # data
    dl[[i]] <- dadm[isin,]
    attr(dl[[i]],"model") <- model
    attr(dl[[i]],"p_names") <- p_names
    attr(dl[[i]],"sampled_p_names") <- sampled_p_names
    attr(dl[[i]],"designs") <- sub_design(designs,isin)
    attr(dl[[i]],"expand") <- expand[isin1]-min(expand[isin1]) + 1
    attr(dl[[i]],"s_expand") <- NULL
    
    attr(dl[[i]],"ok_dadm_winner") <- ok_dadm_winner[isin]
    attr(dl[[i]],"ok_dadm_looser") <- ok_dadm_looser[isin]
    
    attr(dl[[i]],"ok_da_winner") <- ok_da_winner[isin1]
    attr(dl[[i]],"ok_da_looser") <- ok_da_looser[isin1]
    
    attr(dl[[i]],"unique_nort") <- unique_nort[isin] 
    attr(dl[[i]],"unique_nortR") <- unique_nortR[isin] 

    isinlR1 <- slR1==i
    attr(dl[[i]],"expand_nort") <-  expand_nort[isinlR1]- min( expand_nort[isinlR1]) + 1
    attr(dl[[i]],"expand_nortR") <- expand_nortR[isinlR1]-min(expand_nortR[isinlR1]) + 1
    
    attr(dl[[i]],"ok_trials") <- ok_trials[isin2]
    attr(dl[[i]],"expand_winner") <- expand_winner[isin2]-min(expand_winner[isin2]) + 1
    
    if (!is.null(attr(dadm,"expand_uc")))
      attr(dl[[i]],"expand_uc") <- as.numeric(factor(expand_uc[isin2]))
    if (!is.null(attr(dadm,"expand_lc")))
      attr(dl[[i]],"expand_lc") <- as.numeric(factor(expand_lc[isin2]))
    
    if (!is.null(attr(dadm,"LT")))
      attr(dl[[i]],"LT") <- attr(dadm,"LT")[names(attr(dadm,"LT"))==i]
    if (!is.null(attr(dadm,"UT")))
      attr(dl[[i]],"UT") <- attr(dadm,"UT")[names(attr(dadm,"UT"))==i]
    if (!is.null(attr(dadm,"LC")))
      attr(dl[[i]],"LC") <- attr(dadm,"LC")[names(attr(dadm,"LC"))==i]
    if (!is.null(attr(dadm,"UC")))
      attr(dl[[i]],"UC") <- attr(dadm,"UC")[names(attr(dadm,"UC"))==i]
  }
  dl
} 


extractDadms <- function(dadms, names = 1:length(dadms)){
  N_models <- length(dadms)
  pars <- attr(dadms[[1]], "sampled_p_names")
  prior <- attr(dadms[[1]], "prior")
  ll_func <- attr(dadms[[1]], "model")$log_likelihood
  # subjects <- unique(unlist(sapply(dadms, FUN = function(x){return(unique(x$subjects))})))
  subjects <- unique(factor(sapply(dadms, FUN = function(x) levels(x$subjects))))
  # dadm_list <- vector("lis", length = length(subjects))
  # dadm_list[as.numeric(subjects)] <- dm_list(dadms[[1]])
  # dadm_list[as.numeric(subjects)] <- dm_list(dadms[[1]])
  dadm_list <- dm_list(dadms[[1]])
  if(N_models > 1){
    k <- 1
    pars <- paste(names[1], pars, sep = "-")
    dadm_list[as.character(which(!subjects %in% unique(dadms[[1]]$subjects)))] <- NA
    for(dadm in dadms[-1]){
      k <- k + 1
      tmp_list <- vector("list", length = length(subjects))
      tmp_list[as.numeric(unique(dadm$subjects))] <- dm_list(dadm)
      dadm_list <- mapply(list, dadm_list, tmp_list, SIMPLIFY = F)
      curr_pars <- attr(dadm, "sampled_p_names")
      pars <- c(pars, paste(names[k], curr_pars, sep = "-"))
      prior$theta_mu_mean <- c(prior$theta_mu_mean, attr(dadm, "prior")$theta_mu_mean) 
      if(is.matrix(prior$theta_mu_var)){
        prior$theta_mu_var <- adiag(prior$theta_mu_var, attr(dadm, "prior")$theta_mu_var) 
      } else{
        prior$theta_mu_var <- c(prior$theta_mu_var, attr(dadm, "prior")$theta_mu_var) 
      }
    }
    ll_func <- jointLL
  }
  return(list(ll_func = ll_func, pars = pars, prior = prior, 
              dadm_list = dadm_list, subjects = subjects))
}


# type=c("standard","diagonal","blocked","factor","factorRegression","single")[1]
# n_chains=3; rt_resolution=0.02
# prior_list = NULL;par_groups=NULL;n_factors=NULL;constraintMat = NULL;covariates=NULL
# data_list=data[,-4]; design_list=design_B_MT;model_list=NULL
make_samplers <- function(data_list,design_list,model_list=NULL,
  type=c("standard","diagonal","blocked","factor","factorRegression","single")[1],
  n_chains=3,rt_resolution=0.02,
  prior_list = NULL,
  par_groups=NULL,
  n_factors=NULL,constraintMat = NULL,covariates=NULL)
  
{
  if (!(type %in% c("standard","diagonal","blocked","factor","factorRegression","single")))
    stop("type must be one of: standard,diagonal,blocked,factor,factorRegression,single")
  if (class(data_list)=="data.frame") data_list <- list(data_list)
  # Sort subject together and add unique trial within subject integer
  data_list <- lapply(data_list,function(d){
    if (!is.factor(d$subjects)) d$subjects <- factor(d$subjects)
    d <- d[order(d$subjects),]
    add_trials(d)
  })
  if (!is.null(names(design_list)[1]) && names(design_list)[1]=="Flist") 
    design_list <- list(design_list)
  if (length(design_list)!=length(data_list))
      design_list <- rep(design_list,length(data_list))
  if (is.null(model_list)) model_list <- lapply(design_list,function(x){x$model})
  if (any(unlist(lapply(model_list,is.null)))) 
    stop("Must supply model_list if model is not in all design_list components")
  if (!is.null(names(model_list)[1]) && names(model_list)[1]=="type") 
    model_list <- list(model_list)
  if (length(model_list)!=length(data_list))
      model_list <- rep(model_list,length(data_list))
  if (!is.null(names(prior_list)) && any(names(prior_list)=="theta_mu_mean"))
    prior_list <- list(prior_list)
  if (length(prior_list)!=length(data_list))
      prior_list <- rep(prior_list,length(data_list))
  dadm_list <- vector(mode="list",length=length(data_list))
  rt_resolution <- rep(rt_resolution,length.out=length(data_list))
  for (i in 1:length(dadm_list)) {
    message("Processing data set ",i)
    if (!is.null(design_list[[i]]$Ffunctions)) data_list[[i]] <- 
        cbind.data.frame(data_list[[i]],data.frame(lapply(
          design_list[[i]]$Ffunctions,function(f){f(data_list[[i]])})))
    dadm_list[[i]] <- design_model(data=data_list[[i]],design=design_list[[i]],
      model=model_list[[i]],rt_resolution=rt_resolution[i],prior=prior_list[[i]])
  }
  if (type == "standard") {
    source("samplers/pmwg/variants/standard.R")
    out <- pmwgs(dadm_list)
    out$source <- "samplers/pmwg/variants/standard.R"
  } else if(type == "diagonal"){
    source("samplers/pmwg/variants/diag.R")
    out <- pmwgs(dadm_list)
    out$source <- "samplers/pmwg/variants/diag.R"
  } else if (type == "blocked") {
    if (is.null(par_groups)) stop("Must specify par_groups for blocked type")
    source("samplers/pmwg/variants/blocked.R")
    out <- pmwgs(dadm_list,par_groups=par_groups)
    out$source <- "samplers/pmwg/variants/blocked.R"
  } else if (type=="factor") {
    if (is.null(n_factors)) stop("Must specify n_factors for factor type")
    source("samplers/pmwg/variants/factor.R") 
    out <- pmwgs(dadm_list,n_factors=n_factors,constraintMat=constraintMat)
    out$source <- "samplers/pmwg/variants/factor.R"
  } else if (type=="factorRegression") {
    if (is.null(n_factors)) stop("Must specify n_factors for factorRegression type")
    if (is.null(covariates)) stop("Must specify covariates for factorRegression type")
    source("samplers/pmwg/variants/factorRegr.R")
    out <- pmwgs(dadm_list,n_factors=n_factors,constraintMat=constraintMat,
                 covariates=covariates)
    out$source <- "samplers/pmwg/variants/factorRegr.R"
  } else if (type == "single") {
    source("samplers/pmwg/variants/single.R")
    out <- pmwgs(dadm_list)
    out$source <- "samplers/pmwg/variants/single.R"
  }
  # replicate chains
  dadm_lists <- rep(list(out),n_chains)
  # For post predict
  attr(dadm_lists,"data_list") <- data_list
  attr(dadm_lists,"design_list") <- design_list
  attr(dadm_lists,"model_list") <- model_list
  dadm_lists
}
