source("samplers/pmwg/variants/diag.R")

add_info_diag_groups <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  sampler <- add_info_diag(sampler, prior, ...)
  sampler$par_idx <- list(...)$par_idx
  sampler$parGroups <- unique(as.vector(par_idx))
  sampler$group_idx <- aggregate(group ~ subject, sampler$data, mean)[,2]
  return(sampler)
}

sample_store_diag_groups <- function(data, par_names, iters = 1, stage = "init", integrate = T, ...) {
  n_pars <- length(par_names)
  groups <- unique(data$group)
  n_groups <- length(groups)
  samples <- sample_store_standard(data, par_names, iters, stage, integrate, ...)
  samples$theta_mu = array(NA_real_,dim = c(n_pars, n_groups, iters), dimnames = list(par_names, groups, NULL))
  return(samples)
}

get_startpoints_diag_groups<- function(pmwgs, start_mu, start_var){
  startpoints <- get_startpoints_diag(pmwgs, start_mu, start_var)
  startpoints$tmu <- replicate(max(pmwgs$group_idx), startpoints$tmu)
  startpoints$group_idx <- pmwgs$group_idx
  return(startpoints)
}

fill_samples_diag_groups <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  # group level
  samples$theta_mu[,, j] <- group_level$tmu
  samples$theta_var[, , j] <- group_level$tvar
  samples$last_theta_var_inv <- group_level$tvinv
  samples$a_half[,j] <- group_level$a_half
  
  # Random effects
  samples <- fill_samples_RE(samples, proposals, epsilon,j, n_pars)
  return(samples)
}

gibbs_step_diag <- function(sampler, alpha){
  # Gibbs step for diagonal only
  # Get single iter versions, tmu = theta_mu, tvar = theta_var
  last <- last_sample_standard(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  last$tvinv <- diag(last$tvinv)
  
  #Mu
  var_mu = 1.0 / (sampler$n_subjects * last$tvinv + prior$theta_mu_invar)
  mean_mu = var_mu * ((apply(alpha, 1, sum) * last$tvinv + prior$theta_mu_mean * prior$theta_mu_invar))
  tmu <- rnorm(sampler$n_pars, mean_mu, sd = sqrt(var_mu))
  names(tmu) <- sampler$par_names
  
  if(!is.null(hyper$std_shape)){
    # InvGamma alternative (probably inferior) prior
    shape = hyper$std_shape + sampler$n_subjects / 2
    rate = hyper$std_rate + rowSums( (alpha-tmu)^2 ) / 2
    tvinv = rgamma(n=sampler$n_pars, shape=shape, rate=rate)
    tvar = 1/tvinv
    a_half <- NULL
  } else {
    tvinv = rgamma(n=sampler$n_pars, shape=hyper$v_half/2 + sampler$n_subjects/2, rate=hyper$v_half/last$a_half + 
                     rowSums( (alpha-tmu)^2 ) / 2)
    tvar = 1/tvinv
    #Contrary to standard pmwg I use shape, rate for IG()
    a_half <- 1 / rgamma(n = sampler$n_pars, shape = (hyper$v_half + sampler$n_pars) / 2,
                         rate = hyper$v_half * tvinv + 1/(hyper$A_half^2))
  }
  return(list(tmu = tmu, tvar = diag(tvar), tvinv = diag(tvinv), a_half = a_half, alpha = alpha))
}

gibbs_step_diag_groups <- function(sampler, alpha){
  # Gibbs step for group level with multiple means, 
  # mu = theta_mu, var = theta_var
  last <- last_sample_diag_groups(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  last$tvinv <- diag(last$tvinv)
  
  
  #extract previous values (for ease of reading)
  n_subjects <- sampler$n_subjects
  n_pars <- sampler$n_pars
  
  tmu <- last$tmu
  alphatilde <- alpha
  #Update mu
  for(parGroup in sampler$parGroups){
    par_idx <- unique(which(sampler$par_idx == parGroup, arr.ind = T)[,1])
    parGroup_num <- as.numeric(strsplit(parGroup,",")[[1]])
    group_idx <- which(sampler$group_idx %in% parGroup_num)
    var_mu = 1 / (length(group_idx) * last$tvinv[par_idx] + prior$theta_mu_invar[par_idx])
    mean_mu = var_mu * ((rowSums(alpha[par_idx, group_idx, drop = F]) * last$tvinv[par_idx] + prior$theta_mu_mean[par_idx] * prior$theta_mu_invar[par_idx]))
    tmu[par_idx,parGroup_num] <- rnorm(length(par_idx), mean_mu, sd = sqrt(var_mu))
  }
  for (group in unique(sampler$group_idx)){
    group_idx <- which(sampler$group_idx == group)
    alphatilde[,group_idx] <- sweep(alphatilde[,group_idx], 1, tmu[,group]) # This will break if groups != 1:length(groups)
  }
  tvinv = rgamma(n=sampler$n_pars, shape=hyper$v_half/2 + sampler$n_subjects/2, rate=hyper$v_half/last$a_half + 
                   rowSums( (alphatilde)^2 ) / 2)
  tvar = 1/tvinv
  a_half <- 1 / rgamma(n = sampler$n_pars, shape = (hyper$v_half + sampler$n_pars) / 2,
                       rate = hyper$v_half * tvinv + 1/(hyper$A_half^2))
  return(list(tmu = tmu, tvar = diag(tvar), tvinv = diag(tvinv), a_half = a_half, 
              alpha = alpha, group_idx = sampler$group_idx))
}

last_sample_diag_groups <- function(store) {
  list(
    tmu = store$theta_mu[,, store$idx],
    tvar = store$theta_var[, , store$idx],
    tvinv = store$last_theta_var_inv,
    a_half = store$a_half[, store$idx]
  )
}

get_group_level_diag_groups <- function(parameters, s){
  mu <- parameters$tmu[,parameters$group_idx[s]]
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

get_conditionals_diag_groups <- function(s, samples, n_pars){
  groups <- samples$groups
  iteration <- samples$iteration
  pts2_unwound <- log(apply(samples$theta_var,3,diag))
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu[,groups[s],],pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- stats::var(t(all_samples))
  condmvn <- condMVNorm::condMVN(mean = mu_tilde, sigma = var_tilde,
                                 dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                                 X.given = c(samples$theta_mu[,groups[s],iteration], log(diag(samples$theta_var[,,iteration]))))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

filtered_samples_diag_groups <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[,,filter],
    theta_var = sampler$samples$theta_var[, , filter],
    alpha = sampler$samples$alpha[, , filter],
    iteration = length(filter),
    group = sampler$group_idx
  )
}

variant_funs <- list(
  sample_store = sample_store_diag_groups,
  add_info = add_info_diag_groups,
  get_startpoints = get_startpoints_diag_groups,
  fill_samples = fill_samples_diag_groups,
  gibbs_step = gibbs_step_diag_groups,
  get_group_level = get_group_level_diag_groups,
  get_conditionals = get_conditionals_diag_groups,
  filtered_samples = filtered_samples_diag_groups
)
