source("samplers/pmwg/sampling.R")
source("samplers/pmwg/variants/standard.R")

add_info_diag <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, sampler$n_pars), theta_mu_var = rep(1, sampler$n_pars))
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- 1/prior$theta_mu_var #Inverse of the prior
  
  #Hyper parameters
  attr(sampler, "v_half") <- 2
  attr(sampler, "A_half") <- 1
  sampler$prior <- prior
  return(sampler)
}

get_startpoints_diag <- function(pmwgs, start_mu, start_var){
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$n_pars, sd = 1)
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- diag(1/rgamma(pmwgs$n_pars, 10, 5)) #Bit stupid maybe as startpoint
  start_a_half <- 1 / rgamma(n = pmwgs$n_pars, shape = 2, rate = 1)
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var), a_half = start_a_half))
}

get_conditionals_diag <- function(s, samples, n_pars){
  iteration <- samples$iteration
  pts2_unwound <- log(apply(samples$theta_var,3,diag))
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- var(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[,iteration], log(diag(samples$theta_var[,,iteration]))))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
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

variant_funs <- list(
  sample_store = sample_store_standard,
  add_info = add_info_diag,
  get_startpoints = get_startpoints_diag,
  fill_samples = fill_samples_standard,
  gibbs_step = gibbs_step_diag,
  get_group_level = get_group_level_standard,
  get_conditionals = get_conditionals_diag,
  filtered_samples = filtered_samples_standard
)

