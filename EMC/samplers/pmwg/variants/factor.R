source("samplers/pmwg/sampling.R")
source("samplers/pmwg/variants/standard.R")

add_info_factor <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  args <- list(...)
  n_factors <- args$n_factors
  constraintMat <- args$constraintMat
  n_pars <- sampler$n_pars
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, n_pars), 
                  theta_mu_var = rep(1, n_pars),
                  theta_lambda_var = 1)
  }
  # Things I save rather than re-compute inside the loops.
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- diag(1/prior$theta_mu_var) 
  prior$theta_lambda_invar <-1/prior$theta_lambda_var
  
  if(is.null(constraintMat)){
    constraintMat <- matrix(1, nrow = n_pars, ncol = n_factors)
    constraintMat[upper.tri(constraintMat, diag = T)] <- 0 #Now you can't fix one of the diagonal values to 0
  }
  signFix <- F
  constraintMat <- constraintMat != 0 #For indexing
  if(any(diag(constraintMat) != 0)) signFix <- T
  
  #Hyper parameters
  attr(sampler, "al") <- 1
  attr(sampler, "bl") <- 1/2
  attr(sampler, "nu") <- 2
  attr(sampler, "s2") <- 1/2
  attr(sampler, "signFix") <- signFix
  attr(sampler, "constraintMat") <- constraintMat
  
  sampler$prior <- prior
  sampler$n_factors <- n_factors
  return(sampler)
}

sample_store_factor <- function(data, par_names, iters = 1, stage = "init", integrate = T, ...) {
  n_factors <- list(...)$n_factors
  subject_ids <- unique(data$subject)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    theta_lambda = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, NULL, NULL)),
    lambda_untransf = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, NULL, NULL)),
    theta_sig_err_inv = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    theta_psi_inv = array(NA_real_, dim = c(n_factors, n_factors, iters), dimnames = list(NULL, NULL, NULL)),
    theta_eta = array(NA_real_, dim = c(n_subjects, n_factors, iters), dimnames = list(subject_ids, NULL, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

get_startpoints_factor<- function(pmwgs, start_mu, start_var){
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$n_pars, sd = 1)
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
  start_psi_inv <- diag(1, pmwgs$n_factors)
  start_sig_err_inv <- diag(1, pmwgs$n_pars)
  start_lambda <- matrix(0, nrow = pmwgs$n_pars, ncol = pmwgs$n_factors)
  start_lambda[1:pmwgs$n_factors, 1:pmwgs$n_factors] <- diag(1, pmwgs$n_factors)
  start_eta <- matrix(0, nrow = pmwgs$n_subjects, ncol = pmwgs$n_factors)
  return(list(tmu = start_mu, tvar = start_var, lambda = start_lambda, lambda_untransf = start_lambda,
              sig_err_inv = start_sig_err_inv, psi_inv = start_psi_inv,
              eta = start_eta))
}

fill_samples_factor <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$theta_lambda[,,j] <- group_level$lambda
  samples$lambda_untransf[,,j] <- group_level$lambda_untransf
  samples$theta_sig_err_inv[,,j] <- group_level$sig_err_inv
  samples$theta_psi_inv[,,j] <- group_level$psi_inv
  samples$theta_eta[,,j] <- group_level$eta
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}

gibbs_step_factor <- function(sampler, alpha){
  # Gibbs step for group means with parameter expanded factor analysis from Ghosh & Dunson 2009
  # mu = theta_mu, var = theta_var
  last <- last_sample_factor(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  
  #extract previous values (for ease of reading)
  
  alpha <- t(alpha)
  n_subjects <- sampler$n_subjects
  n_pars <- sampler$n_pars
  n_factors <- sampler$n_factors
  constraintMat <- hyper$constraintMat
  
  eta <- matrix(last$eta, n_subjects, n_factors)
  psi_inv <- matrix(last$psi_inv, n_factors)
  sig_err_inv <- last$sig_err_inv
  lambda <- matrix(last$lambda, n_pars, n_factors)
  mu <- last$mu
  
  #Update mu
  mu_sig <- solve(n_subjects * sig_err_inv + prior$theta_mu_invar)
  mu_mu <- mu_sig %*% (sig_err_inv %*% colSums(alpha - eta %*% t(lambda)) + prior$theta_mu_invar%*% prior$theta_mu_mean)
  mu <- rmvnorm(1, mu_mu, mu_sig)
  colnames(mu) <- colnames(alpha)
  # calculate mean-centered observations
  alphatilde <- sweep(alpha, 2, mu)
  
  #Update eta, I do this one first since I don't want to save eta
  eta_sig <- solve(psi_inv + t(lambda) %*% sig_err_inv %*% lambda)
  eta_mu <- eta_sig %*% t(lambda) %*% sig_err_inv %*% t(alphatilde)
  eta[,] <- t(apply(eta_mu, 2, FUN = function(x){rmvnorm(1, x, eta_sig)}))
  
  #Update sig_err
  sig_err_inv <- diag(rgamma(n_pars,shape=(hyper$nu+n_subjects)/2, rate=(hyper$nu*hyper$s2+ colSums((alphatilde - eta %*% t(lambda))^2))/2))
  
  #Update lambda
  for (j in 1:n_pars) {
    constraint <- constraintMat[j,] #T if item is not constraint (bit confusing tbh)
    if(any(constraint)){ #Don't do this if there are no free entries in lambda
      etaS <- eta[,constraint]
      lambda_sig <- solve(sig_err_inv[j,j] * t(etaS) %*% etaS + prior$theta_lambda_invar * diag(1,sum(constraint)))
      lambda_mu <- (lambda_sig * sig_err_inv[j,j]) %*% (t(etaS) %*% alphatilde[,j])
      lambda[j,constraint] <- rmvnorm(1,lambda_mu,lambda_sig)
    }
  }
  
  #Update psi_inv
  psi_inv[,] <- diag(rgamma(n_factors ,shape=(hyper$al+n_subjects)/2,rate=hyper$bl+colSums(eta^2)/2), n_factors)
  
  lambda_orig <- lambda
  #If the diagonals of lambda aren't constrained to be 1, we should fix the signs
  if(hyper$signFix){
    for(l in 1:n_factors){
      mult <- ifelse(lambda[l, l] < 0, -1, 1) #definitely a more clever function for this 
      lambda_orig[,l] <- mult * lambda[, l]
    }
  }
  
  var <- lambda_orig %*% solve(psi_inv) %*% t(lambda_orig) + diag(1/diag((sig_err_inv)))
  lambda_orig <- lambda_orig %*% matrix(diag(sqrt(1/diag(psi_inv)), n_factors), nrow = n_factors)
  return(list(tmu = mu, tvar = var, lambda_untransf = lambda, lambda = lambda_orig, eta = eta, 
              sig_err_inv = sig_err_inv, psi_inv = psi_inv, alpha = t(alpha)))
}

last_sample_factor <- function(store) {
  list(
    mu = store$theta_mu[, store$idx],
    eta = store$theta_eta[,,store$idx],
    lambda = store$lambda_untransf[,,store$idx],
    psi_inv = store$theta_psi_inv[,,store$idx],
    sig_err_inv = store$theta_sig_err_inv[,,store$idx]
  )
}

get_conditionals_factor <- function(s, samples, n_pars){
  iteration <- samples$iteration
  sig_err <- log(apply(samples$theta_sig_err_inv,3,diag))
  psi <- log(apply(samples$theta_psi_inv,3,diag))
  eta <- matrix(samples$theta_eta[s,,], nrow = samples$n_factors)
  lambda <- apply(samples$lambda_untransf, 3, unwind_lambda, samples$constraintMat, samples$n_factors)
  theta_mu <- samples$theta_mu 
  all_samples <- rbind(samples$alpha[, s,],theta_mu, eta, sig_err, psi, lambda)#, sig_err, psi, lambda)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(theta_mu[,iteration],
                                 samples$theta_eta[s,,iteration],
                                 log(diag(samples$theta_sig_err_inv[,, iteration])),
                                 log(apply(samples$theta_psi_inv[,,iteration, drop = F], 3, diag)),
                                 unwind_lambda(samples$lambda_untransf[,, iteration], samples$constraintMat, samples$n_factors)))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

unwind_lambda <- function(lambda, constraintMat, n_factors, reverse = F){
  if(reverse){
    out <- matrix(0, n_randeffect, n_factors)
    out[constraintMat] <- lambda
  } else{
    out <- as.numeric(lambda[constraintMat])
  }
  return(out)
}

filtered_samples_factor <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter],
    lambda_untransf = sampler$samples$lambda_untransf[, , filter, drop = F],
    theta_psi_inv = sampler$samples$theta_psi_inv[, , filter, drop = F],
    theta_sig_err_inv = sampler$samples$theta_sig_err_inv[, , filter],
    theta_eta = sampler$samples$theta_eta[, , filter, drop = F],
    theta_var = sampler$samples$theta_var[,,filter],
    alpha = sampler$samples$alpha[, , filter],
    constraintMat = attributes(sampler)$constraintMat,
    n_factors = sampler$n_factors,
    iteration = length(filter)
  )
}

variant_funs <- list(
  sample_store = sample_store_factor,
  add_info = add_info_factor,
  get_startpoints = get_startpoints_factor,
  fill_samples = fill_samples_factor,
  gibbs_step = gibbs_step_factor,
  get_group_level = get_group_level_standard,
  get_conditionals = get_conditionals_factor,
  filtered_samples = filtered_samples_factor
)

