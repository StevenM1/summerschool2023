source("samplers/pmwg/variants/factor.R")

add_info_factor_groups <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  sampler <- add_info_factor(sampler, prior, ...)
  sampler$par_idx <- list(...)$par_idx
  sampler$parGroups <- unique(as.vector(par_idx))
  sampler$group_idx <- aggregate(group ~ subject, sampler$data, mean)[,2]
  return(sampler)
}

sample_store_factor_groups <- function(data, par_names, iters = 1, stage = "init", integrate = T, ...) {
  n_pars <- length(par_names)
  groups <- unique(data$group)
  n_groups <- length(groups)
  samples <- sample_store_factor(data, par_names, iters, stage, integrate, ...)
  samples$theta_mu = array(NA_real_,dim = c(n_pars, n_groups, iters), dimnames = list(par_names, groups, NULL))
  return(samples)
}

get_startpoints_factor_groups<- function(pmwgs, start_mu, start_var){
  startpoints <- get_startpoints_factor(pmwgs, start_mu, start_var)
  startpoints$tmu <- replicate(max(pmwgs$group_idx), startpoints$tmu)
  startpoints$group_idx <- pmwgs$group_idx
  return(startpoints)
}

fill_samples_factor_groups <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  # Factor side
  samples$theta_lambda[,,j] <- group_level$lambda
  samples$lambda_untransf[,,j] <- group_level$lambda_untransf
  samples$theta_sig_err_inv[,,j] <- group_level$sig_err_inv
  samples$theta_psi_inv[,,j] <- group_level$psi_inv
  samples$theta_eta[,,j] <- group_level$eta
  
  # group level mean and implied variance
  samples$theta_mu[,, j] <- group_level$tmu
  samples$theta_var[, , j] <- group_level$tvar
  
  # Random effects
  samples <- fill_samples_RE(samples, proposals, epsilon,j, n_pars)
  return(samples)
}

gibbs_step_factor_groups <- function(sampler, alpha){
  # Gibbs step for group level with multiple means, 
  # based on parameter expanded factor analysis from Ghosh & Dunson 2009
  # mu = theta_mu, var = theta_var
  last <- last_sample_factor_groups(sampler$samples)
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
  alphatilde <- alpha
  #Update mu
  for(parGroup in sampler$parGroups){
    par_idx <- unique(which(sampler$par_idx == parGroup, arr.ind = T)[,1])
    parGroup_num <- as.numeric(strsplit(parGroup,",")[[1]])
    group_idx <- which(sampler$group_idx %in% parGroup_num)
    mu_sig <- solve(length(group_idx) * sig_err_inv[par_idx, par_idx, drop = F] + prior$theta_mu_invar[par_idx, par_idx, drop = F])
    mu_mu <- mu_sig %*% (sig_err_inv[par_idx, par_idx, drop = F] %*% colSums(alpha[group_idx,par_idx, drop = F]
                                                                             - eta[group_idx,,drop = F] %*% t(lambda[par_idx,, drop = F])) 
                         + prior$theta_mu_invar[par_idx, par_idx,drop = F] %*% prior$theta_mu_mean[par_idx])
    mu[par_idx,parGroup_num] <- rmvnorm(1, mu_mu, mu_sig)
  }
  for (group in unique(sampler$group_idx)){
    group_idx <- which(sampler$group_idx == group)
    alphatilde[group_idx,] <- sweep(alphatilde[group_idx,], 2, mu[,group]) # This will break if groups != 1:length(groups)
  }
  
  
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
              sig_err_inv = sig_err_inv, psi_inv = psi_inv, alpha = t(alpha), group_idx = sampler$group_idx))
}

last_sample_factor_groups <- function(store) {
  list(
    mu = store$theta_mu[,, store$idx],
    eta = store$theta_eta[,,store$idx],
    lambda = store$lambda_untransf[, , store$idx],
    psi_inv = store$theta_psi_inv[,,store$idx],
    sig_err_inv = store$theta_sig_err_inv[,,store$idx]
  )
}

get_group_level_factor_groups <- function(parameters, s){
  mu <- parameters$tmu[,parameters$group_idx[s]]
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

get_conditionals_factor <- function(s, samples, n_pars){
  groups <- samples$groups
  iteration <- samples$iteration
  sig_err <- log(apply(samples$theta_sig_err_inv,3,diag))
  psi <- log(apply(samples$theta_psi_inv,3,diag))
  eta <- matrix(samples$theta_eta[s,,], nrow = samples$n_factors)
  lambda <- apply(samples$lambda_untransf, 3, unwind_lambda, samples$constraintMat, samples$n_factors)
  theta_mu <- samples$theta_mu[,groups[s],] 
  all_samples <- rbind(samples$alpha[, s,],theta_mu, eta, sig_err, psi, lambda)#, sig_err, psi, lambda)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(theta_mu[,iteration],
                                 samples$theta_mu[,groups[s],iteration],
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

get_conditionals_factor_groups <- function(s, samples, n_pars){
  groups <- samples$groups
  iteration <- samples$iteration
  pts2_unwound <- apply(samples$theta_var,3,unwind)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu[,groups[s],],pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  var_tilde <- stats::var(t(all_samples))
  condmvn <- condMVNorm::condMVN(mean = mu_tilde, sigma = var_tilde,
                                 dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                                 X.given = c(samples$theta_mu[,groups[s],iteration], unwind(samples$theta_var[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
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
    iteration = length(filter),
    group = sampler$group_idx
  )
}

variant_funs <- list(
  sample_store = sample_store_factor_groups,
  add_info = add_info_factor_groups,
  get_startpoints = get_startpoints_factor_groups,
  fill_samples = fill_samples_factor_groups,
  gibbs_step = gibbs_step_factor_groups,
  get_group_level = get_group_level_factor_groups,
  get_conditionals = get_conditionals_factor_groups,
  filtered_samples = filtered_samples_factor_groups
)
