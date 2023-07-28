source("samplers/IS2/IS2.R")

unwind <- function(x,reverse=FALSE, diag = TRUE) {
  if (reverse) {
    if(diag){
      out <- diag(exp(x), nrow = length(x))
    } else{
      out <- exp(x)
    }
  } else {
    out <- log(diag(x))
  }
  return(out)
}

group_dist_diag = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  n_randeffect <- info$n_randeffect
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)] 
  param.theta.sig2 <- unwind(param.theta.sig.unwound, reverse = TRUE)
  if (sample){
    return(mvtnorm::rmvnorm(n_samples, param.theta.mu,param.theta.sig2))
  }else{
    logw_second<-max(-5000*info$n_randeffect, mvtnorm::dmvnorm(random_effect, param.theta.mu,param.theta.sig2,log=TRUE))
    return(logw_second)
  }
}

prior_dist_diag = function(parameters, info){
  n_randeffect <- info$n_randeffect
  prior <- info$prior
  hyper <- info$hyper
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind(param.theta.sig.unwound, reverse = TRUE, diag = FALSE)
  param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
  log_prior_mu=mvtnorm::dmvnorm(param.theta.mu, mean = prior$theta_mu_mean, sigma = prior$theta_mu_var, log =TRUE)
  log_prior_sigma = sum(logdinvGamma(param.theta.sig2, shape = hyper$v_half/2, rate = hyper$v_half/param.a))
  log_prior_a = sum(logdinvGamma(param.a,shape = 1/2,rate=1/(hyper$A_half^2)))
  # These are Jacobian corrections for the transformations on these
  logw_den2 <- -sum(log(param.a))
  logw_den3 <- -sum(log(param.theta.sig2))
  return(log_prior_mu + log_prior_sigma + log_prior_a - logw_den3 - logw_den2)
}

logdinvGamma <- function(x, shape, rate){
  alpha <- shape
  beta <- 1/rate
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 
                                                        1) * log(x) - (beta/x)
  return(pmax(log.density, -50)) #Roughly equal to 1e-22 on real scale
}

variant_funs <- list(
  get_all_pars = get_all_pars_standard,
  prior_dist = prior_dist_diag,
  group_dist = group_dist_diag
)

