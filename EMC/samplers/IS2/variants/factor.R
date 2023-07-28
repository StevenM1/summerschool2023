source("samplers/IS2/IS2.R")

get_all_pars_factor <- function(samples, idx, info){
  n_subjects <- samples$n_subjects
  n_iter = length(samples$samples$stage[idx])
  # Extract relevant objects
  alpha <- samples$samples$alpha[,,idx]
  theta_mu <- samples$samples$theta_mu[,idx]
  lambda <- samples$samples$lambda_untransf[,,idx, drop = F]
  psi_inv <- samples$samples$theta_psi_inv[,,idx, drop = F]
  sig_err_inv <- samples$samples$theta_sig_err_inv[,,idx]
  
  constraintMat <- info$hyper$constraintMat
  n_factors <- samples$n_factors
  
  lambda.unwound <- apply(lambda,3,unwind_lambda, constraintMat)
  sig_err_inv.diag <- log(apply(sig_err_inv, 3, diag))
  psi_inv.diag <- matrix(log(apply(psi_inv, 3, diag)), nrow = n_factors)
  
  # Set up
  n_params<- nrow(alpha) + nrow(theta_mu) + nrow(sig_err_inv.diag) + nrow(psi_inv.diag) + nrow(lambda.unwound)
  all_samples=array(dim=c(n_subjects,n_params,n_iter))
  mu_tilde=array(dim = c(n_subjects,n_params))
  var_tilde=array(dim = c(n_subjects,n_params,n_params))
  
  for (j in 1:n_subjects){
    all_samples[j,,] = rbind(alpha[,j,],theta_mu[,],sig_err_inv.diag[,],psi_inv.diag[,],lambda.unwound[,])
    # calculate the mean for re, mu and sigma
    mu_tilde[j,] =apply(all_samples[j,,],1,mean)
    # calculate the covariance matrix for random effects, mu and sigma
    var_tilde[j,,] = cov(t(all_samples[j,,]))
  }
  
  for(i in 1:n_subjects){ #RJI_change: this bit makes sure that the sigma tilde is pos def
    if(!corpcor::is.positive.definite(var_tilde[i,,], tol=1e-8)){
      var_tilde[i,,]<-corpcor::make.positive.definite(var_tilde[i,,], tol=1e-6)
    }
  }
  
  X <- cbind(t(theta_mu),t(sig_err_inv.diag),t(psi_inv.diag), t(lambda.unwound)) 
  info$n_params <- n_params
  info$n_factors <- n_factors
  info$given.ind <- (info$n_randeffect+1):n_params
  info$X.given_ind <- 1:length(info$given.ind)
  return(list(X = X, mu_tilde = mu_tilde, var_tilde = var_tilde, info = info))
}

unwind_lambda <- function(lambda, constraintMat, n_factors = NULL, n_randeffect = NULL, reverse = F){
  if(reverse){
    out <- matrix(0, n_randeffect, n_factors)
    out[constraintMat] <- lambda
  } else{
    out <- as.numeric(lambda[constraintMat])
  }
  return(out)
}

group_dist_factor = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  n_randeffect <- info$n_randeffect
  n_factors <- info$n_factors
  param.theta_mu <- parameters[1:n_randeffect]
  param.sig_err_inv <- exp(parameters[(n_randeffect+1):(n_randeffect + n_randeffect)])
  param.psi_inv <- exp(parameters[(n_randeffect+n_randeffect+1):(n_randeffect + n_randeffect+ n_factors)])
  param.lambda.unwound <- parameters[(n_randeffect+n_randeffect+n_factors+1):length(parameters)] 
  param.lambda <- unwind_lambda(param.lambda.unwound, info$hyper$constraintMat, n_factors, n_randeffect, reverse = T)
  param.var <- param.lambda %*% diag(1/param.psi_inv, length(param.psi_inv)) %*% t(param.lambda) + diag(1/param.sig_err_inv)
  if (sample){
    return(mvtnorm::rmvnorm(n_samples, param.theta_mu, param.var))
  }else{
    logw_second<-max(-5000*info$n_randeffect, mvtnorm::dmvnorm(random_effect, param.theta_mu,param.var,log=TRUE))
    return(logw_second)
  }
}

prior_dist_factor = function(parameters, info){ 
  n_randeffect <- info$n_randeffect
  n_factors <- info$n_factors
  prior <- info$prior
  hyper <- info$hyper
  #Extract and when necessary transform back
  param.theta_mu <- parameters[1:n_randeffect]
  param.sig_err_inv <- exp(parameters[(n_randeffect+1):(n_randeffect + n_randeffect)])
  param.psi_inv <- exp(parameters[(n_randeffect+n_randeffect+1):(n_randeffect + n_randeffect+ n_factors)])
  param.lambda.unwound <- parameters[(n_randeffect+n_randeffect+n_factors+1):length(parameters)] 
  
  log_prior_mu=sum(dnorm(param.theta_mu, mean = prior$theta_mu_mean, sd = sqrt(prior$theta_mu_var), log =TRUE))
  log_prior_sig_err_inv = sum(pmax(-100, dgamma(param.sig_err_inv, shape = hyper$nu/2, rate = (hyper$s2*hyper$nu)/2, log=TRUE)))
  log_prior_psi_inv = sum(pmax(-100, dgamma(param.psi_inv, shape = hyper$al, rate = hyper$bl, log=TRUE)))
  log_prior_lambda=sum(dnorm(param.lambda.unwound, mean = 0, sd = sqrt(prior$theta_lambda_var), log =TRUE))
  
  jac_sig_err_inv <- -sum(log(param.sig_err_inv)) # Jacobian determinant of transformation of log of the sig_err_inv  
  jac_psi_inv <- -sum(log(param.psi_inv)) # Jacobian determinant of transformation of log of the psi_inv
  # Jacobians are actually part of the denominator (dnorm(prop_theta)) since transformations of the data (rather than parameters),
  # warrant a jacobian added. But we add the jacobians here for ease of calculations. 
  return(log_prior_mu + log_prior_sig_err_inv + log_prior_psi_inv + log_prior_lambda - jac_psi_inv - jac_sig_err_inv)
}

variant_funs <- list(
  get_all_pars = get_all_pars_factor,
  prior_dist = prior_dist_factor,
  group_dist = group_dist_factor
)
