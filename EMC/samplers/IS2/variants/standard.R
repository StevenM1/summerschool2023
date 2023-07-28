source("samplers/IS2/IS2.R")

get_all_pars_standard <- function(samples, idx, info){
  n_subjects <- samples$n_subjects
  n_iter = length(samples$samples$stage[idx])
  # Exctract relevant objects
  alpha <- samples$samples$alpha[,,idx]
  theta_mu <- samples$samples$theta_mu[,idx]
  theta_var <- samples$samples$theta_var[,,idx]
  a_half <- log(samples$samples$a_half[,idx])
  theta_var.unwound = apply(theta_var,3,unwind)
  # Set up
  n_params<- samples$n_pars+nrow(theta_var.unwound)+samples$n_pars
  all_samples=array(dim=c(n_subjects,n_params,n_iter))
  mu_tilde=array(dim = c(n_subjects,n_params))
  var_tilde=array(dim = c(n_subjects,n_params,n_params))
  
  for (j in 1:n_subjects){
    all_samples[j,,] = rbind(alpha[,j,],theta_mu[,],theta_var.unwound[,])
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
  X <- cbind(t(theta_mu),t(theta_var.unwound),t(a_half))
  info$n_params <- n_params
  info$given.ind <- (info$n_randeffect+1):n_params
  info$X.given_ind <- 1:(n_params-info$n_randeffect)
  return(list(X = X, mu_tilde = mu_tilde, var_tilde = var_tilde, info = info))
}

robust_diwish = function (W, v, S) { #RJI_change: this function is to protect against weird proposals in the diwish function, where sometimes matrices weren't pos def
  if (!is.matrix(S)) S <- matrix(S)
  if (!is.matrix(W)) W <- matrix(W)
  p <- nrow(S)
  gammapart <- sum(lgamma((v + 1 - 1:p)/2))
  ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
  if (corpcor::is.positive.definite(W, tol=1e-8)){
    cholW<-base::chol(W)
  }else{
    return(1e-10)
  }
  if (corpcor::is.positive.definite(S, tol=1e-8)){
    cholS <- base::chol(S)
  }else{
    return(1e-10)
  }
  halflogdetS <- sum(log(diag(cholS)))
  halflogdetW <- sum(log(diag(cholW)))
  invW <- chol2inv(cholW)
  exptrace <- sum(S * invW)
  lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * exptrace
  lpdf <- lnum - ldenom
  out <- exp(lpdf)
  if(!is.finite(out)) return(1e-100)
  if(out < 1e-10) return(1e-100)
  return(exp(lpdf))
}

unwind=function(x,reverse=FALSE) {
  
  if (reverse) {
    n=sqrt(2*length(x)+0.25)-0.5 ## Dim of matrix.
    out=array(0,dim=c(n,n))
    out[lower.tri(out,diag=TRUE)]=x
    diag(out)=exp(diag(out))
    out=out%*%t(out)
  } else {
    y=t(base::chol(x))
    diag(y)=log(diag(y))
    out=y[lower.tri(y,diag=TRUE)]
  }
  return(out)
}


group_dist_standard = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
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

prior_dist_standard = function(parameters, info){
  n_randeffect <- info$n_randeffect
  prior <- info$prior
  hyper <- info$hyper
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind(param.theta.sig.unwound, reverse = TRUE)
  param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
  log_prior_mu=mvtnorm::dmvnorm(param.theta.mu, mean = prior$theta_mu_mean, sigma = prior$theta_mu_var, log =TRUE)
  log_prior_sigma = log(robust_diwish(param.theta.sig2, v=hyper$v_half+ n_randeffect-1, S = 2*hyper$v_half*diag(1/param.a))) 
  log_prior_a = sum(logdinvGamma(param.a,shape = 1/2,rate=1/(hyper$A_half^2)))
  # These are Jacobian corrections for the transformations on these
  logw_den2 <- -sum(log(param.a))
  logw_den3 <- -(log(2^n_randeffect)+sum((n_randeffect:1+1)*log(diag(param.theta.sig2)))) 
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
  prior_dist = prior_dist_standard,
  group_dist = group_dist_standard
)

