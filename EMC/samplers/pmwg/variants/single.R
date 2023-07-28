library(mvtnorm) ## For the multivariate normal.
library(MASS) ## For matrix inverse.
library(MCMCpack) #For inverse wishart
library(lme4)
library(parallel)
library(matrixcalc)

source("samplers/pmwg/messaging.R")

pmwgs <- function(dadm, pars = NULL, ll_func = NULL, prior = NULL) {
  ###Gets/sets priors, creates pmwgs object and stores essentials
  if(is.data.frame(dadm)) dadm <- list(dadm)
  dadm <- extractDadms(dadm)
  if(is.null(pars)) pars <- dadm$pars
  if(is.null(ll_func)) ll_func <- dadm$ll_func
  if(is.null(prior)) prior <- dadm$prior
  dadm_list <-dadm$dadm_list
  
  # Descriptives
  n_pars <- length(pars)
  subjects <- as.numeric(unique(dadm$subjects))
  n_subjects <- length(subjects)
  
  # Storage for the samples.
  samples <- sample_store(pars, subjects)
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, n_pars), theta_mu_var = diag(rep(1, n_pars)))
  }
  if(is.null(dim(prior$theta_mu_var))){
    prior$theta_mu_var <- diag(prior$theta_mu_var, n_pars)
  }
  sampler <- list(
    data = dadm_list,
    par_names = pars,
    n_pars = n_pars,
    n_subjects = n_subjects,
    subjects = subjects,
    prior = prior,
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler
}

sample_store <- function(par_names, subject_ids, iters = 1, stage = "init", ...) {
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  list(
    epsilon = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    origin = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL))
  )
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 verbose = FALSE, particles = 1000, n_cores = 1, epsilon = NULL) {
  # Create and fill initial random effects for each subject
  likelihoods <- array(NA_real_, dim = c(pmwgs$n_subjects))
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals, pmwgs$prior$theta_mu_mean, pmwgs$prior$theta_mu_var, 
                          particles, pmwgs, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals, pmwgs$prior$theta_mu_mean, pmwgs$prior$theta_mu_var,
                        particles, pmwgs = pmwgs)
  }
  
  if(is.null(epsilon)) epsilon <- rep(set_epsilon(pmwgs$n_pars, verbose), pmwgs$n_subjects)
  if(length(epsilon) == 1) epsilon <- rep(epsilon, pmwgs$n_subjects)
  
  proposals <- simplify2array(proposals)
  pmwgs$init <- TRUE
  pmwgs$samples$alpha[, , 1] <- do.call(cbind, proposals[1,])
  pmwgs$samples$subj_ll[, 1] <- unlist(proposals[2,])
  pmwgs$samples$origin[,1] <- 2
  pmwgs$samples$idx <- 1
  pmwgs$samples$epsilon[,1] <- epsilon
  return(pmwgs)
}

start_proposals <- function(s, start_mu, start_var, n_particles, pmwgs){
  n_particles <- n_particles
  proposals <- particle_draws(n_particles, start_mu, start_var)
  colnames(proposals) <- pmwgs$par_names # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,dadm = pmwgs$data[[which(pmwgs$subjects == s)]])
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}

new_particle_single <- function (s, data, num_particles, prior_mu, prior_var, prev_sample,
                                 eff_mu = NULL, eff_var = NULL, mix_proportion = c(0.5, 0.5, 0), 
                                 likelihood_func = NULL, epsilon = NULL, subjects) 
{
  num_particles <- num_particles[s]
  eff_mu <- prev_sample[,s] #eff_mu[, s]
  eff_var <- eff_var[,,s]
  subj_mu <- prev_sample[,s]
  particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
  cumuNumbers <- cumsum(particle_numbers)
  pop_particles <- particle_draws(particle_numbers[1], prior_mu, 
                                  prior_var)
  ind_particles <- particle_draws(particle_numbers[2], subj_mu, 
                                  prior_var * epsilon[s]^2)
  if(mix_proportion[3] == 0){
    eff_particles <- NULL
  } else{
    eff_particles <- particle_draws(particle_numbers[3], subj_mu, eff_var*epsilon[s]^2)
  }
  proposals <- rbind(pop_particles, ind_particles, eff_particles)
  colnames(proposals) <- names(subj_mu)
  proposals[1, ] <- subj_mu
  lw <- apply(proposals, 1, likelihood_func, dadm = data[[which(subjects == s)]])
  lp <- mvtnorm::dmvnorm(x = proposals, mean = prior_mu, sigma = prior_var, 
                         log = TRUE)
  prop_density <- mvtnorm::dmvnorm(x = proposals, mean = subj_mu, 
                                   sigma = prior_var * (epsilon[s]^2))
  if (mix_proportion[3] == 0) {
    eff_density <- 0
  }
  else {
    eff_density <- mvtnorm::dmvnorm(x = proposals, mean = subj_mu, sigma = eff_var*epsilon[s]^2)
  }
  lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * 
                                             prop_density) + (mix_proportion[3] * eff_density))
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  idx <- sample(x = num_particles, size = 1, prob = weights)
  origin <- min(which(idx <= cumuNumbers))
  return(list(proposal = proposals[idx, ], ll = lw[idx], origin = origin))
}

run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 1000,
                      verbose = TRUE,
                      n_cores = 1,
                      epsilon = NULL,
                      pstar = NULL,
                      n_window = 50,
                      epsilon_upper_bound = 20,
                      mix = NULL, 
                      ...) {
  # Set defaults for NULL values
  if(is.null(mix)) mix <- set_mix(stage, verbose)
  # Set stable (fixed) new_sample argument for this run
  n_pars <- length(pmwgs$par_names)
  # Display stage to screen
  
  alphaStar=-qnorm(pstar/2) #Idk about this one
  n0=round(5/(pstar*(1-pstar))) #Also not questioning this math for now
  if(is.null(epsilon)){
    epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
  }  
  if(length(epsilon) == 1){
    epsilon <- rep(epsilon, pmwgs$n_subjects)
  }
  if(length(particles) == 1){
    particles <- rep(particles, pmwgs$n_subjects)
  }
  
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  if (verbose) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx
  
  prior_mu <- pmwgs$prior$theta_mu_mean
  prior_var <- pmwgs$prior$theta_mu_var
  data <- pmwgs$data
  subjects <- pmwgs$subjects
  n_subjects <- pmwgs$n_subjects
  eff_mu <- NULL
  eff_var <- replicate(n_subjects, prior_var)
  # Main iteration loop
  for (i in 1:iter) {
    if (verbose) {
      update_progress_bar(pb, i, extra = mean(accept_rate(pmwgs)))
    }
    j <- start_iter + i
    if(j == n_window + 1){
      mix <- set_mix("sample", verbose = F)
      epsilon <- epsilon * 2
    }
    if(j > n_window & (j %% 10 == 0 || i == 1 || j == n_window)){
      for(k in 1:n_subjects){
        # I outline quicker ways to do this below, but for now this is easier for testing
        covMat   <- cov(t(pmwgs$samples$alpha[,k,(j-1-n_window):(j-1)]))
        if(!is.positive.semi.definite(covMat)){
          save(covMat, file = "wrongCov.RData")
          covMat <- as.matrix(nearPD(covMat)$mat)
        }
        eff_var[,,k] <- covMat
      }
      # eff_mu <- apply(pmwgs$samples$alpha[,,(j-1-n_window):(j-1)], 1:2, mean)
    }
    
    prev_sample <- matrix(pmwgs$samples$alpha[,,j-1], nrow = n_pars, ncol = n_subjects)
    rownames(prev_sample) <- pmwgs$par_names
    proposals=mclapply(X=1:pmwgs$n_subjects,FUN = new_particle_single, data, particles, prior_mu, prior_var,
                       prev_sample, eff_mu, eff_var, mix, pmwgs$ll_func, epsilon, subjects, mc.cores = n_cores)
    proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))
    alpha <- matrix(proposals[1:n_pars,], nrow = n_pars, ncol = n_subjects)
    ll <- proposals[n_pars + 1,]
    origin <- proposals[n_pars + 2,]
    
    pmwgs$samples$alpha[, , j] <- alpha
    pmwgs$samples$idx <- j
    pmwgs$samples$subj_ll[, j] <- ll
    pmwgs$samples$origin[,j] <- origin
    
    if(!is.null(pstar)){
      if(j > (n_window + n0)){
        acc <-  pmwgs$samples$alpha[1,,j] != pmwgs$samples$alpha[1,,(j-1)]
        epsilon<-pmin(update.epsilon(epsilon^2, acc, pstar, j, pmwgs$n_pars, alphaStar), 15)
      }
    }
    
    pmwgs$samples$epsilon[,j] <- epsilon
  }
  if (verbose) close(pb)
  return(pmwgs)
}

particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  }
  return(rmvnorm(n, mu, covar))
}

update.epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c=((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta=log(sqrt(epsilon2))
  Theta=Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}

set_mix <- function(stage, verbose) {
  if (stage == "sample") {
    mix <- c(0.3, 0.2, 0.5)
  } else {
    mix <- c(0.3, 0.7, 0)
  }
  if(verbose) message(sprintf("mix has been set to c(%s) based on the stage being run",  paste(mix, collapse = ", ")))
  return(mix)
}

set_epsilon <- function(n_pars, verbose) {
  if (n_pars > 15) {
    epsilon <- 0.1
  } else if (n_pars > 10) {
    epsilon <- 0.3
  } else {
    epsilon <- 0.5
  }
  if(verbose) message(sprintf("Epsilon has been set to %.1f based on number of parameters",epsilon))
  return(epsilon)
}

numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
  numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_proportion)
  if (mix_proportion[3] == 0) {
    numbers[3] <- 0
    numbers[2] <- num_particles - numbers[1]
  } else {
    numbers[3] <- num_particles - sum(numbers)
  }
  return(numbers)
}

extend_sampler <- function(sampler, n_samples, stage) {
  # This function takes the sampler and extends it along the intended number of
  # iterations, to ensure that we're not constantly increasing our sampled object
  # by 1. 
  old <- sampler$samples
  par_names <- sampler$par_names
  subject_ids <- sampler$subjects
  start <- old$idx + 1
  end <- old$idx + n_samples
  for(obj_name in names(sampler$samples)){
    obj <- sampler$samples[[obj_name]]
    dimensions <- length(dim(obj))
    if(dimensions > 1 & !obj_name %in% c("sum_samples")){ #Ok not the cleanest
      new_dim <- c(rep(0, (dimensions -1)), n_samples)
      extension <- array(NA_real_, dim = dim(obj) +  new_dim, dimnames = dimnames(obj))
      if(dimensions == 2){ #There must be a cleaner way to to do this
        extension[, -(start:end)] <- obj
      } else{
        extension[,, -(start:end)] <- obj
      }
      sampler$samples[[obj_name]] <- extension
    }
  }
  sampler$samples$stage <- c(old$stage, rep(stage, n_samples))
  sampler
}