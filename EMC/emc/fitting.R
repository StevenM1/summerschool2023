#### Fitting automation
require(parallel)

require(abind)
run_stages <- function(sampler,iter=c(300,0,0),
                       verbose=FALSE,verbose_run_stage=FALSE,
                       max_adapt_trys=2,particles=NA,particle_factor=100, p_accept= NULL, n_cores=1,
                       epsilon = NULL, start_mu = NULL, start_var = NULL, mix=NULL,  
                       pdist_update_n=50,min_unique=200,epsilon_upper_bound=15,
                       eff_var = NULL, eff_mu = NULL) 
  # Initializes (if needed) then runs burn, adapt and sample if iter is not 
  # NA where iter[1] = burn, iter[2] = adapt, iter[3] = sample
  # Adapt stage is run repeatedly up to max_adapt_trys times. 
  # Number of particles is set at particle_factor*sqrt(number of parameters) if
  # particles is NA
{
  
  if (is.na(particles)) 
    particles <- round(particle_factor*sqrt(length(sampler$par_names)))
  if (!sampler$init) {
    sampler <- init(sampler, n_cores = n_cores, 
                    epsilon = epsilon, start_mu = start_mu, start_var = start_var) 
  }
  if (all(iter==0)) return(sampler)
  if ( iter[1] != 0 ) {
    if (verbose) message("Running burn stage")
    sampler <- run_stage(sampler, stage = "burn",iter = iter[1], particles = particles, 
                         n_cores = n_cores, pstar = p_accept, epsilon = epsilon, verbose = verbose_run_stage,
                         min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound, mix=mix)
    if (all(iter[2:3]==0)) return(sampler)
  }
  if (iter[2] != 0) {
    if (verbose) message("Running adapt stage")
    sampler <- run_stage(sampler, stage = "adapt",iter = iter[2], particles = particles, 
                         n_cores = n_cores, pstar = p_accept, epsilon = epsilon, verbose = verbose_run_stage,
                         min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound, mix=mix)
    if (iter[3]==0) return(sampler)
  }
  if (verbose) message("Running sample stage")
  sampler <- run_stage(sampler, stage = "sample", iter = iter[3], epsilon = epsilon,
                       pdist_update_n=pdist_update_n,particles = particles, n_cores = n_cores, 
                       pstar = p_accept, verbose = verbose_run_stage, mix=mix, eff_mu = eff_mu,
                       eff_var = eff_var,
                       min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound)  
  sampler
}



# iter=c(500,0,0);
#   verbose=TRUE;verbose_run_stage=FALSE;mix=NULL;
#   max_adapt_trys=10;particles=NA;particle_factor=100; p_accept= 0.7;
#   cores_per_chain=1;cores_for_chains=NULL;min_unique=200;epsilon_upper_bound=2;
#   epsilon = NULL; start_mu = NULL; start_var = NULL;pdist_update_n=50
# verbose=TRUE; iter=c(3,0,0)
run_chains <- function(samplers,iter=c(300,0,0),
                       verbose=TRUE,verbose_run_stage=FALSE,mix=NULL,
                       max_adapt_trys=10,particles=NA,particle_factor=100, p_accept= 0.7, 
                       cores_per_chain=1,cores_for_chains=NULL,min_unique=200,epsilon_upper_bound=15,
                       epsilon = NULL, start_mu = NULL, start_var = NULL,pdist_update_n=50) 
  # applies run stages over chains preserving list attributes
{
  source(samplers[[1]]$source)
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  run_try <- 0
  repeat {
    samplers_new <- mclapply(samplers,run_stages,iter=iter,
                             verbose=verbose,verbose_run_stage=verbose_run_stage,
                             max_adapt_trys=max_adapt_trys,particles=particles,particle_factor=particle_factor, 
                             p_accept=p_accept, min_unique=min_unique, mix=mix,
                             epsilon = epsilon, start_mu = start_mu, start_var = start_var,
                             pdist_update_n=pdist_update_n,epsilon_upper_bound=epsilon_upper_bound,
                             n_cores=cores_per_chain, mc.cores = cores_for_chains)
    if (class(samplers_new)=="try-error" || 
        any(lapply(samplers_new,class)=="try-error")) {
      save(samplers,iter,particles,particle_factor,p_accept,pdist_update_n,
           epsilon_upper_bound,min_unique,file="fail_run_stage.RData")
      run_try <- run_try + 1
      if (verbose) message("run_stage try error", run_try)
      if (run_try > 10)
        stop("Fail after 10 run_stage try errors, see saved samplers in fail_run_stage.RData")
    } else break
  }
  samplers <- samplers_new
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  samplers
}


run_gd <- function(samplers,iter=NA,max_trys=100,verbose=FALSE,burn=TRUE,
                   max_gd=1.1,thorough=TRUE, mapped=FALSE, shorten = TRUE,
                   epsilon = NULL, verbose_run_stage = FALSE,
                   particles=NA,particle_factor=50, p_accept=NULL,
                   cores_per_chain=1,cores_for_chains=NA,mix=NULL,
                   min_es=NULL,min_iter=NULL,max_iter=NULL) 
  # Repeatedly runs burn or sample to get subject average multivariate 
  # gelman.diag of alpha samples < max_gd (if !through) or if all univariate
  # psrf for every subject and parameter and multivariate < max_gd and
  # if min_es specified at least that effective size in worst case and if 
  # min_iter specified at least that many iterations. If max_iter specified
  # pulls out after max_iter.
  # Trys adding iter (if NA 1/3 initial length) of the length and keeps all or 
  # extra minus first n_remove (initially iter, then grows by iter/2 on each
  # expansion based on mean gd of alphas (mpsrf alone or with psrf depending on thorough)
  # until success or max_trys. Verbose prints out progress after each try
  # Cores used = cores_per_chain*cores_for_chains, latter set to number of chains by default
{
  
  enough_samples <- function(samplers,min_es,min_iter,max_iter,filter="burn") {
    if (!is.null(max_iter) && (sum(samplers[[1]]$samples$stage==filter) >= max_iter)) return(TRUE)
    if (is.null(min_iter)) ok <- TRUE else
      ok <- (sum(samplers[[1]]$samples$stage==filter)) > min_iter
    if (!is.null(min_es)) {
      es_min <- min(es_pmwg(as_mcmc.list(samplers,selection="alpha",filter=filter)))
      ok <- ok & (es_min > min_es)
      attr(ok,"es") <- es_min
    }
    ok
  }
  
  if (thorough) {
    max_name <- ", Max alpha mpsrf/psrf = "
    message("Exit on max alpha psrf/mpsrf < ",max_gd)
  } else {
    max_name <- ", Max alpha msrf = "
    message("Exit on max alpha mpsrf < ",max_gd)
  }
  n_chains <- length(samplers)
  n <- chain_n(samplers)
  if (is.na(cores_for_chains)) cores_for_chains <- n_chains
  if (is.na(iter)) iter <- round(n[,1][1]/3)
  n_remove <- iter
  # Iterate until criterion
  trys <- 0
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  gd <- gd_pmwg(as_mcmc.list(samplers,filter="burn"),!thorough,FALSE,
                filter="burn",mapped=mapped)
  if (all(is.finite(gd))) gd_diff <- apply(gd, 1, max) - 1.5*max_gd else gd_diff <- NA
  repeat {
    run_try <- 0
    repeat {
      new_particle_n <- adaptive_particles(gd_diff, max_gd, particle_factor, particles)
      particle_factor <- new_particle_n$particle_factor
      particles <- new_particle_n$particles
      if(!any(is.na(gd_diff)) & any(gd_diff > .5)) samplers <- check_stuck(samplers) # Maybe a dumb heuristic
      samplers_new <- mclapply(samplers,run_stages,iter=c(iter,0,0),
                               n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                               particles=particles,particle_factor=particle_factor,epsilon=epsilon,
                               verbose=FALSE,verbose_run_stage=verbose_run_stage,
                               mc.cores=cores_for_chains)
      if (!class(samplers_new)=="list" || !all(unlist(lapply(samplers_new,class))=="pmwgs")) 
        gd <- matrix(Inf) else
        gd <- gd_pmwg(as_mcmc.list(samplers_new,filter="burn"),!thorough,FALSE,
                    filter="burn",mapped=mapped)
      if (all(is.finite(gd))) break else {
        run_try <- run_try + 1
        message("gelman diag try error ", run_try)
        if (run_try > max_trys) stop("Gave up after ",max_trys," gelman diag errors.")
      }
    }
    samplers <- samplers_new
    trys <- trys+1
    if (shorten) {
      samplers_short <- lapply(samplers,remove_iterations,select=n_remove,filter="burn")
      if (!class(samplers_short)=="list" || !all(unlist(lapply(samplers_short,class))=="pmwgs")) 
        gd_short <- matrix(Inf) else
        gd_short <- gd_pmwg(as_mcmc.list(samplers_short,filter="burn"),!thorough,FALSE,
                          filter="burn",mapped=mapped)
      if (mean(gd_short) < mean(gd)) {
        gd <- gd_short
        samplers <- samplers_short
        n_remove <- iter
      } else {
        n_remove <- round(n_remove + iter/2)
      }
    }
    enough <- enough_samples(samplers,min_es,min_iter,max_iter,filter=filter)
    if (is.null(attr(enough,"es"))) es_message <- NULL else
      es_message <- paste(", Effective samples =",round(attr(enough,"es")))
    if (all(is.finite(gd))) {
      gd_diff <- (gd[,ncol(gd)] - 2*max_gd)
      ok_gd <- all(gd < max_gd)
      shorten <- !ok_gd
    } else ok_gd <- FALSE
    if (trys == max_trys || (ok_gd & enough)) {
      if (verbose) {
        message("\nFinal multivariate gelman.diag per participant")
        message("Iterations = ",samplers[[1]]$samples$idx,", Mean mpsrf= ",
                round(mean(gd),3),max_name,round(max(gd),3))
      }
      attr(samplers,"data_list") <- data_list
      attr(samplers,"design_list") <- design_list
      attr(samplers,"model_list") <- model_list
      if (ok_gd) attr(samplers,"burnt") <- max_gd else attr(samplers,"burnt") <- NA 
      return(samplers)
    }
    if (verbose) {
      chain_n(samplers)[,"burn"][1]
      message(trys,": Iterations burn = ",chain_n(samplers)[,"burn"][1],", Mean mpsrf= ",
              round(mean(gd),3),max_name,round(max(gd),3),es_message)
    }
  }
}

# ndiscard=80;nstart=120;
# particles=NA; particle_factor = 50; start_mu = NULL; start_var = NULL;
# mix = NULL; verbose=TRUE;verbose_run_stage=FALSE;
# max_gd_trys=100;max_gd=1.1;
# thorough=TRUE;mapped=FALSE; step_size = NA;
# min_es=NULL;min_iter=NULL;max_iter=NULL;
# epsilon = NULL; epsilon_upper_bound=15; p_accept=0.7;
# cores_per_chain=10;cores_for_chains=NULL
auto_burn <- function(samplers,ndiscard=100,nstart=300,
                      particles=NA, particle_factor = 50, start_mu = NULL, start_var = NULL,
                      mix = NULL, verbose=TRUE,verbose_run_stage=FALSE,
                      max_gd_trys=100,max_gd=1.2,
                      thorough=TRUE,mapped=FALSE, step_size = NA,
                      min_es=NULL,min_iter=NULL,max_iter=NULL,
                      epsilon = NULL, epsilon_upper_bound=15, p_accept=0.7,
                      cores_per_chain=1,cores_for_chains=NULL)
  # Takes a pmwgs chains list, initializes it (see run_stages), if !burn adapts
  # and runs burn or sample until gd criterion satisfied (see run_gd for details)
{
  source(samplers[[1]]$source)
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  if (is.null(samplers[[1]]$samples$idx) || samplers[[1]]$samples$idx==1) {
    discard_message <- paste("(first ",ndiscard," then removed)")
    message("Getting initial ",ndiscard + nstart," samples ",discard_message)
    run_try <- 0
    samplers_new <- mclapply(samplers,run_stages,iter=c(nstart + ndiscard,0,0),
                             n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                             particles=particles,particle_factor=particle_factor,epsilon=epsilon,
                             verbose=FALSE,verbose_run_stage=verbose_run_stage,
                             mc.cores=cores_for_chains)
    samplers <- samplers_new
    if(ndiscard!= 0) samplers <- lapply(samplers,remove_iterations,select=ndiscard+1)
    message("Finished initial run")
    attr(samplers,"data_list") <- data_list
    attr(samplers,"design_list") <- design_list
    attr(samplers,"model_list") <- model_list
  } 
  if (max_gd_trys==0) return(samplers)
  message("Beginning iterations to achieve Rhat < ",max_gd)
  run_gd(samplers, max_trys=max_gd_trys,verbose=verbose,
         max_gd=max_gd,thorough=thorough, mapped=mapped, p_accept = p_accept,
         epsilon=epsilon, particles=particles,particle_factor=particle_factor, min_es=min_es,
         min_iter=min_iter, 
         max_iter=max_iter, mix=mix, iter = step_size, verbose_run_stage = verbose_run_stage,
         cores_per_chain=cores_per_chain,cores_for_chains=cores_for_chains)
}


# min_particles = 50; max_particles = 500;min_factor = 25; max_factor = 100; percent_up = 10; percent_down = 5
adaptive_particles <- function(gd_diff, max_gd, particle_factor = NA, particles = NA, 
                               min_particles = 50, max_particles = 500, 
                               min_factor = 25, max_factor = 100,
                               percent_up = 10, percent_down = 5){
  # This function adaptively tunes the particles per participant,
  # so that we can lower the number of particles is we're closer to gd_criterion,
  # percent_up and down are relative to the max. Percent up is scaled by sqrt(gd_diff)
  if (any(is.na(gd_diff))) {
    if (is.na(particles)) 
      return(list(particles=100,particle_factor = particle_factor)) else
      return(list(particles=particles,particle_factor = particle_factor))
  }
  if(is.na(particles)){
    gd_diff[gd_diff > 0] <- sqrt(gd_diff[gd_diff > 0])*(percent_up/100)*max_factor
    gd_diff[gd_diff < 0] <- -(percent_down/100)*max_factor
    particle_factor <- pmin(pmax(min_factor, particle_factor + gd_diff), max_factor)
  } else{
    gd_diff[gd_diff > 0] <- sqrt(gd_diff[gd_diff > 0])*(percent_up/100)*max_particles
    gd_diff[gd_diff < 0] <- -(percent_down/100)*min_particles
    particles <- pmin(pmax(min_particles, particles + gd_diff), max_particles)
  }
  return(list(particles = round(particles), particle_factor = particle_factor))
}

auto_adapt <- function(samplers,max_trys=25,epsilon = NULL, 
                       particles=NA,particle_factor=40, p_accept=.7,
                       cores_per_chain=1,cores_for_chains=NULL,mix=NULL,
                       n_cores_conditional = 1, min_es=NULL,min_unique = 200, 
                       step_size = 25, thin = NULL,
                       verbose=TRUE,verbose_run_stage = FALSE){
  if(verbose) message("Running adapt stage")
  source(samplers[[1]]$source)
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  burnt <- attr(samplers,"burnt")
  trys <- 0
  samplers_new <- mclapply(samplers,run_stages,iter=c(0,min_unique/(length(samplers)*p_accept) - step_size,0),
                           n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                           particles=particles,particle_factor=particle_factor,epsilon=epsilon,
                           verbose=FALSE,verbose_run_stage=verbose_run_stage,
                           mc.cores=cores_for_chains)
  repeat {
    trys <- trys + 1
    samplers_new <- mclapply(samplers_new,run_stages,iter=c(0,step_size,0),
                             n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                             particles=particles,particle_factor=particle_factor,epsilon=epsilon,
                             verbose=FALSE,verbose_run_stage=verbose_run_stage,
                             mc.cores=cores_for_chains)
    test_samples <- lapply(samplers_new, extract_samples, stage = "adapt", thin = thin, samplers_new[[1]]$samples$iteration, thin_eff_only = F)
    keys <- unique(unlist(lapply(test_samples, names)))
    test_samples <- setNames(do.call(mapply, c(abind, lapply(test_samples, '[', keys))), keys)
    test_samples$iteration <- sum(test_samples$iteration)
    # Only need information like n_pars & n_subjects, so only need to pass the first chain
    adapted <- test_adapted(samplers_new[[1]], test_samples, min_unique, n_cores_conditional, verbose_run_stage)
    if(trys > max_trys | adapted) break
  }
  samplers <- samplers_new
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  attr(samplers,"burnt") <- burnt
  attr(samplers, "adapted") <- adapted
  return(samplers)
}


auto_sample <- function(samplers,iter=NA,verbose=TRUE,
                        epsilon = NULL, particles=NA,particle_factor=25, p_accept=.7,
                        cores_per_chain=1,cores_for_chains=NULL,mix=NULL,
                        n_cores_conditional = 1, step_size = 50, thin = NULL,
                        verbose_run_stage = FALSE)
  # Automatically run the sample stage  
{
  if(!attr(samplers, "adapted")){
    warning("samplers should be adapted before you can run sample stage")
    return(samplers)
  }
  if(verbose) message("Running sample stage")
  source(samplers[[1]]$source)
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  burnt <- attr(samplers,"burnt")
  adapted <- attr(samplers,"adapted")
  samplers_new <- samplers
  n_steps <- ceiling(iter/step_size)
  for(step in 1:n_steps){
    if(step == n_steps){
      step_size <- ifelse(iter %% step_size == 0, step_size, iter %% step_size)
    } 
    test_samples <- lapply(samplers_new, extract_samples, 
                           stage = c("adapt", "sample"), thin = thin, 50*trys, thin_eff_only = FALSE)
    keys <- unique(unlist(lapply(test_samples, names)))
    test_samples <- setNames(do.call(mapply, c(abind, lapply(test_samples, '[', keys))), keys)
    test_samples$iteration <- sum(test_samples$iteration)
    conditionals=mclapply(X = 1:samplers_new[[1]]$n_subjects,
                          FUN = variant_funs$get_conditionals,samples = test_samples, 
                          samplers_new[[1]]$n_pars, mc.cores = n_cores_conditional)
    conditionals <- array(unlist(conditionals), dim = c(samplers_new[[1]]$n_pars, 
                                                        samplers_new[[1]]$n_pars + 1, samplers_new[[1]]$n_subjects))
    eff_mu <- conditionals[,1,] #First column is the means
    eff_var <- conditionals[,2:(samplers_new[[1]]$n_pars+1),] #Other columns are the variances
    samplers_new <- mclapply(samplers_new,run_stages,iter=c(0,0,step_size),
                             n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                             particles=particles,particle_factor=particle_factor,epsilon=epsilon,
                             verbose=FALSE,verbose_run_stage=verbose_run_stage, eff_mu = eff_mu,
                             eff_var = eff_var,
                             mc.cores=cores_for_chains)
  }
  samplers <- samplers_new
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  attr(samplers,"burnt") <- burnt
  attr(samplers, "adapted") <- adapted
  return(samplers)
}

run_IS2 <- function(samples, filter = "sample", subfilter = 0, IS_samples = 1000, 
                    stepsize_particles = 500, max_particles = 5000, n_cores = 1, df = 5){
  variant <- basename(samples[[1]]$source)
  source(paste0("samplers/IS2/variants/", variant))
  samples_merged <- merge_samples(samples)
  IS2(samples_merged, filter, subfilter = subfilter, IS_samples, stepsize_particles, max_particles, n_cores, df)
}


test_adapted <- function(sampler, test_samples, min_unique, n_cores_conditional = 1, 
                         verbose = FALSE)
{
  # Only need to check uniqueness for one parameter
  first_par <- test_samples$alpha[1, , ]
  # Split the matrix into a list of vectors by subject
  # Needed for the case where every sample is unique for all subjects
  first_par_list <- split(first_par, seq(NROW(first_par)))
  # Get unique pars (new accepted particles) and check length for
  # all subjects is greater than unq_vals
  n_unique_sub <- lapply(lapply(first_par_list, unique), length)
  n_pars <- sampler$n_pars
  if (all(n_unique_sub > min_unique)) {
    if(verbose){
      message("Enough unique values detected: ", min_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
      mclapply(X = 1:sampler$n_subjects,FUN = variant_funs$get_conditionals,samples = test_samples, 
               n_pars, mc.cores = n_cores_conditional)
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
      if(verbose){
        message("A problem was encountered creating proposal distribution")
        message("Increasing required unique values and continuing adaptation")
      }
      return(FALSE)
    }
    else {
      if(verbose) message("Successfully adapted after ", test_samples$iteration, "iterations - stopping adaptation")
      return(TRUE)
    }
  } else{
    return(FALSE) # Not enough unique particles found
  }
}

check_stuck <- function(samples,filter=c("burn","sample")[1], # dont use adapt
                        start=1,last=TRUE,n=90, # n, from start or last n
                        flat=0.1,dfac=3) 
  # flat: criterion on % change in median first to last 1/3 relative to last 1/3 Let me know what you think.
  
  # dfac: bad if  median(best)-median(chain) > dfac*IQR(best)
  # returns numeric(0) if none stuck, otherwise idices of stuck chains
  
{
  samplell <- function(sampler,filter,subfilter)
    sampler$samples$subj_ll[,sampler$samples$stage==filter][,subfilter]
  
  ns <- unlist(lapply(samples,function(x){sum(x$samples$stage==filter)}))
  if (!all(ns[1]==ns[-1])) stop("Filtered chains not of same length")
  if (n>ns[1]) stop("n is larger than available samples")
  if (last) start <- ns[1]-n + 1
  iter <- 1:n
  lls <- do.call(abind, list(lapply(samples,samplell,filter=filter,subfilter=iter+start-1), along = 3))
  first <- 1:(round(n/3))
  last <- round(2*n/3):n
  first <- apply(lls[,first,],c(1,3),median)
  last <- apply(lls[,last,],c(1,3),median)
  isFlat <- 100*abs((last-first)/last) < flat
  # if(any(isFlat)){
  diff <- apply(last, 1, FUN = function(x) return(max(x) - min(x)))
  IQRs <- apply(lls, c(1,3), FUN = IQR)
  best <- max.col(last)
  worst <- max.col(-last)
  n_subs <- samples[[1]]$n_subjects
  bad_subs <- which(diff > dfac*IQRs[matrix(c(1:n_subs, best), nrow = n_subs)]) # Yay R tricks
  if(any(bad_subs)) samples <- fix_stuck(samples, bad_subs, best, worst)
  # }
  
  return(samples)
}

std_error_IS2 <- function(IS_samples, n_bootstrap = 50000){
  log_marglik_boot= array(dim = n_bootstrap)
  for (i in 1:n_bootstrap){
    log_weight_boot = sample(IS_samples, length(IS_samples), replace = TRUE) #resample with replacement from the lw
    log_marglik_boot[i] <- median(log_weight_boot)
  }
  return(sd(log_marglik_boot))
} 


fix_stuck <- function(samples, bad_subs, best, worst){
  # This function takes the bad subjects and replaces the last entry of the 
  # worst chain with the best chain.
  idx <- samples[[1]]$samples$idx
  for(i in 1:length(bad_subs)){
    samples[[worst[i]]]$samples$alpha[,i,idx] <- samples[[best[i]]]$samples$alpha[,i,idx]
  }
  return(samples)
}


run_emc <- function(file_name,nsample=1000, ...) 
  # Combined auto fitting functions, getting and saving samples to disk.
  #    NB: samples must be first object in file_name file on disk
  # OR if file_name is not a character vector file_name is treated as a samplers 
  #    object and results of fitting returned by the function.
{
  if (is.character(file_name)) {
    sname <- load(file_name)
    if (is.null(attr(get(sname[1]),"burnt")) || is.na(attr(get(sname[1]),"burnt"))) {
      assign(sname[1],auto_burn(get(sname[1]), ...))
      save(list=sname,file=file_name)
    }
    if (is.null(attr(get(sname[1]),"adapted")) && !is.na(attr(get(sname[1]),"burnt"))) {
      assign(sname[1],auto_adapt(get(sname[1]), ...))
      save(list=sname,file=file_name)
    }
    if (!is.null(attr(get(sname[1]),"adapted")) && attr(get(sname[1]),"adapted")) {
      assign(sname[1],auto_sample(get(sname[1]),iter=nsample, ...))
      save(list=sname,file=file_name)
    }
  } else { # file_name is actually a samplers object
    if (is.null(attr(file_name,"burnt")) || is.na(attr(file_name,"burnt")))
      file_name <- auto_burn(file_name, ...)
    if (is.null(attr(file_name,"adapted")) && !is.na(attr(file_name,"burnt")))
      file_name <- auto_adapt(file_name, ...)
    if (!is.null(attr(file_name,"adapted")) && attr(file_name,"adapted"))
      file_name <- auto_sample(file_name,iter=nsample, ...)
    return(file_name)
  }
}

# all_par <- c(particles=NA, particle_factor = 50, mix = NULL,
# epsilon = NULL, epsilon_upper_bound=15, p_accept=0.7,cores_per_chain=1,
# cores_for_chains=NULL,verbose=TRUE,verbose_run_stage = FALSE)
# 
# # step_size = NA (in burn, but = 25 in adapt and 50 sample)
# 
# burn_adapt <- c(min_es=NULL)
# 
# adapt_sample_par <- c(n_cores_conditional = 1,thin = NULL)
# 
# auto_burn_par <- c(ndiscard=100,nstart=300,start_mu = NULL, start_var = NULL,
# max_gd_trys=100,max_gd=1.2,thorough=TRUE,mapped=FALSE, min_iter=NULL,max_iter=NULL)
# 
# auto_adapt_par <- c(max_trys=25, min_unique = 200)
# 
# auto_sample_par <- c(iter=NA)
  