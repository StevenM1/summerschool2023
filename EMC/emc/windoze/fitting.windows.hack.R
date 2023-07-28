#### Fitting automation
require(parallel)

# iter=c(300,NA,NA); verbose=FALSE;verbose_run_stage=FALSE;
#   max_adapt_trys=10;particles=NA;particle_factor=100; p_accept= NULL; n_cores=1;
#   epsilon = NULL; start_mu = NULL; start_var = NULL; mix=NULL;  
#   pdist_update_n=50;min_unique=200;epsilon_upper_bound=2

run_stages <- function(sampler,iter=c(300,NA,NA),
                       verbose=FALSE,verbose_run_stage=FALSE,
  max_adapt_trys=10,particles=NA,particle_factor=100, p_accept= NULL, n_cores=1,
  epsilon = NULL, start_mu = NULL, start_var = NULL, mix=NULL,  
  pdist_update_n=50,min_unique=200,epsilon_upper_bound=2) 
  # Initializes (if needed) then runs burn, adapt and sample if iter is not 
  # NA where iter[1] = burn, iter[2] = adapt, iter[3] = sample
  # Adapt stage is run repeatedly up to max_adapt_trys times. 
  # Number of particles is set at particle_factor*sqrt(number of parameters) if
  # particles is NA
{
  
  if (is.na(particles)) 
      particles <- round(particle_factor*sqrt(length(sampler$par_names)))
  if (!sampler$init) {
      sampler <- init(sampler, n_cores = n_cores, particles = particles, 
        epsilon = epsilon, start_mu = start_mu, start_var = start_var) 
  }
  if (all(is.na(iter))) return(sampler)
  if ( !is.na(iter[1]) ) {
    if (verbose) message("Running burn stage")
      sampler <- run_stage(sampler, stage = "burn",iter = iter[1], particles = particles, 
        n_cores = n_cores, pstar = p_accept, epsilon = epsilon, verbose = verbose_run_stage,
        min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound, mix=mix)
    if (all(is.na(iter[2:3]))) return(sampler)
  }
  if (!is.na(iter[2])) {
    if (verbose) message("Running adapt stage")
    trys <- 0
    idx <- sampler$samples$idx
    repeat {
       sampler <- run_stage(sampler, stage = "adapt", iter = iter[2], epsilon = epsilon, 
          particles = particles, n_cores = n_cores, verbose = verbose_run_stage, pstar = p_accept,
          min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound, mix=mix)
      trys <- trys + 1
      if (verbose) cat(paste0(trys," "))
      if (trys > max_adapt_trys) cat("\nAdaptation failed\n")
      if (trys > max_adapt_trys | sampler$samples$idx-idx < iter[2]) break
      idx <- sampler$samples$idx
    }
    if (trys > max_adapt_trys) 
      attr(sampler,"adapted") <- paste(" failed after",max_adapt_trys*iter[2],"iterations") else
      attr(sampler,"adapted") <- sum(sampler$samples$stage=="adapt")
    if (is.na(iter[3])) return(sampler)
  }
  if (verbose) message("Running sample stage")
  sampler <- run_stage(sampler, stage = "sample", iter = iter[3], epsilon = epsilon,
        pdist_update_n=pdist_update_n,particles = particles, n_cores = n_cores, 
        pstar = p_accept, verbose = verbose_run_stage, mix=mix,
        min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound)  
  sampler
}



# iter=c(NA,1000,500);verbose=FALSE;
#   max_adapt_trys=10;particles=NA;particle_factor=100; p_accept= 0.7;
#   cores_per_chain=1;cores_for_chains=NULL;
#   epsilon = NULL; start_mu = NULL; start_var = NULL;
#   pdist_update_n=50; epsilon_upper_bound=2
run_chains <- function(samplers,iter=c(300,NA,NA),
  verbose=TRUE,verbose_run_stage=FALSE,mix=NULL,
  max_adapt_trys=10,particles=NA,particle_factor=100, p_accept= 0.7, 
  cores_per_chain=1,cores_for_chains=NULL,min_unique=200,epsilon_upper_bound=2,
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
    samplers_new <- lapply(samplers,run_stages,iter=iter,
      verbose=verbose,verbose_run_stage=verbose_run_stage,
      max_adapt_trys=max_adapt_trys,particles=particles,particle_factor=particle_factor, 
      p_accept=p_accept, min_unique=min_unique, mix=mix,
      epsilon = epsilon, start_mu = start_mu, start_var = start_var,
      pdist_update_n=pdist_update_n,epsilon_upper_bound=epsilon_upper_bound,
      n_cores=cores_per_chain)
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



# iter=NA;max_trys=50;verbose=FALSE;burn=FALSE
# max_gd=1.1;thorough=TRUE;epsilon = NULL;particles=NA;particle_factor=100
# p_accept=NULL;cores_per_chain=1;cores_for_chains=NA

run_gd <- function(samplers,iter=NA,max_trys=100,verbose=FALSE,burn=FALSE,
                   max_gd=1.1,thorough=TRUE, epsilon = NULL,pdist_update_n=50, 
                   particles=NA,particle_factor=100, p_accept=NULL,
                   cores_per_chain=1,cores_for_chains=NA,mix=NULL,
                   min_unique=200,epsilon_upper_bound=2,
                   min_es=NULL,min_iter=NULL) 
  # Repeatedly runs burn or sample to get subject average multivariate 
  # gelman.diag of alpha samples < max_gd (if !through) or if all univariate
  # psrf for every subject and parameter and multivariate < max_gd and
  # if min_es specified at least that effective size in worst case and if 
  # min_iter specified at least that many iterations.
  # Trys adding iter (if NA 1/3 initial length) of the length and keeps all or 
  # extra minus first n_remove (initially iter, then grows by iter/2 on each
  # expansion based on mean gd of alphas (mpsrf alone or with psrf depending on thorough)
  # until success or max_trys. Verbose prints out progress after each try
  # Cores used = cores_per_chain*cores_for_chains, latter set to number of chains by default
{
  
  enough_samples <- function(samplers,min_es,min_iter,filter="burn") {
    if (is.null(min_iter)) ok <- TRUE else
      ok <- samplers[[1]]$samples$idx > min_iter
    if (!is.null(min_es)) {
      es_min <- min(es_pmwg(as_mcmc.list(samplers,selection="alpha",filter=filter)))
      ok <- ok & (es_min > min_es)
      attr(ok,"es") <- es_min
    }
    ok
  }
  
  if (!burn & !check_adapt(samplers,FALSE))
    stop("Can only use burn=FALSE if samplers adpated")
  if (burn) filter <- "burn" else filter <- "sample"
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  if (!is.list(samplers) || !all(unlist(lapply(samplers,class))== "pmwgs")) 
      stop("samplers must be a list of pmwgs objects")
  if (thorough) {
    max_name <- ", Max alpha mpsrf/psrf = "
    message("Exit on max alpha psrf/mpsrf < ",max_gd)
  } else {
    max_name <- ", Max alpha msrf = "
    message("Exit on max alpha mpsrf < ",max_gd)
  }
  n_chains <- length(samplers)
  if (is.na(cores_for_chains)) cores_for_chains <- n_chains
  if (burn) stage <- "burn" else stage <- "sample"
  n <- chain_n(samplers)
  if (burn) {
    if (sum(n[,-1])>0) 
      stop("Can only use burn=TRUE if samplers only have burn samples")
    idxs <- n[,1]
  } else {
    if (sum(n[,3])==0) 
      stop("Can only use burn=FALSE if there are sample stage samples")
    idxs <- n[,3]
  }
  if (!all(idxs[1]==idxs[-1]))
    stop("Must have same number of iterations for each chain")
  if (is.na(iter)) iter <- round(idxs[1]/3)
  n_remove <- iter
  # Iterate until criterion
  trys <- 0
  if (is.na(particles))
    particles <- round(particle_factor*sqrt(length(samplers[[1]]$par_names)))
  shorten <- FALSE
  repeat {
    run_try <- 0
    repeat {
      samplers_new <- lapply(samplers,run_stage,stage=stage,iter=iter, mix=mix,
        epsilon = epsilon, particles=particles,pstar=p_accept,verbose=FALSE,
        pdist_update_n=pdist_update_n,n_cores=cores_per_chain,
        min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound)
      if (class(samplers_new)=="try-error" || 
        any(lapply(samplers_new,class)=="try-error")) {
        save(samplers,stage,iter,particles,p_accept,pdist_update_n,
             min_unique,epsilon_upper_bound,epsilon_upper_bound,file="fail_run_stage.RData")
        run_try <- run_try + 1
        if (verbose) message("run_stage try error", run_try)
        if (run_try > 10)
          stop("Fail after 10 run_stage try errors, see saved samplers in fail_run_stage.RData")
      } else {
        gd <- gd_pmwg(as_mcmc.list(samplers_new,filter=filter),!thorough,FALSE,filter=filter)
        if (is.finite(gd[[1]])) break else
          if (verbose) message("gelman diag try error", run_try)
      }
    }
    samplers <- samplers_new
    trys <- trys+1
    if (shorten) {
      samplers_short <- lapply(samplers,remove_iterations,select=n_remove,filter=filter)
      gd_short <- gd_pmwg(as_mcmc.list(samplers_short,filter=filter),!thorough,FALSE,filter=filter)
      if (mean(gd_short) < mean(gd)) {
        gd <- gd_short
        samplers <- samplers_short
        n_remove <- iter
      } else {
        n_remove <- round(n_remove + iter/2)
      }
    }
    enough <- enough_samples(samplers,min_es,min_iter,filter=filter)
    if (is.null(attr(enough,"es"))) es_message <- NULL else
      es_message <- paste(", Effective samples =",round(attr(enough,"es")))
    ok_gd <- all(gd < max_gd)
    shorten <- !ok_gd
    if (trys > max_trys || (ok_gd & enough)) {
      if (verbose) {
        message("Final multivariate gelman.diag per participant")
        message("\nIterations = ",samplers[[1]]$samples$idx,", Mean mpsrf= ",
                round(mean(gd),3),max_name,round(max(gd),3))
      }
      attr(samplers,"data_list") <- data_list
      attr(samplers,"design_list") <- design_list
      attr(samplers,"model_list") <- model_list
      return(samplers)
    }
    if (verbose) {
        chain_n(samplers)[,filter][1]
        message(trys,": Iterations (",filter,") = ",chain_n(samplers)[,filter][1],", Mean mpsrf= ",
                round(mean(gd),3),max_name,round(max(gd),3),es_message)
    }
  }
}


# burn=TRUE;ndiscard=200;nstart=300;nadapt=1000;
#     discard_start=TRUE;start_particles=NA;
#     start_mu = NULL; start_var = NULL;
#     start_mix=NULL;single_start_mix=c(.5,.5);
#     verbose=TRUE;verbose_run_stage=FALSE;
#     max_gd_trys=100;max_gd=1.1;thorough=TRUE;min_es=NULL;min_iter=NULL;
#     epsilon = NULL; epsilon_upper_bound=2;
#     particles=NA;particle_factor=100; sample_particle_factor=100; p_accept=0.7;
#     pdist_update_n=50;min_unique=200; mix=NULL;
#     cores_per_chain=1;cores_for_chains=NULL
# 
# ndiscard=90;nstart=500;min_es=500;cores_per_chain=10

auto_burn <- function(samplers,
    burn=TRUE,ndiscard=200,nstart=300,nadapt=1000,
    discard_start=TRUE,start_particles=NA,
    start_mu = NULL, start_var = NULL,
    start_mix=NULL,single_start_mix=c(.5,.5),
    verbose=TRUE,verbose_run_stage=FALSE,
    max_gd_trys=100,max_gd=1.1,thorough=TRUE,min_es=NULL,min_iter=NULL,
    epsilon = NULL, epsilon_upper_bound=2, 
    particles=NA,particle_factor=100, sample_particle_factor=100, p_accept=0.7,
    pdist_update_n=50,min_unique=200, mix=NULL,
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
    if (ndiscard==0) discard_start <- FALSE
    if (!discard_start) discard_message <- NULL else
      discard_message <- paste("(first ",ndiscard," then removed)")
    message("Getting initial ",ndiscard + nstart," samples ",discard_message)
    run_try <- 0
    if (samplers[[1]]$source=="pmwg/variants/single.R") start_mix <- single_start_mix
    if (ndiscard>0) {
      repeat {
        samplers_new <- lapply(samplers,run_stages,iter=c(ndiscard,NA,NA),
          n_cores=cores_per_chain,p_accept = p_accept,
          epsilon_upper_bound=epsilon_upper_bound, mix=start_mix,
          particles=start_particles,particle_factor=particle_factor,epsilon=epsilon,
          start_mu = start_mu, start_var = start_var, 
          verbose=FALSE,verbose_run_stage=verbose_run_stage,
          pdist_update_n=pdist_update_n,min_unique=min_unique)
        if (class(samplers_new)=="try-error" || 
          any(lapply(samplers_new,class)=="try-error")) {
          save(samplers,ndiscard,particles,particle_factor,p_accept,pdist_update_n,
              min_unique,epsilon_upper_bound,file="fail_run_stage.RData")
          run_try <- run_try + 1
          if (verbose) message("discard run_stage try error", run_try)
          if (run_try > 10)
            stop("Fail after 10 discard run_stage try errors, see saved samplers in fail_run_stage.RData")
        } else break
      }
      samplers <- samplers_new
    }
    run_try <- 0
    repeat {
      samplers_new <- lapply(samplers,run_stages,iter=c(nstart,NA,NA),
        n_cores=cores_per_chain,p_accept = p_accept,
        epsilon_upper_bound=epsilon_upper_bound, mix=mix,
        particles=particles,particle_factor=particle_factor,epsilon=epsilon,
        verbose=FALSE,verbose_run_stage=verbose_run_stage,
        pdist_update_n=pdist_update_n,min_unique=min_unique)
      if (class(samplers_new)=="try-error" || 
        any(lapply(samplers_new,class)=="try-error")) {
        save(samplers,nstart,particles,particle_factor,p_accept,pdist_update_n,
            min_unique,epsilon_upper_bound,file="fail_run_stage.RData")
        run_try <- run_try + 1
        if (verbose) message("burn run_stage try error", run_try)
        if (run_try > 10)
          stop("Fail after 10 burn run_stage try errors, see saved samplers in fail_run_stage.RData")
      } else break
    }
    samplers <- samplers_new
    if (discard_start) samplers <- lapply(samplers,remove_iterations,select=nstart+1)
    message("Finished initial run")
    attr(samplers,"data_list") <- data_list
    attr(samplers,"design_list") <- design_list
    attr(samplers,"model_list") <- model_list
  } 
  if (!burn) { # adapt
    if (!check_adapt(samplers,verbose=FALSE)) {
      message("Running adapt stage")
      run_try <- 0
      repeat {
        samplers_new <- lapply(samplers,run_stages,iter=c(NA,nadapt,NA),
          n_cores=cores_per_chain,p_accept = p_accept,
          epsilon_upper_bound=epsilon_upper_bound, mix=mix,
          particles=particles,particle_factor=particle_factor,epsilon=epsilon,
          verbose=FALSE,verbose_run_stage=verbose_run_stage,
          pdist_update_n=pdist_update_n,min_unique=min_unique)
        if (class(samplers_new)=="try-error" || 
          any(lapply(samplers_new,class)=="try-error")) {
          save(samplers,nadapt,particles,particle_factor,p_accept,pdist_update_n,
              min_unique,epsilon_upper_bound,file="fail_run_stage.RData")
          run_try <- run_try + 1
          if (verbose) message("adapt run_stage try error", run_try)
          if (run_try > 10)
            stop("Fail after 10 adapt run_stage try errors, see saved samplers in fail_run_stage.RData")
        } else break
      }
      samplers <- samplers_new
      if (check_adapt(samplers,verbose=FALSE)) message("Adaptation sucessful")
      attr(samplers,"data_list") <- data_list
      attr(samplers,"design_list") <- design_list
      attr(samplers,"model_list") <- model_list
    }
    if (check_adapt(samplers,verbose=FALSE)) {
      nsample <- chain_n(samplers)[,"sample"]
      if (any(nsample<nstart)) {
        message("Running start sample stage")
        nstart <- pmin(nstart,nstart-min(nsample))
        particle_factor <- sample_particle_factor
        run_try <- 0
        repeat {
          samplers_new <- lapply(samplers,run_stages,iter=c(NA,NA,nstart),
            n_cores=cores_per_chain,p_accept = p_accept,
            epsilon_upper_bound=epsilon_upper_bound, mix=mix,
            particles=particles,particle_factor=particle_factor,epsilon=epsilon,
            verbose=FALSE,verbose_run_stage=verbose_run_stage,
            pdist_update_n=pdist_update_n,min_unique=min_unique)
          if (class(samplers_new)=="try-error" || 
            any(lapply(samplers_new,class)=="try-error")) {
            save(samplers,nadapt,particles,particle_factor,p_accept,pdist_update_n,
                min_unique,epsilon_upper_bound,file="fail_run_stage.RData")
            run_try <- run_try + 1
            if (verbose) message("start sample run_stage try error", run_try)
            if (run_try > 10)
              stop("Fail after 10 start sample run_stage try errors, see saved samplers in fail_run_stage.RData")
          } else break
        }
        samplers <- samplers_new
        attr(samplers,"data_list") <- data_list
        attr(samplers,"design_list") <- design_list
        attr(samplers,"model_list") <- model_list
      }
    } else {
      message("Adaptation failed")
      attr(samplers,"data_list") <- data_list
      attr(samplers,"design_list") <- design_list
      attr(samplers,"model_list") <- model_list
      return(samplers)
    } 
  }
  if (max_gd_trys==0) return(samplers)
  message("Beginning iterations to achieve Rhat < ",max_gd)
  if (!burn) particle_factor <- sample_particle_factor
  run_gd(samplers,burn=burn,max_trys=max_gd_trys,verbose=verbose,
          max_gd=max_gd,thorough=thorough, p_accept = p_accept,
          pdist_update_n=pdist_update_n,min_unique=min_unique,
          epsilon=epsilon, particles=particles,particle_factor=particle_factor,
          min_es=min_es,min_iter=min_iter, mix=mix,
          cores_per_chain=cores_per_chain,cores_for_chains=cores_for_chains)
 }
