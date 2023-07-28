source("samplers/pmwg/sampling.R")

sample_store_nested <- function(data, par_names, iters = 1, stage = "init", integrate = T,...) {
  args <- list(...)
  pop_args <- args$pop_lvl
  group_args <- args$group_lvl
  add_info_to_groups <- function(g, group_samples, group_data, pars){
    group_data <- group_data[[g]]
    samples <- group_samples[[g]]
    subjects <- unique(group_data$subject)
    out <- list(
      par_names = pars,
      subjects = subjects,
      n_pars = length(pars),
      n_subjects = length(subjects),
      samples = samples
    )
    return(out)
  }
  # Pop level
  if(is.null(pop_args$variant)){
    source("pmwg/variants/standard.R")
  } else{
    source(pop_args$variant)
  }
  pop_args$data <- data
  pop_args$par_names <- par_names
  pop_samples <- do.call(variant_funs$sample_store, pop_args)
  # Group level
  if(is.null(group_args$variant)){
    source("pmwg/variants/standard.R")
  } else{
    source(group_args$variant)
  }
  n_groups <- length(unique(data$group))
  # Set up group args
  group_args$data <- data
  group_args$par_names <- par_names
  group_args$integrate <- F
  # Create a single samples object then replicate
  group_sample <- do.call(variant_funs$sample_store, group_args)
  group_samples <- rep(list(group_sample), n_groups)
  # Add the relevant info back
  group_data <- split(data, data$group)
  group_samples <- list(group_samples = lapply(1:n_groups,add_info_to_groups, group_samples, group_data, par_names))
  # Make sure we get variant_funs back
  source("pmwg/variants/nested.R")
  return(c(pop_samples, group_samples))
}

add_info_nested <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  args <- list(...)
  pop_args <- args$pop_lvl
  group_args <- args$group_lvl
  sampler$par_idx <- args$par_idx
  # Pop level
  if(is.null(pop_args$variant)){
    source("pmwg/variants/standard.R")
  } else{
    source(pop_args$variant)
  }
  sampler$pop_funs <- variant_funs
  pop_args$sampler <- sampler
  sampler <- do.call(sampler$pop_funs$add_info, pop_args)
  # Group level
  if(is.null(group_args$variant)){
    source("pmwg/variants/standard.R")
  } else{
    source(group_args$variant)
  }
  
  groups <- unique(sampler$data$group)
  sampler$group_idx <- aggregate(group ~ subject, sampler$data, mean)[,2]
  sampler$n_groups <- length(groups)
  sampler$group_funs <- variant_funs
  for(i in 1:sampler$n_groups){
    group_args$sampler <- sampler$samples$group_samples[[i]]
    sampler$samples$group_samples[[i]] <- do.call(sampler$group_funs$add_info, group_args)
  }
  # Make sure we get the right functions back again
  source("pmwg/variants/nested.R")
  return(sampler)
}

get_startpoints_nested <- function(pmwgs, start_mu, start_var){
  # No custom startpoints between groups for you! (who cares and too much work)
  start_points_pop <- pmwgs$pop_funs$get_startpoints(pmwgs, start_mu, start_var)
  start_points_group <- replicate(pmwgs$n_groups, list(pmwgs$group_funs$get_startpoints(pmwgs, start_mu, start_var)))
  # Make sure they're in the right format for the get_group_level_nested func
  return(list(pop_lvl = start_points_pop, pop_fill = pmwgs$pop_funs$fill_samples, 
              group_lvl = start_points_group, group_fill = pmwgs$group_funs$fill_samples,
              group_idx = pmwgs$group_idx))
}

get_group_level_nested <- function(parameters, s){
  mu <- parameters$group_lvl[[parameters$group_idx[s]]]$tmu
  var <- parameters$group_lvl[[parameters$group_idx[s]]]$tvar
  return(list(mu = mu, var = var))
}

fill_samples_nested <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples <- group_level$pop_fill(samples, group_level$pop_lvl, proposals, epsilon, j = j, n_pars)
  for(g in 1:length(samples$group_samples)){ # 1 line Loops look ugly :(
    samples$group_samples[[g]]$samples <- group_level$group_fill(samples$group_samples[[g]]$samples, group_level$group_lvl[[g]], proposals = NULL, epsilon, j = j, n_pars)
  }
  return(samples)
}

gibbs_step_nested <- function(sampler, alpha){
  # Pfwoah this is clean ;) 
  # Pop level, first create new 'alpha' (actually group level means)
  alpha_groups <- do.call(cbind, lapply(sampler$samples$group_samples, 
                                        FUN = function(x) return(x$samples$theta_mu[,sampler$samples$idx])))
  # The number of observations is now the number of groups
  sampler$n_subjects <- sampler$n_groups 
  pop_lvl <- sampler$pop_funs$gibbs_step(sampler, alpha_groups)
  # Group level
  group_idx <- sampler$group_idx
  group_lvl <- vector(mode = "list", length = sampler$n_groups)
  for(i in 1:sampler$n_groups){
    # First idx the random effects
    alpha_current <- alpha[,sampler$group_idx == i]
    group_samples <- sampler$samples$group_samples[[i]]
    group_samples$samples$idx <- sampler$samples$idx
    group_lvl[[i]] <- sampler$group_funs$gibbs_step(group_samples, alpha_current)
  }
  return(list(pop_lvl = pop_lvl, pop_fill = sampler$pop_funs$fill_samples, 
              group_lvl = group_lvl, group_fill = sampler$group_funs$fill_samples,
              group_idx = sampler$group_idx, alpha = alpha))
}

get_conditionals_nested <- function(s, samples, n_pars){
  samples$theta_mu <- samples$theta_mu[[samples$group[s]]]
  samples$theta_var <- samples$theta_var[[samples$group[s]]]
  out <- samples$cond_fun(s, samples, n_pars)
  return(list(eff_mu = out$eff_mu, eff_var = out$eff_var))
}

filtered_samples_nested <- function(sampler, filter){
  theta_mu <- lapply(sampler$samples$group_samples, 
                     FUN = function(x) return(x$samples$theta_mu[,filter]))
  theta_var <- lapply(sampler$samples$group_samples, 
                     FUN = function(x) return(x$samples$theta_var[,,filter]))
  out <- list(
    theta_mu = theta_mu,
    theta_var = theta_var,
    alpha = sampler$samples$alpha[, , filter],
    iteration = length(filter),
    group = sampler$group_idx,
    cond_fun = sampler$group_funs$get_conditionals
  )
}

variant_funs <- list(
  sample_store = sample_store_nested,
  add_info = add_info_nested,
  get_startpoints = get_startpoints_nested,
  get_group_level = get_group_level_nested,
  fill_samples = fill_samples_nested,
  gibbs_step = gibbs_step_nested,
  filtered_samples = filtered_samples_nested,
  get_conditionals = get_conditionals_nested
)