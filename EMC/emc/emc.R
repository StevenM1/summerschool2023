require(coda)        # EMC
require(rtdists)     # EMC
require(magic)       # PMwG
require(MASS)        # PMwG
require(abind)       # PMwG
require(MCMCpack)    # PMwG IS2
require(corpcor)     # PMwG IS2
require(condMVNorm)  # PMwG IS2
require(parallel)    # PMwG IS2
require(mvtnorm)     # PMwG IS2


source("emc/objects.R")    # pmwg object manipulation + mcmc creation functions
source("emc/statistics.R") # chain statistics
source("emc/plotting.R")   # chain and data plotting
source("emc/data.R")       # data generation 
source("emc/map.R")        # parameter mapping
source("emc/likelihood.R") # likelihoods
source("emc/fitting.R")    # fitting automation
source("emc/design.R")     # make design matrix and view parameters & mapping
source("emc/priors.R")     # sampling from priors
source("emc/pmwg.R")       # interface with pmwg
source("emc/joint_functions.R")       # handles joint models
