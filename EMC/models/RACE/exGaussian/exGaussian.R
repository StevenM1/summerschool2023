# Normal, mu and log(sigma) parameterization

source("models/RACE/exGaussian/exG.R")

exGaussian <- list(
  type="RACE",
  p_types=c("mu","sigma","tau"),
  # Transform to natural scale
  Ntransform=function(x) {
    exp(x)
  },
  # mapped parameter transform
  Mtransform = function(pars) 
    # transform parameters except v back to real line 
    # pars is a matrix output by map_p_vector  
  {
    exGaussian$Ntransform(pars)
  },
  # p_vector transform
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rexGaussian(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dexGaussian(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pexGaussian(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

