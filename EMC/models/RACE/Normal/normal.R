# Normal, mu and log(sigma) parameterization

source("models/RACE/Normal/norm.R")

normal <- list(
  type="RACE",
  p_types=c("mean","sd"),
  # Transform to natural scale
  Ntransform=function(x) {

    get_p_types <- function(nams) 
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
    
    if (!is.matrix(x)) {
      nams <- get_p_types(names(x))
      x[nams == "sd"] <- exp(x[nams == "sd"]) 
    } else {
      nams <- get_p_types(dimnames(x)[[2]])
      x[,nams == "sd"] <- exp(x[,nams == "sd"])
    }
    x
  },
  # mapped parameter transform
  Mtransform = function(pars) 
    # transform parameters except v back to real line 
    # pars is a matrix output by map_p_vector  
  {
    normal$Ntransform(pars)
  },
  # p_vector transform
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rNORMAL(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dNORMAL(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pNORMAL(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

