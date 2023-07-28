# LNR mu and sigma parameterization

source("models/RACE/LNR/lnr.R")

lnrMS <- list(
  type="RACE",
  p_types=c("m","s","t0"),
  # Transform to natural scale
  Ntransform=function(x) {

    get_p_types <- function(nams) 
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
    
    if (!is.matrix(x)) {
      nams <- get_p_types(names(x))
      x[nams != "m"] <- exp(x[nams != "m"]) 
    } else {
      nams <- get_p_types(dimnames(x)[[2]])
      x[,nams != "m"] <- exp(x[,nams != "m"])
    }
    x
  },
  # mapped parameter transform
  Mtransform = function(pars) 
    # transform parameters except v back to real line 
    # pars is a matrix output by map_p_vector  
  {
    lnrMS$Ntransform(pars)
  },
  # p_vector transform scaling parameter by s=1 assumed in lnr.R
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rLNR(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dLNR(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pLNR(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race_missing(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

