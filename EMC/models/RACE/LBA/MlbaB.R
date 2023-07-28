# lba_B parameterization with sv=1 scaling

source("models/RACE/LBA/lba.R")

MlbaB <- list(
  type="RACE",
  p_types=c("v","sv","B","A","t0"),
  # Transform to natural scale
  Ntransform=function(x) {

    get_p_types <- function(nams) 
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))

    if (!is.matrix(x)) {
      nams <- get_p_types(names(x))
      x[nams != "v"] <- exp(x[nams != "v"]) 
    } else {
      nams <- get_p_types(dimnames(x)[[2]])
      x[,nams != "v"] <- exp(x[,nams != "v"])
    }
    x
  },
  # mapped parameter transform
  Mtransform = function(pars) 
    # transform parameters except v back to real line and add b
    # pars is a matrix output by map_p_vector  
  {
    pars <- lbaB$Ntransform(pars)
    pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
    pars
  },
  # p_vector transform, sets sv as a scaling parameter
  transform = function(p) p,
  # Random function for racing accumulator
  rfun=function(lR,pars) rLBA(lR,pars,posdrift=TRUE),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race_missing(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

