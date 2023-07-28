# Lognromal race

#### distribution functions


dLNR <- function(rt,pars)
  # density for single accumulator
{
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
  out[ok] <- dlnorm(rt[ok],meanlog=pars[ok,"m"],sdlog=pars[ok,"s"])
  out
}



pLNR <- function(rt,pars)
  # cumulative density for single accumulator
{
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
  out[ok] <- plnorm(rt[ok],meanlog=pars[ok,"m"],sdlog=pars[ok,"s"])
  out
  
}

#### random


rLNR <- function(lR,pars,p_types=c("m","s","t0")) 
  # lR is an empty latent response factor lR with one level for each accumulator. 
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in 
  # contiguous rows.
  #
  # test
  # pars=cbind(m=c(-1,0),s=c(1,1),t0=c(.3,.3)); lR=factor(c(1,2))
{
  if (!all(p_types %in% dimnames(pars)[[2]])) 
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt <- matrix(rlnorm(dim(pars)[1],meanlog=pars[,"m"],sdlog=pars[,"s"]),
               nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(pars[,"t0"],nrow=length(levels(lR)))[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}

