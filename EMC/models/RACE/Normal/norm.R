# Normal race

#### distribution functions


dNORMAL <- function(rt,pars)
  # density for single accumulator
{
  dnorm(rt,mean=pars[,"mean"],sd=pars[,"sd"])
}



pNORMAL <- function(rt,pars)
  # cumulative density for single accumulator
{
  pnorm(rt,mean=pars[,"mean"],sd=pars[,"sd"])
}

#### random


rNORMAL <- function(lR,pars,p_types=c("mean","sd")) 
  # lR is an empty latent response factor lR with one level for each accumulator. 
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in 
  # contiguous rows.
  #
  # test
  # pars=cbind(mean=c(-1,0),sd=c(1,1)); lR=factor(c(1))
{
  if (!all(p_types %in% dimnames(pars)[[2]])) 
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt <- matrix(rnorm(dim(pars)[1],mean=pars[,"mean"],sd=pars[,"sd"]),
               nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  rt <- dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}

