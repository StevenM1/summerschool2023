# Probit discrete choice

#### distribution functions


pPROBIT <- function(lt,ut,pars)
  # probability between lt and ut
{
  pnorm(ut,mean=pars[,"mean"],sd=pars[,"sd"]) - pnorm(lt,mean=pars[,"mean"],sd=pars[,"sd"])
}

#### random


rPROBIT <- function(lR,pars,p_types=c("mean","sd","threshold"),lt=-Inf) 
  # lR is an empty latent response factor lR with one level for response. 
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in 
  # contiguous rows.
  
{
  if (!all(p_types %in% dimnames(pars)[[2]])) 
    stop("pars must have columns ",paste(p_types,collapse = " "))
  nr <- length(levels(lR)) # Number of responses
  n <- dim(pars)[1]/nr     # Number of simulated trials
  first <- seq(1,dim(pars)[1]-nr+1,length.out=n)  # pick out mean and sd
  threshold <- matrix(pars[,"threshold"],nrow=nr) # format thresholds
  threshold[dim(threshold)[1],] <- Inf
  pmat <- rbind(rnorm(n,pars[first,"mean"],pars[first,"sd"]), # sample normal
                rep(lt,dim(threshold)[2]),threshold)          # lt ... ut
  pmat[dim(pmat)[1],] <- Inf
  R <- factor(apply(pmat,2,function(x){.bincode(x[1],x[-1])}),
              levels=1:length(levels(lR)),labels=levels(lR))
  cbind.data.frame(R=R,rt=NA)
}

