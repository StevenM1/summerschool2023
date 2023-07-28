# Unbounded DDM parameterization on log and probit scales mapped to rtdists
# with s=1 scaling


source("models/DDM/DDM/ddm.R")


ddmTZDt0natural <- list(
  type="DDM",
  p_types=c("v","a","sv","t0","st0","s","Z","SZ","DP"),
  # The "TZD" parameterization defined relative to the "rtdists" package is:
  # natural scale
  #   v = rtdists rate v (positive favors upper)
  # log scale 
  #   t0 > 0: lower bound of non-decision time 
  #   st0 > 0: rtdists width of non-decision time distribution 
  #   a > 0: rtdists upper threshold, a
  #   sv > 0: rtdists v standard deviation sv
  #   s > 0: rtdists moment-to-moment standard deviation, s
  # probit scale
  #   0 < Z < 1: rtdists start point z = Z*a 
  #   0 < SZ < 1: rtdists start-point variability, sz = 2*SZ*min(c(a*Z,a*(1-Z)) 
  #   0 < DP < 1: rtdists d = t0(upper)-t0(lower) = (2*DP-1)*t0  # 
  # 
  # Transform to natural scale
  Ntransform=function(x) {
    
    get_p_types <- function(nams) 
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
    
    if (!is.matrix(x)) {
      nams <- get_p_types(names(x))
      islog <- nams %in%  c("a","sv","st0","s")
      isprobit <- nams %in%  c("Z","SZ","DP")
      x[islog] <- exp(x[islog]) 
      x[isprobit] <- pnorm(x[isprobit])
    } else { 
      nams <- get_p_types(dimnames(x)[[2]])
      islog <- nams %in%  c("a","sv","st0","s")
      isprobit <- nams %in%  c("Z","SZ","DP")
      x[,islog] <- exp(x[,islog])
      x[,isprobit] <- pnorm(x[,isprobit])  
    }
    x
  },
  # mapped parameter transform
  Mtransform = function(pars) 
    # pars is a matrix output by map_p_vector  
  {
    pars <- ddmTZDt0natural$Ntransform(pars)
    pars <- cbind(pars,z=pars[,"a"]*pars[,"Z"],
      sz = 2*pars[,"SZ"]*pars[,"a"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min))
    pars <- cbind(pars, d = pars[,"t0"]*(2*pars[,"DP"]-1))
    attr(pars,"ok") <- 
      !( pars[,"t0"] < 0 | abs(pars[,"v"])>5 | pars[,"a"]>2 | pars[,"sv"]>2 | 
         pars[,"SZ"]>.75 | pars[,"st0"]>.2) 
    if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .01
    if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .01  
    pars
  },
  # p_vector transform, sets s as a scaling parameter
  transform = function(p) p,
  # Random function
  rfun=function(lR,pars) rDDM(lR,pars,precision=3),
  # Density function (PDF)
  dfun=function(rt,R,pars) dDDM(rt,R,pars,precision=3),
  # Probability function (CDF)
  pfun=function(rt,R,pars) pDDM(rt,R,pars,precision=3),
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_ddm(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)
