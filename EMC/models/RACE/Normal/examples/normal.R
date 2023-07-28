rm(list=ls())
library(mvtnorm)
source("emc/emc.R")
source("models/RACE/Normal/normal.R")


# Two choice normal
designNORMAL <- make_design(Flist=list(mean ~ lM, sd ~ 1),
  Ffactors=list(subjects=1,S=1:2),Rlevels=1:2,matchfun=function(d)d$S==d$lR,
  Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert,lR=contr.helmert),
  model=normal)

p_vector <- sampled_p_vector(designNORMAL)
p_vector[1:2] <- c(0,1)
p_vector[3] <- c(log(1)) # sd

dataNORMAL <- make_data(p_vector,design=designNORMAL,trials=1000)
plot_defective_density(dataNORMAL,layout=c(1,2))

dadmNORMAL <- design_model(dataNORMAL,designNORMAL)
par(mfrow=c(1,3))
profile_pmwg(pname="mean",p=p_vector,p_min=-.5,p_max=.5,dadm=dadmNORMAL)
profile_pmwg(pname="mean_lM1",p=p_vector,p_min=0,p_max=2,dadm=dadmNORMAL)
profile_pmwg(pname="sd",p=p_vector,p_min=log(.5),p_max=log(1.5),dadm=dadmNORMAL)


# This example shows EMC works in a model without choice (Rlevels=1, no matchfun)
# Make a reasonably complex design. 
designNORMAL <- make_design(Flist=list(mean ~ A*B, sd ~ B),
  Ffactors=list(subjects=50,A=1:3,B=1:2),Rlevels=1,model=normal)

p_vector <- sampled_p_vector(designNORMAL)
p_vector[1:6] <- c(1,1,2,2,0,0) # Increasing in A and B, no interaction
p_vector[7:8] <- c(log(1,2)) # more variable for B2 than B1 


dataNORMAL <- make_data(p_vector,design=designNORMAL,trials=100)
plot_defective_density(dataNORMAL,layout=c(2,3),rt="topleft")

dadmNORMAL <- design_model(dataNORMAL,designNORMAL)
par(mfrow=c(2,4))
for (i in names(p_vector))
  profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.5,p_max=p_vector[i]+.5,dadm=dadmNORMAL)

samplers <- make_samplers(dataNORMAL,designNORMAL,type="standard",rt_resolution=.001)

