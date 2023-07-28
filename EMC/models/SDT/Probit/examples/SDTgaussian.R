rm(list=ls())
source("emc/emc.R")
source("models/SDT/Probit/SDTgaussian.R")
# This matchfun useful for SDT models, assumes binary noise/signal factor and 
# an even number of confidence ratings
matchfun <- function(d) as.numeric(d$S) == (as.numeric(d$lR)>length(levels(d$lR))/2)


# Two choice probit

# Note that for models of type="SDT" (such as probit) lR is automatically given
# a special "contr.increasing" contrast, and any terms involving the threshold
# for the last level of lR are set to a large positive constant (really this
# should be Inf for the upper integration limit of the top response category,
# but Inf causes problems so the Inf is introduced in the probit model object).

# We adopt the conventional factor S = noise vs. signal nomenclature and 
# scaling where noise mean = 0 and log(sd) = 0 (sd=1)

#### Equal variance binary choice ----

designPROBIT <- make_design(Flist=list(mean ~ S, sd ~ S,threshold ~ lR),
  Ffactors=list(subjects=1,S=1:2),Rlevels=1:2,matchfun=function(d)d$S==d$lR,
  constants=c(mean=0,sd=0),model=probit)

p_vector <- sampled_p_vector(designPROBIT)
p_vector[1] <- 2         # signal mean
p_vector[2] <- log(1)    # signal sd 
p_vector[3] <- 1         # threshold, natural scale

# For type="SDT" models this function strips out the top level of lR as there
# are only r-1 thresholds to estimates (r = number of response levels,
# internally there are still r threshold paramters but the last one is a
# constant
mapped_par(p_vector,designPROBIT) 

# Simulate a large amount of data, bit slow here as designed for n-choice case
dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=100000)
# Check proportion correct
crct=dataPROBIT$R==dataPROBIT$S
tapply(dataPROBIT$R==dataPROBIT$S,dataPROBIT$S,mean)
# Check probability of a hit
tapply(as.numeric(dataPROBIT$R)-1,dataPROBIT$S,mean)
# Should be about the same as this
c(pnorm(1,mean=2,sd=1),1-pnorm(1,mean=2,sd=1))


# profiles
dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=10000)
dadmPROBIT <- design_model(data=dataPROBIT,design=designPROBIT)
par(mfrow=c(1,3))
for (i in names(p_vector))
  print(profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.5,p_max=p_vector[i]+.5,dadm=dadmPROBIT))

dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=1000)
samplers <- make_samplers(dataPROBIT,designPROBIT,type="single")
save(samplers,file="probitIndividualBinary.RData")
# runSingleProbitBinary.R to get 1000 samples
print(load("probitIndividualBinary.RData"))

plot_chains(samples,subfilter=100,filter="burn") # Thoroughly converged by 400
gd_pmwg(samples,subfilter=100,filter="burn")  # 1.01
plot_acfs(samples,subfilter=100,filter="burn")
iat_pmwg(samples,subfilter=100,filter="burn") # ~ 50% yield 
#   mean_S2 sd_S2 threshold threshold_lR2 threshold_lR3 threshold_lR4 threshold_lR5
# 1    2.09  2.22      2.04          1.74          2.07          2.31           2.5

# Excellent recovery, prior completely dominated
tabs <- plot_density(samples,filter="burn",subfilter=100,layout=c(1,3),pars=p_vector)

ppBinary <- post_predict(samples,filter="burn",subfilter=100)

pc <- function(d) 100*mean(d$S==d$R)
# Drilling down on accuracy we see the biggest misfit is in speed for right
# responses by > 5%.
tab <- plot_fit(data=dataPROBIT,pp=ppBinary,stat=pc,
                stat_name="Accuracy (%)",layout=c(1,2))


##### ROC example, 3 level confidence ----



designPROBIT <- make_design(Flist=list(mean ~ S, sd ~ S,threshold ~ lR),
  Ffactors=list(subjects=1,S=1:2),Rlevels=1:6, matchfun=matchfun,
  constants=c(mean=0,sd=0),model=probit)

# Signal > noise variance typical of recognition memory
p_vector <- sampled_p_vector(designPROBIT)
p_vector[1] <- 1            # signal mean
p_vector[2] <- log(1.25)    # signal sd = 1.25 (treatment coding)
p_vector[3] <- -.5          # first threshold untransformed
# other thresholds exponentiated so > 0 then added to previous
# this is done in p_vector transform function
p_vector[4:7] <- log(rep(.5,4))  

# Thresholds evenly spaced with 0.5 gap
mapped_par(p_vector,designPROBIT) 

# Make some data and plot ROCs
dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=10000)

# 0.8 = 1/sd(signal) slope evident 
par(mfrow=c(1,2))
plot_roc(dataPROBIT)
plot_roc(dataPROBIT,zROC=TRUE,qfun=qnorm)

# profiles and sampling
dadmPROBIT <- design_model(data=dataPROBIT,design=designPROBIT)

par(mfrow=c(2,4))
for (i in names(p_vector))
  print(profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.5,p_max=p_vector[i]+.5,dadm=dadmPROBIT))

samplers <- make_samplers(dataPROBIT,designPROBIT,type="single")
save(samplers,file="probitIndividual.RData")
# runSingleProbit.R to get 1000 samples
print(load("probitIndividual.RData"))

plot_chains(samples,subfilter=400) # Thoroughly converged by 400
gd_pmwg(samples,subfilter=400)  # 1.01
plot_acfs(samples,subfilter=400,layout=c(2,4))
iat_pmwg(samples,subfilter=400) # ~ 50% yield 
#   mean_S2 sd_S2 threshold threshold_lR2 threshold_lR3 threshold_lR4 threshold_lR5
# 1    2.09  2.22      2.04          1.74          2.07          2.31           2.5

# Excellent recovery, prior completely dominated
tabs <- plot_density(samples,subfilter=400,layout=c(2,4),pars=p_vector)
#       mean_S2 sd_S2 threshold threshold_lR2 threshold_lR3 threshold_lR4 threshold_lR5
# true    1.000 0.223    -0.500        -0.693        -0.693        -0.693        -0.693
# 2.5%    1.002 0.216    -0.502        -0.714        -0.697        -0.712        -0.701
# 50%     1.015 0.228    -0.494        -0.703        -0.685        -0.700        -0.687
# 97.5%   1.030 0.238    -0.485        -0.690        -0.673        -0.686        -0.672

### Some further examples of more flexible designs ----

##### 3 level confidence, factor A increases d' in a mirror pattern ----

designPROBIT <- make_design(Flist=list(mean ~ A*S, sd ~ S,threshold ~ lR),
  Ffactors=list(subjects=1,S=c("new","old"),A=c("faces","words")),Rlevels=1:6, matchfun=matchfun,
  constants=c(mean=0,sd=0),model=probit)

p_vector <- sampled_p_vector(designPROBIT)
p_vector[1:3] <- c(-.5,1,1)  # mirror effect A=1 0,1, A=2 -.5,1.5
p_vector[4] <- log(1.25) # signal sd = 1.25 (treatment coding)
p_vector[5] <- -1        # first threshold untransformed
p_vector[6:9] <- log(rep(.5,4))  

# shift in criterion by 0.5 evident
mapped_par(p_vector,designPROBIT) 


dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=10000)
par(mfrow=c(2,2))
plot_roc(dataPROBIT[dataPROBIT$A=="faces",],main="A=faces")
plot_roc(dataPROBIT[dataPROBIT$A=="faces",],zROC=TRUE,qfun=qnorm,main="A=faces",lim=c(-1.5,2))
plot_roc(dataPROBIT[dataPROBIT$A=="words",],main="A=words")
plot_roc(dataPROBIT[dataPROBIT$A=="words",],zROC=TRUE,qfun=qnorm,main="A=words",lim=c(-1.5,2))

# profiles
dadmPROBIT <- design_model(data=dataPROBIT,design=designPROBIT)
par(mfrow=c(2,5))
for (i in names(p_vector))
  print(profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.25,p_max=p_vector[i]+.25,dadm=dadmPROBIT))


samplers <- make_samplers(dataPROBIT,designPROBIT,type="single")
save(samplers,file="probitIndividualM.RData")
# runSingleProbitM.R to get 1000 samples
print(load("models/SDT/Probit/examples/samples/probitIndividualM.RData"))

plotChains(samples,subfilter=400) # Thoroughly converged by 400
gd_pmwg(samples,subfilter=400)  # 1.01
plot_acfs(samples,subfilter=400,layout=c(2,4))
iat_pmwg(samples,subfilter=400) # ~ 50% yield 
# Excellent recovery, prior completely dominated
tabs <- plot_density(samples,subfilter=400,layout=c(2,5),pars=p_vector)


##### 3 level confidence, factor A shifts threshold up ----

designPROBIT <- make_design(Flist=list(mean ~ S, sd ~ S,threshold ~ A+lR),
  Ffactors=list(subjects=1,S=1:2,A=1:2),Rlevels=1:6, matchfun=matchfun,
  constants=c(mean=0,sd=0),model=probit)

p_vector <- sampled_p_vector(designPROBIT)
p_vector[1] <- 1         # signal mean
p_vector[2] <- log(1.25) # signal sd = 1.25 (treatment coding)
p_vector[3] <- -1        # first threshold untransformed
p_vector[4] <- .5        # A2 shift up
p_vector[5:8] <- log(rep(.5,4))  

# shift in criterion by 0.5 evident
mapped_par(p_vector,designPROBIT) 


dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=10000)
par(mfrow=c(2,2))
plot_roc(dataPROBIT[dataPROBIT$A==1,],main="A=1")
plot_roc(dataPROBIT[dataPROBIT$A==1,],zROC=TRUE,qfun=qnorm,main="A=1",lim=c(-1.5,2))
plot_roc(dataPROBIT[dataPROBIT$A==2,],main="A=2")
plot_roc(dataPROBIT[dataPROBIT$A==2,],zROC=TRUE,qfun=qnorm,main="A=2",lim=c(-1.5,2))

# profiles
dadmPROBIT <- design_model(data=dataPROBIT,design=designPROBIT)
par(mfrow=c(2,4))
for (i in names(p_vector))
  print(profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.5,p_max=p_vector[i]+.5,dadm=dadmPROBIT))


samplers <- make_samplers(dataPROBIT,designPROBIT,type="single")
save(samplers,file="probitIndividualA.RData")
# runSingleProbitA.R to get 1000 samples
print(load("models/SDT/Probit/examples/samples/probitIndividualA.RData"))

plot_chains(samples,subfilter=400) # Thoroughly converged by 400
gd_pmwg(samples,subfilter=400)  # 1.01
plot_acfs(samples,subfilter=400,layout=c(2,4))
iat_pmwg(samples,subfilter=400) # ~ 50% yield 
#   mean_S2 sd_S2 threshold threshold_lR2 threshold_lR3 threshold_lR4 threshold_lR5
# 1    2.09  2.22      2.04          1.74          2.07          2.31           2.5

# Excellent recovery, prior completely dominated
tabs <- plot_density(samples,subfilter=400,layout=c(2,4),pars=p_vector)


# NB: The real data example shows how to set arbitrarily different (but still 
#     increasing) thresholds for different levels of a factor or factors by
#     by "nesting" lR within those factors.

##### Combine last two ----

# Because of shift in threshold with A and mean, must also fix mean of new A=2
# for identifiable

designPROBIT <- make_design(Flist=list(mean ~ A*S, sd ~ S,threshold ~ A+lR),
  Ffactors=list(subjects=1,S=1:2,A=1:2),Rlevels=1:6, matchfun=matchfun,
  constants=c(mean=0,mean_A2=0,sd=0),model=probit)

p_vector <- sampled_p_vector(designPROBIT)
p_vector[1:2] <- c(1,1)  # mirror effect A=1 0,1, A=2 -.5,1.5
p_vector[3] <- log(1.25) # signal sd = 1.25 (treatment coding)
p_vector[4] <- -1        # first threshold untransformed
p_vector[5] <- .5        # A2 shift up
p_vector[6:9] <- log(rep(.5,4))  

# shift in criterion by 0.5 evident
mapped_par(p_vector,designPROBIT) 


dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=10000)
par(mfrow=c(2,2))
plot_roc(dataPROBIT[dataPROBIT$A==1,],main="A=1")
plot_roc(dataPROBIT[dataPROBIT$A==1,],zROC=TRUE,qfun=qnorm,main="A=1",lim=c(-1.5,2))
plot_roc(dataPROBIT[dataPROBIT$A==2,],main="A=2")
plot_roc(dataPROBIT[dataPROBIT$A==2,],zROC=TRUE,qfun=qnorm,main="A=2",lim=c(-1.5,2))

# profiles
dadmPROBIT <- design_model(data=dataPROBIT,design=designPROBIT)
par(mfrow=c(2,5))
for (i in names(p_vector))
  print(profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.25,p_max=p_vector[i]+.25,dadm=dadmPROBIT))


samplers <- make_samplers(dataPROBIT,designPROBIT,type="single")
# save(samplers,file="probitIndividualMA.RData")
# runSingleProbitMA.R to get 1000 samples
print(load("models/SDT/Probit/examples/samples/probitIndividualMA.RData"))

plotChains(samples,subfilter=400) # Thoroughly converged by 400
gd_pmwg(samples,subfilter=400)  # 1.01
plot_acfs(samples,subfilter=400,layout=c(2,5))
iat_pmwg(samples,subfilter=400) # ~ 30-40% yield 
# Excellent recovery, prior completely dominated
tabs <- plotDensity(samples,subfilter=400,layout=c(2,5),pars=p_vector)












