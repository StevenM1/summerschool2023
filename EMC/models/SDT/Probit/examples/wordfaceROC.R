rm(list=ls())
source("emc/emc.R")
source("models/SDT/Probit/SDTgaussian.R")
# This matchfun useful for SDT models, assumes binary noise/signal factor and 
# an even number of confidence ratings
matchfun <- function(d) as.numeric(d$S) == (as.numeric(d$lR)>length(levels(d$lR))/2)

print(load("Data/wordfaceROC.RData"))
# remove RT for ROC fitting
wordfaceROC$rt <- NA

#### Fit a single subject ----
wordfaceROC1 <- wordfaceROC[wordfaceROC$subjects=="104",]
wordfaceROC1$subjects <- factor(as.character(wordfaceROC1$subjects))

# To allow different thresholds for words and faces we must nest the lR factor
# within the FW factor (FW/lR). This ensures that thresholds are always 
# increasing within each level of FW. FW*lR will not enforce the increasing 
# constraint (but FW + lR will). If you wanted to allow different but increasing
# thresholds across combinations of factors, say A and B, use (A*B)/lR etc.
designFW1 <- make_design(Flist=list(mean ~ FW*S, sd ~ FW*S,threshold ~ FW/lR),
  Ffactors=list(subjects=104,S=c("new","old"),FW=c("faces","words")),Rlevels=1:6, 
  matchfun=matchfun,
  constants=c(mean=0,mean_FWwords=0,sd=0,sd_FWwords=0),model=probit)

p_vector <- sampled_p_vector(designFW1)
# constant mean new face = 0 
p_vector[1:2] <- c(1,1) # old faces, old words 
# constant new face sd = 1
p_vector[3:4] <- log(c(1.25,1)) # old sd = 1.25, no item effect on sd 
p_vector[5:6] <- c(-0.5,-0.5)  # first threshold for faces, shift down by 0.5 for words
p_vector[7:14] <- log(rep(c(.5,.75),4)) # .5 threshold spacing for faces, 0.75 spacing for words

# Check mapping
mapped_par(p_vector,designFW1) 

# Simulate data
simFW1 <- make_data(p_vector,design=designFW1,trials=10000)
par(mfrow=c(2,2))
plot_roc(simFW1[simFW1$FW=="faces",],main="Faces")
plot_roc(simFW1[simFW1$FW=="faces",],zROC=TRUE,qfun=qnorm,main="Faces",lim=c(-2.5,2.5))
plot_roc(simFW1[simFW1$FW=="words",],main="Words")
plot_roc(simFW1[simFW1$FW=="words",],zROC=TRUE,qfun=qnorm,main="Words",lim=c(-2.5,2.5))

# profiles
dadmFWsim <- design_model(data=simFW1,design=designFW1)
par(mfrow=c(2,8))
for (i in names(p_vector))
  print(profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.25,p_max=p_vector[i]+.25,dadm=dadmFWsim))

samplers <- make_samplers(simFW1,designFW1,type="single")
save(samplers,file="probitFWsim.RData")
# runSingleProbitFWsim.R to get 1000 samples

samplers <- make_samplers(wordfaceROC1,designFW1,type="single")
# save(samplers,file="probitFW1.RData")
# runSingleProbitFW1.R to get 1000 samples

# Look at simulation, clearly not close at 1000, added another 2000
print(load("models/SDT/Probit/examples/samples/probitFWsim.RData"))
plot_chains(samples,subfilter=200,layout=c(4,4))  # Fully converged
gd_pmwg(samples,subfilter=200)
plotACFs(samples,subfilter=200,layout=c(4,4)) # Quite autocorrelated
iat_pmwg(samples,subfilter=200) 
# Excellent recovery
tabs <- plot_density(samples,subfilter=200,layout=c(4,4),pars=p_vector)

# Look at data, quick convergence 
print(load("models/SDT/Probit/examples/samples/probitFW1.RData")) 
plot_chains(samples,subfilter=100,layout=c(4,4)) # very fast convergence
gd_pmwg(samples,subfilter=100) # 1.01
plotACFs(samples,subfilter=100,layout=c(4,4)) # Quite autocorrelated
iat_pmwg(samples,subfilter=100) 
tabs <- plot_density(samples,subfilter=100,layout=c(4,4)) # prior domination ok

#### Fit full data set ----

designFW <- make_design(Flist=list(mean ~ FW*S, sd ~ FW*S,threshold ~ FW/lR),
  Ffactors=list(subjects=levels(wordfaceROC$subjects),S=levels(wordfaceROC$S),
                FW=levels(wordfaceROC$FW)),Rlevels=1:6, 
  matchfun=matchfun,
  constants=c(mean=0,mean_FWwords=0,sd=0,sd_FWwords=0),model=probit)

samplers <- make_samplers(wordfaceROC,designFW,type="standard")
# Speedup of 60 reflects the fact that likelihood for SDT models acts on
# counts in each category, so only one likelihood computation is requried for
# each category
# save(samplers,file="probitFW.RData")
#### Explore the fitting script ----  

#  Fitting is performed in the script runProbitFW.R (e.g., on a linux based system
# the command line is R CMD BATCH runProbitFW.R & to run it in background)

# We use three procedures that run the "burn", "adapt" and "sample" until they
# produce samples with the required characteristics. By default each chain gets 
# its own cores (this can be set with the cores_for_chains argument) and the 
# cores_per_chain argument above gives each chain 2 cores, so 6 are used in 
# total. Here are the functions called by runProbitFW.R
# 
# The first "burn" stage by default runs 500 iterations, discards the first 200 
# and repeatedly trys adding new iterations (and possibly removing initial)
# iterations until R hat is less than a criterion (1.1 by default) for all
# random effect parameters in all chains. The aim of this stage is to find
# the posterior mode and get chains suitable for the next "adapt" stage.
#
# samples <- auto_burn(samplers,cores_per_chain=4)
#
# The "adapt" stage develops and approximation to the posterior that will make
# sampling more efficient (less autocorrelated) in the final "sample" stage.
#
# samples <- auto_adapt(samples,cores_per_chain=4)
#
# Here in the final sample stage we ask for 5000 iterations per chain.
#
# samples <- auto_sample(samples,iter=5000,cores_per_chain=4)
# 
# Once sampling is completed the script also gets posterior predictive samples
# to enable model fit checks. By default this is based on randomly selecting 
# iterations from the final (sample) stage, and provides posterior predictives 
# for the random effects. Here we use one core pre participant.
# ppWordFace <- post_predict(samples,n_cores=12)

# Lets load in the results and look at them.


print(load("models/SDT/Probit/examples/samples/probitFW.RData")) 
# We can check the state of samplers
chain_n(samples)
# and check if adaptation has worked 
check_adapt(samples)

# In the burn phase all chains appear to have reached the posterior model 
plot_chains(samples,filter="burn",layout=c(3,5),selection="mu") 
plot_chains(samples,filter="burn",layout=c(3,5),selection="variance") 
plot_chains(samples,filter="burn",layout=c(4,4),selection="correlation") 
plot_chains(samples,filter="burn",layout=c(3,5),selection="alpha") 

# Now lets look at the final samples (the default setting of filter)
# We use some thinning (showing only every 10th iteration) to speed up plotting.
# All appear converged.
plot_chains(samples,layout=c(3,5),selection="mu",thin=10) 
plot_chains(samples,layout=c(3,5),selection="variance",thin=10) 
plot_chains(samples,layout=c(4,4),selection="correlation",thin=10) 
plot_chains(samples,layout=c(3,5),selection="alpha",thin=10) 

# pick a selection and look at the convergence diagnostics
selection="mu"; selection="variance"; selection="alpha"; layout=c(3,5)
selection="correlation"; layout=c(4,7)

# In all cases mixed
gd_pmwg(samples,selection=selection)
# Some large inefficiency with sd parameters for some participants
iat_pmwg(samples,selection=selection,summary_alpha=max) 
# Use min to look at worst case
round(es_pmwg(samples,selection=selection,summary_alpha=min))

# Check if priors dominated by posterior (i.e., data updates them)
tabs <- plot_density(samples,layout=layout,selection=selection)
# Good updating for population level, and tabs is a matrix
round(tabs,3)

# As expected for alpha updating is less and influence of population model 
# acting as a prior is more evident. In this case tabs is a list (e.g., here 
# for the first participant)
round(tabs[[1]],3)

#### Fit
# For type=SDT plot_fit requires a factor (by default "S", argument signalFactor) 
# whose first level is noise and second level is signal in order to construct an 
# ROC. Where this 2 level structure does not apply (e.g., different types 
# of signal) subset the factor. 

# Average fit
par(mfrow=c(2,2))
plot_fit(wordfaceROC,ppWordFace,factors=c("FW","S"))
plot_fit(wordfaceROC,ppWordFace,factors=c("FW","S"),zROC=TRUE,qfun=qnorm,lim=c(-1.5,1.25))

# Individual fits, here just ROC, noiser and sometimes points missing (as 
# participants did not use a confidence level) but overall failry good.
par(mfrow=c(2,4))
plot_fit(wordfaceROC,ppWordFace,zROC=TRUE,qfun=qnorm)


### Look at parameter estimates: mu

# We can look at mu parameters in two ways, 1) in terms of the parameters that were
# actually sampled (and hence always transformed to have no bounds) with names
# organized by type (for the probit model the types are "mean", "sd" and 
# "threshold") that are shown by this function:
sp_names <- p_names(samples); sp_names

# 2) in terms of the parameters that are mapped to the model parameterization,
# which may be bounded, and which correspond to the combinations of factors that
# specified in the formula for the parameter type they correspond to:
mp_names <- p_names(samples,mapped=TRUE); mp_names

# Note that in this form some of the parameters can be constants. The derivation
# of their names can be seen by looking at the design, being the unique cells
# for the factors in the formula for each type
mp_design <- p_names(samples,mapped=TRUE,design=TRUE); mp_design

# Lets look at sampled parameters first

# At hyper level the posterior (black lines) shows moderately strong domination 
# of the prior (red lines).
tab_mu <- plot_density(samples,selection="mu",layout=c(2,7))

# Start by looking at mean and sd parameters
round(tab_mu[,sp_names$mean],2)
round(tab_mu[,sp_names$sd],2)

# To understand these, recall that mean for new faces and words (mean & 
# mean_FWwords) are set to zero and corresponding sd (sd sd_FWwords) set to 1
# (so as sd is sampled on log scale sampled values are set to zero).

# Because we did not specify contrasts R's default treatment code was used.
# We can get a list of the contrasts for each parameter type and look at the
# mean and sd design matrices (which are the same as they share formulas).
maps <- get_map(samples)
maps$mean 
# Reading across row 3 (down column Hence, mean_Sold is the study effect for 
# faces (0.69) and mean_FWwords:Sold is the extra study effect for words 
# relative to faces (1.05).
maps$sd
# The same is true for sd (0.18 and 0.22 respectively), but recall the effect is 
# on the log scale.

tab_mu_mapped <- plot_density(samples,selection="mu",layout=c(2,7),mapped=TRUE)

# Need to remove cells set to constant in mp_names to look at estimates.

# For means we can see that mean_words_old = mean_Sold + mean_FWwords:Sold
round(tab_mu_mapped[,mp_names$mean[-c(1:2)]],2)

# Similarly sd_words_old = exp(sd_Sold) and sd_words_old = exp(sd_Sold + sd_FWwords:Sold)
round(tab_mu_mapped[,mp_names$sd[-c(1:2)]],2)

### Look at parameter estimates: variance

tab_var <- plot_density(samples,selection="variance",layout=c(2,7))
# These estimates reflect indivdiual differences
round(tab_var[,sp_names$mean],2)
round(tab_var[,sp_names$sd],2)

### Look at parameter estimates: correlation

tab_cor <- plot_density(samples,selection="correlation",layout=c(4,4))
# There are 91 correlations (14*13/2). Much of the correlation reflects the 
# design matrix structure.


### Look at parameter estimates: alpha


# Individual participant plots show the prior implied by the population model,
# providing an indication of shrinkage effects. 
tab_alpha <- plot_density(samples,selection="alpha",layout=c(2,7))

# Table of parameters is a list, can look at elements as with mu, e.g., 
round(tab_alpha[[subject_names(samples)[1]]][,sp_names$mean],3)

# As for mu can look at mapped parameters
tab_alpha_mapped <- plot_density(samples,selection="alpha",
                                 layout=c(2,7),mapped=TRUE)
# As expected mean_words_old = .867 + .655
round(tab_alpha_mapped[[subject_names(samples)[1]]][,mp_names$mean[-c(1:2)]],2)


#### Testing population parameter estimates ----

# Suppose we want to test if d'(words) > d'(faces) in the population. This is 
# just the same as testing mu parameter mean_FWwords:Sold > 0. We can do with 
# the p_test function, which acts like a t-test (if we want to compare to 
# something other than zero specify that with the mu argument).
p_test(samples,x_name="mean_FWwords:Sold",selection = "mu")

# The attribute of the table is the probability that samples less than zero
# are observed. Just as in a t-test the complimentary probability is available.
p_test(samples,x_name="mean_FWwords:Sold",selection = "mu",alternative="greater")

# We can make the same test (words > faces) for log(sd)
p_test(samples,x_name="sd_FWwords:Sold",selection = "mu")

# We can present the results for sd instead of log(sd) by specifying a fun argument.
# In this case p_name is just used to name the quantity being tested. 
p_test(samples,x_name="sd",selection = "mu",
       x_fun=function(x){exp(x["sd_FWwords:Sold"])})

# We can also use the fun argument to combine different parameters, here adding
# d'(words) - d'(faces) to d'(faces) to get d'(words) then testing if that is
# greater than 1.
p_test(samples,x_name="d\'(words)",selection = "mu",mu=1,
       x_fun=function(x){sum(x[c("mean_Sold","mean_FWwords:Sold")])})

# The same results is obtained by looking at mapped parameters (NB: mapped 
# analyses are only available for mu and alpha, see below).
p_test(samples,mapped=TRUE,x_name="mean_words_old",selection = "mu",mu=1)

# Turning to thresholds we get estimates for the first face and word threshold
round(tab_mu[,5:6],2)

# To test if they credibly differ:
p_test(samples,selection = "mu",x_name="threshold",
       x_fun=function(x){diff(x[c("threshold","threshold_FWwords")])})

# Alternately we can table all elements by entering the two thresholds through
# separate x and y arguments. 
p_test(x=samples,x_name="threshold_FWwords",y=samples,y_name="threshold",
       selection = "mu")

# Tests can also be performed on population variance (individual difference)
# estimates. For example: 
p_test(x=samples,x_name="mean_FWwords:Sold",y=samples,y_name="mean_Sold",
       selection = "variance")

# Finally, tests can be performed on correlations.
p_test(samples,x_name="threshold_FWwords.threshold",selection = "correlation")


#### Testing individual participants
# We could also do this for individual participants by testing the alpha and 
# specifying a subject name (if not specified tests first, here we give the
# first name explicitly).
p_test(samples,x_name="mean_FWwords:Sold",selection = "alpha",
       x_subject=subject_names(samples)[1])

# We could also compare two subjects, here showing the second subject has a
# smaller d' for faces than the first.
p_test(x=samples,y=samples,
       x_name="mean_Sold",selection = "alpha",
       x_subject=subject_names(samples)[2],
       y_subject=subject_names(samples)[1])
# but a bigger increase for faces than the first.
p_test(x=samples,y=samples,
       x_name="mean_FWwords:Sold",selection = "alpha", 
       x_subject=subject_names(samples)[2],
       y_subject=subject_names(samples)[1])

# Functions can also be specified
p_test(x=samples,y=samples,
       x_fun=function(x){exp(x["sd_FWwords:Sold"])},
       y_fun=function(x){exp(x["sd_FWwords:Sold"])},
       x_name="sd",selection = "alpha", 
       x_subject=subject_names(samples)[2],
       y_subject=subject_names(samples)[1])



#### Parameter recovery study ----
print(load("models/SDT/Probit/examples/samples/probitFW.RData")) 

# Create a single simulated data set.

# Can make up new data using the mean of the alphas (median or random also possible) 
new_dat <- post_predict(samples,use_par="mean",n_post=1)
# Or by sampling a parameter vector from the hyper
new_dat_hyper <- post_predict(samples,hyper=TRUE,n_post=1)

# An attribute keeps the pars used to create
round(attr(new_dat,"pars"),2)
round(attr(new_dat_hyper,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="standard")
# save(samplers,file="RecoveryProbitFixed.RData")
# run in runRecoveryProbitFixed.R
print(load("models/SDT/Probit/examples/samples/RecoveryProbitFixed.RData"))
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plot_density(samples,selection="alpha",filter="burn",layout=c(2,7),pars=pars)
# Some shrinkage but not bad
plot_alpha_recovery(tabs,layout=c(2,7))
plot_alpha_recovery(tabs,layout=c(2,7),do_rmse=TRUE,do_coverage=TRUE)


# can we recover these?
samplers <- make_samplers(new_dat_hyper,design,type="standard")
# save(samplers,file="RecoveryProbitRandom.RData")
# run in runRecoveryProbitRandom.R
print(load("models/SDT/Probit/examples/samples/RecoveryProbitRandom.RData"))
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plot_density(samples,selection="alpha",filter="burn",layout=c(2,7),pars=pars)
# Some shrinkage but not bad
plot_alpha_recovery(tabs,layout=c(2,7))
plot_alpha_recovery(tabs,layout=c(2,7),do_rmse=TRUE,do_coverage=TRUE)



