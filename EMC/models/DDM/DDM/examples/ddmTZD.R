#### Setup ----
rm(list=ls())
source("emc/emc.R")
source("models/DDM/DDM/ddmTZD.R")
# NB: The "TZD" parameterization defined relative to the "rtdists" package is:
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
  #   0 < DP < 1: rtdists d = t0(upper)-t0(lower) = (2*DP-1)*t0


#### Format the data to be analyzed ----

print(load("Data/PNAS.RData"))
# Note that this data was censored at 0.25s and 1.5s

dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

#### Explore the data ----

# 3 x 2 design, 19 subjects
lapply(dat,levels)
# $subjects
#  [1] "as1t" "bd6t" "bl1t" "hsft" "hsgt" "kd6t" "kd9t" "kh6t"
#  [9] "kmat" "ku4t" "na1t" "rmbt" "rt2t" "rt3t" "rt5t" "scat"
# [17] "ta5t" "vf1t" "zk1t"
# 
# $E
# [1] "accuracy" "neutral"  "speed"   
# 
# $S
# [1] "left"  "right"
# 
# $R
# [1] "left"  "right"
# 
# $rt
# NULL

# 800 trials per subject
table(dat$subjects)
# as1t bd6t bl1t hsft hsgt kd6t kd9t kh6t kmat ku4t na1t rmbt rt2t 
#  810  849  843  848  849  849  849  837  846  845  848  831  842 
# rt3t rt5t scat ta5t vf1t zk1t 
#  843  691  849  845  838  806
 
# 70-90% accuracy, .35s - .5s mean RT
plot_defective_density(dat,factors=c("E","S"),layout=c(2,3))


#### Set up the design ----

# Test E factor with Accuracy - Neutral &  Accuracy - Speed contrasts
Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))

# Test stimulus factor with intercept and right-left factor. When applied to rate 
# parameter v_S is the traditional DDM "drift rate" parameter. The 
# intercept term is drift bias. If it is zero drift rate is the same for left
# and right, if positive rate bias favors right, when negative left, e.g.,
# intercept v = 1, v_S = 2 implies an upper rate of 3 and lower rate of -1. 
Vmat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,""))  

######## Wiener diffusion model ----

# Fit a Wiener diffusion model where between-trial variability parameters are set
# to zero in "constants" (NB: this is done on the sampled scale, so qnorm (probit)
# for DP and SZ and log for st0 and sv) with the traditional characterization of
# the speed vs. accuracy emphasis factor (E) selectively influencing threshold (a), 
# and stimulus (S) affecting rate (v). In order to make the model identifiable
# we set moment-to-moment variability to the conventional value of 1. 

design_a <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(S=Vmat,E=Emat),
  Flist=list(v~S,a~E,sv~1, t0~1, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),st0=log(0),DP=qnorm(0.5),SZ=qnorm(0),sv=log(0)),
  model=ddmTZD)

# This produces a 7 parameter model
sampled_p_vector(design_a)

# v is the rate intercept, and v_S the traditional drift rate, a is the threshold
# for accuracy, a_Ea-n accuracy-neutral and a_Ea-s accuracy - neutral. t0 is the
# mean non-decision time and Z the proportional start-point bias (unbiased Z=0.5).

# We will fit a "standard" hierarchical model that allows for population 
# correlations among all parameters (7 x (7-1) /2 = 21 correlations), along with 
# population mean (mu) and variance estimates for each parameter. 

# Here we use the default hyper-prior, reproduced here explicitly (the same
# is obtained if the prior argument is omitted from make_samplers). It sets a
# mean of zero and a variance of 1 for all parameters and assumes they are 
# independent. 
prior=list(theta_mu_mean = rep(0, 7), theta_mu_var = diag(rep(5, 7)))

# make_samplers combines the data and design and makes a list of n_chains
# pmwg objects ready for sampling. Here we use the default 3 chains. Each will
# be sampled independently (multiple chains are useful for assessing convergence).
samplers <- make_samplers(dat,design_a,type="standard",
                          prior_list=prior,n_chains=3,rt_resolution=.02)
# make_samplers also compresses the data, grouping together trials with the same 
# parametersand rt's that differ less than the value specified in the rt_resolution 
# argument. Here we use the default which assumes a seconds scale (using 
# seconds is STRONGLY recommended) appropriate to visual stimuli presented on a
# monitor with 50Hz refresh (i.e., refreshing each 0.02s). Uniform random error
# in rt measurement is introduced in such a setup unless timing is coordinated
# with the sync pulse. Response devices such as the mouse or keyboard buttons
# introduce further error, so this setting is quite conservative. Taking 
# measurement resolution can substantially speed up likelihood calculation, the
# key bottleneck in sampling, in this case by a factor of 4

# save(samplers,file="sPNAS_a.RData")

#### Explore the fitting script ----  

#  Fitting is performed in the script sPNAS_a.R (e.g., on a linux based system
# the command line is R CMD BATCH sPNAS_a.R & to run it in background)

# We use three procedures that run the "burn", "adapt" and "sample" until they
# produce samples with the required characteristics. By default each chain gets 
# its own cores (this can be set with the cores_for_chains argument) and the 
# cores_per_chain argument above gives each chain 2 cores, so 6 are used in 
# total. Here are the functions called by sPNAS_a.R
# 
# The first "burn" stage by default runs 500 iterations, discards the first 200 
# and repeatedly trys adding new iterations (and possibly removing initial)
# iterations until R hat is less than a criterion (1.1 by default) for all
# random effect parameters in all chains. The aim of this stage is to find
# the posterior mode and get chains suitable for the next "adapt" stage.
#
# sPNAS_a <- auto_burn(samplers,cores_per_chain=2)
#
# The "adapt" stage develops and approximation to the posterior that will make
# sampling more efficient (less autocorrelated) in the final "sample" stage.
#
# sPNAS_a <- auto_adapt(sPNAS_a,cores_per_chain=2)
#
# Here in the final sample stage we ask for 1000 iterations per chain.
#
# sPNAS_a <- auto_sample(sPNAS_a,iter=1000,cores_per_chain=2)
# 
# Once sampling is completed the script also gets posterior predictive samples
# to enable model fit checks. By default this is based on randomly selecting 
# iterations from the final (sample) stage, and provides posterior predictives 
# for the random effects. Here we use one core pre participant.
# ppPNAS_a <- post_predict(sPNAS_a,n_cores=19)


#### Load Wiener fit results ----
print(load("models/DDM/DDM/examples/samples/sPNAS_a.RData")) 

#### Check convergence ----

# We can check the state of samplers
chain_n(sPNAS_a)

# Lets first look at the burn samples
plot_chains(sPNAS_a,selection="LL",layout=c(4,5),filter="burn")
par(mfrow=c(2,7))
plot_chains(sPNAS_a,selection="alpha",layout=NULL,filter="burn")
# R hat indicates mostly good mixing
gd_pmwg(sPNAS_a,selection="alpha",filter="burn")
# Default shows multivariate version over parameters. An invisible return
# provides full detail, here printed and rounded, again very good
round(gd_pmwg(sPNAS_a,selection="alpha",filter="burn",print_summary = FALSE),2)


# Focus on the sample stage from here (the default setting of filter) 

# RANDOM EFFECTS (i.e., subject level)
# Participant likelihoods all fat flat hairy caterpillars
plot_chains(sPNAS_a,selection="LL",layout=c(4,5))
# Plot random effects (default selection="alpha"), again they look good
par(mfrow=c(2,7)) # one row per participant
plot_chains(sPNAS_a,selection="alpha",layout=NULL)

# MIXING 
# R hat indicates excellent mixing
round(gd_pmwg(sPNAS_a,selection="alpha",print_summary = FALSE),2)

# SAMPLING EFFICIENCY
# Actual number samples = 3 chains x 533 per chain = 1599. 
# The average effective number shows autocorrelation is quite low
round(es_pmwg(sPNAS_a,selection="alpha"))
# Sometimes you might want to look at worst case summary
round(es_pmwg(sPNAS_a,selection="alpha",summary_alpha=min))
# To get per subject details
round(es_pmwg(sPNAS_a,selection="alpha",summary_alpha=NULL))
# Integrated autocorrelation time (IAT) provides a summary of efficiency, a
# value of 1 means perfect efficiency, larger values indicate the approximate 
# factor by which iterations need to be increased to get a nominal value 
# i.e., Effective size ~ True size / IAT. 
iat_pmwg(sPNAS_a,selection="alpha")


# POPULATION EFFECTS
# Similar analyses as above for random effects

# Population mean
plot_chains(sPNAS_a,selection="mu",layout=c(2,4))
round(gd_pmwg(sPNAS_a,selection="mu"),2) 
round(es_pmwg(sPNAS_a,selection="mu"))
iat_pmwg(sPNAS_a,selection="mu")
# Can also print autocorrelation functions for each chain (can also be done for
# alpha)
par(mfrow=c(3,7))
plot_acfs(sPNAS_a,selection="mu",layout=NULL)

# Population variance
plot_chains(sPNAS_a,selection="variance",layout=c(2,4))
round(gd_pmwg(sPNAS_a,selection="variance"),2)
round(es_pmwg(sPNAS_a,selection="variance"))
iat_pmwg(sPNAS_a,selection="variance")

# There are p*(p-1)/2 correlations, where p = number of parameters.  
plot_chains(sPNAS_a,selection="correlation",layout=c(3,7),ylim=c(-1,1))
# All are estimated quite well without strong autocorrelation.
round(gd_pmwg(sPNAS_a,selection="correlation"),2)
round(es_pmwg(sPNAS_a,selection="correlation"))
iat_pmwg(sPNAS_a,selection="correlation")


####  Fit ----

# By default the plot shows results for all subjects, putting everyone on the
# same x scale, which can make it hard to see fit for some or most subjects
plot_fit(dat,ppPNAS_a,layout=c(2,3))

# You could specify your own limits (e.g,. the data range) and move the legend)
plot_fit(dat,ppPNAS_a,layout=c(2,3),xlim=c(.25,1.5),lpos="right")

# Or plot for a single subject, e.g., the first 
plot_fit(dat,ppPNAS_a,layout=c(2,3),subject="as1t")

# This function (note the "s" in plot_fits) does the x scaling per subject
plot_fits(dat,ppPNAS_a,layout=c(2,3),lpos="right")

# Can also show the average over subjects as subjects is like any other factor,
# so just omit it. We see that the fit is OK but has various misses. 
plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))

# We can also use plot_fit to examine specific aspects of fit by supplying a 
# function that calculates a particular statistic for the data. For example to
# look at accuracy we might use (where d is a data frame containing the observed
# data or or posterior predictives) this to get percent correct
pc <- function(d) 100*mean(d$S==d$R)
# Drilling down on accuracy we see the biggest misfit is in speed for right
# responses by > 5%.
plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),
         stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))

# Conversely for mean RT the biggest misses are over-estimation of RT by 50ms
# or more in the accuracy and neutral conditions
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),
         stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))
round(tab,2)

# However for fast responses (10th percentile) there is global under-prediction.
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.375),
         stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")
round(tab,2)

# and for slow responses (90th percentile) global over-prediction .
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
         stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")
round(tab,2)

# This corresponds to global under-estimation of RT variability.
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),
         stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))
round(tab,2)


# Errors tend to be much faster than the models for all conditions.
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
         stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")
round(tab,2)


#### Posterior parameter inference ----

### Population mean (mu) tests

# Here we use plot_density without plotting to get a table of 95% parameter CIs
ciPNAS <- plot_density(sPNAS_a,layout=c(2,4),selection="mu",do_plot=FALSE)
round(ciPNAS,2)

# We will also do testing of the "mapped" parameters, that is, parameters 
# transformed to the scale used by the model and mapped to the design cells. 
# To illustate here are the mapped posterior medians
mapped_par(ciPNAS[2,],design_a)


# Testing can be performed with the p_test function, which is similar to the 
# base R t.test function. By default its assumes selection="mu". 
# If just one "x" is specified it tests whether a single 
# parameter differs from zero, reproducing one of the credible intervals in the
# table above, but with the probability of being less than zero. 

## Rates

# There is a slight drift-rate bias to left (lower) but not credible at 95%
p_test(x=sPNAS_a,x_name="v")
# In this case it might be better to get the probability of greater than zero:
p_test(x=sPNAS_a,x_name="v",alternative = "greater")
# The main effect of stimulus (i.e., the conventional DDM drift rate) is 
# clearly greater than zero.
p_test(sPNAS_a,x_name="v_S")

## Non-decision time

# Recall non-decision time is on a log scale
p_test(sPNAS_a,x_name="t0")
# We can also test transformed parameters by specifying a function. The x_name
# is just used in the table header. 
p_test(sPNAS_a,x_fun=function(x){exp(x["t0"])},x_name="t0(sec)")
# Alternatively p_test can automatically transform to the scale used by the
# model, here producing the same result.
p_test(x=sPNAS_a,x_name="t0",mapped=TRUE)

## Start-point bias 

# Here we both transform and then test again the unbiased level (0.5). There is 
# a very close to (two tailed 95%) credible bias to left (lower), to get more
# resolution we up the default digits.
p_test(sPNAS_a,x_name="Z",x_fun=function(x){pnorm(x["Z"])},
       mu=0.5,digits=4,alternative = "greater")

## Threshold

# Transforming thresholds to the natural scale
p_test(sPNAS_a,x_name="a",x_fun=function(x){exp(x["a"])})
# Accuracy threshold is credibly but not much greater than neutral
p_test(sPNAS_a,x_name="a_Ea-n")
# Accuracy threshold is much greater than speed.
p_test(sPNAS_a,x_name="a_Ea-s")

# Using mapped=TRUE not only transforms to the scale used by the model but also
# to the cells of the design, automatically generating appropriate names which 
# can be seen using this function.
p_names(sPNAS_a,mapped=TRUE)$a
# For example we could test if the threshold for speed differs from zero
p_test(x=sPNAS_a,x_name="a_speed",mapped=TRUE)
# mapped=TRUE also allows arbitrary comparisons between cells to be done. For 
# example, we can compare the difference between neutral and speed.
p_test(x=sPNAS_a,x_name="n-s",mapped=TRUE,x_fun=function(x){
  diff(x[c("a_speed","a_neutral")])})
# The same thing can be achieved by specifying and x and y argument
p_test(x=sPNAS_a,x_name="a_neutral",y=sPNAS_a,y_name="a_speed",mapped=TRUE)
# More generally we can test differences between combinations of cells using a
# fun. For example, is the average of accuracy and neutral different from speed?
p_test(x=sPNAS_a,x_name="an-s",mapped=TRUE,x_fun=function(x){
  sum(x[c("a_accuracy","a_neutral")])/2 - x["a_speed"]})

### Population variability 
# Here we cannot use mapping, as these parameters are about variability 
# on the sampled scale, but we can make tests.

ciPNAS <- plot_density(sPNAS_a,layout=c(2,4),selection="variance",do_plot=FALSE)
round(ciPNAS,3)

# Suppose we wanted to compare standard deviations of the accuracy-neutral and
# accuracy minus speed contrasts
p_test(x=sPNAS_a,x_name="a_Ea-n",x_fun=function(x){sqrt(x["a_Ea-s"])},
  y=sPNAS_a,y_name="a_Ea-s",y_fun=function(x){sqrt(x["a_Ea-n"])},selection="variance")

### Population correlation
# Again mapped cant be tested. There are lots of these so lets look at them first
# overall.

ciPNAS <- plot_density(sPNAS_a,layout=c(2,4),selection="correlation",do_plot=FALSE)
round(ciPNAS,3)

# The strongest correlations are among the a parameters. Do the two 
# intercept-effect correlations differ?
p_test(x=sPNAS_a,x_name="a_Ea-n.a",y=sPNAS_a,y_name="a_Ea-s.a",selection="correlation")
# Their imprecise estimation makes any difference not credible.

# It is also important to note that the intercept and effect contrasts are themselves
# not orthogonal, which will induce a correlation. This can be seen from the dot
# products of their design matrix being non-zero

get_design_matrix(sPNAS_a)$a

# The two effects are orthogonal so this is not a factor in their positive correlation
p_test(sPNAS_a,x_name="a_Ea-s.a_Ea-n",selection="correlation")

# Otherwise the table above reveals no credible correlations.


### Random effect (alpha) tests

# We can test individual participants, with the first one selected by default 
p_test(x=sPNAS_a,x_fun=function(x){exp(x["t0"])},x_name="t0(sec)",
       selection="alpha",digits=3)

# We can also compare two participants, for example bd6t has slightly but 
# credibly slower non-decision time than as1t.
snams <- subject_names(sPNAS_a)
p_test(selection="alpha",digits=3,
  x=sPNAS_a,x_fun=function(x){exp(x["t0"])},x_name="t0(s)",x_subject=snams[2],
  y=sPNAS_a,y_fun=function(x){exp(x["t0"])},y_name="t0(s)",y_subject=snams[1])
# Once again we can use mapping to the same effect
p_test(selection="alpha",digits=3,mapped=TRUE,
  x=sPNAS_a,x_name="t0",x_subject=snams[2],
  y=sPNAS_a,y_name="t0",y_subject=snams[1])

### Extracting a data frame of parameters

# Parameters can be extracted into a data frame for further analysis, by default
# getting mu from the sample stage
head(parameters_data_frame(sPNAS_a))
# If desired constants can be included
head(parameters_data_frame(sPNAS_a,include_constants=TRUE))
# Or mapped parameters
head(parameters_data_frame(sPNAS_a,mapped=TRUE))
# If random effects are selected a subjects factor column is added
head(parameters_data_frame(sPNAS_a,selection="alpha",mapped=TRUE))


########  FULL DDM ----   

# Once trial-to-trial variability parameters are estimated the sampler can 
# sometimes explore regions where the numerical approximation to the DDM's 
# likelihood become inaccurate. In order to avoid this log_likelihood_ddm (in 
# emc/likelihood.R) imposes the following restrictions on these parameters
# (we also restrict v and a to typical regions), returning low likelihoods when
# they are violated. 
#     abs(v)<5 | a<2 | sv<2 | sv>.1 | SZ<.75 | SZ>.01 | st0<.2
# They are usually appropriate for the seconds scale with s=1 fixed. For 
# some data and/or different scalings they may have to be adjusted. Posterior
# estiamtes should be checked to see if estimates are stacking up against these
# bounds, indicating they might need adjustment. Alternatively, poor convergence
# behaviour may indicate the need for further restriction.

# Lets fit the full model analogous to the Wiener model fit previously.
design_a_full <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(S=Vmat,E=Emat),
  Flist=list(v~S,a~E,sv~1, t0~1, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),DP=qnorm(0.5)),
  model=ddmTZD)

# samplers <- make_samplers(dat,design_a_full,type="standard")
# save(samplers,file="sPNAS_a_full.RData")
# Fits are run by sPNAS_a_full.R 

#### Load Full DDM fitting results ----
print(load("models/DDM/DDM/examples/samples/sPNAS_a_full.RData")) 

#### Convergence ----

# As we start to look at more models, and more complex models, it is 
# convenient to be able to run a standard set of checks, as follows. By
# default this function looks at the sample stage (use filter="burn" to look at
# the burn stage) and pauses after printing out results for each type of
# parameter (interactive=FALSE to turn this off), printing out gd, iat and es.
check_run(sPNAS_a_full)
# It also saves a pdf (by default pdf_name"check_run.pdf", width and height 
# arguments can be used to change default page size) of chain plots (by
# default layout = c(3,4)). These results together show that all of the sample
# stage is well converged, and that sv and particularly SZ are not efficiently
# sampled. Clearly if we wished to do inference with SZ we should get more 
# samples to obtain a sufficient effective size. Note that printed outputs are
# sorted so it is easy to pickup extreme cases, and that both the minimum and 
# mean are used to summarize the effective size of random effects.

#### Fit ----

# Individual fits are improved over the Wiener model 
plot_fit(dat,ppPNAS_a_full,layout=c(2,3))

# The average over subjects shows the improvement very clearly. 
plot_fit(dat,ppPNAS_a_full,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))

# Drilling down on accuracy misfit is much reduced, although it is still present
# for speed for right responses. This suggests a model with Z~E may be warrented.
pc <- function(d) 100*mean(d$S==d$R)
plot_fit(dat,ppPNAS_a_full,layout=c(2,3),factors=c("E","S"),
         stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))

# Mean RT is well estimated.  
tab <- plot_fit(dat,ppPNAS_a_full,layout=c(2,3),factors=c("E","S"),
         stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))

# However for fast responses (10th percentile) under-prediction remains except for speed.
tab <- plot_fit(dat,ppPNAS_a_full,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.375),
         stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")

# and for slow responses (90th percentile) clear over-prediction except for speed.
tab <- plot_fit(dat,ppPNAS_a_full,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
         stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")
round(tab,2)

# This corresponds to under-estimation of RT variability for speed and 
# otherwise over-estimation.
tab <- plot_fit(dat,ppPNAS_a_full,layout=c(2,3),factors=c("E","S"),
         stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))

# Error speed is well estimated except for speed where it is over-estimated.
tab <- plot_fit(dat,ppPNAS_a_full,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
         stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")

#### Posterior parameter inference ----

### Population mean (mu) tests

# Priors are dominated, even for sv and SZ 
ciPNAS_full <- plot_density(sPNAS_a_full,layout=c(2,5),selection="mu")

# Comparing with the Wiener we see higher rates and overall reduced thresholds 
round(ciPNAS_full,2)
round(ciPNAS,2)

# Priors also clearly dominated in mapped parameters
ciPNAS_full_mapped <- plot_density(sPNAS_a_full,layout=c(2,5),selection="mu",mapped=TRUE)

## Rates

# There is a slight drift-rate bias to left (lower) but not credible at 95%
p_test(x=sPNAS_a_full,x_name="v")
# The main effect of stimulus is clearly greater than zero.
p_test(sPNAS_a_full,x_name="v_left",mapped=TRUE)
p_test(sPNAS_a_full,x_name="v_right",mapped=TRUE)
# Rate variability is quite substantial (>20% of the mean)
p_test(x=sPNAS_a_full,x_name="sv",mapped=TRUE)

## Non-decision time 

# t0 is now a lower bound so slightly less than in Wiener
p_test(x=sPNAS_a_full,x_name="t0",mapped=TRUE)
# Variability is quite substantial
p_test(x=sPNAS_a_full,x_name="st0",mapped=TRUE)

## Start-point bias

# Small lower (left) bias
p_test(x=sPNAS_a_full,x_name="Z",mapped=TRUE)
# Small level of start-point variability (sz/2 is 11% of .47, i.e., sz = 0.1)
p_test(x=sPNAS_a_full,x_name="SZ",mapped=TRUE)

## Threshold

# Although the overall threshold is reduced the effects are a bit bigger
p_test(sPNAS_a_full,x_name="a",x_fun=function(x){exp(x["a"])})
p_test(sPNAS_a_full,x_name="a_Ea-n")
p_test(sPNAS_a_full,x_name="a_Ea-s")

# For example we could test if the threshold for speed differs from zero
p_test(x=sPNAS_a_full,x_name="a_speed",mapped=TRUE)
# difference between neutral and speed.
p_test(x=sPNAS_a_full,x_name="n-s",mapped=TRUE,x_fun=function(x){
  diff(x[c("a_speed","a_neutral")])})
# Average of accuracy and neutral vs. speed
p_test(x=sPNAS_a_full,x_name="an-s",mapped=TRUE,x_fun=function(x){
  sum(x[c("a_accuracy","a_neutral")])/2 - x["a_speed"]})

### Variance/correlation

ciPNAS <- plot_density(sPNAS_a_full,layout=c(2,4),selection="variance",do_plot=FALSE)
round(ciPNAS,3)

ciPNAS <- plot_density(sPNAS_a_full,layout=c(2,4),selection="correlation",do_plot=FALSE)
round(ciPNAS,3)

#### Parameter recovery study ----
# The full DDM model is notoriously difficult to sample. Here we check how well 
# its parameters can be recovered in the present design. To do so we create a 
# single simulated data set from that alpha means.
new_dat <- post_predict(sPNAS_a_full,use_par="mean",n_post=1)

# can we recover these?
samplers <- make_samplers(new_dat,design_a_full,type="standard")
# save(samplers,file="RecoveryDDMfull.RData")
# run in RecoveryDDMfull.R
print(load("models/DDM/DDM/examples/samples/RecoveryDDMfull.RData"))

tabs <- plot_density(DDMfull,selection="alpha",layout=c(2,5),mapped=TRUE,
                     pars=attributes(attr(DDMfull,"data_list")[[1]])$pars)
# Poor for sv and terrible for SZ
plot_alpha_recovery(tabs,layout=c(2,5))
# Coverage is fairly decent even in sv and SZ
plot_alpha_recovery(tabs,layout=c(2,5),do_rmse=TRUE,do_coverage=TRUE)

######## Full DDM and threshold, rate, and non-decision time effects of emphasis ----

# In this example we also demonstrate "cell" coding, where there is a separate
# estimate for each cell of the E x S design used for rates. This is not easily
# achieved with the linear model language applied to the E and S factors, so 
# we derive a new factor "SE" with 6 levels that combines their levels. To do
# so we use a function that will create the factors from the Ffactors.
se <- function(d) {factor(paste(d$S,d$E,sep="_"),levels=
  c("left_accuracy","right_accuracy","left_neutral","right_neutral","left_speed","right_speed"))}
# NB1: Although usually not necessary it is good practice to explicitly define 
#      the levels of factors created in this way. 
# NB2: Ffunctions can be used to define arbitrary new factors in this way, making
#      model specification very flexible.

# To archive cell coding we have to also remove the intercept when defining the 
# rate equation (0+SE and SE-1 both work)
# NB: Contrasts can be defined for Ffunctions factors. To achieve cell coding
#     we dont need to do this as it will occur with the default contr.treatment
#     contrast, but here we define it explicitly.
design_avt0_full <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(SE=diag(6)),
  Flist=list(v~0+SE,a~E,sv~1, t0~E, st0~1, s~1, Z~1, SZ~1, DP~1),
  Ffunctions = list(SE=se),
  constants=c(s=log(1),DP=qnorm(0.5)),
  model=ddmTZD)

# samplers <- make_samplers(dat,design_avt0_full,type="standard")
# save(samplers,file="sPNAS_avt0_full.RData")

#### Load full DDM with rate and t0 emphasis effects ----
print(load("models/DDM/DDM/examples/samples/sPNAS_avt0_full.RData")) 

#### Check convergence ----

# In this more complicated case there is some initial non-stationarity in over 
# the first 500 iterations of the sample stage. 
check_run(sPNAS_avt0_full)
# Hence we re-run on only the last 1000. This looks good but lower efficiency
# inference.
check_run(sPNAS_avt0_full,subfilter=1000)

####  Fit ----

# post predict did not use first 1000
plot_fit(dat,ppPNAS_avt0_full,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))
plot_fits(dat,ppPNAS_avt0_full,layout=c(2,3),lpos="right")


pc <- function(d) 100*mean(d$S==d$R)
# Excellent fit in all cases.
plot_fit(dat,ppPNAS_avt0_full,layout=c(2,3),factors=c("E","S"),
         stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))

# Mean RT also excellent
tab <- plot_fit(dat,ppPNAS_avt0_full,layout=c(2,3),factors=c("E","S"),
         stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))

# However for fast responses (10th percentile) still some under-prediction except for speed.
tab <- plot_fit(dat,ppPNAS_avt0_full,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.375),
         stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")
round(tab,2)

# but for slow responses (90th percentile) all good.
tab <- plot_fit(dat,ppPNAS_avt0_full,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
         stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")

# and RT variability now also good.
tab <- plot_fit(dat,ppPNAS_avt0_full,layout=c(2,3),factors=c("E","S"),
         stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))


# Also errors speed also good, with some small residual over-estimation in speed.
tab <- plot_fit(dat,ppPNAS_avt0_full,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
         stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")
round(tab,2)

#### Posterior parameter inference ----

### Population mean (mu) tests

# Clearly the default N(0,1) priors are somewhat inappropriate for the large
# positive and negative rates produced by "cell" coding.
ciPNAS_avt0_full <- plot_density(sPNAS_avt0_full,layout=c(3,6),selection="mu",subfilter=500)

# Cell coding directly gives us the traditional drift rate estimates
round(ciPNAS_avt0_full,2)

# and (by definition) the same results when mapped 
ciPNAS_avt0_full_mapped <- plot_density(sPNAS_avt0_full,layout=c(3,6),
                                        selection="mu",mapped=TRUE,subfilter=500)
round(ciPNAS_avt0_full_mapped,2)

### t0 E effect

# Looking at the t0 effect of E there is no evidence of a accuracy-netural difference
p_test(x=sPNAS_avt0_full,x_name="t0_accuracy",y=sPNAS_avt0_full,y_name="t0_neutral",
       subfilter=500,mapped=TRUE)
# But speed is clearly less by ~40ms
p_test(x=sPNAS_avt0_full,subfilter=500,mapped=TRUE,x_name="t0: av(acc/neut)-speed",
       x_fun=function(x){sum(x[c("t0_accuracy","t0_neutral")])/2 - x["t0_speed"]})

# v E effect
# The pattern of E effects on rates suggests lower rates for speed, perhaps 
# indicative of reduced selective attention and hence reduced discriminative
# quality of the evidence being processed. 

fun <- function(x)
  mean(abs(x[c("v_left_accuracy","v_left_neutral","v_right_accuracy",
               "v_right_neutral")])) - mean(abs(x[c("v_left_speed","v_right_speed")]))

p_test(x=sPNAS_avt0_full,subfilter=500,mapped=TRUE,digits=3,
       x_name="v: av(acc/neut)-speed",x_fun=fun)

# Although the effect is (barely) credible, one might also suspect over-fitting.
# In this case it would be useful to also fit models with selective influence of 
# E on a and t0 vs. a and v. This is left as an exercise.  


#### Model selection ----

# Which model gives the best trade off between goodness of fit and simplicity?

# Information criteria (IC) approaches are based purely on posterior samples. 
# This function produces two types, DIC and BPIC (the latter gives greater 
# weight to simplicity), where smaller is better. 

# It also provides estimates of components of DIC and BPIC :
# 1) EffectiveN: the "effective" number of parameters (which is usually less  
#    than the actual number in hierarchical models)
# 2) meanD: the mean posterior deviance (which can also be used for model
#    selection, but imposing only a weak complexity penalty)
# 3) Dmean: the deviance of the posterior parameter mean, which is a measure
#    of goodness of fit
# 4) minD: an estimates of the minimum posterior deviance, the smaller of Dmean 
#    and the deviance for all samples (note that by default DIC and BPIC uses 
#    this value, base results on Dmean set use_best_fit=FALSE)
pmwg_IC(sPNAS_a)
pmwg_IC(sPNAS_a_full)
pmwg_IC(sPNAS_avt0_full,subfilter=500)
# NB: The nonsensical negative effective parameter count for the Wiener model is
#     conventionally considered not to be a problem, but underlines that 
#     EffectiveN is best interpreted only in  a relative sense. The other values
#     are much more sensible and relative to the actual number of 19*10 = 190 
#     and 19*16=304 alpha parameters (random effects being relevant here as 
#     all of these calculations are based on the alpha posterior likelihoods),
#     and are consistent with hierarchical shrinkage effects.

# We see that the Wiener diffusion wins hugely in DIC/BPIC/meanD, despite its 
# much worse goodness of fit in Dmean and Dmin (also evident in its strong 
# qualitative failures shown graphically above). Overall this pattern suggests
# that the Wiener model be rejected.

# The following function calculates model weights, quantities on the unit 
# interval that under some further assumptions correspond to the probability 
# of a model being the "true" model. Here clearly the avt0 model wins strongly.
# NB: subfilter can be either a single digit, used for all, or a list of the
#     same length as the list of objects.
compare_IC(list(avt0=sPNAS_avt0_full,a=sPNAS_a_full),subfilter=list(500,0))

# This is also true on a per-subject basis except for one participant.
compare_ICs(list(avt0=sPNAS_avt0_full,a=sPNAS_a_full),subfilter=list(0,500))

# However, one might question whether we need both v and t0, so lets fit
# simpler models that drop one or the other.

design_at0_full <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Flist=list(v~S,a~E,sv~1, t0~E, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),DP=qnorm(0.5)),
  model=ddmTZD)

# samplers <- make_samplers(dat,design_at0_full,type="standard")
# save(samplers,file="sPNAS_at0_full.RData")

se <- function(d) {factor(paste(d$S,d$E,sep="_"),levels=
  c("left_accuracy","right_accuracy","left_neutral","right_neutral","left_speed","right_speed"))}
design_av_full <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(SE=diag(6)),
  Flist=list(v~0+SE,a~E,sv~1, t0~1, st0~1, s~1, Z~1, SZ~1, DP~1),
  Ffunctions = list(SE=se),
  constants=c(s=log(1),DP=qnorm(0.5)),
  model=ddmTZD)

# samplers <- make_samplers(dat,design_av_full,type="standard")
# save(samplers,file="sPNAS_av_full.RData")

print(load("models/DDM/DDM/examples/samples/sPNAS_at0_full.RData"))
print(load("models/DDM/DDM/examples/samples/sPNAS_av_full.RData"))

# SZ was slow to converge in mu, needed to 750 for at0 and 500 for av, so ran 
# enough samples to get 1000 left 1000 left without these. 
# All looks good with a_neutral now much more efficient
check_run(sPNAS_at0_full,subfilter=750)
check_run(sPNAS_av_full,subfilter=500,layout=c(3,5))

# avt0 still wins overall
compare_IC(list(avt0=sPNAS_avt0_full,at0=sPNAS_at0_full,av=sPNAS_av_full,a=sPNAS_a_full),
            subfilter=list(500,750,500,0))
# But now there are 5 more equivocal cases, largely favoring the at0 model.
compare_ICs(list(avt0=sPNAS_avt0_full,at0=sPNAS_at0_full,av=sPNAS_av_full,a=sPNAS_a_full),
            subfilter=list(500,750,500,0))

# Also fit a version of the winning model not using cell coding for v
design_avt0_full_nocell <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(SE=diag(6)),
  Flist=list(v~S*E,a~E,sv~1, t0~E, st0~1, s~1, Z~1, SZ~1, DP~1),
  Ffunctions = list(SE=se),
  constants=c(s=log(1),DP=qnorm(0.5)),
  model=ddmTZD)

# samplers <- make_samplers(dat,design_av_full_nocell,type="standard")
# save(samplers,file="sPNAS_av_full_nocell.RData")

# We see it makes very little difference to the DIC and BPIC
print(load("models/DDM/DDM/examples/samples/sPNAS_avt0_full_nocell.RData")) 
compare_IC(list(avt0nocell=sPNAS_avt0_full_nocell,avt0=sPNAS_avt0_full,at0=sPNAS_at0_full,av=sPNAS_av_full,a=sPNAS_a_full),
            subfilter=list(500,500,750,500,0))

