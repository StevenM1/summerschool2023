rm(list=ls())
source("emc/emc.R")
source("models/RACE/LBA/lbaB.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)
# NB: This data has been truncated at 0.25s and 1.5s

# Average rate = intercept, and rate d = difference (match-mismatch) contrast
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))


# Here we fit a series of models

# Only B affected by E
design_B <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,sv~1,B~E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)
# samplers <- make_samplers(dat,design_B,type="standard",rt_resolution=.02)
# save(samplers,file="sPNAS_B.RData")

# B, v and t0 affected by E
design_Bvt0 <- make_design(
  Ffactors=list(subjects=unique(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~E*lM,sv~1,B~lR*E,A~1,t0~E),
  constants=c(sv=log(1)),
  model=lbaB)
# samplers <- make_samplers(dat,design_Bvt0,type="standard",rt_resolution=.02)
# save(samplers,file="sPNAS_Bvt0.RData")

# Up to now we set sv=1 to establish a scale. Now we fit a series of models 
# where sv differs between matching and mismatching accumulators, with its 
# intercept fixed at 1 to establish scaling. Typically allowing this difference
# greatly improves fit and results in sv_TRUE (i.e., match) > sv_FALSE.

# First, allowing E to affect B and t0
design_Bt0_sv <- make_design(
  Ffactors=list(subjects=unique(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,sv~lM,B~lR*E,A~1,t0~E),
  constants=c(sv=log(1)),
  model=lbaB)
# samplers <- make_samplers(dat,design_Bt0_sv,type="standard",rt_resolution=.02)
# save(samplers,file="sPNAS_Bt0_sv.RData")

# Now E affects B and v
design_Bv_sv <- make_design(
  Ffactors=list(subjects=unique(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~E*lM,sv~lM,B~lR*E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)
# samplers <- make_samplers(dat,design_Bv_sv,type="standard",rt_resolution=.02)
# save(samplers,file="sPNAS_Bv_sv.RData")

# Now all three
design_Bvt0_sv <- make_design(
  Ffactors=list(subjects=unique(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~E*lM,sv~lM,B~lR*E,A~1,t0~E),
  constants=c(sv=log(1)),
  model=lbaB)
# samplers <- make_samplers(dat,design_Bvt0_sv,type="standard",rt_resolution=.02)
# save(samplers,file="sPNAS_Bvt0_sv.RData")

# The last model fails to converge (or rather some parameter were so strongly 
# autocorrelated it would have taken an impractically large run to produce 
# anything useful). Looking at estimates, most issues were around the a_n 
# contrast, which might not be surprising given this had a weak manifest
# effect. Here we fit this model again but dropping the a_n contrast for the
# E effects on v and t0, instead estimating an intercept and the difference
# between the average of accuracy and neutral and speed, a 14 parameter model. 
E_AVan_s_mat <- matrix(c(1/4,1/4,-1/2),nrow=3)
dimnames(E_AVan_s_mat) <- list(NULL,c("an-s"))

design_Bvt0_sv_NOa_n <- make_design(
  Ffactors=list(subjects=unique(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(v=list(E=E_AVan_s_mat,lM=ADmat),sv=list(lM=ADmat),
             B=list(E=Emat,lR=ADmat),t0=list(E=E_AVan_s_mat)),
  Flist=list(v~E*lM,sv~lM,B~lR*E,A~1,t0~E),
  constants=c(sv=log(1)),
  model=lbaB)
# samplers <- make_samplers(dat,design_Bvt0_sv_NOa_n,type="standard",rt_resolution=.02)
# save(samplers,file="sPNAS_Bvt0_sv_NOa_n.RData")

# Fits then run in batch mode. For the sv=1 models we used the same sort of 
# script as for the DDM. For the others we used a slightly smarter script that
# only fits the next stage if the previous one worked, e.g., for the Bvsv model:
#
# Here it tests if burn in had been done previously 
#
# if (is.null(attr(samplers,"burnt")) || is.na(attr(samplers,"burnt"))) {
#   sPNAS_Bv_sv <- auto_burn(samplers,cores_per_chain=4)
#   save(sPNAS_Bv_sv,file="sPNAS_Bv_sv.RData")
# }
#
# Here it checks if adaptation was already done
#
# if (is.null(attr(sPNAS_Bv_sv,"adapted")) && !is.na(attr(sPNAS_Bv_sv,"burnt"))) {
#   sPNAS_Bv_sv <- auto_adapt(sPNAS_Bv_sv,cores_per_chain=4)
#   save(sPNAS_Bv_sv,file="sPNAS_Bv_sv.RData")
# }
#
# Here it only moves on to sampling of adaptation was successful
#
# if (!is.null(attr(sPNAS_Bv_sv,"adapted")) && attr(sPNAS_Bv_sv,"adapted")) {
#   sPNAS_Bv_sv <- auto_sample(sPNAS_Bv_sv,iter=2000,cores_per_chain=4)
#   save(sPNAS_Bv_sv,file="sPNAS_Bv_sv.RData")
# }
#
# This was run later after choosing this as the best model and knowing that
# the first 2000 iterations of the sample run had to be discareded, giving
# 5000*3 = 15000 good samples to work with.
#
# sPNAS_Bv_sv <- auto_sample(sPNAS_Bv_sv,iter=4000,cores_per_chain=10)
# save(sPNAS_Bv_sv,file="sPNAS_Bv_sv.RData")
# ppPNAS_Bv_sv <- post_predict(sPNAS_Bv_sv,n_cores=19,subfilter=2000)
# save(ppPNAS_Bv_sv,sPNAS_Bv_sv,file="sPNAS_Bv_sv.RData")

#### Load model results ----
print(load("models/RACE/LBA/examples/samples/sPNAS_B.RData")) 
print(load("models/RACE/LBA/examples/samples/sPNAS_Bvt0.RData"))
print(load("models/RACE/LBA/examples/samples/sPNAS_Bt0_sv.RData"))
print(load("models/RACE/LBA/examples/samples/sPNAS_Bv_sv.RData"))
print(load("models/RACE/LBA/examples/samples/sPNAS_Bvt0_sv.RData")) # Failed model, forced to adapt and sample 
print(load("models/RACE/LBA/examples/samples/sPNAS_Bvt0_sv_NOa_n.RData"))


#### Check convergence ----

# Aim to get 1000 good samples as a basis for model comparison (but note this
# is not enough for inference on some parameters)
check_run(sPNAS_B)
# Slow to converge, but there by 4500
check_run(sPNAS_Bvt0,subfilter = 4500,layout=c(3,3))
# Low efficiency for t0_Ea-n, t0, and B_lRd:Ea-n
check_run(sPNAS_Bt0_sv,subfilter=2000,layout=c(3,5))
# B_lRd:Ea-n, B_Ea-n, v_Ea-n, v_Ea-n:lMd all inefficient. Some slow-wave
# oscillations, definitely need longer series if using the parameter estimates
# (even though likely stationary after 2000)
check_run(sPNAS_Bv_sv,subfilter=2000,layout=c(3,5))
# After extensive burn (although not achieving Rhat < 1.2), adaptation was quick 
# but chains were very poor. Decided not to pursue further.
check_run(sPNAS_Bvt0_sv,subfilter=1000)
# Burned in quickly and converged well after 1000, B_lRd:Ea-n and especially 
# t0_Ean-s quite inefficient
check_run(sPNAS_Bvt0_sv_NOa_n,subfilter=1000,layout=c(3,5)) 

#### Model selection ----

# Comparing all 5 models, clearly need sv (NB: use only 1000 from Bvsv after
# convergence, more added for parameter inference as this is the winning model)
compare_IC(list(B=sPNAS_B,Bvt0=sPNAS_Bvt0,Bt0sv=sPNAS_Bt0_sv,Bvsv=sPNAS_Bv_sv,
  Bvt0sv=sPNAS_Bvt0_sv_NOa_n),subfilter=list(0,4500,2000,2001:3000,1000))
ICs <- compare_ICs(list(B=sPNAS_B, Bvt0=sPNAS_Bvt0, Bt0sv=sPNAS_Bt0_sv, 
  Bvsv=sPNAS_Bv_sv, Bvt0sv=sPNAS_Bvt0_sv_NOa_n),subfilter=list(0,4500,2000,2001:3000,1000))
#
# The latter function returns a list of matrices for each participant, so can use
# that to count winners as follows. 
table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$DIC)]})))
table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$BPIC)]})))
# Mostly Bvsv winning, second Bvt0sv with no Ea_n, and one clear Bt0sv. Could
# try the simpler E contrast with the Bvsv model, but leave that as an exercise.

# Comparing with the best DDM (16 parameter) model to the best (15 parameter)
# LBA model the latter has a slight edge. This is born out in testing fit,
# both models fit well but the small overestimation of error RT in the speed
# condition evident in the DDM is no longer evident.
source("models/DDM/DDM/ddmTZD.R")
compare_IC(list(avt0=sPNAS_avt0_full,Bvsv=sPNAS_Bv_sv),subfilter=list(500,2000))

# For the remaining analysis add 4000 iterations to Bvsv model

####  Fit of winning (Bvsv) model ----

# post predict did not use first 2000, so inference based on 5000*3 samples
plot_fit(dat,ppPNAS_Bv_sv,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))
plot_fits(dat,ppPNAS_Bv_sv,layout=c(2,3),lpos="right")

# No evidence of any misfit
pc <- function(d) 100*mean(d$S==d$R)
plot_fit(dat,ppPNAS_Bv_sv,layout=c(2,3),factors=c("E","S"),
         stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
tab <- plot_fit(dat,ppPNAS_Bv_sv,layout=c(2,3),factors=c("E","S"),
         stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))
tab <- plot_fit(dat,ppPNAS_Bv_sv,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.4),
         stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")
tab <- plot_fit(dat,ppPNAS_Bv_sv,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
         stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")
tab <- plot_fit(dat,ppPNAS_Bv_sv,layout=c(2,3),factors=c("E","S"),
         stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))
tab <- plot_fit(dat,ppPNAS_Bv_sv,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
         stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")

#### Posterior parameter inference ----

### Population mean (mu) tests

# Priors all well dominated except some rate parameters
ciPNAS_Bv_sv <- plot_density(sPNAS_Bv_sv,layout=c(3,6),selection="mu",subfilter=2000)
# On the natural scale it is evident this is because the mismatch (FALSE) rates
# are least well updated, due to the fairly low error rates (errors give the
# most information about FALSE rates).
ciPNAS_Bv_sv_mapped <- plot_density(sPNAS_Bv_sv,layout=c(3,6),
                                        selection="mu",mapped=TRUE)
# Looking at parameters both with and without mapping
round(ciPNAS_Bv_sv,2)
round(ciPNAS_Bv_sv_mapped,2)

# For simpler estimates:
# 1) As is usually found with the LBA sv_true > sv_false 
# 2) t0 is quite small (100ms)
# 3) In this case start point noise (.29) is small relative to B (which has 
#    intercept exp(-.16) = 0.85, so about 1/4 of b).

### B effects

# Recall the map used with 
get_map(sPNAS_Bv_sv)$B

# B_lRd tests threshold right - threshold left, although not quite credible
# it indicates slightly higher right thresholds (i.e., a bias to respond left)
p_test(x=sPNAS_Bv_sv,x_name="B_lRd",subfilter=2000)

# B_Ea-n and B_Ea-s measure differences in response caution (i.e., thresholds
# averaged over left and right accumulators), accuracy-neutral and accuracy-speed 
# respectively. Caution for accuracy is clearly higher than speed, but not 
# credibly greater than neutral.
p_test(x=sPNAS_Bv_sv,x_name="B_Ea-n",subfilter=2000)
p_test(x=sPNAS_Bv_sv,x_name="B_Ea-s",subfilter=2000)
#
# Here we construct a test showing the processing speed
# advantage is greater for speed than neutral
p_test(x=sPNAS_Bv_sv,mapped=TRUE,x_name="average B: neutral-speed",
  x_fun=function(x){mean(x[c("B_left_neutral","B_right_neutral")]) - 
                    mean(x[c("B_left_speed","B_right_speed")])})

# No evidence of a difference between accuracy and neutral thresholds
p_test(x=sPNAS_Bv_sv,x_name="B_Ea-n",subfilter=2000)

# But thresholds clearly higher for speed (by about 0.2 on the natural scale)
p_test(x=sPNAS_Bv_sv,x_name="B_Ea-s",subfilter=2000)

# The remaining terms test interactions with bias (i.e., lR), with evidence of
# a small but just credibly stronger bias to respond left (i.e., a lower threshold
# for the left accumulator) for speed than accuracy.
p_test(x=sPNAS_Bv_sv,x_name="B_lRd:Ea-s",subfilter=1500,digits=2)

### v effects

# Again recall the map used with 
get_map(sPNAS_Bv_sv)$v

# v_Ea-n v_Ea-s indicate that processing rate (the average of matching and 
# mismatching rates) is least in accuracy, slightly greater in neutral and 
# the highest in speed.
p_test(x=sPNAS_Bv_sv,x_name="v_Ea-n",subfilter=2000)
p_test(x=sPNAS_Bv_sv,x_name="v_Ea-s",subfilter=2000)

# Here we show that processing is faster in the speed condition.
p_test(x=sPNAS_Bv_sv,mapped=TRUE,x_name="average v: speed-neutral",
  x_fun=function(x){mean(x[c("v_speed_TRUE","v_speed_FALSE")]) - 
                    mean(x[c("v_neutral_TRUE","v_neutral_FALSE")])})

# v_lMd tests the quality of selective attention, rate match - rate mismatch
p_test(x=sPNAS_Bv_sv,x_name="v_lMd",subfilter=2000)

# v_Ea-n tests if quality neutral - quality accuracy (i.e, 
# (2.55 - -0.14) - (2.71  - 0.24) = 0.22)
p_test(x=sPNAS_Bv_sv,x_name="v_Ea-n:lMd",subfilter=2000)

# The difference is ~5 times larger for accuracy - speed
p_test(x=sPNAS_Bv_sv,x_name="v_Ea-s:lMd",subfilter=2000)

# The neutral - speed difference is also highly credible, and
# can be tested in the mapped form
p_test(x=sPNAS_Bv_sv,mapped=TRUE,x_name="quality: accuracy-neutral",
  x_fun=function(x){diff(x[c("v_neutral_FALSE","v_neutral_TRUE")]) - 
                    diff(x[c("v_speed_FALSE","v_speed_TRUE")])})
# Or from the sampled parameters (i.e., 1.12 - .22)
p_test(x=sPNAS_Bv_sv,x_name="v_Ea-s:lMd",y=sPNAS_Bv_sv,y_name="v_Ea-n:lMd",
       subfilter=2000,digits=3)


