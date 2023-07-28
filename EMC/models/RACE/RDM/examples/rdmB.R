rm(list=ls())
source("emc/emc.R")
source("models/RACE/RDM/rdmB.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Average and difference matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))

#### Models ----

design_B <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~1),
  model=rdmB)

rdm_B <- make_samplers(dat,design_B,type="standard",rt_resolution=.02)
save(rdm_B,file="rdmPNAS_B.RData")


design_Bvt0<- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~E),
  model=rdmB)

rdm_Bvt0<- make_samplers(dat,design_Bvt0,type="standard",rt_resolution=.02)
save(rdm_Bvt0,file="rdmPNAS_Bvt0.RData")

design_Bt0<- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~E),
  model=rdmB)

rdm_Bt0 <- make_samplers(dat,design_Bt0,type="standard",rt_resolution=.02)
save(rdm_Bt0,file="rdmPNAS_Bt0.RData")

design_Bv<- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~1),
  model=rdmB)

rdm_Bv <- make_samplers(dat,design_Bv,type="standard",rt_resolution=.02)
save(rdm_Bv,file="rdmPNAS_Bv.RData")

# All of these models are run in  a single file, run_rdm.R using the run_emc
# function. This function is given the name of the file containing the samplers
# object and any arguments to its constituent functions (auto_burn, auto_adapt
# and auto_sample). For example the first model was fit with 8 cores per chain
# and so 24 cores in total given the default of 3 chains. 
#
# run_emc("rdmPNAS_B.RData",cores_per_chain=8)
#
# If any stage failsrun_check stops, saving the progress so far. The only 
# argument unique to this function is nsample, which determines the number of 
# iterations performed if the sample stage is reached (by default 1000).
#
# After the convergence checking reported below extra samples were obtained for
# some models that converged slowly, so all had 1000 converged samples. This was
# done by adding another run_emc call with nsamples set appropriately (run_emc
# recognizes that burn and adapt stages are complete and so just adds on the
# extra samples at the end). 
#
# Finally for the selected model extra samples were added to obtain 5000 
# iterations, so posterior inference was reliable, and posterior predictives 
# were obtained. (NB: if the file is run several times some lines must be 
# commented out or unwanted sample stage iteration will be added)

#### Check convergence ----

print(load("models/RACE/RDM/examples/samples/rdmPNAS_B.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bt0.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bv.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bvt0.RData"))

check_run(rdm_B)
check_run(rdm_Bt0,subfilter=500)
check_run(rdm_Bv,subfilter=1500)
check_run(rdm_Bvt0,subfilter=1500)

#### Model selection ----

# Bvt0 wins with Bt0 second
compare_IC(list(B=rdm_B,Bt0=rdm_Bt0,Bv=rdm_Bv,Bvt0=rdm_Bvt0),
           subfilter=list(0,500,1500,1500))
# But at the subject level Bv is second
ICs <- compare_ICs(list(B=rdm_B,Bt0=rdm_Bt0,Bv=rdm_Bv,Bvt0=rdm_Bvt0),
                   subfilter=list(0,500,1500,1500:2500))
table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$DIC)]})))
table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$BPIC)]})))

# Comparing with the best DDM (16 parameter) model to the best (15 parameter)
# LBA model and best (16 parameter) RDM the latter comes third. 
source("models/DDM/DDM/ddmTZD.R")
print(load("models/DDM/DDM/examples/samples/sPNAS_avt0_full.RData")) 
source("models/RACE/LBA/lbaB.R")
print(load("models/RACE/LBA/examples/samples/sPNAS_Bv_sv.RData"))
compare_IC(list(DDM_avt0=sPNAS_avt0_full,LBA_Bvsv=sPNAS_Bv_sv,RDM_Bvt0=rdm_Bvt0),
           subfilter=list(500,2000,1500))

# For the remaining analysis add 4000 iterations to Bvt0 model

####  Fit of winning (Bvsv) model ----

# post predict did not use first 1500, so inference based on 5000*3 samples
plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))
plot_fits(dat,pprdm_Bvt0,layout=c(2,3),lpos="right")

# Good fit, slight under-estimation of 10th percentile and over-estimation of 
# error RT in speed
pc <- function(d) 100*mean(d$S==d$R)
plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),
         stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),
         stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.4),
         stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
         stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),
         stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
         stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")

#### Posterior parameter inference ----

### Population mean (mu) tests

# Priors all well dominated except some rate parameters
cirdm_Bvt0 <- plot_density(rdm_Bvt0,layout=c(3,6),selection="mu",subfilter=1500)
# On the natural scale it is evident this is because the mismatch (FALSE) rates
# are least well updated, due to the fairly low error rates (errors give the
# most information about FALSE rates).
cirdm_Bvt0_mapped <- plot_density(rdm_Bvt0,layout=c(3,6),
                                        selection="mu",mapped=TRUE)
# Looking at parameters both with and without mapping
round(cirdm_Bvt0,2)
round(cirdm_Bvt0_mapped,2)

# For simpler estimates:
# 1) t0 is longer than the LBA
# 2) Start point noise slightly larger relative to B than for the LBA.

### B effects

# Recall the map used with 
get_map(rdm_Bvt0)$B

# B_lRd tests threshold right - threshold left, as for the LBA, although not 
# quite credible it indicates slightly higher right thresholds (i.e., a bias to 
# respond left)
p_test(x=rdm_Bvt0,x_name="B_lRd",subfilter=1500,digits=3)

# B_Ea-n and B_Ea-s measure differences in response caution (i.e., thresholds
# averaged over left and right accumulators), accuracy-neutral and accuracy-speed 
# respectively. Caution for accuracy is clearly higher than speed, but not 
# credibly greater than neutral.
p_test(x=rdm_Bvt0,x_name="B_Ea-n",subfilter=1500)
p_test(x=rdm_Bvt0,x_name="B_Ea-s",subfilter=1500)

# Here we construct a test on the natural scale showing caution is greater for 
# neutral than speed
p_test(x=rdm_Bvt0,mapped=TRUE,x_name="average B: neutral-speed",
  x_fun=function(x){mean(x[c("B_left_neutral","B_right_neutral")]) - 
                    mean(x[c("B_left_speed","B_right_speed")])})

# The remaining terms test interactions with bias (i.e., lR), with evidence of
# a small but credibly stronger bias to respond left (i.e., a lower threshold
# for the left accumulator) for speed than accuracy.
p_test(x=rdm_Bvt0,x_name="B_lRd:Ea-s",subfilter=1500,digits=2)

### v effects

# Again recall the map used with 
get_map(rdm_Bvt0)$v

# v_Ea-n v_Ea-s indicate that processing rate (the average of matching and 
# mismatching rates) is less in the accuracy condition than in neutral or speed.
p_test(x=rdm_Bvt0,x_name="v_Ea-n",subfilter=1500)
p_test(x=rdm_Bvt0,x_name="v_Ea-s",subfilter=1500)

# However, neutral and speed do not credibly differ
p_test(x=rdm_Bvt0,mapped=TRUE,x_name="average v: speed-neutral",
  x_fun=function(x){mean(x[c("v_TRUE_speed","v_FALSE_speed")]) - 
                    mean(x[c("v_TRUE_neutral","v_FALSE_neutral")])})

# v_lMd tests the quality of selective attention, rate match - rate mismatch
p_test(x=rdm_Bvt0,x_name="v_lMd",subfilter=1500)

# v_Ea-n indicates quality does not differ credibly between accuracy and neutral.
p_test(x=rdm_Bvt0,x_name="v_lMd:Ea-n",subfilter=1500)

# In contrast there is a strong difference for accuracy - speed
p_test(x=rdm_Bvt0,x_name="v_lMd:Ea-s",subfilter=1500)

# The neutral - speed difference is also highly credible, here is is tested
# in the mapped form
p_test(x=rdm_Bvt0,mapped=TRUE,x_name="quality: accuracy-neutral",
  x_fun=function(x){diff(x[c("v_FALSE_neutral","v_TRUE_neutral")]) - 
                    diff(x[c("v_FALSE_speed","v_TRUE_speed")])})



