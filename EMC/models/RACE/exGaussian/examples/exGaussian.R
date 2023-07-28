rm(list=ls())
source("emc/emc.R")
source("models/RACE/exGaussian/exGaussian.R")

# PNAS data
print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Often descriptive analyses with the ExGaussian are done only on correct 
# responses. 
dat <- dat[dat$S==dat$R,]

# Fit to individual design cells, this will be explained in later lessons, for
# now just note this enables us to fit an 18 parameter model, with separate mu,
# sigma and tau estimates for each design cell. 
se <- function(d) {factor(paste(d$S,d$E,sep="_"),levels=
  c("left_accuracy","right_accuracy","left_neutral","right_neutral","left_speed","right_speed"))}

designEXG <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(SE=diag(6)),
  Flist=list(mu~0+SE,sigma~0+SE,tau~0+SE),
  Ffunctions = list(SE=se),
  model=exGaussian)

# samplers <- make_samplers(dat,designEXG,type="standard",rt_resolution=.02)
# save(samplers,file="sPNAS_exg.RData")

print(load("models/RACE/exGaussian/examples/samples/sPNAS_exg.RData"))
check_run(samplers,layout=c(3,6),width=12,subfilter=500)


#### Check model fit ----

# post predict did not use first 500, so inference based on 2000*3 samples
plot_fit(dat,pp,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))
plot_fits(dat,pp,layout=c(2,3),lpos="right")

# Very accurate fits
tab <- plot_fit(dat,pp,layout=c(2,3),factors=c("E","S"),
         stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))
tab <- plot_fit(dat,pp,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.4),
         stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")
tab <- plot_fit(dat,pp,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
         stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")
tab <- plot_fit(dat,pp,layout=c(2,3),factors=c("E","S"),
         stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))


### Examine population mean (mu) parameters

# Priors all well dominated 
ci <- plot_density(samplers,layout=c(3,6),selection="mu",subfilter=500)
# Now on the natural scale
ci_mapped <- plot_density(samplers,layout=c(3,6),subfilter=500,
                                        selection="mu",mapped=TRUE)

# Looking at parameters on the natural scale
round(ci,2)
round(ci_mapped,2)

# skew: accuracy > neutral > speed
p_test(x=samplers,mapped=TRUE,x_name="left skew: accuracy - speed",
  x_fun=function(x){x["tau_left_accuracy"]-x["tau_left_speed"]})

p_test(x=samplers,mapped=TRUE,x_name="right skew: accuracy - speed",
  x_fun=function(x){x["tau_right_accuracy"]-x["tau_right_speed"]})

p_test(x=samplers,mapped=TRUE,x_name="left skew: accuracy - neutral",
  x_fun=function(x){x["tau_left_accuracy"]-x["tau_left_neutral"]})

p_test(x=samplers,mapped=TRUE,x_name="right skew: accuracy - neutral",
  x_fun=function(x){x["tau_right_accuracy"]-x["tau_right_neutral"]})

# Mean RT: accuracy > neutral > speed
p_test(x=samplers,mapped=TRUE,x_name="left mean: accuracy - speed",
  x_fun=function(x){sum(x[c("mu_left_accuracy","tau_left_accuracy")]) - 
                    sum(x[c("mu_left_speed","tau_left_speed")])})

p_test(x=samplers,mapped=TRUE,x_name="right mean: accuracy - speed",
  x_fun=function(x){sum(x[c("mu_right_accuracy","tau_right_accuracy")]) - 
                    mean(x[c("mu_right_speed","tau_right_speed")])})

p_test(x=samplers,mapped=TRUE,x_name="left mean: accuracy - neutral",
  x_fun=function(x){sum(x[c("mu_left_accuracy","tau_left_accuracy")]) - 
                    sum(x[c("mu_left_neutral","tau_left_neutral")])})

p_test(x=samplers,mapped=TRUE,x_name="right mean: accuracy - neutral",
  x_fun=function(x){sum(x[c("mu_right_accuracy","tau_right_accuracy")]) - 
                    sum(x[c("mu_right_neutral","tau_right_neutral")])})

# SD R: accuracy > neutral > speed
p_test(x=samplers,mapped=TRUE,x_name="left SD: accuracy - speed",
  x_fun=function(x){sqrt(sum(x[c("mu_left_accuracy","tau_left_accuracy")]^2)) - 
                    sqrt(sum(x[c("mu_left_speed","tau_left_speed")]^2))})

p_test(x=samplers,mapped=TRUE,x_name="right SD: accuracy - speed",
  x_fun=function(x){sqrt(sum(x[c("mu_right_accuracy","tau_right_accuracy")]^2)) - 
                    sqrt(sum(x[c("mu_right_speed","tau_right_speed")]^2))})

p_test(x=samplers,mapped=TRUE,x_name="left SD: accuracy - neutral",
  x_fun=function(x){sqrt(sum(x[c("mu_left_accuracy","tau_left_accuracy")]^2)) - 
                    sqrt(sum(x[c("mu_left_neutral","tau_left_neutral")]^2))})

p_test(x=samplers,mapped=TRUE,x_name="right SD: accuracy - neutral",
  x_fun=function(x){sqrt(sum(x[c("mu_right_accuracy","tau_right_accuracy")]^2)) - 
                    sqrt(sum(x[c("mu_right_neutral","tau_right_neutral")]^2))})



