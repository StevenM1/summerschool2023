rm(list=ls())
source("emc/emc.R")

#### DDM ----

print(load("models/DDM/DDM/examples/samples/is2_ddm_a.RData"))
print(load("models/DDM/DDM/examples/samples/is2_ddm_a_full.RData"))
print(load("models/DDM/DDM/examples/samples/is2_ddm_at0_full.RData"))
# following is a replicate to check how much they differ
print(load("models/DDM/DDM/examples/samples/is2_ddm_av_full_1.RData"))
is2_ddm_av_full_1 <- is2_ddm_av_full
print(load("models/DDM/DDM/examples/samples/is2_ddm_av_full.RData"))
print(load("models/DDM/DDM/examples/samples/is2_ddm_avt0_full.RData"))
# Same model as before but no cell coding
print(load("models/DDM/DDM/examples/samples/is2_ddm_avt0_full_nocell.RData"))


#               DIC  wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# avt0nocell -17039 0.126 -16806 0.885        233 -17271 -17449 -17504
# avt0       -17043 0.874 -16802 0.115        240 -17283 -17468 -17523
# at0        -16948 0.000 -16745 0.000        203 -17151 -17314 -17354
# av         -16410 0.000 -16210 0.000        199 -16609 -16777 -16808
# a          -16343 0.000 -16176 0.000        167 -16510 728446 -16677
compare_MLL(list(a=is2_ddm_a,afull=is2_ddm_a_full,at0full=is2_ddm_at0_full,
                 av=is2_ddm_av_full,avt0full=is2_ddm_avt0_full,
                 avt0fullnocell=is2_sPNAS_avt0_full_nocell))
      # avt0full        at0full avt0fullnocell             av              a          afull 
      #     0.65           0.24           0.10           0.01           0.00           0.00

# Cell coding vs. no cell coding now makes quite a difference. 
   
avrep <- list(av=is2_ddm_av_full,av=is2_ddm_av_full_1)
compare_MLL(avrep)    
#  av  av 
# 0.5 0.5 
lapply(avrep,median)
lapply(avrep,mean)
lapply(avrep,sd)
lapply(avrep,IQR)


#### LBA ----

print(load("models/RACE/LBA/examples/samples/is2_lba_B.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bvt0.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bt0_sv.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bv_sv.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bvt0_sv_NOa_n.RData"))

#           DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# B      -15243    0 -15100 0.000        143 -15386 -15485 -15529
# Bvt0   -15946    0 -15711 0.000        236 -16182 -16336 -16418
# Bt0sv  -16870    0 -16667 0.000        203 -17073 -17206 -17276
# Bvsv   -17200    1 -16986 0.999        214 -17414 -17569 -17627
# Bvt0sv -17182    0 -16973 0.001        209 -17391 -17542 -17601
compare_MLL(list(B=is2_lba_B,Bvt0=is2_lba_Bvt0,Bt0sv=is2_lba_Bt0_sv,Bvsv=is2_sPNAS_Bv_sv,Bvt0sv=is2_lba_Bvt0_sv_NOa_n))
# Bvt0sv   Bvsv  Bt0sv   Bvt0      B 
#   0.60   0.32   0.08   0.00   0.00 

#### RDM ----

print(load("models/RACE/RDM/examples/samples/is2_rdm_B.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bv.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bvt0.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bt0.RData"))

#         DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# B    -16514    0 -16333     0        180 -16694 -16844 -16874
# Bt0  -16511    0 -16309     0        202 -16714 -16863 -16916
# Bv   -16628    0 -16420     1        208 -16836 -16999 -17044
# Bvt0 -16648    1 -16403     0        245 -16893 -17054 -17137
compare_MLL(list(B=is2_rdm_B,Bv=is2_rdm_Bv,Bt0=is2_rdm_Bt0,Bvt0=is2_rdm_Bvt0))
#    B   Bv  Bt0 Bvt0 
# 0.82 0.08 0.08 0.01  

# !!! Big discrepancy !!!

#### LNR ----

print(load("models/RACE/LNR/examples/samples/is2_lnr_mu.RData"))
print(load("models/RACE/LNR/examples/samples/is2_lnr_mu_slM.RData"))
print(load("models/RACE/LNR/examples/samples/is2_lnr_mu_slM_t0.RData"))
print(load("models/RACE/LNR/examples/samples/is2_lnr_mu_sElM.RData"))
print(load("models/RACE/LNR/examples/samples/is2_mu_sElM_t0.RData"))

#            DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# mu      -15080    0 -14869     0        211 -15290 -15451 -15501
# mulM    -16840    0 -16633     0        207 -17047 -17217 -17254
# mulMt0  -17102    0 -16879     1        223 -17325 -17509 -17548
# muElM   -17053    0 -16809     0        245 -17298 -17501 -17542
# muElMt0 -17129    1 -16846     0        283 -17411 -17609 -17694
compare_MLL(list(mu=is2_lnr_mu,mulM=is2_lnr_mu_slM,mulMt0=is2_lnr_mu_slM_t0,
                 muElm=is2_lnr_mu_sElM,muElMt0=is2_lnr_mu_sElM_t0))
 # mulMt0    mulM   muElm muElMt0      mu 
 #   0.56    0.41    0.03    0.00    0.00 
    
#### DIC/BPIC best ----

print(load("is2/is2_ddm_avt0_full.RData")) 
print(load("is2/is2_ddm_avt0_full_nocell.RData"))
print(load("is2/is2_lba_Bv_sv.RData"))
print(load("is2/is2_mu_sElM_t0.RData"))
print(load("is2/is2_rdm_Bvt0.RData"))


mll <- list(ddm=is2_ddm_avt0_full,ddm_cell=is2_sPNAS_avt0_full_nocell,
            lba=is2_sPNAS_Bv_sv,rdm=is2_rdm_Bvt0,lnr=is2_lnr_mu_sElM_t0)

compare_MLL(mll)
#  ddm  lba  rdm  lnr 
# 0.50 0.45 0.04 0.02 
    
#### BF best

compare_MLL(list(ddm=is2_ddm_avt0_full,lba=is2_lba_Bvt0_sv_NOa_n,
                 rdm=is2_rdm_B,lnr=is2_lnr_mu_slM_t0))
#  ddm  lba  rdm  lnr 
# 0.50 0.45 0.04 0.02 

#### No cell DDM 

compare_MLL(list(ddm=is2_sPNAS_avt0_full_nocell,lba=is2_lba_Bvt0_sv_NOa_n,
                 rdm=is2_rdm_B,lnr=is2_lnr_mu_slM_t0))
#  lba  ddm  rdm  lnr 
# 0.71 0.22 0.05 0.02 

### Simulation study ----

model_p <- function(mll,nboot=100000) 
  # mll is a list of vectors of marginal log-likelihoods for a set of models
  # picks a vector of mlls for each model in the list randomly with replacement
  # nboot times, calculates model probabilities and averages, the default
  # nboot seems good for stable results at 2 decimal places.
{
  pmp <- function(x) 
    # posterior model probability for a vector of marginal log-likelihoods  
  {
    x <- exp(x-max(x))
    x/sum(x)
  }
  
  sort(apply(apply(do.call(rbind,lapply(mll,function(x){
    x[sample(length(x),nboot,replace=TRUE)]})),2,pmp),1,mean),decreasing=TRUE)
}

round(model_p(mll[-2]),2)
round(model_p(mll[-1]),2)




round(model_p(mll[1:2]),2)
round(model_p(mll[c(1,3)]),2)
round(model_p(mll[c(2,3)]),2)
round(model_p(mll[c(2,4)]),2)

sizeCheck <- function(mll1,mll2) {
  r40 <- numeric(24)
  for (i in 1:length(r40)) {
    indx <- (1+(i-1)*40):(i*40)
    r40[i] <- model_p(list(mll1[[1]][indx],mll2[[1]][indx]))[1]
  }
  print(round(c(range(r40),sd(r40)),2)) 
  r80 <- numeric(12)
  for (i in 1:length(r80)) {
    indx <- (1+(i-1)*80):(i*80)
    r80[i] <- model_p(list(mll1[[1]][indx],mll2[[1]][indx]))[1]
  }
  print(round(c(range(r80),sd(r80)),2)) 
  r120 <- numeric(8)
  for (i in 1:length(r120)) {
    indx <- (1+(i-1)*120):(i*120)
    r120[i] <- model_p(list(mll1[[1]][indx],mll2[[1]][indx]))[1]
  }
  print(round(c(range(r120),sd(r120)),2)) 
  r160 <- numeric(6)
  for (i in 1:length(r160)) {
    indx <- (1+(i-1)*160):(i*160)
    r160[i] <- model_p(list(mll1[[1]][indx],mll2[[1]][indx]))[1]
  }
  print(round(c(range(r160),sd(r160)),2)) 
}

sizeCheck(mll[1],mll[2])
sizeCheck(mll[1],mll[3])
sizeCheck(mll[2],mll[3])
sizeCheck(mll[2],mll[4])






