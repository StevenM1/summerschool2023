rm(list=ls())
source("emc/emc.R")
source("models/RACE/RDM/rdmB.R")
run_emc("rdm_B_MT.RData",cores_per_chain=10)
# # Add samples to get 1000 converged.
# run_emc("rdm_B_MT.RData",cores_per_chain=10,nsample=1000)
# # Add samples to selected model
# run_emc("rdm_B_MT.RData",cores_per_chain=10,nsample=4000)
# # Get posterior predictives
# load("rdm_B_MT.RData")
# pprdm_Bvt0 <- post_predict(rdm_Bvt0,n_cores=19,subfilter=1500)
# save(rdm_B_MT,pprdm_Bvt0,file="rdm_B_MT.RData")

