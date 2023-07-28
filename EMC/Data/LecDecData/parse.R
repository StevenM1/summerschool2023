rm(list=ls())
expt1 <- read.delim("SpeedAccData.txt",header=F)
expt2 <- read.delim("PropData.txt",header=F)
names(expt1) <- c("s","blk","practice","E","stim","F","R","RT","censor")
round(prop.table(table(expt1$censor)),3)
# 0     1 
# 0.979 0.021 
names(expt2) <- c("s","blk","practice","P","stim","F","R","RT","censor")
round(prop.table(table(expt2$censor)),3)
#     0     1 
# 0.983 0.017
expt1 <- expt1[expt1$censor!=1,-c(3,9)]
expt2 <- expt2[expt2$censor!=1,-c(3,9)]
expt1$s=factor(expt1$s)
expt2$s=factor(expt2$s)
expt1$E=factor(expt1$E,labels=c("accuracy","speed"))
expt2$P=factor(expt2$P,labels=c("25","75"))
expt1$S <- expt1$F<4
expt1$S <- factor(expt1$S,labels=c("nw","w"))
expt2$S <- expt2$F<4
expt2$S <- factor(expt2$S,labels=c("nw","w"))
expt1$R <- factor(expt1$R,labels=c("NW","W"))
expt2$R <- factor(expt2$R,labels=c("NW","W"))
expt1$F <- ((expt1$F-1) %% 3)+1
expt2$F <- ((expt2$F-1) %% 3)+1
expt1$F <- factor(expt1$F,labels=c("hf","lf","vlf"))
expt2$F <- factor(expt2$F,labels=c("hf","lf","vlf"))
expt1$C <- toupper(expt1$S)==toupper(expt1$R)
expt2$C <- toupper(expt2$S)==toupper(expt2$R)
# add in combined NW and W factor (4 levels, lump NW together)
expt1$W <- as.character(expt1$F)
expt2$W <- as.character(expt2$F)
expt1$W[expt1$S=="nw"] <- "nw"
expt2$W[expt2$S=="nw"] <- "nw"
expt1$W <- factor(expt1$W,c("hf","lf","vlf","nw"))
expt2$W <- factor(expt2$W,c("hf","lf","vlf","nw"))
save(expt1,expt2,file="ejldt-dat.RData")