

#1) #obey the original method
#there are no general answers 

# This file performs the indirect comparison methods on the simulated data

#rm(list=ls())

# setwd("C:/Users/Antonio/Desktop/population_adjustment_simstudy") 
load(file = "survival_settings.RData") # load simulation settings
source('02_functions.R') # load MAIC functions
source('02_new.R')
set.seed(444)

#options(mc.cores = parallel::detectCores())
scenarios <- nrow(pc) #nrow(pc) # number of simulation scenarios
#replicates <- c(10, 29) # hashtag this when running FullFact (should be 1000)



# for data manipulation 
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)) {install.packages("tidyr"); library(tidyr)}
# for detect cores
if(!require(parallel)) {install.packages("parallel"); library(parallel)}
# for parallel cluster
if(!require(doSNOW)) {install.packages("doSNOW"); library(doSNOW)}
# to conduct standard and weighted cox regressions
if(!require(survival)) {install.packages("survival"); library(survival)}


#if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
#if (!require(tidyr)) {install.packages("tidyr"); library(tidyr)}
#if (!require(multinma)) {install.packages("multinma"); library(multinma)}
#if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
#if (!require(survival)) {install.packages("survival"); library(survival)}
#if (!require(survey)) {install.packages("survey"); library(survey)}
#if (!require(parallel)) {install.packages("parallel"); library(parallel)}
#if (!require(doParallel)) {install.packages("doParallel"); library(doParallel)}

#-------------------------------------------------------------------------------
# #All the code related to parallel clusters here and in the loop at the end where the functions are executed.
# for parallel cluster
#if(!require(doSNOW)) {install.packages("doSNOW"); library(doSNOW)}
# to conduct standard and weighted cox regressions


# load simulated patient-level (AC) and aggregate-level (BC) datasets for all scenarios
IPD.AC.all2 <- ALD.BC.all2 <- vector(mode = "list", scenarios)

#scenarios
for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits = 2), 
                    "b_EM", round(pc$b_EM[i], digits = 2),"meanX_AC", pc$meanX_AC[i], "corX", pc$corX[i]) 
  load(paste0("Data_var/","2IPD_AC_", file.id, ".RData"))
  load(paste0("Data_var/","2ALD_BC_", file.id, ".RData"))  
  IPD.AC.all2[[i]] <- IPD.AC2
  ALD.BC.all2[[i]] <- ALD.BC2

}




### Standard indirect treatment comparison (Bucher method)
bucher.wrapper <- function(data.AC, data.BC) {
  data.AC$trt <- factor(data.AC$trt, levels=c("C","A"))
  d.AC.bucher <- summary(coxph(Surv(time, status)~trt, data=data.AC))$coef[1]
  var.d.AC.bucher <- vcov(coxph(Surv(time, status)~trt, data=data.AC))[[1]]
  d.BC.bucher <- with(data.BC,logHR_B) 
  var.d.BC.bucher <- with(data.BC, var_logHR_B)
  d.AB.bucher <- unname(d.AC.bucher - d.BC.bucher) # average treatment effect
  var.d.AB.bucher <- unname(var.d.AC.bucher+var.d.BC.bucher) # variance of treatment effect
  list(d.AB.bucher, var.d.AB.bucher)
}  

### Matching-adjusted indirect comparison (MAIC)
maic.wrapper <- function(data.AC, data.BC, ems) { # ems indexes the position (columns) of the effect modifiers
  data.AC$trt <- factor(data.AC$trt, levels=c("C","A"))
  AC.ems <- data.AC[,1+ems] # column 1 of IPD is treatment indicator 
  maic.weights <- maic(A.X=AC.ems, B.summary=data.BC[ems]) # maic weights through method of moments
  maic.aess <- approx.ess(maic.weights) # approximate effective sample size
  # fit weighted Cox proportional hazards model using robust=TRUE
  tryCatch({outcome.fit.maic <- coxph(Surv(time, status)~trt, robust=TRUE, weights=maic.weights, data=data.AC)
  d.AC.maic <- summary(outcome.fit.maic)$coef[1]
  var.d.AC.maic <- vcov(outcome.fit.maic)[[1]] # sandwich-like variance estimator for A vs. C
  d.BC.maic <- with(data.BC,logHR_B)
  var.d.BC.maic <- with(data.BC, var_logHR_B)
  d.AB.maic <- d.AC.maic - d.BC.maic # ATE for A vs. B
  var.d.AB.maic <- var.d.AC.maic + var.d.BC.maic # Variance of A vs. B treatment effect
  list(d.AB.maic, var.d.AB.maic, maic.aess)}, error=function(e){message(paste(i))})
}

maic.wrapper2 <- function(data.AC, data.BC, ems) { # ems indexes the position (columns) of the effect modifiers
  data.AC$trt <- factor(data.AC$trt, levels = c("C","A"))
  AC.ems <- data.AC[,1 + ems] # column 1 of IPD is treatment indicator 
  maic.weights2 <- maic2(A.X = AC.ems, B.summary = data.BC[ems],sd.summary=data.BC[ems+4],data_IPD = data.AC,data_ALD = data.BC) # maic weights through method of moments
  maic.aess2 <- approx.ess(maic.weights2) # approximate effective sample size
  # fit weighted Cox proportional hazards model using robust=TRUE
  
  #============================================================================
  #Here in the future a function could be added for bootstrapping SEs.
  tryCatch({outcome.fit.maic <- coxph(Surv(time, status)~trt, robust=TRUE, weights=maic.weights2, data=data.AC)
  d.AC.maic <- summary(outcome.fit.maic)$coef[1]
  var.d.AC.maic <- vcov(outcome.fit.maic)[[1]] # sandwich-like variance estimator for A vs. C
  d.BC.maic <- with(data.BC,logHR_B)
  var.d.BC.maic <- with(data.BC, var_logHR_B)
  d.AB.maic <- d.AC.maic - d.BC.maic # ATE for A vs. B
  var.d.AB.maic <- var.d.AC.maic + var.d.BC.maic # Variance of A vs. B treatment effect
  list(d.AB.maic, var.d.AB.maic, maic.aess2)}, error=function(e){message(paste(i))})
  
  
}



### Simulated treatment comparison (STC) - original "plug-in" approach
stc.wrapper <- function(data.AC, data.BC, pvs, ems) {
  data.AC$trt <- factor(data.AC$trt, levels=c("C","A"))
  AC.chars <- data.AC[,1+pvs] # column 1 of IPD is treatment indicator
  pure.pvs <- setdiff(pvs, ems) # these are purely prognostic variables (are not effect modifiers)
  # these are not centered but the effect modifiers (both interaction and prognostic terms) are
  # fit outcome regresion model with IPD effect modifiers centered at the mean BC values
  stc.coxph <- coxph(as.formula(paste0("Surv(time,status)~", paste0(colnames(AC.chars)[pure.pvs],
                                                                    collapse="+"),
                                       "+",  paste0("trt*I(", colnames(AC.chars)[ems],
                                                    "-", deparse(substitute(data.BC)),
                                                    "$mean.", colnames(AC.chars)[ems], ")",
                                                    collapse="+"))), data=data.AC)
  d.AC.stc <- coef(stc.coxph)["trtA"]
  var.d.AC.stc <- vcov(stc.coxph)["trtA", "trtA"]
  d.BC.stc <- with(data.BC,logHR_B)
  var.d.BC.stc <- with(data.BC, var_logHR_B)
  d.AB.stc <- d.AC.stc - d.BC.stc # A vs. B treatment effect
  var.d.AB.stc <- var.d.AC.stc + var.d.BC.stc # A vs. B variance
  list(d.AB.stc, var.d.AB.stc)  
}  

# set up cluster for parallel computing
num.cores <- detectCores()
cluster <- makeCluster(num.cores, type="SOCK", outfile="")
registerDoSNOW(cluster)
# progress bar
pb <- txtProgressBar(max=replicates, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# combine lists in parallelisation
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


# run indirect comparison methods for all replicates/scenarios in parallel
for(i in 1:scenarios) {
  IPD.AC2 <- IPD.AC.all2[[i]]
  ALD.BC2 <- ALD.BC.all2[[i]]
 
  pvs_i <- pvs
  ems_i <- ems
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits=2), 
                    "b_EM", round(pc$b_EM[i], digits=2),"meanX_AC", pc$meanX_AC[i], "corX", pc$corX[i]) 
  bucher.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                            .init=list(list(), list()), .options.snow=opts,
                            .packages=c("dplyr","tidyr", "survival")) %dopar% {
                              results <- bucher.wrapper(IPD.AC2[[j]], ALD.BC2[[j]])
                              return(results)
                            }
  
  close(pb)
  means <- unlist(bucher.results[[1]])
  variances <- unlist(bucher.results[[2]])
  save(means, file=paste0("weight_results/Bucher/means_", file.id, ".RData"))
  save(variances, file=paste0("weight_results/Bucher/variances_", file.id, ".RData"))  
  maic.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                          .init=list(list(), list(), list()), .options.snow=opts,
                          .packages=c("dplyr","tidyr", "survey")) %dopar% {
                            results <- maic.wrapper(IPD.AC2[[j]], ALD.BC2[[j]], ems=ems_i)
                            return(results)
                          }
  close(pb)
  means <- unlist(maic.results[[1]])
  variances <- unlist(maic.results[[2]])
  approx.ess.maic <- unlist(maic.results[[3]])
  save(means, file=paste0("weight_results/MAIC/means_", file.id, ".RData"))
  save(variances, file=paste0("weight_results/MAIC/variances_", file.id, ".RData"))  
  save(approx.ess.maic, file=paste0("weight_results/MAIC/aess_", file.id, ".RData")) 
  stc.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                         .init=list(list(), list()), .options.snow=opts, 
                         .packages=c("survival")) %dopar% {
                         results <- stc.wrapper(IPD.AC2[[j]], ALD.BC2[[j]], pvs=pvs_i, ems=ems_i)
                            return(results) 
                         }
  close(pb)
  means <- unlist(stc.results[[1]])
  variances <- unlist(stc.results[[2]])
  save(means, file=paste0("weight_results/STC/means_", file.id, ".RData"))
  save(variances, file=paste0("weight_results/STC/variances_", file.id, ".RData"))  
  
  maic.results2 <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                          .init=list(list(), list(), list()), .options.snow=opts,
                          .packages=c("dplyr","tidyr", "survey")) %dopar% {
                            results2 <- maic.wrapper2(IPD.AC2[[j]], ALD.BC2[[j]], ems=ems_i)
                            return(results2)
                          }
  close(pb)
  means2 <- unlist(maic.results2[[1]])
  variances2 <- unlist(maic.results2[[2]])
  approx.ess.maic2 <- unlist(maic.results2[[3]])
  save(means2, file=paste0("weight_results/MAIC2/means_", file.id, ".RData"))
  save(variances2, file=paste0("weight_results/MAIC2/variances_", file.id, ".RData"))  
  save(approx.ess.maic2, file=paste0("weight_results/MAIC2/aess_", file.id, ".RData")) 
  
  
}  



stopCluster(cluster)


#singe result for MIAC2 version1
for(i in 1:scenarios) {
  IPD.AC2 <- IPD.AC.all2[[i]]
  ALD.BC2 <- ALD.BC.all2[[i]]
  
  pvs_i <- pvs
  ems_i <- ems
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits=2), 
                    "b_EM", round(pc$b_EM[i], digits=2),"meanX_AC", pc$meanX_AC[i], "corX", pc$corX[i]) 
  
  maic.results2 <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                           .init=list(list(), list(), list()), .options.snow=opts,
                           .packages=c("dplyr","tidyr", "survey")) %dopar% {
                             results2 <- maic.wrapper2(IPD.AC2[[j]], ALD.BC2[[j]], ems=ems_i)
                             return(results2)
                           }
  close(pb)
  means2 <- unlist(maic.results2[[1]])
  variances2 <- unlist(maic.results2[[2]])
  approx.ess.maic2 <- unlist(maic.results2[[3]])
  save(means2, file=paste0("weight_results/MAIC2/means_", file.id, ".RData"))
  save(variances2, file=paste0("weight_results/MAIC2/variances_", file.id, ".RData"))  
  save(approx.ess.maic2, file=paste0("weight_results/MAIC2/aess_", file.id, ".RData")) 
}

