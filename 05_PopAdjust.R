#1) #obey the original method
#there are no general answers 

# This file performs the indirect comparison methods on the simulated data

#rm(list=ls())

# setwd("C:/Users/Antonio/Desktop/population_adjustment_simstudy") 
load(file = "survival_settings.RData") # load simulation settings
source('02_functions.R') # load MAIC functions
#source('02_new.R')
set.seed(444)

#options(mc.cores = parallel::detectCores())
scenarios <- nrow(pc) #nrow(pc) # number of simulation scenarios
#replicates <- c(10, 29) # hashtag this when running FullFact (should be 1000)

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
IPD.AC.all2 <- IPD.BC.all2 <- ALD.BC.all2 <- vector(mode = "list", scenarios)

#scenarios
for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits = 2), 
                    "b_EM", round(pc$b_EM[i], digits = 2),"meanX_AC", pc$meanX_AC[i], "corX", pc$corX[i]) 
  load(paste0("Data_var/","2IPD_AC_", file.id, ".RData"))
  load(paste0("Data_var/","2IPD_BC_", file.id, ".RData")) 
  load(paste0("Data_var/","2ALD_BC_", file.id, ".RData"))  
  IPD.AC.all2[[i]] <- IPD.AC2
  ALD.BC.all2[[i]] <- ALD.BC2
  IPD.BC.all2[[i]] <- IPD.BC2
}


### Standard indirect treatment comparison (Bucher method)
bucher.wrapper <- function(data.AC, data.BC) {
  data.AC$trt <- factor(data.AC$trt, levels = c("C","A"))
  d.AC.bucher <- summary(coxph(Surv(time, status)~trt, data = data.AC))$coef[1]
  var.d.AC.bucher <- vcov(coxph(Surv(time, status)~trt, data = data.AC))[[1]]
  d.BC.bucher <- with(data.BC,logHR_B) 
  var.d.BC.bucher <- with(data.BC, var_logHR_B)
  d.AB.bucher <- unname(d.AC.bucher - d.BC.bucher) # average treatment effect
  var.d.AB.bucher <- unname(var.d.AC.bucher + var.d.BC.bucher) # variance of treatment effect
  list(d.AB.bucher, var.d.AB.bucher)
}  

### Matching-adjusted indirect comparison (MAIC)
maic.wrapper <- function(data.AC, data.BC, ems) { # ems indexes the position (columns) of the effect modifiers
  data.AC$trt <- factor(data.AC$trt, levels = c("C","A"))
  AC.ems <- data.AC[,1 + ems] # column 1 of IPD is treatment indicator 
  maic.weights <- maic(A.X = AC.ems, B.summary = data.BC[ems]) # maic weights through method of moments
  
  maic.weights.adj <- (maic.weights+e-10) *(1/max(maic.weights+e-10))
  
  data.AC.design <- survey::svydesign(~0, data = data.AC)
  outcome.fit.maic <- survey::svycoxph(Surv(time, status)~trt, design = data.AC.design, weights = maic.weights.adj, data = data.AC)
  
  maic.aess <- approx.ess(maic.weights) # approximate effective sample size
  # fit weighted Cox proportional hazards model using robust=TRUE
  
  #============================================================================
  #Here in the future a function could be added for bootstrapping SEs.
  
  
  # outcome.fit.maic <- coxph(Surv(time, status)~trt, robust = TRUE, weights = maic.weights, data = data.AC)
  d.AC.maic <- summary(outcome.fit.maic)$coef[1]
  
  var.d.AC.maic <- vcov(outcome.fit.maic)[[1]] # sandwich-like variance estimator for A vs. C
  d.BC.maic <- with(data.BC,logHR_B)
  var.d.BC.maic <- with(data.BC, var_logHR_B)
  d.AB.maic <- d.AC.maic - d.BC.maic # ATE for A vs. B
  var.d.AB.maic <- var.d.AC.maic + var.d.BC.maic # Variance of A vs. B treatment effect
  list(d.AB.maic, var.d.AB.maic, maic.aess)
}
### Simulated treatment comparison (STC) - original "plug-in" approach
stc.wrapper <- function(data.AC, data.BC, pvs, ems) {
   data.AC$trt <- factor(data.AC$trt, levels = c("C","A"))
   AC.chars <- data.AC[, 1 + pvs] # column 1 of IPD is treatment indicator
   pure.pvs <- setdiff(pvs, ems) # these are purely prognostic variables (are not effect modifiers)
   # these are not centered but the effect modifiers (both interaction and prognostic terms) are
   # fit outcome regresion model with IPD effect modifiers centered at the mean BC values
   stc.coxph <- coxph(as.formula(paste0("Surv(time,status)~", paste0(colnames(AC.chars)[pure.pvs],
                                                                      collapse = "+"),
                                         "+",  paste0("trt*I(", colnames(AC.chars)[ems],
                                                   "-", deparse(substitute(data.BC)),
                                                    "$mean.", colnames(AC.chars)[ems], ")",
                                                     collapse = "+"))), data = data.AC)
   d.AC.stc <- coef(stc.coxph)["trtA"]
   var.d.AC.stc <- vcov(stc.coxph)["trtA", "trtA"]
   d.BC.stc <- with(data.BC,logHR_B)
   var.d.BC.stc <- with(data.BC, var_logHR_B)
   d.AB.stc <- d.AC.stc - d.BC.stc # A vs. B treatment effect
   var.d.AB.stc <- var.d.AC.stc + var.d.BC.stc # A vs. B variance
   list(d.AB.stc, var.d.AB.stc)  

} 

#Added the code below to run the MLNMR
### ML-NMR using the multinma package
#saveAllWarnings <- function(expr) {
#  withCallingHandlers(expr, 
#                      warning=function(w) {
#                        warn <- conditionMessage(w)
#                        invokeRestart("muffleWarning")
#                      })
#}
#------------------------------------------------------------------------------

##mlnmr.wrapper <- function(data.AC, data.IPD.BC){
  ##data.AC$trtclass <- case_match(data.AC$trt,
                                ## "C" ~ "Placebo",
                               ##  "A" ~ "Active")
  ##data.AC$study <- "AC"
  
  #Used the IPD and thus assuming the Guyot method for recreating IPD is equally accurate.
 ##data.IPD.BC <- data.IPD.BC
 ## data.IPD.BC$trtclass <- case_match(data.IPD.BC$trt,
                             ##        "C" ~ "Placebo",
                             ##        "B" ~ "Active")
  
 ## data.IPD.BC$study <- "BC"
  #Create Covariate means per trt arm
 ## AgDCov <- data.frame()
 ## AgDCov <- data.IPD.BC %>% group_by(trt) %>% 
  ##  summarise(X1 = mean(X1), X2 = mean(X2), X3 = mean(X3), X4 = mean(X4)) %>% 
  ##  mutate(sdX1 = sd(X1), sdX2 = sd(X2), sdX3 = sd(X3), sdX4 = sd(X4))
  
 ## AgDCov$study <- "BC"
  
###  ndmm_net <- combine_network(
###   set_ipd(data.AC,
        ##    study = study,
        #    trt = trt,
         #   trt_class = trtclass,
         #   Surv = Surv(time / 12, status)),
  #  set_agd_surv(data.IPD.BC,
              #   study = study,
              #   trt = trt,
              #   trt_class = trtclass,
               #  Surv = Surv(time / 12, status),
              #   covariates = AgDCov)
 # )
  
 # ndmm_net <- add_integration(ndmm_net,
                         #     X1 = distr(qnorm, mean = X1, sd = sdX1),
                         #     X2 = distr(qnorm, mean = X2, sd = sdX2),
                         #     X3 = distr(qnorm, mean = X3, sd = sdX4),
                         #     X4 = distr(qnorm, mean = X4, sd = sdX4))
  
  #fit <- nma(ndmm_net,
            ## trt_effects = "fixed",
           #  regression = ~(X1 + X2 + X3 + X4)*.trt,
           # likelihood = "weibull",
           #  prior_intercept = normal(0, 100),
           #  prior_trt = normal(0, 100),
           #  prior_reg = normal(0, 100),
           #  prior_aux = half_normal(1),
           #  QR = TRUE, 
           #  iter = 2000)
  #warnings_list <- warnings()
  #64 integration-points
  #loghr <- relative_effects(fit, all_contrasts = TRUE)
  #rhat <- loghr$summary$Rhat
  #ness <- loghr$summary$n_eff
  #rhat_issues <- ifelse(any(rhat>1.05), TRUE, FALSE)
  #d.AB.MLNMR <- log(1/exp(loghr$summary$mean[3]))
 # var.d.AB.MLNMR <- log(1/exp(loghr$summary$sd[3]))^2
 
 # 
  #list(d.AB.MLNMR,var.d.AB.MLNMR, rhat_issues, warnings_list)
#}

#-------------------------------------------------------------------------------
#Removed anything relating to the setting up the cluster and added the MLNMR loop.

#--------------------------------------------------------------------------
#Replaced the cluster setup etc with doparallel etc.
# set up cluster for parallel computing based on the doParallel package
num.cores <- detectCores()
doParallel::registerDoParallel()

#cluster <- makeCluster(num.cores, type="SOCK", outfile="")
#registerDoSNOW(cluster)
# progress bar
pb <- txtProgressBar(max = replicates, style = 3)
 progress <- function(n) setTxtProgressBar(pb, n)
 opts <- list(progress = progress)

# combine lists in parallelisation
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# run indirect comparison methods for all replicates/scenarios in parallel
#scenarios
start_ <- Sys.time()
for (i in 1:scenarios) {
  IPD.AC2 <- IPD.AC.all2[[i]]
  ALD.BC2 <- ALD.BC.all2[[i]]
  IPD.BC2 <- IPD.BC.all2[[i]]
  
  pvs_i <- pvs
  ems_i <- ems
 
 #add this between simulation and "N_AC" 
 #replfrom,"-",replto
   
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits = 2), 
                    "b_EM", round(pc$b_EM[i], digits = 2),"meanX_AC", pc$meanX_AC[i], "corX", pc$corX[i]) 

  
   bucher.results <- foreach(j = 1:1000, .combine = 'comb', .multicombine = TRUE, .init = list(list(), list()),  .packages = c("dplyr","tidyr", "survival")) %dopar% {
                              results <- bucher.wrapper(IPD.AC2[[j]], ALD.BC2[[j]])
                              return(results)
                            }
  close(pb)
  means <- unlist(bucher.results[[1]])
  variances <- unlist(bucher.results[[2]])
  save(means, file = paste0("Results/Bucher/means_", file.id, ".RData"))
  save(variances, file = paste0("Results/Bucher/variances_", file.id, ".RData"))  
  maic.results <- foreach(j = 1:1000, .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list()),  .packages = c("dplyr","tidyr", "survey")) %dopar% {
                            results <- maic.wrapper(IPD.AC2[[j]], ALD.BC2[[j]], ems = ems_i)
                            return(results)
                          }
 close(pb)
  means <- unlist(maic.results[[1]])
  variances <- unlist(maic.results[[2]])
  approx.ess.maic <- unlist(maic.results[[3]])
  save(means, file = paste0("Results/MAIC/means_", file.id, ".RData"))
  save(variances, file = paste0("Results/MAIC/variances_", file.id, ".RData"))  
  save(approx.ess.maic, file = paste0("Results/MAIC/aess_", file.id, ".RData")) 

  
  stc.results <- foreach(j = 1:1000, .combine = 'comb', .multicombine = TRUE, .init = list(list(), list()),  .packages = c("survival")) %dopar% {
    results <- stc.wrapper(IPD.AC2[[j]], ALD.BC2[[j]], pvs = pvs_i, ems = ems_i)
    return(results) 
  }
   close(pb)
  means <- unlist(stc.results[[1]])
  variances <- unlist(stc.results[[2]])
  save(means, file = paste0("Results/STC/means_", file.id, ".RData"))
  save(variances, file = paste0("Results/STC/variances_", file.id, ".RData")) 
  
  # mlnmr.results <- foreach(j = replfrom:replto, .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(),list(),list()), .packages = c("dplyr", "tidyr", "multinma", "parallel")) %dopar% {
  #                          results <- mlnmr.wrapper(IPD.AC[[j]], IPD.BC[[j]])
  #                          return(results) 
  #                        }
  # close(pb)
  # means <- unlist(mlnmr.results[[1]])
  # variances <- unlist(mlnmr.results[[2]])
  # rhat <- unlist(mlnmr.results[[3]])
  # warn <- unlist(mlnmr.results[[4]])
  # #replicates[1]:replicates[2]
  # save(means, file = paste0("Results/MLNMR/means_", file.id, ".RData"))
  # save(variances, file = paste0("Results/MLNMR/variances_", file.id, ".RData"))
  # save(rhat, file = paste0("Results/MLNMR/rhat_", file.id, ".RData"))
  # save(warn, file = paste0("Results/MLNMR/warnings_", file.id, ".RData"))
  
}  

#stopCluster(cluster)
end_ <- Sys.time()
end_ - start_

