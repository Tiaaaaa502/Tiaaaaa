## This file processes the weight_results of the simulation study and computes
## and graphs the relevant performance metrics. 
replicates<-1000
library(dplyr)
library(tidyr)
library(ggplot2)

#load(file="survival_settings.RData")
source("02_Functions.R") # load functions to compute performance measures and for plotting
source("02_new.R")
# TRUE VALUE OF TARGET ESTIMAND IS ZERO
#-------------------------------------------------------------------------------
#ADJUSTED THE DELTA.AB from ZERO to the below
Delta.AB <- log(0.25)-log(0.25)

#Q: Assume this should be on the log scale -------------------------------------data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAABDklEQVR4Xp2UMaoCMRRFFcHCxsrO0tba/i/ADfwV/A2ICFYW7kAY/CSZOGYgnVjbW9taugFbG8f3OpO5A3lzIFh4Dg4ON53VelPp3NWOsW7R+SLVq1FVVbchfhtz/JV6EO99vyF+Ket+pB4ky/ywIX5qXUylHmRfFGMUm9w9+DupB/k/HCY4Lu9K+ZHUgyhbzlBM55pl54HUgxhTzlGscneil92TehBljn9xyMfYctfGq0HSMo748BO38QJ4PCS9UcSjk3oBPBqSXijisUm9AB6L5tGAiEcm9QJ4JDwWFMVDS/EC+FdJuqEovipSvAD+v0i6oCi+7FK8AH7T2roCRfF1neLVIGmLIvpctvFiPtxv6UGYMompAAAAAElFTkSuQmCC

# container of means and variances for all replicates in each scenario  
maic.means.list <- vector(mode="list", scenarios)
maic.variances.list <- vector(mode="list", scenarios)
maic.means.list2 <- vector(mode="list", scenarios)
maic.variances.list2 <- vector(mode="list", scenarios)

stc.means.list <- vector(mode="list", scenarios)
stc.variances.list <- vector(mode="list", scenarios)
bucher.means.list <- vector(mode="list", scenarios)
bucher.variances.list <- vector(mode="list", scenarios)
#mlnmr.means.list <- vector(mode="list", scenarios)
#mlnmr.variances.list <- vector(mode="list", scenarios)

# average treatment effect container 
maic.ate <- rep(NA, scenarios) 
maic.ate2 <- rep(NA, scenarios) 
stc.ate <- rep(NA, scenarios)
bucher.ate <- rep(NA, scenarios)
#mlnmr.ate <- rep(NA, scenarios)

maic.ate.mcse <- rep(NA, scenarios) # monte carlo standard error
maic.ate.mcse2 <- rep(NA, scenarios)
stc.ate.mcse <- rep(NA, scenarios)
bucher.ate.mcse <- rep(NA, scenarios)
#mlnmr.ate.mcse <- rep(NA, scenarios)

# lower bound of confidence interval container
maic.lci <- vector(mode="list", scenarios)
maic.lci2 <- vector(mode="list", scenarios)
stc.lci <- vector(mode="list", scenarios)
bucher.lci <- vector(mode="list", scenarios)
#mlnmr.lci <- vector(mode="list", scenarios)

maic.lci.mean <- rep(NA, scenarios)
maic.lci.mean2 <- rep(NA, scenarios)
stc.lci.mean <- rep(NA, scenarios)
bucher.lci.mean <- rep(NA, scenarios)
#mlnmr.lci.mean <- rep(NA, scenarios)

maic.lci.mcse <- rep(NA, scenarios)
maic.lci.mcse2 <- rep(NA, scenarios)
stc.lci.mcse <- rep(NA, scenarios)
bucher.lci.mcse <- rep(NA, scenarios)
#mlnmr.lci.mcse <- rep(NA, scenarios)

# upper bound of confidence interval container
maic.uci <- vector(mode="list", scenarios)
maic.uci2 <- vector(mode="list", scenarios)
stc.uci <- vector(mode="list", scenarios)
bucher.uci <- vector(mode="list", scenarios)
#mlnmr.uci <- vector(mode="list", scenarios)

maic.uci.mean <- rep(NA, scenarios)
maic.uci.mean2 <- rep(NA, scenarios)
stc.uci.mean <- rep(NA, scenarios)
bucher.uci.mean <- rep(NA, scenarios)
#mlnmr.uci.mean <- rep(NA, scenarios)

maic.uci.mcse <- rep(NA, scenarios)
maic.uci.mcse2 <- rep(NA, scenarios)
stc.uci.mcse <- rep(NA, scenarios)
bucher.uci.mcse <- rep(NA, scenarios)
#mlnmr.uci.mcse <- rep(NA, scenarios)

# Confidence interval widths
maic.ciwidth <- vector(mode="list", scenarios) 
maic.ciwidth2 <- vector(mode="list", scenarios)
stc.ciwidth <- vector(mode="list", scenarios)
bucher.ciwidth <- vector(mode="list", scenarios)
#mlnmr.ciwidth <- vector(mode="list", scenarios)

# Bias containers
maic.bias <- rep(NA, scenarios) 
maic.bias2 <- rep(NA, scenarios) 
stc.bias <- rep(NA, scenarios)
bucher.bias <- rep(NA, scenarios)
#mlnmr.bias <- rep(NA, scenarios)

maic.bias.mcse <- rep(NA, scenarios) 
maic.bias.mcse2 <- rep(NA, scenarios) 
stc.bias.mcse <- rep(NA, scenarios)
bucher.bias.mcse <- rep(NA, scenarios)
#mlnmr.bias.mcse <- rep(NA, scenarios)

# mean absolute error (MAE) containers
maic.abs.err <- vector(mode="list", scenarios)
maic.abs.err2 <- vector(mode="list", scenarios)
stc.abs.err <- vector(mode="list", scenarios)
bucher.abs.err <- vector(mode="list", scenarios)
#mlnmr.abs.err <- vector(mode="list", scenarios)

maic.mae <- rep(NA, scenarios)
maic.mae2 <- rep(NA, scenarios)
stc.mae <- rep(NA, scenarios)
bucher.mae <- rep(NA, scenarios)
#mlnmr.mae <- rep(NA, scenarios)

maic.mae.mcse <- rep(NA, scenarios)
maic.mae.mcse2 <- rep(NA, scenarios)
stc.mae.mcse <- rep(NA, scenarios)
bucher.mae.mcse <- rep(NA, scenarios)
#mlnmr.mae.mcse <- rep(NA, scenarios)

# mean square error (MSE) containers
maic.mse <- rep(NA, scenarios)
maic.mse2 <- rep(NA, scenarios)
stc.mse <- rep(NA, scenarios)
bucher.mse <- rep(NA, scenarios)
#mlnmr.mse <- rep(NA, scenarios)

maic.mse.mcse <- rep(NA, scenarios)
maic.mse.mcse2 <- rep(NA, scenarios)
stc.mse.mcse <- rep(NA, scenarios)
bucher.mse.mcse <- rep(NA, scenarios)
#mlnmr.mse.mcse <- rep(NA, scenarios)

# variability ratio containers
maic.vr <- rep(NA, scenarios)
maic.vr2 <- rep(NA, scenarios)
stc.vr <- rep(NA, scenarios)
bucher.vr <- rep(NA, scenarios)
#mlnmr.vr <- rep(NA, scenarios)

maic.vr.mcse <- rep(NA, scenarios)
maic.vr.mcse2 <- rep(NA, scenarios)
stc.vr.mcse <- rep(NA, scenarios)
bucher.vr.mcse <- rep(NA, scenarios)
#mlnmr.vr.mcse <- rep(NA, scenarios)

# empirical standard error (EmpSE) containers
maic.empse <- rep(NA, scenarios)
maic.empse2 <- rep(NA, scenarios)
stc.empse <- rep(NA, scenarios)
bucher.empse <- rep(NA, scenarios)
#mlnmr.empse <- rep(NA, scenarios)

maic.empse.mcse <- rep(NA, scenarios)
maic.empse.mcse2 <- rep(NA, scenarios)
stc.empse.mcse <- rep(NA, scenarios)
bucher.empse.mcse <- rep(NA, scenarios)
#mlnmr.empse.mcse <- rep(NA, scenarios)

# coverage rate (%) of 95% confidence intervals
maic.cov <- rep(NA, scenarios)
maic.cov2 <- rep(NA, scenarios)
stc.cov <- rep(NA, scenarios)
bucher.cov <- rep(NA, scenarios)
#mlnmr.cov <- rep(NA, scenarios)

maic.cov.mcse <- rep(NA, scenarios)
maic.cov.mcse2 <- rep(NA, scenarios)
stc.cov.mcse <- rep(NA, scenarios)
bucher.cov.mcse <- rep(NA, scenarios)
#mlnmr.cov.mcse <- rep(NA, scenarios)

# % of replicates worse than Bucher for each scenario  
maic.error.worse.than <- rep(NA,scenarios)
maic.error.worse.than.mcse <- rep(NA,scenarios)
maic.error.worse.than2 <- rep(NA,scenarios)
maic.error.worse.than.mcse2 <- rep(NA,scenarios)
stc.error.worse.than <- rep(NA,scenarios)
stc.error.worse.than.mcse <- rep(NA,scenarios)
##mlnmr.error.worse.than <- rep(NA,scenarios)
#mlnmr.error.worse.than.mcse <- rep(NA,scenarios)

# standardized biases
maic.std.bias <- rep(NA,scenarios) 
maic.std.bias2 <- rep(NA,scenarios) 
stc.std.bias <- rep(NA,scenarios)
bucher.std.bias <- rep(NA,scenarios)
#mlnmr.std.bias <- rep(NA,scenarios)

# approximate effective sample sizes
maic.aess.mean <- rep(NA,scenarios)
maic.aess.mcse <- rep(NA,scenarios)
maic.aess.mean2 <- rep(NA,scenarios)
maic.aess.mcse2 <- rep(NA,scenarios)

# number of replicates for which MAIC weighted regression has separability issues 
maic.separability.list <- rep(NA, scenarios)
# assume MAIC analysis stage has separability issues if abs(treatment effect) > max.effect
max.effect <- 1 # 

# number of replicates for which MAIC weighted regression has separability issues 
maic.separability.list2 <- rep(NA, scenarios)

# number of replicates for which STC does not converge 
stc.notconverges.list <- rep(NA, scenarios)
# assume that STC has not converged if variance greater than this number
max.variance <- 5

# number of replicates for which MLNMR does not converge 
#mlnmr.notconverges.list <- rep(NA, scenarios)
# assume that MLNMR has not converged if variance greater than this number
#max.variance.mlnmr <- 5  

#Individual Bias
#mlnmr_bias_nonsum <-rep(NA,scenarios)
stc_bias_nonsum <-rep(NA,scenarios) 
maic_bias_nonsum <- rep(NA,scenarios)
maic_bias_nonsum2 <- rep(NA,scenarios)
bucher_bias_nonsum <-rep(NA,scenarios)

# table storing parameter settings and performance measures for each scenario
scenarios.df <- data.frame()



#ADD REPLICATES SNIPPET TO THIS

#ADD replfrom,"-",replto,

for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits=2), 
                    "b_EM", round(pc$b_EM[i], digits=2),
                    "meanX_AC", pc$meanX_AC[i], "corX", pc$corX[i])   
  # MAIC
  load(paste0("weight_results/MAIC/means_", file.id, ".RData"))
  load(paste0("weight_results/MAIC/variances_", file.id, ".RData"))
  # number of replicates for which MAIC weighted regression has separability issues 
  maic.separability.list[i] <- sum(abs(means[!is.na(means)])>max.effect)
  # discard these replicates (1 of 162000 in this sim. study)
  maic.no.wei.pro <- !is.na(means) & !is.na(variances) #& variances<5 #& variances<1000
  means <- means[maic.no.wei.pro]
  variances <- variances[maic.no.wei.pro]
  maic.no.sep <- abs(means)<max.effect
  #means[!maic.no.sep] <- mean(means[abs(means)<max.effect & !is.na(means)])
  #variances[!maic.no.sep] <- mean(variances[variances<100])
  means <- means[maic.no.sep]
  variances <- variances[maic.no.sep]
  #means <- means[variances<5]
  
  #means <- means[!is.na(means)]
  

  
  
  ### MATCHING-ADJUSTED INDIRECT COMPARISON
  maic.means.list[[i]] <- means
  maic.variances.list[[i]] <- variances
  maic.bias[i] <- bias(maic.means.list[[i]], Delta.AB)  
  maic.bias.mcse[i] <- bias.mcse(maic.means.list[[i]])
  maic.mae[i] <- mae(maic.means.list[[i]], Delta.AB)
  maic.abs.err[[i]] <- maic.means.list[[i]] - (Delta.AB)
  maic.mae.mcse[i] <- mcse.estimate(maic.abs.err[[i]])
  maic.mse[i] <- mse(maic.means.list[[i]], Delta.AB) 
  maic.mse.mcse[i] <- mse.mcse(maic.means.list[[i]], Delta.AB) 
  maic.ate[i] <- mean(maic.means.list[[i]])
  maic.ate.mcse[i] <- mcse.estimate(maic.means.list[[i]])
  # construct confidence interval using normal distribution
  maic.lci[[i]] <- maic.means.list[[i]] + qnorm(0.025)*sqrt(maic.variances.list[[i]])
  maic.uci[[i]] <- maic.means.list[[i]] + qnorm(0.975)*sqrt(maic.variances.list[[i]])
  maic.ciwidth[[i]] <- maic.uci[[i]] - maic.lci[[i]]
  maic.lci.mean[i] <- mean(maic.lci[[i]])
  maic.lci.mcse[i] <- mcse.estimate(maic.lci[[i]])
  maic.uci.mean[i] <- mean(maic.uci[[i]])
  maic.uci.mcse[i] <- mcse.estimate(maic.uci[[i]])
  maic.cov[i] <- coverage(maic.lci[[i]], maic.uci[[i]], Delta.AB)
  maic.cov.mcse[i] <- coverage.mcse(maic.cov[i], length(maic.lci[[i]]))
  maic.empse[i] <- empse(maic.means.list[[i]])
  maic.empse.mcse[i] <- empse.mcse(maic.empse[i], length(maic.means.list[[i]]))
  maic.vr[i] <- var.ratio(maic.means.list[[i]], sqrt(maic.variances.list[[i]]))
  maic.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(maic.variances.list[[i]])), 
                                    emp.se=maic.empse[i],
                                    var.avg.se=mcse.estimate(sqrt(maic.variances.list[[i]]))^2,
                                    var.emp.se=maic.empse.mcse[i]^2)
  maic.std.bias[i] <- (maic.bias[i]*100)/maic.empse[i]
  load(paste0("weight_results/MAIC/aess_", file.id, ".RData"))
  maic.aess.mean[i] <- mean(approx.ess.maic)
  maic.aess.mcse[i] <- mcse.estimate(approx.ess.maic) 
  
  
  ### SIMULATED TREATMENT COMPARISON ("plug-in" approach)
  load(paste0("weight_results/STC/means_", file.id, ".RData"))
  load(paste0("weight_results/STC/variances_", file.id, ".RData"))  
  stc.notconverges.list[i] <- sum(variances>max.variance)
  # discard replicates for which STC did not converge (0 of 162000 in this sim. study)
  stc.converges <- variances<max.variance
  means <- means[stc.converges]
  variances <- variances[stc.converges]
  stc.means.list[[i]] <- means
  stc.variances.list[[i]] <- variances
  stc.bias[i] <- bias(stc.means.list[[i]], Delta.AB)  
  stc.bias.mcse[i] <- bias.mcse(stc.means.list[[i]])
  stc.mae[i] <- mae(stc.means.list[[i]], Delta.AB)
  stc.mae.mcse[i] <- mcse.estimate(stc.means.list[[i]])
  stc.abs.err[[i]] <- stc.means.list[[i]] - (Delta.AB)
  stc.mae.mcse[i] <- mcse.estimate(stc.abs.err[[i]])
  stc.mse[i] <- mse(stc.means.list[[i]], Delta.AB) 
  stc.mse.mcse[i] <- mse.mcse(stc.means.list[[i]], Delta.AB) 
  stc.ate[i] <- mean(stc.means.list[[i]]) 
  stc.ate.mcse[i] <- mcse.estimate(stc.means.list[[i]])
  stc.lci[[i]] <- stc.means.list[[i]] + qnorm(0.025)*sqrt(stc.variances.list[[i]])
  stc.uci[[i]] <- stc.means.list[[i]] + qnorm(0.975)*sqrt(stc.variances.list[[i]])
  stc.ciwidth[[i]] <- stc.uci[[i]] - stc.lci[[i]]
  stc.lci.mean[i] <- mean(stc.lci[[i]])
  stc.lci.mcse[i] <- mcse.estimate(stc.lci[[i]])
  stc.uci.mean[i] <- mean(stc.uci[[i]])
  stc.uci.mcse[i] <- mcse.estimate(stc.uci[[i]])
  stc.cov[i] <- coverage(stc.lci[[i]], stc.uci[[i]], Delta.AB)
  stc.cov.mcse[i] <- coverage.mcse(stc.cov[i], length(stc.lci[[i]]))
  stc.empse[i] <- empse(stc.means.list[[i]])
  stc.empse.mcse[i] <- empse.mcse(stc.empse[i], length(stc.means.list[[i]]))
  stc.vr[i] <- var.ratio(stc.means.list[[i]], sqrt(stc.variances.list[[i]]))
  stc.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(stc.variances.list[[i]])), 
                                   emp.se=stc.empse[i],
                                   var.avg.se=mcse.estimate(sqrt(stc.variances.list[[i]]))^2,
                                   var.emp.se=stc.empse.mcse[i]^2)  
  stc.std.bias[i] <- (stc.bias[i]*100)/stc.empse[i]
  
  
  ### BUCHER METHOD (STANDARD INDIRECT COMPARISON)
  load(paste0("weight_results/Bucher/means_", file.id, ".RData"))
  load(paste0("weight_results/Bucher/variances_", file.id, ".RData")) 
  bucher.means.list[[i]] <- means
  bucher.variances.list[[i]] <- variances
  bucher.bias[i] <- bias(bucher.means.list[[i]], Delta.AB)  
  bucher.bias.mcse[i] <- bias.mcse(bucher.means.list[[i]])
  bucher.mae[i] <- mae(bucher.means.list[[i]], Delta.AB)
  bucher.abs.err[[i]] <- bucher.means.list[[i]] - (Delta.AB)
  bucher.mae.mcse[i] <- mcse.estimate(bucher.abs.err[[i]])
  bucher.mse[i] <- mse(bucher.means.list[[i]], Delta.AB) 
  bucher.mse.mcse[i] <- mse.mcse(bucher.means.list[[i]], Delta.AB) 
  bucher.ate[i] <- mean(bucher.means.list[[i]]) 
  bucher.ate.mcse[i] <- mcse.estimate(bucher.means.list[[i]])
  bucher.lci[[i]] <- bucher.means.list[[i]] + qnorm(0.025)*sqrt(bucher.variances.list[[i]])
  bucher.uci[[i]] <- bucher.means.list[[i]] + qnorm(0.975)*sqrt(bucher.variances.list[[i]])
  bucher.ciwidth[[i]] <- bucher.uci[[i]] - bucher.lci[[i]]
  bucher.lci.mean[i] <- mean(bucher.lci[[i]])
  bucher.lci.mcse[i] <- mcse.estimate(bucher.lci[[i]])
  bucher.uci.mean[i] <- mean(bucher.uci[[i]])
  bucher.uci.mcse[i] <- mcse.estimate(bucher.uci[[i]])
  bucher.cov[i] <- coverage(bucher.lci[[i]], bucher.uci[[i]], Delta.AB)
  bucher.cov.mcse[i] <- coverage.mcse(bucher.cov[i], length(bucher.lci[[i]]))
  bucher.empse[i] <- empse(bucher.means.list[[i]])
  bucher.empse.mcse[i] <- empse.mcse(bucher.empse[i], replicates)
  bucher.vr[i] <- var.ratio(bucher.means.list[[i]], sqrt(bucher.variances.list[[i]]))
  bucher.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(bucher.variances.list[[i]])), 
                                      emp.se=bucher.empse[i],
                                      var.avg.se=mcse.estimate(sqrt(bucher.variances.list[[i]]))^2,
                                      var.emp.se=bucher.empse.mcse[i]^2) 
  bucher.std.bias[i] <- (bucher.bias[i]*100)/bucher.empse[i]
  
  
  # # MAIC2
  # load(paste0("weight_results/MAIC2/means_", file.id, ".RData"))
  # load(paste0("weight_results/MAIC2/variances_", file.id, ".RData"))
  # # number of replicates for which MAIC weighted regression has separability issues
  # maic.separability.list2[i] <- sum(abs(means[!is.na(means)])>max.effect)
  # # discard these replicates (1 of 162000 in this sim. study)
  # maic.no.wei.pro2 <- !is.na(means) & !is.na(variances) #& variances<5 #& variances<1000
  # means <- means[maic.no.wei.pro2]
  # variances <- variances[maic.no.wei.pro2]
  # maic.no.sep2 <- abs(means)<max.effect
  # #means[!maic.no.sep] <- mean(means[abs(means)<max.effect & !is.na(means)])
  # #variances[!maic.no.sep] <- mean(variances[variances<100])
  # means <- means[maic.no.sep2]
  # variances <- variances[maic.no.sep2]
  # #means <- means[variances<5]
  # 
  # #means <- means[!is.na(means)]
  # 
  # 
  # ### MATCHING-ADJUSTED INDIRECT COMPARISON 2(for version2, means, variances, approx.ess.maic)
  # maic.means.list2[[i]] <- means
  # maic.variances.list2[[i]] <- variances
  # maic.bias2[i] <- bias(maic.means.list2[[i]], Delta.AB)
  # maic.bias.mcse2[i] <- bias.mcse(maic.means.list2[[i]])
  # maic.mae2[i] <- mae(maic.means.list2[[i]], Delta.AB)
  # maic.abs.err2[[i]] <- maic.means.list2[[i]] - (Delta.AB)
  # maic.mae.mcse2[i] <- mcse.estimate(maic.abs.err2[[i]])
  # maic.mse2[i] <- mse(maic.means.list2[[i]], Delta.AB)
  # maic.mse.mcse2[i] <- mse.mcse(maic.means.list2[[i]], Delta.AB)
  # maic.ate2[i] <- mean(maic.means.list2[[i]])
  # maic.ate.mcse2[i] <- mcse.estimate(maic.means.list2[[i]])
  # # construct confidence interval using normal distribution
  # maic.lci2[[i]] <- maic.means.list2[[i]] + qnorm(0.025)*sqrt(maic.variances.list2[[i]])
  # maic.uci2[[i]] <- maic.means.list2[[i]] + qnorm(0.975)*sqrt(maic.variances.list2[[i]])
  # maic.ciwidth2[[i]] <- maic.uci2[[i]] - maic.lci2[[i]]
  # maic.lci.mean2[i] <- mean(maic.lci2[[i]])
  # maic.lci.mcse2[i] <- mcse.estimate(maic.lci2[[i]])
  # maic.uci.mean2[i] <- mean(maic.uci2[[i]])
  # maic.uci.mcse2[i] <- mcse.estimate(maic.uci2[[i]])
  # maic.cov2[i] <- coverage(maic.lci2[[i]], maic.uci2[[i]], Delta.AB)
  # maic.cov.mcse2[i] <- coverage.mcse(maic.cov2[i], length(maic.lci2[[i]]))
  # maic.empse2[i] <- empse(maic.means.list2[[i]])
  # maic.empse.mcse2[i] <- empse.mcse(maic.empse2[i], length(maic.means.list2[[i]]))
  # maic.vr2[i] <- var.ratio(maic.means.list2[[i]], sqrt(maic.variances.list2[[i]]))
  # maic.vr.mcse2[i] <- var.ratio.mcse(avg.se=mean(sqrt(maic.variances.list2[[i]])),
  #                                   emp.se=maic.empse2[i],
  #                                   var.avg.se=mcse.estimate(sqrt(maic.variances.list2[[i]]))^2,
  #                                   var.emp.se=maic.empse.mcse2[i]^2)
  # maic.std.bias2[i] <- (maic.bias2[i]*100)/maic.empse2[i]
  # load(paste0("weight_results/MAIC2/aess_", file.id, ".RData"))
  # maic.aess.mean2[i] <- mean(approx.ess.maic)
  # maic.aess.mcse2[i] <- mcse.estimate(approx.ess.maic)


  
  # MAIC2 for verison1,3 (means2, variances2, approx.ess.maic2)
  load(paste0("weight_results/MAIC2/means_", file.id, ".RData"))
  load(paste0("weight_results/MAIC2/variances_", file.id, ".RData"))
  # number of replicates for which MAIC weighted regression has separability issues
  maic.separability.list2[i] <- sum(abs(means2[!is.na(means2)])>max.effect)
  # discard these replicates (1 of 162000 in this sim. study)
  maic.no.wei.pro2 <- !is.na(means2) & !is.na(variances2) #& variances2<5 #& variances<1000
  means2 <- means2[maic.no.wei.pro2]
  variances2 <- variances2[maic.no.wei.pro2]
  maic.no.sep2 <- abs(means2)<max.effect
  #means[!maic.no.sep] <- mean(means[abs(means)<max.effect & !is.na(means)])
  #variances[!maic.no.sep] <- mean(variances[variances<100])
  means2 <- means2[maic.no.sep2]
  variances2 <- variances2[maic.no.sep2]
  #means <- means[variances<5]

  #means <- means[!is.na(means)]


  ### MATCHING-ADJUSTED INDIRECT COMPARISON
  maic.means.list2[[i]] <- means2
  maic.variances.list2[[i]] <- variances2
  maic.bias2[i] <- bias(maic.means.list2[[i]], Delta.AB)
  maic.bias.mcse2[i] <- bias.mcse(maic.means.list2[[i]])
  maic.mae2[i] <- mae(maic.means.list2[[i]], Delta.AB)
  maic.abs.err2[[i]] <- maic.means.list2[[i]] - (Delta.AB)
  maic.mae.mcse2[i] <- mcse.estimate(maic.abs.err2[[i]])
  maic.mse2[i] <- mse(maic.means.list2[[i]], Delta.AB)
  maic.mse.mcse2[i] <- mse.mcse(maic.means.list2[[i]], Delta.AB)
  maic.ate2[i] <- mean(maic.means.list2[[i]])
  maic.ate.mcse2[i] <- mcse.estimate(maic.means.list2[[i]])
  # construct confidence interval using normal distribution
  maic.lci2[[i]] <- maic.means.list2[[i]] + qnorm(0.025)*sqrt(maic.variances.list2[[i]])
  maic.uci2[[i]] <- maic.means.list2[[i]] + qnorm(0.975)*sqrt(maic.variances.list2[[i]])
  maic.ciwidth2[[i]] <- maic.uci2[[i]] - maic.lci2[[i]]
  maic.lci.mean2[i] <- mean(maic.lci2[[i]])
  maic.lci.mcse2[i] <- mcse.estimate(maic.lci2[[i]])
  maic.uci.mean2[i] <- mean(maic.uci2[[i]])
  maic.uci.mcse2[i] <- mcse.estimate(maic.uci2[[i]])
  maic.cov2[i] <- coverage(maic.lci2[[i]], maic.uci2[[i]], Delta.AB)
  maic.cov.mcse2[i] <- coverage.mcse(maic.cov2[i], length(maic.lci2[[i]]))
  maic.empse2[i] <- empse(maic.means.list2[[i]])
  maic.empse.mcse2[i] <- empse.mcse(maic.empse2[i], length(maic.means.list2[[i]]))
  maic.vr2[i] <- var.ratio(maic.means.list2[[i]], sqrt(maic.variances.list2[[i]]))
  maic.vr.mcse2[i] <- var.ratio.mcse(avg.se=mean(sqrt(maic.variances.list2[[i]])),
                                     emp.se=maic.empse2[i],
                                     var.avg.se=mcse.estimate(sqrt(maic.variances.list2[[i]]))^2,
                                     var.emp.se=maic.empse.mcse2[i]^2)
  maic.std.bias2[i] <- (maic.bias2[i]*100)/maic.empse2[i]
  load(paste0("weight_results/MAIC2/aess_", file.id, ".RData"))
  maic.aess.mean2[i] <- mean(approx.ess.maic2)
  maic.aess.mcse2[i] <- mcse.estimate(approx.ess.maic2)


  ### MULTI-LEVEL NETWORK META REGRESSION
 # load(paste0("weight_results/MLNMR/means_", file.id, ".RData"))
 # load(paste0("weight_results/MLNMR/variances_", file.id, ".RData"))  
#  mlnmr.notconverges.list[i] <- sum(variances>max.variance.mlnmr)
  # discard replicates for which mlnmr did not converge (0 of 162000 in this sim. study)
 # mlnmr.converges <- variances<max.variance.mlnmr
 # means <- means[mlnmr.converges]
  #variances <- variances[mlnmr.converges]
  #mlnmr.means.list[[i]] <- means
  #mlnmr.variances.list[[i]] <- variances*-1
  #mlnmr.bias[i] <- bias(mlnmr.means.list[[i]], Delta.AB)  
 # mlnmr.bias.mcse[i] <- bias.mcse(mlnmr.means.list[[i]])
 # mlnmr.mae[i] <- mae(mlnmr.means.list[[i]], Delta.AB)
 # mlnmr.mae.mcse[i] <- mcse.estimate(mlnmr.means.list[[i]])
 # mlnmr.abs.err[[i]] <- mlnmr.means.list[[i]] - (Delta.AB)
 # mlnmr.mae.mcse[i] <- mcse.estimate(mlnmr.abs.err[[i]])
 # mlnmr.mse[i] <- mse(mlnmr.means.list[[i]], Delta.AB) 
 # mlnmr.mse.mcse[i] <- mse.mcse(mlnmr.means.list[[i]], Delta.AB) 
 # mlnmr.ate[i] <- mean(mlnmr.means.list[[i]]) 
 # mlnmr.ate.mcse[i] <- mcse.estimate(mlnmr.means.list[[i]])
 # mlnmr.lci[[i]] <- mlnmr.means.list[[i]] + qnorm(0.025)*sqrt(mlnmr.variances.list[[i]])
 # mlnmr.uci[[i]] <- mlnmr.means.list[[i]] + qnorm(0.975)*sqrt(mlnmr.variances.list[[i]])
 # mlnmr.ciwidth[[i]] <- mlnmr.uci[[i]] - mlnmr.lci[[i]]
 # mlnmr.lci.mean[i] <- mean(mlnmr.lci[[i]])
 # mlnmr.lci.mcse[i] <- mcse.estimate(mlnmr.lci[[i]])
 # mlnmr.uci.mean[i] <- mean(mlnmr.uci[[i]])
  ##mlnmr.uci.mcse[i] <- mcse.estimate(mlnmr.uci[[i]])
 # mlnmr.cov[i] <- coverage(mlnmr.lci[[i]], mlnmr.uci[[i]], Delta.AB)
 # mlnmr.cov.mcse[i] <- coverage.mcse(mlnmr.cov[i], length(mlnmr.lci[[i]]))
 # mlnmr.empse[i] <- empse(mlnmr.means.list[[i]])
 # mlnmr.empse.mcse[i] <- empse.mcse(mlnmr.empse[i], length(mlnmr.means.list[[i]]))
 # mlnmr.vr[i] <- var.ratio(mlnmr.means.list[[i]], sqrt(mlnmr.variances.list[[i]]))
 # mlnmr.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(mlnmr.variances.list[[i]])), 
                                #     emp.se=mlnmr.empse[i],
                                #     var.avg.se=mcse.estimate(sqrt(mlnmr.variances.list[[i]]))^2,
                                #     var.emp.se=mlnmr.empse.mcse[i]^2)  
 # mlnmr.std.bias[i] <- (mlnmr.bias[i]*100)/mlnmr.empse[i]
  

  
  
  truth <- Delta.AB # true baseline A vs. B treatment effect
  maic.error.worse.than[i] <- sum(abs(bucher.means.list[[i]][maic.no.sep]-truth)<abs(maic.means.list[[i]]-truth))/sum(maic.no.sep)
  maic.error.worse.than.mcse[i] <- coverage.mcse(maic.error.worse.than[i], sum(maic.no.sep))
  maic.error.worse.than2[i] <- sum(abs(bucher.means.list[[i]][maic.no.sep]-truth)<abs(maic.means.list2[[i]]-truth))/sum(maic.no.sep)
  maic.error.worse.than.mcse2[i] <- coverage.mcse(maic.error.worse.than2[i], sum(maic.no.sep))
  
  stc.error.worse.than[i] <- sum(abs(bucher.means.list[[i]][stc.converges]-truth)<abs(stc.means.list[[i]]-truth))/sum(stc.converges)
  stc.error.worse.than.mcse[i] <- coverage.mcse(stc.error.worse.than[i], sum(stc.converges))  
  # mlnmr.error.worse.than[i] <- sum(abs(bucher.means.list[[i]][mlnmr.converges]-truth)<abs(mlnmr.means.list[[i]]-truth))/sum(mlnmr.converges)
 # mlnmr.error.worse.than.mcse[i] <- coverage.mcse(mlnmr.error.worse.than[i], sum(mlnmr.converges)) 
  
  
  
  
  
  
  tmp.scenarios <- cbind(i, pc$N_AC[i], pc$b_X[i], pc$b_EM[i], pc$meanX_AC[i], pc$corX[i])
  
  maic.tmp.metrics <- cbind(maic.ate[i], maic.ate.mcse[i], maic.lci.mean[i],
                            maic.lci.mcse[i], maic.uci.mean[i], maic.uci.mcse[i],
                            maic.bias[i], maic.bias.mcse[i], maic.mse[i], maic.mse.mcse[i],
                            maic.mae[i], maic.mae.mcse[i], maic.cov[i], maic.cov.mcse[i],
                            maic.empse[i], maic.empse.mcse[i], maic.vr[i], maic.vr.mcse[i],
                            maic.error.worse.than[i], maic.error.worse.than.mcse[i], 
                            maic.std.bias[i], maic.aess.mean[i], maic.aess.mcse[i])
  
  maic2.tmp.metrics <- cbind(maic.ate2[i], maic.ate.mcse2[i], maic.lci.mean2[i],
                            maic.lci.mcse2[i], maic.uci.mean2[i], maic.uci.mcse2[i],
                            maic.bias2[i], maic.bias.mcse2[i], maic.mse2[i], maic.mse.mcse2[i],
                            maic.mae2[i], maic.mae.mcse2[i], maic.cov2[i], maic.cov.mcse2[i],
                            maic.empse2[i], maic.empse.mcse2[i], maic.vr2[i], maic.vr.mcse2[i],
                            maic.error.worse.than2[i], maic.error.worse.than.mcse2[i], 
                            maic.std.bias2[i], maic.aess.mean2[i], maic.aess.mcse2[i])
  
  
  stc.tmp.metrics <- cbind(stc.ate[i], stc.ate.mcse[i], stc.lci.mean[i],
                           stc.lci.mcse[i], stc.uci.mean[i], stc.uci.mcse[i],
                           stc.bias[i], stc.bias.mcse[i], stc.mse[i], stc.mse.mcse[i],
                           stc.mae[i], stc.mae.mcse[i], stc.cov[i], stc.cov.mcse[i],
                           stc.empse[i], stc.empse.mcse[i], stc.vr[i], stc.vr.mcse[i],
                           stc.error.worse.than[i], stc.error.worse.than.mcse[i],
                           stc.std.bias[i], stc.notconverges.list[i])
  
  bucher.tmp.metrics <- cbind(bucher.ate[i], bucher.ate.mcse[i], bucher.lci.mean[i],
                              bucher.lci.mcse[i], bucher.uci.mean[i], bucher.uci.mcse[i],
                              bucher.bias[i], bucher.bias.mcse[i], bucher.mse[i], bucher.mse.mcse[i],
                              bucher.mae[i], bucher.mae.mcse[i], bucher.cov[i], bucher.cov.mcse[i],
                              bucher.empse[i], bucher.empse.mcse[i], bucher.vr[i], bucher.vr.mcse[i], 
                              bucher.std.bias[i]) 
  
  #mlnmr.tmp.metrics <- cbind(mlnmr.ate[i], mlnmr.ate.mcse[i], mlnmr.lci.mean[i],
                            # mlnmr.lci.mcse[i], mlnmr.uci.mean[i], mlnmr.uci.mcse[i],
                            # mlnmr.bias[i], mlnmr.bias.mcse[i], mlnmr.mse[i], mlnmr.mse.mcse[i],
                            # mlnmr.mae[i], mlnmr.mae.mcse[i], mlnmr.cov[i], mlnmr.cov.mcse[i],
                            # mlnmr.empse[i], mlnmr.empse.mcse[i], mlnmr.vr[i], mlnmr.vr.mcse[i],
                            # mlnmr.error.worse.than[i], mlnmr.error.worse.than.mcse[i],
                            # mlnmr.std.bias[i], mlnmr.notconverges.list[i])
  
  tmp.scenarios<- cbind(tmp.scenarios, maic.tmp.metrics, maic2.tmp.metrics, stc.tmp.metrics, bucher.tmp.metrics)
  
  scenarios.df <- rbind(scenarios.df, tmp.scenarios)
  
 # mlnmr_bias_nonsum <- biaspersim(mlnmr.means.list[[i]], Delta.AB)  
  stc_bias_nonsum <- biaspersim(stc.means.list[[i]], Delta.AB)  
  maic_bias_nonsum <- biaspersim(maic.means.list[[i]], Delta.AB)  
  bucher_bias_nonsum <- biaspersim(bucher.means.list[[i]], Delta.AB) 
  maic_bias_nonsum2 <- biaspersim(maic.means.list2[[i]], Delta.AB)
  
  bias_all <-cbind(bucher_bias_nonsum,maic_bias_nonsum,maic_bias_nonsum2,stc_bias_nonsum) 
}  

colnames(scenarios.df) <- c("Scenario", "N_AC", "b_X","b_EM","meanX_AC", "corX",
                            "maic.ate", "maic.ate.mcse", "maic.lci.mean",
                            "maic.lci.mcse", "maic.uci.mean", "maic.uci.mcse",
                            "maic.bias", "maic.bias.mcse", "maic.mse", "maic.mse.mcse",
                            "maic.mae", "maic.mae.mcse", "maic.cov", "maic.cov.mcse",
                            "maic.empse", "maic.empse.mcse", "maic.vr", "maic.vr.mcse",
                            "maic.error.worse.than", "maic.error.worse.than.mcse", 
                            "maic.std.bias", "maic.aess.mean", "maic.aess.mcse", 
                            "maic.ate2", "maic.ate.mcse2", "maic.lci.mean2",
                            "maic.lci.mcse2", "maic.uci.mean2", "maic.uci.mcse2",
                            "maic.bias2", "maic.bias.mcse2", "maic.mse2", "maic.mse.mcse2",
                            "maic.mae2", "maic.mae.mcse2", "maic.cov2", "maic.cov.mcse2",
                            "maic.empse2", "maic.empse.mcse2", "maic.vr2", "maic.vr.mcse2",
                            "maic.error.worse.than2", "maic.error.worse.than.mcse2", 
                            "maic.std.bias2", "maic.aess.mean2", "maic.aess.mcse2", 
                            "stc.ate", "stc.ate.mcse", 
                            "stc.lci.mean", "stc.lci.mcse", "stc.uci.mean", "stc.uci.mcse",
                            "stc.bias", "stc.bias.mcse", "stc.mse", "stc.mse.mcse",
                            "stc.mae", "stc.mae.mcse", "stc.cov", "stc.cov.mcse",
                            "stc.empse", "stc.empse.mcse", "stc.vr", "stc.vr.mcse",
                            "stc.error.worse.than", "stc.error.worse.than.mcse",
                            "stc.std.bias","stc.no.convergence", 
                            "bucher.ate", "bucher.ate.mcse", "bucher.lci.mean",
                            "bucher.lci.mcse", "bucher.uci.mean", "bucher.uci.mcse",
                            "bucher.bias", "bucher.bias.mcse", "bucher.mse", "bucher.mse.mcse",
                            "bucher.mae", "bucher.mae.mcse", "bucher.cov", "bucher.cov.mcse",
                            "bucher.empse", "bucher.empse.mcse", "bucher.vr", 
                            "bucher.vr.mcse", "bucher.std.bias")
# ,   "mlnmr.ate", "mlnmr.ate.mcse",
# "mlnmr.lci.mean", "mlnmr.lci.mcse", "mlnmr.uci.mean", "mlnmr.uci.mcse",
#   "mlnmr.bias", "mlnmr.bias.mcse", "mlnmr.mse", "mlnmr.mse.mcse",
#   "mlnmr.mae", "mlnmr.mae.mcse", "mlnmr.cov", "mlnmr.cov.mcse",
#   "mlnmr.empse", "mlnmr.empse.mcse", "mlnmr.vr", "mlnmr.vr.mcse",
#   "mlnmr.error.worse.than", "mlnmr.error.worse.than.mcse",
#   "mlnmr.std.bias","mlnmr.no.convergence")

write.csv(scenarios.df, paste0("Analysis/scenarios/","marg",".csv"), row.names = FALSE)


################################################################################
###################################SCATTERPLOTS#################################
################################################################################


bias_all <- as.data.frame(bias_all)
longbias_all <- bias_all %>% pivot_longer(cols = `bucher_bias_nonsum`:`stc_bias_nonsum`, names_to = "ITC method", values_to = "Bias")
ggplot(longbias_all, aes(`ITC method`, Bias)) +
  geom_point(shape=1)

  
################################################################################
###############################TOTAL AVERAGE weight_results############################
################################################################################


scenarios.df <- read.csv(file="Analysis/scenarios/marg.csv", header=TRUE, sep=",")
#scenarios.dfmarg <- read.csv(file="Analysis/scenarios/marg.csv", header=TRUE, sep=",")
#incremental <- scenarios.dfcond- scenarios.dfmarg

# scenarios.dfmarg$b_EM <- round(scenarios.df$b_EM, 2)
# scenarios.dfmarg$b_X <- round(scenarios.df$b_X, 2)

scenarios.df$b_EM <- round(scenarios.df$b_EM, 2)
scenarios.df$b_X <- round(scenarios.df$b_X, 2)


################################################################################
###################################NESTED LOOP PLOTS############################
################################################################################
#RUN BELOW ONLY FOR FULL FACTORIAL!!!!

### MAIN MANUSCRIPT PLOTS ###

# reorder simulation weight_results in nested loop plot presentation order (FOR MSE)
nldata <- nestedloop(scenarios.df,
                     varnames=c("b_EM","N_AC","meanX_AC","b_X","corX"),
                     varlabels=
                       c("Effect-modifiying interaction", 
                         "Subjects in trial with patient-level data", "Covariate overlap",
                         "Prognostic variable effect", 
                         "Covariate correlation"),
                     sign=c(1, 1, 1, 1, 1))
pd.mse <- nldata
pd.mse$meanX_AC <- factor(pd.mse$meanX_AC,levels=c(0.15,0.3,0.45),
                          labels=c("poor","moderate","strong"))
pd.mse$corX <- factor(pd.mse$corX, levels=c(0,0.35),labels=c("none","moderate"))
pd.mse$b_EM <- factor(pd.mse$b_EM,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                    "very strong"))
pd.mse$b_X <- factor(pd.mse$b_X,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                  "very strong"))

### Nested loop plot for MSE (TIFF)
pdf("Analysis/mse.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.mse$maic.mse,
     type="n",
     ylim=c(0,0.7), bty="n",
     xlab="Scenario",
     ylab="Mean square error (MSE)",
     las=1, xaxt="n",
     cex.axis=1.15, cex.lab=1.15) 
lines(pd.mse, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.mse, which="r", ymin.refline=0.5, ymax.refline=0.7, cex.ref=1.0)
lines(pd.mse$maic.mse, col="red", lty=2,type="s", lwd=2.5)
lines(pd.mse$stc.mse, col="green", lty=5, type="s", lwd=2.5)
lines(pd.mse$bucher.mse, col="yellow", lty=1, type="s", lwd=2.5)
lines(pd.mse$maic.mse2, col="blue", lty=4, type="s", lwd=2.5)
legend("topright",lwd=c(2.5,2.5,2.5,2.5),col=c("red","green","yellow","blue"), lty=c(2,5,1,4),
       cex=1.2,bty="n",c("MAIC", "STC","Bucher","MAIC2"))
dev.off()

# ### Nested loop plot for MSE (PDF)
# pdf("Analysis/mse.pdf", width=18,height=15,pointsize=20)
# par(pty="m")
# plot(pd.mse$maic.mse,
#      type="n",
#      ylim=c(0,0.6), bty="n",
#      xlab="Scenario",
#      ylab="Mean square error (MSE)",
#      las=1, xaxt="n",
#      cex.axis=1.15, cex.lab=1.15)
# lines(pd.mse, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
# lines(pd.mse, which="r", ymin.refline=0.45, ymax.refline=0.6, cex.ref=0.9) 
# lines(pd.mse$maic.mse, col="red", lty=2,type="s", lwd=2.5)
# lines(pd.mse$stc.mse, col="green", lty=5, type="s", lwd=2.5)
# lines(pd.mse$bucher.mse, col="blue", lty=4, type="s", lwd=2.5)
# lines(pd.mse$maic.mse2, col="orange", lty=3, type="s", lwd=2.5)
# legend("left",lwd=c(2.5,2.5,2.5,2.5),col=c("red","green","blue","orange"), lty=c(2,5,4,3),
#        cex=1.15,bty="n",c("MAIC", "STC","Bucher","MAIC2"))
# dev.off()

# ### Nested loop plot for empirical standard error (TIFF)
# tiff("Analysis/EmpSE.tiff", res=800, width=18, height=15, units='in')
# par(pty="m")
# plot(pd.ese$maic.empse, type="n",ylim=c(0, 0.6), bty="n", xlab="Scenario",
#      ylab="Empirical standard error (ESE)",las=1, xaxt="n",
#      cex.axis=1.5, cex.lab=1.5)
# lines(pd.ese, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
# lines(pd.ese, which="r", ymin.refline=0.45, ymax.refline=0.6, cex.ref=1.5)
# lines(pd.ese$maic.empse, col="red", lty=2,type="s",lwd=2.5)    
# lines(pd.ese$stc.empse, col="green", lty=5, type="s", lwd=2.5)   
# lines(pd.ese$bucher.empse, col="blue", lty=4, type="s", lwd=2.5) 
# lines(pd.ese$bucher.empse, col="orange", lty=3, type="s", lwd=2.5) 
# legend("bottom", lwd=c(2.5,2.5,2.5,2.5), col=c("red","green","blue","orange"), lty=c(2,5,4,3), 
#        cex=1.5, bty="n", c("MAIC", "STC", "Bucher","MIAC2"))
# dev.off()


# reorder simulation weight_results in nested loop plot presentation order (for empirical standard error)
nldata <- nestedloop(scenarios.df,
                     varnames=c("N_AC","b_EM","meanX_AC","b_X","corX"),
                     varlabels=
                       c("Subjects in trial with patient-level data",
                         "Effect-modifiying interaction", 
                         "Covariate overlap",
                         "Prognostic variable effect", 
                         "Covariate correlation"),
                     sign=c(1, 1, 1, 1, 1))
pd.ese <- nldata
pd.ese$meanX_AC <- factor(pd.ese$meanX_AC,levels=c(0.15,0.3,0.45),
                          labels=c("poor","moderate","strong"))
pd.ese$corX <- factor(pd.ese$corX, levels=c(0,0.35),labels=c("none","moderate"))
pd.ese$b_EM <- factor(pd.ese$b_EM,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                    "very strong"))
pd.ese$b_X <- factor(pd.ese$b_X,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                  "very strong"))
### Nested loop plot for empirical standard error (PDF)
pdf("Analysis/EmpSE.pdf",width=18, height=15,pointsize=20)
par(pty="m")
plot(pd.ese$maic.empse, type="n",ylim=c(0.0, 0.9), bty="n", xlab="Scenario",
     ylab="Empirical standard error (ESE)",las=1, xaxt="n",
     cex.axis=1.15, cex.lab=1.15)
lines(pd.ese, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.ese, which="r", ymin.refline=0.70, ymax.refline=0.9, cex.ref=1.0)
lines(pd.ese$maic.empse, col="red", lty=2,type="s",lwd=2.5)    
lines(pd.ese$stc.empse, col="green", lty=5, type="s", lwd=2.5)   
lines(pd.ese$bucher.empse, col="yellow", lty=1, type="s", lwd=2.5) 
lines(pd.ese$maic.empse2, col="blue", lty=4, type="s", lwd=2.5) 
legend("topright", lwd=c(2.5,2.5,2.5,2.5), col=c("red","green","yellow","blue"), lty=c(2,5,1,4), 
       cex=1.2, bty="n", c("MAIC", "STC", "Bucher","MAIC2"))
dev.off()

# ### Nested loop plot for coverage (TIFF)
# tiff("Analysis/coverage.tiff", res=800, width=18, height=15, units='in')
# par(pty="m")
# plot(pd.cover$maic.cov*100, type="n", ylim=c(20, 100), bty="n", xlab="Scenario",
#      ylab="Coverage of 95% confidence intervals (%)", las=1, xaxt="n",
#      cex.axis=1.5, cex.lab=1.5)
# abline(h=95, col="grey") # nominal 
# lines(pd.cover, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
# lines(pd.cover, which="r", ymin.refline=20, ymax.refline=40, cex.ref=1.5)
# lines(pd.cover$maic.cov*100, col="red", lty=2,type="s", lwd=2.5)    
# lines(pd.cover$stc.cov*100, col="green", lty=5, type="s", lwd=2.5) 
# lines(pd.cover$bucher.cov*100, col="blue", lty=4, type="s", lwd=2.5) 
# lines(pd.cover$bucher.cov*100, col="orange", lty=3, type="s", lwd=2.5) 
# legend("right", lwd=c(2.5,2.5,2.5,2.5), col=c("red","green", "blue","orange"), 
#        lty=c(2,5,4,3), cex=1.5, bty="n",c("MAIC", "STC", "Bucher","MAIC2"))
# dev.off()



# reorder simulation weight_results in nested loop plot presentation order (for coverage)
nldata <- nestedloop(scenarios.df,
                     varnames=c("meanX_AC","b_EM","N_AC","b_X","corX"),
                     varlabels=
                       c("Covariate overlap", "Effect-modifiying interaction", 
                         "Subjects in trial with patient-level data",
                         "Prognostic variable effect", 
                         "Covariate correlation"),
                     sign=c(1, 1, 1, 1, 1))
pd.cover <- nldata
pd.cover$meanX_AC <- factor(pd.cover$meanX_AC,levels=c(0.15,0.3,0.45),
                            labels=c("poor","moderate","strong"))
pd.cover$corX <- factor(pd.cover$corX, levels=c(0,0.35),labels=c("none","moderate"))
pd.cover$b_EM <- factor(pd.cover$b_EM,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                        "very strong"))
pd.cover$b_X <- factor(pd.cover$b_X,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                      "very strong"))

nldata <- nestedloop(scenarios.df,
                     varnames=c("b_EM","meanX_AC","b_X","N_AC","corX"),
                     varlabels=c("Effect-modifying interaction", "Covariate overlap",
                                 "Prognostic variable effect", "Subjects in trial with patient-level data",
                                 "Covariate correlation"),
                     sign=c(1, 1, 1, 1, 1))
pd.bias <- nldata # object to be used in nested loop plot
### Nested loop plot for coverage (PDF)
pdf("Analysis/coverage.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.bias$maic.cov*100, type="n", ylim=c(0, 140), bty="n", xlab="Scenario",
     ylab="Coverage of 95% confidence intervals (%)", las=1, xaxt="n",
     cex.axis=1.15, cex.lab=1.15)
abline(h=95, col="grey") # nominal 
lines(pd.cover, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.cover, which="r", ymin.refline=110, ymax.refline=140, cex.ref=1.0)
lines(pd.cover$maic.cov*100, col="magenta", lty=2,type="s", lwd=2.5)    
lines(pd.cover$stc.cov*100, col="green", lty=5, type="s", lwd=2.5) 
lines(pd.cover$bucher.cov*100, col="yellow", lty=1, type="s", lwd=2.5)
lines(pd.cover$maic.cov2*100, col="cyan", lty=4, type="s", lwd=2.5) 
legend("topright", lwd=c(2.5,2.5,2.5,2.5), col=c("magenta","green", "yellow","cyan"), 
       lty=c(2,5,1,4), cex=1.2, bty="n",c("MAIC", "STC", "Bucher","MAIC2"))
dev.off()

# ### Nested loop plot for variability ratio (TIFF)
# tiff("Analysis/variability_ratio.tiff", res=800, width=18, height=15, units='in')
# par(pty="m")
# plot(pd.cover$maic.vr, type="n", ylim=c(0.8,1.2), bty="n", xlab="Scenario", ylab="Variability ratio",
#      las=1, xaxt="n", cex.axis=1.5, cex.lab=1.5)
# abline(h=1, col="grey") # unbiased variance
# lines(pd.cover, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
# lines(pd.cover, which="r", ymin.refline=1.1, ymax.refline=1.2, cex.ref=1.5)
# lines(pd.cover$maic.vr, col="red", lty=2,type="s",lwd=2.5)    
# lines(pd.cover$stc.vr, col="green", lty=5, type="s", lwd=2.5)   
# lines(pd.cover$bucher.vr, col="blue", lty=4, type="s", lwd=2.5) 
# lines(pd.cover$maic.vr2, col="orange", lty=3, type="s", lwd=2.5) 
# legend("bottom", lwd=c(2.5,2.5,2.5,2.5), col=c("red","green","blue","orange"),lty=c(2,5,4,3), 
#        cex=1.5, bty="n", c("MAIC", "STC", "Bucher","MAIC2"))
# dev.off()




### Nested loop plot for variability ratio (PDF)
pdf("Analysis/variability_ratio.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.cover$maic.vr, type="n", ylim=c(0.4,1.5), bty="n", xlab="Scenario", ylab="Variability ratio",
     las=1, xaxt="n", cex.axis=1.15, cex.lab=1.15)
abline(h=1, col="grey") # unbiased variance
lines(pd.cover, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.cover, which="r", ymin.refline=1.2, ymax.refline=1.5, cex.ref=0.9)
lines(pd.cover$maic.vr, col="red", lty=2,type="s",lwd=2.5)    
lines(pd.cover$stc.vr, col="green", lty=5, type="s", lwd=2.5)   
lines(pd.cover$bucher.vr, col="yellow", lty=1, type="s", lwd=2.5) 
lines(pd.cover$maic.vr2, col="blue", lty=4, type="s", lwd=2.5) 
legend("topright", lwd=c(2.5,2.5,2.5,2.5), col=c("red","green","yellow","blue"),lty=c(2,5,1,4), 
       cex=1.2, bty="n", c("MAIC", "STC", "Bucher","MIAC2"))
dev.off()

# ### Nested loop plot for bias (TIFF)
# tiff("Analysis/bias.tiff", res=800, width=18, height=15, units='in')
# par(pty="m")
# plot(pd.bias$maic.bias, type="n", ylim=c(-1,0.6), bty="n", xlab="Scenario", ylab="Bias",
#      las=1, xaxt="n", cex.axis=1.5, cex.lab=1.5)
# abline(h=0, col="grey") # no bias
# lines(pd.bias, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")) # add vertical lines
# lines(pd.bias, which="r",ymin.refline=0.2, ymax.refline=0.6, cex.ref=1.5) # add reference lines
# lines(pd.bias$maic.bias, col="red", lty=2,type="s", lwd=2.5) # performance measures 
# lines(pd.bias$stc.bias, col="green", lty=5, type="s", lwd=2.5)
# lines(pd.bias$bucher.bias, col="blue", lty=4, type="s", lwd=2.5)
# lines(pd.bias$maic.bias2, col="orange", lty=3, type="s", lwd=2.5)
# legend("bottom",lwd=c(2.5,2.5,2.5,2.5),col=c("red","green", "blue","orange"),
#        lty=c(2,5,4,3),cex=1.5,bty="n", c("MAIC", "STC", "Bucher", "MAIC2")) # legend
# dev.off()


# reorder simulation weight_results in order to be presented in nested loop plot (FOR BIAS).  
# consider re-ordering on a case-by-case basis for each plot
# (function is by R?cker, G., Schwarzer, G. Presenting simulation weight_results in a nested loop plot. 
# BMC Med Res Methodol 14, 129 (2014) doi:10.1186/1471-2288-14-129)
nldata <- nestedloop(scenarios.df,
                     varnames=c("b_EM","meanX_AC","b_X","N_AC","corX"),
                     varlabels=c("Effect-modifying interaction", "Covariate overlap",
                                 "Prognostic variable effect", "Subjects in trial with patient-level data",
                                 "Covariate correlation"),
                     sign=c(1, 1, 1, 1, 1))
pd.bias <- nldata # object to be used in nested loop plot
# use labels instead of numeric values for the following factors
pd.bias$meanX_AC <- factor(pd.bias$meanX_AC, levels=c(0.15,0.3,0.45),labels=c("poor","moderate","strong"))
pd.bias$corX <- factor(pd.bias$corX, levels=c(0,0.35),labels=c("none","moderate"))
pd.bias$b_EM <- factor(pd.bias$b_EM,levels=c(0.40,0.69,1.11), labels=c("moderate","strong",
                                                                       "very strong"))
pd.bias$b_X <- factor(pd.bias$b_X, levels=c(0.40,0.69,1.11), labels=c("moderate","strong",
                                                                      "very strong"))

### Nested loop plot for bias (PDF)
pdf("Analysis/bias.pdf", width=18, height=15, pointsize=20)
par(pty="m")
plot(pd.bias$maic.bias, type="n", ylim=c(-0.8,0.6), bty="n", xlab="Scenario", ylab="Bias",
     las=1, xaxt="n", cex.axis=1.15, cex.lab=1.15)
abline(h=0, col="grey") # no bias
lines(pd.bias, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")) # add vertical lines
lines(pd.bias, which="r",ymin.refline=0.2, ymax.refline=0.6, cex.ref=0.9) # add reference lines
lines(pd.bias$maic.bias, col="red", lty=2,type="s", lwd=2.5) # performance measures 
lines(pd.bias$stc.bias, col="green", lty=5, type="s", lwd=2.5)
lines(pd.bias$bucher.bias, col="yellow", lty=1, type="s", lwd=2.5)
lines(pd.bias$maic.bias2, col="blue", lty=4, type="s", lwd=2.5)
legend("topright",lwd=c(2.5,2.5,2.5),col=c("red","green", "yellow","blue"),
       lty=c(2,5,1,4),cex=1.2,bty="n", c("MAIC", "STC", "Bucher","MAIC2")) # legend
dev.off()


#### SUPPLEMENTARY MATERIAL PLOTS ####

### Nested loop plot for MAE
pdf("Analysis/Supplementary_Material/mae.pdf", width=18, height=15, pointsize=20)
par(pty="m")
plot(pd.mse$maic.mae,type="n",ylim=c(0.1, 0.9),bty="n",xlab="Scenario",
     ylab="Mean absolute error (MAE)",las=1, xaxt="n", cex.axis=1.15, cex.lab=1.15)
lines(pd.mse, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.mse, which="r",ymin.refline=0.72, ymax.refline=0.9,cex.ref=0.9)
lines(pd.mse$maic.mae, col="red", lty=2,type="s", lwd=2.5)
lines(pd.mse$stc.mae, col="green", lty=5, type="s", lwd=2.5)
lines(pd.mse$bucher.mae, col="yellow", lty=1, type="s", lwd=2.5)
lines(pd.mse$maic.mae2, col="blue", lty=4, type="s", lwd=2.5)
legend("topright",lwd=c(2.5,2.5,2.5),col=c("red","green","yellow","blue"),lty=c(2,5,1,4),
       cex=1.15,bty="n",c("MAIC", "STC", "Bucher","MAIC2"))
dev.off()

### Nested loop plot of standardized biases
pdf("Analysis/Supplementary_Material/std_biases.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.bias$maic.std.bias, ylim=c(-550,300), bty="n", xlab="Scenario", ylab="Standardized percentage bias",
     las=1, xaxt="n",cex.axis=1.15, cex.lab=1.15)
abline(h=0, col="grey") # no bias
lines(pd.bias, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.bias, which="r", ymin.refline=100, ymax.refline=300, cex.ref=1.0)
lines(pd.bias$maic.std.bias, col="red", lty=2,type="s", lwd=2.5)
lines(pd.bias$stc.std.bias, col="green", lty=5, type="s", lwd=2.5)
lines(pd.bias$bucher.std.bias, col="yellow", lty=1, type="s", lwd=2.5)
lines(pd.bias$maic.std.bias2, col="blue", lty=4, type="s", lwd=2.5)
legend("topright",lwd=c(2.5,2.5,2.5), col=c("red","green","yellow","blue"), lty=c(2,5,1,4),
       cex=1.15, bty="n", c("MAIC", "STC", "Bucher","MAIC2"))
dev.off()

### Nested loop plot of confidence interval width
pdf("Analysis/Supplementary_Material/CI_width.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.mse$maic.uci.mean-pd.mse$maic.lci.mean, type="n", ylim=c(0.25, 4), bty="n", xlab="Scenario",
     ylab="95% confidence interval width", las=1, xaxt="n",cex.axis=1.15, cex.lab=1.15)
lines(pd.mse, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.mse, which="r", ymin.refline=3.3, ymax.refline=4, cex.ref=0.9)
lines(pd.mse$maic.uci.mean-pd.mse$maic.lci.mean, col="red", lty=2,type="s", lwd=2.5)
lines(pd.mse$stc.uci.mean-pd.mse$stc.lci.mean, col="green", lty=5, type="s", lwd=2.5)
lines(pd.mse$bucher.uci.mean-pd.mse$bucher.lci.mean, col="yellow", lty=1, type="s", lwd=2.5)
lines(pd.mse$maic.uci.mean2-pd.mse$maic.lci.mean2, col="blue", lty=4, type="s", lwd=2.5)
legend("topright", lwd=c(2.5,2.5,2.5), col=c("red","green","yellow","blue"), lty=c(2,5,1,4),
       cex=1.15, bty="n", c("MAIC", "STC", "Bucher","MAIC2"))
dev.off()

### Nested loop plot of % replicates worse than Bucher
pdf("Analysis/Supplementary_Material/worse_than_bucher.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot(pd.cover$maic.error.worse.than*100, type="n", ylim=c(0, 100), bty="n", xlab="Scenario",
     ylab="% of point estimates worse than Bucher", las=1, xaxt="n",cex.axis=1.15, cex.lab=1.15)
lines(pd.cover, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.cover, which="r", ymin.refline=75, ymax.refline=100, cex.ref=0.9)
lines(pd.cover$maic.error.worse.than*100, col="red", lty=2,type="s", lwd=2.5)
lines(pd.cover$stc.error.worse.than*100, col="green", lty=5, type="s", lwd=2.5)
lines(pd.cover$maic.error.worse.than2*100, col="blue", lty=4, type="s", lwd=2.5)
legend("topright", lwd=c(2.5, 2.5), col=c("red","green","blue"), lty=c(2,5,4),
       cex=1.15, bty="n", c("MAIC", "STC","MAIC2"))
dev.off()

# reorder simulation weight_results in nested loop plot presentation order (FOR % reduction in ESS)
nldata <- nestedloop(scenarios.df,
                     varnames=c("meanX_AC","corX","b_EM","N_AC","b_X"),
                     varlabels=
                       c("Covariate overlap", "Covariate correlation",
                         "Effect-modifiying interaction", 
                         "Subjects in trial with patient-level data",
                         "Prognostic variable effect"),
                     sign=c(1, 1, 1, 1, 1))
pd.aess <- nldata
pd.aess$meanX_AC <- factor(pd.aess$meanX_AC,levels=c(0.15,0.3,0.45),
                           labels=c("poor","moderate","strong"))
pd.aess$corX <- factor(pd.aess$corX, levels=c(0,0.35),labels=c("none","moderate"))
pd.aess$b_EM <- factor(pd.aess$b_EM,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                      "very strong"))
pd.aess$b_X <- factor(pd.aess$b_X,levels=c(0.40,0.69,1.11),labels=c("moderate","strong",
                                                                    "very strong"))


### Nested loop plot of MAIC % reduction in effective sample size
pdf("Analysis/Supplementary_Material/maic_percentage_reduction_ess.pdf", width=18,height=15,pointsize=20)
par(pty="m")
plot((pd.aess$N_AC-pd.aess$maic.aess.mean)*100/pd.aess$N_AC, type="n", ylim=c(20,125),
     bty="n", xlab="Scenario", ylab="% reduction in effective sample size", las=1, xaxt="n",
     cex.axis=1.15, cex.lab=1.15)
lines(pd.aess, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.aess, which="r", ymin.refline=105, ymax.refline=125, cex.ref=1.0)
lines((pd.aess$N_AC-pd.aess$maic.aess.mean)*100/pd.aess$N_AC, col="red", lty=2,type="s", lwd=2.5)
lines((pd.aess$N_AC-pd.aess$maic.aess.mean2)*100/pd.aess$N_AC, col="blue", lty=4,type="s", lwd=2.5)

legend("bottomleft", lwd=c(2.5,2.5), col=c("red","blue"), lty=c(2,4), cex=1.2, bty="n", c("MAIC","MAIC2"))
dev.off()


# Plot of bias converging over a simulation scenario (simulation scenario 1)
maic.rolling.bias.1 <- cumsum(maic.means.list[[1]])/(1:replicates) # moving averages
maic.rolling.bias.2 <- cumsum(maic.means.list2[[1]])/(1:replicates) # moving averages
stc.rolling.bias.1 <- cumsum(stc.means.list[[1]])/(1:replicates)
bucher.rolling.bias.1 <- cumsum(bucher.means.list[[1]])/(1:replicates)

pdf("Analysis/Supplementary_material/rolling_bias_convergence.pdf", width=18, height=15, pointsize=20)
plot(maic.rolling.bias.1, xlab="Simulation number", ylab="Rolling bias", col="red",
     type="l",ylim=c(-0.3,0.2), lwd=2, lty=2, cex.axis=1.15, cex.lab=1.15)
abline(h=0, col="grey") # no bias
lines(stc.rolling.bias.1, col="green", lwd=2, lty=5)
lines(bucher.rolling.bias.1, col="yellow", lwd=2, lty=1)
lines(maic.rolling.bias.2, col="blue", lwd=2, lty=4)
legend("topright", legend = c("MAIC", "STC", "Bucher","MAIC2"),
       col = c("red", "green","yellow","blue"), bty = "n",
       pt.cex = 2, cex = 1.2, text.col = "black", horiz = F , inset = c(0.1, 0.1),
       lwd=c(2,2,2,2), lty=c(2,5,1,4))
dev.off()

