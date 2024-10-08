# This file generates the simulation study data

#rm(list=ls())

# setwd("C:/Users/Antonio/Desktop/population_adjustment_simstudy") 

#load(file = "survival_settings.RData") # load simulation setup specifics

if (!dir.exists("Data_var/")) {dir.create("Data_var/")}

# package with copula function to simulate the covariates
if (!require(simstudy)) {install.packages("simstudy"); library(simstudy)}
# package for data manipulation
if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
# Cox regression to summarize the outcomes for the ALD trial in terms of log HR and its variance
if (!require(survival)) {install.packages("survival"); library(survival)}

set.seed(555) # set random seed

#Datasets below varied in:
#1.Level correlation between covariates (corX 0 or 0.35)

#-------------------------------------------------------------------------------
#Replaced cens_rate with cens_rateA and censrateB
gen.data2 <- function(no.chars, no.ems, N_AC, N_BC, b_trt_A, b_trt_B, 
                     b_X, b_EM_A, b_EM_B, meanX_AC, meanX_BC, sdX, weib_shape, 
                     weib_inv_scale, cens_rateA, cens_rateB, corX, allocation) {
  R <- matrix(corX, nrow = no.chars, ncol = no.chars) # set correlation matrix
  diag(R) <- rep(1, no.chars)#Diagonal is of course 1
  states.data2 <- .Random.seed # track random-number-generator state
  N_AC_A <- floor(N_AC*allocation) # number of patients under A in AC
  N_AC_C <- N_AC - N_AC_A # number of patients under C in AC
  N_BC_B <- floor(N_BC*allocation) # number of patients under B in BC
  N_BC_C <- N_BC - N_BC_B # number of patients under C in BC
  # simulate correlated (or uncorrelated) continuous covariates from a Gaussian copula
  # with normally-distributed marginals and the specified correlation structure  
  X_AC_A <- as.matrix(genCorGen(n = N_AC_A, nvar = no.chars, params1 = as.numeric(rep(meanX_AC, no.chars)), 
                                params2 = as.numeric(rep(sdX^2, no.chars)), dist = "normal", corstr = "cs", 
                                corMatrix = R, wide = TRUE)[,-1]) # patients under A in trial AC 
  #N is the number of patients simulated, no.chars is 4
  X_AC_C <- as.matrix(genCorGen(n = N_AC_C, nvar = no.chars, params1 = as.numeric(rep(meanX_AC, no.chars)), 
                                params2 = as.numeric(rep(sdX^2, no.chars)), dist = "normal", corstr = "cs", 
                                corMatrix = R, wide = TRUE)[,-1]) # patients under C in trial AC 
  X_BC_B <- as.matrix(genCorGen(n = N_BC_B, nvar = no.chars, params1 = as.numeric(rep(meanX_BC, no.chars)), 
                                params2 = as.numeric(rep(sdX^2, no.chars)), dist = "normal", corstr = "cs", 
                                corMatrix = R, wide = TRUE)[,-1]) # patients under B in trial BC 
  X_BC_C <- as.matrix(genCorGen(n = N_BC_C, nvar = no.chars, params1 = as.numeric(rep(meanX_BC, no.chars)), 
                                params2 = as.numeric(rep(sdX^2, no.chars)), dist = "normal", corstr = "cs", 
                                corMatrix = R, wide = TRUE)[,-1]) # patients under C in trial BC 
  # set log hazards for each patient
  betaX_AC_A <- rep(b_trt_A, N_AC_A) # intercepts and treatment effects
  betaX_AC_C <- rep(0,N_AC_C)
  betaX_BC_B <- rep(b_trt_B, N_BC_B)
  betaX_BC_C <- rep(0,N_BC_C)
  col.names <- NULL
  for (k in 1:no.chars) {
    col.names <- c(col.names, paste0('X', k))
    betaX_AC_A <- betaX_AC_A + b_X*X_AC_A[,k] # prognostic variable effects
    betaX_AC_C <- betaX_AC_C + b_X*X_AC_C[,k]
    betaX_BC_B <- betaX_BC_B + b_X*X_BC_B[,k]
    betaX_BC_C <- betaX_BC_C + b_X*X_BC_C[,k]
    if (k <= no.ems) {
      betaX_AC_A <- betaX_AC_A + b_EM_A*X_AC_A[,k] # effect modification
      betaX_BC_B <- betaX_BC_B + b_EM_B*X_BC_B[,k]
    }
  }
  X_AC = as.data.frame(rbind(X_AC_A, X_AC_C))
  colnames(X_AC) <- col.names
  betaX_AC = c(betaX_AC_A, betaX_AC_C)
  # Generate AC Weibull-Cox-distributed latent event times according to Bender et al. (2005)
  Tlat = -log(runif(n = N_AC))/(weib_inv_scale*exp(betaX_AC))^(1/weib_shape) 
  
  #-------------------------------------------------------------------------------
  #Replaced cens_rate with cens_rateA 
  C = rexp(n = N_AC, rate = cens_rateA) # AC censoring times
  # AC follow-up times and event indicators
  time = pmin(Tlat, C)
  status = as.numeric(Tlat <= C)
  trt <- c(rep("A", N_AC_A), rep("C", N_AC_C)) # treatment assignment
  IPD.AC2 <- as.data.frame(cbind(trt, X_AC, time, status))
  X_BC = as.data.frame(rbind(X_BC_B, X_BC_C))
  colnames(X_BC) <- col.names
  betaX_BC = c(betaX_BC_B, betaX_BC_C)
  
  # Generate BC Weibull-Cox-distributed latent event times
  Tlat = -log(runif(n = N_BC))/(weib_inv_scale*exp(betaX_BC))^(1/weib_shape) 
  #-------------------------------------------------------------------------------
  #Replaced cens_rate with censrateB
  C = rexp(n = N_BC, rate = cens_rateB) # BC censoring times
  # BC follow-up times and event indicators
  time = pmin(Tlat, C)
  status = as.numeric(Tlat <= C)
  trt <- c(rep("B", N_BC_B), rep("C", N_BC_C)) # treatment assignment
  IPD.BC2 <- as.data.frame(cbind(trt, X_BC, time, status))
  IPD.BC2$trt <- factor(IPD.BC2$trt, levels = c("C","B")) 
  # aggregate the data for the BC trial 
  ALD.BC2 <- as.data.frame(cbind(
    # Trial mean covariate stats and summary outcomes (log HR and variance) for BC
    summarise(IPD.BC2, mean.X1 = mean(X1), mean.X2 = mean(X2), mean.X3 = mean(X3), mean.X4 = mean(X4),
              sd.X1=sd(X1),sd.X2=sd(X2),sd.X3=sd(X3),sd.X4=sd(X4),#add standard deviation
              logHR_B = summary(coxph(Surv(time, status)~trt, data  =  IPD.BC2))$coef[1],
              var_logHR_B = vcov(coxph(Surv(time, status)~trt, data  =  IPD.BC2))[[1]],
              HR_B = summary(coxph(Surv(time, status)~trt, data  =  IPD.BC2))$coef[2])))    
  list(IPD.AC2, IPD.BC2, ALD.BC2, states.data2)
}

##for (i in 1:scenarios) {
# print(i) #removing print as it makes code slow
##IPD.AC <- IPD.BC <- ALD.BC <- vector(mode = "list", replicates)
##states.data <- vector(mode = "list", replicates + 1) # random-number-generator states
##for (j in 1:replicates) {
#-------------------------------------------------------------------------------
#Replaced b_trt with b_trtA and b_trtB in the function below
#Replaced cens_rate with cens_rateA
#Added cens_rateB
##gen.datasets <- gen.data(no.chars = no.chars, no.ems = no.ems, N_AC = pc$N_AC[i], N_BC = N_BC,
##   b_trt_A = b_trtA, b_trt_B = b_trtB, b_X = pc$b_X[i], b_EM_A = pc$b_EM[i],
##   b_EM_B = pc$b_EM[i], meanX_AC = pc$meanX_AC[i], meanX_BC = meanX_BC, sdX = sdX,
##   weib_shape = weib_shape, weib_inv_scale = weib_inv_scale, cens_rateA = cens_rateA,
##   cens_rateB  =  cens_rateB, corX = pc$corX[i], allocation = allocation)
##IPD.AC[[j]] <- gen.datasets[[1]]
## IPD.BC[[j]] <- gen.datasets[[2]]
## ALD.BC[[j]] <- gen.datasets[[3]]
##states.data[[j]] <- gen.datasets[[4]]
##}
##states.data[[replicates + 1]] <- .Random.seed # final random-number-generator state
##file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits = 2), 
##   "b_EM", round(pc$b_EM[i], digits = 2), 
##  "meanX_AC", pc$meanX_AC[i], "corX", pc$corX[i]) 
##save(IPD.AC, file = paste0("Data/",effects,simulation,"IPD_AC_", file.id, ".RData"))
##save(IPD.BC, file = paste0("Data/",effects,simulation,"IPD_BC_", file.id, ".RData"))
##save(ALD.BC, file = paste0("Data/",effects,simulation,"ALD_BC_", file.id, ".RData"))
##save(states.data, file = paste0("Data/states_", file.id, ".RData"))
##}  




for (i in 1:scenarios) {
  
  start_time <- Sys.time()  # Capture start time
  
  cat("Running scenario", i, "of", scenarios, "\n")
  
  IPD.AC2 <- IPD.BC2 <- ALD.BC2 <- vector(mode = "list", replicates)
  
  states.data2 <- vector(mode = "list", replicates + 1) # random-number-generator states
  
  for (j in 1:replicates) {
    
    
    #-------------------------------------------------------------------------------
    
    # Replaced b_trt with b_trtA and b_trtB in the function below
    
    # Replaced cens_rate with cens_rateA
    
    # Added cens_rateB
    
    gen.datasets2 <- gen.data2(no.chars = no.chars, no.ems = no.ems, N_AC = pc$N_AC[i], N_BC = N_BC,
                             
                             b_trt_A = b_trtA, b_trt_B = b_trtB, b_X = pc$b_X[i], b_EM_A = pc$b_EM[i],
                             
                             b_EM_B = pc$b_EM[i], meanX_AC = pc$meanX_AC[i], meanX_BC = meanX_BC, sdX = sdX,
                             
                             weib_shape = weib_shape, weib_inv_scale = weib_inv_scale, cens_rateA = cens_rateA,
                             
                             cens_rateB  =  cens_rateB, corX = pc$corX[i], allocation = allocation)
    
    IPD.AC2[[j]] <- gen.datasets2[[1]]
    
    IPD.BC2[[j]] <- gen.datasets2[[2]]
    
    ALD.BC2[[j]] <- gen.datasets2[[3]]
    
    states.data2[[j]] <- gen.datasets2[[4]]
    
  }
  
  states.data2[[replicates + 1]] <- .Random.seed # final random-number-generator state
  
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits = 2), 
                    
                    "b_EM", round(pc$b_EM[i], digits = 2), 
                    
                    "meanX_AC", pc$meanX_AC[i], "corX", pc$corX[i]) 
  
  save(IPD.AC2, file = paste0("Data_var/", "2IPD_AC_", file.id, ".RData"))
  
  save(IPD.BC2, file = paste0("Data_var/", "2IPD_BC_", file.id, ".RData"))
  
  save(ALD.BC2, file = paste0("Data_var/", "2ALD_BC_", file.id, ".RData"))
  
  save(states.data2, file = paste0("Data_var/2states_", file.id, ".RData"))
  
  end_time <- Sys.time()  # Capture end time
  
  cat("Completed scenario", i, "in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")
  
}















