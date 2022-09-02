rm(list = ls())
library(GrandParisPackage)
library(parallel)
source("simulationParameters.R")
seed <- 122
simulated_POD <- read.table(paste0("simulated_data/simul_data_seed", seed, ".txt"),
                            sep = ";", header = T)
observations <- simulated_POD[, "observations"]


# Estimation parameters ---------------------------------------------------

estTheta = T; estSig = F; updateOrders <- rep(T, n)
n_start_points <- 6
decrease_power <- seq(0.5, 1, length.out = n_start_points)
allRes <- mclapply(decrease_power, function(gradient_power){
  set.seed(floor(1000 * gradient_power))
  thetaStart <- trueTheta
  gradientSteps <- get_grad_steps(gradient_power, 8)
  Res <- fastTangOR(observations, times,
                    thetaModel = thetaStart, sigma2 = trueSigma2,
                    updateOrders = updateOrders, gradientSteps = gradientSteps, 
                    all = T, estimateSigma2 = F, randomWalkParam = 1)
}, mc.cores = detectCores() - 1)

thetaEst <- sapply(allRes, function(x) x$Estimates[,1]) %% (2*pi)
colnames(thetaEst) <- paste0("exp", round(decrease_power, 3))
write.table(thetaEst, 
            file = "simulation_results/est_oneObs_oneStart_severalGrads.txt", 
            col.names = F,
            row.names = F, sep = ";")
