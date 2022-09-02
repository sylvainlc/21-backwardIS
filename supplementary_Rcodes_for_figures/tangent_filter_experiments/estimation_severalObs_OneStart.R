rm(list = ls())
library(GrandParisPackage)
library(parallel)
source("simulationParameters.R")

# Estimation parameters ---------------------------------------------------

estTheta = T; estSig = F; updateOrders <- rep(T, n)
gradientSteps <- get_grad_steps(0.6, cst = 8)
allRes <- mclapply(1:length(dir("simulated_data/")), function(i){
  seed <- 100 + i
  set.seed(seed)
  simulated_POD <- read.table(paste0("simulated_data/simul_data_seed", seed, ".txt"),
                              sep = ";", header = T)
  observations <- simulated_POD[, "observations"]
  thetaStart <- trueTheta
  Res <- fastTangOR(observations, times,
                    thetaModel = thetaStart, sigma2 = trueSigma2, particleSize = 500,
                    updateOrders = updateOrders, gradientSteps = gradientSteps, 
                    all = T, estimateSigma2 = F, randomWalkParam = 1)
  Res
}, mc.cores = detectCores())

thetaEst <- sapply(allRes, function(x) x$Estimates[,1]) %% (2*pi)
write.table(thetaEst, file = "simulation_results/est_severalObs_oneStart.txt", col.names = F,
            row.names = F, sep = ";")
