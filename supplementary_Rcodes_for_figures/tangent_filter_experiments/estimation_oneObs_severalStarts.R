rm(list = ls())
library(GrandParisPackage)
library(parallel)
source("simulationParameters.R")
seed <- 122
set.seed(seed)
simulated_POD <- read.table(paste0("simulated_data/simul_data_seed", seed, ".txt"),
                            sep = ";", header = T)
observations <- simulated_POD[, "observations"]


# Estimation parameters ---------------------------------------------------



get_estimation_one_obs_several_start <- function(IS){
  estTheta = TRUE; estSig = FALSE; updateOrders <- rep(TRUE, n)
  gradientSteps <- get_grad_steps(0.6, cst = 8)
  n_start_points <- 50
  set.seed(seed)
  thetaStart <- runif(1, 0, 2 * pi)
  N_tilde <- ifelse(IS, 5, 2)
  allRes <- mclapply(1:n_start_points, function(seed){
    set.seed(seed)
    thetaStart <- runif(1, 0, 2 * pi)
    fastTangOR(observations, times,
               thetaModel = thetaStart, sigma2 = trueSigma2,
               updateOrders = updateOrders, gradientSteps = gradientSteps, 
               all = FALSE, estimateTheta = estTheta, estimateSigma2 = estSig, 
               randomWalkParam = 1, backwardSampleSize = 2, IS = IS)
  },
  mc.cores = detectCores() - 1)
  thetaEst <- sapply(allRes, function(x) x$Estimates[,1]) %% (2*pi)
  out_path <- paste0("simulation_results/est_oneObs_severalStarts_",
                     ifelse(IS, "IS", "AR"), ".txt")
  write.table(thetaEst, file = out_path, col.names = F,
              row.names = F, sep = ";")
  return(NULL)
}

get_estimation_one_obs_several_start(TRUE)
