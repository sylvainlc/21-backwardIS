rm(list = ls())
library(GrandParisPackage)
library(parallel)
source("simulationParameters.R")
n_traj <- 500
mclapply(1:n_traj, 
         function(i){
  seed <- 100 + i
  set.seed(seed)
  simulated_POD <- SINE_simulate(theta = trueTheta, sigma = trueSigma2, x0 = 0,
                                 times = times)
  write.table(simulated_POD, paste0("simulated_data/simul_data_seed", seed, ".txt"),
              col.names = T, row.names = F, sep = ";")
}, mc.cores= detectCores())
